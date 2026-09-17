#!/usr/bin/env bash
# Run one pipeline per sample, fanned out across the GPUs of a single host.
#
#   container/run_batch.sh SRR1 SRR2 SRR3 ...
#   container/run_batch.sh --samples samples.txt
#   cat samples.txt | container/run_batch.sh
#   container/run_batch.sh S1 S2 -- --model 6b --batch-size 64
#   DRY_RUN=1 container/run_batch.sh S1 S2          # print the commands only
#
# This is the no-scheduler counterpart to container/slurm_example.sbatch. On a
# cluster SLURM decides which GPU a job gets; on one big server nothing does,
# so eight people — or eight samples — all land on GPU 0 and OOM each other.
# Pinning CUDA_VISIBLE_DEVICES per slot is the whole job of this script.
#
# Slots are refilled as they free, not in fixed rounds. Samples differ in depth
# by an order of magnitude, so a barrier every N would leave most GPUs idle
# waiting on the deepest one in each round.
#
# Everything after `--` is passed through to `run.sh pipeline`, so this script
# never has to learn the pipeline's options.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"

# Paths are the image's, not the host's: run.sh binds SAE_DATA onto /data.
FQ_DIR="${SAE_FASTQ_DIR:-/data/fastq_rnaseq}"
FQ1="${SAE_FASTQ_SUFFIX1:-_1.fastq.gz}"
FQ2="${SAE_FASTQ_SUFFIX2:-_2.fastq.gz}"
LOGS="${SAE_LOGS:-$REPO/logs}"

SAMPLES=()
PASSTHRU=()
while [ $# -gt 0 ]; do
  case "$1" in
    --samples) shift; [ -f "${1:-}" ] || { echo "error: no such file: ${1:-}" >&2; exit 2; }
               while read -r s; do
                 s="${s%%#*}"; s="$(echo "$s" | tr -d '[:space:]')"
                 [ -n "$s" ] && SAMPLES+=("$s")
               done < "$1"; shift ;;
    --)        shift; PASSTHRU=("$@"); break ;;
    -h|--help) sed -n '2,19p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *)         SAMPLES+=("$1"); shift ;;
  esac
done
# No arguments and a pipe: read the list from stdin.
if [ ${#SAMPLES[@]} -eq 0 ] && [ ! -t 0 ]; then
  while read -r s; do
    s="${s%%#*}"; s="$(echo "$s" | tr -d '[:space:]')"
    [ -n "$s" ] && SAMPLES+=("$s")
  done
fi
[ ${#SAMPLES[@]} -gt 0 ] || { echo "usage: $0 SAMPLE... [-- pipeline options]" >&2; exit 2; }

# Count GPUs without depending on nvidia-smi, which is packaged separately from
# the driver and is missing on plenty of otherwise working hosts.
count_gpus() {
  local n
  n=$(ls -d /dev/nvidia[0-9]* 2>/dev/null | wc -l | tr -d ' ')
  [ "${n:-0}" -gt 0 ] && { echo "$n"; return; }
  if command -v nvidia-smi >/dev/null 2>&1; then
    n=$(nvidia-smi -L 2>/dev/null | wc -l | tr -d ' ')
    [ "${n:-0}" -gt 0 ] && { echo "$n"; return; }
  fi
  echo 0
}
GPUS="${SAE_GPUS:-$(count_gpus)}"
# With no GPU there is still work to do — the CPU stages — but running eight
# assemblies at once on one box is a memory problem, not a speedup.
if [ "$GPUS" -eq 0 ]; then
  GPUS="${SAE_JOBS:-1}"
  PIN=0
  echo "  no GPUs found — running $GPUS at a time on CPU (SAE_JOBS to change)"
else
  PIN=1
fi

mkdir -p "$LOGS"
cat <<BANNER
  samples : ${#SAMPLES[@]}
  slots   : $GPUS$( [ "$PIN" = 1 ] && printf ' (one GPU each)' )
  logs    : $LOGS
  work    : ${SAE_WORK:-$PWD/work}
  extra   : ${PASSTHRU[*]:-(none)}
BANNER

declare -a SLOT_PID SLOT_SAMPLE
FAILED=()

# Launch one sample on a given slot. run.py namespaces its own output by
# sample (<work>/<sample>/<stage>), so every run shares one work root and the
# dashboard sees them all under it.
launch() {
  local gpu="$1" s="$2" log="$LOGS/$2.log"
  local -a cmd=("$HERE/run.sh" pipeline
                --fastq  "$FQ_DIR/${s}${FQ1}"
                --fastq2 "$FQ_DIR/${s}${FQ2}"
                --sample "$s" ${PASSTHRU[@]+"${PASSTHRU[@]}"})
  if [ -n "${DRY_RUN:-}" ]; then
    printf '  would run [gpu %s]:' "$gpu"; printf ' %q' "${cmd[@]}"; printf '\n'
    return 0
  fi
  echo "  [gpu $gpu] $s -> $log"
  if [ "$PIN" = 1 ]; then
    # Apptainer passes the host environment through by default, but the
    # APPTAINERENV_ form is explicit and survives --cleanenv if it is ever added.
    CUDA_VISIBLE_DEVICES="$gpu" APPTAINERENV_CUDA_VISIBLE_DEVICES="$gpu" \
      "${cmd[@]}" >"$log" 2>&1 &
  else
    "${cmd[@]}" >"$log" 2>&1 &
  fi
  SLOT_PID[$gpu]=$!
  SLOT_SAMPLE[$gpu]="$s"
}

# Block until some slot is free, reaping whatever finished there. Returns the
# slot index in $FREE.
FREE=""
wait_for_slot() {
  while :; do
    local g pid
    for g in $(seq 0 $((GPUS - 1))); do
      pid="${SLOT_PID[$g]:-}"
      if [ -z "$pid" ]; then FREE="$g"; return; fi
      if ! kill -0 "$pid" 2>/dev/null; then
        wait "$pid"; local rc=$?
        [ "$rc" -eq 0 ] || FAILED+=("${SLOT_SAMPLE[$g]} (exit $rc)")
        SLOT_PID[$g]=""; FREE="$g"; return
      fi
    done
    sleep 2
  done
}

n=0
for s in "${SAMPLES[@]}"; do
  if [ -n "${DRY_RUN:-}" ]; then
    launch $(( n % GPUS )) "$s"; n=$(( n + 1 ))
    continue
  fi
  wait_for_slot
  launch "$FREE" "$s"
done

# Drain.
for g in $(seq 0 $((GPUS - 1))); do
  pid="${SLOT_PID[$g]:-}"
  [ -n "$pid" ] || continue
  wait "$pid"; rc=$?
  [ "$rc" -eq 0 ] || FAILED+=("${SLOT_SAMPLE[$g]} (exit $rc)")
done

[ -n "${DRY_RUN:-}" ] && exit 0
if [ ${#FAILED[@]-0} -gt 0 ]; then
  echo
  echo "  ${#FAILED[@]} of ${#SAMPLES[@]} failed:"
  printf '    %s\n' ${FAILED[@]+"${FAILED[@]}"}
  echo "  logs in $LOGS"
  exit 1
fi
echo
echo "  all ${#SAMPLES[@]} samples completed"
