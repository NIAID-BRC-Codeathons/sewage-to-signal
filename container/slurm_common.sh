# Shared prelude for the .sbatch scripts. Sourced, never run.
#
# The same file has to work two ways: submitted with `sbatch` on a cluster, and
# run with `bash` on a laptop or a plain VM. Everything SLURM would provide is
# therefore defaulted rather than assumed.

# SLURM copies the batch script to a spool directory before running it, so $0
# is not inside the repo there. The submit directory is, and SLURM records it.
find_repo() {
  local c
  for c in "${REPO:-}" "${SLURM_SUBMIT_DIR:-}" \
           "$(cd "$(dirname "${BASH_SOURCE[1]:-$0}")/.." 2>/dev/null && pwd)" "$PWD"; do
    [ -n "$c" ] && [ -x "$c/container/run.sh" ] && { printf '%s\n' "$c"; return; }
  done
  echo "error: cannot locate the repository." >&2
  echo "  run this from a checkout, or set REPO=/path/to/sewage-to-signal" >&2
  exit 2
}

sae_env() {
  REPO="$(find_repo)" || exit 2
  export REPO

  if [ -n "${SLURM_JOB_ID:-}" ]; then
    SAE_MODE=slurm
    CPUS="${SLURM_CPUS_PER_TASK:-1}"
    JOB="$SLURM_JOB_ID"
  else
    SAE_MODE=local
    CPUS="$( (command -v nproc >/dev/null 2>&1 && nproc) \
             || sysctl -n hw.ncpu 2>/dev/null || echo 4 )"
    JOB="local-$$"
  fi
  export SAE_MODE CPUS JOB OMP_NUM_THREADS="$CPUS"

  # Cluster layout when it exists, repo-local otherwise. Explicit env wins.
  local scratch="/scratch/${USER:-nobody}"
  [ -d "$scratch" ] || scratch=""
  export SAE_HF="${SAE_HF:-${scratch:+$scratch/hf_cache}}"
  export SAE_HF="${SAE_HF:-${HF_HOME:-$HOME/.cache/huggingface}}"
  export SAE_WORK="${SAE_WORK:-${scratch:+$scratch/sae_work}}"
  export SAE_WORK="${SAE_WORK:-$REPO/work}"
  export SAE_DATA="${SAE_DATA:-$REPO/data}"
  if [ -z "${SAE_SIF:-}" ]; then
    if [ -f /shared/images/sae.sif ]; then export SAE_SIF=/shared/images/sae.sif
    else export SAE_SIF="$REPO/container/sae.sif"; fi
  fi

  mkdir -p "$SAE_WORK" "$SAE_HF" 2>/dev/null || true
}

# DRY_RUN=1 prints the command instead of running it — cheap way to check paths
# before spending an allocation, and how these scripts are tested off-cluster.
sae_run() {
  if [ -n "${DRY_RUN:-}" ]; then
    printf '  would run:'; printf ' %q' "$@"; printf '\n'
    return 0
  fi
  "$@"
}

sae_banner() {
  cat <<BANNER
  mode   : $SAE_MODE (job $JOB, $CPUS cpus)
  repo   : $REPO
  image  : $SAE_SIF$( [ -f "$SAE_SIF" ] || printf ' (MISSING — container/build.sh pull)' )
  work   : $SAE_WORK
  data   : $SAE_DATA
  hf     : $SAE_HF
BANNER
}
