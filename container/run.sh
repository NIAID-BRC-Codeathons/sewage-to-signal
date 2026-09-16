#!/usr/bin/env bash
# Run the pipeline from the Apptainer image, with the right bind mounts.
#
#   container/run.sh pipeline --fastq /data/fastq_rnaseq/x_1.fastq.gz \
#                             --fastq2 /data/fastq_rnaseq/x_2.fastq.gz --sample S1
#   container/run.sh setup status        # what data is present on this node
#   container/run.sh setup model-small   # fetch into the /hf and /data binds
#   container/run.sh web                 # progress UI on http://127.0.0.1:8765
#   container/run.sh query --fasta /data/contig.fna --top-k 8
#   container/run.sh manifest            # package versions baked into the image
#   container/run.sh test                # self-check
#   container/run.sh shell               # interactive
#
# There is one image artifact, sae.sif, and one recipe, sae.def. Where
# Apptainer is installed it runs natively; where it is not but Docker is (a
# macOS laptop), the same .sif runs under Apptainer *inside* a Docker
# container. Nothing runs as a plain Docker image, so the deployed runtime is
# the only runtime ever exercised. SAE_RUNTIME forces apptainer|singularity|nested.
#
# Paths inside the container:
#   /data  <- $SAE_DATA   (default: <repo>/data)
#   /hf    <- $SAE_HF     (default: ${HF_HOME:-$HOME/.cache/huggingface})
#   /work  <- $SAE_WORK   (default: $PWD/work)
#   /atlas <- $SAE_ATLAS  (default: <repo>/sae, holds the parquet tables)
#
# The image is read-only, so every writable path is a bind. The HF cache is
# bound rather than baked because it is ~27 GB and shared between runs.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"

SIF="${SAE_SIF:-$HERE/sae.sif}"
DATA="${SAE_DATA:-$REPO/data}"
HF="${SAE_HF:-${HF_HOME:-$HOME/.cache/huggingface}}"
WORK="${SAE_WORK:-$PWD/work}"
ATLAS="${SAE_ATLAS:-$REPO/sae}"
PORT="${SAE_PORT:-8765}"
# Carries Apptainer only; it is not a second recipe for this project.
AP_IMAGE="${SAE_APPTAINER_IMAGE:-quay.io/singularity/singularity:v4.1.0}"
# Provisioning runs inside the image but writes to the binds: the image is
# read-only, so models land in /hf and data in /data. Not a %apprun, because
# it takes setup.sh's own subcommands.
APP_SETUP="/opt/sae/setup.sh"

RUNNER="${SAE_RUNTIME:-}"
if [ -z "$RUNNER" ]; then
  for r in apptainer singularity; do
    command -v "$r" >/dev/null && { RUNNER="$r"; break; }
  done
  [ -z "$RUNNER" ] && command -v docker >/dev/null && RUNNER=nested
fi
if [ -z "$RUNNER" ]; then
  echo "error: need apptainer or singularity, or docker to nest them in." >&2
  echo "  On a cluster this is usually:  module load apptainer" >&2
  exit 127
fi

if [ ! -f "$SIF" ]; then
  cat >&2 <<EOT
error: image not found: $SIF

  fetch it:        container/build.sh pull     (prebuilt, from GHCR)
  or build it:     container/build.sh
  or point at one: SAE_SIF=/shared/images/sae.sif container/run.sh ...
EOT
  exit 1
fi

mkdir -p "$WORK" "$HF"
cmd="${1:-help}"; shift || true

# --nv passes the host NVIDIA driver through. Harmless to omit on CPU nodes,
# but it fails loudly if requested without a driver, so probe first.
NV=""
if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -L >/dev/null 2>&1; then
  NV="--nv"
fi

BINDS=(-B "$DATA:/data" -B "$HF:/hf" -B "$WORK:/work" -B "$ATLAS:/atlas")
# Let a working copy of the code shadow the baked-in one, for iterating without
# a rebuild. Unset by default so a plain run uses exactly what was built.
[ -n "${SAE_CODE:-}" ] && BINDS+=(-B "$SAE_CODE:/opt/sae/sae")

# --- Apptainer nested in Docker --------------------------------------------
# Two layers of bind translation: Docker puts the host paths under /mnt, then
# Apptainer binds those to the paths the image expects. Building and running a
# SIF both need user namespaces, which inside Docker means --privileged.
if [ "$RUNNER" = nested ]; then
  SIF_DIR="$(cd "$(dirname "$SIF")" && pwd)"; SIF_NAME="$(basename "$SIF")"
  D=(--rm --privileged --platform "${SAE_PLATFORM:-linux/amd64}"
     -v "$SIF_DIR:/mnt/sif" -v "$DATA:/mnt/data" -v "$HF:/mnt/hf"
     -v "$WORK:/mnt/work" -v "$ATLAS:/mnt/atlas")
  B=(-B /mnt/data:/data -B /mnt/hf:/hf -B /mnt/work:/work -B /mnt/atlas:/atlas)
  if [ -n "${SAE_CODE:-}" ]; then
    D+=(-v "$SAE_CODE:/mnt/code"); B+=(-B /mnt/code:/opt/sae/sae)
  fi
  [ -n "$NV" ] && D+=(--gpus all)
  [ -t 0 ] && [ -t 1 ] && D+=(-it)
  IMG="/mnt/sif/$SIF_NAME"

  case "$cmd" in
    pipeline) exec docker run "${D[@]}" "$AP_IMAGE" \
                run $NV "${B[@]}" --app pipeline "$IMG" --work /work "$@" ;;
    query)    exec docker run "${D[@]}" "$AP_IMAGE" \
                run $NV "${B[@]}" --app query "$IMG" "$@" ;;
    web)
      # Apptainer shares the Docker container's network namespace, so binding
      # localhost there would be unreachable through -p. Bind all interfaces
      # inside and publish to the host's loopback only.
      echo "  http://127.0.0.1:$PORT" >&2
      exec docker run "${D[@]}" -p "127.0.0.1:$PORT:8765" "$AP_IMAGE" \
        run "${B[@]}" --app web "$IMG" --host 0.0.0.0 --port 8765 --published "$@" ;;
    setup)    exec docker run "${D[@]}" "$AP_IMAGE" \
                exec "${B[@]}" "$IMG" "$APP_SETUP" "$@" ;;
    manifest) exec docker run "${D[@]}" "$AP_IMAGE" run "${B[@]}" --app manifest "$IMG" ;;
    test)     exec docker run "${D[@]}" "$AP_IMAGE" test "$IMG" ;;
    shell)    exec docker run "${D[@]}" "$AP_IMAGE" shell $NV "${B[@]}" "$IMG" ;;
    help|*)
      sed -n '2,27p' "$0" | sed 's/^# \{0,1\}//'
      echo "runtime: apptainer nested in docker ($AP_IMAGE)" ;;
  esac
  exit $?
fi

# --- native apptainer / singularity -----------------------------------------
case "$cmd" in
  pipeline) exec $RUNNER run $NV "${BINDS[@]}" --app pipeline "$SIF" --work /work "$@" ;;
  query)    exec $RUNNER run $NV "${BINDS[@]}" --app query    "$SIF" "$@" ;;
  web)      exec $RUNNER run     "${BINDS[@]}" --app web      "$SIF" --port "$PORT" "$@" ;;
  setup)    exec $RUNNER exec    "${BINDS[@]}" "$SIF" "$APP_SETUP" "$@" ;;
  manifest) exec $RUNNER run     "${BINDS[@]}" --app manifest "$SIF" ;;
  test)     exec $RUNNER test    "$SIF" ;;
  shell)    exec $RUNNER shell   $NV "${BINDS[@]}" "$SIF" ;;
  help|*)   exec $RUNNER run-help "$SIF" ;;
esac
