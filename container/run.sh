#!/usr/bin/env bash
# Run the pipeline inside the Apptainer image with the right bind mounts.
#
#   container/run.sh pipeline --fastq /data/fastq_rnaseq/x_1.fastq.gz \
#                             --fastq2 /data/fastq_rnaseq/x_2.fastq.gz --sample S1
#   container/run.sh query --fasta /data/contig.fna --top-k 8
#   container/run.sh manifest            # package versions baked into the image
#   container/run.sh shell               # interactive
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

RUNNER=""
command -v apptainer >/dev/null && RUNNER=apptainer
[ -z "$RUNNER" ] && command -v singularity >/dev/null && RUNNER=singularity
if [ -z "$RUNNER" ]; then
  echo "error: neither apptainer nor singularity found on PATH." >&2
  echo "  On a cluster this is usually:  module load apptainer" >&2
  exit 127
fi

if [ ! -f "$SIF" ]; then
  cat >&2 <<EOT
error: image not found: $SIF

  build it:        container/build.sh
  or point at one: SAE_SIF=/shared/images/sae.sif container/run.sh ...
EOT
  exit 1
fi

mkdir -p "$WORK" "$HF"

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

cmd="${1:-help}"; shift || true
case "$cmd" in
  pipeline) exec $RUNNER run $NV "${BINDS[@]}" --app pipeline "$SIF" --work /work "$@" ;;
  query)    exec $RUNNER run $NV "${BINDS[@]}" --app query    "$SIF" "$@" ;;
  manifest) exec $RUNNER run     "${BINDS[@]}" --app manifest "$SIF" ;;
  test)     exec $RUNNER test    "$SIF" ;;
  shell)    exec $RUNNER shell   $NV "${BINDS[@]}" "$SIF" ;;
  help|*)   exec $RUNNER run-help "$SIF" ;;
esac
