#!/usr/bin/env bash
# Run the pipeline inside the container image with the right bind mounts.
#
#   container/run.sh pipeline --fastq /data/fastq_rnaseq/x_1.fastq.gz \
#                             --fastq2 /data/fastq_rnaseq/x_2.fastq.gz --sample S1
#   container/run.sh query --fasta /data/contig.fna --top-k 8
#   container/run.sh manifest            # package versions baked into the image
#   container/run.sh test                # self-check
#   container/run.sh shell               # interactive
#
# Runtimes, in preference order: apptainer, singularity, docker. Override with
# SAE_RUNTIME=docker. Apptainer is the cluster path and takes a .sif; docker is
# the local path on machines without Apptainer (notably macOS) and takes a tag.
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
IMAGE="${SAE_IMAGE:-wastewater-sae:latest}"
DATA="${SAE_DATA:-$REPO/data}"
HF="${SAE_HF:-${HF_HOME:-$HOME/.cache/huggingface}}"
WORK="${SAE_WORK:-$PWD/work}"
ATLAS="${SAE_ATLAS:-$REPO/sae}"

# Entrypoints, kept identical to the %apprun sections in sae.def.
PY=/opt/venv/bin/python
APP_PIPELINE="$PY /opt/sae/sae/pipeline/run.py"
APP_QUERY="$PY /opt/sae/sae/sae_testing_script.py"

RUNNER="${SAE_RUNTIME:-}"
if [ -z "$RUNNER" ]; then
  for r in apptainer singularity docker; do
    command -v "$r" >/dev/null && { RUNNER="$r"; break; }
  done
fi
if [ -z "$RUNNER" ]; then
  echo "error: none of apptainer, singularity or docker found on PATH." >&2
  echo "  On a cluster this is usually:  module load apptainer" >&2
  exit 127
fi

mkdir -p "$WORK" "$HF"
cmd="${1:-help}"; shift || true

# --- docker -----------------------------------------------------------------
if [ "$RUNNER" = docker ]; then
  if ! docker image inspect "$IMAGE" >/dev/null 2>&1; then
    cat >&2 <<EOT
error: image not found: $IMAGE

  build it:        container/build.sh docker
  or point at one: SAE_IMAGE=myrepo/sae:tag container/run.sh ...

  note: build.sh stages the context from 'git archive HEAD', so uncommitted
  changes are not in the image. Shadow them with SAE_CODE=$REPO/sae.
EOT
    exit 1
  fi

  ARGS=(--rm --platform "${SAE_PLATFORM:-linux/amd64}"
        -v "$DATA:/data" -v "$HF:/hf" -v "$WORK:/work" -v "$ATLAS:/atlas")
  # Let a working copy shadow the baked-in code, for iterating without a rebuild.
  [ -n "${SAE_CODE:-}" ] && ARGS+=(-v "$SAE_CODE:/opt/sae/sae")
  # Root in the container would leave root-owned files in /work on Linux.
  # Docker Desktop already maps ownership, and --user breaks it there.
  [ "$(uname -s)" = Linux ] && ARGS+=(--user "$(id -u):$(id -g)")
  command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -L >/dev/null 2>&1 \
    && ARGS+=(--gpus all)
  [ -t 0 ] && [ -t 1 ] && ARGS+=(-it)

  case "$cmd" in
    pipeline) exec docker run "${ARGS[@]}" "$IMAGE" $APP_PIPELINE --work /work "$@" ;;
    query)    exec docker run "${ARGS[@]}" "$IMAGE" $APP_QUERY "$@" ;;
    manifest) exec docker run "${ARGS[@]}" "$IMAGE" cat /opt/image-manifest.txt ;;
    shell)    exec docker run "${ARGS[@]}" "$IMAGE" /bin/bash ;;
    test)
      # Mirrors the %test section of sae.def, which apptainer runs natively.
      exec docker run "${ARGS[@]}" "$IMAGE" /bin/bash -ec '
        echo "--- python deps ---"
        '"$PY"' -c "import torch, esm, pyarrow, pyrodigal, pyhmmer, pyfastx; \
print(\"torch\", torch.__version__); print(\"esm\", esm.__version__); \
print(\"pyrodigal\", pyrodigal.__version__)"
        echo "--- external binaries ---"
        for t in megahit mmseqs fastp; do
          command -v $t >/dev/null || { echo "MISSING $t"; exit 1; }
        done
        megahit --version; mmseqs version; fastp --version 2>&1 | head -1
        echo "--- pipeline imports ---"
        cd /opt/sae/sae/pipeline && '"$PY"' -c \
          "import common, s01_qc, s02_assemble, s03_genes, s04_derep, \
s05_prefilter, s06_embed, s07_match, run; print(\"all stages import\")"' ;;
    help|*)
      sed -n '2,22p' "$0" | sed 's/^# \{0,1\}//'
      echo "runtime: docker (image $IMAGE)" ;;
  esac
  exit $?
fi

# --- apptainer / singularity ------------------------------------------------
if [ ! -f "$SIF" ]; then
  cat >&2 <<EOT
error: image not found: $SIF

  build it:        container/build.sh
  or point at one: SAE_SIF=/shared/images/sae.sif container/run.sh ...
EOT
  exit 1
fi

# --nv passes the host NVIDIA driver through. Harmless to omit on CPU nodes,
# but it fails loudly if requested without a driver, so probe first.
NV=""
if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -L >/dev/null 2>&1; then
  NV="--nv"
fi

BINDS=(-B "$DATA:/data" -B "$HF:/hf" -B "$WORK:/work" -B "$ATLAS:/atlas")
[ -n "${SAE_CODE:-}" ] && BINDS+=(-B "$SAE_CODE:/opt/sae/sae")

case "$cmd" in
  pipeline) exec $RUNNER run $NV "${BINDS[@]}" --app pipeline "$SIF" --work /work "$@" ;;
  query)    exec $RUNNER run $NV "${BINDS[@]}" --app query    "$SIF" "$@" ;;
  manifest) exec $RUNNER run     "${BINDS[@]}" --app manifest "$SIF" ;;
  test)     exec $RUNNER test    "$SIF" ;;
  shell)    exec $RUNNER shell   $NV "${BINDS[@]}" "$SIF" ;;
  help|*)   exec $RUNNER run-help "$SIF" ;;
esac
