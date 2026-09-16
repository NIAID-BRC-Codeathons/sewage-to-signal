#!/usr/bin/env bash
# Build sae.sif. Three routes, because clusters differ in what they allow.
#
#   container/build.sh              # pick the best available route
#   container/build.sh apptainer    # native, needs root or --fakeroot
#   container/build.sh docker       # build with Docker, convert to SIF
#   container/build.sh remote       # Sylabs remote builder (needs `apptainer remote login`)
#
# Must run from a Linux host for the native route. Apptainer cannot build on
# macOS; use the docker route there and convert on the cluster, or build on the
# cluster directly.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
SIF="$HERE/sae.sif"
TAG="wastewater-sae:latest"

route="${1:-auto}"
have() { command -v "$1" >/dev/null 2>&1; }

# The repo working tree holds ~28 GB of data and a 1.2 GB venv next to the
# code. Building from it directly would pull all of that into the image, so
# stage a clean export of tracked files only. This also pins the image to a
# commit rather than to whatever is lying around.
stage_context() {
  have git || { echo "error: git required to stage a clean build context" >&2; exit 127; }
  cd "$REPO"
  git rev-parse --is-inside-work-tree >/dev/null 2>&1 \
    || { echo "error: $REPO is not a git repo" >&2; exit 1; }
  SHA="$(git rev-parse --short HEAD)"
  DIRTY=""
  git diff --quiet HEAD -- sae requirements.txt setup.sh data container 2>/dev/null || DIRTY=" (dirty)"
  STAGE="$(mktemp -d)"
  # HEAD, not the working tree: the image should match a commit.
  git archive HEAD | tar -x -C "$STAGE"
  echo "[build] context: $STAGE  from $SHA$DIRTY  ($(du -sh "$STAGE" | cut -f1))"
  [ -n "$DIRTY" ] && echo "[build] WARNING: uncommitted changes are NOT in this image"
  trap 'rm -rf "$STAGE"' EXIT
}

if [ "$route" = auto ]; then
  if have apptainer || have singularity; then route=apptainer
  elif have docker; then route=docker
  else
    echo "error: need apptainer, singularity or docker on PATH." >&2; exit 127
  fi
  echo "[build] route: $route (auto-selected)"
fi

case "$route" in
  apptainer)
    AP=apptainer; have apptainer || AP=singularity
    stage_context
    cd "$STAGE"
    if [ "$(id -u)" = 0 ]; then
      "$AP" build "$SIF" container/sae.def
    else
      echo "[build] not root — trying --fakeroot"
      "$AP" build --fakeroot "$SIF" container/sae.def
    fi
    ;;
  docker)
    have docker || { echo "error: docker not found" >&2; exit 127; }
    stage_context
    cd "$STAGE"
    echo "[build] docker build (linux/amd64)"
    docker build --platform linux/amd64 --label "git.sha=$SHA" \
      -t "$TAG" -f container/Dockerfile .
    if have apptainer || have singularity; then
      AP=apptainer; have apptainer || AP=singularity
      echo "[build] converting to SIF"
      docker save "$TAG" -o "$HERE/sae-docker.tar"
      "$AP" build "$SIF" "docker-archive://$HERE/sae-docker.tar"
      rm -f "$HERE/sae-docker.tar"
    else
      echo "[build] no apptainer here — image built as docker tag '$TAG'."
      echo "        To finish on a host that has apptainer:"
      echo "          docker save $TAG -o sae-docker.tar"
      echo "          # copy sae-docker.tar to that host, then:"
      echo "          apptainer build sae.sif docker-archive://sae-docker.tar"
    fi
    ;;
  remote)
    AP=apptainer; have apptainer || AP=singularity
    stage_context
    cd "$STAGE"
    "$AP" build --remote "$SIF" container/sae.def
    ;;
  *) echo "unknown route: $route" >&2; exit 2 ;;
esac

if [ -f "$SIF" ]; then
  echo "[build] wrote $SIF ($(du -h "$SIF" | cut -f1))"
  echo "[build] self-check:"
  ${AP:-apptainer} test "$SIF" || echo "  (test reported problems — see above)"
fi
