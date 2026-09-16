#!/usr/bin/env bash
# Build sae.sif from container/sae.def. Three routes, because hosts differ in
# what they allow.
#
#   container/build.sh              # pick the best available route
#   container/build.sh pull         # fetch the prebuilt image from GHCR
#   container/build.sh apptainer    # native, needs root or --fakeroot
#   container/build.sh nested       # Apptainer inside Docker, for macOS
#   container/build.sh remote       # Sylabs remote builder (needs `apptainer remote login`)
#
# On a cluster, prefer `pull`. Building needs root or --fakeroot and HPC sites
# commonly disable fakeroot, so pulling may be the only route that works there.
# CI publishes the image on changes to requirements.txt or sae.def.
#
# sae.def is the only recipe. Apptainer cannot build on macOS, so the nested
# route runs Apptainer in a Docker container and produces the same .sif the
# cluster runs - rather than a second, Docker-shaped image that would have to
# be kept in step and would never exercise Apptainer's own semantics.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
SIF="$HERE/sae.sif"
# Carries Apptainer only; it is not a second recipe for this project.
AP_IMAGE="${SAE_APPTAINER_IMAGE:-quay.io/singularity/singularity:v4.1.0}"
# Published by .github/workflows/container.yml. GHCR paths are lowercase.
ORAS="${SAE_ORAS:-oras://ghcr.io/niaid-brc-codeathons/sewage-to-signal/sae:latest}"

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
  elif have docker; then route=nested
  else
    echo "error: need apptainer, singularity or docker on PATH." >&2; exit 127
  fi
  echo "[build] route: $route (auto-selected)"
fi

case "$route" in
  pull)
    if have apptainer || have singularity; then
      AP=apptainer; have apptainer || AP=singularity
      echo "[build] pull $ORAS"
      "$AP" pull --force "$SIF" "$ORAS"
    elif have docker; then
      echo "[build] pull $ORAS (via apptainer in docker)"
      docker run --rm --privileged --platform "${SAE_PLATFORM:-linux/amd64}" \
        -v "$HERE:/out" "$AP_IMAGE" pull --force /out/sae.sif "$ORAS"
    else
      echo "error: need apptainer, singularity or docker on PATH." >&2; exit 127
    fi
    ;;
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
  nested)
    have docker || { echo "error: docker not found" >&2; exit 127; }
    stage_context
    # Building a SIF creates user namespaces, which inside Docker requires
    # --privileged. On Docker Desktop that is confined to the Linux VM.
    echo "[build] apptainer inside docker ($AP_IMAGE)"
    docker run --rm --privileged --platform "${SAE_PLATFORM:-linux/amd64}" \
      -v "$STAGE:/w" -w /w "$AP_IMAGE" build /w/sae.sif container/sae.def
    mv "$STAGE/sae.sif" "$SIF"
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
  if have apptainer || have singularity; then
    ${AP:-apptainer} test "$SIF" || echo "  (test reported problems — see above)"
  else
    docker run --rm --privileged --platform "${SAE_PLATFORM:-linux/amd64}" \
      -v "$HERE:/mnt/sif" "$AP_IMAGE" test "/mnt/sif/$(basename "$SIF")" \
      || echo "  (test reported problems — see above)"
  fi
fi
