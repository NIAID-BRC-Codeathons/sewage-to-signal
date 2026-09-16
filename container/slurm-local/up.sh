#!/usr/bin/env bash
# A single-node SLURM in Docker, so the scheduler path can be exercised on a
# machine that has no cluster.
#
#   container/slurm-local/up.sh            # start it
#   container/slurm-local/up.sh down       # stop it
#   eval "$(container/slurm-local/up.sh env)"   # put the shims on PATH
#
# The point is that the deployment stops changing shape between a laptop and a
# cluster: with sbatch on PATH the web server submits instead of forking, which
# is what it will do on the cluster. The bin/ shims proxy sbatch, squeue, sacct
# and scancel into this container, so nothing above them knows the difference.
#
# The repo is mounted at its own absolute path, so a path that resolves on the
# host resolves identically inside a job. The Docker socket is mounted too, so
# a job can run container/run.sh and reach the same sae.sif the cluster uses.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
NAME="${SAE_SLURM_NAME:-sae-slurm}"
IMAGE="${SAE_SLURM_IMAGE:-sae-slurm-node}"

# The node needs Apptainer as well as SLURM, so that a job here executes the
# same sae.sif a cluster node would. Built locally; see Dockerfile.
if ! docker image inspect "$IMAGE" >/dev/null 2>&1 && [ "$IMAGE" = sae-slurm-node ]; then
  echo "[slurm-local] building the node image (SLURM + Apptainer)"
  docker build --quiet --platform "${SAE_PLATFORM:-linux/amd64}" \
    -t sae-slurm-node "$HERE" || exit 1
fi

case "${1:-up}" in
  env)
    printf 'export PATH=%q:$PATH\n' "$HERE/bin"
    exit 0 ;;
  down)
    docker rm -f "$NAME" >/dev/null 2>&1 && echo "stopped $NAME" || echo "not running"
    exit 0 ;;
  status)
    docker exec "$NAME" sinfo 2>/dev/null || { echo "not running"; exit 1; }
    exit 0 ;;
esac

command -v docker >/dev/null || { echo "error: docker not found" >&2; exit 127; }
if docker inspect "$NAME" >/dev/null 2>&1; then
  echo "already running: $NAME"
else
  echo "[slurm-local] starting $NAME from $IMAGE"
  docker run -d --name "$NAME" --hostname slurmctl --privileged \
    --platform "${SAE_PLATFORM:-linux/amd64}" \
    -v "$REPO:$REPO" \
    -v "${HF_HOME:-$HOME/.cache/huggingface}:${HF_HOME:-$HOME/.cache/huggingface}" \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -w "$REPO" \
    "$IMAGE" /bin/sh -c 'sudo /etc/startup.sh; sleep infinity' >/dev/null || exit 1
fi

echo -n "[slurm-local] waiting for the controller"
for _ in $(seq 60); do
  if docker exec "$NAME" sinfo >/dev/null 2>&1; then echo " ok"; break; fi
  echo -n "."; sleep 2
done
docker exec "$NAME" sinfo || { echo "controller did not come up" >&2; exit 1; }
cat <<EOT

  Put the shims on PATH, then start the server — it will submit rather than fork:

    eval "\$($HERE/up.sh env)"
    ./sae/.venv/bin/python sae/web/server.py

EOT
