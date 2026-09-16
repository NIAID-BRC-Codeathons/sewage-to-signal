#!/bin/bash
#SBATCH --job-name=sae-web
#SBATCH --output=logs/sae-web-%j.out
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
# Adjust partition/account for your cluster:
##SBATCH --partition=interactive
##SBATCH --account=your_account
#
# Run the progress UI on a compute node and reach it from your laptop.
#
#     sbatch container/slurm_web.sbatch
#     tail -f logs/sae-web-<jobid>.out      # prints the ssh command to use
#
# Compute nodes are not reachable from outside the cluster, so the server binds
# loopback on the node and you tunnel to it. Binding loopback is not just
# convention here: this server accepts uploads and launches subprocesses, so on
# a shared cluster a wider bind would hand every other user on the network a
# way to run commands as you. Do not widen it. If you need a dashboard others
# can see, run it --read-only and still tunnel.
#
# Sizing: these defaults suit monitoring and short stages. The server launches
# pipeline runs as child processes *inside this allocation*, so if you intend
# to start real work from the UI, raise --cpus-per-task/--mem/--time and add
# --gres=gpu:1 for s06_embed. For heavy runs prefer submitting
# slurm_example.sbatch per sample and using this UI read-only to watch them.
set -euo pipefail
module load apptainer 2>/dev/null || module load singularity 2>/dev/null || true

REPO=${REPO:-$HOME/sewage-to-signal}

# Shared scratch, not $HOME: the model cache alone is ~27 GB and every job
# reuses it.
export SAE_HF=${SAE_HF:-/scratch/$USER/hf_cache}
export SAE_WORK=${SAE_WORK:-/scratch/$USER/sae_work}
export SAE_DATA=${SAE_DATA:-$REPO/data}
export SAE_SIF=${SAE_SIF:-/shared/images/sae.sif}
mkdir -p logs "$SAE_WORK" "$SAE_HF"

# A shared node may already have something on a fixed port, and several people
# may run this at once, so take a free one rather than assuming 8765.
pick_port() {
  local p
  for _ in $(seq 25); do
    # $RANDOM is a bash builtin; shuf is GNU-only and not everywhere.
    p=$(( 20000 + RANDOM % 40000 ))
    # A refused connection means nothing is listening, so the port is free.
    if ! (exec 3<>"/dev/tcp/127.0.0.1/$p") 2>/dev/null; then
      echo "$p"; return
    fi
  done
  echo 8765
}
export SAE_PORT=${SAE_PORT:-$(pick_port)}

cat <<BANNER

  sae web UI
    node : $(hostname)
    port : $SAE_PORT
    work : $SAE_WORK

  From your laptop, open a tunnel that terminates on the compute node:

    ssh -N -L ${SAE_PORT}:127.0.0.1:${SAE_PORT} -J ${USER}@LOGIN_HOST ${USER}@$(hostname)

  then browse http://127.0.0.1:${SAE_PORT}

  If your cluster refuses ssh straight to compute nodes, forward through the
  login node instead and read the note above about binding before you do:

    ssh -N -L ${SAE_PORT}:$(hostname):${SAE_PORT} ${USER}@LOGIN_HOST

BANNER

# run.sh binds loopback for native Apptainer, which is what the tunnel expects.
exec "$REPO"/container/run.sh web --work /work
