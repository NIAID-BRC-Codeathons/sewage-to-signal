# The local SLURM node, with Apptainer on it.
#
# This is a TEST HARNESS, not a second recipe for the project. It carries no
# part of the pipeline: it exists so a compute node here looks like a compute
# node on a cluster — SLURM plus Apptainer, running the same sae.sif. Without
# Apptainer the node can schedule a job but not execute one, and the local
# setup would only ever test half the path.
ARG SLURM_IMAGE=nathanhess/slurm:full
FROM ${SLURM_IMAGE}

ARG APPTAINER_VERSION=1.5.3
USER root
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates curl \
 && curl -fsSLO "https://github.com/apptainer/apptainer/releases/download/v${APPTAINER_VERSION}/apptainer_${APPTAINER_VERSION}_amd64.deb" \
 && apt-get install -y "./apptainer_${APPTAINER_VERSION}_amd64.deb" \
 && rm -f "apptainer_${APPTAINER_VERSION}_amd64.deb" \
 && rm -rf /var/lib/apt/lists/* \
 && apptainer --version
