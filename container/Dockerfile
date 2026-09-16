# Mirror of sae.def, for two cases sae.def cannot cover:
#   * building without root/fakeroot — build here, convert to SIF:
#       docker build --platform linux/amd64 -t wastewater-sae:latest -f container/Dockerfile .
#       docker save wastewater-sae:latest -o sae-docker.tar
#       apptainer build sae.sif docker-archive://sae-docker.tar
#   * local smoke-testing on a machine without Apptainer.
# Keep the two files in step: same base, same install order, same paths.
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive \
    MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/venv/bin:/opt/conda/bin:/usr/local/bin:/usr/bin:/bin \
    PYTHONNOUSERSITE=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    HF_HOME=/hf \
    UV_PYTHON_INSTALL_DIR=/opt/uv-python \
    SAE_HOME=/opt/sae \
    MPLCONFIGDIR=/tmp

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl bzip2 xz-utils procps gzip \
        coreutils gawk findutils grep sed tar \
    && rm -rf /var/lib/apt/lists/*

# Bioinformatics binaries: mmseqs2 is not in the Ubuntu archive, so bioconda
# provides one resolver for all three.
RUN mkdir -p /opt/conda \
    && curl -fL -o /usr/local/bin/micromamba \
         https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64 \
    && chmod +x /usr/local/bin/micromamba \
    && /usr/local/bin/micromamba install -y -p /opt/conda \
         -c conda-forge -c bioconda megahit=1.2.9 mmseqs2=18.8cc5c fastp=1.3.7 \
    && /usr/local/bin/micromamba clean -y --all

# uv provisions Python 3.12 itself rather than trusting the base image's python3.
COPY requirements.txt /opt/build/requirements.txt
RUN curl -LsSf https://astral.sh/uv/install.sh | env UV_INSTALL_DIR=/usr/local/bin sh \
    && uv python install 3.12 \
    && uv venv --python 3.12 /opt/venv \
    && uv pip install --python /opt/venv/bin/python -r /opt/build/requirements.txt \
    && rm -rf /root/.cache/uv

COPY sae      /opt/sae/sae
COPY setup.sh /opt/sae/setup.sh
COPY data/fetch_fastq.sh data/subsample_rnaseq.sh data/atlas_manifest.tsv /opt/sae/data/

RUN { echo "built: $(date -u +%Y-%m-%dT%H:%M:%SZ)"; \
      echo "python: $(/opt/venv/bin/python -V 2>&1)"; echo; \
echo "micromamba: $(/usr/local/bin/micromamba --version 2>&1)"; \
            echo "--- bioconda ---"; \
      /usr/local/bin/micromamba list -p /opt/conda | awk 'NF>=2 && $1!="Name" && $1 !~ /^[-\xe2]/ {printf "%s==%s\n", $1, $2}'; \
      echo; echo "--- python packages ---"; uv pip freeze --python /opt/venv/bin/python; \
    } > /opt/image-manifest.txt \
    && mkdir -p /work /data /hf && chmod 777 /work /data /hf

# Documentation only; run.sh publishes this to the host's loopback.
EXPOSE 8765

WORKDIR /work
ENTRYPOINT []
CMD ["/bin/bash"]
