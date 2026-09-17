#!/bin/bash
# ==============================================================================
# Metagenomic Profiling Infrastructure Unified Deployment Script
# Targets: Kraken2, Ganon2 (ganon), EsVirtu (esviritu), VAPER, Kaiju, pyHMMER
# Dependencies: Mamba, Singularity/Apptainer, Git, Wget, Tar, Gzip
# ==============================================================================

# Halt execution immediately if a pipeline command exits with a non-zero status.
set -e

# Enforce group-write permissions for all subsequent files and directories created in this script
umask 002

# ==============================================================================
# Phase 1: Directory Initialization & Unified Virtual Environment Configuration
# ==============================================================================
#BASE_DIR="${HOME}/metagenomic_infrastructure"

BASE_DIR="$(cd "$(dirname "$0")" && pwd)/metagenomic_infrastructure"
ENV_NAME="ww_dm_env"
ENV_DIR="${BASE_DIR}/envs"
ENV_TARGET="${ENV_DIR}/${ENV_NAME}"
CONTAINER_DIR="${BASE_DIR}/singularity_containers"
DB_DIR="${BASE_DIR}/databases"
PIPELINE_DIR="${BASE_DIR}/pipelines"

echo "Initializing high-performance directory structure at ${BASE_DIR}..."
mkdir -p "${CONTAINER_DIR}"
mkdir -p "${DB_DIR}/kraken2"
mkdir -p "${DB_DIR}/ganon"
mkdir -p "${DB_DIR}/esviritu"
mkdir -p "${DB_DIR}/kaiju"
mkdir -p "${DB_DIR}/pyhmmer"
mkdir -p "${PIPELINE_DIR}"
mkdir -p "${ENV_DIR}"

# Post-Creation Fix: Set the sticky bit for group ownership inheritance
# and ensure group write permissions are applied to the base directory.
chmod g+s "${BASE_DIR}"
chmod -R g+w "${BASE_DIR}"

echo "Creating unified Mamba virtual environment: ${ENV_NAME}..."
# Note: VAPER is a Nextflow pipeline, not a standalone compiled binary on bioconda.
# We install 'nextflow' into this unified environment to ensure VAPER can be executed.
# 'raptor' is included explicitly to support Ganon2's HIBF indexing.
# Wrapping the mamba creation command with a temporary umask ensures Conda environment files inherit correctly.
(umask 002 && mamba create -y -c conda-forge -c bioconda -p ${ENV_TARGET}\
    kraken2 \
    ganon \
    esviritu \
    kaiju \
    pyhmmer \
    nextflow \
    raptor \
    wget \
    tar \
    gzip)

echo "Unified environment [${ENV_NAME}] created successfully."
echo "Execute 'conda activate ${ENV_NAME}' to initialize."

# ==============================================================================
# Phase 2: Singularity (OCI) Container Acquisition
# ==============================================================================
echo "Pulling Singularity containers from Quay.io Biocontainers registry..."
cd "${CONTAINER_DIR}"

# Translating Docker OCI images to Singularity .sif formats.
# Note: The 'latest' tag is utilized here for acquisition demonstration. In rigid 
# production environments, specific alphanumeric version tags must be hardcoded.
#singularity pull --name kraken2.sif docker://quay.io/biocontainers/kraken2:latest
#singularity pull --name ganon.sif docker://quay.io/biocontainers/ganon:latest
#singularity pull --name esviritu.sif docker://quay.io/biocontainers/esviritu:latest
#singularity pull --name kaiju.sif docker://quay.io/biocontainers/kaiju:latest
#singularity pull --name pyhmmer.sif docker://quay.io/biocontainers/pyhmmer:latest

# MISSING CONTAINER COMMENT: 
# A standalone Singularity container for VAPER is NOT downloaded here. VAPER operates 
# as a Nextflow pipeline which programmatically manages its own container fetching 
# during execution using native ociAutoPull directives based on its internal configs.

# ==============================================================================
# Phase 3: Standard Database Procurement and Pipeline Cloning
# ==============================================================================
echo "Initiating standard database downloads. Ensure sufficient storage and network bandwidth."

# ---------------------------------------------------------
# Kraken2 Standard Database (AWS S3)
# ---------------------------------------------------------
echo "Downloading Kraken2 Standard Database..."
cd "${DB_DIR}/kraken2"
# Downloading the down-sampled 8GB standard DB to prevent 100GB+ RAM exhaustion on standard nodes.
wget -q -O k2_standard_8gb.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240112.tar.gz
tar -xzf k2_standard_8gb.tar.gz
rm k2_standard_8gb.tar.gz
echo "Kraken2 8GB standard database ready."

# ---------------------------------------------------------
# Ganon2 Standard Database Build
# ---------------------------------------------------------
echo "Configuring Ganon standard default database..."
cd "${DB_DIR}/ganon"
# MISSING STANDARD DATABASE TARBALL COMMENT:
# Ganon does not supply static tarballs for its standard databases. Instead, it relies 
# on dynamic API building via the 'ganon build' command to pull the latest RefSeq/GTDB.
# The command below is commented out to prevent script hanging during lengthy API calls.
# Execute this manually inside the activated conda environment.
echo "# EXECUTE MANUALLY: conda run -n ${ENV_NAME} ganon build --db-prefix ganon_default_db"

# ---------------------------------------------------------
# EsViritu Virus Pathogen Database (Zenodo)
# ---------------------------------------------------------
echo "Downloading EsViritu Virus Pathogen Database v3.2.4 from Zenodo..."
cd "${DB_DIR}/esviritu"
wget -q -O esviritu_db_v3.2.4.tar.gz https://zenodo.org/records/17716199/files/esviritu_db_v3.2.4.tar.gz
tar -xzf esviritu_db_v3.2.4.tar.gz
rm esviritu_db_v3.2.4.tar.gz
# Export the recommended environment variable for EsViritu.
echo "export ESVIRITU_DB=${DB_DIR}/esviritu/v3.2.4" >> "${BASE_DIR}/env_vars.sh"
echo "EsViritu standard database ready."

# ---------------------------------------------------------
# VAPER Pipeline Acquisition
# ---------------------------------------------------------
echo "Cloning VAPER Nextflow Pipeline..."
cd "${PIPELINE_DIR}"
git clone https://github.com/DOH-JDJ0303/vaper.git

# MISSING STANDARD DATABASE COMMENT:
# VAPER does not require an external standalone standard database download recipe. 
# The Nextflow repository cloned above comes stocked internally with its required 
# reference sequences for viral assembly workflows.

# ---------------------------------------------------------
# Kaiju Pre-built Database (AWS S3)
# ---------------------------------------------------------
echo "Downloading Kaiju pre-built index (refseq_ref subset - 54 GB RAM required)..."
cd "${DB_DIR}/kaiju"
# Opting for refseq_ref over the massive 219 GB 'nr' database to preserve memory safety.
wget -q -O kaiju_db_refseq_ref.tgz https://kaiju-idx.s3.eu-central-1.amazonaws.com/2024/kaiju_db_refseq_ref_2024-08-14.tgz
tar -xzf kaiju_db_refseq_ref.tgz
rm kaiju_db_refseq_ref.tgz
echo "Kaiju refseq_ref standard database ready."

# ---------------------------------------------------------
# pyHMMER Pfam Database (EBI FTP)
# ---------------------------------------------------------
echo "Downloading Pfam-A standard HMM database for pyHMMER..."
cd "${DB_DIR}/pyhmmer"
wget -q -O Pfam-A.hmm.gz ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
echo "pyHMMER Pfam database ready."

echo "=============================================================================="
echo "Infrastructure deployment complete. All environments, containers, and databases staged."
echo "=============================================================================="

