Here is a `README.md` you can place in your `src` directory to explain the setup script to your group members.

# Metagenomic Profiling Infrastructure Setup

This directory contains the automated deployment script (`tools_setup.sh`) for provisioning our group's shared metagenomic analysis infrastructure.

The script is designed to run on a multi-user, high-performance computing (HPC) environment. It automatically configures group-write permissions, sets up a unified virtual environment, downloads isolated software containers, and provisions the massive reference databases required for sequence classification.

## Supported Tools

* **Kraken2:** Nucleotide taxonomic classification
* **Ganon2:** Highly compressed metagenomic profiling via Hierarchical Interleaved Bloom Filters
* **Kaiju:** Protein-level translated search for highly divergent sequences
* **EsViritu:** Specialized viral pathogen read mapping
* **pyHMMER:** Cython-optimized profile Hidden Markov Model searches
* **VAPER:** Nextflow-orchestrated viral metagenome assembly pipeline

## What the Script Does

1. **Enforces Group Permissions:** Uses `umask 002` and sets the directory sticky bit (`chmod g+s`) so all created folders, environments, and databases inherit group-write access automatically. This ensures no one in the group gets locked out of the shared tools.
2. **Builds a Unified Mamba Environment:** Installs the core binaries and dependencies for all tools into a single shared environment using the `-p` (prefix) path rather than a user-specific name.
3. **Pulls Singularity Containers:** Downloads OCI-compliant Docker images from Quay.io Biocontainers and converts them to `.sif` Singularity formats for reproducible, immutable pipeline execution.
4. **Provisions Standard Databases:**
* Downloads the 8GB memory-safe Kraken2 standard database.
* Downloads the Kaiju `refseq_ref` database (requires 54 GB of RAM to execute).
* Downloads the EsViritu v3.2.4 database tarball from Zenodo.
* Downloads the Pfam-A standard HMM database for pyHMMER.


5. **Clones Pipelines:** Fetches the VAPER Nextflow pipeline repository directly from GitHub.

## Usage Instructions

### 1. Run the Setup Script

Ensure you have `mamba` and `singularity` available in your module path, then execute:bash
bash tools_setup.sh

```

### 2. Activate the Environment
Because this is a shared group environment, you must activate it using the full absolute directory path rather than a Conda name:
```bash
conda activate /path/to/your/basedir/envs/metagenomics_unified_env

```

### 3. Post-Installation Manual Steps

Due to the architectural constraints of certain tools, a few steps must be executed manually after the script finishes:

* **Build the Ganon Database:** Ganon dynamically builds its database via an API rather than downloading a static file. With the environment activated, run:
```bash
ganon build --db-prefix /path/to/your/basedir/databases/ganon/ganon_default_db

```


* **Export EsViritu Path:** Before running EsViritu, you must point the software to its database by sourcing the generated environment file:
```bash
source /path/to/your/basedir/env_vars.sh

```


* **Run VAPER:** VAPER is executed via Nextflow. You do not need to pull containers for it manually; Nextflow will handle dependency fetching natively using `ociAutoPull` directives during your first run.

