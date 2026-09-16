# Wastewater SAE anomaly detection

Encode predicted proteins from wastewater metagenomes as sparse-autoencoder
feature fingerprints, and match them against ESM Atlas clusters. Two entry
points:

* **`sae/sae_testing_script.py`** — one query sequence (protein or nucleotide),
  end to end, with human-readable feature descriptions.
* **`sae/pipeline/`** — seven modular stages from reads to candidate clusters,
  resumable, each writing a provenance manifest.

## Requirements

**`uv` is required.** There is no pip fallback: `uv` provisions the exact
Python (3.12) this project is pinned to, so the environment never depends on
whatever `python3` the host ships. A pip fallback would quietly build against
a different interpreter and report success.

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh    # or: brew install uv
```

| Command | Needed for | |
|---|---|---|
| `uv` | venv + all Python deps | required |
| `curl` | every data fetch | required |
| `awk` | FASTQ record validation and mate pairing | required |
| `gzip` + `gzcat` (or `zcat`) | subsample compression / reading | required |
| `md5` (or `md5sum`) | verifying mate order after subsampling | required |
| `megahit` or `spades.py` | pipeline `s02_assemble` | no pip equivalent |
| `mmseqs` | `s04_derep` clustering below 100% identity | optional |
| `fastp` | `s01_qc`, ~30× faster than the Python path | optional |
| `aws` | faster, resumable atlas download (else `curl`) | optional |

`./setup.sh preflight` verifies the required set; `./setup.sh status` checks
those plus data and the optional binaries. Python packages are pinned in
`requirements.txt`.

The bioinformatics binaries come from bioconda:

```bash
conda install -c bioconda megahit mmseqs2 fastp
```

## Fresh checkout

```bash
git clone <url> && cd codeathon
./setup.sh                 # status: what is present, what to run next
./setup.sh quickstart      # ~2 GB — enough to run the pipeline end to end
```

`./setup.sh` on its own only reports; it never downloads. Every target is
idempotent and named, so nothing pulls 25 GB unless you ask for it.
See **[DATA.md](DATA.md)** for the size and source of every data item.

```bash
./setup.sh model-6b        # +25 GB, needed for feature *descriptions*
./setup.sh reads           # SARS-CoV-2 amplicon FASTQ
./setup.sh atlas           # ESM Atlas tables, 28.6 GB from the public S3 bucket
```

Without an assembler, pipeline stage `s02` blocks with instructions and the
rest still runs from contigs or proteins.

## Run it

```bash
cd sae
.venv/bin/python sae_testing_script.py --sequence MKTAYIAKQRQ...
.venv/bin/python sae_testing_script.py --fasta contig.fna --orfs all

cd pipeline
python run.py --fastq ../../data/fastq_rnaseq/SRR38294894_1.fastq.gz \
              --fastq2 ../../data/fastq_rnaseq/SRR38294894_2.fastq.gz --sample CHI-A
python run.py --contigs contigs.fa --sample S1 --from s03_genes
```

Or drive it from a browser — progress for every run, plus upload-and-launch:

```bash
python sae/web/server.py      # http://127.0.0.1:8765
container/run.sh web          # the same UI from inside the container
```

It reads the manifests each stage already writes, so runs started from the CLI
show up too. See **[sae/web/README.md](sae/web/README.md)**.

## Layout

```
setup.sh                      one entry point for env + data
requirements.txt              pinned Python deps
DATA.md                       provenance and sizes
data/
  fetch_fastq.sh              NCBI amplicon download
  subsample_rnaseq.sh         ENA streaming subsample
  *.csv, manifest.tsv         run metadata (committed)
  atlas_manifest.tsv          S3 keys + exact byte sizes for verification
sae/
  sae_testing_script.py       single-query tool
  README.md                   caveats, benchmarks, model choice
  codeathon_proposal.html     project proposal
  pipeline/                   s01…s07 + run.py, README
  web/                        progress dashboard + run launcher (stdlib only)
chicago_wastewater_sra*.{md,csv,txt}   cohort definition and accession list
```

## Read these before trusting a result

`sae/README.md` and `sae/pipeline/README.md` carry the caveats that matter —
measured throughput, the 2.3% coverage of the feature→cluster bridge, why the
6B model is required for interpretation but not for retrieval, and what the
CASPER subsample can and cannot support.
