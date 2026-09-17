# Sewage-to-Signal: Early Warning for Immune-Escape Mutations

**NIAID-BRCs AI Codeathon 2.0** · September 16–18, 2026 · Argonne National Laboratory

An AI-assisted wastewater surveillance workflow that detects pathogens of
concern, tracks their prevalence, and flags emerging immune-escape and
therapeutic-resistance mutations.

Project page: https://niaid-brc-codeathons.github.io/projects/sewage-to-signal/

---

> **The scope below is a draft pitch, not a plan.**
>
> The proposal is a one-slide pitch from the organizing team. It exists to
> seed a team, not to constrain one. Scope, methods, target organism, and
> success criteria are all still open — expect them to change substantially.
> Turning this into a real plan is the team's first job, and it lands in the
> project charter due August 28, 2026.

## Proposed scope

**Goal.** Detect pathogens of concern in wastewater, track their prevalence,
identify emerging mutations, and assess potential immune-escape or therapeutic
impact.

**Three-day MVP.** From target-capture wastewater sequencing: classify reads
and identify NIAID pathogens of concern; call protein mutations and track their
geographic, temporal and lineage distribution against public repositories; use
literature RAG plus curated databases to associate mutations with immune
escape, therapeutic resistance or altered pathogenicity; emit a
provenance-linked early-warning report.

**Evaluation.** Taxonomic classification accuracy; concordance of mutation
calls with standard pipelines; accuracy of geographic and temporal context;
precision of literature-derived mutation–phenotype associations; and recovery
of known escape or resistance mutations as high-priority signals.

## What is built so far

A different angle on the same problem, and deliberately narrower than the
pitch: encode predicted proteins from wastewater metagenomes as
sparse-autoencoder feature fingerprints and match them against ESM Atlas
clusters, so a protein with no conventional annotation still gets a functional
description.

* **`sae/sae_testing_script.py`** — one query sequence (protein or nucleotide),
  end to end, with human-readable feature descriptions.
* **`sae/pipeline/`** — seven modular stages from reads to candidate clusters,
  resumable, each writing a provenance manifest.
* **`sae/web/`** — progress dashboard, artifact browser, and run launcher.
* **`container/`** — Apptainer image for cluster use.

Where it does **not** yet meet the pitch, and these are the open gaps:

* **No taxonomic classification.** Nothing identifies pathogens of concern from
  reads; there is no Kraken/Bracken step. `s07_match` surfaces an
  `lca_taxonomy` for a matched cluster representative, which is a weak proxy.
* **No mutation calling.** `s05_prefilter` now separates proteins that are
  *fully* explained by a known family from those only partially explained, and
  the partial ones carry their family label — the hook for comparing a protein
  against its own family. The comparison itself is not implemented.
* **No literature RAG** and no phenotype association.
* **No host or rRNA depletion**, which untargeted wastewater RNA-seq needs.

See `sae/pipeline/README.md` for the measured limits behind each of these.

## Leads

* Alexander Taepper
* Andrew Warren

Team members get access through the
[NIAID-BRC-Codeathons](https://github.com/NIAID-BRC-Codeathons) organization.

---

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
| `megahit` or `spades.py` | `s02_assemble` — **only** way in from raw reads | no pip equivalent |
| `mmseqs` | `s04_derep` clustering below 100% identity | optional |
| `fastp` | `s01_qc`, ~30× faster than the Python path | optional |
| `aws` | faster, resumable atlas download (else `curl`) | optional |

`./src/setup.sh preflight` verifies the required set; `./src/setup.sh status` checks
those plus data and the optional binaries. Python packages are pinned in
`requirements.txt`.

The three bioinformatics binaries have no pip equivalent. The container carries
them pinned, which is the least painful route:

```bash
container/build.sh pull          # prebuilt image from GHCR
```

Otherwise `conda install -c bioconda megahit mmseqs2 fastp`. Everything except
assembly runs without them.

## Fresh checkout

```bash
git clone https://github.com/NIAID-BRC-Codeathons/sewage-to-signal.git
cd sewage-to-signal
./src/setup.sh                 # status: what is present, what to run next
./src/setup.sh quickstart      # ~2 GB — venv, deps, feature table, ESMC-300M, reads
./src/setup.sh pfam            # 400 MB — without it s05 is a pass-through
```

`./src/setup.sh` on its own only reports; it never downloads. Every target is
idempotent and named, so nothing pulls 25 GB unless you ask for it.
See **[DATA.md](DATA.md)** for the size and source of every data item.

```bash
./src/setup.sh model-6b        # +25 GB, needed for feature *descriptions*
./src/setup.sh reads           # SARS-CoV-2 amplicon FASTQ
./src/setup.sh atlas           # ESM Atlas tables, 28.6 GB from the public S3 bucket
```

## Run it

Everything below runs from the repo root and uses the venv interpreter
explicitly — a bare `python` will not have the dependencies.

**One sequence, end to end.** Prints the top SAE features with descriptions.

```bash
./sae/.venv/bin/python sae/sae_testing_script.py \
    --sequence MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ
```

**From a genome or contigs.** Works on a bare checkout — no assembler needed,
because you are supplying the contigs.

```bash
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512.2&rettype=fasta&retmode=text" > sars2.fa

./sae/.venv/bin/python sae/pipeline/run.py \
    --contigs sars2.fa --sample SARS2 \
    --hmm data/pfam/Pfam-A.hmm --bit-cutoffs gathering --model 300m
```

That is the validation case, so you can check the answer: s03 calls 9 genes,
and s05 discards 8 of them as fully explained by Pfam, leaving ORF1a — a
polyprotein that matches 17 Pfam domains and is still a quarter unaccounted
for. Results land in `work/SARS2/`.

**From raw reads.** Needs an assembler, so use the container; `/data` and
`/work` are bind mounts onto `data/` and `./work`.

```bash
container/run.sh pipeline \
    --fastq  /data/fastq_rnaseq/SRR38294894_1.fastq.gz \
    --fastq2 /data/fastq_rnaseq/SRR38294894_2.fastq.gz \
    --sample CHI-A \
    --hmm /data/pfam/Pfam-A.hmm --bit-cutoffs gathering --model 300m
```

**From a browser.** Progress for every run, an artifact browser, and
upload-and-launch.

The server runs *outside* the container, on the login node — it only reads
manifests and submits jobs, so it needs no allocation. Every run it starts goes
to the scheduler; nothing is forked from the server.

```bash
# on a cluster login node
./sae/.venv/bin/python sae/web/server.py                  # http://127.0.0.1:8765

# on a laptop: bring a scheduler, then the same command
container/slurm-local/up.sh                               # SLURM + Apptainer in Docker
eval "$(container/slurm-local/up.sh env)"
./sae/.venv/bin/python sae/web/server.py
```

`container/run.sh web` serves the same UI from inside the image, but it can
only browse — there is no `sbatch` in there, so it cannot launch.

It reads the manifests each stage already writes, so runs started from the CLI
show up too. See **[sae/web/README.md](sae/web/README.md)**.

### Which model, and where to run

`--model 300m` is the right default for a laptop. ESMC-6B is ~12 GB of weights
and will be OOM-killed in a 16 GB container; it also needs a GPU to be
practical. 300m gives you retrieval but no feature *descriptions*, since only
the 6B layer-60 SAE has a published description table.

On Apple Silicon, run `s06_embed` **natively** rather than in the container:
the container is `linux/amd64` under emulation with no MPS, so it falls back to
CPU. Measured on the same protein: 1.24 seq/s native on MPS against 0.33 seq/s
emulated. Assemble in the container, embed outside it.

## Layout

```
src/
  setup.sh                    one entry point for env + data
  setup_subtraction.sh        host-read subtraction provisioning
  README_subtraction.md       what the subtraction step does
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
