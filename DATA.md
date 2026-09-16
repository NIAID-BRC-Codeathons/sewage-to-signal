# Data inventory

Nothing in this table is committed — together it is ~30 GB. `./setup.sh status`
tells you what is present; each row names the target that fetches it.

| Item | Size | Target | Source |
|---|---|---|---|
| SAE feature descriptions | 50 MB | `features` | HF dataset `biohub/ESMC-SAE-Features` |
| ESMC-300M + SAE | 1.7 GB | `model-small` | HF `biohub/ESMC-300M`, `…-sae-layer23-k64-codebook16384` |
| ESMC-6B + SAE | 25 GB | `model-6b` | HF `biohub/ESMC-6B`, `…-sae-layer60-k64-codebook16384` |
| SARS-CoV-2 amplicon FASTQ | 2.4 GB | `reads` | NCBI SRA via `data/fetch_fastq.sh` |
| CASPER metagenome subsample | 100 MB | `rnaseq` | ENA stream via `data/subsample_rnaseq.sh` |
| Pfam-A HMMs | 399 MB gz | `pfam` | EBI `ftp.ebi.ac.uk/pub/databases/Pfam/current_release` |
| ESM Atlas cluster tables | 28.6 GB | `atlas` | AWS S3 `esm-protein-atlas` (public, unauthenticated) |

Model weights go to `~/.cache/huggingface`, which is shared between checkouts,
so a second clone of this repo re-downloads nothing.

## Pfam-A

`s05_prefilter` is a **pass-through without an HMM set** — every protein counts
as dark and goes to the GPU stage, which is the expensive path the funnel exists
to avoid. So this is not an optional extra for any real run:

```bash
./setup.sh pfam
python run.py ... --hmm data/pfam/Pfam-A.hmm --bit-cutoffs gathering
```

Pfam-A is the right default because s05 asks *"is this already explained"* —
breadth beats specificity, and a viral-only set such as VOGdb would leave every
bacterial protein dark and flood s06. Pfam also ships curated per-family
gathering thresholds, which `--bit-cutoffs gathering` uses in place of a flat
E-value.

The download is 399 MB compressed and is decompressed on arrival, since pyhmmer
reads it uncompressed. `current_release` moves, so `relnotes.txt` is kept
alongside to record which release was taken.

## Minimum to be useful

```bash
./setup.sh quickstart      # ~2 GB: venv, deps, feature table, ESMC-300M, CASPER subsample
```

That runs the whole pipeline end to end. Add `model-6b` when you need
human-readable feature descriptions — the 6B layer-60 SAE is the only variant
with a published description table, so it is the only one whose feature ids
mean anything. ESMC-300M is fine for retrieval, where both sides only need to
use the *same* SAE, and is 20x cheaper.

## ESM Atlas cluster tables

`sae/representative_proteins.parquet` (653 MB, 7,723,579 rows) and
`sae/sae_clusters/cluster_members_{0..9,a..f}.parquet` (16 x ~1.7 GB, ~811 M
rows) come from the public ESM Atlas bucket:

```bash
./setup.sh atlas
```

Equivalent to, and verified against, the upstream commands:

```bash
aws s3 sync --no-sign-request \
  s3://esm-protein-atlas/v1/clusters/indexes/secondary/cluster_members/ ./sae_clusters/
aws s3 cp --no-sign-request \
  s3://esm-protein-atlas/v1/clusters/data/representative_proteins.parquet ./
```

`setup.sh atlas` adds three things over the raw commands:

* **Exact size verification.** Every object's byte size is pinned in
  `data/atlas_manifest.tsv`, so a truncated transfer is detected rather than
  silently accepted. All 17 local files were confirmed byte-exact against the
  bucket.
* **Resumability at file granularity.** Only files whose size does not match
  are re-fetched, so an interrupted 28 GB download continues where it stopped.
* **No hard `aws` dependency.** It uses `aws s3 cp` when available (faster,
  resumable within a file) and otherwise plain `curl` against the bucket's
  HTTPS endpoint. Both paths are tested and produce identical, verified bytes.

Without these tables only pipeline stage `s07_match` is affected — it has no
clusters to join against. Every other stage runs normally; pass
`--skip-clusters` to `sae_testing_script.py` in the meantime.
