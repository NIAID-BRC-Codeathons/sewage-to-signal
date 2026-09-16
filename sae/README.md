# SAE → ESM Atlas cluster lookup

Takes a query protein sequence, extracts its top sparse-autoencoder (SAE)
features from an ESMC protein language model, and looks up candidate ESM Atlas
clusters in the local parquet tables.

## Quick start

```bash
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python esm pandas pyarrow

# demo sequence, default 6B backbone
.venv/bin/python sae_testing_script.py

# your own sequence
.venv/bin/python sae_testing_script.py --sequence MKTAYIAK... --top-k 10
.venv/bin/python sae_testing_script.py --fasta query.faa --json-out out.json
```

First run downloads ~25 GB (ESMC-6B) + 0.34 GB (SAE) + 50 MB (feature table)
into `~/.cache/huggingface`.

## Nucleotide input

ESMC is a *protein* language model — it has no nucleotide vocabulary, so DNA
must be translated before it can be embedded. The script does that for you:

```bash
# auto-detected and translated
.venv/bin/python sae_testing_script.py --fasta contig.fna

# every ORF, not just the longest
.venv/bin/python sae_testing_script.py --fasta contig.fna --orfs all --min-orf-aa 60
```

Input type is auto-detected (≥90% ACGT/U/N over ≥20 characters) and can be
forced with `--input-type`. Genes are called with **pyrodigal** in metagenomic
mode, which is the right tool for mixed-organism contigs: it handles gene
boundaries, ribosome binding sites, and flags genes running off the contig edge
as partial. `pip install pyrodigal` — it is a pure wheel, no external binary.

If pyrodigal is missing, or calls nothing above the length cutoff, the script
falls back to plain six-frame translation split on stop codons. Start codons
are *not* required there, because metagenomic contigs are routinely fragments
and demanding a Met start would discard genuine partial genes. The fallback is
noisier — on a 2 kb test contig pyrodigal called 3 genes where six-frame
produced 41 overlapping ORFs — so prefer the gene caller.

The backbone is loaded once and reused across ORFs, so `--orfs all` does not
pay the 24 GB load repeatedly.

| Flag | Meaning |
|---|---|
| `--input-type` | `auto` (default), `protein`, `nucleotide` |
| `--orfs` | `longest` (default) or `all` |
| `--min-orf-aa` | minimum ORF length, default 30 |
| `--no-gene-caller` | force six-frame translation |

Caveats: translation uses the standard genetic code (tables 1 and 11 share the
same codon→amino-acid mapping; they differ only in permitted start codons,
which is the gene caller's concern). Organisms with a divergent code —
*Mycoplasma*, ciliates, some mitochondria — will mistranslate. Frameshifts and
sequencing errors in raw reads truncate ORFs at spurious stops, so assembled
contigs give much better results than raw FASTQ reads.

## What the script does

1. **Features.** Runs the query through `biohub/ESMC-6B` with the
   `biohub/ESMC-6B-sae-layer60-k64-codebook16384` SAE attached to layer 60.
   That yields a 16,384-wide sparse activation per residue, 64 non-zero.
   Activations are max-pooled over residues (`<cls>`/`<eos>` excluded) and
   scaled by the SAE's shipped `idf / max` buffers, then the top-K are kept.
   Ranking is IDF-weighted (see below).
2. **Annotation.** Looks those feature ids up in the `biohub/ESMC-SAE-Features`
   dataset, which gives each feature a `summary`, `category`,
   `exemplar_protein_families` and a reliability `threshold`.
3. **Cluster voting.** Each top feature nominates the UniRef90 proteins it
   fires hardest on (`top_100_uniref_ids`). Those accessions are matched
   against `representative_proteins.parquet` via `uniref_match_accession`.
   Clusters are ranked by how many distinct query features nominated them,
   then by summed activation.

## Important caveats

**The local tables have no SAE-feature column.** `representative_proteins.parquet`
carries `protein_hash`, pfam/taxonomy annotations and a `uniref_match_accession`;
`sae_clusters/cluster_members_*.parquet` is just `cluster_rep_protein_hash →
protein_hash`. So features cannot be matched to clusters directly — the script
bridges through UniRef accessions instead.

**That bridge is sparse.** Of the ~1.55 M distinct UniRef accessions named in
the feature table, only 35,569 (**2.3 %**) appear in this local atlas subset. A
query typically recovers a handful of clusters, not an exhaustive list. A
cluster's absence is *not* evidence that it lacks the feature. For full
coverage you would need SAE features computed over the atlas representatives
themselves, which is not in this dataset.

**The SAE ships placeholder normalization buffers.** The layer-60 SAE's `idf`
and `max` buffers are all ones, so the library's own `normalize_sae=True` is a
no-op. The real UniRef90 statistics (`uniref90_idf` 0→5.0,
`uniref90_max_activation` 7.8→50.7) live in the feature table, and the script
applies those *before* top-K selection. This matters: ranking by raw activation
alone puts generic "meaning not known" features on top, because IDF is exactly
what downweights features that fire on nearly every protein. Use `--no-idf` to
rank by raw activation instead.

**Feature ids are SAE-specific.** The published descriptions apply only to the
6B layer-60 k64 codebook-16384 SAE, which is why it is the default. Pointing
`--sae-repo`/`--layer` elsewhere keeps the ids valid but makes the annotations
meaningless, so the script drops them and skips the cluster join.

**Precision.** The 6B backbone loads in bfloat16 by default (25 GB in fp32 does
not fit a 32 GB machine). Without Transformer Engine / flash-attn installed the
library warns that outputs differ numerically; after the final LayerNorm this
is within rounding noise. Use `--dtype float32` on a larger machine.

## Useful flags

| Flag | Meaning |
|---|---|
| `--sequence` / `--fasta` | query input (defaults to a demo sequence) |
| `--top-k` | number of SAE features to keep (default 5) |
| `--backbone`, `--sae-repo`, `--layer` | swap the model/SAE variant |
| `--uniref-per-feature` | how many of each feature's top-100 UniRef hits to use |
| `--max-clusters` | how many clusters to print (default 15) |
| `--member-counts` | count cluster members; scans ~27 GB, slow |
| `--skip-clusters` | stop after feature extraction + annotation |
| `--no-idf` | rank by raw activation instead of IDF-weighted |
| `--device`, `--dtype` | `auto` picks mps/cuda/cpu and bf16 for 6B |
| `--offline` | use only the local HF cache (skips ETag revalidation) |
| `--verbose` | show HF progress bars and esm kernel warnings |
| `--json-out` | write full results as JSON |

## Startup noise

Two things print on every run and neither indicates a problem:

**"Fetching N files"** is *not* a re-download. `huggingface_hub` revalidates the
already-populated cache on each call; the bar completes instantly (11 files at
~4500 it/s). Measured cost is 0.465 s online vs 0.212 s offline — about 0.25 s
of ETag checks against a 24 GB model load. `--offline` skips it entirely.

**The two `ESMC:` kernel warnings** are `logger.warning` calls fired once when
`esm.models.esmc` is imported, reporting that Transformer Engine and
flash-attn/xformers are missing. On Apple Silicon these are CUDA-only packages,
so the warnings cannot be resolved by installing anything — and per the message
itself, after the final LayerNorm the numerical difference is a few ULP. Both
are suppressed by default; `--verbose` brings them back.

## Sanity checks

**The join is biologically consistent.** Feature 0 is described in the feature
table as the Nudix hydrolase catalytic motif. Feeding it into the cluster join
returns clusters annotated "putative NUDIX hydrolase" /
"NUDIX domain-containing protein" — not arbitrary ones.

**Independent features converge on one theme.** On the demo sequence, the
IDF-weighted top features are dominated by membrane-proximal / juxtamembrane
sensor descriptions, and the clusters they independently nominate are sensor
histidine kinases (PF21623, PF14501), Cache sensor domains (PF02743),
adenylate/guanylate cyclase (PF00211) and methyl-accepting chemotaxis PDC
sensor domains (PF22673) — all bacterial membrane signal transduction. Three
different features (4643, 9171, 4651) land on the same theme, and PF02743
appears from two of them.

## Measured throughput (M2-class Mac, MPS, 200 aa sequences)

| Model | batch=1 | batch=8 | batch=32 | best M seq/day |
|---|---|---|---|---|
| ESMC-6B (bf16) | 1.03 seq/s | 1.24 seq/s | 1.17 seq/s | **0.107** |
| ESMC-300M (fp32) | 14.6 seq/s | 23.0 seq/s | 23.1 seq/s | **1.99** |

Two things this says:

**Batching barely helps the 6B on MPS** (+20% at batch=8, and batch=32 is
*slower* than batch=8). At 12.7 GB of bf16 weights against 32 GB of unified
memory, each forward streams the whole weight set through memory — it is
bandwidth-bound, not compute-bound, so bigger batches buy almost nothing and
eventually cost. The 300M does benefit (1.6x, saturating at batch=8).

**The 20x gap between the models tracks the 20x parameter ratio**, as expected
for a bandwidth-bound workload.

Indexing the 7.7 M atlas representatives on this machine: ~72 days at 6B,
~4 days at 300M. GPU figures would be far better but are not measured here.

### Retrieval vs interpretation

These need different things, and conflating them is expensive:

* **Retrieval** — matching a query against atlas clusters by feature
  fingerprint — only requires that *both sides use the same SAE*. Feature ids
  need to be consistent, not meaningful. ESMC-300M works fine and is 20x
  cheaper.
* **Interpretation** — "what does feature 4643 mean?" — requires the 6B
  layer-60 k64 codebook-16384 SAE, the only variant with published
  descriptions.

So a practical split is: index the atlas with 300M for retrieval, and run 6B
only on the query and its top hits, where human-readable feature descriptions
actually matter.

## Data files

| File | Rows | Contents |
|---|---|---|
| `representative_proteins.parquet` | 7,723,579 | one row per cluster representative: pfam, taxonomy, product name, UniRef match |
| `sae_clusters/cluster_members_*.parquet` | ~811 M | `cluster_rep_protein_hash → protein_hash` membership, 16 shards |
