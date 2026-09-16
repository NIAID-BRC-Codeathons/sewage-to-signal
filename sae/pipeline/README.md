# SAE metagenomics pipeline

Modular stages turning wastewater reads into SAE feature profiles and candidate
ESM Atlas clusters. Each stage is independently runnable, writes a provenance
manifest, and is idempotent — re-running resumes rather than recomputes.

```
s01_qc → s02_assemble → s03_genes → s04_derep → s05_prefilter → s06_embed → s07_match
 reads     contigs        proteins     nr.faa     analyze.faa     parquet     clusters
```

## Why this shape

Reads are **151 bp** = a 50 aa peptide. Too short to contain a domain, and one
indel garbles the frame. Embedding reads directly is both meaningless and
impossible at volume — CASPER alone is ~60.7 Tbp ≈ 400 billion reads, which at
the measured 0.107 M seq/day would take ~12,000 years.

So the pipeline is a funnel, and **s05_prefilter is the main lever**. Homology
search is orders of magnitude cheaper than an ESMC-6B forward pass, and a
protein that is *fully explained* by a known family does not need an SAE. The
SAE earns its cost on everything else, which is also the part you actually care
about for surveillance.

### The three classes

A binary known/dark split discards too much. A hit can be statistically
overwhelming and still explain almost none of the protein, so s05 sorts on
**coverage** as well as significance:

| Class | Test | Fate |
|---|---|---|
| `known` | confident hit **and** covers ≥ `--min-coverage` of the sequence | discarded |
| `partial` | a hit exists, but it is weak **or** incomplete | analysed, **with** its family label |
| `dark` | no significant hit | analysed, no prior |

`known.faa`, `partial.faa` and `dark.faa` are all written for audit;
`analyze.faa` (= partial + dark) is what flows to s06. A
`classification.tsv` records `gene_id, category, family, family_acc, evalue,
coverage, n_domains, aa_len` for every protein.

Coverage is the union of *all* significant domain envelopes across *all*
families, so a genuine multi-domain protein is correctly called complete rather
than penalised for matching several models.

The `partial` class is the interesting one. Those proteins carry a family
label, which is the conditioning key for reference-relative comparison — "is
this Spike unusual *for a Spike*" — as opposed to `dark`, which has no prior
and can only be compared against the global atlas.

Worked example, against a 150 aa HMM built from Spike residues 301–450:

```
gene_id        category  family         evalue     coverage  aa_len
NC_045512.2_3  partial   SpikeFragment  7.17e-105  0.1178    1273
```

E=7e-105 is far past any confidence threshold, yet 88% of the protein is
unaccounted for. Under the old binary rule this was `known` and thrown away.

With Pfam-A, prefer `--bit-cutoffs gathering` — Pfam's curated per-family
thresholds are a better significance test than any flat E-value.

## Usage

```bash
# paired (CASPER metagenomes, streamed from ENA as split mates)
python run.py --fastq ../../data/fastq_rnaseq/SRR38294894_1.fastq.gz \
              --fastq2 ../../data/fastq_rnaseq/SRR38294894_2.fastq.gz --sample CHI-A

# single-file (NCBI endpoint returns concatenated reads, not split mates)
python run.py --fastq ../../data/fastq/SRR40033132.fastq.gz --sample S0153
python run.py --contigs contigs.fa  --sample S1 --from s03_genes
python run.py --proteins prot.faa   --sample S1 --from s06_embed
python s03_genes.py contigs.fa --sample S1      # stages run standalone too
```

`--from`/`--to` bound the range; `--force` overrides the idempotency check.

## Stages

| Stage | Tool | Status |
|---|---|---|
| s01_qc | pyfastx (fastp if present) | works, no external binary; paired via `--fastq2` |
| s02_assemble | MEGAHIT / metaSPAdes | **requires external binary** |
| s03_genes | pyrodigal | works, pure wheel |
| s04_derep | exact hash; MMseqs2 if present | works (exact only without MMseqs2) |
| s05_prefilter | pyhmmer (`--hmm`) or pyswrd (`--ref`) | works; three-class triage; pass-through if neither given |
| s06_embed | ESMC + SAE | works |
| s07_match | pyarrow join | works |

`pip install pyfastx pyrodigal pyhmmer pyswrd` covers everything except the
assembler. There is no pure-Python metagenome assembler, so s02 shells out;
on this machine `conda` needs an interactive `conda tos accept` first.

## What s06 does differently from the prototype

* **Length-bucketed batching.** Sorting by length keeps padding waste low.
  Measured on Apple Silicon the 6B gains only ~20% from batching (1.03 → 1.24
  seq/s) because it is memory-bandwidth-bound; on CUDA expect the usual gain.
* **No densification.** The prototype called `.to_dense()`, materialising an
  L × 16,384 matrix per sequence plus a second copy for normalisation. s06
  takes top-K straight off the sparse COO values with a scatter-reduce into a
  batch × codebook buffer — 2 MB at batch=32 regardless of sequence length.
* **Long-format parquet** `(gene_id, feature_id, activation, raw_activation)`
  instead of per-query JSON.

## Known limits

**Sequences are truncated at `--max-len` (default 1022).** SARS-CoV-2 ORF1a is
4405 aa, so it is cut. Long proteins need chunking with overlap; not yet
implemented.

**The cluster join covers ~2.3%.** s07 bridges feature → UniRef → local
cluster because the local tables have no SAE-feature column. The durable fix is
to run s06 over the 7.7 M atlas representatives once and build a real
feature → cluster index. At the measured 300M throughput that is ~4 days on
this Mac, far less on a GPU, and it replaces the 2.3% bridge with full
coverage. **This is the highest-value next step.**

**The `--ref`/pyswrd backend does not run on Apple Silicon.** `pyswrd.search`
raises `RuntimeError: no supported SIMD backend available` via pyopal on arm64,
so `--hmm` is the only working prefilter backend on this machine. When pyswrd
does run, it yields coverage only if its result type carries query bounds;
without them coverage stays `None` and the protein is never discarded.

**Retrieval and interpretation need different models.** Matching a query to
clusters only requires both sides use the *same* SAE — ESMC-300M is fine and
20× cheaper. Only human-readable feature descriptions require the 6B layer-60
SAE, the single variant with a published description table.

## Validation

End-to-end against the SARS-CoV-2 reference (NC_045512.2), given only the
nucleotide genome:

* **s03** recovers Spike at 21563–25384 (1273 aa), matching the canonical
  annotation exactly, plus ORF1a, ORF1b and ORF3a.
* **s06/s07** independently describe that gene as a viral fusion glycoprotein.
  Spike's top features: f12894 "entry/fusion envelope ectodomains", f15507
  "viral envelope/spike/fusion proteins", f6802 "viral envelope/fusion
  glycoproteins", f8286 "fusion peptides, membrane-proximal/TM anchors".

Nothing in the pipeline was told this was a coronavirus.

Paired QC verified on SRR38294894: 1,000,000 pairs in, 998,003 kept (99.8%),
both mates 998,003 reads with identical name order — pairs stay synchronised
because a pair is dropped when *either* mate fails. That took 126 s in pure
Python (~8 k pairs/s); a full CASPER run would be ~10 h, so install `fastp`
before running one at scale.

Throughput on the SARS-CoV-2 run was 0.31 seq/s, below the 1.03 seq/s benchmarked at
200 aa — these proteins run to the 1022 aa cap and attention is quadratic in
length. Budget by residue, not by protein count.

## Data note

`../../data/fastq/` — 12 runs, all **SARS-CoV-2 amplicon**. Assembling them
yields SARS-CoV-2 genes: a good positive control, useless for discovery. These
are single-file (the NCBI endpoint concatenates mates rather than splitting).

`../../data/fastq_rnaseq/` — **SRR38294894**, a CASPER untargeted metagenome
(PRJNA1247874, RNA-Seq, NovaSeq X, influent CHI-A, 2026-02-01), subsampled to
1 M **paired** reads = 302 Mbp, i.e. 0.11% of the 263 Gbp run.

### How much can that subsample actually support?

FracMinHash k-mer spectrum (k=21) over 400 k reads:

| multiplicity | % of distinct k-mers |
|---|---|
| 1 (singleton) | 81.5% |
| 2 | 9.1% |
| 3–5 | 5.0% |
| 6–20 | 2.8% |
| >20 | 1.6% |

Mean multiplicity 6.54; singletons hold only 12.5% of k-mer *occurrences*.

Two populations. 81.5% of the diversity sits near 1x coverage and will not
assemble. The ~4.4% above 6x carries most of the sequence mass — the abundant
fraction, which in untargeted RNA-seq means rRNA plus highly expressed
transcripts. So this subsample is **enough to exercise the pipeline on real
contigs, genes and proteins**, and **not enough for discovery**, because the
novel long tail is exactly what falls below assembly threshold. It is also a
HEAD subsample, so flowcell-position biased.

### rRNA

Untargeted wastewater RNA-seq is typically rRNA-dominated. pyrodigal will call
spurious short ORFs on rRNA, and those would consume GPU time in s06. A
SortMeRNA/SILVA filter between s02 and s03 is worth adding before any real run.
Not yet implemented.
