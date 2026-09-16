# SAE metagenomics pipeline

Modular stages turning wastewater reads into SAE feature profiles and candidate
ESM Atlas clusters. Each stage is independently runnable, writes a provenance
manifest, and is idempotent — re-running resumes rather than recomputes.

```
s01_qc → s02_assemble → s03_genes → s04_derep → s05_prefilter → s06_embed → s07_match
 reads     contigs        proteins     nr.faa     analyze.faa     parquet     clusters
```

## Why this shape

Reads are **151 bp** = a 50 aa peptide: too short to carry a domain, and one
indel garbles the frame. Embedding them directly is also impossible at volume —
CASPER alone is ~60.7 Tbp ≈ 400 billion reads, which at the measured
0.107 M seq/day would take ~12,000 years.

So the pipeline is a funnel, and **s05_prefilter is the main lever**: homology
search is orders of magnitude cheaper than an ESMC forward pass, and a protein
*fully explained* by a known family does not need an SAE.

### The three classes

A binary known/dark split discards too much — a hit can be statistically
overwhelming and still explain almost none of the protein — so s05 sorts on
**coverage** as well as significance:

| Class | Test | Fate |
|---|---|---|
| `known` | confident hit **and** covers ≥ `--min-coverage` | discarded |
| `partial` | a hit exists, but weak **or** incomplete | analysed, **with** its family label |
| `dark` | no significant hit | analysed, no prior |

Against a 150 aa HMM built from Spike residues 301–450:

```
gene_id        category  family         evalue     coverage  aa_len
NC_045512.2_3  partial   SpikeFragment  7.17e-105  0.1178    1273
```

E=7e-105 is past any confidence threshold, yet 88% of the protein is
unaccounted for. A binary rule calls that `known` and throws it away.

Coverage is the union of *all* significant domain envelopes across *all*
families, so a genuine multi-domain protein is called complete rather than
penalised for matching several models.

All three FASTAs are written for audit; `analyze.faa` (partial + dark) flows to
s06, and `classification.tsv` records `gene_id, category, family, family_acc,
evalue, coverage, n_domains, aa_len`. The family label on `partial` is the
conditioning key for asking whether a protein is unusual *for its own family* —
as opposed to `dark`, which has no prior and can only be compared against the
global atlas.

With Pfam-A prefer `--bit-cutoffs gathering`; the curated per-family thresholds
beat any flat E-value.

## Usage

```bash
# paired (CASPER metagenomes, streamed from ENA as split mates)
python run.py --fastq ../../data/fastq_rnaseq/SRR38294894_1.fastq.gz \
              --fastq2 ../../data/fastq_rnaseq/SRR38294894_2.fastq.gz --sample CHI-A

# single-file (the NCBI endpoint concatenates mates rather than splitting)
python run.py --fastq ../../data/fastq/SRR40033132.fastq.gz --sample S0153
python run.py --contigs contigs.fa --sample S1 --from s03_genes
python run.py --proteins prot.faa  --sample S1 --model 300m    # CPU-sized
python s03_genes.py contigs.fa --sample S1       # stages run standalone too
```

`--from`/`--to` bound the range, `--force` overrides the idempotency check.
`../web/` drives the same thing from a browser.

## Stages

| Stage | Tool | Notes |
|---|---|---|
| s01_qc | pyfastx, or fastp if present | paired via `--fastq2` |
| s02_assemble | MEGAHIT / metaSPAdes | **external binary required** |
| s03_genes | pyrodigal | pure wheel |
| s04_derep | exact hash, MMseqs2 if present | exact-only without MMseqs2 |
| s05_prefilter | pyhmmer (`--hmm`) or pyswrd (`--ref`) | pass-through if neither given |
| s06_embed | ESMC + SAE | `--model {6b,300m}` |
| s07_match | pyarrow join | — |

s02 is the only stage needing something `requirements.txt` cannot install.
`../../container/` carries megahit, mmseqs2 and fastp pinned, which is the
simplest way to get them.

### s06 details

Sorting by length before batching keeps padding waste low, though on Apple
Silicon the 6B gains only ~20% from it (1.03 → 1.24 seq/s), being
memory-bandwidth-bound; on CUDA expect the usual gain. Top-K is taken straight
off the sparse COO values with a scatter-reduce into a batch × codebook buffer —
2 MB at batch=32 regardless of sequence length, rather than densifying an
L × 16,384 matrix per sequence. Output is long-format parquet
`(gene_id, feature_id, activation, raw_activation)`.

## Known limits

**Sequences are truncated at `--max-len` (default 1022).** SARS-CoV-2 ORF1a is
4405 aa, so it is cut. Long proteins need chunking with overlap; not implemented.

**The cluster join covers ~2.3%.** s07 bridges feature → UniRef → local cluster
because the local tables carry no SAE-feature column. The durable fix is to run
s06 over the 7.7 M atlas representatives once and build a real feature → cluster
index — ~4 days at the measured 300M throughput on this Mac, far less on a GPU.
**Highest-value next step.**

**Model choice is a memory decision as much as a quality one.** Retrieval needs
only that both sides use the *same* SAE, so ESMC-300M is fine and 20× cheaper;
only human-readable descriptions need the 6B layer-60 SAE, the one variant with
a published description table. ESMC-6B is ~12 GB of bfloat16 weights, so in a
16 GB Docker VM s06 is SIGKILLed during model load before printing anything —
the log simply stops after s05. Use `--model 300m` there (0.33 seq/s measured on
emulated CPU). `--model` keeps backbone, SAE repo and layer consistent;
`--backbone`/`--sae-repo`/`--layer` override individually.

**`--ref`/pyswrd does not run on Apple Silicon.** `pyswrd.search` raises
`RuntimeError: no supported SIMD backend available` via pyopal on arm64, so
`--hmm` is the only working backend there. Where pyswrd does run it yields
coverage only if its result type carries query bounds; without them coverage
stays `None` and the protein is never discarded.

**No rRNA depletion.** Untargeted wastewater RNA-seq is rRNA-dominated,
pyrodigal calls spurious short ORFs on it, and those consume GPU time in s06.
A SortMeRNA/SILVA filter between s02 and s03 is worth adding before any real run.

## Validation

End-to-end against the SARS-CoV-2 reference (NC_045512.2), given only the
nucleotide genome: **s03** recovers Spike at 21563–25384 (1273 aa), matching the
canonical annotation exactly, plus ORF1a, ORF1b and ORF3a. **s06/s07**
independently describe that gene as a viral fusion glycoprotein — f12894
"entry/fusion envelope ectodomains", f15507 "viral envelope/spike/fusion
proteins", f6802 "viral envelope/fusion glycoproteins", f8286 "fusion peptides,
membrane-proximal/TM anchors". Nothing in the pipeline was told this was a
coronavirus.

Paired QC on SRR38294894: 1,000,000 pairs in, 998,003 kept (99.8%), both mates
in identical name order — a pair is dropped when *either* mate fails. That took
126 s in pure Python (~8 k pairs/s), so a full CASPER run would be ~10 h;
install `fastp` before running one at scale.

Throughput on the SARS-CoV-2 run was 0.31 seq/s against 1.03 benchmarked at
200 aa — those proteins run to the 1022 aa cap and attention is quadratic in
length. Budget by residue, not by protein count.

## Data note

`../../data/fastq/` — 12 runs, all **SARS-CoV-2 amplicon**: a good positive
control, useless for discovery.

`../../data/fastq_rnaseq/` — **SRR38294894**, a CASPER untargeted metagenome
(PRJNA1247874, RNA-Seq, NovaSeq X, influent CHI-A, 2026-02-01), subsampled to
1 M **paired** reads = 302 Mbp, 0.11% of the 263 Gbp run.

FracMinHash k-mer spectrum (k=21) over 400 k reads:

| multiplicity | 1 | 2 | 3–5 | 6–20 | >20 |
|---|---|---|---|---|---|
| % of distinct k-mers | 81.5 | 9.1 | 5.0 | 2.8 | 1.6 |

Mean multiplicity 6.54; singletons hold only 12.5% of k-mer *occurrences*. Two
populations: 81.5% of the diversity sits near 1x coverage and will not assemble,
while the ~4.4% above 6x carries most of the sequence mass — in untargeted
RNA-seq, rRNA plus highly expressed transcripts. So the subsample is enough to
exercise the pipeline on real contigs, genes and proteins, and **not enough for
discovery**, because the novel long tail is exactly what falls below assembly
threshold. It is also a HEAD subsample, so flowcell-position biased.
