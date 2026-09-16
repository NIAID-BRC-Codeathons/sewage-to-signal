# Chicago wastewater surveillance — SRA FASTQ runs, last 12 months

**Query window:** 2025-09-16 → 2026-09-16 (sample collection date)
**Retrieved:** 2026-09-16 via NCBI E-utilities (`esearch`/`efetch` on `db=sra`)
**Result:** **381 run accessions** across 3 BioProjects. All Illumina, all paired FASTQ-retrievable.

Files in this directory:
- `chicago_wastewater_sra_accessions.csv` — full table (run, experiment, biosample, bioproject, collection date, site, strategy, platform, Gbp)
- `chicago_wastewater_sra_runs.txt` — bare run accession list, one per line (for `fasterq-dump`/`prefetch`)

## Projects included

| BioProject | Submitter | Description | Strategy | Platform | Runs | Collection range | Volume |
|---|---|---|---|---|---|---|---|
| [PRJNA1247874](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1247874) | SecureBio (CASPER) | Untargeted metagenomic sequencing of Chicago influent wastewater | RNA-Seq | NovaSeq X | 206 | 2025-09-21 → 2026-06-29 | ~60.7 Tbp |
| [PRJNA989260](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA989260) | Illinois WBE / CDC NWSS | NWSS Illinois project, **jurisdiction `CI` (City of Chicago)** only | AMPLICON (SARS-CoV-2) | NextSeq 2000 | 144 | 2025-09-16 → 2026-01-28 | ~72 Gbp |
| [PRJNA1498729](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1498729) | RIPHL | NWSS "Chicago Wastewater Surveillance" | AMPLICON (SARS-CoV-2) | MiSeq i100 | 31 | 2026-07-29 → 2026-09-03 | ~20 Gbp |

### Sampling sites

**CASPER (PRJNA1247874)** — 5 sites, ~41 runs each:
CHI-A, CHI-B, CHI-D1, CHI-D2 (influent WWTP), CHI-C (airport sewer manhole).

**NWSS Chicago (PRJNA989260 `CI` + PRJNA1498729)** — NWSS `collection_site_id`, with sewershed population:

| Site | Population | Runs (989260) | Runs (1498729) |
|---|---|---|---|
| S0153 | 721,207 | 19 | 4 |
| S0154 | 467,536 | 16 | 4 |
| S0155 | 125,995 | 17 | 2 |
| S0093 | 83,455 | 17 | 2 |
| S0095 | 77,799 | 18 | 3 |
| S0090 | 54,362 | 2 | 6 |
| S0120 | 29,370 | 18 | 4 |
| S0092 | 24,099 | 1 | 0 |
| S0094 | 23,475 | 16 | 1 |
| S0086 | 3,816 | 1 | 0 |
| S0108 | n/a | 19 | 5 |

## Coverage gap

The SARS-CoV-2 amplicon stream has a **~6-month hole (Feb–Jul 2026)**: PRJNA989260 Chicago (`CI`) submissions stop at collection date 2026-01-28, and the RIPHL project PRJNA1498729 only begins 2026-07-29. This looks like a handover of the Chicago NWSS sequencing from the Illinois WBE lab to RIPHL, with nothing deposited in between. Only CASPER metagenomes cover that period (through 2026-06-29).

## Deliberately excluded

| Candidate | Runs | Why excluded |
|---|---|---|
| PRJNA1290815 "Egan Nutrient Removal Pilot" (Black & Veatch) | 55 | `geo_loc_name: USA: Chicago`, but source is a *phosphorus removal bioreactor* — treatment-process research, not pathogen surveillance |
| PRJNA957477 (Verily) | 95 | `ww_surv_jurisdiction: IL`, single site, sewershed pop 86,000 — an Illinois plant, not Chicago; no city/site field to confirm otherwise |
| PRJNA989260 jurisdiction `IL`/`il` | 3,402 | Statewide IDPH submissions with no city attribute. Site IDs are jurisdiction-scoped, so `IL:S0094` ≠ `CI:S0094`; could not be attributed to Chicago |
| PRJNA1451225 | 40 | Stream sediment, Tordera, Spain |
| PRJNA904373 | 7 | *Candidozyma auris* clinical WGS (matched only on submitter RIPHL) |

## Reproducing the accession list

```bash
# 1. keyword + release-date sweep (release date >= collection date, so this is a superset)
esearch -db sra -query 'wastewater AND chicago AND 2025/09/16:2026/09/16[PDAT]' \
  | efetch -format runinfo           # 1,054 runs, 4 BioProjects

# 2. parent projects fetched in full and filtered on sample attributes,
#    since keyword matching alone undercounts (PRJNA989260 has 4,033 runs in window)
esearch -db sra -query 'PRJNA989260 AND 2025/09/16:2026/09/16[PDAT]' | efetch -format xml
#    -> keep ww_surv_jurisdiction == "CI"
esearch -db sra -query 'PRJNA1247874 AND 2025/09/16:2026/09/16[PDAT]' | efetch -format xml
#    -> keep sample TITLE containing "Chicago"
```

## Downloading

**ENA does not have these files.** All 381 accessions are *registered* in ENA, but
`fastq_ftp` is empty for every one of them — the mirror lags recent submissions
(13,705 of PRJNA989260's 16,277 runs have files, all older ones). Use NCBI.

`data/fetch_fastq.sh` pulls runs from NCBI's `sra-reads-be/fastq` endpoint, which
needs no toolkit installed:

```bash
./data/fetch_fastq.sh                          # the 4 default samples
./data/fetch_fastq.sh SRR40033050 SRR40032848  # or specific accessions
```

Endpoint caveats, both verified here:

- It returns **one concatenated gzip**, not split `_1`/`_2` mate files. If you need
  mates separated, install sra-tools and use `fasterq-dump --split-3` instead.
- It sends **no `Accept-Ranges` and no `Content-Length`**, so byte-range resume is
  impossible — `curl -C -` fails with error 56 and an interrupted transfer must
  restart that accession from zero. The script therefore downloads to `<acc>.part`,
  verifies (`gzip -t` plus a line count divisible by 4), then atomically renames, and
  records success in `fastq/manifest.tsv` so re-runs skip finished work.

### Downloaded so far

**Amplicon** — `data/fastq/`, 2.5 GB, 12 runs. Metadata in `samples.csv`, verified
sizes and read counts in `manifest.tsv`. Covers all 11 Chicago sewersheds and
6 of the 8 months that have data.

| Run | Project | Collected | Site | Reads |
|---|---|---|---|---|
| SRR40032820 | PRJNA989260 | 2025-09-25 | S0093 | 2,411,044 |
| SRR40032844 | PRJNA989260 | 2025-10-07 | S0092 | 2,123,184 |
| SRR40032848 | PRJNA989260 | 2025-10-08 | S0086 | 98,230 |
| SRR40032940 | PRJNA989260 | 2025-11-10 | S0094 | 2,084,750 |
| SRR40032946 | PRJNA989260 | 2025-11-13 | S0108 | 4,077,050 |
| SRR40032956 | PRJNA989260 | 2025-10-15 | S0095 | 3,324,024 |
| SRR40033114 | PRJNA989260 | 2025-12-17 | S0155 | 2,052,336 |
| SRR40033132 | PRJNA989260 | 2026-01-27 | S0153 | 18,843,838 |
| SRR40033239 | PRJNA989260 | 2026-01-22 | S0120 | 7,905,406 |
| SRR40033244 | PRJNA989260 | 2026-01-21 | S0154 | 6,706,218 |
| SRR40396823 | PRJNA1498729 | 2026-08-20 | S0090 | 4,363,326 |
| SRR40666459 | PRJNA1498729 | 2026-09-03 | S0120 | 8,343,230 |

`SRR40032848` yielded only 98,230 reads (~30x below its neighbors). It is the sole
run from S0086, the smallest sewershed (pop 3,816), and is probably a low-yield or
failed library rather than a usable sample.

S0120 appears under both PRJNA989260 and PRJNA1498729, so the two projects can be
compared at one sewershed across the lab handover.

**Metagenomic (RNA-Seq)** — `data/fastq_rnaseq/`, 103 MB.

| Run | Collected | Site | Pairs kept | Full run size |
|---|---|---|---|---|
| SRR38294894 | 2026-02-01 | CHI-A influent | 1,000,000 | 88.7 GB (44.4 + 44.3) |

Split `_1`/`_2` mates, 151 bp, verified: 0 malformed records, mate names identical.
Collected inside the Feb-Jun 2026 amplicon gap, so it covers a period nothing else
here reaches.

## Subsampling large runs

`data/subsample_rnaseq.sh ACCESSION [PAIRS]` pulls a prefix of a CASPER run from ENA
(which, unlike NCBI, serves these with `Accept-Ranges: bytes` and split mate files).
It streams and stops once enough records are decoded, so ~107 MB crossed the network
for 1M pairs out of an 88.7 GB run.

Two hazards the script encodes, both hit during development:

- **Never pass `--retry` to a `curl` writing to a pipe.** On a dropped connection curl
  restarts from byte 0 and concatenates the restart into the stream, so the
  decompressor emits misaligned records. ENA drops these connections routinely; this
  produced a file whose line count was divisible by 4 and which therefore passed a
  naive check, while 249 of 390,947 records were actually misaligned. The amplicon
  script's `--retry` is safe only because it writes with `-o file`, where curl
  truncates on retry.
- **A stream cut mid-record leaves a partial record.** The script validates each
  record (header `@`, third line `+`, equal sequence/quality lengths) and stops at the
  first malformed one, then pairs mates by read-name intersection rather than by
  position, since the two streams can end at different points.

### Subsample bias

This is a **head subsample, not a random one**. Reads in these files are ordered by
flowcell position: the 1M pairs kept span tiles 1101-1108 only. That is fine for a
pipeline smoke test, read-length/quality inspection, or presence/absence checks, but
it is not valid for quantitative composition estimates, since tile position correlates
with error and optical-duplicate rates. Proper random subsampling needs the whole file
plus `seqtk sample`.

## Caveats

- Filtering is on **collection_date**. The E-utilities date index is release date (`[PDAT]`); collection date is only in sample attributes, so the pipeline over-fetches on `[PDAT]` and filters locally. Since release date always follows collection date, no record collected in the window can be missed this way.
- **CASPER is ~1,000× larger per run** (~295 Gbp/run deep metagenomes) than the amplicon runs (~0.5–0.6 Gbp/run). Downloading all 206 CASPER runs is ~60 Tbp of reads — filter before pulling.
- Chicago attribution for NWSS records rests on `ww_surv_jurisdiction == "CI"`; NCBI carries no city field for these samples. A handful of `CI`-site records are labeled `IL`/`il` inconsistently in the same project.
