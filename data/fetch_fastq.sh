#!/bin/bash
# Download Chicago wastewater FASTQ from NCBI SRA.
#
# Why NCBI and not ENA: ENA registers these recent accessions but has NOT staged
# their FASTQ files (fastq_ftp is empty for all 381 of our in-window runs).
#
# Endpoint caveats:
#   * Returns ONE concatenated gzip, not split _1/_2 mate files.
#     For proper mate separation install sra-tools and use `fasterq-dump --split-3`.
#   * No Accept-Ranges and no Content-Length -> byte-range resume is impossible
#     (`curl -C -` fails with error 56). An interrupted run restarts that accession.
#
# Integrity: download to <acc>.part, verify, then atomically rename. A partial file
# therefore never occupies the final path, and completion is recorded in manifest.tsv
# so re-runs skip finished work without re-probing the server.
set -u
OUT="$(cd "$(dirname "$0")" && pwd)/fastq"
MAN="$OUT/manifest.tsv"
URL=https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq
mkdir -p "$OUT"
[ -f "$MAN" ] || printf 'run_accession\tbytes\treads\n' > "$MAN"
RUNS="${@:-SRR40033132 SRR40033244 SRR40033239 SRR40666459}"

for acc in $RUNS; do
  f="$OUT/$acc.fastq.gz"; p="$OUT/$acc.part"
  if grep -q "^$acc	" "$MAN" && [ -f "$f" ] \
     && [ "$(wc -c < "$f" | tr -d ' ')" = "$(awk -F'\t' -v a=$acc '$1==a{print $2}' "$MAN")" ]; then
    echo "[skip] $acc verified in manifest"; continue
  fi
  echo "[get ] $acc"
  rm -f "$p"
  if ! curl -sSf --retry 5 --retry-delay 5 --retry-all-errors "$URL?acc=$acc" -o "$p"; then
    echo "[FAIL] $acc download error"; rm -f "$p"; continue
  fi
  # A truncated gzip fails -t; a complete FASTQ has a record count divisible by 4.
  if ! gzip -t "$p" 2>/dev/null; then
    echo "[FAIL] $acc corrupt gzip"; rm -f "$p"; continue
  fi
  lines=$(gzcat "$p" | wc -l | tr -d ' ')
  if [ $((lines % 4)) -ne 0 ] || [ "$lines" -eq 0 ]; then
    echo "[FAIL] $acc truncated FASTQ ($lines lines)"; rm -f "$p"; continue
  fi
  mv "$p" "$f"
  bytes=$(wc -c < "$f" | tr -d ' ')
  grep -v "^$acc	" "$MAN" > "$MAN.tmp" && mv "$MAN.tmp" "$MAN"
  printf '%s\t%s\t%s\n' "$acc" "$bytes" "$((lines / 4))" >> "$MAN"
  echo "[ok  ] $acc  $(du -h "$f" | cut -f1)  $((lines / 4)) reads"
done
