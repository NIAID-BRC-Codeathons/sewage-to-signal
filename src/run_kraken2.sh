#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status
set -euo pipefail

# 1. Define Paths
DB_DIR="/codeathon/wastewater/databases/kraken2"
SAMPLE_DIR="./SRR35987665_out_manual"
THREADS=32

# 2. Define Input Files
R1="${SAMPLE_DIR}/SRR35987665_1.fastq.gz"
R2="${SAMPLE_DIR}/SRR35987665_2.fastq.gz"

# 3. Define Output Names
OUTPUT_FILE="${SAMPLE_DIR}/SRR35987665.kraken"
REPORT_FILE="${SAMPLE_DIR}/SRR35987665.report"

# 4. Run High-Performance Kraken 2 Command
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Starting Kraken 2 classification for SRR35987665..."

kraken2 \
  --db "${DB_DIR}" \
  --threads "${THREADS}" \
  --memory-mapping \
  --paired \
  --gzip-compressed \
  --output "${OUTPUT_FILE}" \
  --report "${REPORT_FILE}" \
  "${R1}" \
  "${R2}"

echo "[$(date +'%Y-%m-%d %H:%M:%S')] Kraken 2 processing complete!"
