#!/bin/bash
# Subsample a CASPER (PRJNA1247874) RNA-Seq metagenome without downloading ~90 GB.
#
# Full runs are ~44 GB per mate. We stream each mate from ENA and stop once enough
# reads are decoded: `head`-style early exit closes the pipe and curl ends the transfer.
#
# Two hazards this script exists to avoid:
#  1. NEVER pass --retry to a curl that writes to a PIPE. On a dropped connection curl
#     restarts from byte 0 and concatenates the restart into the stream, so the
#     decompressor sees [partial][beginning again] and silently emits misaligned
#     records. ENA drops these connections routinely, so this is not hypothetical.
#  2. A stream cut mid-record leaves a partial trailing record. We validate every
#     record (header starts @, 3rd line starts +) and stop at the first malformed one.
#
# This is a HEAD subsample, not a random one: reads are ordered by flowcell position,
# so a prefix is biased toward the first tiles. Fine for a pipeline smoke test or
# read-level inspection; NOT valid for quantitative composition, which needs the
# full file and `seqtk sample`.
set -u
ACC="${1:?usage: subsample_rnaseq.sh SRR_ACCESSION [PAIRS]}"
PAIRS="${2:-1000000}"
OUT="$(cd "$(dirname "$0")" && pwd)/fastq_rnaseq"
mkdir -p "$OUT"

if [ ${#ACC} -le 9 ]; then echo "ENA path layout differs for <=9 char accessions"; exit 2; fi
BASE="ftp.sra.ebi.ac.uk/vol1/fastq/${ACC:0:6}/$(printf '%03d' $((10#${ACC:9})))/$ACC"

# Emit only complete, well-formed 4-line records; exit at the first malformed one.
VALIDATE='
  { buf[NR%4==1?1:(NR%4==2?2:(NR%4==3?3:4))] = $0 }
  NR%4==0 {
    if (buf[1] !~ /^@/ || buf[3] !~ /^\+/ || length(buf[2]) != length(buf[4])) exit
    print buf[1]; print buf[2]; print buf[3]; print buf[4]
    if (++n >= want) exit
  }'

for m in 1 2; do
  url="https://$BASE/${ACC}_$m.fastq.gz"
  full=$(curl -sI "$url" | awk 'tolower($1)=="content-length:"{printf "%.1f", $2/1e9}')
  echo "[get ] ${ACC}_$m  streaming from ${full} GB, want ${PAIRS} reads"
  # No --retry: see hazard 1 above. A dropped connection simply ends the stream.
  curl -sS --no-buffer "$url" 2>/dev/null \
    | gzcat 2>/dev/null \
    | awk -v want="$PAIRS" "$VALIDATE" > "$OUT/${ACC}_$m.fastq"
  echo "       -> $(( $(wc -l < "$OUT/${ACC}_$m.fastq") / 4 )) valid records"
done

# Pair by read name intersection, not by position: the two mates can end at different
# points, and equal counts alone would not prove the same reads were kept.
echo "[pair] intersecting read names"
awk 'NR%4==1{print $1}' "$OUT/${ACC}_1.fastq" | sort > "$OUT/.n1"
awk 'NR%4==1{print $1}' "$OUT/${ACC}_2.fastq" | sort > "$OUT/.n2"
comm -12 "$OUT/.n1" "$OUT/.n2" > "$OUT/.both"
echo "       -> $(wc -l < "$OUT/.both" | tr -d ' ') names in both mates"
for m in 1 2; do
  awk 'NR==FNR{k[$1];next} FNR%4==1{keep=($1 in k)} keep' \
    "$OUT/.both" "$OUT/${ACC}_$m.fastq" > "$OUT/${ACC}_$m.paired"
  mv "$OUT/${ACC}_$m.paired" "$OUT/${ACC}_$m.fastq"
  gzip -f "$OUT/${ACC}_$m.fastq"
done
rm -f "$OUT/.n1" "$OUT/.n2" "$OUT/.both"

# Final check: equal counts, identical name order, all records well-formed.
c1=$(gzcat "$OUT/${ACC}_1.fastq.gz" | awk 'END{print NR/4}')
c2=$(gzcat "$OUT/${ACC}_2.fastq.gz" | awk 'END{print NR/4}')
h1=$(gzcat "$OUT/${ACC}_1.fastq.gz" | awk 'NR%4==1{print $1}' | md5)
h2=$(gzcat "$OUT/${ACC}_2.fastq.gz" | awk 'NR%4==1{print $1}' | md5)
bad=$(gzcat "$OUT/${ACC}_1.fastq.gz" | awk 'NR%4==1 && !/^@/{c++} NR%4==3 && !/^\+/{c++} END{print c+0}')
echo "[ok  ] $c1 pairs | names match: $([ "$h1" = "$h2" ] && echo yes || echo NO) | malformed records: $bad"
