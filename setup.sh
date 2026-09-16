#!/usr/bin/env bash
# One entry point for a fresh checkout.
#
#   ./setup.sh              # status report — tells you what is missing and what to run
#   ./setup.sh preflight    # verify required host commands only
#   ./setup.sh env          # create sae/.venv and install pinned deps
#   ./setup.sh features     # SAE feature descriptions        (~50 MB)
#   ./setup.sh model-small  # ESMC-300M + SAE                 (~1.7 GB)
#   ./setup.sh model-6b     # ESMC-6B + SAE                   (~25 GB)
#   ./setup.sh reads        # SARS-CoV-2 amplicon FASTQ       (~2.4 GB)
#   ./setup.sh rnaseq       # CASPER metagenome subsample     (~100 MB)
#   ./setup.sh atlas        # ESM Atlas cluster tables        (~27 GB)
#   ./setup.sh quickstart   # env + features + model-small + rnaseq  (~2 GB)
#   ./setup.sh all          # everything except atlas
#
# Every target is idempotent: it checks for a complete result and skips.
# Nothing here downloads 25 GB unless you name it.
#
# Requires: uv, curl, awk, gzip + gzcat/zcat, md5 or md5sum.
# `./setup.sh status` verifies all of them. There is deliberately no pip
# fallback — see require_uv() below.

set -uo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
VENV="$ROOT/sae/.venv"
PY="$VENV/bin/python"
PYVER="3.12"   # the interpreter this project is pinned to; uv provisions it

# --- formatting -------------------------------------------------------------
if [ -t 1 ]; then B=$'\033[1m'; G=$'\033[32m'; Y=$'\033[33m'; R=$'\033[31m'; D=$'\033[2m'; N=$'\033[0m'
else B=""; G=""; Y=""; R=""; D=""; N=""; fi
say()  { printf '%s\n' "$*"; }
ok()   { printf '  %sok%s    %s\n' "$G" "$N" "$*"; }
miss() { printf '  %s--%s    %s\n' "$Y" "$N" "$*"; }
bad()  { printf '  %sfail%s  %s\n' "$R" "$N" "$*"; }
head2(){ printf '\n%s%s%s\n' "$B" "$*" "$N"; }

# Is a file present and at least N bytes? Catches truncated downloads.
big_enough() { [ -f "$1" ] && [ "$(wc -c <"$1" | tr -d ' ')" -ge "$2" ]; }

human() { awk -v b="$1" 'BEGIN{
  split("B KB MB GB TB",u," "); i=1; while(b>=1024 && i<5){b/=1024;i++}
  printf (i==1?"%d %s":"%.1f %s"), b, u[i]}'; }

# --- venv -------------------------------------------------------------------
need_venv() {
  [ -x "$PY" ] || { bad "no venv — run: ./setup.sh env"; exit 1; }
}

# uv is required, with no pip fallback. uv provisions Python 3.12 itself, so
# the environment never depends on whatever python3 the host happens to ship
# (this machine's is 3.14). A pip fallback would build a venv against that
# interpreter instead and report success, which is worse than stopping.
require_uv() {
  command -v uv >/dev/null 2>&1 && return 0
  cat >&2 <<'EOT'

error: uv not found on PATH, and it is required.

  install one of:
      curl -LsSf https://astral.sh/uv/install.sh | sh
      brew install uv
      pipx install uv

  then re-run:
      ./setup.sh env

  There is no pip fallback by design: uv pins the interpreter (Python 3.12)
  this project is validated against. Falling back to the host python3 would
  silently build a different environment and call it success.

EOT
  exit 127
}

t_env() {
  head2 "Environment"
  require_uv

  if [ -x "$PY" ]; then
    local v
    v="$("$PY" -c 'import sys;print("%d.%d"%sys.version_info[:2])' 2>/dev/null || echo "?")"
    if [ "$v" != "$PYVER" ]; then
      bad "venv is Python $v, expected $PYVER"
      say "        remove it and re-run:  rm -rf $VENV && ./setup.sh env"
      exit 1
    fi
    ok "venv present (Python $v)"
  else
    uv venv --python "$PYVER" "$VENV" >/dev/null \
      || { bad "uv venv failed (could not provision Python $PYVER)"; exit 1; }
    ok "created venv (Python $PYVER)"
  fi

  uv pip install --python "$PY" -r "$ROOT/requirements.txt" -q \
    || { bad "dependency install failed — see output above"; exit 1; }
  ok "dependencies installed from requirements.txt"
}

# --- Hugging Face assets ----------------------------------------------------
# These land in ~/.cache/huggingface and are shared across checkouts, so a
# second clone costs nothing.
hf_snapshot() {
  need_venv
  "$PY" - "$@" <<'PYEOF'
import sys
from huggingface_hub import snapshot_download
repo, kind = sys.argv[1], (sys.argv[2] if len(sys.argv) > 2 else "model")
p = snapshot_download(repo, repo_type=kind, max_workers=8)
print(f"  -> {p}")
PYEOF
}

hf_present() {
  need_venv
  "$PY" - "$1" "${2:-model}" <<'PYEOF' 2>/dev/null
import sys
from huggingface_hub import snapshot_download
try:
    snapshot_download(sys.argv[1], repo_type=sys.argv[2], local_files_only=True)
except Exception:
    sys.exit(1)
PYEOF
}

t_features() {
  head2 "SAE feature descriptions (~50 MB)"
  if hf_present biohub/ESMC-SAE-Features dataset; then ok "cached"
  else hf_snapshot biohub/ESMC-SAE-Features dataset && ok "downloaded"; fi
}

t_model_small() {
  head2 "ESMC-300M + SAE (~1.7 GB)"
  for r in biohub/ESMC-300M biohub/ESMC-300M-sae-layer23-k64-codebook16384; do
    if hf_present "$r"; then ok "$r cached"; else hf_snapshot "$r" && ok "$r"; fi
  done
}

t_model_6b() {
  head2 "ESMC-6B + SAE (~25 GB)"
  say "  ${D}The only variant with published feature descriptions.${N}"
  for r in biohub/ESMC-6B-sae-layer60-k64-codebook16384 biohub/ESMC-6B; do
    if hf_present "$r"; then ok "$r cached"; else hf_snapshot "$r" && ok "$r"; fi
  done
}

# --- sequencing reads -------------------------------------------------------
t_reads() {
  head2 "SARS-CoV-2 amplicon FASTQ (~2.4 GB)"
  bash "$ROOT/data/fetch_fastq.sh"
}

t_rnaseq() {
  head2 "CASPER metagenome subsample (~100 MB)"
  if big_enough "$ROOT/data/fastq_rnaseq/SRR38294894_1.fastq.gz" 10000000 \
  && big_enough "$ROOT/data/fastq_rnaseq/SRR38294894_2.fastq.gz" 10000000; then
    ok "SRR38294894 present"
  else
    bash "$ROOT/data/subsample_rnaseq.sh" SRR38294894 1000000
  fi
}

# --- ESM Atlas tables -------------------------------------------------------
# Public, unauthenticated bucket. Object names and exact byte sizes are pinned
# in data/atlas_manifest.tsv, so verification is exact rather than a size
# heuristic and a truncated transfer can never be mistaken for a finished one.
S3_HOST="https://esm-protein-atlas.s3.amazonaws.com"
S3_BUCKET="s3://esm-protein-atlas"
MANIFEST="$ROOT/data/atlas_manifest.tsv"

# Print "path<TAB>bytes<TAB>key" for entries whose local size != manifest size.
atlas_pending() {
  awk -F'\t' 'NR>1 && $1!=""' "$MANIFEST" | while IFS=$'\t' read -r path bytes key; do
    local actual=""
    [ -f "$ROOT/$path" ] && actual="$(wc -c <"$ROOT/$path" | tr -d ' ')"
    [ "$actual" = "$bytes" ] || printf '%s\t%s\t%s\n' "$path" "$bytes" "$key"
  done
}

t_atlas() {
  head2 "ESM Atlas cluster tables (~28.6 GB)"
  [ -f "$MANIFEST" ] || { bad "missing $MANIFEST"; exit 1; }

  local total pending
  total=$(awk -F'\t' 'NR>1 && $1!=""' "$MANIFEST" | wc -l | tr -d ' ')
  pending="$(atlas_pending)"
  if [ -z "$pending" ]; then
    ok "all $total files present and byte-exact"
    return 0
  fi

  local n; n=$(printf '%s\n' "$pending" | wc -l | tr -d ' ')
  local transport
  if command -v aws >/dev/null 2>&1; then transport="aws s3 cp"; else transport="curl"; fi
  say "  ${D}$((total - n))/$total already complete; fetching $n via $transport${N}"

  printf '%s\n' "$pending" | while IFS=$'\t' read -r path bytes key; do
    local dest="$ROOT/$path"
    mkdir -p "$(dirname "$dest")"
    say "  get  $path  ($(human "$bytes"))"
    if [ "$transport" = "aws s3 cp" ]; then
      aws s3 cp --no-sign-request --only-show-errors "$S3_BUCKET/$key" "$dest" \
        || { bad "download failed: $path"; exit 1; }
    else
      # .part then rename, so an interrupted transfer never occupies the
      # final path and looks complete to the next run.
      curl -fL -sS --retry 3 --retry-delay 2 -o "$dest.part" "$S3_HOST/$key" \
        || { bad "download failed: $path"; rm -f "$dest.part"; exit 1; }
      mv "$dest.part" "$dest"
    fi
    local actual=""; [ -f "$dest" ] && actual="$(wc -c <"$dest" | tr -d ' ')"
    if [ "$actual" != "$bytes" ]; then
      bad "size mismatch for $path: got $actual, expected $bytes"
      exit 1
    fi
    ok "$path verified"
  done

  if [ -z "$(atlas_pending)" ]; then ok "all $total files present and byte-exact"
  else bad "some files still incomplete — re-run ./setup.sh atlas"; return 1; fi
}

# --- status -----------------------------------------------------------------
# Host commands the scripts in this repo actually invoke.
t_preflight() {
  head2 "Required host commands"
  local fail=0
  command -v uv   >/dev/null && ok "uv      $(uv --version 2>&1 | head -1)" \
    || { bad "uv      REQUIRED — https://astral.sh/uv"; fail=1; }
  command -v curl >/dev/null && ok "curl    all data fetching" \
    || { bad "curl    REQUIRED"; fail=1; }
  command -v awk  >/dev/null && ok "awk     read validation and pairing" \
    || { bad "awk     REQUIRED"; fail=1; }
  command -v gzip >/dev/null && ok "gzip    subsample compression" \
    || { bad "gzip    REQUIRED"; fail=1; }
  if command -v gzcat >/dev/null || command -v zcat >/dev/null; then
    ok "gzcat   decompression (zcat also accepted)"
  else
    bad "gzcat   REQUIRED (or zcat)"; fail=1
  fi
  if command -v md5 >/dev/null || command -v md5sum >/dev/null; then
    ok "md5     mate-order verification (md5sum also accepted)"
  else
    bad "md5     REQUIRED (or md5sum)"; fail=1
  fi
  command -v aws >/dev/null && ok "aws     optional — faster/resumable atlas download" \
    || miss "aws     optional — atlas falls back to curl (same bytes, verified)"
  [ "$fail" = 0 ] || say "\n  ${Y}Install the items marked REQUIRED before continuing.${N}"
  return "$fail"
}

t_status() {
  t_preflight || true
  head2 "Checkout status"
  [ -x "$PY" ] && ok "venv                $("$PY" -V 2>&1)" || miss "venv                ./setup.sh env"

  if [ -x "$PY" ]; then
    hf_present biohub/ESMC-SAE-Features dataset \
      && ok "feature table       cached" || miss "feature table       ./setup.sh features"
    hf_present biohub/ESMC-300M \
      && ok "ESMC-300M           cached" || miss "ESMC-300M           ./setup.sh model-small"
    hf_present biohub/ESMC-6B \
      && ok "ESMC-6B             cached" || miss "ESMC-6B             ./setup.sh model-6b"
  fi

  local n; n=$(ls "$ROOT"/data/fastq/*.fastq.gz 2>/dev/null | wc -l | tr -d ' ')
  [ "$n" -gt 0 ] && ok "amplicon FASTQ      $n run(s)" || miss "amplicon FASTQ      ./setup.sh reads"

  big_enough "$ROOT/data/fastq_rnaseq/SRR38294894_1.fastq.gz" 10000000 \
    && ok "CASPER subsample    present" || miss "CASPER subsample    ./setup.sh rnaseq"

  if [ -f "$MANIFEST" ]; then
    local atotal apend
    atotal=$(awk -F'\t' 'NR>1 && $1!=""' "$MANIFEST" | wc -l | tr -d ' ')
    apend=$(atlas_pending | wc -l | tr -d ' ')
    if [ "$apend" = 0 ]; then ok "atlas tables        $atotal/$atotal byte-exact"
    else miss "atlas tables        $((atotal - apend))/$atotal — ./setup.sh atlas"; fi
  fi

  head2 "External binaries (not installable from pip)"
  for t in megahit spades.py mmseqs fastp; do
    command -v "$t" >/dev/null && ok "$t" || miss "$t  — needed by $( [ "$t" = mmseqs ] && echo s04_derep || { [ "$t" = fastp ] && echo 's01_qc (optional, 30x faster)' || echo s02_assemble; })"
  done
  say ""
  say "  ${D}conda install -c bioconda megahit mmseqs2 fastp${N}"
  say "  ${D}(this machine needs an interactive 'conda tos accept' first)${N}"
  head2 "Next"
  say "  ./setup.sh quickstart     # ~2 GB, enough to run the pipeline end to end"
  say "  cat DATA.md               # provenance and sizes for every item"
}

case "${1:-status}" in
  env)         t_env ;;
  features)    t_features ;;
  model-small) t_model_small ;;
  model-6b)    t_model_6b ;;
  reads)       t_reads ;;
  rnaseq)      t_rnaseq ;;
  atlas)       t_atlas ;;
  quickstart)  t_env; t_features; t_model_small; t_rnaseq; t_status ;;
  all)         t_env; t_features; t_model_small; t_model_6b; t_reads; t_rnaseq; t_status ;;
  preflight)   t_preflight ;;
  status)      t_status ;;
  *) say "unknown target: $1"; sed -n '2,20p' "$0"; exit 2 ;;
esac
