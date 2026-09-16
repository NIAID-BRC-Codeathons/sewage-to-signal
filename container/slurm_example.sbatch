#!/bin/bash
#SBATCH --job-name=sae-pipeline
#SBATCH --output=sae-pipeline-%j.out
# Written to the submit directory on purpose. SLURM opens this file before
# the job body runs, so a subdirectory that does not exist yet yields a job
# that produces no output at all and gives no hint why.
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
# Site-specific — uncomment and adjust for your cluster. These stay commented
# because a gres name that the site does not define is rejected at submission
# ("Invalid generic resource (gres) specification"), not at run time, so an
# active guess here makes the script unsubmittable everywhere else.
##SBATCH --gres=gpu:1           # s06_embed wants a GPU; CPU stages do not
##SBATCH --partition=gpu
##SBATCH --account=your_account

# Example submission. Sizing notes, measured rather than guessed:
#   * s06_embed with ESMC-6B is memory-bandwidth-bound. On Apple M-series it
#     reached 1.03 seq/s at batch 1 and only 1.24 at batch 8; expect a CUDA GPU
#     to do far better, but benchmark before sizing a long run.
#   * Budget by residue, not protein count — attention is quadratic in length.
#     The same model ran 0.31 seq/s on ~1000 aa proteins vs 1.03 at 200 aa.
#   * Assembly is the memory-hungry stage; 64G is a starting point for a single
#     deep metagenome, not a guarantee.
set -euo pipefail
# Works either way:
#     sbatch container/slurm_example.sbatch SAMPLE      # on a cluster
#     bash   container/slurm_example.sbatch SAMPLE      # anywhere else
#     DRY_RUN=1 bash container/slurm_example.sbatch SAMPLE   # just print it
module load apptainer 2>/dev/null || module load singularity 2>/dev/null || true

SAMPLE=${1:?usage: sbatch slurm_example.sbatch SAMPLE_ID}
# Find the prelude explicitly. A failed `source` aborts the shell under
# `set -e` rather than falling through to `||`, and under sbatch the script
# runs from a spool copy, so its own directory is not the repo.
_sae_common=""
for _c in "$(dirname "${BASH_SOURCE[0]:-$0}")" \
          "${REPO:-}/container" \
          "${SLURM_SUBMIT_DIR:-}/container" \
          "$PWD/container" "$PWD"; do
  [ -n "$_c" ] && [ -f "$_c/slurm_common.sh" ] && { _sae_common="$_c/slurm_common.sh"; break; }
done
[ -n "$_sae_common" ] || {
  echo "error: slurm_common.sh not found; submit from a checkout or set REPO=" >&2
  exit 2
}
# shellcheck source=slurm_common.sh
. "$_sae_common"
sae_env
export SAE_WORK="$SAE_WORK/$SAMPLE"
mkdir -p "$SAE_WORK"
sae_banner

sae_run "$REPO"/container/run.sh pipeline \
    --fastq  "/data/fastq_rnaseq/${SAMPLE}_1.fastq.gz" \
    --fastq2 "/data/fastq_rnaseq/${SAMPLE}_2.fastq.gz" \
    --sample "$SAMPLE" \
    --batch-size 16 \
    --hmm /data/pfam/Pfam-A.hmm \
    --bit-cutoffs gathering
