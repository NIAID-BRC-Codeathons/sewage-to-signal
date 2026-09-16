#!/bin/bash
#SBATCH --job-name=sae-pipeline
#SBATCH --output=logs/sae-%j.out
#SBATCH --error=logs/sae-%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --gres=gpu:1            # s06_embed only; drop it to run CPU stages
# Adjust partition/account for your cluster:
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
module load apptainer 2>/dev/null || module load singularity 2>/dev/null || true

REPO=/path/to/codeathon
SAMPLE=${1:?usage: sbatch slurm_example.sbatch SAMPLE_ID}

# Keep the model cache on shared scratch, not $HOME — it is ~27 GB and every
# job reuses it.
export SAE_HF=/scratch/$USER/hf_cache
export SAE_WORK=/scratch/$USER/sae_work/$SAMPLE
export SAE_DATA=$REPO/data
export SAE_SIF=/shared/images/sae.sif
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

mkdir -p logs "$SAE_WORK" "$SAE_HF"

"$REPO"/container/run.sh pipeline \
    --fastq  "/data/fastq_rnaseq/${SAMPLE}_1.fastq.gz" \
    --fastq2 "/data/fastq_rnaseq/${SAMPLE}_2.fastq.gz" \
    --sample "$SAMPLE" \
    --batch-size 16 \
    --hmm /data/Pfam-A.hmm          # prefilter: the main cost lever
