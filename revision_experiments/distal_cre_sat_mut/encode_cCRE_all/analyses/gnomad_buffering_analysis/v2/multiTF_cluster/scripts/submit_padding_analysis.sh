#!/bin/bash
#SBATCH --job-name=pad_emVars
#SBATCH --output=logs/padding_emVars_%A_%a.out
#SBATCH --error=logs/padding_emVars_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --partition=compute
#SBATCH --array=1-22
#SBATCH --mail-user=john.butts@jax.org
#SBATCH --mail-type=END,FAIL

# Runs padding_emVars_analysis.py per chromosome in parallel.
# After all jobs complete, run: python scripts/merge_padding_results.py

cd /pod/2/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/multiTF_cluster

mkdir -p logs results_final/padding_per_chrom

CHROM="chr${SLURM_ARRAY_TASK_ID}"

echo "Starting padding emVars analysis for ${CHROM}"
echo "Date: $(date)"
echo "Node: $(hostname)"

# Activate conda
source ~/.bashrc
conda activate block_calling

python scripts/padding_emVars_analysis.py \
    --chrom ${CHROM} \
    --pad 4 \
    --output results_final/padding_per_chrom/${CHROM}_padding_emVars.tsv

echo "Finished ${CHROM} at: $(date)"
