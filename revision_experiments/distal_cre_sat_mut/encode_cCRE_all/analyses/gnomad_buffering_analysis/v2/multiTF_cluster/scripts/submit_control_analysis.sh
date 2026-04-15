#!/bin/bash
#SBATCH --job-name=ctrl_vars
#SBATCH --output=logs/control_variants_%A_%a.out
#SBATCH --error=logs/control_variants_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --partition=compute
#SBATCH --array=1-22
#SBATCH --mail-user=john.butts@jax.org
#SBATCH --mail-type=END,FAIL

# Runs multiTF_control_analysis.py per chromosome in parallel.
# After all jobs complete, run: python scripts/merge_control_results.py

cd /pod/2/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/multiTF_cluster

mkdir -p logs results_final/control_per_chrom

CHROM="chr${SLURM_ARRAY_TASK_ID}"

echo "Starting control variant analysis for ${CHROM}"
echo "Date: $(date)"
echo "Node: $(hostname)"

# Activate conda
source ~/.bashrc
conda activate block_calling

python scripts/multiTF_control_analysis.py \
    --chrom ${CHROM} \
    --pad 4 \
    --enhancer-source GTEx_MultiTF_allEmVars_HighPIP_gnomAD_v4_annotated_4bp_pad.tsv \
    --enhancer-source openTargets_GWAS_multiTF_allEmVars_highPIP_gnomAD_v4_annotated_4bp_pad.tsv \
    --output results_final/control_per_chrom/${CHROM}_control_variants.tsv

echo "Finished ${CHROM} at: $(date)"
