#!/bin/bash
#SBATCH --job-name=ot_gnomad
#SBATCH --output=logs/ot_gnomad_filter_%A_%a.out
#SBATCH --error=logs/ot_gnomad_filter_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --partition=compute
#SBATCH --array=1-22
#SBATCH --mail-user=john.butts@jax.org
#SBATCH --mail-type=END,FAIL

# Stage 2: per-chrom gnomAD annotation, filtering, and control merging.
# Run AFTER submit_openTargets_build_candidates.sh completes.
#
# Filters applied per chrom:
#   - Lead emVars: kept as-is (in seqlet, PIP >= 0.9)
#   - Shadow emVars: exceeds_lead_skew == True AND has gnomAD v4 AF
#   - Controls: outside seqlets, not emVars, has gnomAD v4 AF
#
# After all array jobs complete, run Stage 3:
#   python scripts/merge_openTargets_seqlet_filtered.py

cd /pod/2/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/multiTF_cluster

mkdir -p logs results_final/ot_filtered_per_chrom

CHROM="chr${SLURM_ARRAY_TASK_ID}"

echo "Stage 2: gnomAD filter for ${CHROM}"
echo "Date: $(date)"
echo "Node: $(hostname)"

source ~/.bashrc
conda activate block_calling

python scripts/openTargets_gnomad_filter.py \
    --chrom ${CHROM} \
    --output results_final/ot_filtered_per_chrom/${CHROM}_ot_filtered.tsv

echo "Finished ${CHROM} at: $(date)"
