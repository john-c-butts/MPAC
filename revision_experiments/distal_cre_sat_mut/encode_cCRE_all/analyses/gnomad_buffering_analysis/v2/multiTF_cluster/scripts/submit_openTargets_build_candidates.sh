#!/bin/bash
#SBATCH --job-name=ot_candidates
#SBATCH --output=logs/ot_build_candidates_%j.out
#SBATCH --error=logs/ot_build_candidates_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --partition=compute
#SBATCH --mail-user=john.butts@jax.org
#SBATCH --mail-type=END,FAIL

# Stage 1: seqlet iteration across all chromosomes (fast, ~5-10 min).
# Outputs per-chrom candidate TSVs to results_final/ot_candidates_per_chrom/
# and phenotype info to results_final/openTargets_phenotype_info.tsv.
#
# After this completes, submit Stage 2:
#   sbatch scripts/submit_openTargets_gnomad_filter.sh

cd /pod/2/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/multiTF_cluster

mkdir -p logs results_final/ot_candidates_per_chrom

echo "Stage 1: building openTargets seqlet candidates"
echo "Date: $(date)"
echo "Node: $(hostname)"

source ~/.bashrc
conda activate block_calling

python scripts/openTargets_build_candidates.py

echo "Stage 1 complete: $(date)"
echo "Next: sbatch scripts/submit_openTargets_gnomad_filter.sh"
