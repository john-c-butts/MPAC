#!/bin/bash
# Submit the full openTargets seqlet-filtered pipeline with SLURM dependencies.
#
# Stage 1: seqlet iteration (single job)     → ot_candidates_per_chrom/
# Stage 2: gnomAD filter (array 1-22)        → ot_filtered_per_chrom/
# Stage 3: merge + phenotype join (single)   → openTargets_GWAS_multiTF_seqletFiltered_with_controls.tsv
#
# Usage (run from multiTF_cluster/):
#   bash scripts/submit_openTargets_seqlet_filtered_pipeline.sh

set -euo pipefail

cd /pod/2/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/multiTF_cluster

mkdir -p logs results_final/ot_candidates_per_chrom results_final/ot_filtered_per_chrom

# ---------------------------------------------------------------------------
# Stage 1: build seqlet candidates (single job)
# ---------------------------------------------------------------------------
STAGE1_ID=$(sbatch --parsable \
    --job-name=ot_s1_candidates \
    --output=logs/ot_s1_candidates_%j.out \
    --error=logs/ot_s1_candidates_%j.err \
    --time=01:00:00 \
    --mem=64G \
    --cpus-per-task=1 \
    --partition=compute \
    --mail-user=john.butts@jax.org \
    --mail-type=FAIL \
    --wrap="source ~/.bashrc && conda activate block_calling && python scripts/openTargets_build_candidates.py")

echo "Stage 1 submitted: job ${STAGE1_ID}"

# ---------------------------------------------------------------------------
# Stage 2: per-chrom gnomAD filter (array, depends on Stage 1)
# ---------------------------------------------------------------------------
STAGE2_ID=$(sbatch --parsable \
    --job-name=ot_s2_gnomad \
    --output=logs/ot_s2_gnomad_%A_%a.out \
    --error=logs/ot_s2_gnomad_%A_%a.err \
    --time=02:00:00 \
    --mem=32G \
    --cpus-per-task=1 \
    --partition=compute \
    --array=1-22 \
    --dependency=afterok:${STAGE1_ID} \
    --mail-user=john.butts@jax.org \
    --mail-type=FAIL \
    --wrap='source ~/.bashrc && conda activate block_calling && CHROM="chr${SLURM_ARRAY_TASK_ID}" && python scripts/openTargets_gnomad_filter.py --chrom ${CHROM} --output results_final/ot_filtered_per_chrom/${CHROM}_ot_filtered.tsv')

echo "Stage 2 submitted: job array ${STAGE2_ID} (depends on ${STAGE1_ID})"

# ---------------------------------------------------------------------------
# Stage 3: merge + phenotype join (single job, depends on all of Stage 2)
# ---------------------------------------------------------------------------
STAGE3_ID=$(sbatch --parsable \
    --job-name=ot_s3_merge \
    --output=logs/ot_s3_merge_%j.out \
    --error=logs/ot_s3_merge_%j.err \
    --time=00:30:00 \
    --mem=16G \
    --cpus-per-task=1 \
    --partition=compute \
    --dependency=afterok:${STAGE2_ID} \
    --mail-user=john.butts@jax.org \
    --mail-type=END,FAIL \
    --wrap="source ~/.bashrc && conda activate block_calling && python scripts/merge_openTargets_seqlet_filtered.py")

echo "Stage 3 submitted: job ${STAGE3_ID} (depends on ${STAGE2_ID})"

echo ""
echo "Pipeline queued. Job chain:"
echo "  Stage 1 (candidates):  ${STAGE1_ID}"
echo "  Stage 2 (gnomAD, x22): ${STAGE2_ID}"
echo "  Stage 3 (merge):        ${STAGE3_ID}"
echo ""
echo "Output: results_final/openTargets_GWAS_multiTF_seqletFiltered_with_controls.tsv"
echo "Monitor with: squeue -u \$USER"
