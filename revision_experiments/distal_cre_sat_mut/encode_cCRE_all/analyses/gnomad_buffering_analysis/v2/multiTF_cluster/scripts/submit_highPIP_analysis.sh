#!/bin/bash
# Submit per-chromosome high-PIP shadow analysis jobs.
# Usage: bash scripts/submit_highPIP_analysis.sh

WORKDIR=/pod/2/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/multiTF_cluster

mkdir -p ${WORKDIR}/results_final/shadow_per_chrom
mkdir -p ${WORKDIR}/logs

echo "Submitting high-PIP shadow analysis jobs for all chromosomes..."

for i in {1..22}; do
    CHROM="chr${i}"
    echo "Submitting ${CHROM}..."
    sbatch \
        --job-name=shadow_${CHROM} \
        --output=${WORKDIR}/logs/shadow_${CHROM}_%j.out \
        --error=${WORKDIR}/logs/shadow_${CHROM}_%j.err \
        --time=04:00:00 \
        --mem=48G \
        --cpus-per-task=1 \
        --partition=compute \
        --wrap="
            cd ${WORKDIR}
            conda init
            conda activate block_calling
            python scripts/highPIP_shadow_analysis.py \
                --chrom ${CHROM} \
                --pip-threshold 0.9 \
                --output results_final/shadow_per_chrom/${CHROM}_shadow_variants.tsv
        "
done

echo ""
echo "All jobs submitted. Check status with: squeue -u \$USER"
echo "After completion, run: python scripts/merge_highPIP_results.py"
