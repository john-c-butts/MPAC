#!/bin/bash
# Submit all chromosome jobs for multi-TF analysis
# Usage: ./submit_all.sh

cd /pod/2/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/multiTF_cluster/slurm

echo "Submitting multi-TF analysis jobs for all chromosomes..."

for i in {1..22}; do
    CHROM="chr${i}"
    echo "Submitting ${CHROM}..."
    sbatch --job-name=multiTF_${CHROM} \
           --output=../logs/multiTF_${CHROM}_%j.out \
           --error=../logs/multiTF_${CHROM}_%j.err \
           job_template.sh ${CHROM}
done

echo ""
echo "All jobs submitted. Check status with: squeue -u $USER"
echo "After completion, run: python scripts/merge_results.py"
