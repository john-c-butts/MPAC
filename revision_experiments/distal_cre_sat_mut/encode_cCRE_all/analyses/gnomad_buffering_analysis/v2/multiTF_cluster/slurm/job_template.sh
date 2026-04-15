#!/bin/bash
#SBATCH --job-name=multiTF_CHROM
#SBATCH --output=../logs/multiTF_CHROM_%j.out
#SBATCH --error=../logs/multiTF_CHROM_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=1
#SBATCH --partition=compute

# Usage: sbatch job_template.sh chr22
# Or use submit_all.sh to submit all chromosomes

CHROM=$1

if [ -z "$CHROM" ]; then
    echo "Error: Chromosome not specified"
    echo "Usage: sbatch job_template.sh chr22"
    exit 1
fi

echo "Starting multi-TF analysis for ${CHROM}"
echo "Date: $(date)"
echo "Node: $(hostname)"

# Set up environment
cd /pod/2/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/multiTF_cluster

# Activate conda environment if needed (uncomment and modify as needed)
# source ~/.bashrc
conda init
conda activate block_calling

# Run the analysis
python scripts/multiTF_analysis.py \
    --chrom ${CHROM} \
    --output results_per_chrom/${CHROM}_multiTF.pkl \
    --max-overlap-bp 2

echo "Finished at: $(date)"
