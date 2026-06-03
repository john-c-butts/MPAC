#!/bin/bash
#SBATCH --job-name=satmut_tensors
#SBATCH --array=1-2
#SBATCH --time=08:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/tMeme012326%a_%j.out
#SBATCH --error=logs/tMeme012326%a_%j.err
#SBATCH --mail-user=john.butts@jax.org
#SBATCH --mail-type=END,FAIL

# Create logs directory if it doesn't exist
mkdir -p logs

# Define input and output directories
INPUT_DIR="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/mpac_preds"
OUTPUT_DIR="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/processed_data/seqlets_calling/seqlets_by_chrom"

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Define input file pattern (adjust to match your naming convention)
INPUT_FILE="${INPUT_DIR}/GRCh38-dELS-chr${SLURM_ARRAY_TASK_ID}-ALL-mpac-017.tsv.gz"

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file $INPUT_FILE not found!"
    exit 1
fi

echo "Processing chromosome: chr${SLURM_ARRAY_TASK_ID}"
echo "Input file: $INPUT_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "SLURM Array Task ID: $SLURM_ARRAY_TASK_ID"

# Run the Python script
python call_tangermeme_rawCorrected.py \
    --input_file "$INPUT_FILE" \
    --output_dir "$OUTPUT_DIR" \
    --chromosome "chr${SLURM_ARRAY_TASK_ID}" \
    --chunksize 1000000 \
    --save_raw_tensors \
    --verify_order

echo "Completed processing for chr${SLURM_ARRAY_TASK_ID}"
