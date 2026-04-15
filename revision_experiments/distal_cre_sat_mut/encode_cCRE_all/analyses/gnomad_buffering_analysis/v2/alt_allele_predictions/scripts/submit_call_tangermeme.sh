#!/bin/bash
#SBATCH --job-name=satmut_tensors
#SBATCH --time=08:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=john.butts@jax.org
#SBATCH --mail-type=END,FAIL

SCRIPT_DIR="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/alt_allele_predictions/scripts"

# Define input file and output directory
INPUT_FILE="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/alt_allele_predictions/processed_data/dELS_10K_satmut_preds_sample_with_EH38E3356274_chr2_86221285_A_G_b38.tsv.gz"
OUTPUT_DIR="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/alt_allele_predictions/processed_data/seqlets"

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file $INPUT_FILE not found!"
    exit 1
fi

echo "Input file: $INPUT_FILE"
echo "Output directory: $OUTPUT_DIR"

# Run the Python script
python ${SCRIPT_DIR}/call_tangermeme_rawCorrected.py \
    --input_file "$INPUT_FILE" \
    --output_dir "$OUTPUT_DIR" \
    --chromosome "dELS_10K_sample" \
    --chunksize 1000000 \
    --save_raw_tensors \
    --verify_order

echo "Completed processing"
