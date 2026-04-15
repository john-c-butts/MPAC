#!/bin/bash
#SBATCH --job-name=dELSseqlets
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=john.butts@jax.org
#SBATCH --mail-type=END,FAIL

SCRIPT_DIR="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/alt_allele_predictions/scripts"

# Define input and output directories
INPUT_DIR="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/alt_allele_predictions/processed_data/seqlets"
OUTPUT_DIR="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/alt_allele_predictions/processed_data/seqlets/annotated_seqlets"
BED_DIR="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/alt_allele_predictions/processed_data/seqlets/bed_files/"

# Create output directories
mkdir -p ${OUTPUT_DIR}
mkdir -p ${BED_DIR}

# Define input files
RAW_TENSOR="${INPUT_DIR}/dELS_10K_sample_raw_tensors.pt"
SEQLET_TENSOR="${INPUT_DIR}/dELS_10K_sample_recursive_seqlets_tensors.pt"
ONEHOT_TENSOR="${INPUT_DIR}/dELS_10K_sample_oneHot_tensors.pt"  # NOTE: must be produced upstream
MOTIFS="/projects/tewhey-lab/buttsj/meme_suite_data/motif_databases/HOCOMOCO/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme"
MOTIF_THRESH="0.01"
MPAC_PREDS="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/alt_allele_predictions/processed_data/dELS_10K_satmut_preds_sample.tsv.gz"

echo "Raw tensors: $RAW_TENSOR"
echo "Seqlet tensors: $SEQLET_TENSOR"
echo "OneHot tensors: $ONEHOT_TENSOR"
echo "Output directory: $OUTPUT_DIR"
echo "Motifs: $MOTIFS"

# Run the Python script
python ${SCRIPT_DIR}/tangermeme_seqlet_call_and_annotate_cli.py \
    --raw_tensor_input "$RAW_TENSOR" \
    --seqlet_formatted_input "$SEQLET_TENSOR" \
    --oneHot_tensor_input "$ONEHOT_TENSOR" \
    --chromosome "dELS_100K_sample" \
    --output_dir "$OUTPUT_DIR" \
    --motif_input "$MOTIFS" \
    --seqlet_threshold "$MOTIF_THRESH" \
    --return_bed_file \
    --sat_mut_preds "$MPAC_PREDS" \
    --bed_output_path "$BED_DIR"

echo "Completed annotating seqlets"
