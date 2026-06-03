#!/bin/bash
#SBATCH --job-name=dELSseqlets
#SBATCH --array=1-22
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/seqlets012326_chr%a_%j.out
#SBATCH --error=logs/seqlets012326_chr%a_%j.err
#SBATCH --mail-user=john.butts@jax.org
#SBATCH --mail-type=END,FAIL

# Create logs directory if it doesn't exist
mkdir -p logs

# Define input and output directories
INPUT_DIR="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/processed_data/seqlets_calling/seqlets_by_chrom"
OUTPUT_DIR="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/processed_data/seqlets_calling/seqlets_by_chrom/annotated_seqlets"
MPAC_DIR="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/mpac_preds"
BED_DIR="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/processed_data/seqlets_calling/seqlets_by_chrom/annotated_seqlets/bed_files/"

# Create output directory
mkdir -p ${OUTPUT_DIR}
# create bed directory
mkdir -p ${BED_DIR}

# Define input file pattern (adjust to match your naming convention)
RAW_TENSOR="${INPUT_DIR}/chr${SLURM_ARRAY_TASK_ID}_raw_tensors.pt"
SEQLET_TENSOR="${INPUT_DIR}/chr${SLURM_ARRAY_TASK_ID}_recursive_seqlets_tensors.pt"
ONEHOT_TENSOR="${INPUT_DIR}/chr${SLURM_ARRAY_TASK_ID}_oneHot_tensors.pt"
MOTIFS="/projects/tewhey-lab/buttsj/meme_suite_data/motif_databases/HOCOMOCO/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme"
MOTIF_THRESH="0.01"
MPAC_PREDS="${MPAC_DIR}/GRCh38-dELS-chr${SLURM_ARRAY_TASK_ID}-ALL-mpac-017.tsv.gz"

echo "Processing chromosome: chr${SLURM_ARRAY_TASK_ID}"
echo "Raw tensors: $RAW_TENSOR"
echo "Seqlet tensors: $SEQLET_TENSOR"
echo "OneHot tensors: $ONEHOT_TENSOR"
echo "Output directory: $OUTPUT_DIR"
echo "Motifs: $MOTIFS"

echo "SLURM Array Task ID: $SLURM_ARRAY_TASK_ID"

# Run the Python script
python tangermeme_seqlet_call_and_annotate_cli.py \
    --raw_tensor_input "$RAW_TENSOR" \
    --seqlet_formatted_input "$SEQLET_TENSOR" \
    --oneHot_tensor_input "$ONEHOT_TENSOR" \
    --chromosome "chr${SLURM_ARRAY_TASK_ID}" \
    --output_dir "$OUTPUT_DIR" \
    --motif_input "$MOTIFS" \
    --seqlet_threshold "$MOTIF_THRESH" \
    --return_bed_file \
    --sat_mut_preds "$MPAC_PREDS" \
    --bed_output_path "$BED_DIR"

echo "Completed annotating seqlets for chr${SLURM_ARRAY_TASK_ID}"
