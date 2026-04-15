#!/bin/bash
#SBATCH --job-name=vep_reg_debug
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --time=72:00:00
#SBATCH --mem=32G
#SBATCH --mail-user=john.butts@jax.org
#SBATCH --mail-type=END,FAIL

module load singularity

VEP_SIF="/projects/tewhey-lab/buttsj/ensembl_vep/ensembl-vep_110.sif"
VEP_DIR="/projects/tewhey-lab/buttsj/vep_data/"
INDIR="../processed_data/vep_vcfs"
OUTDIR="../processed_data/vep_debug_output"

mkdir -p ${OUTDIR}

VCFS=(
    "cosmic_basic_exons_with_pads"
)

FIELDS="Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,BIOTYPE,CANONICAL,NEAREST,DISTANCE,IMPACT,FLAGS"

# Common VEP flags shared across both runs
COMMON_FLAGS="--dir ${VEP_DIR} \
    --cache --offline \
    --format vcf \
    --tab \
    --force_overwrite \
    --nearest transcript \
    --distance 1000 \
    --fork 8 \
    --canonical \
    --flag_pick \
    --fields ${FIELDS}"

for VCF_BASE in "${VCFS[@]}"; do

    INPUT="${INDIR}/${VCF_BASE}.vcf"

    echo "============================================"
    echo "Processing: ${VCF_BASE}"
    echo "Start time: $(date)"
    echo "============================================"

    # ---- Run 1: WITH --regulatory ----
    echo "[$(date)] Starting regulatory run for ${VCF_BASE}..."

    singularity exec ${VEP_SIF} \
        vep ${COMMON_FLAGS} \
            --regulatory \
            --input_file ${INPUT} \
            --output_file ${OUTDIR}/${VCF_BASE}_with_regulatory.tsv

    echo "[$(date)] Finished regulatory run for ${VCF_BASE}"

    # ---- Run 2: WITHOUT --regulatory ----
    echo "[$(date)] Starting no-regulatory run for ${VCF_BASE}..."

    singularity exec ${VEP_SIF} \
        vep ${COMMON_FLAGS} \
            --input_file ${INPUT} \
            --output_file ${OUTDIR}/${VCF_BASE}_no_regulatory.tsv

    echo "[$(date)] Finished no-regulatory run for ${VCF_BASE}"

done

echo "============================================"
echo "All runs complete: $(date)"
echo "============================================"
