#!/bin/bash
#SBATCH --job-name=clinVEP
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --time=72:00:00
#SBATCH --mem=32G
#SBATCH --mail-user=john.butts@jax.org
#SBATCH --mail-type=END,FAIL

module load singularity

VEP_SIF="/projects/tewhey-lab/buttsj/ensembl_vep/ensembl-vep_110.sif"
VEP_DIR="/projects/tewhey-lab/buttsj/vep_data/"
INDIR="../processed_data/chrom_vcfs"
OUTDIR="../processed_data/"

mkdir -p ${OUTDIR}

VCFS=(
    "all_exon_filtered_clinvar_preds"
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
            --no_trim_ids \
            --input_file ${INPUT} \
            --output_file ${OUTDIR}/${VCF_BASE}_with_regulatory.tsv

    echo "[$(date)] Finished regulatory run for ${VCF_BASE}"
    
done

echo "============================================"
echo "All runs complete: $(date)"
echo "============================================"
