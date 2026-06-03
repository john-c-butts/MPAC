#!/bin/bash
# Submit the highPIP eQTL × MPAC pipeline.
#
# Stage 1: SLURM array (chr1–22), one task per chromosome.
# Stage 2: Merge job, runs only after all array tasks succeed.
#
# Cluster: Sumner (JAX)
#   partition=compute, qos=batch
#   Max wall: 3 days | Max mem/node: ~738 GB | Max jobs: 900
#
# Usage:
#   bash submit_highPIP_pipeline.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"

SUSIE_DIR="${BASE_DIR}/raw_data/gtex_v10_fineMapping/SuSiE_fineMapped"
SIGNIF_DIR="${BASE_DIR}/raw_data/gtex_v10_fineMapping/GTEx_Analysis_v10_eQTL_updated"
MPAC_DIR="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/mpac_preds"
PER_CHROM_DIR="${BASE_DIR}/processed_data/per_chrom"
OUT_FILE="${BASE_DIR}/processed_data/highPIP_eQTL_mpac_all.tsv.gz"
LOG_DIR="${BASE_DIR}/processed_data/logs"

mkdir -p "$PER_CHROM_DIR" "$LOG_DIR"

# ---------------------------------------------------------------------------
# Stage 1: per-chromosome array
# ---------------------------------------------------------------------------
ARRAY_JOB_ID=$(sbatch --parsable \
    --job-name=highPIP_per_chrom \
    --partition=compute \
    --qos=batch \
    --array=1-22 \
    --time=4:00:00 \
    --mem=64G \
    --cpus-per-task=1 \
    --output="${LOG_DIR}/per_chrom_%a.out" \
    --error="${LOG_DIR}/per_chrom_%a.err" \
    --wrap="
        CHROM=\"chr\${SLURM_ARRAY_TASK_ID}\"
        python ${SCRIPT_DIR}/highPIP_eQTL_mpac_per_chrom.py \
            --chrom \"\${CHROM}\" \
            --susie_dir  \"${SUSIE_DIR}\" \
            --signif_dir \"${SIGNIF_DIR}\" \
            --mpac_dir   \"${MPAC_DIR}\" \
            --out_dir    \"${PER_CHROM_DIR}\"
    ")

echo "Submitted per-chrom array: job ${ARRAY_JOB_ID} (tasks 1–22)"

# ---------------------------------------------------------------------------
# Stage 2: merge — runs only after all array tasks complete successfully
# ---------------------------------------------------------------------------
MERGE_JOB_ID=$(sbatch --parsable \
    --job-name=highPIP_merge \
    --partition=compute \
    --qos=batch \
    --dependency=afterok:${ARRAY_JOB_ID} \
    --time=1:00:00 \
    --mem=16G \
    --cpus-per-task=1 \
    --output="${LOG_DIR}/merge.out" \
    --error="${LOG_DIR}/merge.err" \
    --wrap="
        python ${SCRIPT_DIR}/merge_highPIP_eQTL_mpac.py \
            --per_chrom_dir \"${PER_CHROM_DIR}\" \
            --out           \"${OUT_FILE}\"
    ")

echo "Submitted merge job:       job ${MERGE_JOB_ID} (runs after ${ARRAY_JOB_ID})"
echo ""
echo "Monitor with: squeue -j ${ARRAY_JOB_ID},${MERGE_JOB_ID}"
echo "Final output: ${OUT_FILE}"
