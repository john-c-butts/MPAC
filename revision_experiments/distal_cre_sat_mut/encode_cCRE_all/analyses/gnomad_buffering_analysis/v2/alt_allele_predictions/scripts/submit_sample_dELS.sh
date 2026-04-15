#!/bin/bash
# Submit the dELS satmut sampling pipeline to SLURM

SCRIPT_DIR="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/alt_allele_predictions/scripts"
PYTHON="python3"

# ── Step 1: collect all IDs and write sampled list ────────────────────────────
echo "Submitting Step 1: collect IDs..."
JOB1=$(sbatch --parsable \
    --job-name=dELS_collect \
    --output="${SCRIPT_DIR}/logs/collect_%j.out" \
    --error="${SCRIPT_DIR}/logs/collect_%j.err" \
    --time=04:00:00 \
    --mem=16G \
    --cpus-per-task=1 \
    --wrap="${PYTHON} ${SCRIPT_DIR}/sample_dELS_satmut_preds.py --step collect")
echo "  Step 1 job ID: ${JOB1}"

# ── Step 2: filter per-chromosome (array job, 0-21 for 22 chromosomes) ────────
echo "Submitting Step 2: filter array job (depends on ${JOB1})..."
JOB2=$(sbatch --parsable \
    --job-name=dELS_filter \
    --output="${SCRIPT_DIR}/logs/filter_%A_%a.out" \
    --error="${SCRIPT_DIR}/logs/filter_%A_%a.err" \
    --time=02:00:00 \
    --mem=4G \
    --cpus-per-task=1 \
    --array=0-21 \
    --dependency=afterok:${JOB1} \
    --wrap="${PYTHON} ${SCRIPT_DIR}/sample_dELS_satmut_preds.py --step filter")
echo "  Step 2 job ID: ${JOB2}"

# ── Step 3: concatenate outputs ───────────────────────────────────────────────
echo "Submitting Step 3: concat (depends on ${JOB2})..."
JOB3=$(sbatch --parsable \
    --job-name=dELS_concat \
    --output="${SCRIPT_DIR}/logs/concat_%j.out" \
    --error="${SCRIPT_DIR}/logs/concat_%j.err" \
    --time=01:00:00 \
    --mem=4G \
    --cpus-per-task=1 \
    --dependency=afterok:${JOB2} \
    --wrap="${PYTHON} ${SCRIPT_DIR}/sample_dELS_satmut_preds.py --step concat")
echo "  Step 3 job ID: ${JOB3}"

echo ""
echo "Pipeline submitted. Monitor with: squeue -u \$USER"
echo "Final output: /projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/alt_allele_predictions/scripts/dELS_10K_satmut_preds_sample.tsv.gz"
