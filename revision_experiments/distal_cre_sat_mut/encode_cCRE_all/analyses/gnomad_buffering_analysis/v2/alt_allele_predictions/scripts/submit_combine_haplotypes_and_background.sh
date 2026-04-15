#!/bin/bash
#SBATCH --job-name=combine_haplos
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=john.butts@jax.org
#SBATCH --mail-type=END,FAIL

SCRIPT_DIR="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/alt_allele_predictions/scripts"

python ${SCRIPT_DIR}/combine_haplotypes_and_background.py
