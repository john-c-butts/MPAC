#!/usr/bin/env bash
#SBATCH --job-name=build_605_altBg_050726
#SBATCH --partition=compute
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/gtex_phenocopy_with_alt_background/logs/build_605_altBg_050726_%j.out
#SBATCH --error=/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/gtex_phenocopy_with_alt_background/logs/build_605_altBg_050726_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=john.butts@jax.org

mkdir -p /projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/gtex_phenocopy_with_alt_background/logs

source /projects/tewhey-lab/buttsj/miniforge3/etc/profile.d/conda.sh
conda activate boda_clone

cd /projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/gtex_phenocopy_with_alt_background

python build_final_dataset_605_altBg_050726.py
