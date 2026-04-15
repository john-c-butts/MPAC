#!/bin/bash
#SBATCH --job-name=promoterAI_seqlets_v2
#SBATCH --partition=high_mem
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --output=/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/promoter_sat_mut_comp/logs/promoterAI_seqlets_v2_%j.out
#SBATCH --error=/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/promoter_sat_mut_comp/logs/promoterAI_seqlets_v2_%j.err
#SBATCH --mail-user=john.butts@jax.org
#SBATCH --mail-type=END,FAIL

V2="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/promoter_sat_mut_comp/v2"

mkdir -p "$(dirname "${V2}/../logs/x")"

cd "${V2}"
python3 v2_promoterAI_seqlet_calling.py
