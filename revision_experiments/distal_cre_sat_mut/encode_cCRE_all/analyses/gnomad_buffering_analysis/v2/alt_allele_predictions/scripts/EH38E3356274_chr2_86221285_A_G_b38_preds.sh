#!/bin/bash -l
#SBATCH --job-name=EH38E3356274
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition=gpu_a100
#SBATCH --qos=gpu_inference
#SBATCH --gres=gpu:1
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=12
#SBATCH --time=6:00:00
#SBATCH --mail-user=john.butts@jax.org
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm_out/slurm_%j_%a.out

# activate boda environment
conda activate boda_clone

#path where all the model artifact tgzs are
model_path="/projects/tewhey-lab/buttsj/boda_ensembl_models/chr2"
#get all model artifact tgzs
models=$(find $model_path -maxdepth 1 -type f -printf "%f\n" | sed "s|^|$model_path/|" | paste -sd " " -)
#get vcf
vcf_file="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/alt_allele_predictions/processed_data/vcfs4satmut/EH38E3356274_chr2_86221285_A_G_b38.vcf"
#make a working directory for just this instance to play in...
mkdir -p "./tmp/tmp/"
cd "./tmp/tmp/"

# some comments for debugging/info
echo "[*] VCF input file ${vcf_file}"
echo "[+] Running!"

python /projects/tewhey-lab/buttsj/boda2/src/vcf_predict_indel.py \
--artifact_path ${models} \
--use_vmap True \
--vcf_file ${vcf_file} \
--fasta_file /projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/alt_allele_predictions/processed_data/alt_genomes/chr2_86221285_A_G_b38.fa \
--output /projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/alt_allele_predictions/mpac_preds/EH38E3356274_chr2_86221285_A_G_b38_preds.vcf \
--window_size 200 \
--relative_start 9 \
--relative_end 180 \
--step_size 10 \
--raw_predictions False \
--strand_reduction mean \
--window_reduction mean \
--use_simple_padding True \
--extract_full_sequences False \
--feature_ids K562 HepG2 SKNSH

echo "[+] Finished!"
