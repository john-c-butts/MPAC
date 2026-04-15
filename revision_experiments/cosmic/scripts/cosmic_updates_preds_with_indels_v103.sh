#!/bin/bash -l
#SBATCH --job-name=CSMCv103
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
#SBATCH --array=1-22
#SBATCH --output=/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/cosmic/scripts/mpac_slurm_output/slurm_%j_%a.out

# activate boda environment
conda activate boda_clone

#path where all the model artifact tgzs are
model_path="/projects/tewhey-lab/buttsj/boda_ensembl_models/chr${SLURM_ARRAY_TASK_ID}"
#get all model artifact tgzs
models=$(find $model_path -maxdepth 1 -type f -printf "%f\n" | sed "s|^|$model_path/|" | paste -sd " " -)
#get vcf
vcf_file="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/cosmic/processed_data/chrom_vcfs/v103/chr${SLURM_ARRAY_TASK_ID}_cosmic_v103_wgs_ncv.vcf"
#make a working directory for just this instance to play in...
mkdir -p "./tmp/${SLURM_ARRAY_TASK_ID}_tmp/"
cd "./tmp/${SLURM_ARRAY_TASK_ID}_tmp/"

# some comments for debugging/info
echo "[*] VCF input file ${vcf_file}"
echo "[+] Running!"

python /projects/tewhey-lab/buttsj/boda2/src/vcf_predict_indel.py \
--artifact_path ${models} \
--use_vmap True \
--vcf_file ${vcf_file} \
--fasta_file /projects/tewhey-lab/buttsj/genomes/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
--output /projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/cosmic/processed_data/mpac_preds/v103/chr${SLURM_ARRAY_TASK_ID}_cosmic_wgs_v103.vcf \
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
