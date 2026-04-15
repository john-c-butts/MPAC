#!/bin/bash
#SBATCH --job-name=mpac_array_test
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition=gpu_v100
#SBATCH --qos=gpu_training
#SBATCH --gres=gpu:1
#SBATCH --mem=128GB
#SBATCH --cpus-per-task=12
#SBATCH --time=96:00:00
#SBATCH --mail-user=john.butts@jax.org
#SBATCH --mail-type=END,FAIL

mod_path="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/bed_files/mpac_array_test/chr22_mods"

python /projects/tewhey-lab/buttsj/boda2/src/vcf_predict.py \
--artifact_path $mod_path/model_artifacts__20230111_190417__993143.tar.gz $mod_path/model_artifacts__20230111_192235__715916.tar.gz $mod_path/model_artifacts__20230111_232551__863644.tar.gz $mod_path/model_artifacts__20230112_000934__941678.tar.gz $mod_path/model_artifacts__20230112_003327__287605.tar.gz $mod_path/model_artifacts__20230112_003539__758935.tar.gz $mod_path/model_artifacts__20230112_014853__150994.tar.gz $mod_path/model_artifacts__20230112_020038__431749.tar.gz $mod_path/model_artifacts__20230112_024037__325298.tar.gz $mod_path/model_artifacts__20230112_182326__436818.tar.gz \
--use_vmap TRUE \
--vcf_file /projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/bed_files/vcf_like/chunked_vcfs/k562_chr22_cre_chunk_8.vcf \
--fasta /projects/tewhey-lab/buttsj/genomes/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
--output /projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/bed_files/mpac_array_test/k562_chr22_test/k562_chr22_cre_chunk_8_mpac_017.vcf \
--relative_start 9 \
--relative_end 181 \
--step_size 10 \
--raw_predictions FALSE \
--strand_reduction mean \
--window_reduction mean \
--feature_ids K562 HepG2 SKNSH

