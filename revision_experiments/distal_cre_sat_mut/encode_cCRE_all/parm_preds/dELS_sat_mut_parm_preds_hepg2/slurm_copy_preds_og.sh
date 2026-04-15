#!/bin/bash

#SBATCH --time=30:00
#SBATCH -n 1
#SBATCH -p compute
#SBATCH -q batch
#SBATCH --mem 4G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=john.butts@jax.org
#SBATCH --job-name=cp_k_parm_all
#SBATCH --array=0-128,130-196,198-338,340-344,346-580,582-774,776-832,834-1395,1397-1469
#SBATCH --output=slurm_out/slurm_%j_%a.out

cp -r chunk${SLURM_ARRAY_TASK_ID} ../final_sat_mut_preds/k562