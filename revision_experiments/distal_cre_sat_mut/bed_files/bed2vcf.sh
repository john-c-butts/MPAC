#!/bin/bash

#SBATCH --job-name=bed2vcf
#SBATCH -N 1
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=64G
#SBATCH --mail-user=john.butts@jax.org
#SBATCH --mail-type=END,FAIL

python bed2vcf.py
