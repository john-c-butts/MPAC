#!/bin/bash
#SBATCH --job-name=openTargetsExplode
#SBATCH --output=../logs/openTargetsExplode_%j.out
#SBATCH --error=../logs/openTargetsExplode_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=john.butts@jax.org
#SBATCH --mail-type=END,FAIL

python openTargets_explode.py
