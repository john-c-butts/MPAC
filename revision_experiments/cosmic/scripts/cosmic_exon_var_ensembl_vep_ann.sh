#!/bin/bash
#SBATCH --job-name=csmcXnVEPall
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --time=72:00:00
#SBATCH --mem=32G
#SBATCH --mail-user=john.butts@jax.org
#SBATCH --mail-type=END,FAIL

module load singularity

singularity exec /projects/tewhey-lab/buttsj/ensembl_vep/vep.sif \
vep --dir /projects/tewhey-lab/buttsj/vep_data/ \
--cache --offline \
--format vcf \
--tab \
--force_overwrite \
--regulatory \
--nearest transcript \
--distance 1000 \
--fork 8 \
--input_file ../processed_data/cosmic.v98.exon.intersecting.vars4ensembl.vep.vcf \
--output_file ../processed_data/cosmic.v98.exon.intersecting.vars4ensembl.vep.annotated.tsv \
