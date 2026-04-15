#!/bin/bash

#SBATCH --job-name=pull_fasta_cres
#SBATCH -N 1
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=64G
#SBATCH --mail-user=john.butts@jax.org
#SBATCH --mail-type=END,FAIL

bedtools getfasta -tab -fi /projects/tewhey-lab/buttsj/genomes/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -bed GRCh38-dELS_only.bed > fasta_seqs/GRCh38-dELS_only.fasta
