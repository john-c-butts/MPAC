#!/bin/bash

#SBATCH --job-name=pull_fasta_cres
#SBATCH -N 1
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=64G
#SBATCH --mail-user=john.butts@jax.org
#SBATCH --mail-type=END,FAIL

bedtools getfasta -tab -fi /projects/tewhey-lab/buttsj/genomes/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -bed processed_data/K562.dELS.only.bed > fasta_seqs/K562.dELS.fasta
bedtools getfasta -tab -fi /projects/tewhey-lab/buttsj/genomes/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -bed processed_data/HepG2.dELS.only.bed > fasta_seqs/HepG2.enhancers.fasta
bedtools getfasta -tab -fi /projects/tewhey-lab/buttsj/genomes/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -bed processed_data/SK-N-SH.dELS.only.bed > fasta_seqs/SK-N-SH.enhancers.fasta
