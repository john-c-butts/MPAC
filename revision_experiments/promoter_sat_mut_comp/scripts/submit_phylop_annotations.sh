#!/bin/bash -l
#SBATCH --job-name=annotate_vcf_w_phyloP
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition=high_mem
#SBATCH --mem=128GB
#SBATCH --time=12:00:00
#SBATCH --mail-user=john.butts@jax.org
#SBATCH --mail-type=END,FAIL

VcfAnnotateFromBigWig -in ../processed_data/mpac.250bp.promoter.preds.with.annotations.parm.act.promoterAI.parm.sat.mut4vcfAnnotateFromBigWig.vcf \
-bw ../raw_data/hg38.phyloP470way.bw \
-name phyloP -mode avg \
-out ../processed_data/mpac.250bp.promoter.preds.with.annotations.parm.act.promoterAI.parm.sat.mut4.phyloP.annotated.vcf
