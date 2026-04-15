#!/bin/bash -l
#SBATCH --job-name=phyloP_241
#SBATCH --partition=compute
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --mail-user=john.butts@jax.org
#SBATCH --mail-type=END,FAIL

VcfAnnotateFromBigWig -in ../processed_data/mpac.250bp.promoter.preds.with.annotations.parm.act.promoterAI.parm.sat.mut4vcfAnnotateFromBigWig.vcf \
-bw ../raw_data/241-mammalian-2020v2.bigWig \
-name phyloP -mode avg \
-out ../processed_data/mpac.250bp.promoter.preds.with.annotations.parm.act.promoterAI.parm.sat.mut4.phyloP.annotated.241.mamm.vcf
