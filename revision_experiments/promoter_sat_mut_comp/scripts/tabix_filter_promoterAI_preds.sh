#!/bin/bash

# script to tabix filter promoterAI preds for 250bp promoters

tabix -R ../raw_data/reformatted_bed_files/gencode.v44.protein.coding.250bp.promoters.autosomes.v2.full.ID.exon.filtered.bed ../processed_data/PrimateAI_and_PromoterAI_scores.hg38.vcf.gz > ../processed_data/PrimateAI_and_PromoterAI_scores.hg38.250bp.filtered.vcf 