#!/bin/bash

# script to tabix filter mpac preds for promoters from promoterAI filtering 

tabix -R ../raw_data/reformatted_bed_files/gencode.v44.protein.coding.250bp.promoters.autosomes.v2.full.ID.exon.filtered.bed ../mpac/processed_data/all.gencode.v44.canonical.protein.coding.1kb.promoters.sat.mut.updated.pos.sorted.vcf.gz > ../mpac/processed_data/all.mpac.preds.tabix.filtered.gencode.250bp.full.ID.112125.vcf