# script to annotate the info column of the promoterAI predictions with gene names #

bgzip -cd ../processed_data/PrimateAI_and_PromoterAI_scores.hg38.vcf.gz | VcfAnnotateFromBed -bed ../raw_data/reformatted_bed_files/gencode.v44.protein.coding.250bp.promoters.autosomes.v2.full.ID.exon.filtered.bed -name gene_name -out ../processed_data/PrimateAI_and_PromoterAI_scores.hg38_geneAnnotated.vcf