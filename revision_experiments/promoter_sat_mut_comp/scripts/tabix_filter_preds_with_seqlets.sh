#!/bin/bash
# script for annotating predictions with seqlet BED files
# why are we doing this?
# so that we can calculate the correlations across models subset for seqlets in addition to just overall correlations

MPAC_PREDS='../mpac/processed_data/all.gencode.v44.canonical.protein.coding.1kb.promoters.sat.mut.updated.pos.vcf.gz'
PROMOTERAI_PREDS='../processed_data/PrimateAI_and_PromoterAI_scores.hg38.vcf.gz'
K562_PARM='/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/promoter_sat_mut_comp/processed_data/k562_parm_preds.vcf.gz'
HEPG2_PARM='/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/promoter_sat_mut_comp/processed_data/hepg2_parm_preds.vcf.gz'
### MPAC ###
# k562
echo "tabix filtering k562 mpac preds for seqlets"
tabix -R ../processed_data/bed_files/mpac_k562_merged_collapsed_seqlets4annotation.bed ${MPAC_PREDS} > ../processed_data/annotated_preds/mpac_k562_seqlet_filtered_preds.tsv

# # hepg2
echo "tabix filtering hepg2 mpac preds for seqlets"
tabix -R ../processed_data/bed_files/mpac_hepg2_merged_collapsed_seqlets4annotation.bed ${MPAC_PREDS} > ../processed_data/annotated_preds/mpac_hepg2_seqlet_filtered_preds.tsv

# # sknsh
echo "tabix filtering sknsh mpac preds for seqlets"
tabix -R ../processed_data/bed_files/mpac_sknsh_merged_collapsed_seqlets4annotation.bed ${MPAC_PREDS} > ../processed_data/annotated_preds/mpac_sknsh_seqlet_filtered_preds.tsv

# ### PromoterAI ###
echo "tabix filtering promoterAI preds for seqlets"
tabix -R ../processed_data/bed_files/promoterAI_merged_collapsed_seqlets4annotation.bed ${PROMOTERAI_PREDS} > ../processed_data/annotated_preds/promoterAI_seqlet_filtered_preds.tsv

# ### PARM ###
# # k562
echo "tabix filtering parm k562 preds for seqlets"
tabix -R ../processed_data/bed_files/PARM_K562_TFs_4_annotation.bed ${K562_PARM} > ../processed_data/annotated_preds/parm_k562_seqlet_filtered_preds.tsv

# # hepg2
echo "tabix filtering parm hepg2 preds for seqlets"
tabix -R ../processed_data/bed_files/PARM_HepG2_TFs_4_annotation.bed ${HEPG2_PARM} > ../processed_data/annotated_preds/parm_hepg2_seqlet_filtered_preds.tsv
