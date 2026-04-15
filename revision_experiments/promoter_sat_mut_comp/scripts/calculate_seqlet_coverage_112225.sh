# script to calculate the seqlet coverage across all promoters #

#!/bin/bash

# store some variables
PATH2BEDS="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/promoter_sat_mut_comp/processed_data/bed_files"
PROMOTER_BED="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/promoter_sat_mut_comp/raw_data/reformatted_bed_files/gencode.v44.protein.coding.250bp.promoters.autosomes.v2.full.ID.exon.filtered.bed"

# k562 calls #
echo "running k562 bedtools coverage: no minimum overlap"
# get percent overlap at CRE level #
bedtools coverage -a ${PROMOTER_BED} -b "${PATH2BEDS}/K562_bedOps_merged_noMin_seqlets_01_112225.bed" > "${PATH2BEDS}/K562_01_noMin_250bp_pro_cover_112225.bed"

# hepg2 calls #
echo "running hepg2 bedtools coverage: no minimum overlap"
# get percent overlap at CRE level #
bedtools coverage -a ${PROMOTER_BED} -b "${PATH2BEDS}/HepG2_bedOps_merged_noMin_seqlets_01_112225.bed" > "${PATH2BEDS}/HepG2_01_noMin_250bp_pro_cover_112225.bed"

# sknsh calls #
echo "running sknsh bedtools coverage: no minimum overlap"
# get percent overlap at CRE level #
bedtools coverage -a ${PROMOTER_BED} -b "${PATH2BEDS}/SKNSH_bedOps_merged_noMin_seqlets_01_112225.bed" > "${PATH2BEDS}/SKNSH_01_noMin_250bp_pro_cover_112225.bed"