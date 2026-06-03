# script to calculate the seqlet coverage across all ENCODE cCREs #

#!/bin/bash

# store some variables
PATH2BEDS="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/processed_data/seqlets_calling/seqlets_by_chrom/annotated_seqlets/bed_files"
DELS_BED="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/processed_data/GRCh38-dELS-only.bed"

# k562 calls #
echo "running k562 bedtools coverage: no minimum overlap"
# get percent overlap at CRE level #
bedtools coverage -a ${DELS_BED} -b "${PATH2BEDS}/all_K562_bedOps_merged_noMin_seqlets_0.01.bed" > "${PATH2BEDS}/all_K562_01_noMin_all_dELS_cover.bed"

# hepg2 calls #
echo "running hepg2 bedtools coverage: no minimum overlap"
# get percent overlap at CRE level #
bedtools coverage -a ${DELS_BED} -b "${PATH2BEDS}/all_HepG2_bedOps_merged_noMin_seqlets_0.01.bed" > "${PATH2BEDS}/all_HepG2_01_noMin_all_dELS_cover.bed"

# sknsh calls #
echo "running sknsh bedtools coverage: no minimum overlap"
# get percent overlap at CRE level #
bedtools coverage -a ${DELS_BED} -b "${PATH2BEDS}/all_SKNSH_bedOps_merged_noMin_seqlets_0.01.bed" > "${PATH2BEDS}/all_SKNSH_01_noMin_all_dELS_cover.bed"