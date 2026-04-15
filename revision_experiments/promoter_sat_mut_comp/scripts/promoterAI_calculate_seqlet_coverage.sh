# script to calculate the seqlet coverage across all promoters #

#!/bin/bash

# store some variables
PATH2BEDS="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/promoter_sat_mut_comp/processed_data/bed_files"
PROMOTER_BED="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/promoter_sat_mut_comp/raw_data/gencode.v44.protein.coding.250bp.promoters.autosomes.v2.bed"

# k562 calls #
echo "running k562 bedtools coverage: no minimum overlap"
# get percent overlap at CRE level #
bedtools coverage -a ${PROMOTER_BED} -b "${PATH2BEDS}/promoterAI_bedOps_merged_noMin_seqlets_01.bed" > "${PATH2BEDS}/promoterAI_01_noMin_250bp_pro_cover.bed"

echo "running k562 bedtools coverage: 2bp overlap"
# get percent overlap at CRE level #
bedtools coverage -a ${PROMOTER_BED} -b "${PATH2BEDS}/promoterAI_bedOps_merged_Min2ovr_seqlets_01.bed" > "${PATH2BEDS}/promoterAI_01_Min2ovr_250bp_pro_cover.bed"

echo "running k562 bedtools coverage: 4bp overlap"
# get percent overlap at CRE level #
bedtools coverage -a ${PROMOTER_BED} -b "${PATH2BEDS}/promoterAI_bedOps_merged_Min4ovr_seqlets_01.bed" > "${PATH2BEDS}/promoterAI_01_Min4ovr_250bp_pro_cover.bed"