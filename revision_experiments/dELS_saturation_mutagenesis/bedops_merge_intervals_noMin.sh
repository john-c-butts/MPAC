#!/bin/bash

# script to redo merge analysis with bedops - we can control minimum overlap here

# store path to bed files as a variable
PATH2BEDS="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/processed_data/seqlets_calling/seqlets_by_chrom/annotated_seqlets/bed_files"

# k562 calls
echo "merging k562 seqlet intervals with no minimum overlap"
bedmap --count --echo-map-range --echo-map-id --delim "\t" ${PATH2BEDS}/all_K562_seqlets_0.01.bed |cut -f2- -|sort-bed --unique - > ${PATH2BEDS}/all_K562_bedOps_merged_noMin_seqlets_0.01.bed

# hepg2 calls
echo "merging hepg2 seqlet intervals with no minimum overlap"
bedmap --count --echo-map-range --echo-map-id --delim "\t" ${PATH2BEDS}/all_HepG2_seqlets_0.01.bed |cut -f2- -|sort-bed --unique - > ${PATH2BEDS}/all_HepG2_bedOps_merged_noMin_seqlets_0.01.bed

# sknsh calls
echo "merging sknsh seqlet intervals with no minimum overlap"
bedmap --count --echo-map-range --echo-map-id --delim "\t" ${PATH2BEDS}/all_SKNSH_seqlets_0.01.bed |cut -f2- -|sort-bed --unique - > ${PATH2BEDS}/all_SKNSH_bedOps_merged_noMin_seqlets_0.01.bed