#!/bin/bash

# script to redo merge analysis with bedops - we can control minimum overlap here

# store path to bed files as a variable
PATH2BEDS="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/promoter_sat_mut_comp/processed_data/bed_files"

# k562 calls
echo "merging promoterAI seqlet intervals with no minimum overlap"
sort-bed ${PATH2BEDS}/promoterAI_annotated_seqlets_01.bed | bedmap --count --echo-map-range --echo-map-id --delim "\t" -|cut -f2- -|sort-bed --unique - > ${PATH2BEDS}/promoterAI_bedOps_merged_noMin_seqlets_01.bed