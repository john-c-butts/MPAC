#!/bin/bash

# script to redo merge analysis with bedops - we can control minimum overlap here
# update 11/22/25: we landed on no minimum overlap (see archive for scripts), just trying to clean up

# store path to bed files as a variable
PATH2BEDS="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/promoter_sat_mut_comp/processed_data/bed_files"

# k562 calls
echo "merging k562 seqlet intervals with no minimum overlap"
sort-bed ${PATH2BEDS}/k562_parm_annotated_seqlets_01_112225.bed | bedmap --count --echo-map-range --echo-map-id --delim "\t" -|cut -f2- -|sort-bed --unique - > ${PATH2BEDS}/parm_k562_tangermeme_bedOps_merged_noMin_seqlets_01_112225.bed

# hepg2 calls
echo "merging hepg2 seqlet intervals with no minimum overlap"
sort-bed ${PATH2BEDS}/hepg2_parm_annotated_seqlets_01_112225.bed | bedmap --count --echo-map-range --echo-map-id --delim "\t" -|cut -f2- -|sort-bed --unique - > ${PATH2BEDS}/parm_hepg2_tangermeme_bedOps_merged_noMin_seqlets_01_112225.bed