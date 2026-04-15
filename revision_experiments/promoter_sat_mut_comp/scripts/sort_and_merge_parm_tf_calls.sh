# script to sort and merge bed file of significant tf calls from PARM

# store path to bed files as a variable
PATH2BEDS="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/promoter_sat_mut_comp/processed_data/bed_files"

# k562 calls
echo "merging significant TF calls from PARM for K562 with no minimum overlap"
sort-bed ${PATH2BEDS}/parm_k562_sig_tf_hits_all.bed| bedmap --count --echo-map-range --echo-map-id --delim "\t" -|cut -f2- -|sort-bed --unique - > ${PATH2BEDS}/PARM_K562_sig_TFs_bedOps_merged_noMin_112225.bed

# hepg2 calls
echo "merging significant TF calls from PARM for HepG2 with no minimum overlap"
sort-bed ${PATH2BEDS}/parm_hepg2_sig_tf_hits_all.bed| bedmap --count --echo-map-range --echo-map-id --delim "\t" -|cut -f2- -|sort-bed --unique - > ${PATH2BEDS}/PARM_HepG2_sig_TFs_bedOps_merged_noMin_112225.bed