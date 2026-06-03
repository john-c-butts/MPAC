# script to concatenate all chromosome divided seqlets into single bed file for calculating coverage #

# import packages
import os

# change to working directory #
os.chdir('/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/processed_data/seqlets_calling/seqlets_by_chrom/annotated_seqlets/bed_files')

# get list of K562 seqlets to cat #
k562_seqs2cat = [i for i in os.listdir() if '_K562_0.01.bed' in i]

# get list of HepG2 seqlets to cat #
hepg2_seqs2cat = [i for i in os.listdir() if '_HepG2_0.01.bed' in i]

# get list of SKNSH seqlets to cat #
sknsh_seqs2cat = [i for i in os.listdir() if '_SKNSH_0.01.bed' in i]

# check that all lengths are 22 (ie all autosomes should be present)
if len(k562_seqs2cat) == 22 and len(hepg2_seqs2cat) == 22 and len(sknsh_seqs2cat) == 22:
    print('all autosomes present')
else:
    print('missing autosome, check bed files')

# build cat and sort commands
k562_cat_cmd = f'cat {' '.join(k562_seqs2cat)} | bedtools sort > all_K562_seqlets_0.01.bed'
hepg2_cat_cmd = f'cat {' '.join(hepg2_seqs2cat)} | bedtools sort > all_HepG2_seqlets_0.01.bed'
sknsh_cat_cmd = f'cat {' '.join(sknsh_seqs2cat)} | bedtools sort > all_SKNSH_seqlets_0.01.bed'

# cat files
print('concatenating K562 seqlet calls')
os.system(k562_cat_cmd)
print('\n')
print('concatenating HepG2 seqlet calls')
os.system(hepg2_cat_cmd)
print('\n')
print('concatenating SKNSH seqlet calls')
os.system(sknsh_cat_cmd)