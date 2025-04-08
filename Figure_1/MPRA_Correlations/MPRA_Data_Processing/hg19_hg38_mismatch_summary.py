#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 13 09:50:23 2023

@author: buttsj
"""

import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as mpl
import operator
from collections import Counter
### Script to try and summarize chromosome/sequence discordance between hg19/hg38 
### and generate modified VCFs for checking with checkVCF.py

# open VCFs passed to Sager and count mismatches
# change to directory of VCFs w/ hg19 IDs
os.chdir('/Users/buttsj/Dropbox (JAX)/Variant_Effects/ukbb_gtex_mpra/ukbb_gtex_mpra_chrom_vcfs/gtex/gtex_hg19_ID_vcfs/')
# open each VCF and store as a dictionary
vcf_dict = {}
for i in [i for i in os.listdir() if i.endswith('.vcf')]:
    vcf_dict.update({i.split('_')[0] : pd.read_csv(i, sep = '\t')})
# iterate through dictionary of VCFs and count chromosome jumps by chromosome
def mismatches_per_chrom (chrom_df_dict):
    mismatch_dict = {}
    for i in chrom_df_dict.keys():
        df = chrom_df_dict.get(i)
        # drop 'chr' from hg38 chromosomes and convert to list
        hg38_chr = [i.split('chr')[-1] for i in df['#chrom']]
        # get hg19 chromosomes
        hg19_chr = [i.split(':')[0] for i in df['id']]
        # compare hg38 and hg19 chromosome numbers
        n_mismatches = [1 if i == j else 0 for i,j in zip(hg38_chr,
                                                          hg19_chr)].count(0)
        # update dictionary with number of mismatches per chromosome
        mismatch_dict.update({i : n_mismatches})
    return mismatch_dict
mismatch_by_chrom = mismatches_per_chrom(vcf_dict)
mismatch_by_chrom = dict(sorted(mismatch_by_chrom.items(), key=operator.itemgetter(1),reverse=True))
# make a DF of mismatches by chromosome for plotting
mismatch_df = pd.DataFrame({'chromosome' : mismatch_by_chrom.keys(),
                            'n chrom changes' : mismatch_by_chrom.values()})

# plot n chromosome changes
mpl.figure(dpi=300)
sns.barplot(data=mismatch_df, 
            x='chromosome',
            y='n chrom changes',
            palette = 'Paired')
mpl.xticks(rotation=90)
mpl.title('N Times Variant Changes Chromosome by Chromosome')
# get the IDs of variants that change chromosome
def get_mismatch_ids (chrom_df_dict):
    mismatch_all = []
    for i in chrom_df_dict.keys():
        df = chrom_df_dict.get(i)
        # drop 'chr' from hg38 chromosomes and convert to list
        hg38_chr = [i.split('chr')[-1] for i in df['#chrom']]
        # get hg19 chromosomes
        hg19_chr = [i.split(':')[0] for i in df['id']]
        # compare hg38 and hg19 chromosome numbers
        mismatch_ids = [k for i,j,k in zip(hg38_chr,
                                           hg19_chr,
                                           df['id'].tolist()) if i != j]
        mismatch_all = mismatch_all + mismatch_ids
    return mismatch_all
chrom_move_ids = get_mismatch_ids(vcf_dict)
# edit chromosome VCFs for checking with checkVCF.py
# changing coordinates to hg19, header names to all caps and add 'format' column
def vcf4check (vcf_dict):
    df_list = []
    for i in vcf_dict.keys():
        # reset columns
        #df = vcf_dict.get(i)
        df_rename = vcf_dict.get(i).rename(columns={'#chrom' : '#CHROM',
                                                    'pos' : 'POS',
                                                    'id' : 'ID',
                                                    'ref' : 'REF',
                                                    'alt' : 'ALT',
                                                    'qual' : 'QUAL',
                                                    'filter' : 'FILTER',
                                                    'info' : 'INFO'})
        # add 'format' column to DF
        df_rename.loc[:,]['FORMAT'] = ['.' for i in range(len(df_rename))]
        # replace CHROM with hg19 values
        df_rename['#CHROM'] = ['chr' + str(i.split(':')[0]) for i in df_rename['ID']]
        # replace POS with hg19 values
        df_rename['POS'] = [i.split(':')[1] for i in df_rename['ID']]
        df_list.append(df_rename)
    hg19_vcf = pd.concat(df_list)
    return hg19_vcf      
hg19_all_chroms = vcf4check(vcf_dict)        
# save hg19 vcf to file for checking with checkVCF.py
hg19_all_chroms.to_csv('/Users/buttsj/Dropbox (JAX)/checkVCF_py/gtex_hg194checkVCF.vcf',
                       sep = '\t',
                       index=False)
# open checkVCF.py output for non-matching REFs
check_vcf_ref = pd.read_csv('/Users/buttsj/Dropbox (JAX)/checkVCF_py/gtex_hg19_check.check.ref',
                            sep = '\t',
                            header=None)
# make counter object of IDs to drop those with multiple entries in chr
id_counter_full = Counter([i for i in hg19_all_chroms['ID']])
id_counter_chr_pos = Counter((':').join([i.split(':')[0], i.split(':')[1]]) for i in hg19_all_chroms['ID'])
# drop duplicate values from chr_pos Counter object
unique_hg19_chr_pos = [i for i in id_counter_chr_pos.keys() if id_counter_chr_pos.get(i) == 1]
# add chr_pos id to hg19_vcf for filtering
hg19_all_chroms.loc[:,]['chr_pos'] = [(':').join([i.split(':')[0],
                                                  i.split(':')[1]]) for i in hg19_all_chroms['ID']]
# filter for only those positions with unique chr_pos IDs
hg19_unique_chr_pos = hg19_all_chroms[hg19_all_chroms['chr_pos'].isin(unique_hg19_chr_pos)]
# parse the unmatched REFs from checKVCF output
mismatch_chr_pos = [(':').join([i.split(':')[0],
                                i.split(':')[1]]) for i in check_vcf_ref[1]]
# mismatch VCF
hg19_mismatch_ref = hg19_unique_chr_pos[hg19_unique_chr_pos['chr_pos'].isin(mismatch_chr_pos)]
# all mismatches appear to have swapped strands? make VCF of only those variants from hg19 ID
# concatenate all vcfs together
all_chroms_vcf = pd.concat(vcf_dict.values())
hg19_mismatch_vcf4check = pd.DataFrame({'#CHROM' : ['chr' + str(i.split(':')[0]) for i in all_chroms_vcf['id']],
                                        'POS' : [i.split(':')[1] for i in all_chroms_vcf['id']],
                                        'ID' : [i for i in all_chroms_vcf['id']],
                                        'REF' : [i.split(':')[2] for i in all_chroms_vcf['id']],
                                        'ALT' : [i.split(':')[3] for i in all_chroms_vcf['id']],
                                        'QUAL' : ['.' for i in range(len(all_chroms_vcf))],
                                        'FILTER' : ['.' for i in range(len(all_chroms_vcf))],
                                        'INFO' : [i for i in all_chroms_vcf['info']],
                                        'FORMAT' : ['.' for i in range(len(all_chroms_vcf))]})
# save hg19 vcf for vcfCheck
hg19_mismatch_vcf4check.to_csv('/Users/buttsj/Dropbox (JAX)/checkVCF_py/hg19_ref_alt_all_chroms_dot.vcf',
                               sep = '\t',
                               index=False)
# no mismatches with corrected vcf! figure out changes of ref/alt between hg19 and hg38
# write list of hg 19 IDs that have changed ref/alt between hg19 and hg38 to file
hg19_hg38_flips = pd.DataFrame({'hg19 id' : [i for i in hg19_mismatch_ref['ID']]})
hg19_hg38_flips.to_csv('/Users/buttsj/Dropbox (JAX)/Variant_Effects/ukbb_gtex_mpra/ukbb_gtex_mpra_chrom_vcfs/gtex/hg19_hg38_ref_alt_flip_ids.txt',
                       sep = '\t',
                       index=False)



