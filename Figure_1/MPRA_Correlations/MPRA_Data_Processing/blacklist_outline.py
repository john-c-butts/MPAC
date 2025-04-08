#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 10:02:49 2023

@author: buttsj
"""

import pandas as pd
import os
from collections import Counter
import matplotlib.pyplot as mpl
from scipy import stats

# open validation MPRA data from Sager
mpra_val = pd.read_csv('/Users/buttsj/Dropbox (JAX)/Variant_Effects/ukbb_gtex_mpra/empirical_ukbb_gtex/MPRA_ALL_v3.txt',
                       sep = ' ',
                       low_memory=False)
# filter validation for only REF alleles
ref_ids = [i for i in mpra_val['IDs'] if 'wC' in i and i.split(':')[-2] == 'R']
mpra_val_ref = mpra_val[mpra_val['IDs'].isin(ref_ids)]
# open my gtex mpra data
gtex_mpra = pd.read_csv('/Users/buttsj/Dropbox (JAX)/Variant_Effects/ukbb_gtex_mpra/gtex.all.chroms.mpra.data.txt',
                        sep = '\t',
                        low_memory=False)
gtex_mpra = gtex_mpra.drop_duplicates(subset='A_log2FC')
gtex_mpra_k562 = gtex_mpra[gtex_mpra['cell_type'] == 'K562']
# correlate Sager's data and my data
# combine validation set and gtex_mpra_k562 data
val_gtex_k562_merge = mpra_val_ref.merge(gtex_mpra_k562, how='inner', left_on='IDs', right_on='ID')
# plot correlation of matched IDs between datasets
mpl.figure(dpi=300)
mpl.scatter(x=val_gtex_k562_merge['K562_mean'],
            y=val_gtex_k562_merge['A_log2FC'],
            s=2,
            marker='.',
            alpha=.4)
mpl.title('MPRA_ALL_v3 vs. GTEX MPRA K562 REF')
mpl.xlabel('MPRA_ALL_v3 MPRA')
mpl.ylabel('GTEX MPRA K562')
# look at highly discordant points between gtex mpra and mpra all
# get abs(delta) between mpra_all and gtex mpra for sorting DF
val_gtex_k562_merge.loc[:,]['delta_act'] = [abs(i - j) for i,j in zip(val_gtex_k562_merge['K562_mean'],
                                                                      val_gtex_k562_merge['A_log2FC'])]
val_gtex_k562_merge_discord = val_gtex_k562_merge.sort_values(by='delta_act', ascending=False).head(10)

# get pearson between both datasets
stats.pearsonr(x=val_gtex_k562_merge['K562_mean'],
               y=val_gtex_k562_merge['A_log2FC'])
# pearsonR .9965
# get number of IDs where mpra values are perfectly matched
[1 if i == j else 0 for i,j in zip(val_gtex_k562_merge['K562_mean'],
                                   val_gtex_k562_merge['A_log2FC'])].count(1)
# get number of positions where both ref and alt in gtex k562 have log2FC SE < 1
len(gtex_mpra_k562[(gtex_mpra_k562['A_log2FC_SE'] < 1) &
                   (gtex_mpra_k562['B_log2FC_SE'] < 1)])
# 178,801
# get IDs of those variants that are being dropped for log2FC SE
k562_log2fc_se_fail_ids = gtex_mpra_k562[(gtex_mpra_k562['A_log2FC_SE'] >= 1) |
                                         (gtex_mpra_k562['B_log2FC_SE'] >= 1)]['ID'].tolist()
# no values are identical between datasets
# open predictions for all of gtex
gtex_preds = pd.read_csv('/Users/buttsj/Dropbox (JAX)/Variant_Effects/ukbb_gtex_mpra/gtex.all.chroms.center.only.vcf',
                         sep = '\t')
# merge predictions and gtex mpra
# match gtex mpra hg38 ID format to predictions chr:pos:ref:alt
gtex_mpra_k562.loc[:,]['hg38_id'] = [(':').join([i.split('_')[0], 
                                                 i.split('_')[1], 
                                                 i.split('_')[2], 
                                                 i.split('_')[3]]) for i in gtex_mpra_k562['variant_hg38']]

gtex_mpra_preds_merge = gtex_preds.merge(gtex_mpra_k562, how='inner', left_on='id', right_on='hg38_id')
# match predictions and all seq mpra

all_preds_merge = gtex_preds.merge()
# missing ~15K variants, all of which are non-snp variants
mpl.figure(dpi=300)
mpl.scatter(y=[float(i) for i in gtex_mpra_preds_merge['k562_ref_pred']], 
            x=gtex_mpra_preds_merge['A_log2FC'],
            s=2,
            marker='.',
            alpha=.4)
mpl.xlim(-5, 17.5)
mpl.ylim(-5, 17.5)
stats.pearsonr(x=[float(i) for i in gtex_mpra_preds_merge['k562_ref_pred']],
               y=gtex_mpra_preds_merge['A_log2FC'])

# calculate correlations by chromosome
# add chrom column for filtering
gtex_mpra_preds_merge.loc[:,]['chr'] = [i.split(':')[0] for i in gtex_mpra_preds_merge['hg38_id']]
pearsonr_chrom_dict = {}
for i in gtex_mpra_preds_merge['chr'].unique():
    chrom_df = gtex_mpra_preds_merge[gtex_mpra_preds_merge['chr'] == i]
    pearsonr_chrom_dict.update({i : stats.pearsonr(x=[float(i) for i in chrom_df['k562_ref_pred']],
                                                   y=chrom_df['A_log2FC'])[0]})
    mpl.figure(dpi=300)
    mpl.scatter(x=[float(i) for i in chrom_df['k562_ref_pred']],
                y=chrom_df['A_log2FC'],
                s=2,
                marker='.',
                alpha=.4)
    mpl.title(i)
# 




# open list of IDs that switch chromosomes between hg19 and hg38
chrom_move_ids = pd.read_csv('/Users/buttsj/Dropbox (JAX)/Variant_Effects/ukbb_gtex_mpra/empirical_ukbb_gtex/chrom_move_ids.txt',
                             sep = '\t')
chrom_move_mpra_preds = gtex_mpra_preds_merge[gtex_mpra_preds_merge['ID'].isin(chrom_move_ids['ID'])]
