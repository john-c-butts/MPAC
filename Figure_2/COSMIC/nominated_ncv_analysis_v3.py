#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 11:22:18 2023

@author: buttsj
"""
#%%
# import packages
import pandas as pd
import os
import matplotlib.pyplot as mpl
import seaborn as sns
from tqdm import tqdm
import numpy as np
pd.set_option('display.max_columns', None)
#%%
# open all predictions for analysis
fuxman_ncv_preds = pd.concat([pd.read_csv(f'/Users/buttsj/Dropbox (JAX)/Variant_Effects/fuxman_bass_ets_mpra/predictions/full_preds/{i}', sep = '\t') for i in os.listdir('/Users/buttsj/Dropbox (JAX)/Variant_Effects/fuxman_bass_ets_mpra/predictions/full_preds/')])
#%%
fuxman_mpra_data = pd.read_excel('41467_2023_36535_MOESM5_ESM.xlsx', sheet_name='MPRAs')
# define function for calling emVars
def call_tfa_bt_emvars (mpra,
                        cell_line):
    # add min fdr for ref/alt activity
    # BH
    mpra['max_bh_log_FDR'] = [max([ref_bf, alt_bf]) for ref_bf, alt_bf in zip(mpra['A.logPadj_BH'],
                                                                              mpra['B.logPadj_BH'])]
    # BF
    mpra['max_bf_log_FDR'] = [max([ref_bf, alt_bf]) for ref_bf, alt_bf in zip(mpra['A.logPadj_BF'],
                                                                              mpra['B.logPadj_BF'])]
    # take higher of two
    mpra['max_fdr'] = [max(bh, bf) for bh, bf in zip(mpra['max_bh_log_FDR'],
                                                     mpra['max_bf_log_FDR'])]
    tfa_vars = mpra[mpra['NCV class'].isin(['TFA-BT NCV'])]
    # filter for cell type
    cell_type_mpra = tfa_vars[tfa_vars['Cell line'] == cell_line]
    # subset for the TFA-BT-NCV and TFA-BT -  unobseved variants
    active_cell = cell_type_mpra[cell_type_mpra['max_fdr'] > -(np.log10(.01))]
    # call emVars
    emvars = active_cell[active_cell['Skew.logFDR'] > -(np.log10(.05))]
    print(f'There are {len(active_cell)} active elements in {cell_line}')
    print(f'There are {len(emvars)} emVars in {cell_line}')
    print(f'The proportion of TFA-BT NCV emVars in {cell_line} is {round(len(emvars) / len(active_cell),2)}')
    return emvars
#%%
# call Jurkat emVars
jurkat_emvars = call_tfa_bt_emvars(fuxman_mpra_data, 'Jurkat')
# call HT-29 emVars
ht_29_emvars = call_tfa_bt_emvars(fuxman_mpra_data, 'HT-29')
# call SK-MEL-28 emVars
sk_mel_28_emvars = call_tfa_bt_emvars(fuxman_mpra_data, 'SK-MEL-28')
# print(len(pd.concat([sk_mel_28_emvars, jurkat_emvars, ht_29_emvars])['NCV (chr:position:ref:alt)'].unique()))
# above line returns 765 emVars consistent with paper
#%%
# generate a list of non-X chromosome emVars to calculate predicted emVar rates
emvars2verify = [i for i in pd.concat([sk_mel_28_emvars,
                                       ht_29_emvars,
                                       jurkat_emvars])['NCV (chr:position:ref:alt)'].unique() if i.split(':')[0] != 'X']
# total of 761 emVars
#%% 
# open bed file with hg19/hg38 coordinates
hg19_38_bed = pd.read_csv('full_fuxman_hg38.bed',
                          sep = '\t',
                          header = None)
# add a column with 38
hg19_38_bed['hg38_id'] = [(':').join([chrom.split('chr')[-1],
                                     str(pos),
                                     var.split(':')[-2],
                                     var.split(':')[-1]]) for chrom, pos, var in zip(hg19_38_bed[0],
                                                                                     hg19_38_bed[1],
                                                                                     hg19_38_bed[3])]
#%%
# make into a dictionary for easily matching IDs
hg38_hg19_dict = dict(zip(hg19_38_bed['hg38_id'],
                          hg19_38_bed[3]))
#%%
# reformat prediction VCFs
# define function to convert predictions to final df
def vcf2df (pred_df):
    # make lists of preds for new column
    # k562
    k_ref = []
    k_alt = []
    k_skew = []
    # hepg2
    h_ref = []
    h_alt = []
    h_skew = []
    # sknsh
    s_ref = []
    s_alt = []
    s_skew = []
    # parse predictions in 'INFO' column
    for i in tqdm(pred_df['INFO']):
        all_preds = i.split(';')
        # k562
        k_ref.append(float(all_preds[0].split('=')[-1]))
        k_alt.append(float(all_preds[3].split('=')[-1]))
        k_skew.append(float(all_preds[6].split('=')[-1]))
        # hepg2
        h_ref.append(float(all_preds[1].split('=')[-1]))
        h_alt.append(float(all_preds[4].split('=')[-1]))
        h_skew.append(float(all_preds[7].split('=')[-1]))
        # sknsh
        s_ref.append(float(all_preds[2].split('=')[-1]))
        s_alt.append(float(all_preds[5].split('=')[-1]))
        s_skew.append(float(all_preds[8].split('=')[-1]))
    
    df = pd.DataFrame({'chrom' : pred_df['chrom'],
                       'pos' : pred_df['pos'],
                       'id' : pred_df['id'],
                       'ref' : pred_df['ref'],
                       'alt' : pred_df['alt'],
                       'k562_ref_pred' : k_ref,
                       'k562_alt_pred' : k_alt,
                       'k562_skew_pred' : k_skew,
                       'hepg2_ref_pred' : h_ref,
                       'hepg2_alt_pred' : h_alt,
                       'hepg2_skew_pred' : h_skew,
                       'sknsh_ref_pred' : s_ref,
                       'sknsh_alt_pred' : s_alt,
                       'sknsh_skew_pred' : s_skew})
    return df
#%%
fuxman_preds = vcf2df(fuxman_ncv_preds)
# add column with ID of variant
fuxman_preds['hg38_id'] = [(':').join([chrom.split('chr')[-1],
                                       str(pos),
                                       ref,
                                       alt]) for chrom, pos, ref, alt in zip(fuxman_preds['chrom'],
                                                                             fuxman_preds['pos'],
                                                                             fuxman_preds['ref'],
                                                                             fuxman_preds['alt'])]
# add hg19 id
fuxman_preds['hg19_id'] = [hg38_hg19_dict.get(i) for i in fuxman_preds['hg38_id']]
#%%
emvar_preds = fuxman_preds[fuxman_preds['hg19_id'].isin(emvars2verify)]
# %%
# calculate the proportion of emVars called by MPAC
emvar_preds['emvar'] = [1 if abs(max([k_skew,
                                 h_skew,
                                 s_skew], key = abs)) > .5 else 0 for k_skew, h_skew, s_skew in zip(emvar_preds['k562_skew_pred'],
                                                                                                   emvar_preds['hepg2_skew_pred'],
                                                                                                   emvar_preds['sknsh_skew_pred'])]
n_emvars = emvar_preds['emvar'].tolist().count(1)
# %%
print(f'The proportion of emVars called by TFA-BT called by MPAC is {round(n_emvars / len(emvar_preds),2)}')
# %%
