#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 19:40:32 2023

@author: buttsj
"""

# script to annotate cosmic predictions
#%%
# import packages
import pandas as pd
import os
import pybedtools
from tqdm import tqdm
#%%
# open predictions
cosmic_preds = pd.read_csv('/Users/buttsj/Dropbox (JAX)/Variant_Effects/cosmic/all.cosmic.v98.autosome.snps.recurrence.017.vcf',
                           sep = '\t')
# open BED file of cosmic predictions
cosmic_bed = pybedtools.BedTool('/Users/buttsj/Dropbox (JAX)/Variant_Effects/cosmic/final_bed_files/all.cosmic.v98.autosome.snps.017.bed')
# open promoter BED files
# 10kb Exon Filtered BED
#tenKB_promoter_bed = pybedtools.BedTool('/Users/buttsj/Dropbox (JAX)/Variant_Effects/cosmic/final_bed_files/gencode.v44.protein.coding.10kb.promoters.autosomes.exon.filtered.bed')
# 5kb Exon Filtered BED
#fiveKB_promoter_bed = pybedtools.BedTool('/Users/buttsj/Dropbox (JAX)/Variant_Effects/cosmic/final_bed_files/gencode.v44.protein.coding.5kb.promoters.autosomes.exon.filtered.bed')
# 1kb Exon Filtered BED
#oneKB_promoter_bed = pybedtools.BedTool('/Users/buttsj/JAX Dropbox/John Butts/Variant_Effects/cosmic/final_bed_files/')
# 500bp Exon Filtered BED
#fiveHundredBP_promoter_bed = pybedtools.BedTool('/Users/buttsj/Dropbox (JAX)/Variant_Effects/cosmic/final_bed_files/gencode.v44.protein.coding.500bp.promoters.autosomes.exon.filtered.bed')
# open 250 bp promoter BED file - 11/26
twoFifty_pro_bed = pybedtools.BedTool('/Users/buttsj/JAX Dropbox/John Butts/Variant_Effects/cosmic/final_bed_files/gencode.v44.protein.coding.canonical.autosomes.0.based.250bp.exon.filtered.bed')
# open ENCODE cCRE BED files
#pELS_bed = pybedtools.BedTool('/Users/buttsj/Dropbox (JAX)/Variant_Effects/encode_ccre_beds/GRCh38-pELS.V4.bed.gz')
#%%# open Muelman DHS tsv
mueleman_dhs = pd.read_csv('/Users/buttsj/Downloads/ENCFF503GCK.tsv',
                           sep = '\t',
                           low_memory=False)
# fiilter muelman dhs and convert to BED file
mueleman_bed = pybedtools.BedTool.from_dataframe(mueleman_dhs.filter(['seqname',
                                                                      'start',
                                                                      'end',
                                                                      'identifier']))
#%%
# define function to add annotation to predictions DF using BED intersect
def add_variate_column (preds_df, bed_file, col_name):
    # COSMIC BED file contains COSV based ID in 'INFO' column - match to DF with that ID #
    # get intersect of input BED file and cosmic BED file
    intersect_bed_ids = [i for i in pd.Series(cosmic_bed.intersect(bed_file, wa=True).to_dataframe()['name'].tolist()).unique()]
    # make a dictionary of intersect cosmic IDs for annotating column
    intersect_dict = dict(zip(intersect_bed_ids, ['+' for i in range(len(intersect_bed_ids))]))
    # add variate to preds_df
    preds_df.loc[:,][col_name] = [intersect_dict.get(i) if i in intersect_dict.keys() else '-' for i in tqdm(preds_df['id'])]
    return preds_df
#%%
# add annotations to prediction df
# 500bp promoter
#cosmic_preds_annotated = add_variate_column(cosmic_preds, fiveHundredBP_promoter_bed, '500bp_promoter')
# 1kb promoter 
#cosmic_preds_annotated = add_variate_column(cosmic_preds_annotated, oneKB_promoter_bed, '1kb_promoter')
# 5kb promoter
#cosmic_preds_annotated = add_variate_column(cosmic_preds_annotated, fiveKB_promoter_bed, '5kb_promoter')
# 10kb promoter
#cosmic_preds_annotated = add_variate_column(cosmic_preds_annotated, tenKB_promoter_bed, '10kb_promoter')
#%%
# add 250 bp promoter column
cosmic_preds_annotated = add_variate_column(cosmic_preds, twoFifty_pro_bed, '250bp_promoter')
# add DHS column
cosmic_preds_annotated = add_variate_column(cosmic_preds_annotated, mueleman_bed, 'DHS_meuleman')
# open list of promoters from Cornell non-coding driver database
cnc_noncoding_promoters = [i for i in pd.read_csv('/Users/buttsj/Dropbox (JAX)/Variant_Effects/cosmic/cnc.noncoding.driver.database.all.promoters.csv')['Gene Name'].unique()]
#%%
# subset promoter bed files on cnc promoters and add variates to DF
# make a list of promoter BED files
promoter_bed_dict = dict(zip([#'10kb_promoter',
                              #'5kb_promoter',
                              #'1kb_promoter',
                              #'500bp_promoter',
                              '250bp_promoter'],[#tenKB_promoter_bed,
                                                 #fiveKB_promoter_bed,
                                                 #oneKB_promoter_bed,
                                                 #fiveHundredBP_promoter_bed,
                                                 twoFifty_pro_bed]))
#%%
# define function to take list of cancer promoters 
def make_cancer_gene_promoter_beds (pro_bed_dict, promoter_list, cancer_gene_list):
    # make a dictionary for storing DF
    promoter_dict = {}
    for i in pro_bed_dict.keys():
        # convert from BED
        bed_df = pro_bed_dict.get(i).to_dataframe()
        # filter bed_df for cancer promoters
        promoter_bed = bed_df[bed_df['name'].isin(promoter_list)]
        # update dictionary with 'promoter' : BED file pair
        promoter_dict.update({f'{cancer_gene_list}_{i}' : pybedtools.BedTool.from_dataframe(promoter_bed)})
    return promoter_dict
#%%
# make a dictionary of cnc promoter beds
cnc_promoter_beds = make_cancer_gene_promoter_beds(promoter_bed_dict, 
                                                   cnc_noncoding_promoters, 
                                                   'cnc_database')
#%%
# add columns for each CNC promoter BED file
for i in cnc_promoter_beds.keys():
    add_variate_column(cosmic_preds_annotated, cnc_promoter_beds.get(i), i)
#%%
# open Weinhold promoter list
weinhold_promoter_list = [i for i in pd.read_csv('/Users/buttsj/JAX Dropbox/John Butts/Variant_Effects/cosmic/Weinhold_2014_Promtoer_Hotspots.csv')['gene_symbol'].unique()]
# make BED files of Weinhold promoters
weinhold_promoter_beds = make_cancer_gene_promoter_beds(promoter_bed_dict, 
                                                        weinhold_promoter_list, 
                                                        'weinhold_2014')
#%%
# add columns for each CNC promoter BED file
for i in weinhold_promoter_beds.keys():
    add_variate_column(cosmic_preds_annotated, weinhold_promoter_beds.get(i), i)

# add binary for recurrence
cosmic_preds_annotated.loc[:,]['recurrent'] = ['+' if i > 1 else '-' for i in cosmic_preds_annotated['recurrence']]
#%%
### annotations as of Oct 9, 2023 ### 
# save annotated predictions to file
#cosmic_preds_annotated.to_csv('/Users/buttsj/Dropbox (JAX)/Variant_Effects/cosmic/all.cosmic.v98.autosome.snps.annotated.017.txt',
#                              sep = '\t',
#                              index = False)
### annotations as of Nov 26, 2023 ### 
# Added 250bp promoter
# save annotated predictions to file
# added meuleman DHS intersection 3/19/2025
cosmic_preds_annotated.to_csv('/Users/buttsj/Dropbox (JAX)/Variant_Effects/cosmic/all.cosmic.v98.autosome.snps.annotated.017.031925.txt',
                              sep = '\t',
                              index = False)

# %%
