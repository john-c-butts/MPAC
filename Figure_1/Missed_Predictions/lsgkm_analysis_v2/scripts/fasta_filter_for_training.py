#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  6 10:42:11 2024

@author: buttsj
"""

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
import random
import os
# open missed predictions text file
missed_preds = pd.read_csv('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/data/all_emvars_missed_predictions_only.tsv',
                           sep = '\t')
# get a random number for adding to path for saving files
rand_int = str(random.randint(1, 100000))
# define function to return a training set and test set of positive and negative sequences for lsgkm training by cell type
def ls_gkm_fastas (missed_preds_df,
                   cell_type,
                   cell_type_list,
                   empirical_emvars,
                   allFasta,
                   percent_training,
                   miss_train_fasta,
                   miss_test_fasta,
                   hit_train_fasta,
                   hit_test_fasta):
    # change to directory
    os.chdir('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2')
    # make a directory for storing fastas
    out_dir = f'{cell_type}_{rand_int}'
    #print(out_dir)
    os.mkdir(out_dir)
    # define path to output
    out_path = f'/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_training_sets/{out_dir}'
    # change to output path
    os.chdir(out_dir)
    # filter for only misses
    miss_ids = missed_preds_df[missed_preds_df['missed_cell_type'].isin(cell_type_list)]['hg19_id'].tolist()
    # add alt IDs
    miss_ids_alt = [idee.replace(':R:', ':A:') for idee in miss_ids]
    # combine for all miss ids
    all_miss_ids = miss_ids + miss_ids_alt
    # get sequences for emvar misses
    # ukbb
    all_miss_seqs = SeqIO.parse(allFasta, 'fasta')
    miss_id = []
    miss_seqs = []
    for seq in tqdm(all_miss_seqs):
        if seq.id in all_miss_ids:
            miss_id.append(seq.id)
            miss_seqs.append(str(seq.seq))
        else:
            continue
    miss_df = pd.DataFrame({'hg19_id' : miss_id,
                                'seq' : miss_seqs}).sort_values(by='hg19_id').drop_duplicates()
    # get all emvar sequences and convert to df
    # get emvar IDs
    # refs
    emvar_ID_ref = list(pd.read_csv(empirical_emvars,
                                    sep = '\t',
                                    low_memory=False)['ID'].unique())
    # alts
    emvar_ID_alt = [idee.replace(':R:', ':A:') for idee in emvar_ID_ref]
    # combine ref alt
    emvar_IDs = emvar_ID_alt + emvar_ID_ref
    # filter for only snps, autosomes, and drop miss IDs
    #emvar_snp_IDs = [idee for idee in emvar_IDs if len(idee.split(':')[2]) == 1 and len(idee.split(':')[3]) == 1 and idee.split(':')[0] != 'X']
    # drop miss IDs
    emvar_hit_IDs = [i for i in emvar_IDs if i not in miss_df['hg19_id'].tolist()]
    # get sequences for emvar hits
    # ukbb
    all_hit_seqs = SeqIO.parse(allFasta, 'fasta')
    hit_id = []
    hit_seqs = []
    for seq in tqdm(all_hit_seqs):
        if seq.id in emvar_hit_IDs:
            hit_id.append(seq.id)
            hit_seqs.append(str(seq.seq))
        else:
            continue
    hit_df = pd.DataFrame({'hg19_id' : hit_id,
                                'seq' : hit_seqs}).sort_values(by='hg19_id').drop_duplicates()
    # divide the hits and misses into hold out and training sets
    # train
    miss_train = miss_df.sample(frac=percent_training)
    miss_test = miss_df[~miss_df['hg19_id'].isin(miss_train['hg19_id'].tolist())]
    # hits
    hit_train = hit_df.sample(len(miss_train))
    # hit test
    hit_test = hit_df[~hit_df['hg19_id'].isin(hit_train['hg19_id'])].sample(len(miss_test))
    # write those dfs to fasta
    # miss training
    with open(f'{miss_train_fasta}.{rand_int}.fa', 'w') as f:
        for k, v in zip(miss_train['hg19_id'],
                        miss_train['seq']):
            f.write(f'>{k}\n{v}\n')
    # miss test
    with open(f'{miss_test_fasta}.{rand_int}.fa', 'w') as f:
        for k, v in zip(miss_test['hg19_id'],
                        miss_test['seq']):
            f.write(f'>{k}\n{v}\n')
    # hits training
    with open(f'{hit_train_fasta}.{rand_int}.fa', 'w') as f:
        for k, v in zip(hit_train['hg19_id'],
                        hit_train['seq']):
            f.write(f'>{k}\n{v}\n')
    # hits test
    with open(f'{hit_test_fasta}.{rand_int}.fa', 'w') as f:
        for k, v in zip(hit_test['hg19_id'],
                        hit_test['seq']):
            f.write(f'>{k}\n{v}\n')
# generate training sets for each cell type
# K562
ls_gkm_fastas(missed_preds,
              'k562',
              ['K562', 'K562:HEPG2', 'K562:SKNSH', 'K562:HEPG2:SKNSH'],
              '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/data/k562_emvars_with_controls.tsv',
              '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/data/all_OL_uncollapsed.fasta',
              .8,
              'k562.miss.80.percent.train',
              'k562.miss.test',
              'k562.hit.80.percent.train',
              'k562.hit.test')
# HepG2
ls_gkm_fastas(missed_preds,
              'hepg2',
              ['HEPG2', 'K562:HEPG2', 'HEPG2:SKNSH', 'K562:HEPG2:SKNSH'],
              '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/data/hepg2_emvars_with_controls.tsv',
              '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/data/all_OL_uncollapsed.fasta',
              .8,
              'hepg2.miss.80.percent.train',
              'hepg2.miss.test',
              'hepg2.hit.80.percent.train',
              'hepg2.hit.test')
# SKNSH
ls_gkm_fastas(missed_preds,
              'sknsh',
              ['SKNSH', 'K562:SKNSH', 'HEPG2:SKNSH', 'K562:HEPG2:SKNSH'],
              '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/data/sknsh_emvars_with_controls.tsv',
              '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/data/all_OL_uncollapsed.fasta',
              .8,
              'sknsh.miss.80.percent.train',
              'sknsh.miss.test',
              'sknsh.hit.80.percent.train',
              'sknsh.hit.test')