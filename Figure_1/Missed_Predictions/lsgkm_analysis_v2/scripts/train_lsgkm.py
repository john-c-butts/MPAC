import pandas as pd
import os

# script to train lsgkm on fastas in directory

# define function
def train_lsgkm (fasta_dir,
                 positive_train,
                 negative_train):
    # change to provided directory
    os.chdir(fasta_dir)
    # make output directory
    folder = fasta_dir.split('/')[-1]
    outdir = f'{folder}_lsgkm_model'
    os.mkdir(outdir)
    # show directory working on
    print(f'training lsgkm on fastas in: {os.getcwd()}')
    print('\n')
    # store path to gkmpredict
    gkmtrain = '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm-svr/bin/gkmtrain'
    # make command
    train_cmd = f'{gkmtrain} {positive_train} {negative_train} {outdir}/{folder}'
    os.system(train_cmd)
# Train lsgkm-svm
# K562
train_lsgkm('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/k562_95033',
            'k562.miss.80.percent.train.95033.fa',
            'k562.hit.80.percent.train.95033.fa')
# HepG2
train_lsgkm('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/hepg2_95033',
            'hepg2.miss.80.percent.train.95033.fa',
            'hepg2.hit.80.percent.train.95033.fa')
# SKNSH
train_lsgkm('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/sknsh_95033',
            'sknsh.miss.80.percent.train.95033.fa',
            'sknsh.hit.80.percent.train.95033.fa')

