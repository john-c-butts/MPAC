import numpy as np
import os
import pandas as pd
from tqdm import tqdm

# script to generate one hot encoded numpy arrays of test fastas and hyp imp scores from gkmexplain

# define function to get one hot encoded vector      
def onehote(seq):
    seq2=list()
    mapping = {"A":[1, 0, 0, 0], "C": [0, 1, 0, 0], "G": [0, 0, 1, 0], "T":[0, 0, 0, 1]}
    for i in seq:
      seq2.append(mapping[i]  if i in mapping.keys() else [0, 0, 0, 0]) 
    return np.array(seq2)

# define function to make a numpy array of hyp imp scores
def hyp_imp_scores2numpy (path2scores):
    # open hypothetical importance scores
    hyp_scores = pd.read_csv(path2scores,
                            sep = '\t',
                            header=None)
    # 
    list_o_lists = []
    for scores in tqdm(hyp_scores[2]):
        tmp = []
        l0 = []
        l1 = []
        l2 = []
        l3 = []
        for split in scores.split(';'):
            split_split = split.split(',')
            l0.append(float(split_split[0]))
            l1.append(float(split_split[1]))
            l2.append(float(split_split[2]))
            l3.append(float(split_split[3]))
        list_o_lists.append(np.asmatrix([l0,l1,l2,l3]))
    hyp_score_array = np.array(list_o_lists)
    return hyp_score_array

# define function to onehot fastas and hyp imp scores
def onehot_and_save (path2scores,
                     path2fasta,
                     path2fasta_dir,
                     onehot_score_out,
                     onehot_fasta_out):
    # one hot encode importance scores
    one_hot_scores = hyp_imp_scores2numpy(path2scores)
    # one hot encode associated fasta file
    # read in the fasta files and one-hot encode
    fasta_seqs = [x.rstrip() for (i,x) in enumerate(open(path2fasta))
              if i%2==1]

    onehot_fasta = np.array([onehote(x).T
                            for x in fasta_seqs])
    # print dimensions of array
    print(one_hot_scores.shape)
    print(onehot_fasta.shape)
    # change to fasta dir for saving
    os.chdir(path2fasta_dir)
    # save arrays to directory
    # hyp scores
    np.savez(onehot_score_out, 
             one_hot_scores)
    # one hot sequences
    np.savez(onehot_fasta_out,
            onehot_fasta)
# convert importance scores and associated fastas and save as numpy arrays
# K562
onehot_and_save('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/k562_95033/k562_95033_lsgkm_model/predictions/k562_95033_hyp_impscores.txt',
                '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/k562_95033/k562.miss.test.95033.fa',
                '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/k562_95033',
                'k562_95033_miss_hyp_impscores',
                'k562_95033_miss_onehot_fastas')
# HepG2
onehot_and_save('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/hepg2_95033/hepg2_95033_lsgkm_model/predictions/hepg2_95033_hyp_impscores.txt',
                '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/hepg2_95033/hepg2.miss.test.95033.fa',
                '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/hepg2_95033',
                'hepg2_95033_miss_hyp_impscores',
                'hepg2_95033_miss_onehot_fastas')
# SKNSH
onehot_and_save('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/sknsh_95033/sknsh_95033_lsgkm_model/predictions/sknsh_95033_hyp_impscores.txt',
                '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/sknsh_95033/sknsh.miss.test.95033.fa',
                '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/sknsh_95033',
                'sknsh_95033_miss_hyp_impscores',
                'sknsh_95033_miss_onehot_fastas')