import pandas as pd
import os

# script to generate predictions from trained lsgkm models

# define function
def lsgkm_predict(path2model,
                  model,
                  pos_test_fasta,
                  neg_test_fasta):
    # change to model directory
    os.chdir(path2model)
    # make directory for model predictions
    os.mkdir('predictions')
    # store path to model predict
    gkmpredict = '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm-svr/bin/gkmpredict'
    # add model
    model_path = f'{path2model}/{model}'
    # store out name as a variable
    # pos
    out_pos = pos_test_fasta.split('/')[-1].split('.fa')[0]
    # neg
    out_neg = neg_test_fasta.split('/')[-1].split('.fa')[0]
    # make positive outfile
    pos_outfile = f'predictions/{out_pos}.predictions.txt'
    # make negative outfile
    neg_outfile = f'predictions/{out_neg}.predictions.txt'
    # make positive predictions cmd
    pos_cmd = f'{gkmpredict} {pos_test_fasta} {model_path} {pos_outfile}'
    # make negative predicitons cmd
    neg_cmd = f'{gkmpredict} {neg_test_fasta} {model_path} {neg_outfile}'
    # run predictions on positive test set
    os.system(pos_cmd)
    # run predictions on negative test set 
    os.system(neg_cmd)
# generate predictions
# K562
lsgkm_predict('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/k562_95033/k562_95033_lsgkm_model',
              'k562_95033.model.txt',
              '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/k562_95033/k562.miss.test.95033.fa',
              '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/k562_95033/k562.hit.test.95033.fa')
# HepG2
lsgkm_predict('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/hepg2_95033/hepg2_95033_lsgkm_model',
              'hepg2_95033.model.txt',
              '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/hepg2_95033/hepg2.miss.test.95033.fa',
              '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/hepg2_95033/hepg2.hit.test.95033.fa')
# SKNSH
lsgkm_predict('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/sknsh_95033/sknsh_95033_lsgkm_model',
              'sknsh_95033.model.txt',
              '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/sknsh_95033/sknsh.miss.test.95033.fa',
              '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/sknsh_95033/sknsh.hit.test.95033.fa')