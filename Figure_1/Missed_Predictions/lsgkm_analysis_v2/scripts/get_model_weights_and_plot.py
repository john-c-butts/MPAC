import matplotlib.pyplot as mpl
import seaborn as sns
import pandas as pd
import os
import numpy as np

# script to generate ecdfs of scores, model weights, and bar plot of top 25 kmers

# define function to generate ecdf plot
def pred_ecdfs(path2preds,
               title,
               outname):
    # change to prediction directory
    os.chdir(path2preds)
    # open hit preds
    hit_preds = [i for i in os.listdir() if 'hit' in i][0]
    # open miss preds
    miss_preds = [i for i in os.listdir() if 'miss' in i][0]
    # open hit preds
    hit_df = pd.read_csv(hit_preds,
                         sep = '\t',
                         header=None)
    # add label
    hit_df['pred_type'] = ['malinois_hit' for i in range(len(hit_df))]
    # open miss preds
    miss_df = pd.read_csv(miss_preds,
                          sep = '\t',
                          header=None)
    # add label
    miss_df['pred_type'] = ['malinois_miss' for i in range(len(miss_df))]
    # combine dfs
    df2plot = pd.concat([hit_df,
                         miss_df])
    # make ecdf
    mpl.figure(dpi=300,
               figsize=(5,7))
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    sns.ecdfplot(data=df2plot,
                 x=1,
                 hue='pred_type',)
    mpl.title(title)
    mpl.savefig(outname)
# generate ecdfs
# K562
pred_ecdfs('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/k562_95033/k562_95033_lsgkm_model/predictions',
           'k562_95033_lsgkm_scores',
           'k562_95033_lsgkm_scores.pdf')
# HepG2
pred_ecdfs('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/hepg2_95033/hepg2_95033_lsgkm_model/predictions',
           'hepg2_95033_lsgkm_scores',
           'hepg2_95033_lsgkm_scores.pdf')
# SKNSH
pred_ecdfs('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/sknsh_95033/sknsh_95033_lsgkm_model/predictions',
           'sknsh_95033_lsgkm_scores',
           'sknsh_95033_lsgkm_scores.pdf')
# define a function to generate all kmers and run
def get_model_weights (path2model,
                       path2preds):
    # store fasta with all 11-mers as a variable
    all_kmers = '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/data/all_11mers.fasta '
    # store path to model predict
    gkmpredict = '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm-svr/bin/gkmpredict'
    # change to predictions directory
    os.chdir(path2preds)
    # get model name from path2model for output
    out_name = path2model.split('/')[-1].split('.')[0]
    outfile = f'{out_name}_model_weights.txt'
    # generate command
    pred_cmd = f'{gkmpredict} {all_kmers} {path2model} {outfile}'
    os.system(pred_cmd)
# K562
get_model_weights('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/k562_95033/k562_95033_lsgkm_model/k562_95033.model.txt',
                  '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/k562_95033/k562_95033_lsgkm_model/predictions')
# HepG2
get_model_weights('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/hepg2_95033/hepg2_95033_lsgkm_model/hepg2_95033.model.txt',
                  '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/hepg2_95033/hepg2_95033_lsgkm_model/predictions')
# SKNSH
get_model_weights('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/sknsh_95033/sknsh_95033_lsgkm_model/sknsh_95033.model.txt',
                  '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/sknsh_95033/sknsh_95033_lsgkm_model/predictions')
# define function to plot top 25 kmers for each model
def plot_top25_kmers (path2preds,
                      path2weights,
                      title,
                      color,
                      outname):
    # change to predictions directory
    os.chdir(path2preds)
    # open model weights
    weights = pd.read_csv(path2weights,
                          sep = '\t',
                          header=None)
    # plot top 25 kmers
    mpl.figure(dpi=300,
               figsize=(5,5))
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    sns.catplot(data=weights[weights[1] > np.quantile(weights[1], .99)].sort_values(by=1, ascending=False).head(25),
                x=0,
                y=1,
                kind='bar',
                color=color)
    mpl.xticks(rotation=90)
    mpl.title(title)
    mpl.tight_layout()
    mpl.savefig(outname)
# plot top 25 kmers
# K562
plot_top25_kmers('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/k562_95033/k562_95033_lsgkm_model/predictions/',
                 '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/k562_95033/k562_95033_lsgkm_model/predictions/k562_95033_model_weights.txt',
                 'lsgkm K562_95033 Top 25 kmers',
                 '#00A79D',
                 'lsgkm K562_95033 Top 25 kmers.pdf')
# HepG2
plot_top25_kmers('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/hepg2_95033/hepg2_95033_lsgkm_model/predictions/',
                 '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/hepg2_95033/hepg2_95033_lsgkm_model/predictions/hepg2_95033_model_weights.txt',
                 'lsgkm HepG2_95033 Top 25 kmers',
                 '#FBB040',
                 'lsgkm HepG2_95033 Top 25 kmers.pdf')
# SKNSH
plot_top25_kmers('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/sknsh_95033/sknsh_95033_lsgkm_model/predictions/',
                 '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/sknsh_95033/sknsh_95033_lsgkm_model/predictions/sknsh_95033_model_weights.txt',
                 'lsgkm SKNSH_95033 Top 25 kmers',
                 '#ED1C24',
                 'lsgkm SKNSH_95033 Top 25 kmers.pdf')
