# script to make heatmaps of enformer 1kg SAD scores vs. DNase ASE skew

# %%
# import libraries
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as mpl
import numpy as np
from scipy import stats

# %%
# define function to open DNase Data, open enformer data, and generate heatmap of correlation of enformer predictions to DNase ASE
# !!! also this defines the comparison set for SEI and MPRA - save list of ids to disk
def enformer_dnase_heatmap (path2dnase,
                            path2enformer,
                            cell_type,
                            heat_out,
                            annot):
    # open dnase data
    data = pd.read_csv(path2dnase,
                       sep = '\t')
    # add refactored es_mean and es_weighed_mean for comparison to
    # es_mean 
    data['refactor_es_mean'] = [-1 * i for i in data['es_mean']]
    # es_weighted_mean
    data['DNase Skew'] = [-1 * i for i in data['es_weighted_mean']]
    # return a list of dnase positives
    dnase_ase_vars = data[data['min_fdr'] <= .05]['id'].tolist()
    # filter full dnase for ase vars
    dnase_ase_data = data[data['min_fdr'] <= .05]
    # open enformer data
    enformer_data = pd.read_csv(path2enformer,
                                sep = '\t')
    # merge data
    all_merge = dnase_ase_data.merge(enformer_data, on = 'id')
    print(len(all_merge))
    ### get list of IDs that are in commmon variants ###
    enformer_ids = all_merge['id'].tolist()
    # get heatmapping
    if cell_type.upper() == 'K562':
        cols2use = ['id',
                    'refactor_es_mean',
                    #'k562_skew_pred',
					'ENCFF899YDP',
					'ENCFF515UNC',
					'ENCFF708UIS',
					'ENCFF413AHU',
					'ENCFF868NHV',
					'ENCFF565YDB',
					'ENCFF971AHO']
        # filter df for corr
        df2corr = all_merge.filter(cols2use)
    elif cell_type.upper() == 'HEPG2':
        cols2use = ['id',
                    'refactor_es_mean',
                    #'hepg2_skew_pred',
			  		'ENCFF136DBS',
					'ENCFF205TKQ',
					'ENCFF577SOF']
        # filter df for corr
        df2corr = all_merge.filter(cols2use)
    # correlate df
    corr = df2corr.corr(numeric_only=True)
    round(corr,2)
    # plot heatmap
    mpl.figure(dpi=300)
    mpl.title(f'Enformer {cell_type.upper()} vs. DNase ASE')
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    sns.heatmap(corr,
            annot=annot,
            vmin=0,
            vmax=1,
            xticklabels=False,
            yticklabels=True)
    mpl.savefig(heat_out)
    mpl.show()
    mpl.close()
    # make scatter plot of best correlation
    # hepg2
    if cell_type.upper() == 'HEPG2':
        pear_corr, pval = stats.pearsonr(all_merge['ENCFF577SOF'],
                                     all_merge['refactor_es_mean'])
        pear_corr = round(pear_corr, 2)
        pval = round(pval, 2)
        pt = (0, 0)
        mpl.figure(dpi=300)
        mpl.rcParams['pdf.fonttype'] = 42
        mpl.rcParams['ps.fonttype'] = 42
        ax = sns.scatterplot(data=all_merge,
                    x= 'refactor_es_mean',
                    y='ENCFF577SOF',
                    color='navy',
                    linewidth=0,
                    s=50)
        ax.axline(pt, 
                  slope=1, 
                  color='black',
                  alpha=.25,
                  linestyle='dashed') 
        mpl.xlabel('DNase Allelic Skew')
        mpl.ylabel('ENCFF577SOF')
        mpl.title(f'Enformer {cell_type.upper()} DNase ASE Correlation | r: {pear_corr}, {pval}')
        mpl.show()
        mpl.close()
    else:
    # k562
        pear_corr, pval = stats.pearsonr(all_merge['ENCFF868NHV'],
                                     all_merge['refactor_es_mean'])
        pear_corr = round(pear_corr, 2)
        pval = round(pval, 2)
        pt = (0, 0)
        mpl.figure(dpi=300)
        mpl.rcParams['pdf.fonttype'] = 42
        mpl.rcParams['ps.fonttype'] = 42
        ax = sns.scatterplot(data=all_merge,
                    x= 'refactor_es_mean',
                    y='ENCFF868NHV',
                    color='navy',
                    linewidth=0,
                    s=50)
        ax.axline(pt, 
                  slope=1, 
                  color='black',
                  alpha=.25,
                  linestyle='dashed') 
        mpl.xlabel('DNase Allelic Skew')
        mpl.ylabel('ENCFF868NHV')
        mpl.title(f'Enformer {cell_type.upper()} DNase ASE Correlation | r: {pear_corr}, {pval}')
        mpl.show()
        mpl.close()
    # plot malinois vs dnase hepg2
    
    return corr, enformer_ids;
# %%
h_corrs, h_ase_enf_ids = enformer_dnase_heatmap ('/Users/buttsj/Dropbox (JAX)/Variant_Effects/Benchmarking_Sets/Vierstra/CAV_017_Preds/hepg2_dnase_cav_merged_017_avg.txt',
                                                 '/Users/buttsj/Dropbox (JAX)/Variant_Effects/Benchmarking_Sets/Vierstra/enformer_comparisons/enformer_precomputed_1kg/all_enformer_1kg_hepg2_dnase_cav_sad.txt',
                                                 'hepg2',
                                                 'hepg2_all_enformer_dnase_cav_correlations.pdf',
                                                 True)
# %%
k_corrs, k_ase_enf_ids = enformer_dnase_heatmap('/Users/buttsj/Dropbox (JAX)/Variant_Effects/Benchmarking_Sets/Vierstra/CAV_017_Preds/k562_dnase_cav_merged_017_avg.txt',
                                                '/Users/buttsj/Dropbox (JAX)/Variant_Effects/Benchmarking_Sets/Vierstra/enformer_comparisons/enformer_precomputed_1kg/all_enformer_1kg_k562_dnase_cav_sad.txt',
                                                'k562',
                                                'k562_all_enformer_dnase_cav_correlations.pdf',
                                                True)
# %%
# define function for filtering correlation df and saving in the style of the other correlation plot dfs
def save_max_corr (corr_df,
                   model,
                   cell_type):
        # drop log2skew row
        no_skew = corr_df[~corr_df.index.isin(['refactor_es_mean','k562_skew_pred'])].filter(['refactor_es_mean'])
        no_skew['feature'] = no_skew.index
        max_skew = no_skew[no_skew['refactor_es_mean'] ==  max(no_skew['refactor_es_mean'])]
        max_skew['model'] = [model]
        max_skew['cell_type'] = [cell_type]
        max_skew['pearson'] = [max_skew['refactor_es_mean'].tolist()[0]]
        return max_skew
# %%
# get max correlation for K562 and save
k_enf_max_corr = save_max_corr(k_corrs, 'Enformer', 'K562')
k_enf_max_corr.to_csv('/Users/buttsj/Dropbox (JAX)/Variant_Effects/Benchmarking_Sets/Vierstra/k562.enformer.max.dnase.ase.corr.tsv',
                      sep = '\t',
                      index = False)
# %%
# get max correlation for HepG2 and save
h_enf_max_corr = save_max_corr(h_corrs, 'Enformer', 'HEPG2')
h_enf_max_corr.to_csv('/Users/buttsj/Dropbox (JAX)/Variant_Effects/Benchmarking_Sets/Vierstra/hepg2.enformer.max.dnase.ase.corr.tsv',
                      sep = '\t',
                      index = False)
# define function for saving malinois correlation to disk
def save_mpra_corr (corr_df,
                    model,
                    cell_type):
    # filter for only dnase corr
    dnase_only = corr_df.filter(['refactor_es_mean'])
    dnase_only['feature'] = dnase_only.index
    # filter for only mpra
    mpra_only = dnase_only[dnase_only['feature'] == f'{cell_type.lower()}_skew_pred']
    mpra_only['cell_type'] = [cell_type.upper()]
    mpra_only['model'] = [model]
    mpra_only['pearson'] = [mpra_only['refactor_es_mean'].tolist()[0]]
    return mpra_only
# %%
# get k562 mpra correlation and save to disk
k_mpra_corr = save_mpra_corr(k_corrs, 'MPRA', 'K562')
k_mpra_corr.to_csv('/Users/buttsj/Dropbox (JAX)/Variant_Effects/Benchmarking_Sets/Vierstra/k562.mpra.skew.dnase.ase.corr.tsv',
                   sep = '\t',
                   index = False)
# %%
h_mpra_corr = save_mpra_corr(h_corrs, 'MPRA', 'HEPG2')
h_mpra_corr.to_csv('/Users/buttsj/Dropbox (JAX)/Variant_Effects/Benchmarking_Sets/Vierstra/hepg2.mpra.skew.dnase.ase.corr.tsv',
                   sep = '\t',
                   index = False)
# %%
# save enformer passing IDs to disk
# k562
pd.DataFrame({'hg38_id' : k_ase_enf_ids}).to_csv('k562.enformer.common.dnase.ase.var.ids.tsv',
                                             sep = '\t',
                                             index = False)
# %%
# save enformer passing IDs to disk
# hepg2
pd.DataFrame({'hg38_id' : h_ase_enf_ids}).to_csv('hepg2.enformer.common.dnase.ase.var.ids.tsv',
                                                 sep = '\t',
                                                 index = False)
# %%
