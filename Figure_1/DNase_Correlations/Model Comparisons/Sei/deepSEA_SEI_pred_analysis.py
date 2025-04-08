# script to generate scatterplots of SEI predictions dnase ASE variants for benchmarking
# %%
# import packages
import pandas as pd
import matplotlib.pyplot as mpl
import seaborn as sns
from scipy import stats
# %%
# open sei predictions and figure out what columns to use for analysis
# hepg2
hepg2_sei = pd.read_csv('/Users/buttsj/Dropbox (JAX)/Variant_Effects/Kipoi_Comparisons/vierstra_cav/kipoi_preds/deepSEA_web_preds/sei_preds/5b3b002f-0002-40d3-aec7-8a789787c482/5b3b002f-0002-40d3-aec7-8a789787c482_hepg2_ase_min_fdr_05_diffs.tsv',
                        sep = '\t')
hepg2_dnase_cols = [i for i in hepg2_sei.keys() if 'HepG2' in i and 'DNase' in i]
# k562
k562_sei = pd.read_csv('/Users/buttsj/Dropbox (JAX)/Variant_Effects/Kipoi_Comparisons/vierstra_cav/kipoi_preds/deepSEA_web_preds/sei_preds/e14814a9-f507-41fa-832d-d249c091e23d/e14814a9-f507-41fa-832d-d249c091e23d_k562_ase_min_fdr_05_diffs.tsv',
                       sep = '\t')
k562_dnase_cols = [i for i in k562_sei.keys() if 'K562' in i and 'DNase' in i]
# %%
# define function to plot scatters and return a df with correlation for summary bar plot
def sei_analysis (sei_preds,
                  sei_cols2use,
                  path2merge_dnase,
                  path2enformer_ids,
                  cell_type,
                  emp_col,
                  plot_color,
                  figsize,
                  annot,
                  annot_x,
                  annot_y):
    # open empirical data
    emp_data = pd.read_csv(path2merge_dnase,
                           sep = '\t')
    # convert id column to chr:pos:ref:alt for matching to empirical
    sei_preds['id'] = [(':').join([chrom,
                                   str(pos),
                                   ref,
                                   alt]) for chrom,
                                             pos,
                                             ref,
                                             alt in zip(sei_preds['chrom'],
                                                        sei_preds['position'],
                                                        sei_preds['ref_allele'],
                                                        sei_preds['alt_allele'])]
    # filter sei preds for relevant dnase columns
    sei_filter = sei_preds.filter(['id'] + sei_cols2use)
    # merge sei preds and empirical data
    merge_df = sei_filter.merge(emp_data,
                                on = 'id')
    # flip sign of DNase data to alt - ref
    merge_df[emp_col] = [-1 * i for i in merge_df[emp_col]]
    # filter for only enformer vars
    merge_df = merge_df[merge_df['id'].isin(pd.read_csv(path2enformer_ids, sep = '\t')['hg38_id'].tolist())]
    print(len(merge_df))
    # make a df for generating a heatmap 
    heat_cols = [emp_col] + sei_cols2use
    heat_df = merge_df.filter(heat_cols)
    heat_corr = heat_df.corr(numeric_only=True)
    round(heat_corr, 2)
    mpl.figure(dpi=300, figsize=figsize)
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    sns.heatmap(heat_corr,
                annot=annot,
                vmin=0,
                vmax=1,
                xticklabels=annot_x,
                yticklabels=annot_y)
    mpl.title(f'{cell_type.upper()} SEI DNase ASE Correlation')
    mpl.savefig(f'sei_corr_plots/{cell_type}_sei_dnase_corr_heatmap_enf_common.pdf')
    # make an empty dictionary for storing correlation info
    correlation_dict = {}
    # iterate through list of columns and get correlations
    for col in sei_cols2use:
        # get output string for col
        col_out = ('_').join(col.split('|'))
        # get correlation
        corr, pval = stats.pearsonr(merge_df[emp_col],
                                    merge_df[col])
        # round correlations to two places
        corr = round(corr, 2)
        pval = round(pval, 2)
        # update correlation_dict
        correlation_dict.update({col : (corr, pval)})
        # plot
        pt = (0, 0)
        mpl.figure(dpi=300)
        mpl.rcParams['pdf.fonttype'] = 42
        mpl.rcParams['ps.fonttype'] = 42
        ax = sns.scatterplot(data=merge_df,
                             x=emp_col,
                             y=col,
                             color=plot_color,
                             linewidth=0,
                             s=50)
        ax.axline(pt, 
                  slope=1, 
                  color='black',
                  alpha=.25,
                  linestyle='dashed') 
        mpl.xlabel('DNase Allelic Skew')
        mpl.ylabel(col)
        mpl.title(f'{cell_type.upper()} DNase ASE Correlation | r: {corr}, {pval}')
        #mpl.savefig(f'/Users/buttsj/Dropbox (JAX)/Variant_Effects/Kipoi_Comparisons/vierstra_cav/sei_corr_plots/{cell_type}_{col_out}_dnase_corr.pdf')
        mpl.show()
        mpl.close()
    return heat_corr, merge_df;
    # make a dataframe of all correlations
    #corr_df = pd.DataFrame({'prediction' : [i for i in correlation_dict.keys()],
    #                        'pearson_r' : [correlation_dict.get(i)[0] for i in correlation_dict.keys()],
    #                        'p_val' : [correlation_dict.get(i)[-1] for i in correlation_dict.keys()],
    #                        'Model' : ['SEI' for i in range(len(correlation_dict.keys()))],
    #                        'cell_type' : [cell_type.upper() for i in range(len(correlation_dict.keys()))]})
    #print(corr_df)
    # save to disk
    #corr_df.to_csv(f'{cell_type}_sei_empirical_dnase_correlation.txt',
    #               sep = '\t',
    #               index = False)
# %%
# hepg2
h_corr, h_merge = sei_analysis(hepg2_sei, 
             hepg2_dnase_cols, 
             '/Users/buttsj/Dropbox (JAX)/Variant_Effects/Kipoi_Comparisons/vierstra_cav/hepg2_dnase_cav_merged_017_avg.txt', 
             '/Users/buttsj/Dropbox (JAX)/Variant_Effects/Benchmarking_Sets/Vierstra/enformer_comparisons/hepg2.enformer.common.dnase.ase.var.ids.tsv',
             'hepg2', 
             'es_mean', 
             'gray',
             (8.5, 5),
             True,
             False,
             True)
# %%
# k562
k_corr, k_merge = sei_analysis(k562_sei,
             k562_dnase_cols,
             '/Users/buttsj/Dropbox (JAX)/Variant_Effects/Kipoi_Comparisons/vierstra_cav/k562_dnase_cav_merged_017_avg.txt',
             '/Users/buttsj/Dropbox (JAX)/Variant_Effects/Benchmarking_Sets/Vierstra/enformer_comparisons/k562.enformer.common.dnase.ase.var.ids.tsv',
             'k562',
             'es_mean',
             'gray',
             (8.5,11),
             False,
             False,
             True)
# %%
# define function to save max correlation to disk for plotting later
def save_max_corr (corr_df,
                   model,
                   cell_type):
        # drop log2skew row
        no_skew = corr_df[corr_df.index != 'es_mean'].filter(['es_mean'])
        no_skew['feature'] = no_skew.index
        max_skew = no_skew[no_skew['es_mean'] ==  max(no_skew['es_mean'])]
        max_skew['model'] = [model]
        max_skew['cell_type'] = [cell_type]
        max_skew['pearson'] = [max_skew['es_mean'].tolist()[0]]
        return max_skew
# %%
k_max_corr = save_max_corr(k_corr, 'SEI', 'K562')
k_max_corr.to_csv('/Users/buttsj/Dropbox (JAX)/Variant_Effects/Benchmarking_Sets/Vierstra/k562.sei.dnase.ase.enformer.vars.corr.tsv',
                  sep = '\t',
                  index = False)
# %%
h_max_corr = save_max_corr(h_corr, 'SEI', 'HEPG2')
h_max_corr.to_csv('/Users/buttsj/Dropbox (JAX)/Variant_Effects/Benchmarking_Sets/Vierstra/hepg2.sei.dnase.ase.enformer.vars.corr.tsv',
                  sep = '\t',
                  index = False)
# %%
# calculate mpac correlation to enformer/sei dnase subset and save for plotting
path2save = '/Users/buttsj/JAX Dropbox/John Butts/Variant_Effects/Benchmarking_Sets/Vierstra/dnase2mpra/'
h_mpac_corr, h_p = stats.pearsonr(h_merge['es_mean'],
                                  h_merge['hepg2_skew_pred'])
h_mpac2save = pd.DataFrame({'feature' : ['hepg2_skew_pred'],
                            'model' : ['MPRA'],
                            'cell_type' : ['HEPG2'],
                            'pearson' : [h_mpac_corr],
                            'test_set' : ['DNase']})
# %%
k_mpac_corr, k_p = stats.pearsonr(k_merge['es_mean'],
                                  k_merge['k562_skew_pred'])
k_mpac2save = pd.DataFrame({'feature' : ['k562_skew_pred'],
                            'model' : ['MPRA'],
                            'cell_type' : ['K562'],
                            'pearson' : [k_mpac_corr],
                            'test_set' : ['DNase']})
# %%
# save mpac correlations to disk
h_mpac2save.to_csv(f'{path2save}/hepg2.mpra.skew.dnase.ase.corr.tsv', sep = '\t', index = False)
# %%
k_mpac2save.to_csv(f'{path2save}/k562.mpra.skew.dnase.ase.corr.tsv', sep = '\t', index = False)
# %%
