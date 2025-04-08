# script to compare precomputed 1KG Enformer preds to MPRA emvars

#%% import packages
import matplotlib.pyplot as mpl
import seaborn as sns
import pandas as pd
from scipy import stats
# %%
# define function to open enformer preds, mpra data and correlate
def enformer_dnase2mpra_corr (path2enformer,
                              path2mpra,
                              cell_type,
                              annot):
    ### open and merge data ###
    # open mpra
    mpra = pd.read_csv(path2mpra,
                       sep = '\t',
                       low_memory=False).drop_duplicates(subset='hg38_id')
    # filter for only skew
    mpra_filtered = mpra.filter(['hg38_id',
                                 'Log2Skew'])
    # open enformer preds
    enformer = pd.read_csv(path2enformer,
                           sep = '\t')
    #enformer_cols = [i for i in enformer.keys() if i != 'hg19_id']
    # drop that hg19 id
    #enformer = enformer.filter(enformer_cols)
    # merge preds on id
    merge_preds = mpra_filtered.merge(enformer, left_on='hg38_id', right_on='hg38_id')
    # print len of merge to have idea of how many precomputed values are there
    print(f'n merge variants: {len(merge_preds)}')
    ### correlation bit ###
    corr = merge_preds.corr(numeric_only=True)
    round(corr, 2)
    # plot heatmap
    mpl.figure(dpi=300)
    mpl.title(f'Enformer {cell_type.upper()} DNase vs. emVars')
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    sns.heatmap(corr,
            annot=annot,
            vmin=0,
            vmax=1,
            cmap='mako')
    return corr, merge_preds

# %%
hepg2_enf_corr, hepg2_merge = enformer_dnase2mpra_corr('/Users/buttsj/Dropbox (JAX)/Variant_Effects/Benchmarking_Sets/Vierstra/dnase2mpra/all_hepg2_emvar_snps_enformer_dnase_preds.tsv',
                                                            '/Users/buttsj/Dropbox (JAX)/Variant_Effects/ukbb_gtex_mpra/final_datasets/hepg2.emvar.merge.017.preds.with.controls.tsv',
                                                            'HepG2',
                                                             True)
#%%
# get p-value of max correlation for reporting
stats.pearsonr(x=hepg2_merge['Log2Skew'],
               y=hepg2_merge['ENCFF577SOF'])
# %%
k562_enf_corr, k562_merge = enformer_dnase2mpra_corr('/Users/buttsj/Dropbox (JAX)/Variant_Effects/Benchmarking_Sets/Vierstra/dnase2mpra/all_k562_emvar_snps_enformer_dnase_preds.tsv',
                                                          '/Users/buttsj/Dropbox (JAX)/Variant_Effects/Benchmarking_Sets/Vierstra/dnase2mpra/all.gtex.traits.mpra.blacklist.filtered.k562.emVars.txt',
                                                          'K562',
                                                          True)
#%%
# get p-value of max correlation for reporting
stats.pearsonr(x=k562_merge['Log2Skew'],
               y=k562_merge['ENCFF868NHV'])
# %%
# define function to save txt file of best correlation from corr df for making plot later
def save_max_corr (corr_df,
                   model,
                   cell_type):
        # drop log2skew row
        no_skew = corr_df[corr_df.index != 'Log2Skew'].filter(['Log2Skew'])
        no_skew['feature'] = no_skew.index
        max_skew = no_skew[no_skew['Log2Skew'] ==  max(no_skew['Log2Skew'])]
        max_skew['model'] = [model]
        max_skew['cell_type'] = [cell_type]
        max_skew['pearson'] = [max_skew['Log2Skew'].tolist()[0]]
        return max_skew
# %%
# get max correlation for K562 and save
k562_max_corr = save_max_corr(k562_enf_corr,
                              'Enformer',
                              'K562')
k562_max_corr.to_csv('k562.max.enformer.dnase.emvar.corr.with.controls.txt',
                     sep = '\t',
                     index = False)
# %%
# get max correlation for HepG2 and save
hepg2_max_corr = save_max_corr(hepg2_enf_corr,
                               'Enformer',
                               'HEPG2')
hepg2_max_corr.to_csv('hepg2.max.enformer.dnase.emvar.corr.with.controls.txt',
                      sep = '\t',
                      index = False)
# %% filter for only common variants from Enformer Predictions and save VCF for doing other predictions
# k562
#k562_enf_emvars.to_csv('/Users/buttsj/Dropbox (JAX)/Variant_Effects/Benchmarking_Sets/Vierstra/dnase2mpra/enformer.common.gtex.traits.mpra.blacklist.filtered.k562.emVars.txt',
#                        sep = '\t',
#                        index = False)
# hepg2
#hepg2_enf_emvars.to_csv('/Users/buttsj/Dropbox (JAX)/Variant_Effects/Benchmarking_Sets/Vierstra/dnase2mpra/enformer.common.gtex.traits.mpra.blacklist.filtered.hepg2.emVars.txt',
#                        sep = '\t',
#                        index = False)
# %% save enformer emvars to disk as vcf and re-run sei predictions
#k562_enf_vcf = pd.DataFrame({'#CHROM' : [i.split(':')[0] for i in k562_enf_emvars['hg38_id']],
#                             'POS' : [i.split(':')[1] for i in k562_enf_emvars['hg38_id']],
#                             'ID' : ['.' for i in range(len(k562_enf_emvars))],
#                             'REF' : [i.split(':')[2] for i in k562_enf_emvars['hg38_id']],
#                             'ALT' : [i.split(':')[-1] for i in k562_enf_emvars['hg38_id']],
#                             'FILTER' : ['.' for i in range(len(k562_enf_emvars))],
#                             'INFO' : ['k562_enformer_common_emvar' for i in range(len(k562_enf_emvars))]})
#hepg2_enf_vcf = pd.DataFrame({'#CHROM' : [i.split(':')[0] for i in hepg2_enf_emvars['hg38_id']],
#                              'POS' : [i.split(':')[1] for i in hepg2_enf_emvars['hg38_id']],
#                              'ID' : ['.' for i in range(len(hepg2_enf_emvars))],
#                              'REF' : [i.split(':')[2] for i in hepg2_enf_emvars['hg38_id']],
#                              'ALT' : [i.split(':')[-1] for i in hepg2_enf_emvars['hg38_id']],
#                              'FILTER' : ['.' for i in range(len(hepg2_enf_emvars))],
#                              'INFO' : ['.' for i in range(len(hepg2_enf_emvars))]})
#k562_enf_vcf.to_csv('k562.enformer.common.var.emvar.hits.vcf',
#                    sep = '\t',
#                    index = False)
#hepg2_enf_vcf.to_csv('hepg2.enformer.common.var.emvar.hits.vcf',
#                     sep = '\t',
#                     index = False)
# %%
