# script to plot correlation of models to DNase and MPRA data

#%% import packages
import matplotlib.pyplot as mpl
import seaborn as sns
import pandas as pd
# %%
# define function to open summary text file - we want to drop that first column that doesn't match
def open_summary(path2summary,
               test_set):
    summary = pd.read_csv(path2summary,
                          sep = '\t',
                          usecols = ['feature',
                                      'model',
                                      'cell_type',
                                      'pearson'])
    # add test set column
    summary['test_set'] = [test_set]
    return summary
# %%
h_enf = open_summary('hepg2.max.enformer.dnase.emvar.corr.txt',
                     'DNase')
# %%
# iterate through all DNase summaries and concatenate into single df
all_dnase = pd.concat([open_summary(i,
                                    'DNase') for i in ['../hepg2.enformer.max.dnase.ase.corr.tsv',
                                                       '../hepg2.sei.dnase.ase.enformer.vars.corr.tsv',
                                                       '/Users/buttsj/JAX Dropbox/John Butts/Variant_Effects/Benchmarking_Sets/Vierstra/dnase2mpra/hepg2.mpra.skew.dnase.ase.corr.tsv',
                                                       '../k562.enformer.max.dnase.ase.corr.tsv',
                                                       '../k562.sei.dnase.ase.enformer.vars.corr.tsv',
                                                       '/Users/buttsj/JAX Dropbox/John Butts/Variant_Effects/Benchmarking_Sets/Vierstra/dnase2mpra/k562.mpra.skew.dnase.ase.corr.tsv']])
# %%
# iterate through all MPRA summaries and concatenate into single df
all_mpra = pd.concat([open_summary(i,
                                   'MPRA_emVar') for i in ['hepg2.max.enformer.dnase.emvar.corr.with.controls.txt',
                                                           'hepg2.max.sei.dnase.emvar.corr.with.controls.txt',
                                                           'hepg2.mpra.pred.enformer.emvar.corr.with.controls.txt',
                                                           'k562.max.enformer.dnase.emvar.corr.with.controls.txt',
                                                           'k562.max.sei.dnase.emvar.corr.with.controls.txt',
                                                           'k562.mpra.pred.enformer.emvar.corr.with.controls.txt']])
# %%
# combine summaries and plot
all_summaries = pd.concat([all_dnase, all_mpra])
# %%
# let's try this out
# define a function to reformat dfs
def sloppy_plot (summary_df,
                 model,
                 cell_type):
    # filter full df
    df = summary_df[(summary_df['model'] == model) &
                    (summary_df['cell_type'] == cell_type)].sort_values(by='test_set')
    reformat_df = pd.DataFrame({'model' : [df['model'].tolist()[0]],
                                'cell_type' : [df['cell_type'].tolist()[0]],
                                'pearson_dnase' : [df['pearson'].tolist()[0]],
                                'pearson_mpra' : [df['pearson'].tolist()[1]]})
    return reformat_df
# %%
enf_h = sloppy_plot(all_summaries,
                   'Enformer',
                    'HEPG2')
enf_k = sloppy_plot(all_summaries,
                    'Enformer',
                    'K562')
sei_h = sloppy_plot(all_summaries,
                    'SEI',
                    'HEPG2')
sei_k = sloppy_plot(all_summaries,
                    'SEI',
                    'K562')
mpra_h = sloppy_plot(all_summaries,
                    'MPRA',
                    'HEPG2')
mpra_k = sloppy_plot(all_summaries,
                    'MPRA',
                    'K562')
all_reformat = pd.concat([enf_h,
                          enf_k,
                          sei_h,
                          sei_k,
                          mpra_h,
                          mpra_k])
# %%
# define a dictionary for coloring plot
palette = {'HEPG2' : '#FBB040',
           'K562' : '#00A79D'}
# %%
mpl.figure(dpi=300)
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.plot([0, 1],
                 [0, 1],
                 linestyle='dashed',
                 alpha=.5,
                 color='k')
sns.scatterplot(data=all_reformat,
                x='pearson_dnase',
                y='pearson_mpra',
                hue='cell_type',
                palette=palette,
                style='model',
                s=75,
                linewidth=0)
mpl.xlim(0,1)
mpl.ylim(0,1)
mpl.xlabel('Pearson r: DNase ASE')
mpl.ylabel('Pearson r: MPRA emVar')
sns.despine()
mpl.legend(loc='lower right')
mpl.title('Model Corelations on DNase ASE and MPRA emVars')
#mpl.savefig('dnase_emvar_model_correlation_scatter_v2.pdf')
# %%
print(all_reformat)
# %%
