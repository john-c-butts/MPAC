# script to filter K562, HepG2, and SK-N-SH for only dELS #

# import packages #
import pandas as pd

# open full bed files #
# k562 #
k562_full = pd.read_csv('K562.enhancers.bed', sep = '\t', header = None)
# hepg2 #
hepg2_full = pd.read_csv('HepG2.enhancers.bed', sep = '\t', header = None)
# sknsh #
sknsh_full = pd.read_csv('SK-N-SH.enhancers.bed', sep = '\t', header = None)
# get list of all 'dELS-related' annotations #
dels2filter = ['dELS', 'dELS,CTCF-bound']
# filter full beds for only dels #
# k562 #
k562_dels = k562_full[k562_full[9].isin(dels2filter)]
# hepg2 #
hepg2_dels = hepg2_full[hepg2_full[9].isin(dels2filter)]
# sknsh #
sknsh_dels = sknsh_full[sknsh_full[9].isin(dels2filter)]
# save dELS only to disk #
# k562 #
k562_dels.to_csv('processed_data/K526.dELS.only.bed', sep = '\t', header = None, index = False)
# hepg2 #
hepg2_dels.to_csv('processed_data/HepG2.dELS.only.bed', sep = '\t', header = None, index = False)
# sknsh #
sknsh_dels.to_csv('processed_data/SK-N-SH.dELS.only.bed', sep = '\t', header = None, index = False)
