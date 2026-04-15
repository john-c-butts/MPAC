# script to filter full cCRE bed file for only dELSs #

# import packages #
import pandas as pd

# open full bed file #
full_bed = pd.read_csv('GRCh38-cCREs.bed', sep = '\t', header = None)

# get all annotations including 'dELS' in there for filtering all intervals #
all_dELS = pd.Series([i for i in full_bed[5] if 'dELS' in i]).unique()

# filter full BED for only dELS #
dELS_only = full_bed[full_bed[5].isin(all_dELS)]

# save dels only BED to disk #
dELS_only.to_csv('GRCh38-dELS_only.bed', sep = '\t', index = False, header = False)
