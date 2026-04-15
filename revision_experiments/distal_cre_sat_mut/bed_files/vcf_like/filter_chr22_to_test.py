# script to filter all sat mut of enahncers for only those on chromosome 22 to test #

# import packages #
import pandas as pd

# open full vcf #
full_vcf = pd.read_csv('all_GRCh38-dELS_vcf_like.tsv', sep = '\t', header = None)
# filter for only chromosome 22 and drop any 'N's
chr22_vcf = full_vcf[(full_vcf[0] == 'chr22') & (full_vcf[4] != 'N')]
# chr22 vcf to disk #
chr22_vcf.to_csv('chr22_GRCh38-dELS_vcf_like.vcf', sep = '\t', index = False, header = False)
