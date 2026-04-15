# script to filter all sat mut of enahncers for only those on chromosome 22 to test #

# import packages #
import pandas as pd
import numpy as np

# open full vcf #
full_vcf = pd.read_csv('K562_enhancers_vcf_like.tsv', sep = '\t', header = None)
# filter for only chromosome 22 and drop any 'N's
chr22_vcf = full_vcf[(full_vcf[0] == 'chr22') & (full_vcf[4] != 'N')]
print(len(chr22_vcf))
# iterate through chrom 22 vcf and write vcfs to disk in 500000 row chunks #
chunk_size = 10000
num_chunks = len(chr22_vcf) // chunk_size + 1

for i, chunk in enumerate(np.array_split(chr22_vcf, num_chunks)):
    chunk.to_csv(f"chunked_vcfs/k562_chr22_cre_chunk_{i}_10k.vcf", sep = '\t', index = False, header = False)
