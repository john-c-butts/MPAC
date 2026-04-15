# script to convert the encode cCRE BED files to VCFs for predictions #

# import packages #
import pandas as pd
from tqdm import tqdm
# open distal CRE bed file #
full_bed = pd.read_csv('../processed_data/EH38E3356274.bed', sep = '\t', header = None)
# drop sex chromosome enhancers #
autosome_bed = full_bed[~full_bed[0].isin(['chrX', 'chrY'])]
# make a dictionary for adding enhancer annotation to fasta tsv #
bed_dict = dict(zip([chrom + ':' + ('-').join([str(start), 
                                               str(stop)]) for chrom, start, stop in zip(autosome_bed[0],
                                                                                         autosome_bed[1],
                                                                                         autosome_bed[2])], # id to match tsv
                     autosome_bed[4]))
# open fasta file
dels_fasta_tsv = pd.read_csv('../processed_data/alt_gtex_fastas/EH38E3356274.fa', sep = '\t', header = None)
# add enhancer id #
dels_fasta_tsv['id'] = [bed_dict.get(i) for i in dels_fasta_tsv[0]]
# iterate through CREs and convert to VCF #
dfs2cat = []
for enhancer, sequence, id in tqdm(zip(dels_fasta_tsv[0],
                                       dels_fasta_tsv[1],
                                       dels_fasta_tsv['id'])):
    # make tmp lists for building vcf-style df #
    chrom_tmp = []
    pos_tmp = []
    ref_tmp = []
    alt_tmp = []
    info_tmp = []
    # get position information #
    # chrom #
    chromosome = enhancer.split(':')[0]
    # start #
    # bed is 0-based so we add 1 to the start, but not the end #
    start = int(enhancer.split('-')[0].split(':')[-1]) + 1
    # end #
    end = int(enhancer.split('-')[-1])
    # iterate through fasta sequence to build df #
    for base in sequence:
        # append chrom 3X #
        chrom_tmp.extend([chromosome, chromosome, chromosome])
        # append id 3X #
        info_tmp.extend([id, id, id])
        # append position 3X and add 1 #
        pos_tmp.extend([start, start, start])
        start +=1
        # append ref 3X #
        ref_tmp.extend([base, base, base])
        # append mutations #
        if base.upper() == 'A':
            alt_tmp.extend(['T','C','G'])
        elif base.upper() == 'T':
            alt_tmp.extend(['A','C','G'])
        elif base.upper() == 'C':
            alt_tmp.extend(['A','T','G'])
        elif base.upper() == 'G':
            alt_tmp.extend(['A','T','C'])
        else:
            alt_tmp.extend(['N','N','N'])
    # build vcf-style df #
    vcf_df = pd.DataFrame({0 : chrom_tmp,
                           1 : pos_tmp,
                           2 : info_tmp,
                           3 : ref_tmp,
                           4 : alt_tmp})
    dfs2cat.append(vcf_df)
# combine all vcf-dfs into single df #
all_vcf_dfs = pd.concat(dfs2cat)
# save concatenated df to disk #
all_vcf_dfs.to_csv('../processed_data/vcfs4satmut/EH38E3356274_chr2_86221285_A_G_b38.vcf', sep = '\t', index = False, header = False)
