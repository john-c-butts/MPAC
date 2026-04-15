# script to reshape the promoterAI preds into a VCF format for easier indexing with tabix #

# import packages
import pandas as pd

# load predictions
promoterAI_raw = pd.read_csv('../raw_data/PrimateAI_and_PromoterAI_scores.hg38.tsv.gz', sep = '\t')

# build VCF
promoterAI_vcf = pd.DataFrame({'#CHROM' : promoterAI_raw['chrom'],
                               'POS' : promoterAI_raw['pos'],
                               'ID' : [(':').join([chrom,
                                                   str(pos),
                                                   ref,
                                                   alt]) for chrom, pos, ref, alt in zip(promoterAI_raw['chrom'],
                                                                                         promoterAI_raw['pos'],
                                                                                         promoterAI_raw['ref'],
                                                                                         promoterAI_raw['alt'])],
                               'REF' : promoterAI_raw['ref'],
                               'ALT' : promoterAI_raw['alt'],
                               'QUAL' : ['.' for i in range(len(promoterAI_raw))],
                               'FILTER' : ['.' for i in range(len(promoterAI_raw))],
                               'INFO' : [f'PAI3D_percentile;{percentile}:PAI3D_thresh;{threshold}:PromoterAI_score;{score}' for percentile, 
                                                                                                                                 threshold, 
                                                                                                                                 score in zip(promoterAI_raw['PAI3D_percentile'],
                                                                                                                                              promoterAI_raw['PAI3D_gene_threshold'],
                                                                                                                                              promoterAI_raw['PromoterAI_score'])]})
promoterAI_vcf.to_csv('../processed_data/PrimateAI_and_PromoterAI_scores.hg38.vcf.gz', sep = '\t', index = False, compression = 'gzip')