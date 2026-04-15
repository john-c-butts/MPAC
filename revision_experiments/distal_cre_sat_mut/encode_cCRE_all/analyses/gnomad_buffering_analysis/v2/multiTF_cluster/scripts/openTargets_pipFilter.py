# import packages
import pandas as pd
import pyarrow.parquet as pq
import os
from tqdm import tqdm
import numpy as np
import gc
import psutil

cols_needed = [
    'studyLocusId', 'studyId', 'leadVariant', 'chromosome', 'position', 'region',
    'beta_leadVariant', 'zScore', 'pValueMantissa_leadVariant', 'pValueExponent_leadVariant',
    'effectAlleleFrequencyFromSource', 'standardError_leadVariant',
    'finemappingMethod', 'credibleSetIndex', 'credibleSetlog10BF',
    'purityMeanR2', 'purityMinR2', 'locusStart', 'locusEnd', 'sampleSize',
    'confidence', 'studyType', 'traitFromSource', 'pubmedId', 'initialSampleSize',
    'nCases', 'nControls', 'nSamples',
    'is95CredibleSet', 'is99CredibleSet', 'logBF', 'posteriorProbability',
    'variantId', 'pValueMantissa', 'pValueExponent', 'beta', 'standardError', 'r2Overall'
]

path2explodedChunks = '/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/multiTF_cluster/exploded_chunks'

# High PIP filter
output_dir = '/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/multiTF_cluster/exploded_chunks/highPIP_09'
os.makedirs(output_dir, exist_ok=True)

parquet_files = sorted([f for f in os.listdir(path2explodedChunks) if f.endswith('.parquet')])

for i, parq in tqdm(enumerate(parquet_files)):
    parqLife = pd.read_parquet(f'{path2explodedChunks}/{parq}', columns=cols_needed)
    
    highPipParq = parqLife[parqLife['posteriorProbability'] > 0.9].copy()
    output_file = f'{output_dir}/openTargetsGWAS_explode_highPIP_chunk_{i:04d}.parquet'
    highPipParq.to_parquet(output_file, index=False)
    
    del parqLife, highPipParq
    gc.collect()

print('done filtering for high PIP (>0.9)')

# Mid PIP filter
output_dir = '/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/multiTF_cluster/exploded_chunks/midPIP_05'
os.makedirs(output_dir, exist_ok=True)

for i, parq in tqdm(enumerate(parquet_files)):
    parqLife = pd.read_parquet(f'{path2explodedChunks}/{parq}', columns=cols_needed)
    
    midPipParq = parqLife[parqLife['posteriorProbability'] > 0.5].copy()
    output_file = f'{output_dir}/openTargetsGWAS_explode_midPIP_chunk_{i:04d}.parquet'
    midPipParq.to_parquet(output_file, index=False)
    
    del parqLife, midPipParq
    gc.collect()

print('done filtering for mid PIP (>0.5)')