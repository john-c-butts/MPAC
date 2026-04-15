### script to explode out the credible set variants from Open Targets ###

# import packages
import pandas as pd
import os
from tqdm import tqdm
import numpy as np
import gc
import psutil

# open the opentargets data
# store the path to the open targets data as a variable (credible_sets)
path2openTargets = '/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/raw_data/openTargets/credible_set'
# open the study info
openTargetsStudyInfo = pd.read_parquet('/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/raw_data/openTargets/study/part-00000-2c4d825c-7fce-4f14-ae0e-0c516c8624e1-c000.snappy.parquet')

# iterate through all open targets data and collapse into a single DF
openTargets2cat = []
for i in tqdm([i for i in os.listdir(path2openTargets) if i.endswith('.parquet')]):
    parqLife = pd.read_parquet(f'{path2openTargets}/{i}')                                                                                                                                                                  
    openTargets2cat.append(parqLife)
openTargetsCredibleSets = pd.concat(openTargets2cat)

# clean up
del openTargets2cat
gc.collect()

# merge study IDs with credible set data
openTargetsMerged = openTargetsCredibleSets.merge(
    openTargetsStudyInfo[['studyId', 'traitFromSource', 'traitFromSourceMappedIds', 'pubmedId', 'initialSampleSize',
                          'cohorts','nCases', 'nControls', 'nSamples', 'diseaseIds']], 
    on='studyId', 
    how='left'
)

# clean up
del openTargetsCredibleSets, openTargetsStudyInfo
gc.collect()

# filter for only gwas
openTargetsGWAS = openTargetsMerged[openTargetsMerged['studyType'] == 'gwas'].copy()
# rename the variantID to leadVariant - will have a duplicate columname
openTargetsGWAS = openTargetsGWAS.rename(columns={'variantId' : 'leadVariant'})

# clean up
del openTargetsMerged
gc.collect()

# chunk and explode out the GWAS credible sets
chunk_size = 25000
output_dir = '/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/multiTF_cluster/exploded_chunks'
os.makedirs(output_dir, exist_ok=True)
dupe_cols = ['pValueMantissa', 'pValueExponent', 'beta', 'standardError']

for i, start in enumerate(range(0, len(openTargetsGWAS), chunk_size)):
    mem = psutil.Process().memory_info().rss / 1e9
    print(f"Chunk {i}, rows {start}-{start + chunk_size}, Memory: {mem:.1f} GB")
    
    output_file = f'{output_dir}/openTargetsGWAS_explode_chunk_{i:04d}.parquet'
    
    if os.path.exists(output_file):
        print(f"Skipping chunk {i}, already exists")
        continue
    
    chunk = openTargetsGWAS.iloc[start:start + chunk_size].copy()
    chunk = chunk.rename(columns={col: f'{col}_leadVariant' for col in dupe_cols})
    chunk['locus'] = chunk['locus'].apply(lambda x: list(x) if isinstance(x, np.ndarray) else x)
    chunk_exploded = chunk.explode('locus').reset_index(drop=True)
    
    mask = chunk_exploded['locus'].notna()
    chunk_valid = chunk_exploded[mask].copy()
    
    chunk_final = pd.concat([
        chunk_valid.drop(columns=['locus']).reset_index(drop=True),
        pd.json_normalize(chunk_valid['locus'].tolist())
    ], axis=1)
    
    chunk_final.to_parquet(output_file, index=False)
    
    del chunk, chunk_exploded, chunk_valid, chunk_final
    gc.collect()

print("Done.")
