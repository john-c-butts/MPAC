#!/usr/bin/env python3
"""
Merge per-chromosome control variant TSVs into a single file.

Usage:
    python scripts/merge_control_results.py
"""

import os
import pandas as pd
from tqdm import tqdm

input_dir = 'results_final/control_per_chrom'
output_path = 'results_final/multiTF_control_variants.tsv'

chrom_files = sorted(
    [f for f in os.listdir(input_dir) if f.endswith('_control_variants.tsv')],
    key=lambda x: int(x.split('_')[0].replace('chr', ''))
)

print(f"Found {len(chrom_files)} per-chromosome files")

dfs = []
for f in tqdm(chrom_files, desc='Merging'):
    df = pd.read_csv(os.path.join(input_dir, f), sep='\t')
    print(f"  {f}: {len(df):,} rows")
    dfs.append(df)

merged = pd.concat(dfs, ignore_index=True)
pre_dedup = len(merged)
merged = merged.drop_duplicates(subset=['variant_id', 'cell_type', 'enhancer_ids'])
print(f"\nTotal: {pre_dedup:,} → {len(merged):,} after dedup")

merged.to_csv(output_path, sep='\t', index=False)
print(f"Saved to {output_path}")
print(f"  Unique variants: {merged['variant_id'].nunique():,}")
for ct, count in merged.groupby('cell_type').size().items():
    print(f"  {ct}: {count:,}")
print(f"  In other seqlet: {merged['in_other_seqlet'].sum():,} "
      f"({100 * merged['in_other_seqlet'].mean():.1f}%)")
