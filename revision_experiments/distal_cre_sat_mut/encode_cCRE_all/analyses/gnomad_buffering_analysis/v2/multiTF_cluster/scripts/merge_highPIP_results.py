#!/usr/bin/env python3
"""
Merge per-chromosome high-PIP shadow variant TSVs into a single file.

Usage:
    python scripts/merge_highPIP_results.py
"""

import os
import pandas as pd
from tqdm import tqdm

input_dir = 'results_final/shadow_per_chrom'
output_path = 'results_final/highPIP_shadow_variants.tsv'

chrom_files = sorted(
    [f for f in os.listdir(input_dir) if f.endswith('_shadow_variants.tsv')],
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
merged = merged.drop_duplicates(subset=['variant_id', 'cell_type', 'enhancer_ids', 'lead_tissue'])
print(f"\nTotal: {pre_dedup:,} → {len(merged):,} after dedup")

# 1. Filter: keep only rows where the lead is a predicted emVar (|lead_abs_skew| > 0.5)
pre_filter = len(merged)
merged = merged[merged['lead_abs_skew'] > 0.5].copy()
print(f"\nLead emVar filter (|lead_abs_skew| > 0.5): {pre_filter:,} → {len(merged):,} rows")

# 2. Add GTEx-format variant ID: chrom_pos_ref_alt_b38
merged['variant_id_gtex'] = (
    merged['chrom'].str.replace('chr', '', regex=False) + '_' +
    merged['pos'].astype(str) + '_' +
    merged['ref'] + '_' +
    merged['alt'] + '_b38'
)

# 3. Add skew fold change relative to lead (signed ratio: variant skew / lead skew)
merged['skew_fc_vs_lead'] = merged['abs_skew'] / merged['lead_abs_skew']

# 4. For exceeds_lead_skew cases, require directional concordance with lead
#    Lead skew sign is inferred from lead_abs_skew and the sign of skew_pred
#    being on the same enhancer. We store the lead's signed skew in the per-chrom
#    files as lead_skew_pred if available, otherwise derive from afc sign as proxy.
#    Here we use the variant's own skew_pred sign vs lead's sign, reconstructed
#    from the lead_variant_id lookup within the same cell_type group.
lead_sign = (
    merged[merged['is_lead_variant']]
    .groupby(['enhancer_ids', 'cell_type', 'lead_tissue'])['skew_pred']
    .first()
    .reset_index()
    .rename(columns={'skew_pred': 'lead_skew_sign_raw'})
)
merged = merged.merge(lead_sign, on=['enhancer_ids', 'cell_type', 'lead_tissue'], how='left')
merged['lead_skew_sign'] = merged['lead_skew_sign_raw'].apply(
    lambda x: 1 if pd.notna(x) and x >= 0 else -1
)
merged.drop(columns='lead_skew_sign_raw', inplace=True)

concordant = (merged['skew_pred'] * merged['lead_skew_sign']) > 0
# Only apply concordance filter to cases; keep all controls
case_mask = merged['exceeds_lead_skew']
pre_concordance = case_mask.sum()
merged = merged[~case_mask | concordant].copy()
print(f"Concordance filter (cases only): {pre_concordance:,} → {merged['exceeds_lead_skew'].sum():,} cases retained")

merged.to_csv(output_path, sep='\t', index=False)
print(f"\nSaved to {output_path}")
print(f"  Unique variants:  {merged['variant_id'].nunique():,}")
print(f"  Unique enhancers: {merged['enhancer_ids'].nunique():,}")
n_exceeds = merged['exceeds_lead_skew'].sum()
print(f"  Exceeds lead (cases): {n_exceeds:,} "
      f"({100 * n_exceeds / len(merged):.1f}%)")
print(f"  Controls:             {(~merged['exceeds_lead_skew']).sum():,}")
print(f"  Per cell type:")
for ct, count in merged.groupby('cell_type').size().items():
    ct_df = merged[merged['cell_type'] == ct]
    print(f"    {ct}: {count:,} ({ct_df['exceeds_lead_skew'].sum():,} exceed lead)")
