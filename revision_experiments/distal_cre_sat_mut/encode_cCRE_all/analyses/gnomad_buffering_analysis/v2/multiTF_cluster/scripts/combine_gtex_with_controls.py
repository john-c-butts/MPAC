#!/usr/bin/env python3
"""
Combine GTEx multiTF allEmVars with matched control variants.

Controls are filtered to only include (enhancer, cell_type) pairs
present in the GTEx file. A 'variant_class' column distinguishes
emVars from controls.

Usage:
    python scripts/combine_gtex_with_controls.py
"""

import pandas as pd

gtex_path = 'GTEx_MultiTF_allEmVars_HighPIP_gnomAD_v4_annotated_4bp_pad.tsv'
control_path = 'results_final/multiTF_control_variants.tsv'
output_path = 'results_final/GTEx_MultiTF_allEmVars_with_controls.tsv'

# --- Load ---
print("Loading GTEx allEmVars...")
gtex = pd.read_csv(gtex_path, sep='\t')
gtex['variant_class'] = 'emVar'
print(f"  {len(gtex):,} rows, {gtex['enhancer_ids'].nunique()} enhancers")

print("Loading control variants...")
ctrl = pd.read_csv(control_path, sep='\t')
print(f"  {len(ctrl):,} total control rows")

# --- Filter controls to GTEx (enhancer, cell_type) pairs ---
gtex_pairs = set(zip(gtex['enhancer_ids'], gtex['cell_type']))
mask = [
    (eid, ct) in gtex_pairs
    for eid, ct in zip(ctrl['enhancer_ids'], ctrl['cell_type'])
]
ctrl_filtered = ctrl[mask].copy()
print(f"  {len(ctrl_filtered):,} controls matching GTEx (enhancer, cell_type) pairs")

# --- Concatenate ---
combined = pd.concat([gtex, ctrl_filtered], ignore_index=True)
print(f"\nCombined: {len(combined):,} rows")
print(f"  variant_class counts:")
for vc, count in combined['variant_class'].value_counts().items():
    print(f"    {vc}: {count:,}")

# --- Save ---
combined.to_csv(output_path, sep='\t', index=False)
print(f"\nSaved to {output_path}")
