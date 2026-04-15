#!/usr/bin/env python3
"""
Stage 2: Per-chrom gnomAD annotation and filtering for openTargets seqlet analysis.

For each chromosome, loads the candidate variants from Stage 1, queries gnomAD v4
allele frequencies, then applies the following filters:

  Lead variants (is_leadVariant == 1):
    - Kept as-is (already confirmed in a seqlet with PIP >= 0.9)

  Non-lead emVars (is_leadVariant == 0):
    - Kept only if BOTH:
        (a) exceeds_lead_skew == True  (stronger predicted effect than lead)
        (b) gnomAD v4 allele frequency is present (variant is in gnomAD)

  Controls:
    - Loaded from results_final/control_per_chrom/{chrom}_control_variants.tsv
    - Filtered to (enhancer_ids, cell_type) pairs present in the filtered emVars
    - Already satisfy: outside seqlets, not emVars (|skew| <= 0.5), have gnomAD AF

Adds a 'variant_class' column: 'lead_emVar', 'shadow_emVar', or 'control'.

Usage:
    python scripts/openTargets_gnomad_filter.py --chrom chr22
    (run from multiTF_cluster/ directory)
"""

import argparse
import os
import subprocess
import tempfile
from io import StringIO

import numpy as np
import pandas as pd
from tqdm import tqdm

GNOMAD_TEMPLATE = (
    '/projects/tewhey-lab/buttsj/Variant_Effects/gnomad/'
    'gnomad_v4/filtered_data/popLevel/{chrom}_gnomAD_v4_pass_popLevel_af.tsv.gz'
)
AF_POP_COLS = [
    'AF', 'AF_afr', 'AF_ami', 'AF_amr', 'AF_asj', 'AF_eas',
    'AF_fin', 'AF_mid', 'AF_nfe', 'AF_remaining', 'AF_sas'
]


def query_gnomad_af(df, chrom):
    """Query gnomAD v4 pop-level AF for all variants in df on this chrom."""
    filepath = GNOMAD_TEMPLATE.format(chrom=chrom)
    if not os.path.exists(filepath):
        print(f"  Warning: gnomAD file not found: {filepath}")
        df['af'] = np.nan
        for col in AF_POP_COLS[1:]:
            df[col] = np.nan
        return df

    positions = df['pos'].unique()

    header_cmd = f"zcat {filepath} | head -1"
    header = subprocess.run(
        header_cmd, shell=True, capture_output=True, text=True
    ).stdout.strip().split('\t')

    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as pf:
        for pos in positions:
            pf.write(f'{chrom}\t{pos}\t\n')
        pattern_file = pf.name

    try:
        grep_cmd = f"zcat {filepath} | grep -F -f {pattern_file}"
        result = subprocess.run(grep_cmd, shell=True, capture_output=True, text=True)
    finally:
        os.unlink(pattern_file)

    if not result.stdout:
        df['af'] = np.nan
        for col in AF_POP_COLS[1:]:
            df[col] = np.nan
        return df

    gnomad_df = pd.read_csv(StringIO(result.stdout), sep='\t', names=header, na_values='.')
    gnomad_df['openTargets_id'] = (
        gnomad_df['CHROM'].str.replace('chr', '') + '_' +
        gnomad_df['POS'].astype(str) + '_' +
        gnomad_df['REF'] + '_' + gnomad_df['ALT']
    )
    gnomad_df = gnomad_df[['openTargets_id'] + AF_POP_COLS].rename(columns={'AF': 'af'})

    df = df.merge(gnomad_df, on='openTargets_id', how='left')
    return df


def main():
    parser = argparse.ArgumentParser(
        description='Stage 2: gnomAD annotation and filtering for openTargets seqlet analysis'
    )
    parser.add_argument('--chrom', required=True, help='Chromosome (e.g., chr22)')
    parser.add_argument(
        '--output', default=None,
        help='Output TSV path. Default: results_final/ot_filtered_per_chrom/{chrom}_ot_filtered.tsv'
    )
    args = parser.parse_args()

    if args.output is None:
        args.output = f'results_final/ot_filtered_per_chrom/{args.chrom}_ot_filtered.tsv'
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    # --- Load candidates ---
    candidates_path = f'results_final/ot_candidates_per_chrom/{args.chrom}_ot_candidates.tsv'
    if not os.path.exists(candidates_path):
        print(f"No candidates file found for {args.chrom}: {candidates_path}")
        return

    print(f"Loading candidates for {args.chrom}...")
    candidates = pd.read_csv(candidates_path, sep='\t')
    print(f"  {len(candidates):,} candidate variants")
    print(f"  Lead variants: {candidates[candidates['is_leadVariant'] == 1]['leadVariant'].nunique()}")

    # --- Query gnomAD v4 ---
    print(f"Querying gnomAD v4 for {args.chrom}...")
    candidates = query_gnomad_af(candidates, args.chrom)
    af_found = candidates['af'].notna().sum()
    print(f"  gnomAD AF found for {af_found:,}/{len(candidates):,} variants")

    # --- Apply filters ---
    lead_mask = candidates['is_leadVariant'] == 1
    shadow_mask = (
        (candidates['is_leadVariant'] == 0) &
        candidates['exceeds_lead_skew'] &
        candidates['af'].notna()
    )

    lead_emVars = candidates[lead_mask].copy()
    shadow_emVars = candidates[shadow_mask].copy()

    lead_emVars['variant_class'] = 'lead_emVar'
    shadow_emVars['variant_class'] = 'shadow_emVar'

    print(f"  Lead emVars retained: {len(lead_emVars):,}")
    print(f"  Shadow emVars retained (exceeds lead + has AF): {len(shadow_emVars):,}")

    filtered_emVars = pd.concat([lead_emVars, shadow_emVars], ignore_index=True)

    if len(filtered_emVars) == 0:
        print("  No filtered emVars — skipping controls and output.")
        return

    # --- Add controls ---
    control_path = f'results_final/control_per_chrom/{args.chrom}_control_variants.tsv'
    if os.path.exists(control_path):
        controls = pd.read_csv(control_path, sep='\t')
        # Filter controls to (enhancer, cell_type) pairs present in filtered emVars
        valid_pairs = set(zip(filtered_emVars['enhancer_ids'], filtered_emVars['cell_type']))
        ctrl_mask = [
            (eid, ct) in valid_pairs
            for eid, ct in zip(controls['enhancer_ids'], controls['cell_type'])
        ]
        controls_filtered = controls[ctrl_mask].copy()
        print(f"  Controls: {len(controls_filtered):,} of {len(controls):,} "
              f"match filtered (enhancer, cell_type) pairs")

        # Align columns: controls may have columns emVars don't and vice versa
        combined = pd.concat([filtered_emVars, controls_filtered], ignore_index=True, sort=False)
    else:
        print(f"  Warning: no control file found at {control_path}")
        combined = filtered_emVars

    print(f"\nCombined output: {len(combined):,} rows")
    for vc, n in combined['variant_class'].value_counts().items():
        print(f"  {vc}: {n:,}")

    combined.to_csv(args.output, sep='\t', index=False)
    print(f"\nSaved to {args.output}")


if __name__ == '__main__':
    main()
