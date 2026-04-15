#!/usr/bin/env python3
"""
Merge per-chromosome multi-TF analysis results into final outputs.

Usage:
    python merge_results.py
    python merge_results.py --input-dir results_per_chrom --output-dir results_final
"""

import argparse
import os
import pickle
import pandas as pd
from glob import glob


def merge_results(input_dir, output_dir):
    """
    Merge all per-chromosome pickle files into final outputs.

    Output structure:
    {
        'variants_all': pd.DataFrame,      # All variants for annotation
        'multiTF_summary': pd.DataFrame,   # Quick reference table
        'multiTF_enhancers': {
            'k562': {enhancer_id: {...}},
            'hepg2': {enhancer_id: {...}},
            'sknsh': {enhancer_id: {...}}
        }
    }
    """
    # Find all per-chromosome pickle files
    pkl_files = sorted(glob(os.path.join(input_dir, 'chr*_multiTF.pkl')))

    if not pkl_files:
        print(f"No pickle files found in {input_dir}")
        return

    print(f"Found {len(pkl_files)} chromosome files to merge:")
    for f in pkl_files:
        print(f"  - {os.path.basename(f)}")

    # Collect results
    all_variants = []
    all_summaries = []
    merged_enhancers = {'k562': {}, 'hepg2': {}, 'sknsh': {}}

    for pkl_file in pkl_files:
        print(f"\nLoading {os.path.basename(pkl_file)}...")
        with open(pkl_file, 'rb') as f:
            data = pickle.load(f)

        chrom = data['chrom']
        n_variants = len(data['variants_all'])
        n_summary = len(data['multiTF_summary'])
        print(f"  {chrom}: {n_variants:,} variants, {n_summary:,} multi-TF enhancers")

        # Collect variants and summaries
        if n_variants > 0:
            all_variants.append(data['variants_all'])
        if n_summary > 0:
            all_summaries.append(data['multiTF_summary'])

        # Merge enhancer data
        for cell_type in ['k562', 'hepg2', 'sknsh']:
            for enh_id, enh_data in data['multiTF_enhancers'][cell_type].items():
                merged_enhancers[cell_type][enh_id] = enh_data

    # Merge dataframes
    print("\nMerging dataframes...")
    variants_df = pd.concat(all_variants, ignore_index=True) if all_variants else pd.DataFrame()
    summary_df = pd.concat(all_summaries, ignore_index=True) if all_summaries else pd.DataFrame()

    # Sort summary by n_emvars
    if len(summary_df) > 0:
        summary_df = summary_df.sort_values('n_emvars', ascending=False).reset_index(drop=True)

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Save TSV files for easy annotation
    print("\nSaving TSV files...")

    variants_path = os.path.join(output_dir, 'multiTF_variants.tsv')
    variants_df.to_csv(variants_path, sep='\t', index=False)
    print(f"  Saved {variants_path} ({len(variants_df):,} rows)")

    summary_path = os.path.join(output_dir, 'multiTF_summary.tsv')
    summary_df.to_csv(summary_path, sep='\t', index=False)
    print(f"  Saved {summary_path} ({len(summary_df):,} rows)")

    # Save full pickle for exploration
    print("\nSaving pickle for exploration...")
    full_data = {
        'variants_all': variants_df,
        'multiTF_summary': summary_df,
        'multiTF_enhancers': merged_enhancers
    }

    pkl_path = os.path.join(output_dir, 'multiTF_analysis.pkl')
    with open(pkl_path, 'wb') as f:
        pickle.dump(full_data, f)
    print(f"  Saved {pkl_path}")

    # Print final summary
    print(f"\n{'='*60}")
    print("MERGE COMPLETE")
    print(f"{'='*60}")
    print(f"Total variants:                  {len(variants_df):,}")
    if len(variants_df) > 0:
        print(f"  - K562:                        {len(variants_df[variants_df['cell_type'] == 'k562']):,}")
        print(f"  - HepG2:                       {len(variants_df[variants_df['cell_type'] == 'hepg2']):,}")
        print(f"  - SKNSH:                       {len(variants_df[variants_df['cell_type'] == 'sknsh']):,}")
        print(f"Variants in multi-TF enhancers:  {variants_df['is_multiTF'].sum():,}")
    print(f"\nMulti-TF enhancers with emVars:")
    for ct in ['k562', 'hepg2', 'sknsh']:
        print(f"  - {ct}: {len(merged_enhancers[ct]):,}")
    print(f"Total in summary table:          {len(summary_df):,}")
    print(f"{'='*60}")


def main():
    parser = argparse.ArgumentParser(description='Merge per-chromosome multi-TF analysis results')
    parser.add_argument('--input-dir', default='results_per_chrom', help='Input directory with per-chrom pickles')
    parser.add_argument('--output-dir', default='results_final', help='Output directory for merged results')
    args = parser.parse_args()

    # Handle relative paths
    base_dir = '/pod/2/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/multiTF_cluster'

    input_dir = args.input_dir if os.path.isabs(args.input_dir) else os.path.join(base_dir, args.input_dir)
    output_dir = args.output_dir if os.path.isabs(args.output_dir) else os.path.join(base_dir, args.output_dir)

    merge_results(input_dir, output_dir)


if __name__ == '__main__':
    main()
