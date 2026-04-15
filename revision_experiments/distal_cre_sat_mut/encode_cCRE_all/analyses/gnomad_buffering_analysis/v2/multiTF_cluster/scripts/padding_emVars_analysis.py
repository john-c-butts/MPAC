#!/usr/bin/env python3
"""
Generate 4bp padding-zone emVars for multiTF enhancers with gnomAD v4 allele frequencies.

For each filtered seqlet in each multiTF enhancer, identifies emVars that fall within
the padding zones flanking the seqlet (but not overlapping it). These variants are
not captured by the exact BED intersection in multiTF_analysis.py but are biologically
relevant as they sit near TF binding sites.

Uses targeted grep by enhancer ID from MPAC files (not full-file loading) for speed.

Can run per-chromosome (for SLURM array parallelism) or all-at-once.

Usage:
    python scripts/padding_emVars_analysis.py --chrom chr9 --output results_final/padding_per_chrom/chr9_padding_emVars.tsv
    python scripts/padding_emVars_analysis.py  # all chromosomes
"""

import argparse
import os
import pickle
import subprocess
import tempfile
from io import StringIO

import pandas as pd
from tqdm import tqdm


def get_padding_regions(chrom_data, cell_types, pad):
    """
    Extract padding regions from filtered seqlets in the pickle data.
    """
    regions = []
    for cell_type in cell_types:
        if cell_type not in chrom_data['multiTF_enhancers']:
            continue
        for enhancer_id, enh_data in chrom_data['multiTF_enhancers'][cell_type].items():
            filtered_seqlets = enh_data['filtered_seqlets']
            for _, seqlet in filtered_seqlets.iterrows():
                regions.append({
                    'cell_type': cell_type,
                    'enhancer_id': enhancer_id,
                    'seqlet_start': seqlet['start'],
                    'seqlet_end': seqlet['end'],
                    'tf_family': seqlet.get('vierstra_cluster', ''),
                    'tf_contrib': seqlet.get('rep_tf_contrib', 0),
                    'pad_left_start': seqlet['start'] - pad,
                    'pad_left_end': seqlet['start'] - 1,
                    'pad_right_start': seqlet['end'] + 1,
                    'pad_right_end': seqlet['end'] + pad,
                })
    return regions


def extract_mpac_by_enhancer(mpac_path, enhancer_ids):
    """
    Use zcat | grep -F to extract MPAC rows for specific enhancer IDs.
    Much faster than grepping by position since there are far fewer enhancer IDs.
    Returns a DataFrame with all rows for those enhancers.
    """
    if not enhancer_ids:
        return pd.DataFrame()

    # Get header
    header_cmd = f"zcat {mpac_path} | head -1"
    header = subprocess.run(
        header_cmd, shell=True, capture_output=True, text=True
    ).stdout.strip().split('\t')

    # Write enhancer ID patterns to temp file
    # Use tab-delimited patterns for specificity: "\tENHANCER_ID\t"
    with tempfile.NamedTemporaryFile(
        mode='w', suffix='.txt', delete=False
    ) as pf:
        for enh_id in enhancer_ids:
            pf.write(f'\t{enh_id}\t\n')
        pattern_file = pf.name

    try:
        grep_cmd = f"zcat {mpac_path} | grep -F -f {pattern_file}"
        result = subprocess.run(
            grep_cmd, shell=True, capture_output=True, text=True
        )
    finally:
        os.unlink(pattern_file)

    if not result.stdout:
        return pd.DataFrame()

    df = pd.read_csv(
        StringIO(result.stdout), sep='\t', names=header, na_values='.'
    )
    return df


def find_padding_emvars(mpac_df, regions, emvar_threshold, existing_var_keys):
    """
    Find emVars in padding zones that are NOT already in multiTF_variants.tsv.
    """
    if len(mpac_df) == 0:
        return []

    padding_emvars = []
    for region in regions:
        cell = region['cell_type']
        enh_id = region['enhancer_id']
        skew_col = f'{cell}_skew_pred'

        # Filter MPAC for this enhancer
        enh_preds = mpac_df[mpac_df['id'] == enh_id]
        if len(enh_preds) == 0:
            continue

        # Left padding zone: [start - pad, start)
        left_mask = (
            (enh_preds['pos'] >= region['pad_left_start']) &
            (enh_preds['pos'] <= region['pad_left_end']) &
            (enh_preds[skew_col].abs() > emvar_threshold)
        )
        # Right padding zone: (end, end + pad]
        right_mask = (
            (enh_preds['pos'] >= region['pad_right_start']) &
            (enh_preds['pos'] <= region['pad_right_end']) &
            (enh_preds[skew_col].abs() > emvar_threshold)
        )

        pad_vars = enh_preds[left_mask | right_mask]
        if len(pad_vars) == 0:
            continue

        for _, var in pad_vars.iterrows():
            var_id = f"{var['chrom']}:{var['pos']}:{var['ref']}:{var['alt']}"
            if (var_id, cell) in existing_var_keys:
                continue
            padding_emvars.append({
                'variant_id': var_id,
                'chrom': var['chrom'],
                'pos': var['pos'],
                'ref': var['ref'],
                'alt': var['alt'],
                'cell_type': cell,
                'skew_pred': var[skew_col],
                'enhancer_ids': enh_id,
                'is_multiTF': True,
                'tf_family': region['tf_family'],
                'tf_contrib': region['tf_contrib'],
                'seqlet_start': region['seqlet_start'],
                'seqlet_end': region['seqlet_end'],
            })

    return padding_emvars


def query_gnomad_af(padding_df, gnomad_template):
    """
    Query gnomAD v4 pop-level allele frequencies for padding emVars.
    """
    af_pop_cols = [
        'AF', 'AF_afr', 'AF_ami', 'AF_amr', 'AF_asj', 'AF_eas',
        'AF_fin', 'AF_mid', 'AF_nfe', 'AF_remaining', 'AF_sas'
    ]
    af_matches = []

    for chrom in tqdm(padding_df['chrom'].unique(), desc='  Querying gnomAD'):
        chrom_vars = padding_df[padding_df['chrom'] == chrom]
        positions = chrom_vars['pos'].unique()

        filepath = gnomad_template.format(chrom=chrom)
        if not os.path.exists(filepath):
            print(f"  Warning: gnomAD file not found for {chrom}: {filepath}")
            continue

        # Get header
        header_cmd = f"zcat {filepath} | head -1"
        header = subprocess.run(
            header_cmd, shell=True, capture_output=True, text=True
        ).stdout.strip().split('\t')

        # Write position patterns to temp file
        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.txt', delete=False
        ) as pf:
            for pos in positions:
                pf.write(f'{chrom}\t{pos}\t\n')
            pattern_file = pf.name

        try:
            grep_cmd = f"zcat {filepath} | grep -F -f {pattern_file}"
            result = subprocess.run(
                grep_cmd, shell=True, capture_output=True, text=True
            )
        finally:
            os.unlink(pattern_file)

        if result.stdout:
            chunk_df = pd.read_csv(
                StringIO(result.stdout), sep='\t', names=header, na_values='.'
            )
            chunk_df['variant_id'] = (
                chunk_df['CHROM'] + ':' +
                chunk_df['POS'].astype(str) + ':' +
                chunk_df['REF'] + ':' +
                chunk_df['ALT']
            )
            af_matches.append(chunk_df[['variant_id'] + af_pop_cols])

    if af_matches:
        all_af = pd.concat(af_matches, ignore_index=True)
        padding_df = padding_df.merge(
            all_af, on='variant_id', how='left'
        ).rename(columns={'AF': 'af'})
    else:
        padding_df['af'] = None
        for col in af_pop_cols[1:]:
            padding_df[col] = None

    return padding_df


def process_chromosome(chrom, pickle_dir, mpac_template, cell_types, pad,
                       emvar_threshold, existing_var_keys):
    """Process a single chromosome and return padding emVars list."""
    pkl_path = os.path.join(pickle_dir, f'{chrom}_multiTF.pkl')
    if not os.path.exists(pkl_path):
        print(f"  Skipping {chrom}: pickle not found")
        return []

    # Load pickle
    with open(pkl_path, 'rb') as f:
        chrom_data = pickle.load(f)

    # Get padding regions from filtered seqlets
    regions = get_padding_regions(chrom_data, cell_types, pad)
    if not regions:
        print(f"  {chrom}: no padding regions")
        return []

    # Collect unique enhancer IDs that have padding zones
    enhancer_ids = set(r['enhancer_id'] for r in regions)
    print(f"  {chrom}: {len(enhancer_ids)} enhancers with "
          f"{len(regions)} seqlet padding zones")

    # Extract MPAC rows for those enhancers via grep
    mpac_path = mpac_template.format(chrom=chrom)
    if not os.path.exists(mpac_path):
        print(f"  Skipping {chrom}: MPAC predictions not found")
        return []

    mpac_df = extract_mpac_by_enhancer(mpac_path, enhancer_ids)
    if len(mpac_df) == 0:
        print(f"  {chrom}: no MPAC rows matched")
        return []

    print(f"  {chrom}: extracted {len(mpac_df)} MPAC rows for {len(enhancer_ids)} enhancers")

    # Find padding-zone emVars
    chrom_padding = find_padding_emvars(
        mpac_df, regions, emvar_threshold, existing_var_keys
    )
    print(f"  {chrom}: {len(chrom_padding)} padding emVars found")

    del mpac_df, chrom_data
    return chrom_padding


def main():
    parser = argparse.ArgumentParser(
        description='Generate padding-zone emVars with gnomAD v4 allele frequencies'
    )
    parser.add_argument(
        '--chrom', default=None,
        help='Process a single chromosome (e.g., chr9). Default: all chromosomes.'
    )
    parser.add_argument(
        '--pad', type=int, default=4,
        help='Padding size in bp (default: 4)'
    )
    parser.add_argument(
        '--output', default=None,
        help='Output TSV path. Default: results_final/padding_emVars_gnomAD_v4.tsv '
             '(or results_final/padding_per_chrom/{chrom}_padding_emVars.tsv for single chrom)'
    )
    args = parser.parse_args()

    # Set default output path
    if args.output is None:
        if args.chrom:
            args.output = f'results_final/padding_per_chrom/{args.chrom}_padding_emVars.tsv'
        else:
            args.output = 'results_final/padding_emVars_gnomAD_v4.tsv'

    # --- Paths ---
    base_path = (
        '/projects/tewhey-lab/buttsj/Variant_Effects/'
        'revision_experiments/distal_cre_sat_mut/encode_cCRE_all'
    )
    pickle_dir = 'results_per_chrom'
    mpac_template = (
        f'{base_path}/mpac_preds/GRCh38-dELS-{{chrom}}-ALL-mpac-017.tsv.gz'
    )
    gnomad_template = (
        '/projects/tewhey-lab/buttsj/Variant_Effects/gnomad/'
        'gnomad_v4/filtered_data/popLevel/{chrom}_gnomAD_v4_pass_popLevel_af.tsv.gz'
    )
    existing_variants_path = 'results_final/multiTF_variants.tsv'

    cell_types = ['k562', 'hepg2', 'sknsh']
    emvar_threshold = 0.5

    # --- Load existing variants to exclude (only current chrom if specified) ---
    existing_var_keys = set()
    if os.path.exists(existing_variants_path):
        print("Loading existing multiTF_variants.tsv...")
        if args.chrom:
            # Only load the current chromosome's variants for fast filtering
            chunks = pd.read_csv(existing_variants_path, sep='\t', chunksize=500000)
            for chunk in chunks:
                chrom_chunk = chunk[chunk['chrom'] == args.chrom]
                existing_var_keys.update(
                    zip(chrom_chunk['variant_id'], chrom_chunk['cell_type'])
                )
            print(f"  {len(existing_var_keys)} existing variant-cell pairs for {args.chrom}")
        else:
            existing_vars = pd.read_csv(existing_variants_path, sep='\t')
            existing_var_keys = set(
                zip(existing_vars['variant_id'], existing_vars['cell_type'])
            )
            print(f"  {len(existing_var_keys)} existing variant-cell pairs to exclude")

    # --- Determine chromosomes to process ---
    if args.chrom:
        chrom_list = [args.chrom]
    else:
        chrom_list = sorted(
            [f.replace('_multiTF.pkl', '')
             for f in os.listdir(pickle_dir)
             if f.endswith('_multiTF.pkl')],
            key=lambda x: int(x.replace('chr', ''))
        )

    # --- Process each chromosome ---
    all_padding_emvars = []

    for chrom in tqdm(chrom_list, desc='Processing chromosomes'):
        chrom_padding = process_chromosome(
            chrom, pickle_dir, mpac_template, cell_types, args.pad,
            emvar_threshold, existing_var_keys
        )
        all_padding_emvars.extend(chrom_padding)

    if not all_padding_emvars:
        print("\nNo padding emVars found.")
        return

    # --- Deduplicate ---
    padding_df = pd.DataFrame(all_padding_emvars)
    pre_dedup = len(padding_df)
    padding_df = padding_df.drop_duplicates(
        subset=['variant_id', 'cell_type', 'enhancer_ids']
    )
    print(f"\nDeduplication: {pre_dedup} → {len(padding_df)} rows")

    # --- Query gnomAD v4 allele frequencies ---
    print("\nQuerying gnomAD v4 allele frequencies...")
    padding_df = query_gnomad_af(padding_df, gnomad_template)

    af_found = padding_df['af'].notna().sum()
    print(f"  AF found for {af_found}/{len(padding_df)} variants")

    # --- Save ---
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    padding_df.to_csv(args.output, sep='\t', index=False)

    print(f"\nSaved {len(padding_df)} padding emVars to {args.output}")
    print(f"  Unique variants: {padding_df['variant_id'].nunique()}")
    print(f"  Per cell type:")
    for ct, count in padding_df.groupby('cell_type').size().items():
        print(f"    {ct}: {count}")


if __name__ == '__main__':
    main()
