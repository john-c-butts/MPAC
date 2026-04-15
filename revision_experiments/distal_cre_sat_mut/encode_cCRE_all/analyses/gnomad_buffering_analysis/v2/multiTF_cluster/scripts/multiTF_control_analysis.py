#!/usr/bin/env python3
"""
Generate control variants for multiTF enhancer analysis.

For each multiTF enhancer, extracts variants that are:
  - Outside the duplicated seqlets (+ padding zone)
  - Not emVars (|skew_pred| <= 0.5)
  - Present in gnomAD v4 (have allele frequency annotation)

Each control variant is annotated with:
  - abs_skew and abs_skew_bin
  - Whether it falls in another (non-duplicated) TF seqlet
  - Percentile rank of its |skew| among all variants in the enhancer
  - lead_abs_skew, max_abs_skew, skew_ratio, exceeds_lead_skew (from source files)

Output columns are compatible with GTEx_MultiTF_allEmVars and
openTargets_GWAS_multiTF_allEmVars files for easy concatenation.

Can run per-chromosome (for SLURM array parallelism).

Usage:
    python scripts/multiTF_control_analysis.py --chrom chr22 \
        --enhancer-source GTEx_MultiTF_allEmVars_HighPIP_gnomAD_v4_annotated_4bp_pad.tsv \
        --enhancer-source openTargets_GWAS_multiTF_allEmVars_highPIP_gnomAD_v4_annotated_4bp_pad.tsv \
        --output results_final/control_per_chrom/chr22_control_variants.tsv
"""

import argparse
import os
import pickle
import subprocess
import tempfile
from io import StringIO

import numpy as np
import pandas as pd
from tqdm import tqdm


# ---------------------------------------------------------------------------
# Reused from padding_emVars_analysis.py
# ---------------------------------------------------------------------------

def extract_mpac_by_enhancer(mpac_path, enhancer_ids):
    """
    Use zcat | grep -F to extract MPAC rows for specific enhancer IDs.
    Returns a DataFrame with all rows for those enhancers.
    """
    if not enhancer_ids:
        return pd.DataFrame()

    header_cmd = f"zcat {mpac_path} | head -1"
    header = subprocess.run(
        header_cmd, shell=True, capture_output=True, text=True
    ).stdout.strip().split('\t')

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


def query_gnomad_af(df, gnomad_template):
    """
    Query gnomAD v4 pop-level allele frequencies.
    Returns df with AF columns merged in.
    """
    af_pop_cols = [
        'AF', 'AF_afr', 'AF_ami', 'AF_amr', 'AF_asj', 'AF_eas',
        'AF_fin', 'AF_mid', 'AF_nfe', 'AF_remaining', 'AF_sas'
    ]
    af_matches = []

    for chrom in tqdm(df['chrom'].unique(), desc='  Querying gnomAD'):
        chrom_vars = df[df['chrom'] == chrom]
        positions = chrom_vars['pos'].unique()

        filepath = gnomad_template.format(chrom=chrom)
        if not os.path.exists(filepath):
            print(f"  Warning: gnomAD file not found for {chrom}: {filepath}")
            continue

        header_cmd = f"zcat {filepath} | head -1"
        header = subprocess.run(
            header_cmd, shell=True, capture_output=True, text=True
        ).stdout.strip().split('\t')

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
        df = df.merge(
            all_af, on='variant_id', how='left'
        ).rename(columns={'AF': 'af'})
    else:
        df['af'] = None
        for col in af_pop_cols[1:]:
            df[col] = None

    return df


# ---------------------------------------------------------------------------
# Control-analysis-specific functions
# ---------------------------------------------------------------------------

def get_intervals(enh_data, pad=4):
    """
    Get exclusion zones and other-seqlet intervals for an enhancer.

    Returns:
        excluded: list of (start, end) tuples — duplicated seqlets + pad
        other_seqlets: list of (start, end, cluster) tuples — non-duplicated seqlets
    """
    filtered = enh_data['filtered_seqlets']
    all_seq = enh_data['all_seqlets']

    # Excluded intervals: duplicated seqlets + padding
    excluded = [
        (row['start'] - pad, row['end'] + pad)
        for _, row in filtered.iterrows()
    ]

    # Identify non-duplicated seqlets by coordinate match
    filtered_coords = set(
        zip(filtered['start'], filtered['end'])
    )
    other_seqlets = [
        (row['start'], row['end'], row.get('vierstra_cluster', ''))
        for _, row in all_seq.iterrows()
        if (row['start'], row['end']) not in filtered_coords
    ]

    return excluded, other_seqlets


def pos_in_intervals(pos, intervals):
    """Check if a position falls within any (start, end) interval."""
    for s, e in intervals:
        if s <= pos <= e:
            return True
    return False


def annotate_other_seqlet(pos, other_seqlet_intervals):
    """
    Check if pos falls in any non-duplicated seqlet.

    Returns:
        (bool, str): in_other_seqlet flag, comma-separated cluster name(s)
    """
    clusters = []
    for s, e, cluster in other_seqlet_intervals:
        if s <= pos <= e:
            if cluster and cluster not in clusters:
                clusters.append(cluster)
    if clusters:
        return True, ','.join(clusters)
    return False, ''


def process_enhancer_controls(enh_preds, enh_data, cell_type,
                              enh_skew_info, emvar_threshold=0.5, pad=4):
    """
    Extract control variants for a single enhancer in a single cell type.

    Args:
        enh_preds: DataFrame of ALL MPAC predictions for this enhancer
        enh_data: dict with 'valid_clusters', 'filtered_seqlets', 'all_seqlets'
        cell_type: k562, hepg2, or sknsh
        enh_skew_info: dict with 'lead_abs_skew', 'max_abs_skew' for this
                       (enhancer, cell_type) from source files
        emvar_threshold: |skew| threshold for emVar definition
        pad: bp padding around duplicated seqlets to exclude

    Returns:
        list of dicts, one per control variant
    """
    skew_col = f'{cell_type}_skew_pred'

    if skew_col not in enh_preds.columns or len(enh_preds) == 0:
        return []

    # Drop rows with NaN skew
    preds = enh_preds.dropna(subset=[skew_col]).copy()
    if len(preds) == 0:
        return []

    # Compute abs_skew for ALL variants (needed for ranking)
    preds['_abs_skew'] = preds[skew_col].abs()
    enhancer_n_variants = len(preds)

    # Build rank series: percentile rank of abs_skew among all enhancer variants
    preds['_rank'] = preds['_abs_skew'].rank(pct=True)

    # Get intervals
    excluded, other_seqlet_intervals = get_intervals(enh_data, pad=pad)

    # Filter: outside excluded intervals
    outside_mask = ~preds['pos'].apply(lambda p: pos_in_intervals(p, excluded))
    preds = preds[outside_mask]

    # Filter: non-emVars
    preds = preds[preds['_abs_skew'] <= emvar_threshold]

    if len(preds) == 0:
        return []

    # Bin absolute skew
    bins = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
    labels = ['0-0.1', '0.1-0.2', '0.2-0.3', '0.3-0.4', '0.4-0.5']
    preds['_bin'] = pd.cut(
        preds['_abs_skew'], bins=bins, labels=labels, include_lowest=True
    )

    # Annotate other-seqlet membership
    other_info = preds['pos'].apply(
        lambda p: annotate_other_seqlet(p, other_seqlet_intervals)
    )
    preds['_in_other'] = other_info.apply(lambda x: x[0])
    preds['_other_cluster'] = other_info.apply(lambda x: x[1])

    # Enhancer-level info
    dup_cluster = ','.join(sorted(enh_data['valid_clusters']))
    lead_abs_skew = enh_skew_info.get('lead_abs_skew', np.nan)
    max_abs_skew = enh_skew_info.get('max_abs_skew', np.nan)

    # Build output rows
    controls = []
    for _, row in preds.iterrows():
        var_id = f"{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}"
        abs_skew = row['_abs_skew']

        controls.append({
            # --- Shared columns (compatible with GTEx/OT files) ---
            'variant_id': var_id,
            'chrom': row['chrom'],
            'pos': row['pos'],
            'ref': row['ref'],
            'alt': row['alt'],
            'cell_type': cell_type,
            'skew_pred': row[skew_col],
            'enhancer_ids': row['id'],
            'is_multiTF': True,
            'tf_family': dup_cluster,
            'tf_contrib': np.nan,
            'tf_instance': np.nan,
            'is_leadVariant': False,
            'abs_skew': abs_skew,
            'lead_abs_skew': lead_abs_skew,
            'max_abs_skew': max_abs_skew,
            'skew_ratio': (max_abs_skew / lead_abs_skew
                           if pd.notna(lead_abs_skew) and lead_abs_skew > 0
                           else np.nan),
            'exceeds_lead_skew': (abs_skew > lead_abs_skew
                                  if pd.notna(lead_abs_skew)
                                  else False),
            # --- Control-specific columns ---
            'variant_class': 'control',
            'abs_skew_bin': row['_bin'],
            'in_other_seqlet': row['_in_other'],
            'other_seqlet_cluster': row['_other_cluster'],
            'enhancer_skew_rank': row['_rank'],
            'enhancer_n_variants': enhancer_n_variants,
        })

    return controls


def load_enhancer_source(source_files):
    """
    Load allowed (enhancer, cell_type) pairs and per-pair skew info
    from one or more source TSV files.

    Only generates controls for (enhancer, cell_type) combinations that
    actually appear in the source analyses.

    Returns:
        allowed_pairs: set of (enhancer_id, cell_type) tuples
        skew_lookup: dict of (enhancer_id, cell_type) -> {lead_abs_skew, max_abs_skew}
    """
    all_dfs = []

    for filepath in source_files:
        df = pd.read_csv(
            filepath, sep='\t',
            usecols=['enhancer_ids', 'cell_type', 'lead_abs_skew', 'max_abs_skew']
        )
        file_enhs = set(df['enhancer_ids'].unique())
        file_pairs = set(zip(df['enhancer_ids'], df['cell_type']))
        print(f"  {os.path.basename(filepath)}: {len(file_enhs)} enhancers, "
              f"{len(file_pairs)} (enhancer, cell_type) pairs")
        all_dfs.append(df)

    # Build per-(enhancer, cell_type) skew lookup
    # Values are consistent within (enhancer, cell_type) per source file,
    # but if an enhancer appears in both GTEx and OT we take the max
    combined = pd.concat(all_dfs, ignore_index=True)
    skew_lookup = {}
    for (enh_id, ct), group in combined.groupby(['enhancer_ids', 'cell_type']):
        skew_lookup[(enh_id, ct)] = {
            'lead_abs_skew': group['lead_abs_skew'].max(),
            'max_abs_skew': group['max_abs_skew'].max(),
        }

    allowed_pairs = set(skew_lookup.keys())
    allowed_enhancers = set(enh for enh, _ in allowed_pairs)
    print(f"  Total: {len(allowed_enhancers)} enhancers, "
          f"{len(allowed_pairs)} (enhancer, cell_type) pairs")
    return allowed_pairs, skew_lookup


def process_chromosome(chrom, pickle_dir, mpac_template, cell_types,
                       emvar_threshold, pad, allowed_pairs, skew_lookup):
    """Process a single chromosome and return control variant list."""
    pkl_path = os.path.join(pickle_dir, f'{chrom}_multiTF.pkl')
    if not os.path.exists(pkl_path):
        print(f"  Skipping {chrom}: pickle not found")
        return []

    with open(pkl_path, 'rb') as f:
        chrom_data = pickle.load(f)

    # Collect unique enhancer IDs across cell types,
    # filtered by allowed (enhancer, cell_type) pairs
    enhancer_ids = set()
    for ct in cell_types:
        if ct in chrom_data['multiTF_enhancers']:
            for eid in chrom_data['multiTF_enhancers'][ct].keys():
                if (eid, ct) in allowed_pairs:
                    enhancer_ids.add(eid)

    if not enhancer_ids:
        print(f"  {chrom}: no matching multiTF enhancers")
        return []

    print(f"  {chrom}: {len(enhancer_ids)} enhancer IDs after filtering")

    # Extract MPAC rows for all enhancers at once
    mpac_path = mpac_template.format(chrom=chrom)
    if not os.path.exists(mpac_path):
        print(f"  Skipping {chrom}: MPAC file not found")
        return []

    mpac_df = extract_mpac_by_enhancer(mpac_path, enhancer_ids)
    if len(mpac_df) == 0:
        print(f"  {chrom}: no MPAC rows matched")
        return []

    print(f"  {chrom}: extracted {len(mpac_df):,} MPAC rows")

    # Process each cell type and enhancer
    all_controls = []
    for ct in cell_types:
        if ct not in chrom_data['multiTF_enhancers']:
            continue

        ct_enhancers = {
            eid: edata
            for eid, edata in chrom_data['multiTF_enhancers'][ct].items()
            if (eid, ct) in allowed_pairs
        }
        ct_count = 0

        for enh_id, enh_data in tqdm(
            ct_enhancers.items(),
            desc=f'    {ct}',
            leave=False
        ):
            enh_preds = mpac_df[mpac_df['id'] == enh_id]
            if len(enh_preds) == 0:
                continue

            # Get lead/max skew from source files
            enh_skew_info = skew_lookup.get(
                (enh_id, ct),
                {'lead_abs_skew': np.nan, 'max_abs_skew': np.nan}
            )

            controls = process_enhancer_controls(
                enh_preds, enh_data, ct, enh_skew_info,
                emvar_threshold=emvar_threshold, pad=pad
            )
            all_controls.extend(controls)
            ct_count += len(controls)

        print(f"    {ct}: {ct_count:,} control variants from "
              f"{len(ct_enhancers):,} enhancers")

    del mpac_df, chrom_data
    return all_controls


def main():
    parser = argparse.ArgumentParser(
        description='Generate control variants for multiTF enhancers'
    )
    parser.add_argument(
        '--chrom', required=True,
        help='Chromosome to process (e.g., chr22)'
    )
    parser.add_argument(
        '--pad', type=int, default=4,
        help='Exclusion padding around duplicated seqlets in bp (default: 4)'
    )
    parser.add_argument(
        '--enhancer-source', action='append', required=True,
        help='TSV file(s) with enhancer_ids column to restrict analysis to. '
             'Can be specified multiple times. Enhancer IDs are unioned across files.'
    )
    parser.add_argument(
        '--output', default=None,
        help='Output TSV path. Default: results_final/control_per_chrom/{chrom}_control_variants.tsv'
    )
    args = parser.parse_args()

    if args.output is None:
        args.output = f'results_final/control_per_chrom/{args.chrom}_control_variants.tsv'

    # --- Load allowed (enhancer, cell_type) pairs and skew info ---
    print("Loading allowed enhancers from source files...")
    allowed_pairs, skew_lookup = load_enhancer_source(args.enhancer_source)

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

    cell_types = ['k562', 'hepg2', 'sknsh']
    emvar_threshold = 0.5

    # --- Process chromosome ---
    print(f"\nProcessing {args.chrom} (pad={args.pad})...")
    all_controls = process_chromosome(
        args.chrom, pickle_dir, mpac_template, cell_types,
        emvar_threshold, args.pad, allowed_pairs, skew_lookup
    )

    if not all_controls:
        print("\nNo control variants found.")
        return

    # --- Build DataFrame and deduplicate ---
    control_df = pd.DataFrame(all_controls)
    pre_dedup = len(control_df)
    control_df = control_df.drop_duplicates(
        subset=['variant_id', 'cell_type', 'enhancer_ids']
    )
    print(f"\nDeduplication: {pre_dedup:,} → {len(control_df):,} rows")

    # --- Query gnomAD v4 allele frequencies ---
    print("\nQuerying gnomAD v4 allele frequencies...")
    control_df = query_gnomad_af(control_df, gnomad_template)

    # --- Filter: keep only variants with AF annotation ---
    pre_filter = len(control_df)
    control_df = control_df[control_df['af'].notna()].reset_index(drop=True)
    print(f"  AF filter: {pre_filter:,} → {len(control_df):,} rows "
          f"({pre_filter - len(control_df):,} dropped without AF)")

    if len(control_df) == 0:
        print("\nNo control variants with gnomAD AF found.")
        return

    # --- Reorder columns for compatibility ---
    shared_cols = [
        'variant_id', 'chrom', 'pos', 'ref', 'alt', 'cell_type',
        'skew_pred', 'enhancer_ids', 'is_multiTF',
        'tf_family', 'tf_contrib', 'tf_instance', 'is_leadVariant',
        'af', 'AF_afr', 'AF_ami', 'AF_amr', 'AF_asj', 'AF_eas',
        'AF_fin', 'AF_mid', 'AF_nfe', 'AF_remaining', 'AF_sas',
        'abs_skew', 'lead_abs_skew', 'max_abs_skew', 'skew_ratio',
        'exceeds_lead_skew',
    ]
    control_cols = [
        'variant_class', 'abs_skew_bin', 'in_other_seqlet',
        'other_seqlet_cluster', 'enhancer_skew_rank', 'enhancer_n_variants',
    ]
    control_df = control_df[shared_cols + control_cols]

    # --- Save ---
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    control_df.to_csv(args.output, sep='\t', index=False)

    print(f"\nSaved {len(control_df):,} control variants to {args.output}")
    print(f"  Unique variants: {control_df['variant_id'].nunique():,}")
    print(f"  Per cell type:")
    for ct, count in control_df.groupby('cell_type').size().items():
        print(f"    {ct}: {count:,}")
    print(f"  In other seqlet: {control_df['in_other_seqlet'].sum():,} "
          f"({100 * control_df['in_other_seqlet'].mean():.1f}%)")
    print(f"  Abs skew bin distribution:")
    for b, count in control_df['abs_skew_bin'].value_counts().sort_index().items():
        print(f"    {b}: {count:,}")


if __name__ == '__main__':
    main()
