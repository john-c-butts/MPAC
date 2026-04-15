#!/usr/bin/env python3
"""
High-PIP shadow variant analysis.

For each dELS enhancer containing a GTEx PIP > 0.9 eQTL variant, identifies
all other variants in that enhancer that have gnomAD v4 AF AND are predicted
to have a larger functional effect (|skew_pred| > lead's |skew_pred|).

These "shadow" variants have larger predicted effects than the confidently
causal eQTL but are nonetheless observed in the population, suggesting the
enhancer is buffered against strong perturbations.

Controls: variants in the same enhancers with gnomAD AF but
|skew_pred| <= lead's. Cases and controls are in the same output file,
distinguished by the 'exceeds_lead_skew' column.

Per (enhancer, cell_type), the lead is defined as the high-PIP variant with
the maximum |skew_pred| in that cell type (conservative: we only flag
'exceeds' if the variant exceeds even the most functionally active lead).

Can run per-chromosome (for SLURM array parallelism).

Usage:
    python scripts/highPIP_shadow_analysis.py --chrom chr22
    python scripts/highPIP_shadow_analysis.py --chrom chr22 \\
        --output results_final/shadow_per_chrom/chr22_shadow_variants.tsv
"""

import argparse
import os
import subprocess
import tempfile
from io import StringIO

import numpy as np
import pandas as pd
from tqdm import tqdm


# ---------------------------------------------------------------------------
# GTEx finemapping
# ---------------------------------------------------------------------------

def load_gtex_high_pip(finemap_dir, chrom, pip_threshold=0.9):
    """
    Load all-tissue GTEx v10 SuSiE finemapping, filter for PIP >= threshold
    on the given chromosome.

    Returns DataFrame with columns:
        chrom, pos, ref, alt, variant_id, pip, af_gtex, cs_id, cs_size,
        afc, afc_se, tissue, phenotype_id, gene_name, biotype
    """
    parquet_files = [
        f for f in os.listdir(finemap_dir)
        if f.endswith('.parquet')
    ]

    dfs = []
    for fname in sorted(parquet_files):
        tissue = fname.replace('.v10.eQTLs.SuSiE_summary.parquet', '')
        df = pd.read_parquet(os.path.join(finemap_dir, fname))

        df = df[df['pip'] >= pip_threshold].copy()
        if len(df) == 0:
            continue

        # Parse variant_id: chr1_64764_C_T_b38
        parts = df['variant_id'].str.split('_', expand=True)
        df['chrom'] = parts[0]
        df['pos'] = parts[1].astype(int)
        df['ref'] = parts[2]
        df['alt'] = parts[3]

        df = df[df['chrom'] == chrom]
        if len(df) == 0:
            continue

        df = df.rename(columns={'af': 'af_gtex'})
        df['tissue'] = tissue
        dfs.append(df)

    if not dfs:
        return pd.DataFrame()

    combined = pd.concat(dfs, ignore_index=True)
    return combined


# ---------------------------------------------------------------------------
# MPAC helpers (mirrors existing scripts)
# ---------------------------------------------------------------------------

def _run_grep(mpac_path, pattern_file):
    cmd = f"zcat {mpac_path} | grep -F -f {pattern_file}"
    return subprocess.run(cmd, shell=True, capture_output=True, text=True)


def _get_mpac_header(mpac_path):
    cmd = f"zcat {mpac_path} | head -1"
    return subprocess.run(
        cmd, shell=True, capture_output=True, text=True
    ).stdout.strip().split('\t')


def extract_mpac_by_position(mpac_path, positions):
    """
    Grep MPAC for rows matching (chrom, pos) pairs.
    positions: iterable of (chrom, pos) tuples.
    """
    if not positions:
        return pd.DataFrame()

    header = _get_mpac_header(mpac_path)

    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as pf:
        for chrom, pos in positions:
            pf.write(f'{chrom}\t{pos}\t\n')
        pattern_file = pf.name

    try:
        result = _run_grep(mpac_path, pattern_file)
    finally:
        os.unlink(pattern_file)

    if not result.stdout:
        return pd.DataFrame()

    return pd.read_csv(StringIO(result.stdout), sep='\t', names=header, na_values='.')


def extract_mpac_by_enhancer(mpac_path, enhancer_ids):
    """
    Grep MPAC for rows matching specific enhancer IDs (the 'id' column).
    """
    if not enhancer_ids:
        return pd.DataFrame()

    header = _get_mpac_header(mpac_path)

    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as pf:
        for enh_id in enhancer_ids:
            pf.write(f'\t{enh_id}\t\n')
        pattern_file = pf.name

    try:
        result = _run_grep(mpac_path, pattern_file)
    finally:
        os.unlink(pattern_file)

    if not result.stdout:
        return pd.DataFrame()

    return pd.read_csv(StringIO(result.stdout), sep='\t', names=header, na_values='.')


def query_gnomad_af(df, gnomad_template):
    """Query gnomAD v4 pop-level allele frequencies. Returns df with AF columns."""
    af_pop_cols = [
        'AF', 'AF_afr', 'AF_ami', 'AF_amr', 'AF_asj', 'AF_eas',
        'AF_fin', 'AF_mid', 'AF_nfe', 'AF_remaining', 'AF_sas'
    ]
    af_matches = []

    for chrom in df['chrom'].unique():
        chrom_vars = df[df['chrom'] == chrom]
        positions = chrom_vars['pos'].unique()

        filepath = gnomad_template.format(chrom=chrom)
        if not os.path.exists(filepath):
            print(f"  Warning: gnomAD file not found for {chrom}: {filepath}")
            continue

        header = _get_mpac_header(filepath)  # same grep pattern works

        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as pf:
            for pos in positions:
                pf.write(f'{chrom}\t{pos}\t\n')
            pattern_file = pf.name

        try:
            result = _run_grep(filepath, pattern_file)
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
        df = df.merge(all_af, on='variant_id', how='left').rename(columns={'AF': 'af'})
    else:
        df['af'] = None
        for col in af_pop_cols[1:]:
            df[col] = None

    return df


# ---------------------------------------------------------------------------
# Core analysis
# ---------------------------------------------------------------------------

def build_lead_lookup(gtex_high_pip, lead_mpac_df, cell_types):
    """
    Match high-PIP GTEx variants to their dELS enhancers via MPAC, then build
    a per-(enhancer, cell_type) lead reference.

    For each (enhancer, cell_type) we retain the high-PIP variant with the
    maximum |skew_pred| in that cell type. This is conservative: we only flag
    a variant as 'exceeds_lead' if it exceeds even the most functionally
    active lead in the enhancer.

    Returns:
        lead_df: DataFrame with one row per (enhancer_id, cell_type, tissue), columns:
            enhancer_id, cell_type, lead_variant_id, lead_abs_skew,
            lead_pip, lead_tissue, phenotype_id, gene_name, biotype,
            afc, afc_se, cs_id, cs_size
    """
    if len(lead_mpac_df) == 0:
        return pd.DataFrame()

    # Build variant_id key for joining
    lead_mpac_df = lead_mpac_df.copy()
    lead_mpac_df['variant_id'] = (
        lead_mpac_df['chrom'] + ':' +
        lead_mpac_df['pos'].astype(str) + ':' +
        lead_mpac_df['ref'] + ':' +
        lead_mpac_df['alt']
    )

    gtex_high_pip = gtex_high_pip.copy()
    gtex_high_pip['variant_id'] = (
        gtex_high_pip['chrom'] + ':' +
        gtex_high_pip['pos'].astype(str) + ':' +
        gtex_high_pip['ref'] + ':' +
        gtex_high_pip['alt']
    )

    # Merge GTEx info with MPAC rows (a variant may be in multiple enhancers)
    mpac_cols = ['variant_id', 'id'] + [f'{ct}_skew_pred' for ct in cell_types]
    gtex_cols = [
        'variant_id', 'pip', 'af_gtex', 'cs_id', 'cs_size',
        'afc', 'afc_se', 'tissue', 'phenotype_id', 'gene_name', 'biotype'
    ]
    merged = gtex_high_pip[gtex_cols].merge(
        lead_mpac_df[mpac_cols],
        on='variant_id',
        how='inner'
    )
    merged = merged.rename(columns={'id': 'enhancer_id'})

    if len(merged) == 0:
        return pd.DataFrame()

    # Melt to long form: one row per (variant, enhancer, cell_type)
    skew_cols = [f'{ct}_skew_pred' for ct in cell_types]
    id_vars = [c for c in merged.columns if c not in skew_cols]
    long = merged.melt(
        id_vars=id_vars,
        value_vars=skew_cols,
        var_name='_skew_col',
        value_name='lead_skew'
    )
    long['cell_type'] = long['_skew_col'].str.replace('_skew_pred', '')
    long = long.drop(columns=['_skew_col'])
    long = long.dropna(subset=['lead_skew'])
    long['lead_abs_skew'] = long['lead_skew'].abs()

    # Per (enhancer, cell_type, tissue): keep the variant with the highest lead_abs_skew
    idx = long.groupby(['enhancer_id', 'cell_type', 'tissue'])['lead_abs_skew'].idxmax()
    lead_df = long.loc[idx].copy()
    lead_df = lead_df.rename(columns={
        'variant_id': 'lead_variant_id',
        'pip': 'lead_pip',
        'tissue': 'lead_tissue',
    })
    lead_df = lead_df.drop(columns=['lead_skew'])

    return lead_df.reset_index(drop=True)


def build_variant_rows(all_mpac_df, lead_df, cell_types):
    """
    For all MPAC variants in the high-PIP enhancers, build output rows
    annotated with lead info and exceeds_lead_skew per (variant, cell_type).

    Returns a DataFrame.
    """
    if len(all_mpac_df) == 0 or len(lead_df) == 0:
        return pd.DataFrame()

    all_mpac_df = all_mpac_df.copy()
    all_mpac_df['variant_id'] = (
        all_mpac_df['chrom'] + ':' +
        all_mpac_df['pos'].astype(str) + ':' +
        all_mpac_df['ref'] + ':' +
        all_mpac_df['alt']
    )

    # Melt MPAC to long form: one row per (variant, cell_type)
    skew_cols = [f'{ct}_skew_pred' for ct in cell_types]
    mpac_id_vars = ['variant_id', 'chrom', 'pos', 'ref', 'alt', 'id']
    mpac_long = all_mpac_df[mpac_id_vars + skew_cols].melt(
        id_vars=mpac_id_vars,
        value_vars=skew_cols,
        var_name='_skew_col',
        value_name='skew_pred'
    )
    mpac_long['cell_type'] = mpac_long['_skew_col'].str.replace('_skew_pred', '')
    mpac_long = mpac_long.drop(columns=['_skew_col'])
    mpac_long = mpac_long.dropna(subset=['skew_pred'])
    mpac_long = mpac_long.rename(columns={'id': 'enhancer_id'})
    mpac_long['abs_skew'] = mpac_long['skew_pred'].abs()

    # Join with lead info
    lead_cols = [
        'enhancer_id', 'cell_type', 'lead_variant_id', 'lead_abs_skew',
        'lead_pip', 'lead_tissue', 'phenotype_id', 'gene_name', 'biotype',
        'afc', 'afc_se', 'cs_id', 'cs_size', 'af_gtex'
    ]
    df = mpac_long.merge(
        lead_df[lead_cols],
        on=['enhancer_id', 'cell_type'],
        how='inner'
    )

    df['exceeds_lead_skew'] = df['abs_skew'] > df['lead_abs_skew']
    df['is_lead_variant'] = df['variant_id'] == df['lead_variant_id']

    # Rename for output compatibility
    df = df.rename(columns={'enhancer_id': 'enhancer_ids'})

    return df


def main():
    parser = argparse.ArgumentParser(
        description='High-PIP shadow variant analysis per chromosome'
    )
    parser.add_argument('--chrom', required=True, help='Chromosome (e.g., chr22)')
    parser.add_argument(
        '--pip-threshold', type=float, default=0.9,
        help='PIP threshold for lead variants (default: 0.9)'
    )
    parser.add_argument('--output', default=None, help='Output TSV path')
    args = parser.parse_args()

    if args.output is None:
        args.output = (
            f'results_final/shadow_per_chrom/{args.chrom}_shadow_variants.tsv'
        )

    # --- Paths ---
    base_path = (
        '/projects/tewhey-lab/buttsj/Variant_Effects/'
        'revision_experiments/distal_cre_sat_mut/encode_cCRE_all'
    )
    finemap_dir = (
        '/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/'
        'promoter_sat_mut_comp/raw_data/gtex_v10_fineMapping/SuSiE_fineMapped'
    )
    mpac_template = (
        f'{base_path}/mpac_preds/GRCh38-dELS-{{chrom}}-ALL-mpac-017.tsv.gz'
    )
    gnomad_template = (
        '/projects/tewhey-lab/buttsj/Variant_Effects/gnomad/'
        'gnomad_v4/filtered_data/popLevel/{chrom}_gnomAD_v4_pass_popLevel_af.tsv.gz'
    )
    cell_types = ['k562', 'hepg2', 'sknsh']

    # --- Step 1: Load GTEx high-PIP variants for this chromosome ---
    print(f"Loading GTEx high-PIP variants for {args.chrom} "
          f"(PIP >= {args.pip_threshold})...")
    gtex_high_pip = load_gtex_high_pip(finemap_dir, args.chrom, args.pip_threshold)
    if len(gtex_high_pip) == 0:
        print(f"  No high-PIP variants found for {args.chrom}. Exiting.")
        return
    print(f"  {len(gtex_high_pip):,} high-PIP variant-tissue rows "
          f"({gtex_high_pip['variant_id'].nunique():,} unique variants, "
          f"{gtex_high_pip['tissue'].nunique()} tissues)")

    # --- Step 2: Find those variants in MPAC (resolves enhancer ID + skew) ---
    mpac_path = mpac_template.format(chrom=args.chrom)
    if not os.path.exists(mpac_path):
        print(f"  MPAC file not found: {mpac_path}. Exiting.")
        return

    print(f"Looking up high-PIP variants in MPAC...")
    positions = set(zip(gtex_high_pip['chrom'], gtex_high_pip['pos']))
    lead_mpac_df = extract_mpac_by_position(mpac_path, positions)
    if len(lead_mpac_df) == 0:
        print(f"  No high-PIP variants found in MPAC for {args.chrom}. Exiting.")
        return
    print(f"  {len(lead_mpac_df):,} MPAC rows for "
          f"{lead_mpac_df['id'].nunique():,} unique enhancers")

    # --- Step 3: Build per-(enhancer, cell_type) lead lookup ---
    print("Building lead lookup...")
    lead_df = build_lead_lookup(gtex_high_pip, lead_mpac_df, cell_types)
    if len(lead_df) == 0:
        print("  No (enhancer, cell_type) leads found. Exiting.")
        return
    enhancer_ids = set(lead_df['enhancer_id'])
    print(f"  {len(enhancer_ids):,} enhancers, "
          f"{len(lead_df):,} (enhancer, cell_type) pairs")

    # --- Step 4: Extract all MPAC variants in those enhancers ---
    print(f"Extracting all MPAC variants for {len(enhancer_ids):,} enhancers...")
    all_mpac_df = extract_mpac_by_enhancer(mpac_path, enhancer_ids)
    if len(all_mpac_df) == 0:
        print("  No MPAC rows extracted. Exiting.")
        return
    print(f"  {len(all_mpac_df):,} variant rows across "
          f"{all_mpac_df['id'].nunique():,} enhancers")

    # --- Step 5: Build annotated variant rows ---
    print("Building annotated variant rows...")
    df = build_variant_rows(all_mpac_df, lead_df, cell_types)
    if len(df) == 0:
        print("  No variant rows built. Exiting.")
        return

    pre_dedup = len(df)
    df = df.drop_duplicates(subset=['variant_id', 'cell_type', 'enhancer_ids', 'lead_tissue'])
    print(f"  {pre_dedup:,} → {len(df):,} rows after dedup")

    # --- Step 6: Query gnomAD AF ---
    print("Querying gnomAD v4 allele frequencies...")
    df = query_gnomad_af(df, gnomad_template)
    af_found = df['af'].notna().sum()
    print(f"  AF found for {af_found:,}/{len(df):,} variants")

    # Keep only variants present in gnomAD
    df = df[df['af'].notna()].reset_index(drop=True)
    print(f"  After AF filter: {len(df):,} rows")

    if len(df) == 0:
        print("No variants with gnomAD AF found. Exiting.")
        return

    # --- Reorder columns ---
    col_order = [
        'variant_id', 'chrom', 'pos', 'ref', 'alt',
        'cell_type', 'skew_pred', 'abs_skew', 'enhancer_ids',
        'af', 'AF_afr', 'AF_ami', 'AF_amr', 'AF_asj', 'AF_eas',
        'AF_fin', 'AF_mid', 'AF_nfe', 'AF_remaining', 'AF_sas',
        'exceeds_lead_skew', 'is_lead_variant',
        'lead_variant_id', 'lead_abs_skew', 'lead_pip',
        'lead_tissue', 'phenotype_id', 'gene_name', 'biotype',
        'af_gtex', 'afc', 'afc_se', 'cs_id', 'cs_size',
    ]
    df = df[[c for c in col_order if c in df.columns]]

    # --- Save ---
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    df.to_csv(args.output, sep='\t', index=False)

    n_exceeds = df['exceeds_lead_skew'].sum()
    n_control = (~df['exceeds_lead_skew']).sum()
    print(f"\nSaved {len(df):,} variants to {args.output}")
    print(f"  Unique variants:  {df['variant_id'].nunique():,}")
    print(f"  Unique enhancers: {df['enhancer_ids'].nunique():,}")
    print(f"  Exceeds lead (cases): {n_exceeds:,} "
          f"({100 * n_exceeds / len(df):.1f}%)")
    print(f"  Controls:             {n_control:,} "
          f"({100 * n_control / len(df):.1f}%)")
    print(f"  Per cell type:")
    for ct, count in df.groupby('cell_type').size().items():
        ct_df = df[df['cell_type'] == ct]
        print(f"    {ct}: {count:,} "
              f"({ct_df['exceeds_lead_skew'].sum():,} exceed lead)")


if __name__ == '__main__':
    main()
