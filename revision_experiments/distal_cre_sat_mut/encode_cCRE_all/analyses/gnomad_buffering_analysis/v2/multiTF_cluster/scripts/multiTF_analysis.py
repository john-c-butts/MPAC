#!/usr/bin/env python3
"""
Multi-TF enhancer analysis script for cluster processing.
Processes one chromosome at a time and saves results to pickle.

Usage:
    python multiTF_analysis.py --chrom chr22 --output results_per_chrom/chr22_multiTF.pkl
"""

import argparse
import os
import pickle
import pandas as pd
import pybedtools
from tqdm import tqdm


def load_vierstra_dict(path_to_excel):
    """Load Vierstra motif to cluster mapping."""
    vierstra_motifs = pd.read_excel(path_to_excel, sheet_name=[0, 1])
    idName_dict = dict(zip(vierstra_motifs[0]['Cluster_ID'], vierstra_motifs[0]['Name']))
    vierstra_motifs[1].loc[:, 'cluster_name'] = [idName_dict.get(i) for i in vierstra_motifs[1]['Cluster_ID']]
    return dict(zip(vierstra_motifs[1]['Motif'], vierstra_motifs[1]['cluster_name']))


def df2bed(emvar_df):
    """Convert emVar DataFrame to BED format (vectorized)."""
    df = emvar_df.copy()
    df['bed_id'] = df['chrom'] + ':' + df['pos'].astype(str) + ':' + df['ref'] + ':' + df['alt'] + '_' + df['id']

    df['bed_start'] = df['pos'] - 1
    df['bed_end'] = df['pos']

    ins_mask = df['ref'].str.len() < df['alt'].str.len()
    df.loc[ins_mask, 'bed_start'] = df.loc[ins_mask, 'pos']
    df.loc[ins_mask, 'bed_end'] = df.loc[ins_mask, 'pos']

    del_mask = df['ref'].str.len() > df['alt'].str.len()
    df.loc[del_mask, 'bed_start'] = df.loc[del_mask, 'pos']
    df.loc[del_mask, 'bed_end'] = df.loc[del_mask, 'pos'] + df.loc[del_mask, 'ref'].str.len() - 1

    bed_df = df[['chrom', 'bed_start', 'bed_end', 'bed_id']].drop_duplicates()
    bed_df.columns = [0, 1, 2, 3]
    return pybedtools.BedTool.from_dataframe(bed_df).sort()


def create_long_form_tf_df(raw_bed_df, vierstra_dict):
    """Create long-form TF dataframe with collapse to representative TF per interval."""
    raw_bed_df = raw_bed_df.copy()
    raw_bed_df[3] = raw_bed_df['name'].astype(str)
    df_split = raw_bed_df.assign(tf_hits_list=raw_bed_df[3].str.split(';'))
    long_df = df_split.explode('tf_hits_list').reset_index(drop=True)

    long_df['full_hit_string'] = long_df['tf_hits_list']
    long_df['hocomoco_tf'] = long_df['full_hit_string'].str.split('_EH').str[0]
    long_df['representative_tf'] = long_df['hocomoco_tf'].str.split('_').str[0]
    long_df['enhancer_id'] = long_df['full_hit_string'].str.split('_').str[2]

    long_df['rep_tf_contrib'] = pd.to_numeric(
        long_df['full_hit_string'].str.split('_').str[-1], errors='coerce'
    ).fillna(0.0)

    long_df['vierstra_cluster'] = long_df['hocomoco_tf'].map(vierstra_dict)
    long_df['activity_class'] = pd.cut(
        long_df['rep_tf_contrib'],
        bins=[-float('inf'), 0, float('inf')],
        labels=['Repressor', 'Activator']
    ).astype(str)
    long_df.loc[long_df['rep_tf_contrib'] == 0, 'activity_class'] = 'Neutral'

    long_df = long_df.drop(columns=['tf_hits_list', 3])

    # Collapse to representative TF per interval (highest |contribution|)
    long_df['abs_contrib'] = long_df['rep_tf_contrib'].abs()
    idx_max = long_df.groupby(['chrom', 'start', 'end'])['abs_contrib'].idxmax()
    collapsed_df = long_df.loc[idx_max].drop(columns=['abs_contrib']).reset_index(drop=True)

    return collapsed_df


def get_multiTF_info(collapsed_seqlets, max_overlap_bp=2):
    """
    Get multi-TF info per enhancer with non-overlap filtering.

    Returns:
        dict: {enhancer_id: {'valid_clusters': [list], 'valid_seqlet_indices': [list]}}
    """
    results = {}

    for (enh_id, cluster), group in collapsed_seqlets.groupby(['enhancer_id', 'vierstra_cluster']):
        if len(group) < 2 or pd.isna(cluster):
            continue

        sorted_group = group.sort_values('start')
        is_valid = True
        prev_end = None
        for _, row in sorted_group.iterrows():
            if prev_end is not None:
                overlap = prev_end - row['start']
                if overlap > max_overlap_bp:
                    is_valid = False
                    break
            prev_end = row['end']

        if is_valid:
            if enh_id not in results:
                results[enh_id] = {'valid_clusters': [], 'valid_seqlet_indices': []}
            results[enh_id]['valid_clusters'].append(cluster)
            results[enh_id]['valid_seqlet_indices'].extend(sorted_group.index.tolist())

    return results


def process_emvar_seqlets(emvar_seqlets_df, mpac_df, cell_type):
    """Add clean variant_id column and merge with predictions."""
    if len(emvar_seqlets_df) == 0:
        return emvar_seqlets_df

    df = emvar_seqlets_df.copy()
    df['variant_id'] = df['blockStarts'].str.split('_').str[0]

    variant_parts = df['variant_id'].str.split(':', expand=True)
    df['var_chrom'] = variant_parts[0]
    df['var_pos'] = variant_parts[1].astype(int)
    df['var_ref'] = variant_parts[2]
    df['var_alt'] = variant_parts[3]

    skew_col = f'{cell_type}_skew_pred'
    mpac_subset = mpac_df[['chrom', 'pos', 'ref', 'alt', skew_col]].copy()
    mpac_subset = mpac_subset.rename(columns={skew_col: 'skew_pred'})

    df = df.merge(
        mpac_subset,
        left_on=['var_chrom', 'var_pos', 'var_ref', 'var_alt'],
        right_on=['chrom', 'pos', 'ref', 'alt'],
        how='left',
        suffixes=('', '_mpac')
    )
    df = df.drop(columns=['chrom_mpac', 'pos_mpac', 'ref_mpac', 'alt_mpac'], errors='ignore')

    return df


def create_variant_dataframe(emvar_seqlets_df, multiTF_enhancer_ids, cell_type):
    """Create variant-level dataframe."""
    if len(emvar_seqlets_df) == 0:
        return pd.DataFrame()

    df = emvar_seqlets_df.copy()

    variant_agg = df.groupby('variant_id').agg({
        'thickEnd': lambda x: ','.join(sorted(set(x))),
        'var_chrom': 'first',
        'var_pos': 'first',
        'var_ref': 'first',
        'var_alt': 'first',
        'skew_pred': 'first'
    }).reset_index()
    variant_agg.columns = ['variant_id', 'enhancer_ids', 'chrom', 'pos', 'ref', 'alt', 'skew_pred']

    def check_multiTF(enh_str):
        return any(e in multiTF_enhancer_ids for e in enh_str.split(','))
    variant_agg['is_multiTF'] = variant_agg['enhancer_ids'].apply(check_multiTF)
    variant_agg['cell_type'] = cell_type

    return variant_agg[['variant_id', 'chrom', 'pos', 'ref', 'alt', 'cell_type',
                        'skew_pred', 'enhancer_ids', 'is_multiTF']]


def multiTF_analysis_single_chrom(chrom, paths, vierstra_dict, max_overlap_bp=2):
    """
    Run multi-TF analysis for a single chromosome.

    Returns nested structure:
    {
        'chrom': str,
        'variants_all': pd.DataFrame,
        'multiTF_summary': pd.DataFrame,
        'multiTF_enhancers': {
            cell_type: {
                enhancer_id: {
                    'valid_clusters': list,
                    'filtered_seqlets': pd.DataFrame,
                    'all_seqlets': pd.DataFrame,
                    'emvars': pd.DataFrame
                }
            }
        }
    }
    """
    print(f"Processing {chrom}...")

    # Read MPAC predictions for this chromosome
    print(f"  Reading MPAC predictions...")
    mpac_path = paths['mpac_template'].format(chrom=chrom)
    allPreds = pd.read_csv(mpac_path, sep='\t')
    print(f"    Loaded {len(allPreds):,} variants")

    # Filter for emVars in each cell type
    print(f"  Filtering for emVars...")
    emvars = {
        'k562': allPreds[allPreds['k562_skew_pred'].abs() > 0.5].copy(),
        'hepg2': allPreds[allPreds['hepg2_skew_pred'].abs() > 0.5].copy(),
        'sknsh': allPreds[allPreds['sknsh_skew_pred'].abs() > 0.5].copy()
    }

    # Convert to BED
    print(f"  Converting emVars to BED...")
    emvar_beds = {ct: df2bed(df) for ct, df in emvars.items()}

    # Load seqlet BED files
    print(f"  Loading seqlet BED files...")
    seqlet_beds = {
        'k562': pybedtools.BedTool(paths['seqlets_k562']),
        'hepg2': pybedtools.BedTool(paths['seqlets_hepg2']),
        'sknsh': pybedtools.BedTool(paths['seqlets_sknsh'])
    }

    # Results storage
    all_variants = []
    all_summaries = []
    multiTF_enhancers = {}

    for cell_type in ['k562', 'hepg2', 'sknsh']:
        print(f"\n  Processing {cell_type}...")
        multiTF_enhancers[cell_type] = {}

        print(f"    Intersecting emVars with seqlets...")
        emvar_seqlets_raw = seqlet_beds[cell_type].intersect(
            emvar_beds[cell_type], wa=True, wb=True
        ).to_dataframe()

        print(f"    Processing emvar_seqlets...")
        emvar_seqlets = process_emvar_seqlets(emvar_seqlets_raw, emvars[cell_type], cell_type)

        print(f"    Creating collapsed seqlets...")
        collapsed_seqlets = create_long_form_tf_df(seqlet_beds[cell_type].to_dataframe(), vierstra_dict)

        print(f"    Finding multi-TF enhancers (max {max_overlap_bp}bp overlap)...")
        multiTF_info = get_multiTF_info(collapsed_seqlets, max_overlap_bp)
        multiTF_enhancer_ids = set(multiTF_info.keys())

        print(f"    Creating variant dataframe...")
        variant_df = create_variant_dataframe(emvar_seqlets, multiTF_enhancer_ids, cell_type)
        all_variants.append(variant_df)

        # Build nested structure for each valid multi-TF enhancer
        print(f"    Building nested enhancer structure...")
        emvar_enhancer_ids = set(emvar_seqlets['thickEnd'].unique()) if len(emvar_seqlets) > 0 else set()
        multiTF_with_emvars = multiTF_enhancer_ids & emvar_enhancer_ids

        summary_rows = []
        for enhancer_id in multiTF_with_emvars:
            info = multiTF_info[enhancer_id]
            valid_clusters = info['valid_clusters']
            valid_indices = info['valid_seqlet_indices']

            # Get all seqlets for this enhancer
            all_seqlets = collapsed_seqlets[collapsed_seqlets['enhancer_id'] == enhancer_id].copy()

            # Get only seqlets from valid clusters
            filtered_seqlets = collapsed_seqlets.loc[
                collapsed_seqlets.index.isin(valid_indices) &
                (collapsed_seqlets['enhancer_id'] == enhancer_id)
            ].copy()

            # Get emvars for this enhancer
            enh_emvars = emvar_seqlets[emvar_seqlets['thickEnd'] == enhancer_id].copy()

            # Store in nested structure
            multiTF_enhancers[cell_type][enhancer_id] = {
                'valid_clusters': valid_clusters,
                'filtered_seqlets': filtered_seqlets,
                'all_seqlets': all_seqlets,
                'emvars': enh_emvars
            }

            # Build summary row
            summary_rows.append({
                'enhancer_id': enhancer_id,
                'cell_type': cell_type,
                'multiTF_clusters': ','.join(sorted(valid_clusters)),
                'n_multiTF_clusters': len(valid_clusters),
                'n_filtered_seqlets': len(filtered_seqlets),
                'n_all_seqlets': len(all_seqlets),
                'n_emvars': len(enh_emvars),
                'chrom': all_seqlets['chrom'].iloc[0] if len(all_seqlets) > 0 else None,
                'start': all_seqlets['start'].min() if len(all_seqlets) > 0 else None,
                'end': all_seqlets['end'].max() if len(all_seqlets) > 0 else None
            })

        summary_df = pd.DataFrame(summary_rows)
        if len(summary_df) > 0:
            summary_df = summary_df.sort_values('n_emvars', ascending=False).reset_index(drop=True)
        all_summaries.append(summary_df)

        print(f"    Multi-TF enhancers (non-overlapping): {len(multiTF_enhancer_ids):,}")
        print(f"    Multi-TF enhancers with emVars: {len(multiTF_with_emvars):,}")
        print(f"    Variants in output: {len(variant_df):,}")

    # Combine results
    results = {
        'chrom': chrom,
        'variants_all': pd.concat(all_variants, ignore_index=True) if all_variants else pd.DataFrame(),
        'multiTF_summary': pd.concat(all_summaries, ignore_index=True) if all_summaries else pd.DataFrame(),
        'multiTF_enhancers': multiTF_enhancers
    }

    return results


def main():
    parser = argparse.ArgumentParser(description='Multi-TF enhancer analysis for a single chromosome')
    parser.add_argument('--chrom', required=True, help='Chromosome to process (e.g., chr22)')
    parser.add_argument('--output', required=True, help='Output pickle file path')
    parser.add_argument('--max-overlap-bp', type=int, default=2, help='Max overlap between seqlets (default: 2)')
    args = parser.parse_args()

    # Define paths
    base_path = '/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all'
    paths = {
        'mpac_template': f'{base_path}/mpac_preds/GRCh38-dELS-{{chrom}}-ALL-mpac-017.tsv.gz',
        'seqlets_k562': f'{base_path}/analyses/gnomad_buffering_analysis/../../processed_data/bed_files/repTF_k562_dELS_seqlets_01.bed',
        'seqlets_hepg2': f'{base_path}/analyses/gnomad_buffering_analysis/../../processed_data/bed_files/repTF_hepg2_dELS_seqlets_01.bed',
        'seqlets_sknsh': f'{base_path}/analyses/gnomad_buffering_analysis/../../processed_data/bed_files/repTF_sknsh_dELS_seqlets_01.bed',
        'vierstra': f'{base_path}/analyses/gnomad_buffering_analysis/motif_annotations.xlsx'
    }

    # Load Vierstra dictionary
    print("Loading Vierstra motif annotations...")
    vierstra_dict = load_vierstra_dict(paths['vierstra'])

    # Run analysis
    results = multiTF_analysis_single_chrom(
        chrom=args.chrom,
        paths=paths,
        vierstra_dict=vierstra_dict,
        max_overlap_bp=args.max_overlap_bp
    )

    # Save results
    print(f"\nSaving results to {args.output}...")
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, 'wb') as f:
        pickle.dump(results, f)

    print("Done!")

    # Print summary
    print(f"\n{'='*50}")
    print(f"Summary for {args.chrom}:")
    print(f"  Total variants: {len(results['variants_all']):,}")
    print(f"  Multi-TF enhancers with emVars:")
    for ct in ['k562', 'hepg2', 'sknsh']:
        print(f"    {ct}: {len(results['multiTF_enhancers'][ct]):,}")
    print(f"{'='*50}")


if __name__ == '__main__':
    main()
