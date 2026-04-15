#!/usr/bin/env python3
"""
Stage 1: Build per-chrom candidate variant files for filtered openTargets analysis.

Loads high-PIP openTargets GWAS credible sets, intersects lead variants with multiTF
emVars, iterates through filtered seqlets (±4bp pad) per enhancer/cell-type, and
outputs per-chrom candidate TSVs for parallel gnomAD annotation in Stage 2.

Filtering applied here:
  - Lead variant PIP >= 0.9 (from exploded_chunks/highPIP_09/ pre-filter)
  - Lead variant must be an emVar in at least one filtered seqlet (±4bp pad)
  - All other seqlet-overlapping emVars are recorded regardless of AF or effect size
    (those filters are applied in Stage 2 after gnomAD annotation)

Also saves results_final/openTargets_phenotype_info.tsv for Stage 3 merging.

Usage:
    python scripts/openTargets_build_candidates.py
    (run from multiTF_cluster/ directory)
"""

import os
import pickle

import numpy as np
import pandas as pd
from tqdm import tqdm


def main():
    path2highPip = 'exploded_chunks/highPIP_09'
    multiTF_variants_path = 'results_final/multiTF_variants.tsv'
    multiTF_summary_path = 'results_final/multiTF_summary.tsv'
    pickle_dir = 'results_per_chrom'
    output_dir = 'results_final/ot_candidates_per_chrom'
    os.makedirs(output_dir, exist_ok=True)

    # --- Load high-PIP openTargets variants ---
    print("Loading high-PIP openTargets variants...")
    chunks = []
    for f in tqdm(sorted(os.listdir(path2highPip))):
        if f.endswith('.parquet'):
            chunks.append(pd.read_parquet(f'{path2highPip}/{f}'))
    highPipAll = pd.concat(chunks, ignore_index=True)
    print(f"  {len(highPipAll):,} high-PIP variant rows loaded")

    # --- Load multiTF emVars ---
    print("Loading multiTF variants...")
    all_emVars = pd.read_csv(multiTF_variants_path, sep='\t')
    multiTF_emVars = all_emVars[all_emVars['is_multiTF']].copy()
    multiTF_emVars['openTargets_id'] = [
        '_'.join([chrom.split('chr')[-1], str(pos), ref, alt])
        for chrom, pos, ref, alt in zip(
            multiTF_emVars['chrom'], multiTF_emVars['pos'],
            multiTF_emVars['ref'], multiTF_emVars['alt']
        )
    ]
    print(f"  {len(multiTF_emVars):,} multiTF emVar rows")

    # --- Filter openTargets to lead variants that are multiTF emVars ---
    ot_emVars = highPipAll[
        highPipAll['leadVariant'].isin(multiTF_emVars['openTargets_id'])
    ]
    n_leads = ot_emVars['leadVariant'].nunique()
    print(f"  {n_leads} unique lead variants overlap multiTF emVars")

    # --- Save phenotype info for Stage 3 ---
    pheno_cols = [
        'leadVariant', 'studyId', 'traitFromSource', 'posteriorProbability',
        'is95CredibleSet', 'finemappingMethod'
    ]
    phenotype_info = ot_emVars[pheno_cols].drop_duplicates()
    pheno_path = 'results_final/openTargets_phenotype_info.tsv'
    phenotype_info.to_csv(pheno_path, sep='\t', index=False)
    print(f"  Saved phenotype info: {len(phenotype_info):,} lead-study pairs → {pheno_path}")

    # --- Build lookup tables ---
    leadVar_enhID = dict(zip(multiTF_emVars['openTargets_id'], multiTF_emVars['enhancer_ids']))
    leadVar_chrom = dict(zip(multiTF_emVars['openTargets_id'], multiTF_emVars['chrom']))
    emvar_groups = multiTF_emVars.groupby(['cell_type', 'enhancer_ids'])

    # --- Group lead variants by chromosome ---
    lead_by_chrom = {}
    for lv in ot_emVars['leadVariant'].unique():
        chrom = leadVar_chrom.get(lv)
        if chrom:
            lead_by_chrom.setdefault(chrom, []).append(lv)
    print(f"  Lead variants span {len(lead_by_chrom)} chromosomes")

    # --- Load multiTF summary for cell-type lookup ---
    multiTF_summary = pd.read_csv(multiTF_summary_path, sep='\t')

    # --- Process each chromosome ---
    total_candidates = 0
    for chrom in tqdm(sorted(lead_by_chrom.keys(), key=lambda x: int(x.replace('chr', ''))),
                      desc='Chromosomes'):
        lead_vars = lead_by_chrom[chrom]
        pkl_path = os.path.join(pickle_dir, f'{chrom}_multiTF.pkl')
        if not os.path.exists(pkl_path):
            print(f"  {chrom}: pickle not found, skipping")
            continue

        with open(pkl_path, 'rb') as f:
            chrom_data = pickle.load(f)

        all_rows = []

        for leadVar in lead_vars:
            enhancer = leadVar_enhID.get(leadVar)
            if not enhancer:
                continue
            enhancer_cells = list(
                multiTF_summary[multiTF_summary['enhancer_id'] == enhancer]['cell_type'].unique()
            )

            for cell in enhancer_cells:
                ct_data = chrom_data['multiTF_enhancers'].get(cell, {})
                if enhancer not in ct_data:
                    continue

                filtered_seqlets = ct_data[enhancer]['filtered_seqlets']

                try:
                    cell_enh_vars = emvar_groups.get_group((cell, enhancer))
                except KeyError:
                    continue

                seqlet_rows = []
                tf_count = 1
                lead_found = False

                for _, seqlet in filtered_seqlets.iterrows():
                    start = seqlet['start']
                    end = seqlet['end']
                    rep_tf = seqlet['vierstra_cluster']
                    contrib = seqlet['rep_tf_contrib']

                    in_seqlet = cell_enh_vars[
                        (cell_enh_vars['pos'] >= start - 4) &
                        (cell_enh_vars['pos'] <= end + 4)
                    ].copy()

                    if len(in_seqlet) == 0:
                        tf_count += 1
                        continue

                    in_seqlet['tf_family'] = rep_tf
                    in_seqlet['tf_contrib'] = contrib
                    in_seqlet['tf_instance'] = tf_count
                    in_seqlet['leadVariant'] = leadVar
                    in_seqlet['is_leadVariant'] = (
                        in_seqlet['openTargets_id'] == leadVar
                    ).astype(int)

                    if in_seqlet['is_leadVariant'].sum() > 0:
                        lead_found = True

                    seqlet_rows.append(in_seqlet)
                    tf_count += 1

                # Only include this (enhancer, cell) if lead was found in a seqlet
                if lead_found and seqlet_rows:
                    all_rows.append(
                        pd.concat(seqlet_rows, ignore_index=True).sort_values('pos')
                    )

        if not all_rows:
            print(f"  {chrom}: no seqlet candidates")
            continue

        candidates = pd.concat(all_rows, ignore_index=True)

        # Compute abs_skew and lead_abs_skew per (enhancer, cell_type, leadVariant)
        candidates['abs_skew'] = candidates['skew_pred'].abs()
        lead_skews = (
            candidates[candidates['is_leadVariant'] == 1]
            .groupby(['enhancer_ids', 'cell_type', 'leadVariant'])['abs_skew']
            .first()
            .rename('lead_abs_skew')
            .reset_index()
        )
        candidates = candidates.merge(
            lead_skews, on=['enhancer_ids', 'cell_type', 'leadVariant'], how='left'
        )
        candidates['exceeds_lead_skew'] = candidates['abs_skew'] > candidates['lead_abs_skew']

        out_path = os.path.join(output_dir, f'{chrom}_ot_candidates.tsv')
        candidates.to_csv(out_path, sep='\t', index=False)
        total_candidates += len(candidates)
        print(f"  {chrom}: {len(candidates):,} candidate variants → {out_path}")

        del chrom_data, all_rows, candidates

    print(f"\nDone. {total_candidates:,} total candidate variants across all chromosomes.")
    print(f"Next step: sbatch scripts/submit_openTargets_gnomad_filter.sh")


if __name__ == '__main__':
    main()
