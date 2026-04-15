#!/usr/bin/env python3
"""
Stage 3: Merge per-chrom filtered openTargets results and add phenotype info.

Merges per-chrom TSVs from Stage 2, then joins in the openTargets phenotype/study
metadata (trait, study ID, PIP, etc.) saved by Stage 1.

Only retains loci where at least one shadow_emVar exists (i.e. a variant with a
stronger predicted effect than the lead that is also in gnomAD). Lead-only loci
(where the lead is the strongest variant in the seqlet) are written to a separate
file for reference.

Controls are carried through only for (enhancer_ids, cell_type) pairs that have
at least one qualifying shadow_emVar after filtering.

Outputs:
  results_final/openTargets_GWAS_multiTF_seqletFiltered_with_controls.tsv  (main)
  results_final/openTargets_GWAS_multiTF_leadOnly.tsv                       (excluded)

Usage:
    python scripts/merge_openTargets_seqlet_filtered.py
    (run from multiTF_cluster/ directory)
"""

import os
import pandas as pd
from tqdm import tqdm


def main():
    filtered_dir = 'results_final/ot_filtered_per_chrom'
    pheno_path = 'results_final/openTargets_phenotype_info.tsv'
    output_path = 'results_final/openTargets_GWAS_multiTF_seqletFiltered_with_controls.tsv'
    leadonly_path = 'results_final/openTargets_GWAS_multiTF_leadOnly.tsv'

    # --- Merge per-chrom filtered TSVs ---
    print("Merging per-chrom filtered results...")
    chrom_files = sorted(
        [f for f in os.listdir(filtered_dir) if f.endswith('_ot_filtered.tsv')],
        key=lambda x: int(x.split('chr')[1].split('_')[0])
    )
    if not chrom_files:
        print(f"No filtered TSVs found in {filtered_dir}")
        return

    dfs = []
    for f in tqdm(chrom_files):
        path = os.path.join(filtered_dir, f)
        df = pd.read_csv(path, sep='\t')
        dfs.append(df)
        print(f"  {f}: {len(df):,} rows")

    merged = pd.concat(dfs, ignore_index=True)
    print(f"\nMerged: {len(merged):,} total rows")
    for vc, n in merged['variant_class'].value_counts().items():
        print(f"  {vc}: {n:,}")

    # --- Merge phenotype/study info onto lead and shadow emVars ---
    print("\nMerging phenotype info...")
    pheno = pd.read_csv(pheno_path, sep='\t')
    print(f"  {len(pheno):,} lead-study pairs, {pheno['leadVariant'].nunique()} unique leads")

    # Only emVar rows have a leadVariant column populated; controls do not
    emvar_mask = merged['variant_class'].isin(['lead_emVar', 'shadow_emVar'])
    emVars = merged[emvar_mask].copy()
    controls = merged[~emvar_mask].copy()

    emVars_pheno = emVars.merge(pheno, on='leadVariant', how='left')
    print(f"  emVar rows before pheno merge: {len(emVars):,}")
    print(f"  emVar rows after pheno merge:  {len(emVars_pheno):,}")

    # --- Drop lead-only loci (no shadow emVar exists) ---
    # Group key: (enhancer_ids, cell_type, leadVariant) — one GWAS association per cell context
    group_key = ['enhancer_ids', 'cell_type', 'leadVariant']
    has_shadow = (
        emVars_pheno[emVars_pheno['variant_class'] == 'shadow_emVar']
        .groupby(group_key)
        .size()
        .gt(0)
        .rename('has_shadow')
        .reset_index()
    )
    emVars_annotated = emVars_pheno.merge(has_shadow, on=group_key, how='left')
    emVars_annotated['has_shadow'] = emVars_annotated['has_shadow'].fillna(False)

    emVars_keep = emVars_annotated[emVars_annotated['has_shadow']].drop(columns='has_shadow')
    emVars_drop = emVars_annotated[~emVars_annotated['has_shadow']].drop(columns='has_shadow')

    n_drop_leads = emVars_drop[emVars_drop['variant_class'] == 'lead_emVar']['leadVariant'].nunique()
    print(f"\nLead-only loci (no shadow emVar): {n_drop_leads} unique lead variants → {leadonly_path}")
    emVars_drop.to_csv(leadonly_path, sep='\t', index=False)

    # --- Filter controls to (enhancer, cell_type) pairs still present after dropping ---
    valid_pairs = set(zip(emVars_keep['enhancer_ids'], emVars_keep['cell_type']))
    ctrl_mask = [
        (eid, ct) in valid_pairs
        for eid, ct in zip(controls['enhancer_ids'], controls['cell_type'])
    ]
    controls_keep = controls[ctrl_mask]
    n_ctrl_dropped = (~pd.Series(ctrl_mask)).sum()
    if n_ctrl_dropped:
        print(f"  Dropped {n_ctrl_dropped:,} control rows from lead-only (enhancer, cell_type) pairs")

    # --- Combine and save main output ---
    final = pd.concat([emVars_keep, controls_keep], ignore_index=True, sort=False)
    print(f"\nFinal combined rows: {len(final):,}")
    print(f"  Unique lead variants: {final['leadVariant'].nunique()}")
    print(f"  Unique enhancers:     {final['enhancer_ids'].nunique()}")
    print(f"  Variant class breakdown:")
    for vc, n in final['variant_class'].value_counts().items():
        print(f"    {vc}: {n:,}")

    final.to_csv(output_path, sep='\t', index=False)
    print(f"\nSaved to {output_path}")


if __name__ == '__main__':
    main()
