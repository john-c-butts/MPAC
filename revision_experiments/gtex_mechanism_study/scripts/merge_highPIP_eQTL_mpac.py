#!/usr/bin/env python3
"""
Merge per-chromosome highPIP eQTL × MPAC results and add global annotations.

Annotations added here (require the full cross-chromosome view):
  - tissue_celltype_concordant per cell type
  - n_gene_associations
  - n_tissues_significant
  - cross_tissue_direction_consistent
  - cross_gene_direction_consistent

Usage:
  python merge_highPIP_eQTL_mpac.py \\
      --per_chrom_dir <path>/processed_data/per_chrom \\
      --out <path>/processed_data/highPIP_eQTL_mpac_all.tsv.gz
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# Biological compatibility: (cell_type → set of GTEx tissue name substrings)
# A lead_tissue matches if any substring appears in it (case-insensitive).
TISSUE_CELLTYPE_MAP = {
    "k562":  {"Whole_Blood", "Cells_EBV-transformed_lymphocytes"},
    "hepg2": {"Liver"},
    "sknsh": {"Brain_", "Nerve_Tibial"},
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--per_chrom_dir", required=True)
    p.add_argument("--out", required=True)
    return p.parse_args()


def load_per_chrom(per_chrom_dir: str) -> pd.DataFrame:
    files = sorted(Path(per_chrom_dir).glob("chr*_highPIP_eQTL_mpac.tsv.gz"))
    if not files:
        sys.exit(f"No per-chrom files found in {per_chrom_dir}")
    print(f"Loading {len(files)} per-chrom files...", file=sys.stderr)
    frames = [pd.read_csv(f, sep="\t") for f in files]
    df = pd.concat(frames, ignore_index=True)
    df = df.drop_duplicates()
    print(f"  {len(df)} total rows after dedup", file=sys.stderr)
    return df


def add_tissue_celltype_concordance(df: pd.DataFrame) -> pd.DataFrame:
    for cell_type, compatible in TISSUE_CELLTYPE_MAP.items():
        col = f"{cell_type}_tissue_concordant"
        mask = pd.Series(False, index=df.index)
        for tissue_substr in compatible:
            mask = mask | df["lead_tissue"].str.contains(tissue_substr, case=False, regex=False)
        df[col] = mask
    return df


def compute_global_stats(df: pd.DataFrame) -> pd.DataFrame:
    """
    Computes per-variant and per-variant-gene summary columns.

    Deduplicates on (variant_id, phenotype_id, lead_tissue) before aggregating
    so that variants in multiple enhancers don't inflate counts or skew
    direction-consistency checks.
    """
    # One row per (variant, gene, tissue) — collapse multi-enhancer duplicates
    vgt = df[["variant_id", "phenotype_id", "lead_tissue", "slope"]].drop_duplicates(
        subset=["variant_id", "phenotype_id", "lead_tissue"]
    )

    # n_gene_associations: distinct genes with pip >= 0.9 per variant across all tissues
    # (all rows in our table already satisfy pip >= 0.9 + significance filters)
    n_genes = (
        vgt.groupby("variant_id")["phenotype_id"]
        .nunique()
        .rename("n_gene_associations")
        .reset_index()
    )

    # n_tissues_significant: distinct tissues per (variant, gene)
    n_tissues = (
        vgt.groupby(["variant_id", "phenotype_id"])["lead_tissue"]
        .nunique()
        .rename("n_tissues_significant")
        .reset_index()
    )

    # cross_tissue_direction_consistent: slope sign consistent across tissues for same (variant, gene)
    def sign_consistent(slopes):
        signs = np.sign(slopes.dropna())
        return bool((signs > 0).all() or (signs < 0).all()) if len(signs) > 0 else pd.NA

    cross_tissue = (
        vgt.groupby(["variant_id", "phenotype_id"])["slope"]
        .agg(sign_consistent)
        .rename("cross_tissue_direction_consistent")
        .reset_index()
    )

    # cross_gene_direction_consistent: slope sign consistent across all genes for a variant
    cross_gene = (
        vgt.groupby("variant_id")["slope"]
        .agg(sign_consistent)
        .rename("cross_gene_direction_consistent")
        .reset_index()
    )

    # Join stats back onto main table
    df = df.merge(n_genes, on="variant_id", how="left")
    df = df.merge(n_tissues, on=["variant_id", "phenotype_id"], how="left")
    df = df.merge(cross_tissue, on=["variant_id", "phenotype_id"], how="left")
    df = df.merge(cross_gene, on="variant_id", how="left")

    return df


FINAL_COL_ORDER = [
    # Row identity
    "variant_id", "chrom", "pos", "ref", "alt",
    "enhancer_id", "lead_tissue", "phenotype_id", "gene_name",
    # GTEx finemapping
    "pip", "afc", "afc_se", "slope", "slope_se",
    "cs_id", "cs_size", "af_gtex",
    # MPAC — K562
    "k562_ref_pred", "k562_alt_pred", "k562_skew_pred",
    # MPAC — HepG2
    "hepg2_ref_pred", "hepg2_alt_pred", "hepg2_skew_pred",
    # MPAC — SK-N-SH
    "sknsh_ref_pred", "sknsh_alt_pred", "sknsh_skew_pred",
    # Annotations
    "k562_tissue_concordant", "hepg2_tissue_concordant", "sknsh_tissue_concordant",
    "n_gene_associations", "n_tissues_significant",
    "cross_tissue_direction_consistent", "cross_gene_direction_consistent",
]


def main():
    args = parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    df = load_per_chrom(args.per_chrom_dir)

    print("Adding tissue-celltype concordance...", file=sys.stderr)
    df = add_tissue_celltype_concordance(df)

    print("Computing global annotations...", file=sys.stderr)
    df = compute_global_stats(df)

    df = df[FINAL_COL_ORDER]
    df = df.sort_values(["chrom", "pos", "lead_tissue", "phenotype_id"]).reset_index(drop=True)

    df.to_csv(args.out, sep="\t", index=False, compression="gzip")
    print(f"Wrote {len(df)} rows to {args.out}", file=sys.stderr)

    # Summary stats
    print(f"\nSummary:", file=sys.stderr)
    print(f"  Unique variants:           {df['variant_id'].nunique()}", file=sys.stderr)
    print(f"  Unique enhancers:          {df['enhancer_id'].nunique()}", file=sys.stderr)
    print(f"  Unique genes:              {df['phenotype_id'].nunique()}", file=sys.stderr)
    print(f"  Unique tissues:            {df['lead_tissue'].nunique()}", file=sys.stderr)
    print(f"  Multi-tissue variants:     {(df.groupby('variant_id')['lead_tissue'].nunique() > 1).sum()}", file=sys.stderr)
    print(f"  Multi-gene variants:       {(df['n_gene_associations'] > 1).sum()}", file=sys.stderr)


if __name__ == "__main__":
    main()
