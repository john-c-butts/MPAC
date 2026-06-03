#!/usr/bin/env python3
"""
Per-chromosome worker for the highPIP eQTL × MPAC pipeline.

For a given chromosome, identifies GTEx v10 variants that are:
  - SuSiE PIP >= 0.9
  - Gene is a significant eGene in that tissue (qval < 0.05)
  - Variant-gene pair present in that tissue's significant pairs
  - Variant overlaps a dELS enhancer in MPAC predictions

Outputs one row per (variant, enhancer, tissue, gene) with GTEx finemapping
stats and all three cell-type MPAC predictions.

Usage (SLURM array, one task per chromosome):
  python highPIP_eQTL_mpac_per_chrom.py \\
      --chrom chr1 \\
      --susie_dir  <path>/SuSiE_fineMapped \\
      --signif_dir <path>/GTEx_Analysis_v10_eQTL_updated \\
      --mpac_dir   <path>/mpac_preds \\
      --out_dir    <path>/processed_data/per_chrom
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd


MPAC_PRED_COLS = [
    "k562_ref_pred", "k562_alt_pred", "k562_skew_pred",
    "hepg2_ref_pred", "hepg2_alt_pred", "hepg2_skew_pred",
    "sknsh_ref_pred", "sknsh_alt_pred", "sknsh_skew_pred",
]

OUTPUT_COLS = [
    "variant_id", "chrom", "pos", "ref", "alt",
    "enhancer_id", "lead_tissue", "phenotype_id", "gene_name",
    "pip", "afc", "afc_se", "slope", "slope_se",
    "cs_id", "cs_size", "af_gtex",
] + MPAC_PRED_COLS


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--chrom", required=True, help="e.g. chr1")
    p.add_argument("--susie_dir", required=True)
    p.add_argument("--signif_dir", required=True)
    p.add_argument("--mpac_dir", required=True)
    p.add_argument("--out_dir", required=True)
    return p.parse_args()


def load_mpac(mpac_dir: str, chrom: str) -> pd.DataFrame:
    path = os.path.join(mpac_dir, f"GRCh38-dELS-{chrom}-ALL-mpac-017.tsv.gz")
    df = pd.read_csv(path, sep="\t", dtype={"pos": np.int64})
    df = df.rename(columns={"id": "enhancer_id"})
    return df


def get_tissues(susie_dir: str) -> list[str]:
    return sorted(
        p.name.split(".")[0]
        for p in Path(susie_dir).glob("*.v10.eQTLs.SuSiE_summary.parquet")
    )


def load_tissue_data(tissue: str, chrom: str, susie_dir: str, signif_dir: str):
    """
    Returns a DataFrame of high-PIP, significance-filtered variant-gene pairs
    for one tissue on one chromosome, with slope from signif_pairs joined in.
    Returns None if no qualifying rows exist.
    """
    susie_path = Path(susie_dir) / f"{tissue}.v10.eQTLs.SuSiE_summary.parquet"
    egenes_path = Path(signif_dir) / f"{tissue}.v10.eGenes.txt.gz"
    signif_path = Path(signif_dir) / f"{tissue}.v10.eQTLs.signif_pairs.parquet"

    for p in (susie_path, egenes_path, signif_path):
        if not p.exists():
            print(f"  [SKIP] missing file: {p}", file=sys.stderr)
            return None

    # SuSiE: pip >= 0.9, current chromosome only
    susie = pd.read_parquet(susie_path)
    susie = susie[
        susie["variant_id"].str.startswith(chrom + "_") & (susie["pip"] >= 0.9)
    ].copy()
    if susie.empty:
        return None
    susie = susie.rename(columns={"af": "af_gtex"})
    susie = susie.drop(columns=["biotype"], errors="ignore")

    # eGenes: keep only genes with qval < 0.05
    egenes = pd.read_csv(egenes_path, sep="\t", usecols=["gene_id", "qval"])
    sig_genes = set(egenes.loc[egenes["qval"] < 0.05, "gene_id"])
    susie = susie[susie["phenotype_id"].isin(sig_genes)]
    if susie.empty:
        return None

    # signif_pairs: presence check + get slope/slope_se; filter to chrom first
    signif = pd.read_parquet(
        signif_path, columns=["gene_id", "variant_id", "slope", "slope_se"]
    )
    signif = signif[signif["variant_id"].str.startswith(chrom + "_")]

    # Inner join: variant-gene pair must be in signif_pairs
    merged = susie.merge(
        signif,
        left_on=["variant_id", "phenotype_id"],
        right_on=["variant_id", "gene_id"],
        how="inner",
    )
    merged = merged.drop(columns=["gene_id"])

    if merged.empty:
        return None

    return merged


def parse_variant_pos_ref_alt(df: pd.DataFrame) -> pd.DataFrame:
    """Add _pos, _ref, _alt columns parsed from variant_id for MPAC joining."""
    parts = df["variant_id"].str.split("_", expand=True)
    df = df.copy()
    df["_pos"] = parts[1].astype(np.int64)
    df["_ref"] = parts[2]
    df["_alt"] = parts[3]
    return df


def main():
    args = parse_args()
    chrom = args.chrom
    os.makedirs(args.out_dir, exist_ok=True)

    print(f"[{chrom}] Loading MPAC predictions...", file=sys.stderr)
    mpac = load_mpac(args.mpac_dir, chrom)

    tissues = get_tissues(args.susie_dir)
    print(f"[{chrom}] Processing {len(tissues)} tissues...", file=sys.stderr)

    tissue_frames = []
    for tissue in tissues:
        df = load_tissue_data(tissue, chrom, args.susie_dir, args.signif_dir)
        if df is None:
            continue
        df["lead_tissue"] = tissue
        tissue_frames.append(df)
        print(f"  {tissue}: {len(df)} high-PIP significant variants", file=sys.stderr)

    if not tissue_frames:
        print(f"[{chrom}] No qualifying variants found.", file=sys.stderr)
        return

    all_variants = pd.concat(tissue_frames, ignore_index=True)

    # Parse pos/ref/alt for MPAC join; filter MPAC to only needed positions first
    all_variants = parse_variant_pos_ref_alt(all_variants)
    needed = all_variants[["_pos", "_ref", "_alt"]].drop_duplicates()
    mpac_filtered = mpac.merge(
        needed, left_on=["pos", "ref", "alt"], right_on=["_pos", "_ref", "_alt"], how="inner"
    ).drop(columns=["_pos", "_ref", "_alt"])

    # Join MPAC: inner join drops variants not in any dELS enhancer
    result = all_variants.merge(
        mpac_filtered,
        left_on=["_pos", "_ref", "_alt"],
        right_on=["pos", "ref", "alt"],
        how="inner",
    )
    result = result.drop(columns=["_pos", "_ref", "_alt"])

    # chrom, pos, ref, alt now come from MPAC (authoritative)
    result = result[OUTPUT_COLS]
    result = result.drop_duplicates()

    out_path = os.path.join(args.out_dir, f"{chrom}_highPIP_eQTL_mpac.tsv.gz")
    result.to_csv(out_path, sep="\t", index=False, compression="gzip")
    print(
        f"[{chrom}] Wrote {len(result)} rows to {out_path}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
