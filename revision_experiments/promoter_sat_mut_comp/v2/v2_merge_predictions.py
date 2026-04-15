"""
v2 Merge Predictions
Consolidates: add_parm_preds2promoter_data.ipynb + add_promoterAI_preds2promtoer_data.ipynb

Starting from MPAC per-variant predictions (with annotations), adds:
  1. PARM baseline activity predictions (per-gene)
  2. PromoterAI variant effect scores (per-variant, from tabix-filtered VCF)
  3. PARM saturation mutagenesis predictions (per-variant, from tabix-filtered VCFs)
     - 95th-percentile emVar thresholds computed from v2 data

Outputs:
  v2/processed_data/mpac.250bp.promoter.preds.with.annotations.parm.act.promoterAI.parm.sat.mut.v2.tsv
  v2/processed_data/mpac.250bp.promoter.preds.with.annotations.parm.act.promoterAI.parm.sat.mut4vcfAnnotateFromBigWig.v2.vcf
"""

import os
import numpy as np
import pandas as pd
from tqdm import tqdm

V2   = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(V2)

# =============================================================================
# STEP 1: LOAD MPAC BASE PREDICTIONS
# =============================================================================

print("Loading MPAC predictions...")
df = pd.read_csv(f"{V2}/processed_data/mpac.250bp.promoter.preds.with.annotations.v2.tsv", sep='\t')
print(f"  {len(df['gene_id'].unique())} promoters, {len(df)} variants")

# =============================================================================
# STEP 2: PROMOTERAI SCORES (per-variant, from tabix-filtered VCF)
# =============================================================================

print("Loading PromoterAI scores...")
pai_vcf = pd.read_csv(
    f"{V2}/processed_data/PrimateAI_and_PromoterAI_scores.hg38.250bp.v2.filtered.vcf",
    sep='\t', header=None,
    usecols=[0, 1, 3, 4, 7],
    names=['chrom', 'pos', 'ref', 'alt', 'info'],
)
pai_vcf['promoterAI_score'] = [
    float(i.split('PromoterAI_score;')[-1])
    for i in tqdm(pai_vcf['info'], desc="  Parsing INFO")
]
pai_vcf['merge_id'] = [
    ':'.join([c, str(p), r, a])
    for c, p, r, a in zip(pai_vcf['chrom'], pai_vcf['pos'], pai_vcf['ref'], pai_vcf['alt'])
]
pai_dict = dict(zip(pai_vcf['merge_id'], pai_vcf['promoterAI_score']))
df['promoterAI_score'] = [pai_dict.get(i) for i in df['merge_id']]
n_missing = df['promoterAI_score'].isna().sum()
print(f"  {len(df) - n_missing} matched, {n_missing} missing")

# =============================================================================
# STEP 3: PARM SAT MUT PREDICTIONS (per-variant, from tabix-filtered VCFs)
# =============================================================================

print("Loading PARM sat mut predictions...")
for cell, col in [('k562', 'parm_k562_pred'), ('hepg2', 'parm_hepg2_pred')]:
    parm_vcf = pd.read_csv(
        f"{V2}/processed_data/all_{cell}_250bp_parm_preds_v2_tabix.vcf",
        sep='\t', header=None,
        names=['chrom', 'pos', 'gene_id', 'ref', 'alt', 'qual', 'filter', 'skew'],
        skiprows=1
    )
    parm_vcf['merge_id'] = [
        ':'.join([c, str(p), r, a])
        for c, p, r, a in zip(parm_vcf['chrom'], parm_vcf['pos'], parm_vcf['ref'], parm_vcf['alt'])
    ]
    # PARM outputs Ref - Alt, so multiply by -1 for Alt - Ref convention
    parm_dict = dict(zip(parm_vcf['merge_id'], parm_vcf['skew']))
    df[col] = [-1 * parm_dict.get(i, np.nan) for i in df['merge_id']]
    print(f"  {cell}: {df[col].notna().sum()} variants matched")

df['parm_mean_pred'] = df[['parm_k562_pred', 'parm_hepg2_pred']].mean(axis=1)

# =============================================================================
# STEP 4: PARM emVar ANNOTATION (95th percentile threshold)
# =============================================================================

print("Computing PARM emVar thresholds...")
all_parm_skews = pd.concat([
    df['parm_k562_pred'].abs(),
    df['parm_hepg2_pred'].abs()
]).dropna()
parm_95 = np.quantile(all_parm_skews, 0.95)
print(f"  95th percentile threshold: {parm_95:.5f}")

df['parm_k562_emVar']  = (df['parm_k562_pred'].abs()  > parm_95).astype(int)
df['parm_hepg2_emVar'] = (df['parm_hepg2_pred'].abs() > parm_95).astype(int)
df['parm_emVar_any']   = ((df['parm_k562_emVar'] + df['parm_hepg2_emVar']) > 0).astype(int)

# =============================================================================
# STEP 5: SAVE MERGED TSV
# =============================================================================

out_tsv = f"{V2}/processed_data/mpac.250bp.promoter.preds.with.annotations.parm.act.promoterAI.parm.sat.mut.v2.tsv"
df.to_csv(out_tsv, sep='\t', index=False)
print(f"Saved merged TSV -> {os.path.basename(out_tsv)}")
print(f"  {len(df['gene_id'].unique())} promoters, {len(df)} variants, {len(df.columns)} columns")

# =============================================================================
# STEP 6: GENERATE VCF FOR VcfAnnotateFromBigWig
# Columns required: CHROM POS ID REF ALT QUAL FILTER INFO
# INFO encodes all numeric annotation columns
# =============================================================================

print("Generating VCF for VcfAnnotateFromBigWig...")
info_cols = [
    'k562_ref_pred', 'k562_alt_pred', 'k562_skew_pred',
    'hepg2_ref_pred', 'hepg2_alt_pred', 'hepg2_skew_pred',
    'sknsh_ref_pred', 'sknsh_alt_pred', 'sknsh_skew_pred',
    'parm_k562_pred', 'parm_hepg2_pred', 'parm_mean_pred',
    'promoterAI_score',
]
info_str = df[info_cols].apply(
    lambda row: ';'.join(f"{c}={v}" for c, v in zip(info_cols, row)), axis=1
)
vcf_df = pd.DataFrame({
    '#CHROM': df['chrom'],
    'POS':    df['pos'],
    'ID':     df['merge_id'],
    'REF':    df['ref'],
    'ALT':    df['alt'],
    'QUAL':   '.',
    'FILTER': '.',
    'INFO':   info_str,
}).sort_values(['#CHROM', 'POS'])

out_vcf = f"{V2}/processed_data/mpac.250bp.promoter.preds.with.annotations.parm.act.promoterAI.parm.sat.mut4vcfAnnotateFromBigWig.v2.vcf"
vcf_df.to_csv(out_vcf, sep='\t', index=False)
print(f"Saved VCF -> {os.path.basename(out_vcf)}")
print("Done.")
