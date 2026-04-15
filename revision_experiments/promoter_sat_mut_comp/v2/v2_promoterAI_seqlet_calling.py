"""
v2 PromoterAI seqlet calling script
Converted from: promoterAI_promoter_filtering_v2.ipynb

Key changes from original:
  - Reads PromoterAI scores directly from v2 filtered VCF (not merged predictions TSV)
  - Assigns gene_id by intersecting positions with archive2 promoter BED
  - Writes all outputs to v2/
"""

import os
import numpy as np
import pandas as pd
import torch
import pickle
from tqdm import tqdm
from tangermeme.seqlet import recursive_seqlets
from tangermeme.annotate import annotate_seqlets
from tangermeme.io import read_meme

V2       = os.path.dirname(os.path.abspath(__file__))  # this script lives in v2/
PROJ     = os.path.dirname(V2)                          # project root
ARCHIVE2 = os.path.join(PROJ, "raw_data/archive2/gencode.v44.protein.coding.250bp.promoters.autosomes.v2.exon.splice.subtracted.bed")
PWM      = "/projects/tewhey-lab/buttsj/meme_suite_data/motif_databases/HOCOMOCO/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme"

# =============================================================================
# LOAD AND PARSE V2 PromoterAI VCF
# INFO format: PAI3D_percentile;{val}:PAI3D_thresh;{val}:PromoterAI_score;{val}
# =============================================================================

print("Loading v2 PromoterAI VCF...")
vcf = pd.read_csv(
    f"{V2}/processed_data/PrimateAI_and_PromoterAI_scores.hg38.250bp.v2.filtered.vcf",
    sep='\t', header=None,
    usecols=[0, 1, 3, 4, 7],
    names=['chrom', 'pos', 'ref', 'alt', 'info'],
    dtype={'chrom': 'category', 'ref': 'category', 'alt': 'category'}
)

print("Parsing PromoterAI scores from INFO column...")
vcf['promoterAI_score'] = [
    float(info.split('PromoterAI_score;')[-1])
    for info in tqdm(vcf['info'], desc="Parsing INFO")
]
vcf.drop(columns=['info'], inplace=True)

# =============================================================================
# ASSIGN GENE_ID FROM ARCHIVE2 BED
# Match each variant position to its promoter interval
# archive2 BED cols: 0=chrom, 1=start, 2=end, 3=gene_name, 4=score, 5=strand, 6=ENSG_ENST_gene
# =============================================================================

print("Loading archive2 promoter BED for gene assignment...")
bed = pd.read_csv(ARCHIVE2, sep='\t', header=None)
bed_intervals = pd.DataFrame({
    'chrom':    bed[0],
    'iv_start': bed[1],
    'iv_end':   bed[2],
    'gene_id':  bed[6],  # full ENSG_ENST_gene from col 6
})

# Assign gene_id via vectorized interval lookup (no cross-join)
print("Assigning gene_ids via interval overlap...")
chrom_dfs = []
for chrom, vcf_chrom in vcf.groupby('chrom', observed=True):
    bed_chrom = bed_intervals[bed_intervals['chrom'] == chrom].sort_values('iv_start').reset_index(drop=True)
    if bed_chrom.empty:
        continue
    starts   = bed_chrom['iv_start'].values
    ends     = bed_chrom['iv_end'].values
    gene_ids = bed_chrom['gene_id'].values
    positions = vcf_chrom['pos'].values
    # For each position find the last interval whose start <= position
    idx = np.searchsorted(starts, positions, side='right') - 1
    in_interval = (idx >= 0) & (positions <= ends[np.maximum(idx, 0)])
    vcf_hit = vcf_chrom.iloc[np.where(in_interval)[0]].copy()
    vcf_hit['gene_id'] = gene_ids[idx[in_interval]]
    chrom_dfs.append(vcf_hit)
merged = pd.concat(chrom_dfs, ignore_index=True)
merged = merged.drop_duplicates(subset=['chrom', 'pos', 'ref', 'alt'])

print(f"  {len(merged['gene_id'].unique())} unique promoters after gene assignment")

# Reformat for tensor generation (cast categoricals back to str)
promoterAI_4_tensors = pd.DataFrame({
    'chrom': merged['chrom'].astype(str),
    'pos':   merged['pos'],
    'id':    merged['gene_id'],
    'ref':   merged['ref'].astype(str),
    'alt':   merged['alt'].astype(str),
    'score': merged['promoterAI_score'],
})

# =============================================================================
# TENSOR GENERATION
# =============================================================================

def predictions_to_tensors(df, target_length=250):
    misfit_genes = []
    base_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    tensor_dict = {}
    for promoter_id, promoter_df in tqdm(df.groupby('id'), desc="Building tensors"):
        position_vectors = []
        for pos, pos_df in promoter_df.groupby('pos'):
            score_vector = [0.0] * 4
            for _, row in pos_df.iterrows():
                score_vector[base_map[row['alt'].upper()]] = row['score']
            position_vectors.append(score_vector)
        t = torch.tensor(position_vectors, dtype=torch.float32).T
        n = t.shape[-1]
        if n == target_length:
            tensor_dict[promoter_id] = t
        elif n < target_length:
            pad = torch.zeros((4, target_length - n), dtype=torch.float32)
            tensor_dict[promoter_id] = torch.cat([t, pad], dim=1)
        else:
            misfit_genes.append(promoter_id)
    return tensor_dict, misfit_genes

print("Generating PromoterAI tensors...")
promoterAI_raw_tensors, misfit_genes = predictions_to_tensors(promoterAI_4_tensors)
print(f"  {len(promoterAI_raw_tensors)} promoters | {len(misfit_genes)} misfits")

gene_name_list = list(promoterAI_raw_tensors.keys())
raw_stack = torch.stack([promoterAI_raw_tensors[i] for i in gene_name_list])

# =============================================================================
# SEQLET / ONESHOT / PLOTLOGO TENSORS
# =============================================================================

def tensors4tangermeme(stacked_raw, gene_list):
    seqlet_dict, plotLogo_dict, oneHot_dict = {}, {}, {}
    for idx, promoter in enumerate(gene_list):
        raw     = stacked_raw[idx]
        oneHot  = 1 * (raw == 0)
        seqlet  = (-1 * raw.sum(dim=0)) / 3
        plotLogo = oneHot * seqlet
        seqlet_dict[promoter]   = seqlet
        plotLogo_dict[promoter] = plotLogo
        oneHot_dict[promoter]   = oneHot
    return seqlet_dict, plotLogo_dict, oneHot_dict

print("Building seqlet/oneHot/plotLogo tensors...")
pAI_seqlets, pAI_plotLogos, pAI_oneHots = tensors4tangermeme(raw_stack, gene_name_list)

pAI_seqlet_stack  = torch.stack([pAI_seqlets[i]   for i in gene_name_list])
pAI_plotLogo_stack = torch.stack([pAI_plotLogos[i] for i in gene_name_list])
pAI_oneHot_stack  = torch.stack([pAI_oneHots[i]   for i in gene_name_list])

# =============================================================================
# SEQLET CALLING AND ANNOTATION
# =============================================================================

print("Calling PromoterAI seqlets (threshold=0.01)...")
promoterAI_seqlet_calls = recursive_seqlets(pAI_seqlet_stack, threshold=0.01)

print("Annotating seqlets...")
motifs = read_meme(PWM)
motif_names = list(motifs.keys())
motif_idxs, motif_pvals = annotate_seqlets(pAI_oneHot_stack, promoterAI_seqlet_calls, PWM)

promoter_idx_dict = {gene_name_list.index(p): p for p in gene_name_list}

full_ann_seqlets = promoterAI_seqlet_calls.copy()
full_ann_seqlets['promoter_id'] = [promoter_idx_dict.get(i) for i in full_ann_seqlets['example_idx']]
full_ann_seqlets['example_idx'] = [motif_names[i] for i in motif_idxs]
full_ann_seqlets['pval']        = [float(p[0]) for p in motif_pvals]

# Save
promoterAI_seqlet_calls.to_csv(
    f"{V2}/promoterAI/processed_data/promoterAI_seqlet_calls_v2.tsv", sep='\t', index=False)
full_ann_seqlets.to_csv(
    f"{V2}/promoterAI/processed_data/promoterAI_full_ann_seqlets_v2.tsv", sep='\t', index=False)
with open(f"{V2}/promoterAI/processed_data/promoterAI_plotLogo_dict_v2.pkl", 'wb') as f:
    pickle.dump(pAI_plotLogos, f)

print("Saved seqlet TSVs -> v2/promoterAI/processed_data/")

# =============================================================================
# SEQLETS TO BED
# =============================================================================

def seqlets2bed(full_preds, full_ann_seqlets):
    promoter_min_pos = full_preds.groupby('id')['pos'].min()
    promoter_chrom   = full_preds.groupby('id')['chrom'].min()
    seqlets = full_ann_seqlets.copy()
    seqlets['promoter_start'] = seqlets['promoter_id'].map(promoter_min_pos)
    seqlets['chrom']          = seqlets['promoter_id'].map(promoter_chrom)
    seqlets.dropna(subset=['promoter_start'], inplace=True)
    seqlets['promoter_start'] = seqlets['promoter_start'].astype(int)
    bed_start = seqlets['promoter_start'] + seqlets['start'] - 1
    bed_end   = seqlets['promoter_start'] + seqlets['end']
    bed_id    = [
        f"{tf}_{pid}_{attr}"
        for tf, pid, attr in zip(seqlets['example_idx'], seqlets['promoter_id'], seqlets['attribution'])
    ]
    return pd.DataFrame({
        'chrom': seqlets['chrom'],
        'start': bed_start,
        'end':   bed_end,
        'id':    bed_id,
    }).sort_values(['chrom', 'start'])

print("Converting PromoterAI seqlets to BED...")
promoterAI_seqlet_bed = seqlets2bed(promoterAI_4_tensors, full_ann_seqlets)
promoterAI_seqlet_bed.to_csv(
    f"{V2}/processed_data/bed_files/promoterAI_annotated_seqlets_01_v2.bed",
    sep='\t', index=False, header=None
)

print(f"  {len(promoterAI_seqlet_bed)} seqlets -> v2/processed_data/bed_files/promoterAI_annotated_seqlets_01_v2.bed")
print("Done.")
