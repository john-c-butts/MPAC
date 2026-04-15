"""
v2 MPAC seqlet calling script
Converted from: mpac_promoter_seqlet_calling_112125.ipynb

Key changes from original:
  - Reads from v2 tabix-filtered VCF (exon+splice subtracted promoters)
  - Reads promoter BED from archive2 (7-col format, ENSG_ENST_gene in col 6)
  - Writes all outputs to v2/
"""

import os
import pandas as pd
import torch
import pickle
from tqdm import tqdm
from tangermeme.seqlet import recursive_seqlets
from tangermeme.annotate import annotate_seqlets
from tangermeme.io import read_meme

V2   = os.path.dirname(os.path.abspath(__file__))  # this script lives in v2/
PROJ = os.path.dirname(V2)                          # project root
PWM  = "/projects/tewhey-lab/buttsj/meme_suite_data/motif_databases/HOCOMOCO/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme"

# =============================================================================
# LOAD PREDICTIONS
# =============================================================================

def vcf2df(pred_df):
    k_ref, k_alt, k_skew = [], [], []
    h_ref, h_alt, h_skew = [], [], []
    s_ref, s_alt, s_skew = [], [], []
    for i in tqdm(pred_df[7], desc="Parsing INFO"):
        all_preds = i.split(';')
        k_ref.append(float(all_preds[0].split('=')[-1]))
        k_alt.append(float(all_preds[3].split('=')[-1]))
        k_skew.append(float(all_preds[6].split('=')[-1]))
        h_ref.append(float(all_preds[1].split('=')[-1]))
        h_alt.append(float(all_preds[4].split('=')[-1]))
        h_skew.append(float(all_preds[7].split('=')[-1]))
        s_ref.append(float(all_preds[2].split('=')[-1]))
        s_alt.append(float(all_preds[5].split('=')[-1]))
        s_skew.append(float(all_preds[8].split('=')[-1]))
    return pd.DataFrame({
        'chrom': pred_df[0],
        'pos':   pred_df[1],
        'id':    pred_df[2],
        'ref':   pred_df[3],
        'alt':   pred_df[4],
        'k562_ref_pred':   k_ref, 'k562_alt_pred':   k_alt, 'k562_skew_pred':   k_skew,
        'hepg2_ref_pred':  h_ref, 'hepg2_alt_pred':  h_alt, 'hepg2_skew_pred':  h_skew,
        'sknsh_ref_pred':  s_ref, 'sknsh_alt_pred':  s_alt, 'sknsh_skew_pred':  s_skew,
        'gene_id':      [i.split('..')[0] for i in pred_df[2]],
    })

print("Loading v2 MPAC predictions...")
all_tabix_mpac_preds = vcf2df(
    pd.read_csv(f"{V2}/processed_data/all.mpac.preds.tabix.filtered.gencode.250bp.v2.vcf",
                sep='\t', header=None)
)
print(f"  {len(all_tabix_mpac_preds['gene_id'].unique())} unique promoters after tabix filter")

# Filter to only ENSG IDs present in the archive2 exon+splice subtracted BED
archive2_bed = pd.read_csv(
    f"{PROJ}/raw_data/archive2/gencode.v44.protein.coding.250bp.promoters.autosomes.v2.exon.splice.subtracted.bed",
    sep='\t', header=None
)
archive2_gene_ids = set(archive2_bed[6])
final_mpac_preds = all_tabix_mpac_preds[
    all_tabix_mpac_preds['gene_id'].isin(archive2_gene_ids)
].copy()
print(f"  {len(final_mpac_preds['gene_id'].unique())} unique promoters after archive2 gene ID filter")

# Add merge_id and activity annotations
final_mpac_preds = final_mpac_preds.copy()
final_mpac_preds['merge_id'] = [
    ':'.join([c, str(p), r, a])
    for c, p, r, a in zip(
        final_mpac_preds['chrom'], final_mpac_preds['pos'],
        final_mpac_preds['ref'],   final_mpac_preds['alt']
    )
]

mean_activity = final_mpac_preds.filter(
    ['k562_ref_pred', 'hepg2_ref_pred', 'sknsh_ref_pred', 'gene_id']
).groupby('gene_id').mean()

for cell, col in [('k562', 'k562_ref_pred'), ('hepg2', 'hepg2_ref_pred'), ('sknsh', 'sknsh_ref_pred')]:
    act_dict = dict(zip(mean_activity.index, mean_activity[col]))
    final_mpac_preds[f'{cell}_mean_ref']    = [act_dict.get(i) for i in final_mpac_preds['gene_id']]
    final_mpac_preds[f'{cell}_mean_active'] = [1 if v >= 1 else 0 for v in final_mpac_preds[f'{cell}_mean_ref']]
    final_mpac_preds[f'{cell}_emVar']       = [1 if abs(v) > 0.5 else 0 for v in final_mpac_preds[f'{cell}_skew_pred']]

final_mpac_preds['mean_active_any'] = [
    1 if sum([k, h, s]) > 0 else 0
    for k, h, s in zip(final_mpac_preds['k562_mean_active'],
                       final_mpac_preds['hepg2_mean_active'],
                       final_mpac_preds['sknsh_mean_active'])
]
final_mpac_preds['emVar_any'] = [
    1 if sum([k, h, s]) > 0 else 0
    for k, h, s in zip(final_mpac_preds['k562_emVar'],
                       final_mpac_preds['hepg2_emVar'],
                       final_mpac_preds['sknsh_emVar'])
]
final_mpac_preds['emVar_all'] = [
    1 if sum([k, h, s]) == 3 else 0
    for k, h, s in zip(final_mpac_preds['k562_emVar'],
                       final_mpac_preds['hepg2_emVar'],
                       final_mpac_preds['sknsh_emVar'])
]

# Add RNA expression
print("Adding RNA expression...")
rna_expression = pd.read_csv(f"{PROJ}/raw_data/rna_celline_filtered.tsv", sep='\t')
for cell, cell_label, col in [
    ('k562',  'K-562',   'k562_rna_tpm'),
    ('hepg2', 'Hep-G2',  'hepg2_rna_tpm'),
    ('sknsh', 'SK-N-SH', 'sknsh_rna_tpm'),
]:
    rna_dict = dict(zip(
        rna_expression[rna_expression['Cell line'] == cell_label]['Gene'],
        rna_expression[rna_expression['Cell line'] == cell_label]['TPM']
    ))
    final_mpac_preds[col] = [rna_dict.get(i.split('.')[0]) for i in final_mpac_preds['gene_id']]

final_mpac_preds.to_csv(
    f"{V2}/processed_data/mpac.250bp.promoter.preds.with.annotations.v2.tsv",
    sep='\t', index=False
)
print(f"  Saved annotations -> v2/processed_data/mpac.250bp.promoter.preds.with.annotations.v2.tsv")

# =============================================================================
# TENSOR GENERATION
# =============================================================================

def predictions_to_tensors_multi_cell(df, target_length=250):
    misfit_genes = []
    cell_types = sorted(set(col.split('_')[0] for col in df.columns if '_skew_pred' in col))
    print(f"  Cell types detected: {', '.join(cell_types)}")
    base_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    tensor_dict = {cell: {} for cell in cell_types}
    for gene_id, gene_df in tqdm(df.groupby('gene_id'), desc="Building tensors"):
        position_vectors_by_cell = {cell: [] for cell in cell_types}
        for pos, pos_df in gene_df.groupby('pos'):
            score_vectors = {cell: [0.0] * 4 for cell in cell_types}
            for _, row in pos_df.iterrows():
                idx = base_map[row['alt'].upper()]
                for cell in cell_types:
                    score_vectors[cell][idx] = row[f'{cell}_skew_pred']
            for cell in cell_types:
                position_vectors_by_cell[cell].append(score_vectors[cell])
        for cell in cell_types:
            t = torch.tensor(position_vectors_by_cell[cell], dtype=torch.float32).T
            n = t.shape[-1]
            if n == target_length:
                tensor_dict[cell][gene_id] = t
            elif n < target_length:
                pad = torch.zeros((4, target_length - n), dtype=torch.float32)
                tensor_dict[cell][gene_id] = torch.cat([t, pad], dim=1)
            else:
                misfit_genes.append(gene_id)
    return tensor_dict, misfit_genes

print("Generating tensors...")
mpac_rawTensors, misfit_ids = predictions_to_tensors_multi_cell(final_mpac_preds)
print(f"  Misfit genes: {len(misfit_ids)}")

promoter_list = list(mpac_rawTensors['k562'].keys())
k_raw = torch.stack([mpac_rawTensors['k562'][i]  for i in promoter_list])
h_raw = torch.stack([mpac_rawTensors['hepg2'][i] for i in promoter_list])
s_raw = torch.stack([mpac_rawTensors['sknsh'][i] for i in promoter_list])

def tensors4tangermeme(stacked_raw, gene_list):
    seqlet_dict, plotLogo_dict, oneHot_dict = {}, {}, {}
    for idx, promoter in enumerate(gene_list):
        raw = stacked_raw[idx]
        oneHot  = 1 * (raw == 0)
        seqlet  = (-1 * raw.sum(dim=0)) / 3
        plotLogo = oneHot * seqlet
        seqlet_dict[promoter]  = seqlet
        plotLogo_dict[promoter] = plotLogo
        oneHot_dict[promoter]  = oneHot
    return seqlet_dict, plotLogo_dict, oneHot_dict

def stack_tensors(seqlet_dict, plotLogo_dict, oneHot_dict, promoters):
    return (
        torch.stack([seqlet_dict[p]  for p in promoters]),
        torch.stack([plotLogo_dict[p] for p in promoters]),
        torch.stack([oneHot_dict[p]  for p in promoters]),
    )

print("Building seqlet/oneHot/plotLogo tensors...")
k_seqlets, k_plotLogos, k_oneHots = tensors4tangermeme(k_raw, promoter_list)
h_seqlets, h_plotLogos, h_oneHots = tensors4tangermeme(h_raw, promoter_list)
s_seqlets, s_plotLogos, s_oneHots = tensors4tangermeme(s_raw, promoter_list)

k_stackedSeqlets, k_stackedPlotLogos, k_stackedOneHots = stack_tensors(k_seqlets, k_plotLogos, k_oneHots, promoter_list)
h_stackedSeqlets, h_stackedPlotLogos, h_stackedOneHots = stack_tensors(h_seqlets, h_plotLogos, h_oneHots, promoter_list)
s_stackedSeqlets, s_stackedPlotLogos, s_stackedOneHots = stack_tensors(s_seqlets, s_plotLogos, s_oneHots, promoter_list)

# =============================================================================
# SEQLET CALLING
# =============================================================================

print("Calling seqlets (threshold=0.01)...")
k_seqlet_calls = recursive_seqlets(k_stackedSeqlets, threshold=0.01)
h_seqlet_calls = recursive_seqlets(h_stackedSeqlets, threshold=0.01)
s_seqlet_calls = recursive_seqlets(s_stackedSeqlets, threshold=0.01)

# =============================================================================
# SEQLET ANNOTATION
# =============================================================================

def ann_seqlets(motif_path, stacked_oneHots, seqlet_calls, promoters):
    motifs = read_meme(motif_path)
    motif_names = list(motifs.keys())
    motif_idxs, motif_pvals = annotate_seqlets(stacked_oneHots, seqlet_calls, motif_path)
    promoter_idx_dict = {promoters.index(p): p for p in promoters}
    full_ann = seqlet_calls.copy()
    full_ann['promoter_id'] = [promoter_idx_dict.get(i) for i in full_ann['example_idx']]
    full_ann['example_idx'] = [motif_names[i] for i in motif_idxs]
    full_ann['pval']        = [float(p[0]) for p in motif_pvals]
    return full_ann

print("Annotating seqlets...")
k_full_ann = ann_seqlets(PWM, k_stackedOneHots, k_seqlet_calls, promoter_list)
h_full_ann = ann_seqlets(PWM, h_stackedOneHots, h_seqlet_calls, promoter_list)
s_full_ann = ann_seqlets(PWM, s_stackedOneHots, s_seqlet_calls, promoter_list)

# Save annotated seqlet TSVs and plotLogo dicts to v2/mpac/
k_seqlet_calls.to_csv(f"{V2}/mpac/processed_data/mpac_k562_seqlet_calls_v2.tsv",    sep='\t', index=False)
k_full_ann.to_csv(    f"{V2}/mpac/processed_data/mpac_k562_full_ann_seqlets_v2.tsv", sep='\t', index=False)
h_seqlet_calls.to_csv(f"{V2}/mpac/processed_data/mpac_hepg2_seqlet_calls_v2.tsv",   sep='\t', index=False)
h_full_ann.to_csv(    f"{V2}/mpac/processed_data/mpac_hepg2_full_ann_seqlets_v2.tsv",sep='\t', index=False)
s_seqlet_calls.to_csv(f"{V2}/mpac/processed_data/mpac_sknsh_seqlet_calls_v2.tsv",   sep='\t', index=False)
s_full_ann.to_csv(    f"{V2}/mpac/processed_data/mpac_sknsh_full_ann_seqlets_v2.tsv",sep='\t', index=False)

for cell, d in [('k562', k_seqlets), ('hepg2', h_seqlets), ('sknsh', s_seqlets)]:
    with open(f"{V2}/mpac/processed_data/mpac_{cell}_plotLogo_dict_v2.pkl", 'wb') as f:
        pickle.dump(d, f)

print("Saved seqlet TSVs and plotLogo dicts -> v2/mpac/processed_data/")

# =============================================================================
# SEQLETS TO BED
# =============================================================================

def seqlets2bed(full_preds, full_ann_seqlets):
    promoter_min_pos = full_preds.groupby('gene_id')['pos'].min()
    promoter_chrom   = full_preds.groupby('gene_id')['chrom'].min()
    seqlets = full_ann_seqlets.copy()
    seqlets['promoter_start'] = seqlets['promoter_id'].map(promoter_min_pos)
    seqlets['chrom']          = seqlets['promoter_id'].map(promoter_chrom)
    seqlets.dropna(subset=['promoter_start'], inplace=True)
    seqlets['promoter_start'] = seqlets['promoter_start'].astype(int)
    bed_start = seqlets['promoter_start'] + seqlets['start'] - 1
    bed_end   = seqlets['promoter_start'] + seqlets['end']
    bed_id    = [
        f"{tf}_{pid}_{attr}"
        for tf, pid, attr in zip(
            seqlets['example_idx'], seqlets['promoter_id'], seqlets['attribution']
        )
    ]
    return pd.DataFrame({
        'chrom': seqlets['chrom'],
        'start': bed_start,
        'end':   bed_end,
        'id':    bed_id,
    }).sort_values(['chrom', 'start'])

print("Converting seqlets to BED...")
k_bed = seqlets2bed(final_mpac_preds, k_full_ann)
h_bed = seqlets2bed(final_mpac_preds, h_full_ann)
s_bed = seqlets2bed(final_mpac_preds, s_full_ann)

k_bed.to_csv(f"{V2}/processed_data/bed_files/k562_mpac_annotated_seqlets_01_v2.bed",  sep='\t', index=False, header=None)
h_bed.to_csv(f"{V2}/processed_data/bed_files/hepg2_mpac_annotated_seqlets_01_v2.bed", sep='\t', index=False, header=None)
s_bed.to_csv(f"{V2}/processed_data/bed_files/sknsh_mpac_annotated_seqlets_01_v2.bed", sep='\t', index=False, header=None)

print("Done. BED files -> v2/processed_data/bed_files/")
print(f"  K562:  {len(k_bed)} seqlets")
print(f"  HepG2: {len(h_bed)} seqlets")
print(f"  SKNSH: {len(s_bed)} seqlets")
