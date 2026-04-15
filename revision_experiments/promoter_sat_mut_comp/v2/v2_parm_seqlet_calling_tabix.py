"""
v2 PARM seqlet calling - tabix-filter approach
Alternative to v2_parm_seqlet_calling.py for comparison.

Instead of filtering PARM folders by the promoter BED upfront, this script:
  1. Processes ALL PARM prediction folders -> one concatenated VCF per cell type
  2. bgzips + tabix indexes the VCF
  3. tabix filters with the exon+splice subtracted BED (same as MPAC/PromoterAI)
  4. Builds tensors and calls seqlets from the filtered predictions

Outputs use _v2_tabix suffix to distinguish from the folder-filter approach.
"""

import glob
import os
import subprocess
import torch
import pandas as pd
import pickle
from tqdm import tqdm
from tangermeme.seqlet import recursive_seqlets
from tangermeme.annotate import annotate_seqlets
from tangermeme.io import read_meme

V2        = os.path.dirname(os.path.abspath(__file__))
PROJ      = os.path.dirname(V2)
PARM_ROOT = "/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/PARM/preds_out"
ARCHIVE2_SUB = os.path.join(PROJ, "raw_data/archive2/gencode.v44.protein.coding.250bp.promoters.autosomes.v2.exon.splice.subtracted.bed")
PWM       = "/projects/tewhey-lab/buttsj/meme_suite_data/motif_databases/HOCOMOCO/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme"

# =============================================================================
# STEP 1: PARSE ALL PARM FOLDERS INTO A LONG-FORMAT DATAFRAME
# =============================================================================

def PARM2DF(path2folders, cell_type):
    dfs = []
    for folder in tqdm(path2folders, desc=f"Parsing {cell_type} PARM folders"):
        interval = folder.split('/')[-1].split('::')[-1]
        gene_id  = folder.split('/')[-1].split('::')[0]
        chrom    = interval.split(':')[0]
        pos      = int(interval.split('-')[0].split(':')[-1]) + 1
        parm_files = sorted(os.listdir(folder))
        mut_file = pd.read_csv(f"{folder}/{parm_files[-1]}", sep='\t')
        for base in ['A', 'C', 'G', 'T']:
            mut_file.loc[mut_file['Ref'] == base, base] = 0
        mut_file['chrom'] = chrom
        pos_list = []
        for i in range(len(mut_file)):
            pos_list.append(pos)
            pos += 1
        mut_file['pos'] = pos_list
        long = mut_file.melt(
            id_vars=['chrom', 'pos', 'Ref'],
            value_vars=['A', 'C', 'G', 'T'],
            var_name='alt', value_name='skew'
        ).rename(columns={'Ref': 'ref'})
        long.sort_values(['pos', 'alt'], inplace=True)
        final = long[long['ref'] != long['alt']].copy()
        final['gene_id']  = gene_id
        final['merge_id'] = [
            ':'.join([c, str(p), r, a])
            for c, p, r, a in zip(final['chrom'], final['pos'], final['ref'], final['alt'])
        ]
        final.reset_index(drop=True, inplace=True)
        dfs.append(final)
    return pd.concat(dfs)

k562_folders  = sorted(glob.glob(f"{PARM_ROOT}/all_gencode_v44_promoter_sat_mut_parm_preds_k562/*"))
hepg2_folders = sorted(glob.glob(f"{PARM_ROOT}/all_gencode_v44_promoter_sat_mut_parm_preds_hepg2/*"))

print(f"Processing {len(k562_folders)} K562 and {len(hepg2_folders)} HepG2 PARM folders...")
k562_parm_preds_full  = PARM2DF(k562_folders,  'K562')
hepg2_parm_preds_full = PARM2DF(hepg2_folders, 'HepG2')

# =============================================================================
# STEP 2: WRITE VCF, BGZIP, TABIX INDEX
# =============================================================================

def write_and_index_vcf(preds, out_vcf_gz):
    vcf = pd.DataFrame({
        '#CHROM': preds['chrom'],
        'POS':    preds['pos'],
        'ID':     preds['gene_id'],
        'REF':    preds['ref'],
        'ALT':    preds['alt'],
        'QUAL':   '.',
        'FILTER': '.',
        'INFO':   preds['skew'],
    }).sort_values(['#CHROM', 'POS'])

    tmp_vcf = out_vcf_gz.replace('.gz', '')
    vcf.to_csv(tmp_vcf, sep='\t', index=False)
    subprocess.run(['bgzip', '-f', tmp_vcf], check=True)
    subprocess.run(['tabix', '-s1', '-b2', '-e2', '-S1', out_vcf_gz], check=True)
    print(f"  Written and indexed: {out_vcf_gz}")

k562_vcf_gz  = f"{V2}/processed_data/all_k562_parm_preds_full_v2_tabix.vcf.gz"
hepg2_vcf_gz = f"{V2}/processed_data/all_hepg2_parm_preds_full_v2_tabix.vcf.gz"

print("Writing and indexing PARM VCFs...")
write_and_index_vcf(k562_parm_preds_full,  k562_vcf_gz)
write_and_index_vcf(hepg2_parm_preds_full, hepg2_vcf_gz)

# =============================================================================
# STEP 3: TABIX FILTER WITH EXON+SPLICE SUBTRACTED BED
# =============================================================================

def tabix_filter(vcf_gz, bed, out_vcf):
    result = subprocess.run(
        ['tabix', '-R', bed, vcf_gz],
        capture_output=True, text=True, check=True
    )
    # Prepend header from the vcf
    header = subprocess.run(
        ['bgzip', '-cd', vcf_gz],
        capture_output=True, text=True
    ).stdout.split('\n')[0]
    with open(out_vcf, 'w') as f:
        f.write(header + '\n')
        f.write(result.stdout)
    print(f"  Filtered: {out_vcf}")

k562_filtered  = f"{V2}/processed_data/all_k562_250bp_parm_preds_v2_tabix.vcf"
hepg2_filtered = f"{V2}/processed_data/all_hepg2_250bp_parm_preds_v2_tabix.vcf"

print("Tabix filtering PARM predictions...")
tabix_filter(k562_vcf_gz,  ARCHIVE2_SUB, k562_filtered)
tabix_filter(hepg2_vcf_gz, ARCHIVE2_SUB, hepg2_filtered)

# =============================================================================
# STEP 4: READ FILTERED VCFS BACK INTO DATAFRAMES
# =============================================================================

def vcf_to_parm_df(vcf_path):
    df = pd.read_csv(vcf_path, sep='\t',
                     names=['chrom', 'pos', 'gene_id', 'ref', 'alt', 'qual', 'filter', 'skew'],
                     skiprows=1)
    df['merge_id'] = [
        ':'.join([c, str(p), r, a])
        for c, p, r, a in zip(df['chrom'], df['pos'], df['ref'], df['alt'])
    ]
    return df

print("Reading filtered predictions...")
k562_parm_preds  = vcf_to_parm_df(k562_filtered)
hepg2_parm_preds = vcf_to_parm_df(hepg2_filtered)
print(f"  K562:  {len(k562_parm_preds['gene_id'].unique())} promoters after tabix filter")
print(f"  HepG2: {len(hepg2_parm_preds['gene_id'].unique())} promoters after tabix filter")

# Filter to only ENSG IDs present in the archive2 exon+splice subtracted BED
archive2_bed = pd.read_csv(ARCHIVE2_SUB, sep='\t', header=None)
archive2_gene_ids = set(archive2_bed[6])
k562_parm_preds  = k562_parm_preds[k562_parm_preds['gene_id'].isin(archive2_gene_ids)].copy()
hepg2_parm_preds = hepg2_parm_preds[hepg2_parm_preds['gene_id'].isin(archive2_gene_ids)].copy()
print(f"  K562:  {len(k562_parm_preds['gene_id'].unique())} promoters after archive2 gene ID filter")
print(f"  HepG2: {len(hepg2_parm_preds['gene_id'].unique())} promoters after archive2 gene ID filter")

# =============================================================================
# STEP 5: BUILD TENSORS FROM FILTERED PREDICTIONS
# =============================================================================

def parm_df_to_tensors(df, target_length=250):
    base_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    raw_dict, seqlet_dict, plotLogo_dict, oneHot_dict = {}, {}, {}, {}
    misfit = []
    for gene_id, gene_df in tqdm(df.groupby('gene_id'), desc="Building tensors"):
        position_vectors = []
        for pos, pos_df in gene_df.groupby('pos'):
            vec = [0.0] * 4
            for _, row in pos_df.iterrows():
                vec[base_map[row['alt'].upper()]] = row['skew']
            position_vectors.append(vec)
        raw = torch.tensor(position_vectors, dtype=torch.float32).T
        n = raw.shape[-1]
        if n == target_length:
            pass
        elif n < target_length:
            pad = torch.zeros((4, target_length - n), dtype=torch.float32)
            raw = torch.cat([raw, pad], dim=1)
        else:
            misfit.append(gene_id)
            continue
        seqlet  = raw.sum(dim=0) / 3
        oneHot  = 1 * (raw == 0)
        plotLogo = seqlet * oneHot
        raw_dict[gene_id]     = raw
        seqlet_dict[gene_id]  = seqlet
        plotLogo_dict[gene_id] = plotLogo
        oneHot_dict[gene_id]  = oneHot
    return raw_dict, seqlet_dict, plotLogo_dict, oneHot_dict, misfit

print("Building K562 tensors...")
k562_raw, k562_seqlet, k562_plotLogo, k562_oneHot, k562_misfit = parm_df_to_tensors(k562_parm_preds)
print(f"  {len(k562_raw)} promoters | {len(k562_misfit)} misfits")

print("Building HepG2 tensors...")
hepg2_raw, hepg2_seqlet, hepg2_plotLogo, hepg2_oneHot, hepg2_misfit = parm_df_to_tensors(hepg2_parm_preds)
print(f"  {len(hepg2_raw)} promoters | {len(hepg2_misfit)} misfits")

# Use intersection of promoters present in both cell types for consistent stacking
promoter_list = sorted(set(k562_raw.keys()) & set(hepg2_raw.keys()))
print(f"  {len(promoter_list)} promoters in both cell types")

# =============================================================================
# STEP 6: STACK, CALL, ANNOTATE SEQLETS
# =============================================================================

def stack(d, promoters):
    return torch.stack([d[p] for p in promoters])

k_stackedSeqlets  = stack(k562_seqlet,  promoter_list)
k_stackedOneHots  = stack(k562_oneHot,  promoter_list)
h_stackedSeqlets  = stack(hepg2_seqlet, promoter_list)
h_stackedOneHots  = stack(hepg2_oneHot, promoter_list)

print("Calling seqlets (threshold=0.01)...")
k_seqlet_calls = recursive_seqlets(k_stackedSeqlets, threshold=0.01)
h_seqlet_calls = recursive_seqlets(h_stackedSeqlets, threshold=0.01)

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

print(f"  K562  seqlets: {len(k_full_ann)}")
print(f"  HepG2 seqlets: {len(h_full_ann)}")

# Save TSVs
k_seqlet_calls.to_csv(f"{V2}/PARM/processed_data/parm_k562_tangermeme_seqlet_calls_v2_tabix.tsv",    sep='\t', index=False)
k_full_ann.to_csv(    f"{V2}/PARM/processed_data/parm_k562_tangermeme_full_ann_seqlets_v2_tabix.tsv", sep='\t', index=False)
h_seqlet_calls.to_csv(f"{V2}/PARM/processed_data/parm_hepg2_tangermeme_seqlet_calls_v2_tabix.tsv",   sep='\t', index=False)
h_full_ann.to_csv(    f"{V2}/PARM/processed_data/parm_hepg2_tangermeme_full_ann_seqlets_v2_tabix.tsv",sep='\t', index=False)

# =============================================================================
# STEP 7: SEQLETS TO BED
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
        f"{tf}:{rho}:{pid}"
        for tf, rho, pid in zip(seqlets['example_idx'], seqlets['attribution'], seqlets['promoter_id'])
    ]
    return pd.DataFrame({
        'chrom': seqlets['chrom'],
        'start': bed_start,
        'end':   bed_end,
        'id':    bed_id,
    }).sort_values(['chrom', 'start'])

print("Converting seqlets to BED...")
k_bed = seqlets2bed(k562_parm_preds, k_full_ann)
h_bed = seqlets2bed(hepg2_parm_preds, h_full_ann)

k_bed.to_csv(f"{V2}/processed_data/bed_files/k562_parm_annotated_seqlets_01_v2_tabix.bed",  sep='\t', index=False, header=None)
h_bed.to_csv(f"{V2}/processed_data/bed_files/hepg2_parm_annotated_seqlets_01_v2_tabix.bed", sep='\t', index=False, header=None)

print(f"  K562  seqlet BED: {len(k_bed)}")
print(f"  HepG2 seqlet BED: {len(h_bed)}")
print("Done. All outputs use _v2_tabix suffix in v2/PARM/processed_data/ and v2/processed_data/bed_files/")
