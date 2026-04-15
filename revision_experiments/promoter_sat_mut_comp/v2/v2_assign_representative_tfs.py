"""
v2 Assign Representative TFs
Converted from: assign_representative_TFs_all_merged.ipynb

For each bedops-merged seqlet interval, assigns a single representative TF:
  1. Single-TF intervals: use that TF
  2. Multi-TF, same TF name: use that TF
  3. Multi-TF, different names: use the TF with highest absolute attribution

Inputs (from steps 3-5 of run_v2_pipeline.sh):
  - {K562,HepG2,SKNSH}_bedOps_merged_noMin_seqlets_01_v2.bed  (MPAC)
  - parm_{k562,hepg2}_tangermeme_bedOps_merged_noMin_seqlets_01_v2.bed  (PARM)
  - promoterAI_bedOps_merged_noMin_seqlets_01_v2.bed  (PromoterAI)

Outputs (to v2/processed_data/bed_files/):
  - mpac_{k562,hepg2,sknsh}_merged_collapsed_repTFs_v2.bed
  - parm_tangermeme_{k562,hepg2}_merged_collapsed_repTFs_v2.bed
  - promoterAI_tangermeme_merge_collapsed_repTFs_v2.bed
"""

import os
import pandas as pd

V2   = os.path.dirname(os.path.abspath(__file__))
BEDS = os.path.join(V2, "processed_data/bed_files")

# =============================================================================
# REPRESENTATIVE TF ASSIGNMENT
# ID format for MPAC/PromoterAI: {TF}_{TF}_{attribution}  (underscore-delimited)
# ID format for PARM:            {TF}-...:{attribution}    (colon-delimited)
# =============================================================================

def assign_representative_tf(df):
    """MPAC / PromoterAI: id format is TF_SUBFAMILY_attribution"""
    rep_tfs, activity_classes, attributions = [], [], []
    for _, row in df.iterrows():
        hits = str(row[3]).split(';')
        best_tf, best_score = '', 0.0
        if len(hits) == 1:
            parts = hits[0].split('_')
            best_tf = '_'.join(parts[:2])
            try: best_score = float(parts[-1])
            except (ValueError, IndexError): best_score = 0.0
        else:
            max_abs, all_names = -1.0, set()
            for hit in hits:
                parts = hit.split('_')
                tf_name = '_'.join(parts[:2])
                all_names.add(tf_name)
                try:
                    score = float(parts[-1])
                    if abs(score) > max_abs:
                        max_abs, best_tf, best_score = abs(score), tf_name, score
                except (ValueError, IndexError):
                    continue
            if len(all_names) == 1:
                best_tf = list(all_names)[0]
        rep_tfs.append(best_tf)
        attributions.append(best_score)
        activity_classes.append('Activator' if best_score > 0 else ('Repressor' if best_score < 0 else 'Neutral'))
    df = df.copy()
    df['representative_tf'] = rep_tfs
    df['activity_class']    = activity_classes
    df['attribution']       = attributions
    return df


def assign_representative_tf_parm(df):
    """PARM tangermeme: id format is TF-...:attribution"""
    rep_tfs, activity_classes, attributions = [], [], []
    for _, row in df.iterrows():
        hits = str(row[3]).split(';')
        best_tf, best_score = '', 0.0
        if len(hits) == 1:
            parts = hits[0].split(':')
            best_tf = parts[0].split('-')[0]
            try: best_score = float(parts[1])
            except (ValueError, IndexError): best_score = 0.0
        else:
            max_abs, all_names = -1.0, set()
            for hit in hits:
                parts = hit.split(':')
                tf_name = parts[0].split('-')[0]
                all_names.add(tf_name)
                try:
                    score = float(parts[1])
                    if abs(score) > max_abs:
                        max_abs, best_tf, best_score = abs(score), tf_name, score
                except (ValueError, IndexError):
                    continue
            if len(all_names) == 1:
                best_tf = list(all_names)[0]
        rep_tfs.append(best_tf)
        attributions.append(best_score)
        activity_classes.append('Activator' if best_score > 0 else ('Repressor' if best_score < 0 else 'Neutral'))
    df = df.copy()
    df['representative_tf'] = rep_tfs
    df['activity_class']    = activity_classes
    df['attribution']       = attributions
    return df


# =============================================================================
# MPAC
# =============================================================================

print("Processing MPAC seqlets...")
for cell in ['K562', 'HEPG2', 'SKNSH']:
    df = pd.read_csv(f"{BEDS}/{cell}_bedOps_merged_noMin_seqlets_01_v2.bed", sep='\t', header=None)
    out = assign_representative_tf(df)
    outpath = f"{BEDS}/mpac_{cell.lower()}_merged_collapsed_repTFs_v2.bed"
    out.to_csv(outpath, sep='\t', index=False, header=None)
    print(f"  {cell}: {len(out)} intervals -> {os.path.basename(outpath)}")

# =============================================================================
# PARM (tangermeme)
# =============================================================================

print("Processing PARM tangermeme seqlets...")
for cell in ['k562', 'hepg2']:
    df = pd.read_csv(f"{BEDS}/parm_{cell}_tangermeme_bedOps_merged_noMin_seqlets_01_v2.bed", sep='\t', header=None)
    out = assign_representative_tf_parm(df)
    outpath = f"{BEDS}/parm_tangermeme_{cell}_merged_collapsed_repTFs_v2.bed"
    out.to_csv(outpath, sep='\t', index=False, header=None)
    print(f"  {cell}: {len(out)} intervals -> {os.path.basename(outpath)}")

# =============================================================================
# PromoterAI
# =============================================================================

print("Processing PromoterAI seqlets...")
df = pd.read_csv(f"{BEDS}/promoterAI_bedOps_merged_noMin_seqlets_01_v2.bed", sep='\t', header=None)
out = assign_representative_tf(df)
outpath = f"{BEDS}/promoterAI_tangermeme_merge_collapsed_repTFs_v2.bed"
out.to_csv(outpath, sep='\t', index=False, header=None)
print(f"  {len(out)} intervals -> {os.path.basename(outpath)}")

print("Done.")
