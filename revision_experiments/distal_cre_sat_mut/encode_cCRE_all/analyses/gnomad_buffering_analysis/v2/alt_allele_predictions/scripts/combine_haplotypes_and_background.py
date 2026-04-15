#!/usr/bin/env python3

# Combines haplotype predictions for EH38E3356274_chr2_86221285_A_G_b38
# with the 100K background dELS sample (which includes ref_ref) and saves
# the concatenated result for downstream seqlet calling.

import pandas as pd
from tqdm import tqdm

BASE_DIR = "/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/alt_allele_predictions"
MPAC_DIR = f"{BASE_DIR}/mpac_preds"
PROCESSED_DIR = f"{BASE_DIR}/processed_data"


def vcf2df(pred_df):
    k_ref, k_alt, k_skew = [], [], []
    h_ref, h_alt, h_skew = [], [], []
    s_ref, s_alt, s_skew = [], [], []

    for i in tqdm(pred_df['INFO']):
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
        'chrom': pred_df['chrom'],
        'pos': pred_df['pos'],
        'id': pred_df['id'],
        'ref': pred_df['ref'],
        'alt': pred_df['alt'],
        'k562_ref_pred': k_ref,
        'k562_alt_pred': k_alt,
        'k562_skew_pred': k_skew,
        'hepg2_ref_pred': h_ref,
        'hepg2_alt_pred': h_alt,
        'hepg2_skew_pred': h_skew,
        'sknsh_ref_pred': s_ref,
        'sknsh_alt_pred': s_alt,
        'sknsh_skew_pred': s_skew
    })


def main():
    # Load and label haplotype predictions
    print("Loading haplotype predictions...")
    ref_alt = vcf2df(pd.read_csv(f"{MPAC_DIR}/EH38E3356274_chr2_86221285_A_G_b38_preds.vcf", sep='\t'))
    ref_alt['id'] = ref_alt['id'].apply(lambda x: f"{x}_chr2_86221285_A_G_b38")

    all_haplos = pd.concat([ref_alt])
    print(f"  Haplotype rows: {len(all_haplos):,}")

    # Load background sample (includes REF_REF)
    print("Loading 10K background sample...")
    background = pd.read_csv(
        f"{PROCESSED_DIR}/dELS_10K_satmut_preds_sample.tsv.gz",
        sep='\t'
    )
    print(f"  Background rows: {len(background):,}")

    # Concatenate
    combined = pd.concat([all_haplos, background])
    print(f"  Combined rows: {len(combined):,}")

    # Save
    out_path = f"{PROCESSED_DIR}/dELS_10K_satmut_preds_sample_with_EH38E3356274_chr2_86221285_A_G_b38.tsv.gz"
    print(f"Saving to {out_path}...")
    combined.to_csv(out_path, sep='\t', index=False,
                    compression={'method': 'gzip', 'compresslevel': 5})
    print("Done.")


if __name__ == "__main__":
    main()
