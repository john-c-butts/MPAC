#!/usr/bin/env python3
"""
Shadow variant dataset — 605 leads, ALL groups use ALT background, no gnomAD AF filter.

All 605 leads routed through {enhancer_id}_{lead_variant} predictions regardless of
flip status. No gnomAD AF filter.

Bug fix vs. prior altBg script: variant_id_gtex and variant_id for lead rows are
corrected to GTEx orientation (ref=lref, alt=lalt) rather than the reversed orientation
that the alt background prediction file produces at the lead position.

Outputs:
  processed_data/final_shadow_variant_dataset_605_altBg_050726.tsv.gz
  processed_data/build_final_dataset_605_altBg_050726.log

Usage:
    python build_final_dataset_605_altBg_050726.py [--control-threshold 0.5]
"""

import argparse
import datetime
import logging
import os
import sys

import numpy as np
import pandas as pd
from tqdm import tqdm

CONTROL_THRESHOLD = 0.5

BASE = os.path.dirname(os.path.abspath(__file__))
RAW  = os.path.join(BASE, 'raw_data')

POST_QC_PATH   = os.path.join(RAW, '605_highPIP_shadow_variants_allTissue.tsv.zip')
PREDS_PATH     = os.path.join(RAW, 'all_predictions_with_background.tsv.gz')
ALT_SEQLET_DIR = os.path.join(RAW, 'alt_background_seqlets')
OUT_DIR        = os.path.join(BASE, 'processed_data')
OUT_PATH       = os.path.join(OUT_DIR, 'final_shadow_variant_dataset_605_altBg_050726.tsv.gz')
LOG_PATH       = os.path.join(OUT_DIR, 'build_final_dataset_605_altBg_050726.log')

FINAL_COLS = [
    'variant_id', 'variant_id_gtex', 'chrom', 'pos', 'ref', 'alt',
    'skew_pred', 'abs_skew', 'cell_type', 'enhancer_ids',
    'exceeds_lead_skew', 'is_lead_variant',
    'lead_variant_id', 'lead_variant_id_gtex',
    'lead_abs_skew', 'lead_pip', 'lead_tissue', 'lead_skew_sign',
    'phenotype_id', 'gene_name', 'gene_short', 'in_seqlet',
]


def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(message)s',
        handlers=[
            logging.FileHandler(LOG_PATH),
            logging.StreamHandler(sys.stdout),
        ],
    )


# ---------------------------------------------------------------------------
# Seqlet annotation (ALT background)
# ---------------------------------------------------------------------------

def load_seqlet_beds():
    beds = {}
    ct_map = {'k562': 'K562', 'hepg2': 'HepG2', 'sknsh': 'SKNSH'}
    for ct, label in ct_map.items():
        path = os.path.join(ALT_SEQLET_DIR, f'{label}_0.01.bed')
        if os.path.exists(path):
            beds[ct] = pd.read_csv(path, sep='\t', header=None, usecols=[0, 1, 2],
                                    names=['chrom', 'start', 'end'])
        else:
            logging.warning('ALT seqlet BED not found: %s', path)
    return beds


def build_seqlet_index(beds):
    index = {}
    for ct, df in beds.items():
        index[ct] = {}
        for chrom, grp in df.groupby('chrom'):
            index[ct][chrom] = (grp['start'].values, grp['end'].values)
    return index


def check_in_seqlet(chrom, pos, seqlet_index, cell_type):
    if cell_type not in seqlet_index:
        return False
    chrom_idx = seqlet_index[cell_type]
    if chrom not in chrom_idx:
        return False
    starts, ends = chrom_idx[chrom]
    return bool(np.any((starts <= pos - 1) & (pos - 1 < ends)))


def annotate_seqlets(df, seqlet_index):
    results = [
        check_in_seqlet(row['chrom'], row['pos'], seqlet_index, row['cell_type'])
        for _, row in tqdm(df.iterrows(), total=len(df), desc='  Annotating seqlets')
    ]
    df = df.copy()
    df['in_seqlet'] = results
    return df


# ---------------------------------------------------------------------------
# Metadata helper
# ---------------------------------------------------------------------------

def _add_metadata(subset, lead_row, enhancer_id, lead_vtx_gtex,
                  lead_abs_skew, lead_skew_sign, exceeds, is_lead, skew_ref):
    subset = subset.copy()
    subset['cell_type']            = lead_row['cell_type']
    subset['enhancer_ids']         = enhancer_id
    subset['exceeds_lead_skew']    = exceeds
    subset['is_lead_variant']      = is_lead
    subset['lead_variant_id']      = lead_row['lead_variant_id']
    subset['lead_abs_skew']        = lead_abs_skew
    subset['lead_pip']             = lead_row['lead_pip']
    subset['lead_tissue']          = lead_row['lead_tissue']
    subset['lead_skew_sign']       = lead_skew_sign
    subset['lead_variant_id_gtex'] = lead_vtx_gtex
    subset['phenotype_id']         = lead_row['phenotype_id']
    subset['gene_name']            = lead_row.get('gene_name', '')
    subset['gene_short']           = lead_row.get('gene_short', '')
    return subset


# ---------------------------------------------------------------------------
# ALT background row derivation
# ---------------------------------------------------------------------------

def derive_alt_rows(group_lead_rows, alt_preds):
    """
    Derive lead + shadow + control rows from alt background predictions for all leads.

    Lead rows: variant_id_gtex and variant_id are corrected to GTEx orientation
    (ref=lref, alt=lalt) because the alt background prediction file stores the lead
    in reversed orientation (ref=lalt, alt=lref).

    Shadow/control rows: orientation is unchanged (reference genome orientation from
    prediction file, consistent with GTEx for non-lead variants).
    """
    new_rows = []

    for (enhancer_id, lead_vtx_gtex), lead_grp in group_lead_rows.groupby(
            ['enhancer_ids', 'lead_variant_id_gtex']):

        job_id    = f'{enhancer_id}_{lead_vtx_gtex.replace("_b38", "")}'
        enh_preds = alt_preds[alt_preds['id'] == job_id].copy()

        if len(enh_preds) == 0:
            logging.warning('no alt_preds for %s', job_id)
            continue

        enh_preds['variant_id'] = (enh_preds['chrom'] + ':' + enh_preds['pos'].astype(str)
                                   + ':' + enh_preds['ref'] + ':' + enh_preds['alt'])
        enh_preds['variant_id_gtex'] = (enh_preds['chrom'] + '_'
                                        + enh_preds['pos'].astype(str) + '_'
                                        + enh_preds['ref'] + '_'
                                        + enh_preds['alt'] + '_b38')

        fields = lead_vtx_gtex.replace('_b38', '').split('_')
        lchrom, lpos, lref, lalt = fields[0], int(fields[1]), fields[2], fields[3]
        lead_vid_orig = f'{lchrom}:{lpos}:{lref}:{lalt}'
        lead_vid_alt  = f'{lchrom}:{lpos}:{lalt}:{lref}'

        for _, lead_row in lead_grp.iterrows():
            ct       = lead_row['cell_type']
            skew_col = f'{ct}_skew_pred'
            if skew_col not in enh_preds.columns:
                continue

            # Lead appears reversed (ref=lalt, alt=lref) in alt background
            lead_in_alt = enh_preds[
                (enh_preds['pos'] == lpos) &
                (enh_preds['ref'] == lalt) &
                (enh_preds['alt'] == lref)
            ]
            if len(lead_in_alt) == 0:
                logging.warning('lead not in alt_preds: %s %s', enhancer_id, ct)
                continue

            lead_skew      = lead_in_alt.iloc[0][skew_col]
            lead_abs_skew  = abs(lead_skew)
            lead_skew_sign = int(np.sign(lead_skew))

            ct_all = enh_preds[['variant_id', 'variant_id_gtex', 'chrom', 'pos',
                                 'ref', 'alt', skew_col]].rename(
                columns={skew_col: 'skew_pred'}).copy()
            ct_all['abs_skew'] = ct_all['skew_pred'].abs()
            ct_all = ct_all[
                (ct_all['variant_id'] != lead_vid_orig) &
                (ct_all['variant_id'] != lead_vid_alt)
            ]

            # Lead row — corrected to GTEx orientation (ref=lref, alt=lalt)
            lr = pd.DataFrame([{
                'variant_id':      lead_vid_orig,
                'variant_id_gtex': lead_vtx_gtex,
                'chrom':           lchrom,
                'pos':             lpos,
                'ref':             lref,
                'alt':             lalt,
                'skew_pred':       lead_skew,
                'abs_skew':        lead_abs_skew,
            }])
            lr = _add_metadata(lr, lead_row, enhancer_id, lead_vtx_gtex,
                                lead_abs_skew, lead_skew_sign,
                                exceeds=False, is_lead=True, skew_ref=lead_skew)
            new_rows.append(lr)

            shadow_mask  = ((ct_all['abs_skew'] > lead_abs_skew) &
                            (np.sign(ct_all['skew_pred']) == lead_skew_sign))
            control_mask = ct_all['abs_skew'] < CONTROL_THRESHOLD

            for mask, exceeds in [(shadow_mask, True), (control_mask, False)]:
                subset = ct_all[mask].copy()
                if len(subset) == 0:
                    continue
                subset = _add_metadata(subset, lead_row, enhancer_id, lead_vtx_gtex,
                                       lead_abs_skew, lead_skew_sign,
                                       exceeds=exceeds, is_lead=False,
                                       skew_ref=lead_skew)
                new_rows.append(subset)

    return pd.concat(new_rows, ignore_index=True) if new_rows else pd.DataFrame()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    global CONTROL_THRESHOLD
    parser = argparse.ArgumentParser()
    parser.add_argument('--control-threshold', type=float, default=CONTROL_THRESHOLD)
    args   = parser.parse_args()
    CONTROL_THRESHOLD = args.control_threshold

    os.makedirs(OUT_DIR, exist_ok=True)
    setup_logging()

    logging.info('Run started: %s', datetime.datetime.now().isoformat())
    logging.info('control_threshold=%.2f', CONTROL_THRESHOLD)

    # --- Load 605 ---
    logging.info('Step 1: Load 605_highPIP_shadow_variants_allTissue.tsv.zip')
    post_qc = pd.read_csv(POST_QC_PATH, sep='\t')
    logging.info('  %d rows, %d unique lead variants',
                 len(post_qc), post_qc['lead_variant_id_gtex'].nunique())

    # --- Synthesize any missing lead rows ---
    logging.info('Step 2: Check for missing lead rows')
    has_lead      = set(post_qc[post_qc['is_lead_variant']]['lead_variant_id_gtex'])
    all_leads     = set(post_qc['lead_variant_id_gtex'].unique())
    missing_leads = all_leads - has_lead
    if missing_leads:
        logging.info('  Synthesizing %d missing lead rows: %s', len(missing_leads), missing_leads)
        syn_rows = []
        for lvid in missing_leads:
            grp = post_qc[post_qc['lead_variant_id_gtex'] == lvid]
            for ct, ct_grp in grp.groupby('cell_type'):
                ref_row  = ct_grp.iloc[0]
                lead_abs = ref_row['lead_abs_skew']
                lead_sgn = ref_row['lead_skew_sign']
                syn = {col: np.nan for col in post_qc.columns}
                syn.update({
                    'variant_id':           ref_row['lead_variant_id'],
                    'variant_id_gtex':      lvid,
                    'chrom':                lvid.split('_')[0],
                    'pos':                  int(lvid.split('_')[1]),
                    'ref':                  lvid.split('_')[2],
                    'alt':                  lvid.split('_')[3],
                    'cell_type':            ct,
                    'skew_pred':            lead_abs * lead_sgn,
                    'abs_skew':             lead_abs,
                    'enhancer_ids':         ref_row['enhancer_ids'],
                    'exceeds_lead_skew':    False,
                    'is_lead_variant':      True,
                    'lead_variant_id':      ref_row['lead_variant_id'],
                    'lead_abs_skew':        lead_abs,
                    'lead_pip':             ref_row['lead_pip'],
                    'lead_tissue':          ref_row['lead_tissue'],
                    'lead_skew_sign':       lead_sgn,
                    'lead_variant_id_gtex': lvid,
                    'phenotype_id':         ref_row['phenotype_id'],
                    'gene_name':            ref_row.get('gene_name', ''),
                    'gene_short':           ref_row.get('gene_short', ''),
                })
                syn_rows.append(syn)
        post_qc = pd.concat([post_qc, pd.DataFrame(syn_rows)], ignore_index=True)
        logging.info('  After synthesis: %d rows', len(post_qc))
    else:
        logging.info('  No missing lead rows.')

    # --- All leads use ALT background ---
    logging.info('Step 3: Build ALT job IDs for all leads')
    lead_rows   = post_qc[post_qc['is_lead_variant']].copy()
    all_keys    = set(zip(lead_rows['enhancer_ids'], lead_rows['lead_variant_id_gtex']))
    alt_job_ids = {f'{eid}_{lvid.replace("_b38", "")}' for eid, lvid in all_keys}
    logging.info('  %d ALT job IDs for %d (enhancer, lead) groups',
                 len(alt_job_ids), len(all_keys))

    # --- Load predictions ---
    logging.info('Step 4: Load ALT background predictions')
    chunks = []
    for chunk in tqdm(pd.read_csv(PREDS_PATH, sep='\t', chunksize=500_000),
                      desc='  Reading predictions'):
        mask = chunk['id'].isin(alt_job_ids)
        if mask.any():
            chunks.append(chunk[mask])
    alt_preds = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()
    n_found   = alt_preds['id'].nunique()
    n_missing = len(alt_job_ids - set(alt_preds['id']))
    logging.info('  %d rows, %d / %d job IDs found, %d missing',
                 len(alt_preds), n_found, len(alt_job_ids), n_missing)
    if n_missing:
        logging.warning('  Missing: %s', sorted(alt_job_ids - set(alt_preds['id']))[:5])

    # --- Derive rows ---
    logging.info('Step 5: Derive rows from ALT background (all leads)')
    combined = derive_alt_rows(lead_rows, alt_preds)
    if len(combined) == 0:
        logging.error('No rows derived.')
        sys.exit(1)
    n_l = combined['is_lead_variant'].sum()
    n_s = (combined['exceeds_lead_skew'] & ~combined['is_lead_variant']).sum()
    n_c = (~combined['exceeds_lead_skew'] & ~combined['is_lead_variant']).sum()
    logging.info('  %d lead, %d shadow, %d control rows', n_l, n_s, n_c)

    # --- Filter to groups with >= 1 shadow ---
    logging.info('Step 6: Filter to groups with >= 1 shadow')
    shadow_groups = set(zip(
        combined[combined['exceeds_lead_skew']]['enhancer_ids'],
        combined[combined['exceeds_lead_skew']]['lead_variant_id_gtex'],
    ))
    combined = combined[combined.apply(
        lambda r: (r['enhancer_ids'], r['lead_variant_id_gtex']) in shadow_groups, axis=1
    )].copy()
    logging.info('  %d groups retained, %d total rows', len(shadow_groups), len(combined))

    # --- Seqlet annotation ---
    logging.info('Step 7: Seqlet annotation (ALT background)')
    beds = load_seqlet_beds()
    seqlet_index = build_seqlet_index(beds)
    logging.info('  Loaded %d BED files', len(beds))
    combined = annotate_seqlets(combined, seqlet_index)
    n_in = combined['in_seqlet'].sum()
    logging.info('  %d / %d in seqlet (%.1f%%)',
                 n_in, len(combined), 100 * n_in / len(combined))

    # --- Select and order final columns ---
    combined = combined[FINAL_COLS]

    # --- Write ---
    combined.to_csv(OUT_PATH, sep='\t', index=False,
                    compression={'method': 'gzip', 'compresslevel': 5})
    logging.info('Wrote %s', OUT_PATH)

    logging.info('=== Final Summary ===')
    logging.info('Lead rows:    %d', combined['is_lead_variant'].sum())
    logging.info('Shadow rows:  %d', (combined['exceeds_lead_skew'] & ~combined['is_lead_variant']).sum())
    logging.info('Control rows: %d', (~combined['exceeds_lead_skew'] & ~combined['is_lead_variant']).sum())
    logging.info('Unique leads: %d', combined['lead_variant_id_gtex'].nunique())
    logging.info('Run complete: %s', datetime.datetime.now().isoformat())


if __name__ == '__main__':
    main()
