"""
Export every 2x2 table behind the odds ratios in COSMIC_analysis_v10.ipynb to
2x2_or_tables_for_cosmic_analyses.tsv, along with each table's odds ratio, 95% CI
and Fisher's exact test p-value.

Setup and functions replicate COSMIC_analysis_v10.ipynb (cells 1-7, 18, 20, 26).
Run from the scripts/ directory:  python3 export_or_tables.py

Output block format (one block per table):
    OR <title>
    <TAB>col1<TAB>col2<TAB>Odds Ratio<TAB>95% CI Lower<TAB>95% CI Upper<TAB>Fisher p-value
    row1<TAB>a<TAB>b<TAB><OR><TAB><lower><TAB><upper><TAB><p>
    row2<TAB>c<TAB>d
OR = scipy fisher_exact odds ratio (a*d / b*c); CI = scipy
stats.contingency.odds_ratio(table).confidence_interval(0.95), as in the notebook.
"""
import pandas as pd
import numpy as np
from scipy import stats
from tqdm import tqdm

pd.options.mode.chained_assignment = None

OUT_FILE = '2x2_or_tables_for_cosmic_analyses.tsv'

# ---- notebook setup (cells 6, 7) ----
cosmic_preds = pd.read_csv(
    '../processed_data/all.cosmic.v98.autosome.with.10bp.indels.annotated.017.081626.tsv',
    sep='\t', low_memory=False)

cosmic_preds['mean_ref_bin'] = pd.cut(cosmic_preds['mean_ref_act'], [float('-inf'), 1, 2, 3, 4, 5, float('inf')])
cosmic_preds['k562_act_bin'] = pd.cut(cosmic_preds['k562_max_act'], [float('-inf'), 1, 2, 3, 4, 5, float('inf')])
cosmic_preds['hepg2_act_bin'] = pd.cut(cosmic_preds['hepg2_max_act'], [float('-inf'), 1, 2, 3, 4, 5, float('inf')])
cosmic_preds['sknsh_act_bin'] = pd.cut(cosmic_preds['sknsh_max_act'], [float('-inf'), 1, 2, 3, 4, 5, float('inf')])

cosmic_preds = cosmic_preds[(cosmic_preds['vep_coding_related'] == 0) | (cosmic_preds['vep_most_severe'] == 'synonymous_variant')]

print('N after filtering:', len(cosmic_preds))


# ---- notebook functions (cells 2, 3, 4, 5) ----
def or_and_ci(table):
    """Odds ratio and two-sided p (fisher_exact) and 95% CI, as computed in the notebook."""
    oddsr, pval = stats.fisher_exact(table)
    ci = stats.contingency.odds_ratio(table).confidence_interval(0.95)
    return oddsr, ci.low, ci.high, pval


def oddsR_calculator(df, column_of_interest, emvar_col):
    pos = df[column_of_interest] == 1
    a = df.loc[pos, emvar_col].sum()
    b = pos.sum() - a
    c = df.loc[~pos, emvar_col].sum()
    d = (~pos).sum() - c
    return np.array([[a, b], [c, d]])


def oddsR_feature_by_skew_bin_all(preds_df, list_of_bins, feature):
    df = preds_df.copy()
    df['mean_skew_bin'] = pd.cut(df['mean_skew'], list_of_bins)
    total_with_feature = df[feature].sum()
    total_without_feature = len(df) - total_with_feature
    results = []
    for bin_label, bin_df in df.groupby('mean_skew_bin', observed=True):
        a = bin_df[feature].sum()
        b = len(bin_df) - a
        c = total_with_feature - a
        d = total_without_feature - b
        results.append({'Bin': bin_label, 'table': np.array([[a, b], [c, d]])})
    return sorted(results, key=lambda r: r['Bin'])


def reformat_df4act_bin(df):
    max_skew_cell = []
    max_skew_l2fc = []
    for k_skew, h_skew, s_skew in zip(df['k562_skew_pred'], df['hepg2_skew_pred'], df['sknsh_skew_pred']):
        skew_list = [k_skew, h_skew, s_skew]
        max_skew = max(skew_list, key=abs)
        if skew_list[0] == max_skew:
            max_skew_cell.append('k562'); max_skew_l2fc.append(max_skew)
        elif skew_list[1] == max_skew:
            max_skew_cell.append('hepg2'); max_skew_l2fc.append(max_skew)
        else:
            max_skew_cell.append('sknsh'); max_skew_l2fc.append(max_skew)
    df['max_skew_cell'] = max_skew_cell
    k_act_bin = df['k562_act_bin'].tolist()
    h_act_bin = df['hepg2_act_bin'].tolist()
    s_act_bin = df['sknsh_act_bin'].tolist()
    act_bin = []
    for i, j in tqdm(zip(range(len(df)), df['max_skew_cell'])):
        if j == 'k562':
            act_bin.append(k_act_bin[i])
        elif j == 'hepg2':
            act_bin.append(h_act_bin[i])
        else:
            act_bin.append(s_act_bin[i])
    df['max_act_bin'] = act_bin
    df['max_skew'] = max_skew_l2fc
    return df


def oddsR_emvar_by_act_bin_max(preds_df, skew_threshold=0.5):
    df = preds_df.copy()
    df['is_emvar'] = (df['max_skew'].abs() > skew_threshold).astype(int)
    recurrent = df[df['recurrent'] == 1]
    unique = df[df['recurrent'] != 1]
    total_recurrent = len(recurrent)
    total_unique = len(unique)
    results = []
    for bin_label in df['max_act_bin'].dropna().unique():
        bin_r = recurrent[recurrent['max_act_bin'] == bin_label]
        bin_u = unique[unique['max_act_bin'] == bin_label]
        a = bin_r['is_emvar'].sum()
        c = total_recurrent - a
        b = bin_u['is_emvar'].sum()
        d = total_unique - b
        results.append({'Bin': bin_label, 'table': np.array([[a, b], [c, d]])})
    return sorted(results, key=lambda r: r['Bin'])


# ---- output helper ----
out_lines = []


def add_table(title, col1, col2, row1_label, row2_label, table):
    (a, b), (c, d) = table
    oddsr, lo, hi, pval = or_and_ci(np.array([[a, b], [c, d]]))
    out_lines.append(f'OR {title}')
    out_lines.append(f'\t{col1}\t{col2}\tOdds Ratio\t95% CI Lower\t95% CI Upper\tFisher p-value')
    out_lines.append(f'{row1_label}\t{a}\t{b}\t{oddsr:.6g}\t{lo:.6g}\t{hi:.6g}\t{pval:.3e}')
    out_lines.append(f'{row2_label}\t{c}\t{d}')


# =========================================================
# 1. Headline odds ratios (notebook cell 10 / Figure 2D)
# =========================================================
headline = [
    ('recurrent', 'Recurrent', 'Recurrent', 'Not Recurrent'),
    ('250bp_pro', '250bp Promoter', 'Promoter', 'Not Promoter'),
    ('cnc_pro', 'Cornell Non-Coding Cancer (CNC) Promoter', 'CNC Promoter', 'Not CNC Promoter'),
    ('recurrent_250bp_pro', 'Recurrent Promoter', 'Recurrent Promoter', 'Not Recurrent Promoter'),
    ('recurrent_cnc_pro', 'Recurrent CNC Promoter', 'Recurrent CNC Promoter', 'Not Recurrent CNC Promoter'),
]
for col, title, poslab, neglab in headline:
    table = oddsR_calculator(cosmic_preds, col, 'emvar_any')
    add_table(f'emVar: {title}', 'emVar', 'non-emVar', poslab, neglab, table)

# =========================================================
# 2. Standalone OR: active variant is recurrent (cell 26)
#    NOTE: reproduced as written in the notebook; the string '(5.0, inf' has no
#    closing ']' so the (5.0, inf] bin is NOT matched as active.
# =========================================================
cosmic_preds4act_bin = reformat_df4act_bin(cosmic_preds)
active_bins = ['(2.0, 3.0]', '(1.0, 2.0]', '(3.0, 4.0]', '(4.0, 5.0]', '(5.0, inf']
is_active = cosmic_preds4act_bin['max_act_bin'].astype(str).isin(active_bins)
active_bin = cosmic_preds4act_bin[is_active]
inactive_bin = cosmic_preds4act_bin[~is_active]
a = active_bin['recurrent'].sum()
b = len(active_bin) - a
c = inactive_bin['recurrent'].sum()
d = len(inactive_bin) - c
add_table('Recurrent: Active Variant (MPAC Activity > 1)', 'Recurrent', 'Not Recurrent',
          'Active', 'Inactive', np.array([[a, b], [c, d]]))

# =========================================================
# 3. Skew-bin breakdowns (cell 18) - 7 feature sets x skew bins
# =========================================================
skew_bins = [float('-inf'), -1, -.5, -.2, -.05, .05, .2, .5, 1, float('inf')]

skew_specs = [
    (cosmic_preds, 'recurrent', 'Recurrent'),
    (cosmic_preds[cosmic_preds['250bp_pro'] == 0], 'recurrent', 'Recurrent: Non-Promoter'),
    (cosmic_preds[(cosmic_preds['250bp_pro'] == 0) & (cosmic_preds['dhs'] == 1)], 'recurrent', 'Recurrent: In DHS'),
    (cosmic_preds, '250bp_pro', '250bp Promoter'),
    (cosmic_preds, 'cnc_pro', 'CNC Promoter'),
    (cosmic_preds, 'recurrent_250bp_pro', 'Recurrent Promoter'),
    (cosmic_preds, 'recurrent_cnc_pro', 'Recurrent CNC Promoter'),
]

for df_in, feature, feature_ann in skew_specs:
    for r in oddsR_feature_by_skew_bin_all(df_in, skew_bins, feature):
        bin_str = str(r['Bin'])
        add_table(
            f'{feature_ann}: Mean Allelic Skew {bin_str}',
            feature_ann, f'Not {feature_ann}',
            f'Skew in {bin_str}', f'Skew outside {bin_str}',
            r['table']
        )

# =========================================================
# 4. Activity-bin breakdowns (cell 20) - 4 subsets x activity bins
# =========================================================
act_specs = [
    (cosmic_preds4act_bin, 'Recurrent emVar by MPAC Activity'),
    (cosmic_preds4act_bin[cosmic_preds4act_bin['250bp_pro'] == 0], 'Recurrent emVar by MPAC Activity: Non-Promoter'),
    (cosmic_preds4act_bin[(cosmic_preds4act_bin['250bp_pro'] == 0) & (cosmic_preds4act_bin['dhs'] == 1)], 'Recurrent emVar by MPAC Activity: Non-Promoter, In DHS'),
    (cosmic_preds4act_bin[cosmic_preds4act_bin['250bp_pro'] == 1], 'Recurrent emVar by MPAC Activity: Promoter'),
]

for df_in, title_ann in act_specs:
    for r in oddsR_emvar_by_act_bin_max(df_in):
        bin_str = str(r['Bin'])
        add_table(
            f'{title_ann}: MPAC Activity {bin_str}',
            f'emVar in Activity {bin_str}', f'Not emVar in Activity {bin_str}',
            'Recurrent', 'Unique',
            r['table']
        )

with open(OUT_FILE, 'w') as f:
    f.write('\n'.join(out_lines) + '\n')

print('wrote', len(out_lines), 'lines to', OUT_FILE)
print('total tables:', sum(1 for l in out_lines if l.startswith('OR ')))
