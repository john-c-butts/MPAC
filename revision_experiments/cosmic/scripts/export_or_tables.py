"""
Export the 2x2 tables behind the figures in COSMIC_analysis_v10.ipynb to
2x2_or_tables_for_cosmic_analyses.xlsx (one sheet per analysis), with each table's odds
ratio, 95% CI and p-value.

Setup and functions replicate COSMIC_analysis_v10.ipynb. Only tables behind plotted figures
(or printed ORs) are exported: Fig 2 forest plots (cells 13, 15, 17), the printed
active-variant OR (cell 22), emVar sharing (cells 32-38) and the tissue/histology-matched
ORs (cells 50, 55). See the 'Notes' sheet in the workbook for how each statistic is computed
(two-sided Fisher exact p-values; OR estimator differs between Fig 2 and the tissue/histology plots).

Run from the scripts/ directory:  python3 export_or_tables.py
"""
import pandas as pd
import numpy as np
from scipy import stats
from tqdm import tqdm

pd.options.mode.chained_assignment = None

OUT_FILE = '2x2_or_tables_for_cosmic_analyses.xlsx'

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
records = []          # one dict per 2x2 table; written to the workbook at the end
SHEET = None          # current sheet / group, set before each section below
GROUP = None


def add_table(title, col1, col2, row1_label, row2_label, table, stats_override=None):
    """stats_override = (OR, CI lower, CI upper, p) when the plotted value is not the
    default (fisher_exact OR, conditional CI, fisher p)."""
    (a, b), (c, d) = table
    if stats_override is None:
        oddsr, lo, hi, pval = or_and_ci(np.array([[a, b], [c, d]]))
    else:
        oddsr, lo, hi, pval = stats_override
    records.append(dict(kind='or', sheet=SHEET, group=GROUP, title=title, cols=(col1, col2),
                        rows=(row1_label, row2_label), counts=(int(a), int(b), int(c), int(d)),
                        stats=(float(oddsr), float(lo), float(hi), float(pval))))


def add_chi2_table(title, col1, col2, row1_label, row2_label, table):
    """2x2 tested with scipy chi2_contingency (default settings, as in the notebook)."""
    (a, b), (c, d) = table
    chi2, pval = stats.chi2_contingency(pd.DataFrame(table))[:2]
    records.append(dict(kind='chi2', sheet=SHEET, group=GROUP, title=title, cols=(col1, col2),
                        rows=(row1_label, row2_label), counts=(int(a), int(b), int(c), int(d)),
                        stats=(float(chi2), float(pval))))


# =========================================================
# 1. Headline odds ratios (notebook cells 10, 13 / Figure 2D)
# =========================================================
headline = [
    ('recurrent', 'Recurrent', 'Recurrent', 'Not Recurrent'),
    ('250bp_pro', '250bp Promoter', 'Promoter', 'Not Promoter'),
    ('cnc_pro', 'Cornell Non-Coding Cancer (CNC) Promoter', 'CNC Promoter', 'Not CNC Promoter'),
    ('recurrent_250bp_pro', 'Recurrent Promoter', 'Recurrent Promoter', 'Not Recurrent Promoter'),
    ('recurrent_cnc_pro', 'Recurrent CNC Promoter', 'Recurrent CNC Promoter', 'Not Recurrent CNC Promoter'),
]
SHEET = 'Fig 2D Headline ORs'
for col, title, poslab, neglab in headline:
    table = oddsR_calculator(cosmic_preds, col, 'emvar_any')
    add_table(f'emVar: {title}', 'emVar', 'non-emVar', poslab, neglab, table)

cosmic_preds4act_bin = reformat_df4act_bin(cosmic_preds)

# =========================================================
# 1b. Active variant is recurrent (cell 22, printed in the notebook)
# =========================================================
SHEET, GROUP = 'Active Variant OR', None
active_bins = ['(2.0, 3.0]', '(1.0, 2.0]', '(3.0, 4.0]', '(4.0, 5.0]', '(5.0, inf]']
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
# 2. Skew-bin breakdowns (cells 14, 15, 17: Fig 2E, 2F) - 7 feature sets x skew bins
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

SHEET = 'Fig 2E-F Skew Bins'
for df_in, feature, feature_ann in skew_specs:
    GROUP = feature_ann
    for r in oddsR_feature_by_skew_bin_all(df_in, skew_bins, feature):
        bin_str = str(r['Bin'])
        add_table(
            f'{feature_ann}: Mean Allelic Skew {bin_str}',
            feature_ann, f'Not {feature_ann}',
            f'Skew in {bin_str}', f'Skew outside {bin_str}',
            r['table']
        )

# =========================================================
# 3. Activity-bin breakdowns (cells 16, 17: Fig 2G; promoter subset not plotted) - 4 subsets x activity bins
# =========================================================
act_specs = [
    (cosmic_preds4act_bin, 'Recurrent emVar by MPAC Activity'),
    (cosmic_preds4act_bin[cosmic_preds4act_bin['250bp_pro'] == 0], 'Recurrent emVar by MPAC Activity: Non-Promoter'),
    (cosmic_preds4act_bin[(cosmic_preds4act_bin['250bp_pro'] == 0) & (cosmic_preds4act_bin['dhs'] == 1)], 'Recurrent emVar by MPAC Activity: Non-Promoter, In DHS'),
]

SHEET = 'Fig 2G Activity Bins'
for df_in, title_ann in act_specs:
    GROUP = title_ann
    for r in oddsR_emvar_by_act_bin_max(df_in):
        bin_str = str(r['Bin'])
        add_table(
            f'{title_ann}: MPAC Activity {bin_str}',
            'Recurrent', 'Unique',
            f'emVar in Activity {bin_str}', f'Not emVar in Activity {bin_str}',
            r['table']
        )

# =========================================================
# 4. emVar sharing across cell types: recurrent vs unique (cells 25-38; bar plot cell 32)
# =========================================================
SHARE_COLS = ['k562_only', 'hepg2_only', 'sknsh_only', 'k562_and_hepg2',
              'k562_and_sknsh', 'hepg2_and_sknsh', 'all']
SHARE_KEYS = {(1, 0, 0): 0, (0, 1, 0): 1, (0, 0, 1): 2, (1, 1, 0): 3,
              (1, 0, 1): 4, (0, 1, 1): 5, (1, 1, 1): 6, (0, 0, 0): None}


def emvar_sharing_summary(cosmic_df):
    """cell 25: one row per variant, one-hot over the 7 sharing categories."""
    rows = []
    f = cosmic_df.filter(['id', 'k562_emvar', 'hepg2_emvar', 'sknsh_emvar', 'recurrent', 'emvar_any'])
    for vid, k, h, s_, rec, anyv in zip(f['id'], f['k562_emvar'], f['hepg2_emvar'],
                                        f['sknsh_emvar'], f['recurrent'], f['emvar_any']):
        key = (k, h, s_)
        if key not in SHARE_KEYS:      # same as notebook: unmatched tuples are dropped
            continue
        onehot = [0] * 7
        if SHARE_KEYS[key] is not None:
            onehot[SHARE_KEYS[key]] = 1
        rows.append([vid] + onehot + [rec, anyv])
    return pd.DataFrame(rows, columns=['id'] + SHARE_COLS + ['recurrent', 'emvar_any'])


SHEET, GROUP = 'emVar Sharing', None
allemVarsSummary = emvar_sharing_summary(cosmic_preds[cosmic_preds['emvar_any'] == 1])
uniqueemVarsSummary = emvar_sharing_summary(cosmic_preds[(cosmic_preds['recurrent'] == 0) &
                                                         (cosmic_preds['emvar_any'] == 1)].copy())
recurrentemVarsSummary = emvar_sharing_summary(cosmic_preds[(cosmic_preds['recurrent'] == 1) &
                                                            (cosmic_preds['emvar_any'] == 1)].copy())

# cell 29-30: cell-type-specific ('only') annotations used by the tissue/histology ORs below
for ct in ['k562', 'hepg2', 'sknsh']:
    only_ids = set(allemVarsSummary.loc[allemVarsSummary[f'{ct}_only'] == 1, 'id'])
    cosmic_preds.loc[:, f'{ct}_emvar_only'] = [1 if i in only_ids else 0 for i in cosmic_preds['id']]

# cell 34: counts per sharing category, unique vs recurrent
mergeCount = pd.DataFrame({
    'Unique Variants': uniqueemVarsSummary[SHARE_COLS].sum(),
    'Recurrent Variants': recurrentemVarsSummary[SHARE_COLS].sum(),
})
u_tot, r_tot = mergeCount['Unique Variants'].sum(), mergeCount['Recurrent Variants'].sum()

# cell 35: emVars shared across all 3 cell types
u_all, r_all = mergeCount.loc['all', 'Unique Variants'], mergeCount.loc['all', 'Recurrent Variants']
add_chi2_table('emVar sharing: shared across all 3 cell types',
               'Unique emVars', 'Recurrent emVars',
               'Active in all 3 cell types', 'Not active in all 3 cell types',
               np.array([[u_all, r_all], [u_tot - u_all, r_tot - r_all]]))
# cell 37: emVars shared across >= 2 cell types
u_any, r_any = mergeCount['Unique Variants'].iloc[3:].sum(), mergeCount['Recurrent Variants'].iloc[3:].sum()
add_chi2_table('emVar sharing: shared across >= 2 cell types',
               'Unique emVars', 'Recurrent emVars',
               'Active in >= 2 cell types', 'Active in 1 cell type',
               np.array([[u_any, r_any], [u_tot - u_any, r_tot - r_any]]))

# =========================================================
# 5. Cell-matched recurrent vs unique ORs (cells 40-56; plots cells 50, 55)
#    Layout: columns Recurrent / Unique, rows emVar / Not emVar.
#    Plotted values: OR = scipy contingency.odds_ratio (conditional MLE), CI = its 95% CI,
#    p = Fisher exact (two-sided). The baseline point is plotted as the fisher_exact OR (cell 47).
# =========================================================
def cm_table(rec_df, uniq_df, emvar_col):
    n_r, n_u = len(rec_df), len(uniq_df)
    r_e = len(rec_df[rec_df[emvar_col] == 1])
    u_e = len(uniq_df[uniq_df[emvar_col] == 1])
    return np.array([[r_e, u_e], [n_r - r_e, n_u - u_e]])


def plotted_stats(table):
    res = stats.contingency.odds_ratio(table)
    ci = res.confidence_interval(confidence_level=0.95)
    return res.statistic, ci.low, ci.high, stats.fisher_exact(table)[1]


def add_cm_table(title, table, label, stats_override):
    add_table(title, 'Recurrent', 'Unique', f'emVar ({label})', f'Not emVar ({label})',
              table, stats_override)


# baseline (cell 40): all variants, any emVar
base_tab = cm_table(cosmic_preds[cosmic_preds['recurrent'] == 1],
                    cosmic_preds[cosmic_preds['recurrent'] == 0], 'emvar_any')
base_or = stats.fisher_exact(base_tab)[0]
_, lo_, hi_, p_ = plotted_stats(base_tab)
for SHEET in ['Tissue Matched', 'Histology Matched']:   # baseline is plotted in both figures
    GROUP = 'Baseline'
    add_cm_table('All Variants (baseline)', base_tab, 'any cell type', (base_or, lo_, hi_, p_))

CT_LABEL = {'k562': 'K562', 'hepg2': 'HepG2', 'sknsh': 'SKNSH'}
groups = [
    ('Primary Site', 'Haematopoietic and Lymphoid Tissue', 'primary_site', ['haematopoietic_and_lymphoid_tissue']),
    ('Primary Site', 'Liver', 'primary_site', ['liver']),
    ('Primary Site', 'Central Nervous System', 'primary_site', ['central_nervous_system']),
    ('Histology', 'Myeloid Leukaemia', 'histology_subtype_1',
     ['acute_myeloid_leukaemia', 'acute_myeloid_leukaemia_therapy_related', 'chronic_myeloid_leukaemia']),
    ('Histology', 'Hepatocellular Carcinoma', 'histology_subtype_1', ['hepatocellular_carcinoma']),
    ('Histology', "Contains 'Neuro'", 'primary_histology',
     [i for i in cosmic_preds['primary_histology'] if 'neuro' in i]),
]
for kind, name, col, values in groups:
    SHEET = 'Tissue Matched' if kind == 'Primary Site' else 'Histology Matched'
    GROUP = f'{kind} matched: {name}'
    sub = cosmic_preds[cosmic_preds[col].isin(values)]
    rec_df, uniq_df = sub[sub['recurrent'] == 1], sub[sub['recurrent'] == 0]
    specs = [('Any cell type', 'emvar_any')]
    specs += [(f'{CT_LABEL[ct]} non-specific', f'{ct}_emvar') for ct in CT_LABEL]
    specs += [(f'{CT_LABEL[ct]} cell-specific', f'{ct}_emvar_only') for ct in CT_LABEL]
    for lab, emvar_col in specs:
        t = cm_table(rec_df, uniq_df, emvar_col)
        add_cm_table(lab, t, lab, plotted_stats(t))

# =========================================================
# Write workbook (one sheet per analysis)
# =========================================================
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side

NOTES = [
    ('2x2 tables for COSMIC analyses (COSMIC_analysis_v10.ipynb)', 'title'),
    ('Only tables behind plotted figures, plus the printed active-variant OR, are included.', None),
    ('', None),
    ('Sheet', 'head'),
    ('Fig 2D Headline ORs: odds a variant is an emVar, by feature (Fig 2D forest plot).', None),
    ('Fig 2E-F Skew Bins: feature odds ratios by mean allelic skew bin (Fig 2E promoter panels, Fig 2F).', None),
    ('Fig 2G Activity Bins: recurrent vs unique emVar odds ratio by MPAC activity bin (Fig 2G).', None),
    ('Active Variant OR: odds that an active variant (MPAC activity > 1) is recurrent (printed in notebook).', None),
    ('emVar Sharing: chi-square tests of emVar sharing across cell types, recurrent vs unique (bar plot).', None),
    ('Tissue Matched / Histology Matched: recurrent vs unique emVar ORs for cell-type matched tumours.', None),
    ('', None),
    ('Statistics', 'head'),
    ('Odds ratio, 95% CI and p-value are as plotted/printed in the notebook.', None),
    ('Fisher p-values are two-sided (scipy.stats.fisher_exact default).', None),
    ('Fig 2 and Active Variant OR: OR = fisher_exact sample odds ratio (a*d / b*c); CI = scipy contingency.odds_ratio 95% CI.', None),
    ('Tissue/Histology sheets: OR and CI = scipy contingency.odds_ratio (conditional MLE); baseline plotted as the fisher_exact OR.', None),
    ('Chi-square: scipy chi2_contingency defaults (Yates continuity correction for 2x2).', None),
    ('Activity-bin tables: columns are Recurrent / Unique; rows are emVar in bin / not emVar in bin (all other bins).', None),
]

thin = Side(style='thin', color='999999')
TITLE_FILL = PatternFill('solid', fgColor='D9E1F2')
GROUP_FILL = PatternFill('solid', fgColor='1F3864')
SHEET_ORDER = ['Fig 2D Headline ORs', 'Fig 2E-F Skew Bins', 'Fig 2G Activity Bins',
               'Active Variant OR', 'emVar Sharing', 'Tissue Matched', 'Histology Matched']

wb = Workbook()
ws = wb.active
ws.title = 'Notes'
for i, (text, style) in enumerate(NOTES, start=1):
    ws.cell(row=i, column=1, value=text)
    if style == 'title':
        ws.cell(row=i, column=1).font = Font(bold=True, size=14)
    elif style == 'head':
        ws.cell(row=i, column=1).font = Font(bold=True)
ws.column_dimensions['A'].width = 130

for sheet in SHEET_ORDER:
    recs = [r for r in records if r['sheet'] == sheet]
    ws = wb.create_sheet(sheet)
    row, last_group = 1, object()
    for r in recs:
        if r['group'] != last_group:
            last_group = r['group']
            if last_group is not None:
                for col in range(1, 8):
                    ws.cell(row=row, column=col).fill = GROUP_FILL
                c = ws.cell(row=row, column=1, value=last_group)
                c.font = Font(bold=True, size=12, color='FFFFFF')
                row += 2
        # title bar
        for col in range(1, 8):
            ws.cell(row=row, column=col).fill = TITLE_FILL
        ws.cell(row=row, column=1, value=r['title']).font = Font(bold=True)
        row += 1
        # header
        if r['kind'] == 'or':
            heads = ['', r['cols'][0], r['cols'][1], 'Odds Ratio', '95% CI Lower', '95% CI Upper', 'Fisher p-value (two-sided)']
        else:
            heads = ['', r['cols'][0], r['cols'][1], 'Chi-square statistic', 'Chi-square p-value']
        for j, h in enumerate(heads, start=1):
            c = ws.cell(row=row, column=j, value=h or None)
            c.font = Font(bold=True)
            c.alignment = Alignment(horizontal='center', wrap_text=True, vertical='center')
            c.border = Border(bottom=thin)
        row += 1
        a_, b_, c_, d_ = r['counts']
        for lab, x, y, first in [(r['rows'][0], a_, b_, True), (r['rows'][1], c_, d_, False)]:
            ws.cell(row=row, column=1, value=lab).font = Font(bold=True)
            ws.cell(row=row, column=2, value=x).number_format = '#,##0'
            ws.cell(row=row, column=3, value=y).number_format = '#,##0'
            if first:
                for j, (v, fmt) in enumerate(zip(r['stats'], ['0.000', '0.000', '0.000', '0.00E+00'] if r['kind'] == 'or' else ['0.000', '0.00E+00']), start=4):
                    ws.cell(row=row, column=j, value=v).number_format = fmt
            row += 1
        row += 1   # blank row between tables
    ws.column_dimensions['A'].width = 62
    for col in 'BCDEFG':
        ws.column_dimensions[col].width = 20
    ws.freeze_panes = 'A1'

wb.save(OUT_FILE)
print('wrote', len(records), 'tables to', OUT_FILE)
for sh in SHEET_ORDER:
    print(f'  {sh}: {sum(1 for r in records if r["sheet"] == sh)}')
