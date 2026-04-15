# ============================================================================
# ADDITIONAL ANALYSES - Add these cells to your notebook
# ============================================================================
# Fixes the violin plot error and adds:
# 1. Absolute skew quantile analysis
# 2. Odds ratio analysis for conserved bases (phyloP >= 2.27)
# ============================================================================

# %% Cell: FIXED PHYLOP DISTRIBUTION VISUALIZATION
# =============================================================================
# PHYLOP DISTRIBUTION VISUALIZATION (FIXED)
# =============================================================================

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: Overall phyloP distribution
ax1 = axes[0, 0]
ax1.hist(sample_data['phyloP'].clip(-5, 10), bins=100, density=True, alpha=0.7, color='steelblue')
ax1.axvline(0, color='red', linestyle='--', linewidth=2, label='phyloP = 0')
ax1.axvline(2.27, color='green', linestyle='--', linewidth=2, label='phyloP = 2.27 (conserved)')
ax1.axvline(sample_data['phyloP'].median(), color='orange', linestyle='-', linewidth=2,
            label=f'Median = {sample_data["phyloP"].median():.2f}')
ax1.set_xlabel('phyloP Score')
ax1.set_ylabel('Density')
ax1.set_title('Overall phyloP Distribution')
ax1.legend()

# Plot 2: phyloP distribution by extreme skew bins
ax2 = axes[0, 1]
for bin_idx, color, label in [(0, 'indianred', 'Strong negative skew (<-1.5)'),
                               (5, 'gray', 'Neutral skew'),
                               (10, 'teal', 'Strong positive skew (>1.5)')]:
    subset = sample_data[sample_data['fixed_bin'] == bin_idx]['phyloP'].clip(-5, 10)
    if len(subset) > 0:
        ax2.hist(subset, bins=50, density=True, alpha=0.5, label=label, color=color)
ax2.axvline(0, color='black', linestyle='--', linewidth=1)
ax2.axvline(2.27, color='green', linestyle='--', linewidth=1, label='Conserved threshold')
ax2.set_xlabel('phyloP Score')
ax2.set_ylabel('Density')
ax2.set_title('phyloP Distribution by Skew Bin')
ax2.legend()

# Plot 3: Violin plot - FIXED to handle empty bins
ax3 = axes[1, 0]
subsample = sample_data.sample(n=min(100000, len(sample_data)), random_state=42)

# Filter out empty bins and track which bins have data
violin_data = []
violin_positions = []
violin_labels = []
for i in range(11):
    bin_data = subsample[subsample['fixed_bin'] == i]['phyloP'].clip(-5, 10).values
    if len(bin_data) > 0:  # Only include bins with data
        violin_data.append(bin_data)
        violin_positions.append(i)
        violin_labels.append(bin_labels[i])

if violin_data:  # Only plot if we have data
    parts = ax3.violinplot(violin_data, positions=violin_positions, showmeans=True, showmedians=True)
    ax3.set_xticks(violin_positions)
    ax3.set_xticklabels(violin_labels, rotation=45, ha='right')
ax3.set_xlabel('Predicted Skew Bin')
ax3.set_ylabel('phyloP Score')
ax3.set_title('phyloP Distribution (Violin Plot)')
ax3.axhline(0, color='gray', linestyle='--', alpha=0.5)
ax3.axhline(2.27, color='green', linestyle='--', alpha=0.5)

# Plot 4: 2D histogram
ax4 = axes[1, 1]
h = ax4.hist2d(sample_data['mean_skew_pred'].clip(-2, 2),
               sample_data['phyloP'].clip(-5, 10),
               bins=50, cmap='viridis', norm=plt.matplotlib.colors.LogNorm())
plt.colorbar(h[3], ax=ax4, label='Count (log scale)')
ax4.axhline(0, color='red', linestyle='--', linewidth=1, alpha=0.7)
ax4.axhline(2.27, color='green', linestyle='--', linewidth=1, alpha=0.7)
ax4.axvline(0, color='red', linestyle='--', linewidth=1, alpha=0.7)
ax4.set_xlabel('Mean Skew Prediction')
ax4.set_ylabel('phyloP Score')
ax4.set_title('Joint Distribution (2D Histogram)')

for ax in axes.flat:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.show()

print('\nphyloP Summary Statistics:')
print(f"  Overall mean: {sample_data['phyloP'].mean():.4f}")
print(f"  Overall median: {sample_data['phyloP'].median():.4f}")
print(f"  % Conserved (phyloP > 0): {100 * (sample_data['phyloP'] > 0).mean():.1f}%")
print(f"  % Strongly conserved (phyloP >= 2.27): {100 * (sample_data['phyloP'] >= 2.27).mean():.1f}%")
print(f"  % Accelerated (phyloP < 0): {100 * (sample_data['phyloP'] < 0).mean():.1f}%")


# %% Cell: ABSOLUTE SKEW QUANTILE ANALYSIS
# =============================================================================
# ABSOLUTE SKEW QUANTILE ANALYSIS
# =============================================================================
# Use |skew| to examine effect magnitude regardless of direction

sample_data['abs_skew'] = sample_data['mean_skew_pred'].abs()

# Create quantile bins based on absolute skew
n_bins = 10
sample_data['abs_quantile_bin'] = pd.qcut(sample_data['abs_skew'], q=n_bins, labels=False, duplicates='drop')

# Get quantile boundaries
abs_quantile_boundaries = sample_data['abs_skew'].quantile(np.linspace(0, 1, n_bins + 1))
print("Absolute skew quantile boundaries:")
print(abs_quantile_boundaries.round(3).to_string())

# Calculate stats
abs_quantile_stats = sample_data.groupby('abs_quantile_bin')['phyloP'].agg(['mean', 'std', 'count']).reset_index()
abs_quantile_stats['sem'] = abs_quantile_stats['std'] / np.sqrt(abs_quantile_stats['count'])

# Plot
fig, ax = plt.subplots(figsize=(10, 6))
x = range(len(abs_quantile_stats))
ax.errorbar(x, abs_quantile_stats['mean'], yerr=abs_quantile_stats['sem']*1.96,
            marker='o', color='darkviolet', linewidth=2, markersize=8, capsize=3)

q_labels = [f'{abs_quantile_boundaries.iloc[i]:.2f}' for i in range(len(abs_quantile_boundaries)-1)]
ax.set_xticks(x)
ax.set_xticklabels(q_labels, rotation=45, ha='right')
ax.set_xlabel('|Skew| (Quantile Bin Lower Bound)')
ax.set_ylabel('Mean phyloP Score')
ax.set_title('Conservation vs Absolute Skew Magnitude (Quantile Bins)')
ax.axhline(0, color='gray', linestyle='--', alpha=0.5)
ax.grid(axis='y', linestyle=':', alpha=0.6)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.show()

print(f"\nAbsolute skew analysis:")
print(f"  Lowest |skew| bin mean phyloP: {abs_quantile_stats['mean'].iloc[0]:.4f}")
print(f"  Highest |skew| bin mean phyloP: {abs_quantile_stats['mean'].iloc[-1]:.4f}")


# %% Cell: ODDS RATIO ANALYSIS FOR CONSERVED BASES
# =============================================================================
# ODDS RATIO ANALYSIS: Probability of being conserved (phyloP >= 2.27)
# =============================================================================
from scipy.stats import fisher_exact

# Define conserved threshold
CONSERVED_THRESHOLD = 2.27

# Add conserved flag
sample_data['is_conserved'] = sample_data['phyloP'] >= CONSERVED_THRESHOLD

# Calculate odds ratio relative to neutral bin (bin 5)
# Odds = P(conserved) / P(not conserved) = n_conserved / n_not_conserved

def calculate_odds_ratio_vs_reference(df, bin_col, ref_bin=5):
    """Calculate odds ratios for each bin vs a reference bin."""
    results = []

    # Get reference bin counts
    ref_data = df[df[bin_col] == ref_bin]
    ref_conserved = ref_data['is_conserved'].sum()
    ref_not_conserved = len(ref_data) - ref_conserved

    for bin_idx in sorted(df[bin_col].dropna().unique()):
        bin_data = df[df[bin_col] == bin_idx]
        bin_conserved = bin_data['is_conserved'].sum()
        bin_not_conserved = len(bin_data) - bin_conserved

        # Create 2x2 contingency table
        # [[bin_conserved, bin_not_conserved], [ref_conserved, ref_not_conserved]]
        table = [[bin_conserved, bin_not_conserved],
                 [ref_conserved, ref_not_conserved]]

        # Fisher's exact test
        odds_ratio, p_value = fisher_exact(table)

        # Calculate 95% CI for log(OR) using Woolf's method
        # SE(log(OR)) = sqrt(1/a + 1/b + 1/c + 1/d)
        a, b, c, d = bin_conserved, bin_not_conserved, ref_conserved, ref_not_conserved
        if all(x > 0 for x in [a, b, c, d]):
            se_log_or = np.sqrt(1/a + 1/b + 1/c + 1/d)
            log_or = np.log(odds_ratio)
            ci_lower = np.exp(log_or - 1.96 * se_log_or)
            ci_upper = np.exp(log_or + 1.96 * se_log_or)
        else:
            ci_lower, ci_upper = np.nan, np.nan

        pct_conserved = 100 * bin_conserved / len(bin_data) if len(bin_data) > 0 else 0

        results.append({
            'bin': bin_idx,
            'n_total': len(bin_data),
            'n_conserved': bin_conserved,
            'pct_conserved': pct_conserved,
            'odds_ratio': odds_ratio,
            'or_ci_lower': ci_lower,
            'or_ci_upper': ci_upper,
            'p_value': p_value
        })

    return pd.DataFrame(results)

# Calculate for fixed bins (reference = neutral bin 5)
or_fixed = calculate_odds_ratio_vs_reference(sample_data, 'fixed_bin', ref_bin=5)
or_fixed['bin_label'] = [bin_labels[int(i)] for i in or_fixed['bin']]

print("Odds Ratios for Conserved Bases (phyloP >= 2.27) vs Neutral Bin:")
print(or_fixed[['bin_label', 'n_total', 'pct_conserved', 'odds_ratio', 'or_ci_lower', 'or_ci_upper', 'p_value']].to_string(index=False))


# %% Cell: ODDS RATIO VISUALIZATION
# =============================================================================
# ODDS RATIO VISUALIZATION
# =============================================================================

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Plot 1: Odds ratio forest plot
ax1 = axes[0]
y_pos = range(len(or_fixed))

# Plot odds ratios with confidence intervals
ax1.errorbar(or_fixed['odds_ratio'], y_pos,
             xerr=[or_fixed['odds_ratio'] - or_fixed['or_ci_lower'],
                   or_fixed['or_ci_upper'] - or_fixed['odds_ratio']],
             fmt='o', color='darkblue', capsize=3, markersize=8)

ax1.axvline(1, color='red', linestyle='--', linewidth=2, label='OR = 1 (no difference)')
ax1.set_yticks(y_pos)
ax1.set_yticklabels(or_fixed['bin_label'])
ax1.set_xlabel('Odds Ratio (vs Neutral Bin)')
ax1.set_ylabel('Predicted Skew Bin')
ax1.set_title(f'Odds of Conservation (phyloP >= {CONSERVED_THRESHOLD})')
ax1.set_xscale('log')
ax1.grid(axis='x', linestyle=':', alpha=0.6)
ax1.legend()

# Plot 2: Percent conserved by bin
ax2 = axes[1]
colors = ['indianred' if i < 5 else 'teal' if i > 5 else 'gray' for i in range(11)]
ax2.bar(range(len(or_fixed)), or_fixed['pct_conserved'], color=colors, alpha=0.7)
ax2.axhline(or_fixed[or_fixed['bin'] == 5]['pct_conserved'].values[0],
            color='gray', linestyle='--', linewidth=2, label='Neutral bin')
ax2.set_xticks(range(len(or_fixed)))
ax2.set_xticklabels(or_fixed['bin_label'], rotation=45, ha='right')
ax2.set_xlabel('Predicted Skew Bin')
ax2.set_ylabel(f'% Conserved (phyloP >= {CONSERVED_THRESHOLD})')
ax2.set_title('Percentage of Conserved Bases by Skew Bin')
ax2.legend()

for ax in axes:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.show()

# Summary
print(f"\nSummary:")
print(f"  Neutral bin % conserved: {or_fixed[or_fixed['bin']==5]['pct_conserved'].values[0]:.2f}%")
print(f"  Strong negative skew (<-1.5) OR: {or_fixed[or_fixed['bin']==0]['odds_ratio'].values[0]:.3f}")
print(f"  Strong positive skew (>1.5) OR: {or_fixed[or_fixed['bin']==10]['odds_ratio'].values[0]:.3f}")


# %% Cell: ODDS RATIO FOR ABSOLUTE SKEW BINS
# =============================================================================
# ODDS RATIO FOR ABSOLUTE SKEW QUANTILE BINS
# =============================================================================

# Calculate OR for absolute skew bins (reference = lowest |skew| bin)
or_abs = calculate_odds_ratio_vs_reference(sample_data, 'abs_quantile_bin', ref_bin=0)

# Create labels from boundaries
or_abs['bin_label'] = [f'{abs_quantile_boundaries.iloc[int(i)]:.2f}-{abs_quantile_boundaries.iloc[int(i)+1]:.2f}'
                        for i in or_abs['bin']]

print("\nOdds Ratios by Absolute Skew (vs Lowest |Skew| Bin):")
print(or_abs[['bin_label', 'n_total', 'pct_conserved', 'odds_ratio', 'p_value']].to_string(index=False))

# Plot
fig, ax = plt.subplots(figsize=(10, 6))

ax.errorbar(range(len(or_abs)), or_abs['odds_ratio'],
            yerr=[or_abs['odds_ratio'] - or_abs['or_ci_lower'],
                  or_abs['or_ci_upper'] - or_abs['odds_ratio']],
            fmt='o-', color='darkviolet', capsize=3, markersize=8, linewidth=2)

ax.axhline(1, color='red', linestyle='--', linewidth=2, label='OR = 1')
ax.set_xticks(range(len(or_abs)))
ax.set_xticklabels([f'{abs_quantile_boundaries.iloc[i]:.2f}' for i in range(len(or_abs))],
                   rotation=45, ha='right')
ax.set_xlabel('|Skew| (Quantile Bin Lower Bound)')
ax.set_ylabel('Odds Ratio (vs Lowest |Skew| Bin)')
ax.set_title(f'Odds of Conservation (phyloP >= {CONSERVED_THRESHOLD}) by |Skew|')
ax.grid(axis='y', linestyle=':', alpha=0.6)
ax.legend()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.show()

print(f"\nKey finding:")
print(f"  Variants with highest |skew| have OR = {or_abs['odds_ratio'].iloc[-1]:.3f} vs lowest |skew|")
print(f"  (95% CI: {or_abs['or_ci_lower'].iloc[-1]:.3f} - {or_abs['or_ci_upper'].iloc[-1]:.3f})")
