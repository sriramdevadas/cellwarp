"""
Reanalyze within-species negative controls at type counts matched
to cross-species comparisons, with bootstrap confidence intervals.

Cross-species reference points:
  - Primary (human-mouse): n=35 types, obs/null=0.522
  - Macaque: n=20 types
  - Mouse lemur: n=15 types

Biology: within-species tissue-pair Procrustes comparisons test whether
cell type geometry is preserved across tissues from the same atlas.
If the cross-species signal (obs/null=0.522) is distinguishable from
the within-species baseline, it supports evolutionary transformation
rather than batch artifact.

Math: obs/null ratio = observed Procrustes distance / median permutation
null distance. Lower values = more geometric coherence. We bootstrap
the mean obs/null across qualifying pairs to get a 95% CI.
"""

import json
import numpy as np
import pandas as pd

np.random.seed(42)

# --- Load data ---
within = pd.read_csv(
    "analysis/expanded_negative_controls/within_species_pairs.csv"
)

print(f"Total within-species pairs: {len(within)}")
print(f"  Human: {(within['species'] == 'human').sum()}")
print(f"  Mouse: {(within['species'] == 'mouse').sum()}")
print()

# --- Distribution of n_types ---
print("=== Cell type count distribution ===")
print(within['n_types'].describe())
print()
print("n_types value counts:")
print(within['n_types'].value_counts().sort_index())
print()

# Cross-species reference
CROSS_SPECIES_OBS_NULL = 0.522

# --- Check n >= 15 threshold ---
n15 = within[within['n_types'] >= 15]
print(f"Pairs with n_types >= 15: {len(n15)}")
print(f"Pairs with n_types >= 12: {len(within[within['n_types'] >= 12])}")
print(f"Pairs with n_types >= 10: {len(within[within['n_types'] >= 10])}")
print(f"Pairs with n_types >= 8:  {len(within[within['n_types'] >= 8])}")
print()

# --- Bootstrap CI function ---
def bootstrap_ci(data, n_boot=10000, ci=0.95, seed=42):
    """Compute bootstrap CI on the mean of data.

    Parameters
    ----------
    data : array-like
        Sample values.
    n_boot : int
        Number of bootstrap resamples.
    ci : float
        Confidence level (0-1).
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    dict with mean, median, ci_lower, ci_upper, n
    """
    rng = np.random.RandomState(seed)
    arr = np.array(data)
    n = len(arr)
    boot_means = np.array([
        arr[rng.randint(0, n, size=n)].mean()
        for _ in range(n_boot)
    ])
    alpha = (1 - ci) / 2
    ci_lower = np.percentile(boot_means, alpha * 100)
    ci_upper = np.percentile(boot_means, (1 - alpha) * 100)
    return {
        "n": int(n),
        "mean": float(arr.mean()),
        "median": float(np.median(arr)),
        "std": float(arr.std(ddof=1)) if n > 1 else 0.0,
        "ci_lower": float(ci_lower),
        "ci_upper": float(ci_upper),
        "boot_mean_of_means": float(boot_means.mean()),
    }


# --- Analyze at multiple thresholds ---
results = {}

thresholds = [15, 12, 10, 8, 6]
for thresh in thresholds:
    subset = within[within['n_types'] >= thresh]
    if len(subset) == 0:
        results[f"n_ge_{thresh}"] = {
            "threshold": thresh,
            "n_pairs": 0,
            "note": "No within-species pairs qualify at this threshold"
        }
        continue

    ratios = subset['obs_to_null_ratio'].values
    stats = bootstrap_ci(ratios, n_boot=10000)

    # Fraction with obs/null < cross-species 0.522
    frac_lower = float((ratios < CROSS_SPECIES_OBS_NULL).sum() / len(ratios))

    entry = {
        "threshold": thresh,
        "n_pairs": int(len(subset)),
        "mean_obs_null": stats["mean"],
        "median_obs_null": stats["median"],
        "std_obs_null": stats["std"],
        "ci_lower_95": stats["ci_lower"],
        "ci_upper_95": stats["ci_upper"],
        "cross_species_obs_null": CROSS_SPECIES_OBS_NULL,
        "fraction_lower_than_cross_species": frac_lower,
        "cross_species_within_ci": stats["ci_lower"] <= CROSS_SPECIES_OBS_NULL <= stats["ci_upper"],
        "pairs": subset[['pair_id', 'n_types', 'obs_to_null_ratio']].to_dict('records'),
    }
    results[f"n_ge_{thresh}"] = entry

    print(f"\n=== Threshold: n_types >= {thresh} ===")
    print(f"  Qualifying pairs: {len(subset)}")
    print(f"  Mean obs/null:   {stats['mean']:.4f}")
    print(f"  Median obs/null: {stats['median']:.4f}")
    print(f"  Std obs/null:    {stats['std']:.4f}")
    print(f"  95% Bootstrap CI on mean: [{stats['ci_lower']:.4f}, {stats['ci_upper']:.4f}]")
    print(f"  Cross-species obs/null:   {CROSS_SPECIES_OBS_NULL}")
    print(f"  Cross-species within CI:  {entry['cross_species_within_ci']}")
    print(f"  Fraction < cross-species: {frac_lower:.3f} ({int(frac_lower*len(subset))}/{len(subset)})")

# --- Key finding: largest available within-species pairs ---
print("\n\n=== KEY ANALYSIS: Largest within-species pairs (n >= 10) ===")
large_pairs = within[within['n_types'] >= 10].sort_values('n_types', ascending=False)
if len(large_pairs) > 0:
    for _, row in large_pairs.iterrows():
        print(f"  {row['pair_id']}: n={row['n_types']}, obs/null={row['obs_to_null_ratio']:.4f}, p={row['p_value']:.6f}")
else:
    print("  (none)")

# --- Effect size: all pairs vs cross-species ---
print("\n\n=== Effect size analysis (all 24 pairs) ===")
all_ratios = within['obs_to_null_ratio'].values
all_stats = bootstrap_ci(all_ratios, n_boot=10000)
print(f"  All pairs mean obs/null:  {all_stats['mean']:.4f}")
print(f"  All pairs median obs/null: {all_stats['median']:.4f}")
print(f"  95% Bootstrap CI on mean: [{all_stats['ci_lower']:.4f}, {all_stats['ci_upper']:.4f}]")
print(f"  Cross-species obs/null:    {CROSS_SPECIES_OBS_NULL}")
print(f"  Cross-species within CI:   {all_stats['ci_lower'] <= CROSS_SPECIES_OBS_NULL <= all_stats['ci_upper']}")
frac_all = (all_ratios < CROSS_SPECIES_OBS_NULL).sum() / len(all_ratios)
print(f"  Fraction < cross-species:  {frac_all:.3f} ({int(frac_all*len(all_ratios))}/{len(all_ratios)})")

# Cohen's d
pooled_std = all_stats['std']
if pooled_std > 0:
    cohens_d = (CROSS_SPECIES_OBS_NULL - all_stats['mean']) / pooled_std
    print(f"  Cohen's d (cross - within)/sd: {cohens_d:.3f}")

# --- Type-count correlation ---
from scipy import stats as sp_stats
corr, corr_p = sp_stats.spearmanr(within['n_types'], within['obs_to_null_ratio'])
print(f"\n  Spearman correlation (n_types vs obs/null): r={corr:.3f}, p={corr_p:.4f}")
print(f"  → {'Significant' if corr_p < 0.05 else 'Not significant'} relationship between type count and obs/null ratio")

# --- Build final results dict ---
final_results = {
    "analysis": "within_species_matched_control_reanalysis",
    "date": "2026-04-06",
    "cross_species_reference": {
        "primary_n_types": 35,
        "primary_obs_null": CROSS_SPECIES_OBS_NULL,
        "macaque_n_types": 20,
        "mouse_lemur_n_types": 15,
    },
    "within_species_summary": {
        "total_pairs": int(len(within)),
        "n_types_range": [int(within['n_types'].min()), int(within['n_types'].max())],
        "n_types_median": float(within['n_types'].median()),
        "pairs_at_n_ge_15": 0,
        "pairs_at_n_ge_12": int(len(within[within['n_types'] >= 12])),
        "pairs_at_n_ge_10": int(len(within[within['n_types'] >= 10])),
        "pairs_at_n_ge_8": int(len(within[within['n_types'] >= 8])),
    },
    "all_pairs_stats": {
        "n": int(len(all_ratios)),
        "mean_obs_null": float(all_stats['mean']),
        "median_obs_null": float(all_stats['median']),
        "std_obs_null": float(all_stats['std']),
        "ci_lower_95": float(all_stats['ci_lower']),
        "ci_upper_95": float(all_stats['ci_upper']),
        "cross_species_within_ci": bool(all_stats['ci_lower'] <= CROSS_SPECIES_OBS_NULL <= all_stats['ci_upper']),
        "fraction_lower_than_cross_species": float(frac_all),
        "cohens_d": float(cohens_d) if pooled_std > 0 else None,
        "spearman_n_types_vs_ratio": {"r": float(corr), "p": float(corr_p)},
    },
    "threshold_results": results,
    "conclusion": "",
}

# --- Determine conclusion ---
ci_contains_cross = all_stats['ci_lower'] <= CROSS_SPECIES_OBS_NULL <= all_stats['ci_upper']
if ci_contains_cross:
    final_results["conclusion"] = (
        f"The cross-species obs/null ratio ({CROSS_SPECIES_OBS_NULL}) falls WITHIN the 95% "
        f"bootstrap CI of the within-species mean [{all_stats['ci_lower']:.4f}, "
        f"{all_stats['ci_upper']:.4f}]. At the available type counts (n=6-12), "
        f"within-species and cross-species coherence are NOT clearly distinguishable. "
        f"However, no within-species pairs reach n>=15, so a direct type-count-matched "
        f"comparison is not possible. The confound is that within-species pairs have "
        f"fewer types (max 12) than cross-species (n=15-35), and type count may "
        f"influence the obs/null ratio."
    )
else:
    direction = "above" if CROSS_SPECIES_OBS_NULL > all_stats['ci_upper'] else "below"
    final_results["conclusion"] = (
        f"The cross-species obs/null ratio ({CROSS_SPECIES_OBS_NULL}) falls {direction} "
        f"the 95% bootstrap CI of the within-species mean [{all_stats['ci_lower']:.4f}, "
        f"{all_stats['ci_upper']:.4f}], indicating distinguishable signal. "
        f"However, no within-species pairs reach n>=15, so a direct type-count-matched "
        f"comparison is not possible."
    )

print(f"\n\n=== CONCLUSION ===")
print(final_results["conclusion"])

# --- Save ---
with open("analysis/within_species_matched/matched_control_results.json", "w") as f:
    json.dump(final_results, f, indent=2)
print("\nSaved: analysis/within_species_matched/matched_control_results.json")
