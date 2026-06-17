#!/usr/bin/env python3
"""
Diagnostic: Expression Level vs Procrustes Rigidity

Tests whether cell-type rigidity (Procrustes residual magnitude) is a
scale artifact of expression level.

Biology: If high-expressing cell types systematically show larger residuals
simply because they have more "room" for deformation, then our rigidity
ranking is confounded and needs normalization.

Math: Computes Spearman rank correlation between several expression-level
metrics and Procrustes residual magnitude across 35 cell types.

Outputs:
  - Console table of all correlations
  - Scatter plot for strongest correlation
  - Normalized residual ranking comparison
"""

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import os

# -------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------
DATA_DIR = "output/phase2/scaled_35types"
OUT_DIR = "output/phase2/diagnostics/expression_level_vs_rigidity"
os.makedirs(OUT_DIR, exist_ok=True)

# -------------------------------------------------------------------
# 1. Load data
# -------------------------------------------------------------------
print("=" * 70)
print("DIAGNOSTIC: Expression Level vs Procrustes Rigidity")
print("=" * 70)

centroids_h = pd.read_csv(f"{DATA_DIR}/centroids_human_35.csv", index_col=0)
centroids_m = pd.read_csv(f"{DATA_DIR}/centroids_mouse_35.csv", index_col=0)
residuals = pd.read_csv(f"{DATA_DIR}/residuals_ranked.csv")

print(f"Human centroids: {centroids_h.shape}")
print(f"Mouse centroids: {centroids_m.shape}")
print(f"Residuals: {residuals.shape[0]} cell types")
print()

# Align cell types
cell_types = residuals["cell_type"].values
assert set(cell_types) == set(centroids_h.index), "Cell type mismatch human"
assert set(cell_types) == set(centroids_m.index), "Cell type mismatch mouse"

# -------------------------------------------------------------------
# 2. Compute expression-level metrics per cell type
# -------------------------------------------------------------------
print("-" * 70)
print("Step 1: Expression-level metrics per cell type")
print("-" * 70)

metrics = pd.DataFrame({"cell_type": cell_types})
metrics = metrics.set_index("cell_type")

# Residual magnitude
metrics["residual"] = residuals.set_index("cell_type").loc[cell_types, "residual_magnitude"].values

for species, centroids, prefix in [
    ("human", centroids_h, "h"),
    ("mouse", centroids_m, "m"),
]:
    vals = centroids.loc[cell_types].values  # (35, 16959)

    # Mean expression across all genes
    metrics[f"{prefix}_mean_expr"] = vals.mean(axis=1)

    # Median expression
    metrics[f"{prefix}_median_expr"] = np.median(vals, axis=1)

    # Detection rate: fraction of genes > 0.1
    metrics[f"{prefix}_detect_rate"] = (vals > 0.1).mean(axis=1)

    # Top-gene concentration: fraction of total expression in top 100 genes
    row_sums = vals.sum(axis=1, keepdims=True)
    # For each row, sort descending and take top 100
    sorted_vals = np.sort(vals, axis=1)[:, ::-1]
    top100_sum = sorted_vals[:, :100].sum(axis=1)
    metrics[f"{prefix}_top100_frac"] = top100_sum / row_sums.flatten()

# Cross-species averages
metrics["avg_mean_expr"] = (metrics["h_mean_expr"] + metrics["m_mean_expr"]) / 2
metrics["avg_detect_rate"] = (metrics["h_detect_rate"] + metrics["m_detect_rate"]) / 2
metrics["avg_top100_frac"] = (metrics["h_top100_frac"] + metrics["m_top100_frac"]) / 2

# Print summary table
print(f"\n{'Cell Type':<45} {'Resid':>7} {'H_mean':>8} {'M_mean':>8} {'H_det':>7} {'M_det':>7} {'H_top100':>8} {'M_top100':>8}")
print("-" * 120)
for ct in cell_types:
    r = metrics.loc[ct]
    print(f"{ct:<45} {r['residual']:7.2f} {r['h_mean_expr']:8.4f} {r['m_mean_expr']:8.4f} "
          f"{r['h_detect_rate']:7.3f} {r['m_detect_rate']:7.3f} "
          f"{r['h_top100_frac']:8.3f} {r['m_top100_frac']:8.3f}")

# -------------------------------------------------------------------
# 3. Spearman correlations with residual magnitude
# -------------------------------------------------------------------
print("\n" + "-" * 70)
print("Step 2: Spearman correlations with residual magnitude")
print("-" * 70)

corr_pairs = [
    ("Human mean expression", "h_mean_expr"),
    ("Mouse mean expression", "m_mean_expr"),
    ("Avg (H+M) mean expression", "avg_mean_expr"),
    ("Human median expression", "h_median_expr"),
    ("Mouse median expression", "m_median_expr"),
    ("Human detection rate (>0.1)", "h_detect_rate"),
    ("Mouse detection rate (>0.1)", "m_detect_rate"),
    ("Avg detection rate", "avg_detect_rate"),
    ("Human top-100 gene concentration", "h_top100_frac"),
    ("Mouse top-100 gene concentration", "m_top100_frac"),
    ("Avg top-100 gene concentration", "avg_top100_frac"),
]

results = []
print(f"\n{'Metric':<40} {'Spearman ρ':>12} {'p-value':>12} {'Signif':>8}")
print("-" * 75)

for label, col in corr_pairs:
    rho, pval = stats.spearmanr(metrics["residual"], metrics[col])
    flag = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else ""
    print(f"{label:<40} {rho:12.4f} {pval:12.6f} {flag:>8}")
    results.append({"metric": label, "column": col, "rho": rho, "pval": pval})

results_df = pd.DataFrame(results)

# Find strongest correlation (by absolute rho)
strongest = results_df.loc[results_df["rho"].abs().idxmax()]
print(f"\nStrongest correlation: {strongest['metric']}")
print(f"  Spearman ρ = {strongest['rho']:.4f}, p = {strongest['pval']:.6f}")

# -------------------------------------------------------------------
# 4. Plot scatter for strongest correlation
# -------------------------------------------------------------------
print("\n" + "-" * 70)
print("Step 3: Scatter plot for strongest correlation")
print("-" * 70)

fig, ax = plt.subplots(figsize=(10, 8))
x = metrics[strongest["column"]]
y = metrics["residual"]

ax.scatter(x, y, s=60, alpha=0.7, edgecolors="black", linewidths=0.5, c="steelblue")

# Label points
for ct in cell_types:
    ax.annotate(ct, (metrics.loc[ct, strongest["column"]], metrics.loc[ct, "residual"]),
                fontsize=6, alpha=0.7, ha="left", va="bottom",
                xytext=(3, 3), textcoords="offset points")

# Regression line
slope, intercept, r_val, p_val, std_err = stats.linregress(x, y)
x_line = np.linspace(x.min(), x.max(), 100)
ax.plot(x_line, slope * x_line + intercept, "r--", alpha=0.5, linewidth=1.5)

ax.set_xlabel(strongest["metric"], fontsize=12)
ax.set_ylabel("Procrustes Residual Magnitude", fontsize=12)
ax.set_title(f"Expression Level vs Rigidity Diagnostic\n"
             f"Spearman ρ = {strongest['rho']:.3f}, p = {strongest['pval']:.4f}",
             fontsize=13)
ax.grid(True, alpha=0.3)

plt.tight_layout()
scatter_path = f"{OUT_DIR}/strongest_correlation_scatter.png"
fig.savefig(scatter_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {scatter_path}")
print(f"Plot shows: Scatter of {strongest['metric']} (x) vs Procrustes residual (y)")
print(f"  for 35 cell types with labeled points and regression line.")

# Also save a multi-panel figure with top 4 correlations
top4 = results_df.reindex(results_df["rho"].abs().sort_values(ascending=False).index).head(4)

fig, axes = plt.subplots(2, 2, figsize=(14, 12))
for idx, (_, row) in enumerate(top4.iterrows()):
    ax = axes[idx // 2, idx % 2]
    x = metrics[row["column"]]
    y = metrics["residual"]
    ax.scatter(x, y, s=40, alpha=0.7, edgecolors="black", linewidths=0.5, c="steelblue")
    slope, intercept, _, _, _ = stats.linregress(x, y)
    x_line = np.linspace(x.min(), x.max(), 100)
    ax.plot(x_line, slope * x_line + intercept, "r--", alpha=0.5)
    ax.set_xlabel(row["metric"], fontsize=10)
    ax.set_ylabel("Residual Magnitude", fontsize=10)
    ax.set_title(f"ρ = {row['rho']:.3f}, p = {row['pval']:.4f}", fontsize=11)
    ax.grid(True, alpha=0.3)

plt.suptitle("Expression Level vs Rigidity — Top 4 Correlations", fontsize=14, y=1.01)
plt.tight_layout()
multi_path = f"{OUT_DIR}/top4_correlations.png"
fig.savefig(multi_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {multi_path}")

# -------------------------------------------------------------------
# 5. DECISION-064 normalization: residual / mean expression
# -------------------------------------------------------------------
print("\n" + "-" * 70)
print("Step 4: DECISION-064 normalized residual ranking")
print("-" * 70)

metrics["norm_residual_by_avg_mean"] = metrics["residual"] / metrics["avg_mean_expr"]

# Rank both
metrics["rank_raw"] = metrics["residual"].rank(ascending=False)
metrics["rank_normalized"] = metrics["norm_residual_by_avg_mean"].rank(ascending=False)

# Spearman between raw and normalized rankings
rho_ranks, p_ranks = stats.spearmanr(metrics["rank_raw"], metrics["rank_normalized"])
print(f"\nSpearman ρ between raw and normalized rankings: {rho_ranks:.4f} (p = {p_ranks:.6f})")

if rho_ranks > 0.9:
    print("→ Rankings are highly stable — normalization barely changes ordering")
elif rho_ranks > 0.7:
    print("→ Rankings are moderately stable — some reordering but broad pattern holds")
else:
    print("→ Rankings change substantially — expression level is a major confound!")

# Show side-by-side ranking
print(f"\n{'Cell Type':<45} {'Raw Rank':>10} {'Norm Rank':>10} {'Δ Rank':>8} {'Resid':>8} {'Norm Resid':>10}")
print("-" * 95)
sorted_metrics = metrics.sort_values("rank_raw")
for ct in sorted_metrics.index:
    r = sorted_metrics.loc[ct]
    delta = int(r["rank_raw"] - r["rank_normalized"])
    print(f"{ct:<45} {int(r['rank_raw']):>10} {int(r['rank_normalized']):>10} {delta:>+8} "
          f"{r['residual']:8.2f} {r['norm_residual_by_avg_mean']:10.4f}")

# Top 5 biggest rank changes
print("\nLargest rank changes after normalization:")
metrics["rank_delta"] = (metrics["rank_raw"] - metrics["rank_normalized"]).abs()
biggest_changes = metrics.sort_values("rank_delta", ascending=False).head(5)
for ct in biggest_changes.index:
    r = biggest_changes.loc[ct]
    print(f"  {ct}: raw rank {int(r['rank_raw'])} → normalized rank {int(r['rank_normalized'])} "
          f"(Δ = {int(r['rank_raw'] - r['rank_normalized']):+d})")

# -------------------------------------------------------------------
# 6. Save results
# -------------------------------------------------------------------
print("\n" + "-" * 70)
print("Step 5: Saving results")
print("-" * 70)

# Save correlations table
results_df.to_csv(f"{OUT_DIR}/correlations.csv", index=False)
print(f"Saved: {OUT_DIR}/correlations.csv")

# Save metrics table
metrics.to_csv(f"{OUT_DIR}/expression_metrics.csv")
print(f"Saved: {OUT_DIR}/expression_metrics.csv")

# Save normalized rankings
norm_ranking = metrics[["residual", "rank_raw", "norm_residual_by_avg_mean",
                         "rank_normalized", "rank_delta", "avg_mean_expr"]].sort_values("rank_raw")
norm_ranking.to_csv(f"{OUT_DIR}/normalized_rankings.csv")
print(f"Saved: {OUT_DIR}/normalized_rankings.csv")

# -------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

any_strong = (results_df["rho"].abs() > 0.5).any()
print(f"\nAny correlation |ρ| > 0.5? {'YES — potential confound' if any_strong else 'NO — expression level is not a confound'}")

if any_strong:
    strong = results_df[results_df["rho"].abs() > 0.5]
    print("  Strong correlations found:")
    for _, row in strong.iterrows():
        print(f"    {row['metric']}: ρ = {row['rho']:.4f}, p = {row['pval']:.6f}")

print(f"\nRaw vs normalized ranking stability: ρ = {rho_ranks:.4f}")
print(f"Interpretation: {'Rankings stable' if rho_ranks > 0.9 else 'Rankings changed' if rho_ranks < 0.7 else 'Moderate stability'}")
print("=" * 70)
