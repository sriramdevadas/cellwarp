#!/usr/bin/env python3
"""
Within-type variance vs Procrustes residual diagnostic.

Biology: If high within-type variance (noisy centroid) explains large Procrustes
residuals, then residual magnitude partly reflects centroid estimation error rather
than genuine cross-species program divergence. This diagnostic quantifies that risk.

Math: For each cell type in each species, compute:
    within_type_variance = (1/n) * sum_i ||x_i - centroid||^2
where x_i is a cell's expression vector in the full 16,959-gene ortholog space.
Then correlate (Spearman) the averaged variance with Procrustes residual magnitude.

Interpretation:
    rho > 0.6 → residuals partly reflect centroid noise — weaker claim
    rho < 0.3 → rigidity score measures genuine program divergence — stronger claim

Inputs:
    data/phase2_scaled/{human,mouse}_scaled.h5ad  (35 cell types, post-norm)
    output/phase2/scaled_35types/residuals_ranked.csv (Procrustes residual magnitudes)

Outputs:
    output/phase2/variance_diagnostic/variance_by_celltype.csv
    output/phase2/variance_diagnostic/variance_vs_residual.png
    output/phase2/variance_diagnostic/diagnostic_results.json
"""

import json
import os
import sys

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import spearmanr


def compute_within_type_variance(adata):
    """Compute mean squared distance from centroid for each cell type.

    For each cell type, computes ||x_i - mu||^2 for every cell i, then averages.
    Works in full gene space (all genes in adata.var).

    Uses the identity: mean(||x_i - mu||^2) = mean(||x_i||^2) - ||mu||^2
    to avoid materializing the full dense (cells x genes) deviation matrix.

    Returns:
        dict: cell_type -> within-type variance (scalar)
    """
    variances = {}
    for ct in sorted(adata.obs["cell_type"].unique()):
        mask = adata.obs["cell_type"] == ct
        X_ct = adata[mask].X  # sparse (n_cells x n_genes)

        n_cells = X_ct.shape[0]

        # Centroid (dense 1D)
        if sp.issparse(X_ct):
            centroid = np.asarray(X_ct.mean(axis=0)).ravel()
            # mean(||x_i||^2): element-wise square then mean
            mean_sq_norms = np.asarray(X_ct.power(2).mean(axis=0)).ravel().sum()
        else:
            centroid = X_ct.mean(axis=0)
            mean_sq_norms = np.mean(np.sum(X_ct ** 2, axis=1))

        centroid_sq_norm = np.sum(centroid ** 2)
        variance = mean_sq_norms - centroid_sq_norm

        variances[ct] = float(variance)
        print(f"  {ct}: n={n_cells}, within-type variance={variance:.2f}")

    return variances


def main():
    # Output directory
    out_dir = "output/phase2/variance_diagnostic"
    os.makedirs(out_dir, exist_ok=True)

    # Load data
    print("Loading human scaled data...")
    human = ad.read_h5ad("data/phase2_scaled/human_scaled.h5ad")
    print("Loading mouse scaled data...")
    mouse = ad.read_h5ad("data/phase2_scaled/mouse_scaled.h5ad")

    # Compute within-type variance
    print("\n--- Human within-type variance (16,959-gene space) ---")
    human_var = compute_within_type_variance(human)
    print("\n--- Mouse within-type variance (16,959-gene space) ---")
    mouse_var = compute_within_type_variance(mouse)

    # Build combined dataframe
    cell_types = sorted(human_var.keys())
    df_var = pd.DataFrame({
        "cell_type": cell_types,
        "human_variance": [human_var[ct] for ct in cell_types],
        "mouse_variance": [mouse_var[ct] for ct in cell_types],
    })
    df_var["mean_variance"] = (df_var["human_variance"] + df_var["mouse_variance"]) / 2

    # Check if human and mouse variances are substantially different
    rho_hm, p_hm = spearmanr(df_var["human_variance"], df_var["mouse_variance"])
    print(f"\nHuman vs mouse variance correlation: Spearman rho={rho_hm:.3f}, p={p_hm:.2e}")
    if rho_hm < 0.5:
        print("  WARNING: Human and mouse variances differ substantially — reporting separately")
        use_mean = False
    else:
        print("  Human and mouse variances are concordant — using mean for correlation")
        use_mean = True

    # Load Procrustes residual magnitudes
    residuals = pd.read_csv("output/phase2/scaled_35types/residuals_ranked.csv")
    df_merged = df_var.merge(residuals[["cell_type", "residual_magnitude"]], on="cell_type")

    # Spearman correlations
    var_col = "mean_variance" if use_mean else None
    results = {}

    # Always compute all three
    for label, col in [("mean", "mean_variance"), ("human", "human_variance"), ("mouse", "mouse_variance")]:
        rho, pval = spearmanr(df_merged[col], df_merged["residual_magnitude"])
        results[f"spearman_rho_{label}"] = float(rho)
        results[f"spearman_p_{label}"] = float(pval)

    primary_rho = results["spearman_rho_mean"]
    primary_p = results["spearman_p_mean"]

    # Interpretation
    if abs(primary_rho) > 0.6:
        interpretation = (
            f"rho={primary_rho:.3f} (|rho|>0.6): Within-type variance substantially correlates "
            "with Procrustes residual magnitude. Residuals partly reflect centroid estimation "
            "noise — the rigidity score claim should be tempered."
        )
    elif abs(primary_rho) > 0.3:
        interpretation = (
            f"rho={primary_rho:.3f} (0.3<|rho|<0.6): Moderate correlation. Within-type variance "
            "contributes somewhat to residual magnitude but does not dominate. Rigidity score "
            "reflects a mix of genuine divergence and centroid noise."
        )
    else:
        interpretation = (
            f"rho={primary_rho:.3f} (|rho|<0.3): No substantial correlation. Within-type variance "
            "does NOT explain Procrustes residual magnitude. Rigidity score measures genuine "
            "cross-species program divergence — stronger claim."
        )

    results["interpretation"] = interpretation
    results["human_mouse_variance_rho"] = float(rho_hm)
    results["human_mouse_variance_p"] = float(p_hm)
    results["n_cell_types"] = len(df_merged)

    # Print summary
    print("\n" + "=" * 70)
    print("WITHIN-TYPE VARIANCE DIAGNOSTIC RESULTS")
    print("=" * 70)
    print(f"Cell types analyzed: {len(df_merged)}")
    print(f"Human-mouse variance concordance: rho={rho_hm:.3f}, p={p_hm:.2e}")
    print(f"\nSpearman correlation (mean variance vs residual magnitude):")
    print(f"  rho = {primary_rho:.3f}, p = {primary_p:.4f}")
    print(f"\nSpearman correlation (human variance vs residual magnitude):")
    print(f"  rho = {results['spearman_rho_human']:.3f}, p = {results['spearman_p_human']:.4f}")
    print(f"\nSpearman correlation (mouse variance vs residual magnitude):")
    print(f"  rho = {results['spearman_rho_mouse']:.3f}, p = {results['spearman_p_mouse']:.4f}")
    print(f"\nInterpretation: {interpretation}")
    print("=" * 70)

    # Save CSV
    df_merged_sorted = df_merged.sort_values("residual_magnitude", ascending=False)
    df_merged_sorted.to_csv(f"{out_dir}/variance_by_celltype.csv", index=False)
    print(f"\nSaved: {out_dir}/variance_by_celltype.csv")

    # Save JSON
    with open(f"{out_dir}/diagnostic_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"Saved: {out_dir}/diagnostic_results.json")

    # Scatter plot
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.scatter(df_merged["mean_variance"], df_merged["residual_magnitude"],
               s=60, alpha=0.7, edgecolors="black", linewidth=0.5)

    # Label each point
    for _, row in df_merged.iterrows():
        label = row["cell_type"]
        # Shorten long names for readability
        if len(label) > 25:
            label = label[:22] + "..."
        ax.annotate(label, (row["mean_variance"], row["residual_magnitude"]),
                    fontsize=7, ha="left", va="bottom",
                    xytext=(4, 4), textcoords="offset points")

    ax.set_xlabel("Within-type variance (mean sq. dist. from centroid, 16,959-gene space)", fontsize=11)
    ax.set_ylabel("Procrustes residual magnitude", fontsize=11)
    ax.set_title(
        f"Within-type Variance vs Procrustes Residual\n"
        f"Spearman rho={primary_rho:.3f}, p={primary_p:.4f}  (n={len(df_merged)} cell types)",
        fontsize=12
    )

    # Add regression line for visual reference
    z = np.polyfit(df_merged["mean_variance"], df_merged["residual_magnitude"], 1)
    x_line = np.linspace(df_merged["mean_variance"].min(), df_merged["mean_variance"].max(), 100)
    ax.plot(x_line, np.polyval(z, x_line), "r--", alpha=0.5, label="Linear fit")
    ax.legend(fontsize=10)

    plt.tight_layout()
    fig.savefig(f"{out_dir}/variance_vs_residual.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_dir}/variance_vs_residual.png")

    # Print top-5 and bottom-5 for context
    print("\n--- Top 5 highest within-type variance ---")
    for _, row in df_merged_sorted.head(5).iterrows():
        print(f"  {row['cell_type']}: var={row['mean_variance']:.1f}, residual={row['residual_magnitude']:.2f}")
    print("\n--- Top 5 lowest within-type variance ---")
    for _, row in df_merged_sorted.tail(5).iterrows():
        print(f"  {row['cell_type']}: var={row['mean_variance']:.1f}, residual={row['residual_magnitude']:.2f}")

    # Text description of plot
    print("\n--- Plot description ---")
    print(f"Scatter plot with 35 labeled points. X-axis: within-type variance (mean of "
          f"human and mouse). Y-axis: Procrustes residual magnitude from 35-type analysis. "
          f"Red dashed regression line overlaid. Correlation is rho={primary_rho:.3f} "
          f"(p={primary_p:.4f}). ", end="")
    if abs(primary_rho) > 0.3:
        print("Points show a visible trend — higher variance cell types tend to have "
              "larger residuals.")
    else:
        print("Points are scattered with no clear trend — variance does not predict residual size.")


if __name__ == "__main__":
    main()
