#!/usr/bin/env python3
"""
Inter-individual variance vs Procrustes residual diagnostic.

Biology: Tests whether cross-species divergence (Procrustes residual) tracks
inter-individual variation within each species. If high-variance cell types also
show large residuals, the same genes may be "plastic" both within and between
species — individual variation and evolutionary variation share the same axis.
If independent, different molecular programs underpin each.

Math: For each cell type in each species:
    1. Compute per-donor centroid: mu_d = (1/n_d) * sum_i x_i for donor d
    2. Compute population centroid: mu = (1/D) * sum_d mu_d  (unweighted by donor)
    3. Inter-donor variance: V = (1/D) * sum_d ||mu_d - mu||^2
       This measures how much the cell type's average expression profile shifts
       from individual to individual, independent of within-individual cell noise.

    Then correlate V with Procrustes residual magnitude across 35 cell types (Spearman).

Interpretation:
    rho > 0.6  → Individual and evolutionary variation share the same axis.
                  Same genes are plastic within and between species.
    rho < 0.3  → Inter-individual and inter-species variation are independent.
                  Different genes drive each.
    Either result is interesting and publishable.

Inputs:
    data/phase2_scaled/{human,mouse}_scaled.h5ad  (35 cell types, post-norm)
    output/phase2/scaled_35types/residuals_ranked.csv (Procrustes residual magnitudes)

Outputs:
    output/phase2/diagnostics/interdonor_variance/
        interdonor_variance_by_celltype.csv
        interdonor_vs_residual.png
        diagnostic_results.json
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

MIN_CELLS_PER_DONOR = 10


def compute_interdonor_variance(adata):
    """Compute inter-donor variance for each cell type.

    For each cell type:
        1. Group cells by donor_id
        2. Keep donors with >= MIN_CELLS_PER_DONOR cells
        3. Compute per-donor centroid (mean expression)
        4. Compute inter-donor variance: mean squared distance of donor
           centroids from the population centroid (unweighted mean of
           donor centroids)

    Returns:
        dict: cell_type -> {
            'interdonor_variance': float,
            'n_donors': int,
            'n_donors_passing': int,
            'mean_cells_per_donor': float,
            'total_cells': int
        }
    """
    results = {}
    for ct in sorted(adata.obs["cell_type"].unique()):
        mask = adata.obs["cell_type"] == ct
        adata_ct = adata[mask]
        n_total = adata_ct.n_obs

        # Group by donor
        donor_counts = adata_ct.obs["donor_id"].value_counts()
        passing_donors = donor_counts[donor_counts >= MIN_CELLS_PER_DONOR].index.tolist()

        if len(passing_donors) < 2:
            print(f"  {ct}: SKIP — only {len(passing_donors)} donor(s) with "
                  f">={MIN_CELLS_PER_DONOR} cells (need ≥2)")
            continue

        # Compute per-donor centroids
        donor_centroids = []
        donor_cell_counts = []
        for donor in passing_donors:
            donor_mask = adata_ct.obs["donor_id"] == donor
            X_donor = adata_ct[donor_mask].X
            if sp.issparse(X_donor):
                centroid = np.asarray(X_donor.mean(axis=0)).ravel()
            else:
                centroid = np.mean(X_donor, axis=0)
            donor_centroids.append(centroid)
            donor_cell_counts.append(int(donor_mask.sum()))

        # Stack centroids: (n_donors x n_genes)
        centroids_matrix = np.array(donor_centroids)

        # Population centroid: unweighted mean of donor centroids
        pop_centroid = centroids_matrix.mean(axis=0)

        # Inter-donor variance: mean squared distance from population centroid
        deviations = centroids_matrix - pop_centroid  # (n_donors x n_genes)
        sq_distances = np.sum(deviations ** 2, axis=1)  # (n_donors,)
        interdonor_var = float(np.mean(sq_distances))

        results[ct] = {
            "interdonor_variance": interdonor_var,
            "n_donors": len(donor_counts),
            "n_donors_passing": len(passing_donors),
            "mean_cells_per_donor": float(np.mean(donor_cell_counts)),
            "total_cells": n_total,
        }

        print(f"  {ct}: n_donors={len(passing_donors)}/{len(donor_counts)}, "
              f"mean_cells/donor={np.mean(donor_cell_counts):.0f}, "
              f"inter-donor var={interdonor_var:.2f}")

    return results


def main():
    out_dir = "output/phase2/diagnostics/interdonor_variance"
    os.makedirs(out_dir, exist_ok=True)

    # Load data
    print("Loading human scaled data...")
    human = ad.read_h5ad("data/phase2_scaled/human_scaled.h5ad")
    print("Loading mouse scaled data...")
    mouse = ad.read_h5ad("data/phase2_scaled/mouse_scaled.h5ad")

    # Compute inter-donor variance
    print(f"\n--- Human inter-donor variance (min {MIN_CELLS_PER_DONOR} cells/donor) ---")
    human_results = compute_interdonor_variance(human)
    print(f"\n--- Mouse inter-donor variance (min {MIN_CELLS_PER_DONOR} cells/donor) ---")
    mouse_results = compute_interdonor_variance(mouse)

    # Find cell types present in both
    shared_types = sorted(set(human_results.keys()) & set(mouse_results.keys()))
    print(f"\nCell types with ≥2 qualifying donors in BOTH species: {len(shared_types)}")

    if len(shared_types) < 5:
        print("ERROR: Too few cell types with sufficient donor coverage. Cannot proceed.")
        sys.exit(1)

    # Build combined dataframe
    rows = []
    for ct in shared_types:
        rows.append({
            "cell_type": ct,
            "human_interdonor_var": human_results[ct]["interdonor_variance"],
            "mouse_interdonor_var": mouse_results[ct]["interdonor_variance"],
            "human_n_donors": human_results[ct]["n_donors_passing"],
            "mouse_n_donors": mouse_results[ct]["n_donors_passing"],
            "human_mean_cells_per_donor": human_results[ct]["mean_cells_per_donor"],
            "mouse_mean_cells_per_donor": mouse_results[ct]["mean_cells_per_donor"],
        })
    df = pd.DataFrame(rows)
    df["mean_interdonor_var"] = (df["human_interdonor_var"] + df["mouse_interdonor_var"]) / 2

    # Check concordance between human and mouse inter-donor variance
    rho_hm, p_hm = spearmanr(df["human_interdonor_var"], df["mouse_interdonor_var"])
    print(f"\nHuman vs mouse inter-donor variance concordance: "
          f"Spearman rho={rho_hm:.3f}, p={p_hm:.2e}")

    # Load Procrustes residual magnitudes
    residuals = pd.read_csv("output/phase2/scaled_35types/residuals_ranked.csv")
    df_merged = df.merge(residuals[["cell_type", "residual_magnitude"]], on="cell_type")
    print(f"Cell types with both inter-donor variance and residuals: {len(df_merged)}")

    # Spearman correlations
    results = {
        "n_cell_types": len(df_merged),
        "min_cells_per_donor": MIN_CELLS_PER_DONOR,
        "human_mouse_interdonor_var_concordance_rho": float(rho_hm),
        "human_mouse_interdonor_var_concordance_p": float(p_hm),
    }

    for label, col in [("mean", "mean_interdonor_var"),
                       ("human", "human_interdonor_var"),
                       ("mouse", "mouse_interdonor_var")]:
        rho, pval = spearmanr(df_merged[col], df_merged["residual_magnitude"])
        results[f"spearman_rho_{label}"] = float(rho)
        results[f"spearman_p_{label}"] = float(pval)

    primary_rho = results["spearman_rho_mean"]
    primary_p = results["spearman_p_mean"]

    # Interpretation
    if abs(primary_rho) > 0.6:
        interpretation = (
            f"rho={primary_rho:.3f} (|rho|>0.6): Inter-individual variance strongly correlates "
            "with Procrustes residual magnitude. Cell types with high donor-to-donor "
            "variation also diverge more cross-species. Individual variation and "
            "evolutionary variation share the same axis — suggests the same genes are "
            "plastic both within and between species."
        )
    elif abs(primary_rho) > 0.3:
        interpretation = (
            f"rho={primary_rho:.3f} (0.3<|rho|<0.6): Moderate correlation between "
            "inter-individual variance and cross-species divergence. Some overlap between "
            "the axes of individual variation and evolutionary divergence, but neither "
            "dominates. A partial shared plasticity model."
        )
    else:
        interpretation = (
            f"rho={primary_rho:.3f} (|rho|<0.3): No substantial correlation. Inter-individual "
            "and inter-species variation are independent. Different molecular programs "
            "underpin within-species donor variability vs cross-species divergence."
        )

    results["interpretation"] = interpretation

    # Per-species donor coverage summary
    results["human_total_donors"] = int(human.obs["donor_id"].nunique())
    results["mouse_total_donors"] = int(mouse.obs["donor_id"].nunique())
    results["human_median_qualifying_donors"] = float(df_merged["human_n_donors"].median())
    results["mouse_median_qualifying_donors"] = float(df_merged["mouse_n_donors"].median())

    # Print summary
    print("\n" + "=" * 70)
    print("INTER-DONOR VARIANCE vs PROCRUSTES RESIDUAL DIAGNOSTIC")
    print("=" * 70)
    print(f"Cell types analyzed: {len(df_merged)}")
    print(f"Donor filter: ≥{MIN_CELLS_PER_DONOR} cells per donor")
    print(f"Human donors: {results['human_total_donors']} total, "
          f"median {results['human_median_qualifying_donors']:.0f} qualifying per cell type")
    print(f"Mouse donors: {results['mouse_total_donors']} total, "
          f"median {results['mouse_median_qualifying_donors']:.0f} qualifying per cell type")
    print(f"\nHuman-mouse inter-donor var concordance: rho={rho_hm:.3f}, p={p_hm:.2e}")
    print(f"\nSpearman correlation (mean inter-donor var vs residual):")
    print(f"  rho = {primary_rho:.3f}, p = {primary_p:.4f}")
    print(f"Spearman correlation (human inter-donor var vs residual):")
    print(f"  rho = {results['spearman_rho_human']:.3f}, p = {results['spearman_p_human']:.4f}")
    print(f"Spearman correlation (mouse inter-donor var vs residual):")
    print(f"  rho = {results['spearman_rho_mouse']:.3f}, p = {results['spearman_p_mouse']:.4f}")
    print(f"\nInterpretation: {interpretation}")
    print("=" * 70)

    # Save CSV (sorted by residual descending)
    df_out = df_merged.sort_values("residual_magnitude", ascending=False)
    df_out.to_csv(f"{out_dir}/interdonor_variance_by_celltype.csv", index=False)
    print(f"\nSaved: {out_dir}/interdonor_variance_by_celltype.csv")

    # Save JSON
    with open(f"{out_dir}/diagnostic_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"Saved: {out_dir}/diagnostic_results.json")

    # --- Scatter plot: 2x2 layout ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Panel A: Mean inter-donor var vs residual
    ax = axes[0, 0]
    ax.scatter(df_merged["mean_interdonor_var"], df_merged["residual_magnitude"],
               s=60, alpha=0.7, edgecolors="black", linewidth=0.5, c="steelblue")
    for _, row in df_merged.iterrows():
        label = row["cell_type"]
        if len(label) > 22:
            label = label[:19] + "..."
        ax.annotate(label, (row["mean_interdonor_var"], row["residual_magnitude"]),
                    fontsize=6, ha="left", va="bottom",
                    xytext=(3, 3), textcoords="offset points")
    z = np.polyfit(df_merged["mean_interdonor_var"], df_merged["residual_magnitude"], 1)
    x_line = np.linspace(df_merged["mean_interdonor_var"].min(),
                         df_merged["mean_interdonor_var"].max(), 100)
    ax.plot(x_line, np.polyval(z, x_line), "r--", alpha=0.5)
    ax.set_xlabel("Mean inter-donor variance", fontsize=10)
    ax.set_ylabel("Procrustes residual magnitude", fontsize=10)
    ax.set_title(f"A) Mean inter-donor var vs residual\n"
                 f"Spearman rho={primary_rho:.3f}, p={primary_p:.4f}",
                 fontsize=11)

    # Panel B: Human inter-donor var vs residual
    ax = axes[0, 1]
    ax.scatter(df_merged["human_interdonor_var"], df_merged["residual_magnitude"],
               s=60, alpha=0.7, edgecolors="black", linewidth=0.5, c="coral")
    z = np.polyfit(df_merged["human_interdonor_var"], df_merged["residual_magnitude"], 1)
    x_line = np.linspace(df_merged["human_interdonor_var"].min(),
                         df_merged["human_interdonor_var"].max(), 100)
    ax.plot(x_line, np.polyval(z, x_line), "r--", alpha=0.5)
    ax.set_xlabel("Human inter-donor variance", fontsize=10)
    ax.set_ylabel("Procrustes residual magnitude", fontsize=10)
    ax.set_title(f"B) Human inter-donor var vs residual\n"
                 f"rho={results['spearman_rho_human']:.3f}, "
                 f"p={results['spearman_p_human']:.4f}", fontsize=11)

    # Panel C: Mouse inter-donor var vs residual
    ax = axes[1, 0]
    ax.scatter(df_merged["mouse_interdonor_var"], df_merged["residual_magnitude"],
               s=60, alpha=0.7, edgecolors="black", linewidth=0.5, c="mediumseagreen")
    z = np.polyfit(df_merged["mouse_interdonor_var"], df_merged["residual_magnitude"], 1)
    x_line = np.linspace(df_merged["mouse_interdonor_var"].min(),
                         df_merged["mouse_interdonor_var"].max(), 100)
    ax.plot(x_line, np.polyval(z, x_line), "r--", alpha=0.5)
    ax.set_xlabel("Mouse inter-donor variance", fontsize=10)
    ax.set_ylabel("Procrustes residual magnitude", fontsize=10)
    ax.set_title(f"C) Mouse inter-donor var vs residual\n"
                 f"rho={results['spearman_rho_mouse']:.3f}, "
                 f"p={results['spearman_p_mouse']:.4f}", fontsize=11)

    # Panel D: Human vs mouse inter-donor var concordance
    ax = axes[1, 1]
    ax.scatter(df_merged["human_interdonor_var"], df_merged["mouse_interdonor_var"],
               s=60, alpha=0.7, edgecolors="black", linewidth=0.5, c="mediumpurple")
    for _, row in df_merged.iterrows():
        label = row["cell_type"]
        if len(label) > 22:
            label = label[:19] + "..."
        ax.annotate(label, (row["human_interdonor_var"], row["mouse_interdonor_var"]),
                    fontsize=6, ha="left", va="bottom",
                    xytext=(3, 3), textcoords="offset points")
    z = np.polyfit(df_merged["human_interdonor_var"], df_merged["mouse_interdonor_var"], 1)
    x_line = np.linspace(df_merged["human_interdonor_var"].min(),
                         df_merged["human_interdonor_var"].max(), 100)
    ax.plot(x_line, np.polyval(z, x_line), "r--", alpha=0.5)
    # Diagonal reference line
    lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
            max(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lims, lims, "k:", alpha=0.3, label="y=x")
    ax.set_xlabel("Human inter-donor variance", fontsize=10)
    ax.set_ylabel("Mouse inter-donor variance", fontsize=10)
    ax.set_title(f"D) Human vs mouse inter-donor var\n"
                 f"rho={rho_hm:.3f}, p={p_hm:.2e}", fontsize=11)
    ax.legend(fontsize=8)

    plt.suptitle("Inter-individual Variance vs Procrustes Residual Diagnostic\n"
                 f"(n={len(df_merged)} cell types, donor filter ≥{MIN_CELLS_PER_DONOR} cells)",
                 fontsize=13, y=1.02)
    plt.tight_layout()
    fig.savefig(f"{out_dir}/interdonor_vs_residual.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_dir}/interdonor_vs_residual.png")

    # --- Print ranked tables ---
    print("\n--- Top 5 highest inter-donor variance (mean of species) ---")
    top5 = df_out.sort_values("mean_interdonor_var", ascending=False).head(5)
    for _, row in top5.iterrows():
        print(f"  {row['cell_type']}: inter-donor var={row['mean_interdonor_var']:.2f}, "
              f"residual={row['residual_magnitude']:.2f}")

    print("\n--- Top 5 lowest inter-donor variance (mean of species) ---")
    bot5 = df_out.sort_values("mean_interdonor_var", ascending=True).head(5)
    for _, row in bot5.iterrows():
        print(f"  {row['cell_type']}: inter-donor var={row['mean_interdonor_var']:.2f}, "
              f"residual={row['residual_magnitude']:.2f}")

    # --- Plot description ---
    print("\n--- Plot description ---")
    print(f"2x2 scatter plot panel. A) Mean inter-donor variance (avg of human and mouse) "
          f"vs Procrustes residual magnitude for {len(df_merged)} cell types. "
          f"B) Human-only inter-donor variance vs residual. "
          f"C) Mouse-only inter-donor variance vs residual. "
          f"D) Human vs mouse inter-donor variance concordance with y=x reference line. "
          f"Primary result: rho={primary_rho:.3f} (p={primary_p:.4f}). ", end="")
    if abs(primary_rho) > 0.6:
        print("Strong positive trend — cell types with high individual variation also "
              "diverge most between species.")
    elif abs(primary_rho) > 0.3:
        print("Moderate trend visible — partial overlap between individual and "
              "evolutionary plasticity axes.")
    else:
        print("No clear trend — inter-individual and inter-species variation appear "
              "to be governed by independent programs.")


if __name__ == "__main__":
    main()
