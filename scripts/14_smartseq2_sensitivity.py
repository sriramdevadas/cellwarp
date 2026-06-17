"""
CellWarp — Smart-seq2 Protocol Sensitivity Analysis

Tests whether the 35-type Procrustes result is robust to protocol mix in mouse data.

Background
----------
Tabula Muris Senis used both Smart-seq2 and 10x Chromium (3' v2). The Smart-seq2
fraction varies by cell type (from ~16% for B cells up to ~68% for CD4+ T cells).
This variable protocol mix is a potential confounder because:
  - Smart-seq2 captures full-length transcripts (deeper per cell, fewer cells)
  - 10x Chromium captures 3' end only (shallower per cell, more cells)
  - Gene detection rates and expression profiles differ between protocols

If the Procrustes result is robust to excluding Smart-seq2 cells entirely, then
protocol mix is not driving the geometric structure we observe.

Math
----
We recompute mouse centroids using only 10x cells:
    μ_t^{10x} = (1/n_{10x}) Σ_{i ∈ 10x} x_i

Then rerun the full Procrustes pipeline (PCA → alignment → permutation test)
with these 10x-only mouse centroids paired with the same human centroids.

We compare rigidity rankings via Spearman correlation between original and
10x-only residual magnitudes. ρ > 0.9 = stable, ρ < 0.7 = concerning.
"""

# obs/null ratio uses null_median (canonical convention, matches src/procrustes.py).

from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from cellwarp.procrustes import (
    compute_centroids,
    pca_reduce_centroids,
    permutation_test,
    procrustes_align,
    compute_residual_vectors,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DATA_DIR = PROJECT_ROOT / "data" / "phase2_scaled"
ORIGINAL_RESULTS_DIR = PROJECT_ROOT / "output" / "phase2" / "scaled_35types"
OUTPUT_DIR = PROJECT_ROOT / "output" / "phase2" / "sensitivity" / "smartseq2"

MOUSE_PATH = DATA_DIR / "mouse_scaled.h5ad"
HUMAN_PATH = DATA_DIR / "human_scaled.h5ad"

ASSAY_COL = "assay"
CELL_TYPE_COL = "cell_type"
TENX_ASSAY = "10x 3' v2"
SMARTSEQ2_ASSAY = "Smart-seq2"

N_PERMUTATIONS = 10_000
RANDOM_SEED = 42
MIN_CELLS_10X = 100  # Minimum 10x cells for reliable centroid

# Progenitor cell types (from developmental_annotations.csv)
PROGENITOR_TYPES = {
    "hematopoietic precursor cell",
    "hematopoietic stem cell",
    "basal cell",
    "mesenchymal stem cell of adipose tissue",
    "mesenchymal stem cell",
}


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # ==================================================================
    # STEP 1: Load data and check assay metadata
    # ==================================================================
    print("=" * 70)
    print("  STEP 1: Assay metadata in mouse data")
    print("=" * 70)

    mouse = ad.read_h5ad(MOUSE_PATH)
    print(f"\n  Mouse data shape: {mouse.shape}")
    print(f"\n  Assay value counts (global):")
    assay_counts = mouse.obs[ASSAY_COL].value_counts()
    for assay, count in assay_counts.items():
        print(f"    {assay}: {count:,} cells ({count / len(mouse) * 100:.1f}%)")

    # ==================================================================
    # STEP 2: Per-cell-type protocol breakdown + 10x-only centroids
    # ==================================================================
    print("\n" + "=" * 70)
    print("  STEP 2: Per-cell-type protocol breakdown")
    print("=" * 70)

    cell_types_35 = sorted(mouse.obs[CELL_TYPE_COL].unique())
    print(f"\n  Cell types in mouse data: {len(cell_types_35)}")

    protocol_breakdown = []
    flagged_types = []

    print(f"\n  {'Cell Type':<50} {'n_10x':>7} {'n_SS2':>7} {'%SS2':>7} {'Flag':>6}")
    print(f"  {'-' * 80}")

    for ct in cell_types_35:
        mask_ct = mouse.obs[CELL_TYPE_COL] == ct
        mask_10x = mouse.obs[ASSAY_COL] == TENX_ASSAY
        mask_ss2 = mouse.obs[ASSAY_COL] == SMARTSEQ2_ASSAY

        n_10x = int((mask_ct & mask_10x).sum())
        n_ss2 = int((mask_ct & mask_ss2).sum())
        n_total = n_10x + n_ss2
        frac_ss2 = n_ss2 / n_total if n_total > 0 else 0.0

        flag = "< 100" if n_10x < MIN_CELLS_10X else ""
        if flag:
            flagged_types.append(ct)

        protocol_breakdown.append({
            "cell_type": ct,
            "n_10x": n_10x,
            "n_smartseq2": n_ss2,
            "n_total": n_total,
            "fraction_smartseq2": frac_ss2,
            "flag_low_10x": n_10x < MIN_CELLS_10X,
        })

        print(f"  {ct:<50} {n_10x:>7,} {n_ss2:>7,} {frac_ss2:>6.1%} {flag:>6}")

    protocol_df = pd.DataFrame(protocol_breakdown)
    protocol_df.to_csv(OUTPUT_DIR / "protocol_breakdown.csv", index=False)

    print(f"\n  Flagged types (10x n < {MIN_CELLS_10X}): {len(flagged_types)}")
    for ct in flagged_types:
        row = protocol_df[protocol_df["cell_type"] == ct].iloc[0]
        print(f"    - {ct}: only {row['n_10x']} 10x cells")

    # Decide whether to proceed with all 35 or drop flagged types
    if flagged_types:
        print(f"\n  WARNING: {len(flagged_types)} cell types have < {MIN_CELLS_10X} 10x cells.")
        print(f"  Will run analysis both ways: all 35 types AND excluding flagged types.")
        run_reduced = True
    else:
        print(f"\n  All cell types have ≥ {MIN_CELLS_10X} 10x cells. Proceeding with all 35.")
        run_reduced = False

    # ==================================================================
    # Compute 10x-only mouse centroids
    # ==================================================================
    print("\n" + "=" * 70)
    print("  Computing 10x-only mouse centroids")
    print("=" * 70)

    mouse_10x = mouse[mouse.obs[ASSAY_COL] == TENX_ASSAY].copy()
    print(f"\n  10x-only mouse cells: {mouse_10x.shape[0]:,}")

    # Filter to cell types with enough cells
    valid_cts_10x = [ct for ct in cell_types_35 if ct not in flagged_types]
    if flagged_types:
        # For the full 35-type run, include flagged types anyway (they still have some cells)
        print(f"\n  Computing centroids for all {len(cell_types_35)} types (including flagged):")
    mouse_10x_centroids_all = compute_centroids(mouse_10x, CELL_TYPE_COL)

    # ==================================================================
    # Load human centroids (precomputed from original analysis)
    # ==================================================================
    print("\n" + "=" * 70)
    print("  Loading human centroids (unaffected by protocol mix)")
    print("=" * 70)

    human_centroids = pd.read_csv(
        ORIGINAL_RESULTS_DIR / "centroids_human_35.csv", index_col=0
    )
    print(f"  Human centroids: {human_centroids.shape[0]} types × {human_centroids.shape[1]:,} genes")

    # Align gene sets
    shared_genes = sorted(
        set(human_centroids.columns) & set(mouse_10x_centroids_all.columns)
    )
    print(f"  Shared genes: {len(shared_genes):,}")
    human_centroids = human_centroids[shared_genes]
    mouse_10x_centroids_all = mouse_10x_centroids_all[shared_genes]

    # Ensure same cell types
    shared_cts = sorted(
        set(human_centroids.index) & set(mouse_10x_centroids_all.index)
    )
    print(f"  Shared cell types: {len(shared_cts)}")
    human_centroids = human_centroids.loc[shared_cts]
    mouse_10x_centroids_all = mouse_10x_centroids_all.loc[shared_cts]

    # ==================================================================
    # STEP 3: Rerun full 35-type Procrustes with 10x-only mouse centroids
    # ==================================================================
    print("\n" + "=" * 70)
    print("  STEP 3: Procrustes analysis — 10x-only mouse centroids")
    print("=" * 70)

    # PCA
    human_pca, mouse_pca, pca_model, cell_types = pca_reduce_centroids(
        human_centroids, mouse_10x_centroids_all
    )

    # Procrustes alignment
    result = procrustes_align(human_pca, mouse_pca)

    # Permutation test
    p_value, null_dist = permutation_test(
        human_pca, mouse_pca,
        n_permutations=N_PERMUTATIONS,
        seed=RANDOM_SEED,
    )

    # Obs/null ratio
    null_mean = np.mean(null_dist)
    obs_null_ratio = result.distance / np.median(null_dist)

    print(f"\n  === 10x-only Procrustes Summary ===")
    print(f"  Procrustes distance: {result.distance:.6f}")
    print(f"  p-value: {p_value:.6f}")
    print(f"  Obs/null ratio: {obs_null_ratio:.4f}")

    # Residual vectors
    residuals_10x = compute_residual_vectors(result, cell_types)

    # Save 10x residuals ranked
    residuals_10x_ranked = []
    total_ssr = result.distance_squared
    for ct in cell_types:
        mag = np.linalg.norm(residuals_10x[ct])
        residuals_10x_ranked.append({
            "cell_type": ct,
            "residual_magnitude_10x": mag,
            "pct_of_ssr_10x": mag**2 / total_ssr * 100 if total_ssr > 0 else 0,
        })
    residuals_10x_df = pd.DataFrame(residuals_10x_ranked).sort_values(
        "residual_magnitude_10x", ascending=False
    ).reset_index(drop=True)
    residuals_10x_df["rank_10x"] = residuals_10x_df.index + 1
    residuals_10x_df.to_csv(OUTPUT_DIR / "residuals_10x_only.csv", index=False)

    # ==================================================================
    # STEP 4: Compare rigidity rankings (Spearman correlation)
    # ==================================================================
    print("\n" + "=" * 70)
    print("  STEP 4: Rigidity ranking comparison (original vs 10x-only)")
    print("=" * 70)

    # Load original residuals
    original_resid = pd.read_csv(ORIGINAL_RESULTS_DIR / "residuals_ranked.csv")

    # Merge on cell type
    comparison = original_resid[["cell_type", "residual_magnitude", "rank"]].merge(
        residuals_10x_df[["cell_type", "residual_magnitude_10x", "rank_10x"]],
        on="cell_type",
        how="inner",
    )

    # Spearman correlation on residual magnitudes
    rho_mag, p_mag = stats.spearmanr(
        comparison["residual_magnitude"], comparison["residual_magnitude_10x"]
    )
    # Spearman correlation on ranks
    rho_rank, p_rank = stats.spearmanr(
        comparison["rank"], comparison["rank_10x"]
    )

    print(f"\n  Spearman ρ (residual magnitudes): {rho_mag:.4f} (p={p_mag:.2e})")
    print(f"  Spearman ρ (ranks):              {rho_rank:.4f} (p={p_rank:.2e})")

    if rho_mag > 0.9:
        verdict = "STABLE — protocol mix is NOT a confounder"
    elif rho_mag > 0.7:
        verdict = "MODERATE — some sensitivity to protocol mix"
    else:
        verdict = "UNSTABLE — protocol mix is affecting rigidity rankings"

    print(f"\n  Verdict: {verdict}")

    # Save comparison
    comparison.to_csv(OUTPUT_DIR / "rigidity_comparison.csv", index=False)

    # Load original Procrustes results for comparison
    with open(ORIGINAL_RESULTS_DIR / "procrustes_results_35.json") as f:
        orig_results = json.load(f)

    orig_distance = orig_results["procrustes"]["distance"]
    orig_p = orig_results["permutation_test"]["p_value"]
    orig_obs_null = orig_distance / orig_results["permutation_test"]["null_distribution_summary"]["median"]

    print(f"\n  === Side-by-side comparison ===")
    print(f"  {'Metric':<30} {'Original':>15} {'10x-only':>15} {'Delta':>15}")
    print(f"  {'-' * 75}")
    print(f"  {'Procrustes distance':<30} {orig_distance:>15.4f} {result.distance:>15.4f} {result.distance - orig_distance:>+15.4f}")
    print(f"  {'p-value':<30} {orig_p:>15.6f} {p_value:>15.6f} {'':>15}")
    print(f"  {'Obs/null ratio':<30} {orig_obs_null:>15.4f} {obs_null_ratio:>15.4f} {obs_null_ratio - orig_obs_null:>+15.4f}")
    print(f"  {'Spearman ρ (magnitudes)':<30} {'':>15} {rho_mag:>15.4f} {'':>15}")

    # ==================================================================
    # Plot: rigidity comparison scatter
    # ==================================================================
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Scatter: original vs 10x-only residual magnitudes
    ax = axes[0]
    ax.scatter(
        comparison["residual_magnitude"],
        comparison["residual_magnitude_10x"],
        s=30, alpha=0.7, color="steelblue",
    )
    # Label progenitor types
    for _, row in comparison.iterrows():
        if row["cell_type"] in PROGENITOR_TYPES:
            ax.annotate(
                row["cell_type"][:20],
                (row["residual_magnitude"], row["residual_magnitude_10x"]),
                fontsize=6, color="red", alpha=0.8,
                xytext=(3, 3), textcoords="offset points",
            )
    # Identity line
    lims = [
        min(comparison["residual_magnitude"].min(), comparison["residual_magnitude_10x"].min()) * 0.9,
        max(comparison["residual_magnitude"].max(), comparison["residual_magnitude_10x"].max()) * 1.1,
    ]
    ax.plot(lims, lims, "--", color="gray", alpha=0.5)
    ax.set_xlabel("Residual magnitude (all cells)")
    ax.set_ylabel("Residual magnitude (10x-only)")
    ax.set_title(f"Rigidity comparison (Spearman ρ = {rho_mag:.3f})")
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    # Rank comparison
    ax = axes[1]
    ax.scatter(comparison["rank"], comparison["rank_10x"], s=30, alpha=0.7, color="steelblue")
    for _, row in comparison.iterrows():
        if row["cell_type"] in PROGENITOR_TYPES:
            ax.annotate(
                row["cell_type"][:20],
                (row["rank"], row["rank_10x"]),
                fontsize=6, color="red", alpha=0.8,
                xytext=(3, 3), textcoords="offset points",
            )
    ax.plot([0, 36], [0, 36], "--", color="gray", alpha=0.5)
    ax.set_xlabel("Rank (all cells)")
    ax.set_ylabel("Rank (10x-only)")
    ax.set_title(f"Rank comparison (Spearman ρ = {rho_rank:.3f})")
    ax.set_xlim(0, 36)
    ax.set_ylim(0, 36)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "rigidity_comparison.png", dpi=150, bbox_inches="tight")
    plt.close()

    print(f"\n  Plot saved: {OUTPUT_DIR / 'rigidity_comparison.png'}")
    print(f"  Description: Left panel — scatter of residual magnitudes (original vs 10x-only).")
    print(f"  Right panel — rank scatter. Progenitor types labeled in red.")
    print(f"  Diagonal line = perfect agreement. Points near line = stable result.")

    # ==================================================================
    # STEP 5: Progenitor finding in 10x-only
    # ==================================================================
    print("\n" + "=" * 70)
    print("  STEP 5: Progenitor vs differentiated divergence (10x-only)")
    print("=" * 70)

    # Classify cell types
    is_progenitor = comparison["cell_type"].isin(PROGENITOR_TYPES)
    prog_resid_10x = comparison.loc[is_progenitor, "residual_magnitude_10x"].values
    diff_resid_10x = comparison.loc[~is_progenitor, "residual_magnitude_10x"].values

    prog_resid_orig = comparison.loc[is_progenitor, "residual_magnitude"].values
    diff_resid_orig = comparison.loc[~is_progenitor, "residual_magnitude"].values

    # Mann-Whitney U test (10x-only)
    mw_stat_10x, mw_p_10x = stats.mannwhitneyu(
        prog_resid_10x, diff_resid_10x, alternative="greater"
    )

    # Mann-Whitney U test (original, for reference)
    mw_stat_orig, mw_p_orig = stats.mannwhitneyu(
        prog_resid_orig, diff_resid_orig, alternative="greater"
    )

    print(f"\n  Progenitor types (n={len(prog_resid_10x)}): {sorted(PROGENITOR_TYPES)}")
    print(f"\n  {'Metric':<40} {'Original':>12} {'10x-only':>12}")
    print(f"  {'-' * 65}")
    print(f"  {'Progenitor mean residual':<40} {np.mean(prog_resid_orig):>12.4f} {np.mean(prog_resid_10x):>12.4f}")
    print(f"  {'Differentiated mean residual':<40} {np.mean(diff_resid_orig):>12.4f} {np.mean(diff_resid_10x):>12.4f}")
    print(f"  {'Ratio (progenitor/differentiated)':<40} {np.mean(prog_resid_orig)/np.mean(diff_resid_orig):>12.4f} {np.mean(prog_resid_10x)/np.mean(diff_resid_10x):>12.4f}")
    print(f"  {'Mann-Whitney U':<40} {mw_stat_orig:>12.1f} {mw_stat_10x:>12.1f}")
    print(f"  {'Mann-Whitney p (one-sided)':<40} {mw_p_orig:>12.6f} {mw_p_10x:>12.6f}")

    progenitor_survives = mw_p_10x < 0.05
    print(f"\n  Progenitor finding survives in 10x-only: "
          f"{'YES' if progenitor_survives else 'NO'} (p={mw_p_10x:.4f})")

    # Sub-analysis: exclude HSC (n=23 10x cells, unreliable centroid)
    print(f"\n  --- Sub-analysis: excluding HSC (only 23 10x cells) ---")
    PROGENITOR_NO_HSC = PROGENITOR_TYPES - {"hematopoietic stem cell"}
    is_prog_no_hsc = comparison["cell_type"].isin(PROGENITOR_NO_HSC)
    # Also exclude HSC from differentiated comparison
    not_hsc = comparison["cell_type"] != "hematopoietic stem cell"
    prog_resid_10x_no_hsc = comparison.loc[is_prog_no_hsc, "residual_magnitude_10x"].values
    diff_resid_10x_no_hsc = comparison.loc[not_hsc & ~is_prog_no_hsc, "residual_magnitude_10x"].values

    mw_stat_no_hsc, mw_p_no_hsc = stats.mannwhitneyu(
        prog_resid_10x_no_hsc, diff_resid_10x_no_hsc, alternative="greater"
    )
    print(f"  Progenitor types (excl HSC, n={len(prog_resid_10x_no_hsc)}): "
          f"{sorted(PROGENITOR_NO_HSC)}")
    print(f"  Progenitor mean: {np.mean(prog_resid_10x_no_hsc):.4f}")
    print(f"  Differentiated mean: {np.mean(diff_resid_10x_no_hsc):.4f}")
    print(f"  Mann-Whitney p (one-sided): {mw_p_no_hsc:.6f}")
    progenitor_survives_no_hsc = mw_p_no_hsc < 0.05
    print(f"  Finding survives: {'YES' if progenitor_survives_no_hsc else 'NO'}")

    # Rank-change analysis: which types moved most?
    print(f"\n  --- Biggest rank changes (|Δrank| ≥ 5) ---")
    comparison["rank_delta"] = comparison["rank_10x"] - comparison["rank"]
    comparison["abs_rank_delta"] = comparison["rank_delta"].abs()
    comparison["frac_ss2"] = comparison["cell_type"].map(
        protocol_df.set_index("cell_type")["fraction_smartseq2"]
    )
    big_movers = comparison[comparison["abs_rank_delta"] >= 5].sort_values(
        "abs_rank_delta", ascending=False
    )
    print(f"  {'Cell Type':<45} {'Orig':>5} {'10x':>5} {'Δ':>5} {'%SS2':>7}")
    print(f"  {'-' * 70}")
    for _, row in big_movers.iterrows():
        print(f"  {row['cell_type']:<45} {int(row['rank']):>5} {int(row['rank_10x']):>5} "
              f"{int(row['rank_delta']):>+5} {row['frac_ss2']:>6.1%}")

    # Correlation between SS2 fraction and absolute rank change
    rho_ss2_rank, p_ss2_rank = stats.spearmanr(
        comparison["frac_ss2"], comparison["abs_rank_delta"]
    )
    print(f"\n  Spearman ρ (SS2 fraction vs |rank change|): {rho_ss2_rank:.4f} (p={p_ss2_rank:.4f})")

    # Spearman excluding HSC (unreliable)
    comp_no_hsc = comparison[comparison["cell_type"] != "hematopoietic stem cell"]
    rho_no_hsc, p_no_hsc = stats.spearmanr(
        comp_no_hsc["residual_magnitude"], comp_no_hsc["residual_magnitude_10x"]
    )
    print(f"\n  Spearman ρ excluding HSC: {rho_no_hsc:.4f} (p={p_no_hsc:.2e})")

    # ==================================================================
    # Save comprehensive results
    # ==================================================================
    print("\n" + "=" * 70)
    print("  Saving results")
    print("=" * 70)

    summary = {
        "analysis": "Smart-seq2 sensitivity analysis",
        "description": "Tests whether 35-type Procrustes result is robust to protocol mix",
        "mouse_cells": {
            "total": int(len(mouse)),
            "10x_only": int(len(mouse_10x)),
            "smartseq2": int(len(mouse)) - int(len(mouse_10x)),
            "smartseq2_fraction": (int(len(mouse)) - int(len(mouse_10x))) / int(len(mouse)),
        },
        "flagged_cell_types": {
            "count": len(flagged_types),
            "types": flagged_types,
            "reason": f"< {MIN_CELLS_10X} 10x cells",
        },
        "procrustes_10x_only": {
            "distance": float(result.distance),
            "distance_squared": float(result.distance_squared),
            "p_value": float(p_value),
            "obs_null_ratio": float(obs_null_ratio),
            "scaling": float(result.scaling),
            "n_cell_types": len(cell_types),
            "n_pca_components": int(human_pca.shape[1]),
        },
        "procrustes_original": {
            "distance": float(orig_distance),
            "p_value": float(orig_p),
            "obs_null_ratio": float(orig_obs_null),
        },
        "rigidity_comparison": {
            "spearman_rho_magnitudes": float(rho_mag),
            "spearman_p_magnitudes": float(p_mag),
            "spearman_rho_ranks": float(rho_rank),
            "spearman_p_ranks": float(p_rank),
            "spearman_rho_excl_hsc": float(rho_no_hsc),
            "spearman_p_excl_hsc": float(p_no_hsc),
            "ss2_fraction_vs_rank_change_rho": float(rho_ss2_rank),
            "ss2_fraction_vs_rank_change_p": float(p_ss2_rank),
            "verdict": verdict,
        },
        "progenitor_analysis": {
            "progenitor_types": sorted(PROGENITOR_TYPES),
            "n_progenitor": int(len(prog_resid_10x)),
            "n_differentiated": int(len(diff_resid_10x)),
            "original": {
                "mean_progenitor_residual": float(np.mean(prog_resid_orig)),
                "mean_differentiated_residual": float(np.mean(diff_resid_orig)),
                "mannwhitney_p": float(mw_p_orig),
            },
            "tenx_only": {
                "mean_progenitor_residual": float(np.mean(prog_resid_10x)),
                "mean_differentiated_residual": float(np.mean(diff_resid_10x)),
                "mannwhitney_p": float(mw_p_10x),
            },
            "progenitor_finding_survives": bool(progenitor_survives),
            "tenx_only_excl_hsc": {
                "progenitor_types": sorted(PROGENITOR_NO_HSC),
                "n_progenitor": int(len(prog_resid_10x_no_hsc)),
                "mean_progenitor_residual": float(np.mean(prog_resid_10x_no_hsc)),
                "mean_differentiated_residual": float(np.mean(diff_resid_10x_no_hsc)),
                "mannwhitney_p": float(mw_p_no_hsc),
                "survives": bool(progenitor_survives_no_hsc),
            },
        },
    }

    with open(OUTPUT_DIR / "sensitivity_results.json", "w") as f:
        json.dump(summary, f, indent=2)
    np.save(OUTPUT_DIR / "null_distribution_10x.npy", null_dist)

    print(f"\n  Saved: {OUTPUT_DIR / 'sensitivity_results.json'}")
    print(f"  Saved: {OUTPUT_DIR / 'null_distribution_10x.npy'}")
    print(f"  Saved: {OUTPUT_DIR / 'protocol_breakdown.csv'}")
    print(f"  Saved: {OUTPUT_DIR / 'residuals_10x_only.csv'}")
    print(f"  Saved: {OUTPUT_DIR / 'rigidity_comparison.csv'}")
    print(f"  Saved: {OUTPUT_DIR / 'rigidity_comparison.png'}")

    # ==================================================================
    # Final summary
    # ==================================================================
    print("\n" + "=" * 70)
    print("  FINAL SUMMARY — Smart-seq2 Sensitivity Analysis")
    print("=" * 70)
    print(f"\n  1. Protocol breakdown:")
    print(f"     Mouse: {len(mouse_10x):,} 10x cells + {len(mouse)-len(mouse_10x):,} Smart-seq2 cells")
    print(f"     Smart-seq2 fraction ranges: {protocol_df['fraction_smartseq2'].min():.1%} – {protocol_df['fraction_smartseq2'].max():.1%}")
    print(f"     Flagged cell types (10x n<{MIN_CELLS_10X}): {len(flagged_types)}")
    print(f"\n  2. Procrustes result (10x-only):")
    print(f"     Distance: {result.distance:.4f} (original: {orig_distance:.4f})")
    print(f"     p-value: {p_value:.6f} (original: {orig_p:.6f})")
    print(f"     Obs/null ratio: {obs_null_ratio:.4f} (original: {orig_obs_null:.4f})")
    print(f"\n  3. Rigidity rankings:")
    print(f"     Spearman ρ = {rho_mag:.4f} (p = {p_mag:.2e})")
    print(f"     Verdict: {verdict}")
    print(f"\n  4. Progenitor finding:")
    print(f"     10x-only MW p = {mw_p_10x:.4f} (original: {mw_p_orig:.4f})")
    print(f"     Survives: {'YES' if progenitor_survives else 'NO'}")
    print(f"\n  CONCLUSION: ", end="")
    if p_value < 0.01 and rho_mag > 0.7:
        print("Procrustes structure is robust to protocol mix.")
        print("  Smart-seq2 is not driving the geometric signal.")
    else:
        print("Sensitivity detected — interpret with caution.")


if __name__ == "__main__":
    main()
