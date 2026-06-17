#!/usr/bin/env python3
"""
PCA Dimensionality Sensitivity Analysis (T1-C from validation agenda)

Biology
-------
Our Procrustes pipeline projects 35 cell-type centroids into a PCA subspace
before alignment. The default threshold (95% variance → 33 components)
is arbitrary. If the geometric signal is robust, the Procrustes distance
should remain significant and the per-cell-type residual ranking should
be stable across a range of dimensionalities.

Math
----
For each PCA cutoff k ∈ {10, 20, 33, 50}:
    1. Fit PCA(n_components=k) on the 70 combined centroids (35 human + 35 mouse).
    2. Run OPA alignment (mouse → human) in k-dimensional space.
    3. Run 10,000-iteration permutation test.
    4. Extract per-cell-type residual magnitudes ‖r_i‖.
    5. Compute Spearman ρ of residual rankings vs the 33-component baseline.

Pass criterion: p < 0.01 AND ranking ρ ≥ 0.70 at ALL cutoffs.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.decomposition import PCA

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from cellwarp.procrustes import (
    RANDOM_SEED,
    _procrustes_distance,
    procrustes_align,
    compute_residual_vectors,
)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

INPUT_DIR = PROJECT_ROOT / "output" / "phase2" / "scaled_35types"
OUTPUT_DIR = PROJECT_ROOT / "output" / "phase2" / "diagnostics" / "pca_sensitivity"

PCA_CUTOFFS = [10, 20, 33, 50]
BASELINE_K = 33
N_PERMUTATIONS = 10_000


def load_centroids() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load pre-computed centroids from CSVs."""
    human = pd.read_csv(INPUT_DIR / "centroids_human_35.csv", index_col=0)
    mouse = pd.read_csv(INPUT_DIR / "centroids_mouse_35.csv", index_col=0)
    cell_types = sorted(human.index.tolist())
    assert sorted(mouse.index.tolist()) == cell_types
    return human.loc[cell_types], mouse.loc[cell_types]


def run_at_cutoff(
    human_centroids: pd.DataFrame,
    mouse_centroids: pd.DataFrame,
    n_components: int,
) -> dict:
    """
    Run full Procrustes pipeline at a given PCA dimensionality.

    Returns dict with keys: n_components, variance_explained, distance,
    p_value, residual_magnitudes (dict cell_type→float), cell_types.
    """
    cell_types = human_centroids.index.tolist()
    n_types = len(cell_types)

    # Stack and fit PCA with fixed component count
    combined = np.vstack([human_centroids.values, mouse_centroids.values])  # (70, G)
    pca = PCA(n_components=n_components, svd_solver="full", random_state=RANDOM_SEED)
    combined_pca = pca.fit_transform(combined)  # (70, k)

    cumvar = float(np.sum(pca.explained_variance_ratio_))

    human_pca = combined_pca[:n_types]
    mouse_pca = combined_pca[n_types:]

    # Procrustes alignment
    result = procrustes_align(human_pca, mouse_pca)

    # Permutation test (silent)
    rng = np.random.RandomState(RANDOM_SEED)
    observed = result.distance
    null_distances = np.zeros(N_PERMUTATIONS)
    for i in range(N_PERMUTATIONS):
        perm = rng.permutation(n_types)
        null_distances[i] = _procrustes_distance(human_pca, mouse_pca[perm])
    n_leq = int(np.sum(null_distances <= observed))
    p_value = (n_leq + 1) / (N_PERMUTATIONS + 1)

    # Residual magnitudes
    residuals = {}
    for i, ct in enumerate(cell_types):
        r = result.aligned_target[i] - result.centered_reference[i]
        residuals[ct] = float(np.linalg.norm(r))

    return {
        "n_components": n_components,
        "variance_explained": cumvar,
        "distance": float(observed),
        "p_value": float(p_value),
        "residual_magnitudes": residuals,
        "cell_types": cell_types,
        "null_distances": null_distances,
    }


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 72)
    print("PCA DIMENSIONALITY SENSITIVITY ANALYSIS (T1-C)")
    print("=" * 72)

    # Load centroids
    human, mouse = load_centroids()
    n_types, n_genes = human.shape
    print(f"\nCentroids: {n_types} cell types × {n_genes:,} genes")
    print(f"PCA cutoffs to test: {PCA_CUTOFFS}")
    print(f"Permutations per cutoff: {N_PERMUTATIONS:,}")

    # Run at each cutoff
    results = {}
    for k in PCA_CUTOFFS:
        print(f"\n{'─' * 72}")
        print(f"Running PCA k={k} ...")
        results[k] = run_at_cutoff(human, mouse, k)
        r = results[k]
        print(f"  Variance explained: {r['variance_explained'] * 100:.1f}%")
        print(f"  Procrustes distance: {r['distance']:.4f}")
        print(f"  p-value: {r['p_value']:.6f}")

    # ---------------------------------------------------------------------------
    # Compute Spearman ρ vs baseline
    # ---------------------------------------------------------------------------
    baseline = results[BASELINE_K]
    baseline_ranking = pd.Series(baseline["residual_magnitudes"]).rank(ascending=False)

    rho_table = {}
    for k in PCA_CUTOFFS:
        r = results[k]
        ranking = pd.Series(r["residual_magnitudes"]).rank(ascending=False)
        rho, pval = spearmanr(baseline_ranking, ranking)
        rho_table[k] = {"rho": rho, "rho_pval": pval}

    # ---------------------------------------------------------------------------
    # Summary table
    # ---------------------------------------------------------------------------
    print(f"\n{'=' * 72}")
    print("SUMMARY TABLE")
    print(f"{'=' * 72}")
    header = f"{'Components':>12} | {'Var Expl':>10} | {'Distance':>10} | {'p-value':>10} | {'ρ vs k=33':>10} | {'Pass?':>6}"
    print(header)
    print("-" * len(header))

    all_pass = True
    rows = []
    for k in PCA_CUTOFFS:
        r = results[k]
        rho = rho_table[k]["rho"]
        p = r["p_value"]
        var_pct = r["variance_explained"] * 100
        dist = r["distance"]
        passes = p < 0.01 and rho >= 0.70
        if not passes:
            all_pass = False
        row = {
            "components": k,
            "variance_explained": f"{var_pct:.1f}%",
            "distance": f"{dist:.4f}",
            "p_value": f"{p:.6f}",
            "rho_vs_baseline": f"{rho:.4f}",
            "pass": passes,
        }
        rows.append(row)
        flag = "PASS" if passes else "FAIL"
        print(
            f"{k:>12} | {var_pct:>9.1f}% | {dist:>10.4f} | {p:>10.6f} | {rho:>10.4f} | {flag:>6}"
        )

    print(f"\n{'=' * 72}")
    verdict = "PASS" if all_pass else "FAIL"
    print(f"OVERALL VERDICT: {verdict}")
    print(f"Criterion: p < 0.01 AND ρ ≥ 0.70 at all cutoffs")
    print(f"{'=' * 72}")

    # ---------------------------------------------------------------------------
    # Save summary table as CSV
    # ---------------------------------------------------------------------------
    summary_df = pd.DataFrame(rows)
    summary_df.to_csv(OUTPUT_DIR / "pca_sensitivity_summary.csv", index=False)

    # Save per-cutoff residual magnitudes
    residual_df = pd.DataFrame(
        {f"k={k}": results[k]["residual_magnitudes"] for k in PCA_CUTOFFS}
    )
    residual_df.index.name = "cell_type"
    residual_df.to_csv(OUTPUT_DIR / "residual_magnitudes_by_k.csv")

    # Save full results as JSON
    json_results = {}
    for k in PCA_CUTOFFS:
        r = results[k]
        json_results[str(k)] = {
            "n_components": r["n_components"],
            "variance_explained": r["variance_explained"],
            "distance": r["distance"],
            "p_value": r["p_value"],
            "spearman_rho_vs_baseline": rho_table[k]["rho"],
            "spearman_rho_pval": rho_table[k]["rho_pval"],
            "residual_magnitudes": r["residual_magnitudes"],
            "null_distribution_summary": {
                "mean": float(np.mean(r["null_distances"])),
                "median": float(np.median(r["null_distances"])),
                "min": float(np.min(r["null_distances"])),
                "max": float(np.max(r["null_distances"])),
            },
        }
    json_results["verdict"] = verdict
    json_results["criterion"] = "p < 0.01 AND rho >= 0.70 at all cutoffs"
    json_results["baseline_k"] = BASELINE_K

    with open(OUTPUT_DIR / "pca_sensitivity_results.json", "w") as f:
        json.dump(json_results, f, indent=2)

    # ---------------------------------------------------------------------------
    # Scatter plots: residual ranking at each cutoff vs baseline (k=33)
    # ---------------------------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    non_baseline = [k for k in PCA_CUTOFFS if k != BASELINE_K]

    for ax, k in zip(axes, non_baseline):
        baseline_mags = pd.Series(baseline["residual_magnitudes"])
        test_mags = pd.Series(results[k]["residual_magnitudes"])

        # Rank (largest residual = rank 1)
        base_rank = baseline_mags.rank(ascending=False)
        test_rank = test_mags.rank(ascending=False)

        rho = rho_table[k]["rho"]
        ax.scatter(base_rank, test_rank, alpha=0.6, s=30, color="steelblue")

        # Identity line
        lims = [0.5, 35.5]
        ax.plot(lims, lims, "k--", alpha=0.3, linewidth=1)

        # Label most diverged and most conserved
        for ct in baseline_mags.nlargest(3).index:
            short = ct[:20] + "..." if len(ct) > 20 else ct
            ax.annotate(
                short,
                (base_rank[ct], test_rank[ct]),
                fontsize=6,
                alpha=0.7,
                textcoords="offset points",
                xytext=(4, 4),
            )
        for ct in baseline_mags.nsmallest(3).index:
            short = ct[:20] + "..." if len(ct) > 20 else ct
            ax.annotate(
                short,
                (base_rank[ct], test_rank[ct]),
                fontsize=6,
                alpha=0.7,
                textcoords="offset points",
                xytext=(4, -8),
            )

        ax.set_xlabel(f"Rank at k=33 (baseline)")
        ax.set_ylabel(f"Rank at k={k}")
        ax.set_title(f"k={k} vs baseline (ρ={rho:.3f})")
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_aspect("equal")

    fig.suptitle(
        "PCA Dimensionality Sensitivity: Residual Ranking Stability",
        fontsize=13,
        fontweight="bold",
    )
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "ranking_scatter_vs_baseline.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nScatter plot saved: {OUTPUT_DIR / 'ranking_scatter_vs_baseline.png'}")

    # ---------------------------------------------------------------------------
    # Distance vs components plot
    # ---------------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ks = [r["n_components"] for r in [results[k] for k in PCA_CUTOFFS]]
    dists = [r["distance"] for r in [results[k] for k in PCA_CUTOFFS]]
    pvals = [r["p_value"] for r in [results[k] for k in PCA_CUTOFFS]]
    varexps = [r["variance_explained"] * 100 for r in [results[k] for k in PCA_CUTOFFS]]

    ax1.plot(ks, dists, "o-", color="steelblue", linewidth=2, markersize=8)
    ax1.set_xlabel("PCA Components")
    ax1.set_ylabel("Procrustes Distance")
    ax1.set_title("Procrustes Distance vs PCA Dimensionality")
    ax1.axhline(y=dists[PCA_CUTOFFS.index(BASELINE_K)], color="gray", linestyle="--", alpha=0.5, label=f"Baseline k={BASELINE_K}")
    # Add variance explained as secondary axis
    ax1b = ax1.twinx()
    ax1b.plot(ks, varexps, "s--", color="coral", linewidth=1.5, markersize=6, alpha=0.7)
    ax1b.set_ylabel("Variance Explained (%)", color="coral")
    ax1b.tick_params(axis="y", labelcolor="coral")
    ax1.legend()

    ax2.semilogy(ks, pvals, "o-", color="darkred", linewidth=2, markersize=8)
    ax2.axhline(y=0.01, color="gray", linestyle="--", alpha=0.5, label="p=0.01 threshold")
    ax2.set_xlabel("PCA Components")
    ax2.set_ylabel("p-value (log scale)")
    ax2.set_title("Permutation Test p-value vs PCA Dimensionality")
    ax2.legend()

    fig.suptitle(
        "PCA Sensitivity: Distance and Significance",
        fontsize=13,
        fontweight="bold",
    )
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "distance_pvalue_vs_k.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Distance/p-value plot saved: {OUTPUT_DIR / 'distance_pvalue_vs_k.png'}")

    # ---------------------------------------------------------------------------
    # Text description of plots
    # ---------------------------------------------------------------------------
    print(f"\n{'=' * 72}")
    print("PLOT DESCRIPTIONS")
    print(f"{'=' * 72}")
    print(
        "\n1. ranking_scatter_vs_baseline.png: Three scatter plots showing the "
        "residual rank of each cell type at k=10, k=20, and k=50 vs the "
        "baseline k=33. Points near the diagonal indicate stable rankings. "
        "Spearman ρ is shown in each panel title."
    )
    print(
        "\n2. distance_pvalue_vs_k.png: Left panel shows Procrustes distance "
        "(blue circles) and variance explained (coral squares) vs PCA components. "
        "Right panel shows permutation p-value (log scale) vs components with "
        "the p=0.01 threshold line."
    )

    print(f"\nAll outputs saved to: {OUTPUT_DIR}/")
    return all_pass


if __name__ == "__main__":
    passed = main()
    sys.exit(0 if passed else 1)
