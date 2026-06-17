#!/usr/bin/env python3
"""
PCA Dimensionality Sensitivity Analysis — T1-C (v2, stricter criterion)

Biology
-------
Our Procrustes pipeline projects 35 cell-type centroids into a PCA subspace
before alignment. The default threshold (95% variance → 33 components) is
arbitrary. If the geometric signal is robust, the Procrustes distance should
remain significant and the per-cell-type residual (rigidity) ranking should
be stable across a range of dimensionalities.

Math
----
For each PCA cutoff k ∈ {10, 20, 33, 50}:
    1. Fit PCA(n_components=k) on the 70 combined centroids (35 human + 35 mouse).
    2. Run OPA alignment (mouse → human) in k-dimensional space.
    3. Run 10,000-iteration permutation test (seed=42, same as main analysis).
    4. Extract per-cell-type residual magnitudes ‖r_i‖ → rigidity ranking.
    5. Compute Spearman ρ of rigidity ranking vs the 95%-variance reference (k=33).

The reference run uses n_components=0.95 (scikit-learn variance threshold) to
confirm it resolves to exactly 33 components. All 4 fixed-k runs are compared
against this reference.

Pass criterion: p < 0.01 at ALL 4 cutoffs, rigidity ranking Spearman ρ > 0.75
vs reference at ALL cutoffs.
"""

# obs/null ratio uses null_median (canonical convention, matches src/procrustes.py).

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
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
)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

INPUT_DIR = PROJECT_ROOT / "output" / "phase2" / "scaled_35types"
OUTPUT_DIR = PROJECT_ROOT / "output" / "validation" / "pca_sensitivity"

PCA_CUTOFFS = [10, 20, 33, 50]
N_PERMUTATIONS = 10_000
RHO_THRESHOLD = 0.75
P_THRESHOLD = 0.01


def load_centroids() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load pre-computed 35-type centroids from CSVs."""
    human = pd.read_csv(INPUT_DIR / "centroids_human_35.csv", index_col=0)
    mouse = pd.read_csv(INPUT_DIR / "centroids_mouse_35.csv", index_col=0)
    cell_types = sorted(human.index.tolist())
    assert sorted(mouse.index.tolist()) == cell_types, "Cell type mismatch"
    return human.loc[cell_types], mouse.loc[cell_types]


def run_reference(
    human_centroids: pd.DataFrame,
    mouse_centroids: pd.DataFrame,
) -> dict:
    """
    Run PCA with n_components=0.95 (variance threshold) to determine how many
    components the 95% rule resolves to. This is the reference run.
    """
    cell_types = human_centroids.index.tolist()
    n_types = len(cell_types)

    combined = np.vstack([human_centroids.values, mouse_centroids.values])
    pca = PCA(n_components=0.95, svd_solver="full", random_state=RANDOM_SEED)
    combined_pca = pca.fit_transform(combined)

    n_components = pca.n_components_
    cumvar = float(np.sum(pca.explained_variance_ratio_))
    print(f"  95% variance threshold → {n_components} components ({cumvar*100:.1f}%)")

    human_pca = combined_pca[:n_types]
    mouse_pca = combined_pca[n_types:]

    result = procrustes_align(human_pca, mouse_pca)

    # Permutation test
    rng = np.random.RandomState(RANDOM_SEED)
    observed = result.distance
    null_distances = np.zeros(N_PERMUTATIONS)
    for i in range(N_PERMUTATIONS):
        perm = rng.permutation(n_types)
        null_distances[i] = _procrustes_distance(human_pca, mouse_pca[perm])
    n_leq = int(np.sum(null_distances <= observed))
    p_value = (n_leq + 1) / (N_PERMUTATIONS + 1)

    null_mean = float(np.mean(null_distances))
    obs_null_ratio = float(observed / np.median(null_distances))

    # Residual magnitudes
    residuals = {}
    for i, ct in enumerate(cell_types):
        r = result.aligned_target[i] - result.centered_reference[i]
        residuals[ct] = float(np.linalg.norm(r))

    return {
        "label": f"ref (95% var → {n_components})",
        "n_components": int(n_components),
        "variance_explained": cumvar,
        "distance": float(observed),
        "p_value": float(p_value),
        "obs_null_ratio": obs_null_ratio,
        "residual_magnitudes": residuals,
        "cell_types": cell_types,
        "null_distances": null_distances,
    }


def run_at_cutoff(
    human_centroids: pd.DataFrame,
    mouse_centroids: pd.DataFrame,
    n_components: int,
) -> dict:
    """
    Run full Procrustes pipeline at a given PCA dimensionality.

    Returns dict with: n_components, variance_explained, distance, p_value,
    obs_null_ratio, residual_magnitudes, cell_types, null_distances.
    """
    cell_types = human_centroids.index.tolist()
    n_types = len(cell_types)

    combined = np.vstack([human_centroids.values, mouse_centroids.values])
    pca = PCA(n_components=n_components, svd_solver="full", random_state=RANDOM_SEED)
    combined_pca = pca.fit_transform(combined)

    cumvar = float(np.sum(pca.explained_variance_ratio_))

    human_pca = combined_pca[:n_types]
    mouse_pca = combined_pca[n_types:]

    # Procrustes alignment
    result = procrustes_align(human_pca, mouse_pca)

    # Permutation test (same seed as main analysis)
    rng = np.random.RandomState(RANDOM_SEED)
    observed = result.distance
    null_distances = np.zeros(N_PERMUTATIONS)
    for i in range(N_PERMUTATIONS):
        perm = rng.permutation(n_types)
        null_distances[i] = _procrustes_distance(human_pca, mouse_pca[perm])
    n_leq = int(np.sum(null_distances <= observed))
    p_value = (n_leq + 1) / (N_PERMUTATIONS + 1)

    null_mean = float(np.mean(null_distances))
    obs_null_ratio = float(observed / np.median(null_distances))

    # Residual magnitudes
    residuals = {}
    for i, ct in enumerate(cell_types):
        r = result.aligned_target[i] - result.centered_reference[i]
        residuals[ct] = float(np.linalg.norm(r))

    return {
        "label": f"k={n_components}",
        "n_components": n_components,
        "variance_explained": cumvar,
        "distance": float(observed),
        "p_value": float(p_value),
        "obs_null_ratio": obs_null_ratio,
        "residual_magnitudes": residuals,
        "cell_types": cell_types,
        "null_distances": null_distances,
    }


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 72)
    print("PCA DIMENSIONALITY SENSITIVITY ANALYSIS (T1-C v2)")
    print(f"Pass criterion: p < {P_THRESHOLD} AND ρ > {RHO_THRESHOLD} at all cutoffs")
    print("=" * 72)

    # Load centroids
    human, mouse = load_centroids()
    n_types, n_genes = human.shape
    print(f"\nCentroids: {n_types} cell types × {n_genes:,} genes")
    print(f"PCA cutoffs: {PCA_CUTOFFS}")
    print(f"Permutations per cutoff: {N_PERMUTATIONS:,}")
    print(f"Random seed: {RANDOM_SEED}")

    # -----------------------------------------------------------------------
    # Reference run: 95% variance threshold
    # -----------------------------------------------------------------------
    print(f"\n{'─' * 72}")
    print("Reference run: PCA with 95% variance threshold ...")
    ref = run_reference(human, mouse)
    ref_k = ref["n_components"]
    print(f"  Distance: {ref['distance']:.4f}")
    print(f"  p-value: {ref['p_value']:.6f}")
    print(f"  obs/null ratio: {ref['obs_null_ratio']:.4f}")

    # -----------------------------------------------------------------------
    # Run at each fixed cutoff
    # -----------------------------------------------------------------------
    results = {"ref": ref}
    for k in PCA_CUTOFFS:
        print(f"\n{'─' * 72}")
        print(f"Running PCA k={k} ...")
        results[k] = run_at_cutoff(human, mouse, k)
        r = results[k]
        print(f"  Variance explained: {r['variance_explained'] * 100:.1f}%")
        print(f"  Procrustes distance: {r['distance']:.4f}")
        print(f"  p-value: {r['p_value']:.6f}")
        print(f"  obs/null ratio: {r['obs_null_ratio']:.4f}")

    # -----------------------------------------------------------------------
    # Compute Spearman ρ vs reference
    # -----------------------------------------------------------------------
    ref_ranking = pd.Series(ref["residual_magnitudes"]).rank(ascending=False)

    rho_table = {}
    for key in ["ref"] + PCA_CUTOFFS:
        r = results[key]
        ranking = pd.Series(r["residual_magnitudes"]).rank(ascending=False)
        rho, pval = spearmanr(ref_ranking, ranking)
        rho_table[key] = {"rho": rho, "rho_pval": pval}

    # -----------------------------------------------------------------------
    # Per-cell-type rigidity ranking at each cutoff
    # -----------------------------------------------------------------------
    rigidity_rankings = {}
    for key in ["ref"] + PCA_CUTOFFS:
        r = results[key]
        mags = pd.Series(r["residual_magnitudes"])
        # Low residual = high rigidity = rank 1
        rigidity_rankings[key] = mags.rank(ascending=True).astype(int)

    rigidity_df = pd.DataFrame(rigidity_rankings)
    rigidity_df.columns = [
        results[k]["label"] if k != "ref" else f"ref (k={ref_k})"
        for k in ["ref"] + PCA_CUTOFFS
    ]
    rigidity_df.index.name = "cell_type"
    rigidity_df.to_csv(OUTPUT_DIR / "rigidity_rankings_by_k.csv")

    # -----------------------------------------------------------------------
    # Summary table
    # -----------------------------------------------------------------------
    print(f"\n{'=' * 72}")
    print("SUMMARY TABLE")
    print(f"{'=' * 72}")
    header = (
        f"{'Label':>18} | {'Comp':>5} | {'Var%':>7} | {'Distance':>10} "
        f"| {'obs/null':>9} | {'p-value':>10} | {'ρ vs ref':>8} | {'Pass?':>6}"
    )
    print(header)
    print("-" * len(header))

    all_pass = True
    rows = []
    all_keys = ["ref"] + PCA_CUTOFFS

    for key in all_keys:
        r = results[key]
        rho = rho_table[key]["rho"]
        p = r["p_value"]
        var_pct = r["variance_explained"] * 100
        dist = r["distance"]
        obs_null = r["obs_null_ratio"]
        is_ref = key == "ref"

        if is_ref:
            passes = True  # reference always passes
            label = f"ref (95% → {ref_k})"
        else:
            passes = p < P_THRESHOLD and rho > RHO_THRESHOLD
            label = f"k={key}"
            if not passes:
                all_pass = False

        row = {
            "label": label,
            "components": r["n_components"],
            "variance_explained_pct": round(var_pct, 1),
            "distance": round(dist, 4),
            "obs_null_ratio": round(obs_null, 4),
            "p_value": p,
            "rho_vs_reference": round(rho, 4),
            "pass": passes,
        }
        rows.append(row)
        flag = "PASS" if passes else "FAIL"
        print(
            f"{label:>18} | {r['n_components']:>5} | {var_pct:>6.1f}% | "
            f"{dist:>10.4f} | {obs_null:>9.4f} | {p:>10.6f} | "
            f"{rho:>8.4f} | {flag:>6}"
        )

    print(f"\n{'=' * 72}")
    verdict = "PASS" if all_pass else "FAIL"
    print(f"OVERALL VERDICT: {verdict}")
    print(f"Criterion: p < {P_THRESHOLD} AND ρ > {RHO_THRESHOLD} at all 4 cutoffs")
    if not all_pass:
        failing = [
            r["label"] for r in rows
            if not r["pass"] and r["label"] != f"ref (95% → {ref_k})"
        ]
        print(f"FAILING cutoffs: {', '.join(failing)}")
    print(f"{'=' * 72}")

    # -----------------------------------------------------------------------
    # Save summary CSV
    # -----------------------------------------------------------------------
    summary_df = pd.DataFrame(rows)
    summary_df.to_csv(OUTPUT_DIR / "pca_sensitivity_summary.csv", index=False)

    # Save per-cutoff residual magnitudes
    residual_df = pd.DataFrame(
        {
            results[k]["label"] if k != "ref" else f"ref (k={ref_k})": results[k]["residual_magnitudes"]
            for k in all_keys
        }
    )
    residual_df.index.name = "cell_type"
    residual_df.to_csv(OUTPUT_DIR / "residual_magnitudes_by_k.csv")

    # Save full results JSON
    json_results = {}
    for key in all_keys:
        r = results[key]
        json_key = "ref" if key == "ref" else str(key)
        json_results[json_key] = {
            "label": r["label"] if key != "ref" else f"ref (95% var → {ref_k})",
            "n_components": r["n_components"],
            "variance_explained": r["variance_explained"],
            "distance": r["distance"],
            "obs_null_ratio": r["obs_null_ratio"],
            "p_value": r["p_value"],
            "spearman_rho_vs_reference": rho_table[key]["rho"],
            "spearman_rho_pval": rho_table[key]["rho_pval"],
            "residual_magnitudes": r["residual_magnitudes"],
            "null_distribution_summary": {
                "mean": float(np.mean(r["null_distances"])),
                "median": float(np.median(r["null_distances"])),
                "min": float(np.min(r["null_distances"])),
                "max": float(np.max(r["null_distances"])),
            },
        }
    json_results["verdict"] = verdict
    json_results["criterion"] = f"p < {P_THRESHOLD} AND rho > {RHO_THRESHOLD} at all cutoffs"
    json_results["reference_k"] = ref_k
    json_results["reference_variance_pct"] = round(ref["variance_explained"] * 100, 1)

    with open(OUTPUT_DIR / "pca_sensitivity_results.json", "w") as f:
        json.dump(json_results, f, indent=2)

    # -----------------------------------------------------------------------
    # Plot 1: Ranking scatter vs reference
    # -----------------------------------------------------------------------
    non_ref_keys = [k for k in PCA_CUTOFFS if k != ref_k]
    fig, axes = plt.subplots(1, len(non_ref_keys), figsize=(5 * len(non_ref_keys), 5))
    if len(non_ref_keys) == 1:
        axes = [axes]

    ref_mags = pd.Series(ref["residual_magnitudes"])
    ref_rank = ref_mags.rank(ascending=False)

    for ax, k in zip(axes, non_ref_keys):
        test_mags = pd.Series(results[k]["residual_magnitudes"])
        test_rank = test_mags.rank(ascending=False)
        rho = rho_table[k]["rho"]
        passes = rho > RHO_THRESHOLD

        ax.scatter(ref_rank, test_rank, alpha=0.6, s=30,
                   color="steelblue" if passes else "tomato")

        # Identity line
        lims = [0.5, 35.5]
        ax.plot(lims, lims, "k--", alpha=0.3, linewidth=1)

        # Label extremes
        for ct in ref_mags.nlargest(3).index:
            short = ct[:18] + ".." if len(ct) > 18 else ct
            ax.annotate(short, (ref_rank[ct], test_rank[ct]),
                        fontsize=5.5, alpha=0.7,
                        textcoords="offset points", xytext=(4, 4))
        for ct in ref_mags.nsmallest(3).index:
            short = ct[:18] + ".." if len(ct) > 18 else ct
            ax.annotate(short, (ref_rank[ct], test_rank[ct]),
                        fontsize=5.5, alpha=0.7,
                        textcoords="offset points", xytext=(4, -8))

        flag = "PASS" if passes else "FAIL"
        ax.set_xlabel(f"Rank at ref (k={ref_k})")
        ax.set_ylabel(f"Rank at k={k}")
        ax.set_title(f"k={k} vs ref (ρ={rho:.3f}, {flag})")
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_aspect("equal")

    fig.suptitle(
        f"PCA Sensitivity: Rigidity Ranking Stability (threshold ρ>{RHO_THRESHOLD})",
        fontsize=13, fontweight="bold",
    )
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "ranking_scatter_vs_reference.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nScatter plot saved: {OUTPUT_DIR / 'ranking_scatter_vs_reference.png'}")

    # -----------------------------------------------------------------------
    # Plot 2: Distance, obs/null, and p-value vs components
    # -----------------------------------------------------------------------
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 5))

    ks = [results[k]["n_components"] for k in PCA_CUTOFFS]
    dists = [results[k]["distance"] for k in PCA_CUTOFFS]
    obs_nulls = [results[k]["obs_null_ratio"] for k in PCA_CUTOFFS]
    pvals = [results[k]["p_value"] for k in PCA_CUTOFFS]
    varexps = [results[k]["variance_explained"] * 100 for k in PCA_CUTOFFS]

    # Panel 1: Distance + variance
    ax1.plot(ks, dists, "o-", color="steelblue", linewidth=2, markersize=8)
    ax1.set_xlabel("PCA Components")
    ax1.set_ylabel("Procrustes Distance", color="steelblue")
    ax1.set_title("Distance vs Dimensionality")
    ax1.axhline(y=ref["distance"], color="gray", linestyle="--", alpha=0.5,
                label=f"ref k={ref_k}")
    ax1b = ax1.twinx()
    ax1b.plot(ks, varexps, "s--", color="coral", linewidth=1.5, markersize=6, alpha=0.7)
    ax1b.set_ylabel("Variance Explained (%)", color="coral")
    ax1b.tick_params(axis="y", labelcolor="coral")
    ax1.legend(loc="lower right")

    # Panel 2: obs/null ratio
    ax2.plot(ks, obs_nulls, "o-", color="teal", linewidth=2, markersize=8)
    ax2.axhline(y=ref["obs_null_ratio"], color="gray", linestyle="--", alpha=0.5,
                label=f"ref k={ref_k}")
    ax2.set_xlabel("PCA Components")
    ax2.set_ylabel("obs/null ratio")
    ax2.set_title("Obs/Null Ratio vs Dimensionality")
    ax2.legend()

    # Panel 3: p-value
    ax3.semilogy(ks, pvals, "o-", color="darkred", linewidth=2, markersize=8)
    ax3.axhline(y=P_THRESHOLD, color="gray", linestyle="--", alpha=0.5,
                label=f"p={P_THRESHOLD} threshold")
    ax3.set_xlabel("PCA Components")
    ax3.set_ylabel("p-value (log scale)")
    ax3.set_title("Permutation p-value vs Dimensionality")
    ax3.legend()

    fig.suptitle(
        "PCA Sensitivity: Distance, Signal Strength, and Significance",
        fontsize=13, fontweight="bold",
    )
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "distance_obsnull_pvalue_vs_k.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Distance/obs-null/p-value plot saved: {OUTPUT_DIR / 'distance_obsnull_pvalue_vs_k.png'}")

    # -----------------------------------------------------------------------
    # Text descriptions
    # -----------------------------------------------------------------------
    print(f"\n{'=' * 72}")
    print("PLOT DESCRIPTIONS")
    print(f"{'=' * 72}")
    print(
        "\n1. ranking_scatter_vs_reference.png: Three scatter plots showing the "
        "residual (rigidity) rank of each cell type at k=10, k=20, and k=50 vs "
        f"the reference k={ref_k}. Points near the diagonal = stable rankings. "
        f"Blue = passes ρ>{RHO_THRESHOLD}, red = fails. Spearman ρ in each title."
    )
    print(
        "\n2. distance_obsnull_pvalue_vs_k.png: Three panels. Left: Procrustes "
        "distance (blue) + variance explained (coral) vs components. Center: "
        "obs/null ratio vs components (lower = stronger signal). Right: "
        "permutation p-value (log scale) vs components with threshold line."
    )

    print(f"\nAll outputs saved to: {OUTPUT_DIR}/")
    return all_pass


if __name__ == "__main__":
    passed = main()
    sys.exit(0 if passed else 1)
