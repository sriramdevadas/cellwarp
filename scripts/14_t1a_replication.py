#!/usr/bin/env python3
"""
CellWarp — T1-A Replication: Procrustes Pipeline

Runs the identical Procrustes analysis from the primary result on independent data:
  - Mouse: Mouse Cell Atlas (Han et al. 2018, BGI, Microwell-seq)
  - Human: Pooled non-Tabula CELLxGENE Census adult healthy data

This is the direct answer to the Tabula batch effect critique. If the cross-species
geometric structure replicates in completely independent data, the signal is real.

Biology
-------
The primary finding: cross-species cell type centroids align geometrically
(Procrustes p=0.0001 for 35 types). The replication tests whether this structure
is atlas-specific or a genuine biological property of homologous cell types.

Steps
-----
  1. Load MCA and HCA replication data
  2. Compute centroids (mean expression per cell type, no per-donor pooling)
  3. PCA on combined centroids (95% variance threshold)
  4. Procrustes alignment (MCA → HCA, same direction as primary)
  5. Permutation test (10,000 iterations, seed=42)
  6. Rigidity ranking correlation with primary Tabula result
  7. T1-B negative control (obs/null ratio comparison)
  8. Generate all outputs and plots

Output
------
  output/validation/t1a_replication/t1a_results.json
  output/validation/t1a_replication/ranking_comparison.csv
  output/validation/t1a_replication/t1a_rigidity_scatter.png
  output/validation/t1a_replication/null_distribution.png

Usage:
    python scripts/14_t1a_replication.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as scipy_stats

# ---------------------------------------------------------------------------
# Project setup
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from cellwarp.procrustes import (
    compute_centroids,
    compute_residual_vectors,
    map_residuals_to_genes,
    pca_reduce_centroids,
    permutation_test,
    procrustes_align,
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

MCA_PATH = PROJECT_ROOT / "data" / "replication" / "mca_t1a.h5ad"
HCA_PATH = PROJECT_ROOT / "data" / "replication" / "hca_t1a.h5ad"
OUTPUT_DIR = PROJECT_ROOT / "output" / "validation" / "t1a_replication"
PRIMARY_RESULTS_PATH = (
    PROJECT_ROOT / "output" / "phase2" / "scaled_35types" / "procrustes_results_35.json"
)
PRIMARY_6TYPE_RESULTS_PATH = (
    PROJECT_ROOT / "output" / "phase2" / "procrustes_results.json"
)

# The original 6 cell types from Phase 2
ORIGINAL_6_TYPES = [
    "B cell",
    "CD4-positive, alpha-beta T cell",
    "CD8-positive, alpha-beta T cell",
    "endothelial cell",
    "hepatocyte",
    "macrophage",
]

RANDOM_SEED = 42
N_PERMUTATIONS = 10_000
CELL_TYPE_COL = "our_cell_type_label"


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def load_replication_data() -> tuple[ad.AnnData, ad.AnnData, list[str]]:
    """
    Load MCA and HCA replication data and determine the intersection of cell types.

    Returns:
        Tuple of (mca_adata, hca_adata, shared_cell_types)
    """
    print("=" * 70)
    print("STEP 1: LOADING REPLICATION DATA")
    print("=" * 70)

    # Load MCA (mouse)
    print(f"\n  Loading MCA: {MCA_PATH}")
    mca = ad.read_h5ad(MCA_PATH)
    print(f"  MCA: {mca.n_obs:,} cells × {mca.n_vars:,} genes")

    mca_types = set(mca.obs[CELL_TYPE_COL].dropna().unique())
    print(f"  MCA cell types ({len(mca_types)}):")
    for ct in sorted(mca_types):
        n = int((mca.obs[CELL_TYPE_COL] == ct).sum())
        print(f"    {ct:<50} {n:>6,}")

    # Load HCA (human)
    print(f"\n  Loading HCA: {HCA_PATH}")
    hca = ad.read_h5ad(HCA_PATH)
    print(f"  HCA: {hca.n_obs:,} cells × {hca.n_vars:,} genes")

    hca_types = set(hca.obs[CELL_TYPE_COL].dropna().unique())
    print(f"  HCA cell types ({len(hca_types)}):")
    for ct in sorted(hca_types):
        n = int((hca.obs[CELL_TYPE_COL] == ct).sum())
        print(f"    {ct:<50} {n:>6,}")

    # Intersection
    shared_types = sorted(mca_types & hca_types)
    print(f"\n  Shared cell types: {len(shared_types)}")
    for ct in shared_types:
        n_mca = int((mca.obs[CELL_TYPE_COL] == ct).sum())
        n_hca = int((hca.obs[CELL_TYPE_COL] == ct).sum())
        in_original = "***" if ct in ORIGINAL_6_TYPES else ""
        print(f"    {ct:<50} MCA={n_mca:>6,}  HCA={n_hca:>6,}  {in_original}")

    # Stop condition: <15 types
    if len(shared_types) < 15:
        print(f"\n  STOP: Only {len(shared_types)} intersection types (< 15 minimum)")
        print("  Types lost:")
        for ct in sorted(mca_types | hca_types):
            if ct not in shared_types:
                where = "MCA only" if ct in mca_types else "HCA only"
                print(f"    {ct}: {where}")
        sys.exit(1)

    # Filter both datasets to shared types
    mca = mca[mca.obs[CELL_TYPE_COL].isin(shared_types)].copy()
    hca = hca[hca.obs[CELL_TYPE_COL].isin(shared_types)].copy()

    # Verify gene spaces match
    assert list(mca.var_names) == list(hca.var_names), (
        f"Gene space mismatch: MCA has {mca.n_vars} genes, HCA has {hca.n_vars}"
    )

    return mca, hca, shared_types


def load_primary_results() -> dict:
    """
    Load the primary 35-type Procrustes results for comparison.

    Returns:
        Dict with primary result data.
    """
    print(f"\n  Loading primary 35-type results: {PRIMARY_RESULTS_PATH}")
    with open(PRIMARY_RESULTS_PATH) as f:
        primary = json.load(f)

    print(f"  Primary: {len(primary['cell_types'])} types, "
          f"distance={primary['procrustes']['distance']:.4f}, "
          f"p={primary['permutation_test']['p_value']:.6f}")

    # Also load 6-type results
    print(f"  Loading primary 6-type results: {PRIMARY_6TYPE_RESULTS_PATH}")
    with open(PRIMARY_6TYPE_RESULTS_PATH) as f:
        primary_6 = json.load(f)

    primary["six_type"] = primary_6
    return primary


# ---------------------------------------------------------------------------
# Procrustes replication
# ---------------------------------------------------------------------------


def run_procrustes_replication(
    mca: ad.AnnData,
    hca: ad.AnnData,
    shared_types: list[str],
) -> dict:
    """
    Run the identical Procrustes pipeline on replication data.

    Steps:
      a) Compute per-type centroids (grand mean, no per-donor pooling)
      b) PCA on combined centroids (95% variance threshold)
      c) Procrustes alignment (MCA → HCA = mouse → human)
      d) Permutation test (10,000 iterations)
      e) Compute residual vectors and map to genes

    Args:
        mca: Mouse Cell Atlas AnnData.
        hca: Pooled HCA AnnData.
        shared_types: Sorted list of cell types in both datasets.

    Returns:
        Dict with all replication results.
    """
    print("\n" + "=" * 70)
    print("STEP 2: PROCRUSTES REPLICATION PIPELINE")
    print("=" * 70)

    # --- Step 2a: Compute centroids ---
    print("\n  Computing MCA (mouse) centroids...")
    mca_centroids = compute_centroids(mca, CELL_TYPE_COL)

    print("\n  Computing HCA (human) centroids...")
    hca_centroids = compute_centroids(hca, CELL_TYPE_COL)

    # Ensure same cell types in same order
    assert sorted(mca_centroids.index.tolist()) == shared_types
    assert sorted(hca_centroids.index.tolist()) == shared_types

    # --- Step 2b: PCA on combined centroids ---
    print("\n  PCA on combined centroids...")
    human_pca, mouse_pca, pca_model, cell_types = pca_reduce_centroids(
        hca_centroids, mca_centroids, variance_threshold=0.95
    )

    n_components = pca_model.n_components_
    cumvar = np.cumsum(pca_model.explained_variance_ratio_)[-1]

    # --- Step 2c: Procrustes alignment ---
    print("\n  Procrustes alignment (MCA → HCA)...")
    result = procrustes_align(human_pca, mouse_pca)

    # --- Step 2d: Permutation test ---
    print("\n  Permutation test...")
    p_value, null_distribution = permutation_test(
        human_pca, mouse_pca,
        n_permutations=N_PERMUTATIONS,
        seed=RANDOM_SEED,
    )

    # --- Step 2e: Residual vectors ---
    print("\n  Computing residual deformation vectors...")
    residuals = compute_residual_vectors(result, cell_types)

    # Map residuals to genes
    gene_names = list(hca.var_names)
    top_genes = map_residuals_to_genes(residuals, pca_model, gene_names)

    # Compute obs/null ratio
    null_median = float(np.median(null_distribution))
    obs_null_ratio = float(result.distance / null_median)

    # Compile results
    replication_results = {
        "procrustes": {
            "distance": float(result.distance),
            "distance_squared": float(result.distance_squared),
            "scaling": float(result.scaling),
            "rotation_det": float(np.linalg.det(result.rotation)),
        },
        "permutation_test": {
            "p_value": float(p_value),
            "n_permutations": N_PERMUTATIONS,
            "null_distribution_summary": {
                "mean": float(np.mean(null_distribution)),
                "median": null_median,
                "std": float(np.std(null_distribution)),
                "min": float(np.min(null_distribution)),
                "max": float(np.max(null_distribution)),
            },
        },
        "pca": {
            "n_components": n_components,
            "cumulative_variance": float(cumvar),
            "n_genes_input": int(hca.n_vars),
        },
        "obs_null_ratio": obs_null_ratio,
        "n_types": len(cell_types),
        "cell_types": cell_types,
        "residuals": {
            ct: {
                "magnitude": float(np.linalg.norm(residuals[ct])),
            }
            for ct in cell_types
        },
        "top_genes": {
            ct: top_genes[ct][["gene", "loading", "abs_loading", "rank"]]
            .to_dict(orient="records")
            for ct in cell_types
        },
        "null_distribution": null_distribution.tolist(),
    }

    return replication_results


# ---------------------------------------------------------------------------
# Ranking comparison
# ---------------------------------------------------------------------------


def compute_ranking_comparison(
    replication_results: dict,
    primary: dict,
) -> tuple[pd.DataFrame, dict]:
    """
    Compare rigidity rankings between T1-A replication and primary Tabula result.

    Matches cell types by exact string name. Computes Spearman correlation
    of per-type residual magnitudes.

    Args:
        replication_results: T1-A Procrustes output.
        primary: Primary 35-type Procrustes output.

    Returns:
        Tuple of (ranking_df, comparison_stats).
    """
    print("\n" + "=" * 70)
    print("STEP 3: RIGIDITY RANKING COMPARISON")
    print("=" * 70)

    # Extract residual magnitudes
    t1a_residuals = {
        ct: info["magnitude"]
        for ct, info in replication_results["residuals"].items()
    }
    tabula_residuals = {
        ct: info["magnitude"]
        for ct, info in primary["residuals"].items()
    }

    # Find matched types
    matched_types = sorted(set(t1a_residuals.keys()) & set(tabula_residuals.keys()))
    print(f"\n  Matched cell types: {len(matched_types)}")

    if len(matched_types) < 3:
        print("  WARNING: Fewer than 3 matched types — correlation meaningless")
        return pd.DataFrame(), {"rho": float("nan"), "p_value": float("nan")}

    # Build comparison table
    rows = []
    for ct in matched_types:
        rows.append({
            "cell_type": ct,
            "tabula_residual": tabula_residuals[ct],
            "t1a_residual": t1a_residuals[ct],
            "included_in_original_6": ct in ORIGINAL_6_TYPES,
        })

    df = pd.DataFrame(rows)

    # Compute ranks (rank 1 = smallest residual = most rigid)
    df["tabula_rank"] = df["tabula_residual"].rank()
    df["t1a_rank"] = df["t1a_residual"].rank()

    # Spearman correlation
    rho, p_val = scipy_stats.spearmanr(
        df["tabula_residual"].values,
        df["t1a_residual"].values,
    )

    print(f"\n  Spearman correlation of residual magnitudes:")
    print(f"    rho = {rho:.4f}")
    print(f"    p   = {p_val:.6f}")
    print(f"    N   = {len(matched_types)} matched types")
    print(f"    Target: rho >= 0.70")
    print(f"    Verdict: {'GO' if rho >= 0.70 else 'NO-GO'}")

    # Print per-type comparison
    print(f"\n  {'Cell Type':<50} {'Tabula Rank':>12} {'T1-A Rank':>10} "
          f"{'Tabula Resid':>13} {'T1-A Resid':>11} {'Orig 6':>7}")
    print("  " + "-" * 105)
    for _, row in df.sort_values("tabula_rank").iterrows():
        orig = "*" if row["included_in_original_6"] else ""
        print(
            f"  {row['cell_type']:<50} "
            f"{row['tabula_rank']:>12.0f} "
            f"{row['t1a_rank']:>10.0f} "
            f"{row['tabula_residual']:>13.4f} "
            f"{row['t1a_residual']:>11.4f} "
            f"{orig:>7}"
        )

    comparison_stats = {
        "rho": float(rho),
        "p_value": float(p_val),
        "n_matched": len(matched_types),
        "target_rho": 0.70,
        "pass": bool(rho >= 0.70),
    }

    return df, comparison_stats


# ---------------------------------------------------------------------------
# T1-B Negative Control
# ---------------------------------------------------------------------------


def t1b_negative_control(
    replication_results: dict,
    primary: dict,
) -> dict:
    """
    T1-B: Compare obs/null ratios across atlases.

    If T1-A obs/null is in the same range as Tabula, the signal is consistent
    across atlases — not Tabula-specific.

    Args:
        replication_results: T1-A results.
        primary: Primary results (35-type and 6-type).

    Returns:
        Dict with T1-B comparison stats.
    """
    print("\n" + "=" * 70)
    print("STEP 4: T1-B NEGATIVE CONTROL (OBS/NULL RATIO COMPARISON)")
    print("=" * 70)

    # T1-A obs/null ratio
    t1a_ratio = replication_results["obs_null_ratio"]

    # Primary 35-type obs/null ratio
    tabula_35_distance = primary["procrustes"]["distance"]
    tabula_35_null_median = primary["permutation_test"]["null_distribution_summary"]["median"]
    tabula_35_ratio = tabula_35_distance / tabula_35_null_median

    # Primary 6-type obs/null ratio
    tabula_6_distance = primary["six_type"]["procrustes"]["distance"]
    tabula_6_null_median = primary["six_type"]["permutation_test"]["null_distribution_summary"]["median"]
    tabula_6_ratio = tabula_6_distance / tabula_6_null_median

    print(f"\n  Obs/Null ratios (lower = stronger signal):")
    print(f"    T1-A replication ({replication_results['n_types']} types): {t1a_ratio:.4f}")
    print(f"    Tabula 35-type:                                  {tabula_35_ratio:.4f}")
    print(f"    Tabula 6-type:                                   {tabula_6_ratio:.4f}")
    print(f"\n  Interpretation:")
    if t1a_ratio < tabula_35_ratio * 1.5:
        print("    T1-A obs/null ratio is in the SAME RANGE as Tabula —")
        print("    signal strength CONSISTENT across independent atlases.")
    else:
        print("    T1-A obs/null ratio is HIGHER than Tabula —")
        print("    weaker signal in replication data (expected: cross-study noise inflates centroid distances).")

    return {
        "t1a_obs_null": float(t1a_ratio),
        "t1a_n_types": replication_results["n_types"],
        "tabula_35_obs_null": float(tabula_35_ratio),
        "tabula_6_obs_null": float(tabula_6_ratio),
    }


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def plot_rigidity_scatter(
    ranking_df: pd.DataFrame,
    comparison_stats: dict,
    output_path: Path,
) -> None:
    """
    Scatter plot of Tabula rank vs T1-A rank with regression line.

    Original 6 types are marked with a different color/marker.
    Annotated with rho and p-value.
    """
    fig, ax = plt.subplots(figsize=(10, 10))

    # Split into original 6 and expansion types
    orig6 = ranking_df[ranking_df["included_in_original_6"]]
    expansion = ranking_df[~ranking_df["included_in_original_6"]]

    # Plot expansion types
    if len(expansion) > 0:
        ax.scatter(
            expansion["tabula_rank"],
            expansion["t1a_rank"],
            c="steelblue",
            s=80,
            alpha=0.7,
            edgecolors="navy",
            linewidths=0.5,
            label="Expansion types",
            zorder=3,
        )

    # Plot original 6 types
    if len(orig6) > 0:
        ax.scatter(
            orig6["tabula_rank"],
            orig6["t1a_rank"],
            c="crimson",
            s=120,
            alpha=0.9,
            edgecolors="darkred",
            linewidths=1.0,
            marker="D",
            label="Original 6 types",
            zorder=4,
        )

    # Label all points
    for _, row in ranking_df.iterrows():
        # Truncate long names
        label = row["cell_type"]
        if len(label) > 25:
            label = label[:22] + "..."
        ax.annotate(
            label,
            (row["tabula_rank"], row["t1a_rank"]),
            fontsize=6,
            ha="left",
            va="bottom",
            xytext=(4, 4),
            textcoords="offset points",
        )

    # Regression line
    n = len(ranking_df)
    x = ranking_df["tabula_rank"].values
    y = ranking_df["t1a_rank"].values
    if n >= 3:
        slope, intercept = np.polyfit(x, y, 1)
        x_line = np.linspace(x.min(), x.max(), 100)
        ax.plot(x_line, slope * x_line + intercept, "k--", alpha=0.3, linewidth=1)

    # Identity line
    max_rank = max(x.max(), y.max())
    ax.plot([0, max_rank + 1], [0, max_rank + 1], "gray", alpha=0.2, linewidth=1,
            linestyle=":", label="Identity line")

    # Annotation
    rho = comparison_stats["rho"]
    p_val = comparison_stats["p_value"]
    ax.text(
        0.05, 0.95,
        f"Spearman rho = {rho:.3f}\np = {p_val:.4f}\nN = {n} types\n"
        f"Target: rho >= 0.70 {'PASS' if rho >= 0.70 else 'FAIL'}",
        transform=ax.transAxes,
        fontsize=11,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="wheat", alpha=0.8),
    )

    ax.set_xlabel("Tabula Rigidity Rank (35-type)", fontsize=12)
    ax.set_ylabel("T1-A Rigidity Rank (replication)", fontsize=12)
    ax.set_title("T1-A Replication: Rigidity Ranking Comparison", fontsize=14)
    ax.legend(loc="lower right", fontsize=10)
    ax.set_aspect("equal", adjustable="box")

    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")

    # Text description
    print(f"\n  Plot description: Scatter plot of cell type rigidity ranks.")
    print(f"  X-axis: rank from primary 35-type Tabula analysis")
    print(f"  Y-axis: rank from T1-A replication (MCA × HCA)")
    print(f"  Red diamonds: original 6 cell types. Blue circles: expansion types.")
    print(f"  Spearman rho={rho:.3f}, p={p_val:.4f}")


def plot_null_distribution(
    replication_results: dict,
    output_path: Path,
) -> None:
    """
    Histogram of T1-A permutation null distribution with observed distance.
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    null_dist = np.array(replication_results["null_distribution"])
    obs_distance = replication_results["procrustes"]["distance"]
    p_value = replication_results["permutation_test"]["p_value"]

    ax.hist(
        null_dist,
        bins=80,
        color="steelblue",
        alpha=0.7,
        edgecolor="navy",
        linewidth=0.3,
        density=True,
        label="Null distribution",
    )

    ax.axvline(
        obs_distance,
        color="crimson",
        linewidth=2.5,
        linestyle="--",
        label=f"Observed (d={obs_distance:.2f})",
    )

    ax.text(
        0.05, 0.95,
        f"T1-A Replication\n"
        f"Observed: {obs_distance:.2f}\n"
        f"Null median: {np.median(null_dist):.2f}\n"
        f"p = {p_value:.6f}\n"
        f"N = {replication_results['n_types']} types, "
        f"{replication_results['pca']['n_components']} PCs",
        transform=ax.transAxes,
        fontsize=11,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="lightyellow", alpha=0.8),
    )

    ax.set_xlabel("Procrustes Distance", fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    ax.set_title(
        "T1-A Replication: Permutation Null Distribution", fontsize=14
    )
    ax.legend(loc="upper right", fontsize=10)

    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")

    # Text description
    print(f"\n  Plot description: Histogram of {N_PERMUTATIONS:,} permuted Procrustes")
    print(f"  distances (blue) with observed distance (red dashed) at {obs_distance:.2f}.")
    print(f"  Null median={np.median(null_dist):.2f}. p={p_value:.6f}.")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------


def main() -> None:
    """Run the complete T1-A replication pipeline."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # ── Step 1: Load data ──
    mca, hca, shared_types = load_replication_data()

    # ── Load primary results for comparison ──
    primary = load_primary_results()

    # ── Step 2: Procrustes replication ──
    replication_results = run_procrustes_replication(mca, hca, shared_types)

    # ── Step 3: Ranking comparison ──
    ranking_df, comparison_stats = compute_ranking_comparison(
        replication_results, primary
    )

    # ── Step 4: T1-B negative control ──
    t1b_stats = t1b_negative_control(replication_results, primary)

    # ── Step 5: Save outputs ──
    print("\n" + "=" * 70)
    print("STEP 5: SAVING OUTPUTS")
    print("=" * 70)

    # Load rescue stats if available
    rescue_stats = {}
    rescue_path = PROJECT_ROOT / "data" / "replication" / "mca_download_stats.json"
    if rescue_path.exists():
        with open(rescue_path) as f:
            mca_stats = json.load(f)
        rescue_stats = mca_stats.get("rescue_stats", {})

    # Count per-type cells (ensure native int for JSON serialization)
    mca_per_type = {
        ct: int((mca.obs[CELL_TYPE_COL] == ct).sum())
        for ct in shared_types
    }
    hca_per_type = {
        ct: int((hca.obs[CELL_TYPE_COL] == ct).sum())
        for ct in shared_types
    }

    # Helper for JSON serialization of numpy types
    class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, (np.integer,)):
                return int(obj)
            if isinstance(obj, (np.floating,)):
                return float(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return super().default(obj)

    # Original 6 in replication
    orig6_in_replication = [ct for ct in ORIGINAL_6_TYPES if ct in shared_types]

    # Compile final results JSON
    results_json = {
        "t1a_procrustes": {
            "distance": replication_results["procrustes"]["distance"],
            "p_value": replication_results["permutation_test"]["p_value"],
            "obs_null_ratio": replication_results["obs_null_ratio"],
            "scaling": replication_results["procrustes"]["scaling"],
            "n_types": replication_results["n_types"],
            "n_components": replication_results["pca"]["n_components"],
            "cumulative_variance": replication_results["pca"]["cumulative_variance"],
        },
        "ranking_comparison": {
            "rho": comparison_stats["rho"],
            "p_value": comparison_stats["p_value"],
            "n_matched": comparison_stats["n_matched"],
            "target_rho": 0.70,
            "pass": comparison_stats["pass"],
        },
        "t1b_obs_null_comparison": t1b_stats,
        "per_type_residuals": {
            ct: {
                "t1a_residual": replication_results["residuals"][ct]["magnitude"],
                "mca_cells": mca_per_type.get(ct, 0),
                "hca_cells": hca_per_type.get(ct, 0),
            }
            for ct in shared_types
        },
        "original_6_in_replication": {
            "count": len(orig6_in_replication),
            "types": orig6_in_replication,
            "missing": [ct for ct in ORIGINAL_6_TYPES if ct not in shared_types],
        },
        "cd4_rescue": rescue_stats,
        "mca_hepatocyte_cells": mca_per_type.get("hepatocyte", 0),
        "excluded_tissues": [],  # Will be populated from download stats
    }

    # Add excluded tissues from MCA download
    if rescue_path.exists():
        results_json["excluded_tissues"] = mca_stats.get(
            "download_stats", {}
        ).get("excluded_tissues", [])

    # Save results JSON
    results_path = OUTPUT_DIR / "t1a_results.json"
    with open(results_path, "w") as f:
        json.dump(results_json, f, indent=2, cls=NumpyEncoder)
    print(f"  Results JSON: {results_path}")
    np.save(OUTPUT_DIR / "null_distribution.npy", np.array(replication_results["null_distribution"]))

    # Save ranking comparison CSV
    ranking_path = OUTPUT_DIR / "ranking_comparison.csv"
    if len(ranking_df) > 0:
        ranking_df.to_csv(ranking_path, index=False)
        print(f"  Ranking CSV: {ranking_path}")

    # Generate plots
    print("\n  Generating plots...")
    if len(ranking_df) > 0:
        plot_rigidity_scatter(
            ranking_df, comparison_stats,
            OUTPUT_DIR / "t1a_rigidity_scatter.png",
        )

    plot_null_distribution(
        replication_results,
        OUTPUT_DIR / "null_distribution.png",
    )

    # ── Final printed summary ──
    print("\n" + "=" * 70)
    print("T1-A REPLICATION — FINAL SUMMARY")
    print("=" * 70)

    p_val = replication_results["permutation_test"]["p_value"]
    distance = replication_results["procrustes"]["distance"]
    obs_null = replication_results["obs_null_ratio"]
    n_types = replication_results["n_types"]
    n_comps = replication_results["pca"]["n_components"]

    print(f"\n  1. T1-A Procrustes: distance={distance:.4f}, "
          f"p={p_val:.6f}, obs/null={obs_null:.4f} "
          f"({n_types} types, {n_comps} PCA components)")

    rho = comparison_stats["rho"]
    rho_p = comparison_stats["p_value"]
    n_matched = comparison_stats["n_matched"]
    print(f"\n  2. Rigidity ranking correlation: "
          f"rho={rho:.4f}, p={rho_p:.6f}, N={n_matched} matched types")

    go_nogo = "GO" if rho >= 0.70 else "NO-GO"
    print(f"\n  3. GO / NO-GO vs target rho >= 0.70: [{go_nogo}]")

    print(f"\n  4. Original 6 types in replication: "
          f"{len(orig6_in_replication)}/6 "
          f"({', '.join(orig6_in_replication)})")
    missing_6 = [ct for ct in ORIGINAL_6_TYPES if ct not in shared_types]
    if missing_6:
        print(f"     Missing: {', '.join(missing_6)}")

    cd4_n = rescue_stats.get("n_classified_cd4", "N/A")
    cd8_n = rescue_stats.get("n_classified_cd8", "N/A")
    ambig_n = rescue_stats.get("n_ambiguous_both", "N/A")
    neither_n = rescue_stats.get("n_neither", "N/A")
    print(f"\n  5. CD4+ T rescue: {cd4_n} classified CD4+, "
          f"{cd8_n} CD8+, {ambig_n} ambiguous, {neither_n} neither")

    hepatocyte_n = mca_per_type.get("hepatocyte", 0)
    hep_status = "PASS" if hepatocyte_n >= 200 else "BORDERLINE"
    print(f"\n  6. MCA hepatocyte cell count: {hepatocyte_n} ({hep_status} at >=200 threshold)")

    excluded = results_json.get("excluded_tissues", [])
    print(f"\n  7. Tissues excluded from MCA due to gene coverage: "
          f"{len(excluded)}")
    for name, reason in excluded:
        print(f"     - {name}: {reason}")

    print(f"\n  8. T1-B obs/null ratio: "
          f"T1-A={t1b_stats['t1a_obs_null']:.4f} vs "
          f"Tabula 35-type={t1b_stats['tabula_35_obs_null']:.4f} vs "
          f"Tabula 6-type={t1b_stats['tabula_6_obs_null']:.4f}")

    print("\n" + "=" * 70)
    print("DONE")
    print("=" * 70)


if __name__ == "__main__":
    main()
