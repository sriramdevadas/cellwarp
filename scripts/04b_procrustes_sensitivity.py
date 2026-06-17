#!/usr/bin/env python3
"""
CellWarp — Phase 2 Procrustes Sensitivity Analysis

Reruns Procrustes analysis with cell type subsets to test whether the significant
alignment (p=0.0035 on all 6 types) is driven primarily by macrophages and/or
B cells, which together account for 84.6% of the total Procrustes SSR.

Analyses:
    A. All 6 cell types (reference — loads existing results)
    B. 5 cell types: exclude macrophage (60.7% of SSR)
    C. 4 cell types: exclude macrophage + B cell (84.6% of SSR combined)

Biology
-------
If the Procrustes alignment remains significant after removing the dominant
residual contributors, the geometric correspondence between species is robust
and not driven by one or two outlier cell types. If significance is lost, the
"D'Arcy Thompson grid" may only apply to a subset of cell types.

Inputs:
    output/phase2/centroids_human.csv
    output/phase2/centroids_mouse.csv

Outputs:
    output/phase2/sensitivity/5types_no_macrophage.json
    output/phase2/sensitivity/4types_no_macrophage_bcell.json
    output/phase2/sensitivity/comparison_table.csv

Usage:
    python scripts/04b_procrustes_sensitivity.py
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import numpy as np
import pandas as pd

from cellwarp.procrustes import (
    RANDOM_SEED,
    compute_residual_vectors,
    map_residuals_to_genes,
    pca_reduce_centroids,
    permutation_test,
    procrustes_align,
    save_results,
)


def run_subset_analysis(
    human_centroids: pd.DataFrame,
    mouse_centroids: pd.DataFrame,
    gene_names: list[str],
    subset_name: str,
    cell_types_to_keep: list[str],
    output_path: Path,
) -> dict:
    """Run full Procrustes pipeline on a subset of cell types."""
    print(f"\n  Filtering to {len(cell_types_to_keep)} cell types:")
    for ct in cell_types_to_keep:
        print(f"    - {ct}")

    h_sub = human_centroids.loc[cell_types_to_keep]
    m_sub = mouse_centroids.loc[cell_types_to_keep]

    # PCA on subset
    human_pca, mouse_pca, pca_model, cell_types = pca_reduce_centroids(
        h_sub, m_sub
    )

    # Procrustes
    result = procrustes_align(human_pca, mouse_pca)

    # Permutation test
    p_value, null_dist = permutation_test(human_pca, mouse_pca)

    # Residuals
    residuals = compute_residual_vectors(result, cell_types)

    # Gene mapping
    top_genes = map_residuals_to_genes(residuals, pca_model, gene_names)

    # Save
    pca_info = {
        "n_components": int(pca_model.n_components_),
        "variance_explained_per_component": (
            pca_model.explained_variance_ratio_.tolist()
        ),
        "cumulative_variance_explained": float(
            sum(pca_model.explained_variance_ratio_)
        ),
        "n_genes_input": len(gene_names),
        "subset_name": subset_name,
        "cell_types_excluded": sorted(
            set(human_centroids.index) - set(cell_types_to_keep)
        ),
    }

    save_results(
        result=result,
        p_value=p_value,
        null_distribution=null_dist,
        residuals=residuals,
        top_genes=top_genes,
        cell_types=cell_types,
        pca_info=pca_info,
        output_path=output_path,
    )

    return {
        "subset": subset_name,
        "n_types": len(cell_types_to_keep),
        "n_pca": int(pca_model.n_components_),
        "variance_explained": float(sum(pca_model.explained_variance_ratio_)) * 100,
        "distance": float(result.distance),
        "ssr": float(result.distance_squared),
        "scaling": float(result.scaling),
        "p_value": float(p_value),
        "significant_001": p_value < 0.01,
        "null_median": float(np.median(null_dist)),
    }


def main() -> None:
    output_dir = Path("./output/phase2/sensitivity")
    output_dir.mkdir(parents=True, exist_ok=True)

    t_start = time.time()

    # ==================================================================
    # Load centroids (from main analysis)
    # ==================================================================
    print("=" * 70)
    print("PHASE 2 — Procrustes Sensitivity Analysis")
    print("=" * 70)

    centroids_dir = Path("./output/phase2")
    human_centroids = pd.read_csv(centroids_dir / "centroids_human.csv", index_col=0)
    mouse_centroids = pd.read_csv(centroids_dir / "centroids_mouse.csv", index_col=0)

    print(f"\n  Loaded centroids: {human_centroids.shape[0]} types × "
          f"{human_centroids.shape[1]:,} genes")

    # Load gene names from main results for consistent reporting
    import anndata as ad
    human = ad.read_h5ad("./data/phase1/human_qc.h5ad")
    gene_names = human.var["feature_name"].tolist()
    del human  # free memory

    all_types = sorted(human_centroids.index.tolist())
    print(f"  Cell types: {all_types}")

    # ==================================================================
    # Load reference results (all 6 types)
    # ==================================================================
    print("\n" + "=" * 70)
    print("REFERENCE: All 6 cell types (from main analysis)")
    print("=" * 70)

    with open(centroids_dir / "procrustes_results.json") as f:
        ref_results = json.load(f)

    ref_summary = {
        "subset": "All 6 types",
        "n_types": 6,
        "n_pca": ref_results["pca"]["n_components"],
        "variance_explained": ref_results["pca"]["cumulative_variance_explained"] * 100,
        "distance": ref_results["procrustes"]["distance"],
        "ssr": ref_results["procrustes"]["distance_squared"],
        "scaling": ref_results["procrustes"]["scaling"],
        "p_value": ref_results["permutation_test"]["p_value"],
        "significant_001": ref_results["permutation_test"]["p_value"] < 0.01,
        "null_median": ref_results["permutation_test"]["null_distribution_summary"]["median"],
    }

    print(f"\n  p-value: {ref_summary['p_value']:.6f}")
    print(f"  Procrustes distance: {ref_summary['distance']:.4f}")

    # ==================================================================
    # Analysis B: 5 types (no macrophage)
    # ==================================================================
    print("\n" + "=" * 70)
    print("SENSITIVITY B: 5 cell types (exclude macrophage)")
    print("=" * 70)

    types_no_macro = [ct for ct in all_types if ct != "macrophage"]
    summary_b = run_subset_analysis(
        human_centroids, mouse_centroids, gene_names,
        subset_name="5 types (no macrophage)",
        cell_types_to_keep=types_no_macro,
        output_path=output_dir / "5types_no_macrophage.json",
    )

    # ==================================================================
    # Analysis C: 4 types (no macrophage, no B cell)
    # ==================================================================
    print("\n" + "=" * 70)
    print("SENSITIVITY C: 4 cell types (exclude macrophage + B cell)")
    print("=" * 70)

    types_no_macro_bcell = [
        ct for ct in all_types
        if ct not in ("macrophage", "B cell")
    ]
    summary_c = run_subset_analysis(
        human_centroids, mouse_centroids, gene_names,
        subset_name="4 types (no macrophage, B cell)",
        cell_types_to_keep=types_no_macro_bcell,
        output_path=output_dir / "4types_no_macrophage_bcell.json",
    )

    # ==================================================================
    # Comparison table
    # ==================================================================
    print("\n" + "=" * 70)
    print("COMPARISON TABLE")
    print("=" * 70)

    summaries = [ref_summary, summary_b, summary_c]
    comparison_df = pd.DataFrame(summaries)
    comparison_df = comparison_df[[
        "subset", "n_types", "n_pca", "variance_explained",
        "distance", "ssr", "scaling", "p_value", "significant_001",
        "null_median",
    ]]

    # Save CSV
    comparison_df.to_csv(output_dir / "comparison_table.csv", index=False)

    # Print formatted table
    print(f"\n  {'Subset':<38} {'Types':>5} {'PCs':>4} {'Var%':>6} "
          f"{'Distance':>10} {'p-value':>10} {'Sig?':>5} {'Null med':>10}")
    print(f"  {'-' * 95}")

    for s in summaries:
        sig = "YES" if s["significant_001"] else "NO"
        print(
            f"  {s['subset']:<38} {s['n_types']:>5} {s['n_pca']:>4} "
            f"{s['variance_explained']:>5.1f}% {s['distance']:>10.4f} "
            f"{s['p_value']:>10.6f} {sig:>5} {s['null_median']:>10.4f}"
        )

    print(f"\n  Key observations:")
    if summary_b["significant_001"]:
        print(f"  - Removing macrophage: p={summary_b['p_value']:.4f} — "
              f"STILL SIGNIFICANT. Alignment not driven by macrophage alone.")
    else:
        print(f"  - Removing macrophage: p={summary_b['p_value']:.4f} — "
              f"NOT SIGNIFICANT. Macrophage may be driving the signal.")

    if summary_c["significant_001"]:
        print(f"  - Removing macrophage + B cell: p={summary_c['p_value']:.4f} — "
              f"STILL SIGNIFICANT. Robust geometric structure.")
    else:
        print(f"  - Removing macrophage + B cell: p={summary_c['p_value']:.4f} — "
              f"NOT SIGNIFICANT. Signal concentrated in macrophage/B cell.")

    # Note on interpretability
    print(f"\n  Note: With 4 cell types, only 4!=24 unique permutations exist.")
    print(f"  Minimum achievable p-value with 4 types: "
          f"{1 / (10_000 + 1):.6f} (Monte Carlo) or {1 / 24:.4f} (exact).")

    t_total = time.time() - t_start
    print(f"\n  Output saved to: {output_dir}/")
    print(f"  Total runtime: {t_total:.1f}s")


if __name__ == "__main__":
    main()
