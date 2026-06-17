#!/usr/bin/env python3
"""
CellWarp — ISSUE-092 Diagnosis: Sun2023 Endothelial Tissue Stratification

Diagnoses whether the Sun2023 residual inversion (endothelial + hepatocyte dominate
Procrustes residuals, opposite to primary result) is an annotation granularity artifact.

Hypothesis: Sun2023 pools endothelial cells from 5 tissues (liver 3283, kidney 1330,
aorta 1188, lung 451, intestine 387) into one centroid, while Tabula Sapiens pools from
different tissues (myometrium 265, adipose 328, muscle 228, pancreas 130, liver 101,
aorta 128). This tissue composition mismatch may inflate the endothelial residual.

Biology
-------
Endothelial cells are highly tissue-specialized. Liver sinusoidal endothelial cells
(LSECs) are fenestrated and express distinct scavenger receptors; kidney glomerular
endothelial cells have unique filtration features; aortic endothelial cells face
hemodynamic shear stress. A pooled centroid across these subtypes will differ from
a centroid pooled from a different tissue mixture, even if the individual subtypes
are conserved across species.

Math
----
For each tissue with ≥200 endothelial cells, compute a separate centroid in the
16,959-gene ortholog space. PCA-reduce alongside the 15 Procrustes type centroids
and Tabula's endothelial centroid. Compute Euclidean distances in PCA space.
Replace pooled endothelial with the best tissue-matched centroid and re-run Procrustes.
If endothelial SSR drops from current level to ≤20%, tissue pooling explains the
inversion.

Output
------
  output/validation/sun2023_issue092_diagnosis/
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from scipy.stats import spearmanr
from sklearn.decomposition import PCA

from cellwarp.procrustes import (
    compute_centroids,
    pca_reduce_centroids,
    procrustes_align,
    permutation_test,
    compute_residual_vectors,
    RANDOM_SEED,
)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

OUTPUT_DIR = Path("output/validation/sun2023_issue092_diagnosis")
SUN2023_H5AD = Path("data/replication/sun2023/sun2023_yc.h5ad")
TABULA_HUMAN_PATH = Path("data/phase1/human_qc.h5ad")
TABULA_CENTROIDS_PATH = Path("output/phase2/scaled_35types/centroids_human_35.csv")
PRIMARY_RESULTS_PATH = Path("output/phase2/scaled_35types/procrustes_results_35.json")
RESIDUALS_RANKED_PATH = Path("output/phase2/scaled_35types/residuals_ranked.csv")
EXPANDED_RESULTS_PATH = Path("output/validation/sun2023_replication_expanded/sun2023_expanded.json")

N_PERMUTATIONS = 10_000
VARIANCE_THRESHOLD = 0.95
MIN_ENDO_CELLS = 200
SSR_THRESHOLD = 0.20  # Task 3 threshold: endothelial must drop to ≤20% SSR


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load expanded result for reference
    with open(EXPANDED_RESULTS_PATH) as f:
        expanded = json.load(f)
    expanded_types = expanded["procrustes"]["cell_types"]
    expanded_residuals = {
        ct: expanded["procrustes"]["per_type_residuals"][ct]["magnitude"]
        for ct in expanded_types
    }
    expanded_total_ssr = sum(v**2 for v in expanded_residuals.values())
    expanded_endo_ssr_pct = expanded_residuals["endothelial cell"]**2 / expanded_total_ssr
    expanded_hep_ssr_pct = expanded_residuals["hepatocyte"]**2 / expanded_total_ssr

    print("=" * 70)
    print("ISSUE-092 DIAGNOSIS: Sun2023 Endothelial Tissue Stratification")
    print("=" * 70)
    print(f"\nExpanded result (n=15) reference:")
    print(f"  Endothelial residual: {expanded_residuals['endothelial cell']:.3f} "
          f"({expanded_endo_ssr_pct:.1%} SSR)")
    print(f"  Hepatocyte residual:  {expanded_residuals['hepatocyte']:.3f} "
          f"({expanded_hep_ssr_pct:.1%} SSR)")

    # ==================================================================
    # TASK 0: Load data and verify tissue metadata
    # ==================================================================
    print("\n" + "=" * 70)
    print("TASK 0: Load Sun2023 h5ad and verify tissue metadata")
    print("=" * 70)

    import anndata as ad

    sun = ad.read_h5ad(SUN2023_H5AD)
    print(f"  Sun2023 shape: {sun.n_obs:,} cells × {sun.n_vars:,} genes")
    print(f"  Obs columns: {list(sun.obs.columns)}")
    print(f"  Tissue column: {'PRESENT' if 'tissue' in sun.obs.columns else 'MISSING'}")

    if "tissue" not in sun.obs.columns:
        print("\n  *** STOP: No tissue column. Cannot proceed. ***")
        return

    # Endothelial tissue breakdown
    endo_mask = sun.obs["cell_type"] == "endothelial cell"
    endo_tissue = sun.obs.loc[endo_mask, "tissue"].value_counts()
    print(f"\n  Endothelial cells total: {endo_mask.sum():,}")
    print("  By tissue:")
    for tissue, n in endo_tissue.items():
        status = "≥200 PASS" if n >= MIN_ENDO_CELLS else "SKIP"
        print(f"    {tissue:<20} {n:>5,}  ({status})")

    # Hepatocyte tissue breakdown
    hep_mask = sun.obs["cell_type"] == "hepatocyte"
    hep_tissue = sun.obs.loc[hep_mask, "tissue"].value_counts()
    print(f"\n  Hepatocyte cells total: {hep_mask.sum():,}")
    print("  By tissue:")
    for tissue, n in hep_tissue.items():
        print(f"    {tissue:<20} {n:>5,}")

    # Load Tabula human data for endothelial tissue composition
    print("\n  Loading Tabula human for endothelial tissue comparison...")
    tabula_h = ad.read_h5ad(TABULA_HUMAN_PATH)
    tabula_endo_mask = tabula_h.obs["cell_type"] == "endothelial cell"
    tabula_endo_tissue = tabula_h.obs.loc[tabula_endo_mask, "tissue"].value_counts()
    print(f"  Tabula endothelial cells total: {tabula_endo_mask.sum():,}")
    print(f"  Top tissues:")
    for tissue, n in tabula_endo_tissue.head(10).items():
        print(f"    {tissue:<35} {n:>5,}")

    # ==================================================================
    # TASK 1: Endothelial tissue stratification
    # ==================================================================
    print("\n" + "=" * 70)
    print("TASK 1: Endothelial Tissue Stratification")
    print("=" * 70)

    # Step 1a: Split by tissue, report counts
    eligible_tissues = [t for t, n in endo_tissue.items() if n >= MIN_ENDO_CELLS]
    print(f"\n  Step 1a: Tissues with ≥{MIN_ENDO_CELLS} endothelial cells: "
          f"{len(eligible_tissues)}")
    for t in eligible_tissues:
        print(f"    {t}: {endo_tissue[t]:,} cells")

    # Step 1b: Compute per-tissue endothelial centroids in gene space
    print(f"\n  Step 1b: Computing per-tissue endothelial centroids...")

    # Sun2023 data is already normalized (CPM+log1p) in human Ensembl ID space
    X_sun = sun.X
    if sp.issparse(X_sun):
        X_sun = X_sun.toarray()

    tissue_centroids = {}
    for tissue in eligible_tissues:
        mask = (sun.obs["cell_type"] == "endothelial cell") & (sun.obs["tissue"] == tissue)
        n_cells = mask.sum()
        centroid = np.mean(X_sun[mask.values], axis=0)
        tissue_centroids[tissue] = centroid
        print(f"    {tissue:<20} n={n_cells:>5,}  centroid norm={np.linalg.norm(centroid):.3f}")

    # Pooled centroid (what the expanded analysis used)
    pooled_mask = sun.obs["cell_type"] == "endothelial cell"
    pooled_centroid = np.mean(X_sun[pooled_mask.values], axis=0)
    print(f"    {'POOLED':<20} n={pooled_mask.sum():>5,}  centroid norm={np.linalg.norm(pooled_centroid):.3f}")

    # Load Tabula human endothelial centroid
    tabula_centroids = pd.read_csv(TABULA_CENTROIDS_PATH, index_col=0)
    tabula_endo = tabula_centroids.loc["endothelial cell"].values
    print(f"    Tabula endothelial centroid norm={np.linalg.norm(tabula_endo):.3f}")

    # Step 1c: PCA-reduce and compute distances
    print(f"\n  Step 1c: Computing distances in shared PCA space...")

    # Build the centroid matrix: all 15 expanded types (Sun2023) + 15 Tabula types
    # Plus extra tissue-specific endothelial centroids for distance calculation
    sun_type_centroids = {}
    for ct in expanded_types:
        ct_mask = sun.obs["cell_type"] == ct
        if ct_mask.sum() > 0:
            sun_type_centroids[ct] = np.mean(X_sun[ct_mask.values], axis=0)

    # Combine Sun2023 type centroids + Tabula type centroids for PCA fitting
    shared_types = sorted(set(sun_type_centroids.keys()) & set(tabula_centroids.index))
    n_types = len(shared_types)

    sun_centroid_matrix = np.array([sun_type_centroids[ct] for ct in shared_types])
    tabula_centroid_matrix = np.array([tabula_centroids.loc[ct].values for ct in shared_types])

    # Stack for PCA fitting (same as Procrustes pipeline)
    combined_for_pca = np.vstack([tabula_centroid_matrix, sun_centroid_matrix])
    print(f"    PCA fitting on {combined_for_pca.shape[0]} centroids × {combined_for_pca.shape[1]} genes")

    # Center and PCA
    center = combined_for_pca.mean(axis=0)
    pca = PCA(n_components=min(combined_for_pca.shape[0] - 1, combined_for_pca.shape[1]))
    pca.fit(combined_for_pca - center)

    # Determine components at 95% variance
    cumvar = np.cumsum(pca.explained_variance_ratio_)
    n_comp = int(np.searchsorted(cumvar, VARIANCE_THRESHOLD) + 1)
    n_comp = min(n_comp, pca.n_components_)
    print(f"    PCA: {n_comp} components at {VARIANCE_THRESHOLD:.0%} variance "
          f"(cumulative: {cumvar[n_comp-1]:.3f})")

    # Project tissue-specific endothelial centroids into PCA space
    tabula_endo_pca = pca.transform((tabula_endo - center).reshape(1, -1))[:, :n_comp][0]
    pooled_pca = pca.transform((pooled_centroid - center).reshape(1, -1))[:, :n_comp][0]

    tissue_pca = {}
    for tissue, centroid in tissue_centroids.items():
        tissue_pca[tissue] = pca.transform((centroid - center).reshape(1, -1))[:, :n_comp][0]

    # Compute Euclidean distances to Tabula endothelial in PCA space
    print(f"\n  Euclidean distances to Tabula endothelial centroid (PCA, k={n_comp}):")
    distances = {}
    for tissue, pca_vec in tissue_pca.items():
        d = np.linalg.norm(pca_vec - tabula_endo_pca)
        distances[tissue] = d

    pooled_dist = np.linalg.norm(pooled_pca - tabula_endo_pca)
    distances["POOLED"] = pooled_dist

    # Sort by distance
    sorted_distances = sorted(distances.items(), key=lambda x: x[1])
    for tissue, d in sorted_distances:
        marker = " <<<" if d == min(distances.values()) else ""
        print(f"    {tissue:<20} distance = {d:>8.3f}{marker}")

    best_tissue = sorted_distances[0][0]
    best_distance = sorted_distances[0][1]
    print(f"\n  CLOSEST tissue match: {best_tissue} (d={best_distance:.3f})")
    print(f"  POOLED distance:      {pooled_dist:.3f}")
    print(f"  Ratio best/pooled:    {best_distance/pooled_dist:.3f}")

    # Also compute gene-space distances for comparison
    print(f"\n  Gene-space Euclidean distances (full 16,959-dim):")
    gene_distances = {}
    for tissue, centroid in tissue_centroids.items():
        d = np.linalg.norm(centroid - tabula_endo)
        gene_distances[tissue] = d
    gene_distances["POOLED"] = np.linalg.norm(pooled_centroid - tabula_endo)
    for tissue, d in sorted(gene_distances.items(), key=lambda x: x[1]):
        print(f"    {tissue:<20} distance = {d:>8.3f}")

    # ==================================================================
    # Step 1d: Recompute Procrustes with tissue-matched endothelial
    # ==================================================================
    print(f"\n  Step 1d: Recomputing Procrustes with tissue-matched endothelial...")

    if best_tissue == "POOLED":
        print("    Best tissue IS the pooled centroid — no improvement possible.")
        tissue_matched_improvement = False
    else:
        # Replace Sun2023 pooled endothelial centroid with best tissue-matched
        sun_centroids_tissuematched = sun_centroid_matrix.copy()
        endo_idx = shared_types.index("endothelial cell")

        # Compute best-tissue centroid in gene space
        best_tissue_centroid = tissue_centroids[best_tissue]
        sun_centroids_tissuematched[endo_idx] = best_tissue_centroid

        # Refit PCA on modified centroids
        combined_modified = np.vstack([tabula_centroid_matrix, sun_centroids_tissuematched])
        center_mod = combined_modified.mean(axis=0)
        pca_mod = PCA(n_components=min(combined_modified.shape[0] - 1, combined_modified.shape[1]))
        pca_mod.fit(combined_modified - center_mod)
        cumvar_mod = np.cumsum(pca_mod.explained_variance_ratio_)
        n_comp_mod = int(np.searchsorted(cumvar_mod, VARIANCE_THRESHOLD) + 1)
        n_comp_mod = min(n_comp_mod, pca_mod.n_components_)

        # Project into PCA space
        human_pca_mod = pca_mod.transform(tabula_centroid_matrix - center_mod)[:, :n_comp_mod]
        mouse_pca_mod = pca_mod.transform(sun_centroids_tissuematched - center_mod)[:, :n_comp_mod]

        print(f"    Modified PCA: {n_comp_mod} components at {VARIANCE_THRESHOLD:.0%} variance")

        # Run Procrustes
        result_mod = procrustes_align(human_pca_mod, mouse_pca_mod)

        # Compute residuals
        residuals_mod = compute_residual_vectors(result_mod, shared_types)
        residual_mags_mod = {ct: float(np.linalg.norm(residuals_mod[ct])) for ct in shared_types}

        total_ssr_mod = sum(v**2 for v in residual_mags_mod.values())
        endo_ssr_mod = residual_mags_mod["endothelial cell"]**2 / total_ssr_mod
        hep_ssr_mod = residual_mags_mod["hepatocyte"]**2 / total_ssr_mod

        print(f"\n    Tissue-matched Procrustes (n={n_types}, endothelial={best_tissue}):")
        print(f"      Distance: {result_mod.distance:.3f} (was {expanded['procrustes']['distance']:.3f})")
        print(f"      Scaling:  {result_mod.scaling:.3f} (was {expanded['procrustes']['scaling']:.3f})")
        print(f"\n    Residual comparison:")
        print(f"      {'Cell type':<50} {'Pooled':>8} {'Tissue-matched':>14} {'Δ':>8}")
        print(f"      {'-'*80}")
        for ct in sorted(residual_mags_mod, key=residual_mags_mod.get, reverse=True):
            pooled_mag = expanded_residuals.get(ct, 0)
            mod_mag = residual_mags_mod[ct]
            delta = mod_mag - pooled_mag
            pooled_pct = pooled_mag**2 / expanded_total_ssr * 100
            mod_pct = mod_mag**2 / total_ssr_mod * 100
            print(f"      {ct:<50} {pooled_mag:>7.3f} ({pooled_pct:>4.1f}%) "
                  f"{mod_mag:>7.3f} ({mod_pct:>4.1f}%) {delta:>+7.3f}")

        print(f"\n    Endothelial SSR: {expanded_endo_ssr_pct:.1%} → {endo_ssr_mod:.1%}")
        print(f"    Hepatocyte SSR:  {expanded_hep_ssr_pct:.1%} → {hep_ssr_mod:.1%}")
        print(f"    SSR threshold met (≤{SSR_THRESHOLD:.0%})? "
              f"{'YES' if endo_ssr_mod <= SSR_THRESHOLD else 'NO'}")

        tissue_matched_improvement = endo_ssr_mod <= SSR_THRESHOLD

    # ==================================================================
    # TASK 2: Hepatocyte check
    # ==================================================================
    print("\n" + "=" * 70)
    print("TASK 2: Hepatocyte Check")
    print("=" * 70)

    n_hep_tissues = len(hep_tissue[hep_tissue > 0])
    if n_hep_tissues == 1 and hep_tissue.index[0] == "liver":
        print(f"\n  CONFIRMED: All {hep_mask.sum():,} hepatocytes from liver only.")
        print(f"  Hepatocyte {expanded_hep_ssr_pct:.1%} SSR is NOT a tissue pooling artifact.")
        print(f"  Alternative explanations:")
        print(f"    1. Alb-rescue annotation bias (script 16/17: liver cells with Alb>0 → hepatocyte)")
        print(f"    2. Genuine biological difference in hepatocyte identity between")
        print(f"       Sun2023 (C57BL/6J 2-month male, 10x 3' v3) and Tabula Sapiens")
        print(f"    3. Sun2023 captures hepatocyte subtypes (periportal/pericentral zonation)")
        print(f"       that shift the centroid relative to Tabula's more uniform sampling")
        hep_explanation = "NOT pooling artifact — single tissue (liver). Likely annotation or biology."
    else:
        print(f"\n  UNEXPECTED: Hepatocytes span {n_hep_tissues} tissues:")
        for tissue, n in hep_tissue.items():
            if n > 0:
                print(f"    {tissue}: {n}")
        hep_explanation = f"Hepatocytes from {n_hep_tissues} tissues — needs stratification."

    # ==================================================================
    # TASK 3: Revised Procrustes (if threshold met)
    # ==================================================================
    print("\n" + "=" * 70)
    print("TASK 3: Revised Procrustes with Tissue-Matched Endothelial")
    print("=" * 70)

    if best_tissue == "POOLED" or not tissue_matched_improvement:
        print(f"\n  Endothelial SSR {'equals pooled' if best_tissue == 'POOLED' else 'stays >' + str(int(SSR_THRESHOLD*100)) + '%'} "
              f"after tissue matching.")
        print(f"  CONCLUSION: Annotation granularity does NOT explain the endothelial inversion.")
        print(f"  The residual inversion is likely biological or protocol-driven.")
        print(f"  Flagging to advisor. ISSUE-092 status: UNRESOLVED.")
        revised_verdict = None
        revised_rho = None
        revised_rho_p = None
    else:
        print(f"\n  Endothelial SSR dropped to {endo_ssr_mod:.1%} (≤{SSR_THRESHOLD:.0%}). "
              f"Running full Procrustes with permutation test...")

        # Full permutation test on tissue-matched Procrustes
        p_mod, null_mod = permutation_test(
            human_pca_mod, mouse_pca_mod, N_PERMUTATIONS, RANDOM_SEED
        )
        obs_null_mod = result_mod.distance / np.median(null_mod)

        print(f"\n  Revised Procrustes results (tissue-matched endothelial):")
        print(f"    p = {p_mod:.4f}")
        print(f"    obs/null = {obs_null_mod:.3f}")
        print(f"    scaling = {result_mod.scaling:.3f}")
        print(f"    distance = {result_mod.distance:.3f}")
        print(f"    PCA components = {n_comp_mod}")

        # Rigidity ranking correlation with primary
        primary_residuals_df = pd.read_csv(RESIDUALS_RANKED_PATH)
        primary_mag_dict = dict(
            zip(primary_residuals_df["cell_type"], primary_residuals_df["residual_magnitude"])
        )
        matched = sorted(set(residual_mags_mod.keys()) & set(primary_mag_dict.keys()))
        if len(matched) >= 4:
            sun_mags = [residual_mags_mod[ct] for ct in matched]
            primary_mags = [primary_mag_dict[ct] for ct in matched]
            rho_mod, rho_p_mod = spearmanr(sun_mags, primary_mags)
        else:
            rho_mod, rho_p_mod = float("nan"), float("nan")

        print(f"    Spearman ρ = {rho_mod:.3f}, p = {rho_p_mod:.4f} (n={len(matched)})")

        # Compare to original expanded
        print(f"\n  Comparison to pooled expanded (n=15):")
        print(f"    {'Metric':<25} {'Pooled':>12} {'Tissue-matched':>14} {'Change':>10}")
        print(f"    {'-'*61}")
        print(f"    {'p-value':<25} {expanded['procrustes']['p_value']:>12.4f} {p_mod:>14.4f}")
        print(f"    {'obs/null':<25} {expanded['procrustes']['obs_null_ratio']:>12.3f} {obs_null_mod:>14.3f} "
              f"{obs_null_mod - expanded['procrustes']['obs_null_ratio']:>+10.3f}")
        print(f"    {'scaling':<25} {expanded['procrustes']['scaling']:>12.3f} {result_mod.scaling:>14.3f}")
        print(f"    {'Spearman ρ':<25} {expanded['rigidity_ranking']['rho']:>12.3f} {rho_mod:>14.3f} "
              f"{rho_mod - expanded['rigidity_ranking']['rho']:>+10.3f}")
        print(f"    {'ρ p-value':<25} {expanded['rigidity_ranking']['p_value']:>12.4f} {rho_p_mod:>14.4f}")

        revised_verdict = "PASS" if rho_mod >= 0.50 and rho_p_mod < 0.05 else "PARTIAL"
        revised_rho = float(rho_mod) if not np.isnan(rho_mod) else None
        revised_rho_p = float(rho_p_mod) if not np.isnan(rho_p_mod) else None

        # Save null distribution
        np.save(OUTPUT_DIR / "null_distribution_tissuematched.npy", null_mod)

        # Plot: null distribution
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        ax = axes[0]
        ax.hist(null_mod, bins=50, alpha=0.7, color="teal", edgecolor="white")
        ax.axvline(result_mod.distance, color="red", linewidth=2,
                   label=f"Observed (d={result_mod.distance:.2f})")
        ax.set_title(f"Tissue-matched Procrustes ({n_types} types)\n"
                     f"p={p_mod:.4f}, obs/null={obs_null_mod:.3f}")
        ax.set_xlabel("Procrustes distance")
        ax.set_ylabel("Count")
        ax.legend()

        ax = axes[1]
        # Scatter: tissue-matched vs primary residuals
        if len(matched) >= 4:
            ax.scatter(
                [primary_mag_dict[ct] for ct in matched],
                [residual_mags_mod[ct] for ct in matched],
                s=60, c="teal", edgecolors="darkslategray", linewidths=0.5, zorder=3,
            )
            for ct in matched:
                short = ct[:22] + "..." if len(ct) > 22 else ct
                ax.annotate(short, (primary_mag_dict[ct], residual_mags_mod[ct]),
                            fontsize=6, ha="left", va="bottom",
                            xytext=(4, 4), textcoords="offset points")
            lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
                    max(ax.get_xlim()[1], ax.get_ylim()[1])]
            ax.plot(lims, lims, "k--", alpha=0.3, zorder=1)
            ax.set_xlim(lims)
            ax.set_ylim(lims)
            ax.set_xlabel("Primary 35-type residual magnitude")
            ax.set_ylabel("Tissue-matched residual magnitude")
            ax.set_title(f"Rigidity: tissue-matched vs primary\n"
                         f"ρ={rho_mod:.3f}, p={rho_p_mod:.4f}")

        plt.suptitle("ISSUE-092: Tissue-Matched Endothelial Procrustes", fontweight="bold")
        plt.tight_layout()
        fig.savefig(OUTPUT_DIR / "procrustes_tissuematched.png", dpi=150, bbox_inches="tight")
        plt.close()

    # ==================================================================
    # Plot: Endothelial tissue distances
    # ==================================================================
    fig, ax = plt.subplots(figsize=(8, 5))
    tissue_names = [t for t, _ in sorted_distances]
    dist_vals = [d for _, d in sorted_distances]
    colors = ["teal" if t != "POOLED" else "salmon" for t in tissue_names]
    bars = ax.barh(range(len(tissue_names)), dist_vals, color=colors, edgecolor="black", linewidth=0.5)
    ax.set_yticks(range(len(tissue_names)))
    ax.set_yticklabels([f"{t} ({endo_tissue.get(t, endo_mask.sum()):,})" if t != "POOLED"
                        else f"POOLED ({endo_mask.sum():,})" for t in tissue_names])
    ax.set_xlabel(f"Euclidean distance to Tabula endothelial (PCA, k={n_comp})")
    ax.set_title("Sun2023 Endothelial: Distance to Tabula by Tissue of Origin")
    for bar, val in zip(bars, dist_vals):
        ax.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height()/2,
                f"{val:.1f}", va="center", fontsize=9)
    ax.invert_yaxis()
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "endothelial_tissue_distances.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\n  Plot: endothelial_tissue_distances.png")

    # ==================================================================
    # Plot: SSR comparison (pooled vs tissue-matched)
    # ==================================================================
    if best_tissue != "POOLED" and tissue_matched_improvement is not None:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # Pooled SSR
        ax = axes[0]
        sorted_pooled = sorted(expanded_residuals.items(), key=lambda x: x[1]**2, reverse=True)
        ct_names_p = [ct[:25] for ct, _ in sorted_pooled]
        ssr_pcts_p = [v**2 / expanded_total_ssr * 100 for _, v in sorted_pooled]
        ax.barh(range(len(ct_names_p)), ssr_pcts_p, color="coral", edgecolor="black", linewidth=0.5)
        ax.set_yticks(range(len(ct_names_p)))
        ax.set_yticklabels(ct_names_p, fontsize=8)
        ax.set_xlabel("% of total SSR")
        ax.set_title("Pooled Endothelial (Original)")
        ax.invert_yaxis()
        ax.axvline(20, color="red", linestyle="--", alpha=0.5, label="20% threshold")
        ax.legend(fontsize=8)

        # Tissue-matched SSR
        if tissue_matched_improvement or endo_ssr_mod is not None:
            ax = axes[1]
            sorted_mod = sorted(residual_mags_mod.items(), key=lambda x: x[1]**2, reverse=True)
            ct_names_m = [ct[:25] for ct, _ in sorted_mod]
            ssr_pcts_m = [v**2 / total_ssr_mod * 100 for _, v in sorted_mod]
            ax.barh(range(len(ct_names_m)), ssr_pcts_m, color="teal", edgecolor="black", linewidth=0.5)
            ax.set_yticks(range(len(ct_names_m)))
            ax.set_yticklabels(ct_names_m, fontsize=8)
            ax.set_xlabel("% of total SSR")
            ax.set_title(f"Tissue-Matched ({best_tissue})")
            ax.invert_yaxis()
            ax.axvline(20, color="red", linestyle="--", alpha=0.5, label="20% threshold")
            ax.legend(fontsize=8)

        plt.suptitle("ISSUE-092: SSR Redistribution After Tissue Matching", fontweight="bold")
        plt.tight_layout()
        fig.savefig(OUTPUT_DIR / "ssr_comparison.png", dpi=150, bbox_inches="tight")
        plt.close()
        print(f"  Plot: ssr_comparison.png")

    # ==================================================================
    # Save results JSON
    # ==================================================================
    print("\n" + "=" * 70)
    print("Saving results")
    print("=" * 70)

    results = {
        "diagnostic": "ISSUE-092: Sun2023 endothelial tissue stratification diagnosis",
        "date": "2026-03-15",
        "task_0_metadata": {
            "tissue_column": "PRESENT",
            "endothelial_by_tissue": {t: int(n) for t, n in endo_tissue.items()},
            "endothelial_tissues_eligible": eligible_tissues,
            "hepatocyte_by_tissue": {t: int(n) for t, n in hep_tissue.items()},
            "tabula_endothelial_top_tissues": {
                t: int(n) for t, n in tabula_endo_tissue.head(10).items()
            },
        },
        "task_1_endothelial_stratification": {
            "pca_components": n_comp,
            "distances_to_tabula_endothelial_pca": {
                t: float(d) for t, d in sorted_distances
            },
            "distances_to_tabula_endothelial_genespace": {
                t: float(d) for t, d in sorted(gene_distances.items(), key=lambda x: x[1])
            },
            "best_tissue_match": best_tissue,
            "best_distance": float(best_distance),
            "pooled_distance": float(pooled_dist),
            "distance_ratio_best_over_pooled": float(best_distance / pooled_dist),
        },
        "task_1d_tissue_matched_procrustes": {
            "endothelial_ssr_pooled_pct": float(expanded_endo_ssr_pct),
            "endothelial_ssr_tissuematched_pct": float(endo_ssr_mod) if best_tissue != "POOLED" else None,
            "ssr_threshold": SSR_THRESHOLD,
            "threshold_met": tissue_matched_improvement if best_tissue != "POOLED" else False,
            "distance": float(result_mod.distance) if best_tissue != "POOLED" else None,
            "scaling": float(result_mod.scaling) if best_tissue != "POOLED" else None,
            "per_type_residuals": {
                ct: {"magnitude_pooled": expanded_residuals.get(ct, 0),
                     "magnitude_tissuematched": residual_mags_mod.get(ct, 0)}
                for ct in shared_types
            } if best_tissue != "POOLED" else None,
        },
        "task_2_hepatocyte": {
            "single_tissue": n_hep_tissues == 1,
            "tissue": hep_tissue.index[0] if n_hep_tissues >= 1 else None,
            "n_cells": int(hep_mask.sum()),
            "ssr_pct": float(expanded_hep_ssr_pct),
            "explanation": hep_explanation,
        },
        "task_3_revised_procrustes": {
            "ran": tissue_matched_improvement if best_tissue != "POOLED" else False,
        },
        "conclusion": {},
    }

    # Add Task 3 results if run
    if revised_verdict is not None:
        results["task_3_revised_procrustes"].update({
            "p_value": float(p_mod),
            "obs_null_ratio": float(obs_null_mod),
            "scaling": float(result_mod.scaling),
            "distance": float(result_mod.distance),
            "pca_components": n_comp_mod,
            "rho_vs_primary": revised_rho,
            "rho_p_value": revised_rho_p,
            "n_matched_types": len(matched),
            "verdict": revised_verdict,
        })

    # Determine ISSUE-092 status
    if tissue_matched_improvement and revised_verdict == "PASS":
        issue_status = "RESOLVED"
        conclusion = (
            f"Tissue pooling explains endothelial inversion. Using {best_tissue} "
            f"endothelial centroid drops SSR from {expanded_endo_ssr_pct:.1%} to "
            f"{endo_ssr_mod:.1%}. Revised Procrustes: ρ={revised_rho:.3f}, p={revised_rho_p:.4f}."
        )
    elif tissue_matched_improvement and revised_verdict == "PARTIAL":
        issue_status = "PARTIALLY RESOLVED"
        conclusion = (
            f"Tissue pooling partially explains endothelial inversion. Using {best_tissue} "
            f"endothelial centroid drops SSR from {expanded_endo_ssr_pct:.1%} to "
            f"{endo_ssr_mod:.1%}. But rigidity ranking still does not replicate "
            f"(ρ={revised_rho:.3f}, p={revised_rho_p:.4f}). Hepatocyte SSR "
            f"({hep_ssr_mod:.1%}) also remains elevated — not a pooling artifact."
        )
    elif best_tissue != "POOLED" and not tissue_matched_improvement:
        issue_status = "UNRESOLVED"
        conclusion = (
            f"Tissue matching does NOT explain the endothelial inversion. Best tissue "
            f"({best_tissue}) reduces distance but endothelial SSR remains "
            f"{endo_ssr_mod:.1%} (>{SSR_THRESHOLD:.0%} threshold). "
            f"The residual inversion is biological or protocol-driven, not an "
            f"annotation granularity artifact."
        )
    else:
        issue_status = "UNRESOLVED"
        conclusion = "Pooled centroid was already closest. No improvement possible."

    results["conclusion"] = {
        "issue_092_status": issue_status,
        "summary": conclusion,
        "endothelial_inversion_explained_by_pooling": tissue_matched_improvement if best_tissue != "POOLED" else False,
        "hepatocyte_inversion_explained_by_pooling": False,
    }

    with open(OUTPUT_DIR / "issue092_diagnosis.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"  Results: {OUTPUT_DIR / 'issue092_diagnosis.json'}")

    # ==================================================================
    # SUMMARY
    # ==================================================================
    print("\n" + "=" * 70)
    print("ISSUE-092 DIAGNOSIS SUMMARY")
    print("=" * 70)
    print(f"\n  Tissue metadata: AVAILABLE")
    print(f"\n  Endothelial tissue breakdown (Sun2023):")
    for t in eligible_tissues:
        print(f"    {t:<20} {endo_tissue[t]:>5,} cells")
    print(f"\n  Tabula endothelial top tissues: myometrium (265), adipose (328), "
          f"muscle (228), pancreas (130)")
    print(f"  Sun2023 endothelial top tissues: liver (3283), kidney (1330), aorta (1188)")
    print(f"  >>> SEVERE tissue composition mismatch <<<")
    print(f"\n  Best tissue match: {best_tissue} (distance {best_distance:.1f} vs "
          f"pooled {pooled_dist:.1f})")
    print(f"\n  Endothelial SSR: {expanded_endo_ssr_pct:.1%} → "
          f"{endo_ssr_mod:.1%}" if best_tissue != "POOLED" else
          f"  Endothelial SSR: {expanded_endo_ssr_pct:.1%} (no change)")
    print(f"  Hepatocyte SSR: {expanded_hep_ssr_pct:.1%} (single tissue, not pooling artifact)")
    print(f"\n  ISSUE-092 STATUS: {issue_status}")
    print(f"  {conclusion}")


if __name__ == "__main__":
    main()
