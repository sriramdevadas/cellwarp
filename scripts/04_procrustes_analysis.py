#!/usr/bin/env python3
"""
CellWarp — Phase 2 Procrustes Analysis Script

Loads QC'd datasets from Phase 1, computes cell type centroids in gene space,
reduces dimensionality via PCA, performs Procrustes alignment (mouse → human),
runs a permutation test for significance, and maps residual deformation vectors
back to gene space to identify cell-type-specific divergence genes.

Biology
-------
Each cell type occupies a position in high-dimensional gene expression space.
By computing the centroid (mean expression vector) per cell type, we reduce each
species to a constellation of 6 points. Procrustes analysis asks: is there a
rigid transformation (rotation + scaling) that maps the mouse constellation onto
the human constellation? If yes, cell types have maintained their relative
geometric relationships across 90 million years of evolution — a D'Arcy Thompson
"grid" in transcriptomic space.

Residual vectors (what the global transformation can't explain) reveal cell-type-
specific evolutionary divergence: genes that changed in one cell type more than
the global pattern predicts.

Pipeline
--------
1. Load normalized data (human_qc.h5ad, mouse_qc.h5ad)
2. Compute mean expression per cell type per species (6 centroids × 16,959 genes)
3. PCA on combined 12 centroids (retain ≥95% variance)
4. Procrustes alignment: center, rotate (SVD), scale mouse → human
5. Permutation test: 10,000 shuffled pairings (p < 0.01 for Phase 2 gate)
6. Compute residual deformation vectors per cell type
7. Map residuals back to gene space (top 20 genes per cell type)
8. Save all results to JSON + intermediate files

Inputs:
    data/phase1/human_qc.h5ad
    data/phase1/mouse_qc.h5ad

Outputs:
    output/phase2/centroids_human.csv       — Mean expression per cell type (human)
    output/phase2/centroids_mouse.csv       — Mean expression per cell type (mouse)
    output/phase2/pca_centroids.npz         — PCA-reduced centroids
    output/phase2/null_distribution.npy     — 10,000 permuted Procrustes distances
    output/phase2/procrustes_results.json   — All results (main output)

Usage:
    python scripts/04_procrustes_analysis.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import anndata as ad
import numpy as np

from cellwarp.procrustes import (
    compute_centroids,
    compute_residual_vectors,
    map_residuals_to_genes,
    pca_reduce_centroids,
    permutation_test,
    procrustes_align,
    save_results,
)


def main() -> None:
    data_dir = Path("./data/phase1")
    output_dir = Path("./output/phase2")
    output_dir.mkdir(parents=True, exist_ok=True)

    t_start = time.time()

    # ==================================================================
    # Load data
    # ==================================================================
    print("=" * 70)
    print("PHASE 2 — Procrustes Analysis")
    print("=" * 70)

    print("\n  Loading normalized datasets...")
    human = ad.read_h5ad(data_dir / "human_qc.h5ad")
    mouse = ad.read_h5ad(data_dir / "mouse_qc.h5ad")
    print(f"  Human: {human.n_obs:,} cells × {human.n_vars:,} genes")
    print(f"  Mouse: {mouse.n_obs:,} cells × {mouse.n_vars:,} genes")

    # Verify gene spaces are aligned
    assert list(human.var_names) == list(mouse.var_names), (
        "Gene spaces are not aligned! Human and mouse must share var_names."
    )
    n_genes = human.n_vars
    print(f"  Shared gene space confirmed: {n_genes:,} genes")

    # ==================================================================
    # Step 1: Compute centroids
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 1: Compute mean expression centroids")
    print("=" * 70)

    centroids_human_path = output_dir / "centroids_human.csv"
    centroids_mouse_path = output_dir / "centroids_mouse.csv"

    if centroids_human_path.exists() and centroids_mouse_path.exists():
        import pandas as pd

        print(f"  Checkpoint found: loading saved centroids")
        human_centroids = pd.read_csv(centroids_human_path, index_col=0)
        mouse_centroids = pd.read_csv(centroids_mouse_path, index_col=0)
        print(
            f"  Human centroids: {human_centroids.shape[0]} types × "
            f"{human_centroids.shape[1]:,} genes"
        )
        print(
            f"  Mouse centroids: {mouse_centroids.shape[0]} types × "
            f"{mouse_centroids.shape[1]:,} genes"
        )
    else:
        print("\n  Human centroids:")
        human_centroids = compute_centroids(human, "cell_type")

        print("\n  Mouse centroids:")
        mouse_centroids = compute_centroids(mouse, "cell_type")

        # Save intermediate (checkpoint)
        human_centroids.to_csv(centroids_human_path)
        mouse_centroids.to_csv(centroids_mouse_path)
        print(f"\n  Saved centroids to {output_dir}/centroids_{{human,mouse}}.csv")

    # ==================================================================
    # Step 2: PCA on combined centroids
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 2: PCA dimensionality reduction (95% variance)")
    print("=" * 70)

    human_pca, mouse_pca, pca_model, cell_types = pca_reduce_centroids(
        human_centroids, mouse_centroids
    )

    # Save intermediate PCA results
    np.savez(
        output_dir / "pca_centroids.npz",
        human=human_pca,
        mouse=mouse_pca,
        cell_types=np.array(cell_types),
        explained_variance_ratio=pca_model.explained_variance_ratio_,
    )
    print(f"  Saved PCA centroids to {output_dir}/pca_centroids.npz")

    # ==================================================================
    # Step 3: Procrustes alignment
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 3: Procrustes alignment (mouse → human)")
    print("=" * 70)

    result = procrustes_align(human_pca, mouse_pca)

    # ==================================================================
    # Step 4: Permutation test
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 4: Permutation test (10,000 iterations)")
    print("=" * 70)

    t_perm_start = time.time()
    p_value, null_dist = permutation_test(human_pca, mouse_pca)
    t_perm = time.time() - t_perm_start
    print(f"  Permutation test completed in {t_perm:.1f}s")

    # Save null distribution
    np.save(output_dir / "null_distribution.npy", null_dist)
    print(f"  Saved null distribution to {output_dir}/null_distribution.npy")

    # ==================================================================
    # Step 5: Residual deformation vectors
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 5: Residual deformation vectors")
    print("=" * 70)

    residuals = compute_residual_vectors(result, cell_types)

    # ==================================================================
    # Step 6: Map residuals to gene space
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 6: Map residuals back to gene space (top 20 genes)")
    print("=" * 70)

    # Use human gene symbols for interpretable output
    gene_names = human.var["feature_name"].tolist()
    top_genes = map_residuals_to_genes(residuals, pca_model, gene_names)

    # ==================================================================
    # Step 7: Save all results
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 7: Save results")
    print("=" * 70)

    pca_info = {
        "n_components": int(pca_model.n_components_),
        "variance_explained_per_component": (
            pca_model.explained_variance_ratio_.tolist()
        ),
        "cumulative_variance_explained": float(
            sum(pca_model.explained_variance_ratio_)
        ),
        "n_genes_input": n_genes,
    }

    save_results(
        result=result,
        p_value=p_value,
        null_distribution=null_dist,
        residuals=residuals,
        top_genes=top_genes,
        cell_types=cell_types,
        pca_info=pca_info,
        output_path=output_dir / "procrustes_results.json",
    )

    # ==================================================================
    # Summary
    # ==================================================================
    t_total = time.time() - t_start

    print("\n" + "=" * 70)
    print("PROCRUSTES ANALYSIS COMPLETE — SUMMARY")
    print("=" * 70)

    print(f"\n  Input: {n_genes:,} shared ortholog genes, 6 cell types, 2 species")
    print(f"  PCA components retained: {pca_model.n_components_}")
    print(
        f"  Cumulative variance explained: "
        f"{sum(pca_model.explained_variance_ratio_) * 100:.1f}%"
    )

    print(f"\n  Procrustes distance (SSR): {result.distance_squared:.6f}")
    print(f"  Procrustes distance (√SSR): {result.distance:.6f}")
    print(f"  Optimal scaling factor: {result.scaling:.6f}")

    print(f"\n  Permutation test p-value: {p_value:.6f}")
    gate_p = "PASS" if p_value < 0.01 else "FAIL"
    print(f"  Phase 2 gate (p < 0.01): {gate_p}")

    print(f"\n  Residual magnitudes:")
    for ct in cell_types:
        mag = np.linalg.norm(residuals[ct])
        print(f"    {ct:<45} {mag:.6f}")

    print(f"\n  Top divergence gene per cell type:")
    for ct in cell_types:
        top = top_genes[ct].iloc[0]
        direction = "↑mouse" if top["loading"] > 0 else "↓mouse"
        print(f"    {ct:<45} {top['gene']:<12} ({direction})")

    print(f"\n  Phase 2 Gate Check (partial):")
    print(f"    Procrustes p < 0.01: {gate_p} (p = {p_value:.6f})")
    print(f"    Interpretable residuals ≥3/6: REQUIRES BIOLOGICAL REVIEW")
    print(f"    Negative control: NOT YET RUN (separate script)")

    print(f"\n  Output files:")
    for f in sorted(output_dir.iterdir()):
        size = f.stat().st_size
        if size > 1024 * 1024:
            size_str = f"{size / 1024 / 1024:.1f} MB"
        elif size > 1024:
            size_str = f"{size / 1024:.1f} KB"
        else:
            size_str = f"{size} B"
        print(f"    {f.name:<35} {size_str:>10}")

    print(f"\n  Total runtime: {t_total:.1f}s")


if __name__ == "__main__":
    main()
