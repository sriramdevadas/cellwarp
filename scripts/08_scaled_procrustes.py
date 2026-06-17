#!/usr/bin/env python3
"""
CellWarp — Scaled Procrustes Analysis (35 Cell Types)

Scales up the Phase 2 Procrustes pipeline from 6 to 35 cell types — every type
with ≥500 cells in both Tabula Sapiens (human) and Tabula Muris Senis (mouse).

Biology
-------
With 6 cell types we demonstrated that cross-species divergence follows a coherent
geometric transformation (Procrustes p=0.0035). Scaling to 35 types tests whether
this structure holds at higher resolution — more cell types means more landmarks
in the D'Arcy Thompson grid, giving both a more powerful statistical test and finer-
grained residual vectors to identify cell-type-specific evolutionary divergence.

With 35 points instead of 6, the permutation test draws from 35! ≈ 1.03×10^40
unique permutations — vastly more statistical power than 6! = 720.

Pipeline
--------
1. Download data for 35 cell types from CZ CELLxGENE Census (≤2000 cells/type/species)
2. Align to shared 16,959 ortholog gene space (same as Phase 1)
3. Normalize (counts per 10k + log1p)
4. Compute centroids (35 × 16,959 per species)
5. PCA on combined 70 centroids (95% variance threshold)
6. Procrustes alignment (mouse → human)
7. Permutation test (10,000 iterations)
8. Compute residual deformation vectors
9. Map residuals to gene space
10. Generate comparison tables (6-type vs 35-type)

Checkpoints saved at each stage to survive crashes.

Inputs:
    data/phase1/orthologs_human_mouse.csv   — cached ortholog mapping
    output/phase2/cell_type_inventory_passing.csv — cell type inventory

Outputs (all in output/phase2/scaled_35types/):
    human_scaled.h5ad               — normalized human data (35 types, ≤2000 cells each)
    mouse_scaled.h5ad               — normalized mouse data (35 types)
    centroids_human_35.csv          — human centroids (35 × 16,959)
    centroids_mouse_35.csv          — mouse centroids (35 × 16,959)
    pca_centroids_35.npz            — PCA-reduced centroids
    null_distribution_35.npy        — 10,000 permuted Procrustes distances
    procrustes_results_35.json      — all results
    comparison_6v35.csv             — 6-type vs 35-type comparison table
    residuals_ranked.csv            — all 35 types ranked by residual magnitude
    top_genes_diverged.csv          — top 10 genes for 5 most diverged types
    top_genes_conserved.csv         — top 10 genes for 5 most conserved types

Usage:
    python scripts/08_scaled_procrustes.py
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from cellwarp.data_loader import (
    HUMAN_COLLECTION,
    HUMAN_ORGANISM,
    MAX_CELLS_PER_TYPE,
    MOUSE_COLLECTION,
    MOUSE_ORGANISM,
    OBS_COLUMNS,
    RANDOM_SEED,
    VAR_COLUMNS,
    build_obs_value_filter,
    download_species_data,
    fetch_orthologs,
    filter_to_shared_orthologs,
    get_dataset_ids_for_collection,
    save_h5ad_atomic,
    subsample_per_cell_type,
)
from cellwarp.procrustes import (
    N_PERMUTATIONS,
    PCA_VARIANCE_THRESHOLD,
    compute_centroids,
    compute_residual_vectors,
    map_residuals_to_genes,
    pca_reduce_centroids,
    permutation_test,
    procrustes_align,
    save_results,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MIN_CELLS_SCALED = 500  # ≥500 cells in both species
OUTPUT_DIR = Path("./output/phase2/scaled_35types")
DATA_DIR = Path("./data/phase2_scaled")

# Existing Phase 1 data paths
ORTHOLOGS_CACHE = Path("./data/phase1/orthologs_human_mouse.csv")
INVENTORY_PATH = Path("./output/phase2/cell_type_inventory_passing.csv")

# Existing 6-type results for comparison
RESULTS_6TYPE = Path("./output/phase2/procrustes_results.json")


# ---------------------------------------------------------------------------
# Helper: get cell types passing the 500-cell gate
# ---------------------------------------------------------------------------

def get_passing_cell_types(inventory_path: Path, min_cells: int = MIN_CELLS_SCALED) -> list[str]:
    """
    Read the inventory CSV and return cell type names with ≥min_cells in both species.

    Args:
        inventory_path: Path to cell_type_inventory_passing.csv
        min_cells: Minimum cell count in both species (default 500).

    Returns:
        Sorted list of Cell Ontology names passing the gate.
    """
    df = pd.read_csv(inventory_path)
    passing = df[df["min_count"] >= min_cells]
    cell_types = sorted(passing["cell_type"].tolist())
    print(f"  Cell types with ≥{min_cells} cells in both species: {len(cell_types)}")
    for ct in cell_types:
        row = passing[passing["cell_type"] == ct].iloc[0]
        print(f"    {ct:<50} H={int(row['human_count']):>7,}  M={int(row['mouse_count']):>7,}")
    return cell_types


# ---------------------------------------------------------------------------
# Stage 1: Download and prepare data
# ---------------------------------------------------------------------------

def stage1_download(cell_types: list[str]) -> tuple[ad.AnnData, ad.AnnData]:
    """
    Download, subsample, align, and normalize data for all cell types.

    Checkpoints:
        DATA_DIR / human_scaled.h5ad  (after normalization)
        DATA_DIR / mouse_scaled.h5ad
        DATA_DIR / human_raw.h5ad     (before normalization)
        DATA_DIR / mouse_raw.h5ad

    Returns:
        Tuple of (human_norm, mouse_norm) AnnData objects.
    """
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    human_norm_path = DATA_DIR / "human_scaled.h5ad"
    mouse_norm_path = DATA_DIR / "mouse_scaled.h5ad"

    # Check for final checkpoint
    if human_norm_path.exists() and mouse_norm_path.exists():
        print("  Loading normalized checkpoint...")
        human = ad.read_h5ad(human_norm_path)
        mouse = ad.read_h5ad(mouse_norm_path)
        print(f"  Human: {human.n_obs:,} cells × {human.n_vars:,} genes")
        print(f"  Mouse: {mouse.n_obs:,} cells × {mouse.n_vars:,} genes")
        _print_cell_type_counts(human, mouse)
        return human, mouse

    # Check for raw aligned checkpoint (pre-normalization)
    human_raw_path = DATA_DIR / "human_raw_aligned.h5ad"
    mouse_raw_path = DATA_DIR / "mouse_raw_aligned.h5ad"

    if human_raw_path.exists() and mouse_raw_path.exists():
        print("  Loading raw aligned checkpoint...")
        human_raw = ad.read_h5ad(human_raw_path)
        mouse_raw = ad.read_h5ad(mouse_raw_path)
    else:
        # Download from Census
        human_raw, mouse_raw = _download_from_census(cell_types)

        # Align to shared ortholog gene space
        print("\n  Aligning to shared ortholog gene space...")
        orthologs = fetch_orthologs(cache_path=ORTHOLOGS_CACHE)
        human_raw, mouse_raw = filter_to_shared_orthologs(human_raw, mouse_raw, orthologs)

        # Save raw aligned checkpoint
        print("  Saving raw aligned checkpoint...")
        save_h5ad_atomic(human_raw, human_raw_path)
        save_h5ad_atomic(mouse_raw, mouse_raw_path)

    # Normalize: counts per 10k + log1p
    print("\n  Normalizing (counts per 10k + log1p)...")
    sc.pp.normalize_total(human_raw, target_sum=10_000)
    sc.pp.log1p(human_raw)
    sc.pp.normalize_total(mouse_raw, target_sum=10_000)
    sc.pp.log1p(mouse_raw)

    # Save normalized checkpoint
    print("  Saving normalized checkpoint...")
    save_h5ad_atomic(human_raw, human_norm_path)
    save_h5ad_atomic(mouse_raw, mouse_norm_path)

    _print_cell_type_counts(human_raw, mouse_raw)
    return human_raw, mouse_raw


def _download_from_census(cell_types: list[str]) -> tuple[ad.AnnData, ad.AnnData]:
    """Download raw data from Census for all cell types."""
    import cellxgene_census

    print("  Opening Census connection...")
    with cellxgene_census.open_soma(census_version="2025-11-08") as census:
        # Get dataset IDs
        human_ds = get_dataset_ids_for_collection(census, HUMAN_COLLECTION)
        mouse_ds = get_dataset_ids_for_collection(census, MOUSE_COLLECTION)

        # Build filter — use cell type names directly from inventory (they are
        # already Census Cell Ontology names)
        human_filter = build_obs_value_filter(cell_types, human_ds)
        mouse_filter = build_obs_value_filter(cell_types, mouse_ds)

        # Download
        print("\n  Downloading human data (35 types)...")
        human_raw = download_species_data(census, HUMAN_ORGANISM, human_filter)

        print("\n  Downloading mouse data (35 types)...")
        mouse_raw = download_species_data(census, MOUSE_ORGANISM, mouse_filter)

    # Subsample to ≤2000 cells per type
    print("\n  Subsampling...")
    human_raw = subsample_per_cell_type(human_raw, MAX_CELLS_PER_TYPE)
    mouse_raw = subsample_per_cell_type(mouse_raw, MAX_CELLS_PER_TYPE)

    return human_raw, mouse_raw


def _print_cell_type_counts(human: ad.AnnData, mouse: ad.AnnData) -> None:
    """Print per-cell-type counts for both species."""
    h_counts = human.obs["cell_type"].value_counts()
    m_counts = mouse.obs["cell_type"].value_counts()
    all_types = sorted(set(h_counts.index) | set(m_counts.index))

    print(f"\n  {'Cell Type':<50} {'Human':>7} {'Mouse':>7}")
    print(f"  {'-' * 66}")
    for ct in all_types:
        h = h_counts.get(ct, 0)
        m = m_counts.get(ct, 0)
        print(f"  {ct:<50} {h:>7,} {m:>7,}")
    print(f"  {'-' * 66}")
    print(f"  {'TOTAL':<50} {human.n_obs:>7,} {mouse.n_obs:>7,}")


# ---------------------------------------------------------------------------
# Stage 2: Centroids
# ---------------------------------------------------------------------------

def stage2_centroids(
    human: ad.AnnData, mouse: ad.AnnData
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Compute mean expression centroids per cell type.

    Checkpoint: OUTPUT_DIR / centroids_{human,mouse}_35.csv
    """
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    h_path = OUTPUT_DIR / "centroids_human_35.csv"
    m_path = OUTPUT_DIR / "centroids_mouse_35.csv"

    if h_path.exists() and m_path.exists():
        print("  Loading centroid checkpoint...")
        h_cent = pd.read_csv(h_path, index_col=0)
        m_cent = pd.read_csv(m_path, index_col=0)
        print(f"  Human centroids: {h_cent.shape[0]} types × {h_cent.shape[1]:,} genes")
        print(f"  Mouse centroids: {m_cent.shape[0]} types × {m_cent.shape[1]:,} genes")
        return h_cent, m_cent

    print("\n  Computing human centroids...")
    h_cent = compute_centroids(human)
    print("\n  Computing mouse centroids...")
    m_cent = compute_centroids(mouse)

    # Save checkpoint
    h_cent.to_csv(h_path)
    m_cent.to_csv(m_path)
    print(f"  Saved centroids to {OUTPUT_DIR}/")

    return h_cent, m_cent


# ---------------------------------------------------------------------------
# Stage 3: PCA
# ---------------------------------------------------------------------------

def stage3_pca(
    h_cent: pd.DataFrame, m_cent: pd.DataFrame
) -> tuple[np.ndarray, np.ndarray, object, list[str]]:
    """
    PCA on combined centroids (95% variance threshold).

    With 35 types × 2 species = 70 centroids, max components = 69.

    Checkpoint: OUTPUT_DIR / pca_centroids_35.npz
    """
    human_pca, mouse_pca, pca_model, cell_types = pca_reduce_centroids(
        h_cent, m_cent, PCA_VARIANCE_THRESHOLD
    )

    np.savez(
        OUTPUT_DIR / "pca_centroids_35.npz",
        human=human_pca,
        mouse=mouse_pca,
        cell_types=np.array(cell_types),
        explained_variance_ratio=pca_model.explained_variance_ratio_,
    )
    print(f"  Saved PCA centroids: {OUTPUT_DIR / 'pca_centroids_35.npz'}")

    return human_pca, mouse_pca, pca_model, cell_types


# ---------------------------------------------------------------------------
# Stage 4: Procrustes + permutation test
# ---------------------------------------------------------------------------

def stage4_procrustes(
    human_pca: np.ndarray, mouse_pca: np.ndarray
) -> tuple[object, float, np.ndarray]:
    """
    Procrustes alignment + permutation test (10,000 iterations).

    Checkpoint: OUTPUT_DIR / null_distribution_35.npy
    """
    result = procrustes_align(human_pca, mouse_pca)

    print("\n  Running permutation test (10,000 iterations)...")
    t0 = time.time()
    p_value, null_dist = permutation_test(human_pca, mouse_pca, N_PERMUTATIONS)
    print(f"  Permutation test: {time.time() - t0:.1f}s")

    np.save(OUTPUT_DIR / "null_distribution_35.npy", null_dist)
    return result, p_value, null_dist


# ---------------------------------------------------------------------------
# Stage 5: Residuals + gene mapping
# ---------------------------------------------------------------------------

def stage5_residuals(
    result: object,
    cell_types: list[str],
    pca_model: object,
    gene_names: list[str],
) -> tuple[dict, dict]:
    """Compute residuals and map back to gene space (top 20 genes per type)."""
    residuals = compute_residual_vectors(result, cell_types)
    top_genes = map_residuals_to_genes(residuals, pca_model, gene_names, n_top=20)
    return residuals, top_genes


# ---------------------------------------------------------------------------
# Stage 6: Comparison and summary tables
# ---------------------------------------------------------------------------

def stage6_comparison(
    result: object,
    p_value: float,
    null_dist: np.ndarray,
    residuals: dict,
    top_genes: dict,
    cell_types: list[str],
    pca_model: object,
    n_genes: int,
) -> None:
    """
    Generate comparison tables, ranked residuals, and top gene lists.

    Outputs:
        comparison_6v35.csv     — side-by-side 6-type vs 35-type
        residuals_ranked.csv    — all 35 types by residual magnitude
        top_genes_diverged.csv  — top 10 genes × 5 most diverged types
        top_genes_conserved.csv — top 10 genes × 5 most conserved types
    """
    # ---- 6-type vs 35-type comparison ----
    print("\n" + "=" * 70)
    print("COMPARISON: 6-TYPE vs 35-TYPE RESULTS")
    print("=" * 70)

    old = None
    if RESULTS_6TYPE.exists():
        with open(RESULTS_6TYPE) as f:
            old = json.load(f)

    rows = []
    metrics = [
        ("Cell types", "6", str(len(cell_types))),
        ("PCA components",
         str(old["pca"]["n_components"]) if old else "N/A",
         str(pca_model.n_components_)),
        ("Cumulative variance",
         f"{old['pca']['cumulative_variance_explained']*100:.1f}%" if old else "N/A",
         f"{sum(pca_model.explained_variance_ratio_)*100:.1f}%"),
        ("Procrustes distance",
         f"{old['procrustes']['distance']:.4f}" if old else "N/A",
         f"{result.distance:.4f}"),
        ("SSR",
         f"{old['procrustes']['distance_squared']:.4f}" if old else "N/A",
         f"{result.distance_squared:.4f}"),
        ("Scaling factor",
         f"{old['procrustes']['scaling']:.4f}" if old else "N/A",
         f"{result.scaling:.4f}"),
        ("Permutation p-value",
         f"{old['permutation_test']['p_value']:.6f}" if old else "N/A",
         f"{p_value:.6f}"),
        ("Significant (p<0.01)",
         str(old["permutation_test"]["significant_at_001"]) if old else "N/A",
         str(p_value < 0.01)),
    ]

    print(f"\n  {'Metric':<30} {'6-type':>15} {'35-type':>15}")
    print(f"  {'-' * 62}")
    for name, v6, v35 in metrics:
        print(f"  {name:<30} {v6:>15} {v35:>15}")
        rows.append({"metric": name, "6_type": v6, "35_type": v35})

    comp_df = pd.DataFrame(rows)
    comp_df.to_csv(OUTPUT_DIR / "comparison_6v35.csv", index=False)

    # ---- Ranked residuals ----
    print("\n" + "=" * 70)
    print("ALL 35 CELL TYPES — RANKED BY RESIDUAL MAGNITUDE")
    print("=" * 70)

    total_ssr = result.distance_squared
    ranked = []
    for ct in cell_types:
        mag = float(np.linalg.norm(residuals[ct]))
        pct = (mag**2 / total_ssr * 100) if total_ssr > 0 else 0.0
        ranked.append({"cell_type": ct, "residual_magnitude": mag, "pct_of_ssr": pct})

    ranked_df = pd.DataFrame(ranked).sort_values("residual_magnitude", ascending=False)
    ranked_df = ranked_df.reset_index(drop=True)
    ranked_df.index = ranked_df.index + 1
    ranked_df.index.name = "rank"

    print(f"\n  {'Rank':<5} {'Cell Type':<50} {'‖residual‖':>12} {'% SSR':>8}")
    print(f"  {'-' * 77}")
    for i, row in ranked_df.iterrows():
        print(f"  {i:<5} {row['cell_type']:<50} {row['residual_magnitude']:>12.4f} {row['pct_of_ssr']:>7.1f}%")

    ranked_df.to_csv(OUTPUT_DIR / "residuals_ranked.csv")

    # ---- Top 5 most diverged ----
    most_diverged = ranked_df.head(5)["cell_type"].tolist()
    most_conserved = ranked_df.tail(5)["cell_type"].tolist()

    print("\n" + "=" * 70)
    print("TOP 10 GENES — 5 MOST DIVERGED CELL TYPES")
    print("=" * 70)

    diverged_rows = []
    for ct in most_diverged:
        genes_df = top_genes[ct].head(10)
        print(f"\n  {ct} (‖r‖ = {np.linalg.norm(residuals[ct]):.4f}):")
        for _, g in genes_df.iterrows():
            direction = "↑mouse" if g["loading"] > 0 else "↓mouse"
            print(f"    {int(g['rank']):>3}. {g['gene']:<15} {g['loading']:>+.4f} ({direction})")
            diverged_rows.append({
                "cell_type": ct,
                "rank": int(g["rank"]),
                "gene": g["gene"],
                "loading": float(g["loading"]),
                "direction": direction,
            })

    pd.DataFrame(diverged_rows).to_csv(OUTPUT_DIR / "top_genes_diverged.csv", index=False)

    # ---- Top 5 most conserved ----
    print("\n" + "=" * 70)
    print("TOP 10 GENES — 5 MOST CONSERVED CELL TYPES")
    print("=" * 70)

    conserved_rows = []
    for ct in most_conserved:
        genes_df = top_genes[ct].head(10)
        print(f"\n  {ct} (‖r‖ = {np.linalg.norm(residuals[ct]):.4f}):")
        for _, g in genes_df.iterrows():
            direction = "↑mouse" if g["loading"] > 0 else "↓mouse"
            print(f"    {int(g['rank']):>3}. {g['gene']:<15} {g['loading']:>+.4f} ({direction})")
            conserved_rows.append({
                "cell_type": ct,
                "rank": int(g["rank"]),
                "gene": g["gene"],
                "loading": float(g["loading"]),
                "direction": direction,
            })

    pd.DataFrame(conserved_rows).to_csv(OUTPUT_DIR / "top_genes_conserved.csv", index=False)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    t_start = time.time()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("SCALED PROCRUSTES ANALYSIS — 35 CELL TYPES")
    print("=" * 70)

    # ---- Identify cell types ----
    print("\n" + "=" * 70)
    print("STAGE 0: Identify passing cell types (≥500 cells both species)")
    print("=" * 70)
    cell_types = get_passing_cell_types(INVENTORY_PATH, MIN_CELLS_SCALED)

    # ---- Stage 1: Download + normalize ----
    print("\n" + "=" * 70)
    print("STAGE 1: Download, align, and normalize data")
    print("=" * 70)
    human, mouse = stage1_download(cell_types)

    # Verify gene alignment
    assert list(human.var_names) == list(mouse.var_names), "Gene space mismatch!"
    n_genes = human.n_vars
    gene_names = human.var["feature_name"].tolist()
    print(f"  Gene space verified: {n_genes:,} shared ortholog genes")

    # ---- Stage 2: Centroids ----
    print("\n" + "=" * 70)
    print("STAGE 2: Compute centroids (35 types × 16,959 genes)")
    print("=" * 70)
    h_cent, m_cent = stage2_centroids(human, mouse)

    # ---- Stage 3: PCA ----
    print("\n" + "=" * 70)
    print("STAGE 3: PCA dimensionality reduction (95% variance)")
    print("=" * 70)
    human_pca, mouse_pca, pca_model, ct_order = stage3_pca(h_cent, m_cent)

    # ---- Stage 4: Procrustes + permutation ----
    print("\n" + "=" * 70)
    print("STAGE 4: Procrustes alignment + permutation test")
    print("=" * 70)
    proc_result, p_value, null_dist = stage4_procrustes(human_pca, mouse_pca)

    # ---- Stage 5: Residuals ----
    print("\n" + "=" * 70)
    print("STAGE 5: Residual vectors → gene space")
    print("=" * 70)
    residuals, top_genes = stage5_residuals(proc_result, ct_order, pca_model, gene_names)

    # ---- Save full results JSON ----
    pca_info = {
        "n_components": int(pca_model.n_components_),
        "variance_explained_per_component": pca_model.explained_variance_ratio_.tolist(),
        "cumulative_variance_explained": float(sum(pca_model.explained_variance_ratio_)),
        "n_genes_input": n_genes,
    }

    save_results(
        result=proc_result,
        p_value=p_value,
        null_distribution=null_dist,
        residuals=residuals,
        top_genes=top_genes,
        cell_types=ct_order,
        pca_info=pca_info,
        output_path=OUTPUT_DIR / "procrustes_results_35.json",
    )

    # ---- Stage 6: Comparison tables ----
    stage6_comparison(
        proc_result, p_value, null_dist, residuals, top_genes,
        ct_order, pca_model, n_genes,
    )

    # ---- Final summary ----
    t_total = time.time() - t_start
    print("\n" + "=" * 70)
    print("SCALED PROCRUSTES ANALYSIS — COMPLETE")
    print("=" * 70)
    print(f"\n  35 cell types × {n_genes:,} genes × 2 species")
    print(f"  PCA: {pca_model.n_components_} components ({sum(pca_model.explained_variance_ratio_)*100:.1f}% variance)")
    print(f"  Procrustes distance: {proc_result.distance:.4f}  (SSR: {proc_result.distance_squared:.4f})")
    print(f"  Scaling: {proc_result.scaling:.4f}")
    print(f"  Permutation p-value: {p_value:.6f}")
    print(f"  Significant (p<0.01): {'YES' if p_value < 0.01 else 'NO'}")
    print(f"  Total runtime: {t_total:.1f}s ({t_total/60:.1f} min)")
    print(f"  Output: {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
