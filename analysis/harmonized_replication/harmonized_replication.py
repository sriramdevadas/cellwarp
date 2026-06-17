#!/usr/bin/env python3
"""
CellWarp — Harmonized Centroid Replication Analysis (FW-024)

Tests whether the negative Spearman correlation (ρ = −0.39) between primary
and CellHint per-type residual rankings is driven by confounds (ontology
granularity, tissue composition, cell count asymmetry) or reflects genuine
disagreement.

Biology
-------
Per-type Procrustes residuals measure how much each cell type deviates from
the global cross-species geometric transformation. The ranking of residuals
(which types are most/least conserved) should replicate if the signal is
biological rather than dataset-specific. The original comparison found
ρ = −0.39 (p = 0.16), suggesting potential non-replication. However, the
C4 investigation identified tissue composition as a significant confound
(CellHint tissue count vs |rank diff|: ρ = −0.53, p = 0.044).

This script applies progressive harmonization to equalize conditions:
  (a) Ontology matching — exclude types with granularity mismatches
  (b) Tissue matching — restrict primary human cells to tissues overlapping
      with CellHint
  (c) Cell count capping — downsample to equalize counts across atlases

Math
----
For each harmonization level, we recompute:
  1. Human centroids (mean expression per type) from restricted cell sets
  2. Joint PCA on (n_types × 2) centroids at 95% variance threshold
  3. Procrustes alignment (mouse → human) minimizing ‖X − sYR‖²
  4. Per-type residuals r_i = aligned_mouse_i − human_i
  5. Spearman ρ between primary and CellHint residual magnitude rankings

Pipeline
--------
  Step 0: Load data (primary h5ad, CellHint centroids, mouse centroids)
  Step 1: Ontology matching — flag/exclude overlapping categories
  Step 2: Tissue matching — map CellHint tissues to primary tissue_general
  Step 3: Cell count capping — downsample primary to match CellHint
  Step 4: Recompute Procrustes for primary (harmonized) and CellHint
  Step 5: Correlate residual rankings
  Step 6: Sensitivity analysis across seeds and harmonization levels

Outputs (all in analysis/harmonized_replication/):
  harmonization_mapping.csv, tissue_restriction_table.csv,
  harmonized_residuals_primary.csv, harmonized_residuals_cellhint.csv,
  correlation_results.json, sensitivity_analysis.csv,
  FW024_results_summary.md
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import spearmanr
import anndata as ad

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from cellwarp.procrustes import (
    pca_reduce_centroids,
    procrustes_align,
    compute_residual_vectors,
    _procrustes_distance,
    permutation_test,
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
OUT_DIR = Path(__file__).resolve().parent
FIG_DIR = PROJECT_ROOT / "figures" / "supplementary"

HUMAN_H5AD = PROJECT_ROOT / "data" / "phase2_scaled" / "human_scaled.h5ad"
MOUSE_H5AD = PROJECT_ROOT / "data" / "phase2_scaled" / "mouse_scaled.h5ad"
MOUSE_CENTROIDS_PATH = PROJECT_ROOT / "output" / "phase2" / "scaled_35types" / "centroids_mouse_35.csv"
CELLHINT_CENTROIDS_PATH = PROJECT_ROOT / "output" / "validation" / "cellhint_replication" / "centroids_cellhint.csv"
CELLHINT_INVENTORY_PATH = PROJECT_ROOT / "output" / "validation" / "cellhint_replication" / "cellhint_inventory.json"
PRIMARY_RESIDUALS_PATH = PROJECT_ROOT / "output" / "phase2" / "scaled_35types" / "residuals_ranked.csv"
PRIMARY_RESULTS_PATH = PROJECT_ROOT / "output" / "phase2" / "scaled_35types" / "procrustes_results_35.json"

SEEDS = [42, 1, 2, 3, 4]

# ---------------------------------------------------------------------------
# CellHint tissue → Primary tissue_general mapping
# ---------------------------------------------------------------------------
# CellHint uses 9 TS-independent tissues. We map each to the tissue_general
# categories used in the primary Tabula Sapiens dataset.
CELLHINT_TO_PRIMARY_TISSUES = {
    "Heart": ["heart"],
    "Blood": ["blood"],
    "Lung": ["lung"],
    "Intestine": ["small intestine", "large intestine", "colon"],
    "Liver": ["liver"],
    "Kidney": ["kidney"],
    "Lymph_node": ["lymph node"],
    "Pancreas": ["pancreas"],
    # Hippocampus has no direct match in Tabula Sapiens
    "Hippocampus": [],
}

# Reverse map: primary tissue_general → CellHint tissue name
PRIMARY_TO_CELLHINT = {}
for ch_tissue, primary_tissues in CELLHINT_TO_PRIMARY_TISSUES.items():
    for pt in primary_tissues:
        PRIMARY_TO_CELLHINT[pt] = ch_tissue


def _print(*args, **kwargs):
    kwargs.setdefault("flush", True)
    print(*args, **kwargs)


# ---------------------------------------------------------------------------
# Step 0: Load data
# ---------------------------------------------------------------------------

def load_all_data():
    """Load primary h5ad, CellHint centroids, mouse centroids, inventory."""
    _print("Loading data...")

    # CellHint inventory (per-tissue cell counts)
    with open(CELLHINT_INVENTORY_PATH) as f:
        ch_inventory = json.load(f)

    # Build CellHint per-tissue counts from inventory audit
    ch_tissue_counts = {}  # {cell_type: {tissue: count}}
    for entry in ch_inventory["cell_count_audit"]:
        ct = entry["cell_type"]
        if entry["status"] != "PASS":
            continue
        # Parse tissue string like "Liver(71,243)"
        tissues_str = entry["tissues"]
        tissue_dict = {}
        for part in tissues_str.split("), "):
            part = part.strip().rstrip(")")
            if "(" not in part:
                continue
            tissue_name, count_str = part.split("(", 1)
            tissue_dict[tissue_name.strip()] = int(count_str.replace(",", ""))
        ch_tissue_counts[ct] = tissue_dict

    # CellHint computation counts (used for centroids)
    ch_comp_counts = {}
    for entry in ch_inventory["cell_count_audit"]:
        if entry["status"] == "PASS":
            ch_comp_counts[entry["cell_type"]] = entry["n_cells_computation"]

    # Load CellHint centroids
    ch_centroids = pd.read_csv(CELLHINT_CENTROIDS_PATH, index_col=0)
    _print(f"  CellHint centroids: {ch_centroids.shape}")

    # Load mouse centroids
    mouse_centroids = pd.read_csv(MOUSE_CENTROIDS_PATH, index_col=0)
    _print(f"  Mouse centroids: {mouse_centroids.shape}")

    # Load primary h5ad (backed for memory efficiency)
    _print(f"  Loading primary human h5ad...")
    human_adata = ad.read_h5ad(HUMAN_H5AD)
    _print(f"  Human h5ad: {human_adata.shape}")

    # Load primary residuals for unharmonized comparison
    primary_residuals = pd.read_csv(PRIMARY_RESIDUALS_PATH)
    _print(f"  Primary residuals: {len(primary_residuals)} types")

    return {
        "ch_centroids": ch_centroids,
        "ch_tissue_counts": ch_tissue_counts,
        "ch_comp_counts": ch_comp_counts,
        "mouse_centroids": mouse_centroids,
        "human_adata": human_adata,
        "primary_residuals": primary_residuals,
    }


# ---------------------------------------------------------------------------
# Step 1: Ontology matching
# ---------------------------------------------------------------------------

def ontology_matching(data):
    """Identify cell types at matching annotation granularity.

    The 15 shared types between CellHint and primary 35-type analysis are:
    B cell, CD4+ T, CD8+ T, T cell, endothelial, epithelial, fibroblast,
    hepatocyte, macrophage, monocyte, myeloid DC, NK cell, neutrophil,
    plasma cell, smooth muscle cell.

    Potential granularity issues:
    - "T cell" overlaps with "CD4+" and "CD8+" T cells. In CellHint, "T cell"
      catches gamma-delta T, MAIT, and generic lymphocytes. In primary (TS),
      "T cell" catches cells annotated as generic T/lymphocyte. These are
      DIFFERENT populations from CD4+/CD8+, but what gets classified as
      "generic T" vs specific subtype varies between annotation systems.
      → Flag but include by default; test exclusion in sensitivity.
    """
    ch_types = set(data["ch_centroids"].index)
    mouse_types = set(data["mouse_centroids"].index)
    primary_types = set(data["human_adata"].obs["cell_type"].unique())

    # Shared types across all three
    shared_types = sorted(ch_types & mouse_types & primary_types)
    _print(f"\n  Shared types (CellHint ∩ Mouse ∩ Primary): {len(shared_types)}")

    mapping_rows = []
    for ct in sorted(ch_types | mouse_types | primary_types):
        in_ch = ct in ch_types
        in_mouse = ct in mouse_types
        in_primary = ct in primary_types

        # Determine inclusion/exclusion
        if ct in shared_types:
            if ct == "T cell":
                status = "included_flagged"
                reason = ("Overlaps with CD4+/CD8+ T cells; what gets classified "
                         "as generic T vs subtype varies between annotations. "
                         "Tested with and without in sensitivity analysis.")
            else:
                status = "included"
                reason = "Present in all three datasets at matching granularity"
        else:
            status = "excluded"
            parts = []
            if not in_ch:
                parts.append("absent from CellHint")
            if not in_mouse:
                parts.append("absent from mouse 35-type set")
            if not in_primary:
                parts.append("absent from primary human")
            reason = "; ".join(parts)

        mapping_rows.append({
            "cell_type": ct,
            "in_cellhint": in_ch,
            "in_mouse_35": in_mouse,
            "in_primary_human": in_primary,
            "status": status,
            "reason": reason,
        })

    mapping_df = pd.DataFrame(mapping_rows)
    mapping_df.to_csv(OUT_DIR / "harmonization_mapping.csv", index=False)
    _print(f"  Saved harmonization_mapping.csv ({len(mapping_df)} types)")

    included = [r["cell_type"] for r in mapping_rows if r["status"].startswith("included")]
    _print(f"  Included types: {len(included)}")
    for ct in included:
        flag = " [FLAGGED: T cell overlap]" if ct == "T cell" else ""
        _print(f"    {ct}{flag}")

    return included, mapping_df


# ---------------------------------------------------------------------------
# Step 2: Tissue matching
# ---------------------------------------------------------------------------

def tissue_matching(data, included_types):
    """For each cell type, find tissue intersection between CellHint and primary.

    CellHint uses 9 tissues. We map these to primary tissue_general names
    and identify which tissues each cell type appears in for both atlases.
    """
    _print("\n  Tissue matching...")

    human_adata = data["human_adata"]
    ch_tissue_counts = data["ch_tissue_counts"]

    tissue_rows = []
    type_tissue_map = {}  # {cell_type: list of primary tissue_general to keep}

    for ct in included_types:
        # CellHint tissues for this type
        ch_tissues = set(ch_tissue_counts.get(ct, {}).keys())
        ch_tissues_with_cells = {t for t in ch_tissues
                                 if ch_tissue_counts[ct].get(t, 0) > 0}

        # Map CellHint tissues to primary tissue_general
        primary_tissues_from_ch = set()
        for ch_t in ch_tissues_with_cells:
            mapped = CELLHINT_TO_PRIMARY_TISSUES.get(ch_t, [])
            primary_tissues_from_ch.update(mapped)

        # Primary tissues for this type (with actual cells)
        ct_mask = human_adata.obs["cell_type"] == ct
        if ct_mask.sum() == 0:
            tissue_rows.append({
                "cell_type": ct,
                "cellhint_tissues": ", ".join(sorted(ch_tissues_with_cells)),
                "primary_tissues": "",
                "intersection_tissues": "",
                "n_cellhint_tissues": len(ch_tissues_with_cells),
                "n_primary_tissues": 0,
                "n_intersection": 0,
                "status": "excluded_no_primary_cells",
            })
            continue

        primary_tissue_counts = human_adata.obs.loc[ct_mask, "tissue_general"].value_counts()
        primary_tissues_with_cells = set(primary_tissue_counts[primary_tissue_counts > 0].index)

        # Intersection: primary tissues that map to CellHint tissues
        intersection = primary_tissues_with_cells & primary_tissues_from_ch

        # Count cells in intersection tissues
        n_cells_intersection = 0
        if intersection:
            intersection_mask = ct_mask & human_adata.obs["tissue_general"].isin(intersection)
            n_cells_intersection = intersection_mask.sum()

        status = "keep" if n_cells_intersection >= 50 else "excluded_low_overlap"

        tissue_rows.append({
            "cell_type": ct,
            "cellhint_tissues": ", ".join(sorted(ch_tissues_with_cells)),
            "primary_tissues_all": ", ".join(sorted(primary_tissues_with_cells)),
            "primary_tissues_matched": ", ".join(sorted(intersection)),
            "n_cellhint_tissues": len(ch_tissues_with_cells),
            "n_primary_tissues_all": len(primary_tissues_with_cells),
            "n_intersection": len(intersection),
            "n_primary_cells_all": int(ct_mask.sum()),
            "n_primary_cells_matched": n_cells_intersection,
            "n_cellhint_cells": data["ch_comp_counts"].get(ct, 0),
            "status": status,
        })

        if status == "keep":
            type_tissue_map[ct] = sorted(intersection)

    tissue_df = pd.DataFrame(tissue_rows)
    tissue_df.to_csv(OUT_DIR / "tissue_restriction_table.csv", index=False)
    _print(f"  Saved tissue_restriction_table.csv")

    kept_types = sorted(type_tissue_map.keys())
    excluded = [ct for ct in included_types if ct not in type_tissue_map]

    _print(f"\n  Types with tissue overlap: {len(kept_types)}")
    for ct in kept_types:
        tissues = type_tissue_map[ct]
        n_matched = tissue_df.loc[tissue_df["cell_type"] == ct, "n_primary_cells_matched"].values[0]
        _print(f"    {ct:<45} {len(tissues)} tissues, {n_matched} cells")

    if excluded:
        _print(f"  Excluded (no/low tissue overlap): {len(excluded)}")
        for ct in excluded:
            _print(f"    {ct}")

    return kept_types, type_tissue_map, tissue_df


# ---------------------------------------------------------------------------
# Step 3: Recompute centroids with harmonization
# ---------------------------------------------------------------------------

def compute_harmonized_centroids(
    human_adata,
    kept_types,
    type_tissue_map,
    ch_comp_counts,
    seed=42,
    apply_tissue_restriction=True,
    apply_count_capping=True,
):
    """Recompute primary human centroids from tissue-restricted, count-capped cells.

    Args:
        human_adata: Full primary human AnnData
        kept_types: Cell types to include
        type_tissue_map: {cell_type: [tissue_general values]}
        ch_comp_counts: {cell_type: CellHint computation cell count}
        seed: Random seed for downsampling
        apply_tissue_restriction: If True, restrict to intersection tissues
        apply_count_capping: If True, cap at min(n_primary, n_cellhint)
    """
    rng = np.random.RandomState(seed)
    gene_ids = human_adata.var_names.tolist()

    centroids = {}
    cell_counts = {}

    for ct in kept_types:
        # Start with all cells of this type
        ct_mask = human_adata.obs["cell_type"] == ct

        if apply_tissue_restriction and ct in type_tissue_map:
            # Restrict to intersection tissues
            tissue_mask = human_adata.obs["tissue_general"].isin(type_tissue_map[ct])
            ct_mask = ct_mask & tissue_mask

        indices = np.where(ct_mask)[0]
        n_available = len(indices)

        if n_available == 0:
            _print(f"    WARNING: {ct} has 0 cells after tissue restriction, skipping")
            continue

        if apply_count_capping:
            n_cellhint = ch_comp_counts.get(ct, n_available)
            n_cap = min(n_available, n_cellhint)
            if n_cap < n_available:
                indices = rng.choice(indices, n_cap, replace=False)

        # Compute centroid
        X = human_adata.X[indices]
        if sp.issparse(X):
            mean_vec = np.asarray(X.mean(axis=0)).flatten()
        else:
            mean_vec = np.mean(X, axis=0)

        centroids[ct] = mean_vec
        cell_counts[ct] = len(indices)

    centroid_df = pd.DataFrame(centroids, index=gene_ids).T
    centroid_df.index.name = "cell_type"

    return centroid_df, cell_counts


# ---------------------------------------------------------------------------
# Step 4: Run Procrustes pipeline
# ---------------------------------------------------------------------------

def run_procrustes_pipeline(human_centroids, mouse_centroids, cell_types, label=""):
    """Run PCA → Procrustes → residuals pipeline.

    Returns dict with residual magnitudes, Procrustes distance, p-value.
    """
    # Ensure matching types
    shared = sorted(set(human_centroids.index) & set(mouse_centroids.index) & set(cell_types))

    if len(shared) < 4:
        _print(f"  [{label}] Only {len(shared)} shared types, skipping")
        return None

    h_sub = human_centroids.loc[shared]
    m_sub = mouse_centroids.loc[shared]

    # PCA
    human_pca, mouse_pca, pca_model, ct_order = pca_reduce_centroids(
        h_sub, m_sub, variance_threshold=0.95
    )

    # Procrustes alignment
    result = procrustes_align(human_pca, mouse_pca)

    # Permutation test (reduced iterations for speed)
    p_val, null_dist = permutation_test(
        human_pca, mouse_pca, n_permutations=10_000, seed=42
    )

    # Per-type residuals
    residuals = {}
    for i, ct in enumerate(ct_order):
        r = result.aligned_target[i] - result.centered_reference[i]
        residuals[ct] = float(np.linalg.norm(r))

    return {
        "cell_types": ct_order,
        "residuals": residuals,
        "distance": float(result.distance),
        "p_value": float(p_val),
        "n_pca": int(pca_model.n_components_),
        "scaling": float(result.scaling),
    }


# ---------------------------------------------------------------------------
# Step 5: Correlate rankings
# ---------------------------------------------------------------------------

def correlate_rankings(primary_residuals, cellhint_residuals, label=""):
    """Compute Spearman ρ between two sets of per-type residual magnitudes."""
    shared = sorted(set(primary_residuals.keys()) & set(cellhint_residuals.keys()))
    if len(shared) < 4:
        return {"rho": np.nan, "p": np.nan, "n": len(shared)}

    p_vals = [primary_residuals[ct] for ct in shared]
    c_vals = [cellhint_residuals[ct] for ct in shared]

    rho, p = spearmanr(p_vals, c_vals)

    # Ranks for reporting
    p_ranks = pd.Series(p_vals, index=shared).rank(ascending=False)
    c_ranks = pd.Series(c_vals, index=shared).rank(ascending=False)

    per_type = []
    for ct in shared:
        per_type.append({
            "cell_type": ct,
            "primary_residual": primary_residuals[ct],
            "cellhint_residual": cellhint_residuals[ct],
            "primary_rank": int(p_ranks[ct]),
            "cellhint_rank": int(c_ranks[ct]),
            "rank_diff": int(p_ranks[ct] - c_ranks[ct]),
        })

    return {
        "rho": float(rho),
        "p": float(p),
        "n": len(shared),
        "label": label,
        "per_type": per_type,
    }


# ---------------------------------------------------------------------------
# Step 6: Sensitivity analysis
# ---------------------------------------------------------------------------

def sensitivity_analysis(data, kept_types, type_tissue_map, mouse_centroids):
    """Test correlation under progressive harmonization and multiple seeds.

    Levels:
      0. Unharmonized (original 15 shared types, original centroids)
      1. Ontology only (exclude T cell, use original centroids)
      2. Ontology + tissue matching (restricted centroids)
      3. Full (ontology + tissue + count capping)

    For each level, also test 5 random seeds for stability.
    """
    _print("\n" + "=" * 70)
    _print("SENSITIVITY ANALYSIS")
    _print("=" * 70)

    human_adata = data["human_adata"]
    ch_centroids = data["ch_centroids"]
    ch_comp_counts = data["ch_comp_counts"]

    # Types without T cell
    types_no_tcell = [ct for ct in kept_types if ct != "T cell"]

    results_rows = []

    # ── Level 0: Unharmonized ──
    _print("\n--- Level 0: Unharmonized (original) ---")
    # Use pre-computed centroids from primary analysis
    human_centroids_orig = pd.read_csv(
        PROJECT_ROOT / "output" / "phase2" / "scaled_35types" / "centroids_human_35.csv",
        index_col=0,
    )

    all_15 = sorted(set(ch_centroids.index) & set(mouse_centroids.index)
                    & set(human_centroids_orig.index))

    proc_primary_orig = run_procrustes_pipeline(
        human_centroids_orig, mouse_centroids, all_15, "primary-unharmonized"
    )
    proc_ch_orig = run_procrustes_pipeline(
        ch_centroids, mouse_centroids, all_15, "cellhint-unharmonized"
    )

    if proc_primary_orig and proc_ch_orig:
        corr = correlate_rankings(
            proc_primary_orig["residuals"], proc_ch_orig["residuals"], "unharmonized"
        )
        results_rows.append({
            "level": "0_unharmonized",
            "seed": 42,
            "n_types": corr["n"],
            "rho": corr["rho"],
            "p_value": corr["p"],
            "description": f"Original 15 shared types, no restrictions",
        })
        _print(f"  ρ = {corr['rho']:.4f}, p = {corr['p']:.4f}, n = {corr['n']}")

        # Persist Level 0 per-type residuals (matched 15-type Procrustes).
        # Enables Paper 1 Figure S4A and Table S3 matched-baseline column
        # without re-executing upstream (FW-024 inventory correction).
        _level0_primary_rows = [
            {"cell_type": ct, "residual_magnitude": proc_primary_orig["residuals"][ct]}
            for ct in proc_primary_orig["cell_types"]
        ]
        _level0_primary_df = pd.DataFrame(_level0_primary_rows)
        _level0_primary_df["rank"] = (
            _level0_primary_df["residual_magnitude"].rank(ascending=False).astype(int)
        )
        _level0_primary_df = _level0_primary_df.sort_values("rank")
        _level0_primary_df.to_csv(
            OUT_DIR / "harmonized_residuals_primary_level0.csv", index=False
        )

        _level0_cellhint_rows = [
            {"cell_type": ct, "residual_magnitude": proc_ch_orig["residuals"][ct]}
            for ct in proc_ch_orig["cell_types"]
        ]
        _level0_cellhint_df = pd.DataFrame(_level0_cellhint_rows)
        _level0_cellhint_df["rank"] = (
            _level0_cellhint_df["residual_magnitude"].rank(ascending=False).astype(int)
        )
        _level0_cellhint_df = _level0_cellhint_df.sort_values("rank")
        _level0_cellhint_df.to_csv(
            OUT_DIR / "harmonized_residuals_cellhint_level0.csv", index=False
        )
        _print(
            "  Saved harmonized_residuals_primary_level0.csv, "
            "harmonized_residuals_cellhint_level0.csv"
        )

    # ── Level 1: Ontology only (exclude T cell) ──
    _print("\n--- Level 1: Ontology only (exclude T cell) ---")
    all_14 = [ct for ct in all_15 if ct != "T cell"]

    proc_primary_ontology = run_procrustes_pipeline(
        human_centroids_orig, mouse_centroids, all_14, "primary-ontology"
    )
    proc_ch_ontology = run_procrustes_pipeline(
        ch_centroids, mouse_centroids, all_14, "cellhint-ontology"
    )

    if proc_primary_ontology and proc_ch_ontology:
        corr = correlate_rankings(
            proc_primary_ontology["residuals"], proc_ch_ontology["residuals"], "ontology"
        )
        results_rows.append({
            "level": "1_ontology_only",
            "seed": 42,
            "n_types": corr["n"],
            "rho": corr["rho"],
            "p_value": corr["p"],
            "description": "Exclude T cell (overlaps CD4+/CD8+)",
        })
        _print(f"  ρ = {corr['rho']:.4f}, p = {corr['p']:.4f}, n = {corr['n']}")

    # ── Level 2: Ontology + tissue matching ──
    _print("\n--- Level 2: Ontology + tissue matching ---")
    for seed in SEEDS:
        _print(f"  Seed {seed}...")
        h_cent, h_counts = compute_harmonized_centroids(
            human_adata, types_no_tcell, type_tissue_map, ch_comp_counts,
            seed=seed, apply_tissue_restriction=True, apply_count_capping=False,
        )
        types_available = sorted(set(h_cent.index) & set(mouse_centroids.index)
                                  & set(ch_centroids.index))

        proc_p = run_procrustes_pipeline(h_cent, mouse_centroids, types_available,
                                          f"primary-tissue-s{seed}")
        proc_c = run_procrustes_pipeline(ch_centroids, mouse_centroids, types_available,
                                          f"cellhint-tissue-s{seed}")

        if proc_p and proc_c:
            corr = correlate_rankings(proc_p["residuals"], proc_c["residuals"],
                                       f"ontology+tissue (seed={seed})")
            results_rows.append({
                "level": "2_ontology_tissue",
                "seed": seed,
                "n_types": corr["n"],
                "rho": corr["rho"],
                "p_value": corr["p"],
                "description": f"Ontology + tissue restriction, seed={seed}",
            })
            _print(f"    ρ = {corr['rho']:.4f}, p = {corr['p']:.4f}, n = {corr['n']}")

    # ── Level 3: Full harmonization (ontology + tissue + count capping) ──
    _print("\n--- Level 3: Full harmonization ---")
    full_results = {}
    for seed in SEEDS:
        _print(f"  Seed {seed}...")
        h_cent, h_counts = compute_harmonized_centroids(
            human_adata, types_no_tcell, type_tissue_map, ch_comp_counts,
            seed=seed, apply_tissue_restriction=True, apply_count_capping=True,
        )
        types_available = sorted(set(h_cent.index) & set(mouse_centroids.index)
                                  & set(ch_centroids.index))

        proc_p = run_procrustes_pipeline(h_cent, mouse_centroids, types_available,
                                          f"primary-full-s{seed}")
        proc_c = run_procrustes_pipeline(ch_centroids, mouse_centroids, types_available,
                                          f"cellhint-full-s{seed}")

        if proc_p and proc_c:
            corr = correlate_rankings(proc_p["residuals"], proc_c["residuals"],
                                       f"full (seed={seed})")
            results_rows.append({
                "level": "3_full_harmonized",
                "seed": seed,
                "n_types": corr["n"],
                "rho": corr["rho"],
                "p_value": corr["p"],
                "description": f"Ontology + tissue + count cap, seed={seed}",
            })
            full_results[seed] = {
                "corr": corr,
                "proc_primary": proc_p,
                "proc_cellhint": proc_c,
                "human_centroids": h_cent,
                "cell_counts": h_counts,
            }
            _print(f"    ρ = {corr['rho']:.4f}, p = {corr['p']:.4f}, n = {corr['n']}")

    # ── Level 2b: Tissue only (no ontology exclusion) ──
    _print("\n--- Level 2b: Tissue matching only (keep T cell) ---")
    h_cent, h_counts = compute_harmonized_centroids(
        human_adata, kept_types, type_tissue_map, ch_comp_counts,
        seed=42, apply_tissue_restriction=True, apply_count_capping=False,
    )
    types_available = sorted(set(h_cent.index) & set(mouse_centroids.index)
                              & set(ch_centroids.index))
    proc_p = run_procrustes_pipeline(h_cent, mouse_centroids, types_available,
                                      "primary-tissue-withT")
    proc_c = run_procrustes_pipeline(ch_centroids, mouse_centroids, types_available,
                                      "cellhint-tissue-withT")
    if proc_p and proc_c:
        corr = correlate_rankings(proc_p["residuals"], proc_c["residuals"],
                                   "tissue_only_with_T")
        results_rows.append({
            "level": "2b_tissue_only",
            "seed": 42,
            "n_types": corr["n"],
            "rho": corr["rho"],
            "p_value": corr["p"],
            "description": "Tissue restriction only (keep T cell, no count cap)",
        })
        _print(f"  ρ = {corr['rho']:.4f}, p = {corr['p']:.4f}, n = {corr['n']}")

    # ── Level 1b: Count capping only (no tissue restriction) ──
    _print("\n--- Level 1b: Count capping only ---")
    h_cent, h_counts = compute_harmonized_centroids(
        human_adata, types_no_tcell, type_tissue_map, ch_comp_counts,
        seed=42, apply_tissue_restriction=False, apply_count_capping=True,
    )
    types_available = sorted(set(h_cent.index) & set(mouse_centroids.index)
                              & set(ch_centroids.index))
    proc_p = run_procrustes_pipeline(h_cent, mouse_centroids, types_available,
                                      "primary-countcap")
    proc_c = run_procrustes_pipeline(ch_centroids, mouse_centroids, types_available,
                                      "cellhint-countcap")
    if proc_p and proc_c:
        corr = correlate_rankings(proc_p["residuals"], proc_c["residuals"],
                                   "count_cap_only")
        results_rows.append({
            "level": "1b_count_cap_only",
            "seed": 42,
            "n_types": corr["n"],
            "rho": corr["rho"],
            "p_value": corr["p"],
            "description": "Ontology + count capping only (no tissue restriction)",
        })
        _print(f"  ρ = {corr['rho']:.4f}, p = {corr['p']:.4f}, n = {corr['n']}")

    sensitivity_df = pd.DataFrame(results_rows)
    sensitivity_df.to_csv(OUT_DIR / "sensitivity_analysis.csv", index=False)
    _print(f"\n  Saved sensitivity_analysis.csv ({len(sensitivity_df)} rows)")

    return sensitivity_df, full_results


# ---------------------------------------------------------------------------
# Step 7: Generate outputs
# ---------------------------------------------------------------------------

def save_outputs(full_results, sensitivity_df, mapping_df, tissue_df):
    """Save all output files."""
    _print("\n" + "=" * 70)
    _print("SAVING OUTPUTS")
    _print("=" * 70)

    # Use seed=42 as canonical result
    if 42 not in full_results:
        _print("  WARNING: seed=42 not in full_results, using first available")
        seed_key = list(full_results.keys())[0]
    else:
        seed_key = 42

    canonical = full_results[seed_key]
    corr = canonical["corr"]
    proc_p = canonical["proc_primary"]
    proc_c = canonical["proc_cellhint"]

    # Save harmonized residuals
    primary_resid_rows = []
    for ct in proc_p["cell_types"]:
        primary_resid_rows.append({
            "cell_type": ct,
            "residual_magnitude": proc_p["residuals"][ct],
            "rank": 0,  # will fill
        })
    primary_resid_df = pd.DataFrame(primary_resid_rows)
    primary_resid_df["rank"] = primary_resid_df["residual_magnitude"].rank(ascending=False).astype(int)
    primary_resid_df = primary_resid_df.sort_values("rank")
    primary_resid_df.to_csv(OUT_DIR / "harmonized_residuals_primary.csv", index=False)

    cellhint_resid_rows = []
    for ct in proc_c["cell_types"]:
        cellhint_resid_rows.append({
            "cell_type": ct,
            "residual_magnitude": proc_c["residuals"][ct],
            "rank": 0,
        })
    cellhint_resid_df = pd.DataFrame(cellhint_resid_rows)
    cellhint_resid_df["rank"] = cellhint_resid_df["residual_magnitude"].rank(ascending=False).astype(int)
    cellhint_resid_df = cellhint_resid_df.sort_values("rank")
    cellhint_resid_df.to_csv(OUT_DIR / "harmonized_residuals_cellhint.csv", index=False)

    _print(f"  Saved harmonized_residuals_primary.csv")
    _print(f"  Saved harmonized_residuals_cellhint.csv")

    # Correlation results JSON
    seed_results = {}
    for seed, res in full_results.items():
        seed_results[str(seed)] = {
            "rho": res["corr"]["rho"],
            "p_value": res["corr"]["p"],
            "n_types": res["corr"]["n"],
        }

    # Compute bootstrap CI on rho across seeds
    seed_rhos = [r["corr"]["rho"] for r in full_results.values()]
    rho_mean = float(np.mean(seed_rhos))
    rho_std = float(np.std(seed_rhos))
    rho_min = float(np.min(seed_rhos))
    rho_max = float(np.max(seed_rhos))

    # Unharmonized rho for reference
    unharmon_row = sensitivity_df[sensitivity_df["level"] == "0_unharmonized"]
    unharmon_rho = float(unharmon_row["rho"].values[0]) if len(unharmon_row) > 0 else np.nan

    corr_results = {
        "canonical_seed": seed_key,
        "rho": corr["rho"],
        "p_value": corr["p"],
        "n_types": corr["n"],
        "cell_types": corr["per_type"],
        "unharmonized_rho": unharmon_rho,
        "rho_improvement": corr["rho"] - unharmon_rho if not np.isnan(unharmon_rho) else None,
        "per_seed_results": seed_results,
        "seed_stability": {
            "mean_rho": rho_mean,
            "std_rho": rho_std,
            "min_rho": rho_min,
            "max_rho": rho_max,
            "cv": abs(rho_std / rho_mean) if rho_mean != 0 else float("inf"),
        },
        "primary_procrustes": {
            "distance": proc_p["distance"],
            "p_value": proc_p["p_value"],
            "n_pca": proc_p["n_pca"],
            "scaling": proc_p["scaling"],
        },
        "cellhint_procrustes": {
            "distance": proc_c["distance"],
            "p_value": proc_c["p_value"],
            "n_pca": proc_c["n_pca"],
            "scaling": proc_c["scaling"],
        },
    }

    with open(OUT_DIR / "correlation_results.json", "w") as f:
        json.dump(corr_results, f, indent=2)
    _print(f"  Saved correlation_results.json")

    return corr_results, primary_resid_df, cellhint_resid_df


# ---------------------------------------------------------------------------
# Step 8: Publication figure
# ---------------------------------------------------------------------------

def generate_figure(corr_results, sensitivity_df, primary_resid_df, cellhint_resid_df):
    """Generate 3-panel publication-quality figure."""
    _print("\n" + "=" * 70)
    _print("GENERATING FIGURE")
    _print("=" * 70)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))

    # ── Panel A: Harmonized scatter ──
    ax = axes[0]
    per_type = corr_results["cell_types"]
    p_ranks = [d["primary_rank"] for d in per_type]
    c_ranks = [d["cellhint_rank"] for d in per_type]
    labels = [d["cell_type"] for d in per_type]

    ax.scatter(p_ranks, c_ranks, s=90, c="#2563eb", alpha=0.85,
               edgecolors="white", linewidth=0.5, zorder=3)

    # Label each point
    for pr, cr, label in zip(p_ranks, c_ranks, labels):
        # Shorten label
        short = label
        if len(short) > 22:
            short = short[:20] + "…"
        ax.annotate(short, (pr, cr), fontsize=6.5, ha="left", va="bottom",
                    xytext=(3, 3), textcoords="offset points")

    n = corr_results["n_types"]
    ax.plot([1, n], [1, n], "k--", alpha=0.3, linewidth=1)
    ax.set_xlabel("Primary rank (harmonized)", fontsize=10)
    ax.set_ylabel("CellHint rank", fontsize=10)
    rho = corr_results["rho"]
    p_val = corr_results["p_value"]
    ax.set_title(f"A. Harmonized ranking\nρ = {rho:.3f}, p = {p_val:.3f}, n = {n}",
                 fontsize=11, fontweight="bold", loc="left")
    ax.set_xlim(0, n + 1)
    ax.set_ylim(0, n + 1)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.invert_xaxis()

    # ── Panel B: Progressive harmonization bar chart ──
    ax = axes[1]

    # Get one representative per level (seed=42)
    level_order = ["0_unharmonized", "1_ontology_only", "2_ontology_tissue",
                   "3_full_harmonized"]
    level_labels = ["Unharmonized\n(15 types)", "Ontology\n(−T cell)",
                    "Ontology\n+ tissue", "Full\n(+ count cap)"]
    level_colors = ["#94a3b8", "#60a5fa", "#3b82f6", "#1d4ed8"]

    rhos = []
    for lvl in level_order:
        subset = sensitivity_df[(sensitivity_df["level"] == lvl) &
                                (sensitivity_df["seed"] == 42)]
        if len(subset) > 0:
            rhos.append(float(subset["rho"].values[0]))
        else:
            rhos.append(np.nan)

    # Error bars for levels with multiple seeds
    yerr = [0] * len(level_order)
    for i, lvl in enumerate(level_order):
        multi = sensitivity_df[sensitivity_df["level"] == lvl]
        if len(multi) > 1:
            yerr[i] = float(multi["rho"].std())

    bars = ax.bar(range(len(level_order)), rhos, color=level_colors,
                  edgecolor="white", linewidth=0.5, yerr=yerr,
                  capsize=4, error_kw={"linewidth": 1.2})
    ax.set_xticks(range(len(level_order)))
    ax.set_xticklabels(level_labels, fontsize=8)
    ax.set_ylabel("Spearman ρ", fontsize=10)
    ax.axhline(y=0, color="black", linewidth=0.5, linestyle="-")
    ax.set_title("B. Progressive harmonization", fontsize=11,
                 fontweight="bold", loc="left")
    ax.set_ylim(min(min(r for r in rhos if not np.isnan(r)) - 0.15, -0.5),
                max(max(r for r in rhos if not np.isnan(r)) + 0.15, 0.5))

    # Annotate ρ values on bars
    for i, (r, bar) in enumerate(zip(rhos, bars)):
        if not np.isnan(r):
            va = "bottom" if r >= 0 else "top"
            offset = 0.02 if r >= 0 else -0.02
            ax.text(i, r + offset + (yerr[i] if r >= 0 else -yerr[i]),
                    f"{r:.2f}", ha="center", va=va, fontsize=8, fontweight="bold")

    # ── Panel C: Rank change for volatile types ──
    ax = axes[2]

    # Identify the types that drove the original negative correlation
    # Compare unharmonized vs harmonized ranks
    unharmon_level = sensitivity_df[sensitivity_df["level"] == "0_unharmonized"]
    harmon_level = sensitivity_df[(sensitivity_df["level"] == "3_full_harmonized") &
                                  (sensitivity_df["seed"] == 42)]

    # Get per-type rank data from the correlation results
    per_type_data = corr_results["cell_types"]
    if per_type_data:
        cts = [d["cell_type"] for d in per_type_data]
        rank_diffs = [abs(d["rank_diff"]) for d in per_type_data]

        # Sort by absolute rank difference
        sorted_idx = np.argsort(rank_diffs)[::-1]
        cts_sorted = [cts[i] for i in sorted_idx]
        rank_diffs_sorted = [rank_diffs[i] for i in sorted_idx]
        p_ranks_sorted = [per_type_data[i]["primary_rank"] for i in sorted_idx]
        c_ranks_sorted = [per_type_data[i]["cellhint_rank"] for i in sorted_idx]

        # Show all types
        n_show = len(cts_sorted)
        y_pos = range(n_show)

        # Shorten labels
        short_labels = []
        for ct in cts_sorted:
            if "CD4" in ct:
                short_labels.append("CD4+ T")
            elif "CD8" in ct:
                short_labels.append("CD8+ T")
            elif "myeloid dendritic" in ct:
                short_labels.append("myeloid DC")
            elif "natural killer" in ct:
                short_labels.append("NK cell")
            elif "smooth muscle" in ct:
                short_labels.append("smooth muscle")
            else:
                short_labels.append(ct)

        # Horizontal bars showing rank difference
        colors = ["#ef4444" if d > 3 else "#f59e0b" if d > 1 else "#22c55e"
                  for d in rank_diffs_sorted]
        ax.barh(y_pos, rank_diffs_sorted, color=colors, edgecolor="white",
                linewidth=0.5, height=0.7)

        ax.set_yticks(y_pos)
        ax.set_yticklabels(short_labels, fontsize=7.5)
        ax.set_xlabel("|Rank difference|", fontsize=10)
        ax.set_title("C. Per-type rank agreement\n(harmonized)",
                     fontsize=11, fontweight="bold", loc="left")
        ax.invert_yaxis()

        # Add rank annotations
        for i, (pr, cr, rd) in enumerate(zip(p_ranks_sorted, c_ranks_sorted, rank_diffs_sorted)):
            ax.text(rd + 0.1, i, f"P:{pr}→C:{cr}", fontsize=6.5, va="center")

    plt.tight_layout()
    fig_path = FIG_DIR / "harmonized_replication_FW024.png"
    plt.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close()
    _print(f"  Saved figure: {fig_path}")

    return fig_path


# ---------------------------------------------------------------------------
# Step 9: Summary report
# ---------------------------------------------------------------------------

def write_summary(corr_results, sensitivity_df, mapping_df, tissue_df, fig_path):
    """Write human-readable summary."""
    _print("\n  Writing summary report...")

    lines = []
    lines.append("# FW-024: Harmonized Centroid Replication Analysis")
    lines.append("")
    lines.append("## Question")
    lines.append("")
    lines.append("Does the negative Spearman correlation (ρ = −0.39) between primary")
    lines.append("and CellHint per-type residual rankings persist after controlling for")
    lines.append("confounds: ontology granularity, tissue composition, and cell count asymmetry?")
    lines.append("")

    lines.append("## Method")
    lines.append("")
    lines.append("Progressive harmonization of conditions between primary (Tabula Sapiens ×")
    lines.append("Tabula Muris Senis) and CellHint replication analyses:")
    lines.append("")
    lines.append("1. **Ontology matching**: Excluded 'T cell' (overlaps CD4+/CD8+ T cell categories)")
    lines.append("2. **Tissue matching**: Restricted primary human cells to tissues that overlap")
    lines.append("   with CellHint's 9 tissues (Heart, Blood, Lung, Intestine, Liver, Kidney,")
    lines.append("   Lymph node, Pancreas). Hippocampus excluded — no TS counterpart.")
    lines.append("3. **Cell count capping**: Downsampled primary cells to min(n_primary, n_cellhint)")
    lines.append("   per type for balanced centroid estimation.")
    lines.append("4. **Recomputed Procrustes**: PCA → alignment → per-type residuals for both")
    lines.append("   harmonized primary and CellHint analyses.")
    lines.append("")

    lines.append("## Results")
    lines.append("")

    # Sensitivity table
    lines.append("### Progressive Harmonization")
    lines.append("")
    lines.append("| Level | Description | n types | ρ | p-value |")
    lines.append("|-------|-------------|---------|---|---------|")

    for _, row in sensitivity_df[sensitivity_df["seed"] == 42].sort_values("level").iterrows():
        lines.append(f"| {row['level']} | {row['description']} | "
                     f"{int(row['n_types'])} | {row['rho']:.3f} | {row['p_value']:.3f} |")
    lines.append("")

    # Seed stability
    seed_info = corr_results.get("seed_stability", {})
    lines.append("### Seed Stability (Full Harmonization)")
    lines.append("")
    lines.append(f"- Mean ρ across 5 seeds: {seed_info.get('mean_rho', 'N/A'):.4f}")
    lines.append(f"- Std: {seed_info.get('std_rho', 'N/A'):.4f}")
    lines.append(f"- Range: [{seed_info.get('min_rho', 'N/A'):.4f}, {seed_info.get('max_rho', 'N/A'):.4f}]")
    lines.append(f"- CV: {seed_info.get('cv', 'N/A'):.3f}")
    lines.append("")

    # Per-seed detail
    lines.append("| Seed | ρ | p-value |")
    lines.append("|------|---|---------|")
    for seed, info in corr_results.get("per_seed_results", {}).items():
        lines.append(f"| {seed} | {info['rho']:.3f} | {info['p_value']:.3f} |")
    lines.append("")

    # Per-type ranking
    lines.append("### Per-Type Ranking (Harmonized, seed=42)")
    lines.append("")
    lines.append("| Cell Type | Primary Rank | CellHint Rank | Rank Diff |")
    lines.append("|-----------|:------------|:-------------|:----------|")
    per_type = sorted(corr_results["cell_types"], key=lambda d: d["primary_rank"])
    for d in per_type:
        shift = abs(d["rank_diff"])
        flag = " **" if shift >= 4 else ""
        endflag = "**" if shift >= 4 else ""
        lines.append(f"| {d['cell_type']} | {d['primary_rank']} | "
                     f"{d['cellhint_rank']} | {flag}{d['rank_diff']}{endflag} |")
    lines.append("")

    # Interpretation
    unharmon_rho = corr_results.get("unharmonized_rho", np.nan)
    harmon_rho = corr_results["rho"]
    improvement = corr_results.get("rho_improvement", None)

    lines.append("## Interpretation")
    lines.append("")
    lines.append(f"Unharmonized ranking correlation: ρ = {unharmon_rho:.3f}")
    lines.append(f"Fully harmonized correlation: ρ = {harmon_rho:.3f}")
    if improvement is not None:
        direction = "improved" if improvement > 0 else "worsened"
        lines.append(f"Change: Δρ = {improvement:+.3f} ({direction})")
    lines.append("")

    if harmon_rho > 0:
        lines.append("After harmonization, the correlation shifts from negative to positive,")
        lines.append("indicating that the original negative correlation was driven by confounds")
        lines.append("(tissue composition and cell count asymmetry) rather than genuine")
        lines.append("biological disagreement. When conditions are equalized, per-type rigidity")
        lines.append("rankings show directional agreement between primary and CellHint analyses.")
    elif harmon_rho > unharmon_rho:
        lines.append("Harmonization improves the correlation (less negative), suggesting that")
        lines.append("tissue composition and cell count confounds partially explain the original")
        lines.append("negative result. However, substantial residual disagreement remains,")
        lines.append("indicating genuine atlas-specific effects on per-type rankings.")
    else:
        lines.append("Harmonization does not improve the correlation, suggesting that the")
        lines.append("rank disagreement reflects genuine differences in how CellHint and")
        lines.append("Tabula Sapiens capture cross-species variation, not merely tissue")
        lines.append("or count confounds.")

    lines.append("")
    lines.append("## Procrustes Results (Harmonized)")
    lines.append("")
    pp = corr_results.get("primary_procrustes", {})
    cp = corr_results.get("cellhint_procrustes", {})
    lines.append(f"- Primary: distance = {pp.get('distance', 'N/A'):.3f}, "
                 f"p = {pp.get('p_value', 'N/A'):.6f}, "
                 f"PCA components = {pp.get('n_pca', 'N/A')}")
    lines.append(f"- CellHint: distance = {cp.get('distance', 'N/A'):.3f}, "
                 f"p = {cp.get('p_value', 'N/A'):.6f}, "
                 f"PCA components = {cp.get('n_pca', 'N/A')}")
    lines.append("")

    lines.append(f"## Figure")
    lines.append("")
    lines.append(f"![Harmonized replication]({fig_path.relative_to(PROJECT_ROOT)})")
    lines.append("")
    lines.append("---")
    lines.append("*Generated: 2026-04-05 (FW-024)*")

    with open(OUT_DIR / "FW024_results_summary.md", "w") as f:
        f.write("\n".join(lines))
    _print(f"  Saved FW024_results_summary.md")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    _print("=" * 70)
    _print("FW-024: HARMONIZED CENTROID REPLICATION ANALYSIS")
    _print("=" * 70)

    # Step 0: Load data
    _print("\n--- STEP 0: Load data ---")
    data = load_all_data()

    # Step 1: Ontology matching
    _print("\n--- STEP 1: Ontology matching ---")
    included_types, mapping_df = ontology_matching(data)

    # Step 2: Tissue matching
    _print("\n--- STEP 2: Tissue matching ---")
    kept_types, type_tissue_map, tissue_df = tissue_matching(data, included_types)

    # Step 3-6: Sensitivity analysis (includes all harmonization levels)
    sensitivity_df, full_results = sensitivity_analysis(
        data, kept_types, type_tissue_map, data["mouse_centroids"]
    )

    if not full_results:
        _print("\n  *** FATAL: No full harmonization results produced ***")
        return

    # Step 7: Save outputs
    corr_results, primary_resid_df, cellhint_resid_df = save_outputs(
        full_results, sensitivity_df, mapping_df, tissue_df
    )

    # Step 8: Generate figure
    fig_path = generate_figure(
        corr_results, sensitivity_df, primary_resid_df, cellhint_resid_df
    )

    # Step 9: Summary report
    write_summary(corr_results, sensitivity_df, mapping_df, tissue_df, fig_path)

    # Final summary
    _print("\n" + "=" * 70)
    _print("FW-024 COMPLETE")
    _print("=" * 70)
    rho = corr_results["rho"]
    p = corr_results["p_value"]
    unharmon = corr_results.get("unharmonized_rho", float("nan"))
    _print(f"  Unharmonized ρ: {unharmon:.4f}")
    _print(f"  Harmonized ρ:   {rho:.4f} (p = {p:.4f})")
    _print(f"  Δρ:             {rho - unharmon:+.4f}")
    _print(f"  Seed stability: {corr_results['seed_stability']['std_rho']:.4f} (std)")
    _print(f"\n  Outputs: {OUT_DIR}")
    _print(f"  Figure:  {fig_path}")

    # Clean up
    del data["human_adata"]
    import gc; gc.collect()


if __name__ == "__main__":
    main()
