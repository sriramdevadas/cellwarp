#!/usr/bin/env python3
"""
CellWarp — CellHint Human-Side Replication

Tests whether the cross-species geometric signal holds when the HUMAN atlas
is substituted (CellHint for Tabula Sapiens) while holding the mouse side
constant (Tabula Mouse).

Biology
-------
CellHint (Xu et al., Cell 2023) harmonizes cell-type annotations across 38
datasets covering 12 tissues and 3.7M human cells from CELLxGENE. We use 9
tissues that are independent of Tabula Sapiens data: blood, heart,
hippocampus, intestine, kidney, liver, lung, lymph node, pancreas.

Three tissues are excluded because they include Tabula Sapiens cells:
bone marrow, skeletal muscle, spleen.

This is the INVERSE of existing replications (Sun2023, PanSci) which hold
Tabula Human constant and vary the mouse atlas.

Math
----
Same Procrustes pipeline as primary analysis. CellHint human centroids
replace Tabula Human centroids; Tabula Mouse centroids are unchanged.
Shared cell types between CellHint and Tabula Mouse 35-type ontology
form the point set for alignment.

Pipeline
--------
  Steps 1-2 (default):
    1. Download CellHint cells from CELLxGENE Census (9 TS-independent tissues)
    2. Map cell types, QC, normalize, compute centroids, report inventory

  Steps 3-4 (--run-procrustes, after advisor confirms):
    3. Procrustes + permutation test (10,000) vs Tabula Mouse
    4. Rigidity ranking Spearman rho vs primary 35-type

Hard constraints
----------------
  - No Tabula Sapiens data on human side
  - Identical pipeline: CPM + log1p, PCA >= 95% variance, 10,000 permutations
  - 16,959 ortholog gene space (same as primary)
  - >= 500 cell gate per type

Output
------
  output/validation/cellhint_replication/
"""

from __future__ import annotations

import argparse
import gc
import json
import sys
from datetime import date
from pathlib import Path

import time as _time

import anndata as ad
import cellxgene_census
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))


def _print(*args, **kwargs):
    """Print with flush to avoid output buffering when not connected to tty."""
    kwargs.setdefault("flush", True)
    print(*args, **kwargs)

from scipy.stats import spearmanr

from cellwarp.procrustes import (
    pca_reduce_centroids,
    procrustes_align,
    permutation_test,
    compute_residual_vectors,
    RANDOM_SEED,
)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

OUTPUT_DIR = Path("output/validation/cellhint_replication")
ORTHOLOG_PATH = Path("data/phase1/orthologs_human_mouse.csv")
MOUSE_CENTROIDS_PATH = Path("output/phase2/scaled_35types/centroids_mouse_35.csv")
RESIDUALS_RANKED_PATH = Path("output/phase2/scaled_35types/residuals_ranked.csv")
HUMAN_QC_PATH = Path("data/phase1/human_qc.h5ad")

N_PERMUTATIONS = 10_000
VARIANCE_THRESHOLD = 0.95
MIN_CELLS = 500
MAX_CELLS_PER_TYPE = 5_000  # Per-type subsample cap per tissue (for memory)

# 9 TS-independent CellHint tissue dataset IDs (Census IDs)
# Collection ID (Discover): 854c0855-23ad-4362-8b77-6b1639e7a9fc
# Census collection_name: "Automatic cell-type harmonization and integration
#   across Human Cell Atlas datasets"
# NOTE: Census dataset_ids differ from Discover API dataset ids.
TS_INDEPENDENT_DATASETS = {
    "Heart": "364bd0c7-f7fd-48ed-99c1-ae26872b1042",
    "Hippocampus": "74014ef8-d2d0-4cbc-8ba2-037f30753ffd",
    "Blood": "d86edd6a-4b5d-437a-ad80-1fd976a5e23a",
    "Lung": "493a8b60-d676-44d1-b022-d14c1ad0b36c",
    "Intestine": "cbec7853-d996-4493-bbf8-5d82857a4e51",
    "Liver": "1d89d081-f43d-401a-be76-bd6ef9259f4e",
    "Kidney": "f95d8919-1f2a-405f-8776-bfecc0ab0f3f",
    "Lymph_node": "ec062e17-ed4b-41a5-b13a-4712b6de5543",
    "Pancreas": "31f657dc-1875-4c4b-a5ca-ce63b3ef3a82",
}

# Excluded tissues (contain Tabula Sapiens data):
#   Skeletal_muscle: 15d374d6-0dfd-4d8e-ade7-81f73dc921ee
#   Spleen:          72f4798d-ba94-48ec-be01-9bb8a1f11f56
#   Bone_marrow:     a20e2f2a-e171-450a-8ddc-8935c2977ad6


# ---------------------------------------------------------------------------
# Cell type mapping
# ---------------------------------------------------------------------------


def map_cellhint_to_ontology(cell_type: str) -> str | None:
    """Map CellHint Cell Ontology annotation to CellWarp ontology.

    CellHint datasets use specific Cell Ontology terms (e.g., "centrilobular
    region hepatocyte", "effector memory CD8-positive, alpha-beta T cell").
    We map these to broader ontology categories that match the Tabula Mouse
    35-type centroid labels.

    Uses specificity-ordered keyword matching to prevent substring collisions
    (e.g., "T follicular helper" matched before generic "T cell", "neuron"
    before "interneuron" could collide with T cell patterns).

    Returns None for cell types not in the CellWarp target ontology.
    """
    ct = cell_type.lower()

    # === HEPATOCYTE — key target type ===
    if "hepatocyte" in ct:
        return "hepatocyte"

    # === CHOLANGIOCYTE ===
    if "cholangiocyte" in ct:
        return "cholangiocyte"

    # === NEURON (before T cell to avoid "interneuron" collision) ===
    if any(x in ct for x in ["neuron", "pyramidal", "granule cell"]):
        return "neuron"

    # === CARDIAC MUSCLE ===
    if any(x in ct for x in ["cardiac muscle", "cardiac myocyte", "cardiomyocyte"]):
        return "cardiac muscle cell"

    # === ASTROCYTE ===
    if "astrocyte" in ct:
        return "astrocyte"

    # === OLIGODENDROCYTE (exclude precursor) ===
    if "oligodendrocyte" in ct and "precursor" not in ct:
        return "oligodendrocyte"

    # === IMMUNE: specific subtypes before generic ===

    # Plasma cell / plasmablast
    if "plasma cell" in ct or "plasmablast" in ct:
        return "plasma cell"

    # NK cells
    if "natural killer" in ct:
        return "natural killer cell"

    # Mast cell
    if "mast cell" in ct:
        return "mast cell"

    # Dendritic cells → myeloid dendritic cell
    if "dendritic cell" in ct:
        return "myeloid dendritic cell"

    # Macrophage (all subtypes including Kupffer and microglia)
    if any(x in ct for x in ["macrophage", "kupffer", "microglial", "microglia"]):
        return "macrophage"

    # Neutrophil
    if "neutrophil" in ct:
        return "neutrophil"

    # Monocyte (all subtypes → single "monocyte")
    if "monocyte" in ct:
        return "monocyte"

    # Erythrocyte
    if ct == "erythrocyte":
        return "erythrocyte"

    # Granulocyte / basophil / eosinophil
    if any(x in ct for x in ["granulocyte", "basophil", "eosinophil"]):
        return "granulocyte"

    # === B CELLS (before T cells to prevent substring collision) ===
    if any(x in ct for x in [
        "b cell", "memory b", "naive b", "follicular b",
        "germinal center b", "class switched", "unswitched memory",
    ]):
        return "B cell"

    # === CD4+ T CELLS (specific subtypes) ===
    if any(x in ct for x in [
        "cd4-positive", "cd4+", "t follicular helper",
        "t-helper", "regulatory t cell", "activated cd4",
    ]):
        return "CD4-positive, alpha-beta T cell"

    # === CD8+ T CELLS ===
    if any(x in ct for x in ["cd8-positive", "cd8+"]):
        return "CD8-positive, alpha-beta T cell"

    # === GENERIC T CELL ===
    if ct == "t cell" or any(x in ct for x in [
        "gamma-delta t", "mucosal invariant t",
    ]):
        return "T cell"
    if ct == "lymphocyte":
        return "T cell"

    # === ENDOTHELIAL (all subtypes including lymphatic) ===
    if any(x in ct for x in ["endothelial", "vasa recta"]):
        return "endothelial cell"

    # === PERICYTE ===
    if "pericyte" in ct:
        return "pericyte"

    # === SMOOTH MUSCLE ===
    if any(x in ct for x in ["smooth muscle", "mural cell"]):
        return "smooth muscle cell"
    if ct == "muscle cell":
        return "smooth muscle cell"

    # === FIBROBLAST (including myofibroblast) ===
    if any(x in ct for x in ["fibroblast", "myofibroblast"]):
        return "fibroblast"

    # === EPITHELIAL (many subtypes) ===
    if any(x in ct for x in [
        "epithelial", "enterocyte", "goblet", "club cell",
        "basal cell", "pulmonary alveolar", "podocyte",
        "tuft cell", "collecting duct", "loop of henle",
        "connecting tubule",
    ]):
        return "epithelial cell"

    # Not mapped: adipocyte, mesothelial, stem cell, progenitors,
    # megakaryocyte, erythroid lineage, ependymal, glial, innate lymphoid,
    # enteroendocrine, paneth, ionocyte, hepatic stellate, skeletal muscle,
    # pancreatic endocrine (A/D/PP/beta), NKT, stromal, etc.
    return None


# ---------------------------------------------------------------------------
# STEP 1-2: Download and inventory
# ---------------------------------------------------------------------------


def run_inventory():
    """Download CellHint data, map cell types, compute centroids, report.

    Downloads from CELLxGENE Census per tissue dataset, maps cell type
    annotations to CellWarp ontology, applies QC (>= 200 genes, <= 20% mito),
    normalizes (CPM + log1p), restricts to 16,959 ortholog gene space,
    and accumulates per-type expression sums for centroid computation.

    Subsamples to <= 5,000 cells per type per tissue to manage memory
    on a laptop. Reports full cell counts from metadata.
    """
    print("=" * 70)
    print("CellHint Human-Side Replication — Steps 1-2: Data & Centroids")
    print("=" * 70)

    # --- Load reference data ---
    print("\nLoading reference data...")
    ortho = pd.read_csv(ORTHOLOG_PATH)
    ortholog_human_ids = set(ortho["human_ensembl_id"])

    mouse_centroids = pd.read_csv(MOUSE_CENTROIDS_PATH, index_col=0)
    full_gene_set = list(mouse_centroids.columns)
    n_full = len(full_gene_set)
    full_gene_idx = {g: i for i, g in enumerate(full_gene_set)}
    print(f"  Ortholog pairs: {len(ortholog_human_ids):,}")
    print(f"  Full gene set: {n_full:,} genes")
    print(f"  Mouse 35-type centroids: {mouse_centroids.shape}")

    # --- Load TS donor list for exclusion verification ---
    print("\nLoading Tabula Sapiens donor list for overlap check...")
    ts_donors = set()
    if HUMAN_QC_PATH.exists():
        human_qc = ad.read_h5ad(HUMAN_QC_PATH, backed="r")
        if "donor_id" in human_qc.obs.columns:
            ts_donors = set(human_qc.obs["donor_id"].unique())
            print(f"  TS donors: {len(ts_donors)}")
        human_qc.file.close()
        del human_qc
    else:
        print(f"  WARNING: {HUMAN_QC_PATH} not found, cannot verify donor overlap")

    # --- Open Census ---
    print(f"\n{'=' * 70}")
    print("STEP 1: Download CellHint from Census")
    print("=" * 70)

    print("\nOpening CELLxGENE Census...")
    census = cellxgene_census.open_soma(census_version="2025-11-08")

    # Verify datasets exist
    print("Checking dataset availability...")
    datasets_df = (
        census["census_info"]["datasets"].read().concat().to_pandas()
    )
    census_ids = set(datasets_df["dataset_id"].values)

    found_datasets = {}
    for tissue, ds_id in TS_INDEPENDENT_DATASETS.items():
        if ds_id in census_ids:
            row = datasets_df[datasets_df["dataset_id"] == ds_id].iloc[0]
            title = row.get("dataset_title", "N/A")
            print(f"  OK  {tissue:<15} {ds_id[:16]}... ({title})")
            found_datasets[tissue] = ds_id
        else:
            print(f"  MISS {tissue:<15} {ds_id[:16]}... NOT IN CENSUS")

    if not found_datasets:
        print("\n  *** FATAL: No CellHint datasets found in Census ***")
        census.close()
        return

    # --- Get ortholog gene soma_joinids for efficient download ---
    print("\nLooking up ortholog gene soma_joinids in Census var table...")
    var_df = (
        census["census_data"]["homo_sapiens"]["ms"]["RNA"]
        .var.read(column_names=["soma_joinid", "feature_id"])
        .concat()
        .to_pandas()
    )
    ortho_var = var_df[var_df["feature_id"].isin(set(full_gene_set))]
    ortholog_soma_ids = ortho_var["soma_joinid"].values.tolist()
    ortho_fid_to_soma = dict(
        zip(ortho_var["feature_id"], ortho_var["soma_joinid"])
    )
    print(
        f"  Ortholog genes in Census: {len(ortholog_soma_ids):,}/{n_full:,} "
        f"({len(ortholog_soma_ids)/n_full:.1%})"
    )
    del var_df
    gc.collect()

    # --- Process each tissue ---
    type_expr_sums: dict[str, np.ndarray] = {}
    type_cell_counts: dict[str, int] = {}
    type_full_counts: dict[str, int] = {}  # Pre-subsample counts
    type_tissue_breakdown: dict[str, dict[str, int]] = {}
    type_tissue_full: dict[str, dict[str, int]] = {}
    tissue_stats: dict[str, dict] = {}
    ts_donor_overlaps: dict[str, list[str]] = {}

    for tissue, ds_id in found_datasets.items():
        print(f"\n--- {tissue} ---")

        # 1. Get obs metadata (fast, small download)
        try:
            obs_df = cellxgene_census.get_obs(
                census,
                "Homo sapiens",
                value_filter=f"dataset_id == '{ds_id}'",
                column_names=["soma_joinid", "cell_type", "tissue", "donor_id"],
            )
        except Exception as e:
            print(f"  ERROR getting metadata: {e}")
            continue

        print(f"  Total cells: {len(obs_df):,}")

        # 2. Check donor overlap with Tabula Sapiens
        if ts_donors and "donor_id" in obs_df.columns:
            ch_donors = set(obs_df["donor_id"].dropna().unique())
            overlap = ts_donors & ch_donors
            if overlap:
                print(
                    f"  *** WARNING: {len(overlap)} donors overlap with TS: "
                    f"{sorted(overlap)[:5]} ***"
                )
                ts_donor_overlaps[tissue] = sorted(overlap)
            else:
                print(f"  Donor overlap with TS: 0 (clean)")

        # 3. Map cell types
        obs_df["our_type"] = obs_df["cell_type"].apply(map_cellhint_to_ontology)
        mapped = obs_df[obs_df["our_type"].notna()].copy()
        unmapped = obs_df[obs_df["our_type"].isna()]

        print(f"  Mapped: {len(mapped):,}, Unmapped: {len(unmapped):,}")

        # Report unmapped types (filter out 0-count Categorical levels)
        if len(unmapped) > 0:
            unmapped_counts = unmapped["cell_type"].value_counts()
            for uct, n in unmapped_counts[unmapped_counts > 0].items():
                print(f"    [skip] {uct}: {n:,}")

        if len(mapped) == 0:
            print(f"  No mapped cells, skipping")
            continue

        # Record full counts (pre-subsample)
        full_type_counts = mapped["our_type"].value_counts()
        for ct, n in full_type_counts.items():
            if ct not in type_full_counts:
                type_full_counts[ct] = 0
                type_tissue_full[ct] = {}
            type_full_counts[ct] += n
            type_tissue_full[ct][tissue] = int(n)

        # 4. Download expression (all cells for dataset × ortholog genes)
        # Use obs_value_filter (sequential scan) NOT obs_coords (random access)
        # — obs_coords with specific soma_joinids is extremely slow on Census.
        print(f"  Downloading expression (obs_value_filter + var_coords)...")
        t0 = _time.time()
        try:
            adata = cellxgene_census.get_anndata(
                census,
                "Homo sapiens",
                obs_value_filter=f"dataset_id == '{ds_id}'",
                var_coords=ortholog_soma_ids,
                obs_column_names=["cell_type", "donor_id"],
                var_column_names=["feature_id", "feature_name"],
            )
        except Exception as e:
            print(f"  ERROR downloading expression: {e}")
            continue

        dt = _time.time() - t0
        print(f"  Downloaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes ({dt:.0f}s)")

        # Map cell types and keep only mapped cells
        adata.obs["our_type"] = adata.obs["cell_type"].apply(
            map_cellhint_to_ontology
        )
        adata = adata[adata.obs["our_type"].notna()].copy()
        print(f"  After type filter: {adata.n_obs:,} cells")

        # Subsample per type locally (for memory during QC/normalization)
        rng = np.random.RandomState(RANDOM_SEED)
        keep_idx = []
        for ct in adata.obs["our_type"].unique():
            ct_mask = np.where(adata.obs["our_type"] == ct)[0]
            if len(ct_mask) > MAX_CELLS_PER_TYPE:
                chosen = rng.choice(ct_mask, MAX_CELLS_PER_TYPE, replace=False)
                keep_idx.extend(chosen)
            else:
                keep_idx.extend(ct_mask)
        keep_idx.sort()
        adata = adata[keep_idx].copy()
        print(f"  After subsample (≤{MAX_CELLS_PER_TYPE}/type): {adata.n_obs:,}")

        # 6. QC: >= 200 genes, <= 20% mito
        n_pre_qc = adata.n_obs
        sc.pp.filter_cells(adata, min_genes=200)
        if adata.n_obs > 0:
            adata.var["mt"] = adata.var["feature_name"].str.startswith("MT-")
            sc.pp.calculate_qc_metrics(
                adata, qc_vars=["mt"], inplace=True
            )
            adata = adata[adata.obs["pct_counts_mt"] <= 20].copy()
        n_post_qc = adata.n_obs
        print(f"  QC: {n_pre_qc:,} → {n_post_qc:,} ({n_pre_qc - n_post_qc:,} removed)")

        if adata.n_obs == 0:
            print(f"  No cells passed QC, skipping")
            del adata
            gc.collect()
            continue

        # 7. Normalize: CPM + log1p (identical to primary pipeline)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        # 8. Per-type accumulation into full gene space
        gene_ids = adata.var["feature_id"].values
        col_to_full = []
        for j, gid in enumerate(gene_ids):
            if gid in full_gene_idx:
                col_to_full.append((j, full_gene_idx[gid]))

        n_ortho_present = len(col_to_full)
        tissue_stats[tissue] = {
            "total_cells": int(len(obs_df)),
            "mapped_cells": int(len(mapped)),
            "qc_cells": int(adata.n_obs),
            "ortho_genes": n_ortho_present,
            "zero_fill": round(1 - n_ortho_present / n_full, 4),
        }
        print(
            f"  Ortholog overlap: {n_ortho_present:,}/{n_full:,} "
            f"(zero-fill: {1 - n_ortho_present/n_full:.1%})"
        )

        type_counts = adata.obs["our_type"].value_counts()
        print(f"  Types after QC:")
        for ct, n in type_counts.items():
            ct_mask = adata.obs["our_type"] == ct
            ct_X = adata[ct_mask].X
            if sp.issparse(ct_X):
                ct_sum = np.array(ct_X.sum(axis=0)).flatten()
            else:
                ct_sum = np.array(ct_X.sum(axis=0)).flatten()

            ct_full = np.zeros(n_full, dtype=np.float64)
            for j_adata, j_full in col_to_full:
                ct_full[j_full] = ct_sum[j_adata]

            if ct not in type_expr_sums:
                type_expr_sums[ct] = np.zeros(n_full, dtype=np.float64)
                type_cell_counts[ct] = 0
                type_tissue_breakdown[ct] = {}
            type_expr_sums[ct] += ct_full
            type_cell_counts[ct] += int(n)
            type_tissue_breakdown[ct][tissue] = int(n)

            print(f"    {ct:<50} {n:>8,}")

        del adata
        gc.collect()

    census.close()
    print("\nCensus connection closed.")

    # ==================================================================
    # STEP 2: Compute centroids and report
    # ==================================================================
    print(f"\n{'=' * 70}")
    print("STEP 2: Centroids & Cell Count Audit")
    print("=" * 70)

    if not type_cell_counts:
        print("  *** FATAL: No cell type data accumulated ***")
        return

    # Cell count audit (use FULL counts for reporting, computation counts for centroid)
    print(f"\n  Cell count audit (>= {MIN_CELLS} gate):")
    print(
        f"  {'Cell type':<50} {'Full':>8} {'Comp':>8}  {'Status'}"
    )
    print(f"  {'-' * 80}")
    audit_rows = []
    for ct in sorted(type_full_counts, key=type_full_counts.get, reverse=True):
        n_full_ct = type_full_counts[ct]
        n_comp = type_cell_counts.get(ct, 0)
        status = "PASS" if n_full_ct >= MIN_CELLS else (
            "BORDERLINE" if n_full_ct >= 200 else "FAIL"
        )
        tissues = ", ".join(
            f"{t}({c:,})" for t, c in
            sorted(
                type_tissue_full.get(ct, {}).items(),
                key=lambda x: -x[1],
            )
        )
        audit_rows.append({
            "cell_type": ct,
            "n_cells_full": n_full_ct,
            "n_cells_computation": n_comp,
            "status": status,
            "tissues": tissues,
        })
        flag = " *** HEPATOCYTE ***" if ct == "hepatocyte" else ""
        print(
            f"  {ct:<50} {n_full_ct:>8,} {n_comp:>8,}  {status}{flag}"
        )

    usable_types = [r["cell_type"] for r in audit_rows if r["status"] == "PASS"]
    print(f"\n  Types passing >= {MIN_CELLS} gate: {len(usable_types)}")

    hep_count = type_full_counts.get("hepatocyte", 0)
    hep_comp = type_cell_counts.get("hepatocyte", 0)
    print(
        f"  Hepatocyte confirmed: {'Y' if hep_count > 0 else 'N'} "
        f"({hep_count:,} full, {hep_comp:,} computation)"
    )

    # Compute centroids from accumulated sums
    centroids = {}
    for ct in usable_types:
        if ct in type_expr_sums and type_cell_counts.get(ct, 0) > 0:
            centroids[ct] = type_expr_sums[ct] / type_cell_counts[ct]

    centroid_df = pd.DataFrame(centroids, index=full_gene_set).T
    centroid_df.index.name = "cell_type"
    print(
        f"\n  CellHint centroid matrix: "
        f"{centroid_df.shape[0]} types x {centroid_df.shape[1]:,} genes"
    )

    # Save centroids
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    centroid_df.to_csv(OUTPUT_DIR / "centroids_cellhint.csv")
    print(f"  Saved: {OUTPUT_DIR / 'centroids_cellhint.csv'}")

    # Find shared types with Tabula Mouse 35-type
    shared_types = sorted(
        set(centroid_df.index) & set(mouse_centroids.index)
    )
    cellhint_only = sorted(
        set(centroid_df.index) - set(mouse_centroids.index)
    )

    print(f"\n  Shared with Tabula Mouse 35-type: {len(shared_types)}")
    for ct in shared_types:
        n = type_full_counts.get(ct, 0)
        print(f"    {ct:<50} {n:>8,}")
    if cellhint_only:
        print(f"\n  CellHint-only (no mouse match): {len(cellhint_only)}")
        for ct in cellhint_only:
            n = type_full_counts.get(ct, 0)
            print(f"    {ct:<50} {n:>8,}")

    # TS donor overlap summary
    if ts_donor_overlaps:
        print(f"\n  *** TS DONOR OVERLAP DETECTED ***")
        for tissue, donors in ts_donor_overlaps.items():
            print(f"    {tissue}: {donors}")
    else:
        print(f"\n  TS donor overlap: NONE (all tissues clean)")

    # Save inventory metadata
    meta = {
        "step": "inventory",
        "date": str(date.today()),
        "tissues_found": list(found_datasets.keys()),
        "tissues_requested": list(TS_INDEPENDENT_DATASETS.keys()),
        "total_mapped_cells": sum(type_full_counts.values()),
        "total_computation_cells": sum(type_cell_counts.values()),
        "hepatocyte_count": hep_count,
        "n_types_above_500": len(usable_types),
        "n_shared_with_mouse_35": len(shared_types),
        "shared_types": shared_types,
        "cellhint_only_types": cellhint_only,
        "cell_count_audit": audit_rows,
        "tissue_stats": tissue_stats,
        "ts_donor_overlaps": ts_donor_overlaps,
    }
    with open(OUTPUT_DIR / "cellhint_inventory.json", "w") as f:
        json.dump(meta, f, indent=2)

    # Final summary
    print(f"\n{'=' * 70}")
    print("INVENTORY COMPLETE — POST FOR ADVISOR REVIEW")
    print("=" * 70)
    print(f"  Total cells (mapped):       {sum(type_full_counts.values()):,}")
    print(f"  Hepatocyte confirmed:       {'Y' if hep_count > 0 else 'N'} ({hep_count:,})")
    print(f"  Types >= 500 cells:         {len(usable_types)}")
    print(f"  Shared with Tabula Mouse:   {len(shared_types)}")
    print(f"  TS donor overlap:           {'YES — INVESTIGATE' if ts_donor_overlaps else 'NONE'}")
    print(f"\n  Centroids: {OUTPUT_DIR / 'centroids_cellhint.csv'}")
    print(f"  Inventory: {OUTPUT_DIR / 'cellhint_inventory.json'}")
    print(f"\n  After advisor confirmation:")
    print(f"  python scripts/33_cellhint_replication.py --run-procrustes")


# ---------------------------------------------------------------------------
# STEP 3-4: Procrustes
# ---------------------------------------------------------------------------


def run_procrustes():
    """Run Procrustes on saved CellHint centroids vs Tabula Mouse.

    Loads pre-computed CellHint human centroids and Tabula Mouse 35-type
    centroids, finds shared cell types, runs PCA + Procrustes alignment +
    10,000-permutation test, computes per-type residuals and rigidity
    ranking correlation against the primary 35-type analysis.
    """
    print("=" * 70)
    print("CellHint Replication — Steps 3-4: Procrustes")
    print("=" * 70)

    # Load centroids
    centroid_path = OUTPUT_DIR / "centroids_cellhint.csv"
    if not centroid_path.exists():
        print(f"  *** FATAL: {centroid_path} not found ***")
        print(f"  Run Steps 1-2 first (without --run-procrustes)")
        return

    ch_centroids = pd.read_csv(centroid_path, index_col=0)
    mouse_centroids = pd.read_csv(MOUSE_CENTROIDS_PATH, index_col=0)

    print(
        f"\n  CellHint: {ch_centroids.shape[0]} types x "
        f"{ch_centroids.shape[1]:,} genes"
    )
    print(
        f"  Mouse:    {mouse_centroids.shape[0]} types x "
        f"{mouse_centroids.shape[1]:,} genes"
    )

    # Shared types
    shared_types = sorted(
        set(ch_centroids.index) & set(mouse_centroids.index)
    )
    n_shared = len(shared_types)
    print(f"  Shared types: {n_shared}")
    for ct in shared_types:
        print(f"    {ct}")

    if n_shared < 4:
        print("\n  *** STOP: < 4 shared types. Cannot run Procrustes. ***")
        return

    ch_sub = ch_centroids.loc[shared_types]
    mouse_sub = mouse_centroids.loc[shared_types]
    assert list(ch_sub.columns) == list(mouse_sub.columns), "Gene space mismatch"

    # --- PCA ---
    print(f"\n--- PCA on combined centroids ---")
    human_pca, mouse_pca, pca_model, types_list = pca_reduce_centroids(
        ch_sub, mouse_sub, variance_threshold=VARIANCE_THRESHOLD
    )

    # --- Procrustes ---
    print(f"\n--- Procrustes alignment ---")
    result = procrustes_align(human_pca, mouse_pca)

    # --- Permutation test ---
    print(f"\n--- Permutation test ({N_PERMUTATIONS:,} iterations) ---")
    p_val, null_dist = permutation_test(
        human_pca, mouse_pca, N_PERMUTATIONS, RANDOM_SEED
    )

    obs_null = result.distance / np.median(null_dist)

    # --- Per-type residuals ---
    print(f"\n--- Per-type residuals ---")
    residuals = compute_residual_vectors(result, types_list)
    residual_mags = {
        ct: float(np.linalg.norm(residuals[ct])) for ct in types_list
    }

    sorted_types = sorted(residual_mags, key=residual_mags.get, reverse=True)
    total_ssr = sum(v**2 for v in residual_mags.values())
    print(f"\n  Residual ranking (n={n_shared}):")
    for i, ct in enumerate(sorted_types, 1):
        pct = residual_mags[ct] ** 2 / total_ssr * 100
        print(
            f"    {i:>2}. {ct:<50} {residual_mags[ct]:>8.3f}  ({pct:.1f}% SSR)"
        )

    # --- Rigidity ranking correlation ---
    print(f"\n{'=' * 70}")
    print("Rigidity Ranking vs Primary 35-type")
    print("=" * 70)

    primary_df = pd.read_csv(RESIDUALS_RANKED_PATH)
    primary_dict = dict(
        zip(primary_df["cell_type"], primary_df["residual_magnitude"])
    )

    matched = sorted(set(residual_mags.keys()) & set(primary_dict.keys()))
    n_matched = len(matched)

    if n_matched >= 4:
        ch_mags = [residual_mags[ct] for ct in matched]
        primary_mags = [primary_dict[ct] for ct in matched]
        rho, rho_p = spearmanr(ch_mags, primary_mags)
    else:
        rho, rho_p = float("nan"), float("nan")

    print(f"\n  Spearman rho = {rho:.3f}, p = {rho_p:.4f} (n={n_matched})")

    # Per-type rank comparison
    if n_matched >= 4:
        ch_rank = pd.Series(
            [residual_mags[ct] for ct in matched], index=matched
        ).rank(ascending=False)
        primary_rank = pd.Series(
            [primary_dict[ct] for ct in matched], index=matched
        ).rank(ascending=False)
        print(
            f"\n  {'Cell type':<50} {'CH rank':>8} {'Prim rank':>10} "
            f"{'Delta':>6}"
        )
        print(f"  {'-' * 76}")
        for ct in matched:
            d = abs(int(ch_rank[ct]) - int(primary_rank[ct]))
            print(
                f"  {ct:<50} {int(ch_rank[ct]):>8} "
                f"{int(primary_rank[ct]):>10} {d:>6}"
            )

    # --- Scaling factor check ---
    scaling = result.scaling
    scaling_flag = ""
    if scaling < 0.8 or scaling > 1.8:
        scaling_flag = " *** OUTSIDE 0.8-1.8 RANGE — FLAG FOR ADVISOR ***"

    # --- Save results ---
    results = {
        "diagnostic": "CellHint human-side replication",
        "date": str(date.today()),
        "dataset": {
            "name": "CellHint (Xu et al., Cell 2023)",
            "source": "CELLxGENE Census (TS-independent subset)",
            "n_tissues": 9,
            "ts_independent": True,
        },
        "procrustes": {
            "n_types": n_shared,
            "cell_types": shared_types,
            "p_value": float(p_val),
            "distance": float(result.distance),
            "obs_null_ratio": float(obs_null),
            "scaling": float(scaling),
            "rotation_det": float(np.linalg.det(result.rotation)),
            "null_median": float(np.median(null_dist)),
            "pca_components": int(pca_model.n_components_),
            "per_type_residuals": {
                ct: {"magnitude": residual_mags[ct]} for ct in types_list
            },
        },
        "rigidity_ranking": {
            "rho": float(rho) if not np.isnan(rho) else None,
            "p_value": float(rho_p) if not np.isnan(rho_p) else None,
            "n_matched": n_matched,
            "matched_types": matched,
        },
        "comparison": {
            "primary_obs_null": 0.522,
            "sun2023_obs_null": 0.554,
            "pansci_obs_null": 0.552,
            "cellhint_obs_null": float(obs_null),
        },
        "random_seed": RANDOM_SEED,
    }

    with open(OUTPUT_DIR / "cellhint_replication.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    np.save(OUTPUT_DIR / "null_distribution.npy", null_dist)

    if n_matched >= 4:
        rank_df = pd.DataFrame(
            {
                "cell_type": matched,
                "cellhint_residual": [residual_mags[ct] for ct in matched],
                "primary_residual": [primary_dict[ct] for ct in matched],
            }
        )
        rank_df.to_csv(OUTPUT_DIR / "ranking_comparison.csv", index=False)

    # --- Plots ---
    # 1. Null distribution
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(
        null_dist, bins=50, alpha=0.7, color="darkcyan", edgecolor="white"
    )
    ax.axvline(
        result.distance,
        color="red",
        linewidth=2,
        label=f"Observed (d={result.distance:.2f})",
    )
    ax.set_title(
        f"CellHint ({n_shared} types) vs Tabula Mouse\n"
        f"p={p_val:.4f}, obs/null={obs_null:.3f}"
    )
    ax.set_xlabel("Procrustes distance")
    ax.set_ylabel("Count")
    ax.legend()
    plt.tight_layout()
    fig.savefig(
        OUTPUT_DIR / "null_distribution.png", dpi=150, bbox_inches="tight"
    )
    plt.close()

    # 2. Rigidity scatter
    if n_matched >= 4:
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.scatter(
            [primary_dict[ct] for ct in matched],
            [residual_mags[ct] for ct in matched],
            s=60,
            c="darkcyan",
            edgecolors="teal",
            linewidths=0.5,
            zorder=3,
        )
        for ct in matched:
            short = ct[:22] + "..." if len(ct) > 22 else ct
            ax.annotate(
                short,
                (primary_dict[ct], residual_mags[ct]),
                fontsize=6,
                ha="left",
                va="bottom",
                xytext=(4, 4),
                textcoords="offset points",
            )
        lims = [
            min(ax.get_xlim()[0], ax.get_ylim()[0]),
            max(ax.get_xlim()[1], ax.get_ylim()[1]),
        ]
        ax.plot(lims, lims, "k--", alpha=0.3, zorder=1)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_xlabel("Primary 35-type residual magnitude")
        ax.set_ylabel("CellHint residual magnitude")
        ax.set_title(
            f"Rigidity: CellHint vs Primary\n"
            f"rho={rho:.3f}, p={rho_p:.4f}, n={n_matched}"
        )
        plt.tight_layout()
        fig.savefig(
            OUTPUT_DIR / "rigidity_scatter.png", dpi=150, bbox_inches="tight"
        )
        plt.close()

    print(
        f"\n  Saved: {OUTPUT_DIR / 'cellhint_replication.json'}, "
        f"null_distribution.npy, plots"
    )

    # --- FINAL SUMMARY ---
    print(f"\n{'=' * 70}")
    print("STEP 4: RESULTS")
    print("=" * 70)
    print(f"\n  obs/null = {obs_null:.3f}")
    print(f"  p = {p_val:.6f}")
    print(f"  n cell types = {n_shared}")
    print(f"  scaling factor = {scaling:.3f}{scaling_flag}")
    print(f"  rigidity rho = {rho:.3f} (p={rho_p:.4f}, n={n_matched})")
    print(f"\n  Comparison (obs/null ratios):")
    print(f"    Primary  (TS human  vs TS mouse):     0.522")
    print(f"    Sun2023  (TS human  vs Sun mouse):    0.554")
    print(f"    PanSci   (TS human  vs PanSci mouse): 0.552")
    print(f"    CellHint (CH human  vs TS mouse):     {obs_null:.3f}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="CellHint human-side replication for CellWarp Procrustes"
    )
    parser.add_argument(
        "--run-procrustes",
        action="store_true",
        help="Run Procrustes (Steps 3-4) using saved centroids from Step 2",
    )
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if args.run_procrustes:
        run_procrustes()
    else:
        run_inventory()


if __name__ == "__main__":
    main()
