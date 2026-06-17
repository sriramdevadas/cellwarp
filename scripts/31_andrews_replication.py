#!/usr/bin/env python3
"""
CellWarp — Fourth Replication: Andrews et al. 2024 (PSC Liver Atlas, healthy controls)

Downloads human single-cell data from Andrews et al. 2024 (CELLxGENE collection
0c8a364b-97b5-4cc8-a593-23c38c6f0ac5), applies pre-registered label mappings
(DECISION-132), then runs identical Procrustes pipeline against Tabula Muris Senis
mouse centroids.

Biology
-------
Andrews et al. profiled the immunological landscape of healthy and PSC human
liver using 10x Chromium scRNA-seq. We use only the healthy donor data (23 adult
donors, caudate lobe of liver). All 6 original CellWarp cell types are present
but require pre-registered label mappings for hepatocyte zonal subtypes and
liver endothelial subtypes.

Pipeline (identical to three existing replications)
---------------------------------------------------
1. Download from Census: collection-specific, adult, normal, 10x, primary
2. Apply DECISION-132 label mappings
3. Normalize: counts per 10k + log1p
4. Compute centroids per type (mean expression, 16,959 ortholog genes)
5. Load existing TMS mouse centroids
6. Joint PCA (95% variance retained)
7. Procrustes alignment (mouse → human)
8. Permutation test (10,000 iterations)

Output
------
  data/replication/andrews_t1a.h5ad
  output/validation/andrews_replication/

Usage:
    python scripts/31_andrews_replication.py
"""

# obs/null ratio uses null_median (canonical convention, matches src/procrustes.py).

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import anndata as ad
import cellxgene_census
import numpy as np
import pandas as pd
import scipy.sparse as sp

# ---------------------------------------------------------------------------
# Project paths
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

DATA_DIR = PROJECT_ROOT / "data" / "replication"
OUTPUT_DIR = PROJECT_ROOT / "output" / "validation" / "andrews_replication"
H5AD_PATH = DATA_DIR / "andrews_t1a.h5ad"
ORTHOLOG_PATH = PROJECT_ROOT / "data" / "phase1" / "orthologs_human_mouse.csv"
PRIMARY_HUMAN_QC = PROJECT_ROOT / "data" / "phase1" / "human_qc.h5ad"

# Existing TMS mouse centroids from 35-type primary analysis
MOUSE_CENTROIDS_PATH = (
    PROJECT_ROOT / "output" / "phase2" / "scaled_35types" / "centroids_mouse_35.csv"
)

# ---------------------------------------------------------------------------
# Census settings
# ---------------------------------------------------------------------------

CENSUS_VERSION = "2025-11-08"
ORGANISM = "Homo sapiens"
COLLECTION_ID = "0c8a364b-97b5-4cc8-a593-23c38c6f0ac5"
RANDOM_SEED = 42
MAX_CELLS_PER_TYPE = 5_000  # Cap per type for computational tractability

# 10x Chromium assays
CHROMIUM_ASSAYS = {
    "10x 3' v1", "10x 3' v2", "10x 3' v3",
    "10x 3' transcription profiling",
    "10x 5' v1", "10x 5' v2",
    "10x 5' transcription profiling",
}

# ---------------------------------------------------------------------------
# DECISION-132: Pre-registered label mappings
# ---------------------------------------------------------------------------

# Census cell_type labels → CellWarp unified label
LABEL_MAPPINGS = {
    # Hepatocyte: pool zonal subtypes
    "hepatocyte": "hepatocyte",
    "periportal region hepatocyte": "hepatocyte",
    "centrilobular region hepatocyte": "hepatocyte",
    "midzonal region hepatocyte": "hepatocyte",
    # Endothelial: pool liver endothelial subtypes
    "endothelial cell of pericentral hepatic sinusoid": "endothelial cell",
    "endothelial cell of periportal hepatic sinusoid": "endothelial cell",
    "vein endothelial cell": "endothelial cell",
    # Exact matches (no pooling)
    "CD8-positive, alpha-beta T cell": "CD8-positive, alpha-beta T cell",
    "CD4-positive, alpha-beta T cell": "CD4-positive, alpha-beta T cell",
    "B cell": "B cell",
    "macrophage": "macrophage",
}

# Final 6 types expected after mapping
TARGET_TYPES = [
    "hepatocyte",
    "CD8-positive, alpha-beta T cell",
    "endothelial cell",
    "CD4-positive, alpha-beta T cell",
    "B cell",
    "macrophage",
]

MIN_CELLS_GATE = 500

OBS_COLUMNS = [
    "cell_type", "tissue", "tissue_general", "assay", "sex",
    "dataset_id", "donor_id", "cell_type_ontology_term_id",
    "is_primary_data", "disease", "development_stage",
]

VAR_COLUMNS = ["feature_id", "feature_name"]


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def load_ortholog_gene_order() -> list[str]:
    """Load the canonical 16,959 human Ensembl gene IDs from primary data."""
    if PRIMARY_HUMAN_QC.exists():
        primary = ad.read_h5ad(PRIMARY_HUMAN_QC, backed="r")
        gene_order = sorted(primary.var_names.tolist())
        primary.file.close()
        print(f"  Loaded {len(gene_order):,} canonical genes from primary data")
    else:
        orthologs = pd.read_csv(ORTHOLOG_PATH)
        gene_order = sorted(orthologs["human_ensembl_id"].values)
        print(f"  Loaded {len(gene_order):,} ortholog genes from {ORTHOLOG_PATH}")
    return gene_order


def get_collection_dataset_ids(census) -> set[str]:
    """Get all dataset IDs belonging to the Andrews collection."""
    datasets_df = census["census_info"]["datasets"].read().concat().to_pandas()
    coll_mask = datasets_df["collection_id"] == COLLECTION_ID
    ds_ids = set(datasets_df.loc[coll_mask, "dataset_id"].tolist())
    print(f"  Andrews collection: {len(ds_ids)} datasets")
    for _, row in datasets_df.loc[coll_mask].iterrows():
        title = row.get("dataset_title", "untitled")
        print(f"    - {row['dataset_id'][:12]}...: {title[:60]}")
    return ds_ids


def download_andrews_data(census, gene_order: list[str]) -> ad.AnnData:
    """
    Download all target cell types from Andrews collection.

    Applies: is_primary_data, disease=normal, adult, 10x Chromium,
    then DECISION-132 label mappings.
    """
    rng = np.random.default_rng(RANDOM_SEED)
    coll_ds_ids = get_collection_dataset_ids(census)

    # Get all source labels we need to query
    source_labels = list(LABEL_MAPPINGS.keys())

    all_adatas = []
    download_stats = {}

    for source_label in source_labels:
        mapped_type = LABEL_MAPPINGS[source_label]
        print(f"\n  Querying: '{source_label}' → '{mapped_type}'")

        # Query Census for this exact label
        safe_label = source_label.replace("'", "\\'")
        obs_df = cellxgene_census.get_obs(
            census, ORGANISM,
            value_filter=(
                f"cell_type == '{safe_label}' "
                f"and is_primary_data == True "
                f"and disease == 'normal'"
            ),
            column_names=OBS_COLUMNS,
        )
        print(f"    Raw Census: {len(obs_df):,} cells")

        if len(obs_df) == 0:
            print(f"    SKIP: 0 cells")
            continue

        # Filter to Andrews collection only
        obs_df = obs_df[obs_df["dataset_id"].isin(coll_ds_ids)]
        print(f"    Andrews collection: {len(obs_df):,} cells")

        if len(obs_df) == 0:
            print(f"    SKIP: 0 cells in Andrews collection")
            continue

        # Adult filter
        ds = obs_df["development_stage"].astype(str).fillna("")
        age_match = ds.str.extract(r"(\d+)-year-old", expand=False)
        has_age = age_match.notna()
        age_numeric = pd.to_numeric(age_match, errors="coerce")
        is_adult_age = has_age & (age_numeric >= 18)
        is_adult_keyword = ds.str.contains("adult", case=False, na=False)
        is_developmental = ds.str.contains(
            "fetal|embryo|newborn|infant|neonatal|child|Carnegie|juvenile",
            case=False, na=False,
        )
        adult_mask = (is_adult_age | is_adult_keyword) & ~is_developmental
        obs_df = obs_df[adult_mask]
        print(f"    After adult filter: {len(obs_df):,} cells")

        if len(obs_df) == 0:
            continue

        # 10x Chromium filter
        obs_df = obs_df[obs_df["assay"].isin(CHROMIUM_ASSAYS)]
        print(f"    After 10x filter: {len(obs_df):,} cells")

        if len(obs_df) == 0:
            continue

        # Download expression data
        soma_joinids = obs_df.index.values
        print(f"    Downloading expression for {len(soma_joinids):,} cells...")

        adata = cellxgene_census.get_anndata(
            census=census,
            organism=ORGANISM,
            obs_coords=soma_joinids,
            var_value_filter=None,
            obs_column_names=OBS_COLUMNS,
            var_column_names=VAR_COLUMNS,
        )
        print(f"    Downloaded: {adata.n_obs:,} × {adata.n_vars:,}")

        # Filter to ortholog gene space
        gene_order_set = set(gene_order)
        if "feature_id" in adata.var.columns:
            gene_mask = adata.var["feature_id"].isin(gene_order_set)
            adata = adata[:, gene_mask].copy()
            adata.var.index = adata.var["feature_id"].values
        else:
            gene_mask = adata.var_names.isin(gene_order_set)
            adata = adata[:, gene_mask].copy()

        print(f"    Ortholog space: {adata.n_obs:,} × {adata.n_vars:,}")

        # Normalize: counts per 10k + log1p
        X = adata.X
        if sp.issparse(X):
            X = X.toarray()
        X = X.astype(np.float64)
        cell_totals = X.sum(axis=1, keepdims=True)
        cell_totals[cell_totals == 0] = 1
        X_norm = np.log1p(X / cell_totals * 10000).astype(np.float32)

        # Reindex to full 16,959 gene space
        current_genes = list(adata.var_names)
        gene_to_col = {g: i for i, g in enumerate(current_genes)}
        full_matrix = np.zeros((adata.n_obs, len(gene_order)), dtype=np.float32)
        for i, gene in enumerate(gene_order):
            if gene in gene_to_col:
                full_matrix[:, i] = X_norm[:, gene_to_col[gene]]

        # Build output AnnData with mapped label
        var_df = pd.DataFrame(index=gene_order)
        var_df.index.name = "gene_id"
        obs_out = adata.obs.copy()
        obs_out["original_cell_type"] = source_label
        obs_out["cell_type_mapped"] = mapped_type

        adata_out = ad.AnnData(
            X=sp.csr_matrix(full_matrix),
            obs=obs_out,
            var=var_df,
        )
        all_adatas.append(adata_out)

        download_stats[source_label] = {
            "mapped_to": mapped_type,
            "n_cells": int(adata_out.n_obs),
            "n_donors": int(obs_df["donor_id"].nunique()),
        }

    # Concatenate
    if not all_adatas:
        raise RuntimeError("No cells downloaded!")

    combined = ad.concat(all_adatas, join="outer")
    print(f"\n  Combined: {combined.n_obs:,} cells × {combined.n_vars:,} genes")

    # Verify gene order
    assert list(combined.var_names) == gene_order, "Gene order mismatch"

    return combined, download_stats


def compute_centroids(adata: ad.AnnData) -> pd.DataFrame:
    """Compute mean expression per mapped cell type."""
    cell_types = sorted(adata.obs["cell_type_mapped"].unique())
    gene_ids = adata.var_names.tolist()

    centroids = {}
    for ct in cell_types:
        mask = adata.obs["cell_type_mapped"] == ct
        mean_vec = np.asarray(adata[mask].X.mean(axis=0)).flatten()
        n_cells = mask.sum()
        n_donors = adata.obs.loc[mask, "donor_id"].nunique()
        centroids[ct] = mean_vec
        print(f"  {ct:<45} {n_cells:>6,} cells  {n_donors:>3} donors")

    df = pd.DataFrame(centroids, index=gene_ids).T
    df.index.name = "cell_type"
    return df


def main():
    t0 = time.time()
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("CellWarp — Fourth Replication: Andrews et al. 2024")
    print("DECISION-132: Pre-registered label mappings")
    print("=" * 80)

    # -----------------------------------------------------------------------
    # Step 1: Load reference data
    # -----------------------------------------------------------------------
    print("\n[Step 1] Loading reference data...")
    gene_order = load_ortholog_gene_order()

    # Load existing TMS mouse centroids
    print(f"  Loading TMS mouse centroids: {MOUSE_CENTROIDS_PATH}")
    mouse_centroids_full = pd.read_csv(MOUSE_CENTROIDS_PATH, index_col=0)
    print(f"  Mouse centroids: {mouse_centroids_full.shape[0]} types × {mouse_centroids_full.shape[1]} genes")

    # -----------------------------------------------------------------------
    # Step 2: Download Andrews data
    # -----------------------------------------------------------------------
    if H5AD_PATH.exists():
        print(f"\n[Step 2] Loading cached Andrews data: {H5AD_PATH}")
        adata = ad.read_h5ad(H5AD_PATH)
        download_stats = None
        print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    else:
        print(f"\n[Step 2] Downloading Andrews data from Census...")
        with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
            adata, download_stats = download_andrews_data(census, gene_order)

        # Save
        print(f"\n  Saving to {H5AD_PATH}...")
        tmp = H5AD_PATH.with_suffix(".h5ad.tmp")
        adata.write_h5ad(tmp)
        tmp.rename(H5AD_PATH)
        print(f"  Saved.")

        if download_stats:
            stats_path = DATA_DIR / "andrews_download_stats.json"
            with open(stats_path, "w") as f:
                json.dump(download_stats, f, indent=2)

    # -----------------------------------------------------------------------
    # Step 3: Cell count gate check
    # -----------------------------------------------------------------------
    print("\n[Step 3] Cell count gate check...")
    print(f"\n  {'Type':<45} {'Cells':>8} {'Donors':>8} {'Status':>8}")
    print(f"  {'-'*71}")

    gate_pass = True
    for ct in TARGET_TYPES:
        mask = adata.obs["cell_type_mapped"] == ct
        n_cells = mask.sum()
        n_donors = adata.obs.loc[mask, "donor_id"].nunique() if mask.any() else 0
        status = "PASS" if n_cells >= MIN_CELLS_GATE else "FAIL"
        if n_cells < MIN_CELLS_GATE:
            gate_pass = False
        print(f"  {ct:<45} {n_cells:>8,} {n_donors:>8} {status:>8}")

    if not gate_pass:
        print("\n  ABORT: Not all types pass ≥500 gate.")
        sys.exit(1)

    print(f"\n  ALL 6 TYPES PASS ≥{MIN_CELLS_GATE} gate.")

    # -----------------------------------------------------------------------
    # Step 4: Compute Andrews human centroids
    # -----------------------------------------------------------------------
    print("\n[Step 4] Computing Andrews human centroids...")
    andrews_centroids = compute_centroids(adata)

    # Save centroids
    centroids_path = OUTPUT_DIR / "centroids_andrews_human.csv"
    andrews_centroids.to_csv(centroids_path)
    print(f"  Saved: {centroids_path}")

    # -----------------------------------------------------------------------
    # Step 5: Match cell types with TMS mouse centroids
    # -----------------------------------------------------------------------
    print("\n[Step 5] Matching cell types with TMS mouse centroids...")
    andrews_types = set(andrews_centroids.index)
    mouse_types = set(mouse_centroids_full.index)
    shared_types = sorted(andrews_types & mouse_types)

    print(f"  Andrews types: {len(andrews_types)}")
    print(f"  Mouse 35-type set: {len(mouse_types)}")
    print(f"  Shared: {len(shared_types)}")
    for ct in shared_types:
        print(f"    - {ct}")

    if len(shared_types) < 6:
        print(f"\n  WARNING: Only {len(shared_types)} shared types (expected 6).")
        # Check which types are missing
        for ct in TARGET_TYPES:
            if ct not in shared_types:
                print(f"    MISSING: {ct}")
                print(f"      In Andrews: {ct in andrews_types}")
                print(f"      In TMS mouse: {ct in mouse_types}")

    # Subset to shared types
    human_cent = andrews_centroids.loc[shared_types]
    mouse_cent = mouse_centroids_full.loc[shared_types]

    # Ensure gene columns match
    shared_genes = sorted(set(human_cent.columns) & set(mouse_cent.columns))
    print(f"  Shared genes: {len(shared_genes):,}")
    human_cent = human_cent[shared_genes]
    mouse_cent = mouse_cent[shared_genes]

    # -----------------------------------------------------------------------
    # Step 6: PCA + Procrustes + Permutation test
    # -----------------------------------------------------------------------
    print("\n[Step 6] PCA on combined centroids...")
    from cellwarp.procrustes import pca_reduce_centroids, procrustes_align, permutation_test

    human_pca, mouse_pca, pca_model, cell_types = pca_reduce_centroids(
        human_cent, mouse_cent, variance_threshold=0.95
    )

    print("\n[Step 7] Procrustes alignment (TMS mouse → Andrews human)...")
    result = procrustes_align(human_pca, mouse_pca)

    print("\n[Step 8] Permutation test (10,000 iterations)...")
    p_value, null_dist = permutation_test(
        human_pca, mouse_pca, n_permutations=10_000, seed=RANDOM_SEED
    )

    # Compute obs/null ratio
    obs_distance = result.distance
    null_mean = np.mean(null_dist)
    obs_null_ratio = obs_distance / np.median(null_dist)

    # Per-type residuals
    residuals = human_pca - result.aligned_target  # In centered coords: X_c - Y_aligned
    per_type_ssr = np.sum(residuals**2, axis=1)
    total_ssr = np.sum(per_type_ssr)
    per_type_frac = per_type_ssr / total_ssr

    # -----------------------------------------------------------------------
    # Step 7: Report
    # -----------------------------------------------------------------------
    elapsed = time.time() - t0

    print("\n" + "=" * 80)
    print("FOURTH REPLICATION RESULT — Andrews et al. 2024")
    print("=" * 80)

    print(f"\n  CELL COUNTS (after all filters):")
    for ct in shared_types:
        mask = adata.obs["cell_type_mapped"] == ct
        n = mask.sum()
        nd = adata.obs.loc[mask, "donor_id"].nunique()
        print(f"    {ct:<45} {n:>6,} cells  {nd:>3} donors")

    print(f"\n  PROCRUSTES RESULT:")
    print(f"    N types:         {len(shared_types)}")
    print(f"    PCA components:  {pca_model.n_components_}")
    print(f"    Variance:        {np.sum(pca_model.explained_variance_ratio_)*100:.1f}%")
    print(f"    Distance:        {obs_distance:.4f}")
    print(f"    Scaling:         {result.scaling:.4f}")
    print(f"    p-value:         {p_value:.4f}")
    print(f"    Obs/null ratio:  {obs_null_ratio:.4f}")
    print(f"    Null mean:       {null_mean:.4f}")

    print(f"\n  PER-TYPE RESIDUALS:")
    for i, ct in enumerate(cell_types):
        print(f"    {ct:<45} {per_type_frac[i]*100:>6.1f}% SSR")

    print(f"\n  COMPARISON WITH EXISTING REPLICATIONS:")
    print(f"    Tabula (35 types):   obs/null = 0.522,  p = 0.0001")
    print(f"    Sun2023 (15 types):  obs/null = 0.554,  p = 0.0001")
    print(f"    PanSci (16 types):   obs/null = 0.552,  p = 0.0001")
    print(f"    Andrews ({len(shared_types)} types):  obs/null = {obs_null_ratio:.3f},  p = {p_value:.4f}")

    print(f"\n  Runtime: {elapsed:.1f}s")

    # Save results
    results = {
        "dataset": "Andrews et al. 2024",
        "collection_id": COLLECTION_ID,
        "n_types": len(shared_types),
        "types": shared_types,
        "n_pca_components": int(pca_model.n_components_),
        "variance_explained": float(np.sum(pca_model.explained_variance_ratio_)),
        "procrustes_distance": float(obs_distance),
        "scaling": float(result.scaling),
        "p_value": float(p_value),
        "obs_null_ratio": float(obs_null_ratio),
        "null_mean": float(null_mean),
        "null_median": float(np.median(null_dist)),
        "per_type_ssr_fraction": {ct: float(per_type_frac[i]) for i, ct in enumerate(cell_types)},
        "comparison": {
            "tabula_35": {"obs_null": 0.522, "p": 0.0001},
            "sun2023_15": {"obs_null": 0.554, "p": 0.0001},
            "pansci_16": {"obs_null": 0.552, "p": 0.0001},
        },
    }
    results_path = OUTPUT_DIR / "andrews_replication_results.json"
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved: {results_path}")

    # Save null distribution
    null_path = OUTPUT_DIR / "null_distribution.npy"
    np.save(null_path, null_dist)

    # Save PCA
    pca_path = OUTPUT_DIR / "pca_andrews.npz"
    np.savez(pca_path, human=human_pca, mouse=mouse_pca,
             components=pca_model.components_,
             explained_variance_ratio=pca_model.explained_variance_ratio_,
             mean=pca_model.mean_)


if __name__ == "__main__":
    main()
