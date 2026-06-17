#!/usr/bin/env python3
"""
CellWarp — T1-A Replication: HCA Human Download (Pooled Census)

Downloads human single-cell data from CZ CELLxGENE Census, excluding all Tabula
Sapiens datasets. Pools cells from 157+ independent collections to build a maximally
independent human atlas for Procrustes replication.

Biology
-------
By pooling data from many independent labs and technologies (10x 3', 10x 5',
Smart-seq2, etc.), we maximize independence from the primary Tabula Sapiens
result. Cross-study variance inflates centroid noise, biasing AGAINST replication
— making any positive result conservative (DECISION-092).

Filters
-------
  - Organism: Homo sapiens
  - is_primary_data: True
  - disease: "normal"
  - development_stage: contains "adult"
  - EXCLUDE all Tabula Sapiens datasets (by dataset_id exclusion)
  - Download cap: <=5,000 cells per type (except hepatocyte: ALL)

Output
------
  data/replication/hca_t1a.h5ad — pooled human data in 16,959 ortholog gene space

Usage:
    python scripts/13_t1a_hca_download.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import cellxgene_census
import numpy as np
import pandas as pd
import scipy.sparse as sp

# ---------------------------------------------------------------------------
# Project setup
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

DATA_DIR = PROJECT_ROOT / "data" / "replication"
OUTPUT_PATH = DATA_DIR / "hca_t1a.h5ad"
ORTHOLOG_PATH = PROJECT_ROOT / "data" / "phase1" / "orthologs_human_mouse.csv"

# Census settings
CENSUS_VERSION = "2025-11-08"
ORGANISM = "Homo sapiens"
RANDOM_SEED = 42
MAX_CELLS_PER_TYPE = 5_000  # DECISION-092: larger cap for cross-study noise
MIN_CELLS_GATE = 200  # DECISION-090: replication gate

# Cell types in the 23-type MCA × HCA intersection
# These use Census cell_type label strings (Cell Ontology terms)
TARGET_TYPES = [
    "B cell",
    "CD4-positive, alpha-beta T cell",
    "CD8-positive, alpha-beta T cell",
    "T cell",
    "basal cell",
    "endothelial cell",
    "enterocyte of epithelium of large intestine",
    "epithelial cell",
    "fibroblast",
    "granulocyte",
    "hematopoietic precursor cell",
    "hepatocyte",
    "luminal epithelial cell of mammary gland",
    "macrophage",
    "mesenchymal stem cell",
    "monocyte",
    "myeloid dendritic cell",
    "natural killer cell",
    "neutrophil",
    "pancreatic acinar cell",
    "pancreatic ductal cell",
    "smooth muscle cell",
    "stromal cell",
]

# Types that get ALL cells (no cap) — hepatocyte is borderline on MCA side
UNCAPPED_TYPES = {"hepatocyte"}

# Metadata columns to keep from Census
OBS_COLUMNS = [
    "cell_type",
    "tissue",
    "tissue_general",
    "assay",
    "sex",
    "dataset_id",
    "donor_id",
    "cell_type_ontology_term_id",
    "is_primary_data",
    "disease",
    "development_stage",
]

VAR_COLUMNS = ["feature_id", "feature_name"]


def load_ortholog_gene_order() -> list[str]:
    """
    Load the canonical sorted list of 16,959 human Ensembl IDs.

    Uses the primary analysis gene set (16,959 genes present in both Tabula
    datasets) rather than the full 17,187-gene ortholog table.

    Returns:
        Sorted list of human Ensembl gene IDs.
    """
    # Use primary data's gene set as canonical
    primary_human_path = PROJECT_ROOT / "data" / "phase1" / "human_qc.h5ad"
    if primary_human_path.exists():
        primary = ad.read_h5ad(primary_human_path, backed="r")
        gene_order = sorted(primary.var_names.tolist())
        primary.file.close()
        print(f"  Loaded {len(gene_order):,} canonical genes from primary data")
    else:
        # Fallback
        orthologs = pd.read_csv(ORTHOLOG_PATH)
        gene_order = sorted(orthologs["human_ensembl_id"].values)
        print(f"  Loaded {len(gene_order):,} ortholog genes from {ORTHOLOG_PATH}")
    return gene_order


def get_tabula_sapiens_dataset_ids(census) -> set[str]:
    """
    Find all Tabula Sapiens dataset IDs in Census for exclusion.

    Returns:
        Set of dataset_id strings to exclude.
    """
    datasets_df = census["census_info"]["datasets"].read().concat().to_pandas()

    # Find Tabula Sapiens datasets
    ts_mask = datasets_df["collection_name"].str.contains(
        "Tabula Sapiens", case=False, na=False
    )
    ts_ids = set(datasets_df.loc[ts_mask, "dataset_id"].tolist())

    print(f"  Found {len(ts_ids)} Tabula Sapiens datasets to exclude:")
    for _, row in datasets_df.loc[ts_mask].iterrows():
        title = row.get("dataset_title", "untitled")
        print(f"    - {row['dataset_id'][:12]}...: {title}")

    return ts_ids


def download_cell_type(
    census,
    cell_type: str,
    ts_dataset_ids: set[str],
    gene_order: list[str],
    rng: np.random.Generator,
) -> ad.AnnData | None:
    """
    Download cells for one cell type from Census, excluding Tabula Sapiens.

    Applies filters: is_primary_data, disease=normal, adult development_stage.
    Subsamples to MAX_CELLS_PER_TYPE unless type is in UNCAPPED_TYPES.

    Args:
        census: Open Census SOMA connection.
        cell_type: Census cell_type label string.
        ts_dataset_ids: Tabula Sapiens dataset IDs to exclude.
        gene_order: Canonical list of human Ensembl IDs.
        rng: NumPy random generator for subsampling.

    Returns:
        AnnData in ortholog gene space, or None if no cells found.
    """
    print(f"\n  Downloading: {cell_type}")

    # Build filter — note: Census API doesn't support NOT IN directly,
    # so we filter after download
    value_filter = (
        f"cell_type == '{cell_type}' "
        f"and is_primary_data == True "
        f"and disease == 'normal'"
    )

    # First get obs to check counts and filter development_stage + Tabula Sapiens
    print(f"    Querying obs metadata...")
    obs_df = cellxgene_census.get_obs(
        census,
        ORGANISM,
        value_filter=value_filter,
        column_names=OBS_COLUMNS,
    )

    if len(obs_df) == 0:
        print(f"    WARNING: 0 cells found for '{cell_type}' — trying alternatives")
        return None

    print(f"    Total Census cells (pre-filter): {len(obs_df):,}")

    # Filter to adult only
    # Census uses age-based ontology terms like "30-year-old stage"
    # not the word "adult". Filter: keep >= 18-year-old OR "adult" keyword.
    # Exclude explicitly developmental stages.
    ds = obs_df["development_stage"].astype(str).fillna("")
    # Extract numeric age where possible
    age_match = ds.str.extract(r"(\d+)-year-old", expand=False)
    has_age = age_match.notna()
    age_numeric = pd.to_numeric(age_match, errors="coerce")
    is_adult_age = has_age & (age_numeric >= 18)
    is_adult_keyword = ds.str.contains("adult", case=False, na=False)
    # Exclude fetal/embryonic/newborn
    is_developmental = ds.str.contains(
        "fetal|embryo|newborn|infant|neonatal|child|Carnegie",
        case=False, na=False,
    )
    adult_mask = (is_adult_age | is_adult_keyword) & ~is_developmental
    obs_df = obs_df[adult_mask]
    print(f"    After adult filter (>=18 years): {len(obs_df):,}")

    # Exclude Tabula Sapiens
    ts_mask = obs_df["dataset_id"].isin(ts_dataset_ids)
    n_ts = int(ts_mask.sum())
    obs_df = obs_df[~ts_mask]
    print(f"    After Tabula Sapiens exclusion ({n_ts:,} removed): {len(obs_df):,}")

    if len(obs_df) == 0:
        print(f"    WARNING: No cells remain after filtering for '{cell_type}'")
        return None

    if len(obs_df) < MIN_CELLS_GATE:
        print(
            f"    WARNING: Only {len(obs_df)} cells < {MIN_CELLS_GATE} gate "
            f"for '{cell_type}'"
        )

    # Determine if we need to subsample
    max_cells = None if cell_type in UNCAPPED_TYPES else MAX_CELLS_PER_TYPE
    if max_cells and len(obs_df) > max_cells:
        sample_indices = rng.choice(len(obs_df), size=max_cells, replace=False)
        sample_indices.sort()
        obs_df = obs_df.iloc[sample_indices]
        print(f"    Subsampled to {len(obs_df):,} cells (cap={max_cells:,})")

    # Get the soma_joinid values for selected cells — stored in the index
    soma_joinids = obs_df.index.values

    # Download expression data for selected cells
    print(f"    Downloading expression matrix...")
    adata = cellxgene_census.get_anndata(
        census=census,
        organism=ORGANISM,
        obs_coords=soma_joinids,
        var_value_filter=None,
        obs_column_names=OBS_COLUMNS,
        var_column_names=VAR_COLUMNS,
    )

    print(f"    Downloaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # Log unique datasets
    if "dataset_id" in adata.obs.columns:
        n_datasets = adata.obs["dataset_id"].nunique()
        print(f"    From {n_datasets} unique datasets")

    # Filter to ortholog gene space
    gene_order_set = set(gene_order)
    if "feature_id" in adata.var.columns:
        gene_mask = adata.var["feature_id"].isin(gene_order_set)
        adata = adata[:, gene_mask].copy()
        adata.var.index = adata.var["feature_id"].values
    else:
        gene_mask = adata.var_names.isin(gene_order_set)
        adata = adata[:, gene_mask].copy()

    print(f"    After ortholog filter: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # Normalize: counts per 10k + log1p
    X = adata.X
    if sp.issparse(X):
        X = X.toarray()
    X = X.astype(np.float64)

    cell_totals = X.sum(axis=1, keepdims=True)
    cell_totals[cell_totals == 0] = 1
    X_norm = np.log1p(X / cell_totals * 10000).astype(np.float32)

    # Reindex to full 16,959 gene space (fill missing with 0)
    current_genes = list(adata.var_names)
    gene_to_col = {g: i for i, g in enumerate(current_genes)}

    full_matrix = np.zeros((adata.n_obs, len(gene_order)), dtype=np.float32)
    for i, gene in enumerate(gene_order):
        if gene in gene_to_col:
            full_matrix[:, i] = X_norm[:, gene_to_col[gene]]

    # Build output AnnData
    var_df = pd.DataFrame(index=gene_order)
    var_df.index.name = "gene_id"

    obs_out = adata.obs.copy()
    obs_out["our_cell_type_label"] = cell_type
    obs_out["cell_type_label"] = cell_type
    obs_out["census_collection"] = obs_out.get("dataset_id", "")

    adata_out = ad.AnnData(
        X=sp.csr_matrix(full_matrix),
        obs=obs_out,
        var=var_df,
    )

    return adata_out


def main() -> None:
    """Main entry point for HCA Census download."""
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    # Check for existing output
    if OUTPUT_PATH.exists():
        print(f"  Output already exists: {OUTPUT_PATH}")
        adata = ad.read_h5ad(OUTPUT_PATH)
        print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")
        print("  Skipping download. Delete the file to re-download.")
        return

    # Load reference data
    print("=" * 70)
    print("LOADING REFERENCE DATA")
    print("=" * 70)

    gene_order = load_ortholog_gene_order()
    rng = np.random.default_rng(RANDOM_SEED)

    # Connect to Census
    print("\n" + "=" * 70)
    print("CONNECTING TO CELLxGENE CENSUS")
    print("=" * 70)
    print(f"  Version: {CENSUS_VERSION}")

    all_adatas = []
    download_stats = {}

    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        # Get Tabula Sapiens dataset IDs for exclusion
        ts_ids = get_tabula_sapiens_dataset_ids(census)

        # Download each cell type
        print("\n" + "=" * 70)
        print(f"DOWNLOADING {len(TARGET_TYPES)} CELL TYPES")
        print("=" * 70)

        for i, cell_type in enumerate(TARGET_TYPES):
            print(f"\n  [{i+1}/{len(TARGET_TYPES)}] {cell_type}")

            adata_ct = download_cell_type(
                census, cell_type, ts_ids, gene_order, rng
            )

            if adata_ct is not None and adata_ct.n_obs > 0:
                all_adatas.append(adata_ct)
                download_stats[cell_type] = {
                    "n_cells": int(adata_ct.n_obs),
                    "n_datasets": int(
                        adata_ct.obs["dataset_id"].nunique()
                        if "dataset_id" in adata_ct.obs.columns
                        else 0
                    ),
                    "status": "PASS" if adata_ct.n_obs >= MIN_CELLS_GATE else "BORDERLINE",
                }
            else:
                download_stats[cell_type] = {
                    "n_cells": 0,
                    "n_datasets": 0,
                    "status": "ABSENT",
                }
                print(f"    WARNING: No cells obtained for {cell_type}")

    # Census connection closed
    print("\n  Census connection closed.")

    if not all_adatas:
        print("\nSTOP: No cell types yielded any cells.")
        sys.exit(1)

    # Concatenate all cell types
    print("\n" + "=" * 70)
    print("CONCATENATING AND SAVING")
    print("=" * 70)

    combined = ad.concat(all_adatas, join="outer")
    print(f"  Combined: {combined.n_obs:,} cells × {combined.n_vars:,} genes")

    # Verify gene order matches canonical
    assert list(combined.var_names) == gene_order, "Gene order mismatch in combined data"

    # Save output
    print(f"\n  Saving to {OUTPUT_PATH}...")
    tmp_path = OUTPUT_PATH.with_suffix(".h5ad.tmp")
    combined.write_h5ad(tmp_path)
    tmp_path.rename(OUTPUT_PATH)
    print(f"  Saved: {OUTPUT_PATH}")

    # Save download stats
    stats_path = DATA_DIR / "hca_download_stats.json"
    full_stats = {
        "total_cells": int(combined.n_obs),
        "total_genes": int(combined.n_vars),
        "n_types": len(all_adatas),
        "per_type": download_stats,
        "census_version": CENSUS_VERSION,
        "tabula_sapiens_excluded": len(ts_ids),
        "max_cells_per_type": MAX_CELLS_PER_TYPE,
        "random_seed": RANDOM_SEED,
    }
    with open(stats_path, "w") as f:
        json.dump(full_stats, f, indent=2)
    print(f"  Stats saved: {stats_path}")

    # Summary
    print("\n" + "=" * 70)
    print("HCA DOWNLOAD COMPLETE — SUMMARY")
    print("=" * 70)
    print(f"  Total cells: {combined.n_obs:,}")
    print(f"  Total genes: {combined.n_vars:,}")
    print(f"  Cell types: {len(all_adatas)}")
    print(f"  Tabula Sapiens datasets excluded: {len(ts_ids)}")

    print(f"\n  {'Cell Type':<50} {'Count':>8} {'Datasets':>12} {'Status':>10}")
    print("  " + "-" * 82)
    for ct in sorted(download_stats.keys()):
        info = download_stats[ct]
        print(
            f"  {ct:<50} {info['n_cells']:>8,} "
            f"{info['n_datasets']:>12} "
            f"{info['status']:>10}"
        )


if __name__ == "__main__":
    main()
