"""
CellWarp — Data Loader Module

Downloads single-cell RNA-seq data from CZ CELLxGENE Census and aligns human/mouse
gene spaces via Ensembl BioMart 1:1 ortholog mapping.

Biology
-------
We download from two reference atlases:
  - Tabula Sapiens (human): multi-organ single-cell atlas profiling ~500k cells
  - Tabula Muris Senis (mouse): aging mouse atlas profiling ~350k cells

For each atlas we extract five homologous cell types (hepatocytes, CD8+ T cells,
endothelial cells, cardiomyocytes, pancreatic beta cells) chosen for tissue diversity
and well-characterized biology.

Math
----
To compare expression programs across species, both matrices must share the same
column (gene) space. We construct this via 1:1 ortholog mapping from Ensembl BioMart:
genes with exactly one copy in each species (~15-16k pairs). After filtering, both
AnnData objects have identical var dimensions — a prerequisite for Procrustes analysis
in Phase 2. Subsampling to ≤2000 cells/type prevents computational imbalance without
sacrificing statistical power (Phase 1 gate requires ≥500 cells/type/species).
"""

from __future__ import annotations

import time
from pathlib import Path

import anndata as ad
import cellxgene_census
import numpy as np
import pandas as pd
import tiledbsoma as soma

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DATA_DIR = Path("./data/phase1")

# Cell type mapping: our project name → candidate CZ CELLxGENE Cell Ontology terms.
# Some cell types may use different names across Census versions, so we list
# alternatives and discover which one exists at runtime.
CELL_TYPE_MAP: dict[str, list[str]] = {
    "Hepatocytes": ["hepatocyte"],
    "CD8+ T cells": ["CD8-positive, alpha-beta T cell"],
    "Endothelial cells": ["endothelial cell"],
    "CD4+ T cells": ["CD4-positive, alpha-beta T cell"],
    "B cells": ["B cell"],
    "Macrophages": ["macrophage"],
}
# Dropped: Cardiomyocytes (absent from human Tabula Sapiens),
# Pancreatic beta cells (only 102 human cells, below 500 gate).
# See STRATEGY.md 2026-03-12, ISSUES 001-002.

HUMAN_COLLECTION = "Tabula Sapiens"
MOUSE_COLLECTION = "Tabula Muris Senis"

HUMAN_ORGANISM = "Homo sapiens"
MOUSE_ORGANISM = "Mus musculus"

MAX_CELLS_PER_TYPE = 2_000
MIN_CELLS_PER_TYPE = 500   # Phase 1 gate criterion
MIN_SHARED_GENES = 12_000  # Phase 1 gate criterion
RANDOM_SEED = 42

# Metadata columns to keep from Census obs
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
]

# Gene metadata columns to keep from Census var
VAR_COLUMNS = ["feature_id", "feature_name"]

# BioMart retry settings
BIOMART_MAX_RETRIES = 3
BIOMART_RETRY_DELAY_S = 30


# ---------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------


def fetch_orthologs(cache_path: Path | None = None) -> pd.DataFrame:
    """
    Query Ensembl BioMart for human-mouse 1:1 orthologs.

    Biology: Orthologs are genes in different species descended from a common
    ancestor. 1:1 orthologs have exactly one copy in each species, meaning the
    gene was neither duplicated nor lost since human-mouse divergence (~90 Mya).
    These give clean, unambiguous cross-species gene mappings.

    Math: The ortholog table defines the column mapping between human and mouse
    expression matrices. After filtering, both matrices have the same number of
    columns in a 1:1 correspondence — required for Procrustes analysis.

    Args:
        cache_path: If provided and exists, load from CSV cache. If provided and
            file does not exist, query BioMart and save result to this path.

    Returns:
        DataFrame with columns: human_ensembl_id, human_gene_name,
        mouse_ensembl_id, mouse_gene_name, orthology_type.
    """
    if cache_path and cache_path.exists():
        print(f"  Loading cached orthologs from {cache_path}")
        df = pd.read_csv(cache_path)
        print(f"  Loaded {len(df)} 1:1 ortholog pairs from cache")
        return df

    from pybiomart import Dataset

    print("  Querying Ensembl BioMart for human-mouse orthologs...")
    print("  (This typically takes 1-2 minutes)")

    result = None
    for attempt in range(1, BIOMART_MAX_RETRIES + 1):
        try:
            dataset = Dataset(
                name="hsapiens_gene_ensembl",
                host="http://www.ensembl.org",
            )
            result = dataset.query(
                attributes=[
                    "ensembl_gene_id",
                    "external_gene_name",
                    "mmusculus_homolog_ensembl_gene",
                    "mmusculus_homolog_associated_gene_name",
                    "mmusculus_homolog_orthology_type",
                ],
                use_attr_names=True,
            )
            break
        except Exception as e:
            if attempt < BIOMART_MAX_RETRIES:
                print(f"  Attempt {attempt}/{BIOMART_MAX_RETRIES} failed: {e}")
                print(f"  Retrying in {BIOMART_RETRY_DELAY_S}s...")
                time.sleep(BIOMART_RETRY_DELAY_S)
            else:
                raise RuntimeError(
                    f"BioMart query failed after {BIOMART_MAX_RETRIES} attempts. "
                    f"Last error: {e}. Check your internet connection and "
                    f"https://www.ensembl.org status."
                ) from e

    assert result is not None

    # Standardize column names
    result.columns = [
        "human_ensembl_id",
        "human_gene_name",
        "mouse_ensembl_id",
        "mouse_gene_name",
        "orthology_type",
    ]

    # Filter to 1:1 orthologs only
    result = result[result["orthology_type"] == "ortholog_one2one"].copy()

    # Drop rows with missing IDs
    result = result.dropna(subset=["human_ensembl_id", "mouse_ensembl_id"])
    result = result[
        (result["human_ensembl_id"].str.len() > 0)
        & (result["mouse_ensembl_id"].str.len() > 0)
    ]
    result = result.drop_duplicates(subset=["human_ensembl_id", "mouse_ensembl_id"])
    result = result.reset_index(drop=True)

    print(f"  Fetched {len(result)} 1:1 ortholog pairs from Ensembl BioMart")

    if cache_path:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        tmp = cache_path.with_suffix(".csv.tmp")
        result.to_csv(tmp, index=False)
        tmp.rename(cache_path)
        print(f"  Saved to {cache_path}")

    return result


def get_dataset_ids_for_collection(
    census: soma.Collection,
    collection_name: str,
) -> list[str]:
    """
    Look up dataset_ids belonging to a specific collection in the Census.

    Biology: CELLxGENE hosts hundreds of datasets from different labs. We restrict
    to specific atlas collections (Tabula Sapiens / Tabula Muris Senis) to avoid
    mixing data from different experimental protocols, which would introduce
    batch effects.

    Args:
        census: Open Census SOMA collection (from cellxgene_census.open_soma).
        collection_name: Collection to search for (e.g., "Tabula Sapiens").

    Returns:
        List of dataset_id strings belonging to this collection.

    Raises:
        ValueError: If no datasets found for the given collection name.
    """
    datasets_df = census["census_info"]["datasets"].read().concat().to_pandas()

    # Exact match first, fall back to substring match
    mask = datasets_df["collection_name"] == collection_name
    if mask.sum() == 0:
        mask = datasets_df["collection_name"].str.contains(
            collection_name, case=False, na=False
        )

    if mask.sum() == 0:
        available = datasets_df["collection_name"].unique()
        raise ValueError(
            f"Collection '{collection_name}' not found in Census. "
            f"Available collections ({len(available)} total) include: "
            f"{sorted(available)[:10]}"
        )

    matched = datasets_df.loc[mask]
    dataset_ids = matched["dataset_id"].tolist()
    actual_name = matched["collection_name"].iloc[0]

    print(f"  Collection '{actual_name}': {len(dataset_ids)} dataset(s)")
    for _, row in matched.iterrows():
        title = row.get("dataset_title", "untitled")
        print(f"    - {row['dataset_id'][:12]}... : {title}")

    return dataset_ids


def discover_cell_type_names(
    census: soma.Collection,
    organism: str,
    dataset_ids: list[str],
    candidates: dict[str, list[str]],
) -> dict[str, str]:
    """
    Verify which cell type names exist in the Census for a given organism and dataset.

    Biology: Cell Ontology terms can vary between Census versions. This function
    queries actual obs metadata to confirm exact string matches before committing
    to a large expression download.

    Args:
        census: Open Census SOMA collection.
        organism: "Homo sapiens" or "Mus musculus".
        dataset_ids: Dataset IDs to restrict the search to.
        candidates: Dict mapping project cell type name → list of candidate
            Census names to try (e.g., {"Pancreatic beta cells": ["type B
            pancreatic cell", "pancreatic beta cell"]}).

    Returns:
        Dict mapping project name → confirmed Census name (only for types found).

    Raises:
        ValueError: If none of the candidate names are found for any cell type.
    """
    confirmed: dict[str, str] = {}
    dataset_filter = _dataset_id_filter(dataset_ids)

    for project_name, census_names in candidates.items():
        found = False
        for candidate in census_names:
            value_filter = (
                f"cell_type == '{candidate}' and is_primary_data == True"
            )
            if dataset_filter:
                value_filter += f" and {dataset_filter}"

            obs_df = cellxgene_census.get_obs(
                census,
                organism,
                value_filter=value_filter,
                column_names=["cell_type"],
            )
            count = len(obs_df)
            if count > 0:
                print(f"  [FOUND] {project_name} → '{candidate}': {count:,} cells")
                confirmed[project_name] = candidate
                found = True
                break
            else:
                print(f"  [NOT FOUND] {project_name} → '{candidate}'")

        if not found:
            print(f"  [MISSING] {project_name}: none of {census_names} found")

    if len(confirmed) == 0:
        raise ValueError(
            f"No cell types found in {organism} Census. "
            "Check cell type names and dataset IDs."
        )

    return confirmed


def build_obs_value_filter(
    cell_type_names: list[str],
    dataset_ids: list[str] | None = None,
) -> str:
    """
    Construct the obs_value_filter string for the Census API.

    Combines cell type, primary data, and healthy tissue filters into a single
    SOMA query expression.

    Args:
        cell_type_names: List of Cell Ontology names to include.
        dataset_ids: Optional dataset_ids to restrict to.

    Returns:
        Filter string for cellxgene_census.get_anndata(obs_value_filter=...).
    """
    if not cell_type_names:
        raise ValueError("cell_type_names must be non-empty")

    # Build cell type filter
    names_str = ", ".join(f"'{n}'" for n in cell_type_names)
    parts = [
        f"cell_type in [{names_str}]",
        "is_primary_data == True",
        "disease == 'normal'",
    ]

    # Add dataset filter if provided
    ds_filter = _dataset_id_filter(dataset_ids)
    if ds_filter:
        parts.append(ds_filter)

    return " and ".join(parts)


def download_species_data(
    census: soma.Collection,
    organism: str,
    obs_value_filter: str,
) -> ad.AnnData:
    """
    Download expression data from Census for one species.

    Biology: Downloads a cells × genes expression matrix (raw UMI counts) from
    CZ CELLxGENE Census. Each row is a single cell, each column is a gene.
    The matrix is sparse (most genes are not expressed in any given cell).

    Math: The returned AnnData.X contains raw integer UMI counts. These are NOT
    normalized yet — normalization happens in script 02. Raw counts are saved
    as the checkpoint to preserve maximum information.

    Args:
        census: Open Census SOMA collection.
        organism: "Homo sapiens" or "Mus musculus".
        obs_value_filter: Filter string from build_obs_value_filter().

    Returns:
        AnnData with raw counts, cell metadata in .obs, gene metadata in .var.
    """
    print(f"  Downloading {organism} data...")
    print(f"  Filter: {obs_value_filter[:120]}{'...' if len(obs_value_filter) > 120 else ''}")

    adata = cellxgene_census.get_anndata(
        census=census,
        organism=organism,
        obs_value_filter=obs_value_filter,
        obs_column_names=OBS_COLUMNS,
        var_column_names=VAR_COLUMNS,
    )

    print(f"  Downloaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # Print per-cell-type breakdown (skip zero-count categorical entries)
    if "cell_type" in adata.obs.columns:
        counts = adata.obs["cell_type"].value_counts()
        for ct, n in counts.items():
            if n > 0:
                print(f"    {ct}: {n:,} cells")

    return adata


def subsample_per_cell_type(
    adata: ad.AnnData,
    max_cells: int = MAX_CELLS_PER_TYPE,
    cell_type_column: str = "cell_type",
    random_seed: int = RANDOM_SEED,
) -> ad.AnnData:
    """
    Randomly subsample cells so no cell type exceeds max_cells.

    Math: Uniform random sampling without replacement within each cell type.
    Using a fixed seed ensures reproducibility across runs.

    Biology: Cell types with fewer than max_cells are kept intact — we need
    ALL available cells for rare types (e.g., pancreatic beta cells).

    Args:
        adata: Full AnnData from download_species_data().
        max_cells: Maximum cells to retain per cell type.
        cell_type_column: Column in .obs containing cell type labels.
        random_seed: For reproducibility.

    Returns:
        Subsampled AnnData (copy, not a view).
    """
    rng = np.random.default_rng(random_seed)
    indices: list[int] = []

    print(f"  Subsampling to ≤{max_cells:,} cells per type (seed={random_seed})...")

    for ct in sorted(adata.obs[cell_type_column].unique()):
        ct_mask = adata.obs[cell_type_column] == ct
        ct_indices = np.where(ct_mask)[0]
        n_total = len(ct_indices)

        if n_total > max_cells:
            selected = rng.choice(ct_indices, size=max_cells, replace=False)
            selected.sort()
            indices.extend(selected.tolist())
            print(f"    {ct}: {n_total:,} → {max_cells:,} (subsampled)")
        else:
            indices.extend(ct_indices.tolist())
            print(f"    {ct}: {n_total:,} (kept all)")

    return adata[indices].copy()


def filter_to_shared_orthologs(
    human_adata: ad.AnnData,
    mouse_adata: ad.AnnData,
    orthologs: pd.DataFrame,
) -> tuple[ad.AnnData, ad.AnnData]:
    """
    Filter both species' AnnData to shared 1:1 ortholog genes only.

    Biology: Not all 1:1 orthologs are present in both datasets — some genes
    may have been filtered during atlas construction or may not be expressed.
    We intersect three sets: BioMart orthologs, human Census genes, and mouse
    Census genes.

    Math: After this step, human_adata.shape[1] == mouse_adata.shape[1], and
    column i in both matrices refers to the same ortholog pair. The mouse matrix
    is re-indexed to use human Ensembl gene IDs so both matrices live in the
    same feature space — a prerequisite for Procrustes analysis.

    Args:
        human_adata: Human AnnData (var must contain 'feature_id' with Ensembl IDs).
        mouse_adata: Mouse AnnData (var must contain 'feature_id' with Ensembl IDs).
        orthologs: DataFrame from fetch_orthologs().

    Returns:
        Tuple (human_filtered, mouse_filtered) with identical var dimensions,
        both indexed by human Ensembl gene IDs.
    """
    # Validate inputs
    for label, adata in [("Human", human_adata), ("Mouse", mouse_adata)]:
        if "feature_id" not in adata.var.columns:
            raise ValueError(
                f"{label} AnnData missing 'feature_id' in var columns. "
                f"Available: {list(adata.var.columns)}"
            )

    # Extract gene ID sets
    human_gene_ids = set(human_adata.var["feature_id"].values)
    mouse_gene_ids = set(mouse_adata.var["feature_id"].values)

    print(f"  Human genes in data: {len(human_gene_ids):,}")
    print(f"  Mouse genes in data: {len(mouse_gene_ids):,}")
    print(f"  BioMart 1:1 ortholog pairs: {len(orthologs):,}")

    # Find orthologs present in BOTH datasets
    shared = orthologs[
        orthologs["human_ensembl_id"].isin(human_gene_ids)
        & orthologs["mouse_ensembl_id"].isin(mouse_gene_ids)
    ].copy()

    # Deduplicate: keep one mapping per human gene and per mouse gene
    shared = shared.drop_duplicates(subset=["human_ensembl_id"])
    shared = shared.drop_duplicates(subset=["mouse_ensembl_id"])
    shared = shared.reset_index(drop=True)

    n_shared = len(shared)
    print(f"  Shared 1:1 orthologs present in both datasets: {n_shared:,}")

    if n_shared < MIN_SHARED_GENES:
        print(
            f"  WARNING: {n_shared:,} < {MIN_SHARED_GENES:,} minimum. "
            f"Phase 1 gene gate will FAIL."
        )

    # Build lookup for mouse → human mapping
    mouse_to_human_id = dict(
        zip(shared["mouse_ensembl_id"], shared["human_ensembl_id"])
    )
    mouse_to_human_name = dict(
        zip(shared["mouse_ensembl_id"], shared["human_gene_name"])
    )

    # Filter human AnnData to shared genes
    human_mask = human_adata.var["feature_id"].isin(shared["human_ensembl_id"].values)
    human_filtered = human_adata[:, human_mask].copy()

    # Set human var index to feature_id for alignment
    human_filtered.var.index = human_filtered.var["feature_id"].values

    # Sort by feature_id for consistent ordering
    gene_order = sorted(shared["human_ensembl_id"].values)
    human_filtered = human_filtered[:, gene_order].copy()

    # Filter mouse AnnData to shared genes
    mouse_mask = mouse_adata.var["feature_id"].isin(shared["mouse_ensembl_id"].values)
    mouse_filtered = mouse_adata[:, mouse_mask].copy()

    # Save original mouse gene info before re-indexing
    mouse_filtered.var["original_mouse_feature_id"] = (
        mouse_filtered.var["feature_id"].values
    )
    mouse_filtered.var["original_mouse_feature_name"] = (
        mouse_filtered.var["feature_name"].values
    )

    # Re-index mouse var to human Ensembl IDs
    mouse_filtered.var["human_ensembl_id"] = [
        mouse_to_human_id[mid] for mid in mouse_filtered.var["feature_id"].values
    ]
    mouse_filtered.var["human_gene_name"] = [
        mouse_to_human_name[mid] for mid in mouse_filtered.var["feature_id"].values
    ]
    mouse_filtered.var.index = mouse_filtered.var["human_ensembl_id"].values

    # Sort mouse to match human gene order
    mouse_filtered = mouse_filtered[:, gene_order].copy()

    # Sanity check
    assert human_filtered.n_vars == mouse_filtered.n_vars, (
        f"Gene count mismatch after alignment: "
        f"human={human_filtered.n_vars}, mouse={mouse_filtered.n_vars}"
    )
    assert list(human_filtered.var.index) == list(mouse_filtered.var.index), (
        "Gene order mismatch after alignment"
    )

    print(
        f"  Final aligned shape: human {human_filtered.shape}, "
        f"mouse {mouse_filtered.shape}"
    )

    return human_filtered, mouse_filtered


def print_download_summary(
    human_adata: ad.AnnData,
    mouse_adata: ad.AnnData,
    cell_type_column: str = "cell_type",
) -> dict:
    """
    Print a human-readable summary and evaluate Phase 1 data gate criteria.

    Args:
        human_adata: Human AnnData after ortholog filtering.
        mouse_adata: Mouse AnnData after ortholog filtering.
        cell_type_column: Column name for cell types in .obs.

    Returns:
        Dict with summary statistics for programmatic use.
    """
    human_counts = human_adata.obs[cell_type_column].value_counts().to_dict()
    mouse_counts = mouse_adata.obs[cell_type_column].value_counts().to_dict()
    shared_genes = human_adata.n_vars

    # Collect all cell types
    all_types = sorted(set(list(human_counts.keys()) + list(mouse_counts.keys())))

    # Print table
    print(f"\n{'Cell Type':<40} {'Human':>8} {'Mouse':>8}")
    print("-" * 58)
    for ct in all_types:
        h = human_counts.get(ct, 0)
        m = mouse_counts.get(ct, 0)
        print(f"{ct:<40} {h:>8,} {m:>8,}")
    print("-" * 58)
    print(f"{'TOTAL':<40} {human_adata.n_obs:>8,} {mouse_adata.n_obs:>8,}")
    print(f"\nShared 1:1 ortholog genes: {shared_genes:,}")

    # Evaluate gate criteria
    gate_cells_pass = all(
        human_counts.get(ct, 0) >= MIN_CELLS_PER_TYPE
        and mouse_counts.get(ct, 0) >= MIN_CELLS_PER_TYPE
        for ct in all_types
    )
    gate_genes_pass = shared_genes >= MIN_SHARED_GENES

    print(f"\nGate check — cells ≥{MIN_CELLS_PER_TYPE}/type/species: "
          f"{'PASS' if gate_cells_pass else 'FAIL'}")
    if not gate_cells_pass:
        for ct in all_types:
            h = human_counts.get(ct, 0)
            m = mouse_counts.get(ct, 0)
            if h < MIN_CELLS_PER_TYPE or m < MIN_CELLS_PER_TYPE:
                print(f"  FAIL: {ct} — human={h:,}, mouse={m:,}")

    print(f"Gate check — shared genes ≥{MIN_SHARED_GENES:,}: "
          f"{'PASS' if gate_genes_pass else 'FAIL'} ({shared_genes:,})")

    return {
        "human_cells_per_type": human_counts,
        "mouse_cells_per_type": mouse_counts,
        "shared_genes": shared_genes,
        "gate_cells_pass": gate_cells_pass,
        "gate_genes_pass": gate_genes_pass,
    }


def save_h5ad_atomic(adata: ad.AnnData, path: Path) -> None:
    """
    Save AnnData to .h5ad using atomic write (tmp file + rename).

    Prevents corrupted checkpoints if the script is interrupted mid-write.

    Args:
        adata: AnnData to save.
        path: Final output path.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".h5ad.tmp")
    adata.write_h5ad(tmp)
    tmp.rename(path)
    print(f"  Saved: {path}")


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _dataset_id_filter(dataset_ids: list[str] | None) -> str:
    """Build a SOMA filter clause for dataset_ids, or empty string if None."""
    if not dataset_ids:
        return ""
    ids_str = ", ".join(f"'{d}'" for d in dataset_ids)
    return f"dataset_id in [{ids_str}]"
