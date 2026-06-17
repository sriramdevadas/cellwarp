"""
CellWarp — Cancer Data Loader Module

Downloads colorectal cancer (CRC) and matched normal colon tissue from CZ CELLxGENE
Census, harmonizes cell type labels to a coarse shared vocabulary, and prepares data
for Procrustes analysis of cancer-as-geometric-deformation.

Biology
-------
Colorectal cancer reshapes the cellular ecosystem of the colon. Tumor tissue acquires
new cell states (e.g., cancer-associated fibroblasts, exhausted T cells, tumor-associated
macrophages) and shifts existing cell types toward altered expression programs. By
applying the same Procrustes framework used for cross-species comparison, we can ask
whether normal-to-tumor transformation follows a coherent geometric deformation in
gene expression space — and whether the "deformation vectors" point to known cancer
biology.

Math
----
We download cells from two conditions (normal colon, CRC tumor) into the same 16,959
shared-ortholog gene space used in the cross-species pipeline. This enables direct
comparison of cancer deformation vectors to species deformation vectors. Coarse cell
type labels serve as landmarks for Procrustes alignment, analogous to homologous cell
types across species.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DATA_DIR = Path("./data/cancer")
OUTPUT_DIR = Path("./output/cancer")

ORGANISM = "Homo sapiens"

# Maximum cells to download per coarse cell type per condition
MAX_CELLS_PER_TYPE = 3_000
# Minimum cells required in BOTH conditions for a cell type to be included
MIN_CELLS_PER_TYPE = 500

RANDOM_SEED = 42

# Metadata columns to retain from Census
OBS_COLUMNS = [
    "cell_type",
    "tissue",
    "tissue_general",
    "disease",
    "donor_id",
    "dataset_id",
    "assay",
    "sex",
    "development_stage",
    "is_primary_data",
]

# Gene metadata columns
VAR_COLUMNS = ["feature_id", "feature_name"]

# ---------------------------------------------------------------------------
# Coarse cell type mapping
# ---------------------------------------------------------------------------

# Maps raw CELLxGENE cell_type labels to 10 coarse categories.
# Keys are lowercase substrings to match against raw labels (case-insensitive).
# Order matters: first match wins. More specific patterns come first.
COARSE_RULES: list[tuple[str, list[str]]] = [
    # IMPORTANT: Mast cell and Epithelial cell MUST come before T cell.
    # "t cell" is a substring of "goblet cell" and "mast cell", so T cell's
    # keyword would incorrectly match these if checked first.
    ("Mast cell", [
        "mast cell",
    ]),
    ("Epithelial cell", [
        "epithelial", "enterocyte", "colonocyte", "goblet", "paneth",
        "tuft", "enteroendocrine", "transit amplifying", "stem cell of colon",
        "intestinal crypt stem cell", "colon epithelial",
    ]),
    # Myofibroblast before T cell — "myofibroblast cell" contains substring
    # "t cell", so must match fibroblast first (ISSUE-054 fix)
    ("Fibroblast", ["myofibroblast"]),
    ("T cell", [
        "cd4", "cd8", "t cell", "regulatory t", "nkt", "gamma-delta",
        "thymocyte", "treg",
    ]),
    ("B cell", [
        "b cell", "plasma cell", "plasmablast", "plasma blast",
    ]),
    ("Macrophage", [
        "macrophage", "monocyte", "dendritic", "langerhans", "kupffer",
        "microglia", "histiocyte",
    ]),
    ("NK cell", [
        "natural killer", "nk cell",
    ]),
    ("Fibroblast", [
        "fibroblast", "stromal", "caf", "cancer associated fibroblast",
    ]),
    ("Endothelial cell", [
        "endothelial",
    ]),
    ("Smooth muscle cell", [
        "smooth muscle",
    ]),
]


# ---------------------------------------------------------------------------
# Census inventory functions
# ---------------------------------------------------------------------------


def query_tissue_values(
    census: Any,
    tissue_keywords: list[str] | None = None,
) -> pd.DataFrame:
    """
    Query all distinct tissue and tissue_general values in Census for Homo sapiens,
    optionally filtered to those containing specific keywords.

    Biology: CELLxGENE annotates each cell with both a specific tissue label
    (e.g., "transverse colon") and a general tissue category (e.g., "large intestine").
    We need to discover the exact strings used before building download filters.

    Args:
        census: Open Census SOMA collection.
        tissue_keywords: If provided, only return rows where tissue or
            tissue_general contains one of these substrings (case-insensitive).

    Returns:
        DataFrame with columns: tissue, tissue_general, cell_count.
    """
    import cellxgene_census

    # Query all tissue values for human primary data
    obs_df = cellxgene_census.get_obs(
        census,
        ORGANISM,
        value_filter="is_primary_data == True",
        column_names=["tissue", "tissue_general"],
    )

    # Count cells per tissue / tissue_general combination
    counts = (
        obs_df.groupby(["tissue_general", "tissue"])
        .size()
        .reset_index(name="cell_count")
        .sort_values("cell_count", ascending=False)
    )

    if tissue_keywords:
        mask = pd.Series(False, index=counts.index)
        for kw in tissue_keywords:
            kw_lower = kw.lower()
            mask |= counts["tissue"].str.lower().str.contains(kw_lower, na=False)
            mask |= counts["tissue_general"].str.lower().str.contains(
                kw_lower, na=False
            )
        counts = counts[mask]

    return counts.reset_index(drop=True)


def query_disease_values(
    census: Any,
    tissue_general_values: list[str],
) -> pd.DataFrame:
    """
    Query all distinct disease values for cells in specified tissue_general categories.

    Biology: CELLxGENE disease annotations use standardized ontology terms.
    We need to discover which disease labels correspond to colorectal cancer
    vs normal tissue before filtering.

    Args:
        census: Open Census SOMA collection.
        tissue_general_values: List of tissue_general values to query
            (e.g., ["large intestine"]).

    Returns:
        DataFrame with columns: disease, cell_count, sorted by count descending.
    """
    import cellxgene_census

    tissue_str = ", ".join(f"'{v}'" for v in tissue_general_values)
    value_filter = (
        f"is_primary_data == True and "
        f"tissue_general in [{tissue_str}]"
    )

    obs_df = cellxgene_census.get_obs(
        census,
        ORGANISM,
        value_filter=value_filter,
        column_names=["disease"],
    )

    counts = (
        obs_df["disease"]
        .value_counts()
        .reset_index()
    )
    counts.columns = ["disease", "cell_count"]

    return counts


def query_cell_type_inventory(
    census: Any,
    tissue_general_values: list[str],
    disease_normal: str,
    disease_tumor: str,
) -> pd.DataFrame:
    """
    Query cell type labels and counts for both normal and tumor conditions.

    Biology: Different conditions may have different cell type compositions.
    Tumors typically gain cancer-associated fibroblasts, exhausted T cells, and
    tumor-associated macrophages, while losing normal epithelial subtypes.

    Args:
        census: Open Census SOMA collection.
        tissue_general_values: Tissue categories to include.
        disease_normal: Disease label for healthy tissue (e.g., "normal").
        disease_tumor: Disease label for CRC (e.g., "colorectal cancer").

    Returns:
        DataFrame with columns: cell_type, normal_count, tumor_count.
    """
    import cellxgene_census

    tissue_str = ", ".join(f"'{v}'" for v in tissue_general_values)
    results = {}

    for condition, disease_label in [
        ("normal", disease_normal),
        ("tumor", disease_tumor),
    ]:
        value_filter = (
            f"is_primary_data == True and "
            f"tissue_general in [{tissue_str}] and "
            f"disease == '{disease_label}'"
        )
        obs_df = cellxgene_census.get_obs(
            census,
            ORGANISM,
            value_filter=value_filter,
            column_names=["cell_type"],
        )
        counts = obs_df["cell_type"].value_counts().to_dict()
        results[condition] = counts

    # Merge into a single DataFrame
    all_types = sorted(
        set(list(results["normal"].keys()) + list(results["tumor"].keys()))
    )
    rows = []
    for ct in all_types:
        rows.append({
            "cell_type": ct,
            "normal_count": results["normal"].get(ct, 0),
            "tumor_count": results["tumor"].get(ct, 0),
        })

    df = pd.DataFrame(rows)
    df = df.sort_values("cell_type").reset_index(drop=True)
    return df


def save_census_inventory(
    tissue_df: pd.DataFrame,
    disease_df: pd.DataFrame,
    cell_type_df: pd.DataFrame,
    output_path: Path,
) -> None:
    """
    Save the full census inventory to a human-readable text file.

    Args:
        tissue_df: Tissue values from query_tissue_values().
        disease_df: Disease values from query_disease_values().
        cell_type_df: Cell type inventory from query_cell_type_inventory().
        output_path: Path to save the inventory file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        f.write("=" * 80 + "\n")
        f.write("CellWarp — CELLxGENE Census Inventory for Colon/CRC Analysis\n")
        f.write("=" * 80 + "\n\n")

        f.write("TISSUE VALUES (colon/intestine related)\n")
        f.write("-" * 80 + "\n")
        f.write(f"{'tissue_general':<30} {'tissue':<35} {'cells':>10}\n")
        f.write("-" * 80 + "\n")
        for _, row in tissue_df.iterrows():
            f.write(
                f"{str(row['tissue_general']):<30} "
                f"{str(row['tissue']):<35} "
                f"{row['cell_count']:>10,}\n"
            )

        f.write(f"\n\nDISEASE VALUES (for colon tissue)\n")
        f.write("-" * 80 + "\n")
        f.write(f"{'disease':<50} {'cells':>10}\n")
        f.write("-" * 80 + "\n")
        for _, row in disease_df.iterrows():
            f.write(f"{str(row['disease']):<50} {row['cell_count']:>10,}\n")

        f.write(f"\n\nCELL TYPE INVENTORY (normal vs tumor)\n")
        f.write("-" * 80 + "\n")
        f.write(f"{'cell_type':<45} {'normal':>10} {'tumor':>10}\n")
        f.write("-" * 80 + "\n")
        for _, row in cell_type_df.iterrows():
            f.write(
                f"{str(row['cell_type']):<45} "
                f"{row['normal_count']:>10,} "
                f"{row['tumor_count']:>10,}\n"
            )
        f.write("-" * 80 + "\n")
        f.write(
            f"{'TOTAL':<45} "
            f"{cell_type_df['normal_count'].sum():>10,} "
            f"{cell_type_df['tumor_count'].sum():>10,}\n"
        )

    print(f"  Saved census inventory: {output_path}")


# ---------------------------------------------------------------------------
# Coarse cell type mapping
# ---------------------------------------------------------------------------


def build_coarse_mapping(
    cell_type_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Map raw CELLxGENE cell_type labels to 10 coarse categories using substring rules.

    Biology: CELLxGENE uses fine-grained Cell Ontology labels (e.g., "CD8-positive,
    alpha-beta cytokine secreting effector T cell"). For Procrustes analysis we need
    matched landmarks between conditions, so we collapse to coarse categories that
    are shared between normal and tumor tissue.

    The mapping uses ordered substring matching: for each raw label, we check
    COARSE_RULES in order and assign the first match. Labels that match no rule
    are assigned "Other".

    Args:
        cell_type_df: DataFrame with columns: cell_type, normal_count, tumor_count
            (from query_cell_type_inventory).

    Returns:
        DataFrame with columns: raw_label, coarse_label, normal_count, tumor_count,
        sorted by coarse_label then raw_label.
    """
    rows = []
    for _, row in cell_type_df.iterrows():
        raw = row["cell_type"]
        coarse = _classify_cell_type(raw)
        rows.append({
            "raw_label": raw,
            "coarse_label": coarse,
            "normal_count": row["normal_count"],
            "tumor_count": row["tumor_count"],
        })

    df = pd.DataFrame(rows)
    df = df.sort_values(["coarse_label", "raw_label"]).reset_index(drop=True)
    return df


def print_coarse_mapping_table(mapping_df: pd.DataFrame) -> pd.DataFrame:
    """
    Print the raw→coarse mapping as a table and return aggregated coarse counts.

    Displays both the per-raw-label mapping and the aggregated coarse-level counts,
    including which coarse types pass the >=500 gate in both conditions.

    Args:
        mapping_df: DataFrame from build_coarse_mapping().

    Returns:
        DataFrame with aggregated coarse counts (coarse_label, normal_count,
        tumor_count, passes_gate).
    """
    print("\n" + "=" * 90)
    print("COARSE CELL TYPE MAPPING — raw_label → coarse_label")
    print("=" * 90)
    print(f"{'raw_label':<55} {'coarse_label':<20} {'normal':>7} {'tumor':>7}")
    print("-" * 90)

    current_coarse = None
    for _, row in mapping_df.iterrows():
        if row["coarse_label"] != current_coarse:
            if current_coarse is not None:
                print()  # Blank line between groups
            current_coarse = row["coarse_label"]
        print(
            f"{str(row['raw_label']):<55} "
            f"{row['coarse_label']:<20} "
            f"{row['normal_count']:>7,} "
            f"{row['tumor_count']:>7,}"
        )

    # Aggregate to coarse level
    agg = (
        mapping_df.groupby("coarse_label")[["normal_count", "tumor_count"]]
        .sum()
        .reset_index()
    )
    agg["passes_gate"] = (
        (agg["normal_count"] >= MIN_CELLS_PER_TYPE)
        & (agg["tumor_count"] >= MIN_CELLS_PER_TYPE)
    )

    print("\n\n" + "=" * 70)
    print("AGGREGATED COARSE COUNTS")
    print("=" * 70)
    print(f"{'coarse_label':<25} {'normal':>10} {'tumor':>10} {'gate':>8}")
    print("-" * 70)
    for _, row in agg.sort_values("coarse_label").iterrows():
        gate_str = "PASS" if row["passes_gate"] else "FAIL"
        print(
            f"{row['coarse_label']:<25} "
            f"{row['normal_count']:>10,} "
            f"{row['tumor_count']:>10,} "
            f"{'  ' + gate_str:>8}"
        )

    n_pass = agg["passes_gate"].sum()
    n_total = len(agg)
    print("-" * 70)
    print(f"  {n_pass}/{n_total} coarse types pass >=500 gate in both conditions")

    return agg


# ---------------------------------------------------------------------------
# Download functions
# ---------------------------------------------------------------------------


def download_condition(
    census: Any,
    tissue_general_values: list[str],
    disease_label: str,
    coarse_mapping: dict[str, str],
    valid_coarse_types: list[str],
    ortholog_gene_ids: list[str],
    max_cells_per_type: int = MAX_CELLS_PER_TYPE,
    random_seed: int = RANDOM_SEED,
) -> ad.AnnData:
    """
    Download cells for one condition using chunked streaming to avoid OOM.

    Instead of materializing the full expression matrix for all matching cells
    (which can be hundreds of thousands of cells x ~60k genes and exceed M-series
    MacBook memory), this function:

      1. Fetches obs metadata ONLY (no X matrix) to discover cell IDs per coarse type.
      2. Subsamples soma_joinids per coarse type to ≤max_cells_per_type.
      3. Fetches expression data ONE coarse type at a time using obs_coords.
      4. Concatenates the per-type AnnData objects.

    Biology: We subsample uniformly across fine-grained subtypes within each coarse
    category, preserving subtype diversity. Fetching per-type keeps peak memory
    proportional to max_cells_per_type × n_genes rather than total_cells × n_genes.

    Math: Cells are subsampled uniformly at random (seeded) within each coarse type.
    The ortholog gene filter restricts to the same 16,959-gene space used in the
    cross-species pipeline.

    Args:
        census: Open Census SOMA collection.
        tissue_general_values: Tissue categories to include.
        disease_label: Disease value to filter on.
        coarse_mapping: Dict mapping raw cell_type → coarse_label.
        valid_coarse_types: Coarse types that passed the >=500 gate.
        ortholog_gene_ids: Human Ensembl gene IDs from ortholog table (the shared
            gene space).
        max_cells_per_type: Maximum cells per coarse type.
        random_seed: For reproducible subsampling.

    Returns:
        AnnData filtered to ortholog genes, with 'coarse_cell_type' in obs,
        subsampled to max_cells_per_type per coarse type.
    """
    import cellxgene_census

    tissue_str = ", ".join(f"'{v}'" for v in tissue_general_values)
    value_filter = (
        f"is_primary_data == True and "
        f"tissue_general in [{tissue_str}] and "
        f"disease == '{disease_label}'"
    )

    # ------------------------------------------------------------------
    # Phase 1: Fetch obs metadata only (lightweight — no expression data)
    # ------------------------------------------------------------------
    print(f"  Fetching obs metadata for disease='{disease_label}'...")
    _log_memory("before metadata fetch")

    obs_df = cellxgene_census.get_obs(
        census,
        ORGANISM,
        value_filter=value_filter,
        column_names=list(dict.fromkeys(["soma_joinid"] + OBS_COLUMNS)),
    )
    print(f"  Metadata fetched: {len(obs_df):,} cells (no expression matrix)")
    _log_memory("after metadata fetch")

    # Apply coarse mapping
    obs_df["coarse_cell_type"] = obs_df["cell_type"].map(coarse_mapping)
    unmapped = obs_df["coarse_cell_type"].isna().sum()
    if unmapped > 0:
        obs_df["coarse_cell_type"] = obs_df["coarse_cell_type"].fillna("Other")
        print(f"  WARNING: {unmapped:,} cells had no coarse mapping -> assigned 'Other'")

    # Filter to valid coarse types
    obs_df = obs_df[obs_df["coarse_cell_type"].isin(valid_coarse_types)].copy()
    print(f"  After coarse type filter: {len(obs_df):,} cells")

    # ------------------------------------------------------------------
    # Phase 2: Subsample soma_joinids per coarse type (metadata only)
    # ------------------------------------------------------------------
    rng = np.random.default_rng(random_seed)
    selected_ids_by_type: dict[str, np.ndarray] = {}

    print(f"  Subsampling to <={max_cells_per_type:,} cells per coarse type (seed={random_seed})...")
    for ct in sorted(valid_coarse_types):
        ct_joinids = obs_df.loc[obs_df["coarse_cell_type"] == ct, "soma_joinid"].values
        n_total = len(ct_joinids)

        if n_total == 0:
            print(f"    {ct}: 0 cells (skipped)")
            continue

        n_select = min(n_total, max_cells_per_type)
        if n_total > max_cells_per_type:
            selected = rng.choice(ct_joinids, size=n_select, replace=False)
            print(f"    {ct}: {n_total:,} -> {n_select:,} (subsampled)")
        else:
            selected = ct_joinids
            print(f"    {ct}: {n_total:,} (kept all)")

        selected_ids_by_type[ct] = np.sort(selected)

    # ------------------------------------------------------------------
    # Phase 3: Fetch expression data ONE coarse type at a time
    # ------------------------------------------------------------------
    print(f"\n  Downloading expression data per coarse type...")
    chunks: list[ad.AnnData] = []

    for ct in sorted(selected_ids_by_type.keys()):
        joinids = selected_ids_by_type[ct]
        _log_memory(f"before fetching {ct}")
        print(f"    Fetching {ct}: {len(joinids):,} cells...")

        chunk = cellxgene_census.get_anndata(
            census=census,
            organism=ORGANISM,
            obs_coords=joinids,
            obs_column_names=OBS_COLUMNS,
            var_column_names=VAR_COLUMNS,
        )
        # Tag with coarse label before filtering genes
        chunk.obs["coarse_cell_type"] = ct
        chunks.append(chunk)

        print(f"    {ct}: {chunk.n_obs:,} cells x {chunk.n_vars:,} genes")
        _log_memory(f"after fetching {ct}")

    # ------------------------------------------------------------------
    # Phase 4: Concatenate and filter to ortholog gene space
    # ------------------------------------------------------------------
    print(f"\n  Concatenating {len(chunks)} chunks...")
    adata = ad.concat(chunks, merge="same")
    del chunks  # free memory
    print(f"  Combined: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    _log_memory("after concat")

    # Filter to ortholog gene space
    gene_mask = adata.var["feature_id"].isin(ortholog_gene_ids)
    n_genes_before = adata.n_vars
    adata = adata[:, gene_mask].copy()
    print(
        f"  Filtered to ortholog genes: {n_genes_before:,} -> {adata.n_vars:,} genes"
    )

    # Set var index to feature_id (human Ensembl IDs) for consistency
    adata.var.index = adata.var["feature_id"].values

    # Sort to match the ortholog gene order
    shared_ids = sorted(set(adata.var.index) & set(ortholog_gene_ids))
    adata = adata[:, shared_ids].copy()

    _log_memory("final")
    return adata


# ---------------------------------------------------------------------------
# Normalization
# ---------------------------------------------------------------------------


def normalize_cancer_data(
    adata: ad.AnnData,
    min_genes: int = 200,
    target_sum: float = 10_000,
) -> ad.AnnData:
    """
    Apply identical normalization to the cross-species pipeline (src/qc.py).

    Steps:
      1. Filter cells with < min_genes detected genes
      2. Normalize counts per target_sum (10,000) + log1p

    We do NOT compute HVGs or PCA here — centroid computation will use the full
    ortholog gene space, consistent with Phase 2 (DECISION-023).

    Biology: Normalization corrects for sequencing depth differences between cells
    so expression levels are comparable. The same normalization as the cross-species
    pipeline ensures cancer deformation vectors are directly comparable to species
    deformation vectors.

    Math: x_norm = x / total_counts * 10,000; x_log = log(1 + x_norm).

    Args:
        adata: Raw count AnnData in ortholog gene space.
        min_genes: Minimum detected genes per cell.
        target_sum: Normalization target count depth.

    Returns:
        Normalized AnnData (filtered, counts per 10k, log1p).
    """
    n_before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=min_genes)
    n_after = adata.n_obs
    n_removed = n_before - n_after
    print(
        f"  QC filter (<{min_genes} genes): "
        f"{n_before:,} -> {n_after:,} ({n_removed:,} removed)"
    )

    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    print(f"  Normalized to {target_sum:,} counts per cell + log1p")

    return adata


# ---------------------------------------------------------------------------
# Summary and reporting
# ---------------------------------------------------------------------------


def print_final_summary(
    normal_adata: ad.AnnData,
    tumor_adata: ad.AnnData,
    dropped_types: list[str],
) -> dict:
    """
    Print a comprehensive summary of the downloaded cancer data.

    Reports:
      - Cells per coarse cell type (normal vs tumor)
      - Dropped cell types
      - Donor counts per condition (flags if < 3)
      - Assay types (flags mismatches between conditions)

    Args:
        normal_adata: Normalized normal colon AnnData.
        tumor_adata: Normalized tumor AnnData.
        dropped_types: Coarse types that failed the >=500 gate.

    Returns:
        Dict with summary statistics for programmatic use.
    """
    print("\n" + "=" * 70)
    print("FINAL DOWNLOAD SUMMARY")
    print("=" * 70)

    # Cell counts per coarse type
    normal_counts = normal_adata.obs["coarse_cell_type"].value_counts().to_dict()
    tumor_counts = tumor_adata.obs["coarse_cell_type"].value_counts().to_dict()
    all_types = sorted(set(list(normal_counts.keys()) + list(tumor_counts.keys())))

    print(f"\n{'Coarse Cell Type':<25} {'Normal':>10} {'Tumor':>10}")
    print("-" * 50)
    for ct in all_types:
        n = normal_counts.get(ct, 0)
        t = tumor_counts.get(ct, 0)
        print(f"{ct:<25} {n:>10,} {t:>10,}")
    print("-" * 50)
    print(
        f"{'TOTAL':<25} {normal_adata.n_obs:>10,} {tumor_adata.n_obs:>10,}"
    )
    print(f"Genes: {normal_adata.n_vars:,}")

    # Dropped types
    if dropped_types:
        print(f"\nDropped cell types (failed >=500 gate):")
        for ct in dropped_types:
            print(f"  - {ct}")
    else:
        print(f"\nNo cell types dropped (all passed >=500 gate).")

    # Donor counts
    print(f"\nDonor counts:")
    normal_donors = normal_adata.obs["donor_id"].nunique()
    tumor_donors = tumor_adata.obs["donor_id"].nunique()
    flag_n = " *** WARNING: < 3 donors!" if normal_donors < 3 else ""
    flag_t = " *** WARNING: < 3 donors!" if tumor_donors < 3 else ""
    print(f"  Normal: {normal_donors} unique donors{flag_n}")
    print(f"  Tumor:  {tumor_donors} unique donors{flag_t}")

    # Per-type donor counts
    print(f"\n  Per-type donor counts (normal / tumor):")
    for ct in all_types:
        n_donors = (
            normal_adata.obs.loc[
                normal_adata.obs["coarse_cell_type"] == ct, "donor_id"
            ].nunique()
            if ct in normal_counts
            else 0
        )
        t_donors = (
            tumor_adata.obs.loc[
                tumor_adata.obs["coarse_cell_type"] == ct, "donor_id"
            ].nunique()
            if ct in tumor_counts
            else 0
        )
        flag = ""
        if n_donors < 3 or t_donors < 3:
            flag = " ***"
        print(f"    {ct:<25} {n_donors:>5} / {t_donors:>5}{flag}")

    # Assay types
    print(f"\nAssay types:")
    normal_assays = sorted(normal_adata.obs["assay"].unique())
    tumor_assays = sorted(tumor_adata.obs["assay"].unique())
    print(f"  Normal: {normal_assays}")
    print(f"  Tumor:  {tumor_assays}")
    if set(normal_assays) != set(tumor_assays):
        print(
            "  *** WARNING: Assay types MISMATCH between conditions. "
            "This is a potential confound."
        )
        only_normal = set(normal_assays) - set(tumor_assays)
        only_tumor = set(tumor_assays) - set(normal_assays)
        if only_normal:
            print(f"    Only in normal: {sorted(only_normal)}")
        if only_tumor:
            print(f"    Only in tumor: {sorted(only_tumor)}")
    else:
        print("  Assay types match between conditions.")

    # Dataset sources
    print(f"\nDataset sources:")
    for condition, adata in [("Normal", normal_adata), ("Tumor", tumor_adata)]:
        datasets = adata.obs["dataset_id"].unique()
        print(f"  {condition}: {len(datasets)} dataset(s)")
        for ds in sorted(datasets):
            n = (adata.obs["dataset_id"] == ds).sum()
            print(f"    {ds[:16]}... : {n:,} cells")

    return {
        "normal_counts": normal_counts,
        "tumor_counts": tumor_counts,
        "normal_donors": normal_donors,
        "tumor_donors": tumor_donors,
        "normal_assays": normal_assays,
        "tumor_assays": tumor_assays,
        "dropped_types": dropped_types,
        "normal_total": normal_adata.n_obs,
        "tumor_total": tumor_adata.n_obs,
        "n_genes": normal_adata.n_vars,
    }


# ---------------------------------------------------------------------------
# File I/O
# ---------------------------------------------------------------------------


def save_h5ad_atomic(adata: ad.AnnData, path: Path) -> None:
    """
    Save AnnData to .h5ad using atomic write (tmp file + rename).

    Prevents corrupted files if the script is interrupted mid-write.

    Args:
        adata: AnnData to save.
        path: Final output path.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".h5ad.tmp")
    adata.write_h5ad(tmp)
    tmp.rename(path)
    print(f"  Saved: {path}")


def load_ortholog_gene_ids(ortholog_path: Path) -> list[str]:
    """
    Load human Ensembl gene IDs from the ortholog mapping table.

    These define the shared gene space used in both the cross-species pipeline
    and the cancer analysis, ensuring deformation vectors are directly comparable.

    Args:
        ortholog_path: Path to orthologs_human_mouse.csv.

    Returns:
        Sorted list of human Ensembl gene IDs.
    """
    df = pd.read_csv(ortholog_path)
    gene_ids = sorted(df["human_ensembl_id"].unique().tolist())
    print(f"  Loaded {len(gene_ids):,} human ortholog gene IDs from {ortholog_path}")
    return gene_ids


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _classify_cell_type(raw_label: str) -> str:
    """
    Classify a raw cell type label into a coarse category using substring matching.

    Args:
        raw_label: Raw cell_type string from CELLxGENE.

    Returns:
        Coarse category string (one of the 10 categories or "Other").
    """
    label_lower = raw_label.lower()
    for coarse_name, keywords in COARSE_RULES:
        for kw in keywords:
            if kw in label_lower:
                return coarse_name
    return "Other"


def _log_memory(label: str) -> None:
    """
    Print current process RSS memory usage if psutil is available.

    Silently skips if psutil is not installed — memory logging is optional
    diagnostics, not a hard dependency.

    Args:
        label: Short description of the checkpoint (e.g., "after fetching T cell").
    """
    try:
        import psutil

        rss_mb = psutil.Process().memory_info().rss / (1024 * 1024)
        print(f"    [mem] {label}: {rss_mb:.0f} MB RSS")
    except ImportError:
        pass
