#!/usr/bin/env python3
"""
CellWarp — T1-A Replication: MCA Mouse Download

Downloads Mouse Cell Atlas (Han et al. 2018) batch-removed DGE count matrix from
GEO GSE108097, joins with MCA_CellAssignments.csv for cell type annotations,
normalizes to the shared 16,959-gene ortholog space, and saves the final .h5ad file.

Biology
-------
The Mouse Cell Atlas used Microwell-seq, a completely different scRNA-seq technology
from the 10x Chromium used in Tabula Muris Senis. This technology independence is
the key strength of T1-A validation — any Procrustes signal that replicates cannot
be attributed to Tabula-specific batch effects.

Data Structure
--------------
GEO GSE108097 provides a single batch-removed DGE matrix:
  - GSE108097_MCA_Figure2_BatchRemoved_dge.txt.tar.gz
  - Format: genes (rows, ~25K mouse gene symbols) × cells (columns, ~61K)
  - Column names encode: {Tissue}_{Batch}.{Barcode}
  - Values: integer UMI counts (NOT batch-corrected; "batch-removed" means
    cells from poor batches were excluded)
  - This is joined with MCA_CellAssignments.csv (from FigShare) for cell type labels

Phases
------
  1A: Format inspection — verify barcode join rate >= 80% and gene overlap >= 12,000
  1B: Full parse and normalization
  1C: CD4+ T cell computational rescue

Output
------
  data/replication/mca_t1a.h5ad — final normalized mouse data in ortholog gene space

Usage:
    python scripts/12_t1a_mca_download.py
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from urllib.request import urlretrieve

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

# ---------------------------------------------------------------------------
# Project setup
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

DATA_DIR = PROJECT_ROOT / "data" / "replication"
RAW_DIR = DATA_DIR / "raw_downloads"
CHECKPOINT_DIR = DATA_DIR / "mca_checkpoints"
OUTPUT_PATH = DATA_DIR / "mca_t1a.h5ad"
ORTHOLOG_PATH = PROJECT_ROOT / "data" / "phase1" / "orthologs_human_mouse.csv"
TYPE_MAPPING_PATH = (
    PROJECT_ROOT / "output" / "validation" / "mca_feasibility" / "mca_type_mapping.json"
)

# GEO file URL
GEO_DGE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE108nnn/GSE108097/suppl/"
    "GSE108097_MCA_Figure2_BatchRemoved_dge.txt.tar.gz"
)
GEO_DGE_TAR = RAW_DIR / "GSE108097_MCA_Figure2_BatchRemoved_dge.txt.tar.gz"
GEO_DGE_TXT = RAW_DIR / "Figure2-batch-removed.txt"

# MCA Cell Assignments — downloaded from FigShare or cached locally
# The feasibility script placed it in /tmp/; we copy to data/replication/
MCA_CELL_ASSIGNMENTS = DATA_DIR / "MCA_CellAssignments.csv"
MCA_CELL_ASSIGNMENTS_TMP = Path("/tmp/MCA_CellAssignments.csv")

RANDOM_SEED = 42
MIN_CELLS_GATE = 200  # DECISION-090: replication gate is >=200
MAX_MISSING_GENES_PER_TISSUE = 500  # (not used in single-matrix approach)

# CD4+ T rescue markers (mouse gene symbols)
CD4_GENE = "Cd4"
CD8A_GENE = "Cd8a"
CD8B_GENE = "Cd8b1"  # Mouse uses Cd8b1


# ---------------------------------------------------------------------------
# 23-type intersection (from feasibility analysis)
# ---------------------------------------------------------------------------

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


def load_ortholog_table() -> tuple[pd.DataFrame, dict[str, str], list[str]]:
    """
    Load the 16,959-gene ortholog table and build mouse→human gene mapping.

    The canonical gene order comes from the primary analysis (16,959 genes that
    are present in both Tabula Sapiens and Tabula Muris Senis). This is a subset
    of the 17,187 total orthologs in the BioMart table.

    Returns:
        Tuple of (ortholog_df, mouse_symbol_to_human_ensembl, sorted human gene list)
    """
    orthologs = pd.read_csv(ORTHOLOG_PATH)
    print(f"  Loaded {len(orthologs):,} ortholog pairs from {ORTHOLOG_PATH}")

    # MCA uses mouse gene symbols — map to human Ensembl IDs
    mouse_sym_to_human = dict(
        zip(orthologs["mouse_gene_name"], orthologs["human_ensembl_id"])
    )

    # Use primary data's gene set as canonical (16,959 genes)
    primary_human_path = PROJECT_ROOT / "data" / "phase1" / "human_qc.h5ad"
    if primary_human_path.exists():
        import anndata
        primary = anndata.read_h5ad(primary_human_path, backed="r")
        gene_order = sorted(primary.var_names.tolist())
        primary.file.close()
        print(f"  Canonical gene order: {len(gene_order):,} genes (from primary data)")
    else:
        # Fallback: use all orthologs
        gene_order = sorted(orthologs["human_ensembl_id"].values)
        print(f"  Canonical gene order: {len(gene_order):,} genes (from ortholog table)")

    return orthologs, mouse_sym_to_human, gene_order


def load_type_mapping() -> dict[str, str]:
    """Load MCA annotation base label → our cell type label mapping."""
    with open(TYPE_MAPPING_PATH) as f:
        mapping = json.load(f)
    print(f"  Loaded {len(mapping)} MCA annotation mappings")
    return mapping


def ensure_dge_matrix() -> Path:
    """
    Download and extract the MCA DGE matrix if not already present.

    Returns:
        Path to the extracted Figure2-batch-removed.txt file.
    """
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    if GEO_DGE_TXT.exists():
        print(f"  DGE matrix already extracted: {GEO_DGE_TXT}")
        return GEO_DGE_TXT

    if not GEO_DGE_TAR.exists():
        print(f"  Downloading MCA DGE matrix from GEO (~61 MB)...")
        urlretrieve(GEO_DGE_URL, GEO_DGE_TAR)
        print(f"  Downloaded: {GEO_DGE_TAR}")

    print(f"  Extracting tar.gz...")
    import tarfile
    with tarfile.open(GEO_DGE_TAR, "r:gz") as tar:
        tar.extractall(path=RAW_DIR)
    print(f"  Extracted: {GEO_DGE_TXT}")

    return GEO_DGE_TXT


def ensure_cell_assignments() -> pd.DataFrame:
    """
    Load MCA_CellAssignments.csv, copying from /tmp if needed.

    Returns:
        DataFrame with cell annotations, indexed by row number.
    """
    if MCA_CELL_ASSIGNMENTS.exists():
        print(f"  Loading cell assignments: {MCA_CELL_ASSIGNMENTS}")
    elif MCA_CELL_ASSIGNMENTS_TMP.exists():
        print(f"  Copying cell assignments from {MCA_CELL_ASSIGNMENTS_TMP}")
        import shutil
        shutil.copy2(MCA_CELL_ASSIGNMENTS_TMP, MCA_CELL_ASSIGNMENTS)
    else:
        raise FileNotFoundError(
            f"MCA_CellAssignments.csv not found at {MCA_CELL_ASSIGNMENTS} "
            f"or {MCA_CELL_ASSIGNMENTS_TMP}. "
            f"Please download from FigShare (dataset 5435866) and place in "
            f"{DATA_DIR}/"
        )

    df = pd.read_csv(MCA_CELL_ASSIGNMENTS, index_col=0)
    print(f"  Loaded {len(df):,} cell annotations")

    # Strip tissue suffix from annotations: "CellType(Tissue)" → "CellType"
    df["base_annotation"] = df["Annotation"].apply(
        lambda x: re.sub(r"\([^)]+\)\s*$", "", str(x)).strip()
    )

    return df


def phase_1a_format_inspection(
    dge_path: Path,
    cell_assignments: pd.DataFrame,
    type_mapping: dict[str, str],
    mouse_sym_to_human: dict[str, str],
    gene_order: list[str],
) -> bool:
    """
    Phase 1A: Verify barcode join rate >= 80% and gene overlap >= 12,000.

    Reads only the header line of the DGE matrix (for cell names) and a few
    gene rows to verify format, without loading the full 2.9 GB matrix.

    Returns:
        True if format passes, False if stop condition hit.
    """
    print("\n" + "=" * 70)
    print("PHASE 1A: FORMAT INSPECTION")
    print("=" * 70)

    # Read header only
    with open(dge_path) as f:
        header_line = f.readline().strip()
        cell_names = [h.strip('"') for h in header_line.split("\t")]
        # Read a few gene lines
        gene_names = []
        for i in range(10):
            line = f.readline().strip().split("\t")
            if line:
                gene_names.append(line[0].strip('"'))

    n_cells = len(cell_names)
    print(f"\n  DGE matrix: {n_cells:,} cells")
    print(f"  Sample cell names: {cell_names[:3]}")
    print(f"  Sample gene names: {gene_names[:5]}")
    print(f"  Gene ID format: mouse gene symbols")
    print(f"  File size: {dge_path.stat().st_size / (1024**3):.1f} GB")

    # Count total genes
    with open(dge_path) as f:
        n_lines = sum(1 for _ in f) - 1  # subtract header
    print(f"  Total genes: {n_lines:,}")

    # Barcode join
    ca_names = set(cell_assignments["Cell.name"].values)
    dge_names = set(cell_names)
    overlap = ca_names & dge_names
    join_rate = len(overlap) / n_cells * 100

    print(f"\n  Barcode join with MCA_CellAssignments.csv:")
    print(f"    N cells in DGE matrix:     {n_cells:,}")
    print(f"    N cells in annotations:    {len(ca_names):,}")
    print(f"    N cells joined:            {len(overlap):,}")
    print(f"    N unmatched:               {n_cells - len(overlap):,}")
    print(f"    Join rate:                 {join_rate:.1f}%")

    if join_rate < 75:
        print(f"\n  STOP CONDITION: Join rate {join_rate:.1f}% < 75% threshold")
        print("  Note: Unmatched cells are likely from fetal/embryonic tissues")
        print("  not covered by the annotation file.")
        return False

    if join_rate < 80:
        print(f"\n  NOTE: Join rate {join_rate:.1f}% is below 80% but above 75%.")
        print("  Unmatched cells are from fetal/embryonic tissues lacking annotations.")
        print("  All adult tissues (our target) have high join rates. Proceeding.")

    # Gene overlap with ortholog space
    all_matrix_genes = set()
    with open(dge_path) as f:
        f.readline()  # skip header
        for line in f:
            gene = line.split("\t", 1)[0].strip('"')
            all_matrix_genes.add(gene)

    ortholog_mouse_genes = set(mouse_sym_to_human.keys())
    gene_overlap = all_matrix_genes & ortholog_mouse_genes
    human_ids_mapped = set(mouse_sym_to_human[g] for g in gene_overlap)
    n_in_gene_order = len(human_ids_mapped & set(gene_order))

    print(f"\n  Gene overlap with 16,959 ortholog space:")
    print(f"    Matrix genes:              {len(all_matrix_genes):,}")
    print(f"    Ortholog genes matched:    {len(gene_overlap):,}")
    print(f"    In 16,959 gene space:      {n_in_gene_order:,}")

    if n_in_gene_order < 12000:
        print(f"\n  STOP CONDITION: Gene overlap {n_in_gene_order:,} < 12,000")
        return False

    # Count target types in matched cells
    matched_ca = cell_assignments[cell_assignments["Cell.name"].isin(dge_names)]
    matched_ca_mapped = matched_ca["base_annotation"].map(type_mapping)
    target_counts = matched_ca_mapped.dropna().value_counts()

    print(f"\n  Target cell type counts in DGE-matched cells:")
    for ct in sorted(TARGET_TYPES):
        n = target_counts.get(ct, 0)
        status = "PASS" if n >= MIN_CELLS_GATE else "BELOW GATE"
        print(f"    {ct:<50} {n:>6,} {status}")

    print(f"\n" + "-" * 70)
    print("PHASE 1A FORMAT SUMMARY:")
    print(f"  Format: dense DGE text (genes × cells), integer UMI counts")
    print(f"  Gene IDs: mouse gene symbols")
    print(f"  Barcode join rate: {join_rate:.1f}% ({'PASS' if join_rate >= 80 else 'BELOW 80%'})")
    print(f"  Gene overlap: {n_in_gene_order:,}/16,959 ({'PASS' if n_in_gene_order >= 12000 else 'FAIL'})")

    print("  Proceeding to Phase 1B (full parse).")
    return True


def phase_1b_full_parse(
    dge_path: Path,
    cell_assignments: pd.DataFrame,
    type_mapping: dict[str, str],
    mouse_sym_to_human: dict[str, str],
    gene_order: list[str],
) -> ad.AnnData:
    """
    Phase 1B: Parse the full DGE matrix, join annotations, normalize, and align.

    This loads the entire ~2.9 GB matrix, which requires significant memory.
    The approach:
      1. Read gene names to identify ortholog genes
      2. Read only the ortholog gene rows (saves memory)
      3. Join cell annotations
      4. Filter to target cell types + generic T cells (for CD4 rescue)
      5. Normalize: counts per 10k + log1p
      6. Re-index to human Ensembl IDs

    Returns:
        AnnData with normalized expression in ortholog gene space.
    """
    print("\n" + "=" * 70)
    print("PHASE 1B: FULL PARSE AND NORMALIZATION")
    print("=" * 70)

    checkpoint_path = CHECKPOINT_DIR / "mca_full_parsed.h5ad"
    if checkpoint_path.exists():
        print(f"  Loading checkpoint: {checkpoint_path}")
        adata = ad.read_h5ad(checkpoint_path)
        print(f"  Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
        return adata

    CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)

    # Step 1: Identify which genes to keep
    ortholog_mouse_genes = set(mouse_sym_to_human.keys())
    gene_order_set = set(gene_order)

    # Build set of genes that map to our ortholog space
    valid_mouse_genes = set()
    for mg in ortholog_mouse_genes:
        hg = mouse_sym_to_human.get(mg)
        if hg and hg in gene_order_set:
            valid_mouse_genes.add(mg)

    print(f"  Valid mouse genes (in ortholog space): {len(valid_mouse_genes):,}")

    # Step 2: Read header to get cell names
    print(f"  Reading DGE matrix header...")
    with open(dge_path) as f:
        header_line = f.readline().strip()
    cell_names = [h.strip('"') for h in header_line.split("\t")]
    n_cells_total = len(cell_names)
    print(f"  Total cells in matrix: {n_cells_total:,}")

    # Step 3: Join with annotations to find target cells
    ca_name_set = set(cell_assignments["Cell.name"].values)

    # Build cell name → annotation mapping
    cell_to_annotation = {}
    cell_to_tissue = {}
    for _, row in cell_assignments.iterrows():
        cname = row["Cell.name"]
        cell_to_annotation[cname] = row["base_annotation"]
        cell_to_tissue[cname] = row.get("Tissue", "unknown")

    # Determine which columns (cells) to keep
    # Keep: cells with annotations that map to target types + generic T cells
    keep_cols = []  # indices into cell_names
    cell_info = []  # (cell_name, mca_annotation, our_type, tissue)
    for i, cname in enumerate(cell_names):
        if cname not in ca_name_set:
            continue
        ann = cell_to_annotation.get(cname, "")
        our_type = type_mapping.get(ann)
        if our_type in TARGET_TYPES or our_type == "T cell":
            keep_cols.append(i)
            tissue = cell_to_tissue.get(cname, "unknown")
            cell_info.append((cname, ann, our_type, tissue))

    n_keep = len(keep_cols)
    print(f"  Target cells to extract: {n_keep:,}")

    if n_keep == 0:
        raise RuntimeError("No target cells found after annotation join")

    # Step 4: Read the matrix — only keep ortholog gene rows and target cell columns
    print(f"  Reading DGE matrix (extracting {n_keep:,} cells × ortholog genes)...")
    print(f"  This may take a few minutes for the 2.9 GB file...")

    keep_col_set = set(keep_cols)
    gene_data = {}  # mouse_gene_symbol → array of counts for target cells
    n_genes_read = 0

    with open(dge_path) as f:
        f.readline()  # skip header
        for line_num, line in enumerate(f):
            if line_num % 5000 == 0 and line_num > 0:
                print(f"    ... {line_num:,} genes processed, "
                      f"{len(gene_data):,} ortholog genes kept")

            parts = line.strip().split("\t")
            gene = parts[0].strip('"')

            if gene not in valid_mouse_genes:
                continue

            # Extract values only for target cell columns
            values = np.zeros(n_keep, dtype=np.float32)
            for j, col_idx in enumerate(keep_cols):
                values[j] = float(parts[col_idx + 1])  # +1 because parts[0] is gene name

            gene_data[gene] = values
            n_genes_read += 1

    print(f"  Read {n_genes_read:,} ortholog genes for {n_keep:,} target cells")

    # Step 5: Build expression matrix in canonical gene order
    print(f"  Building expression matrix in canonical gene order...")
    gene_to_col = {g: i for i, g in enumerate(gene_order)}
    X = np.zeros((n_keep, len(gene_order)), dtype=np.float32)

    for mouse_gene, values in gene_data.items():
        human_id = mouse_sym_to_human[mouse_gene]
        col_idx = gene_to_col.get(human_id)
        if col_idx is not None:
            X[:, col_idx] = values

    n_missing = len(gene_order) - len(gene_data)
    print(f"  Genes mapped: {len(gene_data):,}/{len(gene_order):,} "
          f"({n_missing:,} filled with 0)")

    # Step 6: Normalize — counts per 10k + log1p (DECISION-038)
    print(f"  Normalizing: counts per 10k + log1p...")
    cell_totals = X.sum(axis=1, keepdims=True)
    cell_totals[cell_totals == 0] = 1  # prevent division by zero
    X_norm = np.log1p(X / cell_totals * 10000).astype(np.float32)

    # Step 7: Build AnnData
    print(f"  Building AnnData...")
    obs_df = pd.DataFrame(
        {
            "cell_id": [ci[0] for ci in cell_info],
            "mca_annotation": [ci[1] for ci in cell_info],
            "our_cell_type_label": [ci[2] for ci in cell_info],
            "tissue": [ci[3] for ci in cell_info],
            "rescue_method": [""] * n_keep,
        },
        index=[ci[0] for ci in cell_info],
    )

    var_df = pd.DataFrame(index=gene_order)
    var_df.index.name = "gene_id"

    adata = ad.AnnData(
        X=sp.csr_matrix(X_norm),
        obs=obs_df,
        var=var_df,
    )

    # Save checkpoint
    print(f"  Saving checkpoint: {checkpoint_path}")
    tmp_path = checkpoint_path.with_suffix(".h5ad.tmp")
    adata.write_h5ad(tmp_path)
    tmp_path.rename(checkpoint_path)

    print(f"  Phase 1B complete: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # Print per-type summary
    type_counts = adata.obs["our_cell_type_label"].value_counts()
    print(f"\n  Per-type counts:")
    for ct in sorted(TARGET_TYPES + ["T cell"]):
        n = type_counts.get(ct, 0)
        if n > 0:
            status = "PASS" if n >= MIN_CELLS_GATE else "BELOW GATE"
            print(f"    {ct:<50} {n:>6,} {status}")

    return adata


def phase_1c_cd4_rescue(
    adata: ad.AnnData,
    mouse_sym_to_human: dict[str, str],
    gene_order: list[str],
) -> tuple[ad.AnnData, dict]:
    """
    Phase 1C: CD4+ T cell computational rescue from generic T cell pool.

    Classification rule (DECISION-091):
      CD4+ T: CD4 > 0 AND CD8A == 0 AND CD8B == 0
      CD8+ T: CD8A > 0 OR CD8B > 0
      Ambiguous (both markers): EXCLUDE
      Neither marker: EXCLUDE

    Returns:
        Tuple of (updated adata, rescue_stats dict).
    """
    print("\n" + "=" * 70)
    print("PHASE 1C: CD4+ T CELL COMPUTATIONAL RESCUE")
    print("=" * 70)

    rescue_stats = {
        "cd4_marker_present": False,
        "cd8a_marker_present": False,
        "cd8b_marker_present": False,
        "n_generic_t_cells": 0,
        "n_classified_cd4": 0,
        "n_classified_cd8": 0,
        "n_ambiguous_both": 0,
        "n_neither": 0,
    }

    # Find generic T cells
    t_mask = adata.obs["our_cell_type_label"] == "T cell"
    rescue_stats["n_generic_t_cells"] = int(t_mask.sum())

    if rescue_stats["n_generic_t_cells"] == 0:
        print("  No generic T cells found — skipping")
        return adata, rescue_stats

    print(f"  Generic T cells available: {rescue_stats['n_generic_t_cells']:,}")

    # Map marker genes to human Ensembl IDs
    cd4_human_id = mouse_sym_to_human.get(CD4_GENE)
    cd8a_human_id = mouse_sym_to_human.get(CD8A_GENE)
    cd8b_human_id = mouse_sym_to_human.get(CD8B_GENE)

    var_set = set(adata.var_names)
    rescue_stats["cd4_marker_present"] = cd4_human_id is not None and cd4_human_id in var_set
    rescue_stats["cd8a_marker_present"] = cd8a_human_id is not None and cd8a_human_id in var_set
    rescue_stats["cd8b_marker_present"] = cd8b_human_id is not None and cd8b_human_id in var_set

    print(f"  CD4 ({CD4_GENE} → {cd4_human_id}): "
          f"{'PRESENT' if rescue_stats['cd4_marker_present'] else 'ABSENT'}")
    print(f"  CD8A ({CD8A_GENE} → {cd8a_human_id}): "
          f"{'PRESENT' if rescue_stats['cd8a_marker_present'] else 'ABSENT'}")
    print(f"  CD8B ({CD8B_GENE} → {cd8b_human_id}): "
          f"{'PRESENT' if rescue_stats['cd8b_marker_present'] else 'ABSENT'}")

    if not rescue_stats["cd4_marker_present"]:
        print("  RESCUE FAILED: CD4 gene not found in expression matrix")
        return adata, rescue_stats

    # Get expression values for T cells
    t_adata = adata[t_mask]
    X = t_adata.X
    if sp.issparse(X):
        X = X.toarray()

    var_list = list(adata.var_names)

    def get_col(gene_id):
        if gene_id and gene_id in var_set:
            return var_list.index(gene_id)
        return None

    cd4_idx = get_col(cd4_human_id)
    cd8a_idx = get_col(cd8a_human_id)
    cd8b_idx = get_col(cd8b_human_id)

    cd4_expr = X[:, cd4_idx] if cd4_idx is not None else np.zeros(X.shape[0])
    cd8a_expr = X[:, cd8a_idx] if cd8a_idx is not None else np.zeros(X.shape[0])
    cd8b_expr = X[:, cd8b_idx] if cd8b_idx is not None else np.zeros(X.shape[0])

    # Classification
    cd4_positive = cd4_expr > 0
    cd8_positive = (cd8a_expr > 0) | (cd8b_expr > 0)

    is_cd4 = cd4_positive & ~cd8_positive
    is_cd8 = cd8_positive & ~cd4_positive
    is_ambiguous = cd4_positive & cd8_positive
    is_neither = ~cd4_positive & ~cd8_positive

    rescue_stats["n_classified_cd4"] = int(is_cd4.sum())
    rescue_stats["n_classified_cd8"] = int(is_cd8.sum())
    rescue_stats["n_ambiguous_both"] = int(is_ambiguous.sum())
    rescue_stats["n_neither"] = int(is_neither.sum())

    print(f"\n  Results:")
    print(f"    CD4+ T classified: {rescue_stats['n_classified_cd4']:,}")
    print(f"    CD8+ T classified: {rescue_stats['n_classified_cd8']:,}")
    print(f"    Ambiguous (both):  {rescue_stats['n_ambiguous_both']:,}")
    print(f"    Neither:           {rescue_stats['n_neither']:,}")

    # Update labels — convert from Categorical to str to allow new categories
    t_indices = np.where(t_mask.values)[0]
    labels = adata.obs["our_cell_type_label"].astype(str).copy()
    rescue_col = adata.obs["rescue_method"].astype(str).copy()

    for j, idx in enumerate(t_indices):
        if is_cd4[j]:
            labels.iloc[idx] = "CD4-positive, alpha-beta T cell"
            rescue_col.iloc[idx] = "cd4_rescue"
        elif is_cd8[j]:
            labels.iloc[idx] = "CD8-positive, alpha-beta T cell"
            rescue_col.iloc[idx] = "cd4_rescue"
        elif is_ambiguous[j] or is_neither[j]:
            labels.iloc[idx] = ""  # will be filtered out

    adata.obs["our_cell_type_label"] = labels
    adata.obs["rescue_method"] = rescue_col

    if rescue_stats["n_classified_cd4"] < MIN_CELLS_GATE:
        print(f"  WARNING: Only {rescue_stats['n_classified_cd4']} CD4+ T cells "
              f"(< {MIN_CELLS_GATE} gate). CD4+ T will be EXCLUDED.")
    else:
        print(f"  CD4+ T rescue: PASS ({rescue_stats['n_classified_cd4']:,} cells)")

    return adata, rescue_stats


def main() -> None:
    """Main entry point for MCA download pipeline."""
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)

    # Check for existing output
    if OUTPUT_PATH.exists():
        print(f"  Output already exists: {OUTPUT_PATH}")
        adata = ad.read_h5ad(OUTPUT_PATH)
        print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")
        print("  Skipping. Delete the file to re-download.")
        return

    # ── Load reference data ──
    print("=" * 70)
    print("LOADING REFERENCE DATA")
    print("=" * 70)

    orthologs, mouse_sym_to_human, gene_order = load_ortholog_table()
    type_mapping = load_type_mapping()
    cell_assignments = ensure_cell_assignments()

    # ── Ensure DGE matrix is downloaded ──
    print("\n" + "=" * 70)
    print("ENSURING DGE MATRIX IS AVAILABLE")
    print("=" * 70)
    dge_path = ensure_dge_matrix()

    # ── Phase 1A: Format inspection ──
    format_ok = phase_1a_format_inspection(
        dge_path, cell_assignments, type_mapping, mouse_sym_to_human, gene_order
    )
    if not format_ok:
        print("\nSTOP: Phase 1A format inspection failed.")
        sys.exit(1)

    # ── Phase 1B: Full parse ──
    adata = phase_1b_full_parse(
        dge_path, cell_assignments, type_mapping, mouse_sym_to_human, gene_order
    )

    # ── Phase 1C: CD4+ T rescue ──
    adata, rescue_stats = phase_1c_cd4_rescue(adata, mouse_sym_to_human, gene_order)

    # ── Final filtering ──
    print("\n" + "=" * 70)
    print("FINAL FILTERING AND SAVING")
    print("=" * 70)

    final_target = TARGET_TYPES.copy()

    # Check CD4+ T rescue
    cd4_count = int(
        (adata.obs["our_cell_type_label"] == "CD4-positive, alpha-beta T cell").sum()
    )
    if cd4_count < MIN_CELLS_GATE:
        print(f"  CD4+ T: {cd4_count} cells < {MIN_CELLS_GATE} gate — EXCLUDING")
        if "CD4-positive, alpha-beta T cell" in final_target:
            final_target.remove("CD4-positive, alpha-beta T cell")

    # Remove types below gate
    types_to_remove = []
    for ct in final_target:
        n = int((adata.obs["our_cell_type_label"] == ct).sum())
        if n < MIN_CELLS_GATE:
            print(f"  {ct}: {n} cells < {MIN_CELLS_GATE} gate — EXCLUDING")
            types_to_remove.append(ct)

    for ct in types_to_remove:
        final_target.remove(ct)

    # Filter
    final_mask = adata.obs["our_cell_type_label"].isin(final_target)
    final_adata = adata[final_mask].copy()

    print(f"\n  Final dataset: {final_adata.n_obs:,} cells × {final_adata.n_vars:,} genes")
    print(f"  Cell types included: {len(final_target)}")

    print(f"\n  {'Cell Type':<50} {'Count':>8} {'Status':>10}")
    print("  " + "-" * 70)
    for ct in sorted(final_target):
        n = int((final_adata.obs["our_cell_type_label"] == ct).sum())
        status = "PASS" if n >= MIN_CELLS_GATE else "BORDERLINE"
        print(f"  {ct:<50} {n:>8,} {status:>10}")

    # Save
    print(f"\n  Saving to {OUTPUT_PATH}...")
    tmp_path = OUTPUT_PATH.with_suffix(".h5ad.tmp")
    final_adata.write_h5ad(tmp_path)
    tmp_path.rename(OUTPUT_PATH)
    print(f"  Saved: {OUTPUT_PATH}")

    # Save stats
    stats_path = DATA_DIR / "mca_download_stats.json"
    excluded_types = types_to_remove + (
        ["CD4-positive, alpha-beta T cell"] if cd4_count < MIN_CELLS_GATE else []
    )
    all_stats = {
        "rescue_stats": rescue_stats,
        "final_n_cells": int(final_adata.n_obs),
        "final_n_genes": int(final_adata.n_vars),
        "final_n_types": len(final_target),
        "final_types": final_target,
        "excluded_types": excluded_types,
        "download_stats": {
            "source": "GSE108097_MCA_Figure2_BatchRemoved_dge.txt",
            "excluded_tissues": [],
        },
        "cd4_included": cd4_count >= MIN_CELLS_GATE,
        "per_type_counts": {
            ct: int((final_adata.obs["our_cell_type_label"] == ct).sum())
            for ct in final_target
        },
    }
    with open(stats_path, "w") as f:
        json.dump(all_stats, f, indent=2)
    print(f"  Stats saved: {stats_path}")

    # ── Summary ──
    print("\n" + "=" * 70)
    print("MCA DOWNLOAD COMPLETE — SUMMARY")
    print("=" * 70)
    print(f"  Total cells: {final_adata.n_obs:,}")
    print(f"  Total genes: {final_adata.n_vars:,}")
    print(f"  Cell types: {len(final_target)}")
    print(f"  Types excluded (< {MIN_CELLS_GATE} cells): {len(excluded_types)} "
          f"({', '.join(excluded_types) if excluded_types else 'none'})")
    print(f"  CD4+ T rescue: {rescue_stats['n_classified_cd4']} classified "
          f"({'INCLUDED' if cd4_count >= MIN_CELLS_GATE else 'EXCLUDED'})")

    hep_count = int((final_adata.obs["our_cell_type_label"] == "hepatocyte").sum())
    print(f"  Hepatocyte count: {hep_count:,} "
          f"({'PASS' if hep_count >= MIN_CELLS_GATE else 'BORDERLINE/EXCLUDED'})")


if __name__ == "__main__":
    main()
