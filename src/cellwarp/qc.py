"""
CellWarp — Quality Control and Normalization Module

Filters low-quality cells, normalizes expression, identifies highly variable
genes, and computes dimensionality reductions (PCA, UMAP) for each species
independently.

Biology
-------
Single-cell RNA-seq data is noisy. Cells with very few detected genes are likely
empty droplets, damaged, or doublets. We remove these before normalization.

Normalization (counts-per-10k + log1p) corrects for differences in sequencing
depth between cells so that expression levels are comparable. Highly variable
genes (HVGs) capture biological variation while discarding genes that are
uniformly expressed (housekeeping) or dominated by technical noise.

Math
----
- Counts per 10k: x_norm = x / total_counts * 10_000 (each cell sums to 10k)
- log1p: x_log = log(1 + x_norm) (compresses dynamic range, stabilizes variance)
- HVGs: genes with highest variance-to-mean ratio after normalization
- PCA: linear dimensionality reduction to 50 components (captures >95% variance)
- UMAP: nonlinear 2D projection for visualization only (not used in Procrustes)
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for saving plots
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import seaborn as sns


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MIN_GENES_PER_CELL = 200     # Cells with fewer detected genes are removed
TARGET_SUM = 10_000          # Normalize to this count depth
N_TOP_GENES = 2_000          # Number of HVGs to select
N_PCS = 50                   # PCA components
UMAP_RANDOM_STATE = 42       # Reproducible UMAP
N_NEIGHBORS = 15             # For UMAP neighbor graph


# ---------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------


def filter_low_quality_cells(
    adata: ad.AnnData,
    min_genes: int = MIN_GENES_PER_CELL,
) -> ad.AnnData:
    """
    Remove cells with fewer than min_genes detected genes.

    Biology: Cells with very few detected genes are likely empty droplets,
    damaged cells, or multiplets. Removing them prevents noise from
    distorting downstream centroid computation.

    Args:
        adata: Raw AnnData (expects integer counts in .X).
        min_genes: Minimum number of genes with nonzero expression.

    Returns:
        Filtered AnnData (copy).
    """
    n_before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=min_genes)
    n_after = adata.n_obs
    n_removed = n_before - n_after
    print(f"  Filtered cells with <{min_genes} genes: "
          f"{n_before:,} → {n_after:,} ({n_removed:,} removed, "
          f"{n_removed / n_before * 100:.1f}%)")
    return adata


def compute_qc_metrics(adata: ad.AnnData) -> None:
    """
    Compute standard QC metrics and store them in adata.obs / adata.var.

    Adds: n_genes_by_counts, total_counts, n_cells_by_counts, mean_counts, etc.

    Args:
        adata: AnnData to annotate (modified in place).
    """
    sc.pp.calculate_qc_metrics(adata, inplace=True)


def normalize_and_log(
    adata: ad.AnnData,
    target_sum: float = TARGET_SUM,
) -> ad.AnnData:
    """
    Normalize to counts-per-10k and apply log1p transform.

    Math: For each cell, divide gene counts by total counts, multiply by
    target_sum, then take log(1 + x). This makes expression values
    comparable across cells with different sequencing depths and compresses
    the dynamic range for downstream linear methods (PCA).

    Biology: Raw UMI counts vary wildly between cells due to technical factors
    (capture efficiency, sequencing depth). Normalization removes this
    confound so that expression differences reflect biology, not library size.

    Args:
        adata: AnnData with raw counts. A copy of .X is saved to .raw first.
        target_sum: Target total counts per cell after normalization.

    Returns:
        Normalized AnnData (modified in place, also returned for chaining).
    """
    # Save raw counts for later use
    adata.raw = adata

    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)

    print(f"  Normalized to {target_sum:,} counts per cell + log1p")
    return adata


def select_hvgs(
    adata: ad.AnnData,
    n_top_genes: int = N_TOP_GENES,
) -> ad.AnnData:
    """
    Identify highly variable genes for dimensionality reduction.

    Biology: Most genes are either uniformly expressed (housekeeping) or
    dominated by technical noise. HVGs capture biological variation —
    the genes that actually differ between cell types.

    Math: Uses the Seurat method which computes per-gene mean and dispersion
    on log-normalized data, bins genes by mean, and selects genes with highest
    normalized dispersion within each bin.

    Args:
        adata: Normalized, log-transformed AnnData.
        n_top_genes: Number of HVGs to select.

    Returns:
        AnnData with .var["highly_variable"] column set (modified in place).
    """
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor="seurat",
    )

    n_hvg = adata.var["highly_variable"].sum()
    print(f"  Selected {n_hvg:,} highly variable genes (requested {n_top_genes:,})")
    return adata


def run_pca(
    adata: ad.AnnData,
    n_comps: int = N_PCS,
) -> ad.AnnData:
    """
    Run PCA on highly variable genes.

    Math: Principal component analysis finds the linear combinations of genes
    (components) that capture maximum variance. We use 50 components which
    typically explains >95% of variance, reducing from ~17k genes to 50
    dimensions while preserving most biological signal.

    Args:
        adata: AnnData with HVGs marked. PCA runs on HVG subset only.
        n_comps: Number of principal components to compute.

    Returns:
        AnnData with .obsm["X_pca"] and .uns["pca"] set.
    """
    sc.tl.pca(adata, n_comps=n_comps, mask_var="highly_variable")
    var_explained = adata.uns["pca"]["variance_ratio"].sum() * 100
    print(f"  PCA: {n_comps} components, {var_explained:.1f}% variance explained")
    return adata


def run_umap(
    adata: ad.AnnData,
    n_neighbors: int = N_NEIGHBORS,
    random_state: int = UMAP_RANDOM_STATE,
) -> ad.AnnData:
    """
    Compute neighbor graph and UMAP embedding for visualization.

    Math: UMAP is a nonlinear dimensionality reduction that preserves local
    neighborhood structure. It's used ONLY for visualization — our Procrustes
    analysis in Phase 2 uses PCA, not UMAP.

    Biology: If cell types form distinct clusters in UMAP space, it confirms
    the expression profiles are meaningfully different and our cell type
    annotations are consistent.

    Args:
        adata: AnnData with PCA computed.
        n_neighbors: Number of neighbors for the KNN graph.
        random_state: For reproducible UMAP layout.

    Returns:
        AnnData with .obsm["X_umap"] set.
    """
    sc.pp.neighbors(adata, n_pcs=N_PCS, n_neighbors=n_neighbors)
    sc.tl.umap(adata, random_state=random_state)
    print(f"  UMAP computed (n_neighbors={n_neighbors})")
    return adata


def plot_umap_by_cell_type(
    adata: ad.AnnData,
    species_label: str,
    output_path: Path,
    cell_type_column: str = "cell_type",
) -> None:
    """
    Generate and save a UMAP plot colored by cell type.

    Args:
        adata: AnnData with UMAP computed.
        species_label: "Human" or "Mouse" — used in plot title.
        output_path: Path to save the figure (PNG).
        cell_type_column: Column in .obs for coloring.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(
        adata,
        color=cell_type_column,
        ax=ax,
        show=False,
        title=f"{species_label} — UMAP colored by cell type",
        legend_loc="right margin",
        frameon=True,
    )
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"  Saved UMAP plot: {output_path}")


def print_qc_summary(
    adata: ad.AnnData,
    species_label: str,
    cell_type_column: str = "cell_type",
) -> dict:
    """
    Print per-cell-type QC statistics.

    Args:
        adata: AnnData after QC and normalization.
        species_label: "Human" or "Mouse".
        cell_type_column: Column in .obs for grouping.

    Returns:
        Dict with summary statistics.
    """
    print(f"\n  {species_label} QC Summary:")
    print(f"  {'Cell Type':<40} {'Cells':>8} {'Med Genes':>10} {'Med Counts':>12}")
    print(f"  {'-' * 72}")

    stats = {}
    for ct in sorted(adata.obs[cell_type_column].unique()):
        mask = adata.obs[cell_type_column] == ct
        ct_data = adata.obs.loc[mask]
        n_cells = mask.sum()
        med_genes = int(ct_data["n_genes_by_counts"].median())
        med_counts = int(ct_data["total_counts"].median())
        stats[ct] = {
            "n_cells": int(n_cells),
            "median_genes": med_genes,
            "median_counts": med_counts,
        }
        print(f"  {ct:<40} {n_cells:>8,} {med_genes:>10,} {med_counts:>12,}")

    print(f"  {'-' * 72}")
    total = adata.n_obs
    overall_med_genes = int(adata.obs["n_genes_by_counts"].median())
    overall_med_counts = int(adata.obs["total_counts"].median())
    print(f"  {'TOTAL':<40} {total:>8,} {overall_med_genes:>10,} {overall_med_counts:>12,}")

    return stats


def run_full_qc_pipeline(
    adata: ad.AnnData,
    species_label: str,
    output_dir: Path,
) -> ad.AnnData:
    """
    Run the complete QC pipeline on one species' AnnData.

    Steps: filter → QC metrics → normalize → HVGs → PCA → UMAP → plot.

    Args:
        adata: Raw AnnData (aligned to shared ortholog space).
        species_label: "Human" or "Mouse".
        output_dir: Directory for UMAP plots.

    Returns:
        Processed AnnData with normalization, HVGs, PCA, and UMAP.
    """
    print(f"  Input: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # Step 1: Filter low-quality cells
    adata = filter_low_quality_cells(adata)

    # Step 2: QC metrics (before normalization, on raw counts)
    compute_qc_metrics(adata)

    # Step 3: Normalize
    normalize_and_log(adata)

    # Step 4: HVGs
    select_hvgs(adata)

    # Step 5: PCA
    run_pca(adata)

    # Step 6: UMAP
    run_umap(adata)

    # Step 7: Plot
    plot_path = output_dir / f"umap_{species_label.lower()}.png"
    plot_umap_by_cell_type(adata, species_label, plot_path)

    # Step 8: Summary
    print_qc_summary(adata, species_label)

    return adata
