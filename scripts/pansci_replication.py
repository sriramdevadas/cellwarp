#!/usr/bin/env python3
"""
CellWarp — PanSci Replication: Rigidity Ranking Validation (T1-A)

Runs the Procrustes pipeline on PanSci (Cao lab, Science 2025) mouse data
against Tabula Sapiens human data. Primary goal: rigidity ranking validation
(Spearman ρ vs primary 35-type ranking). Secondary: geometry replication.

Biology
-------
PanSci profiles 21.8M nuclei from 13 mouse organs using EasySci snRNA-seq
(combinatorial indexing, NOT 10x). We use 6-month WT (wild-type) controls.
Median ~1,040 UMI/cell — lower than 10x but vastly more cells. Protocol
mismatch is the known risk (DECISION-097).

Pipeline
--------
  1. Load per-tissue MTX count matrices (Lung, Liver, Colon)
  2. Filter: WT genotype, 06_months age
  3. Map PanSci annotations to our 35-type ontology
  4. CD4/CD8 T cell split via Cd4/Cd8a marker genes
  5. QC: min 200 genes, max 20% mito (verify median genes — Task 0c)
  6. Normalize: CPM + log1p (identical to primary)
  7. Restrict to 16,959 ortholog gene space
  8. Compute centroids (tissue-matched endothelial from lung)
  9. Procrustes + permutation test (10,000 iterations)
  10. Rigidity ranking Spearman ρ vs primary
  11. Protocol sensitivity check (zero-fill vs residual)

Output
------
  output/validation/pansci_replication/

Hard constraints
----------------
  - No reference transfer from Tabula
  - No endothelial pooling across tissues (ISSUE-092)
  - ≥500 cell gate per type (PanSci has millions)
"""

from __future__ import annotations

import gc
import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.io as sio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from scipy.stats import spearmanr

from cellwarp.procrustes import (
    compute_centroids,
    pca_reduce_centroids,
    procrustes_align,
    permutation_test,
    compute_residual_vectors,
    RANDOM_SEED,
)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

OUTPUT_DIR = Path("output/validation/pansci_replication")
PANSCI_DIR = Path("data/replication/pansci")
ORTHOLOG_PATH = Path("data/phase1/orthologs_human_mouse.csv")
TABULA_CENTROIDS_PATH = Path("output/phase2/scaled_35types/centroids_human_35.csv")
PRIMARY_RESULTS_PATH = Path("output/phase2/scaled_35types/procrustes_results_35.json")
RESIDUALS_RANKED_PATH = Path("output/phase2/scaled_35types/residuals_ranked.csv")

N_PERMUTATIONS = 10_000
VARIANCE_THRESHOLD = 0.95
MIN_CELLS = 500  # Primary standard — PanSci has millions
ENDO_TISSUE = "lung"  # Per DECISION-101: lung endothelial best match

# Tissues to load
TISSUES_TO_LOAD = ["lung", "liver", "colon"]

# Organ suffixes for stripping from PanSci type names
ORGAN_SUFFIXES = [
    "-Kidney", "-Lung", "-Liver", "-Heart", "-Muscle", "-Stomach",
    "-BAT", "-iWAT", "-gWAT", "-Ileum", "-Colon", "-Jejunum", "-Duodenum",
]


def strip_organ_suffix(name: str) -> str:
    """Remove tissue suffix from PanSci cell type name."""
    for suffix in ORGAN_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def map_pansci_to_ontology(base_name: str) -> str | None:
    """Map a PanSci base cell type name to our 35-type ontology.

    Uses specificity-ordered keyword matching. Most specific patterns first
    to avoid substring collisions (e.g., 'goblet cells' matching 't cell').
    """
    n = base_name.lower()

    # Hepatocyte (unique, unambiguous)
    if "hepatocyte" in n:
        return "hepatocyte"
    if "hepatic stellate" in n or "cholangiocyte" in n:
        return None

    # Goblet — BEFORE T cell to prevent 'goblet cells' → 't cell' collision
    if "goblet" in n:
        return "large intestine goblet cell"

    # Enterocyte
    if "enterocyte" in n:
        return "enterocyte of epithelium of large intestine"

    # Immune — specific subtypes first
    if "plasma cell" in n or "plasma cells" in n:
        return "plasma cell"
    if "alveolar macrophage" in n or "interstitial macrophage" in n or "kupffer" in n:
        return "macrophage"
    if "dendritic cell" in n or "dendritic cells" in n:
        return "myeloid dendritic cell"
    if "neutrophil" in n or "basophil" in n or "eosinophil" in n:
        return "granulocyte"
    if "monocyte" in n:
        return "monocyte"
    if "mast cell" in n:
        return None

    # B cells
    if "_b cell" in n or "_b cells" in n or "cycling b" in n or "resting b" in n:
        return "B cell"
    if n.startswith("b cell") or n.startswith("b cells"):
        return "B cell"

    # T cells — require word boundary
    if "_t cell" in n or "_t cells" in n:
        return "T cell"
    if "lymphoid cell" in n and "b cell" not in n and "_b " not in n:
        return "T cell"

    # NK cells
    if "natural killer" in n or "nk cell" in n:
        return "natural killer cell"

    # Macrophage (generic)
    if "macrophage" in n:
        return "macrophage"

    # Myeloid (generic — after specific subtypes)
    if "myeloid cell" in n or "myeloid cells" in n:
        return "myeloid leukocyte"

    # Structural
    if "endothelial" in n and "lymphatic" not in n:
        return "endothelial cell"
    if "cardiac fibroblast" in n:
        return "fibroblast of cardiac tissue"
    if "fibroblast" in n or "fibro-adipogenic" in n or "fibro–adipogenic" in n:
        return "fibroblast"
    if "mural cell" in n or "pericyte" in n or "smooth muscle" in n:
        return "smooth muscle cell"

    # Epithelial — specific first
    if "urothelial" in n:
        return "bladder urothelial cell"
    if "basal cell" in n:
        return "basal cell"
    if "acinar cell" in n or "acinar cells" in n:
        return "pancreatic acinar cell"
    if "epithelial" in n:
        return "epithelial cell"

    # Tuft cells — DO NOT map (prevent 'tuft cells' → 't cell' collision)
    if "tuft" in n:
        return None

    return None


def load_tissue_data(tissue: str) -> sp.csr_matrix | None:
    """Load a single tissue's count matrix, metadata, and gene info.

    Returns sparse count matrix (cells × genes) or None if files missing.
    """
    mtx_path = PANSCI_DIR / f"{tissue}_genecount.mtx.gz"
    meta_path = PANSCI_DIR / f"{tissue}_df_cell.csv.gz"
    gene_path = PANSCI_DIR / f"{tissue}_df_gene.csv.gz"

    if not mtx_path.exists():
        print(f"  WARNING: {mtx_path} not found — skipping {tissue}")
        return None, None, None

    print(f"  Loading {tissue} MTX ({mtx_path.stat().st_size / 1e9:.1f} GB)...")
    # MTX files are gene × cell, so we transpose
    mtx = sio.mmread(gzip.open(mtx_path, "rb"))
    mtx = sp.csr_matrix(mtx.T)  # cells × genes
    print(f"    Matrix: {mtx.shape[0]:,} cells × {mtx.shape[1]:,} genes")

    with gzip.open(meta_path, "rt") as f:
        meta = pd.read_csv(f)
    print(f"    Metadata: {len(meta):,} rows")

    with gzip.open(gene_path, "rt") as f:
        genes = pd.read_csv(f)
    print(f"    Genes: {len(genes):,}")

    # Verify dimensions match
    assert mtx.shape[0] == len(meta), \
        f"Cell count mismatch: MTX has {mtx.shape[0]}, metadata has {len(meta)}"
    assert mtx.shape[1] == len(genes), \
        f"Gene count mismatch: MTX has {mtx.shape[1]}, genes has {len(genes)}"

    return mtx, meta, genes


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("PanSci Replication — Rigidity Ranking Validation (T1-A)")
    print("=" * 70)

    # ==================================================================
    # STEP 1: Load ortholog mapping and Tabula centroids
    # ==================================================================
    print("\nSTEP 1: Load reference data")
    print("=" * 70)

    ortho = pd.read_csv(ORTHOLOG_PATH)
    mouse_to_human = dict(zip(ortho["mouse_gene_name"], ortho["human_ensembl_id"]))
    print(f"  Ortholog map: {len(mouse_to_human):,} mouse → human pairs")

    tabula_centroids = pd.read_csv(TABULA_CENTROIDS_PATH, index_col=0)
    print(f"  Tabula centroids: {tabula_centroids.shape[0]} types × {tabula_centroids.shape[1]} genes")

    full_gene_set = list(tabula_centroids.columns)  # Human Ensembl IDs
    n_full = len(full_gene_set)

    # ==================================================================
    # STEP 2: Load and process each tissue
    # ==================================================================
    print(f"\nSTEP 2: Load and process tissues ({', '.join(TISSUES_TO_LOAD)})")
    print("=" * 70)

    # Collect per-type data across tissues
    type_expr_sums = {}  # {cell_type: cumulative sum vector (human gene space)}
    type_cell_counts = {}  # {cell_type: total cells}
    type_tissue_breakdown = {}  # {cell_type: {tissue: count}}
    endo_tissue_expr = {}  # {tissue: {sum, count}} for endothelial tissue matching
    all_per_cell_genes = []  # For median genes check (Task 0c verification)
    all_per_cell_zero_fill = {}  # {cell_type: list of per-cell zero-fill fractions}
    zero_fill_rates = {}  # per-tissue zero-fill in ortholog space

    # CD4/CD8 T cell tracking
    cd4_expr_sum = None
    cd8_expr_sum = None
    cd4_count = 0
    cd8_count = 0
    t_unsplit_sum = None
    t_unsplit_count = 0

    for tissue in TISSUES_TO_LOAD:
        print(f"\n--- Processing {tissue} ---")
        mtx, meta, genes = load_tissue_data(tissue)
        if mtx is None:
            continue

        gene_names = list(genes["gene_name"])

        # Filter to WT, 06_months
        mask_wt6m = (meta["genotype"] == "WT") & (meta["age_group"] == "06_months")
        n_wt6m = mask_wt6m.sum()
        print(f"  WT 6-month filter: {meta.shape[0]:,} → {n_wt6m:,} cells")

        # Apply filter to matrix
        indices = np.where(mask_wt6m.values)[0]
        mtx_filtered = mtx[indices]
        meta_filtered = meta.iloc[indices].copy()

        # Free full matrix
        del mtx
        gc.collect()

        # QC: min 200 genes, max 20% mito
        genes_per_cell = np.diff(mtx_filtered.indptr)
        # Actually compute genes detected (non-zero entries per row)
        genes_detected = np.array((mtx_filtered > 0).sum(axis=1)).flatten()
        total_counts = np.array(mtx_filtered.sum(axis=1)).flatten()

        # Mitochondrial genes
        mt_mask = np.array([g.startswith("mt-") for g in gene_names])
        if mt_mask.sum() > 0:
            mt_counts = np.array(mtx_filtered[:, mt_mask].sum(axis=1)).flatten()
            pct_mt = mt_counts / np.maximum(total_counts, 1) * 100
        else:
            pct_mt = np.zeros(mtx_filtered.shape[0])

        qc_pass = (genes_detected >= 200) & (pct_mt <= 20)
        n_qc = qc_pass.sum()
        print(f"  QC (≥200 genes, ≤20% mito): {n_wt6m:,} → {n_qc:,} cells "
              f"({n_wt6m - n_qc:,} removed)")
        print(f"  Median genes/cell: {np.median(genes_detected[qc_pass]):.0f}, "
              f"median UMI: {np.median(total_counts[qc_pass]):.0f}")

        all_per_cell_genes.extend(genes_detected[qc_pass].tolist())

        # Apply QC filter
        qc_indices = np.where(qc_pass)[0]
        mtx_qc = mtx_filtered[qc_indices]
        meta_qc = meta_filtered.iloc[qc_indices].copy()

        del mtx_filtered
        gc.collect()

        # Map cell types
        meta_qc["base_type"] = meta_qc["main_cell_type_organ"].apply(strip_organ_suffix)
        meta_qc["our_type"] = meta_qc["base_type"].apply(map_pansci_to_ontology)

        # CD4/CD8 T cell split
        t_mask = meta_qc["our_type"] == "T cell"
        n_t = t_mask.sum()
        if n_t > 0:
            cd4_idx_gene = gene_names.index("Cd4") if "Cd4" in gene_names else None
            cd8a_idx_gene = gene_names.index("Cd8a") if "Cd8a" in gene_names else None

            if cd4_idx_gene is not None and cd8a_idx_gene is not None:
                t_indices = np.where(t_mask.values)[0]
                t_matrix = mtx_qc[t_indices]

                cd4_detected = np.array(
                    (t_matrix[:, cd4_idx_gene] > 0).toarray()
                ).flatten()
                cd8_detected = np.array(
                    (t_matrix[:, cd8a_idx_gene] > 0).toarray()
                ).flatten()

                # CD4+/CD8- → CD4+ T; CD8+/CD4- → CD8+ T; rest → T cell (generic)
                is_cd4 = cd4_detected & ~cd8_detected
                is_cd8 = cd8_detected & ~cd4_detected

                n_cd4 = is_cd4.sum()
                n_cd8 = is_cd8.sum()
                n_generic = n_t - n_cd4 - n_cd8

                # Update type labels
                t_labels = meta_qc.loc[t_mask, "our_type"].copy()
                t_obs_idx = t_labels.index
                new_labels = pd.Series("T cell", index=t_obs_idx)
                new_labels.iloc[is_cd4] = "CD4-positive, alpha-beta T cell"
                new_labels.iloc[is_cd8] = "CD8-positive, alpha-beta T cell"
                meta_qc.loc[t_mask, "our_type"] = new_labels.values

                print(f"  T cell CD4/CD8 split: {n_cd4:,} CD4+, {n_cd8:,} CD8+, "
                      f"{n_generic:,} generic ({n_cd4/n_t:.1%} / {n_cd8/n_t:.1%} split rate)")

        # Report type counts for this tissue
        type_counts = meta_qc[meta_qc["our_type"].notna()]["our_type"].value_counts()
        print(f"  Mapped types in {tissue}:")
        for ct, n in type_counts.items():
            print(f"    {ct:<50} {n:>8,}")

        # Normalize: CPM + log1p
        print(f"  Normalizing (CPM + log1p)...")
        # Per-cell normalization
        row_sums = np.array(mtx_qc.sum(axis=1)).flatten()
        # Avoid division by zero
        row_sums[row_sums == 0] = 1
        # Scale to 10,000 counts
        scaling_factors = 1e4 / row_sums

        # Apply to sparse matrix efficiently
        mtx_norm = mtx_qc.multiply(scaling_factors[:, np.newaxis])
        # log1p — convert to CSR to support indexing
        mtx_norm = mtx_norm.log1p()
        if not sp.issparse(mtx_norm) or not isinstance(mtx_norm, sp.csr_matrix):
            mtx_norm = sp.csr_matrix(mtx_norm)

        del mtx_qc
        gc.collect()

        # Map to ortholog gene space
        # Handle duplicate gene names: take first occurrence only
        seen_genes = set()
        shared_mouse_genes = []
        shared_indices = []
        human_ids = []
        for i, g in enumerate(gene_names):
            if g in mouse_to_human and g not in seen_genes:
                shared_mouse_genes.append(g)
                shared_indices.append(i)
                human_ids.append(mouse_to_human[g])
                seen_genes.add(g)

        # Build fast lookup: human_id → full_gene_set index
        full_gene_idx = {g: i for i, g in enumerate(full_gene_set)}

        # Track zero-fill for this tissue
        present_human_ids = set(human_ids)
        missing_human_ids = [g for g in full_gene_set if g not in present_human_ids]
        tissue_zero_fill = len(missing_human_ids) / n_full
        zero_fill_rates[tissue] = tissue_zero_fill
        print(f"  Ortholog overlap: {len(shared_mouse_genes):,}/{len(mouse_to_human):,} "
              f"(zero-fill: {tissue_zero_fill:.1%})")

        # Precompute mapping from shared index → full gene set index
        shared_to_full = []
        for j_shared, hid in enumerate(human_ids):
            if hid in full_gene_idx:
                shared_to_full.append((j_shared, full_gene_idx[hid]))
            # else: gene present in PanSci orthologs but not in Tabula gene set

        # Subset normalized matrix to ortholog genes only (sparse, memory efficient)
        mtx_ortho = mtx_norm[:, shared_indices]  # cells × n_shared_orthologs

        # Build per-type centroids directly from sparse sums (no dense expansion)
        for ct in type_counts.index:
            ct_mask = meta_qc["our_type"] == ct
            ct_indices = np.where(ct_mask.values)[0]
            n_ct = len(ct_indices)

            if n_ct == 0:
                continue

            # Sum expression across cells for this type (sparse → 1D array)
            ct_sparse = mtx_ortho[ct_indices]
            ct_sum_shared = np.array(ct_sparse.sum(axis=0)).flatten()  # length = n_shared

            # Map summed expression to full gene space
            ct_sum_full = np.zeros(n_full, dtype=np.float64)
            for j_shared, j_full in shared_to_full:
                ct_sum_full[j_full] = ct_sum_shared[j_shared]

            # Handle endothelial tissue matching
            if ct == "endothelial cell":
                if tissue not in endo_tissue_expr:
                    endo_tissue_expr[tissue] = {"sum": np.zeros(n_full), "count": 0}
                endo_tissue_expr[tissue]["sum"] += ct_sum_full
                endo_tissue_expr[tissue]["count"] += n_ct

            # Accumulate global sums
            if ct not in type_expr_sums:
                type_expr_sums[ct] = np.zeros(n_full, dtype=np.float64)
                type_cell_counts[ct] = 0
                type_tissue_breakdown[ct] = {}
            type_expr_sums[ct] += ct_sum_full
            type_cell_counts[ct] += n_ct
            type_tissue_breakdown[ct][tissue] = n_ct

            # Per-type zero-fill: fraction of ortholog genes with zero expression
            # Sample up to 1000 cells for efficiency
            if ct not in all_per_cell_zero_fill:
                all_per_cell_zero_fill[ct] = []
            sample_n = min(n_ct, 1000)
            sample_idx = ct_indices[:sample_n]
            sample_sparse = mtx_ortho[sample_idx]
            # Map to full space: count how many of n_full genes are zero
            n_nonzero_shared = np.array((sample_sparse > 0).sum(axis=1)).flatten()
            # Genes in full space that are nonzero = genes in shared space that are nonzero
            # (missing genes are always zero)
            zero_frac = 1.0 - n_nonzero_shared / n_full
            all_per_cell_zero_fill[ct].extend(zero_frac.tolist())

            del ct_sparse
            gc.collect()

        del mtx_ortho

        del mtx_norm
        gc.collect()
        print(f"  {tissue} processing complete.")

    # ==================================================================
    # STEP 3: Compute centroids
    # ==================================================================
    print(f"\n{'=' * 70}")
    print("STEP 3: Compute centroids")
    print("=" * 70)

    # Task 0c verification: median genes per cell
    median_genes = np.median(all_per_cell_genes)
    print(f"\n  Task 0c VERIFICATION: median genes/cell = {median_genes:.0f}")
    if median_genes < 500:
        print(f"  *** FATAL: median genes {median_genes:.0f} < 500. Protocol mismatch. ***")
        print(f"  *** Logging but continuing — this should have been caught earlier. ***")
    elif median_genes < 1000:
        print(f"  WARNING: median genes {median_genes:.0f} in 500-1000 range (low)")
    else:
        print(f"  PASS: median genes {median_genes:.0f} ≥ 1000")

    # Cell count audit
    print(f"\n  Cell count audit (≥{MIN_CELLS} gate):")
    audit_rows = []
    usable_types = []
    for ct, n in sorted(type_cell_counts.items(), key=lambda x: -x[1]):
        if n >= MIN_CELLS:
            status = "PASS"
            usable_types.append(ct)
        elif n >= 200:
            status = "BORDERLINE"
        else:
            status = "FAIL"
        tissues = ", ".join(f"{t}({c:,})" for t, c in
                           sorted(type_tissue_breakdown[ct].items(), key=lambda x: -x[1]))
        audit_rows.append({"cell_type": ct, "n_cells": n, "status": status, "tissues": tissues})
        print(f"  {ct:<50} {n:>8,}  {status}")

    print(f"\n  Types passing ≥{MIN_CELLS} gate: {len(usable_types)}")

    # For endothelial: use lung only (DECISION-101, ISSUE-092)
    if "endothelial cell" in type_cell_counts:
        print(f"\n  Endothelial tissue selection (ISSUE-092/DECISION-101):")
        for tissue, data in endo_tissue_expr.items():
            print(f"    {tissue:<15} {data['count']:>8,} cells")
        if ENDO_TISSUE in endo_tissue_expr:
            print(f"  Using {ENDO_TISSUE} endothelial ({endo_tissue_expr[ENDO_TISSUE]['count']:,} cells)")
        else:
            print(f"  WARNING: {ENDO_TISSUE} not available, using pooled")

    # Compute centroids
    centroids = {}
    for ct in usable_types:
        if ct == "endothelial cell" and ENDO_TISSUE in endo_tissue_expr:
            # Use tissue-matched endothelial
            data = endo_tissue_expr[ENDO_TISSUE]
            centroids[ct] = data["sum"] / data["count"]
            print(f"  {ct}: centroid from {ENDO_TISSUE} ({data['count']:,} cells)")
        else:
            centroids[ct] = type_expr_sums[ct] / type_cell_counts[ct]

    centroid_df = pd.DataFrame(centroids, index=full_gene_set).T
    print(f"\n  PanSci centroid matrix: {centroid_df.shape[0]} types × {centroid_df.shape[1]} genes")

    # Find shared types with Tabula
    shared_types = sorted(set(centroid_df.index) & set(tabula_centroids.index))
    n_shared = len(shared_types)
    print(f"  Shared types with Tabula 35-type: {n_shared}")
    for ct in shared_types:
        print(f"    {ct}")

    if n_shared < 4:
        print("\n  *** STOP: Fewer than 4 shared types. Cannot run Procrustes. ***")
        return

    # Restrict centroids to shared types
    pansci_sub = centroid_df.loc[shared_types]
    tabula_sub = tabula_centroids.loc[shared_types]

    # Save centroids for V1 validation (added 2026-03-21)
    centroid_save_path = Path("data/centroids/pansci_16type_centroids.csv")
    centroid_save_path.parent.mkdir(parents=True, exist_ok=True)
    pansci_sub.to_csv(centroid_save_path)
    print(f"\n  Centroids saved to {centroid_save_path} ({pansci_sub.shape[0]} types × {pansci_sub.shape[1]} genes)")

    # ==================================================================
    # STEP 4: Procrustes
    # ==================================================================
    print(f"\n{'=' * 70}")
    print(f"STEP 4: Procrustes — PanSci ({n_shared} types) → Tabula human")
    print("=" * 70)

    print("\n--- PCA on combined centroids ---")
    human_pca, mouse_pca, pca_model, types_list = pca_reduce_centroids(
        tabula_sub, pansci_sub, variance_threshold=VARIANCE_THRESHOLD
    )

    print("\n--- Procrustes alignment ---")
    result = procrustes_align(human_pca, mouse_pca)

    print(f"\n--- Permutation test ({N_PERMUTATIONS:,} iterations) ---")
    p_val, null_dist = permutation_test(human_pca, mouse_pca, N_PERMUTATIONS, RANDOM_SEED)

    obs_null = result.distance / np.median(null_dist)

    print("\n--- Per-type residuals ---")
    residuals = compute_residual_vectors(result, types_list)
    residual_mags = {ct: float(np.linalg.norm(residuals[ct])) for ct in types_list}

    # Sort by magnitude
    sorted_types = sorted(residual_mags, key=residual_mags.get, reverse=True)
    total_ssr = sum(v ** 2 for v in residual_mags.values())
    print(f"\n  Residual ranking (n={n_shared}):")
    for i, ct in enumerate(sorted_types, 1):
        pct = residual_mags[ct] ** 2 / total_ssr * 100
        print(f"    {i:>2}. {ct:<50} {residual_mags[ct]:>8.3f}  ({pct:.1f}% SSR)")

    # Zero-fill rate
    # Use the gene overlap from the first tissue (all tissues use same gene file)
    overall_zero_fill = np.mean(list(zero_fill_rates.values()))
    print(f"\n  Zero-fill rate: {overall_zero_fill:.1%}")

    # ==================================================================
    # STEP 5: Rigidity ranking correlation (PRIMARY GOAL)
    # ==================================================================
    print(f"\n{'=' * 70}")
    print("STEP 5: Rigidity Ranking Correlation (PRIMARY GOAL)")
    print("=" * 70)

    primary_residuals_df = pd.read_csv(RESIDUALS_RANKED_PATH)
    primary_mag_dict = dict(
        zip(primary_residuals_df["cell_type"], primary_residuals_df["residual_magnitude"])
    )

    matched_types = sorted(set(residual_mags.keys()) & set(primary_mag_dict.keys()))
    n_matched = len(matched_types)

    if n_matched >= 4:
        pansci_mags = [residual_mags[ct] for ct in matched_types]
        primary_mags = [primary_mag_dict[ct] for ct in matched_types]
        rho, rho_p = spearmanr(pansci_mags, primary_mags)
    else:
        rho, rho_p = float("nan"), float("nan")

    print(f"\n  Spearman ρ = {rho:.3f}, p = {rho_p:.4f} (n={n_matched})")

    # Per-type rank comparison
    if n_matched >= 4:
        pansci_rank = pd.Series(pansci_mags, index=matched_types).rank(ascending=False)
        primary_rank = pd.Series(primary_mags, index=matched_types).rank(ascending=False)
        print(f"\n  {'Cell type':<50} {'PanSci rank':>11} {'Primary rank':>12} {'Δ rank':>8}")
        print(f"  {'-' * 81}")
        rank_discrepancies = []
        for ct in matched_types:
            d = abs(int(pansci_rank[ct]) - int(primary_rank[ct]))
            rank_discrepancies.append((ct, d))
            print(f"  {ct:<50} {int(pansci_rank[ct]):>11} {int(primary_rank[ct]):>12} "
                  f"{d:>8}")
        # Largest discrepancies
        rank_discrepancies.sort(key=lambda x: -x[1])
        print(f"\n  Largest rank discrepancies:")
        for ct, d in rank_discrepancies[:5]:
            print(f"    {ct}: Δ={d}")

    # Verdict
    print(f"\n  VERDICT (pre-specified logic):")
    if n_matched >= 15 and rho >= 0.50 and rho_p < 0.05:
        ranking_verdict = "PASS"
        detail = (f"n={n_matched} ≥ 15, ρ={rho:.3f} ≥ 0.50, p={rho_p:.4f} < 0.05. "
                  "Ranking replicates. Flag to advisor.")
    elif n_matched >= 15 and (rho < 0.50 or rho_p >= 0.05):
        ranking_verdict = "FAIL"
        detail = (f"n={n_matched} ≥ 15, but ρ={rho:.3f} {'< 0.50' if rho < 0.50 else '≥ 0.50'}, "
                  f"p={rho_p:.4f}. Third consecutive ranking null. Flag to advisor.")
    else:
        ranking_verdict = "UNDERPOWERED"
        detail = (f"n={n_matched} < 15 after annotation. Do not interpret ρ. "
                  "Flag to advisor.")

    print(f"  {ranking_verdict}: {detail}")

    # ==================================================================
    # STEP 6: Protocol sensitivity check
    # ==================================================================
    print(f"\n{'=' * 70}")
    print("STEP 6: Protocol Sensitivity Check")
    print("=" * 70)

    # Compute per-type mean zero-fill fraction
    type_zero_fill = {}
    for ct in matched_types:
        if ct in all_per_cell_zero_fill and len(all_per_cell_zero_fill[ct]) > 0:
            type_zero_fill[ct] = np.mean(all_per_cell_zero_fill[ct])

    if len(type_zero_fill) >= 4:
        zf_types = sorted(type_zero_fill.keys())
        zf_vals = [type_zero_fill[ct] for ct in zf_types]
        res_vals = [residual_mags[ct] for ct in zf_types]
        rho_zf, p_zf = spearmanr(zf_vals, res_vals)
        print(f"\n  Per-type zero-fill vs Procrustes residual:")
        print(f"  Spearman ρ = {rho_zf:.3f}, p = {p_zf:.4f} (n={len(zf_types)})")
        if abs(rho_zf) > 0.40:
            print(f"  WARNING: |ρ| > 0.40 — low-UMI dropout is inflating residuals "
                  f"for specific types. Ranking may be confounded by protocol.")
        else:
            print(f"  OK: |ρ| ≤ 0.40 — dropout does not systematically inflate "
                  f"residuals for specific types.")
        print(f"\n  Per-type zero-fill fractions:")
        for ct in sorted(type_zero_fill, key=type_zero_fill.get):
            print(f"    {ct:<50} {type_zero_fill[ct]:.3f}")
    else:
        rho_zf, p_zf = float("nan"), float("nan")
        print(f"  Insufficient types for zero-fill correlation")

    # ==================================================================
    # STEP 7: Save outputs
    # ==================================================================
    print(f"\n{'=' * 70}")
    print("Saving outputs")
    print("=" * 70)

    results = {
        "diagnostic": "PanSci replication — rigidity ranking validation (T1-A)",
        "date": "2026-03-15",
        "dataset": {
            "name": "PanSci (Cao lab, Science 2025)",
            "accession": "GEO GSE247719",
            "protocol": "EasySci snRNA-seq (combinatorial indexing)",
            "condition": "WT, 06_months (young adult)",
            "tissues_loaded": TISSUES_TO_LOAD,
            "n_tissues": len(TISSUES_TO_LOAD),
            "median_genes_per_cell": float(median_genes),
            "zero_fill_rate": float(overall_zero_fill),
            "endothelial_tissue": ENDO_TISSUE,
        },
        "procrustes": {
            "n_types": n_shared,
            "cell_types": shared_types,
            "p_value": float(p_val),
            "distance": float(result.distance),
            "obs_null_ratio": float(obs_null),
            "scaling": float(result.scaling),
            "null_median": float(np.median(null_dist)),
            "pca_components": int(pca_model.n_components_),
            "per_type_residuals": {
                ct: {"magnitude": residual_mags[ct]} for ct in types_list
            },
        },
        "rigidity_ranking": {
            "rho": float(rho) if not np.isnan(rho) else None,
            "p_value": float(rho_p) if not np.isnan(rho_p) else None,
            "n_matched_types": n_matched,
            "matched_types": matched_types,
            "verdict": ranking_verdict,
            "detail": detail,
        },
        "protocol_sensitivity": {
            "zerofill_vs_residual_rho": float(rho_zf) if not np.isnan(rho_zf) else None,
            "zerofill_vs_residual_p": float(p_zf) if not np.isnan(p_zf) else None,
            "per_type_zero_fill": {ct: float(v) for ct, v in type_zero_fill.items()}
            if type_zero_fill
            else None,
        },
        "cell_count_audit": audit_rows,
        "random_seed": RANDOM_SEED,
    }

    with open(OUTPUT_DIR / "pansci_replication.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"  Results: {OUTPUT_DIR / 'pansci_replication.json'}")

    np.save(OUTPUT_DIR / "null_distribution.npy", null_dist)

    # Save ranking comparison CSV
    if n_matched >= 4:
        rank_df = pd.DataFrame({
            "cell_type": matched_types,
            "pansci_residual": [residual_mags[ct] for ct in matched_types],
            "primary_residual": [primary_mag_dict[ct] for ct in matched_types],
        })
        rank_df.to_csv(OUTPUT_DIR / "ranking_comparison.csv", index=False)

    # --- Plots ---
    # 1. Null distribution
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(null_dist, bins=50, alpha=0.7, color="darkcyan", edgecolor="white")
    ax.axvline(result.distance, color="red", linewidth=2,
               label=f"Observed (d={result.distance:.2f})")
    ax.set_title(f"PanSci ({n_shared} types) vs Tabula human\n"
                 f"p={p_val:.4f}, obs/null={obs_null:.3f}, ρ={rho:.3f}")
    ax.set_xlabel("Procrustes distance")
    ax.set_ylabel("Count")
    ax.legend()
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "null_distribution.png", dpi=150, bbox_inches="tight")
    plt.close()

    # 2. Rigidity scatter
    if n_matched >= 4:
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.scatter(
            [primary_mag_dict[ct] for ct in matched_types],
            [residual_mags[ct] for ct in matched_types],
            s=60, c="darkcyan", edgecolors="teal", linewidths=0.5, zorder=3,
        )
        for ct in matched_types:
            short = ct[:22] + "..." if len(ct) > 22 else ct
            ax.annotate(short, (primary_mag_dict[ct], residual_mags[ct]),
                        fontsize=6, ha="left", va="bottom",
                        xytext=(4, 4), textcoords="offset points")
        lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
                max(ax.get_xlim()[1], ax.get_ylim()[1])]
        ax.plot(lims, lims, "k--", alpha=0.3, zorder=1)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_xlabel("Primary 35-type residual magnitude")
        ax.set_ylabel("PanSci residual magnitude")
        ax.set_title(f"Rigidity: PanSci vs Primary\nρ={rho:.3f}, p={rho_p:.4f}, n={n_matched}")
        plt.tight_layout()
        fig.savefig(OUTPUT_DIR / "rigidity_scatter.png", dpi=150, bbox_inches="tight")
        plt.close()

    # 3. Protocol sensitivity
    if len(type_zero_fill) >= 4:
        fig, ax = plt.subplots(figsize=(8, 6))
        zf_ct = sorted(type_zero_fill.keys())
        ax.scatter(
            [type_zero_fill[ct] for ct in zf_ct],
            [residual_mags[ct] for ct in zf_ct],
            s=60, c="coral", edgecolors="darkred", linewidths=0.5, zorder=3,
        )
        for ct in zf_ct:
            short = ct[:20] + "..." if len(ct) > 20 else ct
            ax.annotate(short, (type_zero_fill[ct], residual_mags[ct]),
                        fontsize=6, ha="left", va="bottom",
                        xytext=(4, 4), textcoords="offset points")
        ax.set_xlabel("Mean per-cell zero-fill fraction (ortholog space)")
        ax.set_ylabel("Procrustes residual magnitude")
        ax.set_title(f"Protocol Sensitivity: Zero-fill vs Residual\n"
                     f"ρ={rho_zf:.3f}, p={p_zf:.4f}")
        plt.tight_layout()
        fig.savefig(OUTPUT_DIR / "protocol_sensitivity.png", dpi=150, bbox_inches="tight")
        plt.close()

    print(f"  Plots: null_distribution.png, rigidity_scatter.png, protocol_sensitivity.png")

    # ==================================================================
    # FINAL SUMMARY
    # ==================================================================
    print(f"\n{'=' * 70}")
    print("FINAL SUMMARY")
    print("=" * 70)

    print(f"\n  Dataset: PanSci (EasySci snRNA-seq, {', '.join(TISSUES_TO_LOAD)})")
    print(f"  Filter: WT, 06_months")
    print(f"  Median genes/cell: {median_genes:.0f}")
    print(f"  Zero-fill rate: {overall_zero_fill:.1%}")
    print(f"  Endothelial: {ENDO_TISSUE} only (ISSUE-092)")
    print(f"\n  Procrustes: n={n_shared}, p={p_val:.4f}, obs/null={obs_null:.3f}, "
          f"scaling={result.scaling:.3f}")
    print(f"  Rigidity ranking: ρ={rho:.3f}, p={rho_p:.4f} (n={n_matched})")
    print(f"  Protocol sensitivity: zero-fill vs residual ρ={rho_zf:.3f}")
    print(f"\n  VERDICT: {ranking_verdict}")
    print(f"  {detail}")


if __name__ == "__main__":
    main()
