#!/usr/bin/env python3
"""
CellWarp — T1-A Expanded Replication: Sun et al. 2023 (BORDERLINE types added)

Extends the Sun2023 replication from 6 types (scripts/16) to include the 8 BORDERLINE
types from DECISION-098 pre-download gate. Uses IDENTICAL annotation logic (Leiden
clustering + canonical marker gene scoring) on the SAME raw Cell Ranger data — no
new downloads. The goal is to increase n to ≥10 for a powered Spearman ρ test.

Biology
-------
The original run (scripts/16) found geometry replicates at n=6 (p=0.0099, obs/null=0.724)
but rigidity ranking ρ=−0.314 is uninterpretable at n=6 (need |ρ|≥0.707 for p<0.05
with n=6). By adding BORDERLINE types (adventitial, HSC, HPC, MSC, myeloid leukocyte,
NKT, stromal, generic T cell), we can reach n≥10 where |ρ|≥0.564 suffices for p<0.05.

Two-pass annotation preserves original type assignments: BORDERLINE markers only score
clusters that were "unassigned" in the first pass.

Steps
-----
  1. Load Cell Ranger outputs for 8 YC scRNA-seq tissues (same as scripts/16)
  2. QC: min 200 genes, max 20% mitochondrial (same as scripts/16)
  3. Two-pass annotation:
     a. First pass: original 14 marker types (identical to scripts/16)
     b. Second pass: 8 BORDERLINE types on "unassigned" clusters only
  4. Cell count audit with ≥200 gate (DECISION-090)
  5. Normalize: CPM + log1p, restrict to 16,959 ortholog space
  6. Compute Sun2023 expanded centroids
  7. Load Tabula human 35-type centroids (from primary analysis)
  8. Procrustes + permutation test (10,000 iterations)
  9. Spearman ρ vs primary 35-type rigidity ranking
  10. T1-B negative control (obs/null comparison)
  11. Verdict per decision tree

Output
------
  output/validation/sun2023_replication_expanded/

Hard constraints
----------------
  - Does NOT modify scripts/16 or output/validation/sun2023_replication/
  - Annotation is marker-gene-based (no reference transfer from Tabula)
  - Cell gate ≥200 per DECISION-090
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import anndata as ad
import scanpy as sc
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

OUTPUT_DIR = Path("output/validation/sun2023_replication_expanded")
EXTRACTED_DIR = Path("data/replication/sun2023/extracted")
ORTHOLOG_PATH = Path("data/phase1/orthologs_human_mouse.csv")
TABULA_HUMAN_CENTROIDS_PATH = Path(
    "output/phase2/scaled_35types/centroids_human_35.csv"
)
PRIMARY_RESULTS_PATH = Path(
    "output/phase2/scaled_35types/procrustes_results_35.json"
)
RESIDUALS_RANKED_PATH = Path(
    "output/phase2/scaled_35types/residuals_ranked.csv"
)

N_PERMUTATIONS = 10_000
VARIANCE_THRESHOLD = 0.95
MIN_CELLS = 200  # DECISION-090: ≥200 for replication datasets

# Tissues to load (YC scRNA-seq only — identical to scripts/16)
TISSUES = {
    "YC-Liver": "liver",
    "YC-Spleen": "spleen",
    "YC-BM": "bone_marrow",
    "YC-Lung": "lung",
    "YC-Kidney": "kidney",
    "YC-Blood": "blood",
    "YC-Artery": "aorta",
    "YC-Intestine": "intestine",
}

# ---------------------------------------------------------------------------
# PASS 1 markers — identical to scripts/16_sun2023_replication.py
# ---------------------------------------------------------------------------
MARKERS_PASS1 = {
    "hepatocyte": {"pos": ["Alb", "Ttr", "Apoa1"], "neg": []},
    "B cell": {"pos": ["Cd79a", "Cd79b", "Ms4a1"], "neg": ["Cd3e"]},
    "T cell": {"pos": ["Cd3e", "Cd3d"], "neg": ["Cd79a"]},
    "endothelial cell": {"pos": ["Pecam1", "Cdh5"], "neg": []},
    "macrophage": {"pos": ["Adgre1", "Csf1r", "C1qa"], "neg": ["Cd3e", "Cd79a"]},
    "monocyte": {"pos": ["Ly6c2", "Csf1r"], "neg": ["Cd3e", "Cd79a", "Adgre1"]},
    "natural killer cell": {
        "pos": ["Nkg7", "Klrb1c", "Gzma"],
        "neg": ["Cd3e"],
    },
    "neutrophil": {
        "pos": ["S100a8", "S100a9", "Ly6g"],
        "neg": ["Cd3e", "Cd79a"],
    },
    "epithelial cell": {
        "pos": ["Epcam", "Krt8", "Krt18"],
        "neg": ["Pecam1"],
    },
    "fibroblast": {
        "pos": ["Col1a1", "Col1a2", "Dcn"],
        "neg": ["Pecam1", "Epcam"],
    },
    "smooth muscle cell": {
        "pos": ["Acta2", "Myh11", "Tagln"],
        "neg": ["Pecam1"],
    },
    "plasma cell": {"pos": ["Sdc1", "Xbp1", "Jchain"], "neg": ["Ms4a1"]},
    "myeloid dendritic cell": {
        "pos": ["Flt3", "Itgax", "H2-Aa"],
        "neg": ["Cd3e", "Cd79a", "Adgre1"],
    },
    "granulocyte": {
        "pos": ["S100a8", "S100a9"],
        "neg": ["Cd3e", "Cd79a", "Csf1r"],
    },
}

# ---------------------------------------------------------------------------
# PASS 2 markers — 8 BORDERLINE types from DECISION-098
# Canonical markers from immunology/stem cell biology textbooks.
# Only applied to clusters that were "unassigned" in pass 1.
# ---------------------------------------------------------------------------
MARKERS_PASS2 = {
    "adventitial cell": {
        "pos": ["Pi16", "Pdgfra", "Ly6a"],
        "neg": ["Acta2", "Pecam1"],
    },
    "hematopoietic precursor cell": {
        "pos": ["Kit", "Cd34", "Flt3"],
        "neg": ["Cd3e", "Cd79a", "Adgre1", "Ly6g"],
    },
    "hematopoietic stem cell": {
        "pos": ["Kit", "Ly6a", "Hlf"],
        "neg": ["Cd3e", "Cd79a", "Cd34"],
    },
    "mature NK T cell": {
        "pos": ["Cd3e", "Nkg7", "Klrb1c"],
        "neg": ["Cd79a"],
    },
    "mesenchymal stem cell": {
        "pos": ["Lepr", "Cxcl12", "Pdgfra"],
        "neg": ["Cd3e", "Cd79a", "Pecam1"],
    },
    "myeloid leukocyte": {
        "pos": ["Lyz2", "Cd68", "Csf1r"],
        "neg": ["Cd3e", "Cd79a"],
    },
    "stromal cell": {
        "pos": ["Vim", "Col1a1", "Pdgfra"],
        "neg": ["Pecam1", "Epcam", "Cd3e"],
    },
}
# "T cell" (generic) is already captured in pass 1 — cells that are T cell
# but neither CD4+ nor CD8+ remain as "T cell" after the refinement step.


def load_tissue(tissue_dir: str, tissue_label: str) -> ad.AnnData:
    """Load a 10x Cell Ranger matrix and add tissue metadata."""
    path = EXTRACTED_DIR / tissue_dir
    adata = sc.read_10x_mtx(path, var_names="gene_symbols", cache=False)
    adata.obs["tissue"] = tissue_label
    adata.obs["source_dir"] = tissue_dir
    adata.var_names_make_unique()
    return adata


def qc_filter(adata: ad.AnnData) -> ad.AnnData:
    """Apply standard QC: min genes, max mitochondrial fraction."""
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    n_before = adata.n_obs
    adata = adata[adata.obs["n_genes_by_counts"] >= 200].copy()
    adata = adata[adata.obs["pct_counts_mt"] <= 20].copy()
    n_after = adata.n_obs
    print(f"  QC: {n_before:,} → {n_after:,} cells ({n_before - n_after:,} removed)")
    return adata


def score_markers(adata: ad.AnnData, markers: dict) -> pd.DataFrame:
    """Score each cell for each cell type using marker genes.

    Score = mean(positive markers present) - 0.5 * mean(negative markers present).
    Uses raw counts (>0 detection) for robustness with sparse data.
    Identical logic to scripts/16.
    """
    X = adata.X
    if sp.issparse(X):
        X = X.toarray()
    gene_names = list(adata.var_names)
    gene_idx = {g: i for i, g in enumerate(gene_names)}

    scores = {}
    for ct, mrk in markers.items():
        pos_genes = [g for g in mrk["pos"] if g in gene_idx]
        neg_genes = [g for g in mrk["neg"] if g in gene_idx]

        if not pos_genes:
            scores[ct] = np.zeros(adata.n_obs)
            continue

        pos_detected = np.zeros(adata.n_obs)
        for g in pos_genes:
            pos_detected += (X[:, gene_idx[g]] > 0).astype(float)
        pos_score = pos_detected / len(pos_genes)

        neg_detected = np.zeros(adata.n_obs)
        for g in neg_genes:
            neg_detected += (X[:, gene_idx[g]] > 0).astype(float)
        neg_score = neg_detected / max(len(neg_genes), 1)

        scores[ct] = pos_score - 0.5 * neg_score

    return pd.DataFrame(scores, index=adata.obs_names)


def annotate_cells(adata: ad.AnnData) -> ad.AnnData:
    """Two-pass annotation using marker gene scoring.

    Pass 1: Score clusters with original markers (identical to scripts/16).
    Pass 2: For clusters labeled 'unassigned', rescore with BORDERLINE markers.
    """
    print("  Annotating cell types via two-pass marker gene scoring...")

    # --- Clustering (identical to scripts/16) ---
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)
    sc.pp.highly_variable_genes(adata_norm, n_top_genes=2000, flavor="seurat")
    sc.pp.pca(adata_norm, n_comps=30, use_highly_variable=True)
    sc.pp.neighbors(adata_norm, n_neighbors=15)
    sc.tl.leiden(adata_norm, resolution=1.0, random_state=RANDOM_SEED)

    adata.obs["leiden"] = adata_norm.obs["leiden"].values
    cluster_labels = adata.obs["leiden"]
    clusters = sorted(cluster_labels.unique(), key=int)

    # --- Pass 1: Original markers ---
    marker_scores_p1 = score_markers(adata, MARKERS_PASS1)
    cluster_assignments = {}
    for cl in clusters:
        mask = cluster_labels == cl
        mean_scores = marker_scores_p1.loc[mask].mean()
        best_type = mean_scores.idxmax()
        best_score = mean_scores[best_type]
        if best_score >= 0.3:
            cluster_assignments[cl] = best_type
        else:
            cluster_assignments[cl] = "unassigned"

    adata.obs["cell_type"] = cluster_labels.map(cluster_assignments).values
    adata.obs["cell_type"] = adata.obs["cell_type"].astype(str)

    n_unassigned_clusters = sum(1 for v in cluster_assignments.values() if v == "unassigned")
    n_unassigned_cells = (adata.obs["cell_type"] == "unassigned").sum()
    print(f"  Pass 1: {len(clusters)} clusters, {n_unassigned_clusters} unassigned "
          f"({n_unassigned_cells:,} cells)")

    # --- Pass 2: BORDERLINE markers on unassigned clusters ---
    if n_unassigned_clusters > 0:
        marker_scores_p2 = score_markers(adata, MARKERS_PASS2)
        for cl in clusters:
            if cluster_assignments[cl] != "unassigned":
                continue
            mask = cluster_labels == cl
            mean_scores = marker_scores_p2.loc[mask].mean()
            best_type = mean_scores.idxmax()
            best_score = mean_scores[best_type]
            if best_score >= 0.3:
                cluster_assignments[cl] = best_type
                print(f"    Cluster {cl}: {mask.sum():,} cells → {best_type} "
                      f"(score={best_score:.3f})")

        adata.obs["cell_type"] = cluster_labels.map(cluster_assignments).values
        adata.obs["cell_type"] = adata.obs["cell_type"].astype(str)

    n_still_unassigned = (adata.obs["cell_type"] == "unassigned").sum()
    print(f"  Pass 2: {n_still_unassigned:,} cells still unassigned")

    # --- Tissue-aware hepatocyte rescue (identical to scripts/16) ---
    gene_names = list(adata.var_names)
    gene_idx = {g: i for i, g in enumerate(gene_names)}
    X_all = adata.X
    if sp.issparse(X_all):
        X_dense = X_all.toarray()
    else:
        X_dense = X_all

    alb_idx = gene_idx.get("Alb")
    if alb_idx is not None and "tissue" in adata.obs.columns:
        liver_mask = adata.obs["tissue"] == "liver"
        alb_expr = X_dense[:, alb_idx] > 0
        hep_mask = liver_mask.values & alb_expr
        n_hep = hep_mask.sum()
        if n_hep > 0:
            adata.obs.loc[hep_mask, "cell_type"] = "hepatocyte"
            print(f"  Hepatocyte rescue: {n_hep} liver cells with Alb > 0 → hepatocyte")

    # --- CD4/CD8 T cell refinement (identical to scripts/16) ---
    t_mask = adata.obs["cell_type"] == "T cell"
    if t_mask.sum() > 0:
        cd4_idx = gene_idx.get("Cd4")
        cd8a_idx = gene_idx.get("Cd8a")
        if cd4_idx is not None and cd8a_idx is not None:
            t_indices = np.where(t_mask.values)[0]
            cd4_expr = X_dense[t_indices, cd4_idx] > 0
            cd8_expr = X_dense[t_indices, cd8a_idx] > 0
            t_labels = np.array(["T cell"] * len(t_indices), dtype=object)
            t_labels[cd4_expr & ~cd8_expr] = "CD4-positive, alpha-beta T cell"
            t_labels[cd8_expr & ~cd4_expr] = "CD8-positive, alpha-beta T cell"
            adata.obs.loc[t_mask, "cell_type"] = t_labels

    # Print summary
    type_counts = adata.obs["cell_type"].value_counts()
    print(f"\n  Cell type annotation summary ({len(clusters)} clusters → "
          f"{len(type_counts)} types):")
    for ct, n in type_counts.items():
        print(f"    {ct:<45} {n:>6,} cells")

    return adata


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # ==================================================================
    # STEP 1: Load all YC scRNA-seq tissues
    # ==================================================================
    print("=" * 70)
    print("STEP 1: Load Cell Ranger matrices (8 YC scRNA-seq tissues)")
    print("=" * 70)

    tissue_adatas = []
    for tissue_dir, tissue_label in TISSUES.items():
        print(f"\n  Loading {tissue_dir} ({tissue_label})...")
        adata = load_tissue(tissue_dir, tissue_label)
        print(f"    Raw: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
        adata = qc_filter(adata)
        tissue_adatas.append(adata)

    print("\n  Concatenating all tissues...")
    combined = ad.concat(tissue_adatas, join="inner")
    combined.obs_names_make_unique()
    print(f"  Combined: {combined.n_obs:,} cells × {combined.n_vars:,} genes")
    print(f"  Tissues: {combined.obs['tissue'].value_counts().to_dict()}")

    # ==================================================================
    # STEP 2: Annotate cell types (two-pass)
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 2: Annotate cell types (two-pass: original + BORDERLINE)")
    print("=" * 70)

    combined = annotate_cells(combined)

    # ==================================================================
    # STEP 3: Cell count audit (expanded)
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 3: Cell count audit (expanded type set)")
    print("=" * 70)

    # All 35 primary types for completeness
    all_35_types = [
        "B cell", "CD4-positive, alpha-beta T cell",
        "CD8-positive, alpha-beta T cell", "T cell",
        "adventitial cell", "basal cell", "bladder urothelial cell",
        "classical monocyte", "endothelial cell",
        "enterocyte of epithelium of large intestine",
        "epithelial cell", "fibroblast", "fibroblast of cardiac tissue",
        "granulocyte", "hematopoietic precursor cell",
        "hematopoietic stem cell", "hepatocyte", "intermediate monocyte",
        "large intestine goblet cell",
        "luminal epithelial cell of mammary gland",
        "macrophage", "mature NK T cell", "mesenchymal stem cell",
        "mesenchymal stem cell of adipose tissue", "monocyte",
        "myeloid dendritic cell", "myeloid leukocyte",
        "natural killer cell", "neutrophil", "non-classical monocyte",
        "pancreatic acinar cell", "pancreatic ductal cell",
        "plasma cell", "smooth muscle cell", "stromal cell",
    ]

    audit_rows = []
    for ct in all_35_types:
        n = int((combined.obs["cell_type"] == ct).sum())
        if n >= 500:
            status = "PASS"
        elif n >= MIN_CELLS:
            status = "BORDERLINE"
        else:
            status = "FAIL"
        audit_rows.append({"cell_type": ct, "n_cells_YC": n, "status": status})
        if n > 0:
            print(f"  {ct:<50} {n:>6,}  {status}")

    audit_df = pd.DataFrame(audit_rows)
    audit_df.to_csv(OUTPUT_DIR / "cell_count_audit_expanded.csv", index=False)

    n_pass = (audit_df["status"] == "PASS").sum()
    n_border = (audit_df["status"] == "BORDERLINE").sum()
    n_fail = (audit_df["status"] == "FAIL").sum()
    n_usable = n_pass + n_border
    print(f"\n  PASS (≥500): {n_pass}, BORDERLINE (200-499): {n_border}, "
          f"FAIL (<200): {n_fail}")
    print(f"  Usable types (≥{MIN_CELLS}): {n_usable}")

    # Filter to types with ≥200 cells
    usable_types = audit_df[audit_df["status"].isin(["PASS", "BORDERLINE"])][
        "cell_type"
    ].tolist()

    # Only keep types that also exist in the Tabula 35-type centroids
    tabula_h_centroids = pd.read_csv(TABULA_HUMAN_CENTROIDS_PATH, index_col=0)
    tabula_types = set(tabula_h_centroids.index)
    usable_in_tabula = [ct for ct in usable_types if ct in tabula_types]
    dropped = [ct for ct in usable_types if ct not in tabula_types]
    if dropped:
        print(f"\n  Dropped (not in Tabula 35-type): {dropped}")
    usable_types = usable_in_tabula
    print(f"  Types in both Sun2023 (≥{MIN_CELLS}) AND Tabula 35-type: {len(usable_types)}")
    for ct in usable_types:
        n = int((combined.obs["cell_type"] == ct).sum())
        print(f"    {ct:<50} {n:>6,} cells")

    # Filter combined to usable types
    combined = combined[combined.obs["cell_type"].isin(usable_types)].copy()
    print(f"\n  Filtered dataset: {combined.n_obs:,} cells")

    # ==================================================================
    # STEP 4: Normalize and restrict to ortholog space
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 4: Normalize (CPM + log1p) and restrict to ortholog space")
    print("=" * 70)

    ortho = pd.read_csv(ORTHOLOG_PATH)
    mouse_to_human = dict(zip(ortho["mouse_gene_name"], ortho["human_ensembl_id"]))
    print(f"  Ortholog map: {len(mouse_to_human):,} mouse→human gene pairs")

    sun_genes = set(combined.var_names)
    shared_mouse_genes = sorted(sun_genes & set(mouse_to_human.keys()))
    print(f"  Sun2023 genes: {len(sun_genes):,}")
    print(f"  Shared with orthologs: {len(shared_mouse_genes):,} / {len(mouse_to_human):,}")

    combined = combined[:, shared_mouse_genes].copy()

    # Normalize: CPM + log1p
    sc.pp.normalize_total(combined, target_sum=1e4)
    sc.pp.log1p(combined)
    print(f"  Normalized: {combined.n_obs:,} cells × {combined.n_vars:,} genes")

    # Map gene names to human Ensembl IDs
    human_ids = [mouse_to_human[g] for g in combined.var_names]
    combined.var_names = pd.Index(human_ids)
    combined.var_names_make_unique()

    # Match gene space to Tabula 35-type centroids
    full_gene_set = list(tabula_h_centroids.columns)
    n_full = len(full_gene_set)
    existing_genes = set(combined.var_names)
    missing_genes = [g for g in full_gene_set if g not in existing_genes]
    print(f"  Genes in full space: {n_full}")
    print(f"  Present: {len(existing_genes)}, Missing (zero-filled): {len(missing_genes)}")

    if missing_genes:
        X_old = combined.X
        if sp.issparse(X_old):
            X_old = X_old.toarray()
        X_full = np.zeros((combined.n_obs, n_full), dtype=np.float32)
        for j_new, g in enumerate(combined.var_names):
            if g in full_gene_set:
                j_full = full_gene_set.index(g)
                X_full[:, j_full] = X_old[:, j_new]
        obs_df = combined.obs.copy()
        combined = ad.AnnData(
            X=sp.csr_matrix(X_full),
            obs=obs_df,
            var=pd.DataFrame(index=full_gene_set),
        )
    else:
        combined = combined[:, full_gene_set].copy()

    print(f"  Final Sun2023 shape: {combined.n_obs:,} cells × {combined.n_vars:,} genes")
    zero_fill_rate = len(missing_genes) / n_full
    print(f"  Zero-fill rate: {zero_fill_rate:.1%}")

    # ==================================================================
    # STEP 5: Compute centroids
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 5: Compute centroids")
    print("=" * 70)

    print("\n--- Sun2023 mouse centroids (expanded) ---")
    sun_centroids = compute_centroids(combined, "cell_type")

    # Tabula human centroids — restrict to shared types
    shared_types = sorted(set(sun_centroids.index) & set(tabula_h_centroids.index))
    n_shared = len(shared_types)
    print(f"\n  Shared types for Procrustes: {n_shared}")
    for ct in shared_types:
        sun_n = int((combined.obs["cell_type"] == ct).sum())
        print(f"    {ct:<50} Sun2023: {sun_n:>5,} cells")

    if n_shared < 4:
        print("\n  *** STOP: Fewer than 4 shared types. Cannot run Procrustes. ***")
        return

    sun_centroids = sun_centroids.loc[shared_types]
    tabula_h_sub = tabula_h_centroids.loc[shared_types]

    # Save centroids for V1 validation (added 2026-03-21)
    centroid_save_path = Path("data/centroids/sun2023_15type_centroids.csv")
    centroid_save_path.parent.mkdir(parents=True, exist_ok=True)
    sun_centroids.to_csv(centroid_save_path)
    print(f"\n  Centroids saved to {centroid_save_path} ({sun_centroids.shape[0]} types × {sun_centroids.shape[1]} genes)")

    # ==================================================================
    # STEP 6: Procrustes — Sun2023 expanded → Tabula human
    # ==================================================================
    print("\n" + "=" * 70)
    print(f"STEP 6: Procrustes — Sun2023 expanded ({n_shared} types) → Tabula human")
    print("=" * 70)

    print("\n--- PCA on combined centroids ---")
    human_pca, mouse_pca, pca_model, types_list = pca_reduce_centroids(
        tabula_h_sub, sun_centroids, variance_threshold=VARIANCE_THRESHOLD
    )

    print("\n--- Procrustes: Sun2023 mouse → Tabula human ---")
    result = procrustes_align(human_pca, mouse_pca)

    print(f"\n--- Permutation test ({N_PERMUTATIONS:,} iterations) ---")
    p_val, null_dist = permutation_test(
        human_pca, mouse_pca, N_PERMUTATIONS, RANDOM_SEED
    )

    obs_null = result.distance / np.median(null_dist)

    print("\n--- Per-type residuals ---")
    residuals = compute_residual_vectors(result, types_list)
    residual_mags = {ct: float(np.linalg.norm(residuals[ct])) for ct in types_list}

    # Sort by magnitude
    sorted_types = sorted(residual_mags, key=residual_mags.get, reverse=True)
    total_ssr = sum(v**2 for v in residual_mags.values())
    print(f"\n  Residual ranking (n={n_shared}):")
    for i, ct in enumerate(sorted_types, 1):
        pct = residual_mags[ct] ** 2 / total_ssr * 100
        print(f"    {i:>2}. {ct:<50} {residual_mags[ct]:>8.3f}  ({pct:.1f}% SSR)")

    # ==================================================================
    # STEP 7: Rigidity ranking correlation
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 7: Rigidity ranking correlation vs primary 35-type")
    print("=" * 70)

    primary_residuals = pd.read_csv(RESIDUALS_RANKED_PATH)
    primary_mag_dict = dict(
        zip(primary_residuals["cell_type"], primary_residuals["residual_magnitude"])
    )

    matched_types = sorted(set(residual_mags.keys()) & set(primary_mag_dict.keys()))
    n_matched = len(matched_types)
    print(f"  Matched types for ranking: {n_matched}")

    if n_matched >= 4:
        sun_mags = [residual_mags[ct] for ct in matched_types]
        primary_mags = [primary_mag_dict[ct] for ct in matched_types]
        rho, rho_p = spearmanr(sun_mags, primary_mags)
        print(f"  Spearman ρ = {rho:.3f}, p = {rho_p:.4f} (n={n_matched})")

        # Print per-type rank comparison
        sun_rank = pd.Series(sun_mags, index=matched_types).rank(ascending=False)
        primary_rank = pd.Series(primary_mags, index=matched_types).rank(ascending=False)
        print(f"\n  {'Cell type':<50} {'Sun rank':>8} {'Primary rank':>12}")
        print(f"  {'-'*70}")
        for ct in matched_types:
            print(f"  {ct:<50} {int(sun_rank[ct]):>8} {int(primary_rank[ct]):>12}")
    else:
        rho, rho_p = float("nan"), float("nan")
        print(f"  Too few matched types ({n_matched}) for ranking correlation")

    # ==================================================================
    # STEP 8: T1-B Negative Control
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 8: T1-B NEGATIVE CONTROL (obs/null ratio comparison)")
    print("=" * 70)

    # Load primary 35-type results
    with open(PRIMARY_RESULTS_PATH) as f:
        primary_35 = json.load(f)

    tabula_35_distance = primary_35["procrustes"]["distance"]
    tabula_35_null_median = primary_35["permutation_test"][
        "null_distribution_summary"
    ]["median"]
    tabula_35_ratio = tabula_35_distance / tabula_35_null_median

    # Permuted baseline: median(null) / median(null) = 1.0 by construction
    permuted_ratio = 1.0

    print(f"\n  Obs/Null ratios (lower = stronger signal):")
    print(f"    Sun2023 expanded ({n_shared} types):     {obs_null:.4f}")
    print(f"    Tabula 35-type (primary):              {tabula_35_ratio:.4f}")
    print(f"    Permuted baseline (random pairings):   {permuted_ratio:.4f}")
    print(f"\n  T1-B test: Sun2023 obs/null ({obs_null:.3f}) < "
          f"permuted ({permuted_ratio:.3f})? "
          f"{'PASS' if obs_null < permuted_ratio else 'FAIL'}")
    print(f"  Signal consistency: Sun2023 ratio / Tabula ratio = "
          f"{obs_null / tabula_35_ratio:.2f}")

    if obs_null < tabula_35_ratio * 1.5:
        t1b_verdict = "PASS — obs/null in same range as Tabula primary"
    else:
        t1b_verdict = (
            "WEAKER — obs/null higher than Tabula (cross-study noise expected)"
        )
    print(f"  T1-B verdict: {t1b_verdict}")

    # ==================================================================
    # STEP 9: Verdict per decision tree
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 9: VERDICT")
    print("=" * 70)

    if n_shared >= 10 and rho >= 0.50 and rho_p < 0.05:
        verdict = "PASS"
        detail = (
            f"n={n_shared} ≥ 10, ρ={rho:.3f} ≥ 0.50, p={rho_p:.4f} < 0.05. "
            f"T1-A resolved."
        )
    elif n_shared >= 10 and (rho < 0.50 or rho_p >= 0.05):
        verdict = "PARTIAL"
        detail = (
            f"n={n_shared} ≥ 10, but ρ={rho:.3f} "
            f"{'< 0.50' if rho < 0.50 else '≥ 0.50'}, "
            f"p={rho_p:.4f} "
            f"{'≥ 0.05' if rho_p >= 0.05 else '< 0.05'}. "
            f"Ranking unstable. Flag to advisor."
        )
    else:
        verdict = "UNDERPOWERED"
        detail = (
            f"n={n_shared} < 10 after expansion. "
            f"Insufficient types for powered Spearman. Flag to advisor."
        )

    print(f"\n  VERDICT: {verdict}")
    print(f"  {detail}")
    print(f"\n  Procrustes: p={p_val:.4f}, obs/null={obs_null:.3f}, "
          f"scaling={result.scaling:.3f}")
    print(f"  Rigidity: ρ={rho:.3f}, p={rho_p:.4f}, n={n_matched}")
    print(f"  T1-B: {t1b_verdict}")

    # ==================================================================
    # STEP 10: Save outputs
    # ==================================================================
    print("\n" + "=" * 70)
    print("Saving outputs")
    print("=" * 70)

    results = {
        "diagnostic": "Sun2023 expanded 10x replication — T1-A with BORDERLINE types",
        "date": "2026-03-15",
        "dataset": {
            "name": "Sun et al. 2023 (Innovation, CAS Beijing)",
            "protocol": "10x Chromium 3' v3",
            "condition": "YC (young sedentary control)",
            "tissues_loaded": list(TISSUES.keys()),
            "n_tissues": len(TISSUES),
            "total_cells_post_qc": int(combined.n_obs),
            "n_genes": int(combined.n_vars),
            "zero_fill_rate": float(zero_fill_rate),
            "annotation_method": "Two-pass: (1) original 14 markers, (2) 8 BORDERLINE markers on unassigned clusters",
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
        },
        "t1b_negative_control": {
            "sun2023_expanded_obs_null": float(obs_null),
            "tabula_35_obs_null": float(tabula_35_ratio),
            "permuted_baseline": float(permuted_ratio),
            "sun2023_lower_than_permuted": bool(obs_null < permuted_ratio),
            "verdict": t1b_verdict,
        },
        "interpretation": {
            "verdict": verdict,
            "detail": detail,
        },
        "cell_count_audit": audit_rows,
        "random_seed": RANDOM_SEED,
    }

    with open(OUTPUT_DIR / "sun2023_expanded.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"  Results: {OUTPUT_DIR / 'sun2023_expanded.json'}")

    np.save(OUTPUT_DIR / "null_distribution.npy", null_dist)

    # Save ranking comparison as CSV
    if n_matched >= 4:
        rank_df = pd.DataFrame({
            "cell_type": matched_types,
            "sun2023_residual": [residual_mags[ct] for ct in matched_types],
            "primary_residual": [primary_mag_dict[ct] for ct in matched_types],
        })
        rank_df.to_csv(OUTPUT_DIR / "ranking_comparison.csv", index=False)

    # --- Plots ---
    # 1. Null distribution
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(null_dist, bins=50, alpha=0.7, color="forestgreen", edgecolor="white")
    ax.axvline(
        result.distance, color="red", linewidth=2,
        label=f"Observed (d={result.distance:.2f})"
    )
    ax.set_title(
        f"Sun2023 expanded ({n_shared} types) vs Tabula human\n"
        f"p={p_val:.4f}, obs/null={obs_null:.3f}, ρ={rho:.3f}"
    )
    ax.set_xlabel("Procrustes distance")
    ax.set_ylabel("Count")
    ax.legend()
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "null_distribution.png", dpi=150, bbox_inches="tight")
    plt.close()

    # 2. Rigidity scatter: Sun2023 vs primary
    if n_matched >= 4:
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.scatter(
            [primary_mag_dict[ct] for ct in matched_types],
            [residual_mags[ct] for ct in matched_types],
            s=60, c="steelblue", edgecolors="navy", linewidths=0.5, zorder=3,
        )
        for ct in matched_types:
            short_name = ct[:25] + "..." if len(ct) > 25 else ct
            ax.annotate(
                short_name,
                (primary_mag_dict[ct], residual_mags[ct]),
                fontsize=7, ha="left", va="bottom",
                xytext=(4, 4), textcoords="offset points",
            )
        ax.set_xlabel("Primary 35-type residual magnitude")
        ax.set_ylabel("Sun2023 expanded residual magnitude")
        ax.set_title(
            f"Rigidity ranking: Sun2023 vs Primary\n"
            f"Spearman ρ={rho:.3f}, p={rho_p:.4f}, n={n_matched}"
        )
        # Add diagonal reference
        lims = [
            min(ax.get_xlim()[0], ax.get_ylim()[0]),
            max(ax.get_xlim()[1], ax.get_ylim()[1]),
        ]
        ax.plot(lims, lims, "k--", alpha=0.3, zorder=1)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        plt.tight_layout()
        fig.savefig(
            OUTPUT_DIR / "rigidity_scatter.png", dpi=150, bbox_inches="tight"
        )
        plt.close()

    # 3. T1-B bar chart
    fig, ax = plt.subplots(figsize=(6, 4))
    bars = ax.bar(
        ["Sun2023\nexpanded", "Tabula\n35-type", "Permuted\n(baseline)"],
        [obs_null, tabula_35_ratio, permuted_ratio],
        color=["forestgreen", "steelblue", "gray"],
        edgecolor="black", linewidth=0.5,
    )
    ax.set_ylabel("Obs/Null ratio")
    ax.set_title("T1-B Negative Control: Obs/Null Comparison")
    ax.axhline(1.0, color="red", linestyle="--", alpha=0.5, label="Random baseline")
    ax.legend()
    for bar, val in zip(bars, [obs_null, tabula_35_ratio, permuted_ratio]):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                f"{val:.3f}", ha="center", va="bottom", fontsize=9)
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "t1b_comparison.png", dpi=150, bbox_inches="tight")
    plt.close()

    print(f"  Plots: null_distribution.png, rigidity_scatter.png, t1b_comparison.png")

    print(f"\n{'=' * 70}")
    print(f"FINAL VERDICT: {verdict}")
    print(f"{'=' * 70}")
    print(f"  {detail}")
    print(f"  T1-B: {t1b_verdict}")


if __name__ == "__main__":
    main()
