#!/usr/bin/env python3
"""
CellWarp — T1-A Replication: Sun et al. 2023 (10x Chromium 3' v3)

Runs the identical Procrustes pipeline on Sun et al. 2023 young sedentary control
mouse data (10x Chromium, CAS Beijing) against Tabula Sapiens human data.

This is the direct test of whether the cross-species geometric signal replicates
with an independent 10x-based mouse atlas.

Biology
-------
Sun et al. profiled 14 mouse tissues with 10x Chromium 3' v3 (9 scRNA-seq tissues)
in a 2×2 exercise × age study. We use the young sedentary control (YC) arm only.
Cell type annotation uses canonical marker genes — standard immunology/biology,
not reference transfer from Tabula (preserving independence).

Steps
-----
  1. Load Cell Ranger outputs for 8 YC scRNA-seq tissues
  2. QC: min 200 genes, max 20% mitochondrial
  3. Annotate cell types using canonical marker genes (cluster-based)
  4. Map to our ontology, produce cell count audit
  5. Normalize: CPM + log1p (identical to primary pipeline)
  6. Restrict to 16,959 shared 1:1 orthologs
  7. Compute centroids per cell type
  8. Run Procrustes + permutation test (10,000 iterations)
  9. Run same pipeline on Tabula data restricted to same type set
  10. Compare: p-value, obs/null, scaling, rigidity ranking correlation

Output
------
  output/validation/sun2023_replication/
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

OUTPUT_DIR = Path("output/validation/sun2023_replication")
EXTRACTED_DIR = Path("data/replication/sun2023/extracted")
TABULA_HUMAN_PATH = Path("data/phase1/human_qc.h5ad")
TABULA_MOUSE_PATH = Path("data/phase1/mouse_qc.h5ad")
ORTHOLOG_PATH = Path("data/phase1/orthologs_human_mouse.csv")

N_PERMUTATIONS = 10_000
VARIANCE_THRESHOLD = 0.95

# Tissues to load (YC scRNA-seq only)
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

# Canonical marker genes for broad cell type annotation (mouse gene symbols)
# These are textbook markers — not derived from any atlas
MARKERS = {
    "hepatocyte": {"pos": ["Alb", "Ttr", "Apoa1"], "neg": []},
    "B cell": {"pos": ["Cd79a", "Cd79b", "Ms4a1"], "neg": ["Cd3e"]},
    "T cell": {"pos": ["Cd3e", "Cd3d"], "neg": ["Cd79a"]},
    "endothelial cell": {"pos": ["Pecam1", "Cdh5"], "neg": []},
    "macrophage": {"pos": ["Adgre1", "Csf1r", "C1qa"], "neg": ["Cd3e", "Cd79a"]},
    "monocyte": {"pos": ["Ly6c2", "Csf1r"], "neg": ["Cd3e", "Cd79a", "Adgre1"]},
    "natural killer cell": {"pos": ["Nkg7", "Klrb1c", "Gzma"], "neg": ["Cd3e"]},
    "neutrophil": {"pos": ["S100a8", "S100a9", "Ly6g"], "neg": ["Cd3e", "Cd79a"]},
    "epithelial cell": {"pos": ["Epcam", "Krt8", "Krt18"], "neg": ["Pecam1"]},
    "fibroblast": {"pos": ["Col1a1", "Col1a2", "Dcn"], "neg": ["Pecam1", "Epcam"]},
    "smooth muscle cell": {"pos": ["Acta2", "Myh11", "Tagln"], "neg": ["Pecam1"]},
    "plasma cell": {"pos": ["Sdc1", "Xbp1", "Jchain"], "neg": ["Ms4a1"]},
    "myeloid dendritic cell": {"pos": ["Flt3", "Itgax", "H2-Aa"], "neg": ["Cd3e", "Cd79a", "Adgre1"]},
    "granulocyte": {"pos": ["S100a8", "S100a9"], "neg": ["Cd3e", "Cd79a", "Csf1r"]},
}


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
    # Mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    n_before = adata.n_obs
    # Filter
    adata = adata[adata.obs["n_genes_by_counts"] >= 200].copy()
    adata = adata[adata.obs["pct_counts_mt"] <= 20].copy()
    n_after = adata.n_obs
    print(f"  QC: {n_before:,} → {n_after:,} cells ({n_before - n_after:,} removed)")
    return adata


def score_markers(adata: ad.AnnData) -> pd.DataFrame:
    """Score each cell for each cell type using marker genes.

    Returns a DataFrame of shape (n_cells, n_types) with marker scores.
    Score = mean(positive markers present) - mean(negative markers present).
    Uses raw counts (>0 detection) for robustness with sparse data.
    """
    # Work with the sparse matrix directly
    X = adata.X
    if sp.issparse(X):
        X = X.toarray()
    gene_names = list(adata.var_names)
    gene_idx = {g: i for i, g in enumerate(gene_names)}

    scores = {}
    for ct, markers in MARKERS.items():
        pos_genes = [g for g in markers["pos"] if g in gene_idx]
        neg_genes = [g for g in markers["neg"] if g in gene_idx]

        if not pos_genes:
            scores[ct] = np.zeros(adata.n_obs)
            continue

        # Positive: fraction of positive markers detected (>0)
        pos_detected = np.zeros(adata.n_obs)
        for g in pos_genes:
            pos_detected += (X[:, gene_idx[g]] > 0).astype(float)
        pos_score = pos_detected / len(pos_genes)

        # Negative: penalize if negative markers detected
        neg_detected = np.zeros(adata.n_obs)
        for g in neg_genes:
            neg_detected += (X[:, gene_idx[g]] > 0).astype(float)
        neg_score = neg_detected / max(len(neg_genes), 1)

        scores[ct] = pos_score - 0.5 * neg_score

    return pd.DataFrame(scores, index=adata.obs_names)


def annotate_cells(adata: ad.AnnData) -> ad.AnnData:
    """Annotate cells using marker gene scoring.

    Two-pass approach:
    1. Cluster cells (Leiden)
    2. Score each cluster by marker gene enrichment
    3. Assign the dominant cell type to each cluster
    """
    print("  Annotating cell types via marker gene scoring...")

    # Normalize a copy for clustering (don't modify raw counts yet)
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)
    sc.pp.highly_variable_genes(adata_norm, n_top_genes=2000, flavor="seurat")
    sc.pp.pca(adata_norm, n_comps=30, use_highly_variable=True)
    sc.pp.neighbors(adata_norm, n_neighbors=15)
    sc.tl.leiden(adata_norm, resolution=1.0, random_state=RANDOM_SEED)

    # Transfer cluster labels back to original
    adata.obs["leiden"] = adata_norm.obs["leiden"].values

    # Score markers on raw counts
    marker_scores = score_markers(adata)

    # For each cluster, compute mean score per cell type
    cluster_labels = adata.obs["leiden"]
    clusters = sorted(cluster_labels.unique(), key=int)

    cluster_assignments = {}
    for cl in clusters:
        mask = cluster_labels == cl
        mean_scores = marker_scores.loc[mask].mean()
        best_type = mean_scores.idxmax()
        best_score = mean_scores[best_type]

        # Only assign if score is meaningful (>0.3 = at least 30% of pos markers detected)
        if best_score >= 0.3:
            cluster_assignments[cl] = best_type
        else:
            cluster_assignments[cl] = "unassigned"

    # Assign cell types
    adata.obs["cell_type"] = cluster_labels.map(cluster_assignments).values
    adata.obs["cell_type"] = adata.obs["cell_type"].astype(str)

    # Tissue-aware hepatocyte rescue: liver cells with Alb > 0
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

    # CD4/CD8 T cell refinement within T cell clusters
    t_mask = adata.obs["cell_type"] == "T cell"
    if t_mask.sum() > 0:
        cd4_idx = gene_idx.get("Cd4")
        cd8a_idx = gene_idx.get("Cd8a")

        if cd4_idx is not None and cd8a_idx is not None:
            t_indices = np.where(t_mask.values)[0]
            cd4_expr = X_dense[t_indices, cd4_idx] > 0
            cd8_expr = X_dense[t_indices, cd8a_idx] > 0

            # Use object dtype to avoid numpy string truncation
            t_labels = np.array(["T cell"] * len(t_indices), dtype=object)
            t_labels[cd4_expr & ~cd8_expr] = "CD4-positive, alpha-beta T cell"
            t_labels[cd8_expr & ~cd4_expr] = "CD8-positive, alpha-beta T cell"

            adata.obs.loc[t_mask, "cell_type"] = t_labels

    # Print summary
    type_counts = adata.obs["cell_type"].value_counts()
    print(f"  Cell type annotation summary ({len(clusters)} clusters → {len(type_counts)} types):")
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

    # Concatenate all tissues
    print("\n  Concatenating all tissues...")
    combined = ad.concat(tissue_adatas, join="inner")
    combined.obs_names_make_unique()
    print(f"  Combined: {combined.n_obs:,} cells × {combined.n_vars:,} genes")
    print(f"  Tissues: {combined.obs['tissue'].value_counts().to_dict()}")

    # ==================================================================
    # STEP 2: Annotate cell types
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 2: Annotate cell types using canonical marker genes")
    print("=" * 70)

    combined = annotate_cells(combined)

    # ==================================================================
    # STEP 3: Map to our ontology and cell count audit
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 3: Cell count audit")
    print("=" * 70)

    # Our target types (the ones we expect to find)
    target_types = [
        "B cell",
        "CD4-positive, alpha-beta T cell",
        "CD8-positive, alpha-beta T cell",
        "endothelial cell",
        "hepatocyte",
        "macrophage",
        "monocyte",
        "natural killer cell",
        "neutrophil",
        "epithelial cell",
        "fibroblast",
        "smooth muscle cell",
        "plasma cell",
        "myeloid dendritic cell",
        "granulocyte",
    ]

    audit_rows = []
    for ct in target_types:
        n = (combined.obs["cell_type"] == ct).sum()
        if n >= 500:
            status = "PASS"
        elif n >= 200:
            status = "BORDERLINE"
        else:
            status = "FAIL"
        audit_rows.append({"cell_type": ct, "n_cells_YC": int(n), "status": status})
        print(f"  {ct:<45} {n:>6,}  {status}")

    audit_df = pd.DataFrame(audit_rows)
    audit_df.to_csv(OUTPUT_DIR / "cell_count_audit.csv", index=False)

    n_pass = (audit_df["status"] == "PASS").sum()
    n_border = (audit_df["status"] == "BORDERLINE").sum()
    n_fail = (audit_df["status"] == "FAIL").sum()
    print(f"\n  PASS: {n_pass}, BORDERLINE: {n_border}, FAIL: {n_fail}")

    # Filter to types with ≥200 cells
    usable_types = audit_df[audit_df["status"].isin(["PASS", "BORDERLINE"])]["cell_type"].tolist()
    print(f"  Usable types (≥200 cells): {len(usable_types)}")

    if len(usable_types) < 10:
        print("\n  *** STOP: Fewer than 10 usable types. Flagged, not patched. ***")
        with open(OUTPUT_DIR / "STOP_insufficient_types.txt", "w") as f:
            f.write(f"Only {len(usable_types)} types with ≥200 cells.\n")
            f.write(audit_df.to_string())
        return

    # Filter combined to usable types
    combined = combined[combined.obs["cell_type"].isin(usable_types)].copy()
    print(f"  Filtered to usable types: {combined.n_obs:,} cells")

    # ==================================================================
    # STEP 4: Normalize and restrict to ortholog space
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 4: Normalize (CPM + log1p) and restrict to ortholog space")
    print("=" * 70)

    # Load ortholog mapping
    ortho = pd.read_csv(ORTHOLOG_PATH)
    mouse_to_human = dict(zip(ortho["mouse_gene_name"], ortho["human_ensembl_id"]))
    print(f"  Ortholog map: {len(mouse_to_human):,} mouse→human gene pairs")

    # Find overlap between Sun2023 genes and orthologs
    sun_genes = set(combined.var_names)
    shared_mouse_genes = sorted(sun_genes & set(mouse_to_human.keys()))
    print(f"  Sun2023 genes: {len(sun_genes):,}")
    print(f"  Shared with orthologs: {len(shared_mouse_genes):,} / {len(mouse_to_human):,}")

    # Subset to shared ortholog genes
    combined = combined[:, shared_mouse_genes].copy()

    # Normalize: CPM + log1p
    sc.pp.normalize_total(combined, target_sum=1e4)
    sc.pp.log1p(combined)
    print(f"  Normalized: {combined.n_obs:,} cells × {combined.n_vars:,} genes")

    # Map gene names to human Ensembl IDs (to match Tabula Sapiens gene space)
    human_ids = [mouse_to_human[g] for g in combined.var_names]
    combined.var_names = pd.Index(human_ids)
    combined.var_names_make_unique()

    # Pad missing orthologs with zeros to get full 16,959 gene space
    tabula_h = ad.read_h5ad(TABULA_HUMAN_PATH)
    full_gene_set = list(tabula_h.var_names)
    n_full = len(full_gene_set)

    existing_genes = set(combined.var_names)
    missing_genes = [g for g in full_gene_set if g not in existing_genes]
    print(f"  Genes in full space: {n_full}")
    print(f"  Present: {len(existing_genes)}, Missing (zero-filled): {len(missing_genes)}")

    if missing_genes:
        # Create zero columns for missing genes by building full matrix
        existing_idx = [full_gene_set.index(g) for g in combined.var_names if g in full_gene_set]
        X_old = combined.X
        if sp.issparse(X_old):
            X_old = X_old.toarray()

        X_full = np.zeros((combined.n_obs, n_full), dtype=np.float32)
        for j_new, g in enumerate(combined.var_names):
            if g in full_gene_set:
                j_full = full_gene_set.index(g)
                X_full[:, j_full] = X_old[:, j_new]

        # Preserve all obs columns
        obs_df = combined.obs.copy()
        combined = ad.AnnData(
            X=sp.csr_matrix(X_full),
            obs=obs_df,
            var=pd.DataFrame(index=full_gene_set),
        )
    else:
        # Reorder to match Tabula gene order
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

    # Sun2023 mouse centroids
    print("\n--- Sun2023 mouse centroids ---")
    sun_centroids = compute_centroids(combined, "cell_type")

    # Tabula Sapiens human centroids (restricted to shared types)
    print("\n--- Tabula Sapiens human centroids (restricted to shared types) ---")
    tabula_h_sub = tabula_h[tabula_h.obs["cell_type"].isin(usable_types)].copy()
    tabula_h_centroids = compute_centroids(tabula_h_sub, "cell_type")

    # Find intersection of types
    shared_types = sorted(set(sun_centroids.index) & set(tabula_h_centroids.index))
    n_shared = len(shared_types)
    print(f"\n  Shared types for Procrustes: {n_shared}")
    for ct in shared_types:
        print(f"    {ct}")

    if n_shared < 4:
        print("\n  *** STOP: Fewer than 4 shared types. Cannot run Procrustes. ***")
        return

    # Restrict centroids to shared types
    sun_centroids = sun_centroids.loc[shared_types]
    tabula_h_centroids = tabula_h_centroids.loc[shared_types]

    # ==================================================================
    # STEP 6: Procrustes — Sun2023 mouse → Tabula human
    # ==================================================================
    print("\n" + "=" * 70)
    print("COMPARISON A: Sun2023 mouse vs Tabula Sapiens human")
    print("=" * 70)

    print("\n--- PCA on combined centroids ---")
    human_pca_a, mouse_pca_a, pca_a, types_a = pca_reduce_centroids(
        tabula_h_centroids, sun_centroids, variance_threshold=VARIANCE_THRESHOLD
    )

    print("\n--- Procrustes: Sun2023 mouse → Tabula human ---")
    result_a = procrustes_align(human_pca_a, mouse_pca_a)

    print("\n--- Permutation test (10,000 iterations) ---")
    p_a, null_a = permutation_test(human_pca_a, mouse_pca_a, N_PERMUTATIONS, RANDOM_SEED)

    obs_null_a = result_a.distance / np.median(null_a)

    print("\n--- Per-type residuals ---")
    residuals_a = compute_residual_vectors(result_a, types_a)
    residual_mags_a = {ct: float(np.linalg.norm(residuals_a[ct])) for ct in types_a}

    # ==================================================================
    # STEP 7: Primary comparison — Tabula mouse vs Tabula human (same types)
    # ==================================================================
    print("\n" + "=" * 70)
    print("COMPARISON B: Tabula mouse vs Tabula human (same type set)")
    print("=" * 70)

    tabula_m = ad.read_h5ad(TABULA_MOUSE_PATH)
    # Restrict to shared types that exist in Tabula mouse
    tabula_m_types = set(tabula_m.obs["cell_type"].unique())
    primary_shared = sorted(set(shared_types) & tabula_m_types)
    print(f"\n  Types in all three datasets: {len(primary_shared)}")
    for ct in primary_shared:
        print(f"    {ct}")

    if len(primary_shared) < 4:
        print("\n  *** Cannot run comparable primary — Tabula mouse has too few shared types ***")
        # Still save Sun2023 results
        primary_shared = shared_types  # Use what we have

    # Tabula mouse centroids
    print("\n--- Tabula Muris Senis mouse centroids ---")
    tabula_m_sub = tabula_m[tabula_m.obs["cell_type"].isin(primary_shared)].copy()
    tabula_m_centroids = compute_centroids(tabula_m_sub, "cell_type")

    # Tabula human centroids (re-restrict to primary_shared)
    tabula_h_sub2 = tabula_h[tabula_h.obs["cell_type"].isin(primary_shared)].copy()
    tabula_h_centroids2 = compute_centroids(tabula_h_sub2, "cell_type")

    # Re-restrict Sun centroids too
    sun_centroids2 = sun_centroids.loc[[ct for ct in primary_shared if ct in sun_centroids.index]]

    print("\n--- PCA on combined Tabula centroids ---")
    human_pca_b, mouse_pca_b, pca_b, types_b = pca_reduce_centroids(
        tabula_h_centroids2, tabula_m_centroids, variance_threshold=VARIANCE_THRESHOLD
    )

    print("\n--- Procrustes: Tabula mouse → Tabula human ---")
    result_b = procrustes_align(human_pca_b, mouse_pca_b)

    print("\n--- Permutation test (10,000 iterations) ---")
    p_b, null_b = permutation_test(human_pca_b, mouse_pca_b, N_PERMUTATIONS, RANDOM_SEED)

    obs_null_b = result_b.distance / np.median(null_b)

    print("\n--- Per-type residuals ---")
    residuals_b = compute_residual_vectors(result_b, types_b)
    residual_mags_b = {ct: float(np.linalg.norm(residuals_b[ct])) for ct in types_b}

    # ==================================================================
    # STEP 8: Rigidity ranking correlation
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 8: Rigidity ranking correlation")
    print("=" * 70)

    # Load primary 35-type results for ranking comparison
    primary_results_path = Path("output/phase2/scaled_35types/procrustes_results_35.json")
    if primary_results_path.exists():
        with open(primary_results_path) as f:
            primary_35 = json.load(f)
        primary_residuals = {
            ct: np.linalg.norm(primary_35["residuals"][ct]["vector_pca"])
            for ct in primary_35["cell_types"]
            if ct in primary_35["residuals"]
        }
    else:
        # Fallback: use Tabula 6-type result
        primary_residuals = residual_mags_b

    # Match types
    matched_types = sorted(set(residual_mags_a.keys()) & set(primary_residuals.keys()))
    if len(matched_types) >= 4:
        sun_ranks = [residual_mags_a[ct] for ct in matched_types]
        primary_ranks = [primary_residuals[ct] for ct in matched_types]
        rho, rho_p = spearmanr(sun_ranks, primary_ranks)
        print(f"\n  Rigidity ranking Spearman ρ = {rho:.3f}, p = {rho_p:.4f}")
        print(f"  Matched types: {len(matched_types)}")
    else:
        rho, rho_p = float("nan"), float("nan")
        print(f"\n  Too few matched types ({len(matched_types)}) for ranking correlation")

    # ==================================================================
    # STEP 9: Summary and interpretation
    # ==================================================================
    print("\n" + "=" * 70)
    print("SUMMARY COMPARISON TABLE")
    print("=" * 70)

    n_a = len(types_a)
    n_b = len(types_b)

    header = f"{'Comparison':<50} {'n':>3} {'p':>8} {'obs/null':>9} {'scaling':>8} {'ρ':>6}"
    print(header)
    print("-" * len(header))
    print(
        f"{'Primary (Tabula mouse vs Tabula human)':<50} {n_b:>3} {p_b:>8.4f} {obs_null_b:>9.3f} "
        f"{result_b.scaling:>8.3f} {'—':>6}"
    )
    print(
        f"{'Sun2023 YC mouse vs Tabula human':<50} {n_a:>3} {p_a:>8.4f} {obs_null_a:>9.3f} "
        f"{result_a.scaling:>8.3f} {rho:>6.3f}"
    )

    # Interpretation
    print("\n" + "=" * 70)
    print("INTERPRETATION (strict logic, no deviation)")
    print("=" * 70)

    if obs_null_a < 0.80 and p_a < 0.05 and rho > 0.50:
        verdict = "REPLICATION CONFIRMED"
        detail = (
            "obs/null < 0.80 AND p < 0.05 AND rigidity ρ > 0.50. "
            "T1-A existential validation closed."
        )
    elif obs_null_a < 0.80 and p_a < 0.05 and rho <= 0.50:
        verdict = "PARTIAL"
        detail = (
            "obs/null < 0.80 AND p < 0.05 BUT rigidity ρ ≤ 0.50. "
            "Geometry replicates, ranking unstable. Flag to advisor."
        )
    elif obs_null_a > 0.90 or p_a > 0.10:
        verdict = "FAILURE"
        detail = (
            f"obs/null > 0.90 OR p > 0.10. "
            f"Scaling={result_a.scaling:.3f}, zero-fill={zero_fill_rate:.1%}. "
            "Flag to advisor before interpretation."
        )
    else:
        verdict = "AMBIGUOUS"
        detail = (
            f"obs/null={obs_null_a:.3f} (0.80-0.90), p={p_a:.4f} (0.05-0.10). "
            "Flag to advisor, do not interpret."
        )

    print(f"\n  VERDICT: {verdict}")
    print(f"  {detail}")

    # ==================================================================
    # STEP 10: Save outputs
    # ==================================================================
    print("\n" + "=" * 70)
    print("Saving outputs")
    print("=" * 70)

    results = {
        "diagnostic": "Sun2023 10x replication — T1-A independent atlas validation",
        "date": "2026-03-15",
        "dataset": {
            "name": "Sun et al. 2023 (Innovation, CAS Beijing)",
            "protocol": "10x Chromium 3' v3",
            "condition": "YC (young sedentary control, 2-month C57BL/6J male)",
            "tissues_loaded": list(TISSUES.keys()),
            "n_tissues": len(TISSUES),
            "total_cells_post_qc": int(combined.n_obs),
            "n_genes": int(combined.n_vars),
            "gene_overlap_with_orthologs": len(shared_mouse_genes),
            "zero_fill_rate": float(zero_fill_rate),
            "access_method": "OMIX002605 via HTTPS (cncb.ac.cn), .tar.gz Cell Ranger outputs",
        },
        "comparison_a": {
            "name": "Sun2023 YC mouse vs Tabula Sapiens human",
            "n_types": n_a,
            "cell_types": types_a,
            "p_value": float(p_a),
            "distance": float(result_a.distance),
            "obs_null_ratio": float(obs_null_a),
            "scaling": float(result_a.scaling),
            "null_median": float(np.median(null_a)),
            "pca_components": int(pca_a.n_components_),
            "per_type_residuals": {
                ct: {"magnitude": residual_mags_a[ct]} for ct in types_a
            },
        },
        "comparison_b": {
            "name": "Tabula mouse vs Tabula human (same type set)",
            "n_types": n_b,
            "cell_types": types_b,
            "p_value": float(p_b),
            "distance": float(result_b.distance),
            "obs_null_ratio": float(obs_null_b),
            "scaling": float(result_b.scaling),
            "null_median": float(np.median(null_b)),
            "pca_components": int(pca_b.n_components_),
            "per_type_residuals": {
                ct: {"magnitude": residual_mags_b[ct]} for ct in types_b
            },
        },
        "rigidity_ranking": {
            "rho": float(rho) if not np.isnan(rho) else None,
            "p_value": float(rho_p) if not np.isnan(rho_p) else None,
            "n_matched_types": len(matched_types),
            "matched_types": matched_types,
        },
        "interpretation": {
            "verdict": verdict,
            "detail": detail,
        },
        "cell_count_audit": audit_rows,
        "random_seed": RANDOM_SEED,
    }

    with open(OUTPUT_DIR / "sun2023_replication.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"  Results: {OUTPUT_DIR / 'sun2023_replication.json'}")

    np.save(OUTPUT_DIR / "null_a_sun2023.npy", null_a)
    # NOTE: When the cell-type triple-intersection reduces to the same 6 types
    # as the sibling hca_centroid_comparison workflow's comparison_b, the saved
    # null_b_tabula.npy will be byte-identical with the sibling's
    # null_b_tabula_hvm.npy. Both are diagnostic preserves; neither feeds
    # Table 1 or manuscript claims. See R21 MD5-16 investigation.
    np.save(OUTPUT_DIR / "null_b_tabula.npy", null_b)

    # Plot: side-by-side null distributions
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    ax.hist(null_a, bins=50, alpha=0.7, color="forestgreen", edgecolor="white")
    ax.axvline(result_a.distance, color="red", linewidth=2, label=f"Observed (d={result_a.distance:.2f})")
    ax.set_title(f"Sun2023 mouse vs Tabula human ({n_a} types)\np={p_a:.4f}, obs/null={obs_null_a:.3f}")
    ax.set_xlabel("Procrustes distance")
    ax.set_ylabel("Count")
    ax.legend()

    ax = axes[1]
    ax.hist(null_b, bins=50, alpha=0.7, color="darkorange", edgecolor="white")
    ax.axvline(result_b.distance, color="red", linewidth=2, label=f"Observed (d={result_b.distance:.2f})")
    ax.set_title(f"Tabula mouse vs Tabula human ({n_b} types)\np={p_b:.4f}, obs/null={obs_null_b:.3f}")
    ax.set_xlabel("Procrustes distance")
    ax.set_ylabel("Count")
    ax.legend()

    plt.suptitle(f"Sun2023 10x Replication — VERDICT: {verdict}", fontsize=13, fontweight="bold")
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "null_distributions.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Plot: null_distributions.png")

    # Save combined h5ad for reproducibility
    combined.write_h5ad(Path("data/replication/sun2023/sun2023_yc.h5ad"))
    print(f"  Saved: data/replication/sun2023/sun2023_yc.h5ad")

    print(f"\n{'=' * 70}")
    print(f"FINAL VERDICT: {verdict}")
    print(f"{'=' * 70}")
    print(f"  {detail}")


if __name__ == "__main__":
    main()
