#!/usr/bin/env python3
"""
CellWarp — HSC Independent Validation: Progenitor Finding Retest

Validates whether the progenitor divergence finding recovers when the unreliable
Tabula Muris Senis HSC centroid (only 23 10x cells, 98.9% Smart-seq2) is replaced
with an independent mouse HSC centroid from other Census datasets.

Biology
-------
Hematopoietic stem cells (HSCs) sit at the apex of the blood lineage hierarchy.
In Tabula Muris Senis, HSCs are almost entirely profiled by Smart-seq2 (full-length),
not 10x Chromium (3' UMI). This protocol mismatch makes the HSC centroid unreliable
when restricting to 10x-only data (only 23 cells remain).

Census inventory reveals 0 10x mouse HSC cells outside TMS. The only independent
source is 64,160 sci-RNA-seq3 cells from Cao et al. (mouse embryonic development).
sci-RNA-seq3 is UMI-based 3'-end capture — methodologically closer to 10x than
Smart-seq2 — but the cells are embryonic, not adult. This is an imperfect but
informative validation: if progenitor divergence holds even with embryonic HSC from
a completely different study, the signal is unlikely to be a TMS protocol artifact.

Math
----
1. Query Census for all mouse HSCs, identify independent sources.
2. Download sci-RNA-seq3 HSC expression, filter to 16,959 ortholog gene space.
3. Normalize: counts per 10k + log1p (same as pipeline).
4. Compute new HSC centroid: μ_HSC^{indep} = (1/n) Σ x_i
5. Replace old HSC centroid in mouse centroid matrix.
6. Re-run full Procrustes pipeline (PCA → alignment → permutation → residuals).
7. Mann-Whitney U test: progenitor vs differentiated residuals.
"""

# WARNING: This script uses null_mean (arithmetic mean) instead of null_median
# for the obs/null ratio. This is intentional for diagnostic purposes.
# This script is NOT part of the canonical CellWarp pipeline.
# Do not use these obs/null values as canonical results.
# Canonical pipeline: src/procrustes.py (uses median throughout).

from __future__ import annotations

import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import anndata as ad
import cellxgene_census
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

from cellwarp.data_loader import (
    OBS_COLUMNS,
    RANDOM_SEED,
    VAR_COLUMNS,
    fetch_orthologs,
    save_h5ad_atomic,
)
from cellwarp.procrustes import (
    compute_centroids,
    compute_residual_vectors,
    pca_reduce_centroids,
    permutation_test,
    procrustes_align,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ORTHOLOGS_CACHE = PROJECT_ROOT / "data" / "phase1" / "orthologs_human_mouse.csv"
ORIGINAL_RESULTS_DIR = PROJECT_ROOT / "output" / "phase2" / "scaled_35types"
SENSITIVITY_DIR = PROJECT_ROOT / "output" / "phase2" / "sensitivity" / "smartseq2"
OUTPUT_DIR = PROJECT_ROOT / "output" / "phase2" / "sensitivity" / "hsc_validation"
DATA_DIR = PROJECT_ROOT / "data" / "phase1"

HSC_OUTPUT_PATH = DATA_DIR / "mouse_hsc_independent.h5ad"

# Mouse HSC cell type name in Census
HSC_CELL_TYPE = "hematopoietic stem cell"

# Assays considered "UMI-based 3' capture" (protocol-compatible with 10x)
UMI_3PRIME_ASSAYS = ["10x 3' v2", "10x 3' v3", "sci-RNA-seq3"]
TENX_ASSAYS = ["10x 3' v2", "10x 3' v3"]

MAX_HSC_CELLS = 2_000
MIN_HSC_CELLS = 200
N_PERMUTATIONS = 10_000

# Progenitor cell types (from developmental_annotations.csv)
PROGENITOR_TYPES = {
    "hematopoietic precursor cell",
    "hematopoietic stem cell",
    "basal cell",
    "mesenchymal stem cell of adipose tissue",
    "mesenchymal stem cell",
}


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # ==================================================================
    # STEP 1: Inventory query — find independent mouse HSC
    # ==================================================================
    print("=" * 70)
    print("  STEP 1: Census inventory — mouse HSC (non-TMS)")
    print("=" * 70)

    with cellxgene_census.open_soma(census_version="2025-11-08") as census:
        # Get TMS dataset IDs to exclude
        datasets_df = census["census_info"]["datasets"].read().concat().to_pandas()
        tms_mask = datasets_df["collection_name"].str.contains(
            "Tabula Muris", case=False, na=False
        )
        tms_dataset_ids = set(datasets_df.loc[tms_mask, "dataset_id"].tolist())
        print(f"\n  Tabula Muris Senis dataset IDs to exclude: {len(tms_dataset_ids)}")

        # Query all mouse HSC (assay names with apostrophes break SOMA parser,
        # so we filter post-hoc)
        value_filter = (
            f"cell_type == '{HSC_CELL_TYPE}' "
            f"and is_primary_data == True "
            f"and disease == 'normal'"
        )

        print(f"\n  Query filter: {value_filter}")
        print("  Querying Census obs metadata...")

        obs_df = cellxgene_census.get_obs(
            census,
            "Mus musculus",
            value_filter=value_filter,
            column_names=[
                "soma_joinid",
                "cell_type", "tissue", "tissue_general", "assay",
                "dataset_id", "donor_id", "is_primary_data", "disease",
            ],
        )

        print(f"\n  Total mouse HSC (primary, normal): {len(obs_df):,}")

        # Global assay breakdown
        print(f"\n  --- Assay breakdown (all datasets) ---")
        for assay, count in obs_df["assay"].value_counts().items():
            if count > 0:
                print(f"    {assay}: {count:,}")

        # 10x-only check
        is_10x = obs_df["assay"].isin(TENX_ASSAYS)
        obs_10x = obs_df[is_10x].copy()
        non_tms_10x = obs_10x[~obs_10x["dataset_id"].isin(tms_dataset_ids)]
        print(f"\n  10x Chromium HSC (total): {len(obs_10x):,}")
        print(f"  10x Chromium HSC (non-TMS): {len(non_tms_10x):,}")

        if len(non_tms_10x) == 0:
            print(f"\n  *** NO 10x mouse HSC data outside TMS ***")
            print(f"  All {len(obs_10x)} 10x HSC cells belong to Tabula Muris Senis.")
            print(f"  Falling back to sci-RNA-seq3 (UMI-based 3' capture, closest to 10x)")

        # Non-TMS breakdown
        obs_non_tms = obs_df[~obs_df["dataset_id"].isin(tms_dataset_ids)].copy()
        print(f"\n  Non-TMS HSC (all assays): {len(obs_non_tms):,}")

        print(f"\n  --- Non-TMS assay breakdown ---")
        for assay, count in obs_non_tms["assay"].value_counts().items():
            if count > 0:
                print(f"    {assay}: {count:,}")

        print(f"\n  --- Non-TMS tissue breakdown ---")
        for tissue, count in obs_non_tms["tissue"].value_counts().items():
            if count > 0:
                print(f"    {tissue}: {count:,}")

        # Non-TMS dataset breakdown
        print(f"\n  --- Non-TMS datasets ---")
        for ds_id, count in obs_non_tms["dataset_id"].value_counts().head(10).items():
            if count > 0:
                ds_row = datasets_df[datasets_df["dataset_id"] == ds_id]
                if len(ds_row) > 0:
                    coll = ds_row.iloc[0].get("collection_name", "unknown")[:60]
                    title = ds_row.iloc[0].get("dataset_title", "unknown")[:60]
                else:
                    coll = title = "unknown"
                print(f"    {count:>6,} | {coll}")
                print(f"           {title}")

        # Filter to UMI-based 3' assays from non-TMS datasets
        is_umi = obs_non_tms["assay"].isin(UMI_3PRIME_ASSAYS)
        obs_target = obs_non_tms[is_umi].copy()
        n_available = len(obs_target)

        print(f"\n  UMI-based 3' HSC (non-TMS): {n_available:,}")
        if n_available > 0:
            selected_assay = obs_target["assay"].value_counts().index[0]
            selected_tissue = obs_target["tissue"].value_counts().index[0]
            print(f"  Assay: {selected_assay}")
            print(f"  Tissue: {selected_tissue}")

        # Save inventory
        inventory = {
            "query": value_filter,
            "total_hsc": int(len(obs_df)),
            "assay_breakdown": {
                str(k): int(v)
                for k, v in obs_df["assay"].value_counts().items() if v > 0
            },
            "tenx_hsc_total": int(len(obs_10x)),
            "tenx_hsc_non_tms": int(len(non_tms_10x)),
            "non_tms_total": int(len(obs_non_tms)),
            "non_tms_umi_3prime": int(n_available),
            "conclusion": (
                "No 10x HSC outside TMS. Using sci-RNA-seq3 (UMI-based 3' capture) "
                "from Cao et al. embryonic development dataset as closest protocol match."
                if len(non_tms_10x) == 0 and n_available > 0
                else f"Found {len(non_tms_10x)} non-TMS 10x HSC cells"
            ),
        }
        with open(OUTPUT_DIR / "hsc_inventory.json", "w") as f:
            json.dump(inventory, f, indent=2, default=str)

        # ==================================================================
        # STEP 2: Download HSC expression data
        # ==================================================================
        if n_available < MIN_HSC_CELLS:
            print(f"\n  STOP: Only {n_available} non-TMS UMI-based HSC cells found.")
            print(f"  Need ≥{MIN_HSC_CELLS}. Cannot validate.")
            summary = {
                "status": "INSUFFICIENT_DATA",
                "cells_found": n_available,
                "minimum_required": MIN_HSC_CELLS,
            }
            with open(OUTPUT_DIR / "validation_results.json", "w") as f:
                json.dump(summary, f, indent=2)
            return

        print(f"\n{'=' * 70}")
        print(f"  STEP 2: Download HSC expression ({n_available:,} cells available)")
        print(f"{'=' * 70}")

        # Pre-subsample to MAX_HSC_CELLS before downloading expression
        # (downloading 64k cells is too slow; we only need 2000)
        if n_available > MAX_HSC_CELLS:
            rng = np.random.default_rng(RANDOM_SEED)
            sample_idx = rng.choice(n_available, size=MAX_HSC_CELLS, replace=False)
            obs_target = obs_target.iloc[sample_idx]
            print(f"  Pre-subsampled obs to {len(obs_target):,} cells")

        # Use soma_joinid to download only the pre-selected cells
        target_joinids = obs_target["soma_joinid"].tolist()

        print(f"  Downloading expression for {len(obs_target):,} cells "
              f"(soma_joinid range: {min(target_joinids)}-{max(target_joinids)})...")
        hsc_adata = cellxgene_census.get_anndata(
            census=census,
            organism="Mus musculus",
            obs_coords=target_joinids,
            obs_column_names=OBS_COLUMNS,
            var_column_names=VAR_COLUMNS,
        )
        print(f"  Downloaded: {hsc_adata.n_obs:,} cells × {hsc_adata.n_vars:,} genes")

        # Verify cell types
        if "cell_type" in hsc_adata.obs.columns:
            ct_counts = hsc_adata.obs["cell_type"].value_counts()
            for ct_name, ct_count in ct_counts.items():
                if ct_count > 0:
                    print(f"    {ct_name}: {ct_count:,}")

        # Post-hoc filter to UMI-based 3' assays and HSC cell type (safety)
        if "assay" in hsc_adata.obs.columns:
            is_valid = (
                hsc_adata.obs["assay"].isin(UMI_3PRIME_ASSAYS) &
                (hsc_adata.obs["cell_type"] == HSC_CELL_TYPE)
            )
            n_before = hsc_adata.n_obs
            hsc_adata = hsc_adata[is_valid].copy()
            if n_before != hsc_adata.n_obs:
                print(f"  After filter: {hsc_adata.n_obs:,} "
                      f"(removed {n_before - hsc_adata.n_obs})")

    # Filter to ortholog gene space
    print(f"\n  Filtering to ortholog gene space...")
    orthologs = fetch_orthologs(cache_path=ORTHOLOGS_CACHE)

    # Map mouse Ensembl IDs to human Ensembl IDs
    mouse_to_human = dict(
        zip(orthologs["mouse_ensembl_id"], orthologs["human_ensembl_id"])
    )
    mouse_to_human_name = dict(
        zip(orthologs["mouse_ensembl_id"], orthologs["human_gene_name"])
    )

    # Get the target gene space (human Ensembl IDs used in our pipeline)
    mouse_scaled = ad.read_h5ad(
        PROJECT_ROOT / "data" / "phase2_scaled" / "mouse_scaled.h5ad",
        backed="r",
    )
    target_genes = mouse_scaled.var_names.tolist()  # human Ensembl IDs
    target_genes_set = set(target_genes)
    mouse_scaled.file.close()

    # Map HSC gene IDs (mouse Ensembl) to human Ensembl IDs
    hsc_gene_ids = hsc_adata.var["feature_id"].values
    hsc_human_ids = [mouse_to_human.get(gid, None) for gid in hsc_gene_ids]

    # Find genes in our ortholog space
    keep_mask = np.array([
        hid is not None and hid in target_genes_set
        for hid in hsc_human_ids
    ])
    n_matched = int(keep_mask.sum())
    print(f"  HSC genes matching ortholog space: {n_matched:,} / {len(hsc_gene_ids):,}")

    hsc_filtered = hsc_adata[:, keep_mask].copy()

    # Re-index to human Ensembl IDs
    matched_human_ids = [hsc_human_ids[i] for i in range(len(hsc_human_ids)) if keep_mask[i]]
    hsc_filtered.var["original_mouse_feature_id"] = hsc_filtered.var["feature_id"].values
    hsc_filtered.var["original_mouse_feature_name"] = hsc_filtered.var["feature_name"].values
    hsc_filtered.var["human_ensembl_id"] = matched_human_ids
    hsc_filtered.var["human_gene_name"] = [
        mouse_to_human_name.get(mid, "")
        for mid in hsc_filtered.var["feature_id"].values
    ]
    hsc_filtered.var.index = matched_human_ids

    # Intersect with target gene space
    shared_genes = sorted(target_genes_set & set(matched_human_ids))
    n_missing = len(target_genes) - len(shared_genes)
    print(f"  Shared genes with pipeline: {len(shared_genes):,} / {len(target_genes):,}")
    print(f"  Missing genes: {n_missing:,}")

    hsc_filtered = hsc_filtered[:, shared_genes].copy()

    # Normalize: counts per 10k + log1p
    print(f"  Normalizing (counts per 10k + log1p)...")
    sc.pp.normalize_total(hsc_filtered, target_sum=10_000)
    sc.pp.log1p(hsc_filtered)

    # Save
    save_h5ad_atomic(hsc_filtered, HSC_OUTPUT_PATH)
    print(f"  Saved: {HSC_OUTPUT_PATH}")
    print(f"  Final shape: {hsc_filtered.shape}")

    # ==================================================================
    # STEP 3: Compute new HSC centroid and retest
    # ==================================================================
    print(f"\n{'=' * 70}")
    print(f"  STEP 3: Recompute HSC centroid and retest progenitor finding")
    print(f"{'=' * 70}")

    # Compute new HSC centroid
    new_hsc_centroid = np.asarray(hsc_filtered.X.mean(axis=0)).flatten()
    print(f"\n  New HSC centroid computed from {hsc_filtered.n_obs:,} independent cells")
    print(f"  Source: sci-RNA-seq3 (Cao et al., embryonic development)")
    print(f"  Centroid shape: {new_hsc_centroid.shape}")
    print(f"  Centroid mean: {new_hsc_centroid.mean():.6f}, max: {new_hsc_centroid.max():.4f}")

    # Load existing centroids
    human_centroids = pd.read_csv(
        ORIGINAL_RESULTS_DIR / "centroids_human_35.csv", index_col=0
    )
    mouse_centroids_orig = pd.read_csv(
        ORIGINAL_RESULTS_DIR / "centroids_mouse_35.csv", index_col=0
    )

    print(f"  Human centroids: {human_centroids.shape}")
    print(f"  Mouse centroids (original): {mouse_centroids_orig.shape}")

    # Build new HSC centroid series aligned to full gene space
    new_hsc_series = pd.Series(0.0, index=mouse_centroids_orig.columns)
    new_hsc_series[shared_genes] = new_hsc_centroid

    # --- Run A: Replace HSC in ORIGINAL (all-protocol) centroids ---
    print(f"\n  --- Analysis A: Replace HSC in original (all-protocol) centroids ---")

    mouse_centroids_a = mouse_centroids_orig.copy()
    mouse_centroids_a.loc[HSC_CELL_TYPE] = new_hsc_series

    shared_cols = sorted(set(human_centroids.columns) & set(mouse_centroids_a.columns))
    h_a = human_centroids[shared_cols]
    m_a = mouse_centroids_a[shared_cols]
    cts_a = sorted(set(h_a.index) & set(m_a.index))
    h_a = h_a.loc[cts_a]
    m_a = m_a.loc[cts_a]

    print(f"  Cell types: {len(cts_a)}, Genes: {len(shared_cols):,}")

    # PCA + Procrustes
    h_pca_a, m_pca_a, _, ct_list_a = pca_reduce_centroids(h_a, m_a)
    result_a = procrustes_align(h_pca_a, m_pca_a)
    p_a, null_a = permutation_test(h_pca_a, m_pca_a, n_permutations=N_PERMUTATIONS, seed=RANDOM_SEED)
    ratio_a = result_a.distance / np.mean(null_a)

    print(f"\n  Procrustes distance: {result_a.distance:.4f}")
    print(f"  p-value: {p_a:.6f}")
    print(f"  Obs/null ratio: {ratio_a:.4f}")

    # Residuals
    resid_a = compute_residual_vectors(result_a, ct_list_a)
    resid_df_a = _build_residual_df(resid_a, ct_list_a)

    hsc_row_a = resid_df_a[resid_df_a["cell_type"] == HSC_CELL_TYPE].iloc[0]
    print(f"\n  HSC residual: {hsc_row_a['residual_magnitude']:.4f} (rank {hsc_row_a['rank']})")

    # Mann-Whitney
    mw_p_a, prog_mean_a, diff_mean_a = _progenitor_test(resid_df_a, PROGENITOR_TYPES)
    finding_a = mw_p_a < 0.05
    print(f"\n  Progenitor mean: {prog_mean_a:.4f}, Diff mean: {diff_mean_a:.4f}")
    print(f"  Ratio: {prog_mean_a / diff_mean_a:.4f}")
    print(f"  Mann-Whitney p: {mw_p_a:.6f} → {'SIGNIFICANT' if finding_a else 'not significant'}")

    # --- Run B: Replace HSC in 10x-only centroids ---
    print(f"\n  --- Analysis B: Replace HSC in 10x-only mouse centroids ---")

    mouse_scaled_path = PROJECT_ROOT / "data" / "phase2_scaled" / "mouse_scaled.h5ad"
    mouse_all = ad.read_h5ad(mouse_scaled_path)
    mouse_10x = mouse_all[mouse_all.obs["assay"] == "10x 3' v2"].copy()
    del mouse_all

    print(f"\n  Computing 10x-only TMS centroids ({mouse_10x.n_obs:,} cells)...")
    mouse_10x_centroids = compute_centroids(mouse_10x, "cell_type")
    del mouse_10x

    # Replace HSC with independent centroid
    mouse_10x_centroids.loc[HSC_CELL_TYPE] = new_hsc_series[mouse_10x_centroids.columns]

    shared_cols_b = sorted(set(human_centroids.columns) & set(mouse_10x_centroids.columns))
    h_b = human_centroids[shared_cols_b]
    m_b = mouse_10x_centroids[shared_cols_b]
    cts_b = sorted(set(h_b.index) & set(m_b.index))
    h_b = h_b.loc[cts_b]
    m_b = m_b.loc[cts_b]

    print(f"  Cell types: {len(cts_b)}, Genes: {len(shared_cols_b):,}")

    h_pca_b, m_pca_b, _, ct_list_b = pca_reduce_centroids(h_b, m_b)
    result_b = procrustes_align(h_pca_b, m_pca_b)
    p_b, null_b = permutation_test(h_pca_b, m_pca_b, n_permutations=N_PERMUTATIONS, seed=RANDOM_SEED)
    ratio_b = result_b.distance / np.mean(null_b)

    print(f"\n  Procrustes distance: {result_b.distance:.4f}")
    print(f"  p-value: {p_b:.6f}")
    print(f"  Obs/null ratio: {ratio_b:.4f}")

    resid_b = compute_residual_vectors(result_b, ct_list_b)
    resid_df_b = _build_residual_df(resid_b, ct_list_b)

    hsc_row_b = resid_df_b[resid_df_b["cell_type"] == HSC_CELL_TYPE].iloc[0]
    print(f"\n  HSC residual: {hsc_row_b['residual_magnitude']:.4f} (rank {hsc_row_b['rank']})")

    mw_p_b, prog_mean_b, diff_mean_b = _progenitor_test(resid_df_b, PROGENITOR_TYPES)
    finding_b = mw_p_b < 0.05
    print(f"\n  Progenitor mean: {prog_mean_b:.4f}, Diff mean: {diff_mean_b:.4f}")
    print(f"  Ratio: {prog_mean_b / diff_mean_b:.4f}")
    print(f"  Mann-Whitney p: {mw_p_b:.6f} → {'SIGNIFICANT' if finding_b else 'not significant'}")

    # --- Run C: 10x-only, HSC excluded (4 progenitors) ---
    print(f"\n  --- Analysis C: 10x-only, HSC excluded (4 progenitors) ---")
    PROG_NO_HSC = PROGENITOR_TYPES - {HSC_CELL_TYPE}
    mw_p_c, prog_mean_c, diff_mean_c = _progenitor_test(
        resid_df_b, PROG_NO_HSC, exclude_from_diff={HSC_CELL_TYPE}
    )
    finding_c = mw_p_c < 0.05
    print(f"  Progenitor mean (excl HSC): {prog_mean_c:.4f}")
    print(f"  Differentiated mean: {diff_mean_c:.4f}")
    print(f"  Mann-Whitney p: {mw_p_c:.6f} → {'SIGNIFICANT' if finding_c else 'not significant'}")

    # ==================================================================
    # STEP 4: Save results and report
    # ==================================================================
    print(f"\n{'=' * 70}")
    print(f"  STEP 4: Save results and report")
    print(f"{'=' * 70}")

    # Load comparison p-values
    with open(SENSITIVITY_DIR / "sensitivity_results.json") as f:
        ss2 = json.load(f)
    orig_mw_p = ss2["progenitor_analysis"]["original"]["mannwhitney_p"]
    tenx_mw_p = ss2["progenitor_analysis"]["tenx_only"]["mannwhitney_p"]

    resid_df_a.to_csv(OUTPUT_DIR / "residuals_hsc_replaced.csv", index=False)
    resid_df_b.to_csv(OUTPUT_DIR / "residuals_10x_hsc_replaced.csv", index=False)

    summary = {
        "analysis": "HSC independent validation — progenitor finding retest",
        "description": (
            "No 10x mouse HSC data exists outside TMS. Used sci-RNA-seq3 "
            "(UMI-based 3' capture, Cao et al. embryonic development, 64k cells) "
            "as closest protocol match. Caveats: (1) sci-RNA-seq3 not 10x, "
            "(2) embryonic not adult HSCs."
        ),
        "hsc_data": {
            "source": "CZ CELLxGENE Census — Cao et al. embryonic development",
            "assay": "sci-RNA-seq3",
            "tissue": "embryo",
            "n_cells_available": int(n_available),
            "n_cells_used": int(hsc_filtered.n_obs),
            "n_genes_matched": len(shared_genes),
            "n_genes_pipeline": len(target_genes),
            "n_genes_missing": n_missing,
            "10x_non_tms_available": 0,
            "protocol_note": (
                "sci-RNA-seq3 is UMI-based 3'-end capture like 10x. "
                "Both use UMI deduplication and have comparable per-cell depth. "
                "Key difference from Smart-seq2: no full-length bias."
            ),
        },
        "analysis_a_hsc_in_original": {
            "description": "Replace HSC in original (all-protocol) mouse centroids",
            "procrustes_distance": float(result_a.distance),
            "p_value": float(p_a),
            "obs_null_ratio": float(ratio_a),
            "hsc_residual": float(hsc_row_a["residual_magnitude"]),
            "hsc_rank": int(hsc_row_a["rank"]),
            "progenitor_mean_residual": float(prog_mean_a),
            "differentiated_mean_residual": float(diff_mean_a),
            "ratio": float(prog_mean_a / diff_mean_a),
            "mannwhitney_p": float(mw_p_a),
            "significant": bool(finding_a),
        },
        "analysis_b_hsc_in_10x_only": {
            "description": "Replace HSC in 10x-only TMS centroids + independent HSC",
            "procrustes_distance": float(result_b.distance),
            "p_value": float(p_b),
            "obs_null_ratio": float(ratio_b),
            "hsc_residual": float(hsc_row_b["residual_magnitude"]),
            "hsc_rank": int(hsc_row_b["rank"]),
            "progenitor_mean_residual": float(prog_mean_b),
            "differentiated_mean_residual": float(diff_mean_b),
            "ratio": float(prog_mean_b / diff_mean_b),
            "mannwhitney_p": float(mw_p_b),
            "significant": bool(finding_b),
        },
        "analysis_c_excl_hsc": {
            "description": "10x-only centroids, HSC excluded (4 progenitors)",
            "progenitor_mean_residual": float(prog_mean_c),
            "differentiated_mean_residual": float(diff_mean_c),
            "mannwhitney_p": float(mw_p_c),
            "significant": bool(finding_c),
        },
        "p_value_progression": {
            "original_all_protocol": float(orig_mw_p),
            "tenx_only_tms_23_hsc": float(tenx_mw_p),
            "original_plus_indep_hsc": float(mw_p_a),
            "tenx_only_plus_indep_hsc": float(mw_p_b),
            "tenx_only_excl_hsc_4prog": float(mw_p_c),
        },
    }

    with open(OUTPUT_DIR / "validation_results.json", "w") as f:
        json.dump(summary, f, indent=2)
    np.save(OUTPUT_DIR / "null_distribution_a.npy", null_a)
    np.save(OUTPUT_DIR / "null_distribution_b.npy", null_b)

    # Final report
    print(f"\n{'=' * 70}")
    print(f"  FINAL REPORT — HSC Independent Validation")
    print(f"{'=' * 70}")

    print(f"\n  Data source: sci-RNA-seq3 embryonic HSC (Cao et al.)")
    print(f"  Cells used: {hsc_filtered.n_obs:,}")
    print(f"  Genes matched: {len(shared_genes):,} / {len(target_genes):,}")
    print(f"  Caveats: embryonic (not adult), sci-RNA-seq3 (not 10x)")

    print(f"\n  === p-value progression ===")
    print(f"  {'Analysis':<55} {'MW p':>10} {'Sig?':>6}")
    print(f"  {'-' * 73}")
    print(f"  {'Original (all protocol, TMS HSC)':<55} {orig_mw_p:>10.4f} {'*' if orig_mw_p < 0.05 else '':>6}")
    print(f"  {'10x-only TMS (HSC=23 cells, unreliable)':<55} {tenx_mw_p:>10.4f} {'*' if tenx_mw_p < 0.05 else '':>6}")
    print(f"  {'Original + independent embryonic HSC':<55} {mw_p_a:>10.4f} {'*' if finding_a else '':>6}")
    print(f"  {'10x-only + independent embryonic HSC':<55} {mw_p_b:>10.4f} {'*' if finding_b else '':>6}")
    print(f"  {'10x-only, HSC excluded (4 progenitors)':<55} {mw_p_c:>10.4f} {'*' if finding_c else '':>6}")

    print(f"\n  === Interpretation ===")
    if finding_b:
        print(f"  RESULT: Progenitor finding RECOVERS (p={mw_p_b:.4f}).")
        print(f"  The p=0.119 in 10x-only was driven by the unreliable 23-cell HSC centroid.")
        print(f"  With {hsc_filtered.n_obs:,} independent embryonic HSC: signal recovers.")
        print(f"  Conclusion: progenitor divergence is REAL, obscured by HSC protocol confound.")
        print(f"  Caveat: validation uses embryonic HSC (sci-RNA-seq3), not adult 10x.")
    else:
        print(f"  RESULT: Progenitor finding remains non-significant (p={mw_p_b:.4f}).")
        print(f"  Even with {hsc_filtered.n_obs:,} independent HSC cells, effect does not reach p<0.05.")
        print(f"  Conclusion: progenitor divergence is a trend, not a robust finding.")
        print(f"  Recommend framing as directional observation, not central claim.")

    print(f"\n  Saved to: {OUTPUT_DIR}/")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _build_residual_df(
    residuals: dict[str, np.ndarray], cell_types: list[str]
) -> pd.DataFrame:
    """Build ranked residual magnitude DataFrame."""
    data = []
    for ct in cell_types:
        data.append({
            "cell_type": ct,
            "residual_magnitude": float(np.linalg.norm(residuals[ct])),
        })
    df = pd.DataFrame(data).sort_values("residual_magnitude", ascending=False)
    df = df.reset_index(drop=True)
    df["rank"] = df.index + 1
    return df


def _progenitor_test(
    resid_df: pd.DataFrame,
    progenitor_types: set[str],
    exclude_from_diff: set[str] | None = None,
) -> tuple[float, float, float]:
    """Run Mann-Whitney U on progenitor vs differentiated residuals."""
    is_prog = resid_df["cell_type"].isin(progenitor_types)
    prog = resid_df.loc[is_prog, "residual_magnitude"].values

    if exclude_from_diff:
        excl = resid_df["cell_type"].isin(exclude_from_diff)
        diff = resid_df.loc[~is_prog & ~excl, "residual_magnitude"].values
    else:
        diff = resid_df.loc[~is_prog, "residual_magnitude"].values

    _, p = stats.mannwhitneyu(prog, diff, alternative="greater")
    return float(p), float(prog.mean()), float(diff.mean())


if __name__ == "__main__":
    main()
