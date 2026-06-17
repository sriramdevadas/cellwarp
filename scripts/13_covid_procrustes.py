#!/usr/bin/env python3
"""
CellWarp — COVID-19 Disease-Axis Procrustes Pipeline

Replicates the cross-axis rigidity correlation (evolutionary rigidity predicts
disease deformation resistance) using COVID-19 as the primary multi-disease
replication target. Cancer showed ρ=0.407 at n=13 (NS); COVID-19 at n=20 has
power to detect ρ≥0.45 at p<0.05.

Biology
-------
COVID-19 is a systemic disease affecting 12+ tissues. We compare COVID-19
patient cells to tissue-matched normal controls across 20 coherent cell types.
The cross-axis hypothesis predicts that cell types with low evolutionary
Procrustes residual (rigid in cross-species space) should also resist
disease-driven geometric deformation (low COVID deformation score).

Math
----
Identical to scripts/12_cancer_scaled.py:
1. Per-donor centroids → condition-level centroids (controls donor imbalance)
2. Joint PCA (95% variance threshold)
3. Procrustes alignment (normal → COVID)
4. 10,000-iteration permutation test
5. Per-cell-type residual magnitudes
6. Spearman rank correlation with cross-species residuals

Key difference from cancer: COVID cells span multiple tissues per cell type.
We require tissue-matched normal controls (CRITICAL design requirement).

Steps
-----
1. Load cross-species 35-type residuals
2. Download COVID-19 and tissue-matched normal cells from Census
3. Run Procrustes pipeline (centroids → PCA → alignment → permutation)
4. Cross-axis Spearman with sensitivity analyses
5. Severity metadata check (stretch goal — document only)
"""

from __future__ import annotations

import json
import signal
import sys
import tempfile
from collections import defaultdict
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from cellwarp.procrustes import (
    RANDOM_SEED,
    N_PERMUTATIONS,
    PCA_VARIANCE_THRESHOLD,
    N_TOP_GENES,
    procrustes_align,
    permutation_test,
    compute_residual_vectors,
    map_residuals_to_genes,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DATA_DIR = Path("data/disease_replication")
OUTPUT_DIR = Path("output/disease_replication/covid")
CROSS_SPECIES_RESULTS = Path("output/phase2/scaled_35types/procrustes_results_35.json")
RESIDUALS_RANKED = Path("output/phase2/scaled_35types/residuals_ranked.csv")
ORTHOLOG_CSV = Path("data/phase1/orthologs_human_mouse.csv")

# CELLxGENE disease labels
DISEASE_COVID = "COVID-19"
DISEASE_NORMAL = "normal"

# Cell count gates
MIN_CELLS_GATE = 500
MAX_CELLS_PER_TYPE = 3_000

# The 20 cell types to include (from inventory: all PASS)
INCLUDE_TYPES = {
    "stromal cell",
    "epithelial cell",
    "hematopoietic stem cell",
    "basal cell",
    "plasma cell",
    "CD4-positive, alpha-beta T cell",
    "classical monocyte",
    "macrophage",
    "B cell",
    "myeloid dendritic cell",
    "monocyte",
    "natural killer cell",
    "intermediate monocyte",
    "mature NK T cell",
    "fibroblast",
    "smooth muscle cell",
    "endothelial cell",
    "non-classical monocyte",
    "CD8-positive, alpha-beta T cell",
    "neutrophil",
}

# Explicitly excluded
EXCLUDE_TYPES = {
    "granulocyte",              # MARGINAL, unreliable
    "hematopoietic precursor cell",  # MARGINAL, biologically complex in COVID
}

# Census query timeout (seconds) — 600s for expression data downloads
# (COVID spans many tissues; individual downloads are faster than combined)
CENSUS_QUERY_TIMEOUT = 600

# COVID-to-cross-species type mapping for Spearman correlation
# All are exact 1:1 matches to 35-type set except mature NKT cell
COVID_TO_XS_MAP = {
    "stromal cell": (["stromal cell"], "exact"),
    "epithelial cell": (["epithelial cell"], "exact"),
    "hematopoietic stem cell": (["hematopoietic stem cell"], "exact"),
    "basal cell": (["basal cell"], "exact"),
    "plasma cell": (["plasma cell"], "exact"),
    "CD4-positive, alpha-beta T cell": (["CD4-positive, alpha-beta T cell"], "exact"),
    "classical monocyte": (["classical monocyte"], "exact"),
    "macrophage": (["macrophage"], "exact"),
    "B cell": (["B cell"], "exact"),
    "myeloid dendritic cell": (["myeloid dendritic cell"], "exact"),
    "monocyte": (["monocyte"], "exact"),
    "natural killer cell": (["natural killer cell"], "exact"),
    "intermediate monocyte": (["intermediate monocyte"], "exact"),
    "mature NK T cell": (["natural killer cell", "CD8-positive, alpha-beta T cell"], "ambiguous"),
    "fibroblast": (["fibroblast"], "exact"),
    "smooth muscle cell": (["smooth muscle cell"], "exact"),
    "endothelial cell": (["endothelial cell"], "exact"),
    "non-classical monocyte": (["non-classical monocyte"], "exact"),
    "CD8-positive, alpha-beta T cell": (["CD8-positive, alpha-beta T cell"], "exact"),
    "neutrophil": (["neutrophil"], "exact"),
}


@contextmanager
def census_timeout(seconds: int = CENSUS_QUERY_TIMEOUT):
    """Raise TimeoutError if a Census API call exceeds *seconds*."""
    def _handler(signum, frame):
        raise TimeoutError(f"Census query timed out after {seconds}s")
    prev = signal.signal(signal.SIGALRM, _handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, prev)


# ---------------------------------------------------------------------------
# Coarse mapping: Census fine labels → 35-type set
# ---------------------------------------------------------------------------
# Reuses the mapping from the COVID inventory script (disease_inventory_covid.py).
# Census labels for COVID span many sub-types; we map them to the 20 target types.

COARSE_MAP = {
    # Direct matches
    "B cell": "B cell",
    "CD4-positive, alpha-beta T cell": "CD4-positive, alpha-beta T cell",
    "CD8-positive, alpha-beta T cell": "CD8-positive, alpha-beta T cell",
    "classical monocyte": "classical monocyte",
    "non-classical monocyte": "non-classical monocyte",
    "intermediate monocyte": "intermediate monocyte",
    "monocyte": "monocyte",
    "macrophage": "macrophage",
    "natural killer cell": "natural killer cell",
    "mature NK T cell": "mature NK T cell",
    "plasma cell": "plasma cell",
    "endothelial cell": "endothelial cell",
    "fibroblast": "fibroblast",
    "epithelial cell": "epithelial cell",
    "myeloid dendritic cell": "myeloid dendritic cell",
    "neutrophil": "neutrophil",
    "hematopoietic stem cell": "hematopoietic stem cell",
    "smooth muscle cell": "smooth muscle cell",
    "stromal cell": "stromal cell",
    "basal cell": "basal cell",
    # Coarsening — biologically justified
    "memory B cell": "B cell",
    "naive B cell": "B cell",
    "plasmablast": "plasma cell",
    "CD14-positive, CD16-negative classical monocyte": "classical monocyte",
    "CD16-positive, CD56-dim natural killer cell": "natural killer cell",
    "CD56-bright natural killer cell": "natural killer cell",
    "effector memory CD4-positive, alpha-beta T cell": "CD4-positive, alpha-beta T cell",
    "central memory CD4-positive, alpha-beta T cell": "CD4-positive, alpha-beta T cell",
    "naive thymus-derived CD4-positive, alpha-beta T cell": "CD4-positive, alpha-beta T cell",
    "regulatory T cell": "CD4-positive, alpha-beta T cell",
    "effector memory CD8-positive, alpha-beta T cell": "CD8-positive, alpha-beta T cell",
    "central memory CD8-positive, alpha-beta T cell": "CD8-positive, alpha-beta T cell",
    "naive thymus-derived CD8-positive, alpha-beta T cell": "CD8-positive, alpha-beta T cell",
    "conventional dendritic cell": "myeloid dendritic cell",
    "dendritic cell": "myeloid dendritic cell",
    "alveolar macrophage": "macrophage",
    "lung macrophage": "macrophage",
    "Kupffer cell": "macrophage",
    "lung endothelial cell": "endothelial cell",
    "endothelial cell of artery": "endothelial cell",
    "endothelial cell of vascular tree": "endothelial cell",
    "capillary endothelial cell": "endothelial cell",
    "vein endothelial cell": "endothelial cell",
    "blood vessel endothelial cell": "endothelial cell",
    "glomerular endothelial cell": "endothelial cell",
    "endothelial cell of lymphatic vessel": "endothelial cell",
    "epithelial cell of lung": "epithelial cell",
    "type I pneumocyte": "epithelial cell",
    "type II pneumocyte": "epithelial cell",
    "club cell": "epithelial cell",
    "ciliated cell": "epithelial cell",
    "goblet cell": "epithelial cell",
    "respiratory goblet cell": "epithelial cell",
    "bronchial epithelial cell": "epithelial cell",
    "respiratory epithelial cell": "epithelial cell",
    "nasal epithelial cell": "epithelial cell",
    "mucus secreting cell": "epithelial cell",
    "secretory cell": "epithelial cell",
    "lung ciliated cell": "epithelial cell",
    "pulmonary ionocyte": "epithelial cell",
    "vascular associated smooth muscle cell": "smooth muscle cell",
    "bronchial smooth muscle cell": "smooth muscle cell",
    "myofibroblast cell": "fibroblast",
    "lung fibroblast": "fibroblast",
    "fibroblast of lung": "fibroblast",
}


def classify_coarse(raw_label: str) -> str | None:
    """Map a raw Census cell type label to a 20-type target.
    Returns None if unmapped or excluded."""
    target = COARSE_MAP.get(raw_label)
    if target is None:
        return None
    if target not in INCLUDE_TYPES:
        return None
    return target


# ===================================================================
# STEP 1 — Load cross-species residuals
# ===================================================================


def load_cross_species_residuals() -> dict[str, float]:
    """Load the 35-type cross-species residual magnitudes."""
    print("=" * 70)
    print("STEP 1 — Load 35-Type Cross-Species Residuals")
    print("=" * 70)

    residuals_df = pd.read_csv(RESIDUALS_RANKED)
    xs_residuals = dict(zip(residuals_df["cell_type"], residuals_df["residual_magnitude"]))

    print(f"\n  35 cross-species cell types loaded.")
    print(f"  20 target COVID types and their XS residuals:")
    for ct in sorted(INCLUDE_TYPES):
        mag = xs_residuals.get(ct, None)
        if mag is not None:
            print(f"    {ct:<50} {mag:.4f}")
        else:
            print(f"    {ct:<50} NOT IN 35-TYPE SET")

    return xs_residuals


# ===================================================================
# STEP 2 — Download COVID and tissue-matched normal data
# ===================================================================


def download_covid_data(xs_residuals: dict[str, float]):
    """
    Download COVID-19 and tissue-matched normal cells from CELLxGENE Census.

    CRITICAL: For each cell type, verifies that COVID and normal cells come
    from the SAME tissues. Prints tissue source for both conditions.

    Returns: covid_adata, normal_adata
    """
    import anndata as ad
    import cellxgene_census
    import scanpy as sc
    from cellwarp.cancer_loader import load_ortholog_gene_ids, save_h5ad_atomic

    print("\n" + "=" * 70)
    print("STEP 2 — Download COVID-19 and Tissue-Matched Normal Data")
    print("=" * 70)

    # Load ortholog gene IDs
    ortholog_ids = load_ortholog_gene_ids(ORTHOLOG_CSV)
    print(f"  Ortholog gene space: {len(ortholog_ids):,} genes")

    OBS_COLUMNS = [
        "cell_type", "tissue", "tissue_general", "disease",
        "donor_id", "dataset_id", "assay", "sex", "development_stage",
        "is_primary_data",
    ]
    VAR_COLUMNS = ["feature_id", "feature_name"]

    rng = np.random.default_rng(RANDOM_SEED)

    print(f"\n  Opening Census (2025-11-08)...")
    census = cellxgene_census.open_soma(census_version="2025-11-08")

    # -------------------------------------------------------------------
    # Phase 1: Fetch metadata for COVID to determine tissues per type
    # -------------------------------------------------------------------
    print(f"\n  [Phase 1] Fetching COVID metadata to determine tissue sources...")
    with census_timeout(300):
        obs_covid = cellxgene_census.get_obs(
            census, "Homo sapiens",
            value_filter=(
                "is_primary_data == True and "
                f"disease == '{DISEASE_COVID}'"
            ),
            column_names=list(dict.fromkeys(["soma_joinid"] + OBS_COLUMNS)),
        )
    print(f"    Total COVID cells: {len(obs_covid):,}")

    # Apply coarse mapping
    obs_covid["target_type"] = obs_covid["cell_type"].map(
        lambda x: classify_coarse(x)
    )
    obs_covid = obs_covid[obs_covid["target_type"].notna()].copy()
    print(f"    After mapping to 20 target types: {len(obs_covid):,} cells")

    # Build tissue sources per target type
    # Filter to tissues with ≥50 COVID cells to avoid noise and keep queries fast
    type_tissues: dict[str, list[str]] = {}
    print(f"\n  TISSUE SOURCES PER CELL TYPE (COVID, ≥50 cells/tissue):")
    print(f"  {'Cell Type':<45} {'Tissues'}")
    print(f"  {'-' * 90}")
    for ct in sorted(INCLUDE_TYPES):
        ct_mask = obs_covid["target_type"] == ct
        if ct_mask.sum() == 0:
            print(f"  {ct:<45} NO COVID DATA")
            continue
        tissues = obs_covid.loc[ct_mask, "tissue"].value_counts()
        # Filter: only include tissues with real counts (not categorical zeros)
        tissues = tissues[tissues > 0]
        tissue_list = [t for t, c in tissues.items() if c >= 50]
        if not tissue_list:
            # Fall back to ≥10 if no tissue has ≥50
            tissue_list = [t for t, c in tissues.items() if c >= 10]
        type_tissues[ct] = tissue_list
        tissue_str = ", ".join(f"{t}({c:,})" for t, c in tissues.head(5).items() if c > 0)
        n_tissues = len(tissue_list)
        print(f"  {ct:<45} {n_tissues} tissues: {tissue_str}")

    # -------------------------------------------------------------------
    # Phase 2: Download per (raw_type, tissue) with value_filter
    # -------------------------------------------------------------------
    # Strategy: for each target type, iterate raw labels × tissues from
    # smallest to largest. Download each (raw_type, tissue) combo using
    # value_filter (fast indexed lookup). Accumulate until we have 3,000
    # cells per condition. This avoids downloading 579K cells and avoids
    # the scattered-joinid problem that makes obs_coords slow.
    print(f"\n  [Phase 2] Downloading per (raw_type, tissue) — smallest first...")
    print(f"  (≤{MAX_CELLS_PER_TYPE:,} cells/type, tissue-matched normals)")

    covid_chunks = []
    normal_chunks = []
    tissue_match_report = []

    for ct in sorted(INCLUDE_TYPES):
        if ct not in type_tissues or len(type_tissues[ct]) == 0:
            print(f"\n  SKIP {ct}: no COVID tissue data")
            continue

        tissues = type_tissues[ct]
        raw_labels = sorted([k for k, v in COARSE_MAP.items() if v == ct])

        # Build per-(raw, tissue) cell counts from metadata
        ct_meta = obs_covid[obs_covid["target_type"] == ct]
        tissue_raw_counts = []
        for tissue in tissues:
            for raw_label in raw_labels:
                n = ((ct_meta["tissue"] == tissue) &
                     (ct_meta["cell_type"] == raw_label)).sum()
                if n > 0:
                    tissue_raw_counts.append((tissue, raw_label, n))

        # Sort smallest first — download manageable chunks first
        tissue_raw_counts.sort(key=lambda x: x[2])

        print(f"\n  {ct} ({len(raw_labels)} raw, {len(tissues)} tissues, "
              f"{len(tissue_raw_counts)} combos):")

        # Download COVID cells
        ct_covid_chunks = []
        ct_covid_total = 0
        ct_covid_tissues = set()
        for tissue, raw_label, expected_n in tissue_raw_counts:
            if ct_covid_total >= MAX_CELLS_PER_TYPE:
                break
            raw_escaped = raw_label.replace("'", "\\'")
            vf = (
                f"is_primary_data == True and "
                f"disease == '{DISEASE_COVID}' and "
                f"tissue == '{tissue}' and "
                f"cell_type == '{raw_escaped}'"
            )
            print(f"    COVID {raw_label[:30]:30} in {tissue[:25]:25} (~{expected_n:,})...",
                  end="", flush=True)
            try:
                with census_timeout():
                    chunk = cellxgene_census.get_anndata(
                        census=census,
                        organism="Homo sapiens",
                        obs_value_filter=vf,
                        obs_column_names=OBS_COLUMNS,
                        var_column_names=VAR_COLUMNS,
                    )
                ct_covid_chunks.append(chunk)
                ct_covid_total += chunk.n_obs
                ct_covid_tissues.add(tissue)
                print(f" {chunk.n_obs:,}", flush=True)
            except TimeoutError:
                print(f" TIMEOUT", flush=True)
                try:
                    census.close()
                except Exception:
                    pass
                census = cellxgene_census.open_soma(census_version="2025-11-08")

        if not ct_covid_chunks:
            print(f"    → No COVID data for {ct}")
            continue

        # Download normal cells from SAME tissues
        ct_normal_chunks = []
        ct_normal_total = 0
        ct_normal_tissues = set()
        for tissue in sorted(ct_covid_tissues):
            if ct_normal_total >= MAX_CELLS_PER_TYPE:
                break
            for raw_label in raw_labels:
                if ct_normal_total >= MAX_CELLS_PER_TYPE:
                    break
                raw_escaped = raw_label.replace("'", "\\'")
                vf = (
                    f"is_primary_data == True and "
                    f"disease == '{DISEASE_NORMAL}' and "
                    f"tissue == '{tissue}' and "
                    f"cell_type == '{raw_escaped}'"
                )
                print(f"    Normal {raw_label[:30]:30} in {tissue[:25]:25}...",
                      end="", flush=True)
                try:
                    with census_timeout():
                        chunk = cellxgene_census.get_anndata(
                            census=census,
                            organism="Homo sapiens",
                            obs_value_filter=vf,
                            obs_column_names=OBS_COLUMNS,
                            var_column_names=VAR_COLUMNS,
                        )
                    if chunk.n_obs > 0:
                        ct_normal_chunks.append(chunk)
                        ct_normal_total += chunk.n_obs
                        ct_normal_tissues.add(tissue)
                    print(f" {chunk.n_obs:,}", flush=True)
                except TimeoutError:
                    print(f" TIMEOUT", flush=True)
                    try:
                        census.close()
                    except Exception:
                        pass
                    census = cellxgene_census.open_soma(census_version="2025-11-08")

        if not ct_normal_chunks:
            print(f"    → No normal data for {ct}")
            continue

        # Concatenate
        covid_chunk = (ad.concat(ct_covid_chunks, merge="same")
                       if len(ct_covid_chunks) > 1 else ct_covid_chunks[0])
        normal_chunk = (ad.concat(ct_normal_chunks, merge="same")
                        if len(ct_normal_chunks) > 1 else ct_normal_chunks[0])
        del ct_covid_chunks, ct_normal_chunks

        # Subsample to MAX_CELLS_PER_TYPE
        if covid_chunk.n_obs > MAX_CELLS_PER_TYPE:
            idx = rng.choice(covid_chunk.n_obs, size=MAX_CELLS_PER_TYPE, replace=False)
            idx.sort()
            covid_chunk = covid_chunk[idx].copy()
        if normal_chunk.n_obs > MAX_CELLS_PER_TYPE:
            idx = rng.choice(normal_chunk.n_obs, size=MAX_CELLS_PER_TYPE, replace=False)
            idx.sort()
            normal_chunk = normal_chunk[idx].copy()

        overlap = ct_covid_tissues & ct_normal_tissues
        tissue_match_report.append({
            "cell_type": ct,
            "covid_tissues": sorted(ct_covid_tissues),
            "normal_tissues": sorted(ct_normal_tissues),
            "overlap_tissues": sorted(overlap),
            "fully_matched": len(ct_covid_tissues - overlap) == 0,
            "covid_cells": covid_chunk.n_obs,
            "normal_cells": normal_chunk.n_obs,
        })

        print(f"    → COVID={covid_chunk.n_obs:,}, Normal={normal_chunk.n_obs:,}")
        print(f"    Tissue match: {sorted(overlap)}")

        covid_chunk.obs["coarse_cell_type"] = ct
        normal_chunk.obs["coarse_cell_type"] = ct
        covid_chunks.append(covid_chunk)
        normal_chunks.append(normal_chunk)

    census.close()

    # -------------------------------------------------------------------
    # Phase 3: Concatenate, filter to ortholog space, normalize
    # -------------------------------------------------------------------
    print(f"\n  [Phase 3] Concatenating and normalizing...")

    covid_adata = ad.concat(covid_chunks, merge="same")
    normal_adata = ad.concat(normal_chunks, merge="same")
    del covid_chunks, normal_chunks

    print(f"    COVID combined: {covid_adata.n_obs:,} cells × {covid_adata.n_vars:,} genes")
    print(f"    Normal combined: {normal_adata.n_obs:,} cells × {normal_adata.n_vars:,} genes")

    # Filter to ortholog space
    for label, adata in [("COVID", covid_adata), ("Normal", normal_adata)]:
        gene_mask = adata.var["feature_id"].isin(ortholog_ids)
        adata_filt = adata[:, gene_mask].copy()
        adata_filt.var.index = adata_filt.var["feature_id"].values
        shared_ids = sorted(set(adata_filt.var.index) & set(ortholog_ids))
        adata_filt = adata_filt[:, shared_ids].copy()
        # QC
        n_before = adata_filt.n_obs
        sc.pp.filter_cells(adata_filt, min_genes=200)
        n_after = adata_filt.n_obs
        if n_before != n_after:
            print(f"    {label} QC: {n_before:,} → {n_after:,} ({n_before - n_after:,} removed)")
        # Normalize
        sc.pp.normalize_total(adata_filt, target_sum=10_000)
        sc.pp.log1p(adata_filt)
        print(f"    {label}: {adata_filt.n_obs:,} cells × {adata_filt.n_vars:,} genes (normalized)")
        if label == "COVID":
            covid_adata = adata_filt
        else:
            normal_adata = adata_filt

    # -------------------------------------------------------------------
    # Print tissue-match verification table
    # -------------------------------------------------------------------
    print(f"\n  TISSUE-MATCH VERIFICATION TABLE:")
    print(f"  {'Cell Type':<40} {'COVID Cells':>10} {'Normal Cells':>12} "
          f"{'COVID Tissues':>40} {'Normal Tissues':>40} {'Match':>6}")
    print(f"  {'-' * 155}")
    for tm in tissue_match_report:
        ct_str = ", ".join(tm["covid_tissues"][:3])
        nt_str = ", ".join(tm["normal_tissues"][:3])
        if len(tm["covid_tissues"]) > 3:
            ct_str += "..."
        if len(tm["normal_tissues"]) > 3:
            nt_str += "..."
        match = "YES" if tm["fully_matched"] else "PARTIAL"
        print(f"  {tm['cell_type']:<40} {tm['covid_cells']:>10,} {tm['normal_cells']:>12,} "
              f"{ct_str:>40} {nt_str:>40} {match:>6}")

    # Save with atomic writes
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    from cellwarp.cancer_loader import save_h5ad_atomic
    save_h5ad_atomic(covid_adata, DATA_DIR / "covid_condition.h5ad")
    save_h5ad_atomic(normal_adata, DATA_DIR / "covid_normal.h5ad")

    # Save tissue match report
    pd.DataFrame(tissue_match_report).to_csv(
        OUTPUT_DIR / "tissue_match_report.csv", index=False
    )

    # Print checkpoint summary
    print(f"\n  DOWNLOAD SUMMARY:")
    print(f"  {'Cell Type':<40} {'COVID':>8} {'Normal':>8}")
    print(f"  {'-' * 60}")
    all_types = sorted(set(
        list(covid_adata.obs["coarse_cell_type"].unique()) +
        list(normal_adata.obs["coarse_cell_type"].unique())
    ))
    for ct in all_types:
        nc = (covid_adata.obs["coarse_cell_type"] == ct).sum()
        nn = (normal_adata.obs["coarse_cell_type"] == ct).sum()
        print(f"  {ct:<40} {nc:>8,} {nn:>8,}")
    print(f"  {'-' * 60}")
    print(f"  {'TOTAL':<40} {covid_adata.n_obs:>8,} {normal_adata.n_obs:>8,}")
    print(f"  Genes: {covid_adata.n_vars:,}")
    print(f"  COVID donors: {covid_adata.obs['donor_id'].nunique()}")
    print(f"  Normal donors: {normal_adata.obs['donor_id'].nunique()}")
    print(f"\n  CHECKPOINT: Download complete.")

    return covid_adata, normal_adata


# ===================================================================
# STEP 3 — Procrustes pipeline
# ===================================================================


def _donor_centroids(
    adata, cell_types: list[str], gene_names: list[str], label: str
) -> pd.DataFrame:
    """Compute mean-of-donor-means centroids.

    Biology: Averaging per-donor first controls for donor count imbalance.
    A donor with 2,000 cells contributes the same as one with 50 cells.

    Math: μ_{t,c} = (1/D) Σ_d (1/n_d) Σ_i x_i
    """
    centroids = {}
    for ct in cell_types:
        mask = adata.obs["coarse_cell_type"] == ct
        ct_data = adata[mask]
        donors = ct_data.obs["donor_id"].unique()
        donor_means = []
        for d in donors:
            d_mask = ct_data.obs["donor_id"] == d
            d_cells = ct_data[d_mask]
            if d_cells.n_obs > 0:
                donor_means.append(np.asarray(d_cells.X.mean(axis=0)).flatten())
        centroids[ct] = np.mean(donor_means, axis=0)
        print(f"    {label} {ct:<40} {mask.sum():>6,} cells, {len(donors):>4} donors")
    df = pd.DataFrame(centroids, index=gene_names).T
    df.index.name = "cell_type"
    return df


def run_procrustes(covid_adata, normal_adata) -> tuple:
    """
    Run the full Procrustes pipeline on COVID vs normal data.

    Returns: (cell_types, scores, residuals, result, p_value, null_dist,
              pca_model, gene_names, top_genes, dropped)
    """
    import anndata as ad

    print("\n" + "=" * 70)
    print("STEP 3 — COVID Procrustes Pipeline")
    print("=" * 70)

    shared_genes = list(normal_adata.var.index)

    # Identify valid cell types (≥500 in BOTH conditions)
    covid_counts = covid_adata.obs["coarse_cell_type"].value_counts()
    normal_counts = normal_adata.obs["coarse_cell_type"].value_counts()

    all_types = sorted(
        set(covid_counts.index) & set(normal_counts.index)
    )
    passed = []
    dropped = []
    for ct in all_types:
        n_c = covid_counts.get(ct, 0)
        n_n = normal_counts.get(ct, 0)
        if n_c >= MIN_CELLS_GATE and n_n >= MIN_CELLS_GATE:
            passed.append(ct)
        else:
            dropped.append(ct)

    print(f"\n  Cell types passing ≥{MIN_CELLS_GATE} gate ({len(passed)}):")
    for ct in passed:
        print(f"    {ct:<40} covid={covid_counts[ct]:>6,}  normal={normal_counts[ct]:>6,}")
    if dropped:
        print(f"  Dropped ({len(dropped)}):")
        for ct in dropped:
            c = covid_counts.get(ct, 0)
            n = normal_counts.get(ct, 0)
            print(f"    {ct:<40} covid={c:>6,}  normal={n:>6,}")

    cell_types = passed

    # --- Per-donor centroids ---
    print(f"\n  Computing per-donor centroids...")
    normal_centroids = _donor_centroids(normal_adata, cell_types, shared_genes, "Normal")
    covid_centroids = _donor_centroids(covid_adata, cell_types, shared_genes, "COVID")

    # Save centroids
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    normal_centroids.to_csv(OUTPUT_DIR / "centroids_normal.csv")
    covid_centroids.to_csv(OUTPUT_DIR / "centroids_covid.csv")

    # --- PCA ---
    print(f"\n  Running PCA on combined centroids...")
    normal_mat = normal_centroids.loc[cell_types].values
    covid_mat = covid_centroids.loc[cell_types].values
    combined = np.vstack([normal_mat, covid_mat])

    pca = PCA(
        n_components=PCA_VARIANCE_THRESHOLD,
        svd_solver="full",
        random_state=RANDOM_SEED,
    )
    combined_pca = pca.fit_transform(combined)
    n_types = len(cell_types)
    normal_pca = combined_pca[:n_types]
    covid_pca = combined_pca[n_types:]

    cumvar = np.cumsum(pca.explained_variance_ratio_)
    print(f"    PCA: {pca.n_components_} components, {cumvar[-1]*100:.1f}% variance")

    np.savez(
        OUTPUT_DIR / "pca_covid.npz",
        normal_pca=normal_pca, covid_pca=covid_pca,
        components=pca.components_,
        explained_variance_ratio=pca.explained_variance_ratio_,
        mean=pca.mean_,
    )

    # --- Procrustes alignment (normal → COVID) ---
    print(f"\n  Running Procrustes alignment (normal → COVID)...")
    result = procrustes_align(normal_pca, covid_pca)
    det = np.linalg.det(result.rotation)
    assert abs(det - 1.0) < 1e-6, f"Rotation det={det}, expected +1.0"
    print(f"    Distance: {result.distance:.4f}, Scaling: {result.scaling:.6f}")

    # --- Permutation test ---
    print(f"\n  Running permutation test ({N_PERMUTATIONS:,} iterations)...")
    p_value, null_dist = permutation_test(
        normal_pca, covid_pca, n_permutations=N_PERMUTATIONS, seed=RANDOM_SEED
    )
    np.save(OUTPUT_DIR / "null_distribution.npy", null_dist)
    null_median = float(np.median(null_dist))
    obs_null = result.distance / null_median
    print(f"    p={p_value:.6f}, obs/null={obs_null:.3f}")

    # --- Deformation scores ---
    print(f"\n  Computing deformation scores...")
    residuals = compute_residual_vectors(result, cell_types)
    scores = {ct: float(np.linalg.norm(residuals[ct])) for ct in cell_types}

    ranked = sorted(scores.items(), key=lambda x: x[1], reverse=True)
    print(f"\n  DEFORMATION RANKING:")
    print(f"  {'Rank':<6} {'Cell Type':<40} {'Score':>10}")
    print(f"  {'-' * 58}")
    for i, (ct, score) in enumerate(ranked, 1):
        print(f"  {i:<6} {ct:<40} {score:>10.4f}")

    # --- Top genes ---
    ensembl_to_name = dict(zip(
        covid_adata.var["feature_id"], covid_adata.var["feature_name"]
    ))
    readable_genes = [ensembl_to_name.get(g, g) for g in shared_genes]
    top_genes = map_residuals_to_genes(residuals, pca, readable_genes, n_top=N_TOP_GENES)

    return (cell_types, scores, residuals, result, p_value, null_dist,
            pca, shared_genes, top_genes, dropped)


# ===================================================================
# STEP 4 — Cross-axis Spearman correlation + scatter plot
# ===================================================================


def cross_analysis_spearman(
    covid_scores: dict[str, float],
    xs_residuals: dict[str, float],
) -> dict:
    """
    Match COVID cell types to cross-species 35-type residuals. Compute
    Spearman correlation with multiple sensitivity analyses.

    Biology: Tests whether cell types that are evolutionarily rigid
    (low cross-species Procrustes residual) also resist disease-driven
    deformation (low COVID deformation score).

    Math: Spearman rank correlation between COVID deformation scores
    and cross-species residual magnitudes across n=20 matched types.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    print("\n" + "=" * 70)
    print("STEP 4 — Cross-Analysis Spearman Correlation")
    print("=" * 70)

    matched = []
    print(f"\n  FULL MATCH TABLE:")
    print(
        f"  {'COVID Type':<42} {'Deformation':>11} "
        f"{'XS Rigidity':>11} {'XS Match':<42} {'Quality':>8} {'N DS':>5}"
    )
    print(f"  {'-' * 125}")

    for covid_type in sorted(covid_scores.keys()):
        if covid_type not in COVID_TO_XS_MAP:
            print(
                f"  {covid_type:<42} {covid_scores[covid_type]:>11.4f} "
                f"{'':>11} {'[no mapping]':<42} {'SKIP':>8}"
            )
            continue

        xs_list, quality = COVID_TO_XS_MAP[covid_type]

        # Check all xs types exist
        missing = [x for x in xs_list if x not in xs_residuals]
        if missing:
            print(
                f"  {covid_type:<42} {covid_scores[covid_type]:>11.4f} "
                f"{'':>11} {str(xs_list):<42} {'MISSING':>8}"
            )
            continue

        xs_mag = float(np.mean([xs_residuals[x] for x in xs_list]))
        xs_label = " + ".join(xs_list) if len(xs_list) > 1 else xs_list[0]

        matched.append({
            "covid_type": covid_type,
            "xs_type": xs_label,
            "covid_deformation": covid_scores[covid_type],
            "xs_residual": xs_mag,
            "match_quality": quality,
        })
        print(
            f"  {covid_type:<42} {covid_scores[covid_type]:>11.4f} "
            f"{xs_mag:>11.4f} {xs_label:<42} {quality:>8}"
        )

    n_matched = len(matched)
    print(f"\n  n = {n_matched} matched types")

    if n_matched < 3:
        print(f"  WARNING: Too few matches for correlation.")
        return {"n_matched": n_matched, "insufficient": True}

    # --- Full Spearman (a) ---
    covid_vals = [m["covid_deformation"] for m in matched]
    xs_vals = [m["xs_residual"] for m in matched]
    rho, p_value = stats.spearmanr(covid_vals, xs_vals)

    print(f"\n  (a) FULL CORRELATION (n={n_matched}):")
    print(f"      Spearman ρ = {rho:.4f}")
    print(f"      p-value   = {p_value:.6f}")
    print(f"      Significant at α=0.05: {'YES' if p_value < 0.05 else 'NO'}")
    print(f"      Significant at α=0.01: {'YES' if p_value < 0.01 else 'NO'}")

    # --- Without HSC (b) ---
    matched_no_hsc = [m for m in matched if m["covid_type"] != "hematopoietic stem cell"]
    rho_no_hsc, p_no_hsc = None, None
    if len(matched_no_hsc) >= 3:
        c_nh = [m["covid_deformation"] for m in matched_no_hsc]
        x_nh = [m["xs_residual"] for m in matched_no_hsc]
        rho_no_hsc, p_no_hsc = stats.spearmanr(c_nh, x_nh)
        print(f"\n  (b) WITHOUT HSC (protocol sensitivity, n={len(matched_no_hsc)}):")
        print(f"      Spearman ρ = {rho_no_hsc:.4f}")
        print(f"      p-value   = {p_no_hsc:.6f}")
        print(f"      Significant at α=0.05: {'YES' if p_no_hsc < 0.05 else 'NO'}")
        delta_rho = abs(rho_no_hsc - rho)
        print(f"      |Δρ| = {delta_rho:.4f} ({'sensitive' if delta_rho > 0.15 else 'robust'})")

    # --- Without ambiguous matches (NKT) (c) ---
    matched_no_ambig = [m for m in matched if m["match_quality"] != "ambiguous"]
    rho_no_ambig, p_no_ambig = None, None
    if len(matched_no_ambig) >= 3:
        c_na = [m["covid_deformation"] for m in matched_no_ambig]
        x_na = [m["xs_residual"] for m in matched_no_ambig]
        rho_no_ambig, p_no_ambig = stats.spearmanr(c_na, x_na)
        print(f"\n  (c) WITHOUT AMBIGUOUS MATCHES (n={len(matched_no_ambig)}):")
        print(f"      Spearman ρ = {rho_no_ambig:.4f}")
        print(f"      p-value   = {p_no_ambig:.6f}")
        print(f"      Significant at α=0.05: {'YES' if p_no_ambig < 0.05 else 'NO'}")
        delta_rho = abs(rho_no_ambig - rho)
        print(f"      |Δρ| = {delta_rho:.4f} ({'sensitive' if delta_rho > 0.15 else 'robust'})")

    # --- Without both HSC and ambiguous ---
    matched_strict = [
        m for m in matched
        if m["covid_type"] != "hematopoietic stem cell"
        and m["match_quality"] != "ambiguous"
    ]
    rho_strict, p_strict = None, None
    if len(matched_strict) >= 3:
        c_s = [m["covid_deformation"] for m in matched_strict]
        x_s = [m["xs_residual"] for m in matched_strict]
        rho_strict, p_strict = stats.spearmanr(c_s, x_s)
        print(f"\n  (d) WITHOUT HSC AND AMBIGUOUS (strictest, n={len(matched_strict)}):")
        print(f"      Spearman ρ = {rho_strict:.4f}")
        print(f"      p-value   = {p_strict:.6f}")
        print(f"      Significant at α=0.05: {'YES' if p_strict < 0.05 else 'NO'}")

    # --- Comparison with cancer ---
    print(f"\n  COMPARISON WITH CANCER RESULT:")
    print(f"    Cancer:  ρ = 0.407, p = 0.168, n = 13 (NOT SIGNIFICANT)")
    print(f"    COVID:   ρ = {rho:.4f}, p = {p_value:.6f}, n = {n_matched}")
    if p_value < 0.05:
        print(f"    → COVID ACHIEVES SIGNIFICANCE where cancer did not.")
    else:
        print(f"    → COVID also NOT SIGNIFICANT. Cross-axis correlation may not exist.")

    # --- Interpretation ---
    print(f"\n  INTERPRETATION:")
    if p_value < 0.05 and rho > 0:
        print(
            "  POSITIVE CORRELATION: Evolutionarily flexible cell types are also\n"
            "  more deformed by COVID-19. Shared plasticity axis confirmed."
        )
    elif p_value < 0.05 and rho < 0:
        print(
            "  NEGATIVE CORRELATION: Evolutionarily rigid cell types resist\n"
            "  COVID-driven deformation. Shared constraint on rewiring."
        )
    else:
        print(
            "  NO SIGNIFICANT CORRELATION: Evolutionary rigidity and COVID-19\n"
            "  deformation appear independent at this sample size."
        )

    # --- Scatter plot ---
    fig, ax = plt.subplots(figsize=(10, 8))
    colors = []
    for m in matched:
        if m["match_quality"] == "ambiguous":
            colors.append("orange")
        elif m["covid_type"] == "hematopoietic stem cell":
            colors.append("red")
        else:
            colors.append("steelblue")

    ax.scatter(xs_vals, covid_vals, s=60, alpha=0.8,
               c=colors, edgecolors="black", linewidths=0.5)

    for m in matched:
        label = m["covid_type"]
        if len(label) > 28:
            label = label[:25] + "..."
        ax.annotate(
            label,
            (m["xs_residual"], m["covid_deformation"]),
            fontsize=6.5,
            xytext=(5, 5),
            textcoords="offset points",
        )

    # Regression line
    z = np.polyfit(xs_vals, covid_vals, 1)
    x_line = np.linspace(min(xs_vals) * 0.95, max(xs_vals) * 1.05, 100)
    ax.plot(x_line, np.polyval(z, x_line), "--", color="red", alpha=0.5, linewidth=1)

    ax.set_xlabel("Cross-Species Residual Magnitude (35-type)", fontsize=11)
    ax.set_ylabel("COVID-19 Deformation Score (COVID vs Normal)", fontsize=11)
    ax.set_title(
        f"Evolutionary Rigidity vs COVID-19 Deformation\n"
        f"Spearman ρ={rho:.3f}, p={p_value:.4f}, n={n_matched}",
        fontsize=12,
    )
    ax.grid(True, alpha=0.3)

    # Legend for colors
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='steelblue',
               markersize=8, label='Exact match'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='orange',
               markersize=8, label='Ambiguous (NKT)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='red',
               markersize=8, label='HSC (flagged)'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=9)

    fig.tight_layout()
    plot_path = OUTPUT_DIR / "covid_cross_analysis.png"
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)
    print(f"\n  Saved scatter plot: {plot_path}")
    print(f"  Plot: {n_matched} cell types. X=cross-species residual, Y=COVID deformation.")
    print(f"  Blue=exact match, orange=ambiguous (NKT), red=HSC (flagged).")
    if rho > 0:
        print(f"  Points trend upward — flexible types deform more in both contexts.")
    elif rho < 0:
        print(f"  Points trend downward — rigid types resist deformation in both contexts.")

    # Save match table
    match_df = pd.DataFrame(matched)
    match_df.to_csv(OUTPUT_DIR / "cross_analysis_match_table.csv", index=False)
    print(f"  Saved match table: {OUTPUT_DIR / 'cross_analysis_match_table.csv'}")

    sensitivity = {}
    if rho_no_hsc is not None:
        sensitivity["without_HSC"] = {
            "rho": float(rho_no_hsc), "p": float(p_no_hsc),
            "n": len(matched_no_hsc),
        }
    if rho_no_ambig is not None:
        sensitivity["without_ambiguous"] = {
            "rho": float(rho_no_ambig), "p": float(p_no_ambig),
            "n": len(matched_no_ambig),
        }
    if rho_strict is not None:
        sensitivity["without_HSC_and_ambiguous"] = {
            "rho": float(rho_strict), "p": float(p_strict),
            "n": len(matched_strict),
        }

    return {
        "n_matched": n_matched,
        "matched_types": matched,
        "spearman_rho": float(rho),
        "spearman_p": float(p_value),
        "significant_005": bool(p_value < 0.05),
        "significant_001": bool(p_value < 0.01),
        "sensitivity": sensitivity,
        "cancer_comparison": {
            "cancer_rho": 0.407,
            "cancer_p": 0.168,
            "cancer_n": 13,
        },
    }


# ===================================================================
# STEP 5 — Severity metadata check (stretch goal)
# ===================================================================


def severity_metadata_check():
    """
    Check whether COVID-19 datasets in Census have severity stratification.
    This is a documentation-only step — no implementation.
    """
    import cellxgene_census

    print("\n" + "=" * 70)
    print("STEP 5 — Severity Metadata Check (Stretch Goal)")
    print("=" * 70)

    print(f"\n  Opening Census to query dataset metadata...")
    census = cellxgene_census.open_soma(census_version="2025-11-08")

    # Get COVID dataset IDs and their cell counts
    with census_timeout(120):
        obs_covid = cellxgene_census.get_obs(
            census, "Homo sapiens",
            value_filter=(
                f"is_primary_data == True and "
                f"disease == '{DISEASE_COVID}'"
            ),
            column_names=["dataset_id", "disease"],
        )

    ds_counts = obs_covid["dataset_id"].value_counts()
    print(f"\n  {len(ds_counts)} COVID datasets found.")

    # Check disease labels for severity encoding
    all_disease_labels = obs_covid["disease"].value_counts()
    print(f"\n  Disease labels in COVID data:")
    for label, count in all_disease_labels.items():
        if count > 0:
            print(f"    {label}: {count:,} cells")

    # Try to query Census datasets table for collection info
    print(f"\n  Top 10 COVID datasets by cell count:")
    print(f"  {'Rank':>4} {'Dataset ID':<40} {'Cells':>10}")
    print(f"  {'-' * 58}")
    for i, (did, count) in enumerate(ds_counts.head(10).items()):
        print(f"  {i+1:>4} {did:<40} {count:>10,}")

    census.close()

    severity_keywords = ["mild", "moderate", "severe", "critical",
                         "hospitalized", "ventilat", "asymptomatic", "ICU"]
    has_severity = False
    for label, count in all_disease_labels.items():
        if count > 0:
            for kw in severity_keywords:
                if kw.lower() in label.lower():
                    print(f"\n  FOUND severity: '{label}' contains '{kw}' ({count:,})")
                    has_severity = True

    if not has_severity:
        print(f"\n  No severity info in Census cell-level metadata.")
        print(f"  Severity may exist in original study publications but would")
        print(f"  require per-dataset metadata lookup outside Census API.")

    print(f"\n  NOTE: Severity stratification is a stretch goal for future work.")
    print(f"  Could enable mild-vs-severe deformation gradient analysis.")
    print(f"  Would require querying individual dataset publications.")

    return {
        "severity_in_census": has_severity,
        "n_datasets": len(ds_counts),
        "note": "Severity not in Census cell-level metadata. Future work.",
    }


# ===================================================================
# Save results
# ===================================================================


def save_results(
    result, p_value, null_dist, scores, residuals, top_genes,
    cell_types, pca_model, cross_corr, dropped, severity_info,
):
    """Save comprehensive results JSON."""
    null_median = float(np.median(null_dist))

    results_dict = {
        "procrustes": {
            "distance": float(result.distance),
            "distance_squared": float(result.distance_squared),
            "scaling": float(result.scaling),
            "rotation_det": float(np.linalg.det(result.rotation)),
        },
        "permutation_test": {
            "p_value": float(p_value),
            "n_permutations": len(null_dist),
            "null_median": null_median,
            "obs_null_ratio": float(result.distance / null_median),
            "significant_at_001": bool(p_value < 0.01),
        },
        "pca": {
            "n_components": int(pca_model.n_components_),
            "cumulative_variance": float(np.sum(pca_model.explained_variance_ratio_)),
            "per_component_variance": pca_model.explained_variance_ratio_.tolist(),
        },
        "cell_types": cell_types,
        "dropped_types": dropped,
        "deformation_scores": {ct: float(scores[ct]) for ct in cell_types},
        "deformation_ranking": [
            {"rank": i + 1, "cell_type": ct, "score": float(scores[ct])}
            for i, (ct, _) in enumerate(
                sorted(scores.items(), key=lambda x: x[1], reverse=True)
            )
        ],
        "residuals": {
            ct: {
                "vector_pca": residuals[ct].tolist(),
                "magnitude": float(np.linalg.norm(residuals[ct])),
            }
            for ct in cell_types
        },
        "top_genes_per_cell_type": {
            ct: top_genes[ct][["gene", "loading", "abs_loading", "rank"]]
            .to_dict(orient="records")
            for ct in cell_types
        },
        "cross_analysis_correlation": cross_corr,
        "severity_metadata": severity_info,
        "random_seed": RANDOM_SEED,
    }

    output_path = OUTPUT_DIR / "covid_procrustes_results.json"
    tmp_path = output_path.with_suffix(".tmp")
    with open(tmp_path, "w") as f:
        json.dump(results_dict, f, indent=2)
    tmp_path.rename(output_path)
    print(f"\n  Saved: {output_path}")


# ===================================================================
# Summary
# ===================================================================


def print_summary(
    result, p_value, null_dist, scores, cell_types,
    cross_corr, dropped, pca_model,
):
    """Print analyst-format summary."""
    null_median = float(np.median(null_dist))
    obs_null = result.distance / null_median

    print("\n" + "=" * 70)
    print("COVID-19 PROCRUSTES — SUMMARY")
    print("=" * 70)

    print(f"""
1. WHAT WAS DONE
   COVID-19 Procrustes with {len(cell_types)} tissue-matched cell types.
   Per-donor centroids averaged to control for donor imbalance.
   Dropped: {dropped if dropped else 'none'}.
   Gene space: {pca_model.mean_.shape[0]:,} ortholog genes.

2. KEY NUMBERS
   Procrustes distance:     {result.distance:.4f}
   Scaling factor:          {result.scaling:.6f}
   Permutation p-value:     {p_value:.6f}
   Null median:             {null_median:.4f}
   Obs/null ratio:          {obs_null:.3f}
   PCA components:          {pca_model.n_components_}
   Significant at α=0.01:   {'YES' if p_value < 0.01 else 'NO'}

3. DEFORMATION RANKING""")

    ranked = sorted(scores.items(), key=lambda x: x[1], reverse=True)
    print(f"   {'Rank':<6} {'Cell Type':<40} {'Score':>10}")
    print(f"   {'-' * 58}")
    for i, (ct, score) in enumerate(ranked, 1):
        print(f"   {i:<6} {ct:<40} {score:>10.4f}")

    cc = cross_corr
    print(f"""
4. CROSS-ANALYSIS (CRITICAL)
   Matched types:   {cc.get('n_matched', 'N/A')}
   Spearman ρ:      {cc.get('spearman_rho', 'N/A'):.4f}
   p-value:         {cc.get('spearman_p', 'N/A'):.6f}
   Significant:     {cc.get('significant_005', 'N/A')}

5. COMPARISON: CANCER vs COVID
   Cancer:  ρ = 0.407, p = 0.168, n = 13 (NOT SIGNIFICANT)
   COVID:   ρ = {cc.get('spearman_rho', '?'):.4f}, p = {cc.get('spearman_p', '?'):.6f}, n = {cc.get('n_matched', '?')}""")

    sens = cc.get("sensitivity", {})
    if sens:
        print(f"\n6. SENSITIVITY ANALYSIS")
        if "without_HSC" in sens:
            s = sens["without_HSC"]
            print(f"   Without HSC (protocol caveat):     ρ={s['rho']:.4f}, p={s['p']:.6f}, n={s['n']}")
        if "without_ambiguous" in sens:
            s = sens["without_ambiguous"]
            print(f"   Without ambiguous (NKT):           ρ={s['rho']:.4f}, p={s['p']:.6f}, n={s['n']}")
        if "without_HSC_and_ambiguous" in sens:
            s = sens["without_HSC_and_ambiguous"]
            print(f"   Without HSC + ambiguous (strict):  ρ={s['rho']:.4f}, p={s['p']:.6f}, n={s['n']}")

    if cc.get("spearman_p") is not None:
        p = cc["spearman_p"]
        rho = cc["spearman_rho"]
        if p < 0.05 and rho > 0:
            print(
                "\n   POSITIVE CORRELATION CONFIRMED with adequate power.\n"
                "   Evolutionary flexibility predicts COVID-19 vulnerability."
            )
        elif p < 0.05 and rho < 0:
            print(
                "\n   NEGATIVE CORRELATION CONFIRMED with adequate power.\n"
                "   Evolutionary rigidity predicts COVID-19 resistance."
            )
        elif p >= 0.05:
            print(
                "\n   NOT SIGNIFICANT even with n=20.\n"
                "   Evolutionary rigidity and disease deformation are independent."
            )


# ===================================================================
# Main
# ===================================================================


def main():
    print("\n" + "#" * 70)
    print("# CellWarp — COVID-19 Disease-Axis Procrustes Pipeline")
    print("# Primary replication target for cross-axis rigidity correlation")
    print("#" * 70 + "\n")

    # Step 1 — Load cross-species residuals
    xs_residuals = load_cross_species_residuals()

    # Step 2 — Download data
    covid_adata, normal_adata = download_covid_data(xs_residuals)

    print(f"\n  {'='*70}")
    print(f"  CHECKPOINT: Download complete. Proceeding to Procrustes pipeline.")
    print(f"  {'='*70}\n")

    # Step 3 — Procrustes pipeline
    (cell_types, scores, residuals, result, p_value, null_dist,
     pca_model, gene_names, top_genes, dropped) = run_procrustes(
        covid_adata, normal_adata
    )
    del covid_adata, normal_adata  # free memory

    # Step 4 — Cross-axis Spearman
    cross_corr = cross_analysis_spearman(scores, xs_residuals)

    # Step 5 — Severity metadata check
    severity_info = severity_metadata_check()

    # Save results
    save_results(
        result, p_value, null_dist, scores, residuals, top_genes,
        cell_types, pca_model, cross_corr, dropped, severity_info,
    )

    # Summary
    print_summary(
        result, p_value, null_dist, scores, cell_types,
        cross_corr, dropped, pca_model,
    )

    print(f"\n  All outputs saved to {OUTPUT_DIR}/")
    print("  Done.\n")


if __name__ == "__main__":
    main()
