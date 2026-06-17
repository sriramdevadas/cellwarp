#!/usr/bin/env python3
"""
CellWarp — Phase 2 Human-vs-Human Negative Control

Downloads a second independent human dataset from CZ CELLxGENE Census
(excluding Tabula Sapiens), processes it identically, and runs Procrustes
alignment against our existing Tabula Sapiens data.

Biology
-------
If aligning two independent HUMAN atlases produces a Procrustes distance
and p-value similar to our human-vs-mouse result, then our observed signal
is batch effects (atlas-level technical differences), not biology (species
differences). This is the critical negative control for Phase 2.

Expected outcome if our cross-species signal is real:
  - Human-vs-mouse:  small distance, significant p (already observed: p=0.0035)
  - Human-vs-human:  distance large relative to null, NON-significant p

Pipeline
--------
1. Explore Census metadata → find 2-3 non-TS datasets covering all 6 types
2. Download expression from those datasets only (efficient contiguous access)
3. Subsample to ≤2,000 cells per type
4. Filter to same 16,959 ortholog genes, normalize identically
5. Compute centroids, PCA, Procrustes (Human-Atlas-2 → Tabula Sapiens)
6. Permutation test (10,000 iterations)
7. Compare with human-vs-mouse results

Inputs:
    data/phase1/human_qc.h5ad               — Tabula Sapiens reference
    output/phase2/centroids_human.csv        — TS centroids (precomputed)
    output/phase2/procrustes_results.json    — Human-vs-mouse results

Outputs:
    data/phase2/human2_negctrl.h5ad                      — Second human atlas
    output/phase2/negative_control/centroids_human2.csv   — H2 centroids
    output/phase2/negative_control/negctrl_results.json   — Procrustes results
    output/phase2/negative_control/source_datasets.json   — Data provenance
    output/phase2/negative_control/comparison.csv         — Side-by-side table

Usage:
    python scripts/04c_negative_control.py
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from cellwarp.data_loader import save_h5ad_atomic
from cellwarp.procrustes import (
    RANDOM_SEED,
    compute_centroids,
    compute_residual_vectors,
    map_residuals_to_genes,
    pca_reduce_centroids,
    permutation_test,
    procrustes_align,
    save_results,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CELL_TYPES = [
    "B cell",
    "CD4-positive, alpha-beta T cell",
    "CD8-positive, alpha-beta T cell",
    "endothelial cell",
    "hepatocyte",
    "macrophage",
]
TS_COLLECTION_PATTERN = "Tabula Sapiens"
MAX_CELLS_PER_TYPE = 2_000
MIN_CELLS_FOR_COVERAGE = 100  # Min cells to count as "covered" in dataset
OUTPUT_DIR = Path("./output/phase2/negative_control")
DATA_DIR = Path("./data/phase2")


# ---------------------------------------------------------------------------
# Dataset selection and download
# ---------------------------------------------------------------------------


def find_covering_datasets(
    census,
    datasets_df: pd.DataFrame,
    ts_dataset_ids: set[str],
) -> tuple[list[str], pd.DataFrame, dict]:
    """
    Find a minimal set of non-TS datasets that collectively cover all 6 types.

    Uses a greedy set-cover algorithm: at each step, pick the dataset that
    covers the most currently-uncovered cell types (with ≥MIN_CELLS_FOR_COVERAGE
    cells). Typically finds a cover with 2-3 datasets.

    Returns:
        Tuple of (selected_dataset_ids, counts_matrix, source_info_dict).
    """
    import cellxgene_census

    print("\n  Querying Census metadata for non-TS cells per dataset...")

    # Collect per-dataset per-type cell counts
    all_obs = []
    for ct in CELL_TYPES:
        t0 = time.time()
        obs = cellxgene_census.get_obs(
            census,
            "Homo sapiens",
            value_filter=(
                f"cell_type == '{ct}' "
                f"and is_primary_data == True "
                f"and disease == 'normal'"
            ),
            column_names=["cell_type", "dataset_id"],
        )
        obs = obs[~obs["dataset_id"].isin(ts_dataset_ids)]
        all_obs.append(obs)
        dt = time.time() - t0
        print(
            f"    {ct:<45} {len(obs):>8,} non-TS cells "
            f"({obs['dataset_id'].nunique():>4} datasets) [{dt:.1f}s]"
        )

    combined = pd.concat(all_obs, ignore_index=True)

    # Build dataset × cell_type count matrix
    counts = (
        combined.groupby(["dataset_id", "cell_type"], observed=True)
        .size()
        .unstack(fill_value=0)
    )

    # Greedy set cover
    selected = []
    covered = set()
    remaining = set(CELL_TYPES)
    source_info = {}

    print("\n  Greedy dataset selection:")
    while remaining:
        best_ds = None
        best_new = set()
        best_total = 0

        for ds_id in counts.index:
            if ds_id in selected:
                continue
            new_types = set()
            for ct in remaining:
                if ct in counts.columns and counts.loc[ds_id, ct] >= MIN_CELLS_FOR_COVERAGE:
                    new_types.add(ct)
            total = counts.loc[ds_id].sum()
            if len(new_types) > len(best_new) or (
                len(new_types) == len(best_new) and total > best_total
            ):
                best_ds = ds_id
                best_new = new_types
                best_total = total

        if best_ds is None or len(best_new) == 0:
            print(f"  WARNING: Cannot cover types: {remaining}")
            break

        selected.append(best_ds)
        covered.update(best_new)
        remaining -= best_new

        row = datasets_df[datasets_df["dataset_id"] == best_ds]
        coll = row["collection_name"].iloc[0] if len(row) > 0 else "unknown"

        types_in_ds = {
            ct: int(counts.loc[best_ds, ct])
            for ct in CELL_TYPES
            if ct in counts.columns and counts.loc[best_ds, ct] > 0
        }
        source_info[best_ds] = {
            "collection": coll,
            "types_covered": list(best_new),
            "cell_counts": types_in_ds,
        }

        print(f"    Dataset {len(selected)}: {coll[:60]}")
        print(f"      ID: {best_ds}")
        print(f"      New types: {sorted(best_new)}")
        print(f"      Cell counts: {types_in_ds}")

    print(f"\n  Coverage: {len(covered)}/6 types from {len(selected)} datasets")
    return selected, counts, source_info


def download_second_human_atlas(
    target_gene_ids: set[str],
    target_gene_order: list[str],
) -> ad.AnnData:
    """
    Download a second independent human dataset from CELLxGENE Census.

    Strategy: find 2-3 non-TS datasets covering all 6 types via greedy
    set cover, then download ALL cells of our types from those datasets
    (efficient contiguous block access), then subsample to ≤2,000/type.
    """
    import cellxgene_census

    with cellxgene_census.open_soma(census_version="2025-11-08") as census:
        # Get datasets table and TS IDs
        datasets_df = (
            census["census_info"]["datasets"].read().concat().to_pandas()
        )
        ts_mask = datasets_df["collection_name"].str.contains(
            TS_COLLECTION_PATTERN, case=False, na=False
        )
        ts_dataset_ids = set(datasets_df[ts_mask]["dataset_id"].tolist())
        print(f"  Excluding {len(ts_dataset_ids)} Tabula Sapiens datasets")

        # Find minimal covering datasets
        selected_ds, counts, source_info = find_covering_datasets(
            census, datasets_df, ts_dataset_ids
        )

        # Save source info
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        with open(OUTPUT_DIR / "source_datasets.json", "w") as f:
            json.dump(source_info, f, indent=2, default=str)

        # Build filter for selected datasets
        names_str = ", ".join(f"'{ct}'" for ct in CELL_TYPES)
        ds_str = ", ".join(f"'{d}'" for d in selected_ds)
        obs_filter = (
            f"cell_type in [{names_str}] "
            f"and dataset_id in [{ds_str}] "
            f"and is_primary_data == True "
            f"and disease == 'normal'"
        )

        # Download expression
        print(f"\n  Downloading expression from {len(selected_ds)} datasets...")
        t0 = time.time()
        adata = cellxgene_census.get_anndata(
            census,
            "Homo sapiens",
            obs_value_filter=obs_filter,
            obs_column_names=[
                "cell_type",
                "tissue",
                "assay",
                "dataset_id",
                "donor_id",
                "is_primary_data",
                "disease",
            ],
            var_column_names=["feature_id", "feature_name"],
        )
        dt = time.time() - t0
        print(
            f"  Downloaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes "
            f"[{dt:.1f}s]"
        )

    # Per-type cell counts before subsampling
    print(f"\n  Pre-subsample cell counts:")
    for ct in CELL_TYPES:
        n = (adata.obs["cell_type"] == ct).sum()
        print(f"    {ct:<45} {n:>8,}")

    # Subsample to MAX_CELLS_PER_TYPE
    print(f"\n  Subsampling to ≤{MAX_CELLS_PER_TYPE:,} per type...")
    rng = np.random.default_rng(RANDOM_SEED)
    keep_indices = []
    for ct in sorted(adata.obs["cell_type"].unique()):
        ct_idx = np.where(adata.obs["cell_type"] == ct)[0]
        if len(ct_idx) > MAX_CELLS_PER_TYPE:
            selected = rng.choice(
                ct_idx, size=MAX_CELLS_PER_TYPE, replace=False
            )
            selected.sort()
            keep_indices.extend(selected.tolist())
            print(f"    {ct:<45} {len(ct_idx):>8,} → {MAX_CELLS_PER_TYPE:,}")
        else:
            keep_indices.extend(ct_idx.tolist())
            print(f"    {ct:<45} {len(ct_idx):>8,} (kept all)")

    adata = adata[keep_indices].copy()

    # Filter to our target genes
    print(f"\n  Filtering to {len(target_gene_ids):,} shared ortholog genes...")
    gene_mask = adata.var["feature_id"].isin(target_gene_ids)
    n_matched = gene_mask.sum()
    print(f"  Matched {n_matched:,} / {len(target_gene_ids):,} target genes")

    adata = adata[:, gene_mask].copy()
    adata.var.index = adata.var["feature_id"].values

    # Align to target gene order (genes present in both)
    available = set(adata.var_names)
    common_order = [g for g in target_gene_order if g in available]
    if len(common_order) < len(target_gene_order):
        n_miss = len(target_gene_order) - len(common_order)
        print(f"  WARNING: {n_miss} target genes missing ({len(common_order):,} retained)")
    adata = adata[:, common_order].copy()

    # Normalize identically: counts per 10k + log1p
    print(f"\n  Normalizing: counts per 10k + log1p...")
    adata.raw = adata
    sc.pp.normalize_total(adata, target_sum=10_000)
    sc.pp.log1p(adata)

    # Summary
    print(f"\n  Second human atlas summary:")
    print(f"  {'Cell Type':<45} {'Cells':>8} {'Donors':>8} {'Tissues':>8}")
    print(f"  {'-' * 72}")
    for ct in sorted(adata.obs["cell_type"].unique()):
        mask = adata.obs["cell_type"] == ct
        n = mask.sum()
        nd = adata.obs.loc[mask, "donor_id"].nunique()
        nt = adata.obs.loc[mask, "tissue"].nunique()
        print(f"  {ct:<45} {n:>8,} {nd:>8} {nt:>8}")
    print(f"  {'-' * 72}")
    print(f"  {'TOTAL':<45} {adata.n_obs:>8,}")
    print(f"  Genes: {adata.n_vars:,}")

    return adata


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    t_start = time.time()
    checkpoint_path = DATA_DIR / "human2_negctrl.h5ad"

    print("=" * 70)
    print("PHASE 2 — Human-vs-Human Negative Control")
    print("=" * 70)

    # ==================================================================
    # Load reference data
    # ==================================================================
    print("\n  Loading Tabula Sapiens reference...")
    human_ts = ad.read_h5ad("./data/phase1/human_qc.h5ad")
    target_gene_ids = set(human_ts.var_names)
    target_gene_order = list(human_ts.var_names)
    gene_symbols = human_ts.var["feature_name"].tolist()
    print(
        f"  TS reference: {human_ts.n_obs:,} cells × {human_ts.n_vars:,} genes"
    )

    # Load precomputed TS centroids
    ts_centroids = pd.read_csv(
        "./output/phase2/centroids_human.csv", index_col=0
    )

    # ==================================================================
    # Step 1: Download second human atlas (or load checkpoint)
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 1: Download second human atlas")
    print("=" * 70)

    if checkpoint_path.exists():
        print(f"  Checkpoint found: {checkpoint_path}")
        human2 = ad.read_h5ad(checkpoint_path)
        print(f"  Loaded: {human2.n_obs:,} cells × {human2.n_vars:,} genes")
    else:
        human2 = download_second_human_atlas(target_gene_ids, target_gene_order)
        save_h5ad_atomic(human2, checkpoint_path)

    # Verify gene alignment
    h2_genes = set(human2.var_names)
    ts_genes = set(human_ts.var_names)
    shared = h2_genes & ts_genes
    print(f"\n  Gene alignment: {len(shared):,} shared genes")

    if len(shared) < len(ts_genes):
        n_missing = len(ts_genes) - len(shared)
        print(f"  Restricting to {len(shared):,} shared genes")
        shared_order = [g for g in target_gene_order if g in shared]
        ts_centroids = ts_centroids[shared_order]
        gene_symbols = [
            s for g, s in zip(target_gene_order, gene_symbols) if g in shared
        ]
        human2 = human2[:, shared_order].copy()

    del human_ts  # Free memory

    # ==================================================================
    # Step 2: Compute Human2 centroids
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 2: Compute centroids for second human atlas")
    print("=" * 70)

    h2_centroids = compute_centroids(human2, "cell_type")
    h2_centroids.to_csv(OUTPUT_DIR / "centroids_human2.csv")
    del human2

    # ==================================================================
    # Step 3: PCA on combined TS + Human2 centroids
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 3: PCA on combined TS + Human2 centroids")
    print("=" * 70)

    h1_pca, h2_pca, pca_model, cell_types = pca_reduce_centroids(
        ts_centroids, h2_centroids
    )

    # ==================================================================
    # Step 4: Procrustes alignment (Human2 → TS)
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 4: Procrustes alignment (Human2 → Tabula Sapiens)")
    print("=" * 70)

    result = procrustes_align(h1_pca, h2_pca)

    # ==================================================================
    # Step 5: Permutation test
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 5: Permutation test (10,000 iterations)")
    print("=" * 70)

    p_value, null_dist = permutation_test(h1_pca, h2_pca)

    np.save(OUTPUT_DIR / "null_distribution_negctrl.npy", null_dist)

    # ==================================================================
    # Step 6: Residuals and gene mapping
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 6: Residual deformation vectors")
    print("=" * 70)

    residuals = compute_residual_vectors(result, cell_types)
    top_genes = map_residuals_to_genes(residuals, pca_model, gene_symbols)

    # ==================================================================
    # Step 7: Save results
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 7: Save results")
    print("=" * 70)

    pca_info = {
        "n_components": int(pca_model.n_components_),
        "variance_explained_per_component": (
            pca_model.explained_variance_ratio_.tolist()
        ),
        "cumulative_variance_explained": float(
            sum(pca_model.explained_variance_ratio_)
        ),
        "n_genes_input": len(gene_symbols),
        "comparison": "human_vs_human (TS vs non-TS atlas)",
    }

    save_results(
        result=result,
        p_value=p_value,
        null_distribution=null_dist,
        residuals=residuals,
        top_genes=top_genes,
        cell_types=cell_types,
        pca_info=pca_info,
        output_path=OUTPUT_DIR / "negctrl_results.json",
    )

    # ==================================================================
    # Step 8: Compare with human-vs-mouse results
    # ==================================================================
    print("\n" + "=" * 70)
    print("COMPARISON: Human-vs-Mouse vs Human-vs-Human")
    print("=" * 70)

    with open("./output/phase2/procrustes_results.json") as f:
        hvm = json.load(f)

    comparison = [
        {
            "comparison": "Human vs Mouse",
            "distance": hvm["procrustes"]["distance"],
            "ssr": hvm["procrustes"]["distance_squared"],
            "scaling": hvm["procrustes"]["scaling"],
            "p_value": hvm["permutation_test"]["p_value"],
            "significant_001": hvm["permutation_test"]["p_value"] < 0.01,
            "null_median": hvm["permutation_test"][
                "null_distribution_summary"
            ]["median"],
            "n_pca": hvm["pca"]["n_components"],
            "obs_to_null_ratio": (
                hvm["procrustes"]["distance"]
                / hvm["permutation_test"]["null_distribution_summary"]["median"]
            ),
        },
        {
            "comparison": "Human vs Human (neg. ctrl)",
            "distance": float(result.distance),
            "ssr": float(result.distance_squared),
            "scaling": float(result.scaling),
            "p_value": float(p_value),
            "significant_001": p_value < 0.01,
            "null_median": float(np.median(null_dist)),
            "n_pca": int(pca_model.n_components_),
            "obs_to_null_ratio": (
                float(result.distance) / float(np.median(null_dist))
            ),
        },
    ]

    comparison_df = pd.DataFrame(comparison)
    comparison_df.to_csv(OUTPUT_DIR / "comparison.csv", index=False)

    print(
        f"\n  {'Comparison':<30} {'PCs':>4} {'Distance':>10} "
        f"{'Null Med':>10} {'Obs/Null':>8} {'p-value':>10} {'Sig?':>5}"
    )
    print(f"  {'-' * 82}")

    for row in comparison:
        sig = "YES" if row["significant_001"] else "NO"
        print(
            f"  {row['comparison']:<30} {row['n_pca']:>4} "
            f"{row['distance']:>10.4f} {row['null_median']:>10.4f} "
            f"{row['obs_to_null_ratio']:>8.3f} "
            f"{row['p_value']:>10.6f} {sig:>5}"
        )

    # Interpretation
    hvm_p = hvm["permutation_test"]["p_value"]
    hvh_p = float(p_value)

    print(f"\n  Interpretation:")
    if hvh_p >= 0.01 and hvm_p < 0.01:
        print(
            f"  GOOD: Human-vs-mouse significant (p={hvm_p:.4f}) but "
            f"human-vs-human NOT (p={hvh_p:.4f})."
        )
        print(
            f"  → Cross-species geometric structure is biological, "
            f"not a batch effect."
        )
        print(f"  → Negative control PASSED.")
    elif hvh_p < 0.01 and hvm_p < 0.01:
        hvh_ratio = comparison[1]["obs_to_null_ratio"]
        hvm_ratio = comparison[0]["obs_to_null_ratio"]
        if hvh_ratio > hvm_ratio * 1.5:
            print(
                f"  MIXED: Both significant, but human-vs-human has weaker "
                f"alignment (obs/null={hvh_ratio:.3f} vs {hvm_ratio:.3f})."
            )
            print(
                f"  → Cross-species signal is STRONGER than batch signal. "
                f"Partial pass — document in paper."
            )
        else:
            print(
                f"  CRITICAL: Both significant with similar strength. "
                f"Cross-species signal may be a batch effect."
            )
            print(f"  → Negative control FAILED. Investigate further.")
    elif hvh_p >= 0.01 and hvm_p >= 0.01:
        print(
            f"  NOTE: Neither comparison significant. "
            f"Pipeline may lack power."
        )
    else:
        print(
            f"  UNEXPECTED: Human-vs-human significant but "
            f"human-vs-mouse not."
        )

    t_total = time.time() - t_start

    print(f"\n  Output files:")
    for f in sorted(OUTPUT_DIR.iterdir()):
        size = f.stat().st_size
        if size > 1024 * 1024:
            s = f"{size / 1024 / 1024:.1f} MB"
        elif size > 1024:
            s = f"{size / 1024:.1f} KB"
        else:
            s = f"{size} B"
        print(f"    {f.name:<45} {s:>10}")

    print(f"\n  Total runtime: {t_total:.1f}s")


if __name__ == "__main__":
    main()
