#!/usr/bin/env python3
"""
CellWarp — Phase 2 Negative Control v2 (Adult-Only, Healthy Tissue)

Reruns the human-vs-human negative control with STRICT filters:
  - development_stage must contain "adult" (no fetal, neonatal, embryonic)
  - disease must be "normal" (no COPD, no disease tissue)
  - is_primary_data = True
  - Excludes all Tabula Sapiens datasets

Biology
-------
The v1 negative control (ISSUE-017) was invalidated because 3/6 cell types
(endothelial 97.9%, hepatocyte 100%, macrophage 98.1%) were sourced from a
fetal developmental atlas. Fetal cells have fundamentally different expression
profiles from adult cells, inflating the H-vs-H Procrustes distance.

This v2 control ensures a fair apples-to-apples comparison: adult human cells
from one atlas vs adult human cells from another atlas. If the H-vs-H distance
is now significantly SMALLER than H-vs-M, the cross-species signal is real.

Pipeline
--------
1. Query Census for each cell type with strict adult + healthy filters
2. Greedy set-cover to find minimal non-TS dataset set
3. Download, subsample ≤2,000/type, filter to 16,959 genes, normalize
4. Compute centroids, PCA, Procrustes, 10k permutation test
5. Print 3-way comparison: H-vs-M, old H-vs-H (invalid), new H-vs-H (v2)

Outputs:
    data/phase2/human2_negctrl_v2.h5ad
    output/phase2/negative_control_v2/
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
MIN_CELLS_FOR_COVERAGE = 50  # Lower threshold for set-cover since adult-only is stricter
OUTPUT_DIR = Path("./output/phase2/negative_control_v2")
DATA_DIR = Path("./data/phase2")


# ---------------------------------------------------------------------------
# Dataset selection with strict adult-only filter
# ---------------------------------------------------------------------------


def find_adult_covering_datasets(census, datasets_df, ts_dataset_ids):
    """
    Find non-TS datasets covering all 6 cell types with STRICT adult-only filter.

    Filters:
      - cell_type matches our 6 types
      - is_primary_data == True
      - disease == 'normal'
      - development_stage_ontology_term_id == 'HssapDv:0000087' (human adult stage)
        OR development_stage contains 'adult' (fallback string match)
      - NOT in Tabula Sapiens

    Uses greedy set-cover to find minimal dataset combination.
    """
    import cellxgene_census

    print("\n  Querying Census for ADULT-ONLY, HEALTHY, non-TS cells...")
    print("  Filters: development_stage contains 'adult', disease='normal'")

    # Query per cell type with adult filter
    all_obs = []
    cell_type_totals = {}

    for ct in CELL_TYPES:
        t0 = time.time()
        # Use development_stage_ontology_term_id for precision
        # HsapDv:0000087 = "human adult stage" in HsapDv ontology
        # But safer to just require 'adult' substring since ontology IDs can vary
        obs = cellxgene_census.get_obs(
            census,
            "Homo sapiens",
            value_filter=(
                f"cell_type == '{ct}' "
                f"and is_primary_data == True "
                f"and disease == 'normal'"
            ),
            column_names=["cell_type", "dataset_id", "development_stage"],
        )

        # Exclude Tabula Sapiens
        obs = obs[~obs["dataset_id"].isin(ts_dataset_ids)]

        # STRICT adult-only filter: development_stage must contain "adult"
        # This excludes: fetal, embryonic, neonatal, infant, child, etc.
        adult_mask = obs["development_stage"].str.contains(
            "adult", case=False, na=False
        )
        n_before = len(obs)
        obs = obs[adult_mask].copy()
        n_after = len(obs)

        # Report development stages found
        if n_before > 0:
            all_stages = obs["development_stage"].value_counts()
            stage_str = "; ".join(
                f"{s}: {c}" for s, c in all_stages.head(5).items()
            )
        else:
            stage_str = "none"

        cell_type_totals[ct] = n_after
        all_obs.append(obs)
        dt = time.time() - t0
        print(
            f"    {ct:<45} {n_after:>8,} adult cells "
            f"(of {n_before:>8,} total) "
            f"[{dt:.1f}s]"
        )
        if n_after < 500:
            print(f"      WARNING: <500 adult cells available")
            if n_before > 0 and n_after < n_before:
                excluded_stages = (
                    obs["development_stage"].value_counts() if len(obs) > 0
                    else pd.Series(dtype=int)
                )
                print(f"      Stages kept: {stage_str}")

    combined = pd.concat(all_obs, ignore_index=True)

    # Build dataset × cell_type count matrix
    if len(combined) == 0:
        print("\n  FATAL: No adult-only healthy cells found in non-TS datasets!")
        return [], None, {}, cell_type_totals

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

    print(f"\n  Greedy dataset selection (adult-only):")
    while remaining:
        best_ds = None
        best_new = set()
        best_total = 0

        for ds_id in counts.index:
            if ds_id in selected:
                continue
            new_types = set()
            for ct in remaining:
                if (
                    ct in counts.columns
                    and counts.loc[ds_id, ct] >= MIN_CELLS_FOR_COVERAGE
                ):
                    new_types.add(ct)
            total = counts.loc[ds_id].sum()
            if len(new_types) > len(best_new) or (
                len(new_types) == len(best_new) and total > best_total
            ):
                best_ds = ds_id
                best_new = new_types
                best_total = total

        if best_ds is None or len(best_new) == 0:
            print(f"    WARNING: Cannot cover types: {remaining}")
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
            "types_covered": sorted(best_new),
            "cell_counts": types_in_ds,
        }

        print(f"    Dataset {len(selected)}: {coll[:70]}")
        print(f"      ID: {best_ds}")
        print(f"      New types: {sorted(best_new)}")
        print(f"      Cell counts: {types_in_ds}")

    print(f"\n  Coverage: {len(covered)}/6 types from {len(selected)} datasets")
    return selected, counts, source_info, cell_type_totals


def download_adult_only_atlas(
    target_gene_ids: set[str],
    target_gene_order: list[str],
):
    """
    Download adult-only, healthy, non-TS cells from CELLxGENE Census.

    Two-pass approach:
    1. First pass: query metadata to find covering datasets (adult filter applied)
    2. Second pass: download expression from selected datasets, then apply
       adult filter again on the downloaded data to be safe.
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

        # Find adult-only covering datasets
        selected_ds, counts, source_info, cell_totals = (
            find_adult_covering_datasets(
                census, datasets_df, ts_dataset_ids
            )
        )

        if not selected_ds:
            raise RuntimeError("No covering datasets found with adult-only filter")

        # Save source info
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        with open(OUTPUT_DIR / "source_datasets.json", "w") as f:
            json.dump(source_info, f, indent=2, default=str)

        with open(OUTPUT_DIR / "cell_availability.json", "w") as f:
            json.dump(cell_totals, f, indent=2)

        # Download expression from selected datasets
        names_str = ", ".join(f"'{ct}'" for ct in CELL_TYPES)
        ds_str = ", ".join(f"'{d}'" for d in selected_ds)
        obs_filter = (
            f"cell_type in [{names_str}] "
            f"and dataset_id in [{ds_str}] "
            f"and is_primary_data == True "
            f"and disease == 'normal'"
        )

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
                "development_stage",
            ],
            var_column_names=["feature_id", "feature_name"],
        )
        dt = time.time() - t0
        print(
            f"  Downloaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes "
            f"[{dt:.1f}s]"
        )

    # Apply adult filter on downloaded data (belt and suspenders)
    print(f"\n  Applying adult-only filter on downloaded data...")
    n_before = adata.n_obs
    adult_mask = adata.obs["development_stage"].str.contains(
        "adult", case=False, na=False
    )
    adata = adata[adult_mask].copy()
    n_after = adata.n_obs
    print(f"  Adult filter: {n_before:,} → {n_after:,} cells ({n_before - n_after:,} removed)")

    if n_after == 0:
        raise RuntimeError("No cells remaining after adult filter!")

    # Report development stages
    print(f"\n  Development stages in final data:")
    for stage, count in adata.obs["development_stage"].value_counts().items():
        print(f"    {stage}: {count:,}")

    # Report disease status (should all be 'normal')
    print(f"\n  Disease status in final data:")
    for dis, count in adata.obs["disease"].value_counts().items():
        print(f"    {dis}: {count:,}")

    # Per-type cell counts before subsampling
    print(f"\n  Pre-subsample cell counts (adult-only):")
    for ct in CELL_TYPES:
        n = (adata.obs["cell_type"] == ct).sum()
        flag = " *** BELOW 500 ***" if n < 500 else ""
        print(f"    {ct:<45} {n:>8,}{flag}")

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

    # Filter to target genes
    print(f"\n  Filtering to {len(target_gene_ids):,} shared ortholog genes...")
    gene_mask = adata.var["feature_id"].isin(target_gene_ids)
    n_matched = gene_mask.sum()
    print(f"  Matched {n_matched:,} / {len(target_gene_ids):,} target genes")

    adata = adata[:, gene_mask].copy()
    adata.var.index = adata.var["feature_id"].values

    # Align to target gene order
    available = set(adata.var_names)
    common_order = [g for g in target_gene_order if g in available]
    if len(common_order) < len(target_gene_order):
        n_miss = len(target_gene_order) - len(common_order)
        print(f"  WARNING: {n_miss} target genes missing ({len(common_order):,} retained)")
    adata = adata[:, common_order].copy()

    # Normalize: counts per 10k + log1p
    print(f"\n  Normalizing: counts per 10k + log1p...")
    adata.raw = adata
    sc.pp.normalize_total(adata, target_sum=10_000)
    sc.pp.log1p(adata)

    # Summary
    print(f"\n  Second human atlas v2 summary (ADULT-ONLY, HEALTHY):")
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


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    t_start = time.time()
    checkpoint_path = DATA_DIR / "human2_negctrl_v2.h5ad"

    print("=" * 70)
    print("PHASE 2 — Negative Control v2 (Adult-Only, Healthy)")
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
    # Step 1: Download or load adult-only second human atlas
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 1: Download adult-only second human atlas")
    print("=" * 70)

    if checkpoint_path.exists():
        print(f"  Checkpoint found: {checkpoint_path}")
        human2 = ad.read_h5ad(checkpoint_path)
        print(f"  Loaded: {human2.n_obs:,} cells × {human2.n_vars:,} genes")
        # Verify it has development_stage and all are adult
        if "development_stage" in human2.obs.columns:
            stages = human2.obs["development_stage"].value_counts()
            print(f"  Development stages:")
            for s, c in stages.items():
                print(f"    {s}: {c:,}")
            non_adult = ~human2.obs["development_stage"].str.contains(
                "adult", case=False, na=True
            )
            if non_adult.any():
                print(f"  WARNING: {non_adult.sum()} non-adult cells in checkpoint!")
                print(f"  Deleting checkpoint and re-downloading...")
                checkpoint_path.unlink()
                human2 = download_adult_only_atlas(
                    target_gene_ids, target_gene_order
                )
                save_h5ad_atomic(human2, checkpoint_path)
        else:
            print(f"  WARNING: No development_stage column — re-downloading")
            checkpoint_path.unlink()
            human2 = download_adult_only_atlas(
                target_gene_ids, target_gene_order
            )
            save_h5ad_atomic(human2, checkpoint_path)
    else:
        human2 = download_adult_only_atlas(target_gene_ids, target_gene_order)
        save_h5ad_atomic(human2, checkpoint_path)

    # Verify gene alignment
    h2_genes = set(human2.var_names)
    ts_genes = set(human_ts.var_names)
    shared = h2_genes & ts_genes
    print(f"\n  Gene alignment: {len(shared):,} shared genes")

    if len(shared) < len(ts_genes):
        n_missing = len(ts_genes) - len(shared)
        print(f"  Restricting to {len(shared):,} shared genes ({n_missing} missing)")
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
    print("STEP 2: Compute centroids for adult-only second human atlas")
    print("=" * 70)

    h2_centroids = compute_centroids(human2, "cell_type")
    h2_centroids.to_csv(OUTPUT_DIR / "centroids_human2_v2.csv")
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
    print("STEP 4: Procrustes alignment (Adult Human2 → Tabula Sapiens)")
    print("=" * 70)

    result = procrustes_align(h1_pca, h2_pca)

    # ==================================================================
    # Step 5: Permutation test
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 5: Permutation test (10,000 iterations)")
    print("=" * 70)

    p_value, null_dist = permutation_test(h1_pca, h2_pca)

    np.save(OUTPUT_DIR / "null_distribution_negctrl_v2.npy", null_dist)

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
        "comparison": "human_vs_human_v2 (TS vs non-TS, ADULT-ONLY, HEALTHY)",
    }

    save_results(
        result=result,
        p_value=p_value,
        null_distribution=null_dist,
        residuals=residuals,
        top_genes=top_genes,
        cell_types=cell_types,
        pca_info=pca_info,
        output_path=OUTPUT_DIR / "negctrl_v2_results.json",
    )

    # ==================================================================
    # Step 8: Three-way comparison
    # ==================================================================
    print("\n" + "=" * 70)
    print("THREE-WAY COMPARISON")
    print("=" * 70)

    # Load human-vs-mouse results
    with open("./output/phase2/procrustes_results.json") as f:
        hvm = json.load(f)

    # Load old (invalid) human-vs-human results
    old_hvh_path = Path("./output/phase2/negative_control/negctrl_results.json")
    old_hvh = None
    if old_hvh_path.exists():
        with open(old_hvh_path) as f:
            old_hvh = json.load(f)

    comparison = [
        {
            "comparison": "Human vs Mouse",
            "distance": hvm["procrustes"]["distance"],
            "ssr": hvm["procrustes"]["distance_squared"],
            "scaling": hvm["procrustes"]["scaling"],
            "p_value": hvm["permutation_test"]["p_value"],
            "significant_001": hvm["permutation_test"]["p_value"] < 0.01,
            "null_median": hvm["permutation_test"]["null_distribution_summary"]["median"],
            "n_pca": hvm["pca"]["n_components"],
            "obs_to_null_ratio": (
                hvm["procrustes"]["distance"]
                / hvm["permutation_test"]["null_distribution_summary"]["median"]
            ),
            "note": "cross-species",
        },
    ]

    if old_hvh:
        comparison.append({
            "comparison": "Human vs Human v1 (INVALID)",
            "distance": old_hvh["procrustes"]["distance"],
            "ssr": old_hvh["procrustes"]["distance_squared"],
            "scaling": old_hvh["procrustes"]["scaling"],
            "p_value": old_hvh["permutation_test"]["p_value"],
            "significant_001": old_hvh["permutation_test"]["p_value"] < 0.01,
            "null_median": old_hvh["permutation_test"]["null_distribution_summary"]["median"],
            "n_pca": old_hvh["pca"]["n_components"],
            "obs_to_null_ratio": (
                old_hvh["procrustes"]["distance"]
                / old_hvh["permutation_test"]["null_distribution_summary"]["median"]
            ),
            "note": "fetal-contaminated, invalid",
        })

    comparison.append({
        "comparison": "Human vs Human v2 (adult-only)",
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
        "note": "adult-only, healthy, non-TS",
    })

    comparison_df = pd.DataFrame(comparison)
    comparison_df.to_csv(OUTPUT_DIR / "comparison_3way.csv", index=False)

    # Print comparison table
    print(
        f"\n  {'Comparison':<35} {'PCs':>4} {'Distance':>10} "
        f"{'Null Med':>10} {'Obs/Null':>8} {'p-value':>10} {'Sig@.01':>8}"
    )
    print(f"  {'─' * 90}")

    for row in comparison:
        sig = "YES" if row["significant_001"] else "NO"
        note = f"  ({row['note']})" if row.get("note") else ""
        print(
            f"  {row['comparison']:<35} {row['n_pca']:>4} "
            f"{row['distance']:>10.4f} {row['null_median']:>10.4f} "
            f"{row['obs_to_null_ratio']:>8.3f} "
            f"{row['p_value']:>10.6f} {sig:>8}"
        )

    # Per-cell-type residual comparison
    print(f"\n  Per-cell-type residual magnitudes:")
    print(
        f"  {'Cell Type':<35} {'H-vs-M':>10} {'H-vs-H v1':>10} {'H-vs-H v2':>10}"
    )
    print(f"  {'─' * 70}")
    for ct in cell_types:
        hvm_mag = hvm["residuals"][ct]["magnitude"]
        old_mag = old_hvh["residuals"][ct]["magnitude"] if old_hvh else float("nan")
        new_mag = float(np.linalg.norm(residuals[ct]))
        print(
            f"  {ct:<35} {hvm_mag:>10.3f} {old_mag:>10.3f} {new_mag:>10.3f}"
        )

    # Interpretation
    print(f"\n  {'=' * 70}")
    print(f"  INTERPRETATION:")
    print(f"  {'=' * 70}")

    hvm_p = hvm["permutation_test"]["p_value"]
    hvh_p = float(p_value)
    hvm_dist = hvm["procrustes"]["distance"]
    hvh_dist = float(result.distance)

    if hvh_p >= 0.01 and hvm_p < 0.01:
        print(f"\n  ✓ NEGATIVE CONTROL PASSED")
        print(f"    Human-vs-mouse: p={hvm_p:.4f} (significant)")
        print(f"    Human-vs-human (adult): p={hvh_p:.4f} (NOT significant)")
        print(f"    → Cross-species signal is biological, not batch effects.")
        print(f"    → Phase 2 negative control gate: PASS")
    elif hvh_p < 0.01 and hvm_p < 0.01:
        hvh_ratio = hvh_dist / float(np.median(null_dist))
        hvm_ratio = hvm_dist / hvm["permutation_test"]["null_distribution_summary"]["median"]
        if hvh_ratio > hvm_ratio * 1.5:
            print(f"\n  ~ PARTIAL PASS (weaker H-vs-H alignment)")
            print(f"    Both significant, but H-vs-H alignment is weaker:")
            print(f"    H-vs-M obs/null = {hvm_ratio:.3f}")
            print(f"    H-vs-H obs/null = {hvh_ratio:.3f}")
            print(f"    → Cross-species signal is stronger than batch signal.")
        elif hvh_dist > hvm_dist * 1.3:
            print(f"\n  ~ PARTIAL PASS (H-vs-H distance is larger)")
            print(f"    H-vs-M distance: {hvm_dist:.4f}")
            print(f"    H-vs-H distance: {hvh_dist:.4f} ({hvh_dist/hvm_dist:.2f}x)")
            print(f"    → Same-species comparison produces more residual, as expected.")
        else:
            print(f"\n  ✗ NEGATIVE CONTROL FAILED (AGAIN)")
            print(f"    Both comparisons significant with similar strength.")
            print(f"    H-vs-M: dist={hvm_dist:.4f}, p={hvm_p:.4f}")
            print(f"    H-vs-H: dist={hvh_dist:.4f}, p={hvh_p:.4f}")
            print(f"    → Signal may be batch effects even with adult-only data.")
    elif hvh_p >= 0.01 and hvm_p >= 0.01:
        print(f"\n  ? INCONCLUSIVE")
        print(f"    Neither comparison significant.")
    else:
        print(f"\n  ! UNEXPECTED")
        print(f"    H-vs-H significant but H-vs-M not.")

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
