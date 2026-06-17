"""
Mouse Lemur CellWarp Pipeline
==============================
Runs the full CellWarp Procrustes pipeline on human vs mouse lemur:
  Step 1: Data preparation (filter, ortholog mapping, normalize, centroids)
  Step 2: Procrustes alignment + permutation test (10K)
  Step 3: Per-type residuals + ranking correlation with primary (human-mouse)
  Step 4: Evolutionary distance context

Biology
-------
Microcebus murinus (gray mouse lemur) diverged from humans ~75 Mya — intermediate
between macaque (~25-30 Mya) and mouse (~90 Mya). Adding mouse lemur gives a
fourth data point on the evolutionary distance axis:

  macaque (25 Mya) < mouse lemur (75 Mya) < mouse (90 Mya)

If Procrustes geometric coherence scales with evolutionary distance, we expect:
  obs/null(macaque) < obs/null(lemur) < obs/null(mouse)

The Tabula Microcebus (Ezran et al., Nature 2025) is a 244K-cell multi-tissue
atlas from 4 donors, sharing 15 cell types with our 35 primary analysis types.

Math
----
Same as primary pipeline:
1. Centroids: μ_t = (1/n) Σ x_i per cell type
2. PCA: joint space retaining ≥95% variance
3. Procrustes: min_{R,s,t} ‖X - (sYR + 1t')‖_F^2
4. Permutation test: 10K permutations of type pairings
5. Residuals: r_i = aligned_lemur_i - human_i
6. Ranking correlation: Spearman ρ vs primary analysis

Key design choice: We use pre-computed human 35-type centroids from the primary
analysis rather than re-downloading the full human data. The primary centroids
were computed from Tabula Sapiens data with all 35 types, so we just filter to
the 15 types shared with mouse lemur. Both centroid sets are then restricted to
the shared ortholog gene space before PCA.
"""

import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

# Add src to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))
from cellwarp.procrustes import (
    pca_reduce_centroids,
    procrustes_align,
    permutation_test,
    compute_residual_vectors,
    map_residuals_to_genes,
    CELL_TYPE_COL,
)

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
OUT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_ROOT / "data"

RANDOM_SEED = 42
MAX_CELLS_PER_TYPE = 2_000
MIN_CELLS_PER_TYPE = 500
N_PERMUTATIONS = 10_000

# Pre-computed human 35-type centroids
HUMAN_CENTROIDS_PATH = PROJECT_ROOT / "output/phase2/scaled_35types/centroids_human_35.csv"

# Primary analysis residuals
PRIMARY_RESIDUALS = PROJECT_ROOT / "output/phase2/scaled_35types/residuals_ranked.csv"

# Mouse lemur data
LEMUR_H5AD = DATA_DIR / "mouse_lemur" / "tabula_microcebus_LCA_complete.h5ad"

# Ortholog mapping
LEMUR_ORTHOLOGS = OUT_DIR / "biomart_mouse_lemur_human_orthologs.csv"

# The 15 cell types passing the >=500 gate in mouse lemur data
PASSING_TYPES = [
    "B cell",
    "CD4-positive, alpha-beta T cell",
    "CD8-positive, alpha-beta T cell",
    "T cell",
    "endothelial cell",
    "enterocyte of epithelium of large intestine",
    "fibroblast",
    "macrophage",
    "mature NK T cell",
    "mesenchymal stem cell",
    "monocyte",
    "natural killer cell",
    "neutrophil",
    "pancreatic acinar cell",
    "plasma cell",
]


def step1_data_preparation():
    """Step 1: Load lemur data, compute centroids, load human centroids, align genes."""
    print("=" * 70)
    print("STEP 1: DATA PREPARATION")
    print("=" * 70)

    # ── Load mouse lemur data ─────────────────────────────────────────
    print("\n[1a] Loading mouse lemur data...")
    lemur = ad.read_h5ad(LEMUR_H5AD)
    print(f"  Raw: {lemur.n_obs:,} cells × {lemur.n_vars:,} genes")

    # Filter to 10x 3' v2 only
    print("\n[1b] Filtering to 10x 3' v2 only...")
    pre = lemur.n_obs
    lemur = lemur[lemur.obs["assay"] == "10x 3' v2"].copy()
    print(f"  After 10x filter: {lemur.n_obs:,} cells (removed {pre - lemur.n_obs:,} Smart-seq2)")

    # Filter to passing cell types
    print(f"\n[1c] Filtering to {len(PASSING_TYPES)} passing cell types...")
    lemur = lemur[lemur.obs["cell_type"].isin(PASSING_TYPES)].copy()
    print(f"  After type filter: {lemur.n_obs:,} cells")

    # Check per-type counts after 10x filtering
    ct_counts = lemur.obs["cell_type"].value_counts()
    print(f"\n  Per-type counts (10x only):")
    final_types = []
    for ct in sorted(PASSING_TYPES):
        n = ct_counts.get(ct, 0)
        status = "PASS" if n >= MIN_CELLS_PER_TYPE else "DROPPED"
        if n >= MIN_CELLS_PER_TYPE:
            final_types.append(ct)
        print(f"    {ct:<55} {n:>8,} {status}")

    if len(final_types) < len(PASSING_TYPES):
        print(f"\n  WARNING: {len(PASSING_TYPES) - len(final_types)} types dropped "
              f"after 10x-only filtering")
        lemur = lemur[lemur.obs["cell_type"].isin(final_types)].copy()

    print(f"\n  Final types: {len(final_types)}")

    # ── Subsample ─────────────────────────────────────────────────────
    print(f"\n[1d] Subsampling to ≤{MAX_CELLS_PER_TYPE:,} cells/type...")
    rng = np.random.default_rng(RANDOM_SEED)
    indices = []
    for ct in sorted(final_types):
        ct_idx = np.where(lemur.obs["cell_type"] == ct)[0]
        if len(ct_idx) > MAX_CELLS_PER_TYPE:
            selected = rng.choice(ct_idx, size=MAX_CELLS_PER_TYPE, replace=False)
            selected.sort()
            indices.extend(selected.tolist())
            print(f"    {ct}: {len(ct_idx):,} → {MAX_CELLS_PER_TYPE:,}")
        else:
            indices.extend(ct_idx.tolist())
            print(f"    {ct}: {len(ct_idx):,} (kept all)")
    lemur = lemur[indices].copy()
    print(f"  Subsampled: {lemur.n_obs:,} cells")

    # ── Map orthologs ─────────────────────────────────────────────────
    print(f"\n[1e] Loading ortholog mapping...")
    orthologs = pd.read_csv(LEMUR_ORTHOLOGS)
    print(f"  Ortholog pairs loaded: {len(orthologs):,}")

    lemur_to_human = dict(zip(
        orthologs["lemur_ensembl_id"],
        orthologs["human_ensembl_id"],
    ))

    # ── Load pre-computed human centroids ─────────────────────────────
    print(f"\n[1f] Loading pre-computed human 35-type centroids...")
    human_centroids_full = pd.read_csv(HUMAN_CENTROIDS_PATH, index_col=0)
    print(f"  Human centroids: {human_centroids_full.shape[0]} types × "
          f"{human_centroids_full.shape[1]:,} genes")

    # Filter to matching types
    human_centroids = human_centroids_full.loc[
        human_centroids_full.index.isin(final_types)
    ].copy()
    print(f"  After filtering to {len(final_types)} types: {human_centroids.shape}")

    human_gene_ids = set(human_centroids.columns)
    lemur_gene_ids = set(lemur.var_names)

    # Find usable ortholog pairs
    usable_lemur_ids = []
    usable_human_ids = []
    for lid in lemur_gene_ids:
        hid = lemur_to_human.get(lid)
        if hid and hid in human_gene_ids:
            usable_lemur_ids.append(lid)
            usable_human_ids.append(hid)

    print(f"  Lemur genes in data: {len(lemur_gene_ids):,}")
    print(f"  Human genes in centroids: {len(human_gene_ids):,}")
    print(f"  Usable ortholog pairs: {len(usable_lemur_ids):,}")

    # Sort by human Ensembl ID
    pairs = sorted(zip(usable_human_ids, usable_lemur_ids), key=lambda x: x[0])
    gene_order_human = [p[0] for p in pairs]
    gene_order_lemur = [p[1] for p in pairs]

    # Filter human centroids to shared genes
    human_centroids = human_centroids[gene_order_human]

    # Filter and re-index lemur data to human gene space
    lemur_filt = lemur[:, gene_order_lemur].copy()
    lemur_filt.var["original_lemur_id"] = lemur_filt.var_names.copy()
    lemur_filt.var_names = pd.Index(gene_order_human)

    n_shared_genes = len(gene_order_human)
    print(f"  Aligned gene space: {n_shared_genes:,} genes")

    # ── Normalize lemur data ──────────────────────────────────────────
    print(f"\n[1g] Checking/normalizing lemur data...")
    import scipy.sparse as sp
    if sp.issparse(lemur_filt.X):
        sample = lemur_filt.X[:100].toarray()
    else:
        sample = lemur_filt.X[:100]

    max_val = sample.max()
    mean_val = sample.mean()
    print(f"  Lemur X sample: max={max_val:.2f}, mean={mean_val:.4f}")

    if max_val > 50:
        print("  Lemur appears to be raw counts. Normalizing...")
        sc.pp.normalize_total(lemur_filt, target_sum=1e4)
        sc.pp.log1p(lemur_filt)
    elif max_val > 15:
        print("  Lemur appears CPM-normalized. Log-transforming...")
        sc.pp.log1p(lemur_filt)
    else:
        print("  Lemur appears already log-normalized. No transformation needed.")

    # ── Compute lemur centroids ───────────────────────────────────────
    print(f"\n[1h] Computing lemur centroids...")
    lemur_centroids = {}
    for ct in sorted(final_types):
        mask = lemur_filt.obs["cell_type"] == ct
        if sp.issparse(lemur_filt.X):
            mean_vec = np.asarray(lemur_filt[mask].X.mean(axis=0)).flatten()
        else:
            mean_vec = np.mean(lemur_filt[mask].X, axis=0)
        lemur_centroids[ct] = mean_vec
        n_cells = mask.sum()
        print(f"    {ct:<55} {n_cells:>6,} cells → centroid ({n_shared_genes:,} genes)")

    lemur_centroids_df = pd.DataFrame(lemur_centroids, index=gene_order_human).T
    lemur_centroids_df.index.name = "cell_type"

    # Save cell counts
    l_counts = lemur_filt.obs["cell_type"].value_counts()
    count_rows = []
    for ct in final_types:
        count_rows.append({
            "cell_type": ct,
            "human_cells": "from_precomputed_centroids",
            "lemur_cells": int(l_counts.get(ct, 0)),
        })
    pd.DataFrame(count_rows).to_csv(OUT_DIR / "centroid_cell_counts.csv", index=False)

    return human_centroids, lemur_centroids_df, final_types, n_shared_genes


def step2_procrustes(human_centroids, lemur_centroids, cell_types):
    """Step 2: PCA → Procrustes → permutation test."""
    print("\n" + "=" * 70)
    print("STEP 2: PROCRUSTES ANALYSIS")
    print("=" * 70)

    # ── PCA ────────────────────────────────────���──────────────────────
    print(f"\n[2a] PCA reduction ({len(cell_types)} types × {human_centroids.shape[1]:,} genes)...")
    human_pca, lemur_pca, pca_model, pca_types = pca_reduce_centroids(
        human_centroids, lemur_centroids
    )

    pca_info = {
        "n_components": int(pca_model.n_components_),
        "variance_explained": pca_model.explained_variance_ratio_.tolist(),
        "cumulative_variance": float(np.cumsum(pca_model.explained_variance_ratio_)[-1]),
    }

    # ── Procrustes alignment ──────────────────────────────────────────
    print(f"\n[2b] Procrustes alignment (lemur → human)...")
    result = procrustes_align(human_pca, lemur_pca)

    # ── Permutation test ──────────────────────────────────────────────
    print(f"\n[2c] Permutation test ({N_PERMUTATIONS:,} iterations)...")
    p_value, null_dist = permutation_test(
        human_pca, lemur_pca,
        n_permutations=N_PERMUTATIONS,
        seed=RANDOM_SEED,
    )

    obs_null_ratio = result.distance / np.median(null_dist)
    obs_null_ratio_mean = result.distance / np.mean(null_dist)

    print(f"\n  obs/null ratio (median): {obs_null_ratio:.4f}")
    print(f"  obs/null ratio (mean):   {obs_null_ratio_mean:.4f}")

    return result, p_value, null_dist, pca_model, pca_types, pca_info


def step3_residuals_and_ranking(result, cell_types, pca_model, human_centroids):
    """Step 3: Per-type residuals + ranking correlation."""
    print("\n" + "=" * 70)
    print("STEP 3: PER-TYPE RESIDUALS AND RANKING")
    print("=" * 70)

    # ── Residuals ─────────────────────────────────────────────────────
    print(f"\n[3a] Computing per-type residuals...")
    residuals = compute_residual_vectors(result, cell_types)

    # Rank by magnitude
    resid_mags = {ct: float(np.linalg.norm(residuals[ct])) for ct in cell_types}
    ranked = sorted(resid_mags.items(), key=lambda x: x[1], reverse=True)

    print(f"\n  Residual ranking (most → least divergent):")
    print(f"  {'Rank':>4} {'Cell Type':<55} {'Magnitude':>10}")
    print(f"  " + "-" * 72)
    for i, (ct, mag) in enumerate(ranked, 1):
        print(f"  {i:>4} {ct:<55} {mag:>10.4f}")

    # Save residuals CSV
    resid_rows = []
    for i, (ct, mag) in enumerate(ranked, 1):
        pct = (mag ** 2 / result.distance_squared * 100) if result.distance_squared > 0 else 0
        resid_rows.append({
            "rank": i,
            "cell_type": ct,
            "residual_magnitude": mag,
            "pct_of_ssr": pct,
        })
    resid_df = pd.DataFrame(resid_rows)
    resid_df.to_csv(OUT_DIR / "per_type_residuals.csv", index=False)

    # ── Map residuals to gene space ───────────────────────────────────
    print(f"\n[3b] Mapping residuals to gene space...")
    gene_names = human_centroids.columns.tolist()
    top_genes = map_residuals_to_genes(residuals, pca_model, gene_names, n_top=20)

    # ── Ranking correlation with primary ──────────────────────────────
    print(f"\n[3c] Ranking correlation with primary (human-mouse)...")
    primary_df = pd.read_csv(PRIMARY_RESIDUALS)
    primary_lookup = {row["cell_type"]: row["residual_magnitude"]
                      for _, row in primary_df.iterrows()}

    shared = [ct for ct in cell_types if ct in primary_lookup]
    n_shared = len(shared)
    print(f"  Types shared with primary 35: {n_shared}")

    if n_shared >= 5:
        comp_rows = []
        for ct in shared:
            comp_rows.append({
                "cell_type": ct,
                "primary_residual": primary_lookup[ct],
                "lemur_residual": resid_mags[ct],
            })
        comp_df = pd.DataFrame(comp_rows)
        comp_df["primary_rank"] = comp_df["primary_residual"].rank(ascending=False).astype(int)
        comp_df["lemur_rank"] = comp_df["lemur_residual"].rank(ascending=False).astype(int)
        comp_df["rank_shift"] = abs(comp_df["primary_rank"] - comp_df["lemur_rank"])
        comp_df = comp_df.sort_values("primary_rank")

        rho, rho_p = spearmanr(comp_df["primary_residual"], comp_df["lemur_residual"])

        print(f"  Spearman ρ = {rho:.4f}")
        print(f"  p-value   = {rho_p:.4f}")
        print(f"  Significant at α=0.05: {'YES' if rho_p < 0.05 else 'NO'}")

        print(f"\n  {'Cell Type':<55} {'P_rank':>6} {'L_rank':>6} {'Shift':>6}")
        print(f"  " + "-" * 76)
        for _, row in comp_df.iterrows():
            shift = int(row["rank_shift"])
            flag = " **" if shift >= 5 else ""
            print(f"  {row['cell_type']:<55} "
                  f"{int(row['primary_rank']):>6} {int(row['lemur_rank']):>6} "
                  f"{shift:>6}{flag}")

        comp_df.to_csv(OUT_DIR / "ranking_comparison_with_primary.csv", index=False)
    else:
        rho, rho_p = float("nan"), float("nan")
        print(f"  Too few shared types ({n_shared}) for meaningful correlation")

    # Also correlate with macaque ranking
    macaque_results_path = PROJECT_ROOT / "output/macaque_pipeline/primary_procrustes_results.json"
    mac_rho, mac_rho_p = float("nan"), float("nan")
    mac_shared_n = 0
    if macaque_results_path.exists():
        print(f"\n[3d] Ranking correlation with macaque (human-macaque)...")
        with open(macaque_results_path) as f:
            mac_data = json.load(f)
        mac_residuals = {ct: d["magnitude"] for ct, d in mac_data["residuals"].items()}
        mac_shared = [ct for ct in cell_types if ct in mac_residuals]
        mac_shared_n = len(mac_shared)
        print(f"  Types shared with macaque: {mac_shared_n}")

        if mac_shared_n >= 5:
            mac_comp = []
            for ct in mac_shared:
                mac_comp.append({
                    "cell_type": ct,
                    "macaque_residual": mac_residuals[ct],
                    "lemur_residual": resid_mags[ct],
                })
            mac_df = pd.DataFrame(mac_comp)
            mac_rho, mac_rho_p = spearmanr(mac_df["macaque_residual"], mac_df["lemur_residual"])
            print(f"  Spearman ρ (lemur vs macaque) = {mac_rho:.4f}, p = {mac_rho_p:.4f}")
            mac_df.to_csv(OUT_DIR / "ranking_comparison_with_macaque.csv", index=False)

    ranking_results = {
        "vs_primary": {"rho": float(rho), "p": float(rho_p), "n_types": n_shared},
        "vs_macaque": {"rho": float(mac_rho), "p": float(mac_rho_p), "n_types": mac_shared_n},
    }

    return residuals, resid_mags, top_genes, ranking_results


def step4_evolutionary_context(obs_null_ratio, p_value, n_types):
    """Step 4: Place result in evolutionary distance context."""
    print("\n" + "=" * 70)
    print("STEP 4: EVOLUTIONARY DISTANCE CONTEXT")
    print("=" * 70)

    species_pairs = [
        {
            "pair": "human-macaque",
            "divergence_mya": 25,
            "obs_null": 0.841,
            "p_value": 0.0002,
            "n_types": 20,
            "source": "DECISION-131",
        },
        {
            "pair": "human-mouse_lemur",
            "divergence_mya": 75,
            "obs_null": float(obs_null_ratio),
            "p_value": float(p_value),
            "n_types": n_types,
            "source": "This analysis",
        },
        {
            "pair": "human-mouse",
            "divergence_mya": 90,
            "obs_null": 0.522,
            "p_value": 0.0001,
            "n_types": 35,
            "source": "Primary analysis",
        },
    ]

    print(f"\n  {'Pair':<25} {'Div (Mya)':>10} {'obs/null':>10} {'p-value':>10} {'n_types':>8}")
    print(f"  " + "-" * 68)
    for sp in sorted(species_pairs, key=lambda x: x["divergence_mya"]):
        print(f"  {sp['pair']:<25} {sp['divergence_mya']:>10} "
              f"{sp['obs_null']:>10.4f} {sp['p_value']:>10.4f} {sp['n_types']:>8}")

    # Check ordering
    sorted_by_dist = sorted(species_pairs, key=lambda x: x["divergence_mya"])
    obs_nulls = [sp["obs_null"] for sp in sorted_by_dist]
    distances = [sp["divergence_mya"] for sp in sorted_by_dist]

    monotone_decreasing = all(obs_nulls[i] >= obs_nulls[i+1]
                              for i in range(len(obs_nulls)-1))

    print(f"\n  Evolutionary distance ordering check:")
    print(f"    Macaque (25 Mya) obs/null = {obs_nulls[0]:.4f}")
    print(f"    Mouse lemur (75 Mya) obs/null = {obs_nulls[1]:.4f}")
    print(f"    Mouse (90 Mya) obs/null = {obs_nulls[2]:.4f}")
    print(f"    Monotonically decreasing: {'YES' if monotone_decreasing else 'NO'}")

    if monotone_decreasing:
        print(f"    → Geometric coherence DECREASES with evolutionary distance.")
        print(f"      Three points isn't a regression, but the ordering is informative.")
    else:
        print(f"    → Ordering does NOT follow evolutionary distance.")
        print(f"      Confounds: n_types differs (20, {species_pairs[1]['n_types']}, 35),")
        print(f"      cell type composition varies, different atlas technologies.")

    rho_dist, p_dist = spearmanr(distances, obs_nulls)
    print(f"\n  Spearman ρ(divergence, obs/null) = {rho_dist:.4f} (p = {p_dist:.4f})")
    print(f"  NOTE: n=3 points. This is descriptive, not a statistical test.")

    all_sig = all(sp["p_value"] < 0.01 for sp in species_pairs)
    print(f"\n  All three species pairs significant at α=0.01: {'YES' if all_sig else 'NO'}")
    if all_sig:
        print(f"  → Global geometric coherence confirmed across three evolutionary distances:")
        print(f"    25, 75, and 90 million years of primate-rodent divergence.")

    return species_pairs, monotone_decreasing


def save_outputs(result, p_value, null_dist, pca_info, cell_types,
                 residuals, resid_mags, top_genes, ranking_results,
                 species_pairs, monotone, n_shared_genes):
    """Step 5: Save all outputs."""
    print("\n" + "=" * 70)
    print("STEP 5: SAVING OUTPUTS")
    print("=" * 70)

    obs_null_ratio = result.distance / np.median(null_dist)

    # ── Procrustes results JSON ───────────────────────────────────────
    results_json = {
        "analysis": "human_vs_mouse_lemur",
        "species": {
            "reference": "Homo sapiens (Tabula Sapiens via Census, pre-computed centroids)",
            "target": "Microcebus murinus (Tabula Microcebus, Ezran et al. 2025)",
        },
        "divergence_mya": 75,
        "n_types": len(cell_types),
        "cell_types": cell_types,
        "gene_space": n_shared_genes,
        "procrustes": {
            "distance": float(result.distance),
            "distance_squared": float(result.distance_squared),
            "scaling": float(result.scaling),
            "rotation_det": float(np.linalg.det(result.rotation)),
        },
        "permutation_test": {
            "p_value": float(p_value),
            "obs_null_ratio": float(obs_null_ratio),
            "null_median": float(np.median(null_dist)),
            "null_mean": float(np.mean(null_dist)),
            "null_std": float(np.std(null_dist)),
            "significant_at_001": bool(p_value < 0.01),
            "n_permutations": N_PERMUTATIONS,
        },
        "pca": pca_info,
        "residuals": {
            ct: {
                "magnitude": float(resid_mags[ct]),
                "pct_ssr": float(resid_mags[ct] ** 2 / result.distance_squared * 100)
                if result.distance_squared > 0 else 0.0,
            }
            for ct in cell_types
        },
        "ranking_correlation": ranking_results,
        "evolutionary_context": {
            "species_pairs": species_pairs,
            "monotone_with_distance": monotone,
        },
    }

    with open(OUT_DIR / "procrustes_results.json", "w") as f:
        json.dump(results_json, f, indent=2)
    print(f"  Saved: procrustes_results.json")

    # ── Ranking correlation JSON ──────────────────────────────────────
    with open(OUT_DIR / "ranking_correlation.json", "w") as f:
        json.dump(ranking_results, f, indent=2)
    print(f"  Saved: ranking_correlation.json")

    # ── Cell type mapping CSV ─────────────────────────────────────────
    mapping_rows = []
    for ct in cell_types:
        mapping_rows.append({
            "cell_type": ct,
            "in_primary_35": True,
            "residual_magnitude": resid_mags[ct],
        })
    pd.DataFrame(mapping_rows).to_csv(OUT_DIR / "cell_type_mapping.csv", index=False)
    print(f"  Saved: cell_type_mapping.csv")

    # ── Null distribution ─────────────────────────────────────────────
    np.save(OUT_DIR / "null_distribution.npy", null_dist)
    print(f"  Saved: null_distribution.npy")

    # ── Summary markdown ─────────────────────────────────────��────────
    ranked = sorted(resid_mags.items(), key=lambda x: x[1], reverse=True)

    summary = f"""# Mouse Lemur CellWarp Analysis — Summary

**Date:** 2026-04-05
**Species pair:** Human (*H. sapiens*) vs Mouse lemur (*M. murinus*)
**Evolutionary distance:** ~75 Mya
**Data sources:**
- Human: Tabula Sapiens (pre-computed centroids from primary 35-type analysis)
- Mouse lemur: Tabula Microcebus (Ezran et al., Nature 2025, CELLxGENE Discover)

---

## Global Coherence (Procrustes)

| Metric | Value |
|--------|-------|
| Cell types | {len(cell_types)} |
| Shared ortholog genes | {n_shared_genes:,} |
| Procrustes distance | {result.distance:.4f} |
| Null median distance | {np.median(null_dist):.4f} |
| **obs/null ratio** | **{obs_null_ratio:.4f}** |
| **p-value** | **{p_value:.6f}** |
| Significant (alpha=0.01) | {'YES' if p_value < 0.01 else 'NO'} |
| PCA components | {pca_info['n_components']} |
| Cumulative variance | {pca_info['cumulative_variance']*100:.1f}% |
| Scaling factor | {result.scaling:.6f} |
| Rotation det | {np.linalg.det(result.rotation):+.6f} |

## Comparison Across Species Pairs

| Species Pair | Div (Mya) | obs/null | p-value | n types |
|-------------|-----------|----------|---------|---------|
| Human-macaque | 25 | 0.841 | 0.0002 | 20 |
| **Human-mouse lemur** | **75** | **{obs_null_ratio:.4f}** | **{p_value:.4f}** | **{len(cell_types)}** |
| Human-mouse | 90 | 0.522 | <0.0001 | 35 |

obs/null monotonically decreasing with distance: **{'YES' if monotone else 'NO'}**

## Per-Type Residuals (ranked by divergence)

| Rank | Cell Type | Magnitude | % SSR |
|------|-----------|-----------|-------|
"""

    for i, (ct, mag) in enumerate(ranked, 1):
        pct = mag ** 2 / result.distance_squared * 100 if result.distance_squared > 0 else 0
        summary += f"| {i} | {ct} | {mag:.4f} | {pct:.1f}% |\n"

    summary += f"""
## Ranking Correlation

### vs Primary (human-mouse, 35 types)
- Spearman rho = {ranking_results['vs_primary']['rho']:.4f}
- p-value = {ranking_results['vs_primary']['p']:.4f}
- n shared types = {ranking_results['vs_primary']['n_types']}

### vs Macaque (human-macaque, 20 types)
- Spearman rho = {ranking_results['vs_macaque']['rho']:.4f}
- p-value = {ranking_results['vs_macaque']['p']:.4f}
- n shared types = {ranking_results['vs_macaque']['n_types']}

**Expected:** Ranking correlation is weak/non-significant based on simulation
findings and five prior non-replications (MCA rho=0.120, Sun2023 rho=0.146,
PanSci rho=0.194, CellHint rho=-0.386, macaque rho=0.137).

## Key Observations

1. **Global geometric coherence at 75 Mya**: {'Confirmed' if p_value < 0.01 else 'Not confirmed'}
   — the Procrustes transformation structure {'holds' if p_value < 0.01 else 'does not hold'}
   at an intermediate evolutionary distance.

2. **obs/null ratio = {obs_null_ratio:.3f}**: {'Falls between macaque (0.841) and mouse (0.522), consistent with evolutionary distance scaling.' if 0.522 <= obs_null_ratio <= 0.841 else 'Does not follow expected evolutionary distance scaling.'}

3. **Three-species geometric conservation**: {'All three pairs significant (p < 0.01). Cross-species geometric coherence is a robust signal across 25-90 Mya of primate and rodent evolution.' if p_value < 0.01 else 'Two of three pairs significant.'}

## Data Quality Notes

- Mouse lemur data filtered to 10x 3' v2 only (95% of atlas)
- 15 cell types pass >=500 cell gate (exactly at threshold)
- Hepatocyte has 458 cells (just below gate) — excluded
- 60% immune types, 40% non-immune (better balance than macaque 65%)
- 4 donors only (vs 47 in RIRA macaque data)
- Human centroids from pre-computed primary analysis (Tabula Sapiens, same
  centroids used in the 35-type primary pipeline)
"""

    with open(OUT_DIR / "mouse_lemur_summary.md", "w") as f:
        f.write(summary)
    print(f"  Saved: mouse_lemur_summary.md")


def main():
    print("=" * 70)
    print("CELLWARP: HUMAN vs MOUSE LEMUR PIPELINE")
    print("=" * 70)

    # Step 1: Data preparation
    human_centroids, lemur_centroids, cell_types, n_genes = step1_data_preparation()

    # Step 2: Procrustes
    result, p_value, null_dist, pca_model, pca_types, pca_info = \
        step2_procrustes(human_centroids, lemur_centroids, cell_types)

    # Step 3: Residuals and ranking
    residuals, resid_mags, top_genes, ranking_results = \
        step3_residuals_and_ranking(result, pca_types, pca_model, human_centroids)

    # Step 4: Evolutionary context
    obs_null_ratio = result.distance / np.median(null_dist)
    species_pairs, monotone = step4_evolutionary_context(
        obs_null_ratio, p_value, len(cell_types)
    )

    # Step 5: Save outputs
    save_outputs(
        result, p_value, null_dist, pca_info, pca_types,
        residuals, resid_mags, top_genes, ranking_results,
        species_pairs, monotone, n_genes,
    )

    # ── Final console summary ─────────────────────────────────────────
    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    print(f"  Species: Human vs Mouse lemur (~75 Mya)")
    print(f"  Cell types: {len(cell_types)}")
    print(f"  Shared genes: {n_genes:,}")
    print(f"  obs/null ratio: {obs_null_ratio:.4f}")
    print(f"  p-value: {p_value:.6f}")
    print(f"  Significant: {'YES' if p_value < 0.01 else 'NO'}")
    print(f"  Ranking rho (vs primary): {ranking_results['vs_primary']['rho']:.4f}")
    print(f"  Monotone with distance: {'YES' if monotone else 'NO'}")
    print(f"\n  All outputs saved to: {OUT_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
