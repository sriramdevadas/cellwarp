#!/usr/bin/env python3
"""
CellWarp — Donor-Split Analysis with Shared PCA Space

Reruns the donor-split within-species baseline control using a shared PCA space
for both the within-species and cross-species comparisons within each split.

Fix for the independent-PCA confound
-------------------------------------
The original analysis (script 41) fits PCA separately for each comparison:
  - Within-species PCA on (human_half1 + human_half2) → ~21 components
  - Cross-species PCA on (human_half1 + mouse) → ~29 components

This means the obs/null ratios are computed in different coordinate systems
with different dimensionalities, which could bias the delta.

This script instead fits a single PCA on all three centroid matrices combined:
  (human_half1 + human_half2 + mouse) → one shared space
Then projects all three into the same coordinate system before running
both Procrustes analyses. This makes the obs/null ratios directly comparable.

Same constraints as original:
  - Same 100 splits, same seeds (RANDOM_SEED=42), same donor partitions
  - Same 95% variance threshold for PCA
  - Same 1000 permutations per split
  - Same >=100 cells per half threshold

Inputs:
    data/phase2_scaled/human_scaled.h5ad
    data/phase2_scaled/mouse_scaled.h5ad

Outputs (new files in analysis/donor_split/, existing files untouched):
    donor_split_shared_pca_distributions.csv
    donor_split_shared_pca_results.json
    donor_split_shared_pca_figure.png
    donor_split_comparison.json
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.decomposition import PCA

import anndata as ad

from cellwarp.procrustes import (
    _procrustes_distance,
    PCA_VARIANCE_THRESHOLD,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

RANDOM_SEED = 42
N_SPLITS = 100
N_PERMUTATIONS_PER_SPLIT = 1000
MIN_CELLS_PER_HALF = 100
OUTPUT_DIR = Path("./analysis/donor_split")

HUMAN_DATA_PATH = Path("./data/phase2_scaled/human_scaled.h5ad")
MOUSE_DATA_PATH = Path("./data/phase2_scaled/mouse_scaled.h5ad")

INFEASIBLE_TYPES = [
    "pancreatic acinar cell",
    "luminal epithelial cell of mammary gland",
    "pancreatic ductal cell",
]

CROSS_SPECIES_OBS_NULL_PRIMARY = 0.5223
SELF_COMPARISON_OBS_NULL = 0.033


# ---------------------------------------------------------------------------
# Helpers (copied from script 41 to keep this script self-contained)
# ---------------------------------------------------------------------------


def load_data_h5ad(path: Path) -> ad.AnnData:
    """Load h5ad via h5py to avoid anndata version issues with uns/log1p."""
    import h5py
    import scipy.sparse as sp

    with h5py.File(path, "r") as f:
        X_group = f["X"]
        if "data" in X_group:
            data = X_group["data"][:]
            indices = X_group["indices"][:]
            indptr = X_group["indptr"][:]
            shape = tuple(X_group.attrs.get("shape", X_group.attrs.get("h5sparse_shape")))
            X = sp.csr_matrix((data, indices, indptr), shape=shape)
        else:
            X = X_group[:]

        obs_data = {}
        for key in f["obs"].keys():
            if key == "__categories":
                continue
            item = f["obs"][key]
            if isinstance(item, h5py.Dataset):
                obs_data[key] = item[:]
            elif isinstance(item, h5py.Group):
                if "codes" in item and "categories" in item:
                    codes = item["codes"][:]
                    cats = item["categories"][:]
                    if cats.dtype.kind in ("S", "O"):
                        cats = [c.decode() if isinstance(c, bytes) else c for c in cats]
                    obs_data[key] = pd.Categorical.from_codes(codes, categories=cats)
        obs = pd.DataFrame(obs_data)

        var_data = {}
        for key in f["var"].keys():
            if key == "__categories":
                continue
            item = f["var"][key]
            if isinstance(item, h5py.Dataset):
                d = item[:]
                if d.dtype.kind in ("S", "O"):
                    d = [v.decode() if isinstance(v, bytes) else v for v in d]
                var_data[key] = d
            elif isinstance(item, h5py.Group):
                if "codes" in item and "categories" in item:
                    codes = item["codes"][:]
                    cats = item["categories"][:]
                    if cats.dtype.kind in ("S", "O"):
                        cats = [c.decode() if isinstance(c, bytes) else c for c in cats]
                    var_data[key] = pd.Categorical.from_codes(codes, categories=cats)
        var = pd.DataFrame(var_data)

        if "_index" in obs.columns:
            obs.index = [str(x.decode() if isinstance(x, bytes) else x) for x in obs["_index"]]
            obs = obs.drop(columns=["_index"])
        if "_index" in var.columns:
            var.index = [str(x.decode() if isinstance(x, bytes) else x) for x in var["_index"]]
            var = var.drop(columns=["_index"])

    return ad.AnnData(X=X, obs=obs, var=var)


def compute_centroids_from_subset(adata: ad.AnnData, cell_types: list[str]) -> pd.DataFrame:
    """Compute mean expression centroids for specified cell types."""
    gene_ids = adata.var_names.tolist()
    centroids = {}
    for ct in cell_types:
        mask = adata.obs["cell_type"] == ct
        if mask.sum() == 0:
            continue
        centroids[ct] = np.asarray(adata[mask].X.mean(axis=0)).flatten()
    df = pd.DataFrame(centroids, index=gene_ids).T
    df.index.name = "cell_type"
    return df


def balanced_donor_split(donors: list[str], cells_per_donor: dict[str, int],
                         rng: np.random.RandomState) -> tuple[list[str], list[str]]:
    """Randomly partition donors into two balanced groups (shuffle then greedy)."""
    shuffled = list(donors)
    rng.shuffle(shuffled)
    half1, half2 = [], []
    n1, n2 = 0, 0
    for d in shuffled:
        c = cells_per_donor[d]
        if n1 <= n2:
            half1.append(d)
            n1 += c
        else:
            half2.append(d)
            n2 += c
    return half1, half2


def compute_obs_null_ratio(X: np.ndarray, Y: np.ndarray,
                           n_perms: int, seed: int) -> dict:
    """Compute Procrustes distance, permutation null, and obs/null ratio."""
    n = X.shape[0]
    rng = np.random.RandomState(seed)

    observed = _procrustes_distance(X, Y)

    null_distances = np.zeros(n_perms)
    for i in range(n_perms):
        perm = rng.permutation(n)
        null_distances[i] = _procrustes_distance(X, Y[perm])

    null_mean = np.mean(null_distances)
    null_median = np.median(null_distances)
    obs_null_ratio = float(observed / null_median) if null_median > 0 else float("inf")
    n_leq = int(np.sum(null_distances <= observed))
    p_value = (n_leq + 1) / (n_perms + 1)

    return {
        "observed_distance": float(observed),
        "null_mean": float(null_mean),
        "null_median": float(np.median(null_distances)),
        "null_std": float(np.std(null_distances)),
        "obs_null_ratio": obs_null_ratio,
        "p_value": p_value,
        "n_types": n,
    }


def shared_pca_reduce(h1_centroids: pd.DataFrame,
                      h2_centroids: pd.DataFrame,
                      mouse_centroids: pd.DataFrame,
                      variance_threshold: float = PCA_VARIANCE_THRESHOLD,
                      ) -> tuple[np.ndarray, np.ndarray, np.ndarray, int, list[str]]:
    """
    Fit a single PCA on all three centroid matrices and project each into it.

    Math: Stack all three centroid matrices into a (3n x G) matrix:
        C = [human_half1; human_half2; mouse]
    Fit PCA to C, retaining components for >=variance_threshold cumulative variance.
    Project each subset into this shared PCA space.

    Args:
        h1_centroids: (n, G) DataFrame — human half1 centroids.
        h2_centroids: (n, G) DataFrame — human half2 centroids.
        mouse_centroids: (n, G) DataFrame — mouse centroids.
        variance_threshold: Minimum cumulative variance to retain.

    Returns:
        Tuple of (h1_pca, h2_pca, mouse_pca, n_components, cell_types).
    """
    cell_types = sorted(h1_centroids.index.tolist())
    assert sorted(h2_centroids.index.tolist()) == cell_types
    assert sorted(mouse_centroids.index.tolist()) == cell_types

    h1_mat = h1_centroids.loc[cell_types].values   # (n, G)
    h2_mat = h2_centroids.loc[cell_types].values   # (n, G)
    m_mat = mouse_centroids.loc[cell_types].values  # (n, G)

    # Stack all three: (3n, G)
    combined = np.vstack([h1_mat, h2_mat, m_mat])
    n_types = len(cell_types)

    pca = PCA(
        n_components=variance_threshold,
        svd_solver="full",
        random_state=RANDOM_SEED,
    )
    combined_pca = pca.fit_transform(combined)  # (3n, k)
    k = pca.n_components_

    h1_pca = combined_pca[:n_types]          # (n, k)
    h2_pca = combined_pca[n_types:2*n_types] # (n, k)
    m_pca = combined_pca[2*n_types:]         # (n, k)

    return h1_pca, h2_pca, m_pca, k, cell_types


def main():
    t0 = time.time()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("DONOR-SPLIT ANALYSIS — SHARED PCA SPACE")
    print("=" * 80)

    # ------------------------------------------------------------------
    # Load data (reuse loader from script 41)
    # ------------------------------------------------------------------
    print("\n[1] Loading data...")
    human = load_data_h5ad(HUMAN_DATA_PATH)
    mouse = load_data_h5ad(MOUSE_DATA_PATH)
    print(f"    Human: {human.n_obs:,} cells x {human.n_vars:,} genes")
    print(f"    Mouse: {mouse.n_obs:,} cells x {mouse.n_vars:,} genes")

    all_types = sorted(human.obs["cell_type"].unique())
    feasible_types = [ct for ct in all_types if ct not in INFEASIBLE_TYPES]
    print(f"    Feasible cell types: {len(feasible_types)} / {len(all_types)}")

    # Precompute donor info
    donor_info = {}
    for ct in feasible_types:
        ct_mask = human.obs["cell_type"] == ct
        ct_obs = human.obs[ct_mask]
        cpd = ct_obs.groupby("donor_id").size().to_dict()
        donor_info[ct] = {
            "donors": list(cpd.keys()),
            "cells_per_donor": cpd,
            "total_cells": int(ct_mask.sum()),
        }

    all_donors = sorted(set(
        d for info in donor_info.values() for d in info["donors"]
    ))
    global_cpd = human.obs.groupby("donor_id").size().to_dict()
    print(f"    Total donors (feasible types): {len(all_donors)}")

    # Precompute mouse centroids
    print("\n    Computing mouse centroids (once)...")
    mouse_centroids_all = compute_centroids_from_subset(mouse, feasible_types)
    print(f"    Mouse centroids: {mouse_centroids_all.shape}")

    # ------------------------------------------------------------------
    # Run donor splits with shared PCA
    # ------------------------------------------------------------------
    print(f"\n[2] Running {N_SPLITS} donor splits (SHARED PCA)...")
    print(f"    Permutations per split: {N_PERMUTATIONS_PER_SPLIT}")
    print("-" * 80)

    rng = np.random.RandomState(RANDOM_SEED)
    split_results = []

    for split_i in range(N_SPLITS):
        split_seed = rng.randint(0, 2**31)

        # Same donor partition as original
        split_rng = np.random.RandomState(split_seed)
        half1_donors, half2_donors = balanced_donor_split(
            all_donors, global_cpd, split_rng
        )

        # Compute centroids from each donor half
        half1_types = []
        half1_centroids = {}
        half2_centroids = {}

        gene_ids = human.var_names.tolist()

        for ct in feasible_types:
            info = donor_info[ct]
            cpd = info["cells_per_donor"]

            h1_donors_ct = [d for d in half1_donors if d in cpd]
            h2_donors_ct = [d for d in half2_donors if d in cpd]
            h1_cells = sum(cpd[d] for d in h1_donors_ct)
            h2_cells = sum(cpd[d] for d in h2_donors_ct)

            if h1_cells < MIN_CELLS_PER_HALF or h2_cells < MIN_CELLS_PER_HALF:
                continue

            h1_mask = (human.obs["cell_type"] == ct) & (human.obs["donor_id"].isin(h1_donors_ct))
            h2_mask = (human.obs["cell_type"] == ct) & (human.obs["donor_id"].isin(h2_donors_ct))

            h1_mean = np.asarray(human[h1_mask].X.mean(axis=0)).flatten()
            h2_mean = np.asarray(human[h2_mask].X.mean(axis=0)).flatten()

            half1_centroids[ct] = h1_mean
            half2_centroids[ct] = h2_mean
            half1_types.append(ct)

        if len(half1_types) < 4:
            print(f"  Split {split_i + 1:>3}: SKIPPED (only {len(half1_types)} types)")
            continue

        # Build centroid DataFrames for the shared type set
        h1_df = pd.DataFrame(half1_centroids, index=gene_ids).T
        h1_df.index.name = "cell_type"
        h2_df = pd.DataFrame(half2_centroids, index=gene_ids).T
        h2_df.index.name = "cell_type"

        # Mouse centroids restricted to the same type set
        shared_types = sorted(
            set(half1_types) & set(mouse_centroids_all.index)
        )
        if len(shared_types) < 4:
            print(f"  Split {split_i + 1:>3}: SKIPPED (only {len(shared_types)} shared types)")
            continue

        h1_df = h1_df.loc[shared_types]
        h2_df = h2_df.loc[shared_types]
        mouse_sub = mouse_centroids_all.loc[shared_types]

        # ---- SHARED PCA: fit once on all three ----
        try:
            h1_pca, h2_pca, m_pca, n_components, ct_list = shared_pca_reduce(
                h1_df, h2_df, mouse_sub,
                variance_threshold=PCA_VARIANCE_THRESHOLD,
            )
        except Exception as e:
            print(f"  Split {split_i + 1:>3}: Shared PCA FAILED ({e})")
            continue

        # ---- WITHIN-SPECIES: half1 vs half2 in shared PCA space ----
        ws_result = compute_obs_null_ratio(
            h1_pca, h2_pca, N_PERMUTATIONS_PER_SPLIT, split_seed
        )

        # ---- CROSS-SPECIES: half1 vs mouse in shared PCA space ----
        cs_result = compute_obs_null_ratio(
            h1_pca, m_pca, N_PERMUTATIONS_PER_SPLIT, split_seed
        )

        delta = cs_result["obs_null_ratio"] - ws_result["obs_null_ratio"]

        split_results.append({
            "split": split_i + 1,
            "seed": split_seed,
            "n_half1_donors": len(half1_donors),
            "n_half2_donors": len(half2_donors),
            "n_types": len(shared_types),
            "n_pca_components": n_components,
            "ws_obs_distance": ws_result["observed_distance"],
            "ws_null_mean": ws_result["null_mean"],
            "ws_obs_null_ratio": ws_result["obs_null_ratio"],
            "ws_p_value": ws_result["p_value"],
            "cs_obs_distance": cs_result["observed_distance"],
            "cs_null_mean": cs_result["null_mean"],
            "cs_obs_null_ratio": cs_result["obs_null_ratio"],
            "cs_p_value": cs_result["p_value"],
            "delta_obs_null": delta,
        })

        if (split_i + 1) % 10 == 0 or split_i == 0:
            print(f"  Split {split_i + 1:>3}: "
                  f"ws={ws_result['obs_null_ratio']:.3f} cs={cs_result['obs_null_ratio']:.3f} "
                  f"delta={delta:+.3f} k={n_components} n={len(shared_types)}")

    # ------------------------------------------------------------------
    # Analyze results
    # ------------------------------------------------------------------
    print(f"\n\n[3] RESULTS — SHARED PCA ({len(split_results)} successful splits)")
    print("=" * 80)

    if len(split_results) == 0:
        print("  ERROR: No successful splits. Cannot proceed.")
        return

    df = pd.DataFrame(split_results)

    ws_ratios = df["ws_obs_null_ratio"].values
    cs_ratios = df["cs_obs_null_ratio"].values
    deltas = df["delta_obs_null"].values
    k_values = df["n_pca_components"].values

    print(f"\n  PCA COMPONENTS (shared space):")
    print(f"    Median: {np.median(k_values):.0f}")
    print(f"    Range:  [{np.min(k_values)}, {np.max(k_values)}]")
    print(f"    Mean:   {np.mean(k_values):.1f}")

    print(f"\n  WITHIN-SPECIES (donor-split) obs/null ratio:")
    print(f"    Mean:   {np.mean(ws_ratios):.4f}")
    print(f"    Median: {np.median(ws_ratios):.4f}")
    print(f"    Std:    {np.std(ws_ratios):.4f}")
    print(f"    Range:  [{np.min(ws_ratios):.4f}, {np.max(ws_ratios):.4f}]")
    print(f"    95% CI: [{np.percentile(ws_ratios, 2.5):.4f}, {np.percentile(ws_ratios, 97.5):.4f}]")

    ws_sig = (df["ws_p_value"] < 0.01).sum()
    print(f"    Significant (p<0.01): {ws_sig} / {len(df)} splits")

    print(f"\n  CROSS-SPECIES (matched) obs/null ratio:")
    print(f"    Mean:   {np.mean(cs_ratios):.4f}")
    print(f"    Median: {np.median(cs_ratios):.4f}")
    print(f"    Std:    {np.std(cs_ratios):.4f}")
    print(f"    Range:  [{np.min(cs_ratios):.4f}, {np.max(cs_ratios):.4f}]")
    print(f"    95% CI: [{np.percentile(cs_ratios, 2.5):.4f}, {np.percentile(cs_ratios, 97.5):.4f}]")

    cs_sig = (df["cs_p_value"] < 0.01).sum()
    print(f"    Significant (p<0.01): {cs_sig} / {len(df)} splits")

    print(f"\n  DELTA (cross-species - donor-split):")
    print(f"    Mean:   {np.mean(deltas):+.4f}")
    print(f"    Median: {np.median(deltas):+.4f}")
    print(f"    Std:    {np.std(deltas):.4f}")
    print(f"    95% CI: [{np.percentile(deltas, 2.5):+.4f}, {np.percentile(deltas, 97.5):+.4f}]")

    n_positive = (deltas > 0).sum()
    print(f"    Positive (cs > ws): {n_positive} / {len(deltas)} ({100*n_positive/len(deltas):.1f}%)")

    t_stat, t_pval = stats.ttest_rel(cs_ratios, ws_ratios)
    print(f"    Paired t-test: t={t_stat:.3f}, p={t_pval:.2e}")

    w_stat, w_pval = stats.wilcoxon(deltas)
    print(f"    Wilcoxon signed-rank: W={w_stat:.1f}, p={w_pval:.2e}")

    # Early stop check
    if n_positive / len(deltas) < 0.5 or np.percentile(deltas, 97.5) < 0:
        print("\n  *** ALERT: Delta has flipped sign or become non-significant! ***")
        print("  *** STOP — this changes the manuscript plan. ***")

    consistently_positive = n_positive / len(deltas) > 0.95
    print(f"\n  INTERPRETATION:")
    if consistently_positive and np.percentile(deltas, 2.5) > 0:
        print(f"    Delta is consistently positive (95% CI excludes 0).")
        print(f"    Cross-species is LESS COHERENT than within-species.")
        print(f"    => EVOLUTIONARY SIGNAL CONFIRMED in shared PCA space.")
    elif n_positive / len(deltas) > 0.5:
        print(f"    Delta is mostly positive ({100*n_positive/len(deltas):.0f}% of splits).")
        print(f"    Evidence for evolutionary signal, but some overlap.")
    else:
        print(f"    Delta is not consistently positive.")
        print(f"    Within-species variation may explain the cross-species signal.")

    # ------------------------------------------------------------------
    # Save results (NEW files only — do not overwrite originals)
    # ------------------------------------------------------------------
    df.to_csv(OUTPUT_DIR / "donor_split_shared_pca_distributions.csv", index=False)

    summary = {
        "pca_strategy": "shared",
        "n_splits": len(df),
        "n_splits_attempted": N_SPLITS,
        "n_permutations_per_split": N_PERMUTATIONS_PER_SPLIT,
        "min_cells_per_half": MIN_CELLS_PER_HALF,
        "random_seed": RANDOM_SEED,
        "infeasible_types": INFEASIBLE_TYPES,
        "pca_components": {
            "median": float(np.median(k_values)),
            "min": int(np.min(k_values)),
            "max": int(np.max(k_values)),
            "mean": float(np.mean(k_values)),
        },
        "within_species": {
            "mean_obs_null_ratio": float(np.mean(ws_ratios)),
            "median_obs_null_ratio": float(np.median(ws_ratios)),
            "std_obs_null_ratio": float(np.std(ws_ratios)),
            "ci_95": [float(np.percentile(ws_ratios, 2.5)),
                      float(np.percentile(ws_ratios, 97.5))],
            "pct_significant_p001": float(ws_sig / len(df)),
            "median_n_types": float(df["n_types"].median()),
        },
        "cross_species_matched": {
            "mean_obs_null_ratio": float(np.mean(cs_ratios)),
            "median_obs_null_ratio": float(np.median(cs_ratios)),
            "std_obs_null_ratio": float(np.std(cs_ratios)),
            "ci_95": [float(np.percentile(cs_ratios, 2.5)),
                      float(np.percentile(cs_ratios, 97.5))],
            "pct_significant_p001": float(cs_sig / len(df)),
        },
        "delta": {
            "mean": float(np.mean(deltas)),
            "median": float(np.median(deltas)),
            "std": float(np.std(deltas)),
            "ci_95": [float(np.percentile(deltas, 2.5)),
                      float(np.percentile(deltas, 97.5))],
            "pct_positive": float(n_positive / len(deltas)),
            "paired_ttest_p": float(t_pval),
            "wilcoxon_p": float(w_pval),
        },
        "reference_values": {
            "primary_cross_species_obs_null": CROSS_SPECIES_OBS_NULL_PRIMARY,
            "self_comparison_obs_null": SELF_COMPARISON_OBS_NULL,
        },
        "runtime_seconds": time.time() - t0,
    }

    with open(OUTPUT_DIR / "donor_split_shared_pca_results.json", "w") as f:
        json.dump(summary, f, indent=2)

    # ------------------------------------------------------------------
    # Comparison JSON (shared vs independent PCA)
    # ------------------------------------------------------------------
    indep_path = OUTPUT_DIR / "donor_split_results.json"
    if indep_path.exists():
        with open(indep_path) as f:
            indep = json.load(f)

        # Get median k values from independent PCA distributions
        indep_dist_path = OUTPUT_DIR / "donor_split_distributions.csv"
        indep_df = pd.read_csv(indep_dist_path) if indep_dist_path.exists() else None

        comparison = {
            "independent_pca": {
                "ws_median": indep["within_species"]["median_obs_null_ratio"],
                "cs_median": indep["cross_species_matched"]["median_obs_null_ratio"],
                "delta_median": indep["delta"]["median"],
                "delta_ci_95": indep["delta"]["ci_95"],
                "delta_positive_fraction": indep["delta"]["pct_positive"],
                "median_k_ws": "~21 (varies per split, not recorded per-split)",
                "median_k_cs": "~29 (varies per split, not recorded per-split)",
            },
            "shared_pca": {
                "ws_median": float(np.median(ws_ratios)),
                "cs_median": float(np.median(cs_ratios)),
                "delta_median": float(np.median(deltas)),
                "delta_ci_95": [float(np.percentile(deltas, 2.5)),
                                float(np.percentile(deltas, 97.5))],
                "delta_positive_fraction": float(n_positive / len(deltas)),
                "median_k": float(np.median(k_values)),
            },
            "change": {
                "ws_median_change": float(np.median(ws_ratios)) - indep["within_species"]["median_obs_null_ratio"],
                "cs_median_change": float(np.median(cs_ratios)) - indep["cross_species_matched"]["median_obs_null_ratio"],
                "delta_median_change": float(np.median(deltas)) - indep["delta"]["median"],
            },
        }

        with open(OUTPUT_DIR / "donor_split_comparison.json", "w") as f:
            json.dump(comparison, f, indent=2)
        print(f"\n  Comparison saved to: {OUTPUT_DIR}/donor_split_comparison.json")

    print(f"  Results saved to: {OUTPUT_DIR}/donor_split_shared_pca_results.json")
    print(f"  Distributions saved to: {OUTPUT_DIR}/donor_split_shared_pca_distributions.csv")

    # ------------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------------
    print("\n[4] Generating figure...")
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))

    # Panel A: Within-species histogram
    ax = axes[0]
    ax.hist(ws_ratios, bins=25, color="#4C72B0", alpha=0.8, edgecolor="white")
    ax.axvline(np.median(ws_ratios), color="red", linestyle="--", linewidth=1.5,
               label=f"Median = {np.median(ws_ratios):.3f}")
    ax.set_xlabel("obs/null ratio", fontsize=11)
    ax.set_ylabel("Count", fontsize=11)
    ax.set_title("A. Within-species\n(donor-split, shared PCA)", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)

    # Panel B: Cross-species histogram
    ax = axes[1]
    ax.hist(cs_ratios, bins=25, color="#DD8452", alpha=0.8, edgecolor="white")
    ax.axvline(np.median(cs_ratios), color="red", linestyle="--", linewidth=1.5,
               label=f"Median = {np.median(cs_ratios):.3f}")
    ax.axvline(CROSS_SPECIES_OBS_NULL_PRIMARY, color="black", linestyle=":",
               linewidth=1.5, label=f"Primary (35-type) = {CROSS_SPECIES_OBS_NULL_PRIMARY:.3f}")
    ax.set_xlabel("obs/null ratio", fontsize=11)
    ax.set_title("B. Cross-species\n(matched subsets, shared PCA)", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)

    # Panel C: Three-way comparison
    ax = axes[2]
    positions = [1, 2, 3]
    colors = ["#55A868", "#4C72B0", "#DD8452"]
    labels = ["Self-comparison\n(subsampling)", "Donor-split\n(within-species)", "Cross-species\n(human vs mouse)"]

    ax.scatter([1], [SELF_COMPARISON_OBS_NULL], s=120, color=colors[0],
               zorder=5, edgecolors="black", linewidth=1)

    bp1 = ax.boxplot([ws_ratios], positions=[2], widths=0.5,
                     patch_artist=True, showfliers=False)
    bp1["boxes"][0].set_facecolor(colors[1])
    bp1["boxes"][0].set_alpha(0.7)
    bp1["medians"][0].set_color("black")

    bp2 = ax.boxplot([cs_ratios], positions=[3], widths=0.5,
                     patch_artist=True, showfliers=False)
    bp2["boxes"][0].set_facecolor(colors[2])
    bp2["boxes"][0].set_alpha(0.7)
    bp2["medians"][0].set_color("black")

    ax.axhline(1.0, color="gray", linestyle=":", linewidth=1,
               label="Permutation null (1.0)")

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("obs/null ratio", fontsize=11)
    ax.set_title("C. Coherence hierarchy\n(shared PCA)", fontsize=12, fontweight="bold")
    ax.legend(fontsize=8, loc="upper left")
    ax.set_ylim(bottom=-0.05)

    plt.tight_layout()
    fig_path = OUTPUT_DIR / "donor_split_shared_pca_figure.png"
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Figure saved to: {fig_path}")

    print(f"\n  Figure description:")
    print(f"    Panel A: Within-species (donor-split) obs/null in shared PCA space. Median = {np.median(ws_ratios):.3f}.")
    print(f"    Panel B: Cross-species obs/null in same shared PCA space. Median = {np.median(cs_ratios):.3f}.")
    print(f"    Panel C: Coherence hierarchy: self-comparison ({SELF_COMPARISON_OBS_NULL}) < "
          f"donor-split ({np.median(ws_ratios):.3f}) < cross-species ({np.median(cs_ratios):.3f}) < null (1.0).")
    print(f"    All comparisons in identical {np.median(k_values):.0f}-dimensional PCA space per split.")

    elapsed = time.time() - t0
    print(f"\n  Total runtime: {elapsed:.1f} seconds ({elapsed/60:.1f} minutes)")


if __name__ == "__main__":
    main()
