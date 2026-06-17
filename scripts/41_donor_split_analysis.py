#!/usr/bin/env python3
"""
CellWarp — Donor-Split Within-Species Baseline Control

Tests whether cross-species geometric coherence exceeds same-protocol within-species
coherence, providing evidence for an evolutionary component in the signal.

Logic
-----
Split Tabula Sapiens donors into two non-overlapping groups, compute cell-type
centroids independently on each half, run Procrustes between the two human halves.
Compare the within-species obs/null ratio against the cross-species obs/null ratio.
If cross-species obs/null > donor-split obs/null (i.e., cross-species is less
coherent), there is an evolutionary component beyond what within-species variation
can explain.

Biology
-------
The donor-split control isolates "protocol + stochastic + individual variation" from
"protocol + stochastic + individual + evolutionary variation". If cross-species
Procrustes is significantly less coherent than donor-split Procrustes (higher obs/null),
the difference must come from evolutionary divergence — the only factor that differs
between the two comparisons.

Pipeline (per split iteration)
------------------------------
1. Randomly partition human donors into two balanced groups
2. For each feasible cell type, compute centroids on each donor half
3. Drop cell types with <100 cells in either half
4. PCA on combined centroids (human_half1 + human_half2)
5. Procrustes alignment (half2 → half1)
6. Permutation test (1000 permutations)
7. Record obs/null ratio
8. ALSO: PCA on (human_half1 + mouse), Procrustes, permutation test
   restricted to the same cell type subset → cross-species obs/null

Inputs:
    data/phase2_scaled/human_scaled.h5ad   — normalized human data (35 types)
    data/phase2_scaled/mouse_scaled.h5ad   — normalized mouse data (35 types)

Outputs (all in analysis/donor_split/):
    donor_split_results.json       — full results
    donor_split_figure.png         — three-way comparison figure
    donor_split_distributions.csv  — per-split metrics
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
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

from cellwarp.procrustes import (
    _procrustes_distance,
    compute_centroids,
    pca_reduce_centroids,
    permutation_test,
    procrustes_align,
    PCA_VARIANCE_THRESHOLD,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

RANDOM_SEED = 42
N_SPLITS = 100                  # Number of random donor splits
N_PERMUTATIONS_PER_SPLIT = 1000  # Permutations per split (lighter than 10k)
MIN_CELLS_PER_HALF = 100        # Minimum cells in each half for valid centroid
OUTPUT_DIR = Path("./analysis/donor_split")

HUMAN_DATA_PATH = Path("./data/phase2_scaled/human_scaled.h5ad")
MOUSE_DATA_PATH = Path("./data/phase2_scaled/mouse_scaled.h5ad")

# Known infeasible types (from feasibility analysis: <4 donors)
INFEASIBLE_TYPES = [
    "pancreatic acinar cell",
    "luminal epithelial cell of mammary gland",
    "pancreatic ductal cell",
]

# Reference values from primary analysis
CROSS_SPECIES_OBS_NULL_PRIMARY = 0.5223  # 35-type obs/null ratio
SELF_COMPARISON_OBS_NULL = 0.033          # from negative controls


def load_data_h5ad(path: Path) -> ad.AnnData:
    """Load h5ad, handling anndata version issues with uns/log1p."""
    import h5py
    import scipy.sparse as sp

    with h5py.File(path, "r") as f:
        # Read X matrix
        X_group = f["X"]
        if "data" in X_group:
            # Sparse CSR
            data = X_group["data"][:]
            indices = X_group["indices"][:]
            indptr = X_group["indptr"][:]
            shape = tuple(X_group.attrs.get("shape", X_group.attrs.get("h5sparse_shape")))
            X = sp.csr_matrix((data, indices, indptr), shape=shape)
        else:
            X = X_group[:]

        # Read obs
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
                    obs_data[key] = pd.Categorical.from_codes(
                        codes, categories=cats
                    )
        obs = pd.DataFrame(obs_data)

        # Read var
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
                    var_data[key] = pd.Categorical.from_codes(
                        codes, categories=cats
                    )
        var = pd.DataFrame(var_data)

        # Set index
        if "_index" in obs.columns:
            obs.index = [str(x.decode() if isinstance(x, bytes) else x) for x in obs["_index"]]
            obs = obs.drop(columns=["_index"])
        if "_index" in var.columns:
            var.index = [str(x.decode() if isinstance(x, bytes) else x) for x in var["_index"]]
            var = var.drop(columns=["_index"])

    adata = ad.AnnData(X=X, obs=obs, var=var)
    return adata


def compute_obs_null_ratio(X: np.ndarray, Y: np.ndarray,
                           n_perms: int, seed: int) -> dict:
    """Compute Procrustes distance, permutation null, and obs/null ratio silently."""
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


def compute_centroids_from_subset(adata: ad.AnnData, cell_types: list[str]) -> pd.DataFrame:
    """Compute centroids for specified cell types from an AnnData subset."""
    gene_ids = adata.var_names.tolist()
    centroids = {}
    for ct in cell_types:
        mask = adata.obs["cell_type"] == ct
        n_cells = mask.sum()
        if n_cells == 0:
            continue
        mean_vec = np.asarray(adata[mask].X.mean(axis=0)).flatten()
        centroids[ct] = mean_vec
    df = pd.DataFrame(centroids, index=gene_ids).T
    df.index.name = "cell_type"
    return df


def balanced_donor_split(donors: list[str], cells_per_donor: dict[str, int],
                         rng: np.random.RandomState) -> tuple[list[str], list[str]]:
    """Randomly partition donors into two balanced groups.

    Strategy: shuffle donors, then greedily assign to the group with fewer cells.
    This ensures roughly equal cell counts while being stochastic.
    """
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


def main():
    t0 = time.time()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("DONOR-SPLIT WITHIN-SPECIES BASELINE CONTROL")
    print("=" * 80)

    # ------------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------------
    print("\n[1] Loading data...")
    human = load_data_h5ad(HUMAN_DATA_PATH)
    mouse = load_data_h5ad(MOUSE_DATA_PATH)
    print(f"    Human: {human.n_obs:,} cells x {human.n_vars:,} genes")
    print(f"    Mouse: {mouse.n_obs:,} cells x {mouse.n_vars:,} genes")

    # Identify feasible cell types
    all_types = sorted(human.obs["cell_type"].unique())
    feasible_types = [ct for ct in all_types if ct not in INFEASIBLE_TYPES]
    print(f"    Feasible cell types: {len(feasible_types)} / {len(all_types)}")

    # Precompute: donors per cell type and cell counts
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

    # All unique donors across feasible types
    all_donors = sorted(set(
        d for info in donor_info.values() for d in info["donors"]
    ))
    global_cpd = human.obs.groupby("donor_id").size().to_dict()
    print(f"    Total donors (feasible types): {len(all_donors)}")

    # Also need mouse centroids for cross-species comparison
    # Precompute these once (they don't change across splits)
    print("\n    Computing mouse centroids (once)...")
    mouse_centroids_all = compute_centroids_from_subset(mouse, feasible_types)
    print(f"    Mouse centroids: {mouse_centroids_all.shape}")

    # ------------------------------------------------------------------
    # Run donor splits
    # ------------------------------------------------------------------
    print(f"\n[2] Running {N_SPLITS} donor splits...")
    print(f"    Permutations per split: {N_PERMUTATIONS_PER_SPLIT}")
    print("-" * 80)

    rng = np.random.RandomState(RANDOM_SEED)
    split_results = []

    for split_i in range(N_SPLITS):
        split_seed = rng.randint(0, 2**31)

        # Partition donors globally
        split_rng = np.random.RandomState(split_seed)
        half1_donors, half2_donors = balanced_donor_split(
            all_donors, global_cpd, split_rng
        )

        # For each feasible type, compute centroids from each half
        half1_types = []  # Types that pass the cell count threshold
        half1_centroids = {}
        half2_centroids = {}

        for ct in feasible_types:
            info = donor_info[ct]
            cpd = info["cells_per_donor"]

            # Count cells in each half for this cell type
            h1_donors_ct = [d for d in half1_donors if d in cpd]
            h2_donors_ct = [d for d in half2_donors if d in cpd]
            h1_cells = sum(cpd[d] for d in h1_donors_ct)
            h2_cells = sum(cpd[d] for d in h2_donors_ct)

            if h1_cells < MIN_CELLS_PER_HALF or h2_cells < MIN_CELLS_PER_HALF:
                continue

            # Compute centroids from each half
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

        # Build centroid DataFrames
        gene_ids = human.var_names.tolist()
        h1_df = pd.DataFrame(half1_centroids, index=gene_ids).T
        h1_df.index.name = "cell_type"
        h2_df = pd.DataFrame(half2_centroids, index=gene_ids).T
        h2_df.index.name = "cell_type"

        # ---- WITHIN-SPECIES: half1 vs half2 ----
        try:
            h1_pca, h2_pca, pca_ws, ct_ws = pca_reduce_centroids(
                h1_df, h2_df, variance_threshold=PCA_VARIANCE_THRESHOLD
            )
        except Exception as e:
            print(f"  Split {split_i + 1:>3}: PCA FAILED ({e})")
            continue

        ws_result = compute_obs_null_ratio(
            h1_pca, h2_pca, N_PERMUTATIONS_PER_SPLIT, split_seed
        )

        # ---- CROSS-SPECIES: half1 vs mouse (same type subset) ----
        # Use mouse centroids restricted to the same types
        mouse_sub = mouse_centroids_all.loc[
            [ct for ct in half1_types if ct in mouse_centroids_all.index]
        ]
        h1_sub = h1_df.loc[[ct for ct in half1_types if ct in mouse_sub.index]]

        # Ensure matching types
        shared_types = sorted(set(h1_sub.index) & set(mouse_sub.index))
        if len(shared_types) < 4:
            print(f"  Split {split_i + 1:>3}: SKIPPED cross-species (only {len(shared_types)} shared types)")
            continue

        h1_sub = h1_sub.loc[shared_types]
        mouse_sub = mouse_sub.loc[shared_types]

        try:
            h1_pca_cs, m_pca_cs, pca_cs, ct_cs = pca_reduce_centroids(
                h1_sub, mouse_sub, variance_threshold=PCA_VARIANCE_THRESHOLD
            )
        except Exception as e:
            print(f"  Split {split_i + 1:>3}: Cross-species PCA FAILED ({e})")
            continue

        cs_result = compute_obs_null_ratio(
            h1_pca_cs, m_pca_cs, N_PERMUTATIONS_PER_SPLIT, split_seed
        )

        delta = cs_result["obs_null_ratio"] - ws_result["obs_null_ratio"]

        split_results.append({
            "split": split_i + 1,
            "seed": split_seed,
            "n_half1_donors": len(half1_donors),
            "n_half2_donors": len(half2_donors),
            "n_types_ws": ws_result["n_types"],
            "ws_obs_distance": ws_result["observed_distance"],
            "ws_null_mean": ws_result["null_mean"],
            "ws_obs_null_ratio": ws_result["obs_null_ratio"],
            "ws_p_value": ws_result["p_value"],
            "n_types_cs": cs_result["n_types"],
            "cs_obs_distance": cs_result["observed_distance"],
            "cs_null_mean": cs_result["null_mean"],
            "cs_obs_null_ratio": cs_result["obs_null_ratio"],
            "cs_p_value": cs_result["p_value"],
            "delta_obs_null": delta,
        })

        if (split_i + 1) % 10 == 0 or split_i == 0:
            print(f"  Split {split_i + 1:>3}: "
                  f"ws={ws_result['obs_null_ratio']:.3f} (p={ws_result['p_value']:.4f}, n={ws_result['n_types']}), "
                  f"cs={cs_result['obs_null_ratio']:.3f} (p={cs_result['p_value']:.4f}, n={cs_result['n_types']}), "
                  f"delta={delta:+.3f}")

    # ------------------------------------------------------------------
    # Analyze results
    # ------------------------------------------------------------------
    print(f"\n\n[3] RESULTS ({len(split_results)} successful splits)")
    print("=" * 80)

    if len(split_results) == 0:
        print("  ERROR: No successful splits. Cannot proceed.")
        return

    df = pd.DataFrame(split_results)

    # Within-species statistics
    ws_ratios = df["ws_obs_null_ratio"].values
    cs_ratios = df["cs_obs_null_ratio"].values
    deltas = df["delta_obs_null"].values

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

    # Paired t-test
    t_stat, t_pval = stats.ttest_rel(cs_ratios, ws_ratios)
    print(f"    Paired t-test: t={t_stat:.3f}, p={t_pval:.2e}")

    # Wilcoxon signed-rank test
    w_stat, w_pval = stats.wilcoxon(deltas)
    print(f"    Wilcoxon signed-rank: W={w_stat:.1f}, p={w_pval:.2e}")

    consistently_positive = n_positive / len(deltas) > 0.95
    print(f"\n  INTERPRETATION:")
    if consistently_positive and np.percentile(deltas, 2.5) > 0:
        print(f"    Delta is consistently positive (95% CI excludes 0).")
        print(f"    Cross-species is LESS COHERENT than within-species.")
        print(f"    => EVOLUTIONARY SIGNAL EXISTS beyond within-species variation.")
    elif n_positive / len(deltas) > 0.5:
        print(f"    Delta is mostly positive ({100*n_positive/len(deltas):.0f}% of splits).")
        print(f"    Evidence for evolutionary signal, but some overlap.")
    else:
        print(f"    Delta is not consistently positive.")
        print(f"    Within-species variation may explain the cross-species signal.")

    # ------------------------------------------------------------------
    # Save results
    # ------------------------------------------------------------------
    df.to_csv(OUTPUT_DIR / "donor_split_distributions.csv", index=False)

    summary = {
        "n_splits": len(df),
        "n_splits_attempted": N_SPLITS,
        "n_permutations_per_split": N_PERMUTATIONS_PER_SPLIT,
        "min_cells_per_half": MIN_CELLS_PER_HALF,
        "random_seed": RANDOM_SEED,
        "infeasible_types": INFEASIBLE_TYPES,
        "within_species": {
            "mean_obs_null_ratio": float(np.mean(ws_ratios)),
            "median_obs_null_ratio": float(np.median(ws_ratios)),
            "std_obs_null_ratio": float(np.std(ws_ratios)),
            "ci_95": [float(np.percentile(ws_ratios, 2.5)),
                      float(np.percentile(ws_ratios, 97.5))],
            "pct_significant_p001": float(ws_sig / len(df)),
            "median_n_types": float(df["n_types_ws"].median()),
        },
        "cross_species_matched": {
            "mean_obs_null_ratio": float(np.mean(cs_ratios)),
            "median_obs_null_ratio": float(np.median(cs_ratios)),
            "std_obs_null_ratio": float(np.std(cs_ratios)),
            "ci_95": [float(np.percentile(cs_ratios, 2.5)),
                      float(np.percentile(cs_ratios, 97.5))],
            "pct_significant_p001": float(cs_sig / len(df)),
            "median_n_types": float(df["n_types_cs"].median()),
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

    with open(OUTPUT_DIR / "donor_split_results.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\n  Results saved to: {OUTPUT_DIR}/donor_split_results.json")
    print(f"  Distributions saved to: {OUTPUT_DIR}/donor_split_distributions.csv")

    # ------------------------------------------------------------------
    # Figure: three-way comparison
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
    ax.set_title("A. Within-species\n(donor-split)", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)

    # Panel B: Cross-species histogram
    ax = axes[1]
    ax.hist(cs_ratios, bins=25, color="#DD8452", alpha=0.8, edgecolor="white")
    ax.axvline(np.median(cs_ratios), color="red", linestyle="--", linewidth=1.5,
               label=f"Median = {np.median(cs_ratios):.3f}")
    ax.axvline(CROSS_SPECIES_OBS_NULL_PRIMARY, color="black", linestyle=":",
               linewidth=1.5, label=f"Primary (35-type) = {CROSS_SPECIES_OBS_NULL_PRIMARY:.3f}")
    ax.set_xlabel("obs/null ratio", fontsize=11)
    ax.set_title("B. Cross-species\n(matched subsets)", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)

    # Panel C: Three-way comparison (violin/box)
    ax = axes[2]
    positions = [1, 2, 3]
    colors = ["#55A868", "#4C72B0", "#DD8452"]
    labels = ["Self-comparison\n(subsampling)", "Donor-split\n(within-species)", "Cross-species\n(human vs mouse)"]

    # Self-comparison is a point estimate
    ax.scatter([1], [SELF_COMPARISON_OBS_NULL], s=120, color=colors[0],
               zorder=5, edgecolors="black", linewidth=1)

    # Donor-split as boxplot
    bp1 = ax.boxplot([ws_ratios], positions=[2], widths=0.5,
                     patch_artist=True, showfliers=False)
    bp1["boxes"][0].set_facecolor(colors[1])
    bp1["boxes"][0].set_alpha(0.7)
    bp1["medians"][0].set_color("black")

    # Cross-species as boxplot
    bp2 = ax.boxplot([cs_ratios], positions=[3], widths=0.5,
                     patch_artist=True, showfliers=False)
    bp2["boxes"][0].set_facecolor(colors[2])
    bp2["boxes"][0].set_alpha(0.7)
    bp2["medians"][0].set_color("black")

    # Permutation null reference line
    ax.axhline(1.0, color="gray", linestyle=":", linewidth=1,
               label="Permutation null (1.0)")

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("obs/null ratio", fontsize=11)
    ax.set_title("C. Coherence hierarchy", fontsize=12, fontweight="bold")
    ax.legend(fontsize=8, loc="upper left")
    ax.set_ylim(bottom=-0.05)

    plt.tight_layout()
    fig_path = OUTPUT_DIR / "donor_split_figure.png"
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Figure saved to: {fig_path}")

    # ------------------------------------------------------------------
    # Text description of figure
    # ------------------------------------------------------------------
    print(f"\n  Figure description:")
    print(f"    Panel A: Histogram of within-species (donor-split) obs/null ratios across")
    print(f"    {len(df)} random splits. Median = {np.median(ws_ratios):.3f}.")
    print(f"    Panel B: Histogram of cross-species obs/null ratios (human half1 vs mouse)")
    print(f"    restricted to the same cell type subsets. Median = {np.median(cs_ratios):.3f}.")
    print(f"    Panel C: Three-way comparison showing coherence hierarchy:")
    print(f"    self-comparison ({SELF_COMPARISON_OBS_NULL}) < donor-split ({np.median(ws_ratios):.3f}) < cross-species ({np.median(cs_ratios):.3f}) < null (1.0).")

    elapsed = time.time() - t0
    print(f"\n  Total runtime: {elapsed:.1f} seconds ({elapsed/60:.1f} minutes)")


if __name__ == "__main__":
    main()
