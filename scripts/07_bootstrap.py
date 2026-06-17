#!/usr/bin/env python3
"""
CellWarp — Bootstrap Robustness Test (Phase 3)

Tests the stability of the 35-cell-type Procrustes result by resampling cells
within each cell type and re-running the full pipeline 100 times.

Biology
-------
If the Procrustes alignment is robust, the estimated distance should be stable
across subsamples. High variability would suggest the geometric signal depends
on specific cells rather than reflecting a general property of cell type
expression programs. The coefficient of variation (CV = std/mean) must be < 0.2
for Phase 3 gate passage.

Math
----
For each of B=100 bootstrap iterations:
  1. For each cell type t in each species s: draw 50% of cells uniformly
     at random (with seed = iteration number for reproducibility).
  2. Recompute centroid: μ_t^(b) = (1/n_b) Σ x_i  (n_b = ⌊n/2⌋ cells)
  3. Run PCA on combined centroids (95% variance threshold).
  4. Run Procrustes alignment (mouse → human).
  5. Run permutation test (1,000 iterations for speed).
  6. Record: Procrustes distance d^(b), p-value p^(b).

Gate criterion: CV(d^(1), ..., d^(B)) = std/mean < 0.2

Inputs:
    data/phase2_scaled/human_scaled.h5ad    — normalized 35-type human data
    data/phase2_scaled/mouse_scaled.h5ad    — normalized 35-type mouse data

Outputs (all in output/phase3/bootstrap/):
    bootstrap_results.csv       — distance and p-value per iteration
    bootstrap_summary.json      — mean, std, CV, gate result
    bootstrap_histogram.png     — histogram of distances with mean ± 1 SD

Usage:
    python scripts/07_bootstrap.py
"""

from __future__ import annotations

import contextlib
import io
import json
import sys
import time
import warnings
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cellwarp.procrustes import (
    PCA_VARIANCE_THRESHOLD,
    _procrustes_distance,
    pca_reduce_centroids,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

N_BOOTSTRAP = 100
SUBSAMPLE_FRACTION = 0.5
N_PERM_BOOTSTRAP = 1_000  # Reduced from 10,000 for speed
CELL_TYPE_COL = "cell_type"

# NOTE: Task description says data/phase1/ but those files contain only 6 types.
# The 35-type normalized data is in data/phase2_scaled/.
HUMAN_DATA = Path("data/phase2_scaled/human_scaled.h5ad")
MOUSE_DATA = Path("data/phase2_scaled/mouse_scaled.h5ad")
OUTPUT_DIR = Path("output/phase3/bootstrap")

CV_THRESHOLD = 0.2  # Gate criterion


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def subsample_cells(
    adata: ad.AnnData, fraction: float, seed: int
) -> ad.AnnData:
    """
    Subsample a fraction of cells per cell type with a given random seed.

    Biology: Subsampling tests how sensitive the centroid estimates (and
    downstream Procrustes geometry) are to the specific cells drawn from
    each cell type population. With 50% subsampling, we simulate having
    access to a different biological sample of the same cell types.

    Args:
        adata: Annotated data with cell_type column in .obs.
        fraction: Fraction of cells to retain per type (0-1).
        seed: Random seed for reproducibility.

    Returns:
        AnnData with subsampled cells (copy, does not modify original).
    """
    rng = np.random.RandomState(seed)
    indices = []
    for ct in sorted(adata.obs[CELL_TYPE_COL].unique()):
        ct_idx = np.where(adata.obs[CELL_TYPE_COL].values == ct)[0]
        n_keep = max(1, int(len(ct_idx) * fraction))
        chosen = rng.choice(ct_idx, size=n_keep, replace=False)
        indices.extend(chosen)
    indices = sorted(indices)
    return adata[indices].copy()


def compute_centroids_silent(adata: ad.AnnData) -> pd.DataFrame:
    """
    Compute mean expression per cell type (silent — no printing).

    Math: For cell type t with cells {c_1, ..., c_n}:
        μ_t = (1/n) Σ_{i=1}^{n} x_{c_i}

    Args:
        adata: AnnData with log-normalized expression in .X (may be sparse).

    Returns:
        DataFrame of shape (n_types, n_genes), indexed by cell type name.
    """
    cell_types = sorted(adata.obs[CELL_TYPE_COL].unique())
    centroids = {}
    for ct in cell_types:
        mask = adata.obs[CELL_TYPE_COL] == ct
        mean_vec = np.asarray(adata[mask].X.mean(axis=0)).flatten()
        centroids[ct] = mean_vec
    df = pd.DataFrame(centroids, index=adata.var_names.tolist()).T
    df.index.name = "cell_type"
    return df


def silent_permutation_test(
    X: np.ndarray, Y: np.ndarray, n_perm: int, seed: int
) -> tuple[float, float]:
    """
    Silent permutation test — returns (distance, p_value) only.

    Math: Same as cellwarp.procrustes.permutation_test but without printing
    and returns only the scalar summaries needed for bootstrap recording.

    Args:
        X: Reference points (n × k). Human centroids in PCA space.
        Y: Target points (n × k). Mouse centroids in PCA space.
        n_perm: Number of permutations.
        seed: Random seed.

    Returns:
        Tuple of (observed_distance, p_value).
    """
    n = X.shape[0]
    rng = np.random.RandomState(seed)
    observed = _procrustes_distance(X, Y)

    null_distances = np.zeros(n_perm)
    for i in range(n_perm):
        perm = rng.permutation(n)
        null_distances[i] = _procrustes_distance(X, Y[perm])

    n_leq = int(np.sum(null_distances <= observed))
    p_value = (n_leq + 1) / (n_perm + 1)
    return observed, p_value


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    t_start = time.time()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("PHASE 3 — BOOTSTRAP ROBUSTNESS TEST (35 CELL TYPES)")
    print("=" * 70)

    # ---- Load data ----
    print("\n  Loading 35-type normalized data...")
    human = ad.read_h5ad(HUMAN_DATA)
    mouse = ad.read_h5ad(MOUSE_DATA)
    print(f"  Human: {human.n_obs:,} cells × {human.n_vars:,} genes")
    print(f"  Mouse: {mouse.n_obs:,} cells × {mouse.n_vars:,} genes")
    n_types = len(human.obs[CELL_TYPE_COL].unique())
    print(f"  Cell types: {n_types}")

    # ---- Bootstrap loop ----
    print(
        f"\n  Running {N_BOOTSTRAP} bootstrap iterations "
        f"({SUBSAMPLE_FRACTION * 100:.0f}% subsample, "
        f"{N_PERM_BOOTSTRAP:,} permutations each)..."
    )

    results = []

    # Suppress numpy overflow warnings from high-dimensional det() computation
    # (ISSUE-022: 33-dim determinant overflows float64, sign is still correct)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)

        for i in range(N_BOOTSTRAP):
            t_iter = time.time()

            # 1. Subsample 50% of cells per type per species
            human_sub = subsample_cells(human, SUBSAMPLE_FRACTION, seed=i)
            mouse_sub = subsample_cells(mouse, SUBSAMPLE_FRACTION, seed=i)

            # 2. Compute centroids from subsampled cells
            h_cent = compute_centroids_silent(human_sub)
            m_cent = compute_centroids_silent(mouse_sub)

            # Free memory from subsampled AnnData
            del human_sub, mouse_sub

            # 3. PCA on combined centroids (95% variance threshold, silent)
            with contextlib.redirect_stdout(io.StringIO()):
                h_pca, m_pca, pca_model, cell_types = pca_reduce_centroids(
                    h_cent, m_cent, PCA_VARIANCE_THRESHOLD
                )

            # 4 + 5. Procrustes + permutation test (silent)
            # Use seed=i+1000 to avoid correlation with subsampling seed
            distance, p_value = silent_permutation_test(
                h_pca, m_pca, N_PERM_BOOTSTRAP, seed=i + 1000
            )

            n_pc = pca_model.n_components_
            elapsed = time.time() - t_iter

            results.append({
                "iteration": i,
                "distance": distance,
                "p_value": p_value,
                "n_pca_components": n_pc,
                "runtime_s": elapsed,
            })

            # Progress indicator every 10 iterations
            if (i + 1) % 10 == 0 or i == 0:
                print(
                    f"    [{i + 1:>3}/{N_BOOTSTRAP}] "
                    f"d={distance:.4f}, p={p_value:.4f}, "
                    f"PC={n_pc}, {elapsed:.1f}s"
                )

    # ---- Compute statistics ----
    df = pd.DataFrame(results)
    distances = df["distance"].values
    p_values = df["p_value"].values

    mean_d = float(np.mean(distances))
    std_d = float(np.std(distances))
    cv = std_d / mean_d
    median_d = float(np.median(distances))
    min_d = float(np.min(distances))
    max_d = float(np.max(distances))
    frac_sig = float(np.mean(p_values < 0.01))
    mean_p = float(np.mean(p_values))
    mean_pc = float(np.mean(df["n_pca_components"].values))

    cv_pass = cv < CV_THRESHOLD

    # ---- Save results ----
    df.to_csv(OUTPUT_DIR / "bootstrap_results.csv", index=False)

    summary = {
        "n_bootstrap": N_BOOTSTRAP,
        "subsample_fraction": SUBSAMPLE_FRACTION,
        "n_permutations_per_iteration": N_PERM_BOOTSTRAP,
        "distances": {
            "mean": mean_d,
            "std": std_d,
            "cv": cv,
            "median": median_d,
            "min": min_d,
            "max": max_d,
        },
        "p_values": {
            "mean": mean_p,
            "fraction_significant_001": frac_sig,
        },
        "pca_components": {
            "mean": mean_pc,
        },
        "gate_criterion": {
            "threshold": CV_THRESHOLD,
            "cv": cv,
            "pass": bool(cv_pass),
        },
    }

    with open(OUTPUT_DIR / "bootstrap_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # ---- Plot histogram ----
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(distances, bins=20, color="steelblue", edgecolor="white", alpha=0.8)
    ax.axvline(
        mean_d, color="red", linestyle="-", linewidth=2,
        label=f"Mean = {mean_d:.2f}",
    )
    ax.axvline(
        mean_d - std_d, color="red", linestyle="--", linewidth=1,
        label=f"Mean − 1 SD = {mean_d - std_d:.2f}",
    )
    ax.axvline(
        mean_d + std_d, color="red", linestyle="--", linewidth=1,
        label=f"Mean + 1 SD = {mean_d + std_d:.2f}",
    )
    ax.set_xlabel("Procrustes Distance", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.set_title(
        f"Bootstrap Procrustes Distances (n={N_BOOTSTRAP}, 50% subsample)\n"
        f"CV = {cv:.4f} {'(PASS < 0.2)' if cv_pass else '(FAIL ≥ 0.2)'}",
        fontsize=13,
    )
    ax.legend(fontsize=10)
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "bootstrap_histogram.png", dpi=150)
    plt.close(fig)

    # ---- Summary ----
    total_time = time.time() - t_start

    print("\n" + "=" * 70)
    print("BOOTSTRAP ROBUSTNESS — RESULTS")
    print("=" * 70)
    print(f"\n  Bootstrap iterations: {N_BOOTSTRAP}")
    print(f"  Subsample fraction: {SUBSAMPLE_FRACTION * 100:.0f}%")
    print(f"  Permutations per iteration: {N_PERM_BOOTSTRAP:,}")
    print(f"\n  Procrustes distance distribution:")
    print(f"    Mean:   {mean_d:.4f}")
    print(f"    Std:    {std_d:.4f}")
    print(f"    CV:     {cv:.4f}")
    print(f"    Median: {median_d:.4f}")
    print(f"    Range:  [{min_d:.4f}, {max_d:.4f}]")
    print(f"\n  Permutation p-values:")
    print(f"    Mean p-value: {mean_p:.6f}")
    print(f"    Fraction significant (p < 0.01): {frac_sig * 100:.1f}%")
    print(f"\n  PCA components (mean): {mean_pc:.1f}")
    print(f"\n  GATE CRITERION: CV < {CV_THRESHOLD}")
    print(f"    CV = {cv:.4f}")
    print(f"    Result: {'PASS' if cv_pass else 'FAIL'}")
    print(f"\n  Total runtime: {total_time:.1f}s ({total_time / 60:.1f} min)")
    print(f"  Output: {OUTPUT_DIR}/")

    print(
        f"\n  Plot description: Histogram of {N_BOOTSTRAP} bootstrap "
        f"Procrustes distances."
    )
    if cv_pass:
        print(
            f"  The distribution is narrow (CV={cv:.4f}), centered around "
            f"{mean_d:.2f}."
        )
        print(
            f"  Red lines mark mean ± 1 SD. The low CV indicates the Procrustes"
        )
        print(
            f"  distance is stable across subsamples — the geometric signal "
            f"is robust."
        )
    else:
        print(
            f"  The distribution has CV={cv:.4f}, centered around {mean_d:.2f}."
        )
        print(
            f"  Red lines mark mean ± 1 SD. The high CV indicates the Procrustes"
        )
        print(
            f"  distance varies substantially across subsamples."
        )


if __name__ == "__main__":
    main()
