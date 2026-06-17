#!/usr/bin/env python3
"""
CellWarp — 1M Permutation Test for Independent-PCA Sensitivity Analysis (C1)

Runs the Procrustes permutation test with 1,000,000 label permutations on the
independent-PCA aligned centroids, upgrading the C1 analysis from 10K to 1M
permutations (same approach as C10 in analysis/permutation_1M/).

Biology & Math: Same null model as the primary analysis. H0: the correspondence
between human and mouse cell types is no better than random. For each permutation,
shuffle which mouse centroid maps to which human centroid and compute Procrustes
distance. p-value = (# null ≤ observed + 1) / (B + 1).

Input:  analysis/independent_pca_sensitivity/pca_centroids_independent.npz
Output: analysis/independent_pca_sensitivity/ (updates results JSON + null dist)
"""

from __future__ import annotations

import json
import time
import warnings
from pathlib import Path

import numpy as np
from numpy.linalg import svd

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
N_PERMUTATIONS = 1_000_000
RANDOM_SEED = 42
PROJECT_ROOT = Path(__file__).resolve().parent.parent
INPUT_DIR = PROJECT_ROOT / "analysis" / "independent_pca_sensitivity"
OUTPUT_DIR = INPUT_DIR  # update in place


# ---------------------------------------------------------------------------
# Procrustes distance (self-contained, same math as src/procrustes.py)
# ---------------------------------------------------------------------------
def _procrustes_distance(X: np.ndarray, Y: np.ndarray) -> float:
    """Procrustes distance without reflection, no printing.

    Uses sign from SVD of V @ U.T to avoid overflow in np.linalg.det
    for high-dimensional matrices.
    """
    n, k = X.shape
    X_c = X - X.mean(axis=0)
    Y_c = Y - Y.mean(axis=0)

    M = X_c.T @ Y_c
    U, sigma, Vt = svd(M)
    V = Vt.T

    # Determine sign of det(V @ U.T) via SVD to avoid overflow on large k
    VUt = V @ U.T
    det_val = np.linalg.det(VUt)
    if np.isfinite(det_val):
        sign = np.sign(det_val)
    else:
        eigvals = np.linalg.eigvals(VUt)
        sign = np.sign(np.real(np.prod(eigvals)))

    D_diag = np.ones(k)
    D_diag[-1] = sign

    ss_Y = np.sum(Y_c**2)
    trace_sigma_D = np.sum(sigma * D_diag)
    s = trace_sigma_D / ss_Y

    Y_aligned = s * (Y_c @ (V * D_diag) @ U.T)
    return np.sqrt(np.sum((X_c - Y_aligned) ** 2))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print(f"INDEPENDENT-PCA PERMUTATION TEST — {N_PERMUTATIONS:,} permutations")
    print("  C1 SENSITIVITY ANALYSIS: 35 cell types × 33 PCs (independent PCA)")
    print("=" * 70)

    # Load pre-computed independent-PCA aligned centroids
    npz_path = INPUT_DIR / "pca_centroids_independent.npz"
    data = np.load(npz_path, allow_pickle=True)
    human_pca = data["human"]
    mouse_aligned = data["mouse_aligned"]
    cell_types = data["cell_types"].tolist()

    n_types, n_pcs = human_pca.shape
    print(f"\n  Loaded independent-PCA centroids: {n_types} cell types × {n_pcs} PCs")
    print(f"  Cell types ({n_types}): {cell_types[:5]} ... {cell_types[-3:]}")

    assert n_types == 35, f"Expected 35 cell types, got {n_types}"
    assert n_pcs == 33, f"Expected 33 PCs, got {n_pcs}"

    # Observed Procrustes distance
    observed = _procrustes_distance(human_pca, mouse_aligned)
    print(f"\n  Observed Procrustes distance: {observed:.3f} (expected ~52.716)")
    assert abs(observed - 52.716) < 0.1, (
        f"Observed distance {observed:.3f} does not match expected 52.716"
    )

    # Permutation test
    print(f"\n  Running {N_PERMUTATIONS:,} permutations (seed={RANDOM_SEED})...")
    rng = np.random.RandomState(RANDOM_SEED)
    null_distances = np.zeros(N_PERMUTATIONS)

    t_start = time.time()
    for i in range(N_PERMUTATIONS):
        perm = rng.permutation(n_types)
        null_distances[i] = _procrustes_distance(human_pca, mouse_aligned[perm])
        if (i + 1) % 200_000 == 0:
            elapsed = time.time() - t_start
            rate = (i + 1) / elapsed
            eta = (N_PERMUTATIONS - i - 1) / rate
            print(f"    {i + 1:>10,} / {N_PERMUTATIONS:,}  "
                  f"({elapsed:.1f}s elapsed, ~{eta:.1f}s remaining)")

    runtime = time.time() - t_start

    # p-value (conservative: +1 includes observed itself)
    n_leq = int(np.sum(null_distances <= observed))
    p_value = (n_leq + 1) / (N_PERMUTATIONS + 1)

    null_mean = float(np.mean(null_distances))
    null_median = float(np.median(null_distances))
    null_std = float(np.std(null_distances))
    obs_null_ratio = observed / null_median

    # Results
    print(f"\n  {'=' * 50}")
    print(f"  RESULTS — independent PCA, {N_PERMUTATIONS:,} permutations")
    print(f"  {'=' * 50}")
    print(f"  Observed Procrustes distance: {observed:.3f}")
    print(f"  obs/null ratio: {obs_null_ratio:.4f}")
    print(f"  Permutations: {N_PERMUTATIONS:,}")
    print(f"  Null distribution: mean={null_mean:.3f}, median={null_median:.3f}")
    print(f"  Null std: {null_std:.3f}")
    print(f"  Null range: [{np.min(null_distances):.3f}, "
          f"{np.max(null_distances):.3f}]")
    print(f"  Null 5th percentile: {np.percentile(null_distances, 5):.3f}")
    print(f"  Null 1st percentile: {np.percentile(null_distances, 1):.3f}")
    print(f"  Null 0.1th percentile: {np.percentile(null_distances, 0.1):.3f}")
    print(f"  Null 0.01th percentile: {np.percentile(null_distances, 0.01):.3f}")
    print(f"  Permuted distances ≤ observed: {n_leq:,} / {N_PERMUTATIONS:,}")

    if p_value < 1e-6:
        p_str = f"p < 10⁻⁶ (exact: {p_value:.2e})"
    else:
        p_str = f"p = {p_value:.6e}"
    print(f"  p-value: {p_str}")
    print(f"  Runtime: {runtime:.1f}s ({runtime / 60:.2f} min)")

    # Save null distribution
    np.save(OUTPUT_DIR / "null_distribution_independent_pca.npy", null_distances)
    print(f"\n  Saved null distribution ({N_PERMUTATIONS:,} values)")

    # -----------------------------------------------------------------------
    # Update existing results JSON — preserve everything, update permutation
    # -----------------------------------------------------------------------
    results_path = OUTPUT_DIR / "independent_pca_results.json"
    with open(results_path) as f:
        results = json.load(f)

    # Store prior 10K result for reference
    prior_perm = results["permutation_test"].copy()

    # Update config
    results["config"]["n_permutations"] = N_PERMUTATIONS

    # Update permutation test section
    results["permutation_test"] = {
        "observed_distance": float(observed),
        "p_value": float(p_value),
        "p_value_str": p_str,
        "n_leq_observed": n_leq,
        "n_permutations": N_PERMUTATIONS,
        "null_distribution_summary": {
            "mean": null_mean,
            "median": null_median,
            "std": null_std,
            "min": float(np.min(null_distances)),
            "max": float(np.max(null_distances)),
            "percentile_0_01": float(np.percentile(null_distances, 0.01)),
            "percentile_0_1": float(np.percentile(null_distances, 0.1)),
            "percentile_1": float(np.percentile(null_distances, 1)),
            "percentile_2_5": float(np.percentile(null_distances, 2.5)),
            "percentile_5": float(np.percentile(null_distances, 5)),
            "percentile_95": float(np.percentile(null_distances, 95)),
            "percentile_97_5": float(np.percentile(null_distances, 97.5)),
            "percentile_99": float(np.percentile(null_distances, 99)),
            "percentile_99_9": float(np.percentile(null_distances, 99.9)),
        },
        "obs_null_ratio": float(obs_null_ratio),
        "significant_at_001": bool(p_value < 0.01),
        "prior_10k_result": {
            "n_permutations": prior_perm["n_permutations"],
            "p_value": prior_perm["p_value"],
            "n_leq_observed": prior_perm["n_leq_observed"],
            "obs_null_ratio": prior_perm["obs_null_ratio"],
        },
    }

    # Update comparison section's independent_pca p-value
    results["comparison_to_joint_pca"]["independent_pca"]["p_value"] = float(p_value)

    # Save updated JSON
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"  Updated: {results_path.name}")

    # Summary comparison
    print(f"\n  {'=' * 50}")
    print(f"  10K → 1M COMPARISON")
    print(f"  {'=' * 50}")
    print(f"  {'Metric':<30} {'10K':>15} {'1M':>15}")
    print(f"  {'-'*62}")
    print(f"  {'Null ≤ observed':<30} {'0 / 10,000':>15} {f'{n_leq:,} / 1,000,000':>15}")
    print(f"  {'p-value':<30} {'< 1e-4':>15} {p_str:>15}")
    print(f"  {'Null mean':<30} {prior_perm['null_distribution_summary']['mean']:>15.3f} {null_mean:>15.3f}")
    print(f"  {'Null min':<30} {prior_perm['null_distribution_summary']['min']:>15.3f} {np.min(null_distances):>15.3f}")
    print(f"  {'Gap (null min - obs)':<30} {prior_perm['null_distribution_summary']['min'] - observed:>15.3f} {np.min(null_distances) - observed:>15.3f}")
    print()


if __name__ == "__main__":
    main()
