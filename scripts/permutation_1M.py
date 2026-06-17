#!/usr/bin/env python3
"""
CellWarp — 1M Permutation Test

Runs the PRIMARY Procrustes permutation test with 1,000,000 label permutations
to resolve the p-value beyond the 10K floor (p=0.0001). Uses the 35 cell type ×
33 PC dataset from the scaled analysis (script 08).

Biology & Math: Same as src/procrustes.py permutation_test(). Null hypothesis
H0: the correspondence between human and mouse cell types is no better than
random. For each permutation, shuffle which mouse centroid maps to which human
centroid and compute Procrustes distance. p-value = (# null ≤ observed + 1) /
(B + 1).

Input:  output/phase2/scaled_35types/pca_centroids_35.npz (from 08_scaled_procrustes.py)
Output: analysis/permutation_1M/
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
OUTPUT_DIR = PROJECT_ROOT / "analysis" / "permutation_1M"

# Suppress overflow warnings from det() on 33×33 matrices (harmless —
# we only need the sign, which is correctly computed from the SVD)
warnings.filterwarnings("ignore", message=".*encountered in det.*")


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
    # det(V @ U.T) = det(V) * det(U.T) = product of all singular value signs
    # But since U, V are from SVD they are orthogonal, det = ±1
    # Safest: compute det directly but handle overflow by falling back to sign
    VUt = V @ U.T
    det_val = np.linalg.det(VUt)
    if np.isfinite(det_val):
        sign = np.sign(det_val)
    else:
        # Fallback: det of orthogonal matrix is ±1; use SVD of VUt
        _, s_vu, _ = svd(VUt)
        sign = np.sign(np.prod(s_vu)) * np.sign(det_val) if det_val != 0 else 1.0
        # For orthogonal matrices, use the product of diagonal of the
        # Schur decomposition, but simplest: check if reflection
        sign = -1.0 if np.sum(np.log(np.abs(s_vu))) < -30 else 1.0
        # Actually for orthogonal matrices all singular values = 1,
        # so product = 1. The sign comes from whether it's a proper rotation.
        # Use the reliable method: sign = sign of product of eigenvalues
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
    print(f"PROCRUSTES PERMUTATION TEST — {N_PERMUTATIONS:,} permutations")
    print("  PRIMARY ANALYSIS: 35 cell types × 33 PCs")
    print("=" * 70)

    # Load PCA centroids from scaled 35-type analysis
    pca_path = PROJECT_ROOT / "output" / "phase2" / "scaled_35types" / "pca_centroids_35.npz"
    data = np.load(pca_path)
    human_pca = data["human"]
    mouse_pca = data["mouse"]
    cell_types = data["cell_types"].tolist()

    n_types, n_pcs = human_pca.shape
    print(f"\n  Loaded PCA centroids: {n_types} cell types × {n_pcs} PCs")
    print(f"  Cell types ({n_types}): {cell_types[:5]} ... {cell_types[-3:]}")

    # Verify dimensions match expected primary analysis
    assert n_types == 35, f"Expected 35 cell types, got {n_types}"
    assert n_pcs == 33, f"Expected 33 PCs, got {n_pcs}"

    # Observed Procrustes distance
    observed = _procrustes_distance(human_pca, mouse_pca)
    print(f"\n  Observed Procrustes distance: {observed:.3f} (expected ~61.153)")
    assert abs(observed - 61.153) < 0.1, (
        f"Observed distance {observed:.3f} does not match expected 61.153"
    )

    # Permutation test
    print(f"\n  Running {N_PERMUTATIONS:,} permutations (seed={RANDOM_SEED})...")
    rng = np.random.RandomState(RANDOM_SEED)
    null_distances = np.zeros(N_PERMUTATIONS)

    t_start = time.time()
    for i in range(N_PERMUTATIONS):
        perm = rng.permutation(n_types)
        null_distances[i] = _procrustes_distance(human_pca, mouse_pca[perm])
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

    # Results
    print(f"\n  {'=' * 50}")
    print(f"  RESULTS — 35 cell types, {N_PERMUTATIONS:,} permutations")
    print(f"  {'=' * 50}")
    print(f"  Observed Procrustes distance: {observed:.3f}")
    print(f"  obs/null ratio: {observed / np.median(null_distances):.3f} (expected ~0.522)")
    print(f"  Permutations: {N_PERMUTATIONS:,}")
    print(f"  Null distribution: mean={np.mean(null_distances):.3f}, "
          f"median={np.median(null_distances):.3f}")
    print(f"  Null std: {np.std(null_distances):.3f}")
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
    np.save(OUTPUT_DIR / "null_distribution_1M.npy", null_distances)
    print(f"\n  Saved null distribution: {OUTPUT_DIR / 'null_distribution_1M.npy'}")

    # Save JSON results
    results = {
        "test": "Procrustes permutation test (1M) — PRIMARY 35-type analysis",
        "motivation": "tighter p-value bound at 1,000,000 permutations",
        "observed_procrustes_distance": float(observed),
        "obs_null_ratio": float(observed / np.median(null_distances)),
        "n_permutations": N_PERMUTATIONS,
        "random_seed": RANDOM_SEED,
        "n_cell_types": n_types,
        "n_pca_components": n_pcs,
        "n_null_leq_observed": n_leq,
        "p_value": float(p_value),
        "p_value_str": p_str,
        "null_distribution_summary": {
            "mean": float(np.mean(null_distances)),
            "median": float(np.median(null_distances)),
            "std": float(np.std(null_distances)),
            "min": float(np.min(null_distances)),
            "max": float(np.max(null_distances)),
            "percentile_0_01": float(np.percentile(null_distances, 0.01)),
            "percentile_0_1": float(np.percentile(null_distances, 0.1)),
            "percentile_1": float(np.percentile(null_distances, 1)),
            "percentile_5": float(np.percentile(null_distances, 5)),
            "percentile_95": float(np.percentile(null_distances, 95)),
            "percentile_99": float(np.percentile(null_distances, 99)),
            "percentile_99_9": float(np.percentile(null_distances, 99.9)),
        },
        "runtime_seconds": round(runtime, 1),
        "prior_10k_result": {
            "n_permutations": 10_000,
            "p_value": 0.0001,
            "obs_null_ratio": 0.522,
            "observed_distance": 61.153,
            "null_median": 117.105,
            "note": "From output/phase2/scaled_35types/procrustes_results_35.json",
        },
    }

    with open(OUTPUT_DIR / "results_1M.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"  Saved JSON results: {OUTPUT_DIR / 'results_1M.json'}")


if __name__ == "__main__":
    main()
