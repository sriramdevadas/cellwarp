#!/usr/bin/env python3
"""
CellWarp fast-path reproduction demo — no download required (takes a few minutes).

Loads the deposited PCA centroids and reproduces the PRIMARY headline result:
the 35-type human-mouse Procrustes geometric-coherence obs/null ratio (median
denominator, the canonical convention) and its 1,000,000-permutation p-value.
Prints PASS/FAIL against the published values (obs/null = 0.522, p < 1e-6).

This is the turnkey, one-command, no-network demonstration. It reuses the exact
canonical pipeline functions (cellwarp.procrustes) and the deposited centroids,
so it requires no atlas downloads, no QC, and no re-clustering. It is read-only:
it does NOT overwrite any committed analysis output.

Input  : output/phase2/scaled_35types/pca_centroids_35.npz (deposited)
Usage  : python reproduce/fast_path.py
Runtime: a few minutes (1,000,000 permutations, seed 42)
"""
from __future__ import annotations

import sys
import time
import warnings
from pathlib import Path

import numpy as np

# Suppress the harmless det() RuntimeWarnings on the 33x33 rotation matrix. It is
# orthogonal, so its determinant is ±1 and nothing overflows; the warnings are
# backend-dependent FPU status flags and the sign is correct regardless
# (matches scripts/permutation_1M.py).
warnings.filterwarnings("ignore", message=".*encountered in det.*")

REPO_ROOT = Path(__file__).resolve().parent.parent
# Make cellwarp importable whether or not it was pip-installed.
sys.path.insert(0, str(REPO_ROOT / "src"))
from cellwarp.procrustes import _procrustes_distance, permutation_test  # noqa: E402

CENTROIDS = REPO_ROOT / "output" / "phase2" / "scaled_35types" / "pca_centroids_35.npz"
PUBLISHED_OBS_NULL = 0.522
OBS_NULL_TOL = 0.005
N_PERMUTATIONS = 1_000_000
SEED = 42


def main() -> int:
    print("=" * 68)
    print("  CellWarp fast-path reproduction (no download; the permutation step takes a few minutes)")
    print("=" * 68)

    if not CENTROIDS.exists():
        print(f"  ERROR: deposited centroids not found at {CENTROIDS}")
        return 1

    data = np.load(CENTROIDS)
    human, mouse = data["human"], data["mouse"]
    n_types, n_pcs = human.shape
    print(f"  Loaded deposited centroids: {n_types} cell types x {n_pcs} PCs")
    print(f"    ({CENTROIDS.relative_to(REPO_ROOT)})")

    observed = _procrustes_distance(human, mouse)
    print(f"\n  Observed Procrustes distance: {observed:.3f}")
    print(f"  Running {N_PERMUTATIONS:,} label permutations (seed={SEED}) —"
          f" this can take a few minutes, with no output until it finishes.")

    t0 = time.time()
    p_value, null = permutation_test(human, mouse, n_permutations=N_PERMUTATIONS, seed=SEED)
    runtime = time.time() - t0

    null_median = float(np.median(null))
    obs_null = observed / null_median

    obs_ok = abs(obs_null - PUBLISHED_OBS_NULL) <= OBS_NULL_TOL
    p_ok = p_value < 1e-6

    print("\n" + "-" * 68)
    print("  RESULT (median denominator = canonical convention)")
    print("-" * 68)
    print(f"  obs/null = {obs_null:.4f}   published {PUBLISHED_OBS_NULL:.3f}"
          f"   -> {'PASS' if obs_ok else 'FAIL'} (tol +/-{OBS_NULL_TOL})")
    print(f"  p-value  = {p_value:.2e}   published p < 1e-6"
          f"        -> {'PASS' if p_ok else 'FAIL'}")
    print(f"  runtime  = {runtime:.1f}s")
    print("=" * 68)

    if obs_ok and p_ok:
        print("  PASS - headline reproduced from deposited centroids, no download.")
        print("=" * 68)
        return 0
    print("  FAIL - reproduced values do not match published headline.")
    print("=" * 68)
    return 1


if __name__ == "__main__":
    sys.exit(main())
