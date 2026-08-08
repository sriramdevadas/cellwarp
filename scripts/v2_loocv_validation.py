#!/usr/bin/env python3
"""
V2 — LOOCV Implementation: Independent Validation
===================================================
Verifies the leave-one-out cross-validation result used as the empirical
independence test for the Cardini/PCA independence argument.

NO imports from src/ — Procrustes reimplemented from scratch.

Canonical targets:
  Correct recoveries: 35/35
  Mean error-to-null ratio: 0.42 (tolerance ±0.02)

The null is geometric: mean Euclidean distance from predicted position
to each of the 34 training human centroids. NOT a permutation null.

Author: V2 software validation (independent of src/)
Date: 2026-03-21
"""

import json
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import svd
from sklearn.decomposition import PCA

# ── Constants ──────────────────────────────────────────────────────────
RANDOM_SEED = 42
VARIANCE_THRESHOLD = 0.95
MEAN_RATIO_TOLERANCE = 0.02

# ── File paths ─────────────────────────────────────────────────────────
BASE = Path(__file__).resolve().parent.parent
HUMAN_35 = BASE / "output/phase2/scaled_35types/centroids_human_35.csv"
MOUSE_35 = BASE / "output/phase2/scaled_35types/centroids_mouse_35.csv"
OUTPUT_DIR = BASE / "output/validation/v2_loocv"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ======================================================================
# Independent Procrustes implementation (identical to V1)
# ======================================================================

def procrustes_fit(X: np.ndarray, Y: np.ndarray) -> dict:
    """
    Ordinary Procrustes Analysis: align Y to X via rotation + uniform
    scaling (no reflection). Returns rotation R, scaling s, and
    translations for applying to new points.

    Math:
        1. Center both matrices
        2. Cross-covariance M = X_c^T @ Y_c
        3. SVD: M = U Σ V^T
        4. Rotation R = V D U^T  (D corrects for reflection)
        5. Scaling s = tr(Σ D) / ||Y_c||²_F
    """
    n, k = X.shape
    mu_X = X.mean(axis=0)
    mu_Y = Y.mean(axis=0)
    X_c = X - mu_X
    Y_c = Y - mu_Y

    M = X_c.T @ Y_c
    U, sigma, Vt = svd(M)
    V = Vt.T

    # Reflection correction: ensure proper rotation (det = +1)
    d = np.linalg.det(V @ U.T)
    D_diag = np.ones(k)
    D_diag[-1] = np.sign(d)

    # Optimal rotation
    R = (V * D_diag) @ U.T

    # Optimal scaling
    ss_Y = np.sum(Y_c ** 2)
    s = np.sum(sigma * D_diag) / ss_Y

    return {
        "rotation": R,
        "scaling": s,
        "translation_ref": mu_X,
        "translation_target": mu_Y,
    }


# ======================================================================
# LOOCV procedure
# ======================================================================

def run_loocv(
    human_centroids: pd.DataFrame,
    mouse_centroids: pd.DataFrame,
) -> list[dict]:
    """
    Leave-one-out cross-validation on cell types.

    For each of 35 types held out:
      1. Remove held-out type from both species
      2. PCA on 68 training centroids (95% variance)
      3. Procrustes on 34 training pairs
      4. Project held-out mouse centroid, apply learned transform
      5. error = ||predicted - actual_human||
      6. null = mean(||predicted - each training human||)
      7. ratio = error / null; correct if ratio < 1.0
    """
    cell_types = sorted(human_centroids.index.tolist())
    assert sorted(mouse_centroids.index.tolist()) == cell_types
    n_types = len(cell_types)

    results = []

    # Suppress the RuntimeWarnings from det() on the orthogonal rotation matrix.
    # Nothing overflows -- det = ±1; see scripts/07_bootstrap.py.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)

        for i, held_out in enumerate(cell_types):
            # 1. Remove held-out type
            train_types = [ct for ct in cell_types if ct != held_out]
            h_train = human_centroids.loc[train_types].values  # (34, G)
            m_train = mouse_centroids.loc[train_types].values  # (34, G)

            # 2. PCA on 68 training centroids
            combined = np.vstack([h_train, m_train])  # (68, G)
            pca = PCA(
                n_components=VARIANCE_THRESHOLD,
                svd_solver="full",
                random_state=RANDOM_SEED,
            )
            combined_pca = pca.fit_transform(combined)
            n_train = len(train_types)
            h_pca = combined_pca[:n_train]   # (34, k)
            m_pca = combined_pca[n_train:]   # (34, k)
            n_components = pca.n_components_

            # 3. Procrustes on training pairs
            proc = procrustes_fit(h_pca, m_pca)

            # 4. Project held-out centroids into training PCA space
            held_h_gene = human_centroids.loc[held_out].values.reshape(1, -1)
            held_m_gene = mouse_centroids.loc[held_out].values.reshape(1, -1)
            held_h_pca = pca.transform(held_h_gene).flatten()  # (k,)
            held_m_pca = pca.transform(held_m_gene).flatten()  # (k,)

            # 5. Apply learned transform to held-out mouse centroid
            # predicted = s * (held_m - μ_Y) @ R + μ_X
            held_m_centered = held_m_pca - proc["translation_target"]
            predicted = (
                proc["scaling"] * (held_m_centered @ proc["rotation"])
                + proc["translation_ref"]
            )

            # 6. Error = Euclidean distance from predicted to actual human
            error = float(np.linalg.norm(predicted - held_h_pca))

            # 7. Null = mean distance from predicted to all 34 training
            #    human centroids (geometric comparison, not permutation)
            null_distances = np.array([
                np.linalg.norm(predicted - h_pca[j])
                for j in range(n_train)
            ])
            null_dist = float(np.mean(null_distances))

            # Ratio
            ratio = error / null_dist if null_dist > 0 else float("inf")
            correct = ratio < 1.0

            results.append({
                "fold": i,
                "held_out": held_out,
                "error": error,
                "null_dist": null_dist,
                "ratio": ratio,
                "correct": correct,
                "n_pca_components": int(n_components),
            })

            status = "CORRECT" if correct else "*** FAIL ***"
            print(
                f"  {i + 1:>2}/35  {held_out:<50s} "
                f"error={error:>8.3f}  null={null_dist:>8.3f}  "
                f"ratio={ratio:.4f}  PCA={n_components}  {status}"
            )

    return results


# ======================================================================
# Main
# ======================================================================

def main():
    print("=" * 70)
    print("  V2 — LOOCV IMPLEMENTATION: INDEPENDENT VALIDATION")
    print(f"  Variance threshold: {VARIANCE_THRESHOLD}")
    print(f"  Mean ratio tolerance: ±{MEAN_RATIO_TOLERANCE}")
    print("=" * 70)

    # Load centroids
    print("\nLoading 35-type centroids...")
    h35 = pd.read_csv(HUMAN_35, index_col=0)
    m35 = pd.read_csv(MOUSE_35, index_col=0)
    print(f"  Human: {h35.shape}, Mouse: {m35.shape}")
    assert h35.shape == m35.shape == (35, 16959)

    # Run LOOCV
    print(f"\n{'=' * 70}")
    print("  LOOCV — 35 folds")
    print(f"{'=' * 70}\n")

    t0 = time.time()
    results = run_loocv(h35, m35)
    elapsed = time.time() - t0

    # ── Aggregate ─────────────────────────────────────────────────────
    n_correct = sum(1 for r in results if r["correct"])
    n_total = len(results)
    mean_ratio = float(np.mean([r["ratio"] for r in results]))
    median_ratio = float(np.median([r["ratio"] for r in results]))
    min_ratio = min(r["ratio"] for r in results)
    max_ratio = max(r["ratio"] for r in results)
    best_type = min(results, key=lambda r: r["ratio"])["held_out"]
    worst_type = max(results, key=lambda r: r["ratio"])["held_out"]

    pca_counts = [r["n_pca_components"] for r in results]
    pca_min = min(pca_counts)
    pca_max = max(pca_counts)
    pca_mode = max(set(pca_counts), key=pca_counts.count)

    print(f"\n{'=' * 70}")
    print("  LOOCV RESULTS")
    print(f"{'=' * 70}")
    print(f"  Correct recoveries: {n_correct}/{n_total}")
    print(f"  Mean ratio:   {mean_ratio:.4f}")
    print(f"  Median ratio: {median_ratio:.4f}")
    print(f"  Min ratio:    {min_ratio:.4f} ({best_type})")
    print(f"  Max ratio:    {max_ratio:.4f} ({worst_type})")
    print(f"  PCA components: {pca_min}–{pca_max} (mode={pca_mode})")
    print(f"  Runtime: {elapsed:.1f}s")

    # ── Pass/Fail ─────────────────────────────────────────────────────
    print(f"\n{'=' * 70}")
    print("  PASS/FAIL ASSESSMENT")
    print(f"{'=' * 70}")

    # Check 1: 35/35 correct
    recovery_pass = n_correct == 35
    if recovery_pass:
        print(f"  Correct recoveries: PASS ({n_correct}/35)")
    else:
        print(f"  Correct recoveries: *** FAIL *** ({n_correct}/35)")
        failing = [r for r in results if not r["correct"]]
        for r in failing:
            print(
                f"    FAILING FOLD: {r['held_out']}, "
                f"error={r['error']:.4f}, null={r['null_dist']:.4f}, "
                f"ratio={r['ratio']:.4f}"
            )

    # Check 2: mean ratio within tolerance
    expected_mean = 0.42
    ratio_diff = abs(mean_ratio - expected_mean)
    ratio_pass = ratio_diff <= MEAN_RATIO_TOLERANCE
    if ratio_pass:
        print(
            f"  Mean ratio: PASS "
            f"(expected {expected_mean}, got {mean_ratio:.4f}, "
            f"diff={ratio_diff:.4f})"
        )
    else:
        print(
            f"  Mean ratio: *** FAIL *** "
            f"(expected {expected_mean}, got {mean_ratio:.4f}, "
            f"diff={ratio_diff:.4f})"
        )

    # Overall verdict
    verdict = "PASS" if (recovery_pass and ratio_pass) else "FAIL"
    print(f"\n  VERDICT: {verdict}")

    # ── Save results ──────────────────────────────────────────────────
    output = {
        "validation": "V2 — LOOCV",
        "date": "2026-03-21",
        "verdict": verdict,
        "n_correct": n_correct,
        "n_total": n_total,
        "mean_ratio": mean_ratio,
        "median_ratio": median_ratio,
        "min_ratio": min_ratio,
        "max_ratio": max_ratio,
        "best_type": best_type,
        "worst_type": worst_type,
        "expected_mean_ratio": expected_mean,
        "ratio_diff": ratio_diff,
        "ratio_tolerance": MEAN_RATIO_TOLERANCE,
        "pca_component_range": [pca_min, pca_max],
        "pca_component_mode": pca_mode,
        "runtime_s": elapsed,
        "per_fold": results,
    }

    out_path = OUTPUT_DIR / "v2_loocv_results.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, default=lambda x: int(x) if isinstance(x, np.integer) else float(x) if isinstance(x, np.floating) else x)
    print(f"\n  Results saved to: {out_path}")


if __name__ == "__main__":
    main()
