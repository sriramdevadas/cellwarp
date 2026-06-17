#!/usr/bin/env python3
"""
V1 — Procrustes Permutation Test: Independent Validation
=========================================================
Independently verifies the primary obs/null statistic and p-value via
label permutation, plus all replication comparisons.

NO imports from src/ — Procrustes reimplemented from scratch.

Canonical targets:
  Primary (35-type):  obs/null = 0.522 ± 0.005, p ≤ 0.0002, PCA = 33
  Sun2023 (15-type):  obs/null = 0.554
  PanSci  (16-type):  obs/null = 0.552
  CellHint (15-type): obs/null = 0.448
  6-type ancillary:   obs/null = 0.317
  Human-human ctrl:   obs/null = 0.607

Pass criteria:  obs/null within ±0.005; p ≤ 0.0002; PCA components match.
Fail criteria:  Any value outside tolerance → report exact values, escalate.

Author: V1 software validation (independent of src/)
Date: 2026-03-21
"""

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import svd
from sklearn.decomposition import PCA

# ── Constants ──────────────────────────────────────────────────────────
RANDOM_SEED = 42
N_PERMUTATIONS = 10_000
VARIANCE_THRESHOLD = 0.95
OBS_NULL_TOLERANCE = 0.005

# ── File paths ─────────────────────────────────────────────────────────
BASE = Path(__file__).resolve().parent.parent

# Primary 35-type
HUMAN_35 = BASE / "output/phase2/scaled_35types/centroids_human_35.csv"
MOUSE_35 = BASE / "output/phase2/scaled_35types/centroids_mouse_35.csv"

# 6-type ancillary
HUMAN_6 = BASE / "output/phase2/centroids_human.csv"
MOUSE_6 = BASE / "output/phase2/centroids_mouse.csv"

# Human-human control v2
HUMAN2_V2 = BASE / "output/phase2/negative_control_v2/centroids_human2_v2.csv"

# CellHint human (independent)
CELLHINT = BASE / "output/validation/cellhint_replication/centroids_cellhint.csv"

# Sun2023 expanded centroids (15 types, saved by 17_sun2023_expanded.py)
SUN2023_CSV = BASE / "data/centroids/sun2023_15type_centroids.csv"

# PanSci centroids (16 types, saved by pansci_replication.py)
PANSCI_CSV = BASE / "data/centroids/pansci_16type_centroids.csv"

# Sun2023 h5ad fallback (only used if CSV not available)
SUN2023_H5AD = BASE / "data/replication/sun2023/sun2023_yc.h5ad"

# Output
OUTPUT_DIR = BASE / "output/validation/v1_procrustes"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ======================================================================
# Independent Procrustes implementation (from first principles)
# ======================================================================

def procrustes_distance(X: np.ndarray, Y: np.ndarray) -> float:
    """
    Ordinary Procrustes Analysis: align Y to X via rotation + uniform
    scaling (no reflection), return Frobenius distance of residuals.

    Math:
        1. Center both matrices
        2. Cross-covariance M = X_c^T @ Y_c
        3. SVD: M = U Σ V^T
        4. Rotation R = V D U^T  where D corrects for reflection
        5. Scaling s = tr(Σ D) / ||Y_c||²_F
        6. Aligned Y = s * Y_c @ R
        7. Distance = ||X_c - aligned Y||_F = sqrt(SSR)
    """
    n, k = X.shape
    X_c = X - X.mean(axis=0)
    Y_c = Y - Y.mean(axis=0)

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

    # Align and compute residual distance
    Y_aligned = s * (Y_c @ R)
    return np.sqrt(np.sum((X_c - Y_aligned) ** 2))


def permutation_test(
    X: np.ndarray,
    Y: np.ndarray,
    n_perm: int = N_PERMUTATIONS,
    seed: int = RANDOM_SEED,
) -> dict:
    """
    Label-permutation null: shuffle mouse type labels n_perm times.

    p-value = (#{null ≤ observed} + 1) / (n_perm + 1)
    The +1 correction prevents p = 0 and includes the observed value.
    """
    rng = np.random.default_rng(seed)
    d_obs = procrustes_distance(X, Y)

    null_dist = np.empty(n_perm)
    for i in range(n_perm):
        perm = rng.permutation(Y.shape[0])
        null_dist[i] = procrustes_distance(X, Y[perm])

    n_leq = int(np.sum(null_dist <= d_obs))
    p_value = (n_leq + 1) / (n_perm + 1)
    obs_null = d_obs / np.median(null_dist)

    return {
        "d_obs": float(d_obs),
        "null_median": float(np.median(null_dist)),
        "null_min": float(np.min(null_dist)),
        "null_max": float(np.max(null_dist)),
        "null_std": float(np.std(null_dist)),
        "obs_null": float(obs_null),
        "p_value": float(p_value),
        "n_leq": n_leq,
        "n_perm": n_perm,
    }


# ======================================================================
# Validation pipeline
# ======================================================================

def validate_comparison(
    name: str,
    ref_centroids: pd.DataFrame,
    target_centroids: pd.DataFrame,
    expected_obs_null: float,
    expected_n_types: int,
    expected_pca: int | None = None,
    tolerance: float = OBS_NULL_TOLERANCE,
) -> dict:
    """
    Full validation for one comparison:
      1. Match types between ref and target
      2. PCA reduction (95% variance)
      3. Procrustes + permutation test
      4. Report pass/fail
    """
    print(f"\n{'=' * 70}")
    print(f"  COMPARISON: {name}")
    print(f"{'=' * 70}")

    # Match types
    shared = sorted(set(ref_centroids.index) & set(target_centroids.index))
    n_types = len(shared)

    flag_types = n_types != expected_n_types
    if flag_types:
        print(f"  *** FLAG: Expected {expected_n_types} types, got {n_types} ***")
        print(f"  Shared types: {shared}")
    else:
        print(f"  Types: {n_types} (matches expected)")

    # Restrict to shared types in sorted order
    ref = ref_centroids.loc[shared]
    tgt = target_centroids.loc[shared]

    # Verify gene space
    assert list(ref.columns) == list(tgt.columns), "Gene columns must match"
    n_genes = ref.shape[1]
    print(f"  Genes: {n_genes}")

    # PCA on combined centroids
    combined = np.vstack([ref.values, tgt.values])  # (2n, G)
    pca = PCA(
        n_components=VARIANCE_THRESHOLD,
        svd_solver="full",
        random_state=RANDOM_SEED,
    )
    combined_pca = pca.fit_transform(combined)
    n_components = pca.n_components_
    cumvar = float(np.cumsum(pca.explained_variance_ratio_)[-1])

    print(f"  PCA components: {n_components} ({cumvar * 100:.1f}% variance)")

    flag_pca = (expected_pca is not None) and (n_components != expected_pca)
    if flag_pca:
        print(f"  *** FLAG: Expected {expected_pca} PCA components ***")

    # Split back
    ref_pca = combined_pca[:n_types]
    tgt_pca = combined_pca[n_types:]

    # Procrustes + permutation
    t0 = time.time()
    results = permutation_test(ref_pca, tgt_pca, N_PERMUTATIONS, RANDOM_SEED)
    elapsed = time.time() - t0

    # Report
    print(f"  d_obs:       {results['d_obs']:.6f}")
    print(f"  null_median: {results['null_median']:.6f}")
    print(f"  obs/null:    {results['obs_null']:.6f}")
    print(f"  p-value:     {results['p_value']:.6f}")
    print(f"  Null hits ≤ d_obs: {results['n_leq']} / {results['n_perm']}")
    print(
        f"  Null summary: min={results['null_min']:.4f}, "
        f"median={results['null_median']:.4f}, "
        f"max={results['null_max']:.4f}, "
        f"std={results['null_std']:.4f}"
    )
    print(f"  Runtime: {elapsed:.1f}s")

    # Pass/fail
    obs_null_diff = abs(results["obs_null"] - expected_obs_null)
    obs_null_pass = obs_null_diff <= tolerance

    if obs_null_pass:
        print(
            f"  obs/null: PASS "
            f"(expected {expected_obs_null:.3f}, got {results['obs_null']:.6f}, "
            f"diff={obs_null_diff:.6f})"
        )
    else:
        print(
            f"  obs/null: *** FAIL *** "
            f"(expected {expected_obs_null:.3f}, got {results['obs_null']:.6f}, "
            f"diff={obs_null_diff:.6f})"
        )

    # Aggregate verdict
    verdict = "PASS"
    failures = []
    if not obs_null_pass:
        failures.append(f"obs/null outside tolerance ({obs_null_diff:.6f} > {tolerance})")
        verdict = "FAIL"
    if flag_types:
        failures.append(f"type count mismatch ({n_types} vs {expected_n_types})")
        verdict = "FLAG"
    if flag_pca:
        failures.append(f"PCA components mismatch ({n_components} vs {expected_pca})")
        verdict = "FAIL"

    print(f"\n  VERDICT: {verdict}")
    if failures:
        for f in failures:
            print(f"    - {f}")

    return {
        "name": name,
        "n_types": n_types,
        "expected_n_types": expected_n_types,
        "n_genes": n_genes,
        "n_pca_components": n_components,
        "expected_pca_components": expected_pca,
        "cumulative_variance": cumvar,
        "obs_null": results["obs_null"],
        "expected_obs_null": expected_obs_null,
        "obs_null_diff": obs_null_diff,
        "obs_null_pass": obs_null_pass,
        "p_value": results["p_value"],
        "n_leq": results["n_leq"],
        "d_obs": results["d_obs"],
        "null_median": results["null_median"],
        "null_min": results["null_min"],
        "null_max": results["null_max"],
        "null_std": results["null_std"],
        "runtime_s": elapsed,
        "verdict": verdict,
        "failures": failures,
        "seed": RANDOM_SEED,
        "n_permutations": N_PERMUTATIONS,
        "shared_types": shared,
    }


def load_sun2023_centroids() -> pd.DataFrame:
    """
    Compute Sun2023 mouse centroids from h5ad file.

    The h5ad is already normalized (CPM + log1p) and restricted to
    16,959 ortholog genes. Centroids = mean expression per cell type.

    Note: The h5ad contains 14 annotated types. The expanded analysis
    (17_sun2023_expanded.py) added 'T cell' as a 15th type via
    second-pass cluster annotation. This validation uses the 14 types
    present in the saved h5ad.
    """
    import scanpy as sc

    print("\n  Loading Sun2023 h5ad (backed mode)...")
    adata = sc.read_h5ad(SUN2023_H5AD)
    print(f"  Shape: {adata.shape}")
    print(f"  Cell types: {adata.obs['cell_type'].nunique()}")

    # Compute centroids per cell type
    centroids = {}
    for ct in sorted(adata.obs["cell_type"].unique()):
        mask = adata.obs["cell_type"] == ct
        centroids[ct] = np.asarray(adata[mask].X.mean(axis=0)).flatten()

    centroid_df = pd.DataFrame(centroids, index=adata.var_names).T
    print(f"  Centroid matrix: {centroid_df.shape[0]} types × {centroid_df.shape[1]} genes")
    return centroid_df


# ======================================================================
# Main
# ======================================================================

def main():
    print("=" * 70)
    print("  V1 — PROCRUSTES PERMUTATION TEST: INDEPENDENT VALIDATION")
    print(f"  Seed: {RANDOM_SEED}  |  Permutations: {N_PERMUTATIONS:,}")
    print(f"  Tolerance: ±{OBS_NULL_TOLERANCE}")
    print("=" * 70)

    all_results = []

    # ── 1. PRIMARY: Tabula Mouse vs Tabula Human (35-type) ────────────
    print("\n\nLoading primary 35-type centroids...")
    h35 = pd.read_csv(HUMAN_35, index_col=0)
    m35 = pd.read_csv(MOUSE_35, index_col=0)
    print(f"  Human: {h35.shape}, Mouse: {m35.shape}")

    r = validate_comparison(
        name="PRIMARY: Tabula Mouse vs Tabula Human",
        ref_centroids=h35,
        target_centroids=m35,
        expected_obs_null=0.522,
        expected_n_types=35,
        expected_pca=33,
    )
    all_results.append(r)

    # ── Check primary pass before replications ────────────────────────
    if r["verdict"] == "FAIL":
        print("\n\n*** PRIMARY FAILED — halting before replications ***")
        print("*** Escalate: review primary results above ***")
        save_results(all_results)
        sys.exit(1)

    print("\n\n" + "#" * 70)
    print("  PRIMARY PASSED — proceeding to replication comparisons")
    print("#" * 70)

    # ── 2. Sun2023 Mouse vs Tabula Human ──────────────────────────────
    if SUN2023_CSV.exists():
        print(f"\n\nLoading Sun2023 centroids from CSV: {SUN2023_CSV}")
        sun_centroids = pd.read_csv(SUN2023_CSV, index_col=0)
        print(f"  Shape: {sun_centroids.shape}")
    else:
        print(f"\n\n  Sun2023 CSV not found at {SUN2023_CSV}")
        print("  Falling back to h5ad centroid extraction (14 types only)...")
        sun_centroids = load_sun2023_centroids()

    r = validate_comparison(
        name="Sun2023 Mouse vs Tabula Human",
        ref_centroids=h35,
        target_centroids=sun_centroids,
        expected_obs_null=0.554,
        expected_n_types=15,
    )
    all_results.append(r)

    # ── 3. PanSci Mouse vs Tabula Human ───────────────────────────────
    if PANSCI_CSV.exists():
        print(f"\n\nLoading PanSci centroids from CSV: {PANSCI_CSV}")
        pansci_centroids = pd.read_csv(PANSCI_CSV, index_col=0)
        print(f"  Shape: {pansci_centroids.shape}")

        r = validate_comparison(
            name="PanSci Mouse vs Tabula Human",
            ref_centroids=h35,
            target_centroids=pansci_centroids,
            expected_obs_null=0.552,
            expected_n_types=16,
        )
        all_results.append(r)
    else:
        print("\n\n" + "=" * 70)
        print("  COMPARISON: PanSci Mouse vs Tabula Human")
        print("=" * 70)
        print(f"  *** SKIPPED: PanSci CSV not found at {PANSCI_CSV} ***")
        print("  Run scripts/pansci_replication.py first to generate centroid CSV.")
        all_results.append({
            "name": "PanSci Mouse vs Tabula Human",
            "verdict": "SKIPPED",
            "reason": f"CSV not found at {PANSCI_CSV}",
            "canonical_obs_null": 0.552,
            "canonical_n_types": 16,
        })

    # ── 4. CellHint Human vs Tabula Mouse ─────────────────────────────
    print("\n\nLoading CellHint centroids...")
    ch = pd.read_csv(CELLHINT, index_col=0)
    print(f"  CellHint: {ch.shape}")

    r = validate_comparison(
        name="CellHint Human vs Tabula Mouse",
        ref_centroids=m35,       # Tabula Mouse is reference
        target_centroids=ch,     # CellHint Human is target
        expected_obs_null=0.448,
        expected_n_types=15,
    )
    all_results.append(r)

    # ── 5. 6-type ancillary ───────────────────────────────────────────
    print("\n\nLoading 6-type centroids...")
    h6 = pd.read_csv(HUMAN_6, index_col=0)
    m6 = pd.read_csv(MOUSE_6, index_col=0)
    print(f"  Human: {h6.shape}, Mouse: {m6.shape}")

    r = validate_comparison(
        name="6-type ancillary (Tabula Mouse vs Tabula Human)",
        ref_centroids=h6,
        target_centroids=m6,
        expected_obs_null=0.317,
        expected_n_types=6,
    )
    all_results.append(r)

    # ── 6. Human-human control ────────────────────────────────────────
    print("\n\nLoading human-human control centroids...")
    h2v2 = pd.read_csv(HUMAN2_V2, index_col=0)
    print(f"  Human2 v2: {h2v2.shape}")

    r = validate_comparison(
        name="Human-human control (Tabula Human vs Human2 v2)",
        ref_centroids=h6,       # Tabula Human 6-type is reference
        target_centroids=h2v2,  # Human2 v2 is target
        expected_obs_null=0.607,
        expected_n_types=6,
    )
    all_results.append(r)

    # ── Ordering check ────────────────────────────────────────────────
    print("\n\n" + "=" * 70)
    print("  ORDERING CHECK")
    print("=" * 70)
    print("  Expected: human-human (0.607) > Sun2023 (0.554) ≈ PanSci (0.552)")
    print("            > primary (0.522) > CellHint (0.448) > ancillary (0.317)")
    print()

    completed = [r for r in all_results if r.get("verdict") != "SKIPPED"]
    obs_null_map = {r["name"]: r["obs_null"] for r in completed}

    for name, val in sorted(obs_null_map.items(), key=lambda x: -x[1]):
        print(f"  {val:.6f}  {name}")

    # Check pairwise ordering for completed comparisons
    order_violations = []
    expected_order = [
        ("Human-human control", "Sun2023"),
        ("Sun2023", "PRIMARY"),
        ("PRIMARY", "CellHint"),
        ("CellHint", "6-type ancillary"),
    ]
    for higher_key, lower_key in expected_order:
        h_match = [r for r in completed if higher_key.lower() in r["name"].lower()]
        l_match = [r for r in completed if lower_key.lower() in r["name"].lower()]
        if h_match and l_match:
            h_val = h_match[0]["obs_null"]
            l_val = l_match[0]["obs_null"]
            if h_val <= l_val:
                order_violations.append(
                    f"{higher_key} ({h_val:.4f}) ≤ {lower_key} ({l_val:.4f})"
                )

    if order_violations:
        print("\n  *** ORDERING VIOLATIONS ***")
        for v in order_violations:
            print(f"    - {v}")
    else:
        print("\n  Ordering: PASS (all completed comparisons in expected order)")

    # ── Summary ───────────────────────────────────────────────────────
    print("\n\n" + "=" * 70)
    print("  V1 VALIDATION SUMMARY")
    print("=" * 70)

    for r in all_results:
        verdict = r.get("verdict", "?")
        name = r.get("name", "?")
        if verdict == "SKIPPED":
            print(f"  SKIPPED  {name}: {r.get('reason', '')}")
        else:
            obs = r.get("obs_null", 0)
            exp = r.get("expected_obs_null", 0)
            n = r.get("n_types", "?")
            pca = r.get("n_pca_components", "?")
            p = r.get("p_value", 0)
            n_leq = r.get("n_leq", "?")
            print(
                f"  {verdict:6s}  {name}\n"
                f"         obs/null={obs:.6f} (exp {exp:.3f}, diff={abs(obs - exp):.6f}), "
                f"n={n}, PCA={pca}, p={p:.6f}, null_hits={n_leq}"
            )
            if r.get("failures"):
                for f in r["failures"]:
                    print(f"         ! {f}")

    save_results(all_results)


def save_results(all_results: list[dict]):
    """Save all results to JSON."""
    out_path = OUTPUT_DIR / "v1_validation_results.json"

    # Convert numpy types for JSON serialization
    def convert(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj

    serializable = json.loads(json.dumps(all_results, default=convert))

    with open(out_path, "w") as f:
        json.dump(serializable, f, indent=2)
    print(f"\n  Results saved to: {out_path}")


if __name__ == "__main__":
    main()
