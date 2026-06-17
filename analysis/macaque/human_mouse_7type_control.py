#!/usr/bin/env python3
"""Human-mouse 7-type control for the D1 pilot.

Purpose: apples-to-apples comparison of n=7 Procrustes statistics between
species pairs. If human-mouse 7-type yields similar obs/null and similar
p ≈ 0.01, the macaque p is a property of n=7 permutation-space ceiling,
not a weaker macaque comparison.

Uses committed centroids from the symmetric primary pipeline
(normalize_total(1e4) + log1p on raw counts from Tabula Sapiens and
Tabula Muris Senis). Same 7 cell types as the D1 Qu pilot.
"""
from __future__ import annotations
import json, sys, time
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT / "src"))
from procrustes import pca_reduce_centroids, procrustes_align, permutation_test  # noqa

RANDOM_SEED = 42
N_PERMUTATIONS = 10_000

HUMAN_CSV = PROJECT / "output/phase2/scaled_35types/centroids_human_35.csv"
MOUSE_CSV = PROJECT / "output/phase2/scaled_35types/centroids_mouse_35.csv"

OUT_JSON = PROJECT / "output/macaque_pipeline/human_mouse_7type_control.json"

TYPES_7 = [
    "bladder urothelial cell",
    "endothelial cell",
    "epithelial cell",
    "fibroblast",
    "hepatocyte",
    "smooth muscle cell",
    "stromal cell",
]


def main():
    t0 = time.time()
    print("=" * 70)
    print("Human-mouse 7-type control (matches D1 target types)")
    print("=" * 70)

    human = pd.read_csv(HUMAN_CSV, index_col=0)
    mouse = pd.read_csv(MOUSE_CSV, index_col=0)
    print(f"Human full: {human.shape}, Mouse full: {mouse.shape}")

    # Verify all 7 types exist in both
    missing_h = [t for t in TYPES_7 if t not in human.index]
    missing_m = [t for t in TYPES_7 if t not in mouse.index]
    if missing_h or missing_m:
        raise SystemExit(f"Missing types — human: {missing_h}, mouse: {missing_m}")

    # Subset (use 16,959-gene human-mouse space; both CSVs share the same columns)
    assert list(human.columns) == list(mouse.columns), "Human/mouse gene columns differ"
    gene_space = len(human.columns)
    human_7 = human.loc[TYPES_7].copy()
    mouse_7 = mouse.loc[TYPES_7].copy()
    print(f"Subset: human_7 {human_7.shape}, mouse_7 {mouse_7.shape}")

    # Sanity
    hv, mv = human_7.to_numpy(), mouse_7.to_numpy()
    print(f"Human:  min={hv.min():.3f} max={hv.max():.3f} std={hv.std():.3f}")
    print(f"Mouse:  min={mv.min():.3f} max={mv.max():.3f} std={mv.std():.3f}")
    print(f"std ratio (mouse/human): {mv.std() / max(hv.std(), 1e-12):.3f}")

    # PCA + Procrustes
    Xp, Yp, pca, order = pca_reduce_centroids(human_7, mouse_7)
    result = procrustes_align(Xp, Yp)
    p, null = permutation_test(Xp, Yp, n_permutations=N_PERMUTATIONS, seed=RANDOM_SEED)
    obs_null_med = result.distance / np.median(null)
    obs_null_mean = result.distance / np.mean(null)

    # Per-type residuals
    residuals = []
    for i, typ in enumerate(order):
        r = result.aligned_target[i] - result.centered_reference[i]
        mag = float(np.linalg.norm(r))
        pct = 100.0 * mag**2 / result.distance_squared if result.distance_squared > 0 else 0.0
        residuals.append({"type": typ, "magnitude": mag, "pct_ssr": pct})
    residuals_sorted = sorted(residuals, key=lambda r: -r["magnitude"])

    out = {
        "analysis": "HUMAN_MOUSE_7TYPE_CONTROL",
        "cell_types": order,
        "gene_space": gene_space,
        "procrustes": {
            "distance": float(result.distance),
            "distance_squared": float(result.distance_squared),
            "scaling": float(result.scaling),
            "rotation_det": float(np.linalg.det(result.rotation)),
        },
        "permutation_test": {
            "p_value": float(p),
            "obs_null_ratio_median": float(obs_null_med),
            "obs_null_ratio_mean": float(obs_null_mean),
            "null_median": float(np.median(null)),
            "null_mean": float(np.mean(null)),
            "n_permutations": N_PERMUTATIONS,
            "permuted_le_observed": int(np.sum(null <= result.distance)),
        },
        "pca": {
            "n_components": int(pca.n_components_),
            "variance_explained": pca.explained_variance_ratio_.tolist(),
            "cumulative_variance": float(np.cumsum(pca.explained_variance_ratio_)[-1]),
        },
        "per_type_residuals_ranked": residuals_sorted,
        "seed": RANDOM_SEED,
        "runtime_s": time.time() - t0,
        "preprocessing": "committed symmetric: normalize_total(1e4) + log1p (raw Tabula)",
    }

    def _default(o):
        if isinstance(o, np.bool_): return bool(o)
        if isinstance(o, np.integer): return int(o)
        if isinstance(o, np.floating): return float(o)
        if isinstance(o, np.ndarray): return o.tolist()
        raise TypeError(f"not serialisable: {type(o).__name__}")
    OUT_JSON.write_text(json.dumps(out, indent=2, default=_default))

    print("\n" + "=" * 70)
    print("HUMAN-MOUSE 7-TYPE CONTROL RESULTS")
    print("=" * 70)
    print(f"  gene_space:      {gene_space}")
    print(f"  distance:        {result.distance:.4f}")
    print(f"  scaling:         {result.scaling:.6f}")
    print(f"  obs/null (med):  {obs_null_med:.4f}")
    print(f"  obs/null (mean): {obs_null_mean:.4f}")
    print(f"  p-value:         {p:.6f}")
    print(f"  permuted ≤ obs:  {int(np.sum(null <= result.distance))} / {N_PERMUTATIONS}")
    print(f"  PCA k:           {pca.n_components_}")
    print()
    print("  Per-type residuals (high → low):")
    for d in residuals_sorted:
        print(f"    {d['type']:<30s}  mag={d['magnitude']:.3f}  pct_ssr={d['pct_ssr']:.1f}%")
    print(f"\n  Results: {OUT_JSON}")


if __name__ == "__main__":
    main()
