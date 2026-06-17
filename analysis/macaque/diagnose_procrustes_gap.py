#!/usr/bin/env python3
"""Diagnose why my 13-type Procrustes obs/null (0.497) differs from committed
sensitivity (0.749). Centroid cell counts + 13,927-gene list reproduced exactly.

Hypotheses to test:
  H1: The committed sensitivity inherited the 20-type PCA basis (8 components)
      rather than recomputing fresh PCA on the 13-type stack (11 components).
  H2: The committed pipeline applied a different normalization to RIRA (e.g.,
      log1p without normalize_total, or a different target_sum).
  H3: The committed pipeline used `procrustes_align(X, Y)` with X=macaque,
      Y=human instead of X=human, Y=macaque.

Reuses saved reconstructed centroids; re-runs Procrustes under each scenario.
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

PROJECT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT / "src"))
from procrustes import procrustes_align, permutation_test  # noqa

RANDOM_SEED = 42
N_PERMS = 10_000

HUMAN_CSV = PROJECT / "output/phase2/scaled_35types/centroids_human_35.csv"
MAC_CSV = PROJECT / "output/macaque_pipeline/reconstruction_rira13_centroids.csv"
EXPECTED_OBS_NULL = 0.7488054428
EXPECTED_DIST = 16.5307225745
EXPECTED_SCALING = 0.0534434720

TYPES_13 = [
    "B cell",
    "CD4-positive, alpha-beta T cell",
    "CD8-positive, alpha-beta T cell",
    "T cell",
    "classical monocyte",
    "granulocyte",
    "hematopoietic precursor cell",
    "intermediate monocyte",
    "macrophage",
    "myeloid dendritic cell",
    "myeloid leukocyte",
    "natural killer cell",
    "non-classical monocyte",
]


def joint_pca(X_df, Y_df, n_components):
    """Stack + PCA at the given n_components (int or float variance threshold)."""
    types = sorted(X_df.index.tolist())
    assert sorted(Y_df.index.tolist()) == types
    X = X_df.loc[types].values
    Y = Y_df.loc[types].values
    stacked = np.vstack([X, Y])
    pca = PCA(n_components=n_components, svd_solver="full", random_state=RANDOM_SEED)
    proj = pca.fit_transform(stacked)
    n = len(types)
    return proj[:n], proj[n:], pca.n_components_, pca.explained_variance_ratio_.sum()


def run_procrustes(X_pca, Y_pca, label, swap=False):
    """Run Procrustes + permutation and print concise summary."""
    if swap:
        X_pca, Y_pca = Y_pca, X_pca
    r = procrustes_align(X_pca, Y_pca)
    p, null = permutation_test(X_pca, Y_pca, n_permutations=N_PERMS, seed=RANDOM_SEED)
    obs_null_med = r.distance / np.median(null)
    obs_null_mean = r.distance / np.mean(null)
    delta = abs(obs_null_med - EXPECTED_OBS_NULL)
    tag = "PASS" if delta < 0.01 else "FAIL"
    print(f"\n  >> {label}")
    print(f"     distance={r.distance:.4f}  scaling={r.scaling:.6f}  "
          f"obs/null(med)={obs_null_med:.4f}  obs/null(mean)={obs_null_mean:.4f}  "
          f"p={p:.6f}  [{tag}, Δ={delta:.4f} from 0.749]")
    return r, obs_null_med, p


def main():
    print("Loading centroids…")
    human = pd.read_csv(HUMAN_CSV, index_col=0)
    mac = pd.read_csv(MAC_CSV, index_col=0)

    gene_ids = list(mac.columns)
    human_13 = human.loc[TYPES_13, gene_ids].copy()
    mac_13 = mac.loc[TYPES_13, gene_ids].copy()

    print(f"  human_13: {human_13.shape}  mac_13: {mac_13.shape}")
    print(f"  Expected: obs/null = {EXPECTED_OBS_NULL:.4f}, distance = {EXPECTED_DIST:.4f}, scaling = {EXPECTED_SCALING:.6f}")

    print("\n=== H1: PCA dimensionality ===")
    for nc in [8, 11, 0.95, 0.90, 0.99, 12, 7, 6]:
        Xp, Yp, k, cv = joint_pca(human_13, mac_13, nc)
        run_procrustes(Xp, Yp, f"n_components={nc}  → k={k}, cumvar={cv:.3f}")

    print("\n=== H3: swap X/Y (human↔macaque as ref/target) ===")
    Xp, Yp, k, cv = joint_pca(human_13, mac_13, 0.95)
    run_procrustes(Xp, Yp, f"swap (Y=human, X=macaque), k={k}", swap=True)


if __name__ == "__main__":
    main()
