#!/usr/bin/env python3
"""Re-run the K12 + D1-7 + human-mouse-12 permutation tests, saving the full
10K null arrays to .npy so Figure 5 panels can draw proper histograms.

Inputs:
  output/macaque_pipeline/reconstruction_qu12_centroids.csv
  output/macaque_pipeline/reconstruction_qu7_D1_centroids.csv
  output/phase2/scaled_35types/centroids_human_35.csv
  output/phase2/scaled_35types/centroids_mouse_35.csv

Outputs:
  output/macaque_pipeline/null_distribution_qu12.npy
  output/macaque_pipeline/null_distribution_qu7_D1.npy
  output/macaque_pipeline/null_distribution_hm12.npy
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT / "src"))
from procrustes import pca_reduce_centroids, procrustes_align, permutation_test  # noqa

SEED = 42
N_PERMS = 10_000

TYPES_12 = [
    "B cell", "bladder urothelial cell", "endothelial cell", "epithelial cell",
    "fibroblast", "hepatocyte", "macrophage", "monocyte", "neutrophil",
    "plasma cell", "smooth muscle cell", "stromal cell",
]
TYPES_7 = [
    "bladder urothelial cell", "endothelial cell", "epithelial cell",
    "fibroblast", "hepatocyte", "smooth muscle cell", "stromal cell",
]

OUT = PROJECT / "output/macaque_pipeline"
HUMAN = pd.read_csv(PROJECT / "output/phase2/scaled_35types/centroids_human_35.csv", index_col=0)
MOUSE = pd.read_csv(PROJECT / "output/phase2/scaled_35types/centroids_mouse_35.csv", index_col=0)
MAC12 = pd.read_csv(OUT / "reconstruction_qu12_centroids.csv", index_col=0)
MAC7 = pd.read_csv(OUT / "reconstruction_qu7_D1_centroids.csv", index_col=0)


def run(X_df, Y_df, label, out_npy):
    Xp, Yp, pca, order = pca_reduce_centroids(X_df, Y_df)
    r = procrustes_align(Xp, Yp)
    p, null = permutation_test(Xp, Yp, n_permutations=N_PERMS, seed=SEED)
    np.save(out_npy, null)
    return {
        "label": label,
        "distance": float(r.distance),
        "obs_null_median": float(r.distance / np.median(null)),
        "p_value": float(p),
        "null_median": float(np.median(null)),
        "null_path": str(out_npy),
        "n_types": len(order),
        "k_pca": int(pca.n_components_),
    }


def main():
    gene_ids = list(MAC12.columns)
    summaries = []
    # K12 macaque primary (human subset to 12 types, in 13,927 gene space)
    H12 = HUMAN.loc[TYPES_12, gene_ids].copy()
    summaries.append(run(H12, MAC12, "Qu12_primary", OUT / "null_distribution_qu12.npy"))

    # D1 7-type no-immune
    H7 = HUMAN.loc[TYPES_7, gene_ids].copy()
    summaries.append(run(H7, MAC7, "Qu7_D1_sensitivity", OUT / "null_distribution_qu7_D1.npy"))

    # Human-mouse 12-type control (uses full 16,959 human-mouse space — per prior run it's at 16,959)
    HUMAN_HM = HUMAN.loc[TYPES_12].copy()
    MOUSE_HM = MOUSE.loc[TYPES_12].copy()
    summaries.append(run(HUMAN_HM, MOUSE_HM, "hm12_control", OUT / "null_distribution_hm12.npy"))

    for s in summaries:
        print(s)


if __name__ == "__main__":
    main()
