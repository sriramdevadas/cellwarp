#!/usr/bin/env python3
"""Compute ||X_c||, ||Y_c||, trace(ΣD) for my reconstruction and compare to
what the committed sensitivity's (s=0.0534, d²=273.3) implies.

Tests:
  - raw 13,927-gene norms
  - after joint PCA at multiple thresholds
  - alternative normalizations for Y (target_sum 1e6, no normalize_total,
    log2(CPM+1), raw log-counts, inverse-transform from human_scaled convention)
"""
from __future__ import annotations
import sys
from pathlib import Path

import gzip
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io as sio
import scipy.sparse as sp

PROJECT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT / "src"))
from procrustes import pca_reduce_centroids, procrustes_align, permutation_test  # noqa

RANDOM_SEED = 42
TYPES_13 = [
    "B cell", "CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell",
    "T cell", "classical monocyte", "granulocyte", "hematopoietic precursor cell",
    "intermediate monocyte", "macrophage", "myeloid dendritic cell",
    "myeloid leukocyte", "natural killer cell", "non-classical monocyte",
]

GENE_LIST = PROJECT / "output/macaque_pipeline/reconstruction_rira13_gene_list.csv"
HUMAN_CSV = PROJECT / "output/phase2/scaled_35types/centroids_human_35.csv"
RIRA_CONV = PROJECT / "data/macaque/rira/converted"
RIRA_META = PROJECT / "data/macaque/rira/rira_metadata.csv"


def harmonize(meta):
    imm = meta["RIRA_Immune_v2.cellclass"]; tnk = meta["RIRA_TNK_v2.cellclass"]; mye = meta["RIRA_Myeloid_v3.cellclass"]
    out = pd.Series(pd.NA, index=meta.index, dtype="object")
    out[imm == "Bcell"] = "B cell"
    is_tnk = imm == "T_NK"
    out[is_tnk & (tnk == "CD4+ T Cells")] = "CD4-positive, alpha-beta T cell"
    out[is_tnk & (tnk == "CD8+ T Cells")] = "CD8-positive, alpha-beta T cell"
    out[is_tnk & (tnk == "NK Cells")] = "natural killer cell"
    out[is_tnk & ~tnk.isin(["CD4+ T Cells", "CD8+ T Cells", "NK Cells"])] = "T cell"
    is_mye = imm == "Myeloid"
    out[is_mye & (mye == "CD14+ Monocytes")] = "classical monocyte"
    out[is_mye & (mye == "Inflammatory Monocytes")] = "intermediate monocyte"
    out[is_mye & (mye == "CD16+ Monocytes")] = "non-classical monocyte"
    out[is_mye & mye.isin(["Macrophages", "Alv. mac."])] = "macrophage"
    out[is_mye & mye.isin(["DC", "pDC", "Mature DC"])] = "myeloid dendritic cell"
    out[is_mye & (mye == "Myelocytes")] = "granulocyte"
    out[is_mye & (mye == "Promyelocytes")] = "hematopoietic precursor cell"
    out[is_mye & ~mye.isin(["CD14+ Monocytes", "Inflammatory Monocytes", "CD16+ Monocytes",
                            "Macrophages", "Alv. mac.", "DC", "pDC", "Mature DC",
                            "Myelocytes", "Promyelocytes"])] = "myeloid leukocyte"
    return out


def compute_centroids(adata_sub, types):
    n_genes = adata_sub.n_vars
    out = np.zeros((len(types), n_genes))
    for i, t in enumerate(types):
        m = (adata_sub.obs["target_type"] == t).to_numpy()
        if m.sum() == 0:
            continue
        X = adata_sub.X[m]
        out[i] = np.asarray(X.mean(axis=0)).flatten()
    return out


def report(label, X, Y):
    Xc = X - X.mean(axis=0)
    Yc = Y - Y.mean(axis=0)
    nx, ny = (Xc**2).sum(), (Yc**2).sum()
    M = Xc.T @ Yc
    from numpy.linalg import svd
    U, sig, Vt = svd(M, full_matrices=False)
    V = Vt.T
    d = np.linalg.det(V @ U.T)
    D = np.ones(len(sig)); D[-1] = np.sign(d) if len(sig) else 1
    tr_sD = float((sig * D).sum())
    s = tr_sD / ny if ny else 0.0
    # Distance²
    d_sq = nx - tr_sD**2 / ny if ny else float("inf")
    print(f"  {label:<50}  ||X_c||²={nx:.2f}  ||Y_c||²={ny:.2f}  "
          f"trace(ΣD)={tr_sD:.4f}  s={s:.6f}  d²={d_sq:.2f}  d={np.sqrt(max(d_sq,0)):.3f}")
    return {"nx": nx, "ny": ny, "tr_sD": tr_sD, "s": s, "d_sq": d_sq}


def main():
    print("=" * 70)
    print("Diagnostic: committed sensitivity implies ||Y_c||² >> ||X_c||² (s≈0.053)")
    print("=" * 70)
    print(f"Committed target: d²=273.3, s=0.0534 → ||X_c||²≈d²+s²·||Y_c||²")
    print(f"If ||Y_c||²≈10000 then ||X_c||²≈273+28.5=301.5; if ||Y_c||²≈30000 then ||X_c||²≈273+85.6=358.9")

    # Load saved reconstruction
    mac_csv = PROJECT / "output/macaque_pipeline/reconstruction_rira13_centroids.csv"
    mac = pd.read_csv(mac_csv, index_col=0)
    human = pd.read_csv(HUMAN_CSV, index_col=0)
    gene_ids = list(mac.columns)
    X = human.loc[TYPES_13, gene_ids].values
    Y = mac.loc[TYPES_13, gene_ids].values

    print("\n--- RAW 13927-dim ---")
    report("X=human, Y=macaque (my reconstruction)", X, Y)

    print("\n--- after joint PCA ---")
    for nc in [0.95, 0.99, 0.90, 0.999]:
        from sklearn.decomposition import PCA
        stacked = np.vstack([X, Y])
        pca = PCA(n_components=nc, svd_solver="full", random_state=RANDOM_SEED)
        proj = pca.fit_transform(stacked)
        k = pca.n_components_
        Xp = proj[:13]; Yp = proj[13:]
        report(f"PCA cumvar={nc}, k={k}", Xp, Yp)

    # Now test alternative normalizations. Load matrix + subsample once.
    print("\n--- Alternative normalizations (Y reconstruction) ---")
    gene_df = pd.read_csv(GENE_LIST)
    three_sym = list(gene_df["symbol"])
    three_ensg = list(gene_df["ensg"])

    print("Loading RIRA matrix once…")
    with gzip.open(RIRA_CONV / "matrix.mtx.gz", "rt") as f:
        mtx = sio.mmread(f)
    mtx = sp.csr_matrix(mtx.T)
    with gzip.open(RIRA_CONV / "barcodes.tsv.gz", "rt") as f:
        barcodes = [line.strip() for line in f]
    with gzip.open(RIRA_CONV / "genes.tsv.gz", "rt") as f:
        genes = [line.strip() for line in f]

    meta = pd.read_csv(RIRA_META, usecols=["cell_id", "RIRA_Immune_v2.cellclass",
                                            "RIRA_TNK_v2.cellclass", "RIRA_Myeloid_v3.cellclass"],
                       low_memory=False).set_index("cell_id")
    meta = meta.reindex(barcodes)
    target_type = harmonize(meta)
    kept = target_type.notna().to_numpy()

    # Subsample 2000/type
    rng = np.random.default_rng(RANDOM_SEED)
    sub = np.zeros(len(barcodes), dtype=bool)
    for t in TYPES_13:
        idx = np.where((target_type == t).to_numpy())[0]
        sel = rng.choice(idx, size=min(2000, len(idx)), replace=False) if len(idx) > 2000 else idx
        sub[sel] = True
    keep = kept & sub
    # Gene filter to three-way
    gene_to_idx = {g: i for i, g in enumerate(genes)}
    sym_idx = [gene_to_idx[s] for s in three_sym if s in gene_to_idx]
    raw = mtx[keep][:, sym_idx]
    tt = target_type[keep].to_numpy()

    def build_adata(X_sparse, norm):
        a = ad.AnnData(X=X_sparse, obs=pd.DataFrame({"target_type": tt}, index=np.arange(X_sparse.shape[0]).astype(str)))
        a.var_names = three_sym  # symbols; alignment works since we preserve order
        if norm == "none":
            pass
        elif norm == "log1p_only":
            sc.pp.log1p(a)
        elif norm == "nt1e4_log1p":
            sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
        elif norm == "nt1e6_log1p":
            sc.pp.normalize_total(a, target_sum=1e6); sc.pp.log1p(a)
        elif norm == "log2cpm":
            sc.pp.normalize_total(a, target_sum=1e6); a.X = a.X.log1p()  # natural log OK for comparison
        return a

    for norm in ["nt1e4_log1p", "nt1e6_log1p", "log1p_only"]:
        a = build_adata(raw.copy(), norm)
        Ynew = compute_centroids(a, TYPES_13)
        print(f"\n  Normalization = {norm}")
        print(f"    Y range: min={Ynew.min():.3f} max={Ynew.max():.3f} mean={Ynew.mean():.3f} std={Ynew.std():.3f}")
        info = report(f"    X=human, Y=macaque({norm})", X, Ynew)

        # Also run PCA + Procrustes test
        from sklearn.decomposition import PCA
        stacked = np.vstack([X, Ynew])
        pca = PCA(n_components=0.95, svd_solver="full", random_state=RANDOM_SEED)
        proj = pca.fit_transform(stacked)
        Xp, Yp = proj[:13], proj[13:]
        r = procrustes_align(Xp, Yp)
        p, null = permutation_test(Xp, Yp, n_permutations=1000, seed=RANDOM_SEED)
        obs_null = r.distance / np.median(null)
        print(f"    → PCA(k={pca.n_components_}) Procrustes: d={r.distance:.3f}  s={r.scaling:.4f}  "
              f"obs/null(med)={obs_null:.4f}  [target 0.749]")


if __name__ == "__main__":
    main()
