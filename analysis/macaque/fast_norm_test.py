#!/usr/bin/env python3
"""Fast hypothesis test: what normalization of RIRA .X reproduces committed s≈0.053?
RIRA matrix data is NOT raw counts (non-integer, max 16937, per-cell max ~100).
Options tested: as-is, log1p only, normalize_total(1e4)+log1p.
"""
from __future__ import annotations
import gzip, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
import anndata as ad
import scanpy as sc

PROJECT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT / "src"))
from procrustes import procrustes_align, permutation_test  # noqa
from sklearn.decomposition import PCA

RANDOM_SEED = 42
TYPES_13 = [
    "B cell", "CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell",
    "T cell", "classical monocyte", "granulocyte", "hematopoietic precursor cell",
    "intermediate monocyte", "macrophage", "myeloid dendritic cell",
    "myeloid leukocyte", "natural killer cell", "non-classical monocyte",
]
RIRA_CONV = PROJECT / "data/macaque/rira/converted"
RIRA_META = PROJECT / "data/macaque/rira/rira_metadata.csv"
GENE_LIST = PROJECT / "output/macaque_pipeline/reconstruction_rira13_gene_list.csv"
HUMAN_CSV = PROJECT / "output/phase2/scaled_35types/centroids_human_35.csv"


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


def main():
    t0 = time.time()
    print(f"[{time.time()-t0:.0f}s] Loading gene list + human centroids…", flush=True)
    gene_df = pd.read_csv(GENE_LIST)
    three_sym = list(gene_df["symbol"]); three_ensg = list(gene_df["ensg"])
    human = pd.read_csv(HUMAN_CSV, index_col=0)
    X = human.loc[TYPES_13, three_ensg].values  # (13, 13927)

    print(f"[{time.time()-t0:.0f}s] Loading RIRA matrix…", flush=True)
    with gzip.open(RIRA_CONV / "matrix.mtx.gz", "rt") as f:
        mtx = sio.mmread(f)
    mtx = sp.csr_matrix(mtx.T)
    with gzip.open(RIRA_CONV / "barcodes.tsv.gz", "rt") as f:
        barcodes = [line.strip() for line in f]
    with gzip.open(RIRA_CONV / "genes.tsv.gz", "rt") as f:
        genes = [line.strip() for line in f]
    print(f"[{time.time()-t0:.0f}s] Matrix loaded, shape={mtx.shape}", flush=True)

    meta = pd.read_csv(RIRA_META, usecols=["cell_id", "RIRA_Immune_v2.cellclass",
                                            "RIRA_TNK_v2.cellclass", "RIRA_Myeloid_v3.cellclass"],
                       low_memory=False).set_index("cell_id").reindex(barcodes)
    tgt = harmonize(meta)
    kept = tgt.notna().to_numpy()

    rng = np.random.default_rng(RANDOM_SEED)
    sub = np.zeros(len(barcodes), dtype=bool)
    for t in TYPES_13:
        idx = np.where((tgt == t).to_numpy())[0]
        sel = rng.choice(idx, size=min(2000, len(idx)), replace=False) if len(idx) > 2000 else idx
        sub[sel] = True
    keep = kept & sub
    print(f"[{time.time()-t0:.0f}s] Kept {int(keep.sum()):,} cells after 2000/type subsample", flush=True)

    # Gene filter
    gmap = {g: i for i, g in enumerate(genes)}
    sym_idx = [gmap[s] for s in three_sym if s in gmap]
    raw = mtx[keep][:, sym_idx]  # (n_cells, 13927) sparse
    tt = tgt[keep].to_numpy()
    print(f"[{time.time()-t0:.0f}s] Gene-filtered shape: {raw.shape}", flush=True)

    def centroids(X_mat):
        out = np.zeros((13, X_mat.shape[1]))
        for i, t in enumerate(TYPES_13):
            m = (tt == t)
            if m.sum() == 0: continue
            out[i] = np.asarray(X_mat[m].mean(axis=0)).flatten()
        return out

    def test(Ymat, label):
        Xc = X - X.mean(axis=0); Yc = Ymat - Ymat.mean(axis=0)
        raw_info = f"||X_c||²={((Xc)**2).sum():.2f} ||Y_c||²={((Yc)**2).sum():.2f}"
        # PCA at 95%
        stacked = np.vstack([X, Ymat])
        pca = PCA(n_components=0.95, svd_solver="full", random_state=RANDOM_SEED)
        proj = pca.fit_transform(stacked)
        Xp, Yp = proj[:13], proj[13:]
        r = procrustes_align(Xp, Yp)
        # Use 10000 perms for a solid number
        p, null = permutation_test(Xp, Yp, n_permutations=10000, seed=RANDOM_SEED)
        obs_null = r.distance / np.median(null)
        print(f"\n  >> {label}", flush=True)
        print(f"     Y stats: min={Ymat.min():.3f} max={Ymat.max():.3f} mean={Ymat.mean():.3f} std={Ymat.std():.3f}", flush=True)
        print(f"     Raw: {raw_info}", flush=True)
        print(f"     PCA k={pca.n_components_}, d={r.distance:.3f}, s={r.scaling:.6f}, "
              f"obs/null(med)={obs_null:.4f}  [target 0.7488, Δ={abs(obs_null-0.7488):.4f}]", flush=True)

    # Test 1: as-is (data already normalized)
    print(f"\n[{time.time()-t0:.0f}s] Test 1: use matrix values as-is (no normalization)", flush=True)
    Y = centroids(raw)
    test(Y, "as-is (no normalization)")

    # Test 2: log1p only
    print(f"\n[{time.time()-t0:.0f}s] Test 2: log1p only", flush=True)
    Ylog = centroids(raw.copy().log1p())
    test(Ylog, "log1p only")

    # Test 3: normalize_total + log1p (my original)
    print(f"\n[{time.time()-t0:.0f}s] Test 3: normalize_total(1e4)+log1p", flush=True)
    a = ad.AnnData(X=raw.copy())
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    test(centroids(a.X), "normalize_total(1e4)+log1p (original)")

    # Test 4: normalize_total + log1p with target_sum=1e6
    print(f"\n[{time.time()-t0:.0f}s] Test 4: normalize_total(1e6)+log1p", flush=True)
    a = ad.AnnData(X=raw.copy())
    sc.pp.normalize_total(a, target_sum=1e6); sc.pp.log1p(a)
    test(centroids(a.X), "normalize_total(1e6)+log1p")


if __name__ == "__main__":
    main()
