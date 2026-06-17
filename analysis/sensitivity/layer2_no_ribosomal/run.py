#!/usr/bin/env python3
"""
Layer-2 ribosomal-exclusion sensitivity.

Refits the primary joint PCA and the Layer-2 covariance / Krzanowski-S
analysis after dropping all RPL/RPS genes from the ortholog space.

Methodology mirrors scripts/t3b_ellipsoid_alignment.py exactly, with two
differences:
  1. Input gene space is filtered to drop ribosomal genes via regex
     ^(RPL|RPS|Rpl|Rps)\\d on the gene-symbol side of the ortholog map.
  2. The Procrustes rotation R is recomputed in the new joint-PCA basis
     (not loaded from procrustes_results_35.json, which used the full
     gene space).

Outputs:
  - results.json: pre/post Krzanowski S + null permutations, CPC1 drivers,
    cell counts, PCA component count
  - cpc1_drivers.csv: per-cell-type old vs new rank-1 CPC1 driver
"""
from __future__ import annotations

import json
import re
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.decomposition import PCA

PROJECT = Path(__file__).resolve().parents[3]
OUT_DIR = Path(__file__).resolve().parent
DATA_SCALED = PROJECT / "data" / "phase2_scaled"
PHASE2_DIR = PROJECT / "output" / "phase2" / "scaled_35types"
ELLIPSOID_DIR = PROJECT / "output" / "mechanistic" / "ellipsoid_alignment"
ORTHOLOGS = PROJECT / "data" / "phase1" / "orthologs_human_mouse.csv"

SEED = 42
N_PERM = 10_000
K_VALUES = [1, 3, 5]
PCA_VARIANCE = 0.95
RIBO_RE = re.compile(r"^(RPL|RPS|Rpl|Rps)\d")


def subspace_overlap(V1: np.ndarray, V2: np.ndarray, k: int) -> float:
    A, B = V1[:, :k], V2[:, :k]
    M = A.T @ B
    return float(np.trace(M.T @ M))


def procrustes_R(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """Return the rotation matrix that aligns Y onto X (no reflection)."""
    X_c = X - X.mean(axis=0)
    Y_c = Y - Y.mean(axis=0)
    M = X_c.T @ Y_c
    U, _, Vt = np.linalg.svd(M)
    V = Vt.T
    D = np.eye(M.shape[0])
    if np.linalg.det(V @ U.T) < 0:
        D[-1, -1] = -1.0
    return V @ D @ U.T


def main() -> None:
    t0 = time.time()
    print("=" * 70)
    print("Layer-2 ribosomal-exclusion sensitivity")
    print("=" * 70)

    # A1.1: identify ribosomal genes
    orth = pd.read_csv(ORTHOLOGS)
    ribo_mask = orth["human_gene_name"].astype(str).str.match(RIBO_RE) | \
                orth["mouse_gene_name"].astype(str).str.match(RIBO_RE)
    ribo_ensembl = set(orth.loc[ribo_mask, "human_ensembl_id"])

    hc = pd.read_csv(PHASE2_DIR / "centroids_human_35.csv", index_col=0)
    mc = pd.read_csv(PHASE2_DIR / "centroids_mouse_35.csv", index_col=0)
    saved = np.load(PHASE2_DIR / "pca_centroids_35.npz", allow_pickle=True)
    cell_types = list(saved["cell_types"])
    hc = hc.loc[cell_types]
    mc = mc.loc[cell_types]

    all_genes = list(hc.columns)
    ribo_in_space = sorted(g for g in all_genes if g in ribo_ensembl)
    kept_genes = [g for g in all_genes if g not in ribo_ensembl]
    print(f"\nRibosomal genes excluded: {len(ribo_in_space)} of {len(all_genes)}")
    print(f"Kept gene count: {len(kept_genes)}")

    # A1.2: refit joint PCA at 95% variance on the reduced gene space
    hc_red = hc[kept_genes]
    mc_red = mc[kept_genes]
    combined = np.vstack([hc_red.values, mc_red.values])  # (70, G_red)
    pca = PCA(n_components=PCA_VARIANCE, svd_solver="full", random_state=SEED)
    combined_pca = pca.fit_transform(combined)
    n_pc = pca.n_components_
    cumvar = float(np.sum(pca.explained_variance_ratio_))
    print(f"\nReduced PCA: {n_pc} components, {cumvar*100:.1f}% variance")
    print(f"  (primary analysis: 33 components at 95.2% variance)")

    H_pca = combined_pca[: len(cell_types)]
    M_pca = combined_pca[len(cell_types) :]

    # Recompute rotation in the new basis (centroid-optimal R for Layer 1)
    R_new = procrustes_R(H_pca, M_pca)
    print(f"R recomputed in reduced basis: shape={R_new.shape}, "
          f"det={np.linalg.det(R_new):+.4f}")

    # A1.3: project single cells through new PCA
    print("\nProjecting single cells into reduced PCA basis ...")
    cells = {}
    for species, h5 in [("human", "human_scaled.h5ad"),
                        ("mouse", "mouse_scaled.h5ad")]:
        adata = ad.read_h5ad(DATA_SCALED / h5)
        # Restrict to kept_genes order
        X = adata[:, kept_genes].X
        if sp.issparse(X):
            X = X.toarray()
        X = np.asarray(X, dtype=np.float32)
        X_pca = (X - pca.mean_) @ pca.components_.T
        cells[species] = {
            "X_pca": X_pca,
            "cell_types": adata.obs["cell_type"].values,
        }
        print(f"  {species}: {X_pca.shape[0]:,} cells projected")

    # A1.3 continued: per-type covariance + eigendecomposition
    print("\nComputing per-type covariance + eigendecomposition ...")
    eigen = {}
    for species in ("human", "mouse"):
        ct_arr = cells[species]["cell_types"]
        Xp = cells[species]["X_pca"]
        for ct in cell_types:
            mask = ct_arr == ct
            cells_arr = Xp[mask]
            n = cells_arr.shape[0]
            if n < 5:
                print(f"  ⚠ {species}/{ct}: only {n} cells")
                continue
            centroid = cells_arr.mean(axis=0)
            centered = cells_arr - centroid
            cov = (centered.T @ centered) / (n - 1)
            eigvals, eigvecs = np.linalg.eigh(cov)
            idx = np.argsort(eigvals)[::-1]
            eigen[(species, ct)] = {
                "n": n,
                "cov": cov,
                "eigvals": eigvals[idx],
                "eigvecs": eigvecs[:, idx],
            }

    # A1.6: Krzanowski S statistics, pre and post rotation, with nulls
    print("\nComputing Krzanowski S (pre/post rotation) + null permutations ...")
    V_H = [eigen[("human", ct)]["eigvecs"] for ct in cell_types]
    V_M = [eigen[("mouse", ct)]["eigvecs"] for ct in cell_types]
    V_M_rot = [R_new @ v for v in V_M]

    obs_pre = {k: np.mean([subspace_overlap(V_H[j], V_M[j], k) / k
                           for j in range(len(cell_types))]) for k in K_VALUES}
    obs_post = {k: np.mean([subspace_overlap(V_H[j], V_M_rot[j], k) / k
                            for j in range(len(cell_types))]) for k in K_VALUES}

    rng = np.random.RandomState(SEED)
    null_pre = {k: np.zeros(N_PERM) for k in K_VALUES}
    null_post = {k: np.zeros(N_PERM) for k in K_VALUES}
    t_perm = time.time()
    for i in range(N_PERM):
        perm = rng.permutation(len(cell_types))
        for k in K_VALUES:
            pre_s, post_s = [], []
            for j in range(len(cell_types)):
                pre_s.append(subspace_overlap(V_H[j], V_M[perm[j]], k) / k)
                post_s.append(subspace_overlap(V_H[j], V_M_rot[perm[j]], k) / k)
            null_pre[k][i] = np.mean(pre_s)
            null_post[k][i] = np.mean(post_s)
    print(f"  Permutation runtime: {time.time()-t_perm:.0f}s")

    pre_summary = {}
    post_summary = {}
    for k in K_VALUES:
        p_pre = float((np.sum(null_pre[k] >= obs_pre[k]) + 1) / (N_PERM + 1))
        p_post = float((np.sum(null_post[k] >= obs_post[k]) + 1) / (N_PERM + 1))
        pre_summary[k] = {
            "observed_S": float(obs_pre[k]),
            "null_mean": float(null_pre[k].mean()),
            "null_std": float(null_pre[k].std()),
            "p_value": p_pre,
        }
        post_summary[k] = {
            "observed_S": float(obs_post[k]),
            "null_mean": float(null_post[k].mean()),
            "null_std": float(null_post[k].std()),
            "p_value": p_post,
        }

    # A1.4-A1.5: CPC1 + project back to reduced gene space
    print("\nComputing CPC1 + projecting to reduced gene space ...")
    W = pca.components_  # (n_pc, G_red)
    kept_gene_arr = np.array(kept_genes)

    # Load existing 35-type rank-1 drivers (full-space; ribosomal-dominated)
    with open(ELLIPSOID_DIR / "cpc_genes.json") as f:
        old_cpc = json.load(f)
    old_35 = old_cpc.get("35type", {})

    # Gene symbol mapping
    ens_to_sym = dict(zip(orth["human_ensembl_id"], orth["human_gene_name"]))

    cpc_rows = []
    for ct in cell_types:
        if ("human", ct) not in eigen or ("mouse", ct) not in eigen:
            continue
        h = eigen[("human", ct)]
        m = eigen[("mouse", ct)]
        weighted_cov = h["n"] * h["cov"] + m["n"] * m["cov"]
        evals, evecs = np.linalg.eigh(weighted_cov)
        idx = np.argsort(evals)[::-1]
        cpc1 = evecs[:, idx[0]]
        cpc1_genes = cpc1 @ W  # (G_red,)
        abs_load = np.abs(cpc1_genes)
        top5 = np.argsort(abs_load)[::-1][:5]
        new_top1_id = kept_gene_arr[top5[0]]
        new_top1_sym = ens_to_sym.get(new_top1_id, new_top1_id)
        new_top5 = [(kept_gene_arr[i], ens_to_sym.get(kept_gene_arr[i], kept_gene_arr[i]),
                     float(cpc1_genes[i])) for i in top5]

        old_entry = old_35.get(ct, {}).get("top5_genes", [])
        old_top1_id = (old_entry[0]["gene"] if isinstance(old_entry[0], dict)
                       else old_entry[0][0]) if old_entry else None
        old_top1_sym = ens_to_sym.get(old_top1_id, old_top1_id) if old_top1_id else "—"
        was_ribosomal = bool(old_top1_sym and RIBO_RE.match(str(old_top1_sym)))

        cpc_rows.append({
            "cell_type": ct,
            "old_top1_ensembl": old_top1_id,
            "old_top1_symbol": old_top1_sym,
            "old_top1_was_ribosomal": was_ribosomal,
            "new_top1_ensembl": new_top1_id,
            "new_top1_symbol": new_top1_sym,
            "new_top1_loading": float(cpc1_genes[top5[0]]),
            "new_top5": new_top5,
        })

    cpc_df = pd.DataFrame([{k: v for k, v in r.items() if k != "new_top5"}
                           for r in cpc_rows])
    cpc_df.to_csv(OUT_DIR / "cpc1_drivers.csv", index=False)
    print(f"  Saved {OUT_DIR/'cpc1_drivers.csv'}")

    # Save full results
    results = {
        "metadata": {
            "script": str(Path(__file__).resolve().relative_to(PROJECT)),
            "seed": SEED,
            "n_permutations": N_PERM,
            "k_values": K_VALUES,
            "pca_variance_threshold": PCA_VARIANCE,
            "ribosomal_regex": RIBO_RE.pattern,
            "runtime_sec": time.time() - t0,
        },
        "ribosomal_filter": {
            "n_excluded_genes": len(ribo_in_space),
            "n_genes_in_full_space": len(all_genes),
            "n_kept_genes": len(kept_genes),
        },
        "joint_pca_reduced": {
            "n_components": int(n_pc),
            "cumulative_variance": cumvar,
            "primary_comparison": {
                "n_components": 33,
                "cumulative_variance": 0.952,
            },
        },
        "rotation": {
            "shape": list(R_new.shape),
            "determinant": float(np.linalg.det(R_new)),
        },
        "krzanowski_pre_rotation": pre_summary,
        "krzanowski_post_rotation": post_summary,
        "primary_comparison": {
            "pre_rotation": {"k=5": {"S": 0.483, "null_mean": 0.375, "p": "<1e-4"}},
            "post_rotation": {"k=5": {"S": 0.230, "null_mean": 0.180, "p": "<1e-4"}},
        },
        "cpc1_drivers": cpc_rows,
        "cpc1_summary": {
            "n_originally_ribosomal_top1": int(sum(r["old_top1_was_ribosomal"] for r in cpc_rows)),
            "n_new_top1_still_ribosomal": int(sum(
                bool(RIBO_RE.match(str(r["new_top1_symbol"])))
                for r in cpc_rows)),
        },
    }
    out_path = OUT_DIR / "results.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\n✓ Saved {out_path}")
    print(f"Total runtime: {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
