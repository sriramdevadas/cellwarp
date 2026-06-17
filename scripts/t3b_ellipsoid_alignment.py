#!/usr/bin/env python3
"""
T3-B: Covariance Ellipsoid Alignment Analysis

Tests whether the SHAPE of within-cell-type variation (not just centroid position)
is conserved across species after Procrustes rotation. Implements Krzanowski (1984)
subspace similarity and Common Principal Components (Flury 1988) analysis in the
CellWarp cross-species single-cell framework.

Biology: Each cell type is a cloud of cells in gene expression space. That cloud has
an ellipsoidal shape defined by its covariance matrix. If cell type identity is deeply
constrained, the cloud shape itself should be conserved across species — human and
mouse cells of the same type spread in similar directions.

Math: Subspace overlap S(k) = trace(V_H' V_M V_M' V_H) where V_H, V_M are top-k
eigenvectors of within-type covariance matrices. S(k)/k ∈ [0,1] measures alignment
of the top-k principal axes of variation.

References:
- Krzanowski (1984) "Between-groups comparison of principal components" JASA
- Flury (1988) "Common Principal Components and Related Multivariate Models"
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
from sklearn.decomposition import PCA

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT = Path(__file__).resolve().parent.parent
DATA_PHASE1 = PROJECT / "data" / "phase1"
DATA_SCALED = PROJECT / "data" / "phase2_scaled"
OUTPUT_SCALED = PROJECT / "output" / "phase2" / "scaled_35types"
OUTPUT_DIR = PROJECT / "output" / "mechanistic" / "ellipsoid_alignment"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

SEED = 42
N_PERM = 10_000
K_VALUES = [1, 3, 5]

SIX_TYPES = [
    "B cell",
    "CD4-positive, alpha-beta T cell",
    "CD8-positive, alpha-beta T cell",
    "endothelial cell",
    "hepatocyte",
    "macrophage",
]

# ---------------------------------------------------------------------------
# Step 0: Reconstruct PCA model from saved centroids
# ---------------------------------------------------------------------------
def reconstruct_pca():
    """Refit PCA on 70 combined centroids to get components_ and mean_."""
    print("=" * 70)
    print("STEP 0: Reconstructing PCA model from centroids")
    print("=" * 70)

    saved = np.load(OUTPUT_SCALED / "pca_centroids_35.npz", allow_pickle=True)
    cell_types = list(saved["cell_types"])
    saved_human = saved["human"]  # (35, 33)
    saved_mouse = saved["mouse"]  # (35, 33)

    # Load full-gene centroids
    hc = pd.read_csv(OUTPUT_SCALED / "centroids_human_35.csv", index_col=0)
    mc = pd.read_csv(OUTPUT_SCALED / "centroids_mouse_35.csv", index_col=0)

    # Order by cell_types
    hc = hc.loc[cell_types]
    mc = mc.loc[cell_types]
    gene_names = list(hc.columns)

    combined = np.vstack([hc.values, mc.values])  # (70, 16959)

    pca = PCA(n_components=33, svd_solver="full", random_state=SEED)
    combined_pca = pca.fit_transform(combined)

    # Verify match
    err_h = np.max(np.abs(combined_pca[:35] - saved_human))
    err_m = np.max(np.abs(combined_pca[35:] - saved_mouse))
    print(f"  PCA reconstruction max error: human={err_h:.2e}, mouse={err_m:.2e}")
    assert err_h < 0.01 and err_m < 0.01, "PCA reconstruction mismatch!"

    print(f"  PCA: {pca.n_components_} components, "
          f"{sum(pca.explained_variance_ratio_)*100:.1f}% variance")
    print(f"  Gene space: {len(gene_names)} genes")

    return pca, cell_types, gene_names


# ---------------------------------------------------------------------------
# Step 1: Project individual cells into PCA space
# ---------------------------------------------------------------------------
def project_cells(pca, gene_names, scale="35type"):
    """Project all individual cells into 33-D PCA space."""
    print("\n" + "=" * 70)
    print(f"STEP 1: Projecting individual cells into PCA space ({scale})")
    print("=" * 70)

    if scale == "35type":
        h_path = DATA_SCALED / "human_scaled.h5ad"
        m_path = DATA_SCALED / "mouse_scaled.h5ad"
    else:
        h_path = DATA_PHASE1 / "human_qc.h5ad"
        m_path = DATA_PHASE1 / "mouse_qc.h5ad"

    results = {}
    for species, path in [("human", h_path), ("mouse", m_path)]:
        adata = ad.read_h5ad(path)
        # Ensure gene order matches PCA
        assert set(gene_names).issubset(set(adata.var_names)), \
            f"Gene mismatch for {species}"
        X = adata[:, gene_names].X
        if sp.issparse(X):
            X = X.toarray()
        X = np.asarray(X, dtype=np.float32)

        # Project to PCA space
        X_pca = (X - pca.mean_) @ pca.components_.T  # (n_cells, 33)

        ct_col = adata.obs["cell_type"].values
        if scale == "6type":
            mask = np.isin(ct_col, SIX_TYPES)
            X_pca = X_pca[mask]
            ct_col = ct_col[mask]

        results[species] = {"X_pca": X_pca, "cell_types": ct_col}

        # Report counts
        unique, counts = np.unique(ct_col, return_counts=True)
        print(f"\n  {species}: {len(ct_col)} cells, {len(unique)} types")
        for ct, n in sorted(zip(unique, counts), key=lambda x: -x[1]):
            flag = " ⚠️ <100" if n < 100 else ""
            print(f"    {ct:<50} {n:>6}{flag}")

    return results


# ---------------------------------------------------------------------------
# Step 2-3: Covariance matrices and eigendecomposition
# ---------------------------------------------------------------------------
def compute_covariance_eigen(cell_data, scale_label):
    """Compute within-type covariance matrices and eigendecompose."""
    print("\n" + "=" * 70)
    print(f"STEPS 2-3: Covariance matrices & eigendecomposition ({scale_label})")
    print("=" * 70)

    results = {}
    for species in ["human", "mouse"]:
        X_pca = cell_data[species]["X_pca"]
        ct_col = cell_data[species]["cell_types"]
        unique_cts = sorted(np.unique(ct_col))

        for ct in unique_cts:
            mask = ct_col == ct
            cells = X_pca[mask]  # (n, 33)
            n = cells.shape[0]

            # Center on zero
            centroid = cells.mean(axis=0)
            centered = cells - centroid

            # Covariance matrix
            cov = (centered.T @ centered) / (n - 1)  # (33, 33)
            cond = np.linalg.cond(cov)

            # Eigendecomposition (sorted descending)
            eigvals, eigvecs = np.linalg.eigh(cov)
            idx = np.argsort(eigvals)[::-1]
            eigvals = eigvals[idx]
            eigvecs = eigvecs[:, idx]

            # Variance explained fractions
            total_var = eigvals.sum()
            var_frac = eigvals / total_var if total_var > 0 else eigvals

            key = (species, ct)
            results[key] = {
                "n_cells": n,
                "centroid": centroid,
                "cov": cov,
                "eigvals": eigvals,
                "eigvecs": eigvecs,
                "cond": cond,
                "var_frac": var_frac,
            }

    # Report
    cts = sorted(set(ct for _, ct in results.keys()))
    print(f"\n  {'Cell Type':<50} {'Species':<8} {'N':>6} "
          f"{'Cond':>10} {'Top1%':>7} {'Top3%':>7} {'Top5%':>7}")
    print("  " + "-" * 100)
    for ct in cts:
        for sp_name in ["human", "mouse"]:
            r = results[(sp_name, ct)]
            vf = r["var_frac"]
            cond_str = f"{r['cond']:.0f}" if r['cond'] < 1e6 else f"{r['cond']:.1e}"
            flag = " ⚠️" if r["cond"] > 1000 else ""
            n_flag = " LOW" if r["n_cells"] < 100 else ""
            print(f"  {ct:<50} {sp_name:<8} {r['n_cells']:>6}{n_flag} "
                  f"{cond_str:>10}{flag} {vf[0]*100:>6.1f}% "
                  f"{sum(vf[:3])*100:>6.1f}% {sum(vf[:5])*100:>6.1f}%")

    return results


# ---------------------------------------------------------------------------
# Step 4-5: Subspace overlap (pre and post Procrustes)
# ---------------------------------------------------------------------------
def subspace_overlap(V1, V2, k):
    """Krzanowski (1984) subspace similarity: S(k) = trace(V1' V2 V2' V1)."""
    A = V1[:, :k]
    B = V2[:, :k]
    M = A.T @ B  # (k, k)
    return np.trace(M.T @ M)  # = ||A'B||_F^2


def compute_alignment_scores(eigen_results, rotation_matrix, scale_label):
    """Compute pre- and post-Procrustes ellipsoid alignment."""
    print("\n" + "=" * 70)
    print(f"STEPS 4-5: Ellipsoid alignment pre/post Procrustes ({scale_label})")
    print("=" * 70)

    cts = sorted(set(ct for sp, ct in eigen_results.keys() if sp == "human"))
    R = rotation_matrix  # 33×33

    records = []
    for ct in cts:
        h = eigen_results[("human", ct)]
        m = eigen_results[("mouse", ct)]
        V_H = h["eigvecs"]  # (33, 33)
        V_M = m["eigvecs"]  # (33, 33)

        # Rotate mouse eigenvectors by Procrustes R
        V_M_rot = R @ V_M  # (33, 33)

        for k in K_VALUES:
            s_pre = subspace_overlap(V_H, V_M, k) / k
            s_post = subspace_overlap(V_H, V_M_rot, k) / k
            records.append({
                "cell_type": ct,
                "k": k,
                "S_pre": s_pre,
                "S_post": s_post,
                "improvement": s_post - s_pre,
            })

    df = pd.DataFrame(records)

    # Report
    print(f"\n  {'Cell Type':<50} {'k':>3} {'Pre':>7} {'Post':>7} {'Δ':>7}")
    print("  " + "-" * 80)
    for ct in cts:
        for k in K_VALUES:
            row = df[(df.cell_type == ct) & (df.k == k)].iloc[0]
            print(f"  {ct:<50} {k:>3} {row.S_pre:>7.3f} {row.S_post:>7.3f} "
                  f"{row.improvement:>+7.3f}")

    # Summary
    for k in K_VALUES:
        sub = df[df.k == k]
        print(f"\n  k={k}: mean pre={sub.S_pre.mean():.3f}, "
              f"mean post={sub.S_post.mean():.3f}, "
              f"mean Δ={sub.improvement.mean():+.3f}")

    return df


# ---------------------------------------------------------------------------
# Step 6: Label shuffle permutation test
# ---------------------------------------------------------------------------
def perm_test_label_shuffle(eigen_results, rotation_matrix, n_perm=N_PERM):
    """Shuffle cell type labels — test BOTH pre and post-Procrustes alignment."""
    print("\n" + "=" * 70)
    print(f"STEP 6: Label shuffle permutation test ({n_perm} iterations)")
    print("=" * 70)

    cts = sorted(set(ct for sp, ct in eigen_results.keys() if sp == "human"))
    R = rotation_matrix

    # Pre-extract eigenvectors
    V_H_list = [eigen_results[("human", ct)]["eigvecs"] for ct in cts]
    V_M_list = [eigen_results[("mouse", ct)]["eigvecs"] for ct in cts]
    V_M_rot_list = [R @ v for v in V_M_list]
    n_types = len(cts)

    # Observed: both pre and post-Procrustes
    obs_pre, obs_post = {}, {}
    for k in K_VALUES:
        pre_scores, post_scores = [], []
        for j in range(n_types):
            pre_scores.append(subspace_overlap(V_H_list[j], V_M_list[j], k) / k)
            post_scores.append(subspace_overlap(V_H_list[j], V_M_rot_list[j], k) / k)
        obs_pre[k] = np.mean(pre_scores)
        obs_post[k] = np.mean(post_scores)

    rng = np.random.RandomState(SEED)
    null_pre = {k: np.zeros(n_perm) for k in K_VALUES}
    null_post = {k: np.zeros(n_perm) for k in K_VALUES}

    t0 = time.time()
    for i in range(n_perm):
        perm = rng.permutation(n_types)
        for k in K_VALUES:
            pre_s, post_s = [], []
            for j in range(n_types):
                pre_s.append(
                    subspace_overlap(V_H_list[j], V_M_list[perm[j]], k) / k)
                post_s.append(
                    subspace_overlap(V_H_list[j], V_M_rot_list[perm[j]], k) / k)
            null_pre[k][i] = np.mean(pre_s)
            null_post[k][i] = np.mean(post_s)
    elapsed = time.time() - t0
    print(f"  Runtime: {elapsed:.1f}s")

    results_pre, results_post = {}, {}
    print("\n  PRE-PROCRUSTES (raw alignment):")
    for k in K_VALUES:
        n_geq = int(np.sum(null_pre[k] >= obs_pre[k])) + 1
        p = n_geq / (n_perm + 1)
        results_pre[k] = {"observed": obs_pre[k], "null": null_pre[k], "p_value": p}
        print(f"    k={k}: observed={obs_pre[k]:.4f}, null mean={null_pre[k].mean():.4f}, "
              f"null std={null_pre[k].std():.4f}, p={p:.4f}")

    print("\n  POST-PROCRUSTES (after rotation):")
    for k in K_VALUES:
        n_geq = int(np.sum(null_post[k] >= obs_post[k])) + 1
        p = n_geq / (n_perm + 1)
        results_post[k] = {"observed": obs_post[k], "null": null_post[k], "p_value": p}
        print(f"    k={k}: observed={obs_post[k]:.4f}, null mean={null_post[k].mean():.4f}, "
              f"null std={null_post[k].std():.4f}, p={p:.4f}")

    return results_pre, results_post


# ---------------------------------------------------------------------------
# Step 7: Random rotation permutation test
# ---------------------------------------------------------------------------
def random_orthogonal(n, rng):
    """Sample uniformly from SO(n) via QR of Gaussian matrix."""
    Z = rng.randn(n, n)
    Q, R_ = np.linalg.qr(Z)
    # Make deterministic signs
    d = np.diag(R_)
    Q = Q * np.sign(d)
    # Ensure SO(n) (det=+1)
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    return Q


def perm_test_random_rotation(eigen_results, rotation_matrix, n_perm=N_PERM):
    """Compare Procrustes rotation to random rotations for ellipsoid alignment."""
    print("\n" + "=" * 70)
    print(f"STEP 7: Random rotation permutation test ({n_perm} iterations)")
    print("=" * 70)

    cts = sorted(set(ct for sp, ct in eigen_results.keys() if sp == "human"))
    R = rotation_matrix
    dim = R.shape[0]

    # Observed: post-Procrustes
    obs = {}
    for k in K_VALUES:
        scores = []
        for ct in cts:
            V_H = eigen_results[("human", ct)]["eigvecs"]
            V_M_rot = R @ eigen_results[("mouse", ct)]["eigvecs"]
            scores.append(subspace_overlap(V_H, V_M_rot, k) / k)
        obs[k] = np.mean(scores)

    rng = np.random.RandomState(SEED + 1)
    null = {k: np.zeros(n_perm) for k in K_VALUES}

    V_H_list = [eigen_results[("human", ct)]["eigvecs"] for ct in cts]
    V_M_list = [eigen_results[("mouse", ct)]["eigvecs"] for ct in cts]
    n_types = len(cts)

    t0 = time.time()
    for i in range(n_perm):
        Q = random_orthogonal(dim, rng)
        for k in K_VALUES:
            scores = []
            for j in range(n_types):
                V_M_rot = Q @ V_M_list[j]
                scores.append(subspace_overlap(V_H_list[j], V_M_rot, k) / k)
            null[k][i] = np.mean(scores)
        if (i + 1) % 2000 == 0:
            print(f"  ... {i+1}/{n_perm} ({time.time()-t0:.0f}s)")
    elapsed = time.time() - t0
    print(f"  Runtime: {elapsed:.1f}s")

    results = {}
    for k in K_VALUES:
        n_geq = int(np.sum(null[k] >= obs[k])) + 1
        p = n_geq / (n_perm + 1)
        results[k] = {"observed": obs[k], "null": null[k], "p_value": p}
        print(f"  k={k}: observed={obs[k]:.4f}, null mean={null[k].mean():.4f}, "
              f"null std={null[k].std():.4f}, p={p:.4f}")

    return results


# ---------------------------------------------------------------------------
# Step 8: Correlation with rigidity
# ---------------------------------------------------------------------------
def correlate_with_rigidity(alignment_df, residuals_df, scale_label):
    """Spearman correlation of ellipsoid alignment vs Procrustes residual."""
    print("\n" + "=" * 70)
    print(f"STEP 8: Ellipsoid alignment vs rigidity correlation ({scale_label})")
    print("=" * 70)

    merged = alignment_df.merge(
        residuals_df[["cell_type", "residual_magnitude"]],
        on="cell_type", how="inner"
    )

    records = []
    for k in K_VALUES:
        for metric, label in [("S_pre", "pre"), ("S_post", "post")]:
            sub = merged[merged.k == k]
            rho, p = stats.spearmanr(sub[metric], sub["residual_magnitude"])
            records.append({
                "k": k, "metric": label, "rho": rho, "p": p,
            })

    df = pd.DataFrame(records)
    # FDR correction (Benjamini-Hochberg)
    from scipy.stats import false_discovery_control
    # Manual BH
    ps = df["p"].values
    n = len(ps)
    sorted_idx = np.argsort(ps)
    q = np.zeros(n)
    for rank_i, idx in enumerate(sorted_idx):
        q[idx] = ps[idx] * n / (rank_i + 1)
    # Enforce monotonicity
    for i in range(n - 2, -1, -1):
        q[sorted_idx[i]] = min(q[sorted_idx[i]], q[sorted_idx[i + 1]] if i < n - 1 else 1.0)
    q = np.minimum(q, 1.0)
    df["q_fdr"] = q

    print(f"\n  {'k':>3} {'Metric':<6} {'ρ':>8} {'p':>10} {'q_FDR':>10}")
    print("  " + "-" * 40)
    for _, row in df.iterrows():
        sig = " *" if row.q_fdr < 0.05 else ""
        print(f"  {row.k:>3} {row.metric:<6} {row.rho:>+8.3f} {row.p:>10.4f} "
              f"{row.q_fdr:>10.4f}{sig}")

    return df


# ---------------------------------------------------------------------------
# Step 9: Eigenvalue ratio conservation
# ---------------------------------------------------------------------------
def eigenvalue_conservation(eigen_results, residuals_df, scale_label):
    """Test whether eigenvalue magnitude structure is conserved.

    Uses Pearson correlation on normalized eigenvalue profiles (not Spearman,
    since both vectors are sorted descending — Spearman would always be 1.0).
    Also computes cosine similarity as a secondary metric.
    """
    print("\n" + "=" * 70)
    print(f"STEP 9: Eigenvalue ratio conservation ({scale_label})")
    print("=" * 70)

    cts = sorted(set(ct for sp, ct in eigen_results.keys() if sp == "human"))
    top_k = 10

    records = []
    for ct in cts:
        h_vals = eigen_results[("human", ct)]["eigvals"][:top_k]
        m_vals = eigen_results[("mouse", ct)]["eigvals"][:top_k]
        h_norm = h_vals / h_vals.sum() if h_vals.sum() > 0 else h_vals
        m_norm = m_vals / m_vals.sum() if m_vals.sum() > 0 else m_vals

        # Pearson correlation on normalized eigenvalue profiles
        r_pearson, p_pearson = stats.pearsonr(h_norm, m_norm)

        # Cosine similarity
        cos_sim = float(np.dot(h_norm, m_norm) /
                       (np.linalg.norm(h_norm) * np.linalg.norm(m_norm)))

        # L2 distance between normalized profiles
        l2_dist = float(np.linalg.norm(h_norm - m_norm))

        records.append({
            "cell_type": ct,
            "pearson_r": r_pearson,
            "pearson_p": p_pearson,
            "cosine_sim": cos_sim,
            "l2_distance": l2_dist,
            "human_top1_frac": float(h_norm[0]),
            "mouse_top1_frac": float(m_norm[0]),
            "human_norm": h_norm.tolist(),
            "mouse_norm": m_norm.tolist(),
        })

    df = pd.DataFrame(records)

    print(f"\n  {'Cell Type':<50} {'Pearson r':>10} {'p':>10} {'Cosine':>8} {'L2':>8}")
    print("  " + "-" * 90)
    for _, row in df.iterrows():
        print(f"  {row.cell_type:<50} {row.pearson_r:>+10.3f} {row.pearson_p:>10.4f} "
              f"{row.cosine_sim:>8.3f} {row.l2_distance:>8.3f}")

    print(f"\n  Mean Pearson r: {df.pearson_r.mean():.3f}")
    print(f"  Mean cosine similarity: {df.cosine_sim.mean():.3f}")

    # Correlate eigenvalue conservation with rigidity
    merged = df.merge(
        residuals_df[["cell_type", "residual_magnitude"]],
        on="cell_type", how="inner"
    )
    if len(merged) > 5:
        rho_rig, p_rig = stats.spearmanr(
            merged["pearson_r"], merged["residual_magnitude"]
        )
        print(f"\n  Eigenvalue conservation (Pearson r) vs rigidity: "
              f"ρ={rho_rig:+.3f}, p={p_rig:.4f}")
    else:
        rho_rig, p_rig = np.nan, np.nan

    return df, rho_rig, p_rig


# ---------------------------------------------------------------------------
# Step 10: Krzanowski Common Principal Components
# ---------------------------------------------------------------------------
def common_principal_components(eigen_results, pca_model, gene_names, scale_label):
    """Krzanowski (1984) CPC approximation."""
    print("\n" + "=" * 70)
    print(f"STEP 10: Common Principal Components ({scale_label})")
    print("=" * 70)

    cts = sorted(set(ct for sp, ct in eigen_results.keys() if sp == "human"))
    W = pca_model.components_  # (33, 16959)

    cpc_results = {}
    for ct in cts:
        h = eigen_results[("human", ct)]
        m = eigen_results[("mouse", ct)]

        # Weighted sum of covariance matrices
        n_h = h["n_cells"]
        n_m = m["n_cells"]
        weighted_cov = n_h * h["cov"] + n_m * m["cov"]

        # Eigendecompose weighted sum → CPCs
        eigvals, eigvecs = np.linalg.eigh(weighted_cov)
        idx = np.argsort(eigvals)[::-1]
        eigvals = eigvals[idx]
        eigvecs = eigvecs[:, idx]

        # How much variance does top CPC explain in each species separately?
        cpc1 = eigvecs[:, 0]  # (33,)
        var_h = cpc1 @ h["cov"] @ cpc1
        var_m = cpc1 @ m["cov"] @ cpc1
        total_var_h = np.trace(h["cov"])
        total_var_m = np.trace(m["cov"])
        frac_h = var_h / total_var_h if total_var_h > 0 else 0
        frac_m = var_m / total_var_m if total_var_m > 0 else 0

        # Project CPC back to gene space
        cpc1_genes = cpc1 @ W  # (16959,)
        abs_loadings = np.abs(cpc1_genes)
        top5_idx = np.argsort(abs_loadings)[::-1][:5]
        top5_genes = [(gene_names[i], float(cpc1_genes[i])) for i in top5_idx]

        cpc_results[ct] = {
            "cpc1_var_frac_human": float(frac_h),
            "cpc1_var_frac_mouse": float(frac_m),
            "top5_genes": top5_genes,
            "cpc1_eigval": float(eigvals[0]),
            "total_eigval": float(eigvals.sum()),
        }

        print(f"\n  {ct}:")
        print(f"    CPC1 variance: human={frac_h:.1%}, mouse={frac_m:.1%}")
        print(f"    Top 5 genes: {', '.join(g for g, _ in top5_genes)}")

    return cpc_results


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def make_plots(alignment_df, residuals_df, perm_label_pre, perm_label_post,
               perm_rotation, eigen_results, eigenval_df, scale_label):
    """Generate all figures."""
    print("\n" + "=" * 70)
    print(f"Generating plots ({scale_label})")
    print("=" * 70)

    prefix = OUTPUT_DIR / f"{scale_label}_"

    # Merge with residuals for rigidity rank
    cts_ranked = residuals_df.sort_values("residual_magnitude").cell_type.tolist()
    cts_in_align = sorted(alignment_df.cell_type.unique())

    # Filter to types present in both
    cts_ranked = [c for c in cts_ranked if c in cts_in_align]

    # ---- Plot 1: Heatmap of post-Procrustes S(k)/k ----
    fig, ax = plt.subplots(figsize=(8, max(6, len(cts_ranked) * 0.3)))
    pivot = alignment_df.pivot(index="cell_type", columns="k", values="S_post")
    pivot = pivot.loc[[c for c in cts_ranked if c in pivot.index]]
    im = ax.imshow(pivot.values, aspect="auto", cmap="YlOrRd", vmin=0, vmax=1)
    ax.set_yticks(range(len(pivot)))
    ax.set_yticklabels(pivot.index, fontsize=7)
    ax.set_xticks(range(len(K_VALUES)))
    ax.set_xticklabels([f"k={k}" for k in K_VALUES])
    ax.set_xlabel("Number of axes (k)")
    ax.set_ylabel("Cell type (sorted by rigidity, most rigid at bottom)")
    ax.set_title(f"Post-Procrustes Ellipsoid Alignment S(k)/k\n({scale_label})")
    plt.colorbar(im, ax=ax, label="S(k)/k")
    # Add text
    for i in range(len(pivot)):
        for j in range(len(K_VALUES)):
            val = pivot.values[i, j]
            color = "white" if val > 0.6 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                    fontsize=6, color=color)
    plt.tight_layout()
    fig.savefig(prefix / "" if False else f"{prefix}heatmap.png", dpi=150)
    fig.savefig(str(prefix) + "heatmap.png", dpi=150)
    plt.close(fig)
    print(f"  Saved: {prefix}heatmap.png")

    # ---- Plot 2: Pre vs Post bar chart at k=3 ----
    fig, ax = plt.subplots(figsize=(10, max(5, len(cts_ranked) * 0.25)))
    sub = alignment_df[alignment_df.k == 3].set_index("cell_type")
    sub = sub.loc[[c for c in cts_ranked if c in sub.index]]
    y_pos = np.arange(len(sub))
    ax.barh(y_pos - 0.15, sub.S_pre, 0.3, label="Pre-Procrustes", alpha=0.7)
    ax.barh(y_pos + 0.15, sub.S_post, 0.3, label="Post-Procrustes", alpha=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(sub.index, fontsize=7)
    ax.set_xlabel("S(3)/3 alignment score")
    ax.set_title(f"Ellipsoid Alignment: Pre vs Post Procrustes (k=3)\n({scale_label})")
    ax.legend()
    ax.set_xlim(0, 1)
    plt.tight_layout()
    fig.savefig(str(prefix) + "pre_vs_post_k3.png", dpi=150)
    plt.close(fig)
    print(f"  Saved: {prefix}pre_vs_post_k3.png")

    # ---- Plot 3: Null distributions (3 rows × 3 cols) ----
    fig, axes = plt.subplots(3, 3, figsize=(15, 11))
    for col, k in enumerate(K_VALUES):
        # Row 0: Label shuffle PRE-Procrustes
        ax = axes[0, col]
        null_data = perm_label_pre[k]["null"]
        obs_val = perm_label_pre[k]["observed"]
        ax.hist(null_data, bins=50, alpha=0.7, color="seagreen", density=True)
        ax.axvline(obs_val, color="red", linewidth=2, label=f"Obs={obs_val:.3f}")
        ax.set_title(f"Label Shuffle PRE k={k}\np={perm_label_pre[k]['p_value']:.4f}")
        ax.legend(fontsize=7)
        ax.set_xlabel("Mean S(k)/k")

        # Row 1: Label shuffle POST-Procrustes
        ax = axes[1, col]
        null_data = perm_label_post[k]["null"]
        obs_val = perm_label_post[k]["observed"]
        ax.hist(null_data, bins=50, alpha=0.7, color="steelblue", density=True)
        ax.axvline(obs_val, color="red", linewidth=2, label=f"Obs={obs_val:.3f}")
        ax.set_title(f"Label Shuffle POST k={k}\np={perm_label_post[k]['p_value']:.4f}")
        ax.legend(fontsize=7)
        ax.set_xlabel("Mean S(k)/k")

        # Row 2: Random rotation
        ax = axes[2, col]
        null_data = perm_rotation[k]["null"]
        obs_val = perm_rotation[k]["observed"]
        ax.hist(null_data, bins=50, alpha=0.7, color="darkorange", density=True)
        ax.axvline(obs_val, color="red", linewidth=2, label=f"Obs={obs_val:.3f}")
        ax.set_title(f"Random Rotation k={k}\np={perm_rotation[k]['p_value']:.4f}")
        ax.legend(fontsize=7)
        ax.set_xlabel("Mean S(k)/k")

    axes[0, 0].set_ylabel("Label Shuffle\n(Pre-Procrustes)")
    axes[1, 0].set_ylabel("Label Shuffle\n(Post-Procrustes)")
    axes[2, 0].set_ylabel("Random Rotation\nDensity")
    fig.suptitle(f"Permutation Test Null Distributions ({scale_label})", fontsize=14)
    plt.tight_layout()
    fig.savefig(str(prefix) + "null_distributions.png", dpi=150)
    plt.close(fig)
    print(f"  Saved: {prefix}null_distributions.png")

    # ---- Plot 4: Alignment vs rigidity scatter ----
    merged = alignment_df.merge(
        residuals_df[["cell_type", "residual_magnitude"]],
        on="cell_type", how="inner"
    )
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for col, k in enumerate(K_VALUES):
        ax = axes[col]
        sub = merged[merged.k == k]
        ax.scatter(sub.residual_magnitude, sub.S_post, alpha=0.7, s=40)
        for _, row in sub.iterrows():
            label = row.cell_type[:20]
            ax.annotate(label, (row.residual_magnitude, row.S_post),
                       fontsize=5, alpha=0.7)
        rho, p = stats.spearmanr(sub.residual_magnitude, sub.S_post)
        ax.set_xlabel("Procrustes residual (rigidity)")
        ax.set_ylabel(f"S({k})/{k} post-Procrustes")
        ax.set_title(f"k={k}: ρ={rho:+.3f}, p={p:.3f}")
    fig.suptitle(f"Ellipsoid Alignment vs Rigidity ({scale_label})", fontsize=13)
    plt.tight_layout()
    fig.savefig(str(prefix) + "alignment_vs_rigidity.png", dpi=150)
    plt.close(fig)
    print(f"  Saved: {prefix}alignment_vs_rigidity.png")

    # ---- Plot 5: Eigenvalue profiles (6-type only or top 6 from 35) ----
    cts_for_eigenplot = cts_ranked[:6] if len(cts_ranked) > 6 else cts_ranked
    n_plot = len(cts_for_eigenplot)
    fig, axes = plt.subplots(2, min(3, n_plot), figsize=(12, 7))
    if n_plot <= 3:
        axes = axes.reshape(2, -1)
    top_k = 10
    for i, ct in enumerate(cts_for_eigenplot[:6]):
        row_idx = i // 3
        col_idx = i % 3
        if n_plot <= 3 and row_idx > 0:
            break
        ax = axes[row_idx, col_idx] if n_plot > 3 else axes[0, i]
        h_vals = eigen_results[("human", ct)]["var_frac"][:top_k]
        m_vals = eigen_results[("mouse", ct)]["var_frac"][:top_k]
        x = np.arange(1, top_k + 1)
        ax.bar(x - 0.15, h_vals, 0.3, label="Human", alpha=0.7)
        ax.bar(x + 0.15, m_vals, 0.3, label="Mouse", alpha=0.7)
        ax.set_title(ct[:30], fontsize=9)
        ax.set_xlabel("Eigenvalue rank")
        ax.set_ylabel("Fraction of variance")
        if i == 0:
            ax.legend(fontsize=7)
    fig.suptitle(f"Eigenvalue Profiles: Human vs Mouse ({scale_label})", fontsize=13)
    plt.tight_layout()
    fig.savefig(str(prefix) + "eigenvalue_profiles.png", dpi=150)
    plt.close(fig)
    print(f"  Saved: {prefix}eigenvalue_profiles.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def run_analysis(scale, pca, gene_names, cell_types_35, rotation_matrix,
                 residuals_df):
    """Run full analysis for one scale (6type or 35type)."""
    scale_label = f"{scale}"

    # Project cells
    cell_data = project_cells(pca, gene_names, scale=scale)

    # Covariance + eigen
    eigen_results = compute_covariance_eigen(cell_data, scale_label)

    # Alignment scores
    alignment_df = compute_alignment_scores(eigen_results, rotation_matrix, scale_label)

    # Filter residuals to types present
    cts_present = sorted(set(ct for sp, ct in eigen_results.keys() if sp == "human"))
    res_filtered = residuals_df[residuals_df.cell_type.isin(cts_present)].copy()

    # Permutation tests
    perm_label_pre, perm_label_post = perm_test_label_shuffle(
        eigen_results, rotation_matrix)
    perm_rotation = perm_test_random_rotation(eigen_results, rotation_matrix)

    # Correlation with rigidity
    rig_corr = correlate_with_rigidity(alignment_df, res_filtered, scale_label)

    # Eigenvalue conservation
    eigenval_df, ev_rig_rho, ev_rig_p = eigenvalue_conservation(
        eigen_results, res_filtered, scale_label
    )

    # CPC
    cpc_results = common_principal_components(
        eigen_results, pca, gene_names, scale_label
    )

    # Plots
    make_plots(alignment_df, res_filtered, perm_label_pre, perm_label_post,
               perm_rotation, eigen_results, eigenval_df, scale_label)

    return {
        "alignment_df": alignment_df,
        "perm_label_pre": perm_label_pre,
        "perm_label_post": perm_label_post,
        "perm_rotation": perm_rotation,
        "rig_corr": rig_corr,
        "eigenval_df": eigenval_df,
        "eigenval_rig": (ev_rig_rho, ev_rig_p),
        "cpc_results": cpc_results,
        "eigen_results": eigen_results,
    }


def save_outputs(results_6, results_35):
    """Save all outputs to disk."""
    print("\n" + "=" * 70)
    print("Saving outputs")
    print("=" * 70)

    for label, res in [("6type", results_6), ("35type", results_35)]:
        if res is None:
            continue

        # Alignment scores CSV
        res["alignment_df"].to_csv(
            OUTPUT_DIR / f"{label}_alignment_scores.csv", index=False
        )

        # Eigenvalue conservation CSV
        res["eigenval_df"].to_csv(
            OUTPUT_DIR / f"{label}_eigenvalue_conservation.csv", index=False
        )

        # Rigidity correlation CSV
        res["rig_corr"].to_csv(
            OUTPUT_DIR / f"{label}_rigidity_correlation.csv", index=False
        )

    # Permutation results JSON
    perm_json = {}
    for label, res in [("6type", results_6), ("35type", results_35)]:
        if res is None:
            continue
        perm_json[label] = {}
        for test_name, test_data in [
            ("label_shuffle_pre", res["perm_label_pre"]),
            ("label_shuffle_post", res["perm_label_post"]),
            ("random_rotation", res["perm_rotation"]),
        ]:
            perm_json[label][test_name] = {}
            for k in K_VALUES:
                perm_json[label][test_name][f"k={k}"] = {
                    "observed": float(test_data[k]["observed"]),
                    "null_mean": float(test_data[k]["null"].mean()),
                    "null_std": float(test_data[k]["null"].std()),
                    "p_value": float(test_data[k]["p_value"]),
                }

    with open(OUTPUT_DIR / "permutation_results.json", "w") as f:
        json.dump(perm_json, f, indent=2)

    # CPC genes JSON
    cpc_json = {}
    for label, res in [("6type", results_6), ("35type", results_35)]:
        if res is None:
            continue
        cpc_json[label] = {}
        for ct, data in res["cpc_results"].items():
            cpc_json[label][ct] = {
                "cpc1_var_frac_human": data["cpc1_var_frac_human"],
                "cpc1_var_frac_mouse": data["cpc1_var_frac_mouse"],
                "top5_genes": [
                    {"gene": g, "loading": float(l)} for g, l in data["top5_genes"]
                ],
            }

    with open(OUTPUT_DIR / "cpc_genes.json", "w") as f:
        json.dump(cpc_json, f, indent=2)

    # Summary stats JSON
    summary = {}
    for label, res in [("6type", results_6), ("35type", results_35)]:
        if res is None:
            continue
        adf = res["alignment_df"]
        summary[label] = {
            "n_types": len(adf.cell_type.unique()),
            "mean_alignment": {
                f"k={k}": {
                    "pre": float(adf[adf.k == k].S_pre.mean()),
                    "post": float(adf[adf.k == k].S_post.mean()),
                    "improvement": float(adf[adf.k == k].improvement.mean()),
                }
                for k in K_VALUES
            },
            "eigenval_vs_rigidity": {
                "rho": float(res["eigenval_rig"][0]) if not (isinstance(res["eigenval_rig"][0], float) and np.isnan(res["eigenval_rig"][0])) else None,
                "p": float(res["eigenval_rig"][1]) if not (isinstance(res["eigenval_rig"][1], float) and np.isnan(res["eigenval_rig"][1])) else None,
            },
        }

    with open(OUTPUT_DIR / "summary_stats.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("  All outputs saved to:", OUTPUT_DIR)


def main():
    t_start = time.time()

    # Reconstruct PCA
    pca, cell_types_35, gene_names = reconstruct_pca()

    # Load Procrustes rotation
    with open(OUTPUT_SCALED / "procrustes_results_35.json") as f:
        proc = json.load(f)
    R = np.array(proc["procrustes"]["rotation_matrix"])
    print(f"  Rotation matrix: {R.shape}, det={np.linalg.det(R):.4f}")

    # Load residuals
    residuals_df = pd.read_csv(OUTPUT_SCALED / "residuals_ranked.csv")

    # Run 35-type analysis
    print("\n\n" + "#" * 70)
    print("# 35-TYPE ANALYSIS")
    print("#" * 70)
    results_35 = run_analysis("35type", pca, gene_names, cell_types_35, R, residuals_df)

    # Run 6-type analysis
    print("\n\n" + "#" * 70)
    print("# 6-TYPE ANALYSIS")
    print("#" * 70)
    results_6 = run_analysis("6type", pca, gene_names, cell_types_35, R, residuals_df)

    # Save
    save_outputs(results_6, results_35)

    elapsed = time.time() - t_start
    print(f"\nTotal runtime: {elapsed/60:.1f} minutes")

    # ---- FINAL SUMMARY ----
    print("\n\n" + "=" * 70)
    print("FINAL SUMMARY — T3-B Covariance Ellipsoid Alignment")
    print("=" * 70)

    for label, res in [("35-TYPE", results_35), ("6-TYPE", results_6)]:
        print(f"\n--- {label} ---")
        adf = res["alignment_df"]
        n_types = len(adf.cell_type.unique())
        print(f"  Cell types: {n_types}")

        for k in K_VALUES:
            sub = adf[adf.k == k]
            print(f"\n  k={k}:")
            print(f"    Mean alignment: pre={sub.S_pre.mean():.3f}, "
                  f"post={sub.S_post.mean():.3f}, "
                  f"Δ={sub.improvement.mean():+.3f}")
            print(f"    Label shuffle PRE p={res['perm_label_pre'][k]['p_value']:.4f}")
            print(f"    Label shuffle POST p={res['perm_label_post'][k]['p_value']:.4f}")
            print(f"    Random rotation p={res['perm_rotation'][k]['p_value']:.4f}")

        ev_rho, ev_p = res["eigenval_rig"]
        if np.isnan(ev_rho):
            print(f"\n  Eigenvalue conservation vs rigidity: N/A (too few types)")
        else:
            print(f"\n  Eigenvalue conservation vs rigidity: "
                  f"ρ={ev_rho:+.3f}, p={ev_p:.4f}")

        # Rigidity correlations
        print(f"\n  Alignment vs rigidity (FDR-corrected):")
        for _, row in res["rig_corr"].iterrows():
            sig = " *" if row.q_fdr < 0.05 else ""
            print(f"    k={int(row.k)} {row.metric}: ρ={row.rho:+.3f}, "
                  f"p={row.p:.4f}, q={row.q_fdr:.4f}{sig}")


if __name__ == "__main__":
    main()
