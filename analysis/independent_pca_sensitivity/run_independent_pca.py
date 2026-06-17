#!/usr/bin/env python3
"""
Independent-PCA Sensitivity Analysis for CellWarp

Reviewer concern: Joint PCA on 70 centroids (35 human + 35 mouse) optimizes
the subspace for shared variance, potentially inflating apparent geometric
coherence and biasing per-type residuals.

This script addresses the concern by:
1. Computing PCA SEPARATELY on 35 human and 35 mouse centroids
2. Aligning the two PCA subspaces via Procrustes on loading matrices
   (does NOT use cell-type pairing — only aligns PCA axes in gene space)
3. Projecting both species into the aligned subspace
4. Running Procrustes superimposition + permutation test (10,000 iterations)
5. Computing per-type residuals and comparing rigidity rankings to joint PCA

Math — Subspace Alignment via Procrustes on Loadings
-----------------------------------------------------
Let W_H (k × G) and W_M (k × G) be the PCA loading matrices for human and
mouse respectively, where rows are orthonormal PCA axes in gene space.

We seek rotation Q (k × k) such that Q @ W_M ≈ W_H, i.e., Q maps mouse
PCA axes to their best-matching human PCA axes.

    C = W_H @ W_M.T    (k × k cross-product of loading matrices)
    C = U Σ V^T         (SVD)
    Q = U V^T           (optimal rotation, det = +1 enforced)

The aligned mouse scores become: M_aligned = M_pca @ Q.T

Crucially, this alignment uses ONLY the geometry of the PCA axes in gene
space — it does not use the cell-type correspondence. The pairing information
enters only at the Procrustes superimposition step.

Biology: If the geometric coherence in the joint-PCA analysis is an artifact
of shared embedding, the independent-PCA analysis should show substantially
weaker coherence (higher obs/null ratio, higher p-value). If the signal is
genuine, the independent-PCA results should closely reproduce the joint-PCA
findings.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from numpy.linalg import svd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "output" / "phase2" / "scaled_35types"
OUTPUT_DIR = PROJECT_ROOT / "analysis" / "independent_pca_sensitivity"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

RANDOM_SEED = 42
N_PERMUTATIONS = 10_000
PCA_VARIANCE_THRESHOLD = 0.95
JOINT_PCA_N_COMPONENTS = 33  # from primary analysis


# ---------------------------------------------------------------------------
# Procrustes helpers (from src/procrustes.py, adapted for silent operation)
# ---------------------------------------------------------------------------


def _procrustes_distance(X: np.ndarray, Y: np.ndarray) -> float:
    """Compute Procrustes distance (no reflection, silent)."""
    n, k = X.shape
    X_c = X - X.mean(axis=0)
    Y_c = Y - Y.mean(axis=0)

    M = X_c.T @ Y_c
    U, sigma, Vt = svd(M)
    V = Vt.T

    d = np.linalg.det(V @ U.T)
    D_diag = np.ones(k)
    D_diag[-1] = np.sign(d)

    ss_Y = np.sum(Y_c ** 2)
    trace_sigma_D = np.sum(sigma * D_diag)
    s = trace_sigma_D / ss_Y

    Y_aligned = s * (Y_c @ (V * D_diag) @ U.T)
    return np.sqrt(np.sum((X_c - Y_aligned) ** 2))


def procrustes_align(X: np.ndarray, Y: np.ndarray):
    """Full Procrustes alignment returning all components."""
    n, k = X.shape
    mu_X = X.mean(axis=0)
    mu_Y = Y.mean(axis=0)
    X_c = X - mu_X
    Y_c = Y - mu_Y

    M = X_c.T @ Y_c
    U, sigma, Vt = svd(M)
    V = Vt.T

    d = np.linalg.det(V @ U.T)
    D = np.eye(k)
    D[-1, -1] = np.sign(d)

    R = V @ D @ U.T
    ss_Y = np.sum(Y_c ** 2)
    trace_sigma_D = np.sum(sigma * np.diag(D))
    s = trace_sigma_D / ss_Y

    Y_aligned = s * (Y_c @ R)
    d_sq = np.sum((X_c - Y_aligned) ** 2)

    return {
        "rotation": R,
        "scaling": s,
        "distance": np.sqrt(d_sq),
        "distance_squared": d_sq,
        "aligned_target": Y_aligned,
        "centered_reference": X_c,
        "mu_X": mu_X,
        "mu_Y": mu_Y,
    }


def permutation_test(
    X: np.ndarray,
    Y: np.ndarray,
    n_permutations: int = N_PERMUTATIONS,
    seed: int = RANDOM_SEED,
) -> tuple[float, np.ndarray, float]:
    """Permutation test returning (p_value, null_distribution, observed)."""
    rng = np.random.RandomState(seed)
    observed = _procrustes_distance(X, Y)

    null_distances = np.zeros(n_permutations)
    for i in range(n_permutations):
        perm = rng.permutation(X.shape[0])
        null_distances[i] = _procrustes_distance(X, Y[perm])

    n_leq = int(np.sum(null_distances <= observed))
    p_value = (n_leq + 1) / (n_permutations + 1)

    return p_value, null_distances, observed


# ---------------------------------------------------------------------------
# Subspace alignment
# ---------------------------------------------------------------------------


def align_subspaces_procrustes(
    W_H: np.ndarray, W_M: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """
    Align mouse PCA subspace to human PCA subspace via Procrustes on loadings.

    W_H, W_M: (k × G) loading matrices with orthonormal rows.

    Returns:
        Q: (k × k) rotation matrix mapping mouse PCA axes → human PCA axes
        singular_values: (k,) alignment quality per axis (1.0 = perfect match)
    """
    C = W_H @ W_M.T  # (k × k)
    U, sigma, Vt = svd(C)
    V = Vt.T

    # Enforce proper rotation (det = +1)
    d = np.linalg.det(V @ U.T)
    D = np.eye(len(sigma))
    D[-1, -1] = np.sign(d)

    Q = V @ D @ U.T  # not U @ V.T — we want Q @ W_M ≈ W_H, solved as min ‖W_H - Q W_M‖
    # Actually: we want Q s.t. Q @ W_M ≈ W_H
    # C = W_H @ W_M.T, SVD of C = U Σ V^T
    # Optimal Q = U D V^T (maps mouse loadings to human loadings)
    Q = U @ D @ Vt

    return Q, sigma


def align_subspaces_cca(
    H_pca: np.ndarray, M_pca: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Align subspaces via Canonical Correlation Analysis.

    Unlike Procrustes-on-loadings, CCA uses the cell-type pairing to find
    maximally correlated subspace axes. Included for comparison.

    Returns:
        H_cca: (n × k) human canonical variates
        M_cca: (n × k) mouse canonical variates
        correlations: (k,) canonical correlations
    """
    from sklearn.cross_decomposition import CCA as SkCCA

    n, k = H_pca.shape
    # CCA with k components
    n_components = min(k, n - 1)
    cca = SkCCA(n_components=n_components, max_iter=1000, tol=1e-10)
    H_cca, M_cca = cca.fit_transform(H_pca, M_pca)

    # Compute canonical correlations
    correlations = np.array([
        np.corrcoef(H_cca[:, i], M_cca[:, i])[0, 1]
        for i in range(n_components)
    ])

    return H_cca, M_cca, correlations


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------


def main():
    print("=" * 72)
    print("INDEPENDENT-PCA SENSITIVITY ANALYSIS")
    print("Addresses reviewer concern: joint PCA may inflate geometric coherence")
    print("=" * 72)

    # ------------------------------------------------------------------
    # 1. Load centroid data
    # ------------------------------------------------------------------
    print("\n[1] Loading 35-type centroids from gene space...")

    human_centroids = pd.read_csv(DATA_DIR / "centroids_human_35.csv", index_col=0)
    mouse_centroids = pd.read_csv(DATA_DIR / "centroids_mouse_35.csv", index_col=0)

    cell_types = sorted(human_centroids.index.tolist())
    assert sorted(mouse_centroids.index.tolist()) == cell_types
    n_types = len(cell_types)
    n_genes = human_centroids.shape[1]
    print(f"  {n_types} cell types × {n_genes:,} genes")

    H_gene = human_centroids.loc[cell_types].values  # (35, G)
    M_gene = mouse_centroids.loc[cell_types].values  # (35, G)

    # ------------------------------------------------------------------
    # 2. Load joint-PCA results for comparison
    # ------------------------------------------------------------------
    print("\n[2] Loading joint-PCA results for comparison...")

    joint_residuals = pd.read_csv(DATA_DIR / "residuals_ranked.csv")
    joint_ranking = {
        row["cell_type"]: row["rank"]
        for _, row in joint_residuals.iterrows()
    }
    print(f"  Joint PCA: {len(joint_ranking)} cell types ranked")

    joint_npz = np.load(DATA_DIR / "pca_centroids_35.npz", allow_pickle=True)
    joint_n_components = joint_npz["human"].shape[1]
    print(f"  Joint PCA used {joint_n_components} components")

    with open(DATA_DIR / "procrustes_results_35.json") as f:
        joint_results = json.load(f)
    joint_distance = joint_results["procrustes"]["distance"]
    joint_null_median = joint_results["permutation_test"]["null_distribution_summary"]["median"]
    joint_obs_null = joint_distance / joint_null_median
    joint_p = joint_results["permutation_test"]["p_value"]
    print(f"  Joint PCA: distance={joint_distance:.3f}, null_median={joint_null_median:.3f}")
    print(f"  Joint PCA: obs/null={joint_obs_null:.4f}, p={joint_p}")

    # ------------------------------------------------------------------
    # 3. Compute PCA SEPARATELY on human and mouse
    # ------------------------------------------------------------------
    print("\n[3] Computing PCA separately on each species...")

    # Determine number of components: use same as joint (33), or 95% variance
    # With 35 samples, max possible = 34 components each
    max_components = min(n_types - 1, n_genes)  # 34

    # Fit full PCA first to check variance
    pca_human_full = PCA(n_components=max_components, svd_solver="full",
                         random_state=RANDOM_SEED)
    pca_mouse_full = PCA(n_components=max_components, svd_solver="full",
                         random_state=RANDOM_SEED)

    H_pca_full = pca_human_full.fit_transform(H_gene)
    M_pca_full = pca_mouse_full.fit_transform(M_gene)

    # Cumulative variance with joint_n_components (33)
    cumvar_h = np.cumsum(pca_human_full.explained_variance_ratio_)
    cumvar_m = np.cumsum(pca_mouse_full.explained_variance_ratio_)

    print(f"\n  Human PCA ({max_components} max components):")
    print(f"    Variance at k={joint_n_components}: {cumvar_h[joint_n_components-1]*100:.2f}%")
    print(f"    Variance at k={max_components}: {cumvar_h[-1]*100:.2f}%")

    print(f"\n  Mouse PCA ({max_components} max components):")
    print(f"    Variance at k={joint_n_components}: {cumvar_m[joint_n_components-1]*100:.2f}%")
    print(f"    Variance at k={max_components}: {cumvar_m[-1]*100:.2f}%")

    # Find k for 95% variance threshold in each
    k_95_h = int(np.searchsorted(cumvar_h, PCA_VARIANCE_THRESHOLD) + 1)
    k_95_m = int(np.searchsorted(cumvar_m, PCA_VARIANCE_THRESHOLD) + 1)
    print(f"\n  Components for 95% variance: human={k_95_h}, mouse={k_95_m}")

    # Use joint_n_components (33) for direct comparability
    k = joint_n_components
    print(f"  Using k={k} components (matching joint PCA) for primary analysis")

    H_pca = H_pca_full[:, :k]  # (35, 33)
    M_pca = M_pca_full[:, :k]  # (35, 33)

    W_H = pca_human_full.components_[:k]  # (33, G)
    W_M = pca_mouse_full.components_[:k]  # (33, G)

    print(f"\n  Per-component variance explained:")
    print(f"  {'PC':<5} {'Human':>10} {'Mouse':>10} {'Joint':>10}")
    print(f"  {'-'*40}")
    joint_var = joint_npz.get("explained_variance_ratio", None)
    for i in range(min(10, k)):
        h_var = pca_human_full.explained_variance_ratio_[i] * 100
        m_var = pca_mouse_full.explained_variance_ratio_[i] * 100
        j_str = ""
        if joint_var is not None and i < len(joint_var):
            j_str = f"{joint_var[i]*100:.2f}%"
        print(f"  PC{i+1:<3} {h_var:>9.2f}% {m_var:>9.2f}% {j_str:>10}")
    print(f"  ...")

    # ------------------------------------------------------------------
    # 4. Align subspaces via Procrustes on loading matrices
    # ------------------------------------------------------------------
    print("\n[4] Aligning PCA subspaces via Procrustes on loading matrices...")

    Q, alignment_sigmas = align_subspaces_procrustes(W_H, W_M)

    print(f"  Rotation matrix Q: {Q.shape}, det = {np.linalg.det(Q):+.6f}")
    print(f"  Alignment singular values (top 10):")
    for i in range(min(10, len(alignment_sigmas))):
        print(f"    Axis {i+1}: σ = {alignment_sigmas[i]:.6f}")
    print(f"  Mean alignment quality: {np.mean(alignment_sigmas):.6f}")
    print(f"  Min alignment quality: {np.min(alignment_sigmas):.6f}")

    # Rotate mouse PCA scores to human PCA coordinate system
    M_aligned = M_pca @ Q.T  # (35, 33)

    print(f"\n  Human PCA scores shape: {H_pca.shape}")
    print(f"  Mouse aligned scores shape: {M_aligned.shape}")

    # ------------------------------------------------------------------
    # 5. Procrustes superimposition in independent-PCA space
    # ------------------------------------------------------------------
    print("\n[5] Running Procrustes superimposition in independent-PCA space...")

    result = procrustes_align(H_pca, M_aligned)
    print(f"  Procrustes distance: {result['distance']:.6f}")
    print(f"  Procrustes distance²: {result['distance_squared']:.6f}")
    print(f"  Scaling factor: {result['scaling']:.6f}")
    print(f"  Rotation determinant: {np.linalg.det(result['rotation']):+.6f}")

    # ------------------------------------------------------------------
    # 6. Permutation test (10,000 iterations)
    # ------------------------------------------------------------------
    print(f"\n[6] Running permutation test ({N_PERMUTATIONS:,} iterations)...")
    t0 = time.time()
    p_value, null_dist, observed = permutation_test(H_pca, M_aligned)
    elapsed = time.time() - t0
    print(f"  Completed in {elapsed:.1f}s")

    null_mean = np.mean(null_dist)
    null_median = np.median(null_dist)
    null_std = np.std(null_dist)
    obs_null_ratio = observed / null_median
    n_leq = int(np.sum(null_dist <= observed))

    print(f"\n  Observed distance: {observed:.6f}")
    print(f"  Null distribution: mean={null_mean:.6f}, median={null_median:.6f}, std={null_std:.6f}")
    print(f"  Null range: [{np.min(null_dist):.6f}, {np.max(null_dist):.6f}]")
    print(f"  Permuted distances ≤ observed: {n_leq} / {N_PERMUTATIONS}")
    print(f"  p-value: {p_value:.6f}")
    print(f"  obs/null ratio: {obs_null_ratio:.4f}")
    print(f"  Significant at α=0.01: {'YES' if p_value < 0.01 else 'NO'}")

    # ------------------------------------------------------------------
    # 7. Per-type residuals and rigidity ranking
    # ------------------------------------------------------------------
    print("\n[7] Computing per-type residuals in independent-PCA space...")

    total_ssr = result["distance_squared"]
    indep_residuals = []
    for i, ct in enumerate(cell_types):
        r = result["aligned_target"][i] - result["centered_reference"][i]
        mag = np.linalg.norm(r)
        pct = (mag ** 2 / total_ssr * 100) if total_ssr > 0 else 0.0
        indep_residuals.append({
            "cell_type": ct,
            "residual_magnitude": mag,
            "pct_of_ssr": pct,
            "residual_vector": r.tolist(),
        })

    # Sort by residual magnitude (descending = most divergent first)
    indep_residuals.sort(key=lambda x: x["residual_magnitude"], reverse=True)
    for rank, entry in enumerate(indep_residuals, 1):
        entry["rank"] = rank

    print(f"\n  {'Rank':<5} {'Cell Type':<45} {'Magnitude':>10} {'% SSR':>8}")
    print(f"  {'-'*72}")
    for entry in indep_residuals:
        print(
            f"  {entry['rank']:<5} {entry['cell_type']:<45} "
            f"{entry['residual_magnitude']:>10.4f} {entry['pct_of_ssr']:>7.2f}%"
        )

    # ------------------------------------------------------------------
    # 8. Spearman ρ between joint and independent rigidity rankings
    # ------------------------------------------------------------------
    print("\n[8] Comparing rigidity rankings: joint PCA vs independent PCA...")

    indep_ranking = {e["cell_type"]: e["rank"] for e in indep_residuals}

    # Build paired arrays (same cell type order)
    joint_ranks = []
    indep_ranks = []
    for ct in cell_types:
        joint_ranks.append(joint_ranking[ct])
        indep_ranks.append(indep_ranking[ct])

    rho, rho_p = stats.spearmanr(joint_ranks, indep_ranks)
    print(f"  Spearman ρ = {rho:.4f} (p = {rho_p:.6f})")

    # Check for dramatic rank changes (>10 positions)
    print(f"\n  Rank comparison (|Δ| > 10 flagged):")
    print(f"  {'Cell Type':<45} {'Joint':>6} {'Indep':>6} {'Δ':>6} {'Flag':>6}")
    print(f"  {'-'*72}")
    dramatic_changes = []
    for ct in cell_types:
        j_r = joint_ranking[ct]
        i_r = indep_ranking[ct]
        delta = i_r - j_r
        flag = "***" if abs(delta) > 10 else ""
        if abs(delta) > 10:
            dramatic_changes.append((ct, j_r, i_r, delta))
        print(f"  {ct:<45} {j_r:>6} {i_r:>6} {delta:>+6} {flag:>6}")

    if dramatic_changes:
        print(f"\n  ⚠ {len(dramatic_changes)} cell type(s) changed >10 rank positions:")
        for ct, j_r, i_r, delta in dramatic_changes:
            print(f"    {ct}: joint rank {j_r} → indep rank {i_r} (Δ={delta:+d})")
    else:
        print(f"\n  No cell types changed >10 rank positions.")

    # ------------------------------------------------------------------
    # 9. CCA comparison (secondary analysis)
    # ------------------------------------------------------------------
    print("\n[9] CCA comparison (secondary analysis — uses pairing info)...")

    try:
        H_cca, M_cca, can_corrs = align_subspaces_cca(H_pca, M_pca)
        print(f"  CCA canonical correlations (top 10):")
        for i in range(min(10, len(can_corrs))):
            print(f"    CC{i+1}: {can_corrs[i]:.6f}")
        print(f"  Mean canonical correlation: {np.mean(can_corrs):.6f}")

        # Procrustes in CCA space
        p_cca, null_cca, obs_cca = permutation_test(H_cca, M_cca)
        null_mean_cca = np.mean(null_cca)
        obs_null_cca = obs_cca / np.median(null_cca)
        print(f"\n  CCA-space Procrustes:")
        print(f"    Observed distance: {obs_cca:.6f}")
        print(f"    Null mean: {null_mean_cca:.6f}")
        print(f"    obs/null ratio: {obs_null_cca:.4f}")
        print(f"    p-value: {p_cca:.6f}")

        cca_results = {
            "observed_distance": float(obs_cca),
            "null_mean": float(null_mean_cca),
            "obs_null_ratio": float(obs_null_cca),
            "p_value": float(p_cca),
            "canonical_correlations": can_corrs.tolist(),
        }
    except Exception as e:
        print(f"  CCA failed: {e}")
        cca_results = {"error": str(e)}

    # ------------------------------------------------------------------
    # 10. Save all outputs
    # ------------------------------------------------------------------
    print("\n[10] Saving outputs...")

    # Save residuals ranking
    residuals_df = pd.DataFrame([
        {
            "rank": e["rank"],
            "cell_type": e["cell_type"],
            "residual_magnitude": e["residual_magnitude"],
            "pct_of_ssr": e["pct_of_ssr"],
            "joint_rank": joint_ranking[e["cell_type"]],
            "rank_delta": e["rank"] - joint_ranking[e["cell_type"]],
        }
        for e in indep_residuals
    ])
    residuals_df.to_csv(OUTPUT_DIR / "residuals_ranked_independent_pca.csv", index=False)
    print(f"  Saved: residuals_ranked_independent_pca.csv")

    # Save null distribution
    np.save(OUTPUT_DIR / "null_distribution_independent_pca.npy", null_dist)
    print(f"  Saved: null_distribution_independent_pca.npy")

    # Save PCA centroids in aligned space
    np.savez(
        OUTPUT_DIR / "pca_centroids_independent.npz",
        human=H_pca,
        mouse_aligned=M_aligned,
        cell_types=np.array(cell_types),
        explained_variance_ratio_human=pca_human_full.explained_variance_ratio_[:k],
        explained_variance_ratio_mouse=pca_mouse_full.explained_variance_ratio_[:k],
        alignment_singular_values=alignment_sigmas,
        subspace_rotation_Q=Q,
    )
    print(f"  Saved: pca_centroids_independent.npz")

    # Save full results JSON
    full_results = {
        "method": "Independent PCA with Procrustes subspace alignment",
        "description": (
            "PCA computed separately on 35 human and 35 mouse centroids. "
            "Subspaces aligned via Procrustes on loading matrices (does not "
            "use cell-type pairing). Procrustes superimposition then tested "
            "in the aligned subspace."
        ),
        "config": {
            "n_cell_types": n_types,
            "n_genes": n_genes,
            "n_pca_components": k,
            "n_permutations": N_PERMUTATIONS,
            "random_seed": RANDOM_SEED,
            "pca_variance_threshold": PCA_VARIANCE_THRESHOLD,
        },
        "separate_pca": {
            "human": {
                "n_components": k,
                "cumulative_variance": float(cumvar_h[k - 1]),
                "variance_per_component": pca_human_full.explained_variance_ratio_[:k].tolist(),
            },
            "mouse": {
                "n_components": k,
                "cumulative_variance": float(cumvar_m[k - 1]),
                "variance_per_component": pca_mouse_full.explained_variance_ratio_[:k].tolist(),
            },
            "k_for_95pct_human": int(k_95_h),
            "k_for_95pct_mouse": int(k_95_m),
        },
        "subspace_alignment": {
            "method": "Procrustes on loading matrices",
            "rotation_Q_det": float(np.linalg.det(Q)),
            "singular_values": alignment_sigmas.tolist(),
            "mean_alignment_quality": float(np.mean(alignment_sigmas)),
            "min_alignment_quality": float(np.min(alignment_sigmas)),
        },
        "procrustes": {
            "distance": float(result["distance"]),
            "distance_squared": float(result["distance_squared"]),
            "scaling": float(result["scaling"]),
            "rotation_det": float(np.linalg.det(result["rotation"])),
        },
        "permutation_test": {
            "observed_distance": float(observed),
            "p_value": float(p_value),
            "n_leq_observed": n_leq,
            "n_permutations": N_PERMUTATIONS,
            "null_distribution_summary": {
                "mean": float(null_mean),
                "median": float(null_median),
                "std": float(null_std),
                "min": float(np.min(null_dist)),
                "max": float(np.max(null_dist)),
                "percentile_2_5": float(np.percentile(null_dist, 2.5)),
                "percentile_97_5": float(np.percentile(null_dist, 97.5)),
            },
            "obs_null_ratio": float(obs_null_ratio),
            "significant_at_001": bool(p_value < 0.01),
        },
        "comparison_to_joint_pca": {
            "joint_pca": {
                "distance": float(joint_distance),
                "null_mean": float(joint_null_mean),
                "obs_null_ratio": float(joint_obs_null),
                "p_value": float(joint_p),
                "n_components": joint_n_components,
            },
            "independent_pca": {
                "distance": float(observed),
                "null_mean": float(null_mean),
                "obs_null_ratio": float(obs_null_ratio),
                "p_value": float(p_value),
                "n_components": k,
            },
            "spearman_rho": float(rho),
            "spearman_p": float(rho_p),
            "n_dramatic_rank_changes": len(dramatic_changes),
            "dramatic_rank_changes": [
                {"cell_type": ct, "joint_rank": j, "indep_rank": i, "delta": d}
                for ct, j, i, d in dramatic_changes
            ],
        },
        "cca_secondary": cca_results,
        "residuals": {
            e["cell_type"]: {
                "rank": e["rank"],
                "magnitude": e["residual_magnitude"],
                "pct_of_ssr": e["pct_of_ssr"],
                "joint_rank": joint_ranking[e["cell_type"]],
                "rank_delta": e["rank"] - joint_ranking[e["cell_type"]],
            }
            for e in indep_residuals
        },
    }

    results_path = OUTPUT_DIR / "independent_pca_results.json"
    with open(results_path, "w") as f:
        json.dump(full_results, f, indent=2)
    print(f"  Saved: independent_pca_results.json")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)
    print("SUMMARY: Independent PCA vs Joint PCA")
    print("=" * 72)
    print(f"\n  {'Metric':<30} {'Joint PCA':>15} {'Independent PCA':>15}")
    print(f"  {'-'*62}")
    print(f"  {'Procrustes distance':<30} {joint_distance:>15.3f} {observed:>15.3f}")
    print(f"  {'Null mean':<30} {joint_null_mean:>15.3f} {null_mean:>15.3f}")
    print(f"  {'obs/null ratio':<30} {joint_obs_null:>15.4f} {obs_null_ratio:>15.4f}")
    print(f"  {'p-value':<30} {joint_p:>15.4f} {p_value:>15.4f}")
    print(f"  {'Significant (α=0.01)':<30} {'YES':>15} {'YES' if p_value < 0.01 else 'NO':>15}")
    print(f"  {'PCA components':<30} {joint_n_components:>15} {k:>15}")
    print(f"\n  Spearman ρ (rigidity ranking): {rho:.4f} (p = {rho_p:.6f})")
    print(f"  Rank changes > 10 positions: {len(dramatic_changes)}")

    if p_value < 0.01 and rho > 0.7:
        print(f"\n  CONCLUSION: Independent PCA reproduces the joint-PCA finding.")
        print(f"  The geometric coherence signal is NOT an artifact of joint embedding.")
    elif p_value < 0.01:
        print(f"\n  CONCLUSION: Independent PCA preserves significance but rankings diverge.")
        print(f"  Global coherence is robust; per-type residuals are partially embedding-dependent.")
    else:
        print(f"\n  CONCLUSION: Independent PCA does NOT reproduce the joint-PCA finding.")
        print(f"  The reviewer concern may be warranted — further investigation needed.")

    print()
    return full_results


if __name__ == "__main__":
    results = main()
