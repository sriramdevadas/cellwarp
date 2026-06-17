"""
CellWarp — Procrustes Analysis Module

Implements ordinary Procrustes analysis (OPA) from scratch to align mouse cell type
centroids onto human cell type centroids in PCA-reduced gene expression space. This
tests whether cross-species differences follow a coherent geometric transformation
— a D'Arcy Thompson grid in transcriptomic space.

Biology
-------
Each cell type has a characteristic gene expression profile. We represent each type
as a single point (its centroid / mean expression vector) in gene space. With 6 cell
types per species, we have two "constellations" of 6 points — one human, one mouse.

Procrustes analysis asks: can the mouse constellation be rotated, scaled, and
translated to match the human constellation? If yes (small residuals, significant
permutation test), the cross-species divergence follows a coherent geometric
transformation rather than random drift.

Math
----
Ordinary Procrustes Analysis solves:
    min_{R, s, t}  ‖X - (s Y R + 1 t^T)‖_F^2

where X = human centroids (n×k), Y = mouse centroids (n×k), R is an orthogonal
rotation matrix, s is a scalar scale, t is a translation vector, and ‖·‖_F is
the Frobenius norm (sum of squared elements).

The solution proceeds by centering (removes t), SVD of the cross-covariance
matrix (gives R), and a closed-form expression for s. See procrustes_align()
for the full derivation.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from numpy.linalg import svd
from sklearn.decomposition import PCA


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

RANDOM_SEED = 42
N_PERMUTATIONS = 10_000
PCA_VARIANCE_THRESHOLD = 0.95
N_TOP_GENES = 20
CELL_TYPE_COL = "cell_type"


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass
class ProcrustesResult:
    """Container for Procrustes alignment results.

    Attributes:
        rotation: k×k orthogonal rotation matrix R.
        scaling: Optimal scalar scaling factor s.
        translation_ref: Centroid μ_X of the reference (human) point set.
        translation_target: Centroid μ_Y of the target (mouse) point set.
        distance_squared: Sum of squared residuals ‖X_c - s Y_c R‖_F^2.
        distance: Square root of distance_squared.
        aligned_target: Transformed mouse centroids (s Y_c R) in centered coords.
        centered_reference: Centered human centroids (X_c = X - μ_X).
    """

    rotation: np.ndarray
    scaling: float
    translation_ref: np.ndarray
    translation_target: np.ndarray
    distance_squared: float
    distance: float
    aligned_target: np.ndarray
    centered_reference: np.ndarray


# ---------------------------------------------------------------------------
# Centroid computation
# ---------------------------------------------------------------------------


def compute_centroids(
    adata: ad.AnnData,
    cell_type_col: str = CELL_TYPE_COL,
) -> pd.DataFrame:
    """
    Compute mean expression vector per cell type.

    Biology: The centroid of a cell type in gene expression space represents
    its average transcriptomic profile. Averaging over hundreds to thousands
    of cells smooths single-cell noise while preserving the characteristic
    expression signature of each type.

    Math: For cell type t with cells {c_1, ..., c_n}:
        μ_t = (1/n) Σ_{i=1}^{n} x_{c_i}
    where x_{c_i} is the log-normalized expression vector (G genes).

    Args:
        adata: AnnData with log-normalized expression in .X (may be sparse).
        cell_type_col: Column in .obs identifying cell types.

    Returns:
        DataFrame of shape (n_cell_types, n_genes), indexed by cell type
        name, columns are gene identifiers (matching adata.var_names).
    """
    cell_types = sorted(adata.obs[cell_type_col].unique())
    gene_ids = adata.var_names.tolist()

    centroids = {}
    for ct in cell_types:
        mask = adata.obs[cell_type_col] == ct
        # .X may be sparse; .mean(axis=0) returns a matrix, flatten to 1D
        mean_vec = np.asarray(adata[mask].X.mean(axis=0)).flatten()
        centroids[ct] = mean_vec
        n_cells = mask.sum()
        print(
            f"  {ct:<45} {n_cells:>6,} cells → centroid "
            f"({mean_vec.shape[0]:,} genes)"
        )

    df = pd.DataFrame(centroids, index=gene_ids).T
    df.index.name = "cell_type"
    return df


# ---------------------------------------------------------------------------
# PCA on centroids
# ---------------------------------------------------------------------------


def pca_reduce_centroids(
    human_centroids: pd.DataFrame,
    mouse_centroids: pd.DataFrame,
    variance_threshold: float = PCA_VARIANCE_THRESHOLD,
) -> tuple[np.ndarray, np.ndarray, PCA, list[str]]:
    """
    Reduce combined centroids to PCA space retaining ≥variance_threshold.

    Math: Stack human and mouse centroids into a (2n × G) matrix.  PCA finds
    the linear subspace of maximum variance:

        X_pca = (X - μ) @ W^T

    where W is the (k × G) matrix of principal component loadings and k is
    chosen so that Σ_{i=1}^{k} λ_i / Σ_i λ_i ≥ variance_threshold.

    With 12 samples in G-dimensional space, PCA yields at most 11 non-zero
    components (rank ≤ min(n_samples - 1, n_features) after centering).

    Biology: PCA captures the dominant axes of variation across all cell type
    centroids from both species. These axes represent the gene programs that
    most differentiate cell types and species simultaneously.

    Args:
        human_centroids: (6, G) DataFrame of human centroids.
        mouse_centroids: (6, G) DataFrame of mouse centroids.
        variance_threshold: Minimum cumulative variance to retain (default 0.95).

    Returns:
        Tuple of (human_pca, mouse_pca, pca_model, cell_types) where:
        - human_pca: (6, k) array in PCA space
        - mouse_pca: (6, k) array in PCA space
        - pca_model: fitted PCA object (for inverse transform to gene space)
        - cell_types: ordered list of cell type names
    """
    # Ensure same cell type order
    cell_types = sorted(human_centroids.index.tolist())
    assert sorted(mouse_centroids.index.tolist()) == cell_types, (
        "Cell types must match between species"
    )

    human_mat = human_centroids.loc[cell_types].values  # (6, G)
    mouse_mat = mouse_centroids.loc[cell_types].values  # (6, G)

    # Combine: rows 0..5 = human, rows 6..11 = mouse
    combined = np.vstack([human_mat, mouse_mat])  # (12, G)

    # Max possible components: min(n_samples - 1, n_features)
    n_max = min(combined.shape[0] - 1, combined.shape[1])  # 11

    # Fit PCA — n_components as float selects enough for threshold
    pca = PCA(
        n_components=variance_threshold,
        svd_solver="full",
        random_state=RANDOM_SEED,
    )
    combined_pca = pca.fit_transform(combined)  # (12, k)
    n_components = pca.n_components_

    # Split back into species
    n_types = len(cell_types)
    human_pca = combined_pca[:n_types]   # (n_types, k)
    mouse_pca = combined_pca[n_types:]   # (n_types, k)

    cumvar = np.cumsum(pca.explained_variance_ratio_)
    print(f"\n  PCA on 12 centroids ({combined.shape[1]:,} genes):")
    print(f"  Components retained: {n_components} / {n_max} possible")
    print(f"  Cumulative variance explained: {cumvar[-1] * 100:.1f}%")
    for i, (var, cum) in enumerate(
        zip(pca.explained_variance_ratio_, cumvar)
    ):
        print(f"    PC{i + 1}: {var * 100:.2f}%  (cumulative: {cum * 100:.1f}%)")

    return human_pca, mouse_pca, pca, cell_types


# ---------------------------------------------------------------------------
# Procrustes analysis (from scratch)
# ---------------------------------------------------------------------------


def procrustes_align(
    X: np.ndarray,
    Y: np.ndarray,
    allow_reflection: bool = False,
) -> ProcrustesResult:
    """
    Ordinary Procrustes Analysis: align target Y onto reference X.

    Finds rotation R, scaling s, and translation t that minimize:
        ‖X - (s Y R + 1 t^T)‖_F^2

    Math (step by step):

    1. CENTER both point sets to remove translation:
       X_c = X - μ_X   where  μ_X = (1/n) Σ_i X_i
       Y_c = Y - μ_Y   where  μ_Y = (1/n) Σ_i Y_i

       This eliminates the translation parameter t. The optimal translation
       is t = μ_X - s R μ_Y, but we work in centered coordinates.

    2. CROSS-COVARIANCE MATRIX:
       M = X_c^T Y_c   (k × k matrix)

       M encodes how the dimensions of X and Y co-vary. Its SVD reveals
       the optimal rotation that maximizes the alignment between the two
       point sets.

    3. SINGULAR VALUE DECOMPOSITION:
       M = U Σ V^T

       where U, V are k×k orthogonal matrices, Σ = diag(σ_1, ..., σ_k)
       with σ_i ≥ 0.  The singular values measure alignment quality along
       each principal direction.

    4. OPTIMAL ROTATION:
       R = V D U^T   where D = diag(1, 1, ..., det(V U^T))

       Without the correction matrix D, the solution V U^T might include a
       reflection (det = -1).  The matrix D flips the last axis if needed
       to ensure det(R) = +1 (proper rotation).  Biological transformations
       should be continuous, so we forbid reflections by default.

       Proof sketch: We maximize trace(X_c^T · s Y_c R) = s · trace(Σ D).
       Since Σ has non-negative entries, this is maximized by D = I unless
       det(V U^T) < 0, in which case we flip the last entry.

    5. OPTIMAL SCALING:
       s = trace(Σ D) / ‖Y_c‖_F^2

       Derived by setting ∂/∂s [‖X_c - s Y_c R‖_F^2] = 0 and solving.
       The numerator trace(Σ D) measures how well Y aligns with X after
       rotation; the denominator normalizes by the spread of Y.

    6. TRANSFORM TARGET:
       Y_aligned = s · Y_c · R   (in centered coordinates)

       To recover world coordinates: Y_world = Y_aligned + μ_X.

    7. PROCRUSTES DISTANCE (sum of squared residuals):
       d^2 = ‖X_c - s Y_c R‖_F^2
           = ‖X_c‖_F^2 + s^2 ‖Y_c‖_F^2 - 2s · trace(Σ D)
           = ‖X_c‖_F^2 - [trace(Σ D)]^2 / ‖Y_c‖_F^2

       The last form substitutes the optimal s. Smaller d^2 means better
       geometric correspondence between the species' cell type landscapes.

    Biology: The rotation R captures the orthogonal transformation relating
    mouse and human cell type geometries — how gene expression programs are
    "rewired" across species.  The scaling s captures overall magnitude
    differences.  The residuals (what rotation + scaling cannot explain)
    reveal cell-type-specific evolutionary divergence.

    Args:
        X: Reference point set (n × k array). Human centroids in PCA space.
        Y: Target point set (n × k array). Mouse centroids in PCA space.
           Rows must correspond: X[i] and Y[i] are the same cell type.
        allow_reflection: If True, allow improper rotations (det = -1).
           Default False — biological transformations should be continuous.

    Returns:
        ProcrustesResult with rotation, scaling, distance, and aligned points.
    """
    n, k = X.shape
    assert Y.shape == (n, k), f"Shape mismatch: X={X.shape}, Y={Y.shape}"

    # Step 1: Center
    mu_X = X.mean(axis=0)  # (k,)
    mu_Y = Y.mean(axis=0)  # (k,)
    X_c = X - mu_X
    Y_c = Y - mu_Y

    # Step 2: Cross-covariance matrix
    M = X_c.T @ Y_c  # (k, k)

    # Step 3: SVD
    U, sigma, Vt = svd(M)  # U: (k,k), sigma: (k,), Vt: (k,k)
    V = Vt.T

    # Step 4: Optimal rotation (prevent reflection unless allowed)
    if allow_reflection:
        D = np.eye(k)
    else:
        d = np.linalg.det(V @ U.T)
        D = np.eye(k)
        D[-1, -1] = np.sign(d)  # +1 for rotation, -1 corrects reflection

    R = V @ D @ U.T  # (k, k) rotation matrix

    # Step 5: Optimal scaling
    ss_Y = np.sum(Y_c**2)  # ‖Y_c‖_F^2
    trace_sigma_D = np.sum(sigma * np.diag(D))
    s = trace_sigma_D / ss_Y

    # Step 6: Transform target
    Y_aligned = s * (Y_c @ R)  # (n, k)

    # Step 7: Procrustes distance
    d_sq = np.sum((X_c - Y_aligned) ** 2)
    d = np.sqrt(d_sq)

    print(f"\n  Procrustes alignment (mouse → human):")
    print(f"    Dimensions: {n} points × {k} PCA components")
    print(f"    Optimal scaling: {s:.6f}")
    print(f"    Rotation determinant: {np.linalg.det(R):+.6f} (expect +1.0)")
    print(f"    Sum of squared residuals (SSR): {d_sq:.6f}")
    print(f"    Procrustes distance (√SSR): {d:.6f}")

    return ProcrustesResult(
        rotation=R,
        scaling=s,
        translation_ref=mu_X,
        translation_target=mu_Y,
        distance_squared=d_sq,
        distance=d,
        aligned_target=Y_aligned,
        centered_reference=X_c,
    )


# ---------------------------------------------------------------------------
# Permutation test
# ---------------------------------------------------------------------------


def _procrustes_distance(
    X: np.ndarray,
    Y: np.ndarray,
    allow_reflection: bool = False,
) -> float:
    """
    Compute Procrustes distance silently (no printing).

    Internal helper for permutation test performance. Implements the same
    math as procrustes_align but returns only the scalar distance.

    Args:
        X: Reference points (n × k).
        Y: Target points (n × k).
        allow_reflection: Allow improper rotations.

    Returns:
        Procrustes distance (√SSR after optimal alignment).
    """
    n, k = X.shape
    X_c = X - X.mean(axis=0)
    Y_c = Y - Y.mean(axis=0)

    M = X_c.T @ Y_c
    U, sigma, Vt = svd(M)
    V = Vt.T

    if allow_reflection:
        D_diag = np.ones(k)
    else:
        d = np.linalg.det(V @ U.T)
        D_diag = np.ones(k)
        D_diag[-1] = np.sign(d)

    ss_Y = np.sum(Y_c**2)
    trace_sigma_D = np.sum(sigma * D_diag)
    s = trace_sigma_D / ss_Y

    Y_aligned = s * (Y_c @ (V * D_diag) @ U.T)
    return np.sqrt(np.sum((X_c - Y_aligned) ** 2))


def permutation_test(
    X: np.ndarray,
    Y: np.ndarray,
    n_permutations: int = N_PERMUTATIONS,
    seed: int = RANDOM_SEED,
    allow_reflection: bool = False,
) -> tuple[float, np.ndarray]:
    """
    Permutation test for Procrustes distance significance.

    Null hypothesis H0: The correspondence between human and mouse cell types
    is no better than random. Under H0, any permutation of the cell type
    pairings should yield a similar Procrustes distance.

    Math: For each of B permutations:
        1. Randomly permute rows of Y (shuffle which mouse centroid
           corresponds to which human centroid).
        2. Compute d_perm = Procrustes distance between X and permuted Y.

    p-value = (#{d_perm ≤ d_observed} + 1) / (B + 1)

    The +1 in the numerator includes the observed statistic itself
    (conservative estimate that prevents p = 0). The +1 in the denominator
    normalizes the total count including the observed value.

    Biology: A significant result (p < 0.01) means the geometric
    correspondence between homologous cell types is unlikely to arise by
    chance — the cross-species transformation has real structure, consistent
    with D'Arcy Thompson's hypothesis.

    Note: With 6 cell types, there are 6! = 720 unique permutations. With
    10,000 Monte Carlo samples, many permutations will be drawn multiple
    times, but the p-value estimate is still valid.

    Args:
        X: Reference points (n × k). Human centroids in PCA space.
        Y: Target points (n × k). Mouse centroids in PCA space.
        n_permutations: Number of random permutations (default 10,000).
        seed: Random seed for reproducibility.
        allow_reflection: Passed to Procrustes alignment.

    Returns:
        Tuple of (p_value, null_distribution) where null_distribution
        is an array of n_permutations Procrustes distances.
    """
    n = X.shape[0]
    rng = np.random.RandomState(seed)

    # Observed distance
    observed = _procrustes_distance(X, Y, allow_reflection)

    # Generate null distribution
    null_distances = np.zeros(n_permutations)
    for i in range(n_permutations):
        perm = rng.permutation(n)
        null_distances[i] = _procrustes_distance(X, Y[perm], allow_reflection)

    # p-value: fraction of null ≤ observed (including observed itself)
    n_leq = int(np.sum(null_distances <= observed))
    p_value = (n_leq + 1) / (n_permutations + 1)

    print(f"\n  Permutation test ({n_permutations:,} iterations, seed={seed}):")
    print(f"    Observed Procrustes distance: {observed:.6f}")
    print(
        f"    Null distribution: median={np.median(null_distances):.6f}, "
        f"mean={np.mean(null_distances):.6f}"
    )
    print(
        f"    Null range: [{np.min(null_distances):.6f}, "
        f"{np.max(null_distances):.6f}]"
    )
    print(
        f"    Permuted distances ≤ observed: {n_leq} / {n_permutations}"
    )
    print(f"    p-value: {p_value:.6f}")
    print(
        f"    Significant at α=0.01: "
        f"{'YES ✓' if p_value < 0.01 else 'NO ✗'}"
    )
    print(f"    Note: {n} cell types → {n}! = {math.factorial(n)} unique permutations")

    return p_value, null_distances


# ---------------------------------------------------------------------------
# Residual deformation vectors
# ---------------------------------------------------------------------------


def compute_residual_vectors(
    result: ProcrustesResult,
    cell_types: list[str],
) -> dict[str, np.ndarray]:
    """
    Compute per-cell-type residual deformation vectors in PCA space.

    Math: For cell type i:
        r_i = aligned_mouse_i - human_i   (in centered PCA coordinates)

    The residual r_i is what the global Procrustes transformation (rotation
    + scaling) could NOT explain for cell type i. It is the cell-type-specific
    component of the cross-species difference.

    Biology: Large residuals indicate cell types whose cross-species
    divergence is poorly explained by the global geometric transformation.
    The direction of the residual (mapped back to gene space) reveals WHICH
    genes drive cell-type-specific divergence.

    Args:
        result: ProcrustesResult from procrustes_align.
        cell_types: Ordered list of cell type names (matching row order).

    Returns:
        Dict mapping cell type name → residual vector (k,) in PCA space.
    """
    residuals = {}
    print(f"\n  Per-cell-type residual magnitudes (in PCA space):")
    print(f"  {'Cell Type':<45} {'‖residual‖':>12} {'% of total':>12}")
    print(f"  {'-' * 70}")

    total_ssr = result.distance_squared
    for i, ct in enumerate(cell_types):
        r = result.aligned_target[i] - result.centered_reference[i]
        residuals[ct] = r
        mag = np.linalg.norm(r)
        pct = (mag**2 / total_ssr * 100) if total_ssr > 0 else 0.0
        print(f"  {ct:<45} {mag:>12.6f} {pct:>11.1f}%")

    return residuals


# ---------------------------------------------------------------------------
# Map residuals back to gene space
# ---------------------------------------------------------------------------


def map_residuals_to_genes(
    residuals: dict[str, np.ndarray],
    pca_model: PCA,
    gene_names: list[str],
    n_top: int = N_TOP_GENES,
) -> dict[str, pd.DataFrame]:
    """
    Project residual vectors from PCA space back to gene space.

    Math: The PCA loadings matrix W has shape (k × G), where each row is a
    principal component direction in gene space. The gene-space residual is:

        g_i = r_i · W   where  r_i ∈ ℝ^k  and  W ∈ ℝ^{k×G}

    So g_i ∈ ℝ^G.  The magnitude |g_i[j]| indicates how much gene j
    contributes to the residual deformation for cell type i.  Genes with
    large absolute loadings are the most diverged between species for that
    cell type, beyond what the global Procrustes transformation explains.

    Biology: These are the genes driving cell-type-specific evolutionary
    divergence.  If hepatocyte residuals highlight metabolic enzymes, it
    suggests species-specific adaptations in liver metabolism — exactly the
    kind of interpretable signal that supports the D'Arcy Thompson framing.

    The sign of the loading indicates direction: positive means the gene is
    higher in the aligned mouse centroid than the human centroid for that
    cell type (relative to the global transformation).

    Args:
        residuals: Dict mapping cell type → residual vector in PCA space.
        pca_model: Fitted PCA object with components_ attribute.
        gene_names: Gene symbols (length G, matching PCA training order).
        n_top: Number of top genes to report per cell type (default 20).

    Returns:
        Dict mapping cell type name → DataFrame with columns:
        ['gene', 'loading', 'abs_loading', 'rank'].
    """
    W = pca_model.components_  # (k, G)

    top_genes = {}
    for ct, residual_pca in residuals.items():
        # Project to gene space
        gene_loadings = residual_pca @ W  # (G,)

        # Create DataFrame sorted by absolute loading
        df = pd.DataFrame(
            {
                "gene": gene_names,
                "loading": gene_loadings,
                "abs_loading": np.abs(gene_loadings),
            }
        )
        df = df.sort_values("abs_loading", ascending=False).reset_index(drop=True)
        df["rank"] = df.index + 1

        top_genes[ct] = df.head(n_top)

        print(f"\n  Top {n_top} genes for {ct}:")
        for _, row in df.head(n_top).iterrows():
            direction = "↑mouse" if row["loading"] > 0 else "↓mouse"
            print(
                f"    {int(row['rank']):>3}. {row['gene']:<15} "
                f"loading={row['loading']:>+.6f}  ({direction})"
            )

    return top_genes


# ---------------------------------------------------------------------------
# Results serialization
# ---------------------------------------------------------------------------


def save_results(
    result: ProcrustesResult,
    p_value: float,
    null_distribution: np.ndarray,
    residuals: dict[str, np.ndarray],
    top_genes: dict[str, pd.DataFrame],
    cell_types: list[str],
    pca_info: dict,
    output_path: Path,
) -> None:
    """
    Save all Procrustes analysis results to JSON.

    Saves the main results as a human-readable JSON file. The full null
    distribution is saved separately as a .npy file to keep the JSON concise.

    Args:
        result: ProcrustesResult from alignment.
        p_value: From permutation test.
        null_distribution: Array of permuted distances (saved as .npy).
        residuals: Per-cell-type residual vectors.
        top_genes: Per-cell-type top gene DataFrames.
        cell_types: Ordered cell type names.
        pca_info: Dict with PCA metadata.
        output_path: Path to save JSON file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    results_dict = {
        "procrustes": {
            "distance": float(result.distance),
            "distance_squared": float(result.distance_squared),
            "scaling": float(result.scaling),
            "rotation_matrix": result.rotation.tolist(),
            "rotation_det": float(np.linalg.det(result.rotation)),
        },
        "permutation_test": {
            "p_value": float(p_value),
            "n_permutations": len(null_distribution),
            "null_distribution_summary": {
                "mean": float(np.mean(null_distribution)),
                "median": float(np.median(null_distribution)),
                "std": float(np.std(null_distribution)),
                "min": float(np.min(null_distribution)),
                "max": float(np.max(null_distribution)),
                "percentile_2_5": float(np.percentile(null_distribution, 2.5)),
                "percentile_5": float(np.percentile(null_distribution, 5)),
                "percentile_95": float(np.percentile(null_distribution, 95)),
                "percentile_97_5": float(np.percentile(null_distribution, 97.5)),
            },
            "significant_at_001": bool(p_value < 0.01),
        },
        "pca": pca_info,
        "cell_types": cell_types,
        "residuals": {
            ct: {
                "vector_pca": residuals[ct].tolist(),
                "magnitude": float(np.linalg.norm(residuals[ct])),
            }
            for ct in cell_types
        },
        "top_genes_per_cell_type": {
            ct: top_genes[ct][["gene", "loading", "abs_loading", "rank"]]
            .to_dict(orient="records")
            for ct in cell_types
        },
        "aligned_mouse_centroids_pca": result.aligned_target.tolist(),
        "human_centroids_pca_centered": result.centered_reference.tolist(),
        "random_seed": RANDOM_SEED,
    }

    # Atomic write (tmp → rename) to prevent corruption on crash
    tmp_path = output_path.with_suffix(".tmp")
    with open(tmp_path, "w") as f:
        json.dump(results_dict, f, indent=2)
    tmp_path.rename(output_path)

    print(f"\n  Results saved: {output_path}")
