"""
CellWarp — Liang-Wagner Treeness Module

Implements the treeness test from Liang & Wagner (2015, Nature Communications)
to assess whether cell type transcriptomes have tree-like geometric structure.
Ported from the authors' R implementation (github.com/cloverliang/LCE,
src/TreenessTest.R).

Biology
-------
Liang & Wagner (2015) demonstrated that normal cell type transcriptomes exhibit
tree-like structure — the pairwise distance matrix fits a phylogenetic tree
topology. This geometric property breaks down in cancer. CellWarp's Procrustes
rigidity measures a different geometric property: how geometrically conserved
cell type positions are across species. Comparing these two independent geometric
frameworks reveals whether they capture the same or distinct aspects of cell
identity geometry.

Math
----
Delta statistic for a tetrad {a, b, c, d} (four-point condition):
  Given six pairwise Euclidean distances, form three partition sums:
    H1 = d(a,b) + d(c,d)
    H2 = d(a,c) + d(b,d)
    H3 = d(a,d) + d(b,c)
  Sort: H(1) <= H(2) <= H(3)
  delta = (H(3) - H(2)) / (H(3) - H(1))

  delta = 1: perfect tree (four-point condition satisfied exactly)
  delta = 0: maximally non-tree-like (all three partition sums equal)

Analytic p-value (Holland et al. 2002):
  Under H0 (no tree constraint, multivariate normal), delta has CDF:
    F(delta) = (3/pi) * (arctan((2*delta - 1) / sqrt(3)) + pi/6)
  P-value for tree structure: p = 1 - F(delta).
  Small p indicates evidence of tree-like geometry.

Aggregation via Storey's method (Storey 2002):
  pi0 = estimated proportion of tetrads with no tree structure (true nulls).
  1 - pi0 = fraction of tetrads with genuine tree structure.

References
----------
- Holland BR et al. (2002) "delta plots: a tool for analyzing phylogenetic
  distance data" Mol Biol Evol 19(12):2051-2059
- Liang C & Wagner G (2015) "Functional tradeoffs and the evolution of
  cell type transcriptome similarity" Nat Commun 6:8320
- Storey JD (2002) "A direct approach to false discovery rates"
  J R Stat Soc B 64(3):479-498
"""

from __future__ import annotations

from itertools import combinations

import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.spatial.distance import pdist, squareform


def compute_all_deltas(
    centroids: np.ndarray,
) -> tuple[np.ndarray, list[tuple[int, ...]]]:
    """Compute delta statistic for all C(n, 4) tetrads.

    Precomputes the full pairwise Euclidean distance matrix, then evaluates
    the four-point condition for every possible tetrad.

    Math: For tetrad {a, b, c, d}, the three partition sums are:
      H1 = d(a,b) + d(c,d)
      H2 = d(a,c) + d(b,d)
      H3 = d(a,d) + d(b,c)
    After sorting H(1) <= H(2) <= H(3):
      delta = (H(3) - H(2)) / (H(3) - H(1))

    Args:
        centroids: (n, G) array where n = number of cell types, G = genes.

    Returns:
        (deltas, tetrads) where:
        - deltas: array of length C(n, 4) with delta values in [0, 1]
        - tetrads: list of (i, j, k, l) index tuples
    """
    n = centroids.shape[0]
    dist_vec = pdist(centroids, metric="euclidean")
    dist_matrix = squareform(dist_vec)

    tetrads = list(combinations(range(n), 4))
    deltas = np.empty(len(tetrads))

    for idx, (i, j, k, l) in enumerate(tetrads):
        d_ij = dist_matrix[i, j]
        d_ik = dist_matrix[i, k]
        d_il = dist_matrix[i, l]
        d_jk = dist_matrix[j, k]
        d_jl = dist_matrix[j, l]
        d_kl = dist_matrix[k, l]

        H1 = d_ij + d_kl
        H2 = d_ik + d_jl
        H3 = d_il + d_jk

        vals = sorted([H1, H2, H3])
        denom = vals[2] - vals[0]
        deltas[idx] = (vals[2] - vals[1]) / denom if denom > 1e-15 else 0.0

    return deltas, tetrads


def holland_pvalues(deltas: np.ndarray) -> np.ndarray:
    """Analytic p-values for tree structure (Holland et al. 2002).

    Under H0 (no tree constraint, multivariate normal), delta has CDF:
      F(delta) = (3/pi) * (arctan((2*delta - 1) / sqrt(3)) + pi/6)

    The p-value for tree structure (H1: delta is large) is 1 - F(delta).
    Small p = evidence of tree structure.

    Biology: A small p-value for a tetrad means those four cell types have
    pairwise distances consistent with a phylogenetic tree. Aggregating
    across all tetrads reveals whether the full set of cell types has
    tree-like geometry.

    Args:
        deltas: Array of delta values in [0, 1].

    Returns:
        Array of p-values (same length as deltas). Small p = tree-like.
    """
    cdf = (3 / np.pi) * (np.arctan((2 * deltas - 1) / np.sqrt(3)) + np.pi / 6)
    return 1 - cdf


def estimate_pi0_storey(
    pvalues: np.ndarray,
    lambda_range: np.ndarray | None = None,
) -> float:
    """Estimate pi0 (proportion of true nulls) using Storey's method.

    For a range of lambda thresholds, estimates:
      pi0_hat(lambda) = #{p_i > lambda} / (m * (1 - lambda))

    Fits a natural cubic spline to (lambda, pi0_hat) and evaluates at
    max(lambda). Result clipped to [1/m, 1.0].

    Math: Under the complete null, all p-values are uniform and
    pi0_hat(lambda) = 1 for all lambda. Under enrichment for tree
    structure, p-values concentrate near 0, and pi0_hat decreases
    as lambda increases because fewer p-values exceed the threshold.

    Args:
        pvalues: Array of p-values from holland_pvalues().
        lambda_range: Threshold values (default: 0.05, 0.10, ..., 0.95).

    Returns:
        Estimated pi0 in [0, 1]. Low pi0 means most tetrads show tree
        structure; pi0 near 1 means no tree structure.
    """
    if lambda_range is None:
        lambda_range = np.arange(0.05, 0.96, 0.05)

    m = len(pvalues)
    pi0_hat = np.array(
        [np.sum(pvalues > lam) / (m * (1 - lam)) for lam in lambda_range]
    )

    # Natural cubic spline fit, evaluate at max(lambda)
    try:
        spline = UnivariateSpline(
            lambda_range, pi0_hat, k=3, s=len(lambda_range)
        )
        pi0 = float(spline(max(lambda_range)))
    except Exception:
        # Fallback: median of high-lambda estimates
        pi0 = float(np.median(pi0_hat[-5:]))

    return float(np.clip(pi0, 1 / m, 1.0))


def per_celltype_treeness(
    deltas: np.ndarray,
    tetrads: list[tuple[int, ...]],
    n_types: int,
) -> np.ndarray:
    """Compute per-cell-type treeness score as mean delta over containing tetrads.

    For cell type i, computes mean(delta) across all C(n-1, 3) tetrads that
    include cell type i. This measures how much including this cell type
    contributes to or degrades tree fit.

    Biology: Cell types whose inclusion improves tree fit (higher mean delta)
    occupy positions in expression space that are consistent with a
    phylogenetic branching pattern. Cell types that degrade tree fit occupy
    positions that violate tree geometry — potentially reflecting convergent
    evolution, horizontal gene transfer effects, or unique functional
    specializations.

    Args:
        deltas: Array of delta values for all tetrads.
        tetrads: List of (i, j, k, l) tuples (indices into cell type list).
        n_types: Total number of cell types.

    Returns:
        Array of length n_types with mean delta per cell type.
    """
    # Build index: which tetrads contain each cell type
    type_to_tetrad_idx: dict[int, list[int]] = {i: [] for i in range(n_types)}
    for idx, tet in enumerate(tetrads):
        for t in tet:
            type_to_tetrad_idx[t].append(idx)

    scores = np.empty(n_types)
    for i in range(n_types):
        indices = type_to_tetrad_idx[i]
        scores[i] = np.mean(deltas[indices])

    return scores
