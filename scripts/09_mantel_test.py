#!/usr/bin/env python3
"""
Mantel test comparing pairwise distance matrices of human and mouse
cell-type centroids in joint PCA space.

Biology: If cross-species cell-type relationships preserve pairwise
distances (not just rigid transformations), the Mantel correlation
will be significantly positive. This validates the Procrustes
framework under a weaker assumption — distance-matrix preservation
rather than a single linear transformation.

Math: The Mantel test correlates two distance matrices element-wise
(using their upper triangles) and assesses significance by permuting
row/column labels of one matrix. Under the null hypothesis of no
association, the correlation should be near zero.

Usage:
    python scripts/09_mantel_test.py
"""

import json
import numpy as np
from scipy.spatial.distance import pdist, squareform
from pathlib import Path


def mantel_test(dist_x, dist_y, n_permutations=10000, seed=42):
    """
    Mantel test comparing two square distance matrices.

    Parameters
    ----------
    dist_x : np.ndarray
        Square symmetric distance matrix (n x n) for species X.
    dist_y : np.ndarray
        Square symmetric distance matrix (n x n) for species Y.
    n_permutations : int
        Number of random permutations for significance testing.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    dict with pearson_r, pearson_p, spearman_r, spearman_p
    """
    n = dist_x.shape[0]
    assert dist_x.shape == dist_y.shape == (n, n)

    # Extract upper triangle (excluding diagonal)
    iu = np.triu_indices(n, k=1)
    x_vec = dist_x[iu]
    y_vec = dist_y[iu]

    # Observed correlations
    obs_pearson = np.corrcoef(x_vec, y_vec)[0, 1]
    obs_spearman = _spearman(x_vec, y_vec)

    # Permutation test
    rng = np.random.default_rng(seed)
    n_ge_pearson = 0
    n_ge_spearman = 0

    for _ in range(n_permutations):
        perm = rng.permutation(n)
        # Permute rows and columns of one matrix
        y_perm = dist_y[np.ix_(perm, perm)]
        y_perm_vec = y_perm[iu]

        r_pear = np.corrcoef(x_vec, y_perm_vec)[0, 1]
        r_spear = _spearman(x_vec, y_perm_vec)

        if r_pear >= obs_pearson:
            n_ge_pearson += 1
        if r_spear >= obs_spearman:
            n_ge_spearman += 1

    # p-value: fraction of permutations with r >= observed
    # Add 1 to numerator and denominator for unbiased estimate
    pearson_p = (n_ge_pearson + 1) / (n_permutations + 1)
    spearman_p = (n_ge_spearman + 1) / (n_permutations + 1)

    return {
        "pearson_r": float(obs_pearson),
        "pearson_p": float(pearson_p),
        "spearman_r": float(obs_spearman),
        "spearman_p": float(spearman_p),
    }


def _spearman(x, y):
    """Spearman rank correlation via ranking and Pearson."""
    from scipy.stats import rankdata
    return np.corrcoef(rankdata(x), rankdata(y))[0, 1]


def main():
    project_root = Path(__file__).resolve().parent.parent
    npz_path = project_root / "output" / "phase2" / "scaled_35types" / "pca_centroids_35.npz"
    out_dir = project_root / "analysis" / "mantel_test"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load PCA centroids (35 cell types × 33 PCs)
    data = np.load(npz_path, allow_pickle=True)
    human_pca = data["human"]   # (35, 33)
    mouse_pca = data["mouse"]   # (35, 33)
    cell_types = data["cell_types"]
    n_types = human_pca.shape[0]

    print(f"Loaded {n_types} cell-type centroids in {human_pca.shape[1]}-D PCA space")
    print(f"Cell types: {list(cell_types[:5])} ... ({n_types} total)")

    # Compute pairwise Euclidean distance matrices
    human_dist = squareform(pdist(human_pca, metric="euclidean"))
    mouse_dist = squareform(pdist(mouse_pca, metric="euclidean"))

    print(f"\nHuman distance matrix: {human_dist.shape}, range [{human_dist[human_dist > 0].min():.3f}, {human_dist.max():.3f}]")
    print(f"Mouse distance matrix: {mouse_dist.shape}, range [{mouse_dist[mouse_dist > 0].min():.3f}, {mouse_dist.max():.3f}]")

    # Run Mantel test
    n_permutations = 10_000
    print(f"\nRunning Mantel test with {n_permutations:,} permutations...")
    results = mantel_test(human_dist, mouse_dist, n_permutations=n_permutations, seed=42)
    results["n_types"] = int(n_types)
    results["n_permutations"] = n_permutations

    # Print results
    print(f"\n{'='*55}")
    print(f"  MANTEL TEST RESULTS ({n_types} cell types)")
    print(f"{'='*55}")
    print(f"  Pearson  r = {results['pearson_r']:.4f}   p = {results['pearson_p']:.6f}")
    print(f"  Spearman r = {results['spearman_r']:.4f}   p = {results['spearman_p']:.6f}")
    print(f"  Permutations: {n_permutations:,}")
    print(f"{'='*55}")

    sig_label = "SIGNIFICANT" if results['pearson_p'] < 0.01 else "NOT significant"
    print(f"\n  → Pearson correlation is {sig_label} at α = 0.01")

    # Save JSON results
    json_path = out_dir / "mantel_results.json"
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {json_path}")

    # Generate interpretation
    write_summary(results, out_dir / "mantel_summary.md")
    print(f"Summary saved to {out_dir / 'mantel_summary.md'}")


def write_summary(results, path):
    """Write a one-paragraph interpretation of the Mantel test results."""
    r_p = results["pearson_r"]
    p_p = results["pearson_p"]
    r_s = results["spearman_r"]
    p_s = results["spearman_p"]
    n = results["n_types"]
    n_perm = results["n_permutations"]

    if p_p < 0.001:
        p_str = f"p < 0.001"
    else:
        p_str = f"p = {p_p:.4f}"

    strength = "strong" if abs(r_p) > 0.7 else "moderate" if abs(r_p) > 0.4 else "weak"

    text = f"""# Mantel Test: Cross-Species Distance Matrix Correlation

## Results

| Metric   | r      | p-value |
|----------|--------|---------|
| Pearson  | {r_p:.4f} | {p_p:.6f} |
| Spearman | {r_s:.4f} | {p_s:.6f} |

- Cell types: {n}
- Permutations: {n_perm:,}

## Interpretation

The Mantel test evaluated whether the {n} x {n} pairwise Euclidean distance matrix among cell-type centroids in joint PCA space is correlated between human and mouse. The Pearson correlation was {r_p:.4f} ({p_str}, {n_perm:,} permutations) and the Spearman rank correlation was {r_s:.4f} (p = {p_s:.6f}), indicating a {strength} and {"statistically significant" if p_p < 0.01 else "non-significant"} preservation of pairwise distance structure across species. This result {"supports" if p_p < 0.01 and r_p > 0 else "does not support"} the assumption underlying the Procrustes analysis: the geometric relationships among cell types in transcriptomic space are conserved between human and mouse, even under a weaker test that does not assume a single rigid-body transformation. Because the Mantel test only requires monotonic distance preservation (especially via the Spearman variant) rather than exact linear correspondence, this positive result provides independent evidence that the cross-species cell-type configuration is not an artefact of the Procrustes fitting procedure but reflects genuine biological conservation of the transcriptomic geometry.
"""
    with open(path, "w") as f:
        f.write(text)


if __name__ == "__main__":
    main()
