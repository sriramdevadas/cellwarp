"""
End-to-end integration smoke test.

Generates two fake species with planted geometric signal, runs the full
pipeline (PCA -> centroids -> Procrustes -> permutation), and verifies
the result is significant.

Must complete in <30 seconds.
"""

import numpy as np
import pytest

from cellwarp.procrustes import (
    _procrustes_distance,
    compute_residual_vectors,
    permutation_test,
    procrustes_align,
)


class TestEndToEnd:
    """Full pipeline on synthetic data."""

    def test_full_pipeline(self):
        """
        Generate two fake species: 5 cell types, 100 cells each, 500 genes.
        Plant geometric signal (shared centroid structure + Gaussian noise).
        Run PCA -> centroid -> Procrustes -> permutation.
        Assert obs/null < 1, p < 0.05, outputs well-formed.
        """
        rng = np.random.RandomState(42)
        n_types = 5
        n_cells = 100
        n_genes = 500
        n_pcs = 10
        cell_types = [f"type_{i}" for i in range(n_types)]

        # ── Step 1: Generate synthetic gene expression data ──
        # Shared centroid structure in gene space
        centroids = rng.randn(n_types, n_genes) * 2.0

        # Species 1 (human-like): cells around centroids
        human_data = np.vstack([
            centroids[i] + rng.randn(n_cells, n_genes) * 0.5
            for i in range(n_types)
        ])  # (500, 500)
        human_labels = np.repeat(cell_types, n_cells)

        # Species 2 (mouse-like): same centroids, rotated + noise
        R_true = np.linalg.qr(rng.randn(n_genes, n_genes))[0]
        mouse_data = np.vstack([
            (centroids[i] @ R_true.T) + rng.randn(n_cells, n_genes) * 0.5
            for i in range(n_types)
        ])
        mouse_labels = np.repeat(cell_types, n_cells)

        # ── Step 2: Compute centroids ──
        human_centroids = np.array([
            human_data[human_labels == ct].mean(axis=0) for ct in cell_types
        ])
        mouse_centroids = np.array([
            mouse_data[mouse_labels == ct].mean(axis=0) for ct in cell_types
        ])

        # ── Step 3: Joint PCA ──
        from sklearn.decomposition import PCA
        joint = np.vstack([human_centroids, mouse_centroids])
        pca = PCA(n_components=n_pcs, random_state=42)
        joint_pca = pca.fit_transform(joint)
        X = joint_pca[:n_types]  # human
        Y = joint_pca[n_types:]  # mouse

        # ── Step 4: Procrustes alignment ──
        result = procrustes_align(X, Y)

        assert isinstance(result.rotation, np.ndarray)
        assert result.rotation.shape == (n_pcs, n_pcs)
        assert result.distance >= 0
        assert result.scaling > 0

        # ── Step 5: Permutation test ──
        p_value, null_dist = permutation_test(
            X, Y, n_permutations=500, seed=42
        )
        observed = _procrustes_distance(X, Y)
        obs_null = observed / np.mean(null_dist)

        assert obs_null < 1.0, f"obs/null = {obs_null:.3f}, expected < 1"
        assert p_value < 0.05, f"p = {p_value:.4f}, expected < 0.05"
        assert len(null_dist) == 500

        # ── Step 6: Residual vectors ──
        residuals = compute_residual_vectors(result, cell_types)

        assert len(residuals) == n_types
        for ct in cell_types:
            assert ct in residuals
            assert residuals[ct].shape == (n_pcs,)

    def test_no_signal_pipeline(self):
        """
        Random data (no shared structure) should NOT be significant.
        """
        rng = np.random.RandomState(99)
        n_types, k = 5, 10

        X = rng.randn(n_types, k)
        Y = rng.randn(n_types, k)

        result = procrustes_align(X, Y)
        p_value, null_dist = permutation_test(
            X, Y, n_permutations=500, seed=99
        )

        assert p_value > 0.01, (
            f"p = {p_value:.4f}, random data should not be significant"
        )
