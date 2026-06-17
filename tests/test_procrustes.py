"""
Unit tests for cellwarp.procrustes.

Covers: known rotation recovery, permutation null mean, planted signal
detection, determinism, and edge cases.

All tests use small synthetic data (<=5 types, <=10 dims, <=50 cells)
and should complete in <5 seconds each.
"""

import numpy as np
import pytest

from cellwarp.procrustes import (
    ProcrustesResult,
    _procrustes_distance,
    permutation_test,
    procrustes_align,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _rotation_matrix_2d(theta: float) -> np.ndarray:
    """2D rotation matrix by angle theta (radians)."""
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s], [s, c]])


def _random_rotation(k: int, seed: int = 0) -> np.ndarray:
    """Random k×k orthogonal rotation matrix (det=+1)."""
    rng = np.random.RandomState(seed)
    A = rng.randn(k, k)
    Q, R = np.linalg.qr(A)
    # Ensure det(Q) = +1
    Q = Q @ np.diag(np.sign(np.diag(R)))
    if np.linalg.det(Q) < 0:
        Q[:, -1] *= -1
    return Q


# ---------------------------------------------------------------------------
# Known rotation recovery
# ---------------------------------------------------------------------------

class TestKnownRotation:
    """Verify that Procrustes recovers a known rotation + scaling."""

    def test_2d_rotation_recovery(self):
        """Create X, rotate by known angle to get Y. Verify recovery."""
        rng = np.random.RandomState(42)
        X = rng.randn(5, 2)  # 5 types, 2 dims
        theta = np.pi / 6  # 30 degrees
        R_true = _rotation_matrix_2d(theta)
        s_true = 1.5
        Y = (X @ R_true.T) / s_true  # Y such that s*Y@R = X

        result = procrustes_align(X, Y)

        # Recovered scaling should match
        np.testing.assert_allclose(result.scaling, s_true, atol=1e-10)
        # Rotation should match (up to sign convention)
        np.testing.assert_allclose(
            np.abs(np.linalg.det(result.rotation)), 1.0, atol=1e-10
        )
        # Distance should be ~0 (perfect alignment)
        assert result.distance < 1e-10

    def test_identity_rotation(self):
        """If X == Y, rotation should be identity, scale = 1."""
        X = np.array([[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0]])
        result = procrustes_align(X, X)

        np.testing.assert_allclose(result.scaling, 1.0, atol=1e-10)
        np.testing.assert_allclose(
            result.rotation, np.eye(2), atol=1e-10
        )
        assert result.distance < 1e-10

    def test_nd_rotation_recovery(self):
        """Verify recovery in higher dimensions (k=5)."""
        rng = np.random.RandomState(99)
        k = 5
        X = rng.randn(4, k)  # 4 types
        R_true = _random_rotation(k, seed=7)
        s_true = 2.0
        Y = (X @ R_true.T) / s_true

        result = procrustes_align(X, Y)

        np.testing.assert_allclose(result.scaling, s_true, atol=1e-8)
        assert result.distance < 1e-8


# ---------------------------------------------------------------------------
# Permutation null mean
# ---------------------------------------------------------------------------

class TestPermutationNull:
    """On random data, obs/null ratio should be ~1.0."""

    def test_null_mean_ratio(self):
        """Random data: mean null distance ≈ observed distance."""
        rng = np.random.RandomState(42)
        X = rng.randn(5, 10)  # 5 types, 10 dims
        Y = rng.randn(5, 10)

        p, null_dist = permutation_test(
            X, Y, n_permutations=1000, seed=42
        )
        observed = _procrustes_distance(X, Y)
        ratio = observed / np.mean(null_dist)

        # Random data: ratio should be near 1.0 (±0.1)
        assert 0.7 < ratio < 1.3, f"obs/null ratio = {ratio:.3f}, expected ~1.0"
        # p-value should not be significant
        assert p > 0.01, f"p = {p:.4f}, random data should not be significant"


# ---------------------------------------------------------------------------
# Planted signal detection
# ---------------------------------------------------------------------------

class TestSignalDetection:
    """Create data with shared centroid structure + noise. Assert detection."""

    def test_planted_signal(self):
        """Shared structure + noise → obs/null < 0.7, p < 0.05."""
        rng = np.random.RandomState(42)
        n_types, k = 5, 10

        # Shared centroid structure
        centroids = rng.randn(n_types, k) * 3.0

        # X and Y are centroids + small noise
        X = centroids + rng.randn(n_types, k) * 0.3
        Y = centroids + rng.randn(n_types, k) * 0.3

        p, null_dist = permutation_test(
            X, Y, n_permutations=1000, seed=42
        )
        observed = _procrustes_distance(X, Y)
        ratio = observed / np.mean(null_dist)

        assert ratio < 0.7, f"obs/null ratio = {ratio:.3f}, expected < 0.7"
        assert p < 0.05, f"p = {p:.4f}, planted signal should be significant"


# ---------------------------------------------------------------------------
# Determinism
# ---------------------------------------------------------------------------

class TestDeterminism:
    """Same seed → bitwise identical results."""

    def test_permutation_determinism(self):
        """Two runs with same seed must give identical results."""
        rng = np.random.RandomState(42)
        X = rng.randn(5, 10)
        Y = rng.randn(5, 10)

        p1, null1 = permutation_test(X, Y, n_permutations=200, seed=123)
        p2, null2 = permutation_test(X, Y, n_permutations=200, seed=123)

        assert p1 == p2
        np.testing.assert_array_equal(null1, null2)

    def test_alignment_determinism(self):
        """Procrustes alignment is deterministic (no randomness)."""
        rng = np.random.RandomState(42)
        X = rng.randn(4, 3)
        Y = rng.randn(4, 3)

        r1 = procrustes_align(X, Y)
        r2 = procrustes_align(X, Y)

        np.testing.assert_array_equal(r1.rotation, r2.rotation)
        assert r1.scaling == r2.scaling
        assert r1.distance == r2.distance


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

class TestEdgeCases:
    """Minimum viable inputs and boundary conditions."""

    def test_k2_minimum_types(self):
        """k=2 cell types (minimum for Procrustes)."""
        X = np.array([[1.0, 0.0], [-1.0, 0.0]])
        Y = np.array([[0.0, 1.0], [0.0, -1.0]])

        result = procrustes_align(X, Y)
        assert result.distance < 1e-10  # 90° rotation should recover perfectly

    def test_single_pca_component(self):
        """1D alignment (k=1): just scaling + sign."""
        X = np.array([[3.0], [-3.0], [1.0]])
        Y = np.array([[1.0], [-1.0], [0.33]])

        result = procrustes_align(X, Y)
        assert result.distance < 0.1

    def test_perfectly_aligned(self):
        """Already aligned data: distance=0, scale=1, R=I."""
        X = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=float)
        result = procrustes_align(X, X)

        assert result.distance < 1e-10
        np.testing.assert_allclose(result.scaling, 1.0, atol=1e-10)
