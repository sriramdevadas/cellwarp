"""
Tests for permutation test calibration and power.

Covers: null calibration (uniform p-values), power (planted signal
rejection rate), and seed sensitivity.
"""

import numpy as np
import pytest
from scipy import stats

from cellwarp.procrustes import permutation_test


# ---------------------------------------------------------------------------
# Null calibration
# ---------------------------------------------------------------------------

class TestNullCalibration:
    """Under H0 (no signal), p-values should be approximately uniform."""

    def test_pvalue_uniformity(self):
        """
        Generate 200 null datasets. Run permutation test on each.
        Assert p-values are approximately uniform (KS test, p > 0.01).
        """
        n_datasets = 200
        n_perms = 100
        p_values = []

        for i in range(n_datasets):
            rng = np.random.RandomState(i)
            X = rng.randn(5, 8)
            Y = rng.randn(5, 8)
            p, _ = permutation_test(X, Y, n_permutations=n_perms, seed=i + 1000)
            p_values.append(p)

        p_values = np.array(p_values)

        # KS test against uniform distribution
        ks_stat, ks_p = stats.kstest(p_values, "uniform")
        assert ks_p > 0.01, (
            f"p-values not uniform: KS stat={ks_stat:.3f}, "
            f"KS p={ks_p:.4f}. Mean p={np.mean(p_values):.3f}"
        )


# ---------------------------------------------------------------------------
# Power
# ---------------------------------------------------------------------------

class TestPower:
    """With planted signal, rejection rate should be high."""

    def test_rejection_rate(self):
        """
        Generate 200 signal datasets. Assert rejection rate > 80% at alpha=0.05.
        """
        n_datasets = 200
        n_perms = 100
        alpha = 0.05
        rejections = 0

        for i in range(n_datasets):
            rng = np.random.RandomState(i)
            # Shared structure with noise
            centroids = rng.randn(5, 8) * 3.0
            X = centroids + rng.randn(5, 8) * 0.3
            Y = centroids + rng.randn(5, 8) * 0.3

            p, _ = permutation_test(X, Y, n_permutations=n_perms, seed=i + 2000)
            if p < alpha:
                rejections += 1

        rate = rejections / n_datasets
        assert rate > 0.80, (
            f"Power too low: {rate:.2%} rejections at alpha={alpha}, "
            f"expected > 80%"
        )


# ---------------------------------------------------------------------------
# Seed sensitivity
# ---------------------------------------------------------------------------

class TestSeedSensitivity:
    """p-values should be consistent across different seeds."""

    def test_pvalue_consistency(self):
        """
        Run same analysis with 5 different seeds. Assert p-values are
        similar (within ±0.05 of each other, given 500 permutations).
        """
        rng = np.random.RandomState(42)
        # Create moderately correlated data
        centroids = rng.randn(5, 8) * 2.0
        X = centroids + rng.randn(5, 8) * 0.5
        Y = centroids + rng.randn(5, 8) * 0.5

        p_values = []
        for seed in [10, 20, 30, 40, 50]:
            p, _ = permutation_test(X, Y, n_permutations=500, seed=seed)
            p_values.append(p)

        p_values = np.array(p_values)
        spread = np.max(p_values) - np.min(p_values)
        assert spread < 0.05, (
            f"p-values too variable across seeds: {p_values}, spread={spread:.4f}"
        )
