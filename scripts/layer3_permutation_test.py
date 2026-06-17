#!/usr/bin/env python3
"""
Layer 3 Global Permutation Test — Eigenvalue Conservation

Tests whether the observed mean Pearson r between matched human and mouse
eigenvalue profiles (r=0.953) is significantly greater than expected under
the null hypothesis of no cell-type-specific eigenvalue conservation.

Null model: Shuffle the assignment of mouse eigenvalue profiles to cell type
labels, then recompute mean Pearson r across all 35 pairs. This preserves
marginal distributions but breaks the cell-type-specific pairing.

Matches the Layer 2 methodology (label-shuffle permutation test) from
t3b_ellipsoid_alignment.py.
"""

from __future__ import annotations

import json
import time
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT = Path(__file__).resolve().parent.parent
INPUT_CSV = PROJECT / "output" / "mechanistic" / "ellipsoid_alignment" / "35type_eigenvalue_conservation.csv"
OUTPUT_DIR = PROJECT / "output" / "layer3_permutation"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

SEED = 42
N_PERM = 10_000

# ---------------------------------------------------------------------------
# Step 1: Load eigenvalue profiles and confirm observed statistic
# ---------------------------------------------------------------------------
print("=" * 70)
print("LAYER 3 GLOBAL PERMUTATION TEST")
print("=" * 70)

print(f"\nSource: {INPUT_CSV}")
df = pd.read_csv(INPUT_CSV)
n_types = len(df)

print(f"Cell types loaded: {n_types}")

# Parse the normalized eigenvalue profiles (stored as string lists)
import ast

human_profiles = np.array([ast.literal_eval(row) for row in df["human_norm"]])
mouse_profiles = np.array([ast.literal_eval(row) for row in df["mouse_norm"]])

print(f"Profile shape: {human_profiles.shape} (n_types × top_k eigenvalues)")

# Recompute observed mean Pearson r to verify
observed_rs = []
for i in range(n_types):
    r, _ = stats.pearsonr(human_profiles[i], mouse_profiles[i])
    observed_rs.append(r)

observed_rs = np.array(observed_rs)
observed_mean_r = float(np.mean(observed_rs))

print(f"\nObserved per-type Pearson r range: {observed_rs.min():.3f} – {observed_rs.max():.3f}")
print(f"Observed mean Pearson r: {observed_mean_r:.6f}")

# Confirm match with expected 0.953
if abs(observed_mean_r - 0.953) < 0.001:
    print("  → CONFIRMED: matches expected r=0.953")
    status = "CONFIRMED"
else:
    print(f"  → MISMATCH: expected 0.953, got {observed_mean_r:.3f}")
    status = "MISMATCH"

# Also verify against CSV column
csv_mean_r = float(df["pearson_r"].mean())
print(f"  CSV column mean: {csv_mean_r:.6f} (cross-check)")

# ---------------------------------------------------------------------------
# Step 2: Permutation test — shuffle mouse labels
# ---------------------------------------------------------------------------
print(f"\n{'=' * 70}")
print(f"Permutation test: {N_PERM} iterations")
print("=" * 70)

rng = np.random.RandomState(SEED)
null_mean_rs = np.zeros(N_PERM)

t0 = time.time()
for i in range(N_PERM):
    perm = rng.permutation(n_types)
    perm_rs = []
    for j in range(n_types):
        r, _ = stats.pearsonr(human_profiles[j], mouse_profiles[perm[j]])
        perm_rs.append(r)
    null_mean_rs[i] = np.mean(perm_rs)

elapsed = time.time() - t0
print(f"  Runtime: {elapsed:.1f}s")

# ---------------------------------------------------------------------------
# Step 3: Empirical p-value
# ---------------------------------------------------------------------------
n_geq = int(np.sum(null_mean_rs >= observed_mean_r)) + 1  # +1 for observed
empirical_p = n_geq / (N_PERM + 1)

null_mean = float(np.mean(null_mean_rs))
null_std = float(np.std(null_mean_rs))

print(f"\n{'=' * 70}")
print("RESULTS")
print("=" * 70)
print(f"  Observed mean r:      {observed_mean_r:.6f}  ({status} with 0.953)")
print(f"  Permutation null mean: {null_mean:.6f}")
print(f"  Permutation null SD:   {null_std:.6f}")
print(f"  Empirical p-value:     {empirical_p}")
print(f"  n_permutations:        {N_PERM}")
print(f"  n_types:               {n_types}")
print(f"  Count >= observed:     {n_geq - 1} / {N_PERM}  (+1 for observed)")

# ---------------------------------------------------------------------------
# Step 4: Save output
# ---------------------------------------------------------------------------
results = {
    "observed_mean_r": observed_mean_r,
    "permutation_mean": null_mean,
    "permutation_std": null_std,
    "empirical_p": empirical_p,
    "n_permutations": N_PERM,
    "n_types": n_types,
    "date_run": str(date.today()),
    "per_type_r_min": float(observed_rs.min()),
    "per_type_r_max": float(observed_rs.max()),
    "source_file": str(INPUT_CSV.relative_to(PROJECT)),
    "seed": SEED,
}

out_path = OUTPUT_DIR / "layer3_permutation_results.json"
with open(out_path, "w") as f:
    json.dump(results, f, indent=2)

print(f"\n  Output saved to: {out_path}")
print("\nDone.")
