#!/usr/bin/env python3
"""RIRA .X normalization diagnostic.

Reports:
  1. Value range (min, max, median, p95)
  2. Fraction of exact zeros
  3. Per-cell sum distribution (mean, sd, median, p5, p95)
  4. Per-gene variance distribution (histogram bins)
  5. Presence of negatives
  6. Signature check against common conventions:
     - Raw UMI counts (integers, per-cell sum = library size)
     - normalize_total(target_sum) → log1p (values ≤ log1p(target_sum), sum(exp(x)-1)=target_sum)
     - SCTransform Pearson residuals (signed, std≈1 per gene)
     - SCTransform corrected counts (integer-like, non-fixed lib size)
     - z-scoring (signed, std=1 per gene)

Uses a 20,000-cell random sample for column-wise stats to keep memory
manageable; full-matrix stats for sparse operations that are O(nnz).
"""
from __future__ import annotations
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

RIRA_CONV = Path(__file__).resolve().parent.parent.parent / "data/macaque/rira/converted"
SEED = 42
N_CELL_SAMPLE = 20_000


def main():
    print("Loading RIRA matrix…", flush=True)
    with gzip.open(RIRA_CONV / "matrix.mtx.gz", "rt") as f:
        mtx = sio.mmread(f)
    X = sp.csr_matrix(mtx.T)  # (cells, genes)
    n_cells, n_genes = X.shape
    print(f"Shape: {n_cells:,} cells × {n_genes:,} genes; nnz = {X.nnz:,}", flush=True)

    vals = X.data
    # --- 1. Value range ------------------------------------------------
    print("\n[1] Value range (over non-zero entries)")
    # Also include zeros: min over ALL values
    total = n_cells * n_genes
    zero_count = total - X.nnz
    min_all = 0.0 if zero_count > 0 else vals.min()
    print(f"  min (all):           {min_all:.6f}")
    print(f"  min (non-zero):      {vals.min():.6f}")
    print(f"  max:                 {vals.max():.6f}")
    print(f"  median (non-zero):   {np.median(vals):.6f}")
    print(f"  p95 (non-zero):      {np.percentile(vals, 95):.6f}")
    print(f"  p99 (non-zero):      {np.percentile(vals, 99):.6f}")
    print(f"  p99.9 (non-zero):    {np.percentile(vals, 99.9):.6f}")
    print(f"  non-integer fraction (non-zero): {(vals != vals.astype(int)).mean():.4%}")

    # --- 2. Fraction of zeros ------------------------------------------
    print("\n[2] Zero structure")
    print(f"  total entries:       {total:,}")
    print(f"  non-zero entries:    {X.nnz:,}")
    print(f"  exact zeros:         {zero_count:,}  ({zero_count/total:.4%})")

    # --- 5. Negatives --------------------------------------------------
    print("\n[5] Negative values")
    n_neg = int((vals < 0).sum())
    print(f"  negative values: {n_neg:,}  ({'none' if n_neg == 0 else f'min={vals.min():.4f}'})")

    # --- 3. Per-cell sum distribution ----------------------------------
    print("\n[3] Per-cell sum distribution (all cells)")
    # .sum(axis=1) on sparse returns a matrix
    cell_sums = np.asarray(X.sum(axis=1)).flatten()
    print(f"  mean:    {cell_sums.mean():.3f}")
    print(f"  std:     {cell_sums.std():.3f}")
    print(f"  median:  {np.median(cell_sums):.3f}")
    print(f"  min:     {cell_sums.min():.3f}")
    print(f"  max:     {cell_sums.max():.3f}")
    print(f"  p5:      {np.percentile(cell_sums, 5):.3f}")
    print(f"  p95:     {np.percentile(cell_sums, 95):.3f}")
    # Check if log-normalized: sum(exp(x)-1) per cell should be constant if so.
    # Only feasible on cells with max value <50 (else overflow).
    small_cells = np.where(
        np.asarray(X.max(axis=1).todense()).flatten() < 30
    )[0]
    print(f"\n  Testing log-normalized signature (cells with max<30, n={len(small_cells):,}):")
    if len(small_cells) > 0:
        # Sample 100 such cells
        rng = np.random.default_rng(SEED)
        sel = rng.choice(small_cells, size=min(100, len(small_cells)), replace=False)
        exp_sums = []
        for i in sel:
            row = X.getrow(i).toarray().flatten()
            nz = row[row > 0]
            if len(nz) > 0:
                exp_sums.append(np.expm1(nz).sum())
        exp_sums = np.array(exp_sums)
        print(f"    sum(exp(x)-1): mean={exp_sums.mean():.1f}, median={np.median(exp_sums):.1f}, "
              f"std={exp_sums.std():.1f}")
        print(f"    range: [{exp_sums.min():.1f}, {exp_sums.max():.1f}]")
        # If log1p(norm-total(T)), this should be very close to T uniformly.

    # --- 4. Per-gene variance distribution -----------------------------
    print("\n[4] Per-gene variance distribution (20,000-cell sample)")
    rng = np.random.default_rng(SEED)
    sel_cells = np.sort(rng.choice(n_cells, size=N_CELL_SAMPLE, replace=False))
    X_sub = X[sel_cells]  # (20000, n_genes)
    # Per-gene mean and variance
    mean_g = np.asarray(X_sub.mean(axis=0)).flatten()
    # Var(X) = E[X²] - E[X]²
    X_sq = X_sub.multiply(X_sub)
    mean_sq_g = np.asarray(X_sq.mean(axis=0)).flatten()
    var_g = mean_sq_g - mean_g ** 2
    var_g = np.maximum(var_g, 0)
    std_g = np.sqrt(var_g)
    # Focus on genes with nonzero variance
    nz_mask = var_g > 1e-12
    print(f"  Genes with nonzero variance: {int(nz_mask.sum()):,} / {n_genes:,}")
    print(f"  Per-gene mean:   p5={np.percentile(mean_g[nz_mask], 5):.4f}  "
          f"median={np.median(mean_g[nz_mask]):.4f}  p95={np.percentile(mean_g[nz_mask], 95):.4f}")
    print(f"  Per-gene std:    p5={np.percentile(std_g[nz_mask], 5):.4f}  "
          f"median={np.median(std_g[nz_mask]):.4f}  p95={np.percentile(std_g[nz_mask], 95):.4f}")
    print(f"  Per-gene var:    p5={np.percentile(var_g[nz_mask], 5):.4f}  "
          f"median={np.median(var_g[nz_mask]):.4f}  p95={np.percentile(var_g[nz_mask], 95):.4f}")
    # Fraction of genes with std ∈ [0.8, 1.2] (z-score signature)
    unit_std = ((std_g > 0.8) & (std_g < 1.2)).sum()
    print(f"  Per-gene std ∈ [0.8, 1.2]: {unit_std:,} / {n_genes:,}  ({unit_std / n_genes:.2%})")
    # If z-scored, should be ~100% near 1.0.

    # --- 6. Convention inference ---------------------------------------
    print("\n[6] Convention inference")
    print("  Raw UMI counts:     ", end="")
    is_integer = (vals != vals.astype(int)).mean() < 0.01
    print("YES" if is_integer else f"NO (non-integer fraction = {(vals != vals.astype(int)).mean():.4%})")
    print("  log1p(normalize_total):  ", end="")
    # Max should be ≤ log1p(target_sum). If max > 20, this is ruled out.
    max_val = vals.max()
    if max_val > 20:
        print(f"NO (max={max_val:.1f} > 20; log1p would cap ≤ log1p(1e6) ≈ 13.8)")
    else:
        print(f"PLAUSIBLE (max={max_val:.1f})")
    print("  Pearson residuals:  ", end="")
    print("NO (no negatives)" if n_neg == 0 else "POSSIBLE")
    print("  z-scored:           ", end="")
    print("NO (no negatives)" if n_neg == 0 else f"MAYBE ({unit_std/n_genes:.1%} genes have std≈1)")
    print("  SCTransform corrected counts:  ", end="")
    # Signature: non-negative, non-integer but integer-like, per-cell sum
    # tracks library size roughly, no upper bound
    if n_neg == 0 and max_val > 50:
        print("CONSISTENT (non-negative, non-integer, large per-value max)")
    else:
        print("UNCLEAR")


if __name__ == "__main__":
    main()
