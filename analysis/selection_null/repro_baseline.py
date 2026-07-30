#!/usr/bin/env python3
"""
Baseline reproduction + code-path lock for the selection-aware null.

Reproduces the published Layer-1 numbers FROM SCRATCH in this fresh env, using
the UNMODIFIED published code path (gate_lib.per_gene_corr for C; the conserved
top-quartile selection; cellwarp.procrustes pca_reduce_centroids ->
procrustes_align -> permutation_test for obs/null). This is exactly the machinery
the sigma-derangement null wraps, so reproducing 0.384 / 0.522 here both confirms
the checkout runs end-to-end and locks the baseline the wrapper must hit at
sigma = identity.

Run: ../cellwarp/.venv/bin/python repro_baseline.py
"""
import sys, io, contextlib
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "analysis" / "conserved_contribution"))
import gate_lib as G  # noqa: E402  (unmodified published library)
from cellwarp.procrustes import (  # noqa: E402  (unmodified published pipeline)
    pca_reduce_centroids, procrustes_align, permutation_test,
)

def obs_null_ratio(gene_id_list, h_df, m_df, n_perm=2000):
    """Verbatim re-implementation of run_gate.obs_null_ratio (unmodified pipeline)."""
    cols = [g for g in gene_id_list if g in h_df.columns]
    hc = h_df[cols]; mc = m_df[cols]
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        hp, mp, _, _ = pca_reduce_centroids(hc, mc, 0.95)
        res = procrustes_align(hp, mp)
        _, null = permutation_test(hp, mp, n_perm)
    return float(res.distance), float(np.median(null)), float(res.distance / np.median(null))

def main():
    h, m = G.load_centroids()
    types = list(h.index); genes = list(h.columns)
    H, M = h.values, m.values
    print(f"centroids: human {H.shape} mouse {M.shape} | types={len(types)} genes={len(genes)}")

    # C the gate way (per-gene Pearson across the 35 matched centroids)
    C = G.per_gene_corr(H, M, "pearson")
    valid = ~np.isnan(C)
    q75 = np.quantile(C[valid], 0.75)
    cons_mask = valid & (C >= q75)
    cons_ids = [genes[j] for j in range(len(genes)) if cons_mask[j]]
    print(f"valid C = {int(valid.sum())} (expect 15940) | Q75 = {q75:.4f} (expect ~0.592) "
          f"| n_conserved = {len(cons_ids)} (expect 3985)")

    allg = obs_null_ratio(genes, h, m)
    cons = obs_null_ratio(cons_ids, h, m)
    print(f"ALL-genes obs/null  = {allg[2]:.4f}   (published 0.522)  [d={allg[0]:.3f} null_med={allg[1]:.3f}]")
    print(f"CONSERVED obs/null  = {cons[2]:.4f}   (published 0.384)  [d={cons[0]:.3f} null_med={cons[1]:.3f}]")

    # vectorized C equivalence check (the wrapper uses a fast vectorized C in the
    # N>=1000 loop; selection math is not "the method" -- prove it is identical)
    Hc = H - H.mean(0, keepdims=True); Mc = M - M.mean(0, keepdims=True)
    num = (Hc * Mc).sum(0)
    den = np.sqrt((Hc**2).sum(0) * (Mc**2).sum(0))
    with np.errstate(invalid="ignore", divide="ignore"):
        Cvec = np.where(den > 0, num / den, np.nan)
    same_valid = np.array_equal(np.isnan(C), np.isnan(Cvec))
    max_abs = np.nanmax(np.abs(C - Cvec))
    cons_vec = set(genes[j] for j in range(len(genes)) if (~np.isnan(Cvec[j])) and (Cvec[j] >= np.quantile(Cvec[~np.isnan(Cvec)], 0.75)))
    print(f"vectorized-C check: same NaN mask={same_valid} | max|dC|={max_abs:.2e} | "
          f"conserved-set identical={cons_vec == set(cons_ids)}")

if __name__ == "__main__":
    main()
