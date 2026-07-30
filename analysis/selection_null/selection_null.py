#!/usr/bin/env python3
"""
CellWarp Layer-1 — selection-aware null (correspondence-permutation, selection-propagated).

Pre-registered design (advisor-authored, 2026-06-28). Bounds how much of the
conserved-C quartile obs/null (published 0.384) is MANUFACTURED by selecting genes
on the conservation score C versus EARNED by genuine cross-species cell-identity
geometry.

THE PIPELINE IS UNTOUCHED. This wrapper only (1) re-pairs the 35 cell types under a
permutation sigma, (2) recomputes C_sigma, (3) re-selects the top-quartile genes, and
(4) re-runs the *unmodified* published obs/null:
    cellwarp.procrustes.pca_reduce_centroids -> procrustes_align -> permutation_test
and the per-gene C uses gate_lib's definition (np.corrcoef across the 35 centroids;
here vectorized — proven identical to gate_lib.per_gene_corr to 8e-16 in repro_baseline.py).

Per draw (mode=derangement, the PRIMARY test):
  1. sigma = derangement of the 35 types (no type maps to itself) -> destroys the
     genuine cross-species correspondence while preserving every gene's marginals
     and the entire selection procedure.
  2. M_sigma = mouse centroids with rows permuted by sigma (human row i <-> mouse row sigma(i)).
  3. C_sigma = per-gene Pearson(H[:,j], M_sigma[:,j]) across the 35 paired centroids.
  4. Select top quartile by C_sigma (C_sigma >= Q75).
  5. obs/null_sigma = UNMODIFIED joint-PCA + Procrustes obs/null on those genes under the
     sigma pairing (inner label-shuffle null, n_perm = 2000 — matching the matched-random control).
Repeat N >= 1000 -> the selection-induced obs/null distribution.

mode=labelshuffle: optional cross-check — outer draws are FULL permutations (the headline
label-shuffle family, fixed points allowed) instead of derangements; C is recomputed and
genes reselected inside that loop. Same machinery.

Outputs (in ./outputs/): per-draw CSV, the N x 35 sigma array (.npy), a summary JSON, and a
markdown report. The pre-registered reading is SURFACED, not declared (advisor adjudicates).

Run (from this dir, with the fresh clone's venv):
  ../cellwarp/.venv/bin/python selection_null.py --n 1000 --nperm 2000 --workers 8 --mode derangement
"""
from __future__ import annotations

import argparse
import contextlib
import io
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths + unmodified published code (imported, never edited)
# ---------------------------------------------------------------------------
REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "analysis" / "conserved_contribution"))
import gate_lib as G  # noqa: E402  (unmodified: load_centroids, per_gene_corr, SEED)
from cellwarp.procrustes import (  # noqa: E402  (unmodified pipeline)
    pca_reduce_centroids, procrustes_align, permutation_test,
)

# det() on near-singular kxk matrices inside the *published* pipeline emits benign
# RuntimeWarnings (np.sign handles the inf/nan sign); results match the deposit exactly.
warnings.filterwarnings("ignore", category=RuntimeWarning)

# ---------------------------------------------------------------------------
# Load substrate ONCE at import (so spawned workers get it too). Module-level.
# ---------------------------------------------------------------------------
_h, _m = G.load_centroids()                     # 35 types x 16,959 orthologs, name-aligned
TYPES = list(_h.index)
GENES = np.asarray(list(_h.columns))
H = _h.values.astype(float)                     # (35, 16959) human centroids
M = _m.values.astype(float)                     # (35, 16959) mouse centroids, true pairing = row k <-> row k
N_TYPES = H.shape[0]
IDX = [f"{i:02d}" for i in range(N_TYPES)]       # sorted integer labels -> neutralize the
                                                 # pipeline's internal index re-sort, preserving sigma-pairing
H_DF = pd.DataFrame(H, index=IDX, columns=GENES)  # human never permuted -> constant


def vectorized_C(Hmat: np.ndarray, Mmat: np.ndarray) -> np.ndarray:
    """Per-gene Pearson across the 35 paired centroids (== gate_lib.per_gene_corr, vectorized)."""
    Hc = Hmat - Hmat.mean(0, keepdims=True)
    Mc = Mmat - Mmat.mean(0, keepdims=True)
    num = (Hc * Mc).sum(0)
    den = np.sqrt((Hc ** 2).sum(0) * (Mc ** 2).sum(0))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, num / den, np.nan)


def obs_null_ratio(gene_ids, h_df: pd.DataFrame, m_df: pd.DataFrame, n_perm: int = 2000):
    """Verbatim run_gate.obs_null_ratio: the UNMODIFIED published obs/null on a gene subset."""
    cols = [g for g in gene_ids if g in h_df.columns]
    hc = h_df[cols]; mc = m_df[cols]
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        hp, mp, _, _ = pca_reduce_centroids(hc, mc, 0.95)
        res = procrustes_align(hp, mp)
        _, null = permutation_test(hp, mp, n_perm)   # inner label-shuffle null, seed=42 (published)
    return float(res.distance), float(np.median(null)), float(res.distance / np.median(null))


def select_conserved(C: np.ndarray):
    """Top-quartile selection, exactly as run_gate (Q75 on defined-C genes; C >= Q75)."""
    valid = ~np.isnan(C)
    q75 = float(np.quantile(C[valid], 0.75))
    cons_mask = valid & (C >= q75)
    return cons_mask, q75, int(valid.sum())


def run_sigma(task):
    """One null draw under permutation `perm`. Returns the conserved-set obs/null under sigma."""
    draw_id, perm, n_perm = task
    perm = np.asarray(perm)
    Mp = M[perm]                                  # sigma-paired mouse centroids
    C = vectorized_C(H, Mp)                        # C_sigma
    cons_mask, q75, n_valid = select_conserved(C)
    cons_ids = GENES[cons_mask].tolist()
    m_df = pd.DataFrame(Mp, index=IDX, columns=GENES)
    dist, null_med, ratio = obs_null_ratio(cons_ids, H_DF, m_df, n_perm=n_perm)
    return {
        "draw_id": int(draw_id), "ratio": ratio, "n_conserved": int(cons_mask.sum()),
        "n_valid": n_valid, "q75": q75, "distance": dist, "null_median": null_med,
        "n_fixed_points": int((perm == np.arange(N_TYPES)).sum()),
    }


# ---------------------------------------------------------------------------
# Draw generators
# ---------------------------------------------------------------------------
def random_derangement(rng: np.random.Generator, n: int) -> np.ndarray:
    """A uniform-ish derangement of 0..n-1 (no fixed point) by rejection sampling."""
    while True:
        p = rng.permutation(n)
        if not np.any(p == np.arange(n)):
            return p


def make_draws(mode: str, n_draws: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    if mode == "derangement":
        return np.stack([random_derangement(rng, N_TYPES) for _ in range(n_draws)])
    elif mode == "labelshuffle":          # full permutations (fixed points allowed)
        return np.stack([rng.permutation(N_TYPES) for _ in range(n_draws)])
    raise ValueError(mode)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=1000, help="number of sigma draws (>=1000 pre-registered)")
    ap.add_argument("--nperm", type=int, default=2000, help="inner label-shuffle permutations")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--seed", type=int, default=12345, help="seed for the sigma draws (wrapper-level)")
    ap.add_argument("--mode", choices=["derangement", "labelshuffle"], default="derangement")
    ap.add_argument("--out", type=str, default=str(Path(__file__).resolve().parent / "outputs"))
    args = ap.parse_args()

    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    tag = args.mode
    t0 = time.time()

    # --- Real reference values (recomputed from scratch via the same code path) ---
    C_real = vectorized_C(H, M)                                   # identity pairing = TRUE correspondence
    cons_mask_real, q75_real, n_valid = select_conserved(C_real)
    real_dist, real_nm, real_conserved = obs_null_ratio(GENES[cons_mask_real].tolist(), H_DF,
                                                         pd.DataFrame(M, index=IDX, columns=GENES), args.nperm)
    full_dist, full_nm, full_space = obs_null_ratio(GENES.tolist(), H_DF,
                                                    pd.DataFrame(M, index=IDX, columns=GENES), args.nperm)
    print(f"[real] valid_C={n_valid} Q75={q75_real:.4f} n_conserved={int(cons_mask_real.sum())} "
          f"| conserved obs/null={real_conserved:.4f} (pub 0.384) | full-space={full_space:.4f} (pub 0.522)")

    # --- Generate draws + run the null in parallel ---
    perms = make_draws(args.mode, args.n, args.seed)
    np.save(out / f"sigma_perms_{tag}.npy", perms)
    tasks = [(i, perms[i].tolist(), args.nperm) for i in range(args.n)]

    print(f"[null] mode={args.mode} N={args.n} nperm={args.nperm} workers={args.workers} ...")
    results = []
    import multiprocessing as mp
    ctx = mp.get_context("spawn")
    try:
        from tqdm import tqdm
        pbar = tqdm(total=args.n, smoothing=0.05)
    except Exception:
        pbar = None
    with ctx.Pool(processes=args.workers) as pool:
        for r in pool.imap_unordered(run_sigma, tasks, chunksize=4):
            results.append(r)
            if pbar: pbar.update(1)
            elif len(results) % 50 == 0:
                print(f"  {len(results)}/{args.n}")
    if pbar: pbar.close()

    df = pd.DataFrame(results).sort_values("draw_id").reset_index(drop=True)
    df.to_csv(out / f"sigma_null_draws_{tag}.csv", index=False)

    # --- Pre-registered reading (SURFACED, not declared) ---
    ratios = df["ratio"].to_numpy()
    mean, sd = float(ratios.mean()), float(ratios.std(ddof=1))
    p01, p05, p50 = (float(np.percentile(ratios, q)) for q in (1, 5, 50))
    rmin, rmax = float(ratios.min()), float(ratios.max())
    z = (real_conserved - mean) / sd
    frac_le_real = float((ratios <= real_conserved).mean())     # empirical percentile of the real value
    n_below_real = int((ratios <= real_conserved).sum())

    summary = {
        "design": "correspondence-permutation, selection-propagated (pre-registered 2026-06-28)",
        "mode": args.mode, "n_draws": args.n, "n_perm_inner": args.nperm,
        "seed_sigma": args.seed, "seed_inner_pipeline": int(G.SEED),
        "substrate": {"types": N_TYPES, "genes_total": int(len(GENES)), "valid_C": n_valid,
                      "n_conserved_quartile": int(cons_mask_real.sum()), "Q75_real": q75_real},
        "real": {"conserved_obs_null": real_conserved, "full_space_obs_null": full_space,
                 "published_conserved": 0.384, "published_full_space": 0.522,
                 "conserved_distance": real_dist, "conserved_null_median": real_nm},
        "sigma_null": {"mean": mean, "sd": sd, "p01": p01, "p05": p05, "median": p50,
                       "min": rmin, "max": rmax,
                       "n_conserved_per_draw": {"min": int(df.n_conserved.min()),
                                                "max": int(df.n_conserved.max())}},
        "real_position": {"z": z, "empirical_percentile": frac_le_real,
                          "n_draws_at_or_below_real": n_below_real},
        "preregistered_criteria_SURFACED_not_declared": {
            "PASS_if": "real conserved obs/null < sigma-null 1st percentile (equivalently z <= -3)",
            "FAIL_if": "real falls within / near the sigma-null distribution",
            "real_below_1st_percentile": bool(real_conserved < p01),
            "z_le_minus3": bool(z <= -3.0),
            "diagnostic_sigma_null_center_vs_full_space_0p522": mean - full_space,
            "note": "Engineer surfaces; advisor adjudicates the margin. No verdict declared here.",
        },
        "runtime_sec": round(time.time() - t0, 1),
    }
    with open(out / f"selection_null_summary_{tag}.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("\n" + "=" * 70)
    print(f"SELECTION-AWARE NULL ({args.mode}) — N={args.n}, nperm={args.nperm}")
    print("=" * 70)
    print(f"  real conserved obs/null      = {real_conserved:.4f}   (published 0.384)")
    print(f"  full-space obs/null          = {full_space:.4f}   (published 0.522)")
    print(f"  sigma-null obs/null mean+-sd = {mean:.4f} +/- {sd:.4f}")
    print(f"  sigma-null 1st percentile    = {p01:.4f}   (5th = {p05:.4f}, median = {p50:.4f})")
    print(f"  sigma-null range             = [{rmin:.4f}, {rmax:.4f}]")
    print(f"  real position: z = {z:.2f} | empirical percentile = {frac_le_real*100:.2f}% "
          f"({n_below_real}/{args.n} draws <= real)")
    print(f"  [pre-registered, SURFACED] real < sigma-null 1st pct? {real_conserved < p01} ; "
          f"z <= -3? {z <= -3.0}")
    print(f"  [diagnostic] sigma-null center vs full-space 0.522: {mean:.4f} vs {full_space:.4f} "
          f"(delta {mean-full_space:+.4f})")
    print(f"  runtime {summary['runtime_sec']}s -> outputs/ ({tag})")
    print("  Advisor adjudicates the margin; this run surfaces, does not declare.")


if __name__ == "__main__":
    main()
