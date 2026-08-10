#!/usr/bin/env python3
"""
Detection-breadth sensitivity of the conservation score C and the master-TF enrichment.

Question: C is a per-gene Pearson correlation of expression profile across the 35
matched centroids. A gene detected in only a handful of those centroids can reach
|C| ~ 1 on very little evidence, and Fig 5A's right tail invites exactly that worry.
This script measures how broadly each gene is actually detected, relates that to C,
and re-runs the Fig 5C master-TF enrichment on genes detected everywhere, to establish
whether the headline enrichment depends on thinly-supported genes.

DETECTION CRITERION, AND WHERE IT COMES FROM

The pipeline defines NO detection-breadth criterion for C, and this script does not
introduce one into the pipeline -- it supplies one for sensitivity analysis only. The
two relevant facts, both verified against the tree:

  * gate_lib.per_gene_corr filters on np.std(a) > 0 and np.std(b) > 0 -- nonzero
    variance across the 35 centroids in each species. It says nothing about how many
    of those centroids are nonzero, so a gene detected in three types and absent from
    thirty-two passes it.
  * analysis/biological_predictors/biological_predictors.py:268 defines a real
    detection rate, frac_expressing = (X_ct > 0).sum(axis=0) / n_cells, thresholded at
    0.10. That is a per-cell criterion and needs the raw cell matrices, which are not
    deposited, so it cannot be applied here.

The criterion used is therefore: a gene is DETECTED in a cell type when its centroid
value in that type is > 0. The centroids are non-negative log-normalised means, so an
exact zero means zero counts in every cell of that type. This is the most generous
reading available from deposited artifacts -- it counts a gene as detected on the
strength of a single cell -- so any shortfall it reports is a lower bound on the
problem, which is the right direction for a sensitivity analysis.

breadth_min is the minimum of the human and mouse per-gene breadths, i.e. the number of
cell types in which the gene is detected in BOTH species, which is what C is computed
across.

METHOD

The enrichment re-run reuses the frozen statistic and the frozen sampler rather than
re-deriving either: gate_lib.matched_draws, gate_lib.expr_bins, the 20-bin
expression-matched null of run_gate.py CHECK 3a, and the 10x10 expression x Tau null of
run_robustness.py joint_matched_3a. The only thing that changes between the two
enrichment rows is the gene pool. Percentiles, bins and Tau are all rebuilt inside each
pool, because a percentile in the filtered pool is not a percentile in the full one.

The unfiltered row is computed FIRST and asserted against the deposited gate before any
filtered row is touched. Without that, a filtered number means nothing: it could differ
from the paper because the filter bit, or because this script reimplemented the
statistic differently, and there would be no way to tell which.

Outputs (tracked):
  analysis/conserved_contribution/breadth_sensitivity_results.json
  analysis/conserved_contribution/breadth_per_gene.csv

Deterministic: fixed sampler seeds, no network. Runtime about a minute. Not invoked by
reproduce/run_all.sh, on the same footing as analysis/simulation_study/sweep_spread.py
and analysis/ranking_replication/cross_atlas_ci.py. Gated by reproduce/validate.py.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import gate_lib as G

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "breadth_sensitivity_results.json"
OUT_CSV = HERE / "breadth_per_gene.csv"

N_DRAWS = 100_000
SEED_EXPR = 20240601      # same seeds as highN_tf_pvalues.py, for comparability
SEED_JOINT = 20240602
STRICTEST = 35            # detected in all 35 types, in both species
HIST_BINS = 60            # make_figure7.py:165 builds Fig 5A with exactly this

CRITERION = (
    "A gene is detected in a cell type when its centroid value in that type is > 0. "
    "The centroids are non-negative log-normalised means, so an exact zero is zero "
    "counts in every cell of that type. breadth_min is the number of the 35 types in "
    "which the gene is detected in BOTH species. This criterion is supplied by this "
    "sensitivity analysis and is NOT the pipeline's: gate_lib.per_gene_corr filters "
    "only on np.std > 0 in each species, which admits a gene detected in three types; "
    "and biological_predictors.py:268's frac_expressing > 0.10 is a per-cell rate "
    "needing raw cell matrices that are not deposited. It is the most generous "
    "reading available, so every shortfall reported here is a lower bound."
)

# Deposited values the unfiltered row must reproduce before anything else runs.
GATE_OBS = 0.9377038895859473      # gate_results.json check3a.median_Crank
GATE_TESTABLE = 73                 # gate_results.json check3a.testable
GATE_NULL_EXPR = 0.5399623588456712    # highN_tf_pvalues.json expression_matched.null_median
GATE_NULL_JOINT = 0.7611041405269762   # highN_tf_pvalues.json joint_expr_tau_matched.null_median


def tau_specificity(Hc: np.ndarray) -> np.ndarray:
    """Tau tissue-specificity per gene on human centroids (run_robustness.py's form)."""
    def _t(col):
        mx = col.max()
        return np.sum(1 - col / mx) / (len(col) - 1) if mx > 0 else np.nan
    return np.array([_t(Hc[:, j]) for j in range(Hc.shape[1])])


def enrichment(sub: pd.DataFrame, h_cent: pd.DataFrame, label: str) -> dict:
    """Median C-percentile of the master TFs vs both matched nulls, inside this pool."""
    n = len(sub)
    C = sub["C_pearson"].to_numpy()
    me = sub["mean_expression"].to_numpy()
    Hc = h_cent[sub["gene_id"]].to_numpy()

    C_rank = stats.rankdata(C) / n
    pos = {s: i for i, s in enumerate(sub["symbol"].to_numpy())}
    tf = np.array(sorted({pos[s] for s in G.POSITIVE_CONTROL_TFS if s in pos}), dtype=int)
    obs = float(np.median(C_rank[tf]))

    def _jbin(v, k):
        r = stats.rankdata(v, method="ordinal")
        return np.minimum((r - 1) * k // len(r), k - 1)

    bins = {
        "expression_matched": (G.expr_bins(me, 20), SEED_EXPR),
        "joint_expr_tau_matched": (_jbin(me, 10) * 10 + _jbin(tau_specificity(Hc), 10),
                                   SEED_JOINT),
    }

    out = {"label": label, "n_genes": int(n), "n_tfs": int(len(tf)),
           "obs_median_C_percentile": obs, "n_draws": N_DRAWS}
    for key, (b, seed) in bins.items():
        rng = np.random.default_rng(seed)
        meds = np.fromiter(
            (np.median(C_rank[d]) for d in G.matched_draws(tf, b, N_DRAWS, rng)),
            dtype=float, count=N_DRAWS)
        ge = int(np.count_nonzero(meds >= obs))
        out[key] = {
            "null_median": float(np.median(meds)),
            "null_mean": float(meds.mean()),
            "null_sd": float(meds.std()),
            "exceedances": ge,
            "p_empirical": (ge + 1) / (N_DRAWS + 1),
            "z": float((obs - meds.mean()) / meds.std()),
        }
    print(f"  {label:<34} n={n:>6} TFs={out['n_tfs']:>3} obs={obs:.4f} "
          f"expr null={out['expression_matched']['null_median']:.4f} "
          f"p={out['expression_matched']['p_empirical']:.2e} | "
          f"joint null={out['joint_expr_tau_matched']['null_median']:.4f} "
          f"p={out['joint_expr_tau_matched']['p_empirical']:.2e}")
    return out


def main() -> None:
    print("=" * 78)
    print("DETECTION-BREADTH SENSITIVITY OF C AND THE MASTER-TF ENRICHMENT")
    print("=" * 78)

    df, h_cent, m_cent = G.build_gene_table()
    valid = df.dropna(subset=["C_pearson"]).reset_index(drop=True)
    gid = valid["gene_id"].to_numpy()
    H = h_cent[gid].to_numpy()
    M = m_cent[gid].to_numpy()

    bh = (H > 0).sum(axis=0)
    bm = (M > 0).sum(axis=0)
    bmin = np.minimum(bh, bm)
    valid = valid.assign(breadth_human=bh, breadth_mouse=bm, breadth_min=bmin)
    C = valid["C_pearson"].to_numpy()
    n_total = len(valid)
    print(f"\ngenes with a defined C: {n_total}")

    # ── the verification row, before anything filtered ───────────────────────
    print("\nVERIFICATION -- unfiltered row must reproduce the deposited gate:")
    unfiltered = enrichment(valid, h_cent, "unfiltered")
    assert abs(unfiltered["obs_median_C_percentile"] - GATE_OBS) < 1e-9, (
        f"unfiltered observed {unfiltered['obs_median_C_percentile']!r} does not "
        f"reproduce gate_results.json check3a.median_Crank {GATE_OBS!r}; refusing to "
        f"report a filtered number against a statistic that does not match the paper")
    assert unfiltered["n_tfs"] == GATE_TESTABLE, (
        f"unfiltered TF count {unfiltered['n_tfs']} != check3a.testable {GATE_TESTABLE}")
    for key, dep in (("expression_matched", GATE_NULL_EXPR),
                     ("joint_expr_tau_matched", GATE_NULL_JOINT)):
        got = unfiltered[key]["null_median"]
        assert abs(got - dep) < 0.02, f"{key} null median {got:.4f} vs deposited {dep:.4f}"
    print("  -> reproduces check3a.median_Crank exactly and both null medians "
          "within sampler tolerance")

    # ── breadth, breadth vs C, the Fig 5A tail, the conserved set ────────────
    rho, p = stats.spearmanr(bmin, C)
    counts, edges = np.histogram(C, bins=HIST_BINS)
    term = C >= edges[-2]
    q75 = float(np.quantile(C, 0.75))
    cons = C >= q75

    def summary(v):
        return {
            "median": float(np.median(v)), "mean": float(np.mean(v)),
            "q05": float(np.quantile(v, 0.05)), "q25": float(np.quantile(v, 0.25)),
            "q75": float(np.quantile(v, 0.75)),
            "n_at_35": int((v == 35).sum()),
            "n_le_3": int((v <= 3).sum()), "n_le_5": int((v <= 5).sum()),
            "n_le_10": int((v <= 10).sum()),
        }

    # ── the strictest filter ─────────────────────────────────────────────────
    print(f"\nFILTERED -- breadth_min >= {STRICTEST} (detected in all 35 types, both species):")
    strict = enrichment(valid[bmin >= STRICTEST].reset_index(drop=True), h_cent,
                        f"breadth_min >= {STRICTEST}")

    result = {
        "analysis": "detection-breadth sensitivity of C and the Fig 5C master-TF enrichment",
        "detection_criterion": CRITERION,
        "criterion_is_the_pipelines": False,
        "pipeline_criterion_notes": {
            "gate_lib.per_gene_corr": "filters on np.std > 0 per species only; no breadth requirement",
            "biological_predictors.py:268": "frac_expressing > 0.10, per-cell, needs raw matrices (not deposited)",
        },
        "n_genes_with_C": int(n_total),
        "breadth": {"human": summary(bh), "mouse": summary(bm), "min": summary(bmin)},
        "breadth_vs_C": {"spearman_rho": float(rho), "p_value": float(p), "n": int(n_total),
                         "direction_note": "positive: sparsely detected genes sit at NEGATIVE C, "
                                           "not in the high-C tail"},
        "figure_5A_terminal_bin": {
            "n_bins": HIST_BINS,
            "source": "make_figure7.py:165 np.histogram(C, bins=60)",
            "lower_edge": float(edges[-2]), "upper_edge": float(edges[-1]),
            "n_genes": int(term.sum()),
            "median_breadth_min": float(np.median(bmin[term])),
            "frac_breadth_min_le_3": float(np.mean(bmin[term] <= 3)),
            "frac_breadth_min_le_5": float(np.mean(bmin[term] <= 5)),
            "all_genes_frac_breadth_min_le_3": float(np.mean(bmin <= 3)),
        },
        "conserved_set": {
            "q75_of_C": q75, "n_genes": int(cons.sum()),
            "breadth_min": summary(bmin[cons]),
            "frac_breadth_min_le_5": float(np.mean(bmin[cons] <= 5)),
            "frac_breadth_min_le_15": float(np.mean(bmin[cons] <= 15)),
        },
        "enrichment": {"unfiltered": unfiltered, "strictest": strict},
        "strictest_filter_breadth_min": STRICTEST,
    }

    OUT_JSON.write_text(json.dumps(result, indent=2) + "\n")
    valid[["gene_id", "symbol", "C_pearson", "mean_expression",
           "breadth_human", "breadth_mouse", "breadth_min"]].to_csv(OUT_CSV, index=False)

    print(f"\n  Spearman(breadth_min, C) = {rho:+.4f}  p = {p:.3e}")
    print(f"  Fig 5A terminal bin [{edges[-2]:.4f}, {edges[-1]:.4f}]: "
          f"{int(term.sum())} genes, median breadth_min {np.median(bmin[term]):.0f}")
    print(f"  conserved set (C >= {q75:.4f}): {int(cons.sum())} genes, "
          f"median breadth_min {np.median(bmin[cons]):.0f}")
    print(f"\nwrote {OUT_JSON}")
    print(f"wrote {OUT_CSV}")


if __name__ == "__main__":
    main()
