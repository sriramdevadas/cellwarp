#!/usr/bin/env python3
"""
Out-of-sample retention of the conserved-contribution gene set, measured on C alone.

Splits donors of both species into halves, computes C on one half, selects the top
quartile on that half, and reads the selected genes' mean C on the other half, so
selection and evaluation share no donor.

Produces the three C values and the retained fraction reported in S1 Text, "The
objection is answered instead out of sample": mean C of 0.663 on the other half
against 0.730 for genes selected on the evaluation half and 0.318 for all genes,
so roughly 84% of the apparent gain survives.

Reads the repository copy of docs/supplementary_materials/table_S11_gene_conservation.csv
rather than a working copy, so a divergence between the two would surface here. All
inputs are tracked.

Only the two path constants differ from the version that produced those values; the
code is unchanged.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
S11 = REPO / "docs/supplementary_materials/table_S11_gene_conservation.csv"
OUT = HERE
SEED = 42

R = {}
df = pd.read_csv(S11)
print(f"S11 rows: {len(df)}   columns: {list(df.columns)}")

A = df["donor_split_C_A"]
B = df["donor_split_C_B"]

# ---- (1) corpus sizes; assert non-empty BEFORE any statistic (note 11) -------
nA, nB = int(A.notna().sum()), int(B.notna().sum())
paired = A.notna() & B.notna()
nP = int(paired.sum())
print("\n=== (1) corpus ===")
print(f"  n non-null C_A        : {nA}")
print(f"  n non-null C_B        : {nB}")
print(f"  n non-null in BOTH    : {nP}")
print(f"  n rows total          : {len(df)}")
assert nP > 0, "paired corpus is EMPTY -- every result below would be void"
R["corpus"] = dict(n_rows=len(df), n_A=nA, n_B=nB, n_paired=nP)

a = A[paired].to_numpy(float)
b = B[paired].to_numpy(float)

# ---- (2) Spearman ------------------------------------------------------------
rho, p = stats.spearmanr(a, b)
pear = float(np.corrcoef(a, b)[0, 1])
print("\n=== (2) cross-half agreement (split 0) ===")
print(f"  Spearman rho = {rho:.6f}   p = {p:.3e}   n = {nP}")
print(f"  Pearson  r   = {pear:.6f}")
print(f"  S1 reports median rho = 0.80 ACROSS 100 SPLITS (not split 0); "
      f"donor_stability_results.json cross_half_C_spearman_median = 0.7954227089705543")
R["spearman"] = dict(rho=float(rho), p=float(p), pearson=pear, n=nP)

# ---- (3) quartile thresholds -------------------------------------------------
qA, qB = float(np.quantile(a, 0.75)), float(np.quantile(b, 0.75))
topA = a >= qA
topB = b >= qB
print("\n=== (3) top-quartile thresholds on the paired subset ===")
print(f"  Q75(C_A) = {qA:.6f}   |topA| = {int(topA.sum())}")
print(f"  Q75(C_B) = {qB:.6f}   |topB| = {int(topB.sum())}")
R["quartiles"] = dict(q75_A=qA, q75_B=qB, n_topA=int(topA.sum()), n_topB=int(topB.sum()))

# ---- (4) Jaccard -------------------------------------------------------------
inter = int((topA & topB).sum())
union = int((topA | topB).sum())
jac = inter / union
print("\n=== (4) Jaccard of the two top quartiles ===")
print(f"  |A n B| = {inter}   |A u B| = {union}   Jaccard = {jac:.6f}")
print(f"  S1 reports median Jaccard = 0.58 ACROSS 100 SPLITS; "
      f"donor_stability_results.json conserved_jaccard_median = 0.5797518560552324")
R["jaccard"] = dict(intersection=inter, union=union, jaccard=jac)


def sem(x):
    return float(np.std(x, ddof=1) / np.sqrt(len(x)))


def winners_curse(sel_vec, eval_vec, sel_name, eval_name):
    """Three means on eval_vec: all / selected-by-eval (in-sample) / selected-by-sel (OOS)."""
    q_sel = np.quantile(sel_vec, 0.75)
    q_eval = np.quantile(eval_vec, 0.75)
    m_all = eval_vec
    m_in = eval_vec[eval_vec >= q_eval]
    m_oos = eval_vec[sel_vec >= q_sel]
    out = dict(
        eval_on=eval_name, selected_on=sel_name,
        mean_all=float(m_all.mean()), se_all=sem(m_all), n_all=len(m_all),
        mean_in_sample=float(m_in.mean()), se_in_sample=sem(m_in), n_in_sample=len(m_in),
        mean_out_of_sample=float(m_oos.mean()), se_out_of_sample=sem(m_oos), n_out_of_sample=len(m_oos),
    )
    out["earned_retention__oos_minus_all"] = out["mean_out_of_sample"] - out["mean_all"]
    out["manufactured__in_minus_oos"] = out["mean_in_sample"] - out["mean_out_of_sample"]
    return out


def show(d, title):
    print(f"\n--- {title} ---")
    print(f"  mean C_{d['eval_on']} over ALL paired genes          = {d['mean_all']:.6f} "
          f"+/- {d['se_all']:.6f}   (n={d['n_all']})")
    print(f"  mean C_{d['eval_on']} over top-quartile by C_{d['eval_on']} (IN-sample)  = {d['mean_in_sample']:.6f} "
          f"+/- {d['se_in_sample']:.6f}   (n={d['n_in_sample']})")
    print(f"  mean C_{d['eval_on']} over top-quartile by C_{d['selected_on']} (OUT-of-sample) = {d['mean_out_of_sample']:.6f} "
          f"+/- {d['se_out_of_sample']:.6f}   (n={d['n_out_of_sample']})")
    print(f"  EARNED retention   (OOS - all) = {d['earned_retention__oos_minus_all']:+.6f}")
    print(f"  MANUFACTURED       (IN  - OOS) = {d['manufactured__in_minus_oos']:+.6f}")


# ---- (5) select on A, evaluate on B -----------------------------------------
print("\n=== (5) winner's curse: select on A, evaluate on B ===")
d_AB = winners_curse(a, b, "A", "B")
show(d_AB, "evaluate on half B")
R["select_A_eval_B"] = d_AB

# ---- (6) exchanged ----------------------------------------------------------
print("\n=== (6) exchanged: select on B, evaluate on A ===")
d_BA = winners_curse(b, a, "B", "A")
show(d_BA, "evaluate on half A")
R["select_B_eval_A"] = d_BA

# ---- (7) NEGATIVE CONTROL ---------------------------------------------------
# random quartile of the same size must collapse the OOS mean onto the all-genes mean
print("\n=== (7) NEGATIVE CONTROL: random quartile instead of top-C_A ===")
rng = np.random.default_rng(SEED)
k = int(topA.sum())
draws = []
for _ in range(1000):
    idx = rng.choice(nP, size=k, replace=False)
    draws.append(float(b[idx].mean()))
draws = np.array(draws)
print(f"  k (same size as top quartile) = {k}, 1000 draws, seed {SEED}")
print(f"  mean C_B over ALL paired genes            = {b.mean():.6f}")
print(f"  random-quartile mean C_B: mean over draws = {draws.mean():.6f}")
print(f"                            sd  over draws  = {draws.std(ddof=1):.6f}")
print(f"                            2.5 - 97.5 pct  = [{np.percentile(draws,2.5):.6f}, {np.percentile(draws,97.5):.6f}]")
print(f"  |random mean - all-genes mean|            = {abs(draws.mean()-b.mean()):.3e}")
print(f"  observed top-C_A OOS mean                 = {d_AB['mean_out_of_sample']:.6f}")
z = (d_AB["mean_out_of_sample"] - draws.mean()) / draws.std(ddof=1)
print(f"  z of top-C_A OOS mean vs random-quartile distribution = {z:.2f}")
R["negative_control"] = dict(
    k=k, n_draws=1000, seed=SEED,
    mean_all=float(b.mean()), random_mean=float(draws.mean()), random_sd=float(draws.std(ddof=1)),
    random_p2_5=float(np.percentile(draws, 2.5)), random_p97_5=float(np.percentile(draws, 97.5)),
    abs_gap_random_vs_all=float(abs(draws.mean() - b.mean())),
    observed_oos_mean=d_AB["mean_out_of_sample"], z_vs_random=float(z),
)

(OUT / "block1_form_a_results.json").write_text(json.dumps(R, indent=2))
print(f"\nwrote {OUT/'block1_form_a_results.json'}")
