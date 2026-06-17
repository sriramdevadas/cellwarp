"""
High-resolution empirical p-values for the master-TF positive control (cb M3).

The pre-registered gate (run_gate.py CHECK 3a) and the joint expr x Tau
robustness check (run_robustness.py joint_matched_3a) each ran N = 1,000 matched
draws, so both report the empirical floor p = 1/(1000+1) = 0.000999. The medians
they report (master-TF 0.94; expr-matched background 0.54; expr+Tau-matched
background 0.76) are deterministic and unchanged here; only the p-resolution
improves.

This script re-runs BOTH matched-background permutation tests for the SAME
median-C-percentile statistic at high N, reusing the FROZEN sampler
(gate_lib.matched_draws) and the SAME bin constructions used by the two
producers:
  - expression-matched: 20 equal-frequency mean-expression bins (run_gate.py)
  - joint expr x Tau-matched: 10 x 10 mean-expression x Tau bins (run_robustness.py)

It changes nothing in the frozen gate / robustness outputs; it writes only
highN_tf_pvalues.json (additive). Medians are asserted unchanged against the
deposited values as a guard.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

import gate_lib as G

HERE = Path(__file__).resolve().parent
N = int(sys.argv[1]) if len(sys.argv) > 1 else 1_000_000
CHUNK = 25_000

# ── frozen inputs (identical to run_gate.py / run_robustness.py) ──────────────
df, h_cent, m_cent = G.build_gene_table()
valid = df.dropna(subset=["C_pearson"]).reset_index(drop=True)
C = valid["C_pearson"].values
me = valid["mean_expression"].values
n = len(valid)

# Tau specificity on human centroids (exactly as run_robustness.py / the gate)
Hc = h_cent[valid["gene_id"]].values  # 35 x n
def _tau(col):
    mx = col.max()
    return np.sum(1 - col / mx) / (len(col) - 1) if mx > 0 else np.nan
Tau = np.array([_tau(Hc[:, j]) for j in range(n)])

C_rank = stats.rankdata(C) / n
sym_to_pos = {s: i for i, s in enumerate(valid["symbol"].values)}
tf_pos = np.array(sorted({sym_to_pos[s] for s in G.POSITIVE_CONTROL_TFS
                          if s in sym_to_pos}), dtype=int)
obs = float(np.median(C_rank[tf_pos]))

# expression-matched bins (run_gate.py: G.expr_bins(me, 20))
ebins = G.expr_bins(me, 20)
# joint expr x Tau bins (run_robustness.py)
def _jbin(v, k):
    r = stats.rankdata(v, method="ordinal")
    return np.minimum((r - 1) * k // len(r), k - 1)
jbins = _jbin(me, 10) * 10 + _jbin(Tau, 10)

median_rank = lambda pos: float(np.median(C_rank[pos]))


def high_n(bins, seed, label):
    """Empirical p that a matched-null median C-percentile >= observed, at N draws."""
    rng = np.random.default_rng(seed)
    ge = 0
    done = 0
    sample = []  # capped sample of null medians for mean/sd/median reporting
    while done < N:
        k = min(CHUNK, N - done)
        draws = G.matched_draws(tf_pos, bins, k, rng)
        meds = np.fromiter((np.median(C_rank[d]) for d in draws), float, k)
        ge += int(np.count_nonzero(meds >= obs))
        if len(sample) < 200_000:
            sample.append(meds)
        done += k
    sample = np.concatenate(sample)
    nmean, nsd, nmed = float(sample.mean()), float(sample.std()), float(np.median(sample))
    p_emp = (ge + 1) / (N + 1)
    z = (obs - nmean) / nsd if nsd > 0 else float("inf")
    # normal-approximation p (context for the empirical floor when ge == 0)
    p_norm = float(stats.norm.sf(z))
    print(f"[{label}] N={N:,}  obs_median={obs:.4f}  null_median={nmed:.4f}  "
          f"null_mean={nmean:.4f}  null_sd={nsd:.4f}")
    print(f"    exceedances={ge}  p_empirical={p_emp:.3e}  z={z:.2f}  p_normal_approx={p_norm:.2e}")
    return dict(n_draws=N, obs_median=obs, null_median=nmed, null_mean=nmean,
                null_sd=nsd, exceedances=ge, p_empirical=p_emp, z=float(z),
                p_normal_approx=p_norm)


# Guard: medians must match the deposited values (deterministic).
assert abs(obs - 0.9377038895859473) < 1e-9, f"obs median drifted: {obs}"

res_expr = high_n(ebins, 20240601, "expr-matched (20-bin)")
res_joint = high_n(jbins, 20240602, "joint expr+Tau-matched (10x10-bin)")

# Median guards against the deposited matched medians (resolution only changes p).
assert abs(res_expr["null_median"] - 0.54) < 0.02, res_expr["null_median"]
assert abs(res_joint["null_median"] - 0.76) < 0.02, res_joint["null_median"]

out = dict(
    statistic="median C-percentile of the 73 master TFs vs matched background",
    n_testable_tfs=int(len(tf_pos)),
    expression_matched=res_expr,
    joint_expr_tau_matched=res_joint,
    note=("Medians unchanged from the frozen gate (0.94 / 0.54 / 0.76); only the "
          "empirical p-resolution improves over the N=1000 floor (0.000999). "
          "p_empirical is (exceedances+1)/(N+1); p_normal_approx is the normal-tail "
          "context for the standardized effect z."),
)
with open(HERE / "highN_tf_pvalues.json", "w") as f:
    json.dump(out, f, indent=2)
print("\nSaved", HERE / "highN_tf_pvalues.json")
