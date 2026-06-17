"""
Post-hoc robustness (NOT pre-registered): is the conserved-contribution signal driven by cross-species
conservation, or merely by cell-type SPECIFICITY (peaked genes correlate more easily,
and identity TFs/markers are specific)?

Test: re-run the 3a TF positive control and 3b enrichment against a background matched
jointly on BOTH mean expression AND specificity (Tau). If identity TFs / markers still
beat a specificity+expression-matched background, the signal is conservation, not
specificity per se. Plus the canonical Hartigan dip test for check 1.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import diptest

import gate_lib as G

rng = np.random.default_rng(G.SEED)
HERE = Path(__file__).resolve().parent

df, h_cent, m_cent = G.build_gene_table()
valid = df.dropna(subset=["C_pearson"]).copy().reset_index(drop=True)
C = valid["C_pearson"].values
me = valid["mean_expression"].values
n = len(valid)

# Tau specificity on human centroids (log1p, non-negative). Tau in [0,1].
Hc = h_cent[valid["gene_id"]].values  # 35 x n
def tau(col):
    x = col.copy()
    mx = x.max()
    if mx <= 0:
        return np.nan
    return np.sum(1 - x / mx) / (len(x) - 1)
Tau = np.array([tau(Hc[:, j]) for j in range(n)])
valid["tau"] = Tau

rho_C_tau = stats.spearmanr(C, Tau)[0]
rho_C_expr = stats.spearmanr(C, me)[0]
print(f"rho(C, Tau specificity) = {rho_C_tau:.3f}")
print(f"rho(C, mean_expr)       = {rho_C_expr:.3f}")
print(f"rho(Tau, mean_expr)     = {stats.spearmanr(Tau, me)[0]:.3f}")

# dip test (canonical multimodality)
dstat, dip_p = diptest.diptest(C)
print(f"Hartigan dip test on C: D={dstat:.4f}, p={dip_p:.4f}  "
      f"({'multimodal' if dip_p < 0.05 else 'unimodal / continuum'})")

# ---- joint expression x Tau matched background ----
def jbin(v, k):
    r = stats.rankdata(v, method="ordinal")
    return np.minimum((r - 1) * k // len(r), k - 1)
be = jbin(me, 10); bt = jbin(Tau, 10)
jbins = be * 10 + bt  # 100 joint bins

all_idx = np.arange(n)
C_rank = stats.rankdata(C) / n
q75 = np.quantile(C, 0.75)
cons_idx = all_idx[C >= q75]

sym_to_pos = {s: i for i, s in enumerate(valid["symbol"].values)}
def syms_to_pos(symbols):
    return np.array(sorted({sym_to_pos[s] for s in symbols if s in sym_to_pos}), dtype=int)

tf_pc_pos = syms_to_pos(G.POSITIVE_CONTROL_TFS)
tf_all_pos = syms_to_pos(pd.read_csv(G.TF_ACTIVITY, index_col=0).columns.tolist())
cm_pos = syms_to_pos(pd.read_csv(G.CELLMARKER_H)["gene_symbol"].dropna().unique().tolist())

def matched_pvalue(target_pos, stat_fn, ndraw=1000):
    obs = stat_fn(target_pos)
    draws = G.matched_draws(target_pos, jbins, ndraw, rng)
    null = np.array([stat_fn(d) for d in draws])
    return obs, float(np.median(null)), G.emp_p_greater(obs, null)

# 3a under joint matching
median_rank = lambda pos: float(np.median(C_rank[pos]))
o, nm, p = matched_pvalue(tf_pc_pos, median_rank, 1000)
tau_above_med = float(np.mean(Tau[tf_pc_pos] > np.median(Tau)))
print(f"\n[3a joint expr+Tau matched] identity-TF median C-percentile={o:.3f} "
      f"(matched {nm:.3f}) p={p:.4f}  | frac TFs above median Tau={tau_above_med:.2f}")
res3a = dict(median_Crank=o, matched_median=nm, p=float(p), tf_frac_high_tau=tau_above_med)

# 3b under joint matching
def overlap(set_pos, ref_pos, ndraw=1000):
    ref = set(ref_pos.tolist())
    obs = int(sum(pp in ref for pp in set_pos))
    null = np.array([int(sum(pp in ref for pp in d))
                     for d in G.matched_draws(set_pos, jbins, ndraw, rng)])
    return obs, float(np.mean(null)), G.emp_p_greater(obs, null)
otf, etf, ptf = overlap(cons_idx, tf_all_pos)
ocm, ecm, pcm = overlap(cons_idx, cm_pos)
print(f"[3b joint matched] TF: conserved {otf} vs matched {etf:.1f} ({otf/etf:.2f}x) p={ptf:.4f}")
print(f"[3b joint matched] CellMarker: conserved {ocm} vs matched {ecm:.1f} ({ocm/ecm:.2f}x) p={pcm:.4f}")

out = dict(
    rho_C_tau=float(rho_C_tau), rho_C_expr=float(rho_C_expr),
    rho_tau_expr=float(stats.spearmanr(Tau, me)[0]),
    dip_D=float(dstat), dip_p=float(dip_p), dip_multimodal=bool(dip_p < 0.05),
    joint_matched_3a=res3a,
    joint_matched_3b=dict(
        TF=dict(observed=otf, matched=etf, fold=float(otf/etf), p=float(ptf)),
        CellMarker=dict(observed=ocm, matched=ecm, fold=float(ocm/ecm), p=float(pcm)),
    ),
)
with open(HERE / "robustness_results.json", "w") as f:
    json.dump(out, f, indent=2)
print(f"\nSaved {HERE/'robustness_results.json'}")
