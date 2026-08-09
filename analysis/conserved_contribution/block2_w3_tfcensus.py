#!/usr/bin/env python3
"""
Matched-null enrichment of conserved contribution against a full TF census.

The 73 master transcription factors were mapped to the 35 matched lineages, so they
were selected against the very types that define C. The expression- and
specificity-matched nulls control those two properties but neither controls gene-set
composition. This repeats the enrichment against a full human TF census instead of
the curated master set, so the comparison is against all TFs rather than against a
set chosen alongside C.

Backs the master-TF enrichment reported in Results section 5 and drawn in Fig 5C.

NOT TRACKED: reads a third-party human TF census CSV (DatabaseExtract_v_1.01.csv)
that is not deposited in this repository and must be obtained from its own source.
highN_tf_pvalues.py and highN_tf_pvalues.json in this directory are tracked.

Only the two path constants differ from the version that produced those values; the
code is unchanged.
"""
import ast
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
CC = REPO / "analysis/conserved_contribution"
W3 = HERE
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(CC))
import gate_lib as G  # noqa: E402

N_DRAWS = 50_000            # published run used 1e6; stated, not hidden
LAMBERT = W3 / "DatabaseExtract_v_1.01.csv"

# ---------------------------------------------------------------------------
# Exec the published setup, stopping before the expensive part and the write
# ---------------------------------------------------------------------------
src = (CC / "highN_tf_pvalues.py").read_text()
tree = ast.parse(src)
keep = []
for node in tree.body:
    if isinstance(node, ast.Assign) and any(
            isinstance(t, ast.Name) and t.id == "res_expr" for t in node.targets):
        break
    keep.append(node)
ns = {"__name__": "_frozen_setup", "__file__": str(CC / "highN_tf_pvalues.py")}
sys.argv = ["highN_tf_pvalues.py"]          # so N falls to its 1e6 default, rebound below
exec(compile(ast.Module(body=keep, type_ignores=[]), "highN_tf_pvalues.py", "exec"), ns)
print("published setup exec'd; the deposited-median assert passed "
      "(obs == 0.9377038895859473)")

valid = ns["valid"]; C_rank = ns["C_rank"]; me = ns["me"]; Tau = ns["Tau"]
ebins = ns["ebins"]; jbins = ns["jbins"]; tf_pos = ns["tf_pos"]; obs = ns["obs"]
high_n = ns["high_n"]; n = ns["n"]
ns["N"] = N_DRAWS                            # high_n reads the module global N
print(f"valid genes {n}   master TFs in space {len(tf_pos)}   obs median {obs:.6f}")

R = {"n_draws": N_DRAWS, "n_valid": int(n), "obs_median_73TF": float(obs)}

# ---------------------------------------------------------------------------
# FAITHFULNESS GATE
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print(f"FAITHFULNESS GATE -- reproduce 0.94 / 0.54 / 0.76 (N={N_DRAWS:,})")
print("=" * 78)
g_expr = high_n(ebins, 20240601, "expr-matched (20-bin)")
g_joint = high_n(jbins, 20240602, "joint expr+Tau (10x10)")
pub = json.load(open(CC / "highN_tf_pvalues.json"))
print(f"\n  {'quantity':38s} {'here':>10s} {'deposited':>10s} {'diff':>10s}")
rows = [("obs median C-percentile (73 TFs)", obs, pub["expression_matched"]["obs_median"]),
        ("expr-matched null median", g_expr["null_median"], pub["expression_matched"]["null_median"]),
        ("joint expr+Tau null median", g_joint["null_median"], pub["joint_expr_tau_matched"]["null_median"])]
ok = True
for nm, a, b in rows:
    print(f"  {nm:38s} {a:10.4f} {b:10.4f} {a-b:+10.4f}")
    if abs(a - b) > 0.02:
        ok = False
print(f"  -> GATE {'PASS' if ok else 'FAIL'} (tolerance 0.02 on the null medians; "
      f"they are Monte-Carlo estimates at {N_DRAWS:,} vs the deposited 1,000,000)")
R["gate"] = dict(passed=bool(ok), expr=g_expr, joint=g_joint,
                 deposited_expr_median=pub["expression_matched"]["null_median"],
                 deposited_joint_median=pub["joint_expr_tau_matched"]["null_median"])
if not ok:
    raise SystemExit("ABORT: faithfulness gate failed; nothing below would mean anything")

# ---------------------------------------------------------------------------
# Gene sets
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("GENE SETS")
print("=" * 78)
gid_to_pos = {g: i for i, g in enumerate(valid["gene_id"].values)}
sym_to_pos = {s: i for i, s in enumerate(valid["symbol"].values)}   # published: last wins

lam = pd.read_csv(LAMBERT, low_memory=False)
lam_tf = lam[lam["Is TF?"] == "Yes"]
lam_ens = set(lam_tf["Ensembl ID"].dropna().astype(str))
lam_sym = set(lam_tf["HGNC symbol"].dropna().astype(str))
print(f"  Lambert v_1.01: {len(lam)} rows, Is TF?==Yes: {len(lam_tf)}, "
      f"unique Ensembl IDs: {len(lam_ens)}")

# ARM A -- full TF census, mapped BY ENSEMBL ID (Dispatch 32 Block 4: feature_name
# is non-unique; the centroid axis is Ensembl-keyed, so a symbol join is the one
# avoidable error here)
census_pos = np.array(sorted({gid_to_pos[g] for g in lam_ens if g in gid_to_pos}), int)
print(f"  ARM A full TF census inside the {n}-gene valid space: {len(census_pos)}")
n_master_in_census = len({p for p in tf_pos} & set(census_pos.tolist()))
print(f"    ...of which are among the 73 master TFs: {n_master_in_census} / {len(tf_pos)}")
# sanity: how many would a SYMBOL join have found, for contrast
sym_join = len({sym_to_pos[s] for s in lam_sym if s in sym_to_pos})
print(f"    (a symbol-keyed join would have found {sym_join}; Ensembl join used)")

# ARM B -- marker-matched non-TF background
cm = pd.read_csv(G.CELLMARKER_H)
cm_syms = set(cm["gene_symbol"].dropna().astype(str))
cm_all_pos = {sym_to_pos[s] for s in cm_syms if s in sym_to_pos}
cm_nontf_pos = np.array(sorted(cm_all_pos - set(census_pos.tolist())), int)
print(f"  CellMarker symbols: {len(cm_syms)} unique; in valid space: {len(cm_all_pos)}")
print(f"  ARM B non-TF markers (CellMarker minus Lambert census): {len(cm_nontf_pos)}")
print(f"    removed as TFs: {len(cm_all_pos) - len(cm_nontf_pos)}")

R["sets"] = dict(master_tf=int(len(tf_pos)), tf_census=int(len(census_pos)),
                 master_in_census=int(n_master_in_census),
                 symbol_join_would_find=int(sym_join),
                 cellmarker_in_space=int(len(cm_all_pos)),
                 marker_nontf=int(len(cm_nontf_pos)),
                 lambert_rows=int(len(lam)), lambert_is_tf_yes=int(len(lam_tf)))

# ---------------------------------------------------------------------------
# Per-bin candidate counts (the degeneracy caution)
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("PER-BIN CANDIDATE COUNTS -- can the matched null actually draw?")
print("=" * 78)
def bin_report(label, pos, bins, binname):
    tb = bins[pos]
    need = pd.Series(tb).value_counts()
    avail = pd.Series(bins).value_counts()
    ratio = (need / avail.reindex(need.index)).sort_values(ascending=False)
    tight = int((need / avail.reindex(need.index) > 0.5).sum())
    exhaust = int((need >= avail.reindex(need.index)).sum())
    print(f"  {label} x {binname}: {len(need)} occupied bins, "
          f"need/avail max {ratio.iloc[0]:.3f}, median {ratio.median():.4f}")
    print(f"    bins where the query needs >50% of the bin: {tight}   "
          f"bins where it needs ALL of them (draw is deterministic): {exhaust}")
    return dict(n_bins=int(len(need)), max_ratio=float(ratio.iloc[0]),
                median_ratio=float(ratio.median()), n_gt_half=tight, n_exhausted=exhaust,
                min_avail=int(avail.reindex(need.index).min()))

R["bins"] = {}
for label, pos in (("73 master TFs", tf_pos), ("TF census", census_pos),
                   ("non-TF markers", cm_nontf_pos)):
    R["bins"][label] = dict(
        expr=bin_report(label, pos, ebins, "expr(20)"),
        joint=bin_report(label, pos, jbins, "joint(10x10)"))

# ---------------------------------------------------------------------------
# The three gene sets against both nulls
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print(f"RESULTS -- median C percentile vs both matched nulls (N={N_DRAWS:,} draws)")
print("=" * 78)
res = {}
for label, pos, seed in (("73 master TFs", tf_pos, 20240601),
                         ("full TF census (Lambert v_1.01)", census_pos, 20240701),
                         ("marker-matched non-TF background", cm_nontf_pos, 20240801)):
    o = float(np.median(C_rank[pos]))
    ns["obs"] = o          # high_n compares against the module-global `obs`
    e = high_n(ebins, seed, f"{label} / expr")
    j = high_n(jbins, seed + 1, f"{label} / joint")
    res[label] = dict(n_genes=int(len(pos)), obs_median=o, expr=e, joint=j)
    print(f"\n  {label}  (n={len(pos)})")
    print(f"    observed median C percentile      : {o:.4f}")
    print(f"    expression-matched null median    : {e['null_median']:.4f}   "
          f"p = {e['p_empirical']:.3e}   z = {e['z']:+.2f}")
    print(f"    joint expr+Tau-matched null median: {j['null_median']:.4f}   "
          f"p = {j['p_empirical']:.3e}   z = {j['z']:+.2f}")
ns["obs"] = obs
R["results"] = res

(W3 / "block2_w3_results.json").write_text(json.dumps(R, indent=2, default=float))
print(f"\nwrote {W3/'block2_w3_results.json'}")
