#!/usr/bin/env python3
"""
Matched-n primary baselines for every bar in the replication panel.

Restricts the tracked primary centroids to each replication's own matched type set,
runs the published pipeline unchanged, and reports the primary's obs/null at that n.
Then draws a 1,000-subset distribution at each distinct n, so each bar carries an
interval rather than being read against the primary's 35-type value.

Produces the matched-baseline comparison reported in Results, "each replication is
read against the primary restricted to that replication's own matched types", and
the accompanying S1 Text section on matched-n baselines.

All inputs are tracked.

Only the two path constants differ from the version that produced those values; the
code is unchanged.
"""
import contextlib
import io
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
OUT = HERE
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REPO / "src"))
from cellwarp.procrustes import pca_reduce_centroids, procrustes_align, permutation_test  # noqa

NPERM, SEED, NSUB = 10_000, 42, 1000
SUBSET_SEED = 42
PUBLISHED = 0.5222043226858066


@contextlib.contextmanager
def quiet():
    with contextlib.redirect_stdout(io.StringIO()):
        yield


hc = pd.read_csv(REPO / "output/phase2/scaled_35types/centroids_human_35.csv", index_col=0)
mc = pd.read_csv(REPO / "output/phase2/scaled_35types/centroids_mouse_35.csv", index_col=0)
saved = np.load(REPO / "output/phase2/scaled_35types/pca_centroids_35.npz", allow_pickle=True)
CT35 = [str(c) for c in saved["cell_types"]]
hc, mc = hc.loc[CT35], mc.loc[CT35]
genes = [c for c in hc.columns if c.startswith("ENSG")]
hc, mc = hc[genes], mc[genes]
print(f"primary centroids: {len(CT35)} types x {len(genes)} genes")


def run(types):
    """The published pipeline on a type subset. Returns obs, null median, ratio, p, k."""
    ts = sorted(types)
    with quiet():
        h, m, pca, cts = pca_reduce_centroids(hc.loc[ts], mc.loc[ts], 0.95)
        res = procrustes_align(h, m)
        p, null = permutation_test(h, m, n_permutations=NPERM, seed=SEED)
    k = int(h.shape[1])
    nm = float(np.median(null))
    return dict(n=len(ts), k=k, obs=float(res.distance), null_median=nm,
                obs_null=float(res.distance / nm), p=float(p),
                cumvar=float(np.cumsum(pca.explained_variance_ratio_)[k - 1]))


R = {"nperm": NPERM, "seed": SEED, "n_subsets": NSUB, "subset_seed": SUBSET_SEED,
     "published_target": PUBLISHED, "n_genes": len(genes)}

# ---------------------------------------------------------------- faithfulness gate
print("\n" + "=" * 84)
print("FAITHFULNESS GATE -- published pipeline on all 35 types")
print("=" * 84)
t0 = time.time()
full = run(CT35)
print(f"  n={full['n']}  k={full['k']}  cumvar={full['cumvar']:.6f}")
print(f"  obs={full['obs']:.9f}  null_median={full['null_median']:.9f}")
print(f"  obs/null = {full['obs_null']:.16f}")
print(f"  published= {PUBLISHED:.16f}")
print(f"  diff     = {full['obs_null'] - PUBLISHED:+.3e}")
print(f"  k == 33 : {full['k'] == 33}      ({time.time()-t0:.1f}s per call)")
gate = abs(full["obs_null"] - PUBLISHED) < 1e-9 and full["k"] == 33
print(f"  -> GATE {'PASS' if gate else 'FAIL'}")
R["gate"] = dict(passed=bool(gate), **full)
if not gate:
    raise SystemExit("ABORT: the published headline does not reproduce; nothing below is valid")

# ---------------------------------------------------------------- the seven bars
BARS = [
    ("Primary (TSxTMS)", "analysis/permutation_1M/results_1M.json", None, None),
    ("Sun2023 (10x v3)", "output/validation/sun2023_replication_expanded/sun2023_expanded.json",
     "procrustes", "cell_types"),
    ("PanSci (EasySci)", "output/validation/pansci_replication/pansci_replication.json",
     "procrustes", "cell_types"),
    ("CellHint (Human)", "output/validation/cellhint_replication/cellhint_replication.json",
     "procrustes", "cell_types"),
    ("pan-Census (pooled)", "analysis/census_replication/replication_results.json",
     "permutation_test", "ROOT:cell_types"),
    ("Andrews (liver)", "output/validation/andrews_replication/andrews_replication_results.json",
     None, "types"),
    ("MCAxHCA", "output/validation/t1a_replication/t1a_results.json",
     "t1a_procrustes", "PTR:per_type_residuals"),
]

print("\n" + "=" * 84)
print("PER-BAR MATCHED TYPE SETS")
print("=" * 84)
bars = []
for label, path, sub, tkey in BARS:
    d = json.loads((REPO / path).read_text())
    node = d[sub] if sub else d
    on = node.get("obs_null_ratio", d.get("obs_null_ratio"))
    npub = node.get("n_types", d.get("n_types", d.get("n_cell_types")))
    pv = node.get("p_value", d.get("p_value"))
    if tkey is None:
        tl = list(CT35)
    elif tkey.startswith("ROOT:"):
        tl = list(d[tkey[5:]])
    elif tkey.startswith("PTR:"):
        tl = sorted(d[tkey[4:]].keys())
    else:
        tl = list(node[tkey])
    tl = [str(t) for t in tl]
    assert len(tl) > 0, f"{label}: EMPTY type list"                     # note 11
    inter = sorted(set(tl) & set(CT35))
    assert len(inter) > 0, f"{label}: EMPTY intersection with the primary 35"
    match = (len(tl) == npub)
    print(f"\n  {label}")
    print(f"    deposited obs/null = {on:.10f}   n printed = {npub}   p = {pv:.4e}")
    print(f"    recovered type list length = {len(tl)}   matches printed n: {match}"
          + ("" if match else "   *** MISMATCH ***"))
    print(f"    intersection with the primary 35 = {len(inter)}"
          + ("" if len(inter) == len(tl) else f"   ({len(tl)-len(inter)} not in the primary)"))
    if len(inter) != len(tl):
        print(f"      not in primary: {sorted(set(tl)-set(CT35))}")
    bars.append(dict(label=label, dep_obs_null=float(on), n_printed=int(npub),
                     p=float(pv), types=tl, n_types=len(tl),
                     intersection=inter, n_inter=len(inter), n_matches_printed=bool(match)))

# ---------------------------------------------------------------- matched-n baselines
print("\n" + "=" * 84)
print("MATCHED-n PRIMARY BASELINES (primary centroids restricted to each bar's intersection)")
print("=" * 84)
print(f"  {'bar':22s} {'n':>3s} {'k':>3s} {'obs':>10s} {'null_med':>10s} {'obs/null':>10s} {'p':>11s}")
for b in bars:
    r = run(b["intersection"])
    b["matched_n_primary"] = r
    print(f"  {b['label']:22s} {r['n']:3d} {r['k']:3d} {r['obs']:10.4f} {r['null_median']:10.4f} "
          f"{r['obs_null']:10.6f} {r['p']:11.4e}")
R["bars"] = bars

# ---------------------------------------------------------------- subset distributions
NS = sorted({b["n_inter"] for b in bars} | {6, 12, 15, 16, 17, 22, 35})
print("\n" + "=" * 84)
print(f"RANDOM-SUBSET DISTRIBUTIONS -- {NSUB} subsets of the 35 types per n, seed {SUBSET_SEED}")
print("=" * 84)
print(f"  {'n':>3s} {'draws':>6s} {'median':>9s} {'p5':>9s} {'p95':>9s} {'min':>9s} {'max':>9s} {'med k':>6s} {'sec':>7s}")
dists = {}
for n in NS:
    t0 = time.time()
    rng = np.random.default_rng(SUBSET_SEED)
    ndraw = 1 if n == len(CT35) else NSUB          # n=35: only one possible subset
    vals, ks = [], []
    for _ in range(ndraw):
        sub = list(rng.choice(CT35, size=n, replace=False)) if n < len(CT35) else list(CT35)
        r = run(sub)
        vals.append(r["obs_null"]); ks.append(r["k"])
    v = np.array(vals)
    dists[n] = dict(n=n, n_draws=ndraw, median=float(np.median(v)),
                    p5=float(np.percentile(v, 5)), p95=float(np.percentile(v, 95)),
                    min=float(v.min()), max=float(v.max()),
                    median_k=float(np.median(ks)), values=[float(x) for x in v],
                    sec=round(time.time() - t0, 1))
    d = dists[n]
    print(f"  {n:3d} {ndraw:6d} {d['median']:9.5f} {d['p5']:9.5f} {d['p95']:9.5f} "
          f"{d['min']:9.5f} {d['max']:9.5f} {d['median_k']:6.1f} {d['sec']:7.1f}")
R["distributions"] = {str(k): {kk: vv for kk, vv in v.items() if kk != "values"}
                      for k, v in dists.items()}
R["distribution_values"] = {str(k): v["values"] for k, v in dists.items()}

# n = 35 NEGATIVE CONTROL
print("\n" + "=" * 84)
print("NEGATIVE CONTROL -- n = 35, where the only possible subset is the full set")
print("=" * 84)
d35 = dists[35]
print(f"  draws {d35['n_draws']}   value {d35['median']:.16f}")
print(f"  published                {PUBLISHED:.16f}")
print(f"  diff                     {d35['median']-PUBLISHED:+.3e}")
print(f"  range collapses to a point: {d35['min'] == d35['max']}")
negctl = abs(d35["median"] - PUBLISHED) < 1e-9 and d35["min"] == d35["max"]
print(f"  -> NEGCTL {'ok' if negctl else 'WRONG -- the subsetting path is not the published one; block void'}")
R["negative_control_passed"] = bool(negctl)

# ---------------------------------------------------------------- place each bar
print("\n" + "=" * 84)
print("EACH REPLICATION AGAINST ITS OWN MATCHED-n DISTRIBUTION")
print("=" * 84)
print(f"  {'bar':22s} {'n':>3s} {'deposited':>10s} {'matched-n':>10s} {'subset med':>11s} "
      f"{'[p5, p95]':>20s} {'pctile':>7s} {'z':>7s}")
for b in bars:
    n = b["n_inter"]
    v = np.array(dists[n]["values"])
    dep = b["dep_obs_null"]
    pct = float((v <= dep).mean() * 100)
    z = float((dep - v.mean()) / v.std(ddof=1)) if len(v) > 1 else float("nan")
    b["placement"] = dict(percentile=pct, z=z, subset_median=dists[n]["median"],
                          subset_p5=dists[n]["p5"], subset_p95=dists[n]["p95"])
    print(f"  {b['label']:22s} {n:3d} {dep:10.6f} {b['matched_n_primary']['obs_null']:10.6f} "
          f"{dists[n]['median']:11.6f} [{dists[n]['p5']:8.5f},{dists[n]['p95']:8.5f}] "
          f"{pct:7.1f} {z:+7.2f}")

(OUT / "block2_matched_n_results.json").write_text(json.dumps(R, indent=2))
print(f"\nwrote {OUT/'block2_matched_n_results.json'}")
