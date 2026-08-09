#!/usr/bin/env python3
"""
Whether replicated cell types are the more conserved ones in the primary.

The matched-n primary baselines sit below their own random-subset medians in six of
seven bars. That is not an n effect, since each baseline is at its own n, so the
types the replications happen to match are as a set more conserved in the primary
than a random draw of the same size. This tests the obvious mechanism: that the
matched types are the low-residual ones.

Backs the claim in Results that "each replication's matched types are also, with one
exception, more conserved in the primary than the types it misses, which is why the
matched baselines fall below random type sets of the same size". The statistics this
script computes are not quoted in the submitted text; the sentence states the
relationship without them, and none of the tests here reaches significance.

NOT TRACKED: reads block2_matched_n_results.json, written by block2_matched_n.py in
this directory. That file is not deposited; run block2_matched_n.py first. The
residual ranking it reads is tracked.

Only the three path constants differ from the version that produced those values;
the code is unchanged.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
W5 = HERE
OUT = HERE
OUT.mkdir(parents=True, exist_ok=True)

# ---- per-type Procrustes residuals for the 35 primary types ------------------
res = pd.read_csv(REPO / "output/phase2/scaled_35types/residuals_ranked.csv")
res = res.set_index("cell_type")["residual_magnitude"]
print(f"per-type residuals: {len(res)} types   "
      f"range [{res.min():.4f}, {res.max():.4f}]   median {res.median():.4f}")
assert len(res) == 35

# ---- each replication's matched set, from Block 2's recovered lists ----------
B2 = json.loads((W5 / "block2_matched_n_results.json").read_text())
bars = {b["label"]: b["types"] for b in B2["bars"] if b["label"] != "Primary (TSxTMS)"}
for k, v in bars.items():
    print(f"  {k:22s} n={len(v)}")
assert len(bars) == 6

CT35 = list(res.index)
count = pd.Series({t: sum(t in v for v in bars.values()) for t in CT35})
ever = count > 0
print(f"\ntypes appearing in >=1 replication matched set: {int(ever.sum())} / 35")
print(f"types appearing in NONE                        : {int((~ever).sum())} / 35")
assert ever.sum() > 0 and (~ever).sum() > 0        # note 11: both arms non-empty

R = {"n_types": 35, "n_ever": int(ever.sum()), "n_never": int((~ever).sum())}

# ---- (A) the split the user asked for ---------------------------------------
a, b = res[ever].values, res[~ever].values
u, pu = stats.mannwhitneyu(a, b, alternative="two-sided")
t, pt = stats.ttest_ind(a, b, equal_var=False)
print("\n" + "=" * 78)
print("(A) mean per-type residual, split by 'appears in ANY replication matched set'")
print("=" * 78)
print(f"  appears in >=1 (n={len(a)}): mean {a.mean():.4f}  sd {a.std(ddof=1):.4f}  median {np.median(a):.4f}")
print(f"  appears in none (n={len(b)}): mean {b.mean():.4f}  sd {b.std(ddof=1):.4f}  median {np.median(b):.4f}")
print(f"  difference (replicated - not) = {a.mean()-b.mean():+.4f}")
print(f"  Mann-Whitney U = {u:.1f}   p = {pu:.4f}")
print(f"  Welch t = {t:+.3f}   p = {pt:.4f}")
R["split"] = dict(n_ever=len(a), n_never=len(b), mean_ever=float(a.mean()),
                  mean_never=float(b.mean()), median_ever=float(np.median(a)),
                  median_never=float(np.median(b)), diff=float(a.mean()-b.mean()),
                  mannwhitney_u=float(u), mannwhitney_p=float(pu),
                  welch_t=float(t), welch_p=float(pt))

# ---- (B) dose-response: residual vs how MANY replications match the type ----
rho, prho = stats.spearmanr(count.values, res.values)
print("\n" + "=" * 78)
print("(B) dose-response -- residual vs number of replications matching the type")
print("=" * 78)
print(f"  Spearman rho = {rho:+.4f}   p = {prho:.4f}   n = 35")
print(f"\n  {'#reps':>5s} {'n types':>8s} {'mean resid':>11s} {'median':>9s}")
for c in sorted(count.unique()):
    m = count == c
    print(f"  {c:5d} {int(m.sum()):8d} {res[m].mean():11.4f} {res[m].median():9.4f}")
R["dose_response"] = dict(spearman_rho=float(rho), spearman_p=float(prho),
                          by_count={int(c): dict(n=int((count == c).sum()),
                                                 mean=float(res[count == c].mean()),
                                                 median=float(res[count == c].median()))
                                    for c in sorted(count.unique())})

# ---- (C) per replication, its matched types vs the rest ---------------------
print("\n" + "=" * 78)
print("(C) per replication -- mean residual of its matched types vs the other 35-n")
print("=" * 78)
print(f"  {'replication':22s} {'n':>3s} {'matched mean':>13s} {'rest mean':>10s} {'diff':>8s} {'MW p':>8s}")
per = {}
for lbl, tl in bars.items():
    inn = res[[t for t in CT35 if t in set(tl)]].values
    out = res[[t for t in CT35 if t not in set(tl)]].values
    _, p = stats.mannwhitneyu(inn, out, alternative="two-sided")
    per[lbl] = dict(n=len(inn), matched_mean=float(inn.mean()), rest_mean=float(out.mean()),
                    diff=float(inn.mean()-out.mean()), mw_p=float(p))
    print(f"  {lbl:22s} {len(inn):3d} {inn.mean():13.4f} {out.mean():10.4f} "
          f"{inn.mean()-out.mean():+8.4f} {p:8.4f}")
R["per_replication"] = per

# ---- (D) does it actually explain item 1? -----------------------------------
# item 1's gap = (matched-n baseline) - (random-subset median at the same n)
print("\n" + "=" * 78)
print("(D) does residual rank explain item 1's gap?")
print("=" * 78)
print(f"  {'replication':22s} {'n':>3s} {'matched-n':>10s} {'subset med':>11s} {'item-1 gap':>11s} {'resid diff':>11s}")
rows = []
for bb in B2["bars"]:
    if bb["label"] == "Primary (TSxTMS)":
        continue
    gap = bb["matched_n_primary"]["obs_null"] - bb["placement"]["subset_median"]
    rd = per[bb["label"]]["diff"]
    rows.append((bb["label"], bb["n_inter"], bb["matched_n_primary"]["obs_null"],
                 bb["placement"]["subset_median"], gap, rd))
    print(f"  {bb['label']:22s} {bb['n_inter']:3d} {bb['matched_n_primary']['obs_null']:10.6f} "
          f"{bb['placement']['subset_median']:11.6f} {gap:+11.6f} {rd:+11.4f}")
g = np.array([r[4] for r in rows]); d = np.array([r[5] for r in rows])
rho2, p2 = stats.spearmanr(d, g)
print(f"\n  Spearman(residual difference, item-1 gap) = {rho2:+.4f}   p = {p2:.4f}   n = {len(rows)}")
print(f"  (a POSITIVE rho means: the more low-residual the matched set, the more")
print(f"   negative the gap -- i.e. residual rank explains the item-1 effect)")
R["explains_item1"] = dict(spearman_rho=float(rho2), spearman_p=float(p2), n=len(rows),
                           rows=[dict(label=r[0], n=r[1], matched_n=r[2],
                                      subset_median=r[3], gap=r[4], resid_diff=r[5])
                                 for r in rows])

(OUT / "task2_residual_mechanism.json").write_text(json.dumps(R, indent=2))
print(f"\nwrote {OUT/'task2_residual_mechanism.json'}")
