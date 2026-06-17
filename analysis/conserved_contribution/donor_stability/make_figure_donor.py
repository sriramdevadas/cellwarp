"""Donor-stability figure: per-split distributions (light recompute, no matched-draws)."""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import gate_lib as G

genes = list(np.load(HERE / "agg_human_cap10000.npz", allow_pickle=True)["genes"])
core = pd.read_csv(HERE.parent / "gene_conservation_core.csv").set_index("gene_id").loc[genes].reset_index()
valid = core["C_pearson"].notna().values
sym = core["symbol"].values
sym2pos = {}
[sym2pos.setdefault(s, i) for i, s in enumerate(sym)]
tf_pc = np.array(sorted({sym2pos[s] for s in G.POSITIVE_CONTROL_TFS if s in sym2pos}))
NG = len(genes)
res = json.load(open(HERE / "donor_stability_results.json"))

def Cvec(H, M):
    ok = ~(np.isnan(H).any(1) | np.isnan(M).any(1)); H, M = H[ok], M[ok]
    Hc = H - H.mean(0); Mc = M - M.mean(0)
    den = np.sqrt((Hc**2).sum(0) * (Mc**2).sum(0))
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(den == 0, np.nan, (Hc*Mc).sum(0)/den)
    return r

class Agg:
    def __init__(s, tag, cap=10000):
        d = np.load(HERE/f"agg_{tag}_cap{cap}.npz", allow_pickle=True); s.g = d["gsums"]
        gg = pd.read_csv(HERE/f"agg_{tag}_cap{cap}_groups.csv")
        s.t = gg.type_idx.values; s.dn = gg.donor.astype(str).values; s.c = gg["count"].values.astype(float)
    def donors(s): return np.unique(s.dn)
    def cen(s, ds):
        sel = np.isin(s.dn, np.array(list(ds))); H = np.full((35, NG), np.nan)
        for t in range(35):
            m = sel & (s.t == t); cc = s.c[m].sum()
            if cc > 0: H[t] = s.g[m].sum(0)/cc
        return H
H, M = Agg("human"), Agg("mouse")
hd, md = H.donors(), M.donors()
def tfpct(C):
    v = valid & ~np.isnan(C); idx = np.where(v)[0]; rank = stats.rankdata(C[idx])/len(idx)
    loc = {g: i for i, g in enumerate(idx)}
    return np.median([rank[loc[g]] for g in tf_pc if g in loc])

cross, tfp = [], []
for s in range(100):
    rs = np.random.default_rng(100+s); hp = rs.permutation(hd); mp = rs.permutation(md)
    CA = Cvec(H.cen(set(hp[:len(hp)//2])), M.cen(set(mp[:len(mp)//2])))
    CB = Cvec(H.cen(set(hp[len(hp)//2:])), M.cen(set(mp[len(mp)//2:])))
    m = valid & ~np.isnan(CA) & ~np.isnan(CB)
    cross.append(stats.spearmanr(CA[m], CB[m])[0]); tfp += [tfpct(CA), tfpct(CB)]
# ceiling
csh, csm = np.load(HERE/"agg_human_cs.npz"), np.load(HERE/"agg_mouse_cs.npz")
def csc(cs, cnt, r):
    Hh = cs[r].astype(float).copy(); c = cnt[r].astype(float)
    Hh[c > 0] /= c[c > 0, None]; Hh[c == 0] = np.nan; return Hh
ceil = []
for r in range(20):
    C1 = Cvec(csc(csh["csA"], csh["cntA"], r), csc(csm["csA"], csm["cntA"], r))
    C2 = Cvec(csc(csh["csB"], csh["cntB"], r), csc(csm["csB"], csm["cntB"], r))
    m = valid & ~np.isnan(C1) & ~np.isnan(C2); ceil.append(stats.spearmanr(C1[m], C2[m])[0])

fig, ax = plt.subplots(2, 2, figsize=(12, 9))
a = ax[0, 0]
a.hist(cross, bins=20, color="#4C772B0".replace("#","#") if False else "#4C72B0", alpha=.85, label="donor-split (100)")
a.axvline(np.median(ceil), color="green", ls="--", lw=2, label=f"cell ceiling={np.median(ceil):.2f}")
a.axvline(res["null_shuffle_spearman_95"], color="red", ls="--", lw=2, label=f"null95={res['null_shuffle_spearman_95']:.2f}")
a.axvline(np.median(cross), color="navy", ls="-", lw=1)
a.set_title(f"A. Cross-half C agreement\ndonor-split median={np.median(cross):.2f}, gap to ceiling={np.median(ceil)-np.median(cross):.2f}")
a.set_xlabel("Spearman(C_halfA, C_halfB)"); a.set_xlim(0, 1); a.legend(fontsize=8)

b = ax[0, 1]
b.hist(tfp, bins=20, color="crimson", alpha=.8)
b.axvline(0.75, color="k", ls="--", label="recovery threshold 0.75")
b.axvline(np.median(tfp), color="darkred", lw=2, label=f"median={np.median(tfp):.2f}")
b.set_title(f"B. Master-TF recovery per donor-half (n=200)\n100% halves > 0.75; enrichment sig in 100%")
b.set_xlabel("median TF conservation percentile"); b.set_xlim(0, 1); b.legend(fontsize=8)

c = ax[1, 0]
caps = [500, 2000, 10000]; ch = [res["caps"][str(k)]["cross_half"] for k in caps]
both = [res["caps"][str(k)]["frac_halves_BOTH"] for k in caps]
c.plot(range(3), ch, "o-", color="#4C72B0", label="cross-half C")
c.plot(range(3), both, "s-", color="crimson", label="frac halves passing BOTH")
c.set_xticks(range(3)); c.set_xticklabels([f"cap {k}" for k in caps])
c.set_ylim(0, 1.05); c.set_title("C. Cell-count-cap robustness\n(TF finding cap-invariant -> not power-limited)")
c.set_ylabel("value"); c.legend(fontsize=8)

d = ax[1, 1]
cp = res["cross_protocol"]
labels = ["per-gene C\n(10x vs SS)", "TF pct\n10x", "TF pct\nSmart-seq2", "donor-split\ncross-half C"]
vals = [cp["spearman_C10x_vs_CSS"], cp["tf_pct_10x"], cp["tf_pct_SS"], np.median(cross)]
cols = ["#888", "#2a9d8f", "#264653", "#4C72B0"]
d.bar(labels, vals, color=cols)
for i, v in enumerate(vals): d.text(i, v+.02, f"{v:.2f}", ha="center", fontsize=9)
d.axhline(res["null_shuffle_spearman_95"], color="red", ls=":", label="null95")
d.set_ylim(0, 1.0); d.set_title("D. Cross-protocol (cross-'site' analogue)\nTF finding holds in each protocol (p<=0.004)")
d.set_ylabel("Spearman / percentile"); d.tick_params(axis="x", labelsize=8); d.legend(fontsize=8)

plt.tight_layout(); plt.savefig(HERE/"donor_stability_figure.png", dpi=140, bbox_inches="tight")
print("saved donor_stability_figure.png | cross median", round(float(np.median(cross)),3),
      "| ceiling", round(float(np.median(ceil)),3), "| tfp median", round(float(np.median(tfp)),3))
