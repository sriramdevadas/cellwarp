"""
Build main Figure 2 — the conserved-contribution gene set.

Publication-quality 4-panel figure (Cell Systems style via cellwarp.figure_style):
  7A  distribution of the per-gene conservation score C across the valid orthologs
      (a broad continuum, no discrete module)
  7B  C vs mean expression (rho ~ 0.22) and C vs Tau specificity (rho ~ 0.06):
      the two independence controls
  7C  the master-TF C-percentiles against the expression-matched and the joint
      (expression + Tau)-matched backgrounds
  7D  donor-split reproducibility of per-gene C, with the cell-sampling ceiling
      and the shuffle null marked

Sources (all verified, produced in-place against the deposit's tracked inputs):
  gene_conservation_core.csv, gate_results.json, robustness_results.json,
  donor_stability/donor_stability_results.json, and the gitignored Census
  aggregates in donor_stability/ for the donor-split recomputation.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec

HERE = Path(__file__).resolve().parent
DONOR = HERE / "donor_stability"
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))                 # gate_lib
sys.path.insert(0, str(ROOT / "src"))         # cellwarp.figure_style
import gate_lib as G
from cellwarp import figure_style as fs

fs.apply_style()

# ---------------------------------------------------------------------------
# Load verified data
# ---------------------------------------------------------------------------
df, h_cent, m_cent = G.build_gene_table()
valid = df.dropna(subset=["C_pearson"]).reset_index(drop=True)
C = valid["C_pearson"].values
me = valid["mean_expression"].values
n = len(valid)

Hc = h_cent[valid["gene_id"]].values  # 35 x n
Tau = np.array([np.sum(1 - Hc[:, j] / Hc[:, j].max()) / (Hc.shape[0] - 1)
                if Hc[:, j].max() > 0 else np.nan for j in range(n)])

gate = json.load(open(HERE / "gate_results.json"))
rob = json.load(open(HERE / "robustness_results.json"))
don = json.load(open(DONOR / "donor_stability_results.json"))

q75 = gate["thresholds"]["q75"]
q25 = gate["thresholds"]["q25"]
rho_expr = gate["check2"]["spearman_C_vs_expr"]
rho_tau = rob["rho_C_tau"]
dip_D = rob["dip_D"]

# C percentile rank and master-TF positions
C_rank = stats.rankdata(C) / n
sym2pos = {s: i for i, s in enumerate(valid["symbol"].values)}
tf_pos = np.array(sorted({sym2pos[s] for s in G.POSITIVE_CONTROL_TFS if s in sym2pos}))
tf_obs_median = float(np.median(C_rank[tf_pos]))

# ---------------------------------------------------------------------------
# 7C: matched-null median-percentile distributions (regenerated, seeded)
# ---------------------------------------------------------------------------
bins = G.expr_bins(me, n_bins=20)
def jbin(v, k):
    r = stats.rankdata(v, method="ordinal")
    return np.minimum((r - 1) * k // len(r), k - 1)
jbins = jbin(me, 10) * 10 + jbin(np.nan_to_num(Tau, nan=np.nanmedian(Tau)), 10)

rng = np.random.default_rng(42)
null_e = np.array([np.median(C_rank[d]) for d in G.matched_draws(tf_pos, bins, 1000, rng)])
rng = np.random.default_rng(123)
null_j = np.array([np.median(C_rank[d]) for d in G.matched_draws(tf_pos, jbins, 1000, rng)])

# ---------------------------------------------------------------------------
# 7D: donor-split cross-half C, cell-sampling ceiling (light recompute from npz)
# ---------------------------------------------------------------------------
genes = list(np.load(DONOR / "agg_human_cap10000.npz", allow_pickle=True)["genes"])
core = pd.read_csv(HERE / "gene_conservation_core.csv").set_index("gene_id").loc[genes].reset_index()
valid_d = core["C_pearson"].notna().values
NG = len(genes)

def Cvec(H, M):
    ok = ~(np.isnan(H).any(1) | np.isnan(M).any(1))
    H, M = H[ok], M[ok]
    Hcc = H - H.mean(0); Mcc = M - M.mean(0)
    den = np.sqrt((Hcc ** 2).sum(0) * (Mcc ** 2).sum(0))
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(den == 0, np.nan, (Hcc * Mcc).sum(0) / den)
    return r

class Agg:
    def __init__(s, tag, cap=10000):
        d = np.load(DONOR / f"agg_{tag}_cap{cap}.npz", allow_pickle=True); s.g = d["gsums"]
        gg = pd.read_csv(DONOR / f"agg_{tag}_cap{cap}_groups.csv")
        s.t = gg.type_idx.values; s.dn = gg.donor.astype(str).values
        s.c = gg["count"].values.astype(float)
    def donors(s): return np.unique(s.dn)
    def cen(s, ds):
        sel = np.isin(s.dn, np.array(list(ds))); H = np.full((35, NG), np.nan)
        for t in range(35):
            m = sel & (s.t == t); cc = s.c[m].sum()
            if cc > 0:
                H[t] = s.g[m].sum(0) / cc
        return H

H, M = Agg("human"), Agg("mouse")
hd, md = H.donors(), M.donors()
cross = []
for s in range(100):
    rs = np.random.default_rng(100 + s)
    hp = rs.permutation(hd); mp = rs.permutation(md)
    CA = Cvec(H.cen(set(hp[:len(hp)//2])), M.cen(set(mp[:len(mp)//2])))
    CB = Cvec(H.cen(set(hp[len(hp)//2:])), M.cen(set(mp[len(mp)//2:])))
    m = valid_d & ~np.isnan(CA) & ~np.isnan(CB)
    cross.append(stats.spearmanr(CA[m], CB[m])[0])
cross = np.array(cross)

csh = np.load(DONOR / "agg_human_cs.npz"); csm = np.load(DONOR / "agg_mouse_cs.npz")
def csc(cs, cnt, r):
    Hh = cs[r].astype(float).copy(); c = cnt[r].astype(float)
    Hh[c > 0] /= c[c > 0, None]; Hh[c == 0] = np.nan; return Hh
ceil = []
for r in range(20):
    C1 = Cvec(csc(csh["csA"], csh["cntA"], r), csc(csm["csA"], csm["cntA"], r))
    C2 = Cvec(csc(csh["csB"], csh["cntB"], r), csc(csm["csB"], csm["cntB"], r))
    m = valid_d & ~np.isnan(C1) & ~np.isnan(C2)
    ceil.append(stats.spearmanr(C1[m], C2[m])[0])
ceil_med = float(np.median(ceil))
cross_med = float(np.median(cross))
null95 = don["null_shuffle_spearman_95"]
frac_both = don["donor_split_cap10000"]["frac_halves_BOTH"]

print(f"recomputed: cross-half median={cross_med:.3f} (json {don['donor_split_cap10000']['cross_half_C_spearman_median']:.3f}), "
      f"ceiling={ceil_med:.3f} (json {don['ceiling_cellsplit_spearman_median']:.3f})")
print(f"7C: TF obs median={tf_obs_median:.3f}; null_e median={np.median(null_e):.3f}; null_j median={np.median(null_j):.3f}")

# ===========================================================================
# Figure layout: 2 rows x 6 cols. Top: A | B(expr) | B(Tau). Bottom: C | D.
# ===========================================================================
fig = plt.figure(figsize=(fs.COL2, 5.4))
gs = gridspec.GridSpec(2, 6, figure=fig, hspace=0.62, wspace=1.7,
                       height_ratios=[1, 1])
axA = fig.add_subplot(gs[0, 0:2])
axBe = fig.add_subplot(gs[0, 2:4])
axBt = fig.add_subplot(gs[0, 4:6])
axC = fig.add_subplot(gs[1, 0:3])
axD = fig.add_subplot(gs[1, 3:6])
for ax in (axA, axBe, axBt, axC, axD):
    fs.clean_spine(ax)

# ---- 7A: distribution of C --------------------------------------------------
counts, edges = np.histogram(C, bins=60)
centers = 0.5 * (edges[:-1] + edges[1:])
wdt = edges[1] - edges[0]
col = np.where(centers >= q75, fs.C_TEAL, np.where(centers <= q25, fs.C_LIGHTGRAY, fs.C_BLUE))
axA.bar(centers, counts, width=wdt, color=col, linewidth=0)
axA.axvline(q75, color=fs.C_DARKGRAY, ls="--", lw=0.8)
axA.axvline(q25, color=fs.C_DARKGRAY, ls="--", lw=0.8)
ymax = counts.max()
axA.text(q75 + 0.02, ymax * 0.96, f"conserved\n(Q75={q75:.2f})", fontsize=6,
         color=fs.C_TEAL, va="top", ha="left")
axA.text(q25 - 0.02, ymax * 0.96, f"divergent\n(Q25={q25:.2f})", fontsize=6,
         color=fs.C_GRAY, va="top", ha="right")
axA.set_xlabel("Conservation score $C$")
axA.set_ylabel("Orthologs")
axA.set_xlim(-0.6, 1.02)
axA.set_title(f"$n$ = {n:,} valid orthologs\nbroad continuum (dip $D$ = {dip_D:.3f})",
              fontsize=7, loc="center", pad=4)
fs.add_panel_label(axA, "A", x=-0.20, y=1.30)

# ---- 7B: independence controls (C vs expression; C vs Tau) ------------------
def decile_means(x, y, k=10):
    dec = pd.qcut(pd.Series(x).rank(method="first"), k, labels=False)
    xs = [float(np.mean(x[dec == d])) for d in range(k)]
    ys = [float(np.mean(y[dec == d])) for d in range(k)]
    return np.array(xs), np.array(ys)

axBe.hexbin(me, C, gridsize=40, cmap="Blues", mincnt=1, bins="log", linewidths=0)
dx, dy = decile_means(me, C)
axBe.plot(dx, dy, "o-", color=fs.C_ORANGE, ms=3, lw=1.2, label="decile mean")
axBe.set_xlabel("Mean expression (log1p)")
axBe.set_ylabel("Conservation score $C$")
axBe.set_title(f"vs expression\nSpearman $\\rho$ = {rho_expr:.2f}", fontsize=7, pad=4)
axBe.legend(loc="lower right", fontsize=6, frameon=False)
fs.add_panel_label(axBe, "B", x=-0.34, y=1.30)

axBt.hexbin(Tau, C, gridsize=40, cmap="Purples", mincnt=1, bins="log", linewidths=0)
tx, ty = decile_means(Tau, C)
axBt.plot(tx, ty, "o-", color=fs.C_PURPLE, ms=3, lw=1.2, label="decile mean")
axBt.set_xlabel("Specificity (Tau)")
axBt.set_title(f"vs specificity\nSpearman $\\rho$ = {rho_tau:.2f}", fontsize=7, pad=4)
axBt.set_ylim(axBe.get_ylim())
axBt.legend(loc="lower right", fontsize=6, frameon=False)

# ---- 7C: master-TF percentiles vs matched backgrounds ----------------------
axC.hist(null_e, bins=40, density=True, color=fs.C_GRAY, alpha=0.55,
         label=f"expr-matched null (med {np.median(null_e):.2f})")
axC.hist(null_j, bins=40, density=True, color=fs.C_TEAL, alpha=0.55,
         label=f"expr+Tau-matched null (med {np.median(null_j):.2f})")
axC.axvline(tf_obs_median, color=fs.C_ORANGE, lw=2.0,
            label=f"master TFs (median {tf_obs_median:.2f})")
top = axC.get_ylim()[1]
axC.set_ylim(-0.06 * top, 1.34 * top)            # bottom room for rug; top headroom for legend
# rug of the individual 73 TF percentiles
axC.plot(C_rank[tf_pos], np.full(len(tf_pos), -0.035 * top), "|", color=fs.C_ORANGE,
         ms=6, mew=0.9, clip_on=False)
# value annotation at the observed master-TF median line
axC.text(tf_obs_median - 0.015, 1.30 * top, f"{tf_obs_median:.2f}", color=fs.C_ORANGE,
         fontsize=7, fontweight="bold", ha="right", va="top")
axC.set_xlim(0, 1)
axC.set_xlabel("$C$ percentile rank")
axC.set_ylabel("Density (matched-null medians)")
axC.set_title(f"{len(tf_pos)} master TFs vs matched backgrounds · $p < 10^{{-6}}$ vs both",
              fontsize=7, pad=4)
axC.legend(loc="upper left", fontsize=6, frameon=False)
fs.add_panel_label(axC, "C", x=-0.13, y=1.18)

# ---- 7D: donor-split reproducibility of per-gene C -------------------------
axD.hist(cross, bins=18, color=fs.C_BLUE, alpha=0.85, label="donor-split (n=100)")
axD.axvline(cross_med, color=fs.C_DARKGRAY, lw=1.8,
            label=f"donor-split median = {cross_med:.2f}")
axD.axvline(ceil_med, color=fs.C_TEAL, ls="--", lw=1.4,
            label=f"cell-sampling ceiling = {ceil_med:.2f}")
axD.axvline(null95, color=fs.C_GRAY, ls=":", lw=1.4,
            label=f"shuffle null (95th = {null95:.3f})")
top_d = axD.get_ylim()[1]
axD.set_ylim(0, 1.32 * top_d)                     # headroom for legend + median annotation
axD.text(cross_med - 0.012, 1.28 * top_d, f"{cross_med:.2f}", color=fs.C_DARKGRAY,
         fontsize=7, fontweight="bold", ha="right", va="top")
axD.set_xlim(0, 1)
axD.set_xlabel("Cross-half agreement of $C$  (Spearman $C_A$, $C_B$)")
axD.set_ylabel("Donor splits")
axD.set_title(f"donor-split median = {cross_med:.2f} · master TFs recovered in "
              f"{int(frac_both*200)}/200 halves", fontsize=7, pad=4)
axD.legend(loc="upper left", fontsize=6, frameon=False)
fs.add_panel_label(axD, "D", x=-0.13, y=1.18)

# figure-wide definition of the core quantity
fig.text(0.5, -0.01,
         "$C$ = per-gene Pearson correlation of the expression profile across the 35 matched "
         "cell-type centroids, human vs mouse (computed independently of the Procrustes axes).",
         ha="center", va="top", fontsize=6, style="italic", color=fs.C_DARKGRAY)

# ---------------------------------------------------------------------------
stem = ROOT / "figures" / "main" / "fig7_conserved_contribution"
fig.savefig(str(stem) + ".pdf", format="pdf", dpi=fs.DPI, bbox_inches="tight", pad_inches=0.05)
fig.savefig(str(stem) + ".png", format="png", dpi=fs.DPI, bbox_inches="tight", pad_inches=0.05)
print("saved", stem.with_suffix(".pdf").name, "+ .png (300 dpi)")
