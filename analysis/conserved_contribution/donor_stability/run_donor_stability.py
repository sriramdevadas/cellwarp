"""
Donor-stability gate. Consumes the per-(type,donor,proto) aggregates + cell-split
half-centroids from pull_aggregate.py. C and all controls are FROZEN (gate_lib).
Only donor resampling is new.

Outputs: donor_stability_results.json (+ printed report).
"""
from __future__ import annotations
import io, json, contextlib, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))  # for gate_lib
import gate_lib as G
from cellwarp.procrustes import pca_reduce_centroids, procrustes_align, permutation_test

ROOT = G.ROOT
NSPLIT = 100
NSPLIT_CAP = 30
K_CEIL = 20
NDRAW = 200
rng = np.random.default_rng(7)

# ---- frozen gene annotations / controls (from the conserved-contribution gate) ----
core = pd.read_csv(HERE.parent / "gene_conservation_core.csv")
genes = list(np.load(HERE / "agg_human_cap10000.npz", allow_pickle=True)["genes"])
core = core.set_index("gene_id").loc[genes].reset_index()
NG = len(genes)
mean_expr = core["mean_expression"].values.astype(float)
valid = core["C_pearson"].notna().values & np.isfinite(mean_expr)
# Tau specificity from deposit (frozen)
h_dep, m_dep = G.load_centroids()
Hdep = h_dep[genes].values
Tau = np.array([np.sum(1 - Hdep[:, j] / Hdep[:, j].max()) / (Hdep.shape[0] - 1)
                if Hdep[:, j].max() > 0 else np.nan for j in range(NG)])
# frozen joint expr x Tau bins (the hardened control)
def jbin(v, k):
    r = stats.rankdata(v, method="ordinal")
    return np.minimum((r - 1) * k // len(r), k - 1)
jbins = jbin(mean_expr, 10) * 10 + jbin(np.nan_to_num(Tau, nan=np.nanmedian(Tau)), 10)

sym = core["symbol"].values
sym2pos = {}
for i, s in enumerate(sym):
    sym2pos.setdefault(s, i)
def syms_pos(symbols):
    return np.array(sorted({sym2pos[s] for s in symbols if s in sym2pos}), int)
tf_pc = syms_pos(G.POSITIVE_CONTROL_TFS)
tf_all = syms_pos(pd.read_csv(G.TF_ACTIVITY, index_col=0).columns.tolist())
cm = syms_pos(pd.read_csv(G.CELLMARKER_H)["gene_symbol"].dropna().unique().tolist())

# ---- vectorized per-gene Pearson across types (verified == gate_lib) ----
def Cvec(H, M):
    ok = ~(np.isnan(H).any(1) | np.isnan(M).any(1))
    H, M = H[ok], M[ok]
    Hc = H - H.mean(0); Mc = M - M.mean(0)
    num = (Hc * Mc).sum(0)
    den = np.sqrt((Hc ** 2).sum(0) * (Mc ** 2).sum(0))
    with np.errstate(invalid="ignore", divide="ignore"):
        r = num / den
    r[den == 0] = np.nan
    return r, int(ok.sum())

# ---- aggregate loader / centroid reconstruction ----
class Agg:
    def __init__(self, tag, cap):
        d = np.load(HERE / f"agg_{tag}_cap{cap}.npz", allow_pickle=True)
        self.gsums = d["gsums"]
        g = pd.read_csv(HERE / f"agg_{tag}_cap{cap}_groups.csv")
        self.tidx = g["type_idx"].values
        self.donor = g["donor"].astype(str).values
        self.proto = g["proto"].astype(str).values
        self.count = g["count"].values.astype(float)
        self.types = list(d["types"])
    def all_donors(self):
        return np.unique(self.donor)
    def centroid(self, donor_set=None, proto=None):
        sel = np.ones(len(self.tidx), bool)
        if donor_set is not None:
            sel &= np.isin(self.donor, np.array(list(donor_set)))
        if proto is not None:
            sel &= (self.proto == proto)
        H = np.full((35, NG), np.nan)
        for t in range(35):
            m = sel & (self.tidx == t)
            c = self.count[m].sum()
            if c > 0:
                H[t] = self.gsums[m].sum(0) / c
        return H

def obs_null(H, M, nperm=2000):
    ok = ~(np.isnan(H).any(1) | np.isnan(M).any(1))
    hc = pd.DataFrame(H[ok]); mc = pd.DataFrame(M[ok])
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        hp, mp, _, _ = pca_reduce_centroids(hc, mc, 0.95)
        res = procrustes_align(hp, mp)
        _, null = permutation_test(hp, mp, nperm)
    return float(res.distance / np.median(null))

# ---- the frozen master-TF finding evaluated on a C vector ----
def eval_finding(C, ndraw=NDRAW):
    v = valid & ~np.isnan(C)
    idx = np.where(v)[0]
    Cv = C[idx]
    rank = stats.rankdata(Cv) / len(Cv)
    loc = {g: i for i, g in enumerate(idx)}
    tf_loc = np.array([loc[g] for g in tf_pc if g in loc])
    med_pct = float(np.median(rank[tf_loc])) if len(tf_loc) else np.nan
    frac_above = float(np.mean(Cv[tf_loc] > np.median(Cv))) if len(tf_loc) else np.nan
    # conserved set = top quartile; TF(CollecTRI) enrichment vs joint-matched bg
    q = np.quantile(Cv, 0.75)
    cons = np.where(Cv >= q)[0]
    jb = jbins[idx]
    def overlap(ref_global):
        ref = set(loc[g] for g in ref_global if g in loc)
        obs = len(set(cons) & ref)
        draws = G.matched_draws(cons, jb, ndraw, rng)
        null = np.array([len(set(d.tolist()) & ref) for d in draws])
        return obs, float(np.mean(null)), G.emp_p_greater(obs, null)
    o_tf, e_tf, p_tf = overlap(tf_all)
    o_cm, e_cm, p_cm = overlap(cm)
    return dict(median_tf_pct=med_pct, frac_tf_above_median=frac_above,
                tf_obs=o_tf, tf_exp=e_tf, tf_p=float(p_tf),
                cm_obs=o_cm, cm_exp=e_cm, cm_p=float(p_cm),
                n_valid=len(idx), conserved=set(idx[cons].tolist()))

def cons_set_ids(C):
    v = valid & ~np.isnan(C); idx = np.where(v)[0]; Cv = C[idx]
    return set(idx[Cv >= np.quantile(Cv, 0.75)].tolist())

R = {}
print("=" * 74, "\nDONOR-STABILITY GATE\n", "=" * 74)

# ============ VALIDITY GATE ============
print("\n-- VALIDITY: reproduce deposit C + obs/null=0.522 --")
H = Agg("human", 10000); M = Agg("mouse", 10000)
Hfull = H.centroid(); Mfull = M.centroid()
Cfull, nfull = Cvec(Hfull, Mfull)
# verify vectorized Pearson matches gate_lib on deposit centroids
Cdep_check, _ = Cvec(Hdep, m_dep[genes].values)
gate_dep = core["C_pearson"].values
mask = ~np.isnan(Cdep_check) & ~np.isnan(gate_dep)
R["vec_vs_gatelib_maxdiff"] = float(np.nanmax(np.abs(Cdep_check[mask] - gate_dep[mask])))
mask2 = ~np.isnan(Cfull) & ~np.isnan(gate_dep)
R["validity"] = dict(
    corr_Cnew_vs_deposit=float(stats.spearmanr(Cfull[mask2], gate_dep[mask2])[0]),
    pearson_Cnew_vs_deposit=float(np.corrcoef(Cfull[mask2], gate_dep[mask2])[0, 1]),
    obs_null_full=obs_null(Hfull, Mfull),
    n_genes_valid_full=nfull,
    centroid_corr_human=float(np.nanmedian([np.corrcoef(Hfull[t], Hdep[t])[0, 1] for t in range(35)])),
    vec_vs_gatelib_maxdiff=R["vec_vs_gatelib_maxdiff"],
)
print(f"  vectorized-Pearson vs gate_lib max|diff| = {R['vec_vs_gatelib_maxdiff']:.2e}")
print(f"  corr(C_new, C_deposit): spearman={R['validity']['corr_Cnew_vs_deposit']:.3f} "
      f"pearson={R['validity']['pearson_Cnew_vs_deposit']:.3f}")
print(f"  obs/null (full, fresh pull) = {R['validity']['obs_null_full']:.3f} (anchor 0.522)")
print(f"  median per-type human centroid corr new-vs-deposit = {R['validity']['centroid_corr_human']:.3f}")
fin_full = eval_finding(Cfull, ndraw=1000)
R["finding_full"] = {k: v for k, v in fin_full.items() if k != "conserved"}
print(f"  full-data master-TF: median pct={fin_full['median_tf_pct']:.3f}, "
      f"TF enrich {fin_full['tf_obs']}/{fin_full['tf_exp']:.1f} p={fin_full['tf_p']:.3f}, "
      f"CellMarker {fin_full['cm_obs']}/{fin_full['cm_exp']:.1f} p={fin_full['cm_p']:.3f}")

# ============ DONOR-SPLIT (main, cap 10000) ============
def run_donor_split(H, M, nsplit, ndraw=NDRAW, label="cap10000"):
    hd = H.all_donors(); md = M.all_donors()
    cross, jacc, half_recs = [], [], []
    ndrop = []
    for s in range(nsplit):
        rs = np.random.default_rng(100 + s)
        hp = rs.permutation(hd); mp = rs.permutation(md)
        hA, hB = set(hp[:len(hp)//2]), set(hp[len(hp)//2:])
        mA, mB = set(mp[:len(mp)//2]), set(mp[len(mp)//2:])
        CA, _ = Cvec(H.centroid(hA), M.centroid(mA))
        CB, _ = Cvec(H.centroid(hB), M.centroid(mB))
        m = ~np.isnan(CA) & ~np.isnan(CB) & valid
        cross.append(stats.spearmanr(CA[m], CB[m])[0])
        jacc.append(len(cons_set_ids(CA) & cons_set_ids(CB)) /
                    len(cons_set_ids(CA) | cons_set_ids(CB)))
        if s < (nsplit if ndraw else 0):
            for C in (CA, CB):
                f = eval_finding(C, ndraw=ndraw)
                half_recs.append(dict(median_tf_pct=f["median_tf_pct"], tf_p=f["tf_p"],
                                      cm_p=f["cm_p"], tf_fold=f["tf_obs"]/max(f["tf_exp"], 1e-9)))
    hr = pd.DataFrame(half_recs)
    out = dict(
        n_split=nsplit,
        cross_half_C_spearman_median=float(np.median(cross)),
        cross_half_C_spearman_p2_5=float(np.percentile(cross, 2.5)),
        conserved_jaccard_median=float(np.median(jacc)),
        tf_pct_median=float(hr.median_tf_pct.median()),
        tf_pct_min=float(hr.median_tf_pct.min()),
        frac_halves_tf_recovered=float((hr.median_tf_pct > 0.75).mean()),
        frac_halves_tf_enrich_sig=float((hr.tf_p < 0.05).mean()),
        frac_halves_cm_enrich_sig=float((hr.cm_p < 0.05).mean()),
        frac_halves_BOTH=float(((hr.median_tf_pct > 0.75) & (hr.tf_p < 0.05)).mean()),
        tf_fold_median=float(hr.tf_fold.median()),
    )
    return out, cross

print(f"\n-- DONOR-SPLIT (cap10000, {NSPLIT} splits) --")
ds_main, cross_main = run_donor_split(H, M, NSPLIT)
R["donor_split_cap10000"] = ds_main
for k, v in ds_main.items():
    print(f"  {k}: {v:.4f}" if isinstance(v, float) else f"  {k}: {v}")

# ============ CEILING (cell-bootstrap within same donors) ============
print("\n-- CEILING (cell-split within same donors) --")
csh = np.load(HERE / "agg_human_cs.npz"); csm = np.load(HERE / "agg_mouse_cs.npz")
def cs_centroid(cs, cnt, r):
    H = cs[r].astype(float).copy()
    c = cnt[r].astype(float)
    H[c > 0] /= c[c > 0, None]
    H[c == 0] = np.nan
    return H
ceil = []
for r in range(K_CEIL):
    C1, _ = Cvec(cs_centroid(csh["csA"], csh["cntA"], r), cs_centroid(csm["csA"], csm["cntA"], r))
    C2, _ = Cvec(cs_centroid(csh["csB"], csh["cntB"], r), cs_centroid(csm["csB"], csm["cntB"], r))
    m = ~np.isnan(C1) & ~np.isnan(C2) & valid
    ceil.append(stats.spearmanr(C1[m], C2[m])[0])
R["ceiling_cellsplit_spearman_median"] = float(np.median(ceil))
print(f"  ceiling (same donors, diff cells) cross-half C spearman = {np.median(ceil):.4f}")
R["gap_ceiling_minus_donorsplit"] = float(np.median(ceil) - ds_main["cross_half_C_spearman_median"])
print(f"  GAP ceiling - donor-split = {R['gap_ceiling_minus_donorsplit']:.4f} (donor-specific instability)")

# ============ NULL (shuffle gene labels) ============
nb = []
rs = np.random.default_rng(9)
# use the last split's CA/CB analog: recompute one split
hd = H.all_donors(); md = M.all_donors()
hp = rs.permutation(hd); mp = rs.permutation(md)
CA, _ = Cvec(H.centroid(set(hp[:len(hp)//2])), M.centroid(set(mp[:len(mp)//2])))
CB, _ = Cvec(H.centroid(set(hp[len(hp)//2:])), M.centroid(set(mp[len(mp)//2:])))
m = ~np.isnan(CA) & ~np.isnan(CB) & valid
for s in range(200):
    rr = np.random.default_rng(800 + s)
    perm = rr.permutation(np.where(m)[0])
    nb.append(stats.spearmanr(CA[m], CB[perm])[0])
R["null_shuffle_spearman_95"] = float(np.percentile(nb, 95))
print(f"  NULL shuffle-gene cross-half C 95pct = {R['null_shuffle_spearman_95']:.4f} (chance ~0)")

# ============ CAP / POWER ============
print("\n-- CAP / POWER (cross-half C + TF recovery vs cells) --")
R["caps"] = {}
for cap in (500, 2000, 10000):
    Hc = Agg("human", cap); Mc = Agg("mouse", cap)
    o, _ = run_donor_split(Hc, Mc, NSPLIT_CAP, ndraw=NDRAW, label=f"cap{cap}")
    R["caps"][cap] = dict(cross_half=o["cross_half_C_spearman_median"],
                          tf_pct_median=o["tf_pct_median"],
                          frac_halves_BOTH=o["frac_halves_BOTH"])
    print(f"  cap{cap}: cross-half C={o['cross_half_C_spearman_median']:.3f}, "
          f"TF pct median={o['tf_pct_median']:.3f}, halves passing BOTH={o['frac_halves_BOTH']:.2f}")

# ============ CROSS-PROTOCOL (10x vs Smart-seq2) ============
print("\n-- CROSS-PROTOCOL (10x vs Smart-seq2 within atlas) --")
def proto_C(H, M, proto):
    return Cvec(H.centroid(proto=proto), M.centroid(proto=proto))
C10, n10 = proto_C(H, M, "10x")
Css, nss = proto_C(H, M, "SS")
m = ~np.isnan(C10) & ~np.isnan(Css) & valid
R["cross_protocol"] = dict(
    n_types_10x=n10, n_types_SS=nss,
    spearman_C10x_vs_CSS=float(stats.spearmanr(C10[m], Css[m])[0]),
    tf_pct_10x=eval_finding(C10, ndraw=1000)["median_tf_pct"],
    tf_pct_SS=eval_finding(Css, ndraw=1000)["median_tf_pct"],
    tf_p_10x=eval_finding(C10, ndraw=1000)["tf_p"],
    tf_p_SS=eval_finding(Css, ndraw=1000)["tf_p"],
)
cp = R["cross_protocol"]
print(f"  C types covered: 10x={n10}, SS={nss}")
print(f"  Spearman(C_10x, C_SS) = {cp['spearman_C10x_vs_CSS']:.3f}")
print(f"  TF median pct: 10x={cp['tf_pct_10x']:.3f} (p={cp['tf_p_10x']:.3f}); "
      f"SS={cp['tf_pct_SS']:.3f} (p={cp['tf_p_SS']:.3f})")

json.dump(R, open(HERE / "donor_stability_results.json", "w"), indent=2)
print("\nSaved donor_stability_results.json")
