"""
Conserved-contribution gene-set gate — runs the pre-registered checks.

Outputs:
  gene_conservation_core.csv        per-gene C_pearson/C_spearman/mean_expression
  gate_results.json                 all numeric results
  (printed)                         structured human-readable report
"""
from __future__ import annotations

import io
import json
import contextlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.mixture import GaussianMixture

import gate_lib as G

rng = np.random.default_rng(G.SEED)
HERE = Path(__file__).resolve().parent
R = {}  # results dict

print("=" * 74)
print("CONSERVED-CONTRIBUTION GATE")
print("=" * 74)

# ---------------------------------------------------------------------------
# Core quantity
# ---------------------------------------------------------------------------
df, h_cent, m_cent = G.build_gene_table()
df.to_csv(HERE / "gene_conservation_core.csv", index=False)
valid = df.dropna(subset=["C_pearson"]).copy().reset_index(drop=True)
C = valid["C_pearson"].values
me = valid["mean_expression"].values
n_valid = len(valid)
R["n_genes_total"] = int(len(df))
R["n_valid"] = int(n_valid)
print(f"\nGenes: {len(df)} total, {n_valid} with defined C_pearson "
      f"({len(df)-n_valid} constant-profile excluded)")
print(f"C_pearson: mean={C.mean():.3f} median={np.median(C):.3f} "
      f"sd={C.std():.3f}  frac>0.5={np.mean(C>0.5):.3f} frac>0.8={np.mean(C>0.8):.3f} "
      f"frac<0={np.mean(C<0):.3f}")
R["C_summary"] = dict(mean=float(C.mean()), median=float(np.median(C)), sd=float(C.std()),
                      frac_gt0_5=float(np.mean(C > 0.5)), frac_gt0_8=float(np.mean(C > 0.8)),
                      frac_lt0=float(np.mean(C < 0)))
# sanity: pearson vs spearman vs (circular) procrustes_contribution
gc = pd.read_csv(G.GENE_CONS_TABLE)[["gene_id", "procrustes_contribution"]]
mtmp = valid.merge(gc, on="gene_id", how="left")
R["rho_pearson_spearman"] = float(stats.spearmanr(mtmp.C_pearson, mtmp.C_spearman, nan_policy="omit")[0])
sub = mtmp.dropna(subset=["procrustes_contribution"])
R["rho_pearson_vs_circular_loading"] = float(stats.spearmanr(sub.C_pearson, sub.procrustes_contribution)[0])
print(f"sanity rho(C_pearson, C_spearman)={R['rho_pearson_spearman']:.3f}; "
      f"rho(C_pearson, circular loading)={R['rho_pearson_vs_circular_loading']:.3f}")

# expression-matched bins (computed once on valid genes)
bins = G.expr_bins(me, n_bins=20)
all_idx = np.arange(n_valid)

# conserved / divergent sets
q75, q25 = np.quantile(C, 0.75), np.quantile(C, 0.25)
cons_idx = all_idx[C >= q75]
div_idx = all_idx[C <= q25]
strict_idx = all_idx[C > 0.8]
R["thresholds"] = dict(q75=float(q75), q25=float(q25),
                       n_conserved=int(len(cons_idx)), n_divergent=int(len(div_idx)),
                       n_strict_gt0_8=int(len(strict_idx)))
print(f"\nConserved set (C>=Q75={q75:.3f}): {len(cons_idx)} | "
      f"Divergent (C<=Q25={q25:.3f}): {len(div_idx)} | strict C>0.8: {len(strict_idx)}")

# ---------------------------------------------------------------------------
# CHECK 1 — distribution structure
# ---------------------------------------------------------------------------
print("\n" + "-" * 74 + "\nCHECK 1 — distribution structure")
# bimodality coefficient (Sarle): (skew^2 + 1) / kurtosis(+3) ; >0.555 suggests bimodal
skew = stats.skew(C); kurt = stats.kurtosis(C, fisher=True)
bc = (skew**2 + 1) / (kurt + 3)
# GMM BIC, 1 vs 2 vs 3 components
X = C.reshape(-1, 1)
bics = {}
for k in (1, 2, 3):
    gm = GaussianMixture(n_components=k, covariance_type="full", random_state=G.SEED, n_init=3).fit(X)
    bics[k] = float(gm.bic(X))
best_k = min(bics, key=bics.get)
# dip test if available
dip_p = None
try:
    import diptest
    _, dip_p = diptest.diptest(C)
except Exception:
    dip_p = None
R["check1"] = dict(skew=float(skew), kurtosis=float(kurt), bimodality_coef=float(bc),
                   gmm_bic=bics, gmm_best_k=int(best_k), dip_p=(float(dip_p) if dip_p is not None else None))
print(f"  skew={skew:.3f} kurtosis={kurt:.3f} bimodality_coef={bc:.3f} (>0.555 ~ multimodal)")
print(f"  GMM BIC: {{1:{bics[1]:.0f}, 2:{bics[2]:.0f}, 3:{bics[3]:.0f}}} -> best k={best_k}")
print(f"  Hartigan dip p = {dip_p}")
structured = (bc > 0.555) or (best_k >= 2 and (bics[1] - bics[best_k]) > 10) or (dip_p is not None and dip_p < 0.05)
R["check1"]["structured"] = bool(structured)
print(f"  => {'STRUCTURED (separable)' if structured else 'SMOOTH CONTINUUM (thresholds arbitrary; weakens)'}")

# ---------------------------------------------------------------------------
# CHECK 2 — expression independence
# ---------------------------------------------------------------------------
print("\n" + "-" * 74 + "\nCHECK 2 — expression independence")
rho_expr = stats.spearmanr(C, me)[0]
# decile-of-expression means of C
dec = pd.qcut(pd.Series(me).rank(method="first"), 10, labels=False)
dec_means = [float(C[dec == d].mean()) for d in range(10)]
R["check2"] = dict(spearman_C_vs_expr=float(rho_expr), decile_C_means=dec_means)
if abs(rho_expr) < 0.3:
    band = "PASS (<0.3: not primarily an expression proxy)"
elif abs(rho_expr) < 0.5:
    band = "PARTIAL (0.3-0.5: must clear check3 vs expr-matched bg)"
else:
    band = "FAIL (>=0.5: conservation ~ expression)"
R["check2"]["band"] = band
print(f"  Spearman(C, mean_expr) = {rho_expr:.3f} -> {band}")
print(f"  C by expression decile (low->high): " + " ".join(f"{x:.2f}" for x in dec_means))

# ---------------------------------------------------------------------------
# Build reference gene sets (symbols -> valid-table positional idx)
# ---------------------------------------------------------------------------
sym2ens, _ = G.ortholog_maps()
sym_to_pos = {s: i for i, s in enumerate(valid["symbol"].values)}
ens_to_pos = {g: i for i, g in enumerate(valid["gene_id"].values)}

def syms_to_pos(symbols):
    pos = []
    for s in symbols:
        if s in sym_to_pos:
            pos.append(sym_to_pos[s])
    return np.array(sorted(set(pos)), dtype=int)

# positive-control TFs present & valid
tf_pc_pos = syms_to_pos(G.POSITIVE_CONTROL_TFS)
R["tf_pc_total"] = len(G.POSITIVE_CONTROL_TFS)
R["tf_pc_testable"] = int(len(tf_pc_pos))

# local CollecTRI TF list (296)
tf_all = pd.read_csv(G.TF_ACTIVITY, index_col=0).columns.tolist()
tf_all_pos = syms_to_pos(tf_all)
# CellMarker identity markers
cm = pd.read_csv(G.CELLMARKER_H)
cm_syms = cm["gene_symbol"].dropna().unique().tolist()
cm_pos = syms_to_pos(cm_syms)

# percentile rank of C (0..1)
C_rank = stats.rankdata(C) / n_valid

# ---------------------------------------------------------------------------
# CHECK 3a — TF positive control (vs expression-matched background)
# ---------------------------------------------------------------------------
print("\n" + "-" * 74 + "\nCHECK 3a — identity-TF positive control")
def matched_stat_pvalue(target_pos, stat_fn, n=1000):
    obs = stat_fn(target_pos)
    draws = G.matched_draws(target_pos, bins, n, rng)
    null = np.array([stat_fn(d) for d in draws])
    return obs, null

median_rank = lambda pos: float(np.median(C_rank[pos]))
obs_mr, null_mr = matched_stat_pvalue(tf_pc_pos, median_rank, 1000)
p_3a = G.emp_p_greater(obs_mr, null_mr)
frac_above_median = float(np.mean(C[tf_pc_pos] > np.median(C)))
frac_in_conserved = float(np.mean(np.isin(tf_pc_pos, cons_idx)))
R["check3a"] = dict(testable=int(len(tf_pc_pos)), median_Crank=obs_mr,
                    null_median_Crank=float(np.median(null_mr)), p_value=float(p_3a),
                    frac_above_global_median=frac_above_median,
                    frac_in_conserved_set=frac_in_conserved, enrichment_vs_0_25=float(frac_in_conserved/0.25))
pass_3a = (p_3a < 0.05) and (frac_above_median > 0.5)
R["check3a"]["pass"] = bool(pass_3a)
print(f"  testable identity TFs: {len(tf_pc_pos)}/{len(G.POSITIVE_CONTROL_TFS)}")
print(f"  median C-percentile of TF list = {obs_mr:.3f} (matched-bg median {np.median(null_mr):.3f}), "
      f"empirical p = {p_3a:.4f}")
print(f"  frac above global median C = {frac_above_median:.2f}; "
      f"frac in conserved set = {frac_in_conserved:.2f} (vs 0.25 expected; {frac_in_conserved/0.25:.2f}x)")
print(f"  => 3a {'PASS' if pass_3a else 'FAIL'}")

# ---------------------------------------------------------------------------
# CHECK 3b — coherent enrichment vs expression-matched background
# ---------------------------------------------------------------------------
print("\n" + "-" * 74 + "\nCHECK 3b — coherent enrichment (expr-matched)")
def overlap_enrichment(set_pos, ref_pos, n=1000):
    ref_set = set(ref_pos.tolist())
    obs = int(np.sum([p in ref_set for p in set_pos]))
    draws = G.matched_draws(set_pos, bins, n, rng)
    null = np.array([int(np.sum([p in ref_set for p in d])) for d in draws])
    p = G.emp_p_greater(obs, null)
    return obs, float(np.mean(null)), float(p)

# H3b-i: conserved set enriched for TFs (CollecTRI 296)
obs_tf, exp_tf, p_tf = overlap_enrichment(cons_idx, tf_all_pos, 1000)
# H3b-ii: conserved set enriched for CellMarker identity markers
obs_cm, exp_cm, p_cm = overlap_enrichment(cons_idx, cm_pos, 1000)
# BH-FDR across the two
from statsmodels.stats.multitest import multipletests
ps = [p_tf, p_cm]
rej, padj, _, _ = multipletests(ps, method="fdr_bh")
R["check3b"] = dict(
    H_tf=dict(ref_n=int(len(tf_all_pos)), observed=obs_tf, expected_matched=exp_tf,
              fold=float(obs_tf/exp_tf if exp_tf else np.nan), p=float(p_tf), p_adj=float(padj[0])),
    H_cellmarker=dict(ref_n=int(len(cm_pos)), observed=obs_cm, expected_matched=exp_cm,
                      fold=float(obs_cm/exp_cm if exp_cm else np.nan), p=float(p_cm), p_adj=float(padj[1])),
)
pass_3b = bool(((padj[0] < 0.05 and obs_tf > exp_tf)) or ((padj[1] < 0.05 and obs_cm > exp_cm)))
R["check3b"]["pass"] = pass_3b
print(f"  H3b-i  TF (CollecTRI {len(tf_all_pos)} valid): conserved has {obs_tf} vs matched {exp_tf:.1f} "
      f"({obs_tf/exp_tf:.2f}x) p={p_tf:.4f} padj={padj[0]:.4f}")
print(f"  H3b-ii CellMarker ({len(cm_pos)} valid): conserved has {obs_cm} vs matched {exp_cm:.1f} "
      f"({obs_cm/exp_cm:.2f}x) p={p_cm:.4f} padj={padj[1]:.4f}")
print(f"  => 3b {'PASS' if pass_3b else 'FAIL'}")

# ---------------------------------------------------------------------------
# CHECK 4 — membership stability (type-level fallback)
# ---------------------------------------------------------------------------
print("\n" + "-" * 74 + "\nCHECK 4 — membership stability (type-level; raw cells absent)")
H_full, M_full = h_cent.values, m_cent.values
types = list(h_cent.index)
gene_ids_full = list(h_cent.columns)
base_pos_ids = set(valid.iloc[cons_idx]["gene_id"])

def conserved_ids_from_types(keep_type_idx):
    Hs, Ms = H_full[keep_type_idx], M_full[keep_type_idx]
    c = G.per_gene_corr(Hs, Ms, "pearson")
    ok = ~np.isnan(c)
    ids = np.array(gene_ids_full)[ok]
    cc = c[ok]
    thr = np.quantile(cc, 0.75)
    return set(ids[cc >= thr]), dict(zip(ids, cc))

def jaccard(a, b):
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if (a | b) else 0.0

# leave-1-type-out
j1 = []
for i in range(len(types)):
    keep = [k for k in range(len(types)) if k != i]
    s, _ = conserved_ids_from_types(keep)
    j1.append(jaccard(base_pos_ids, s))
# leave-5-types-out (100 random)
j5 = []
for _ in range(100):
    keep = sorted(rng.choice(len(types), len(types) - 5, replace=False))
    s, _ = conserved_ids_from_types(keep)
    j5.append(jaccard(base_pos_ids, s))
# cell-count restriction: types with min_count >= 2000
inv = pd.read_csv(ROOT_INV := G.ROOT / "output/phase2/cell_type_inventory_passing.csv")
hi = set(inv[inv["min_count"] >= 2000]["cell_type"])
keep_hi = [k for k, t in enumerate(types) if t in hi]
s_hi, _ = conserved_ids_from_types(keep_hi)
j_hi = jaccard(base_pos_ids, s_hi)
# pearson vs spearman membership
csp = valid["C_spearman"].values
sp_cons = set(valid.iloc[all_idx[csp >= np.nanquantile(csp, 0.75)]]["gene_id"])
j_ps = jaccard(base_pos_ids, sp_cons)
R["check4"] = dict(
    leave1_jaccard_median=float(np.median(j1)), leave1_jaccard_min=float(np.min(j1)),
    leave5_jaccard_median=float(np.median(j5)), leave5_jaccard_p05=float(np.percentile(j5, 5)),
    highcount_n_types=int(len(keep_hi)), highcount_jaccard=float(j_hi),
    pearson_vs_spearman_jaccard=float(j_ps),
)
pass_4 = (np.median(j1) >= 0.60) and (np.median(j5) >= 0.50)
R["check4"]["pass"] = bool(pass_4)
print(f"  leave-1-type-out Jaccard: median={np.median(j1):.3f} min={np.min(j1):.3f}")
print(f"  leave-5-type-out Jaccard: median={np.median(j5):.3f} p5={np.percentile(j5,5):.3f}")
print(f"  high-count types only (min>=2000, n={len(keep_hi)}): Jaccard={j_hi:.3f}")
print(f"  Pearson-vs-Spearman conserved-set Jaccard={j_ps:.3f}")
print(f"  => 4 {'PASS' if pass_4 else 'FAIL'} (type-level stability)")

# ---------------------------------------------------------------------------
# SECONDARY — geometry attribution (obs/null) — reuse cellwarp routines
# ---------------------------------------------------------------------------
print("\n" + "-" * 74 + "\nSECONDARY — geometry attribution (caveated; selection-biased)")
from cellwarp.procrustes import pca_reduce_centroids, procrustes_align, permutation_test

def obs_null_ratio(gene_id_list, n_perm=2000):
    cols = [g for g in gene_id_list if g in h_cent.columns]
    hc = h_cent[cols]; mc = m_cent[cols]
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        hp, mp, _, _ = pca_reduce_centroids(hc, mc, 0.95)
        res = procrustes_align(hp, mp)
        _, null = permutation_test(hp, mp, n_perm)
    return float(res.distance), float(np.median(null)), float(res.distance / np.median(null))

cons_ids = valid.iloc[cons_idx]["gene_id"].tolist()
div_ids = valid.iloc[div_idx]["gene_id"].tolist()
# validity gate: all genes -> should reproduce ~0.522
allg_obs, allg_nm, allg_ratio = obs_null_ratio(gene_ids_full)
cons_obs, cons_nm, cons_ratio = obs_null_ratio(cons_ids)
div_obs, div_nm, div_ratio = obs_null_ratio(div_ids)
# expression-matched random sets of size = conserved
mr_ratios = []
for d in G.matched_draws(cons_idx, bins, 20, rng):
    ids = valid.iloc[d]["gene_id"].tolist()
    mr_ratios.append(obs_null_ratio(ids)[2])
R["secondary"] = dict(
    validity_all_genes_ratio=allg_ratio,
    conserved=dict(obs=cons_obs, null_median=cons_nm, ratio=cons_ratio),
    divergent=dict(obs=div_obs, null_median=div_nm, ratio=div_ratio),
    matched_random_ratio_mean=float(np.mean(mr_ratios)),
    matched_random_ratio_sd=float(np.std(mr_ratios)),
    matched_random_ratios=[float(x) for x in mr_ratios],
)
print(f"  validity (all {len(gene_ids_full)} genes) obs/null = {allg_ratio:.3f} (target 0.522)")
print(f"  CONSERVED set obs/null = {cons_ratio:.3f}")
print(f"  DIVERGENT set obs/null = {div_ratio:.3f}")
print(f"  EXPR-MATCHED random (n=20) obs/null = {np.mean(mr_ratios):.3f} +/- {np.std(mr_ratios):.3f}")
print("  [caveat] low conserved obs/null is partly expected by selection; informative contrast")
print("   is conserved vs divergent, and conserved vs expr-matched-random.")

# ---------------------------------------------------------------------------
# VERDICT
# ---------------------------------------------------------------------------
print("\n" + "=" * 74 + "\nVERDICT")
c2_fail = abs(rho_expr) >= 0.5
verdict = "PASS" if (not c2_fail and pass_3a and pass_3b and pass_4) else "FAIL"
if verdict == "PASS" and not structured:
    verdict = "PASS (note: smooth continuum in check 1)"
R["verdict"] = dict(check2_fail=bool(c2_fail), pass_3a=bool(pass_3a), pass_3b=bool(pass_3b),
                    pass_4=bool(pass_4), structured=bool(structured), verdict=verdict)
print(f"  check2 expression-proxy fail: {c2_fail}")
print(f"  3a identity-TF recovery: {'PASS' if pass_3a else 'FAIL'}")
print(f"  3b coherent enrichment:  {'PASS' if pass_3b else 'FAIL'}")
print(f"  4  membership stability: {'PASS' if pass_4 else 'FAIL'}")
print(f"  ==> {verdict}")

with open(HERE / "gate_results.json", "w") as f:
    json.dump(R, f, indent=2)
print(f"\nSaved {HERE/'gate_results.json'}")
