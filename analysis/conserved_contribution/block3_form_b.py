#!/usr/bin/env python3
"""
Out-of-sample test of geometric circularity in the conserved-contribution selection.

Splits donors into halves A and B for both species independently, using
run_donor_stability.py's own construction and seeds, computes C on half A, selects
genes on C_A, and evaluates the published obs/null geometry on half-B centroids.
Selection and evaluation share no donor.

Produces the out-of-sample geometry figures reported in S1 Text, "The objection is
answered instead out of sample": the selected set outperforming the full unselected
gene space in all forty comparisons across twenty splits in both directions, the
margin median 0.121 and minimum 0.092, and the two z values against one hundred
expression-matched random sets drawn on the evaluation half.

NOT TRACKED: reads a capped per-donor aggregate (agg_human_cap10000.npz) that is not
deposited in this repository. gene_conservation_core.csv, run_gate.py and
run_donor_stability.py in this directory are tracked.

Only the two path constants differ from the version that produced those values; the
code is unchanged.
"""
from __future__ import annotations

import ast
import contextlib
import io
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
CC = REPO / "analysis/conserved_contribution"
DONOR = CC / "donor_stability"
OUT = HERE

sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(CC))
import gate_lib as G  # noqa: E402
from cellwarp.procrustes import pca_reduce_centroids, procrustes_align, permutation_test  # noqa: E402

NPERM = 2000
N_MATCHED_DRAWS = 100
N_SPLITS_FULL = 1      # split 0: all five arms
N_SPLITS_TIDF = 20     # splits 0..19: arms T/I/F
SEED_MATCHED = 7       # run_donor_stability.py:25 uses default_rng(7); mirror it


# ---------------------------------------------------------------------------
# AST extraction of the published definitions
# ---------------------------------------------------------------------------
def extract(pyfile: Path, names: set[str], ns: dict) -> dict:
    src = pyfile.read_text()
    tree = ast.parse(src)
    got = {}
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)) and node.name in names:
            seg = ast.get_source_segment(src, node)
            got[node.name] = seg
            exec(compile(ast.Module(body=[node], type_ignores=[]), str(pyfile), "exec"), ns)
    missing = names - set(got)
    if missing:
        raise SystemExit(f"ABORT: {pyfile.name} does not define {sorted(missing)} -- "
                         "the published code does not compose as this check assumes.")
    return got


ns = dict(np=np, pd=pd, stats=stats, io=io, contextlib=contextlib, G=G,
          pca_reduce_centroids=pca_reduce_centroids,
          procrustes_align=procrustes_align, permutation_test=permutation_test)

src_ds = extract(DONOR / "run_donor_stability.py", {"Agg", "Cvec"}, ns)
src_gate = extract(CC / "run_gate.py", {"obs_null_ratio"}, ns)
print("=" * 78)
print("PUBLISHED CODE REUSED VERBATIM (AST-extracted, not retyped)")
print("=" * 78)
for name, seg in {**src_ds, **src_gate}.items():
    print(f"\n----- {name} -----")
    print(seg)
print("=" * 78)

Agg = ns["Agg"]
Cvec = ns["Cvec"]
obs_null_ratio = ns["obs_null_ratio"]

# Agg's methods reference module-level HERE and NG in their original file.
ns["HERE"] = DONOR

# ---------------------------------------------------------------------------
# Frozen gene annotations, exactly as run_donor_stability.py builds them
# ---------------------------------------------------------------------------
core = pd.read_csv(CC / "gene_conservation_core.csv")
genes = list(np.load(DONOR / "agg_human_cap10000.npz", allow_pickle=True)["genes"])
core = core.set_index("gene_id").loc[genes].reset_index()
NG = len(genes)
ns["NG"] = NG
mean_expr = core["mean_expression"].values.astype(float)
deposit_valid = core["C_pearson"].notna().values & np.isfinite(mean_expr)
print(f"\ngenes in aggregate: {NG}")
print(f"deposit-valid genes (C_pearson defined & finite mean_expr): {int(deposit_valid.sum())}")

H = Agg("human", 10000)
M = Agg("mouse", 10000)
hd_all, md_all = H.all_donors(), M.all_donors()
print(f"human donors: {len(hd_all)}   mouse donors: {len(md_all)}")

results = {"config": dict(nperm=NPERM, n_matched_draws=N_MATCHED_DRAWS,
                          n_splits_TIDF=N_SPLITS_TIDF, seed_matched=SEED_MATCHED,
                          n_genes=NG, n_deposit_valid=int(deposit_valid.sum()),
                          n_human_donors=len(hd_all), n_mouse_donors=len(md_all)),
           "splits": []}


def half_frames(hset, mset):
    """(n_ok x NG) half-donor centroids as DataFrames keyed by gene id.

    IMPORTANT composition note. run_gate.py's obs_null_ratio was written against the
    DEPOSIT centroids, which have no NaN rows, so it does not mask incomplete types.
    Half-donor centroids DO have NaN rows: a cell type with no cells in this donor
    half comes back all-NaN from Agg.centroid. run_donor_stability.py's own obs_null
    handles exactly this with

        ok = ~(np.isnan(H).any(1) | np.isnan(M).any(1))

    and Cvec applies the identical mask before computing C. We therefore apply that
    same mask HERE, once per half, before handing the frames to obs_null_ratio --
    so every arm within a split is evaluated on the same set of cell types, and the
    masking rule is the published one rather than a new one.
    """
    Hc = H.centroid(hset)
    Mc = M.centroid(mset)
    ok = ~(np.isnan(Hc).any(1) | np.isnan(Mc).any(1))
    return (pd.DataFrame(Hc[ok], columns=genes),
            pd.DataFrame(Mc[ok], columns=genes),
            int(ok.sum()))


def run_arms(sel_C, eval_h, eval_m, universe_pos, arms, rng, label):
    """Evaluate each arm's gene set on the given (already half-) centroid frames."""
    # bind the published function's module-level centroids to THIS half
    ns["h_cent"] = eval_h
    ns["m_cent"] = eval_m
    out = {}
    for arm_name, gene_ids in arms.items():
        t0 = time.time()
        obs, nullmed, ratio = obs_null_ratio(gene_ids, NPERM)
        out[arm_name] = dict(n_genes=len(gene_ids), obs=obs, null_median=nullmed,
                             ratio=ratio, sec=round(time.time() - t0, 2))
        print(f"    {label} {arm_name:26s} n={len(gene_ids):6d}  obs={obs:9.4f}  "
              f"null_med={nullmed:9.4f}  obs/null={ratio:.4f}  ({out[arm_name]['sec']}s)")
    return out


def one_direction(split, sel_tag, eval_tag, sel_C, eval_C, eval_h, eval_m,
                  universe, full_arms):
    """sel_tag/eval_tag in {'A','B'}; select on sel_C, evaluate geometry on eval_* frames."""
    upos = np.where(universe)[0]
    cs = sel_C[upos]
    ce = eval_C[upos]
    q75_sel, q25_sel = np.quantile(cs, 0.75), np.quantile(cs, 0.25)
    q75_eval = np.quantile(ce, 0.75)

    T_pos = upos[cs >= q75_sel]
    I_pos = upos[ce >= q75_eval]
    D_pos = upos[cs <= q25_sel]

    gid = np.array(genes)
    arms = {
        f"T top-quartile by C_{sel_tag} (OUT-of-sample)": gid[T_pos].tolist(),
        f"I top-quartile by C_{eval_tag} (IN-sample)": gid[I_pos].tolist(),
        "F all valid genes (FLOOR)": gid[upos].tolist(),
    }
    if full_arms:
        arms[f"D bottom-quartile by C_{sel_tag} (divergent)"] = gid[D_pos].tolist()

    rng = np.random.default_rng(SEED_MATCHED)
    res = run_arms(sel_C, eval_h, eval_m, upos, arms, rng,
                   f"split{split} sel={sel_tag}->eval={eval_tag}")

    if full_arms:
        # arm M: expression-matched random, same size as T, published sampler
        me = mean_expr[upos]
        bins = G.expr_bins(me, n_bins=20)
        T_local = np.searchsorted(upos, T_pos)          # positions within the universe
        draws = G.matched_draws(T_local, bins, N_MATCHED_DRAWS, rng)
        mr = []
        t0 = time.time()
        for j, d in enumerate(draws):
            ids = gid[upos[d]].tolist()
            mr.append(obs_null_ratio(ids, NPERM)[2])
        mr = np.array(mr)
        res["M expr-matched random (published sampler)"] = dict(
            n_draws=N_MATCHED_DRAWS, n_genes=int(len(T_pos)),
            mean=float(mr.mean()), sd=float(mr.std(ddof=1)),
            p5=float(np.percentile(mr, 5)), p50=float(np.percentile(mr, 50)),
            p95=float(np.percentile(mr, 95)),
            min=float(mr.min()), max=float(mr.max()),
            ratios=[float(x) for x in mr], sec=round(time.time() - t0, 2))
        m = res["M expr-matched random (published sampler)"]
        print(f"    split{split} sel={sel_tag}->eval={eval_tag} "
              f"M expr-matched random n_draws={N_MATCHED_DRAWS} "
              f"mean={m['mean']:.4f} sd={m['sd']:.4f} p5={m['p5']:.4f} "
              f"[{m['min']:.4f},{m['max']:.4f}] ({m['sec']}s)")
    return res


t_start = time.time()
for split in range(N_SPLITS_TIDF):
    full = split < N_SPLITS_FULL
    rs = np.random.default_rng(100 + split)          # run_donor_stability.py:174
    hp = rs.permutation(hd_all)
    mp = rs.permutation(md_all)
    hA, hB = set(hp[:len(hp) // 2]), set(hp[len(hp) // 2:])
    mA, mB = set(mp[:len(mp) // 2]), set(mp[len(mp) // 2:])

    C_A, nA_types = Cvec(H.centroid(hA), M.centroid(mA))
    C_B, nB_types = Cvec(H.centroid(hB), M.centroid(mB))
    universe = deposit_valid & ~np.isnan(C_A) & ~np.isnan(C_B)

    hA_f, mA_f, nA_ok = half_frames(hA, mA)
    hB_f, mB_f, nB_ok = half_frames(hB, mB)

    print(f"\n=== SPLIT {split}  (arms {'T/I/F/M/D' if full else 'T/I/F'}) ===")
    print(f"  donors: human A={len(hA)} B={len(hB)} | mouse A={len(mA)} B={len(mB)}")
    print(f"  types with complete centroids: A={nA_types} B={nB_types} "
          f"(frames: A={nA_ok} rows, B={nB_ok} rows)")
    print(f"  universe (deposit-valid & C_A & C_B defined) = {int(universe.sum())}")

    rec = dict(split=split,
               n_human_A=len(hA), n_human_B=len(hB), n_mouse_A=len(mA), n_mouse_B=len(mB),
               n_types_A=nA_types, n_types_B=nB_types,
               n_types_frame_A=nA_ok, n_types_frame_B=nB_ok,
               n_universe=int(universe.sum()),
               spearman_CA_CB=float(stats.spearmanr(C_A[universe], C_B[universe])[0]))
    print(f"  spearman(C_A, C_B) on universe = {rec['spearman_CA_CB']:.4f}")

    rec["select_A_eval_B"] = one_direction(split, "A", "B", C_A, C_B, hB_f, mB_f,
                                           universe, full)
    rec["select_B_eval_A"] = one_direction(split, "B", "A", C_B, C_A, hA_f, mA_f,
                                           universe, full)
    results["splits"].append(rec)
    (OUT / "block3_form_b_results.json").write_text(json.dumps(results, indent=2))

results["wall_sec"] = round(time.time() - t_start, 1)
(OUT / "block3_form_b_results.json").write_text(json.dumps(results, indent=2))
print(f"\nwall time: {results['wall_sec']}s")
print(f"wrote {OUT/'block3_form_b_results.json'}")
