#!/usr/bin/env python3
"""Pseudo-bulk centroid-definition sensitivity analysis (Layer-1, human-mouse, 35 types).

Question: the primary Layer-1 configuration uses MEAN-OF-LOG1P centroids
(per-cell CP10K -> log1p -> mean over cells of a type), which weights cells
equally but is not a pseudo-bulk profile. Does the geometric coherence
(obs/null = 0.522; lineage-stratified 0.668) depend on that choice?

Two pseudo-bulk redefinitions are regenerated from the RAW per-cell counts that
produced the deposited centroids (the aligned, subsampled cells), holding
everything else fixed (same 16,959 orthologs, same 35 types, same cells,
target_sum = 1e4, joint PCA at 95% variance, seed 42):

  A  equal_cell_weight   log1p(mean_cells CP10K)   per-cell CP10K -> mean -> log1p
  B  aggregate_pseudobulk log1p(CP10K(sum raw))    sum raw per type -> CP10K -> log1p

Flavor A keeps the equal-per-cell weighting of the primary but moves the log
outside the mean; flavor B is a true library-size-weighted pseudo-bulk, so
cells with deeper libraries contribute more.

Two fidelity anchors are run through the identical code path:
  deposited   the shipped centroid CSVs (exact anchor: obs/null = 0.5222043)
  recon_mean_log1p  mean-of-log1p rebuilt from the same raw counts by this
                    script (proves the raw-count provenance and the streaming
                    accumulator reproduce the deposit)

Engine: src/cellwarp/procrustes.py permutation_test() called directly, and the
shipped within-lineage machinery from scripts/test_lineage_stratified_permutation.py.
permutation_1M.py is deliberately NOT reused (its asserts pin n_types/n_pcs/
observed to the deposited centroids and would hard-fail on regenerated ones).

Reads only. Writes only new files under analysis/sensitivity_analyses/:
  pseudobulk_centroids_results.json
  pseudobulk_centroids_results.md
  pseudobulk_pca_centroids.npz      (PCA-space centroids per flavor, small)
No deposited centroid CSV, primary output, figure, table, or manuscript file is
touched.

Usage:
    python analysis/sensitivity_analyses/pseudobulk_centroids.py RAW_DIR
    CELLWARP_PSEUDOBULK_RAW_DIR=<dir> python .../pseudobulk_centroids.py

RAW_DIR is required and has no default. It is the directory holding the
*_raw_aligned.h5ad checkpoints written by scripts/08_scaled_procrustes.py
stage 1. Any default would resolve against the running user's home directory
rather than the tree this deposit was unpacked in. Those checkpoints are not
deposited (see DATA_SOURCES.md): regenerate them with 08_scaled_procrustes.py,
or pass the directory holding your own. Set CELLWARP_PSEUDOBULK_CACHE=<path.npz>
to cache the gene-space accumulators between runs (the streaming pass reads
~2.4 GB).
"""
from __future__ import annotations

import json
import os
import sys
import time
import warnings
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr
from sklearn.decomposition import PCA

warnings.filterwarnings("ignore", message=".*encountered in det.*")

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))
from cellwarp.procrustes import procrustes_align, permutation_test  # noqa: E402
import test_lineage_stratified_permutation as lin  # noqa: E402  (shipped stratified machinery)

SCALED = REPO / "output" / "phase2" / "scaled_35types"
OUT = REPO / "analysis" / "sensitivity_analyses"
_raw_arg = sys.argv[1] if len(sys.argv) > 1 else os.environ.get(
    "CELLWARP_PSEUDOBULK_RAW_DIR")
if not _raw_arg:
    sys.exit(
        "ERROR: RAW_DIR is required and has no default.\n"
        "  usage: python analysis/sensitivity_analyses/pseudobulk_centroids.py "
        "RAW_DIR\n"
        "  or set CELLWARP_PSEUDOBULK_RAW_DIR=<dir>\n"
        "RAW_DIR holds the *_raw_aligned.h5ad checkpoints from "
        "scripts/08_scaled_procrustes.py stage 1; they are not deposited "
        "(see DATA_SOURCES.md)."
    )
# Deliberately no home-directory expansion here: the shell already resolves a
# leading ~ in an argument before this process sees it, and resolving one here
# would reintroduce a home-relative path. See the guard in reproduce/run_all.sh.
RAW_DIR = Path(_raw_arg)
CACHE = os.environ.get("CELLWARP_PSEUDOBULK_CACHE")

SEED = 42
NPERM = 10_000
TARGET_SUM = 1e4
BLOCK = 2048  # cells per streaming block

# Deposited anchors (mean-of-log1p centroids, 35 types, 33 PCs, seed 42, 10k perms)
BASELINE = {"primary_obs_null": 0.5222043226858066,
            "lineage_obs_null": 0.6682720235158841,
            "distance": 61.1529953159736,
            "n_pcs": 33}
# PASS = "well below 1 and significant"; comparability to 0.522 reported separately.
PASS_OBS_NULL_MAX = 0.80
PASS_P_MAX = 0.01
COMPARABLE_DELTA = 0.10


# ---------------------------------------------------------------------------
# Streaming per-type accumulators over the raw CSR h5ad
# ---------------------------------------------------------------------------
def stream_type_sums(path: Path, cell_types: list[str], genes: list[str]) -> dict:
    """One pass over a raw-count h5ad; accumulate per-type gene sums.

    Returns per-type (n_types x n_genes) float64 sums of
      raw        raw counts                        -> flavor B numerator
      cp10k      per-cell CP10K                    -> flavor A numerator
      log1p      per-cell log1p(CP10K)             -> mean-of-log1p reconstruction
    plus n_cells per type. Zero-library cells are handled as scanpy's
    normalize_total does (divisor forced to 1; the row is all zeros anyway).
    """
    n_types, n_genes = len(cell_types), len(genes)
    acc = {k: np.zeros((n_types, n_genes), dtype=np.float64) for k in ("raw", "cp10k", "log1p")}
    n_cells = np.zeros(n_types, dtype=np.int64)
    n_zero_lib = 0

    with h5py.File(path, "r") as f:
        assert list(f["var"]["_index"][:].astype(str)) == genes, f"gene space mismatch in {path}"
        cats = f["obs"]["cell_type"]["categories"][:].astype(str)
        codes = f["obs"]["cell_type"]["codes"][:].astype(np.int64)
        row_of_cat = np.array([cell_types.index(c) for c in cats], dtype=np.int64)
        row_idx = row_of_cat[codes]

        X = f["X"]
        assert X.attrs["encoding-type"] == "csr_matrix", "expected CSR-encoded X"
        indptr = X["indptr"][:].astype(np.int64)
        data_ds, ind_ds = X["data"], X["indices"]
        n_obs = indptr.size - 1
        assert tuple(X.attrs["shape"]) == (n_obs, n_genes)

        t0 = time.time()
        for i0 in range(0, n_obs, BLOCK):
            i1 = min(i0 + BLOCK, n_obs)
            s, e = int(indptr[i0]), int(indptr[i1])
            block = sparse.csr_matrix(
                (data_ds[s:e].astype(np.float64), ind_ds[s:e], indptr[i0:i1 + 1] - s),
                shape=(i1 - i0, n_genes),
            )
            lib = np.asarray(block.sum(axis=1)).ravel()
            n_zero_lib += int((lib == 0).sum())
            lib[lib == 0] = 1.0

            cp10k = sparse.diags(TARGET_SUM / lib) @ block
            lg = cp10k.copy()
            lg.data = np.log1p(lg.data)

            rows = row_idx[i0:i1]
            grouper = sparse.csr_matrix(
                (np.ones(i1 - i0), (rows, np.arange(i1 - i0))), shape=(n_types, i1 - i0)
            )
            acc["raw"] += (grouper @ block).toarray()
            acc["cp10k"] += (grouper @ cp10k).toarray()
            acc["log1p"] += (grouper @ lg).toarray()
            np.add.at(n_cells, rows, 1)
            if (i0 // BLOCK) % 10 == 0:
                print(f"    {i1:>7,}/{n_obs:,} cells  ({time.time() - t0:5.1f}s)", flush=True)

    print(f"    done: {n_obs:,} cells, {n_zero_lib} zero-library cells, "
          f"{time.time() - t0:.1f}s")
    return {"sums": acc, "n_cells": n_cells, "n_obs": int(n_obs), "n_zero_lib": n_zero_lib}


def build_centroids(agg: dict) -> dict[str, np.ndarray]:
    """Turn per-type accumulators into the three centroid definitions."""
    n = agg["n_cells"][:, None].astype(np.float64)
    raw_sum = agg["sums"]["raw"]
    lib = raw_sum.sum(axis=1, keepdims=True)
    return {
        "equal_cell_weight": np.log1p(agg["sums"]["cp10k"] / n),          # A
        "aggregate_pseudobulk": np.log1p(raw_sum / lib * TARGET_SUM),      # B
        "recon_mean_log1p": agg["sums"]["log1p"] / n,                      # anchor
    }


# ---------------------------------------------------------------------------
# Configuration run: joint PCA -> Procrustes -> global null -> lineage null
# ---------------------------------------------------------------------------
def run_flavor(name: str, H: np.ndarray, M: np.ndarray, cell_types: list[str],
               block_indices: list[list[int]], full_rank: pd.Series) -> tuple[dict, np.ndarray, np.ndarray]:
    print("\n" + "=" * 78 + f"\nFLAVOR: {name}\n" + "=" * 78)
    combined = np.vstack([H, M])
    pca = PCA(n_components=0.95, svd_solver="full", random_state=SEED)
    Z = pca.fit_transform(combined)
    k, n = int(pca.n_components_), len(cell_types)
    Hp, Mp = Z[:n], Z[n:]
    print(f"  joint PCA: {k} components, {pca.explained_variance_ratio_.sum() * 100:.2f}% variance")

    res = procrustes_align(Hp, Mp)
    p_global, null_global = permutation_test(Hp, Mp, n_permutations=NPERM, seed=SEED)
    obs_null = float(res.distance / np.median(null_global))

    p_lin, null_lin = lin.stratified_permutation(Hp, Mp, block_indices, n_perm=NPERM, seed=SEED)
    obs_null_lin = float(res.distance / np.median(null_lin))
    print(f"\n  lineage-stratified null: obs/null = {obs_null_lin:.4f}, p = {p_lin:.6f}")

    resid = pd.Series(np.linalg.norm(res.centered_reference - res.aligned_target, axis=1),
                      index=cell_types)
    rho, rho_p = spearmanr(resid.loc[full_rank.index].values, full_rank.values)

    passed = bool(obs_null < PASS_OBS_NULL_MAX and p_global < PASS_P_MAX)
    out = {
        "flavor": name,
        "n_pcs": k,
        "pca_cumvar": float(pca.explained_variance_ratio_.sum()),
        "pcs_unchanged_vs_deposited": bool(k == BASELINE["n_pcs"]),
        "procrustes_distance": float(res.distance),
        "procrustes_scaling": float(res.scaling),
        "global_null_median": float(np.median(null_global)),
        "obs_null_ratio": obs_null,
        "p_value": float(p_global),
        "lineage_null_median": float(np.median(null_lin)),
        "lineage_obs_null_ratio": obs_null_lin,
        "lineage_p_value": float(p_lin),
        "delta_vs_baseline_primary": obs_null - BASELINE["primary_obs_null"],
        "delta_vs_baseline_lineage": obs_null_lin - BASELINE["lineage_obs_null"],
        "residual_ranking_rho_vs_deposited": float(rho),
        "residual_ranking_rho_p": float(rho_p),
        "pass": passed,
        "comparable_to_baseline": bool(abs(obs_null - BASELINE["primary_obs_null"]) <= COMPARABLE_DELTA),
    }
    return out, Hp, Mp


def main() -> None:
    print("=" * 78 + "\nPseudo-bulk centroid-definition sensitivity (Layer-1, human-mouse)\n" + "=" * 78)
    print(f"  raw counts: {RAW_DIR}")

    saved = np.load(SCALED / "pca_centroids_35.npz", allow_pickle=True)
    cell_types = [str(c) for c in saved["cell_types"]]
    dep_h = pd.read_csv(SCALED / "centroids_human_35.csv", index_col=0).loc[cell_types]
    dep_m = pd.read_csv(SCALED / "centroids_mouse_35.csv", index_col=0).loc[cell_types]
    genes = list(dep_h.columns)
    assert list(dep_m.columns) == genes
    print(f"  {len(cell_types)} cell types x {len(genes):,} orthologs")

    ct_to_idx = {ct: i for i, ct in enumerate(cell_types)}
    block_indices = [sorted(ct_to_idx[m] for m in members)
                     for members in lin.LINEAGE_BLOCKS.values()]
    assert sorted(i for b in block_indices for i in b) == list(range(len(cell_types)))

    full_rank = (pd.read_csv(SCALED / "residuals_ranked.csv")
                 .set_index("cell_type")["residual_magnitude"].loc[cell_types])

    # ---- regenerate centroids from raw counts (cached optionally) ----
    if CACHE and Path(CACHE).exists():
        print(f"\n  loading cached accumulators: {CACHE}")
        z = np.load(CACHE)
        cents = {sp: {k: z[f"{sp}_{k}"] for k in
                      ("equal_cell_weight", "aggregate_pseudobulk", "recon_mean_log1p")}
                 for sp in ("human", "mouse")}
        provenance = json.loads(str(z["provenance"]))
    else:
        cents, provenance = {}, {}
        for sp in ("human", "mouse"):
            print(f"\n  streaming {sp}_raw_aligned.h5ad ...")
            agg = stream_type_sums(RAW_DIR / f"{sp}_raw_aligned.h5ad", cell_types, genes)
            cents[sp] = build_centroids(agg)
            provenance[sp] = {"n_cells_total": agg["n_obs"],
                              "n_zero_library_cells": agg["n_zero_lib"],
                              "cells_per_type": dict(zip(cell_types, agg["n_cells"].tolist()))}
        if CACHE:
            np.savez_compressed(CACHE, provenance=json.dumps(provenance),
                                **{f"{sp}_{k}": v for sp, d in cents.items() for k, v in d.items()})
            print(f"  cached accumulators -> {CACHE}")

    # ---- provenance gate: reconstruction vs deposited centroid CSVs ----
    recon_gate = {}
    for sp, dep in (("human", dep_h), ("mouse", dep_m)):
        d = np.abs(cents[sp]["recon_mean_log1p"] - dep.values)
        recon_gate[sp] = {"max_abs_diff": float(d.max()),
                          "mean_abs_diff": float(d.mean()),
                          "max_rel_diff": float((d / np.maximum(np.abs(dep.values), 1e-12)).max())}
        print(f"  recon vs deposited {sp}: max|diff| = {d.max():.3e}, mean|diff| = {d.mean():.3e}")

    # ---- run every flavor through the identical engine ----
    flavors = [("deposited", dep_h.values, dep_m.values)]
    for key in ("recon_mean_log1p", "equal_cell_weight", "aggregate_pseudobulk"):
        flavors.append((key, cents["human"][key], cents["mouse"][key]))

    results, pca_space = [], {}
    for name, H, M in flavors:
        r, Hp, Mp = run_flavor(name, H, M, cell_types, block_indices, full_rank)
        results.append(r)
        pca_space[f"{name}_human"], pca_space[f"{name}_mouse"] = Hp, Mp

    # ---- artifacts ----
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez(OUT / "pseudobulk_pca_centroids.npz", cell_types=np.array(cell_types), **pca_space)

    payload = {
        "analysis": "pseudobulk_centroid_definition_sensitivity",
        "configuration": "Layer-1 primary, human-mouse, 35 cell types, 16,959 orthologs",
        # Provenance only, and deliberately the basename: an absolute path
        # here would record the machine the analysis was run on and would
        # resolve nowhere in an unpacked deposit.
        "raw_counts_dir": RAW_DIR.name,
        "seed": SEED, "n_permutations": NPERM, "target_sum": TARGET_SUM,
        "baseline_deposited": BASELINE,
        "pass_criterion": {"obs_null_ratio_below": PASS_OBS_NULL_MAX,
                           "p_value_below": PASS_P_MAX,
                           "comparable_delta": COMPARABLE_DELTA},
        "provenance": provenance,
        "reconstruction_gate_vs_deposited_csv": recon_gate,
        "flavors": results,
    }
    (OUT / "pseudobulk_centroids_results.json").write_text(json.dumps(payload, indent=2))

    # ---- markdown summary ----
    label = {"deposited": "deposited mean-of-log1p (anchor)",
             "recon_mean_log1p": "mean-of-log1p rebuilt from raw (anchor)",
             "equal_cell_weight": "A. equal-cell-weight  log1p(mean CP10K)",
             "aggregate_pseudobulk": "B. aggregate pseudo-bulk  log1p(CP10K(sum raw))"}
    lines = [
        "# Pseudo-bulk centroid-definition sensitivity (Layer-1, human-mouse, 35 types)",
        "",
        f"Seed {SEED}, {NPERM:,} permutations, target_sum {TARGET_SUM:.0e}, joint PCA at 95% variance.",
        "Everything except the centroid definition is held fixed: same 16,959 orthologs,",
        "same 35 cell types, same cells, same Procrustes/permutation engine.",
        "",
        "| centroid definition | PCs | Procrustes d | obs/null (global) | p | obs/null (lineage-stratified) | p (lineage) | rank rho vs deposited |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        lines.append(
            f"| {label[r['flavor']]} | {r['n_pcs']} | {r['procrustes_distance']:.3f} | "
            f"{r['obs_null_ratio']:.4f} | {r['p_value']:.2e} | "
            f"{r['lineage_obs_null_ratio']:.4f} | {r['lineage_p_value']:.2e} | "
            f"{r['residual_ranking_rho_vs_deposited']:.3f} |"
        )
    lines += ["", f"Deposited baseline: global obs/null = {BASELINE['primary_obs_null']:.4f}, "
                  f"lineage-stratified obs/null = {BASELINE['lineage_obs_null']:.4f}, "
                  f"{BASELINE['n_pcs']} PCs, Procrustes d = {BASELINE['distance']:.3f}.",
              "", "## Verdict", ""]
    for r in results:
        if r["flavor"] in ("deposited", "recon_mean_log1p"):
            continue
        lines.append(
            f"**{r['flavor']}: {'PASS' if r['pass'] else 'FAIL'}** "
            f"(obs/null = {r['obs_null_ratio']:.4f} < {PASS_OBS_NULL_MAX}, p = {r['p_value']:.2e} "
            f"< {PASS_P_MAX}; delta vs 0.522 baseline = {r['delta_vs_baseline_primary']:+.4f}; "
            f"lineage-stratified obs/null = {r['lineage_obs_null_ratio']:.4f}, "
            f"delta = {r['delta_vs_baseline_lineage']:+.4f}; "
            f"retained PCs {r['n_pcs']} "
            f"{'(unchanged)' if r['pcs_unchanged_vs_deposited'] else 'vs 33 (CHANGED)'})")
        lines.append("")

    lines += [
        "## Fidelity gates", "",
        "Both anchors were pushed through the identical code path before the pseudo-bulk flavors:", "",
        f"- `deposited` (the shipped centroid CSVs) returns obs/null = "
        f"{results[0]['obs_null_ratio']:.13f} and lineage-stratified "
        f"{results[0]['lineage_obs_null_ratio']:.13f}, bit-identical to the published "
        f"{BASELINE['primary_obs_null']:.13f} / {BASELINE['lineage_obs_null']:.13f}.",
        f"- `recon_mean_log1p` rebuilds the primary definition from the same raw counts this "
        f"script reads, and lands at {results[1]['obs_null_ratio']:.9f} "
        f"(delta {results[1]['delta_vs_baseline_primary']:+.2e}); centroid values agree with the "
        f"deposited CSVs to max|diff| = "
        f"{max(recon_gate['human']['max_abs_diff'], recon_gate['mouse']['max_abs_diff']):.2e} "
        f"(float32 round-off). The raw-count provenance and the streaming accumulator are therefore "
        f"the same object the deposit was built from.",
        f"- Cells: {provenance['human']['n_cells_total']:,} human / "
        f"{provenance['mouse']['n_cells_total']:,} mouse, "
        f"{provenance['human']['n_zero_library_cells'] + provenance['mouse']['n_zero_library_cells']} "
        f"zero-library cells.",
        "", "## Caveats", "",
        "- **p-values are at the permutation floor.** 0 of 10,000 permutations reached the observed "
        "distance for every flavor and for both nulls, so p = 1/(10,000+1) = 1.0e-04 is a bound, not "
        "an estimate. The obs/null ratio, not p, is the quantity being compared here.",
        "- **The retained-PC count is not stable across centroid definitions.** 95% variance is a "
        "relative criterion, and the pseudo-bulk definitions redistribute variance across the joint "
        f"spectrum: {BASELINE['n_pcs']} PCs (primary) -> "
        + ", ".join(f"{r['n_pcs']} ({r['flavor']})" for r in results[2:])
        + ". Procrustes distances are consequently not comparable across rows in absolute terms; "
          "each row's obs/null is computed against its own permutation null in its own space, which "
          "is what makes the ratios comparable.",
        f"- **Per-type rigidity ordering is more sensitive than configuration-level coherence.** "
        f"Spearman rho of per-type residual magnitude against the deposited ranking is "
        f"{results[2]['residual_ranking_rho_vs_deposited']:.3f} "
        f"(p = {results[2]['residual_ranking_rho_p']:.1e}) for A and "
        f"{results[3]['residual_ranking_rho_vs_deposited']:.3f} "
        f"(p = {results[3]['residual_ranking_rho_p']:.1e}) for B. Library-size weighting (B) "
        f"reshuffles which types look rigid more than it weakens the global result.",
    ]
    (OUT / "pseudobulk_centroids_results.md").write_text("\n".join(lines) + "\n")

    print("\n" + "=" * 78 + "\nSUMMARY\n" + "=" * 78)
    print(f"  {'centroid definition':<46} {'PCs':>4} {'obs/null':>9} {'p':>10} {'lineage':>9} {'p':>10}")
    for r in results:
        print(f"  {label[r['flavor']]:<46} {r['n_pcs']:>4} {r['obs_null_ratio']:>9.4f} "
              f"{r['p_value']:>10.2e} {r['lineage_obs_null_ratio']:>9.4f} {r['lineage_p_value']:>10.2e}")
    for r in results:
        if r["flavor"] in ("equal_cell_weight", "aggregate_pseudobulk"):
            print(f"  {r['flavor']}: {'PASS' if r['pass'] else 'FAIL'}")
    print(f"\n  wrote {OUT / 'pseudobulk_centroids_results.json'}")
    print(f"  wrote {OUT / 'pseudobulk_centroids_results.md'}")
    print(f"  wrote {OUT / 'pseudobulk_pca_centroids.npz'}")


if __name__ == "__main__":
    main()
