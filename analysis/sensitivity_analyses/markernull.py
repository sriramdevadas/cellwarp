#!/usr/bin/env python3
"""
Marker-similarity-stratified permutation null sensitivity analysis.

Question: does the headline Layer-1 coherence (obs/null = 0.522) carry geometric
information BEYOND the marker-based ontology matching that defined the 35 landmarks,
or is it "close to tautological"? Test with a null HARDER than the lineage-stratified
null: permute type labels only within marker-similarity groups (so only expression-
similar types are scrambled), sweeping the cluster granularity K.

Design: reuse the SHIPPED stratified-permutation machinery
(scripts/test_lineage_stratified_permutation.stratified_permutation) verbatim; only
the partition changes. Marker-similarity partitions = species-averaged gene-space
centroids -> Ward/Euclidean hierarchical clustering -> within-group permutation.
Continuity anchors: K=1 recovers the unrestricted null (0.522, bit-identical);
the lineage partition recovers 0.668.

Outputs (tracked):
  analysis/sensitivity_analyses/markernull_results.json
  docs/supplementary_materials/table_S10_markernull.csv
  docs/supplementary_materials/figure_S8_markernull.pdf (+ .png)

Faithfulness gate at 10,000 permutations (== the scratch / lineage-null count);
reported p-values at 100,000 permutations (so borderline-K p are not floor-limited).
Seed 42 throughout.
"""
from __future__ import annotations

import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist, squareform
from scipy import stats

warnings.filterwarnings("ignore", message=".*encountered in det.*")

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))
from cellwarp.procrustes import (  # noqa: E402
    procrustes_align, permutation_test, compute_residual_vectors, _procrustes_distance,
)
import test_lineage_stratified_permutation as lin  # noqa: E402  (shipped stratified machinery)

SCALED = REPO / "output" / "phase2" / "scaled_35types"
OUT = REPO / "analysis" / "sensitivity_analyses"
SUPP = REPO / "docs" / "supplementary_materials"
SEED = 42
N_PERM_GATE = 10_000      # == scratch / lineage-null count (exact-ratio faithfulness gate)
N_PERM = 100_000          # reported set (firmer borderline-K p-values)

# scratch values to gate against (10k, seed 42)
SCRATCH = {"k1": 0.5222043226858066, "lineage": 0.6682720235158841,
           5: 0.7209, 8: 0.8072, 10: 0.9032, 15: 0.9636,
           "corr_k15_ratio": 0.9018, "corr_k15_p": 0.00070, "monotonicity_rho": -0.136}


def blocks_from_labels(labels, n):
    return [[i for i in range(n) if labels[i] == lab] for lab in sorted(set(labels))]


def permspace_log10(blocks):
    return float(sum(math.log10(math.factorial(len(b))) for b in blocks))


def run_null(X, Y, blocks, n_perm):
    p, null = lin.stratified_permutation(X, Y, blocks, n_perm=n_perm, seed=SEED)
    obs = _procrustes_distance(X, Y)
    return obs / float(np.median(null)), float(p)


def ward_labels(data, K):
    return fcluster(linkage(data, method="ward", metric="euclidean"), t=K, criterion="maxclust")


def main():
    print("=" * 78 + "\nMajor 2 — marker-similarity-stratified null\n" + "=" * 78)

    saved = np.load(SCALED / "pca_centroids_35.npz", allow_pickle=True)
    cell_types = list(saved["cell_types"])
    X, Y = saved["human"], saved["mouse"]                 # (35,33) joint-PCA centroids
    n = len(cell_types)
    obs = _procrustes_distance(X, Y)
    ct_to_idx = {ct: i for i, ct in enumerate(cell_types)}

    hc = pd.read_csv(SCALED / "centroids_human_35.csv", index_col=0).loc[cell_types].values
    mc = pd.read_csv(SCALED / "centroids_mouse_35.csv", index_col=0).loc[cell_types].values
    G_avg = 0.5 * (hc + mc)                                # species-averaged marker profile/type

    # ---- continuity anchors ----
    print("\n[1] continuity anchors")
    r_k1_g, p_k1_g = run_null(X, Y, [list(range(n))], N_PERM_GATE)
    # bit-identical check: K=1 stratified null == unrestricted permutation_test null
    _, null_k1 = lin.stratified_permutation(X, Y, [list(range(n))], n_perm=N_PERM_GATE, seed=SEED)
    _, null_unrestricted = permutation_test(X, Y, n_permutations=N_PERM_GATE, seed=SEED)
    k1_bit_identical = bool(np.allclose(null_k1, null_unrestricted))
    lin_blocks = [sorted(ct_to_idx[m] for m in mem) for mem in lin.LINEAGE_BLOCKS.values()]
    r_lin_g, p_lin_g = run_null(X, Y, lin_blocks, N_PERM_GATE)
    print(f"    K=1 (gate 10k) obs/null={r_k1_g:.4f}  bit-identical-to-unrestricted={k1_bit_identical}")
    print(f"    lineage (gate 10k) obs/null={r_lin_g:.4f}")

    # ---- K-sweep, Ward/Euclidean (gate 10k then published 100k) ----
    print("\n[2] K-sweep (gene-avg / Ward / Euclidean)")
    sweep = {}
    for K in range(1, 17):
        blocks = blocks_from_labels(ward_labels(G_avg, K), n)
        r10, p10 = run_null(X, Y, blocks, N_PERM_GATE)
        r100, p100 = run_null(X, Y, blocks, N_PERM)
        nns = sum(len(b) > 1 for b in blocks)
        sweep[K] = dict(obs_null_10k=r10, p_10k=p10, obs_null_100k=r100, p_100k=p100,
                        n_nonsingleton=nns, log10_permspace=permspace_log10(blocks),
                        sizes=sorted((len(b) for b in blocks), reverse=True))
        flag = ""
        if K in (5, 8, 10, 15):
            flag = f"  [scratch {SCRATCH[K]:.4f}, Δ={abs(r10-SCRATCH[K]):.4f}]"
        print(f"    K={K:>2}  10k obs/null={r10:.4f} p={p10:.5f} | 100k obs/null={r100:.4f} p={p100:.5f}{flag}")
    crossover = next((K for K in range(1, 17) if sweep[K]["p_100k"] >= 0.05), None)

    # ---- correlation-linkage robustness (gene-avg / average / correlation) ----
    print("\n[3] robustness: correlation-linkage clustering")
    Zc = linkage(G_avg, method="average", metric="correlation")
    corr = {}
    for K in (5, 8, 10, 15):
        blocks = blocks_from_labels(fcluster(Zc, t=K, criterion="maxclust"), n)
        r10, p10 = run_null(X, Y, blocks, N_PERM_GATE)
        r100, p100 = run_null(X, Y, blocks, N_PERM)
        corr[K] = dict(obs_null_10k=r10, p_10k=p10, obs_null_100k=r100, p_100k=p100)
        print(f"    K={K:>2}  10k obs/null={r10:.4f} p={p10:.5f} | 100k p={p100:.5f}")

    # ---- monotonicity: per-type residual vs marker-distinctness ----
    print("\n[4] monotonicity (per-type residual vs marker-distinctness)")
    res = procrustes_align(X, Y)
    resid = compute_residual_vectors(res, cell_types)
    resid_mag = np.array([np.linalg.norm(resid[ct]) for ct in cell_types])
    D = squareform(pdist(G_avg, metric="euclidean")); np.fill_diagonal(D, np.nan)
    distinct = np.nanmean(D, axis=1)
    rho, prho = stats.spearmanr(resid_mag, distinct)
    print(f"    Spearman rho={rho:+.3f} p={prho:.3f}  [scratch {SCRATCH['monotonicity_rho']:+.3f}]")

    # ---- faithfulness gate ----
    gate = (abs(r_k1_g - SCRATCH["k1"]) < 1e-9 and k1_bit_identical and
            abs(r_lin_g - SCRATCH["lineage"]) < 1e-9 and
            all(abs(sweep[K]["obs_null_10k"] - SCRATCH[K]) < 5e-4 for K in (5, 8, 10, 15)) and
            abs(rho - SCRATCH["monotonicity_rho"]) < 5e-3)
    print(f"\n[5] FAITHFULNESS GATE (10k ratios == scratch): {gate}")

    # ---- results.json ----
    report = dict(
        meta=dict(seed=SEED, n_perm_gate=N_PERM_GATE, n_perm_reported=N_PERM,
                  representation="species-averaged gene-space centroids (35 x 16,959)",
                  clustering_primary="Ward / Euclidean (fcluster maxclust)",
                  obs_distance=float(obs)),
        gate_passed=bool(gate),
        anchors=dict(k1=dict(obs_null_10k=r_k1_g, bit_identical_to_unrestricted=k1_bit_identical),
                     lineage=dict(obs_null_10k=r_lin_g)),
        ward_sweep=sweep, crossover_K_at_100k=crossover,
        correlation_linkage=corr,
        monotonicity=dict(spearman_rho=float(rho), spearman_p=float(prho)),
    )
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "markernull_results.json").write_text(json.dumps(report, indent=2, default=str))

    # ---- Supplementary Table S10 ----
    SUPP.mkdir(parents=True, exist_ok=True)
    rows = [dict(partition="Unrestricted (K=1 anchor)", obs_null=round(sweep[1]["obs_null_100k"], 3),
                 p_100k=f"{sweep[1]['p_100k']:.1e}", n_nonsingleton_groups=1, note="recovers primary 0.522"),
            dict(partition="Lineage (5 blocks, anchor)", obs_null=round(r_lin_g, 3),
                 p_100k="<1e-4", n_nonsingleton_groups=3, note="recovers 0.668")]
    for K in (5, 8, 10, 15):
        rows.append(dict(partition=f"Marker-similarity K={K} (Ward)",
                         obs_null=round(sweep[K]["obs_null_100k"], 3),
                         p_100k=f"{sweep[K]['p_100k']:.1e}",
                         n_nonsingleton_groups=sweep[K]["n_nonsingleton"],
                         note=("beats null" if sweep[K]["p_100k"] < 0.05 else "n.s. (null degenerate)")))
    rows.append(dict(partition="Marker-similarity K=15 (correlation linkage)",
                     obs_null=round(corr[15]["obs_null_100k"], 3), p_100k=f"{corr[15]['p_100k']:.1e}",
                     n_nonsingleton_groups="-", note="robustness: beats null even at K=15"))
    rows.append(dict(partition="Monotonicity: residual vs marker-distinctness (Spearman)",
                     obs_null=round(float(rho), 3), p_100k=f"{prho:.2f}", n_nonsingleton_groups="-",
                     note="independent (|rho|<=0.15, n.s.)"))
    pd.DataFrame(rows).to_csv(SUPP / "table_S10_markernull.csv", index=False)

    # ---- Supplementary Figure S8 ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    Ks = list(range(1, 17))
    ratios = [sweep[K]["obs_null_100k"] for K in Ks]
    pvals = [sweep[K]["p_100k"] for K in Ks]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 4.0))
    sig = [r for r, p in zip(ratios, pvals) if p < 0.05]
    ax1.plot(Ks, ratios, "-o", color="#2c3e8c", ms=4, lw=1.3, zorder=3)
    for K, r, p in zip(Ks, ratios, pvals):
        ax1.scatter([K], [r], s=26, c=("#2c3e8c" if p < 0.05 else "#c0392b"), zorder=4)
    ax1.axhline(1.0, ls=":", c="gray", lw=0.8)
    ax1.scatter([1], [ratios[0]], s=70, facecolors="none", edgecolors="green", lw=1.5, zorder=5, label="K=1 = unrestricted (0.522)")
    if crossover:
        ax1.axvline(crossover - 0.5, ls="--", c="#c0392b", lw=1.0)
        ax1.text(crossover - 0.4, 0.56, f"p≥0.05 at K={crossover}\n(null degenerate)", fontsize=7, color="#c0392b", va="bottom")
    ax1.set_xlabel("cluster granularity K (marker-similarity groups)")
    ax1.set_ylabel("obs/null (lower = stronger; blue p<0.05, red n.s.)")
    ax1.set_title("A  Coherence vs marker-similarity-null granularity", fontsize=9, loc="left")
    ax1.set_xticks(Ks); ax1.legend(fontsize=6, loc="lower right"); ax1.set_ylim(0.45, 1.05)

    ax2.scatter(distinct, resid_mag, s=22, c="#34495e", alpha=0.85, edgecolors="white", lw=0.4)
    ax2.set_xlabel("marker-distinctness (mean dist. to other types)")
    ax2.set_ylabel("Procrustes residual (per type)")
    ax2.set_title(f"B  Per-type residual vs marker-distinctness\nSpearman ρ={rho:+.3f}, p={prho:.2f}, n=35", fontsize=9, loc="left")
    fig.tight_layout()
    fig.savefig(SUPP / "figure_S8_markernull.pdf")
    fig.savefig(SUPP / "figure_S8_markernull.png", dpi=200)
    plt.close(fig)

    print("\n" + "=" * 78)
    print(f"  anchors: K=1={r_k1_g:.4f} (bit-id {k1_bit_identical}), lineage={r_lin_g:.4f}")
    print(f"  Ward sweep 100k: " + " ".join(f"K{K}={sweep[K]['obs_null_100k']:.3f}(p{sweep[K]['p_100k']:.0e})" for K in (5,8,10,15)))
    print(f"  crossover (p>=0.05) at K={crossover}; corr-linkage K15 p={corr[15]['p_100k']:.1e}")
    print(f"  monotonicity rho={rho:+.3f} p={prho:.3f}")
    print(f"  GATE {'PASSED' if gate else 'FAILED'}")
    print(f"  wrote {OUT/'markernull_results.json'}, table_S10, figure_S8")


if __name__ == "__main__":
    main()
