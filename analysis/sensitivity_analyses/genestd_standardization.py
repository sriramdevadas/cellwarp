#!/usr/bin/env python3
"""
Per-gene standardization sensitivity analysis (Layer-1 + CPC1).

Question: are (i) the Layer-1 cross-species geometric coherence, (ii) the per-type
divergence ranking, and (iii) the CPC1 ribosomal dominance (Table S6: 25/35) artifacts
of the primary centroid normalization (means of log1p CP10K, no per-gene scaling)?

Design: reuse the SHIPPED primary pipeline unchanged (deposited gene-space centroids;
sklearn joint PCA; cellwarp.procrustes alignment + label-permutation null; the
generate_table_S6.py CPC1 computation) and inject ONLY per-gene standardization at
two stages, holding PCA dimensionality at the primary's 33 components:

  base : raw log1p-CP10K centroids                         (faithfulness gate -> 0.522, 25/35)
  A    : per-gene z-score across the 70 CENTROIDS, then PCA-33
  B    : per-gene z-score across CELLS per species (pp.scale-style), then PCA-33

Outputs (tracked):
  analysis/sensitivity_analyses/genestd_results.json
  docs/supplementary_materials/table_S9_genestd_standardization.csv

Producer for Supplementary Table S9. Null: seed 42; faithfulness gate at 10,000
permutations (== primary), reported set at 100,000 permutations.
"""
from __future__ import annotations

import gc
import json
import re
import sys
import time
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
from sklearn.decomposition import PCA

warnings.filterwarnings("ignore", message=".*encountered in det.*")

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
from cellwarp.procrustes import (  # noqa: E402
    procrustes_align, permutation_test, compute_residual_vectors,
)

SCALED = REPO / "output" / "phase2" / "scaled_35types"
DATA = REPO / "data" / "phase2_scaled"
OUT = REPO / "analysis" / "sensitivity_analyses"
SUPP = REPO / "docs" / "supplementary_materials"
SEED = 42
N_PCA = 33
TOP_N = 20
N_PERM_GATE = 10_000     # exact-match gate (== primary's permutation count)
N_PERM = 100_000         # reported set (firmer ratio/p)

SHIPPED_OBS_NULL = 0.5222043226858066
SHIPPED_SPECIFIC = {  # the 10 cell-type-specific types in shipped Table S6
    "stromal cell", "hematopoietic stem cell", "neutrophil", "plasma cell",
    "classical monocyte", "large intestine goblet cell", "monocyte",
    "intermediate monocyte", "mature NK T cell", "hepatocyte",
}


def is_ribosomal(symbol: str) -> bool:
    return bool(re.match(r"^(RPL|RPS|Rpl|Rps)\d", str(symbol)))


def joint_pca(combined, n_components):
    p = PCA(n_components=n_components, svd_solver="full", random_state=SEED)
    Z = p.fit_transform(combined)
    return p, Z[:35], Z[35:]


def layer1(human_pca, mouse_pca, cell_types, n_perm, label):
    """Shipped Procrustes + label-permutation null. Returns obs/null + per-type residuals."""
    res = procrustes_align(human_pca, mouse_pca)
    _, null = permutation_test(human_pca, mouse_pca, n_permutations=n_perm, seed=SEED)
    null_median = float(np.median(null))
    p_lt = float((np.sum(null <= res.distance) + 1) / (n_perm + 1))
    resid = compute_residual_vectors(res, cell_types)
    mags = {ct: float(np.linalg.norm(resid[ct])) for ct in cell_types}
    print(f"    {label:16s} obs/null={res.distance/null_median:.4f}  p={p_lt:.2e}")
    return dict(distance=float(res.distance), null_median=null_median,
                obs_null=res.distance / null_median, p_value=p_lt, resid_mag=mags)


def chunked_gene_mean_std(path, gene_ids, chunk=4000):
    a = ad.read_h5ad(path, backed="r")
    assert list(a.var_names) == list(gene_ids), "gene order mismatch"
    n, G = a.shape
    s = np.zeros(G); s2 = np.zeros(G)
    for i in range(0, n, chunk):
        Xc = a.X[i:min(i + chunk, n)]
        Xc = (Xc.toarray() if sp.issparse(Xc) else np.asarray(Xc)).astype(np.float64)
        s += Xc.sum(0); s2 += (Xc * Xc).sum(0)
    a.file.close()
    mean = s / n
    return mean, np.sqrt(np.clip(s2 / n - mean ** 2, 0.0, None)), n


def covs_for_species(path, gene_ids, transformers, pcas):
    """Per-type within-type covariance (33x33) in each variant's PCA space."""
    print(f"    loading {path.name} (full)...", flush=True)
    t0 = time.time()
    a = ad.read_h5ad(path)[:, gene_ids]
    X_all = a.X.tocsr() if sp.issparse(a.X) else a.X
    ct_col = a.obs["cell_type"].values
    cts = sorted(np.unique(ct_col))
    print(f"      {X_all.shape[0]} cells, {len(cts)} types ({time.time()-t0:.0f}s)", flush=True)
    out = {V: {} for V in pcas}; ncells = {}
    for ct in cts:
        idx = np.where(ct_col == ct)[0]
        Xslice = X_all[idx]
        Xt = (Xslice.toarray() if sp.issparse(Xslice) else np.asarray(Xslice)).astype(np.float64)
        ncells[ct] = Xt.shape[0]
        for V in pcas:
            Xp = pcas[V].transform(transformers[V](Xt))
            cen = Xp - Xp.mean(0)
            out[V][ct] = (cen.T @ cen) / (cen.shape[0] - 1)
    del X_all, a; gc.collect()
    return out, ncells


def cpc1_classify(cov_h, cov_m, n_h, n_m, W, gene_ids, ens2sym):
    """Replicate generate_table_S6.py: weighted cov -> top eigenvector -> gene loadings."""
    weighted = n_h * cov_h + n_m * cov_m
    eigvals, eigvecs = np.linalg.eigh(weighted)
    eigvecs = eigvecs[:, np.argsort(eigvals)[::-1]]
    cpc1 = eigvecs[:, 0]
    cpc1_genes = cpc1 @ W
    top_idx = np.argsort(np.abs(cpc1_genes))[::-1][:TOP_N]
    syms = [ens2sym.get(gene_ids[i], gene_ids[i]) for i in top_idx]
    rank1 = syms[0]
    return dict(classification="ribosomal-dominated" if is_ribosomal(rank1) else "cell-type-specific",
                rank1=rank1, n_ribo_top20=int(sum(is_ribosomal(s) for s in syms)),
                top5=syms[:5])


def main():
    t_start = time.time()
    print("=" * 78 + "\nM1 — per-gene standardization sensitivity (Layer-1 + CPC1)\n" + "=" * 78)

    saved = np.load(SCALED / "pca_centroids_35.npz", allow_pickle=True)
    cell_types = list(saved["cell_types"])
    saved_h, saved_m = saved["human"], saved["mouse"]
    hc = pd.read_csv(SCALED / "centroids_human_35.csv", index_col=0).loc[cell_types]
    mc = pd.read_csv(SCALED / "centroids_mouse_35.csv", index_col=0).loc[cell_types]
    gene_ids = list(hc.columns)
    H, M = hc.values.astype(np.float64), mc.values.astype(np.float64)
    combined = np.vstack([H, M])
    prim = pd.read_csv(SCALED / "residuals_ranked.csv")
    prim_mag = {r.cell_type: r.residual_magnitude for r in prim.itertuples()}
    a_meta = ad.read_h5ad(DATA / "human_scaled.h5ad", backed="r")
    ens2sym = dict(zip(a_meta.var_names, a_meta.var["feature_name"])); a_meta.file.close()

    # ---- cell-level per-gene mean/std (Variant B) ----
    print("\n[1] cell-level per-gene mean/std (chunked)...", flush=True)
    muH, sdH, nH = chunked_gene_mean_std(DATA / "human_scaled.h5ad", gene_ids)
    muM, sdM, nM = chunked_gene_mean_std(DATA / "mouse_scaled.h5ad", gene_ids)
    print(f"    human {nH} cells, mouse {nM} cells")

    # ---- feature spaces + PCA-33 ----
    mu_g, sd_g = combined.mean(0), combined.std(0, ddof=0)
    sdg_s = np.where(sd_g == 0, 1.0, sd_g)
    sdH_s = np.where(sdH == 0, 1.0, sdH); sdM_s = np.where(sdM == 0, 1.0, sdM)
    combinedA = (combined - mu_g) / sdg_s
    combinedB = np.vstack([(H - muH) / sdH_s, (M - muM) / sdM_s])

    print("\n[2] Layer-1 obs/null + ranking (joint PCA-33 per variant)...")
    cts = list(cell_types)
    layer1_out = {}

    # base: primary PCA(0.95)->33 exactly (gate) + 100k
    pca_b, h_b, m_b = joint_pca(combined, 0.95)
    pca_err = max(float(np.max(np.abs(h_b - saved_h))), float(np.max(np.abs(m_b - saved_m))))
    base_gate = layer1(h_b, m_b, cts, N_PERM_GATE, "base@10k")
    base_100k = layer1(h_b, m_b, cts, N_PERM, "base@100k")
    base_rho = float(stats.spearmanr([base_gate["resid_mag"][c] for c in cts],
                                     [prim_mag[c] for c in cts]).statistic)
    layer1_out["base"] = dict(pca_max_err_vs_deposited=pca_err, gate_10k=base_gate,
                              run_100k=base_100k, ranking_rho_vs_primary=base_rho)

    results = {}
    for name, cmb in [("A", combinedA), ("B", combinedB)]:
        pca_v, h_v, m_v = joint_pca(cmb, N_PCA)
        r = layer1(h_v, m_v, cts, N_PERM, f"variant {name}")
        rho = float(stats.spearmanr([r["resid_mag"][c] for c in cts],
                                    [prim_mag[c] for c in cts]).statistic)
        layer1_out[name] = dict(cumvar=float(pca_v.explained_variance_ratio_.sum()),
                                obs_null=r["obs_null"], p_value=r["p_value"],
                                ranking_rho_vs_primary=rho, resid_mag=r["resid_mag"])

    # ---- CPC1 covariance pass ----
    print("\n[3] CPC1 within-type covariances + classification...")
    pcas = {"base": pca_b,
            "A": PCA(N_PCA, svd_solver="full", random_state=SEED).fit(combinedA),
            "B": PCA(N_PCA, svd_solver="full", random_state=SEED).fit(combinedB)}
    Ws = {V: pcas[V].components_ for V in pcas}
    th = {"base": lambda X: X, "A": lambda X: (X - mu_g) / sdg_s, "B": lambda X: (X - muH) / sdH_s}
    tm = {"base": lambda X: X, "A": lambda X: (X - mu_g) / sdg_s, "B": lambda X: (X - muM) / sdM_s}
    cov_h, nmap_h = covs_for_species(DATA / "human_scaled.h5ad", gene_ids, th, pcas)
    cov_m, nmap_m = covs_for_species(DATA / "mouse_scaled.h5ad", gene_ids, tm, pcas)

    cpc1_out = {}
    for V in ["base", "A", "B"]:
        rows = {ct: cpc1_classify(cov_h[V][ct], cov_m[V][ct], nmap_h[ct], nmap_m[ct],
                                  Ws[V], gene_ids, ens2sym) for ct in cell_types}
        ribo = [ct for ct in cell_types if rows[ct]["classification"] == "ribosomal-dominated"]
        cpc1_out[V] = dict(n_ribosomal_dominated=len(ribo),
                           specific_types=sorted(ct for ct in cell_types if ct not in ribo),
                           mean_ribo_in_top20=float(np.mean([rows[ct]["n_ribo_top20"] for ct in cell_types])),
                           rank1=rows)
        print(f"    [{V}] ribosomal-dominated = {len(ribo)}/35 ; mean #ribo in top-20 = {cpc1_out[V]['mean_ribo_in_top20']:.2f}")

    # ---- faithfulness gates ----
    gate_l1 = abs(base_gate["obs_null"] - SHIPPED_OBS_NULL) < 1e-9 and abs(base_rho - 1.0) < 1e-9
    gate_cpc1 = (cpc1_out["base"]["n_ribosomal_dominated"] == 25 and
                 set(cpc1_out["base"]["specific_types"]) == SHIPPED_SPECIFIC)
    print(f"\n[4] FAITHFULNESS GATES: Layer-1 base==shipped: {gate_l1} ; CPC1 base==25/35+set: {gate_cpc1}")
    if not (gate_l1 and gate_cpc1):
        print("  !! GATE FAILED — base does not reproduce shipped; abort flag in results.json")

    # ---- variant-B CPC1 rank-1 marker loadings (cell-type-specific examples) ----
    B_markers = {ct: cpc1_out["B"]["rank1"][ct]["rank1"] for ct in cell_types}

    # ---- write results.json ----
    report = dict(
        meta=dict(seed=SEED, n_pca=N_PCA, n_perm_gate=N_PERM_GATE, n_perm_reported=N_PERM,
                  runtime_sec=time.time() - t_start),
        gates=dict(layer1_base=bool(gate_l1), cpc1_base=bool(gate_cpc1)),
        layer1={k: ({kk: vv for kk, vv in v.items() if kk not in ("resid_mag", "gate_10k", "run_100k")}
                    | ({"obs_null_10k": v["gate_10k"]["obs_null"], "obs_null_100k": v["run_100k"]["obs_null"]}
                       if k == "base" else {}))
                for k, v in layer1_out.items()},
        cpc1={V: {kk: vv for kk, vv in cpc1_out[V].items() if kk != "rank1"} for V in ["base", "A", "B"]},
        variantB_cpc1_rank1=B_markers,
        per_type_resid={k: layer1_out.get(k, {}).get("resid_mag") for k in ("A", "B")},
    )
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "genestd_results.json").write_text(json.dumps(report, indent=2, default=str))

    # ---- Supplementary Table S9 ----
    SUPP.mkdir(parents=True, exist_ok=True)
    markers_of_interest = ["B cell", "plasma cell", "hepatocyte", "endothelial cell",
                           "natural killer cell", "basal cell", "smooth muscle cell"]
    s9 = pd.DataFrame([
        dict(scheme="Primary (no per-gene scaling)", layer1_obs_null=round(base_100k["obs_null"], 3),
             layer1_p=f"<1e-5", ranking_rho_vs_primary=round(base_rho, 3),
             cpc1_ribosomal_dominated_of_35=cpc1_out["base"]["n_ribosomal_dominated"],
             cpc1_mean_ribo_in_top20=round(cpc1_out["base"]["mean_ribo_in_top20"], 2)),
        dict(scheme="A: per-gene z across 70 centroids", layer1_obs_null=round(layer1_out["A"]["obs_null"], 3),
             layer1_p="<1e-5", ranking_rho_vs_primary=round(layer1_out["A"]["ranking_rho_vs_primary"], 3),
             cpc1_ribosomal_dominated_of_35=cpc1_out["A"]["n_ribosomal_dominated"],
             cpc1_mean_ribo_in_top20=round(cpc1_out["A"]["mean_ribo_in_top20"], 2)),
        dict(scheme="B: per-gene z across cells per species (pp.scale)", layer1_obs_null=round(layer1_out["B"]["obs_null"], 3),
             layer1_p="<1e-5", ranking_rho_vs_primary=round(layer1_out["B"]["ranking_rho_vs_primary"], 3),
             cpc1_ribosomal_dominated_of_35=cpc1_out["B"]["n_ribosomal_dominated"],
             cpc1_mean_ribo_in_top20=round(cpc1_out["B"]["mean_ribo_in_top20"], 2)),
    ])
    s9.to_csv(SUPP / "table_S9_genestd_standardization.csv", index=False)
    # Variant-B rank-1 CPC1 markers for a representative panel (for the Table S6 recharacterization)
    s9b = pd.DataFrame([dict(cell_type=ct, variantB_CPC1_rank1_driver=B_markers[ct]) for ct in markers_of_interest])
    s9b.to_csv(SUPP / "table_S9_schemeB_CPC1_markers.csv", index=False)

    print("\n" + "=" * 78)
    print(f"  base   obs/null={base_100k['obs_null']:.4f} (gate10k {base_gate['obs_null']:.4f})  rho={base_rho:.3f}  CPC1 ribo=25/35")
    print(f"  A      obs/null={layer1_out['A']['obs_null']:.4f}  rho={layer1_out['A']['ranking_rho_vs_primary']:.3f}  CPC1 ribo={cpc1_out['A']['n_ribosomal_dominated']}/35")
    print(f"  B      obs/null={layer1_out['B']['obs_null']:.4f}  rho={layer1_out['B']['ranking_rho_vs_primary']:.3f}  CPC1 ribo={cpc1_out['B']['n_ribosomal_dominated']}/35")
    print(f"  Variant-B rank-1 CPC1 markers (panel): " + ", ".join(f"{ct.split()[0]}:{B_markers[ct]}" for ct in markers_of_interest))
    print(f"  wrote {OUT/'genestd_results.json'}")
    print(f"  wrote {SUPP/'table_S9_genestd_standardization.csv'}")
    print(f"  runtime {time.time()-t_start:.0f}s")


if __name__ == "__main__":
    main()
