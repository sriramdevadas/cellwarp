#!/usr/bin/env python3
"""Layer-1 ribosomal/housekeeping exclusion sensitivity analysis.
Drops ribosomal (and ribosomal+housekeeping) genes from the 35-type centroids,
refits joint PCA (95%), recomputes Procrustes obs/null + per-type rigidity ranking
+ 10k-perm null. Reuses the canonical src/cellwarp/procrustes.py. Outputs to
analysis/sensitivity_analyses/. Does NOT touch manuscript/Table S1/figures/validate/CROSSWALK.
"""
import ast, json, re, sys
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.decomposition import PCA
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
from cellwarp.procrustes import procrustes_align, permutation_test

OUT = REPO / "analysis/sensitivity_analyses"
SEED, NPERM = 42, 10000
RIBO_RE = re.compile(r"^(RPL|RPS|Rpl|Rps)\d")

# ---- HOUNKPE housekeeping set (repo canonical), parsed from the script ----
hk_src = (REPO / "scripts/12_housekeeping_ratio.py").read_text()
m = re.search(r"HOUNKPE_HK_GENES\s*=\s*(\{.*?\})", hk_src, re.S)
HK_SYMBOLS = ast.literal_eval(m.group(1))
print(f"HOUNKPE housekeeping symbols: {len(HK_SYMBOLS)}")

# ---- data ----
hc = pd.read_csv(REPO/"output/phase2/scaled_35types/centroids_human_35.csv", index_col=0)
mc = pd.read_csv(REPO/"output/phase2/scaled_35types/centroids_mouse_35.csv", index_col=0)
saved = np.load(REPO/"output/phase2/scaled_35types/pca_centroids_35.npz", allow_pickle=True)
cell_types = [str(c) for c in saved["cell_types"]]
hc = hc.loc[cell_types]; mc = mc.loc[cell_types]
genes = [c for c in hc.columns if c.startswith("ENSG")]
print(f"cell_types={len(cell_types)}  genes={len(genes)}")

orth = pd.read_csv(REPO/"data/phase1/orthologs_human_mouse.csv")
ribo_mask = orth["human_gene_name"].astype(str).str.match(RIBO_RE) | orth["mouse_gene_name"].astype(str).str.match(RIBO_RE)
ribo_ens = set(orth.loc[ribo_mask, "human_ensembl_id"])
sym2ens = dict(zip(orth["human_gene_name"], orth["human_ensembl_id"]))
hk_ens = set(sym2ens[s] for s in HK_SYMBOLS if s in sym2ens)
print(f"ribosomal ENSG in space: {len(ribo_ens & set(genes))} | housekeeping ENSG in space: {len(hk_ens & set(genes))}")

full_rank = pd.read_csv(REPO/"output/phase2/scaled_35types/residuals_ranked.csv").set_index("cell_type")["residual_magnitude"]

def run_variant(name, excluded):
    kept = [g for g in genes if g not in excluded]
    H = hc[kept].values; M = mc[kept].values
    combined = np.vstack([H, M])
    pca = PCA(n_components=0.95, svd_solver="full", random_state=SEED)
    Z = pca.fit_transform(combined)
    n = len(cell_types)
    Hp, Mp = Z[:n], Z[n:]
    res = procrustes_align(Hp, Mp)
    _, null = permutation_test(Hp, Mp, n_permutations=NPERM, seed=SEED)
    obs_null = float(res.distance / np.median(null))
    resid = np.linalg.norm(res.centered_reference - res.aligned_target, axis=1)
    rank_s = pd.Series(resid, index=cell_types)
    common = [c for c in cell_types if c in full_rank.index]
    rho, p = spearmanr(rank_s.loc[common].values, full_rank.loc[common].values)
    ranking = rank_s.sort_values(ascending=False)  # rank1 = most displaced/flexible
    out = {
        "variant": name, "n_excluded": len(set(excluded)&set(genes)), "n_kept": len(kept),
        "pca_components": int(pca.n_components_), "pca_cumvar": float(pca.explained_variance_ratio_.sum()),
        "procrustes_distance": float(res.distance), "scaling_s": float(res.scaling),
        "null_median": float(np.median(null)), "obs_null_ratio": obs_null,
        "ranking_rho_vs_fullspace": float(rho), "ranking_rho_p": float(p),
        "top5_flexible": list(ranking.index[:5]), "top5_rigid": list(ranking.index[-5:][::-1]),
    }
    rank_df = pd.DataFrame({"cell_type": ranking.index, "residual_magnitude": ranking.values,
                            "rank_flexible_to_rigid": range(1, n+1)})
    rank_df.to_csv(OUT/f"layer1_exclusion_ranking_{name}.csv", index=False)
    return out

variants = [
    ("full_sanity", set()),
    ("ribosomal_only", ribo_ens),
    ("ribosomal_plus_housekeeping", ribo_ens | hk_ens),
]
results = [run_variant(n, e) for n, e in variants]
json.dump({"seed": SEED, "n_permutations": NPERM, "headline_obs_null": 0.522,
           "housekeeping_list": "HOUNKPE_HOUSEKEEPING_GENES (MSigDB)", "variants": results},
          open(OUT/"layer1_exclusion_results.json", "w"), indent=2)
print("\n===== M4 SUMMARY =====")
for r in results:
    print(f"  {r['variant']:30s} obs/null={r['obs_null_ratio']:.4f}  (PCA {r['pca_components']}c, {r['pca_cumvar']*100:.1f}%)  "
          f"n_excl={r['n_excluded']:5d}  ranking-rho vs full={r['ranking_rho_vs_fullspace']:.3f} (p={r['ranking_rho_p']:.1e})")
