"""
Census Pan-Census Replication — Step 2 (Local h5ad Files)
==========================================================
Processes locally-downloaded h5ad files from independent Census datasets.
Much faster than Census API for expression data access.
"""

import gc
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))
from cellwarp.procrustes import (
    compute_centroids,
    pca_reduce_centroids,
    procrustes_align,
    permutation_test,
    compute_residual_vectors,
)

OUT_DIR = Path(__file__).parent
DATA_DIR = Path(__file__).resolve().parent.parent.parent
H5AD_DIR = OUT_DIR / "h5ad_cache"

RANDOM_SEED = 42
MAX_CELLS = 2_000
MIN_CELLS = 200

PRIMARY_PATH = DATA_DIR / "output/phase2/scaled_35types/residuals_ranked.csv"
ORTHO_PATH = DATA_DIR / "data/phase1/orthologs_human_mouse.csv"

TARGET_TYPES = [
    "B cell", "CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell",
    "T cell", "basal cell", "endothelial cell", "epithelial cell", "fibroblast",
    "fibroblast of cardiac tissue", "granulocyte", "hepatocyte",
    "macrophage", "mature NK T cell", "mesenchymal stem cell", "monocyte",
    "myeloid leukocyte", "natural killer cell", "neutrophil",
    "pancreatic acinar cell", "pancreatic ductal cell", "plasma cell",
    "smooth muscle cell", "stromal cell",
    # HSC dropped: no independent mouse data outside 43GB dataset
]

MOUSE_FILES = [
    "7b6bab5a-f9c4-4a56-9ed4-3b9079b14867.h5ad",
    "58b01044-c5e5-4b0f-8a2d-6ebf951e01ff.h5ad",
    "4c4cfb38-c2af-4524-8ef8-bbcf1b6e2670.h5ad",
    "a2da8d7b-54a8-47d1-a0d3-aafcd0535f00.h5ad",
    "731e0ae7-e600-470f-a6dc-8c35c28d6c3d.h5ad",
    "047d57f2-4d14-45de-aa98-336c6f583750.h5ad",
    "6c6b4c47-096d-4084-97e7-714ee10c556c.h5ad",
    "49e4ffcc-5444-406d-bdee-577127404ba8.h5ad",
    "25818bf7-e2a7-41ec-8ff2-bc369c0ff4f5.h5ad",
]

HUMAN_FILES = [
    "2adb1f8a-a6b1-4909-8ee8-484814e2d4bf.h5ad",
    "b61a921b-7fa3-4b42-b455-aaaf32447920.h5ad",
    "fd072bc3-2dfb-46f8-b4e3-467cb3223182.h5ad",
    "37b21763-7f0f-41ae-9001-60bad6e2841d.h5ad",
    "ee195b7d-184d-4dfa-9b1c-51a7e601ac11.h5ad",
    "65badd7a-9262-4fd1-9ce2-eb5dc0ca8039.h5ad",
]


def load_and_filter(h5ad_path, target_types):
    """Load h5ad, filter to target cell types + healthy primary data."""
    adata = sc.read_h5ad(h5ad_path, backed="r")

    # Determine cell_type column
    ct_col = "cell_type" if "cell_type" in adata.obs.columns else None
    if ct_col is None:
        for col in adata.obs.columns:
            if "cell_type" in col.lower():
                ct_col = col
                break
    if ct_col is None:
        print(f"    No cell_type column found. Columns: {list(adata.obs.columns)[:10]}")
        return None

    # Filter
    mask = adata.obs[ct_col].isin(target_types)
    if "is_primary_data" in adata.obs.columns:
        mask = mask & (adata.obs["is_primary_data"] == True)
    if "disease" in adata.obs.columns:
        mask = mask & (adata.obs["disease"] == "normal")

    n_match = mask.sum()
    if n_match == 0:
        return None

    result = adata[mask].to_memory()
    # Ensure cell_type column name is standardized
    if ct_col != "cell_type":
        result.obs["cell_type"] = result.obs[ct_col]

    return result


def process_species(files, target_types, label):
    """Load all files for one species, merge, subsample."""
    rng = np.random.default_rng(RANDOM_SEED)
    type_cells = {}

    all_adatas = []
    for fname in files:
        fpath = H5AD_DIR / fname
        if not fpath.exists():
            print(f"  MISSING: {fname}")
            continue

        print(f"  {fname[:12]}...", end=" ", flush=True)
        adata = load_and_filter(fpath, target_types)
        if adata is None or adata.n_obs == 0:
            print("0 matching cells")
            continue

        # Subsample per type (respecting budget already used)
        indices = []
        for ct in target_types:
            ct_mask = adata.obs["cell_type"] == ct
            ct_idx = np.where(ct_mask)[0]
            n_avail = len(ct_idx)
            if n_avail == 0:
                continue

            already = type_cells.get(ct, 0)
            needed = min(n_avail, MAX_CELLS - already)
            if needed <= 0:
                continue

            if n_avail > needed:
                sel = rng.choice(ct_idx, size=needed, replace=False)
                sel.sort()
            else:
                sel = ct_idx

            indices.extend(sel.tolist())
            type_cells[ct] = already + len(sel)

        if indices:
            sub = adata[indices].copy()
            all_adatas.append(sub)
            cts = sub.obs["cell_type"].nunique()
            print(f"{adata.n_obs:,} total, {sub.n_obs:,} selected ({cts} types)")
        else:
            print("0 selected")

        del adata
        gc.collect()

    if not all_adatas:
        raise RuntimeError(f"No {label} data loaded")

    # Merge
    merged = ad.concat(all_adatas, join="inner")
    del all_adatas
    gc.collect()

    print(f"  {label} total: {merged.n_obs:,} x {merged.n_vars:,}")
    print(f"  Per-type counts:")
    for ct in sorted(type_cells.keys()):
        print(f"    {ct:<50} {type_cells[ct]:>5}")

    return merged, type_cells


def filter_orthologs(h, m, ortho):
    """Filter to shared 1:1 orthologs."""
    # Find gene ID column
    h_id_col = "feature_id" if "feature_id" in h.var.columns else h.var.index.name or "index"
    m_id_col = "feature_id" if "feature_id" in m.var.columns else m.var.index.name or "index"

    # Get gene IDs
    if "feature_id" in h.var.columns:
        h_ids = set(h.var["feature_id"])
    else:
        h_ids = set(h.var.index)

    if "feature_id" in m.var.columns:
        m_ids = set(m.var["feature_id"])
    else:
        m_ids = set(m.var.index)

    # Try both Ensembl IDs and gene names
    sh = ortho[ortho["human_ensembl_id"].isin(h_ids) & ortho["mouse_ensembl_id"].isin(m_ids)]

    if len(sh) < 1000:
        # Try gene names
        h_names = set(h.var["feature_name"]) if "feature_name" in h.var.columns else set()
        m_names = set(m.var["feature_name"]) if "feature_name" in m.var.columns else set()
        sh2 = ortho[ortho["human_gene_name"].isin(h_names) & ortho["mouse_gene_name"].isin(m_names)]
        if len(sh2) > len(sh):
            print(f"  Using gene names instead of Ensembl IDs ({len(sh2)} vs {len(sh)})")
            sh = sh2
            h_id_col, m_id_col = "human_gene_name", "mouse_gene_name"
            if "feature_name" in h.var.columns:
                h.var.index = h.var["feature_name"].values
            if "feature_name" in m.var.columns:
                m.var.index = m.var["feature_name"].values
        else:
            h_id_col, m_id_col = "human_ensembl_id", "mouse_ensembl_id"
    else:
        h_id_col, m_id_col = "human_ensembl_id", "mouse_ensembl_id"

    sh = sh.drop_duplicates(h_id_col).drop_duplicates(m_id_col)
    m2h = dict(zip(sh[m_id_col], sh[h_id_col]))

    # Filter and align
    h_gene_set = set(sh[h_id_col])
    m_gene_set = set(sh[m_id_col])

    if "feature_id" in h.var.columns and h_id_col == "human_ensembl_id":
        h_mask = h.var["feature_id"].isin(h_gene_set)
        hf = h[:, h_mask].copy()
        hf.var.index = hf.var["feature_id"].values
    elif "feature_name" in h.var.columns and h_id_col == "human_gene_name":
        h_mask = h.var["feature_name"].isin(h_gene_set)
        hf = h[:, h_mask].copy()
        hf.var.index = hf.var["feature_name"].values
    else:
        hf = h[:, [g for g in h.var.index if g in h_gene_set]].copy()

    if "feature_id" in m.var.columns and m_id_col == "mouse_ensembl_id":
        m_mask = m.var["feature_id"].isin(m_gene_set)
        mf = m[:, m_mask].copy()
        mf.var.index = [m2h[g] for g in mf.var["feature_id"]]
    elif "feature_name" in m.var.columns and m_id_col == "mouse_gene_name":
        m_mask = m.var["feature_name"].isin(m_gene_set)
        mf = m[:, m_mask].copy()
        mf.var.index = [m2h[g] for g in mf.var["feature_name"]]
    else:
        m_genes = [g for g in m.var.index if g in m_gene_set]
        mf = m[:, m_genes].copy()
        mf.var.index = [m2h[g] for g in m_genes]

    # Common genes
    common = sorted(set(hf.var.index) & set(mf.var.index))
    hf = hf[:, common].copy()
    mf = mf[:, common].copy()

    print(f"  Shared orthologs: {len(common):,}")
    return hf, mf, len(common)


def main():
    print("=" * 70)
    print("CENSUS PAN-CENSUS REPLICATION (local h5ad)")
    print("=" * 70)

    primary = pd.read_csv(PRIMARY_PATH)
    plookup = {r["cell_type"]: r["residual_magnitude"] for _, r in primary.iterrows()}
    ortho = pd.read_csv(ORTHO_PATH)
    print(f"Types: {len(TARGET_TYPES)}, Orthologs: {len(ortho):,}")

    # Process mouse
    print("\n[1] Mouse datasets...")
    mouse, m_counts = process_species(MOUSE_FILES, TARGET_TYPES, "Mouse")

    # Process human
    print("\n[2] Human datasets...")
    human, h_counts = process_species(HUMAN_FILES, TARGET_TYPES, "Human")

    # Shared types
    shared = sorted(ct for ct in TARGET_TYPES
                    if m_counts.get(ct, 0) >= MIN_CELLS and h_counts.get(ct, 0) >= MIN_CELLS)
    print(f"\nShared (>={MIN_CELLS}): {len(shared)} types")
    mouse = mouse[mouse.obs["cell_type"].isin(shared)].copy()
    human = human[human.obs["cell_type"].isin(shared)].copy()

    # Orthologs
    print("\n[3] Orthologs...")
    hf, mf, ngenes = filter_orthologs(human, mouse, ortho)
    del human, mouse; gc.collect()

    # Normalize
    print("\n[4] Normalize...")
    sc.pp.normalize_total(hf, target_sum=1e4); sc.pp.log1p(hf)
    sc.pp.normalize_total(mf, target_sum=1e4); sc.pp.log1p(mf)

    final_types = sorted(set(hf.obs["cell_type"].unique()) & set(mf.obs["cell_type"].unique()))
    for ct in final_types:
        h = (hf.obs["cell_type"] == ct).sum()
        m = (mf.obs["cell_type"] == ct).sum()
        print(f"  {ct:<50} H:{h:>5} M:{m:>5}")

    # Pipeline
    print("\n[5] Centroids...")
    hc = compute_centroids(hf); mc = compute_centroids(mf)

    print("\n[6] PCA...")
    hp, mp, pca, cell_types = pca_reduce_centroids(hc, mc)

    print("\n[7] Procrustes...")
    result = procrustes_align(hp, mp)

    print("\n[8] Permutation test (10K)...")
    pval, null = permutation_test(hp, mp)

    print("\n[9] Residuals...")
    resids = compute_residual_vectors(result, cell_types)

    obs_null = result.distance / np.median(null)
    rr = {ct: float(np.linalg.norm(resids[ct])) for ct in cell_types}
    ranked = sorted(rr.items(), key=lambda x: x[1], reverse=True)

    print("\n" + "=" * 70)
    print("RESULTS")
    print("=" * 70)
    for i, (ct, mag) in enumerate(ranked, 1):
        print(f"  {i:>2}. {ct:<50} {mag:.4f}")

    # Ranking correlation
    sp = [ct for ct in cell_types if ct in plookup]
    ns = len(sp)
    comp = pd.DataFrame([{"cell_type": ct, "primary_residual": plookup[ct],
                           "replication_residual": rr[ct]} for ct in sp])
    comp["primary_rank"] = comp["primary_residual"].rank(ascending=False).astype(int)
    comp["replication_rank"] = comp["replication_residual"].rank(ascending=False).astype(int)
    comp["rank_shift"] = abs(comp["primary_rank"] - comp["replication_rank"])
    comp = comp.sort_values("primary_rank")
    rho, rpval = spearmanr(comp["primary_residual"], comp["replication_residual"])

    print(f"\n  Spearman rho = {rho:.4f}, p = {rpval:.4f}, n = {ns}")
    for _, r in comp.iterrows():
        print(f"  {r['cell_type']:<50} P:{int(r['primary_rank']):>2} "
              f"R:{int(r['replication_rank']):>2} shift:{int(r['rank_shift'])}")

    # Save outputs
    print("\n[10] Saving...")
    pd.DataFrame([{"cell_type": ct, "in_primary": ct in plookup,
                    "human_cells": h_counts.get(ct, 0), "mouse_cells": m_counts.get(ct, 0)}
                   for ct in cell_types]).to_csv(OUT_DIR / "cell_type_matching.csv", index=False)

    pd.DataFrame([{"rank": i, "cell_type": ct, "residual_magnitude": mag,
                    "primary_residual": plookup.get(ct)}
                   for i, (ct, mag) in enumerate(ranked, 1)]
    ).to_csv(OUT_DIR / "replication_residuals.csv", index=False)
    comp.to_csv(OUT_DIR / "ranking_comparison.csv", index=False)

    rj = {"dataset": "Pan-Census (local h5ad)", "census_version": "2025-11-08",
          "n_cell_types": len(cell_types), "cell_types": cell_types,
          "n_shared_with_primary": ns, "n_shared_genes": ngenes,
          "human_cells": int(hf.n_obs), "mouse_cells": int(mf.n_obs),
          "procrustes": {"distance": float(result.distance),
                         "distance_squared": float(result.distance_squared),
                         "scaling": float(result.scaling)},
          "permutation_test": {"p_value": float(pval), "n_permutations": len(null),
                               "null_mean": float(np.mean(null)),
                               "null_median": float(np.median(null)),
                               "null_std": float(np.std(null)),
                               "obs_null_ratio": float(obs_null)},
          "ranking_correlation": {"spearman_rho": float(rho), "p_value": float(rpval),
                                  "n_types": ns},
          "per_type_residuals": {ct: float(rr[ct]) for ct in cell_types},
          "comparison": {"Sun2023": {"rho": 0.146, "p": 0.603, "n": 15},
                         "PanSci": {"rho": 0.194, "p": 0.471, "n": 16},
                         "CellHint": {"rho": -0.042, "p": 0.897, "n": 12},
                         "primary": {"obs_null": 0.522, "p": 0.0001, "n_types": 35}}}
    np.save(OUT_DIR / "null_distribution.npy", null)
    with open(OUT_DIR / "replication_results.json", "w") as f:
        json.dump(rj, f, indent=2)

    with open(OUT_DIR / "dataset_selection.md", "w") as f:
        f.write(f"# Census Replication - Dataset Selection\n\n"
                f"## Pan-Census Independent Pool\n"
                f"Downloaded h5ad files from {len(MOUSE_FILES)} mouse + {len(HUMAN_FILES)} human\n"
                f"independent CELLxGENE datasets (2025-11-08 build).\n\n"
                f"Excludes: Tabula Sapiens, Tabula Muris Senis, CellHint, PanSci.\n"
                f"HSC dropped: no independent mouse data outside 43GB embryonic dataset.\n\n"
                f"## Final: {len(cell_types)} types, {ngenes:,} genes, <={MAX_CELLS} cells/type\n")

    with open(OUT_DIR / "census_replication_summary.md", "w") as f:
        f.write(f"# Census Pan-Census Replication\n\n"
                f"## Global Coherence\n"
                f"- obs/null = {obs_null:.4f} (primary: 0.522)\n"
                f"- p = {pval:.6f}\n\n"
                f"## Ranking: rho = {rho:.3f}, p = {rpval:.3f}, n = {ns}\n\n"
                f"| Replication | rho | p | n |\n|--|--|--|--|\n"
                f"| Sun2023 | 0.146 | 0.603 | 15 |\n"
                f"| PanSci | 0.194 | 0.471 | 16 |\n"
                f"| CellHint | -0.042 | 0.897 | 12 |\n"
                f"| **Pan-Census** | **{rho:.3f}** | **{rpval:.3f}** | **{ns}** |\n\n"
                f"## Per-Type\n| Type | P | R | Shift |\n|--|--|--|--|\n")
        for _, r in comp.iterrows():
            f.write(f"| {r['cell_type']} | {int(r['primary_rank'])} | "
                    f"{int(r['replication_rank'])} | {int(r['rank_shift'])} |\n")

    print(f"\n{'='*70}")
    print(f"DONE: {len(cell_types)} types, rho={rho:.3f}, obs/null={obs_null:.3f}")
    print(f"Sun2023=0.146 | PanSci=0.194 | CellHint=-0.042 | Census={rho:.3f}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
