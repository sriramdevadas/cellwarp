#!/usr/bin/env python3
"""R15 Part A.2 — Layer-2 covariance ellipsoid alignment on the PanSci atlas.

Sibling to scripts/t3b_ellipsoid_alignment.py (which runs the primary
TS-vs-Tabula-Muris analysis). This script reuses the shape-agnostic
helpers from the primary and adds three PanSci-specific functions:
    - reconstruct_pca_pansci()  : 16 TS + 16 PanSci centroids -> joint PCA
    - project_cells_pansci()    : TS h5ad + PanSci MTX -> 17-D per-cell coords
    - main_pansci()             : orchestrates C2-C9 + persists outputs

The Layer-1 (centroid Procrustes) result must reproduce the existing PanSci
replication at output/validation/pansci_replication/pansci_replication.json
to machine precision; otherwise the script halts before Layer-2 compute.

Conventions inherited from the primary:
    - PCA basis: trained jointly on 32 centroids (16 TS + 16 PanSci);
      individual cells projected as (X - pca.mean_) @ pca.components_.T
    - Endothelial cell uses PanSci lung tissue only (DECISION-101 /
      pansci_replication.py:533-537), both for centroid and per-cell input.
    - perm_test_label_shuffle is called once with the centroid-optimal R;
      it returns both pre-rotation (V_M raw) and post-rotation (R @ V_M)
      nulls sharing identical permutation indices (matches primary at
      t3b_ellipsoid_alignment.py:316-326). Seed 42, 10,000 perms.

No edits to t3b_ellipsoid_alignment.py or pansci_replication.py.
Outputs go to output/twolayer_pansci_replication/.
"""

from __future__ import annotations

import gc
import gzip
import json
import sys
import time
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
from scipy import stats

# Project root on sys.path so we can import the primary helpers and cellwarp.
PROJECT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT / "src"))
sys.path.insert(0, str(PROJECT / "scripts"))

from cellwarp.procrustes import (  # noqa: E402
    RANDOM_SEED,
    pca_reduce_centroids,
    procrustes_align,
    permutation_test,
)
from t3b_ellipsoid_alignment import (  # noqa: E402
    compute_covariance_eigen,
    subspace_overlap,
    compute_alignment_scores,
    perm_test_label_shuffle,
    common_principal_components,
    eigenvalue_conservation,
    correlate_with_rigidity,
)

# Reuse PanSci pre-processing helpers verbatim (read-only import).
from pansci_replication import (  # noqa: E402
    map_pansci_to_ontology,
    strip_organ_suffix,
    ENDO_TISSUE,
)


# ---------------------------------------------------------------------------
# Paths and constants
# ---------------------------------------------------------------------------
DATA_PHASE1 = PROJECT / "data" / "phase1"
DATA_SCALED = PROJECT / "data" / "phase2_scaled"
ORTHOLOG_PATH = DATA_PHASE1 / "orthologs_human_mouse.csv"
TS_H5AD_PATH = DATA_SCALED / "human_scaled.h5ad"
PANSCI_DIR = PROJECT / "data" / "replication" / "pansci"
PANSCI_CENTROIDS_PATH = PROJECT / "data" / "centroids" / "pansci_16type_centroids.csv"
TS_CENTROIDS_PATH = PROJECT / "output" / "phase2" / "scaled_35types" / "centroids_human_35.csv"
EXISTING_REPL_JSON = PROJECT / "output" / "validation" / "pansci_replication" / "pansci_replication.json"

OUTPUT_DIR = PROJECT / "output" / "twolayer_pansci_replication"
SCRATCH_DIR = OUTPUT_DIR / "scratch"

SEED = 42
N_PERM = 10_000
K_VALUES = [1, 3, 5]
VARIANCE_THRESHOLD = 0.95
TISSUES_TO_LOAD = ["lung", "liver", "colon"]

# Sixteen matched cell types (locked from the existing PanSci replication).
MATCHED_TYPES_16 = [
    "B cell",
    "CD4-positive, alpha-beta T cell",
    "CD8-positive, alpha-beta T cell",
    "T cell",
    "endothelial cell",
    "epithelial cell",
    "fibroblast",
    "granulocyte",
    "hepatocyte",
    "large intestine goblet cell",
    "macrophage",
    "monocyte",
    "myeloid dendritic cell",
    "myeloid leukocyte",
    "plasma cell",
    "smooth muscle cell",
]


# ---------------------------------------------------------------------------
# C2: Reconstruct the joint PCA model and verify Layer-1 reproduction
# ---------------------------------------------------------------------------
def reconstruct_pca_pansci(variance_threshold: float = VARIANCE_THRESHOLD):
    """Refit the joint PCA on 32 centroids (16 TS + 16 PanSci).

    Returns
    -------
    pca : sklearn PCA
        Joint PCA model trained on the 32 centroid matrix.
    human_pca, mouse_pca : (16, k) arrays
        Centroids in the joint PCA space (TS-as-human, PanSci-as-mouse).
    types_list : list[str]
        Cell type ordering used by pca_reduce_centroids (sorted ascending).
    gene_names : list[str]
        16,959 ortholog gene names (human Ensembl IDs).
    """
    print("=" * 70)
    print("C2: Reconstruct joint PCA model from 16 TS + 16 PanSci centroids")
    print("=" * 70)

    ts_full = pd.read_csv(TS_CENTROIDS_PATH, index_col=0)
    pansci_full = pd.read_csv(PANSCI_CENTROIDS_PATH, index_col=0)
    assert list(ts_full.columns) == list(pansci_full.columns), \
        "TS and PanSci centroids must share gene column order"
    gene_names = list(ts_full.columns)
    assert len(gene_names) == 16959, f"Expected 16,959 genes, got {len(gene_names)}"

    missing_ts = [t for t in MATCHED_TYPES_16 if t not in ts_full.index]
    missing_ps = [t for t in MATCHED_TYPES_16 if t not in pansci_full.index]
    assert not missing_ts, f"TS missing matched types: {missing_ts}"
    assert not missing_ps, f"PanSci missing matched types: {missing_ps}"

    ts_sub = ts_full.loc[MATCHED_TYPES_16]
    ps_sub = pansci_full.loc[MATCHED_TYPES_16]

    human_pca, mouse_pca, pca, types_list = pca_reduce_centroids(
        ts_sub, ps_sub, variance_threshold=variance_threshold
    )
    return pca, human_pca, mouse_pca, types_list, gene_names


def verify_layer1_reproduction(human_pca, mouse_pca, types_list):
    """Reproduce the Layer-1 Procrustes result and compare to the existing JSON."""
    print("\n" + "=" * 70)
    print("C2: Layer-1 reproduction sanity check")
    print("=" * 70)

    result = procrustes_align(human_pca, mouse_pca)
    p_val, null_dist = permutation_test(
        human_pca, mouse_pca, N_PERM, RANDOM_SEED
    )
    obs_null = result.distance / np.median(null_dist)
    det_R = float(np.linalg.det(result.rotation))

    with EXISTING_REPL_JSON.open() as fh:
        existing = json.load(fh)
    e = existing["procrustes"]
    expected_obs_null = e["obs_null_ratio"]
    expected_p = e["p_value"]
    expected_dist = e["distance"]
    expected_scaling = e["scaling"]
    expected_n_pc = e["pca_components"]

    # Variance retained is checked by the caller against pca.explained_variance_ratio_.
    drift = {
        "obs_null": abs(obs_null - expected_obs_null) > 1e-9,
        "distance": abs(result.distance - expected_dist) > 1e-9,
        "scaling": abs(result.scaling - expected_scaling) > 1e-9,
        "n_pca": human_pca.shape[1] != expected_n_pc,
        "p_value": p_val > 1e-4,  # primary required p < 1e-4
        "det_R": abs(det_R - 1.0) > 1e-9,
    }
    any_drift = any(drift.values())

    print(f"\n  Layer-1 reproduction:")
    print(f"    obs/null ratio  : {obs_null:.15f}  (expected {expected_obs_null:.15f})")
    print(f"    Procrustes dist : {result.distance:.15f}  (expected {expected_dist:.15f})")
    print(f"    Scaling s       : {result.scaling:.15f}  (expected {expected_scaling:.15f})")
    print(f"    PCA components  : {human_pca.shape[1]}  (expected {expected_n_pc})")
    print(f"    p-value         : {p_val:.6f}  (expected < 1e-4)")
    print(f"    det(R)          : {det_R:+.15f}  (expected +1)")
    if any_drift:
        print(f"\n  *** DRIFT DETECTED ***  fields={drift}")
        raise RuntimeError(f"Layer-1 reproduction DRIFT: {drift}")

    return result, p_val, null_dist, obs_null


# ---------------------------------------------------------------------------
# C3: Per-cell projection helpers
# ---------------------------------------------------------------------------
def _project_ts_cells(pca, gene_names):
    """Project TS h5ad cells (filtered to MATCHED_TYPES_16) into joint PCA space."""
    print("\n  --- TS (human) per-cell projection ---")
    adata = ad.read_h5ad(TS_H5AD_PATH)
    print(f"    TS h5ad: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

    assert set(gene_names).issubset(set(adata.var_names)), \
        "TS h5ad missing some ortholog gene IDs"

    type_mask = adata.obs["cell_type"].isin(MATCHED_TYPES_16).values
    n_kept = int(type_mask.sum())
    print(f"    Filter to 16 matched types: -> {n_kept:,} cells")

    X = adata[type_mask, gene_names].X
    if sp.issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float32)

    X_pca = (X - pca.mean_.astype(np.float32)) @ pca.components_.T.astype(np.float32)
    ct = adata.obs.loc[type_mask, "cell_type"].astype(str).values

    counts = pd.Series(ct).value_counts().to_dict()
    print(f"    Projected: {X_pca.shape[0]:,} cells x {X_pca.shape[1]} components")
    for t in MATCHED_TYPES_16:
        print(f"      TS  {t:<50} n={counts.get(t, 0):,}")
    return {"X_pca": X_pca, "cell_types": ct, "counts": counts}


def _project_pansci_cells(pca, gene_names):
    """Re-run the PanSci pre-processing pipeline cell-by-cell, project to PCA.

    Mirrors pansci_replication.py steps 1-3 (WT/06_months filter, QC, CD4/CD8
    split, CPM+log1p, map to 16,959 ortholog space), but emits per-cell
    projected coordinates rather than centroids. Endothelial cell PanSci data
    is restricted to ENDO_TISSUE ("lung") only.
    """
    print("\n  --- PanSci (mouse) per-cell projection ---")
    ortho = pd.read_csv(ORTHOLOG_PATH)
    mouse_to_human = dict(zip(ortho["mouse_gene_name"], ortho["human_ensembl_id"]))
    full_gene_idx = {g: i for i, g in enumerate(gene_names)}
    n_full = len(gene_names)

    pca_components = pca.components_.astype(np.float64)  # (k, 16959)
    pca_mean = pca.mean_.astype(np.float64)              # (16959,)
    k = pca_components.shape[0]

    type_chunks = {ct: [] for ct in MATCHED_TYPES_16}  # ct -> list of (n_i, k) arrays
    type_tissue_counts = {ct: {} for ct in MATCHED_TYPES_16}

    for tissue in TISSUES_TO_LOAD:
        print(f"\n    [{tissue}] loading MTX ...")
        mtx_path = PANSCI_DIR / f"{tissue}_genecount.mtx.gz"
        meta_path = PANSCI_DIR / f"{tissue}_df_cell.csv.gz"
        gene_path = PANSCI_DIR / f"{tissue}_df_gene.csv.gz"
        t0 = time.time()
        mtx = sio.mmread(gzip.open(mtx_path, "rb"))
        mtx = sp.csr_matrix(mtx.T)  # cells x genes
        with gzip.open(meta_path, "rt") as fh:
            meta = pd.read_csv(fh)
        with gzip.open(gene_path, "rt") as fh:
            genes_df = pd.read_csv(fh)
        gene_names_ps = list(genes_df["gene_name"])
        assert mtx.shape == (len(meta), len(gene_names_ps))
        print(f"      loaded {mtx.shape[0]:,} cells x {mtx.shape[1]:,} genes "
              f"in {time.time()-t0:.0f}s")

        # WT + 06_months filter
        mask = (meta["genotype"] == "WT") & (meta["age_group"] == "06_months")
        idx = np.where(mask.values)[0]
        mtx = mtx[idx]
        meta = meta.iloc[idx].copy()
        print(f"      WT/06m: {len(meta):,} cells")

        # QC: >=200 genes, <=20% mito
        genes_detected = np.array((mtx > 0).sum(axis=1)).flatten()
        totals = np.array(mtx.sum(axis=1)).flatten()
        mt_mask = np.array([g.startswith("mt-") for g in gene_names_ps])
        if mt_mask.sum() > 0:
            mt_counts = np.array(mtx[:, mt_mask].sum(axis=1)).flatten()
            pct_mt = mt_counts / np.maximum(totals, 1) * 100
        else:
            pct_mt = np.zeros(mtx.shape[0])
        qc = (genes_detected >= 200) & (pct_mt <= 20)
        idx2 = np.where(qc)[0]
        mtx = mtx[idx2]
        meta = meta.iloc[idx2].copy()
        print(f"      QC pass: {len(meta):,} cells")

        # Cell type mapping (PanSci -> 35-type ontology) + CD4/CD8 split.
        meta["base_type"] = meta["main_cell_type_organ"].apply(strip_organ_suffix)
        meta["our_type"] = meta["base_type"].apply(map_pansci_to_ontology)

        t_mask = meta["our_type"] == "T cell"
        if t_mask.any():
            cd4_i = gene_names_ps.index("Cd4") if "Cd4" in gene_names_ps else None
            cd8_i = gene_names_ps.index("Cd8a") if "Cd8a" in gene_names_ps else None
            if cd4_i is not None and cd8_i is not None:
                t_idx = np.where(t_mask.values)[0]
                t_mat = mtx[t_idx]
                cd4 = np.array((t_mat[:, cd4_i] > 0).toarray()).flatten()
                cd8 = np.array((t_mat[:, cd8_i] > 0).toarray()).flatten()
                is_cd4 = cd4 & ~cd8
                is_cd8 = cd8 & ~cd4
                t_obs_idx = meta.loc[t_mask].index
                new_labels = pd.Series("T cell", index=t_obs_idx)
                new_labels.iloc[is_cd4] = "CD4-positive, alpha-beta T cell"
                new_labels.iloc[is_cd8] = "CD8-positive, alpha-beta T cell"
                meta.loc[t_mask, "our_type"] = new_labels.values

        # Normalize CPM + log1p (target sum 1e4 per pansci_replication.py:380).
        row_sums = np.array(mtx.sum(axis=1)).flatten()
        row_sums[row_sums == 0] = 1
        scaling_factors = 1e4 / row_sums
        mtx_norm = mtx.multiply(scaling_factors[:, np.newaxis]).log1p()
        if not isinstance(mtx_norm, sp.csr_matrix):
            mtx_norm = sp.csr_matrix(mtx_norm)
        del mtx
        gc.collect()

        # Map PanSci gene index -> joint-PCA full-gene index (first occurrence).
        seen = set()
        shared_indices = []
        shared_full_indices = []
        for i, g in enumerate(gene_names_ps):
            if g in seen:
                continue
            if g not in mouse_to_human:
                continue
            hid = mouse_to_human[g]
            if hid not in full_gene_idx:
                continue
            seen.add(g)
            shared_indices.append(i)
            shared_full_indices.append(full_gene_idx[hid])
        shared_indices = np.asarray(shared_indices, dtype=np.int64)
        shared_full_indices = np.asarray(shared_full_indices, dtype=np.int64)
        print(f"      ortholog overlap: {len(shared_indices):,} genes "
              f"(out of {n_full:,} joint-PCA columns)")

        # Restrict the normalized matrix to shared mouse genes.
        mtx_shared = mtx_norm[:, shared_indices].tocsr()
        del mtx_norm
        gc.collect()

        # Project each type's cells.  We expand only one type-block at a time
        # to keep peak memory bounded: dense block size <= n_ct * 16,959 floats.
        for ct in MATCHED_TYPES_16:
            ct_mask = meta["our_type"] == ct
            if ct == "endothelial cell" and tissue != ENDO_TISSUE:
                # Endothelial PanSci side: lung only (DECISION-101).
                continue
            n_ct = int(ct_mask.sum())
            if n_ct == 0:
                continue
            type_tissue_counts[ct][tissue] = n_ct
            ct_idx = np.where(ct_mask.values)[0]

            # Chunk processing to bound dense memory (~ 1 GB / chunk).
            chunk_size = max(1, int(8 * 1024**3 / (n_full * 8)))  # ~8 GB cap
            chunk_size = min(chunk_size, 5000)  # firm cap regardless
            for start in range(0, n_ct, chunk_size):
                sub = ct_idx[start:start + chunk_size]
                sparse_block = mtx_shared[sub]
                # Expand sparse (n, n_shared) -> dense (n, n_full)
                n_block = sparse_block.shape[0]
                dense_block = np.zeros((n_block, n_full), dtype=np.float64)
                # CSR per-row gather: faster as full toarray then scatter columns
                sparse_dense = sparse_block.toarray()  # (n_block, n_shared)
                dense_block[:, shared_full_indices] = sparse_dense
                del sparse_dense
                # Project: (X - mean) @ components.T
                proj = (dense_block - pca_mean) @ pca_components.T
                del dense_block
                type_chunks[ct].append(proj.astype(np.float32))
            print(f"      [{tissue}] {ct:<45} n={n_ct:>8,}")

        del mtx_shared
        gc.collect()

    # Concatenate per-type chunks across tissues.
    proj_by_type = {}
    counts = {}
    for ct in MATCHED_TYPES_16:
        if not type_chunks[ct]:
            proj_by_type[ct] = np.zeros((0, k), dtype=np.float32)
            counts[ct] = 0
            continue
        arr = np.vstack(type_chunks[ct])
        proj_by_type[ct] = arr
        counts[ct] = arr.shape[0]

    # Build flat arrays in MATCHED_TYPES_16 order.
    X_pca_all = []
    cell_types_all = []
    for ct in MATCHED_TYPES_16:
        arr = proj_by_type[ct]
        if arr.shape[0] == 0:
            continue
        X_pca_all.append(arr)
        cell_types_all.extend([ct] * arr.shape[0])
    X_pca_all = np.vstack(X_pca_all) if X_pca_all else np.zeros((0, k), dtype=np.float32)
    cell_types_all = np.asarray(cell_types_all, dtype=object)

    print(f"\n    PanSci projected: {X_pca_all.shape[0]:,} cells x {k} components")
    for ct in MATCHED_TYPES_16:
        ts_breakdown = ", ".join(
            f"{t}({c:,})" for t, c in type_tissue_counts[ct].items()
        )
        print(f"      PS  {ct:<50} n={counts[ct]:>8,}  [{ts_breakdown}]")

    return {
        "X_pca": X_pca_all,
        "cell_types": cell_types_all,
        "counts": counts,
        "tissue_breakdown": type_tissue_counts,
    }


def project_cells_pansci(pca, gene_names):
    """Return cell_data dict in the shape compute_covariance_eigen expects."""
    print("\n" + "=" * 70)
    print("C3: Per-cell projection into joint PCA space")
    print("=" * 70)

    ts = _project_ts_cells(pca, gene_names)
    ps = _project_pansci_cells(pca, gene_names)

    # compute_covariance_eigen iterates over species names ["human", "mouse"];
    # we map TS -> "human", PanSci -> "mouse".
    cell_data = {
        "human": {"X_pca": ts["X_pca"], "cell_types": ts["cell_types"]},
        "mouse": {"X_pca": ps["X_pca"], "cell_types": ps["cell_types"]},
    }
    counts = {
        "TS":     {ct: int(ts["counts"].get(ct, 0)) for ct in MATCHED_TYPES_16},
        "PanSci": {ct: int(ps["counts"].get(ct, 0)) for ct in MATCHED_TYPES_16},
        "PanSci_tissue_breakdown": {
            ct: {t: int(n) for t, n in ps["tissue_breakdown"][ct].items()}
            for ct in MATCHED_TYPES_16
        },
    }
    return cell_data, counts


# ---------------------------------------------------------------------------
# C8: per-type S(k=5) vs Layer-1 residual correlation
# ---------------------------------------------------------------------------
def per_type_residual_correlation(alignment_df):
    """Spearman rho between per-type S_pre(k=5) and Layer-1 per-type residual."""
    print("\n" + "=" * 70)
    print("C8: per-type S_pre(k=5) vs Layer-1 residual correlation")
    print("=" * 70)

    with EXISTING_REPL_JSON.open() as fh:
        repl = json.load(fh)
    per_type_res = repl["procrustes"]["per_type_residuals"]

    sub = alignment_df[alignment_df.k == 5].set_index("cell_type")
    rows = []
    for ct in MATCHED_TYPES_16:
        if ct not in sub.index or ct not in per_type_res:
            continue
        rows.append({
            "cell_type": ct,
            "S_pre_k5": float(sub.loc[ct, "S_pre"]),
            "S_post_k5": float(sub.loc[ct, "S_post"]),
            "S_pre_k1": float(alignment_df[(alignment_df.k == 1)
                                         & (alignment_df.cell_type == ct)].iloc[0].S_pre),
            "S_post_k1": float(alignment_df[(alignment_df.k == 1)
                                          & (alignment_df.cell_type == ct)].iloc[0].S_post),
            "S_pre_k3": float(alignment_df[(alignment_df.k == 3)
                                         & (alignment_df.cell_type == ct)].iloc[0].S_pre),
            "S_post_k3": float(alignment_df[(alignment_df.k == 3)
                                          & (alignment_df.cell_type == ct)].iloc[0].S_post),
            "layer1_residual_magnitude": float(per_type_res[ct]["magnitude"]),
        })
    df = pd.DataFrame(rows)
    rho, p = stats.spearmanr(df["S_pre_k5"], df["layer1_residual_magnitude"])
    print(f"\n  n = {len(df)}")
    print(f"  Spearman rho = {rho:+.4f}, p = {p:.4f}")
    return df, float(rho), float(p)


# ---------------------------------------------------------------------------
# Output serialisation
# ---------------------------------------------------------------------------
def _layer_block(perm_result):
    """Pack a per-k permutation result into the summary schema."""
    out = {}
    for k in K_VALUES:
        r = perm_result[k]
        null = r["null"]
        out[f"k{k}"] = {
            "S": float(r["observed"]),
            "null_mean": float(null.mean()),
            "null_sd": float(null.std()),
            "p": float(r["p_value"]),
            "margin": float(r["observed"] - null.mean()),
        }
    return out


def save_outputs(
    layer1, alignment_df, perm_pre, perm_post,
    per_type_corr_df, per_type_rho, per_type_p,
    cpc_results, counts, gene_names,
):
    print("\n" + "=" * 70)
    print("Saving outputs to", OUTPUT_DIR)
    print("=" * 70)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)

    # --- per_type_S_pansci.csv ---
    per_type_corr_df.to_csv(OUTPUT_DIR / "per_type_S_pansci.csv", index=False)

    # --- cpc1_drivers_pansci.csv + classification ---
    cpc_rows = []
    ribo_types = []
    cts_types = []
    for ct in MATCHED_TYPES_16:
        if ct not in cpc_results:
            continue
        top5 = cpc_results[ct]["top5_genes"]
        rank1 = top5[0][0] if top5 else ""
        classification = "ribosomal" if _is_ribosomal(rank1) else "cell_type_specific"
        if classification == "ribosomal":
            ribo_types.append(ct)
        else:
            cts_types.append({"type": ct, "rank1_driver": rank1})
        cpc_rows.append({
            "cell_type": ct,
            "rank1_driver": rank1,
            "top5_drivers": ";".join(g for g, _ in top5),
            "classification": classification,
            "cpc1_var_frac_TS": cpc_results[ct]["cpc1_var_frac_human"],
            "cpc1_var_frac_PanSci": cpc_results[ct]["cpc1_var_frac_mouse"],
        })
    cpc_df = pd.DataFrame(cpc_rows)
    cpc_df.to_csv(OUTPUT_DIR / "cpc1_drivers_pansci.csv", index=False)

    # --- null_distributions_pansci.npz ---
    null_arrays = {}
    for k in K_VALUES:
        null_arrays[f"pre_k{k}"] = perm_pre[k]["null"]
        null_arrays[f"post_k{k}"] = perm_post[k]["null"]
    np.savez_compressed(OUTPUT_DIR / "null_distributions_pansci.npz", **null_arrays)

    # --- pansci_layer2_summary.json ---
    summary = {
        "layer1_reproduction": {
            "obs_over_null": float(layer1["obs_null"]),
            "p_value": float(layer1["p_value"]),
            "n_pca_components": int(layer1["n_pc"]),
            "variance_retained": float(layer1["variance_retained"]),
            "procrustes_scale_s": float(layer1["scaling"]),
            "R_determinant": float(layer1["det_R"]),
        },
        "layer2_pre_rotation": _layer_block(perm_pre),
        "layer2_post_rotation": _layer_block(perm_post),
        "per_type_correlation": {
            "spearman_rho": per_type_rho,
            "p_value": per_type_p,
            "n": int(len(per_type_corr_df)),
        },
        "cpc1_classification": {
            "ribosomal_dominated_count": len(ribo_types),
            "ribosomal_dominated_types": ribo_types,
            "cell_type_specific_count": len(cts_types),
            "cell_type_specific": cts_types,
        },
        "cell_counts": {
            "TS_per_type": counts["TS"],
            "PanSci_per_type": counts["PanSci"],
            "PanSci_tissue_breakdown": counts["PanSci_tissue_breakdown"],
        },
        "qc_acknowledgments": {
            "endothelial_tissue_restriction": ENDO_TISSUE,
            "w1_low_proxy_types": [],
            "w2_single_tissue_types": [],
            "_note": "W-1 and W-2 lists are merged separately (not regenerated here).",
        },
        "config": {
            "random_seed": SEED,
            "n_permutations": N_PERM,
            "k_values": K_VALUES,
            "variance_threshold": VARIANCE_THRESHOLD,
        },
    }
    with (OUTPUT_DIR / "pansci_layer2_summary.json").open("w") as fh:
        json.dump(summary, fh, indent=2)

    # --- krzanowski_verification_pansci.md (human summary) ---
    md = []
    md.append("# Krzanowski verification — PanSci Layer-2\n")
    md.append("Generated by `scripts/t3b_ellipsoid_alignment_pansci.py` (R15 Part A.2).\n")
    md.append("## Layer-1 reproduction\n")
    md.append(f"- obs/null ratio: **{layer1['obs_null']:.6f}** "
              f"(existing: {layer1['expected_obs_null']:.6f})")
    md.append(f"- Procrustes distance: {layer1['distance']:.6f}")
    md.append(f"- p-value: {layer1['p_value']:.6f}  (< 1e-4 required)")
    md.append(f"- PCA components: {layer1['n_pc']}, variance retained: {layer1['variance_retained']:.4f}")
    md.append(f"- Procrustes scale s: {layer1['scaling']:.6f}")
    md.append(f"- det(R) = {layer1['det_R']:+.6e}")
    md.append("")
    md.append("## Layer-2 aggregate S (mean over 16 types)\n")
    md.append("| k | S_pre | null_pre_mean | p_pre | S_post | null_post_mean | p_post |")
    md.append("|---|-------|---------------|-------|--------|----------------|--------|")
    for k in K_VALUES:
        pre, post = perm_pre[k], perm_post[k]
        md.append(
            f"| {k} | {pre['observed']:.4f} | {pre['null'].mean():.4f} | {pre['p_value']:.4g} "
            f"| {post['observed']:.4f} | {post['null'].mean():.4f} | {post['p_value']:.4g} |"
        )
    md.append("")
    md.append("## Per-type S vs Layer-1 residual\n")
    md.append(f"- Spearman rho = **{per_type_rho:+.4f}** (p = {per_type_p:.4f}, n = {len(per_type_corr_df)})")
    md.append("")
    md.append("## CPC1 classification\n")
    md.append(f"- Ribosomal-dominated: **{len(ribo_types)} / 16**")
    md.append(f"  - Types: {', '.join(ribo_types) if ribo_types else '(none)'}")
    md.append(f"- Cell-type-specific: **{len(cts_types)} / 16**")
    for entry in cts_types:
        md.append(f"  - {entry['type']}: rank-1 driver = `{entry['rank1_driver']}`")
    (OUTPUT_DIR / "krzanowski_verification_pansci.md").write_text("\n".join(md) + "\n")

    # Sizes
    for name in [
        "krzanowski_verification_pansci.md",
        "per_type_S_pansci.csv",
        "cpc1_drivers_pansci.csv",
        "null_distributions_pansci.npz",
        "pansci_layer2_summary.json",
    ]:
        p = OUTPUT_DIR / name
        print(f"  {p}  ({p.stat().st_size:,} bytes)")


def _is_ribosomal(symbol: str) -> bool:
    import re
    return bool(re.match(r"^(RPL|RPS|Rpl|Rps)\d", symbol or ""))


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------
def main_pansci():
    t0 = time.time()

    # C2: PCA reconstruction + Layer-1 sanity ---------------------------------
    pca, human_pca, mouse_pca, types_list, gene_names = reconstruct_pca_pansci()
    variance_retained = float(np.sum(pca.explained_variance_ratio_))
    n_pc = int(pca.n_components_)

    # Persist R, s, types_list to scratch BEFORE permutation runs to make a
    # mid-run crash auditable.
    result, p_val, null_dist_l1, obs_null = verify_layer1_reproduction(
        human_pca, mouse_pca, types_list
    )
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    np.save(SCRATCH_DIR / "R_centroid_optimal.npy", result.rotation)
    np.save(SCRATCH_DIR / "human_pca_centroids.npy", human_pca)
    np.save(SCRATCH_DIR / "mouse_pca_centroids.npy", mouse_pca)
    with (SCRATCH_DIR / "layer1_meta.json").open("w") as fh:
        json.dump({
            "scaling": float(result.scaling),
            "distance": float(result.distance),
            "obs_null": float(obs_null),
            "p_value": float(p_val),
            "det_R": float(np.linalg.det(result.rotation)),
            "n_pc": n_pc,
            "variance_retained": variance_retained,
            "types_list": types_list,
        }, fh, indent=2)

    with EXISTING_REPL_JSON.open() as fh:
        existing = json.load(fh)
    layer1 = {
        "obs_null": obs_null,
        "expected_obs_null": existing["procrustes"]["obs_null_ratio"],
        "p_value": p_val,
        "distance": float(result.distance),
        "scaling": float(result.scaling),
        "det_R": float(np.linalg.det(result.rotation)),
        "n_pc": n_pc,
        "variance_retained": variance_retained,
    }

    # C3: per-cell projection -------------------------------------------------
    cell_data, counts = project_cells_pansci(pca, gene_names)
    np.savez_compressed(
        SCRATCH_DIR / "per_cell_projections.npz",
        ts_X_pca=cell_data["human"]["X_pca"],
        ts_cell_types=cell_data["human"]["cell_types"].astype(str),
        ps_X_pca=cell_data["mouse"]["X_pca"],
        ps_cell_types=cell_data["mouse"]["cell_types"].astype(str),
    )

    # C4: per-type covariance + eigendecomposition ----------------------------
    eigen_results = compute_covariance_eigen(cell_data, scale_label="pansci16type")

    # C5+C6: Krzanowski S pre and post ---------------------------------------
    alignment_df = compute_alignment_scores(
        eigen_results, result.rotation, scale_label="pansci16type"
    )

    # C7: permutation nulls (one call -> both pre and post share perm indices) -
    perm_pre, perm_post = perm_test_label_shuffle(
        eigen_results, result.rotation, n_perm=N_PERM
    )

    # C8: per-type S vs Layer-1 residual --------------------------------------
    per_type_corr_df, per_type_rho, per_type_p = per_type_residual_correlation(
        alignment_df
    )

    # C9: CPC1 driver scan ----------------------------------------------------
    cpc_results = common_principal_components(
        eigen_results, pca, gene_names, scale_label="pansci16type"
    )

    # Persist all outputs ------------------------------------------------------
    save_outputs(
        layer1, alignment_df, perm_pre, perm_post,
        per_type_corr_df, per_type_rho, per_type_p,
        cpc_results, counts, gene_names,
    )

    elapsed = (time.time() - t0) / 60
    print(f"\nDONE in {elapsed:.1f} min. Review output values before use.")


if __name__ == "__main__":
    main_pansci()
