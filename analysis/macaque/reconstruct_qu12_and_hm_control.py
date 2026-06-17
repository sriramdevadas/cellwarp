#!/usr/bin/env python3
"""Expanded Option K: 12-type human-macaque (Qu-only, symmetric) + 12-type
human-mouse control. Locks in the D1 conservative+moderate expansion with
n=12 = D1's 7 + {B cell, macrophage, plasma cell, neutrophil, monocyte}.

Macaque pipeline: Qu raw 10x → normalize_total(1e4) → log1p → per-type
centroids at ≤ 2,000 cells/type (seed 42). Same convention as D1. No
per-organ cap — see docstring "Per-organ subsampling decision" below.

Human-mouse control: subset committed centroids_{human,mouse}_35.csv
to the same 12 types. Same Procrustes convention.

Per-organ subsampling decision
------------------------------
Two of the new types (neutrophil, plasma cell) are organ-concentrated
in Qu: neutrophil 49% Liver, plasma cell 71% GI. Applying a per-organ
cap would partially flatten the within-Qu tissue composition, but:
  (a) D1 did not use a per-organ cap, so applying it here breaks
      methodological consistency with the already-reported D1 result;
  (b) The committed Tabula Sapiens centroids (centroids_human_35.csv)
      were computed from a 2,000-cell-per-type random subsample of the
      full TS cell pool — i.e., TS's own tissue composition is
      preserved, not organ-capped. TS and Qu tissue distributions
      differ intrinsically (TS has lymphoid/blood-heavy sampling for
      immune types; Qu covers 16 solid tissues), and a per-organ cap
      within Qu does not align Qu to TS — it just flattens Qu.
  (c) The pipeline convention for this project is uniform per-type
      subsample (MAX_CELLS_PER_TYPE), applied identically to all types.

We therefore apply the uniform 2,000-cell-per-type cap and flag the
neutrophil/plasma cell organ concentrations in the report. Readers who
want a per-organ-balanced sensitivity can rerun with a modified cap;
the committed driver supports it via a flag.
"""
from __future__ import annotations
import json, sys, time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

PROJECT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT / "src"))
from procrustes import pca_reduce_centroids, procrustes_align, permutation_test  # noqa

RANDOM_SEED = 42
N_PERMUTATIONS = 10_000
MAX_CELLS_PER_TYPE = 2_000

QU_META_CSV = PROJECT / "data/macaque/qu_2022/qu_metadata.csv"
AGG_H5AD    = PROJECT / "output/macaque_pipeline/qu_dry_run_aggregate.h5ad"
HUMAN_CSV   = PROJECT / "output/phase2/scaled_35types/centroids_human_35.csv"
MOUSE_CSV   = PROJECT / "output/phase2/scaled_35types/centroids_mouse_35.csv"

OUT_DIR     = PROJECT / "output/macaque_pipeline"
OUT_QU_JSON = OUT_DIR / "reconstruction_qu12_results.json"
OUT_HM_JSON = OUT_DIR / "human_mouse_12type_control.json"
OUT_QU_CENTROIDS = OUT_DIR / "reconstruction_qu12_centroids.csv"

# 12 target types — D1 baseline 7 + conservative clean fine-matches 4 + monocyte generic 1
TYPES_12 = [
    "B cell",
    "bladder urothelial cell",
    "endothelial cell",
    "epithelial cell",
    "fibroblast",
    "hepatocyte",
    "macrophage",
    "monocyte",
    "neutrophil",
    "plasma cell",
    "smooth muscle cell",
    "stromal cell",
]


def harmonize_qu_labels(meta: pd.DataFrame) -> pd.Series:
    """Map Qu native author labels → 12 CellWarp target types.

    Mapping table (fine-cluster-first, stromal fallback):
      bladder urothelial cell  ← finalcluster=Epithelial cell, Organ=Bladder
      hepatocyte               ← finalcluster=Epithelial cell, Organ=Liver
      epithelial cell          ← finalcluster=Epithelial cell, Organ ∉ {Bladder, Liver}
      endothelial cell         ← finalcluster=Endothelial cell
      fibroblast               ← finalcluster=Fibroblasts
      smooth muscle cell       ← finalcluster=Smooth muscle cell
      B cell                   ← finalcluster=B cell
      macrophage               ← finalcluster=Macrophage
      plasma cell              ← finalcluster=Plasma cell
      neutrophil               ← finalcluster=Neutrophil
      monocyte                 ← finalcluster=Monocyte (generic CL:0000576; not split into
                                  classical/intermediate/non-classical)
      stromal cell             ← majorcluster=Stromal cell (residual after specific rules)
    """
    fc = meta["finalcluster"].astype(str)
    mc = meta["majorcluster"].astype(str)
    org = meta["Organ"].astype(str)
    out = pd.Series(pd.NA, index=meta.index, dtype="object")

    # Epithelial is split by organ
    out[(fc == "Epithelial cell") & (org == "Bladder")] = "bladder urothelial cell"
    out[(fc == "Epithelial cell") & (org == "Liver")] = "hepatocyte"
    out[(fc == "Epithelial cell") & ~org.isin(["Bladder", "Liver"])] = "epithelial cell"
    # Other fine-label 1:1 mappings
    out[fc == "Endothelial cell"] = "endothelial cell"
    out[fc == "Fibroblasts"] = "fibroblast"
    out[fc == "Smooth muscle cell"] = "smooth muscle cell"
    out[fc == "B cell"] = "B cell"
    out[fc == "Macrophage"] = "macrophage"
    out[fc == "Plasma cell"] = "plasma cell"
    out[fc == "Neutrophil"] = "neutrophil"
    out[fc == "Monocyte"] = "monocyte"
    # Stromal fallback: claim unclaimed majorcluster=Stromal
    remaining = out.isna() & (mc == "Stromal cell")
    out[remaining] = "stromal cell"
    return out


def build_sample_tissue_crosswalk(sample_vals, bcs, adata) -> dict[str, str]:
    """Meta.sample → h5ad.tissue via max barcode overlap (same pattern as D1)."""
    h5ad_bc = adata.obs["barcode"].astype(str).to_numpy()
    h5ad_tt = adata.obs["tissue"].astype(str).to_numpy()
    unique_tissues = sorted(set(h5ad_tt))
    crosswalk = {}
    for s in pd.unique(sample_vals):
        meta_bc_for_sample = set(bcs[sample_vals == s])
        best, best_ov = None, 0
        for tt in unique_tissues:
            tt_bcs = set(h5ad_bc[h5ad_tt == tt])
            ov = len(meta_bc_for_sample & tt_bcs)
            if ov > best_ov:
                best_ov = ov
                best = tt
        crosswalk[s] = best
    return crosswalk


def run_procrustes_pair(X_df, Y_df, label):
    """Run joint PCA + Procrustes + 10k perms; return summary dict."""
    Xp, Yp, pca, order = pca_reduce_centroids(X_df, Y_df)
    r = procrustes_align(Xp, Yp)
    p, null = permutation_test(Xp, Yp, n_permutations=N_PERMUTATIONS, seed=RANDOM_SEED)
    obs_null_med = r.distance / np.median(null)
    obs_null_mean = r.distance / np.mean(null)
    residuals = []
    for i, typ in enumerate(order):
        vec = r.aligned_target[i] - r.centered_reference[i]
        mag = float(np.linalg.norm(vec))
        pct = 100.0 * mag**2 / r.distance_squared if r.distance_squared > 0 else 0.0
        residuals.append({"type": typ, "magnitude": mag, "pct_ssr": pct})
    residuals_sorted = sorted(residuals, key=lambda x: -x["magnitude"])
    return {
        "label": label,
        "order": order,
        "procrustes": {
            "distance": float(r.distance),
            "scaling": float(r.scaling),
            "distance_squared": float(r.distance_squared),
            "rotation_det": float(np.linalg.det(r.rotation)),
        },
        "permutation_test": {
            "p_value": float(p),
            "obs_null_ratio_median": float(obs_null_med),
            "obs_null_ratio_mean": float(obs_null_mean),
            "null_median": float(np.median(null)),
            "null_mean": float(np.mean(null)),
            "n_permutations": N_PERMUTATIONS,
            "permuted_le_observed": int(np.sum(null <= r.distance)),
        },
        "pca": {
            "n_components": int(pca.n_components_),
            "variance_explained": pca.explained_variance_ratio_.tolist(),
            "cumulative_variance": float(np.cumsum(pca.explained_variance_ratio_)[-1]),
        },
        "per_type_residuals_ranked": residuals_sorted,
    }


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("12-type expanded: human-macaque (Qu) + human-mouse control")
    print("=" * 70)

    # ── 1. Load Qu metadata, harmonize, filter
    print("\n[1] Loading Qu metadata + harmonizing labels")
    qu_meta = pd.read_csv(QU_META_CSV, low_memory=False)
    qu_meta["target_type"] = harmonize_qu_labels(qu_meta)

    print("\n  Per-type available cells (Qu, post-harmonize):")
    available = {}
    for t in TYPES_12:
        n = int((qu_meta["target_type"] == t).sum())
        available[t] = n
        print(f"    {t:<30s} {n:>7,}")

    # ── 2. Barcode matching via sample → tissue crosswalk
    print("\n[2] Loading aggregate h5ad + barcode crosswalk")
    adata = ad.read_h5ad(AGG_H5AD)
    sep = "@@_"
    ids = qu_meta["cell_id"].astype(str).to_numpy()
    bcs = np.array([x.split(sep, 1)[1] if sep in x else x for x in ids])
    samples = qu_meta["sample"].astype(str).to_numpy()
    crosswalk = build_sample_tissue_crosswalk(samples, bcs, adata)
    remapped = np.array([crosswalk.get(s, "UNKNOWN") for s in samples])
    meta_keys = np.array([f"{b}__{t}" for b, t in zip(bcs, remapped)])
    h5ad_keys = [f"{bc}__{tt}" for bc, tt in
                 zip(adata.obs["barcode"].astype(str), adata.obs["tissue"].astype(str))]
    key_to_row = {k: i for i, k in enumerate(h5ad_keys)}
    matched_idx, matched_meta = [], []
    for i, k in enumerate(meta_keys):
        if k in key_to_row:
            matched_idx.append(key_to_row[k])
            matched_meta.append(i)
    matched_idx = np.array(matched_idx, dtype=np.int64)
    matched_meta = np.array(matched_meta, dtype=np.int64)
    print(f"  Matched {len(matched_idx):,} / {len(meta_keys):,} metadata cells")

    types_assigned = qu_meta["target_type"].to_numpy()[matched_meta]
    included = set(TYPES_12)
    type_mask = np.array([
        (t is not None) and (not (isinstance(t, float) and np.isnan(t))) and t in included
        for t in types_assigned
    ])
    final_idx = matched_idx[type_mask]
    final_types = types_assigned[type_mask].astype(str)
    adata = adata[final_idx].copy()
    adata.obs["target_type"] = final_types
    print(f"  After type filter: {adata.n_obs:,} cells in 12 target types")

    # ── 3. Subsample
    print(f"\n[3] Subsampling 2,000/type (seed {RANDOM_SEED})")
    rng = np.random.default_rng(RANDOM_SEED)
    keep = np.zeros(adata.n_obs, dtype=bool)
    used = {}
    for t in TYPES_12:
        idx = np.where((adata.obs["target_type"] == t).to_numpy())[0]
        sel = rng.choice(idx, size=min(MAX_CELLS_PER_TYPE, len(idx)), replace=False) if len(idx) > MAX_CELLS_PER_TYPE else idx
        keep[sel] = True
        used[t] = int(sel.size if hasattr(sel, "size") else len(sel))
        print(f"  {t:<30s} {len(idx):>7,} → {used[t]:>5,}")
    adata = adata[keep].copy()

    # ── 4. Normalize + centroids (macaque)
    print("\n[4] normalize_total(1e4) + log1p on 12-type filtered Qu cells")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    mac = np.zeros((12, adata.n_vars), dtype=np.float64)
    for i, t in enumerate(TYPES_12):
        m = (adata.obs["target_type"] == t).to_numpy()
        mac[i] = np.asarray(adata.X[m].mean(axis=0)).flatten()
    mac_df = pd.DataFrame(mac, index=TYPES_12, columns=adata.var_names)
    mac_df.to_csv(OUT_QU_CENTROIDS)

    # ── 5. Human + mouse subsets
    print("\n[5] Subsetting human and mouse centroids")
    human_full = pd.read_csv(HUMAN_CSV, index_col=0)
    mouse_full = pd.read_csv(MOUSE_CSV, index_col=0)
    missing_h = [t for t in TYPES_12 if t not in human_full.index]
    missing_m = [t for t in TYPES_12 if t not in mouse_full.index]
    if missing_h or missing_m:
        raise SystemExit(f"Missing types — human:{missing_h} mouse:{missing_m}")
    gene_ids = list(mac_df.columns)  # 13,927 three-way
    human_12 = human_full.loc[TYPES_12, gene_ids].copy()
    # For H-M control: use human-mouse 2-way space (16,959) since mouse/human
    # are already matched in that space, and to mirror the 7-type control.
    human_12_hm = human_full.loc[TYPES_12].copy()
    mouse_12_hm = mouse_full.loc[TYPES_12].copy()

    # ── 6. Procrustes: human-macaque 12-type (13,927 gene space)
    print("\n[6a] Human-macaque 12-type Procrustes (13,927 three-way gene space)")
    res_mac = run_procrustes_pair(human_12, mac_df, "human-macaque-12")

    # ── 7. Procrustes: human-mouse 12-type (16,959 gene space, matches 7-type control)
    print("\n[6b] Human-mouse 12-type control (16,959 human-mouse gene space)")
    res_mm = run_procrustes_pair(human_12_hm, mouse_12_hm, "human-mouse-12")

    # Also do human-mouse 12-type in the 13,927 space for direct apples-to-apples
    print("\n[6c] Human-mouse 12-type in 13,927-gene three-way space (for parity with macaque)")
    mouse_12_3way = mouse_full.loc[TYPES_12, gene_ids].copy()
    res_mm_3way = run_procrustes_pair(human_12, mouse_12_3way, "human-mouse-12-3way-gene-space")

    # ── 8. Per-organ concentration reporting for the new types
    org_conc = {}
    for t in ["B cell", "macrophage", "plasma cell", "neutrophil", "monocyte"]:
        sub = qu_meta[qu_meta["target_type"] == t]
        cnt = sub["Organ"].value_counts().head(5).to_dict()
        total = len(sub)
        top_frac = max(cnt.values()) / total if total else 0.0
        org_conc[t] = {"total": total, "top5_organs": cnt, "top_organ_fraction": top_frac}

    out_combined = {
        "analysis": "QU12_AND_HM12_CONTROL",
        "types_included": TYPES_12,
        "gene_space_macaque": 13927,
        "gene_space_human_mouse": 16959,
        "preprocessing": {
            "qu": "raw UMI → normalize_total(1e4) → log1p; max 2,000 cells/type (seed 42)",
            "human_tabula_sapiens": "committed centroids_human_35.csv (symmetric: normalize_total+log1p)",
            "mouse_tabula_muris_senis": "committed centroids_mouse_35.csv (symmetric)",
        },
        "per_organ_subsampling": "NOT applied (uniform 2k/type; see driver docstring for rationale)",
        "monocyte_mapping": "Qu Monocyte → TS monocyte (CL:0000576, generic); not split into classical/intermediate/non-classical",
        "per_type_cells": {
            "available_post_harmonize": available,
            "used_after_subsample": used,
        },
        "organ_concentration_flagged_types": org_conc,
        "human_macaque_12": res_mac,
        "human_mouse_12_control": res_mm,
        "human_mouse_12_three_way_gene_space": res_mm_3way,
        "seed": RANDOM_SEED,
        "committed_context": {
            "primary_human_mouse_35type_obs_null": 0.522,
            "primary_human_macaque_20type_obs_null": 0.841,
            "sensitivity_rira_13type_obs_null": 0.7488,
            "d1_pilot_human_macaque_7type_obs_null": 0.7334,
            "d1_pilot_human_macaque_7type_p": 0.0125,
            "human_mouse_7type_control_obs_null": 0.4848,
            "human_mouse_7type_control_p": 0.00060,
        },
    }

    def _default(o):
        if isinstance(o, np.bool_): return bool(o)
        if isinstance(o, np.integer): return int(o)
        if isinstance(o, np.floating): return float(o)
        if isinstance(o, np.ndarray): return o.tolist()
        raise TypeError(f"not serialisable: {type(o).__name__}")
    # Split results into two JSON files per user request pattern
    OUT_QU_JSON.write_text(json.dumps({
        "analysis": "QU12_HUMAN_MACAQUE",
        "types_included": TYPES_12,
        "gene_space": 13927,
        **res_mac,
        "per_type_cells": {"available": available, "used": used},
        "organ_concentration_flagged_types": org_conc,
        "monocyte_mapping_note": "CL:0000576 generic; not subtype-split",
    }, indent=2, default=_default))
    OUT_HM_JSON.write_text(json.dumps({
        "analysis": "HUMAN_MOUSE_12TYPE_CONTROL",
        "types_included": TYPES_12,
        "gene_space": 16959,
        "control_16959": res_mm,
        "control_13927_parity": res_mm_3way,
    }, indent=2, default=_default))

    # Console summary
    def _fmt(r):
        p = r["permutation_test"]
        proc = r["procrustes"]
        return (f"d={proc['distance']:.3f}  s={proc['scaling']:.4f}  "
                f"obs/null(med)={p['obs_null_ratio_median']:.4f}  "
                f"p={p['p_value']:.6f}  "
                f"perm≤obs={p['permuted_le_observed']}/{p['n_permutations']}  "
                f"k={r['pca']['n_components']}")

    print("\n" + "=" * 70)
    print("RESULTS SUMMARY")
    print("=" * 70)
    print(f"  [Qu 12-type (13,927 genes)]   {_fmt(res_mac)}")
    print(f"  [H-M 12-type (16,959 genes)]  {_fmt(res_mm)}")
    print(f"  [H-M 12-type (13,927 genes)]  {_fmt(res_mm_3way)}")

    print("\n  Per-type residuals — Qu (high→low):")
    for d in res_mac["per_type_residuals_ranked"]:
        print(f"    {d['type']:<30s}  mag={d['magnitude']:.3f}  pct_ssr={d['pct_ssr']:.1f}%")
    print("\n  Per-type residuals — H-M 12-type 16,959 (high→low):")
    for d in res_mm["per_type_residuals_ranked"]:
        print(f"    {d['type']:<30s}  mag={d['magnitude']:.3f}  pct_ssr={d['pct_ssr']:.1f}%")

    print("\n  Organ concentration for new types:")
    for t, info in org_conc.items():
        top_org, top_n = next(iter(info["top5_organs"].items()))
        print(f"    {t:<15s}  top={top_org} ({info['top_organ_fraction']:.0%} of {info['total']:,})")


if __name__ == "__main__":
    main()
