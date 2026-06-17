#!/usr/bin/env python3
"""Test 1-Qu-pilot (D1) — Qu-only macaque extension, 7 non-immune types.

Uses Qu et al. 2022 author-provided labels (finalcluster + majorcluster
+ Organ) with the exact mapping specified by the user for this pilot.
Symmetric pipeline from raw: Qu 10x UMIs → normalize_total(1e4) → log1p
→ per-type centroids; human X from committed centroids_human_35.csv.

Gene space: 13,927 three-way ortholog intersection (as used in the
committed primary macaque and in our RIRA-13 reconstruction). Pairwise
human-macaque (~14,558) was considered; we kept 13,927 for consistency
with prior analyses and because the 631-gene difference is unlikely to
materially affect an n=7 Procrustes at ≥95% variance PCA truncation.

Author-label → target mapping (per user spec):
  bladder urothelial cell ← finalcluster=Epithelial cell, Organ=Bladder
  hepatocyte              ← finalcluster=Epithelial cell, Organ=Liver
  endothelial cell        ← finalcluster=Endothelial cell (all organs)
  epithelial cell         ← finalcluster=Epithelial cell (all organs)
  fibroblast              ← finalcluster=Fibroblasts
  smooth muscle cell      ← finalcluster=Smooth muscle cell
  stromal cell            ← majorcluster=Stromal cell
"""
from __future__ import annotations
import json, sys, time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

PROJECT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT / "src"))
from procrustes import pca_reduce_centroids, procrustes_align, permutation_test  # noqa

RANDOM_SEED = 42
N_PERMUTATIONS = 10_000
MAX_CELLS_PER_TYPE = 2_000

QU_META_CSV = PROJECT / "data/macaque/qu_2022/qu_metadata.csv"
AGG_H5AD    = PROJECT / "output/macaque_pipeline/qu_dry_run_aggregate.h5ad"
HUMAN_CENTROIDS_CSV = PROJECT / "output/phase2/scaled_35types/centroids_human_35.csv"

OUT_DIR = PROJECT / "output/macaque_pipeline"
OUT_JSON = OUT_DIR / "reconstruction_qu7_D1_results.json"
OUT_CENTROIDS = OUT_DIR / "reconstruction_qu7_D1_centroids.csv"

QU_TARGET_TYPES = [
    "bladder urothelial cell",
    "endothelial cell",
    "epithelial cell",
    "fibroblast",
    "hepatocyte",
    "smooth muscle cell",
    "stromal cell",
]


def harmonize_qu_author_labels(meta: pd.DataFrame) -> pd.Series:
    """Map Qu native author labels → 7 target types per user-specified table."""
    fc = meta["finalcluster"].astype(str)
    mc = meta["majorcluster"].astype(str)
    org = meta["Organ"].astype(str)
    out = pd.Series(pd.NA, index=meta.index, dtype="object")

    # Most-specific first so that "bladder urothelial" and "hepatocyte"
    # are assigned before "epithelial cell" catches everything else.
    out[(fc == "Epithelial cell") & (org == "Bladder")] = "bladder urothelial cell"
    out[(fc == "Epithelial cell") & (org == "Liver")] = "hepatocyte"
    # "epithelial cell" = Epithelial cells not claimed by bladder/liver
    epi_mask = (fc == "Epithelial cell") & ~org.isin(["Bladder", "Liver"])
    out[epi_mask] = "epithelial cell"
    out[fc == "Endothelial cell"] = "endothelial cell"
    out[fc == "Fibroblasts"] = "fibroblast"
    out[fc == "Smooth muscle cell"] = "smooth muscle cell"
    # Stromal = majorcluster=Stromal, overrides if reached (user said
    # majorcluster=Stromal cell). This includes fibroblasts/endothelial/smooth
    # under majorcluster=Stromal, so we apply it BEFORE the specific finalcluster
    # rules above to mirror user's explicit hierarchy.
    #
    # BUT: the user's explicit mapping places stromal cell AFTER the
    # fibroblast/endothelial/smooth specific rules. Read literally, the
    # mapping is: if row's finalcluster matches a specific rule use that;
    # else if majorcluster=Stromal use stromal. Implement exactly that.
    remaining = out.isna() & (mc == "Stromal cell")
    out[remaining] = "stromal cell"
    return out


def main():
    t = {}
    t0 = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("D1 pilot — Qu-only macaque, 7 author-label-defined types")
    print("=" * 70)

    # ── 1. Load Qu metadata + harmonize labels ────────────────────────────
    print(f"\n[1] Loading Qu metadata CSV ({QU_META_CSV.name})")
    t_s = time.time()
    qu_meta = pd.read_csv(QU_META_CSV, low_memory=False)
    t["load_metadata"] = time.time() - t_s
    print(f"  Shape: {qu_meta.shape}  ({t['load_metadata']:.1f}s)")

    qu_meta["target_type"] = harmonize_qu_author_labels(qu_meta)
    counts = qu_meta["target_type"].value_counts(dropna=False).to_dict()
    print("\n[1a] Qu author-label counts per target type:")
    low_power = []
    for typ in QU_TARGET_TYPES:
        n = counts.get(typ, 0)
        flag = " [LOW-POWER: <50]" if n < 50 else ""
        print(f"  {typ:<30s} {n:>7,}{flag}")
        if n < 50:
            low_power.append(typ)
    print(f"  (unmapped: {counts.get(np.nan, 0):,})")

    included_types = [t for t in QU_TARGET_TYPES if counts.get(t, 0) >= 50]
    if low_power:
        print(f"\n  Excluding {len(low_power)} low-power type(s): {low_power}")
    print(f"  Pilot n_types: {len(included_types)}")

    # ── 2. Parse metadata cell_id → (sample, barcode) for barcode matching
    # Format is "{Organ}@@_{BARCODE-1}" but the first component is the
    # sample tag not the Organ (e.g., "Liver-Normal@@_..."). Verify by
    # checking against 'sample' column.
    print("\n[2] Parsing cell_id format")
    sample_vals = qu_meta["sample"].astype(str).to_numpy()
    ids = qu_meta["cell_id"].astype(str).to_numpy()
    # Split on "@@_"
    sep = "@@_"
    split_ok = np.array([sep in x for x in ids])
    print(f"  cell_ids with '@@_' separator: {split_ok.sum():,} / {len(ids):,}")
    prefixes = np.array([x.split(sep, 1)[0] if sep in x else x for x in ids])
    bcs = np.array([x.split(sep, 1)[1] if sep in x else x for x in ids])
    unique_prefix = pd.Series(prefixes).value_counts()
    print(f"  Unique prefixes: {len(unique_prefix)}")
    # Compare prefix vs 'sample' and 'orig.ident'
    match_sample = (prefixes == sample_vals).mean()
    match_orig = (prefixes == qu_meta["orig.ident"].astype(str).to_numpy()).mean()
    print(f"  prefix matches metadata 'sample': {match_sample:.1%}")
    print(f"  prefix matches metadata 'orig.ident': {match_orig:.1%}")

    # ── 3. Load aggregate h5ad and align via (tissue, barcode)
    print(f"\n[3] Loading aggregate h5ad")
    t_s = time.time()
    adata = ad.read_h5ad(AGG_H5AD)
    t["load_h5ad"] = time.time() - t_s
    print(f"  Shape: {adata.shape}  ({t['load_h5ad']:.1f}s)")

    # Qu metadata uses its own sample/orig.ident naming (e.g., "Liver",
    # "Liver-Normal") that doesn't match the h5ad tissue tags derived from
    # GSM filenames ("Liver1", "Liver2"). We build a crosswalk from sample
    # → h5ad tissue by maximum barcode overlap per metadata.sample value.
    print("\n[3a] Building meta.sample → h5ad.tissue crosswalk by barcode overlap")
    h5ad_barcodes = adata.obs["barcode"].astype(str).to_numpy()
    h5ad_tissues = adata.obs["tissue"].astype(str).to_numpy()
    # Index h5ad by (barcode, tissue)
    h5ad_key_to_row = {f"{bc}__{tt}": i for i, (bc, tt) in
                       enumerate(zip(h5ad_barcodes, h5ad_tissues))}
    unique_tissues = sorted(set(h5ad_tissues))

    sample_to_tissue: dict[str, str] = {}
    for s in pd.unique(sample_vals):
        meta_bc_for_sample = set(bcs[sample_vals == s])
        # For each candidate tissue, count overlap
        best, best_ov = None, 0
        for tt in unique_tissues:
            tt_bcs = set(h5ad_barcodes[h5ad_tissues == tt])
            ov = len(meta_bc_for_sample & tt_bcs)
            if ov > best_ov:
                best_ov = ov
                best = tt
        sample_to_tissue[s] = best
        print(f"  meta.sample={s:<20s} → h5ad.tissue={best}  (overlap {best_ov:,})")

    # Now build metadata keys using the crosswalked tissue
    remapped_tissue = np.array([sample_to_tissue.get(s, "UNKNOWN") for s in sample_vals])
    meta_keys = np.array([f"{b}__{t}" for b, t in zip(bcs, remapped_tissue)])
    matched_idx: list[int] = []
    matched_meta_rows: list[int] = []
    for i, k in enumerate(meta_keys):
        if k in h5ad_key_to_row:
            matched_idx.append(h5ad_key_to_row[k])
            matched_meta_rows.append(i)
    matched_idx = np.array(matched_idx, dtype=np.int64)
    matched_meta_rows = np.array(matched_meta_rows, dtype=np.int64)
    print(f"  Matched {len(matched_idx):,} / {len(meta_keys):,} metadata cells to h5ad")

    # Take only cells that matched AND are in included_types
    target_types_assigned = qu_meta["target_type"].to_numpy()[matched_meta_rows]
    included_set = set(included_types)
    type_mask = np.array([
        (t is not None) and (not (isinstance(t, float) and np.isnan(t))) and t in included_set
        for t in target_types_assigned
    ])
    final_idx = matched_idx[type_mask]
    final_types = target_types_assigned[type_mask].astype(str)
    print(f"  After type filter: {len(final_idx):,} cells")

    adata = adata[final_idx].copy()
    adata.obs["target_type"] = final_types

    # ── 4. Subsample per type
    print(f"\n[4] Subsampling per type to max {MAX_CELLS_PER_TYPE:,} (seed {RANDOM_SEED})")
    rng = np.random.default_rng(RANDOM_SEED)
    sub_mask = np.zeros(adata.n_obs, dtype=bool)
    per_type_available = {}
    per_type_used = {}
    for typ in included_types:
        idx = np.where((adata.obs["target_type"] == typ).to_numpy())[0]
        per_type_available[typ] = int(len(idx))
        sel = rng.choice(idx, size=min(MAX_CELLS_PER_TYPE, len(idx)), replace=False) if len(idx) > MAX_CELLS_PER_TYPE else idx
        sub_mask[sel] = True
        per_type_used[typ] = int(sub_mask[sel].sum())
        print(f"  {typ:<30s} {len(idx):>7,} → {int(sub_mask[sel].sum()):>5,}")
    adata = adata[sub_mask].copy()

    # ── 5. Normalize + centroids
    print(f"\n[5] normalize_total(1e4) + log1p")
    t_s = time.time()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    t["normalize"] = time.time() - t_s

    print(f"\n[6] Compute per-type centroids")
    mac_cents = np.zeros((len(included_types), adata.n_vars), dtype=np.float64)
    for i, typ in enumerate(included_types):
        m = (adata.obs["target_type"] == typ).to_numpy()
        mac_cents[i] = np.asarray(adata.X[m].mean(axis=0)).flatten()
    mac_df = pd.DataFrame(mac_cents, index=included_types, columns=adata.var_names)
    mac_df.to_csv(OUT_CENTROIDS)

    # ── 7. Human subset + PCA + Procrustes
    print(f"\n[7] Subset human centroids + joint PCA + Procrustes")
    human = pd.read_csv(HUMAN_CENTROIDS_CSV, index_col=0)
    gene_ids = list(adata.var_names)
    human_sub = human.loc[included_types, gene_ids].copy()
    hv, mv = human_sub.to_numpy(), mac_df.to_numpy()
    print(f"  Human: min={hv.min():.3f} max={hv.max():.3f} std={hv.std():.3f}")
    print(f"  Macaque: min={mv.min():.3f} max={mv.max():.3f} std={mv.std():.3f}")
    print(f"  std ratio (mac/human): {mv.std() / max(hv.std(), 1e-12):.3f}")

    Xp, Yp, pca, order = pca_reduce_centroids(human_sub, mac_df)
    result = procrustes_align(Xp, Yp)
    p, null = permutation_test(Xp, Yp, n_permutations=N_PERMUTATIONS, seed=RANDOM_SEED)
    obs_null = result.distance / np.median(null)
    obs_null_mean = result.distance / np.mean(null)

    # Per-type residuals
    residuals = {}
    for i, typ in enumerate(order):
        r = result.aligned_target[i] - result.centered_reference[i]
        mag = float(np.linalg.norm(r))
        pct = 100.0 * mag**2 / result.distance_squared if result.distance_squared > 0 else 0.0
        residuals[typ] = {"magnitude": mag, "pct_ssr": pct}
    # Rank by magnitude
    ranked = sorted(residuals.items(), key=lambda kv: -kv[1]["magnitude"])

    out = {
        "analysis": "QU7_D1_PILOT",
        "gene_space_used": adata.n_vars,
        "gene_space_rationale": "13,927 three-way ortholog intersection (consistent with prior analyses)",
        "n_types_requested": 7,
        "n_types_analyzed": len(included_types),
        "included_types": included_types,
        "excluded_low_power": low_power,
        "per_type_available_cells": per_type_available,
        "per_type_used_cells": per_type_used,
        "label_mapping": {
            "bladder urothelial cell": "finalcluster=Epithelial cell AND Organ=Bladder",
            "hepatocyte":               "finalcluster=Epithelial cell AND Organ=Liver",
            "endothelial cell":          "finalcluster=Endothelial cell",
            "epithelial cell":           "finalcluster=Epithelial cell AND Organ NOT IN {Bladder, Liver}",
            "fibroblast":                "finalcluster=Fibroblasts",
            "smooth muscle cell":        "finalcluster=Smooth muscle cell",
            "stromal cell":              "majorcluster=Stromal cell (AFTER specific finalcluster rules)",
        },
        "procrustes": {
            "distance": float(result.distance),
            "distance_squared": float(result.distance_squared),
            "scaling": float(result.scaling),
            "rotation_det": float(np.linalg.det(result.rotation)),
        },
        "permutation_test": {
            "p_value": float(p),
            "obs_null_ratio_median": float(obs_null),
            "obs_null_ratio_mean": float(obs_null_mean),
            "null_median": float(np.median(null)),
            "null_mean": float(np.mean(null)),
            "n_permutations": N_PERMUTATIONS,
        },
        "pca": {
            "n_components": int(pca.n_components_),
            "variance_explained": pca.explained_variance_ratio_.tolist(),
            "cumulative_variance": float(np.cumsum(pca.explained_variance_ratio_)[-1]),
        },
        "per_type_residuals_ranked": [
            {"type": t, "magnitude": residuals[t]["magnitude"], "pct_ssr": residuals[t]["pct_ssr"]}
            for t, _ in ranked
        ],
        "committed_context": {
            "primary_obs_null_n20_asymmetric": 0.841,
            "sensitivity_obs_null_n13_RIRA_only": 0.7488,
            "primary_human_mouse_obs_null_n35": 0.522,
        },
        "seed": RANDOM_SEED,
        "timings_s": t,
        "runtime_total_s": time.time() - t0,
        "preprocessing": "raw UMI → normalize_total(1e4) → log1p (symmetric with human/mouse)",
    }

    def _default(o):
        if isinstance(o, np.bool_): return bool(o)
        if isinstance(o, np.integer): return int(o)
        if isinstance(o, np.floating): return float(o)
        if isinstance(o, np.ndarray): return o.tolist()
        raise TypeError(f"not serialisable: {type(o).__name__}")
    OUT_JSON.write_text(json.dumps(out, indent=2, default=_default))

    dt = time.time() - t0
    print("\n" + "=" * 70)
    print(f"D1 PILOT RESULTS (runtime {dt:.1f}s)")
    print("=" * 70)
    print(f"  n_types:          {len(included_types)}")
    print(f"  gene_space:       {adata.n_vars} (three-way ortholog)")
    print(f"  distance:         {result.distance:.4f}")
    print(f"  scaling:          {result.scaling:.6f}")
    print(f"  obs/null (med):   {obs_null:.4f}")
    print(f"  obs/null (mean):  {obs_null_mean:.4f}")
    print(f"  p-value:          {p:.6f}")
    print(f"  n_permutations:   {N_PERMUTATIONS}")
    print(f"  PCA components:   {pca.n_components_}")
    print()
    print(f"  Per-type residuals (high → low):")
    for t_, d_ in ranked:
        print(f"    {t_:<30s}  mag={d_['magnitude']:.3f}  pct_ssr={d_['pct_ssr']:.1f}%")
    print(f"\n  Results JSON: {OUT_JSON}")


if __name__ == "__main__":
    main()
