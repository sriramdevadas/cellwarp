#!/usr/bin/env python3
"""
Bilateral mt-fraction re-run via CZ CELLxGENE Census raw counts.

An earlier attempt used local mouse_scaled.h5ad, which has
mt-counts stripped (mouse-side mt genes were zeroed during the
ortholog-alignment pipeline). We re-pull TMS raw
counts for the 35 matched mouse cell types from Census 2025-11-08
(pre-ortholog filtering, full per-species gene space) to recover the
mouse mt-fraction signal.

Methodology:
  - Filter: is_primary_data == True, disease == 'normal',
    cell_type ∈ {35 matched names from Table S5}, Tabula Muris
    Senis dataset IDs (MOUSE_COLLECTION).
  - Cap at MAX_CELLS_PER_TYPE = 2,000 per type (RANDOM_SEED = 42)
    — same convention as primary 08_scaled_procrustes pipeline.
  - Per-cell mt-fraction = sum(mt-gene UMI counts) / sum(total UMI
    counts); mt-genes are those with feature_name starting with 'mt-'.
  - Per-cell-type mt-fraction = mean of per-cell mt-fraction over all
    cells of that type.
  - Combined per-type mt-fraction = (human + mouse) / 2.
  - Spearman correlation of combined per-type mt-fraction vs the
    35-type rigidity ranking from analysis/.../residuals_ranked.csv.

Outputs:
  results.json — full numerics
  per_type_table.csv — 35 × {human, mouse, combined, rigidity_rank}
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

PROJECT = Path(__file__).resolve().parents[3]
OUT_DIR = Path(__file__).resolve().parent
PHASE2_DIR = PROJECT / "output" / "phase2" / "scaled_35types"
A3_PATH = PROJECT / "analysis/sensitivity/mt_cellcycle/results.json"

sys.path.insert(0, str(PROJECT / "src"))
from cellwarp.data_loader import (  # noqa: E402
    HUMAN_COLLECTION,
    HUMAN_ORGANISM,
    MOUSE_COLLECTION,
    MOUSE_ORGANISM,
    MAX_CELLS_PER_TYPE,
    RANDOM_SEED,
    build_obs_value_filter,
    download_species_data,
    get_dataset_ids_for_collection,
    subsample_per_cell_type,
)

CENSUS_VERSION = "2025-11-08"
LOCAL_CACHE_MOUSE = OUT_DIR / "tms_35types_raw.h5ad"
LOCAL_CACHE_HUMAN = OUT_DIR / "ts_35types_raw.h5ad"


def mt_fraction_from_adata(adata: ad.AnnData, mt_prefix: str) -> pd.Series:
    """Per-cell mt-fraction over all genes available in adata (full per-species)."""
    feature_names = adata.var.get("feature_name", pd.Series(adata.var_names, index=adata.var_names))
    mt_mask = feature_names.astype(str).str.startswith(mt_prefix)
    n_mt = int(mt_mask.sum())
    print(f"  '{mt_prefix}'-prefix mt genes: {n_mt}")
    mt_indices = np.where(mt_mask.values)[0]
    X = adata.X
    if sp.issparse(X):
        total = np.asarray(X.sum(axis=1)).ravel()
        mt_X = X[:, mt_indices]
        mt_counts = np.asarray(mt_X.sum(axis=1)).ravel()
    else:
        total = X.sum(axis=1)
        mt_counts = X[:, mt_indices].sum(axis=1)
    fraction = np.where(total > 0, mt_counts / total, np.nan)
    return pd.Series(fraction, index=adata.obs.index)


def download_or_load(species: str, cell_types: list[str]) -> ad.AnnData:
    cache = LOCAL_CACHE_MOUSE if species == "mouse" else LOCAL_CACHE_HUMAN
    if cache.exists():
        print(f"  Loading cached {cache.name} ...")
        return ad.read_h5ad(cache)

    import cellxgene_census
    organism = MOUSE_ORGANISM if species == "mouse" else HUMAN_ORGANISM
    collection = MOUSE_COLLECTION if species == "mouse" else HUMAN_COLLECTION
    print(f"  Opening Census ({CENSUS_VERSION}) for {species} ...")
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        ds_ids = get_dataset_ids_for_collection(census, collection)
        obs_filter = build_obs_value_filter(cell_types, ds_ids)
        adata = download_species_data(census, organism, obs_filter)
    adata = subsample_per_cell_type(adata, MAX_CELLS_PER_TYPE,
                                     cell_type_column="cell_type",
                                     random_seed=RANDOM_SEED)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(cache)
    print(f"  Cached → {cache}")
    return adata


def main() -> None:
    t0 = time.time()
    print("=" * 70)
    print("Bilateral mt-fraction re-run via Census 2025-11-08")
    print("=" * 70)

    s5 = pd.read_csv(PROJECT / "docs/supplementary_materials/table_S5.csv")
    # Use the human_cell_type names — these match the Census cell_type strings
    # that the rest of the pipeline uses (Table S5 is the canonical 35-type list).
    cell_types = sorted(s5["human_cell_type"].tolist())
    cell_ontology_ids = sorted(s5["cell_ontology_id"].tolist())
    assert len(cell_types) == 35, f"Expected 35, got {len(cell_types)}"
    print(f"\n35 matched cell types from Table S5")

    # Mouse: re-pull raw from Census to get pre-ortholog full gene space
    print("\nMouse download ...")
    mouse = download_or_load("mouse", cell_types)
    print(f"  mouse: {mouse.n_obs:,} cells × {mouse.n_vars:,} genes "
          f"({mouse.obs['cell_type'].nunique()} types)")

    # Human: re-pull raw similarly (for symmetry; the A3 mt-fraction used the
    # restricted ortholog space, so this gives the cleaner bilateral comparison)
    print("\nHuman download ...")
    human = download_or_load("human", cell_types)
    print(f"  human: {human.n_obs:,} cells × {human.n_vars:,} genes "
          f"({human.obs['cell_type'].nunique()} types)")

    # Compute mt-fraction per cell
    print("\nmt-fraction (mouse, mt-* prefix on feature_name) ...")
    mt_mouse = mt_fraction_from_adata(mouse, "mt-")
    mouse.obs["mt_fraction"] = mt_mouse.values
    mouse_per_type = mouse.obs.groupby("cell_type", observed=True)["mt_fraction"].mean()
    print(f"  mouse per-type mt-fraction range: "
          f"{mouse_per_type.min()*100:.2f}% – {mouse_per_type.max()*100:.2f}%")

    print("\nmt-fraction (human, MT- prefix on feature_name) ...")
    mt_human = mt_fraction_from_adata(human, "MT-")
    human.obs["mt_fraction"] = mt_human.values
    human_per_type = human.obs.groupby("cell_type", observed=True)["mt_fraction"].mean()
    print(f"  human per-type mt-fraction range: "
          f"{human_per_type.min()*100:.2f}% – {human_per_type.max()*100:.2f}%")

    # Combined
    common_types = sorted(set(mouse_per_type.index) & set(human_per_type.index))
    print(f"\nCommon types: {len(common_types)} (expected 35)")

    per_type = pd.DataFrame({
        "cell_type": common_types,
        "human_mt_fraction": [float(human_per_type[t]) for t in common_types],
        "mouse_mt_fraction": [float(mouse_per_type[t]) for t in common_types],
    })
    per_type["combined_mt_fraction"] = (per_type["human_mt_fraction"]
                                         + per_type["mouse_mt_fraction"]) / 2

    # Spearman vs rigidity ranking
    residuals = pd.read_csv(PHASE2_DIR / "residuals_ranked.csv")
    rank_map = dict(zip(residuals["cell_type"], residuals["rank"]))
    per_type["rigidity_rank"] = per_type["cell_type"].map(rank_map)
    per_type = per_type.dropna(subset=["rigidity_rank"]).reset_index(drop=True)
    n = len(per_type)
    rho_c, p_c = stats.spearmanr(per_type["combined_mt_fraction"],
                                 per_type["rigidity_rank"])
    rho_h, p_h = stats.spearmanr(per_type["human_mt_fraction"],
                                 per_type["rigidity_rank"])
    rho_m, p_m = stats.spearmanr(per_type["mouse_mt_fraction"],
                                 per_type["rigidity_rank"])

    print(f"\nSpearman correlation vs 35-type rigidity rank (n={n}):")
    print(f"  Combined: ρ = {rho_c:+.4f}, p = {p_c:.4f}")
    print(f"  Human:    ρ = {rho_h:+.4f}, p = {p_h:.4f}")
    print(f"  Mouse:    ρ = {rho_m:+.4f}, p = {p_m:.4f}")

    per_type.to_csv(OUT_DIR / "per_type_table.csv", index=False)

    results = {
        "metadata": {
            "script": str(Path(__file__).resolve().relative_to(PROJECT)),
            "runtime_sec": time.time() - t0,
            "census_version": CENSUS_VERSION,
            "human_collection": HUMAN_COLLECTION,
            "mouse_collection": MOUSE_COLLECTION,
            "max_cells_per_type": MAX_CELLS_PER_TYPE,
            "random_seed": RANDOM_SEED,
            "methodology": (
                "Raw UMI counts pulled fresh from Census 2025-11-08, "
                "pre-ortholog filtering. mt-fraction = sum(mt-gene UMI "
                "counts) / sum(total UMI counts) per cell, mean per type. "
                "Mouse mt genes identified by feature_name ^mt-, human by "
                "^MT-. Combined = arithmetic mean (human + mouse) / 2."
            ),
        },
        "human_atlas": {
            "n_cells": int(human.n_obs),
            "n_genes_full": int(human.n_vars),
            "n_cell_types": int(human.obs["cell_type"].nunique()),
        },
        "mouse_atlas": {
            "n_cells": int(mouse.n_obs),
            "n_genes_full": int(mouse.n_vars),
            "n_cell_types": int(mouse.obs["cell_type"].nunique()),
        },
        "per_type_table": per_type.to_dict("records"),
        "spearman_vs_rigidity": {
            "combined": {"rho": float(rho_c), "p_value": float(p_c), "n": n},
            "human": {"rho": float(rho_h), "p_value": float(p_h), "n": n},
            "mouse": {"rho": float(rho_m), "p_value": float(p_m), "n": n},
        },
    }
    with open(OUT_DIR / "results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\n✓ Saved {OUT_DIR/'results.json'}")
    print(f"Runtime: {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
