"""Identify top independent datasets per species that cover the most target cell types."""
import cellxgene_census
import pandas as pd
import numpy as np
from pathlib import Path

OUT = Path(__file__).parent
TARGET = sorted(pd.read_csv(OUT / "passing_cell_types.csv")["cell_type"].tolist())

EXCLUDED_NAMES = ["Tabula Sapiens", "Tabula Muris Senis", "CellHint", "PanSci", "EasySci"]

with cellxgene_census.open_soma(census_version="2025-11-08") as census:
    ds_df = census["census_info"]["datasets"].read().concat().to_pandas()

    excluded = set()
    for n in EXCLUDED_NAMES:
        excluded.update(ds_df.loc[ds_df["collection_name"].str.contains(n, case=False, na=False), "dataset_id"])

    for org, label in [("Mus musculus", "MOUSE"), ("Homo sapiens", "HUMAN")]:
        print(f"\n{'='*70}\n{label}\n{'='*70}")
        names = ", ".join(f"'{ct}'" for ct in TARGET)
        obs = cellxgene_census.get_obs(
            census, org,
            value_filter=f"cell_type in [{names}] and is_primary_data == True and disease == 'normal'",
            column_names=["cell_type", "dataset_id"],
        )
        obs = obs[~obs["dataset_id"].isin(excluded)]

        # Per-dataset stats
        ds_stats = []
        for ds_id, grp in obs.groupby("dataset_id", observed=True):
            ct_counts = grp["cell_type"].value_counts()
            passing = [ct for ct in ct_counts.index if ct in TARGET and ct_counts[ct] >= 200]
            ds_stats.append({
                "dataset_id": ds_id,
                "total_cells": len(grp),
                "n_types": len(passing),
                "passing_types": passing,
                "type_list": "; ".join(sorted(passing)),
            })
        ds_stats = pd.DataFrame(ds_stats).sort_values("n_types", ascending=False)

        # Show top datasets
        print(f"\nTop datasets by cell type coverage (>= 200 cells/type):")
        for _, r in ds_stats.head(15).iterrows():
            title = ds_df.loc[ds_df["dataset_id"] == r["dataset_id"], "dataset_title"].values
            title = title[0][:60] if len(title) > 0 else "?"
            print(f"  {r['dataset_id'][:12]}  {r['n_types']:>2} types  "
                  f"{r['total_cells']:>8,} cells  {title}")

        # Greedy set cover: pick datasets to maximize type coverage
        covered = set()
        selected = []
        remaining = ds_stats.copy()
        while len(covered) < len(TARGET) and len(remaining) > 0:
            # Score: number of NEW types covered
            remaining["new_types"] = remaining["passing_types"].apply(
                lambda x: len(set(x) - covered)
            )
            remaining = remaining.sort_values("new_types", ascending=False)
            best = remaining.iloc[0]
            if best["new_types"] == 0:
                break
            selected.append(best)
            covered.update(best["passing_types"])
            remaining = remaining.iloc[1:]

        print(f"\nGreedy set cover: {len(selected)} datasets cover {len(covered)}/{len(TARGET)} types")
        total_cells = 0
        for s in selected:
            title = ds_df.loc[ds_df["dataset_id"] == s["dataset_id"], "dataset_title"].values
            title = title[0][:50] if len(title) > 0 else "?"
            total_cells += s["total_cells"]
            print(f"  {s['dataset_id'][:12]}  +{len(set(s['passing_types']) - (covered - set(s['passing_types'])))} new  "
                  f"{s['total_cells']:>8,} cells  types: {s['type_list'][:80]}")
        print(f"  Total cells: {total_cells:,}")
        missing = set(TARGET) - covered
        if missing:
            print(f"  Missing types: {sorted(missing)}")

print("\nDone.")
