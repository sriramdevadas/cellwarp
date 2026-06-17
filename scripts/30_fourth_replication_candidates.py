#!/usr/bin/env python3
"""
CellWarp — Fourth Replication Dataset Candidate Search

Queries CZ CELLxGENE Census for candidate human datasets that can serve as
a fourth independent replication of the core Procrustes result.

Biology
-------
We need a human atlas independent of Tabula Sapiens (primary), Sun2023
(mouse replication), and PanSci (mouse replication). The fourth dataset
pairs with Tabula Muris Senis mouse to form a new human-mouse Procrustes
test. Must cover the original 6 cell types (hepatocyte, CD8+ T, endothelial,
CD4+ T, B cell, macrophage) with ≥500 cells each.

Criteria
--------
  - Organism: Homo sapiens
  - is_primary_data: True
  - disease: "normal"
  - development_stage: adult (>=18 years)
  - assay: 10x Chromium (3' or 5')
  - NOT Tabula Sapiens collection
  - ≥500 cells per type for the 6 original types
  - Prefer single collection (stronger independence claim than pooled)

Output
------
  Prints ranked candidate collections with cell type coverage details.

Usage:
    python scripts/30_fourth_replication_candidates.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import cellxgene_census
import pandas as pd

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

CENSUS_VERSION = "2025-11-08"
ORGANISM = "Homo sapiens"

# The 6 original cell types (Census cell_type labels)
TARGET_TYPES = [
    "hepatocyte",
    "CD8-positive, alpha-beta T cell",
    "endothelial cell",
    "CD4-positive, alpha-beta T cell",
    "B cell",
    "macrophage",
]

SHORT_NAMES = {
    "hepatocyte": "Hepatocyte",
    "CD8-positive, alpha-beta T cell": "CD8+ T",
    "endothelial cell": "Endothelial",
    "CD4-positive, alpha-beta T cell": "CD4+ T",
    "B cell": "B cell",
    "macrophage": "Macrophage",
}

MIN_CELLS = 500

CHROMIUM_ASSAYS = {
    "10x 3' v1",
    "10x 3' v2",
    "10x 3' v3",
    "10x 3' transcription profiling",
    "10x 5' v1",
    "10x 5' v2",
    "10x 5' transcription profiling",
}


def main():
    print("=" * 80)
    print("CellWarp — Fourth Replication Dataset Candidate Search")
    print("=" * 80)

    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        # Step 1: Get dataset/collection metadata
        print("\n[Step 1] Loading Census dataset metadata...")
        datasets_df = census["census_info"]["datasets"].read().concat().to_pandas()
        print(f"  Total datasets in Census: {len(datasets_df):,}")

        # Identify Tabula Sapiens datasets
        ts_mask = datasets_df["collection_name"].str.contains(
            "Tabula Sapiens", case=False, na=False
        )
        ts_dataset_ids = set(datasets_df.loc[ts_mask, "dataset_id"].tolist())
        print(f"  Tabula Sapiens datasets to exclude: {len(ts_dataset_ids)}")

        # Build dataset_id → collection mapping
        ds_to_coll_id = dict(zip(datasets_df["dataset_id"], datasets_df["collection_id"]))
        ds_to_coll_name = dict(zip(datasets_df["dataset_id"], datasets_df["collection_name"]))

        # Step 2: Query for each target cell type separately
        print("\n[Step 2] Querying Census for each target cell type...")

        all_records = []

        for ct in TARGET_TYPES:
            short = SHORT_NAMES[ct]
            print(f"\n  Querying: {short} ({ct})...")

            obs_df = cellxgene_census.get_obs(
                census,
                ORGANISM,
                value_filter=(
                    f"is_primary_data == True "
                    f"and disease == 'normal' "
                    f"and cell_type == '{ct}'"
                ),
                column_names=[
                    "soma_joinid",
                    "cell_type",
                    "dataset_id",
                    "assay",
                    "development_stage",
                    "tissue",
                    "donor_id",
                ],
            )
            print(f"    Raw: {len(obs_df):,} cells")

            # Adult filter
            adult_mask = obs_df["development_stage"].str.contains(
                "adult|year-old|years old", case=False, na=False
            )
            obs_df = obs_df[adult_mask]
            print(f"    Adult: {len(obs_df):,} cells")

            # Exclude Tabula Sapiens
            obs_df = obs_df[~obs_df["dataset_id"].isin(ts_dataset_ids)]
            print(f"    Non-TS: {len(obs_df):,} cells")

            # 10x Chromium filter
            obs_df = obs_df[obs_df["assay"].isin(CHROMIUM_ASSAYS)]
            print(f"    10x Chromium: {len(obs_df):,} cells")

            if len(obs_df) == 0:
                continue

            # Map to collections
            obs_df = obs_df.copy()
            obs_df["collection_id"] = obs_df["dataset_id"].map(ds_to_coll_id)
            obs_df["collection_name"] = obs_df["dataset_id"].map(ds_to_coll_name)

            # Aggregate per collection
            for coll_id, grp in obs_df.groupby("collection_id"):
                all_records.append({
                    "collection_id": coll_id,
                    "collection_name": grp["collection_name"].iloc[0],
                    "cell_type": ct,
                    "short_name": short,
                    "n_cells": len(grp),
                    "n_donors": grp["donor_id"].nunique(),
                    "tissues": ", ".join(sorted(grp["tissue"].dropna().unique())),
                    "assays": ", ".join(sorted(grp["assay"].unique())),
                    "n_datasets": grp["dataset_id"].nunique(),
                })

        records_df = pd.DataFrame(all_records)

        # Step 3: Score collections
        print("\n\n[Step 3] Scoring collections...")

        # For each collection, count how many types pass ≥500
        collection_scores = []
        for coll_id, grp in records_df.groupby("collection_id"):
            coll_name = grp["collection_name"].iloc[0]
            type_info = {}
            types_passing = 0
            total_cells = 0

            for _, row in grp.iterrows():
                ct = row["cell_type"]
                n = row["n_cells"]
                passing = n >= MIN_CELLS
                type_info[ct] = {
                    "n_cells": n,
                    "n_donors": row["n_donors"],
                    "tissues": row["tissues"],
                    "assays": row["assays"],
                    "pass": passing,
                }
                total_cells += n
                if passing:
                    types_passing += 1

            collection_scores.append({
                "collection_id": coll_id,
                "collection_name": coll_name,
                "types_passing": types_passing,
                "types_total": len(grp),
                "type_info": type_info,
                "total_cells": total_cells,
            })

        # Sort: most types passing first, then by total cells
        collection_scores.sort(key=lambda x: (x["types_passing"], x["total_cells"]), reverse=True)

        # Step 4: Report
        print("\n" + "=" * 80)
        full_pass = [c for c in collection_scores if c["types_passing"] == 6]
        print(f"Collections with ALL 6 types ≥{MIN_CELLS}: {len(full_pass)}")
        print("=" * 80)

        if full_pass:
            for i, cand in enumerate(full_pass[:10]):
                _print_candidate(i + 1, cand, "FULL PASS")
        else:
            print("\nNo single collection covers all 6 types at ≥500 cells.")
            print("Showing top partial-coverage collections:\n")

            shown = 0
            for cand in collection_scores:
                if cand["types_passing"] >= 3 and shown < 15:
                    _print_candidate(shown + 1, cand, "PARTIAL")
                    shown += 1

        # Step 5: Hepatocyte bottleneck analysis
        print("\n" + "=" * 80)
        print("HEPATOCYTE BOTTLENECK ANALYSIS")
        print("=" * 80)
        hep_records = records_df[records_df["cell_type"] == "hepatocyte"].sort_values(
            "n_cells", ascending=False
        )
        if len(hep_records) == 0:
            print("  No hepatocyte records found in any collection!")
        else:
            print(f"\n  Collections with ANY hepatocytes (10x, adult, normal):")
            for _, row in hep_records.iterrows():
                status = "PASS" if row["n_cells"] >= MIN_CELLS else "FAIL"
                # Check what other target types this collection has
                coll_types = records_df[records_df["collection_id"] == row["collection_id"]]
                other_types = [
                    SHORT_NAMES[r["cell_type"]]
                    for _, r in coll_types.iterrows()
                    if r["cell_type"] != "hepatocyte" and r["n_cells"] >= MIN_CELLS
                ]
                print(
                    f"  [{status}] {row['n_cells']:>6,} cells | {row['n_donors']:>3} donors | "
                    f"{row['collection_name'][:50]}"
                )
                print(f"         Collection ID: {row['collection_id']}")
                print(f"         Tissues: {row['tissues']}")
                print(f"         Assays: {row['assays']}")
                if other_types:
                    print(f"         Other target types passing: {', '.join(other_types)}")
                else:
                    print(f"         Other target types passing: NONE")
                print()

        # Step 6: Multi-collection strategy
        print("\n" + "=" * 80)
        print("MULTI-COLLECTION STRATEGY (if single collection impossible)")
        print("=" * 80)

        # Find the smallest set of collections that covers all 6 types
        # Start with hepatocyte collections (bottleneck) and see what else they need
        if len(hep_records) > 0:
            for _, hep_row in hep_records.iterrows():
                if hep_row["n_cells"] < MIN_CELLS:
                    continue
                coll_id = hep_row["collection_id"]
                coll_name = hep_row["collection_name"]
                coll_types = records_df[records_df["collection_id"] == coll_id]
                passing_types = set(
                    coll_types[coll_types["n_cells"] >= MIN_CELLS]["cell_type"].tolist()
                )
                missing = set(TARGET_TYPES) - passing_types
                print(f"\n  Starting from: {coll_name[:60]}")
                print(f"  Has: {', '.join(SHORT_NAMES[t] for t in passing_types)}")
                if missing:
                    print(f"  Missing: {', '.join(SHORT_NAMES[t] for t in missing)}")
                else:
                    print(f"  COVERS ALL 6 — single collection sufficient!")

        # Summary
        print("\n" + "=" * 80)
        print("SUMMARY")
        print("=" * 80)
        print(f"\nPer-type totals (all non-TS collections, 10x Chromium, adult, normal):")
        for ct in TARGET_TYPES:
            ct_recs = records_df[records_df["cell_type"] == ct]
            total = ct_recs["n_cells"].sum()
            n_colls = len(ct_recs)
            print(f"  {SHORT_NAMES[ct]:<15} {total:>10,} cells across {n_colls:>3} collections")


def _print_candidate(rank, cand, label):
    """Print details for one candidate collection."""
    print(f"\n--- {label} Candidate {rank}: {cand['collection_name'][:70]}")
    print(f"    Collection ID: {cand['collection_id']}")
    print(f"    Types passing: {cand['types_passing']}/{len(TARGET_TYPES)}")
    print(f"    Total target cells: {cand['total_cells']:,}")
    print()
    print(f"    {'Type':<15} {'Cells':>8} {'Donors':>7} {'Tissues':<35} {'Status':>8}")
    print(f"    {'-'*75}")
    for ct in TARGET_TYPES:
        short = SHORT_NAMES[ct]
        if ct in cand["type_info"]:
            info = cand["type_info"][ct]
            status = "PASS" if info["pass"] else "FAIL"
            tissues_str = info["tissues"][:35] if isinstance(info["tissues"], str) else "—"
            print(
                f"    {short:<15} {info['n_cells']:>8,} {info['n_donors']:>7} "
                f"{tissues_str:<35} {status:>8}"
            )
        else:
            print(f"    {short:<15} {'0':>8} {'0':>7} {'—':<35} {'ABSENT':>8}")


if __name__ == "__main__":
    main()
