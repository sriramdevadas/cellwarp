#!/usr/bin/env python3
"""
CellWarp — Tabula Sapiens 2.0 Verification

VERIFICATION 1: Query Census for all Tabula Sapiens donor_ids.
Report which TSP donors exist, cells per donor, tissues per donor.

VERIFICATION 2: Check primary CellWarp dataset for donor overlap.

No download. No pipeline. Verification only.
"""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import cellxgene_census
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
CENSUS_VERSION = "2025-11-08"
ORGANISM = "Homo sapiens"

# Primary CellWarp human data paths (check multiple possible locations)
PRIMARY_PATHS = [
    PROJECT_ROOT / "data" / "phase1" / "human_raw.h5ad",
    PROJECT_ROOT / "data" / "phase1" / "human_qc.h5ad",
    PROJECT_ROOT / "data" / "phase1" / "human_aligned.h5ad",
]


def main():
    print("=" * 80)
    print("CellWarp — Tabula Sapiens 2.0 Donor Verification")
    print("=" * 80)

    # ==================================================================
    # VERIFICATION 1: Census Tabula Sapiens donors
    # ==================================================================
    print("\n[VERIFICATION 1] Querying Census for Tabula Sapiens donors...")

    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        # Find all Tabula Sapiens datasets
        datasets_df = census["census_info"]["datasets"].read().concat().to_pandas()
        ts_mask = datasets_df["collection_name"].str.contains(
            "Tabula Sapiens", case=False, na=False
        )
        ts_datasets = datasets_df[ts_mask]
        ts_dataset_ids = set(ts_datasets["dataset_id"].tolist())
        ts_collection_ids = ts_datasets["collection_id"].unique()

        print(f"\n  Tabula Sapiens collections: {len(ts_collection_ids)}")
        for cid in ts_collection_ids:
            cname = ts_datasets[ts_datasets["collection_id"] == cid]["collection_name"].iloc[0]
            n_ds = len(ts_datasets[ts_datasets["collection_id"] == cid])
            print(f"    {cid}: {cname} ({n_ds} datasets)")

        print(f"\n  Total Tabula Sapiens datasets: {len(ts_dataset_ids)}")

        # Query obs for all TS datasets
        print("\n  Querying all Tabula Sapiens cells in Census...")
        all_obs_frames = []
        for ds_id in sorted(ts_dataset_ids):
            try:
                obs = cellxgene_census.get_obs(
                    census, ORGANISM,
                    value_filter=f"dataset_id == '{ds_id}' and is_primary_data == True",
                    column_names=[
                        "soma_joinid", "donor_id", "cell_type", "tissue",
                        "assay", "disease", "development_stage", "dataset_id",
                    ],
                )
                if len(obs) > 0:
                    all_obs_frames.append(obs)
            except Exception as e:
                print(f"    Error querying {ds_id}: {e}")

        if not all_obs_frames:
            print("  ERROR: No Tabula Sapiens cells found in Census!")
            return

        ts_obs = pd.concat(all_obs_frames, ignore_index=True)
        print(f"\n  Total Tabula Sapiens cells in Census: {len(ts_obs):,}")

        # Donor analysis
        print(f"\n  ALL TABULA SAPIENS DONOR IDs IN CENSUS:")
        donor_summary = (
            ts_obs.groupby("donor_id")
            .agg(
                n_cells=("soma_joinid", "count"),
                tissues=("tissue", lambda x: ", ".join(sorted(x.dropna().unique())[:5])),
                n_tissues=("tissue", "nunique"),
                assays=("assay", lambda x: ", ".join(sorted(x.unique()))),
                diseases=("disease", lambda x: ", ".join(sorted(x.unique()))),
            )
            .sort_index()
        )

        print(f"\n  {'Donor ID':<20} {'Cells':>8} {'Tissues':>4} {'Assays':<30} {'Sample tissues'}")
        print(f"  {'-'*100}")
        for donor_id, row in donor_summary.iterrows():
            print(
                f"  {str(donor_id):<20} {row['n_cells']:>8,} {row['n_tissues']:>4} "
                f"{row['assays'][:30]:<30} {row['tissues'][:40]}"
            )

        all_census_donors = set(donor_summary.index.tolist())

        # Check specifically for TSP17-TSP30
        print(f"\n\n  TSP17-TSP30 STATUS IN CENSUS:")
        tsp_17_30 = set()
        for i in range(17, 31):
            donor_id = f"TSP{i}"
            tsp_17_30.add(donor_id)
            if donor_id in all_census_donors:
                row = donor_summary.loc[donor_id]
                print(f"    {donor_id}: PRESENT — {row['n_cells']:,} cells, {row['n_tissues']} tissues")
            else:
                # Check partial matches
                partial = [d for d in all_census_donors if str(d).startswith(f"TSP{i}")]
                if partial:
                    print(f"    {donor_id}: NOT EXACT — but partial matches: {partial}")
                else:
                    print(f"    {donor_id}: ABSENT")

    # ==================================================================
    # VERIFICATION 2: Primary CellWarp dataset donors
    # ==================================================================
    print("\n\n[VERIFICATION 2] Checking primary CellWarp dataset...")

    primary_donors = set()
    primary_path_used = None

    for path in PRIMARY_PATHS:
        if path.exists():
            print(f"\n  Loading: {path}")
            adata = ad.read_h5ad(path, backed="r")
            print(f"  Shape: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

            # Check for donor_id column
            obs_cols = list(adata.obs.columns)
            print(f"  Obs columns: {obs_cols[:15]}")

            donor_col = None
            for candidate in ["donor_id", "donor", "sample_id", "individual"]:
                if candidate in obs_cols:
                    donor_col = candidate
                    break

            if donor_col:
                donors = sorted(adata.obs[donor_col].unique().tolist())
                primary_donors = set(donors)
                primary_path_used = path
                print(f"\n  Donor column: '{donor_col}'")
                print(f"  Unique donors: {len(donors)}")
                print(f"  Donor IDs: {donors}")

                # Per-donor breakdown
                print(f"\n  {'Donor ID':<20} {'Cells':>8}")
                print(f"  {'-'*30}")
                for d in donors:
                    n = (adata.obs[donor_col] == d).sum()
                    print(f"  {str(d):<20} {n:>8,}")
            else:
                print(f"  WARNING: No donor_id column found in {obs_cols}")

            adata.file.close()
            break  # Use first existing file
    else:
        print("  WARNING: No primary human data file found!")
        # Try the scaled 35-type data
        scaled_path = PROJECT_ROOT / "data" / "phase1" / "human_aligned.h5ad"
        if not scaled_path.exists():
            print(f"  Checked: {[str(p) for p in PRIMARY_PATHS]}")

    # ==================================================================
    # OVERLAP ANALYSIS
    # ==================================================================
    print("\n\n[OVERLAP ANALYSIS]")
    print("=" * 80)

    print(f"\n  1. TSP donor IDs in Census: {len(all_census_donors)}")
    for d in sorted(all_census_donors):
        print(f"     {d}")

    print(f"\n  2. TSP donor IDs in primary CellWarp: {len(primary_donors)}")
    for d in sorted(primary_donors):
        print(f"     {d}")

    overlap = all_census_donors & primary_donors
    print(f"\n  3. Overlap: {len(overlap)}")
    for d in sorted(overlap):
        print(f"     {d}")

    # TSP17-30 that are in Census but NOT in primary
    tsp_17_30_in_census = tsp_17_30 & all_census_donors
    tsp_17_30_independent = tsp_17_30_in_census - primary_donors
    print(f"\n  4. TSP17-TSP30 donors in Census but NOT in primary: {len(tsp_17_30_independent)}")
    for d in sorted(tsp_17_30_independent):
        row = donor_summary.loc[d]
        print(f"     {d}: {row['n_cells']:,} cells, {row['n_tissues']} tissues")

    # Also check: any Census donors NOT in primary (regardless of TSP17-30 range)
    all_independent = all_census_donors - primary_donors
    print(f"\n  BONUS: ALL Census TS donors NOT in primary: {len(all_independent)}")
    for d in sorted(all_independent):
        if d in donor_summary.index:
            row = donor_summary.loc[d]
            print(f"     {d}: {row['n_cells']:,} cells, {row['n_tissues']} tissues")


if __name__ == "__main__":
    main()
