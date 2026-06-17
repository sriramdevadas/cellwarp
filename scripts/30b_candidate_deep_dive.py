#!/usr/bin/env python3
"""
CellWarp — Deep dive on top 3 hepatocyte-containing collections.

Checks all cell type labels (especially endothelial variants),
donor demographics, and publication details.

Usage:
    python scripts/30b_candidate_deep_dive.py
"""

from __future__ import annotations

import cellxgene_census
import pandas as pd

CENSUS_VERSION = "2025-11-08"
ORGANISM = "Homo sapiens"

CANDIDATE_COLLECTIONS = {
    "0c8a364b-97b5-4cc8-a593-23c38c6f0ac5": "Liver immunology (5/6 types)",
    "74e10dc4-cbb2-4605-a189-8a1cd8e44d8c": "Spatial proteogenomics liver (4/6 types)",
    "ff69f0ee-fef6-4895-9f48-6c64a68c8289": "Pediatric liver (check age)",
}

CHROMIUM_ASSAYS = {
    "10x 3' v1", "10x 3' v2", "10x 3' v3",
    "10x 3' transcription profiling",
    "10x 5' v1", "10x 5' v2",
    "10x 5' transcription profiling",
}

ENDOTHELIAL_KEYWORDS = ["endothelial", "sinusoidal", "lsec", "vascular"]


def main():
    print("=" * 80)
    print("CellWarp — Candidate Deep Dive (Hepatocyte Collections)")
    print("=" * 80)

    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        # Get dataset metadata
        datasets_df = census["census_info"]["datasets"].read().concat().to_pandas()

        for coll_id, desc in CANDIDATE_COLLECTIONS.items():
            print(f"\n{'='*80}")
            print(f"COLLECTION: {desc}")
            print(f"ID: {coll_id}")
            print(f"{'='*80}")

            coll_datasets = datasets_df[datasets_df["collection_id"] == coll_id]
            if len(coll_datasets) == 0:
                print("  NOT FOUND!")
                continue

            coll_name = coll_datasets["collection_name"].iloc[0]
            print(f"  Full name: {coll_name}")
            print(f"  N datasets: {len(coll_datasets)}")
            for _, ds in coll_datasets.iterrows():
                print(f"    Dataset: {ds.get('dataset_title', 'N/A')}")
                print(f"    ID: {ds['dataset_id']}")

            # Query per dataset_id
            all_frames = []
            for _, ds in coll_datasets.iterrows():
                ds_id = ds["dataset_id"]
                try:
                    obs = cellxgene_census.get_obs(
                        census, ORGANISM,
                        value_filter=(
                            f"dataset_id == '{ds_id}' "
                            f"and is_primary_data == True "
                            f"and disease == 'normal'"
                        ),
                        column_names=[
                            "soma_joinid", "cell_type", "dataset_id", "assay",
                            "development_stage", "tissue", "donor_id", "sex",
                        ],
                    )
                    all_frames.append(obs)
                except Exception as e:
                    print(f"    Error querying {ds_id}: {e}")

            if not all_frames:
                print("  No data retrieved!")
                continue

            coll_obs = pd.concat(all_frames, ignore_index=True)
            print(f"\n  Total cells (normal, primary): {len(coll_obs):,}")

            # Adult filter
            adult_mask = coll_obs["development_stage"].str.contains(
                "adult|year-old|years old", case=False, na=False
            )
            n_adult = adult_mask.sum()
            n_nonadult = len(coll_obs) - n_adult
            print(f"  Adult: {n_adult:,}  |  Non-adult: {n_nonadult:,}")

            # Show development stages
            print(f"\n  DEVELOPMENT STAGES (all cells):")
            for stage, n in coll_obs["development_stage"].value_counts().head(10).items():
                print(f"    {stage}: {n:,}")

            # Work with adult + 10x subset
            coll_10x = coll_obs[adult_mask & coll_obs["assay"].isin(CHROMIUM_ASSAYS)]
            print(f"\n  Adult 10x Chromium cells: {len(coll_10x):,}")

            if len(coll_10x) == 0:
                # Also show non-adult 10x if no adult 10x
                coll_10x_all = coll_obs[coll_obs["assay"].isin(CHROMIUM_ASSAYS)]
                print(f"  (All-age 10x Chromium: {len(coll_10x_all):,})")
                coll_10x = coll_10x_all  # Show all-age for pediatric

            # All cell types
            print(f"\n  ALL CELL TYPES:")
            type_counts = (
                coll_10x.groupby("cell_type")
                .agg(n_cells=("soma_joinid", "count"), n_donors=("donor_id", "nunique"))
                .sort_values("n_cells", ascending=False)
            )
            for ct, row in type_counts.iterrows():
                marker = ""
                if any(kw in str(ct).lower() for kw in ENDOTHELIAL_KEYWORDS):
                    marker = " *** ENDOTHELIAL ***"
                gate = "PASS" if row["n_cells"] >= 500 else "fail"
                print(f"    [{gate:>4}] {ct:<55} {row['n_cells']:>6,} cells  {row['n_donors']:>3} donors{marker}")

            # Tissues
            print(f"\n  TISSUES:")
            for tissue in sorted(coll_10x["tissue"].dropna().unique()):
                print(f"    {tissue}: {len(coll_10x[coll_10x['tissue'] == tissue]):,}")

            # Assays
            print(f"\n  ASSAYS:")
            for assay, n in coll_10x["assay"].value_counts().items():
                print(f"    {assay}: {n:,}")

            # Donors
            print(f"\n  DONORS: {coll_10x['donor_id'].nunique()}")
            print(f"  SEX: {dict(coll_10x['sex'].value_counts())}")


if __name__ == "__main__":
    main()
