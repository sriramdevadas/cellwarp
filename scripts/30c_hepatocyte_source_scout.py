#!/usr/bin/env python3
"""
CellWarp — Hepatocyte Source Scout

Searches Census for non-liver-only hepatocyte sources:
1. Multi-tissue collections that include liver + other organs
2. iPSC-derived hepatocytes / hepatocyte-like cells
3. Organoid-derived hepatocytes
4. Any collection with hepatocytes AND non-liver cell types from non-liver tissues

The goal: find a human dataset where hepatocytes come from liver but the other
5 target types (CD4/CD8 T, B, endothelial, macrophage) come from diverse tissues,
mimicking Tabula Sapiens tissue diversity.

Usage:
    python scripts/30c_hepatocyte_source_scout.py
"""

from __future__ import annotations

import cellxgene_census
import pandas as pd

CENSUS_VERSION = "2025-11-08"
ORGANISM = "Homo sapiens"

CHROMIUM_ASSAYS = {
    "10x 3' v1", "10x 3' v2", "10x 3' v3",
    "10x 3' transcription profiling",
    "10x 5' v1", "10x 5' v2",
    "10x 5' transcription profiling",
}

# All hepatocyte-related labels to search
HEPATOCYTE_LABELS = [
    "hepatocyte",
    "periportal region hepatocyte",
    "centrilobular region hepatocyte",
    "midzonal region hepatocyte",
    "hepatocyte-like cell",  # iPSC-derived
]

# Organoid / iPSC keywords to search in cell_type labels
IPSC_ORGANOID_KEYWORDS = [
    "hepatocyte-like",
    "iPSC",
    "organoid",
    "stem cell-derived",
]

TARGET_TYPES_OTHER = [
    "CD8-positive, alpha-beta T cell",
    "CD4-positive, alpha-beta T cell",
    "B cell",
    "endothelial cell",
    "macrophage",
]


def main():
    print("=" * 80)
    print("CellWarp — Hepatocyte Source Scout")
    print("=" * 80)

    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        datasets_df = census["census_info"]["datasets"].read().concat().to_pandas()
        ds_to_coll_id = dict(zip(datasets_df["dataset_id"], datasets_df["collection_id"]))
        ds_to_coll_name = dict(zip(datasets_df["dataset_id"], datasets_df["collection_name"]))

        # Tabula Sapiens exclusion
        ts_mask = datasets_df["collection_name"].str.contains(
            "Tabula Sapiens", case=False, na=False
        )
        ts_ds_ids = set(datasets_df.loc[ts_mask, "dataset_id"].tolist())

        # =================================================================
        # Part 1: All hepatocyte sources (any label)
        # =================================================================
        print("\n[Part 1] All hepatocyte sources in Census...")
        all_hep_obs = []
        for label in HEPATOCYTE_LABELS:
            print(f"\n  Querying: '{label}'...")
            try:
                obs = cellxgene_census.get_obs(
                    census, ORGANISM,
                    value_filter=(
                        f"cell_type == '{label}' "
                        f"and is_primary_data == True "
                        f"and disease == 'normal'"
                    ),
                    column_names=[
                        "soma_joinid", "cell_type", "dataset_id", "assay",
                        "development_stage", "tissue", "donor_id",
                    ],
                )
                print(f"    Found: {len(obs):,} cells")
                if len(obs) > 0:
                    all_hep_obs.append(obs)
            except Exception as e:
                print(f"    Error: {e}")

        if not all_hep_obs:
            print("\n  No hepatocytes found!")
            return

        hep_df = pd.concat(all_hep_obs, ignore_index=True)
        print(f"\n  Total hepatocyte-label cells: {len(hep_df):,}")

        # Adult filter
        adult_mask = hep_df["development_stage"].str.contains(
            "adult|year-old|years old", case=False, na=False
        )
        hep_df = hep_df[adult_mask]
        print(f"  Adult: {len(hep_df):,}")

        # Exclude Tabula Sapiens
        hep_df = hep_df[~hep_df["dataset_id"].isin(ts_ds_ids)]
        print(f"  Non-TS: {len(hep_df):,}")

        # Map to collections
        hep_df["collection_id"] = hep_df["dataset_id"].map(ds_to_coll_id)
        hep_df["collection_name"] = hep_df["dataset_id"].map(ds_to_coll_name)

        # Show ALL hepatocyte collections with assay and tissue breakdown
        print(f"\n  HEPATOCYTE COLLECTIONS (all assays, adult, normal, non-TS):")
        print(f"  {'Collection':<55} {'Cells':>7} {'Assays':<30} {'Tissues'}")
        print(f"  {'-'*120}")

        for coll_id, grp in hep_df.groupby("collection_id"):
            cname = grp["collection_name"].iloc[0]
            assays = ", ".join(sorted(grp["assay"].unique()))
            tissues = ", ".join(sorted(grp["tissue"].dropna().unique()))
            is_10x = any(a in CHROMIUM_ASSAYS for a in grp["assay"].unique())
            marker = " [10x]" if is_10x else ""
            print(
                f"  {str(cname)[:55]:<55} {len(grp):>7,} "
                f"{assays[:30]:<30} {tissues[:50]}{marker}"
            )

        # =================================================================
        # Part 2: iPSC / organoid hepatocytes
        # =================================================================
        print(f"\n\n[Part 2] iPSC / organoid hepatocyte search...")

        # Search for hepatocyte-like cells
        for kw in IPSC_ORGANOID_KEYWORDS:
            print(f"\n  Searching cell_type containing '{kw}'...")
            try:
                obs = cellxgene_census.get_obs(
                    census, ORGANISM,
                    value_filter=(
                        f"is_primary_data == True "
                        f"and disease == 'normal'"
                    ),
                    column_names=["soma_joinid", "cell_type", "dataset_id",
                                  "assay", "tissue", "development_stage"],
                )
                # This is too expensive — search in the hepatocyte data instead
                print(f"    (skipping full Census scan — too expensive)")
                break
            except Exception:
                break

        # Check tissue labels for organoid/iPSC signals in hepatocyte data
        print(f"\n  Checking tissue labels in hepatocyte data for organoid/iPSC signals:")
        all_tissues = sorted(hep_df["tissue"].dropna().unique())
        for t in all_tissues:
            n = len(hep_df[hep_df["tissue"] == t])
            flag = ""
            t_lower = t.lower()
            if any(kw in t_lower for kw in ["organoid", "ipsc", "vitro", "culture", "cell line"]):
                flag = " *** iPSC/ORGANOID ***"
            print(f"    {t}: {n:,}{flag}")

        # =================================================================
        # Part 3: Multi-tissue collections with hepatocytes
        # =================================================================
        print(f"\n\n[Part 3] Multi-tissue collections with hepatocytes...")
        print(f"  Checking which hepatocyte collections also have non-liver cell types...\n")

        # 10x hepatocyte collections only
        hep_10x = hep_df[hep_df["assay"].isin(CHROMIUM_ASSAYS)]
        hep_colls = hep_10x["collection_id"].unique()

        for coll_id in hep_colls:
            coll_name = ds_to_coll_name.get(
                hep_10x[hep_10x["collection_id"] == coll_id]["dataset_id"].iloc[0],
                "Unknown"
            )
            hep_count = len(hep_10x[hep_10x["collection_id"] == coll_id])

            # Get all dataset_ids for this collection
            coll_ds_ids = set(
                datasets_df[datasets_df["collection_id"] == coll_id]["dataset_id"].tolist()
            )

            print(f"\n  --- {str(coll_name)[:65]}")
            print(f"      Hepatocytes: {hep_count:,} (10x)")
            print(f"      Datasets: {len(coll_ds_ids)}")

            # Check other 5 target types in this collection
            other_type_coverage = {}
            for ct in TARGET_TYPES_OTHER:
                safe_ct = ct.replace("'", "\\'")
                try:
                    ct_obs = cellxgene_census.get_obs(
                        census, ORGANISM,
                        value_filter=(
                            f"cell_type == '{safe_ct}' "
                            f"and is_primary_data == True "
                            f"and disease == 'normal'"
                        ),
                        column_names=["soma_joinid", "dataset_id", "assay",
                                      "tissue", "development_stage"],
                    )
                    # Filter to this collection, adult, 10x
                    ct_obs = ct_obs[ct_obs["dataset_id"].isin(coll_ds_ids)]
                    adult_m = ct_obs["development_stage"].str.contains(
                        "adult|year-old|years old", case=False, na=False
                    )
                    ct_obs = ct_obs[adult_m & ct_obs["assay"].isin(CHROMIUM_ASSAYS)]

                    tissues = ", ".join(sorted(ct_obs["tissue"].dropna().unique())[:3])
                    other_type_coverage[ct] = {
                        "n": len(ct_obs),
                        "tissues": tissues,
                    }
                except Exception as e:
                    other_type_coverage[ct] = {"n": 0, "tissues": f"Error: {e}"}

            # Report
            n_pass = 0
            for ct in TARGET_TYPES_OTHER:
                info = other_type_coverage[ct]
                short = ct.split(",")[0] if "," in ct else ct
                status = "PASS" if info["n"] >= 500 else "fail"
                if info["n"] >= 500:
                    n_pass += 1
                print(
                    f"      {short:<35} {info['n']:>6,}  [{status}]  "
                    f"{info['tissues'][:40]}"
                )

            all_tissues_coll = set()
            for ct in TARGET_TYPES_OTHER:
                if other_type_coverage[ct]["tissues"] and "Error" not in other_type_coverage[ct]["tissues"]:
                    all_tissues_coll.update(other_type_coverage[ct]["tissues"].split(", "))
            # Check hepatocyte tissue
            hep_tissues = hep_10x[hep_10x["collection_id"] == coll_id]["tissue"].dropna().unique()
            all_tissues_coll.update(hep_tissues)

            non_liver = [t for t in all_tissues_coll if "liver" not in t.lower()]
            print(f"      Non-liver tissues present: {', '.join(sorted(non_liver)[:5]) if non_liver else 'NONE'}")
            print(f"      Target types ≥500: {n_pass}/5 (+ hepatocyte = {n_pass+1}/6)")

    print("\n" + "=" * 80)
    print("SCOUT COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
