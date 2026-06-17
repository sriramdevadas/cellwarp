"""
Census Exploration — Step 1
============================
Surveys all available collections in CELLxGENE Census (build 2025-11-08) to
identify independent cross-species datasets for Procrustes replication.

Goal: Find a large set of human + mouse cells from collections NOT used in
the primary analysis (Tabula Sapiens, Tabula Muris Senis) or existing
replications (CellHint, PanSci, Sun2023 [external]).
"""

import cellxgene_census
import pandas as pd
from pathlib import Path

OUT_DIR = Path(__file__).parent

# Collections to exclude (used in primary analysis or existing replications)
EXCLUDED_COLLECTIONS = {
    # Primary analysis
    "Tabula Sapiens",
    # Primary analysis (mouse)
    "Tabula Muris Senis",
    # Existing replications
    "CellHint: harmonised cell-type annotations for interoperable single-cell and spatial genomic atlases",
    # PanSci
    "A single-cell atlas of the mouse brain vasculature",  # not actually PanSci
}

# PanSci dataset IDs — exclude these specifically
PANSCI_DATASET_IDS = set()  # Will be populated from census

# Cell types from the primary 35-type analysis
PRIMARY_CELL_TYPES = [
    "stromal cell", "epithelial cell", "hematopoietic precursor cell",
    "hematopoietic stem cell", "pancreatic acinar cell", "basal cell",
    "T cell", "neutrophil", "fibroblast of cardiac tissue",
    "myeloid leukocyte", "mesenchymal stem cell of adipose tissue",
    "plasma cell", "mesenchymal stem cell",
    "CD4-positive, alpha-beta T cell", "classical monocyte",
    "macrophage", "B cell",
    "luminal epithelial cell of mammary gland",
    "large intestine goblet cell",
    "enterocyte of epithelium of large intestine",
    "myeloid dendritic cell", "monocyte", "natural killer cell",
    "intermediate monocyte", "mature NK T cell", "adventitial cell",
    "granulocyte", "fibroblast", "bladder urothelial cell",
    "pancreatic ductal cell", "smooth muscle cell", "hepatocyte",
    "endothelial cell", "non-classical monocyte",
    "CD8-positive, alpha-beta T cell",
]


def main():
    print("=" * 70)
    print("CENSUS EXPLORATION — Identifying independent replication datasets")
    print("=" * 70)

    with cellxgene_census.open_soma(census_version="2025-11-08") as census:
        # ── Step 1: List all datasets ─────────────────────────────────────
        print("\n[1] Loading dataset catalog...")
        datasets_df = census["census_info"]["datasets"].read().concat().to_pandas()
        print(f"  Total datasets in Census: {len(datasets_df)}")
        print(f"  Total collections: {datasets_df['collection_name'].nunique()}")

        # Save full dataset catalog
        datasets_df.to_csv(OUT_DIR / "census_datasets_full.csv", index=False)

        # ── Step 2: Identify collections with both human and mouse ────────
        print("\n[2] Identifying collections with both human and mouse data...")

        # Get organism info per dataset
        human_datasets = set()
        mouse_datasets = set()

        # Check which datasets have human data
        print("  Querying human cell counts per dataset...")
        human_obs = cellxgene_census.get_obs(
            census, "Homo sapiens",
            value_filter="is_primary_data == True and disease == 'normal'",
            column_names=["dataset_id", "cell_type"],
        )
        human_per_dataset = human_obs.groupby("dataset_id").size().reset_index(name="n_human_cells")
        human_datasets = set(human_per_dataset["dataset_id"])
        print(f"  Human: {len(human_obs):,} healthy primary cells across {len(human_datasets)} datasets")

        # Human cell type counts per dataset
        human_ct_counts = human_obs.groupby(["dataset_id", "cell_type"]).size().reset_index(name="n_cells")

        print("  Querying mouse cell counts per dataset...")
        mouse_obs = cellxgene_census.get_obs(
            census, "Mus musculus",
            value_filter="is_primary_data == True and disease == 'normal'",
            column_names=["dataset_id", "cell_type"],
        )
        mouse_per_dataset = mouse_obs.groupby("dataset_id").size().reset_index(name="n_mouse_cells")
        mouse_datasets = set(mouse_per_dataset["dataset_id"])
        print(f"  Mouse: {len(mouse_obs):,} healthy primary cells across {len(mouse_datasets)} datasets")

        # Mouse cell type counts per dataset
        mouse_ct_counts = mouse_obs.groupby(["dataset_id", "cell_type"]).size().reset_index(name="n_cells")

        # ── Step 3: Map datasets to collections ──────────────────────────
        print("\n[3] Mapping datasets to collections...")

        # Merge dataset info
        ds_info = datasets_df[["dataset_id", "collection_name", "dataset_title"]].copy()

        # Human collections
        human_collections = (
            ds_info[ds_info["dataset_id"].isin(human_datasets)]
            .merge(human_per_dataset, on="dataset_id")
        )
        human_by_collection = (
            human_collections.groupby("collection_name")
            .agg(n_datasets=("dataset_id", "count"), total_human_cells=("n_human_cells", "sum"))
            .reset_index()
        )

        # Mouse collections
        mouse_collections = (
            ds_info[ds_info["dataset_id"].isin(mouse_datasets)]
            .merge(mouse_per_dataset, on="dataset_id")
        )
        mouse_by_collection = (
            mouse_collections.groupby("collection_name")
            .agg(n_datasets=("dataset_id", "count"), total_mouse_cells=("n_mouse_cells", "sum"))
            .reset_index()
        )

        # Find collections with BOTH species
        both = human_by_collection.merge(
            mouse_by_collection, on="collection_name", suffixes=("_human", "_mouse")
        )
        both["total_cells"] = both["total_human_cells"] + both["total_mouse_cells"]
        both = both.sort_values("total_cells", ascending=False)

        print(f"\n  Collections with BOTH human and mouse healthy primary data: {len(both)}")
        print(f"\n  {'Collection':<60} {'Human':>10} {'Mouse':>10} {'Total':>10}")
        print("  " + "-" * 95)
        for _, row in both.head(20).iterrows():
            name = row["collection_name"][:58]
            excluded = " [EXCLUDED]" if any(ex in row["collection_name"] for ex in EXCLUDED_COLLECTIONS) else ""
            print(f"  {name:<60} {row['total_human_cells']:>10,} {row['total_mouse_cells']:>10,} {row['total_cells']:>10,}{excluded}")

        both.to_csv(OUT_DIR / "collections_both_species.csv", index=False)

        # ── Step 4: "Pan-Census" approach ─────────────────────────────────
        # Instead of restricting to a single collection, pool ALL non-excluded data
        print("\n[4] Pan-Census approach: pooling all independent data...")

        # Identify excluded dataset IDs
        excluded_ds_ids = set()
        for ex_name in EXCLUDED_COLLECTIONS:
            mask = datasets_df["collection_name"].str.contains(ex_name, case=False, na=False)
            excluded_ds_ids.update(datasets_df.loc[mask, "dataset_id"].tolist())

        # Also find PanSci datasets
        pansci_mask = datasets_df["collection_name"].str.contains("PanSci|EasySci|Cao", case=False, na=False)
        pansci_ids = set(datasets_df.loc[pansci_mask, "dataset_id"].tolist())
        excluded_ds_ids.update(pansci_ids)
        if pansci_ids:
            print(f"  Excluding PanSci datasets: {len(pansci_ids)}")

        print(f"  Total excluded dataset IDs: {len(excluded_ds_ids)}")

        # Filter to independent data
        human_independent = human_ct_counts[~human_ct_counts["dataset_id"].isin(excluded_ds_ids)]
        mouse_independent = mouse_ct_counts[~mouse_ct_counts["dataset_id"].isin(excluded_ds_ids)]

        # Count cells per cell type (independent data only)
        human_by_ct = human_independent.groupby("cell_type")["n_cells"].sum().reset_index()
        human_by_ct.columns = ["cell_type", "human_cells"]
        mouse_by_ct = mouse_independent.groupby("cell_type")["n_cells"].sum().reset_index()
        mouse_by_ct.columns = ["cell_type", "mouse_cells"]

        # Merge to find types with both species
        ct_both = human_by_ct.merge(mouse_by_ct, on="cell_type", how="inner")
        ct_both["total"] = ct_both["human_cells"] + ct_both["mouse_cells"]
        ct_both["min_species"] = ct_both[["human_cells", "mouse_cells"]].min(axis=1)
        ct_both = ct_both.sort_values("min_species", ascending=False)

        print(f"\n  Cell types with BOTH human and mouse independent data: {ct_both.shape[0]}")

        # Filter to types that overlap with primary analysis
        ct_overlap = ct_both[ct_both["cell_type"].isin(PRIMARY_CELL_TYPES)].copy()
        ct_overlap = ct_overlap.sort_values("min_species", ascending=False)

        print(f"  Of which overlap with primary 35 types: {len(ct_overlap)}")
        print(f"\n  {'Cell Type':<50} {'Human':>8} {'Mouse':>8} {'Min':>8}")
        print("  " + "-" * 78)
        for _, row in ct_overlap.iterrows():
            gate = "✓" if row["min_species"] >= 500 else "✗"
            print(f"  {row['cell_type']:<50} {row['human_cells']:>8,} {row['mouse_cells']:>8,} {row['min_species']:>8,} {gate}")

        # Types passing the 500-cell gate
        ct_passing = ct_overlap[ct_overlap["min_species"] >= 500]
        print(f"\n  Types passing ≥500 cells/species gate: {len(ct_passing)}")

        ct_overlap.to_csv(OUT_DIR / "independent_cell_type_counts.csv", index=False)
        ct_passing.to_csv(OUT_DIR / "passing_cell_types.csv", index=False)

        # ── Step 5: Collection-level breakdown for independent data ───────
        print("\n[5] Top independent collections contributing to passing cell types...")

        passing_types = set(ct_passing["cell_type"])

        # For each passing type, which collections contribute?
        human_passing = human_independent[human_independent["cell_type"].isin(passing_types)]
        mouse_passing = mouse_independent[mouse_independent["cell_type"].isin(passing_types)]

        # Map dataset_id → collection_name
        ds_to_collection = dict(zip(datasets_df["dataset_id"], datasets_df["collection_name"]))

        human_passing = human_passing.copy()
        human_passing["collection"] = human_passing["dataset_id"].map(ds_to_collection)
        mouse_passing = mouse_passing.copy()
        mouse_passing["collection"] = mouse_passing["dataset_id"].map(ds_to_collection)

        # Top human collections
        h_coll = human_passing.groupby("collection")["n_cells"].sum().sort_values(ascending=False)
        print(f"\n  Top human collections (for passing cell types):")
        for coll, n in h_coll.head(10).items():
            print(f"    {coll[:65]:<68} {n:>8,}")

        # Top mouse collections
        m_coll = mouse_passing.groupby("collection")["n_cells"].sum().sort_values(ascending=False)
        print(f"\n  Top mouse collections (for passing cell types):")
        for coll, n in m_coll.head(10).items():
            print(f"    {coll[:65]:<68} {n:>8,}")

        # ── Step 6: Single-collection candidates ─────────────────────────
        print("\n[6] Single-collection cross-species candidates (independent)...")

        # Filter the 'both' table to exclude known collections
        independent_both = both[~both["collection_name"].apply(
            lambda x: any(ex.lower() in x.lower() for ex in EXCLUDED_COLLECTIONS)
        )].copy()

        # For each, count overlapping cell types
        for idx, row in independent_both.iterrows():
            coll_name = row["collection_name"]
            coll_ds_ids = set(datasets_df.loc[
                datasets_df["collection_name"] == coll_name, "dataset_id"
            ])

            h_types = set(human_ct_counts.loc[
                human_ct_counts["dataset_id"].isin(coll_ds_ids), "cell_type"
            ])
            m_types = set(mouse_ct_counts.loc[
                mouse_ct_counts["dataset_id"].isin(coll_ds_ids), "cell_type"
            ])

            shared_with_primary = (h_types & m_types) & set(PRIMARY_CELL_TYPES)
            independent_both.loc[idx, "n_shared_primary_types"] = len(shared_with_primary)
            independent_both.loc[idx, "shared_types"] = "; ".join(sorted(shared_with_primary))

        independent_both = independent_both.sort_values("n_shared_primary_types", ascending=False)
        independent_both.to_csv(OUT_DIR / "independent_cross_species_collections.csv", index=False)

        print(f"\n  {'Collection':<55} {'H cells':>8} {'M cells':>8} {'Shared types':>12}")
        print("  " + "-" * 88)
        for _, row in independent_both.head(15).iterrows():
            n_shared = int(row.get("n_shared_primary_types", 0))
            print(f"  {row['collection_name'][:53]:<55} "
                  f"{row['total_human_cells']:>8,} {row['total_mouse_cells']:>8,} "
                  f"{n_shared:>12}")

    print("\n" + "=" * 70)
    print("EXPLORATION COMPLETE — Files saved to", OUT_DIR)
    print("=" * 70)


if __name__ == "__main__":
    main()
