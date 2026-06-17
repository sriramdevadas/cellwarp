"""
CellWarp — Cell Type Inventory Query

Queries CZ CELLxGENE Census for ALL cell types present in both Tabula Sapiens (human)
and Tabula Muris Senis (mouse) with ≥200 cells per species.

Biology
-------
We want to maximize the number of homologous cell type landmarks for Procrustes
analysis. More landmarks = more statistical power and a more complete picture of
the cross-species geometric transformation. This script inventories what's available
before committing to a download.

Filters: is_primary_data=True, disease="normal" (same as Phase 1 download).

Output: ./output/phase2/cell_type_inventory.csv
"""

from __future__ import annotations

from pathlib import Path

import cellxgene_census
import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

OUTPUT_DIR = Path("./output/phase2")
OUTPUT_CSV = OUTPUT_DIR / "cell_type_inventory.csv"

MIN_CELLS = 200  # Minimum cells per species to include

# We need dataset_ids to restrict to each atlas
HUMAN_COLLECTION = "Tabula Sapiens"
MOUSE_COLLECTION = "Tabula Muris Senis"

HUMAN_ORGANISM = "Homo sapiens"
MOUSE_ORGANISM = "Mus musculus"


def get_dataset_ids(census, collection_name: str) -> list[str]:
    """Look up dataset_ids for a collection."""
    datasets_df = census["census_info"]["datasets"].read().concat().to_pandas()
    mask = datasets_df["collection_name"] == collection_name
    if mask.sum() == 0:
        mask = datasets_df["collection_name"].str.contains(
            collection_name, case=False, na=False
        )
    if mask.sum() == 0:
        raise ValueError(f"Collection '{collection_name}' not found")
    matched = datasets_df.loc[mask]
    ids = matched["dataset_id"].tolist()
    actual = matched["collection_name"].iloc[0]
    print(f"  {actual}: {len(ids)} dataset(s)")
    return ids


def get_cell_type_counts(
    census, organism: str, dataset_ids: list[str]
) -> pd.Series:
    """
    Query Census for cell type counts with our standard filters.

    Returns a Series: index=cell_type, values=count, sorted descending.
    """
    # Build filter
    ids_str = ", ".join(f"'{d}'" for d in dataset_ids)
    value_filter = (
        f"is_primary_data == True and disease == 'normal' "
        f"and dataset_id in [{ids_str}]"
    )

    print(f"  Querying {organism} cell types...")
    obs_df = cellxgene_census.get_obs(
        census,
        organism,
        value_filter=value_filter,
        column_names=["cell_type"],
    )

    counts = obs_df["cell_type"].value_counts()
    # Remove zero-count categorical entries
    counts = counts[counts > 0]
    print(f"  Found {len(counts)} distinct cell types ({counts.sum():,} total cells)")
    return counts.sort_values(ascending=False)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("CellWarp — Cell Type Inventory")
    print("Querying CELLxGENE Census for shared cell types")
    print(f"Minimum cells per species: {MIN_CELLS}")
    print("=" * 70)

    print("\nOpening Census (stable release)...")
    with cellxgene_census.open_soma(census_version="2025-11-08") as census:
        # Get dataset IDs
        print("\nLooking up atlas collections...")
        human_ids = get_dataset_ids(census, HUMAN_COLLECTION)
        mouse_ids = get_dataset_ids(census, MOUSE_COLLECTION)

        # Get cell type counts
        print("\nQuerying cell type counts...")
        human_counts = get_cell_type_counts(census, HUMAN_ORGANISM, human_ids)
        mouse_counts = get_cell_type_counts(census, MOUSE_ORGANISM, mouse_ids)

    # Build inventory: cell types present in BOTH species
    human_types = set(human_counts.index)
    mouse_types = set(mouse_counts.index)
    shared_types = human_types & mouse_types

    print(f"\n{'=' * 70}")
    print(f"Human-only cell types: {len(human_types - mouse_types)}")
    print(f"Mouse-only cell types: {len(mouse_types - human_types)}")
    print(f"Shared cell types (any count): {len(shared_types)}")

    # Build dataframe for shared types
    rows = []
    for ct in sorted(shared_types):
        h = int(human_counts.get(ct, 0))
        m = int(mouse_counts.get(ct, 0))
        rows.append({
            "cell_type": ct,
            "human_count": h,
            "mouse_count": m,
            "min_count": min(h, m),
            "passes_200_gate": min(h, m) >= MIN_CELLS,
        })

    df = pd.DataFrame(rows)
    df = df.sort_values("min_count", ascending=False).reset_index(drop=True)

    # Filter to ≥200 in both species
    passing = df[df["passes_200_gate"]].copy()
    failing = df[~df["passes_200_gate"]].copy()

    # Print passing types
    print(f"\n{'=' * 70}")
    print(f"SHARED CELL TYPES WITH ≥{MIN_CELLS} CELLS IN BOTH SPECIES: {len(passing)}")
    print(f"{'=' * 70}")
    print(f"\n{'Cell Type':<45} {'Human':>8} {'Mouse':>8} {'Min':>8}")
    print("-" * 71)
    for _, row in passing.iterrows():
        ct = row["cell_type"]
        # Mark our current 6 types
        marker = " *" if ct in {
            "hepatocyte", "CD8-positive, alpha-beta T cell",
            "endothelial cell", "CD4-positive, alpha-beta T cell",
            "B cell", "macrophage",
        } else ""
        print(f"{ct:<45} {row['human_count']:>8,} {row['mouse_count']:>8,} {row['min_count']:>8,}{marker}")
    print("-" * 71)
    print(f"{'TOTAL':<45} {passing['human_count'].sum():>8,} {passing['mouse_count'].sum():>8,}")
    print(f"\n* = currently in our Phase 1-2 analysis (6 types)")

    # Print failing types (close misses)
    if len(failing) > 0:
        close_misses = failing[failing["min_count"] >= 50].sort_values(
            "min_count", ascending=False
        )
        if len(close_misses) > 0:
            print(f"\n{'=' * 70}")
            print(f"CLOSE MISSES (50-{MIN_CELLS-1} cells in one species): {len(close_misses)}")
            print(f"{'=' * 70}")
            print(f"\n{'Cell Type':<45} {'Human':>8} {'Mouse':>8} {'Min':>8}")
            print("-" * 71)
            for _, row in close_misses.iterrows():
                print(f"{row['cell_type']:<45} {row['human_count']:>8,} {row['mouse_count']:>8,} {row['min_count']:>8,}")

    # Print summary stats
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print(f"{'=' * 70}")
    print(f"Total distinct cell types in Tabula Sapiens:    {len(human_types)}")
    print(f"Total distinct cell types in Tabula Muris Senis: {len(mouse_types)}")
    print(f"Shared (present in both):                       {len(shared_types)}")
    print(f"Shared with ≥{MIN_CELLS} cells in both:               {len(passing)}")
    print(f"Currently using:                                6")
    print(f"Potential expansion:                            +{len(passing) - 6} types")

    # Also show types that pass ≥500 (our current Phase 1 gate)
    passing_500 = df[df["min_count"] >= 500]
    print(f"Shared with ≥500 cells in both:                {len(passing_500)}")

    # Save full inventory
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"\nSaved full inventory to {OUTPUT_CSV}")

    # Also save just the passing types for easy reference
    passing_csv = OUTPUT_DIR / "cell_type_inventory_passing.csv"
    passing.to_csv(passing_csv, index=False)
    print(f"Saved passing types (≥{MIN_CELLS}) to {passing_csv}")


if __name__ == "__main__":
    main()
