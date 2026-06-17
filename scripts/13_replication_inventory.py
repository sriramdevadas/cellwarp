#!/usr/bin/env python3
"""
CellWarp — Independent Replication Atlas Inventory

Queries CELLxGENE Census to determine whether sufficient non-Tabula data exists
for an independent replication of the 35-type Procrustes pipeline on a completely
different human/mouse atlas pair.

Biology
-------
Independent replication is the gold standard for validating our cross-species
geometric morphometric findings. If the Procrustes deformation pattern replicates
on entirely different atlases (no Tabula Sapiens, no Tabula Muris Senis), the
result cannot be attributed to atlas-specific batch effects or protocol artifacts.

Math
----
No computation — this is a feasibility census. For each of the 35 cell types, we
count available cells in non-Tabula datasets, applying the same ≥500-cell gate
used in the main analysis.

Steps
-----
1. Load the 35 cell type names from the main analysis results.
2. Query Census for human cells (exclude Tabula Sapiens), adult-only.
3. Query Census for mouse cells (exclude Tabula Muris Senis and Tabula Muris).
4. Print inventory table with PASS/MARGINAL/FAIL per cell type.
5. Summarize feasibility and identify gaps.
6. Report top contributing non-Tabula collections per species.

Outputs:
    output/validation/replication_inventory.txt  — Human-readable report
    output/validation/replication_inventory.csv  — Machine-readable counts table
"""

from __future__ import annotations

import json
import sys
import time
from io import StringIO
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import cellxgene_census
import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

RESULTS_JSON = Path("output/phase2/scaled_35types/procrustes_results_35.json")
OUTPUT_DIR = Path("output/validation")

HUMAN_ORGANISM = "Homo sapiens"
MOUSE_ORGANISM = "Mus musculus"

# Collection substrings to EXCLUDE via dataset_id lookup (DECISION-014 pattern).
# We search census_info/datasets for collection_name containing these substrings,
# collect all matching dataset_ids, then exclude obs rows by dataset_id.
# "Tabula Muris" catches both "Tabula Muris" and "Tabula Muris Senis".
EXCLUDE_HUMAN_COLLECTION_SUBSTR = ["Tabula Sapiens"]
EXCLUDE_MOUSE_COLLECTION_SUBSTR = ["Tabula Muris"]

MIN_CELLS_PASS = 500
MIN_CELLS_MARGINAL = 200


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def log(msg: str = "") -> None:
    """Print with immediate flush so output is visible in non-interactive mode."""
    print(msg, flush=True)


def load_cell_types() -> list[str]:
    """Load the 35 cell type names from the main Procrustes results JSON."""
    with open(RESULTS_JSON) as f:
        data = json.load(f)
    cell_types = data["cell_types"]
    log(f"  Loaded {len(cell_types)} cell types from {RESULTS_JSON}")
    return cell_types


def get_excluded_dataset_ids(
    datasets_df: pd.DataFrame,
    collection_substrings: list[str],
) -> set[str]:
    """Find dataset_ids belonging to excluded collections via substring match.

    Census obs metadata does not include collection_name, so we query
    census_info/datasets for collection_name *containing* the substring,
    collect all matching dataset_ids, then use those IDs to filter obs rows.
    This is the DECISION-014 pattern — ensures no Tabula data leaks through
    even if collection_name has variant suffixes or versioning.

    Args:
        datasets_df: Full Census datasets table.
        collection_substrings: Substrings to match against collection_name.
            E.g. "Tabula Muris" catches both "Tabula Muris" and
            "Tabula Muris Senis".

    Returns:
        Set of dataset_id strings to exclude.
    """
    excluded = set()
    for substr in collection_substrings:
        mask = datasets_df["collection_name"].str.contains(
            substr, case=False, na=False
        )
        matched = datasets_df.loc[mask]
        ids = set(matched["dataset_id"].tolist())

        # Print matched collection names for auditability
        matched_names = sorted(matched["collection_name"].unique())
        log(f"    Substring '{substr}' matched {len(ids)} dataset(s) "
            f"from {len(matched_names)} collection(s):")
        for cname in matched_names:
            n = (matched["collection_name"] == cname).sum()
            log(f"      - '{cname}' ({n} datasets)")

        if len(ids) == 0:
            log(f"    WARNING: No datasets matched '{substr}' — "
                f"Tabula data may leak into replication!")

        excluded |= ids
    return excluded


def query_species_obs(
    census,
    organism: str,
    cell_types: list[str],
    excluded_dataset_ids: set[str],
    apply_adult_filter: bool,
) -> pd.DataFrame:
    """Query Census obs for all 35 cell types in a single batched query.

    Uses cell_type in [...] to fetch all types at once (1 query instead of 35),
    then post-filters to exclude Tabula datasets by dataset_id and optionally
    restrict to adult development stage.

    Args:
        census: Open Census SOMA collection.
        organism: "Homo sapiens" or "Mus musculus".
        cell_types: List of cell type names to query.
        excluded_dataset_ids: Dataset IDs to exclude from counts.
        apply_adult_filter: If True, filter to rows where development_stage
            contains 'adult'.

    Returns:
        Filtered DataFrame with columns: cell_type, dataset_id, development_stage.
    """
    # Build batched filter: cell_type in ['type1', 'type2', ...]
    names_str = ", ".join(f"'{ct}'" for ct in cell_types)
    value_filter = (
        f"cell_type in [{names_str}] "
        f"and is_primary_data == True "
        f"and disease == 'normal'"
    )

    log(f"  Querying {organism} (single batched query for {len(cell_types)} types)...")
    t0 = time.time()

    obs = cellxgene_census.get_obs(
        census,
        organism,
        value_filter=value_filter,
        column_names=["cell_type", "dataset_id", "development_stage"],
    )

    dt = time.time() - t0
    n_total = len(obs)
    log(f"  Raw Census result: {n_total:,} cells [{dt:.1f}s]")

    # Exclude Tabula datasets by dataset_id (DECISION-014 pattern)
    obs = obs[~obs["dataset_id"].isin(excluded_dataset_ids)]
    n_after_excl = len(obs)
    log(f"  After Tabula exclusion: {n_after_excl:,} cells "
        f"(removed {n_total - n_after_excl:,})")

    # Adult-only filter
    if apply_adult_filter and "development_stage" in obs.columns:
        n_with_stage = obs["development_stage"].notna().sum()
        if n_with_stage > 0:
            adult_mask = obs["development_stage"].str.contains(
                "adult", case=False, na=False
            )
            n_before = len(obs)
            obs = obs[adult_mask].copy()
            log(f"  After adult filter: {len(obs):,} cells "
                f"(removed {n_before - len(obs):,} non-adult)")
        else:
            log(f"  development_stage not populated — skipping adult filter")

    # Per-type breakdown
    type_counts = obs["cell_type"].value_counts()
    for ct in cell_types:
        n = type_counts.get(ct, 0)
        status_char = "  " if n >= MIN_CELLS_PASS else " *" if n >= MIN_CELLS_MARGINAL else " <<"
        log(f"    {ct:<50} {n:>8,}{status_char}")

    return obs


def classify_status(human_n: int, mouse_n: int) -> str:
    """Classify a cell type as PASS, MARGINAL, or FAIL."""
    if human_n >= MIN_CELLS_PASS and mouse_n >= MIN_CELLS_PASS:
        return "PASS"
    if human_n < MIN_CELLS_MARGINAL or mouse_n < MIN_CELLS_MARGINAL:
        return "FAIL"
    return "MARGINAL"


def gap_reason(human_n: int, mouse_n: int) -> str:
    """Describe why a cell type fails — human gap, mouse gap, or both."""
    reasons = []
    if human_n < MIN_CELLS_PASS:
        reasons.append(f"human={human_n}")
    if mouse_n < MIN_CELLS_PASS:
        reasons.append(f"mouse={mouse_n}")
    return ", ".join(reasons) if reasons else ""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Capture all output to both console and file
    report = StringIO()

    def tee(msg: str = ""):
        """Print to both stdout and the report buffer."""
        print(msg, flush=True)
        report.write(msg + "\n")

    # ------------------------------------------------------------------
    # STEP 1 — Load cell types
    # ------------------------------------------------------------------
    tee("=" * 78)
    tee("STEP 1 — Load 35 Cell Types from Main Analysis")
    tee("=" * 78)
    cell_types = load_cell_types()
    tee(f"\n  Cell types: {len(cell_types)}")
    for i, ct in enumerate(cell_types, 1):
        tee(f"    {i:>2}. {ct}")

    # ------------------------------------------------------------------
    # STEP 2 & 3 — Query Census
    # ------------------------------------------------------------------
    tee("\n" + "=" * 78)
    tee("STEP 2 — Query CELLxGENE Census (census_version='2025-11-08')")
    tee("=" * 78)

    with cellxgene_census.open_soma(census_version="2025-11-08") as census:
        # Load datasets table for collection → dataset_id mapping
        tee("\n  Loading Census datasets table...")
        datasets_df = (
            census["census_info"]["datasets"].read().concat().to_pandas()
        )
        tee(f"  Census contains {len(datasets_df)} datasets from "
            f"{datasets_df['collection_name'].nunique()} collections")

        # Build dataset_id → collection_name lookup
        ds_to_coll = dict(
            zip(datasets_df["dataset_id"], datasets_df["collection_name"])
        )

        # Identify excluded dataset IDs
        tee("\n  --- Human exclusions ---")
        human_excluded = get_excluded_dataset_ids(
            datasets_df, EXCLUDE_HUMAN_COLLECTION_SUBSTR
        )
        tee(f"  Total human-excluded dataset IDs: {len(human_excluded)}")

        tee("\n  --- Mouse exclusions ---")
        mouse_excluded = get_excluded_dataset_ids(
            datasets_df, EXCLUDE_MOUSE_COLLECTION_SUBSTR
        )
        tee(f"  Total mouse-excluded dataset IDs: {len(mouse_excluded)}")

        # ------------------------------------------------------------------
        # STEP 2a — Human (adult-only, exclude Tabula Sapiens)
        # ------------------------------------------------------------------
        tee("\n" + "-" * 78)
        tee("  HUMAN — adult, normal, is_primary_data, non-Tabula Sapiens")
        tee("-" * 78)

        human_obs = query_species_obs(
            census,
            HUMAN_ORGANISM,
            cell_types,
            human_excluded,
            apply_adult_filter=True,
        )

        # ------------------------------------------------------------------
        # STEP 2b — Mouse (exclude Tabula Muris / TMS)
        # ------------------------------------------------------------------
        tee("\n" + "-" * 78)
        tee("  MOUSE — normal, is_primary_data, non-Tabula Muris (Senis)")
        tee("  Note: adult filter applied if development_stage is populated")
        tee("-" * 78)

        mouse_obs = query_species_obs(
            census,
            MOUSE_ORGANISM,
            cell_types,
            mouse_excluded,
            apply_adult_filter=True,
        )

    # Census connection closed — work with collected DataFrames
    tee("\n  Census connection closed.")

    # ------------------------------------------------------------------
    # Derive counts from collected obs
    # ------------------------------------------------------------------
    human_counts = (
        human_obs.groupby("cell_type", observed=True).size().to_dict()
        if len(human_obs) > 0
        else {}
    )
    mouse_counts = (
        mouse_obs.groupby("cell_type", observed=True).size().to_dict()
        if len(mouse_obs) > 0
        else {}
    )

    # ------------------------------------------------------------------
    # STEP 3 — Build results table
    # ------------------------------------------------------------------
    tee("\n" + "=" * 78)
    tee("STEP 3 — Inventory Results")
    tee("=" * 78)

    rows = []
    for ct in cell_types:
        h = human_counts.get(ct, 0)
        m = mouse_counts.get(ct, 0)
        status = classify_status(h, m)
        gap = gap_reason(h, m)
        rows.append({
            "cell_type": ct,
            "human_non_tabula": h,
            "mouse_non_tabula": m,
            "status": status,
            "gap_detail": gap,
        })

    df = pd.DataFrame(rows)

    # Print formatted table
    tee(f"\n  {'Cell Type':<50} {'Human':>8} {'Mouse':>8}  {'Status':<10}")
    tee("  " + "-" * 82)
    for _, row in df.iterrows():
        marker = ""
        if row["status"] == "FAIL":
            marker = " <<"
        elif row["status"] == "MARGINAL":
            marker = " *"
        tee(
            f"  {row['cell_type']:<50} "
            f"{row['human_non_tabula']:>8,} "
            f"{row['mouse_non_tabula']:>8,}  "
            f"{row['status']:<10}{marker}"
        )
    tee("  " + "-" * 82)
    tee(f"  {'TOTAL':<50} "
        f"{df['human_non_tabula'].sum():>8,} "
        f"{df['mouse_non_tabula'].sum():>8,}")

    # ------------------------------------------------------------------
    # STEP 4 — Feasibility summary
    # ------------------------------------------------------------------
    tee("\n" + "=" * 78)
    tee("STEP 4 — Feasibility Summary")
    tee("=" * 78)

    n_pass = (df["status"] == "PASS").sum()
    n_marginal = (df["status"] == "MARGINAL").sum()
    n_fail = (df["status"] == "FAIL").sum()

    tee(f"\n  PASS (≥{MIN_CELLS_PASS} both species):     {n_pass}/35")
    tee(f"  MARGINAL ({MIN_CELLS_MARGINAL}-{MIN_CELLS_PASS - 1} either): "
        f"      {n_marginal}/35")
    tee(f"  FAIL (<{MIN_CELLS_MARGINAL} either):            {n_fail}/35")

    # Detail on failures
    if n_fail > 0:
        tee(f"\n  --- Types that FAIL ---")
        for _, row in df[df["status"] == "FAIL"].iterrows():
            tee(f"    {row['cell_type']:<50} {row['gap_detail']}")

    if n_marginal > 0:
        tee(f"\n  --- Types that are MARGINAL ---")
        for _, row in df[df["status"] == "MARGINAL"].iterrows():
            tee(f"    {row['cell_type']:<50} {row['gap_detail']}")

    # Human-only and mouse-only gaps
    human_gap = df[
        (df["human_non_tabula"] < MIN_CELLS_PASS)
        & (df["mouse_non_tabula"] >= MIN_CELLS_PASS)
    ]
    mouse_gap = df[
        (df["mouse_non_tabula"] < MIN_CELLS_PASS)
        & (df["human_non_tabula"] >= MIN_CELLS_PASS)
    ]
    both_gap = df[
        (df["human_non_tabula"] < MIN_CELLS_PASS)
        & (df["mouse_non_tabula"] < MIN_CELLS_PASS)
    ]

    tee(f"\n  Gap breakdown:")
    tee(f"    Human gap only:  {len(human_gap)} types")
    tee(f"    Mouse gap only:  {len(mouse_gap)} types")
    tee(f"    Both species:    {len(both_gap)} types")

    # Maximum feasible n for independent replication
    tee(f"\n  Maximum feasible n for independent replication:")
    tee(f"    If requiring ≥500 both species: {n_pass} cell types")
    n_pass_200 = ((df["human_non_tabula"] >= MIN_CELLS_MARGINAL)
                  & (df["mouse_non_tabula"] >= MIN_CELLS_MARGINAL)).sum()
    tee(f"    If relaxing to ≥200 both species: {n_pass_200} cell types")

    # ------------------------------------------------------------------
    # STEP 5 — Top contributing collections
    # ------------------------------------------------------------------
    tee("\n" + "=" * 78)
    tee("STEP 5 — Top Contributing Collections (non-Tabula)")
    tee("=" * 78)

    for species, obs_df, label in [
        ("HUMAN", human_obs, "Human"),
        ("MOUSE", mouse_obs, "Mouse"),
    ]:
        tee(f"\n  --- {species}: Top 5 collections by cell count ---")

        if len(obs_df) == 0:
            tee(f"    No non-Tabula cells found for {label}")
            continue

        obs_with_coll = obs_df[["cell_type", "dataset_id"]].copy()
        obs_with_coll["collection_name"] = obs_with_coll["dataset_id"].map(ds_to_coll)

        coll_stats = (
            obs_with_coll.groupby("collection_name", observed=True)
            .agg(
                n_cells=("cell_type", "size"),
                n_cell_types=("cell_type", "nunique"),
                n_datasets=("dataset_id", "nunique"),
            )
            .sort_values("n_cells", ascending=False)
            .reset_index()
        )

        for _, row in coll_stats.head(5).iterrows():
            tee(
                f"    {row['n_cells']:>9,} cells | "
                f"{row['n_cell_types']:>2} types | "
                f"{row['n_datasets']:>2} datasets | "
                f"{row['collection_name']}"
            )

        if len(coll_stats) > 5:
            rest_cells = coll_stats.iloc[5:]["n_cells"].sum()
            rest_n = len(coll_stats) - 5
            tee(f"    {rest_cells:>9,} cells | ... {rest_n} more collections")

        # Coverage concentration
        total = coll_stats["n_cells"].sum()
        if total > 0:
            top1_pct = coll_stats.iloc[0]["n_cells"] / total * 100
            top3_sum = coll_stats.head(3)["n_cells"].sum()
            top3_pct = top3_sum / total * 100
            tee(
                f"\n    Concentration: top-1 = {top1_pct:.1f}%, "
                f"top-3 = {top3_pct:.1f}% "
                f"({len(coll_stats)} total collections)"
            )

    # ------------------------------------------------------------------
    # Save outputs
    # ------------------------------------------------------------------
    tee("\n" + "=" * 78)
    tee("Saving outputs")
    tee("=" * 78)

    # CSV
    csv_path = OUTPUT_DIR / "replication_inventory.csv"
    df.to_csv(csv_path, index=False)
    tee(f"  Saved: {csv_path}")

    # Text report
    txt_path = OUTPUT_DIR / "replication_inventory.txt"
    with open(txt_path, "w") as f:
        f.write(report.getvalue())
    tee(f"  Saved: {txt_path}")

    tee("\n  Done. This is a feasibility check — no data was downloaded.")


if __name__ == "__main__":
    main()
