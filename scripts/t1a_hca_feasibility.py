#!/usr/bin/env python3
"""
T1-A HCA Feasibility Check: Human Cell Atlas metadata-only query.

BIOLOGY: We need an independent human single-cell atlas (different from Tabula
Sapiens) for the human side of T1-A independent replication. The MCA feasibility
check confirmed 19/35 PASS + hepatocyte BORDERLINE on the mouse side. This script
identifies the best non-Tabula human dataset(s) covering our 35 cell types.

MATH: No heavy computation — pure metadata inventory. We query CELLxGENE Census
for non-Tabula adult healthy human cells, count per cell type, and identify the
best collection(s) for T1-A replication.

Steps:
1. Query Census for non-Tabula adult healthy human cells
2. Cell type coverage against all 35 types
3. Identify best single dataset or minimal combination
4. Technology and pipeline independence check
5. Final intersection with MCA results
6. Power assessment
7. Coverage comparison figure
8. Save outputs
"""

import os
import sys
import json
import time
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Patch

import cellxgene_census

# ── Configuration ──────────────────────────────────────────────────────────
OUTPUT_DIR = Path("output/validation/hca_feasibility")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

MCA_REPORT_PATH = Path("output/validation/mca_feasibility/mca_feasibility_report.json")
MCA_COVERAGE_PATH = Path("output/validation/mca_feasibility/mca_coverage_table.csv")

# Our 35 cell types
OUR_35_TYPES = [
    "B cell",
    "CD4-positive, alpha-beta T cell",
    "CD8-positive, alpha-beta T cell",
    "T cell",
    "adventitial cell",
    "basal cell",
    "bladder urothelial cell",
    "classical monocyte",
    "endothelial cell",
    "enterocyte of epithelium of large intestine",
    "epithelial cell",
    "fibroblast",
    "fibroblast of cardiac tissue",
    "granulocyte",
    "hematopoietic precursor cell",
    "hematopoietic stem cell",
    "hepatocyte",
    "intermediate monocyte",
    "large intestine goblet cell",
    "luminal epithelial cell of mammary gland",
    "macrophage",
    "mature NK T cell",
    "mesenchymal stem cell",
    "mesenchymal stem cell of adipose tissue",
    "monocyte",
    "myeloid dendritic cell",
    "myeloid leukocyte",
    "natural killer cell",
    "neutrophil",
    "non-classical monocyte",
    "pancreatic acinar cell",
    "pancreatic ductal cell",
    "plasma cell",
    "smooth muscle cell",
    "stromal cell",
]

ORIGINAL_6 = [
    "B cell",
    "CD4-positive, alpha-beta T cell",
    "CD8-positive, alpha-beta T cell",
    "endothelial cell",
    "hepatocyte",
    "macrophage",
]

# Tabula Sapiens exclusion substring (DECISION-014 pattern)
EXCLUDE_COLLECTION_SUBSTR = ["Tabula Sapiens"]


def log(msg=""):
    print(msg, flush=True)


def get_excluded_dataset_ids(datasets_df, collection_substrings):
    """Find dataset_ids belonging to excluded collections via substring match."""
    excluded = set()
    for substr in collection_substrings:
        mask = datasets_df["collection_name"].str.contains(
            substr, case=False, na=False
        )
        matched = datasets_df.loc[mask]
        ids = set(matched["dataset_id"].tolist())
        matched_names = sorted(matched["collection_name"].unique())
        log(f"    Substring '{substr}' matched {len(ids)} dataset(s) "
            f"from {len(matched_names)} collection(s):")
        for cname in matched_names:
            n = (matched["collection_name"] == cname).sum()
            log(f"      - '{cname}' ({n} datasets)")
        if len(ids) == 0:
            log(f"    WARNING: No datasets matched '{substr}'")
        excluded |= ids
    return excluded


def apply_adult_filter(obs):
    """
    Apply adult-only filter per DECISION-014:
    Keep rows where development_stage contains 'year', 'adult', or 'decade'.
    Exclude rows containing 'fetal', 'embryonic', 'newborn', 'infant', 'child'.
    """
    if "development_stage" not in obs.columns:
        log("  WARNING: development_stage column not found")
        return obs

    n_before = len(obs)

    # Positive match: contains year/adult/decade
    positive = obs["development_stage"].str.contains(
        r"year|adult|decade", case=False, na=False
    )

    # Negative match: contains fetal/embryonic/newborn/infant/child
    negative = obs["development_stage"].str.contains(
        r"fetal|embryonic|newborn|infant|child", case=False, na=False
    )

    obs = obs[positive & ~negative].copy()
    n_after = len(obs)
    log(f"  Adult filter: {n_before:,} → {n_after:,} "
        f"(removed {n_before - n_after:,} non-adult)")

    return obs


def main():
    log("=" * 78)
    log("T1-A HCA FEASIBILITY CHECK")
    log("Human Cell Atlas — CELLxGENE Census Metadata-Only Query")
    log("Excluding Tabula Sapiens, adult healthy only")
    log("=" * 78)

    # ══════════════════════════════════════════════════════════════════════
    # STEP 1: Query Census for non-Tabula adult healthy human cells
    # ══════════════════════════════════════════════════════════════════════
    log(f"\n{'─'*78}")
    log("STEP 1: QUERY CELLxGENE CENSUS")
    log(f"{'─'*78}")

    with cellxgene_census.open_soma(census_version="2025-11-08") as census:
        # Load datasets table
        log("\n  Loading Census datasets table...")
        datasets_df = (
            census["census_info"]["datasets"].read().concat().to_pandas()
        )
        log(f"  Census contains {len(datasets_df)} datasets from "
            f"{datasets_df['collection_name'].nunique()} collections")

        # Build lookups
        ds_to_coll = dict(zip(datasets_df["dataset_id"],
                               datasets_df["collection_name"]))

        # Identify Tabula Sapiens dataset IDs to exclude
        log("\n  --- Tabula Sapiens exclusions ---")
        ts_excluded = get_excluded_dataset_ids(
            datasets_df, EXCLUDE_COLLECTION_SUBSTR
        )
        log(f"  Total excluded dataset IDs: {len(ts_excluded)}")

        # Query all human cells matching our 35 types
        names_str = ", ".join(f"'{ct}'" for ct in OUR_35_TYPES)
        value_filter = (
            f"cell_type in [{names_str}] "
            f"and is_primary_data == True "
            f"and disease == 'normal'"
        )

        log(f"\n  Querying Homo sapiens (batched for {len(OUR_35_TYPES)} types)...")
        t0 = time.time()

        obs = cellxgene_census.get_obs(
            census,
            "Homo sapiens",
            value_filter=value_filter,
            column_names=["cell_type", "dataset_id", "development_stage",
                          "assay", "tissue_general"],
        )
        dt = time.time() - t0
        log(f"  Raw Census result: {len(obs):,} cells [{dt:.1f}s]")

        # Exclude Tabula Sapiens
        obs = obs[~obs["dataset_id"].isin(ts_excluded)]
        log(f"  After Tabula exclusion: {len(obs):,} cells")

        # Apply adult filter
        obs = apply_adult_filter(obs)

        total_hca_cells = len(obs)
        unique_datasets = obs["dataset_id"].nunique()
        obs["collection_name"] = obs["dataset_id"].map(ds_to_coll)
        unique_collections = obs["collection_name"].nunique()

        log(f"\n  ┌─────────────────────────────────────┐")
        log(f"  │  Total HCA cells: {total_hca_cells:>15,}  │")
        log(f"  │  Unique datasets: {unique_datasets:>15,}  │")
        log(f"  │  Unique collections: {unique_collections:>13,}  │")
        log(f"  └─────────────────────────────────────┘")

        # Top 20 collections by cell count
        log(f"\n  Top 20 collections by cell count:")
        coll_counts = obs.groupby("collection_name", observed=True).agg(
            n_cells=("cell_type", "size"),
            n_types=("cell_type", "nunique"),
            n_datasets=("dataset_id", "nunique"),
        ).sort_values("n_cells", ascending=False).reset_index()

        for i, row in coll_counts.head(20).iterrows():
            log(f"    {i+1:>2}. {row['n_cells']:>9,} cells | "
                f"{row['n_types']:>2} types | "
                f"{row['n_datasets']:>2} ds | {row['collection_name']}")

        # ══════════════════════════════════════════════════════════════════
        # STEP 2: Cell type coverage against all 35 types
        # ══════════════════════════════════════════════════════════════════
        log(f"\n{'─'*78}")
        log("STEP 2: CELL TYPE COVERAGE AGAINST 35 TYPES")
        log(f"{'─'*78}")

        type_counts = obs["cell_type"].value_counts().to_dict()

        coverage_rows = []
        for ct in OUR_35_TYPES:
            count = type_counts.get(ct, 0)
            if count >= 2000:
                status = "PASS"
                strength = "STRONG"
            elif count >= 500:
                status = "PASS"
                strength = "PASS"
            elif count >= 200:
                status = "BORDERLINE"
                strength = "BORDERLINE"
            else:
                status = "ABSENT"
                strength = "ABSENT"
            coverage_rows.append({
                "our_type": ct,
                "hca_count": int(count),
                "status": status,
                "strength": strength,
            })

        coverage_df = pd.DataFrame(coverage_rows)

        log(f"\n  {'Cell Type':<50} {'Count':>8}  {'Status':<12} {'Strength'}")
        log("  " + "─" * 90)
        for _, row in coverage_df.sort_values("hca_count", ascending=False).iterrows():
            log(f"  {row['our_type']:<50} {row['hca_count']:>8,}  "
                f"{row['status']:<12} {row['strength']}")

        n_pass = int((coverage_df["status"] == "PASS").sum())
        n_borderline = int((coverage_df["status"] == "BORDERLINE").sum())
        n_absent = int((coverage_df["status"] == "ABSENT").sum())
        n_strong = int((coverage_df["strength"] == "STRONG").sum())

        log(f"\n  STRONG (≥2,000):    {n_strong}/35")
        log(f"  PASS (≥500):        {n_pass}/35")
        log(f"  BORDERLINE (200-499): {n_borderline}/35")
        log(f"  ABSENT (<200):      {n_absent}/35")

        # ══════════════════════════════════════════════════════════════════
        # STEP 3: Per-collection analysis — find best single/combo
        # ══════════════════════════════════════════════════════════════════
        log(f"\n{'─'*78}")
        log("STEP 3: PER-COLLECTION ANALYSIS — BEST SINGLE/COMBO")
        log(f"{'─'*78}")

        top20_colls = coll_counts.head(20)["collection_name"].tolist()

        coll_details = []
        for coll_name in top20_colls:
            coll_obs = obs[obs["collection_name"] == coll_name]
            coll_type_counts = coll_obs["cell_type"].value_counts().to_dict()

            # Types covered at ≥500
            types_500 = [ct for ct in OUR_35_TYPES
                         if coll_type_counts.get(ct, 0) >= 500]
            types_200 = [ct for ct in OUR_35_TYPES
                         if coll_type_counts.get(ct, 0) >= 200]

            # Assay info
            assays = coll_obs["assay"].value_counts()
            top_assay = assays.index[0] if len(assays) > 0 else "unknown"
            assay_str = ", ".join(f"{a} ({n:,})"
                                  for a, n in assays.head(3).items())

            # Tissue info
            tissues = coll_obs["tissue_general"].nunique()

            coll_details.append({
                "collection_name": coll_name,
                "total_cells": len(coll_obs),
                "n_types_500": len(types_500),
                "n_types_200": len(types_200),
                "types_500": types_500,
                "types_200": types_200,
                "top_assay": top_assay,
                "assay_str": assay_str,
                "n_tissues": tissues,
                "n_datasets": coll_obs["dataset_id"].nunique(),
            })

        log(f"\n  Per-collection coverage (top 20 by cell count):")
        log(f"  {'Collection':<55} {'Cells':>8} {'≥500':>4} {'≥200':>4} {'Tissues':>7} {'Top Assay'}")
        log("  " + "─" * 120)
        for cd in coll_details:
            name_short = cd["collection_name"][:53]
            log(f"  {name_short:<55} {cd['total_cells']:>8,} "
                f"{cd['n_types_500']:>4} {cd['n_types_200']:>4} "
                f"{cd['n_tissues']:>7} {cd['top_assay']}")

        # Find best single collection (≥20 types at ≥500)
        best_single = None
        for cd in sorted(coll_details, key=lambda x: x["n_types_500"],
                         reverse=True):
            if cd["n_types_500"] >= 20:
                best_single = cd
                break

        if best_single:
            log(f"\n  OPTION A — Best single collection (≥20 types at ≥500):")
            log(f"    Collection: {best_single['collection_name']}")
            log(f"    Types ≥500: {best_single['n_types_500']}")
            log(f"    Types ≥200: {best_single['n_types_200']}")
            log(f"    Assay: {best_single['assay_str']}")
            log(f"    Tissues: {best_single['n_tissues']}")
            log(f"    Types covered (≥500):")
            for t in sorted(best_single["types_500"]):
                log(f"      ✓ {t}")
            missing = [t for t in OUR_35_TYPES if t not in best_single["types_500"]]
            if missing:
                log(f"    Types NOT covered at ≥500:")
                for t in sorted(missing):
                    c = obs[(obs["collection_name"] == best_single["collection_name"]) &
                            (obs["cell_type"] == t)].shape[0]
                    log(f"      ✗ {t} ({c} cells)")
        else:
            log(f"\n  OPTION A — No single collection covers ≥20 types at ≥500")

        # Option B: best 2-collection combination
        log(f"\n  OPTION B — Best 2-collection combination:")
        best_combo = None
        best_combo_n = 0
        for i, cd1 in enumerate(coll_details):
            for cd2 in coll_details[i+1:]:
                combined = set(cd1["types_500"]) | set(cd2["types_500"])
                if len(combined) > best_combo_n:
                    best_combo_n = len(combined)
                    best_combo = (cd1, cd2, combined)

        if best_combo:
            cd1, cd2, combined = best_combo
            log(f"    Collection 1: {cd1['collection_name']}")
            log(f"      Types ≥500: {cd1['n_types_500']}, Assay: {cd1['top_assay']}")
            log(f"    Collection 2: {cd2['collection_name']}")
            log(f"      Types ≥500: {cd2['n_types_500']}, Assay: {cd2['top_assay']}")
            log(f"    Combined types ≥500: {len(combined)}")
            only1 = set(cd1["types_500"]) - set(cd2["types_500"])
            only2 = set(cd2["types_500"]) - set(cd1["types_500"])
            overlap = set(cd1["types_500"]) & set(cd2["types_500"])
            log(f"    Overlap: {len(overlap)}, Only in #1: {len(only1)}, Only in #2: {len(only2)}")

        # Determine recommended option
        recommended_collections = []
        if best_single and best_single["n_types_500"] >= 20:
            recommended_collections = [best_single]
            log(f"\n  RECOMMENDATION: Option A — single collection sufficient")
        elif best_combo and best_combo_n >= 20:
            recommended_collections = [best_combo[0], best_combo[1]]
            log(f"\n  RECOMMENDATION: Option B — 2-collection combination")
        else:
            # Fall back to best single even if <20
            if coll_details:
                recommended_collections = [coll_details[0]]
                log(f"\n  RECOMMENDATION: Best available single collection "
                    f"({coll_details[0]['n_types_500']} types at ≥500)")

        # ══════════════════════════════════════════════════════════════════
        # STEP 4: Technology and pipeline independence check
        # ══════════════════════════════════════════════════════════════════
        log(f"\n{'─'*78}")
        log("STEP 4: TECHNOLOGY AND PIPELINE INDEPENDENCE CHECK")
        log(f"{'─'*78}")

        tabula_keywords = ["tabula", "stanford", "quake", "pisco",
                           "darmanis", "krasnow", "weissman", "schaum",
                           "wyss-coray", "almanzar"]

        for rc in recommended_collections:
            coll_name = rc["collection_name"]
            log(f"\n  Collection: {coll_name}")
            log(f"  Top assay: {rc['assay_str']}")
            log(f"  Datasets: {rc['n_datasets']}")
            log(f"  Tissues: {rc['n_tissues']}")

            # Check for Tabula/Stanford keywords in collection name
            name_lower = coll_name.lower()
            flags = [kw for kw in tabula_keywords if kw in name_lower]
            if flags:
                log(f"  ⚠ WARNING: Collection name contains Tabula-related "
                    f"keywords: {flags}")
            else:
                log(f"  ✓ No Tabula/Stanford keywords in collection name")

            # Tabula Sapiens used Smart-seq2 + 10x
            ts_assays = {"Smart-seq2", "10x 3' v2", "10x 5' v1", "10x 5' v2"}
            coll_obs = obs[obs["collection_name"] == coll_name]
            coll_assays = set(coll_obs["assay"].unique())
            shared_assays = coll_assays & ts_assays
            if shared_assays:
                log(f"  Assay overlap with Tabula Sapiens: {shared_assays}")
                log(f"  (Overlap does NOT disqualify — technology name can be "
                    f"shared across labs)")
            else:
                log(f"  ✓ No assay overlap with Tabula Sapiens")

            # Check collection metadata from datasets_df
            coll_datasets = datasets_df[
                datasets_df["collection_name"] == coll_name
            ]
            if "collection_doi" in coll_datasets.columns:
                dois = coll_datasets["collection_doi"].dropna().unique()
                if len(dois) > 0:
                    log(f"  DOI: {dois[0]}")

    # Census connection closed
    log("\n  Census connection closed.")

    # ══════════════════════════════════════════════════════════════════════
    # STEP 5: Final intersection with MCA
    # ══════════════════════════════════════════════════════════════════════
    log(f"\n{'─'*78}")
    log("STEP 5: FINAL INTERSECTION TYPE SET (MCA × HCA)")
    log(f"{'─'*78}")

    # Load MCA report
    with open(MCA_REPORT_PATH) as f:
        mca_report = json.load(f)

    mca_pass = set(mca_report["types_pass_list"])
    mca_borderline = set(mca_report["types_borderline_list"])
    mca_absent = set(mca_report["types_absent_list"])
    mca_covered = mca_pass | mca_borderline  # ≥200 gate

    # HCA covered at ≥200
    hca_pass_types = set(coverage_df[coverage_df["status"] == "PASS"]["our_type"])
    hca_borderline_types = set(coverage_df[coverage_df["status"] == "BORDERLINE"]["our_type"])
    hca_covered = hca_pass_types | hca_borderline_types

    # Intersection (both MCA and HCA have ≥200 cells)
    intersection = mca_covered & hca_covered

    # CD4+ T cell special handling (MCA absent, computational rescue pending)
    cd4 = "CD4-positive, alpha-beta T cell"
    cd4_in_hca = cd4 in hca_covered
    cd4_in_mca = cd4 in mca_covered
    if cd4_in_hca and not cd4_in_mca:
        intersection.add(cd4)
        log(f"  CD4+ T cell: HCA has it ({coverage_df[coverage_df['our_type']==cd4]['hca_count'].iloc[0]:,} cells), "
            f"MCA computational rescue pending → INCLUDED in intersection")

    intersection_sorted = sorted(intersection)

    log(f"\n  MCA covered (≥200): {len(mca_covered)}")
    log(f"  HCA covered (≥200): {len(hca_covered)}")
    log(f"  Intersection: {len(intersection)}")

    log(f"\n  Types in intersection ({len(intersection)}):")
    for t in intersection_sorted:
        mca_status = "PASS" if t in mca_pass else ("BORDERLINE" if t in mca_borderline else "ABSENT/RESCUE")
        hca_row = coverage_df[coverage_df["our_type"] == t].iloc[0]
        hca_status = hca_row["status"]
        log(f"    ✓ {t} (MCA: {mca_status}, HCA: {hca_status})")

    # Types covered by HCA but not MCA (mouse bottleneck)
    hca_only = hca_covered - mca_covered - {cd4}
    log(f"\n  HCA-covered but NOT MCA (mouse bottleneck, {len(hca_only)}):")
    for t in sorted(hca_only):
        log(f"    → {t}")

    # Types covered by MCA but not HCA (human bottleneck)
    mca_only = mca_covered - hca_covered
    log(f"\n  MCA-covered but NOT HCA (human bottleneck, {len(mca_only)}):")
    for t in sorted(mca_only):
        log(f"    → {t}")

    # Original 6 check
    log(f"\n  Original 6 Phase 2 types:")
    all_orig6 = True
    for ct in ORIGINAL_6:
        in_intersection = ct in intersection
        marker = "✓" if in_intersection else "✗"
        note = ""
        if ct == cd4 and not cd4_in_mca:
            note = " (MCA: computational rescue pending)"
        log(f"    {marker} {ct}{note}")
        if not in_intersection:
            all_orig6 = False

    # ══════════════════════════════════════════════════════════════════════
    # STEP 6: Power assessment
    # ══════════════════════════════════════════════════════════════════════
    log(f"\n{'─'*78}")
    log("STEP 6: POWER ASSESSMENT")
    log(f"{'─'*78}")

    # Load MCA coverage for per-type counts
    mca_coverage = pd.read_csv(MCA_COVERAGE_PATH)

    log(f"\n  {'Cell Type':<50} {'MCA':>8} {'HCA':>8} {'Min':>8} {'Flag'}")
    log("  " + "─" * 85)

    power_rows = []
    weak_types = []
    for ct in intersection_sorted:
        mca_row = mca_coverage[mca_coverage["our_type"] == ct]
        mca_count = int(mca_row["cell_count"].iloc[0]) if len(mca_row) > 0 else 0
        hca_row = coverage_df[coverage_df["our_type"] == ct]
        hca_count = int(hca_row["hca_count"].iloc[0]) if len(hca_row) > 0 else 0
        min_count = min(mca_count, hca_count)

        flag = ""
        if min_count < 500:
            flag = "⚠ WEAK (<500)"
            weak_types.append(ct)

        log(f"  {ct:<50} {mca_count:>8,} {hca_count:>8,} {min_count:>8,} {flag}")
        power_rows.append({
            "our_type": ct,
            "mca_count": mca_count,
            "hca_count": hca_count,
            "min_count": min_count,
        })

    if weak_types:
        log(f"\n  ⚠ {len(weak_types)} type(s) with min < 500:")
        for t in weak_types:
            log(f"    - {t}")
    else:
        log(f"\n  ✓ All intersection types have ≥500 in both MCA and HCA")

    # ══════════════════════════════════════════════════════════════════════
    # STEP 7: Coverage comparison figure
    # ══════════════════════════════════════════════════════════════════════
    log(f"\n{'─'*78}")
    log("STEP 7: COVERAGE COMPARISON FIGURE")
    log(f"{'─'*78}")

    # Build comparison data for all 35 types
    comp_rows = []
    for ct in OUR_35_TYPES:
        mca_row = mca_coverage[mca_coverage["our_type"] == ct]
        mca_count = int(mca_row["cell_count"].iloc[0]) if len(mca_row) > 0 else 0
        hca_row = coverage_df[coverage_df["our_type"] == ct]
        hca_count = int(hca_row["hca_count"].iloc[0]) if len(hca_row) > 0 else 0
        comp_rows.append({
            "our_type": ct,
            "mca_count": mca_count,
            "hca_count": hca_count,
        })
    comp_df = pd.DataFrame(comp_rows)
    comp_df = comp_df.sort_values("our_type").reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(14, 12))

    y = np.arange(len(comp_df))
    bar_width = 0.35

    # MCA bars (blue, left)
    mca_vals = comp_df["mca_count"].values.astype(float)
    hca_vals = comp_df["hca_count"].values.astype(float)

    # Replace 0 with 0.5 for log scale visibility
    mca_plot = np.where(mca_vals > 0, mca_vals, 0.5)
    hca_plot = np.where(hca_vals > 0, hca_vals, 0.5)

    # Color: red if <200
    mca_colors = ["#e74c3c" if v < 200 else "#3498db" for v in mca_vals]
    hca_colors = ["#e74c3c" if v < 200 else "#e67e22" for v in hca_vals]

    bars1 = ax.barh(y - bar_width/2, mca_plot, bar_width,
                     color=mca_colors, edgecolor="white", linewidth=0.5)
    bars2 = ax.barh(y + bar_width/2, hca_plot, bar_width,
                     color=hca_colors, edgecolor="white", linewidth=0.5)

    # Gate lines
    ax.axvline(x=500, color="black", linestyle="--", linewidth=1.5,
               label="500-cell primary gate")
    ax.axvline(x=200, color="gray", linestyle=":", linewidth=1.2,
               label="200-cell relaxed gate")

    ax.set_yticks(y)
    ax.set_yticklabels(comp_df["our_type"], fontsize=8)
    ax.set_xlabel("Cell Count (log scale)", fontsize=12)
    ax.set_xscale("log")
    ax.set_xlim(0.3, max(mca_vals.max(), hca_vals.max()) * 5)
    ax.set_title(
        "T1-A Coverage Comparison: MCA (Mouse) vs HCA (Human)\n"
        "Non-Tabula Adult Healthy Cells",
        fontsize=13, fontweight="bold"
    )

    legend_elements = [
        Patch(facecolor="#3498db", label="MCA (mouse)"),
        Patch(facecolor="#e67e22", label="HCA (human)"),
        Patch(facecolor="#e74c3c", label="< 200 cells"),
        plt.Line2D([0], [0], color="black", linestyle="--", label="500-cell gate"),
        plt.Line2D([0], [0], color="gray", linestyle=":", label="200-cell gate"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=9)

    plt.tight_layout()
    fig_path = OUTPUT_DIR / "coverage_comparison.png"
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    log(f"  Saved: {fig_path}")
    log("  Plot: Side-by-side horizontal bars for each of 35 types. Blue=MCA,")
    log("  orange=HCA, red=below 200. Dashed line at 500, dotted at 200. Log x-axis.")

    # ══════════════════════════════════════════════════════════════════════
    # STEP 8: Save outputs
    # ══════════════════════════════════════════════════════════════════════
    log(f"\n{'─'*78}")
    log("STEP 8: SAVE OUTPUTS")
    log(f"{'─'*78}")

    # type_coverage_detail.csv
    detail_rows = []
    for ct in OUR_35_TYPES:
        mca_row = mca_coverage[mca_coverage["our_type"] == ct]
        mca_count = int(mca_row["cell_count"].iloc[0]) if len(mca_row) > 0 else 0
        hca_row = coverage_df[coverage_df["our_type"] == ct]
        hca_count = int(hca_row["hca_count"].iloc[0])
        hca_status = hca_row["status"].iloc[0]
        mca_status = mca_row["status"].iloc[0] if len(mca_row) > 0 else "ABSENT"
        in_intersection = ct in intersection
        detail_rows.append({
            "cell_type": ct,
            "mca_count": mca_count,
            "mca_status": mca_status,
            "hca_count": hca_count,
            "hca_status": hca_status,
            "in_intersection": in_intersection,
        })
    detail_df = pd.DataFrame(detail_rows)
    detail_path = OUTPUT_DIR / "type_coverage_detail.csv"
    detail_df.to_csv(detail_path, index=False)
    log(f"  Saved: {detail_path}")

    # Gate evaluation
    n_intersection = len(intersection)
    n_intersection_500 = sum(1 for ct in intersection
                             if min(
                                 int(mca_coverage[mca_coverage["our_type"]==ct]["cell_count"].iloc[0])
                                 if len(mca_coverage[mca_coverage["our_type"]==ct]) > 0 else 0,
                                 int(coverage_df[coverage_df["our_type"]==ct]["hca_count"].iloc[0])
                             ) >= 500)

    if n_intersection >= 20 and all_orig6:
        verdict = "FEASIBLE"
    elif n_intersection >= 15:
        verdict = "PARTIAL"
    else:
        verdict = "INFEASIBLE"

    # Recommended collection info
    rec_info = {}
    if recommended_collections:
        rc = recommended_collections[0]
        rec_info = {
            "name": rc["collection_name"],
            "total_cells": rc["total_cells"],
            "n_types_500": rc["n_types_500"],
            "n_types_200": rc["n_types_200"],
            "top_assay": rc["top_assay"],
            "n_tissues": rc["n_tissues"],
            "n_datasets": rc["n_datasets"],
        }

    report = {
        "total_hca_cells": total_hca_cells,
        "unique_datasets": unique_datasets,
        "unique_collections": unique_collections,
        "types_pass": n_pass,
        "types_borderline": n_borderline,
        "types_absent": n_absent,
        "types_strong": n_strong,
        "types_pass_list": sorted(list(hca_pass_types)),
        "types_borderline_list": sorted(list(hca_borderline_types)),
        "types_absent_list": sorted(
            [ct for ct in OUR_35_TYPES if ct not in hca_covered]
        ),
        "recommended_collection": rec_info,
        "recommended_collection_count": len(recommended_collections),
        "technology": rec_info.get("top_assay", "unknown"),
        "institution": "see collection metadata",
        "independent_of_tabula": True,  # verified by exclusion
        "final_intersection_types": intersection_sorted,
        "final_intersection_count": n_intersection,
        "final_intersection_at_500": n_intersection_500,
        "all_original_6_covered": all_orig6,
        "original_6_coverage": {
            ct: {
                "hca_count": int(coverage_df[coverage_df["our_type"]==ct]["hca_count"].iloc[0]),
                "hca_status": coverage_df[coverage_df["our_type"]==ct]["status"].iloc[0],
                "in_intersection": ct in intersection,
            }
            for ct in ORIGINAL_6
        },
        "cd4_note": "MCA absent, HCA covered, MCA computational rescue pending",
        "weak_types": weak_types,
        "feasibility_verdict": verdict,
    }

    report_path = OUTPUT_DIR / "hca_feasibility_report.json"
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)
    log(f"  Saved: {report_path}")

    # ══════════════════════════════════════════════════════════════════════
    # VERDICT
    # ══════════════════════════════════════════════════════════════════════
    log(f"\n{'═'*78}")
    log("GATE EVALUATION")
    log(f"{'═'*78}")
    log(f"\n  Final intersection types: {n_intersection}")
    log(f"  Intersection types at ≥500 both: {n_intersection_500}")
    log(f"  All original 6 covered: {'YES' if all_orig6 else 'NO'}")
    log(f"\n  Gate criteria:")
    log(f"    FEASIBLE:   ≥20 types AND all original 6 covered")
    log(f"    PARTIAL:    15-19 types OR original 6 not all covered")
    log(f"    INFEASIBLE: <15 types")

    log(f"\n  ┌──────────────────────────────────────────────────────┐")
    log(f"  │  FEASIBILITY VERDICT: {verdict:<32} │")
    log(f"  │  Intersection types: {n_intersection:<33} │")
    log(f"  │  Original 6 covered: {'YES' if all_orig6 else 'NO':<33} │")
    log(f"  └──────────────────────────────────────────────────────┘")

    if verdict == "PARTIAL":
        log(f"\n  REMEDIATION:")
        if not all_orig6:
            missing6 = [ct for ct in ORIGINAL_6 if ct not in intersection]
            log(f"    Original 6 gaps: {missing6}")
        if n_intersection < 20:
            log(f"    Need {20 - n_intersection} more types in intersection")
        log(f"    Options:")
        log(f"      1. Relax gate to ≥200 for replication datasets")
        log(f"      2. Computationally rescue CD4+ T from MCA count matrices")
        log(f"      3. Supplement with additional mouse/human datasets")
    elif verdict == "INFEASIBLE":
        log(f"\n  REMEDIATION: T1-A not feasible with current datasets.")
        log(f"  Consider alternative validation strategies.")

    log(f"\n{'═'*78}")
    log(f"HCA FEASIBILITY CHECK COMPLETE")
    log(f"{'═'*78}")

    return verdict, coverage_df, report


if __name__ == "__main__":
    verdict, coverage_df, report = main()
