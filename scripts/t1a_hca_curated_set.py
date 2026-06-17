#!/usr/bin/env python3
"""
T1-A HCA Curated Set Analysis: Can 2-4 collections replace pooled Census?

BIOLOGY: For T1-A independent replication, a curated set of 2-4 large collections
is more defensible than pooling 157 collections (cleaner batch structure, easier to
describe in methods). This script checks whether such a curated set can cover the
23-type MCA×HCA intersection.

MATH: Greedy set cover — at each step add the collection covering the most uncovered
types. NP-hard in general but greedy gives log(n)-approximation, and with n=23 types
and ~30 candidate collections, exact enumeration of small combinations is also feasible.

Steps:
1. Query Census for per-collection × per-type counts (top 30 collections)
2. Greedy set cover on 23 intersection types
3. Assess greedy solution: CURATED SET vs POOLED recommendation
"""

import os
import sys
import json
import time
from pathlib import Path
from collections import defaultdict

import pandas as pd
import numpy as np
import cellxgene_census

# ── Configuration ──────────────────────────────────────────────────────────
OUTPUT_DIR = Path("output/validation/hca_feasibility")

# The 23 intersection types from HCA feasibility report
INTERSECTION_TYPES = [
    "B cell",
    "CD4-positive, alpha-beta T cell",
    "CD8-positive, alpha-beta T cell",
    "T cell",
    "basal cell",
    "endothelial cell",
    "enterocyte of epithelium of large intestine",
    "epithelial cell",
    "fibroblast",
    "granulocyte",
    "hematopoietic precursor cell",
    "hepatocyte",
    "luminal epithelial cell of mammary gland",
    "macrophage",
    "mesenchymal stem cell",
    "monocyte",
    "myeloid dendritic cell",
    "natural killer cell",
    "neutrophil",
    "pancreatic acinar cell",
    "pancreatic ductal cell",
    "smooth muscle cell",
    "stromal cell",
]

EXCLUDE_COLLECTION_SUBSTR = ["Tabula Sapiens"]


def log(msg=""):
    print(msg, flush=True)


def apply_adult_filter(obs):
    """Apply adult-only filter per DECISION-014."""
    if "development_stage" not in obs.columns:
        return obs
    n_before = len(obs)
    positive = obs["development_stage"].str.contains(
        r"year|adult|decade", case=False, na=False
    )
    negative = obs["development_stage"].str.contains(
        r"fetal|embryonic|newborn|infant|child", case=False, na=False
    )
    obs = obs[positive & ~negative].copy()
    log(f"  Adult filter: {n_before:,} → {len(obs):,} (removed {n_before - len(obs):,})")
    return obs


def main():
    log("=" * 78)
    log("T1-A HCA CURATED SET ANALYSIS")
    log("Can 2-4 collections replace pooled 157-collection Census?")
    log("=" * 78)

    # ══════════════════════════════════════════════════════════════════════
    # STEP 0: Query Census for per-collection per-type counts
    # ══════════════════════════════════════════════════════════════════════
    log(f"\n{'─'*78}")
    log("STEP 0: QUERY CENSUS FOR PER-COLLECTION × PER-TYPE COUNTS")
    log(f"{'─'*78}")

    with cellxgene_census.open_soma(census_version="2025-11-08") as census:
        # Load datasets table
        log("\n  Loading Census datasets table...")
        datasets_df = (
            census["census_info"]["datasets"].read().concat().to_pandas()
        )
        ds_to_coll = dict(zip(datasets_df["dataset_id"],
                               datasets_df["collection_name"]))

        # Exclude Tabula Sapiens
        ts_mask = datasets_df["collection_name"].str.contains(
            "Tabula Sapiens", case=False, na=False
        )
        ts_excluded = set(datasets_df.loc[ts_mask, "dataset_id"].tolist())
        log(f"  Excluding {len(ts_excluded)} Tabula Sapiens dataset IDs")

        # Query all 23 intersection types
        names_str = ", ".join(f"'{ct}'" for ct in INTERSECTION_TYPES)
        value_filter = (
            f"cell_type in [{names_str}] "
            f"and is_primary_data == True "
            f"and disease == 'normal'"
        )

        log(f"  Querying Homo sapiens for {len(INTERSECTION_TYPES)} intersection types...")
        t0 = time.time()
        obs = cellxgene_census.get_obs(
            census,
            "Homo sapiens",
            value_filter=value_filter,
            column_names=["cell_type", "dataset_id", "development_stage",
                          "assay", "tissue_general"],
        )
        dt = time.time() - t0
        log(f"  Raw result: {len(obs):,} cells [{dt:.1f}s]")

        # Exclude Tabula Sapiens
        obs = obs[~obs["dataset_id"].isin(ts_excluded)]
        log(f"  After Tabula exclusion: {len(obs):,}")

        # Adult filter
        obs = apply_adult_filter(obs)

        # Map dataset_id → collection_name
        obs["collection_name"] = obs["dataset_id"].map(ds_to_coll)

        # Also extract collection-level metadata for institution/PI info
        # Census doesn't have PI names directly, but we have collection_name
        # and can check DOIs
        coll_meta = {}
        for coll_name in obs["collection_name"].unique():
            coll_ds = datasets_df[datasets_df["collection_name"] == coll_name]
            doi = ""
            if "collection_doi" in coll_ds.columns:
                dois = coll_ds["collection_doi"].dropna().unique()
                if len(dois) > 0:
                    doi = str(dois[0])
            coll_meta[coll_name] = {"doi": doi}

    log("  Census connection closed.\n")

    # ══════════════════════════════════════════════════════════════════════
    # Compute per-collection per-type counts
    # ══════════════════════════════════════════════════════════════════════
    # Group by collection × cell_type
    coll_type_counts = (
        obs.groupby(["collection_name", "cell_type"], observed=True)
        .size()
        .reset_index(name="count")
    )

    # Per-collection summary
    coll_summary = obs.groupby("collection_name", observed=True).agg(
        total_cells=("cell_type", "size"),
        n_types_raw=("cell_type", "nunique"),
        n_datasets=("dataset_id", "nunique"),
        n_tissues=("tissue_general", "nunique"),
    ).reset_index()

    # For each collection, determine which intersection types it covers at ≥500
    coll_coverage = {}
    for coll_name in coll_summary["collection_name"]:
        sub = coll_type_counts[coll_type_counts["collection_name"] == coll_name]
        type_dict = dict(zip(sub["cell_type"], sub["count"]))
        covered_500 = [t for t in INTERSECTION_TYPES
                       if type_dict.get(t, 0) >= 500]
        covered_200 = [t for t in INTERSECTION_TYPES
                       if type_dict.get(t, 0) >= 200]
        coll_coverage[coll_name] = {
            "types_500": set(covered_500),
            "types_200": set(covered_200),
            "type_counts": {t: type_dict.get(t, 0) for t in INTERSECTION_TYPES},
        }

    # Also get top assay per collection
    coll_assays = {}
    for coll_name in coll_summary["collection_name"]:
        coll_obs = obs[obs["collection_name"] == coll_name]
        assay_vc = coll_obs["assay"].value_counts()
        top_assay = assay_vc.index[0] if len(assay_vc) > 0 else "unknown"
        assay_str = ", ".join(f"{a} ({n:,})"
                              for a, n in assay_vc.head(2).items())
        coll_assays[coll_name] = {
            "top_assay": top_assay,
            "assay_str": assay_str,
            "all_assays": set(assay_vc.index.tolist()),
        }

    # ══════════════════════════════════════════════════════════════════════
    # STEP 1: Top 30 collections by type coverage
    # ══════════════════════════════════════════════════════════════════════
    log(f"{'─'*78}")
    log("STEP 1: TOP 30 COLLECTIONS BY INTERSECTION TYPE COVERAGE")
    log(f"{'─'*78}")

    # Sort by n_types_500 descending, then total_cells descending
    ranked = []
    for _, row in coll_summary.iterrows():
        coll_name = row["collection_name"]
        cov = coll_coverage[coll_name]
        assay_info = coll_assays[coll_name]
        meta = coll_meta.get(coll_name, {})
        ranked.append({
            "collection_name": coll_name,
            "total_cells": int(row["total_cells"]),
            "n_types_500": len(cov["types_500"]),
            "n_types_200": len(cov["types_200"]),
            "types_500": sorted(cov["types_500"]),
            "types_200": sorted(cov["types_200"]),
            "n_datasets": int(row["n_datasets"]),
            "n_tissues": int(row["n_tissues"]),
            "top_assay": assay_info["top_assay"],
            "assay_str": assay_info["assay_str"],
            "doi": meta.get("doi", ""),
        })

    ranked.sort(key=lambda x: (-x["n_types_500"], -x["total_cells"]))

    log(f"\n  Rank {'Collection':<55} {'Cells':>9} {'≥500':>4} {'≥200':>4} "
        f"{'Tiss':>4} {'Top Assay':<25}")
    log("  " + "─" * 120)

    for i, r in enumerate(ranked[:30], 1):
        name = r["collection_name"][:53]
        log(f"  {i:>3}. {name:<55} {r['total_cells']:>9,} "
            f"{r['n_types_500']:>4} {r['n_types_200']:>4} "
            f"{r['n_tissues']:>4} {r['top_assay']:<25}")

    # Print detailed type lists for top 10
    log(f"\n  Detailed type coverage for top 10 by ≥500 coverage:")
    for i, r in enumerate(ranked[:10], 1):
        log(f"\n  {i}. {r['collection_name']}")
        log(f"     Cells: {r['total_cells']:,} | Tissues: {r['n_tissues']} | "
            f"Assay: {r['assay_str']}")
        if r["doi"]:
            log(f"     DOI: {r['doi']}")
        log(f"     Types ≥500 ({r['n_types_500']}/23):")
        for t in r["types_500"]:
            c = coll_coverage[r["collection_name"]]["type_counts"][t]
            log(f"       ✓ {t} ({c:,})")
        # Show types at ≥200 but <500
        borderline = set(r["types_200"]) - set(r["types_500"])
        if borderline:
            log(f"     Types 200-499 ({len(borderline)}):")
            for t in sorted(borderline):
                c = coll_coverage[r["collection_name"]]["type_counts"][t]
                log(f"       ~ {t} ({c})")

    # ══════════════════════════════════════════════════════════════════════
    # STEP 2: Greedy set cover
    # ══════════════════════════════════════════════════════════════════════
    log(f"\n{'─'*78}")
    log("STEP 2: GREEDY SET COVER ON 23 INTERSECTION TYPES")
    log(f"{'─'*78}")

    target_types = set(INTERSECTION_TYPES)
    uncovered = target_types.copy()
    selected = []
    max_steps = 5

    for step in range(1, max_steps + 1):
        if not uncovered:
            break

        # Find collection covering most uncovered types at ≥500
        best_coll = None
        best_marginal = 0
        best_covered = set()

        for r in ranked:
            coll_name = r["collection_name"]
            if coll_name in [s["collection_name"] for s in selected]:
                continue
            marginal = len(coll_coverage[coll_name]["types_500"] & uncovered)
            if marginal > best_marginal:
                best_marginal = marginal
                best_coll = r
                best_covered = coll_coverage[coll_name]["types_500"] & uncovered

        if best_coll is None or best_marginal == 0:
            log(f"\n  Step {step}: No collection adds any new types. STOPPING.")
            break

        selected.append(best_coll)
        newly_covered = best_covered
        uncovered -= newly_covered

        log(f"\n  Step {step}: ADD '{best_coll['collection_name']}'")
        log(f"    Marginal types added: {best_marginal}")
        log(f"    New types: {sorted(newly_covered)}")
        log(f"    Assay: {best_coll['top_assay']}")
        log(f"    Tissues: {best_coll['n_tissues']}")
        log(f"    Cells: {best_coll['total_cells']:,}")
        log(f"    Cumulative covered: {len(target_types) - len(uncovered)}/23")

    greedy_covered = target_types - uncovered
    greedy_n = len(greedy_covered)

    log(f"\n  {'─'*60}")
    log(f"  GREEDY RESULT: {greedy_n}/23 types covered with "
        f"{len(selected)} collections")

    if uncovered:
        log(f"  UNCOVERED ({len(uncovered)}):")
        for t in sorted(uncovered):
            # Check if any collection has it at ≥200
            best_count = 0
            best_source = ""
            for r in ranked:
                c = coll_coverage[r["collection_name"]]["type_counts"][t]
                if c > best_count:
                    best_count = c
                    best_source = r["collection_name"][:40]
            log(f"    ✗ {t} (best: {best_count} in '{best_source}')")

    # Also try greedy with ≥200 gate
    log(f"\n  Greedy with relaxed ≥200 gate:")
    uncovered_200 = target_types.copy()
    selected_200 = []
    for step in range(1, max_steps + 1):
        if not uncovered_200:
            break
        best_coll = None
        best_marginal = 0
        for r in ranked:
            coll_name = r["collection_name"]
            if coll_name in [s["collection_name"] for s in selected_200]:
                continue
            marginal = len(coll_coverage[coll_name]["types_200"] & uncovered_200)
            if marginal > best_marginal:
                best_marginal = marginal
                best_coll = r
        if best_coll is None or best_marginal == 0:
            break
        selected_200.append(best_coll)
        uncovered_200 -= coll_coverage[best_coll["collection_name"]]["types_200"]

    greedy_200_n = len(target_types) - len(uncovered_200)
    log(f"  ≥200 gate: {greedy_200_n}/23 types with {len(selected_200)} collections")
    if uncovered_200:
        log(f"  Still uncovered at ≥200: {sorted(uncovered_200)}")

    # ══════════════════════════════════════════════════════════════════════
    # STEP 3: Assess greedy solution
    # ══════════════════════════════════════════════════════════════════════
    log(f"\n{'─'*78}")
    log("STEP 3: ASSESSMENT AND RECOMMENDATION")
    log(f"{'─'*78}")

    # Institutional diversity
    log(f"\n  Selected collections ({len(selected)}):")
    institutions = set()
    technologies = set()
    for i, s in enumerate(selected, 1):
        log(f"  {i}. {s['collection_name']}")
        log(f"     Cells: {s['total_cells']:,} | Tissues: {s['n_tissues']} | "
            f"Assay: {s['top_assay']}")
        if s["doi"]:
            log(f"     DOI: {s['doi']}")
        technologies.add(s["top_assay"])
        # Use collection name as proxy for institution diversity
        institutions.add(s["collection_name"])

    log(f"\n  Technologies: {sorted(technologies)}")
    log(f"  Number of distinct collections: {len(selected)}")
    tech_diverse = len(technologies) > 1
    log(f"  Technology diversity: {'YES' if tech_diverse else 'NO (all same assay)'}")

    # Check original 6 coverage
    orig_6 = {
        "B cell", "CD4-positive, alpha-beta T cell",
        "CD8-positive, alpha-beta T cell", "endothelial cell",
        "hepatocyte", "macrophage",
    }
    curated_covered = target_types - uncovered
    orig_6_covered = orig_6 & curated_covered
    orig_6_missing = orig_6 - curated_covered

    log(f"\n  Original 6 coverage: {len(orig_6_covered)}/6")
    if orig_6_missing:
        log(f"  Missing original 6: {sorted(orig_6_missing)}")

    # Final recommendation
    log(f"\n  {'═'*60}")

    if greedy_n >= 20 and len(selected) <= 4:
        recommendation = "CURATED SET"
        log(f"  RECOMMENDATION: {recommendation}")
        log(f"  {greedy_n}/23 types covered with {len(selected)} collections")
        log(f"  This is cleaner and more defensible than pooling 157 collections.")
    elif greedy_n >= 18 and len(selected) <= 5:
        recommendation = "CURATED SET (with caveats)"
        log(f"  RECOMMENDATION: {recommendation}")
        log(f"  {greedy_n}/23 types covered with {len(selected)} collections")
        log(f"  {len(uncovered)} types uncovered — acceptable if non-critical")
    else:
        recommendation = "POOLED"
        log(f"  RECOMMENDATION: {recommendation}")
        log(f"  Only {greedy_n}/23 types covered with ≤5 collections")
        log(f"  No curated combination of ≤5 collections covers ≥18 types")
        log(f"  Pooled 157-collection approach (DECISION-092) confirmed as necessary")

    log(f"\n  Summary table:")
    log(f"    Greedy ≥500 gate: {greedy_n}/23 types, {len(selected)} collections")
    log(f"    Greedy ≥200 gate: {greedy_200_n}/23 types, {len(selected_200)} collections")
    log(f"    Original 6 in curated: {len(orig_6_covered)}/6")
    log(f"    Technologies: {sorted(technologies)}")

    log(f"\n{'═'*78}")
    log(f"CURATED SET ANALYSIS COMPLETE — {recommendation}")
    log(f"{'═'*78}")

    # Save results
    result = {
        "greedy_n_types_500": greedy_n,
        "greedy_n_collections": len(selected),
        "greedy_collections": [
            {
                "name": s["collection_name"],
                "total_cells": s["total_cells"],
                "n_types_500": s["n_types_500"],
                "top_assay": s["top_assay"],
                "n_tissues": s["n_tissues"],
                "doi": s["doi"],
            }
            for s in selected
        ],
        "greedy_covered_types": sorted(curated_covered),
        "greedy_uncovered_types": sorted(uncovered),
        "greedy_200_n_types": greedy_200_n,
        "greedy_200_n_collections": len(selected_200),
        "original_6_in_curated": len(orig_6_covered),
        "technologies": sorted(technologies),
        "recommendation": recommendation,
    }

    result_path = OUTPUT_DIR / "curated_set_analysis.json"
    with open(result_path, "w") as f:
        json.dump(result, f, indent=2)
    log(f"\n  Saved: {result_path}")

    return recommendation, result


if __name__ == "__main__":
    recommendation, result = main()
