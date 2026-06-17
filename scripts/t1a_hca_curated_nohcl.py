#!/usr/bin/env python3
"""
Quick sensitivity: greedy set cover WITHOUT the HCL collection (same Guo lab as MCA).

If MCA is the mouse replication dataset and HCL is the human anchor, both come from
Guoji Guo's lab at Zhejiang University using microwell-seq. A shared-lab confound
undermines the independence claim for T1-A replication.

This script re-runs greedy set cover on the curated_set_analysis.json results,
excluding HCL, to see if a lab-independent curated set is feasible.
"""

import json, sys
from pathlib import Path

import cellxgene_census
import pandas as pd
import time

OUTPUT_DIR = Path("output/validation/hca_feasibility")

INTERSECTION_TYPES = [
    "B cell", "CD4-positive, alpha-beta T cell",
    "CD8-positive, alpha-beta T cell", "T cell", "basal cell",
    "endothelial cell", "enterocyte of epithelium of large intestine",
    "epithelial cell", "fibroblast", "granulocyte",
    "hematopoietic precursor cell", "hepatocyte",
    "luminal epithelial cell of mammary gland", "macrophage",
    "mesenchymal stem cell", "monocyte", "myeloid dendritic cell",
    "natural killer cell", "neutrophil", "pancreatic acinar cell",
    "pancreatic ductal cell", "smooth muscle cell", "stromal cell",
]

EXCLUDE_COLLECTION_SUBSTR = ["Tabula Sapiens"]
HCL_NAME = "Construction of a human cell landscape at single-cell level"


def log(msg=""):
    print(msg, flush=True)


def main():
    log("=" * 78)
    log("SENSITIVITY: GREEDY SET COVER WITHOUT HCL (shared Guo lab)")
    log("=" * 78)

    # Re-query Census (same as before but track more collections)
    log("\n  Querying Census...")

    with cellxgene_census.open_soma(census_version="2025-11-08") as census:
        datasets_df = (
            census["census_info"]["datasets"].read().concat().to_pandas()
        )
        ds_to_coll = dict(zip(datasets_df["dataset_id"],
                               datasets_df["collection_name"]))
        ts_mask = datasets_df["collection_name"].str.contains(
            "Tabula Sapiens", case=False, na=False
        )
        ts_excluded = set(datasets_df.loc[ts_mask, "dataset_id"].tolist())

        names_str = ", ".join(f"'{ct}'" for ct in INTERSECTION_TYPES)
        value_filter = (
            f"cell_type in [{names_str}] "
            f"and is_primary_data == True "
            f"and disease == 'normal'"
        )

        t0 = time.time()
        obs = cellxgene_census.get_obs(
            census, "Homo sapiens", value_filter=value_filter,
            column_names=["cell_type", "dataset_id", "development_stage",
                          "assay", "tissue_general"],
        )
        log(f"  Raw: {len(obs):,} cells [{time.time()-t0:.1f}s]")

        obs = obs[~obs["dataset_id"].isin(ts_excluded)]
        positive = obs["development_stage"].str.contains(
            r"year|adult|decade", case=False, na=False
        )
        negative = obs["development_stage"].str.contains(
            r"fetal|embryonic|newborn|infant|child", case=False, na=False
        )
        obs = obs[positive & ~negative].copy()
        obs["collection_name"] = obs["dataset_id"].map(ds_to_coll)
        log(f"  After filters: {len(obs):,} cells")

        # Collection metadata for DOIs
        coll_meta = {}
        for coll_name in obs["collection_name"].unique():
            coll_ds = datasets_df[datasets_df["collection_name"] == coll_name]
            doi = ""
            if "collection_doi" in coll_ds.columns:
                dois = coll_ds["collection_doi"].dropna().unique()
                if len(dois) > 0:
                    doi = str(dois[0])
            coll_meta[coll_name] = {"doi": doi}

    log("  Census closed.\n")

    # Per-collection per-type counts
    coll_type_counts = (
        obs.groupby(["collection_name", "cell_type"], observed=True)
        .size().reset_index(name="count")
    )
    coll_summary = obs.groupby("collection_name", observed=True).agg(
        total_cells=("cell_type", "size"),
        n_datasets=("dataset_id", "nunique"),
        n_tissues=("tissue_general", "nunique"),
    ).reset_index()

    coll_coverage = {}
    coll_assays = {}
    for coll_name in coll_summary["collection_name"]:
        sub = coll_type_counts[coll_type_counts["collection_name"] == coll_name]
        type_dict = dict(zip(sub["cell_type"], sub["count"]))
        coll_coverage[coll_name] = {
            "types_500": set(t for t in INTERSECTION_TYPES
                             if type_dict.get(t, 0) >= 500),
            "types_200": set(t for t in INTERSECTION_TYPES
                             if type_dict.get(t, 0) >= 200),
        }
        coll_obs = obs[obs["collection_name"] == coll_name]
        assay_vc = coll_obs["assay"].value_counts()
        coll_assays[coll_name] = assay_vc.index[0] if len(assay_vc) > 0 else "unknown"

    # Build ranked list
    ranked = []
    for _, row in coll_summary.iterrows():
        cn = row["collection_name"]
        ranked.append({
            "collection_name": cn,
            "total_cells": int(row["total_cells"]),
            "n_types_500": len(coll_coverage[cn]["types_500"]),
            "n_types_200": len(coll_coverage[cn]["types_200"]),
            "n_tissues": int(row["n_tissues"]),
            "top_assay": coll_assays[cn],
            "doi": coll_meta.get(cn, {}).get("doi", ""),
        })
    ranked.sort(key=lambda x: (-x["n_types_500"], -x["total_cells"]))

    # ── Greedy WITHOUT HCL ──────────────────────────────────────────────
    log(f"{'─'*78}")
    log("GREEDY SET COVER EXCLUDING HCL")
    log(f"{'─'*78}")

    target = set(INTERSECTION_TYPES)
    uncovered = target.copy()
    selected = []

    for step in range(1, 8):  # allow up to 7 to see where it converges
        if not uncovered:
            break
        best_coll = None
        best_marginal = 0
        best_newly = set()
        for r in ranked:
            cn = r["collection_name"]
            if cn == HCL_NAME:
                continue
            if cn in [s["collection_name"] for s in selected]:
                continue
            newly = coll_coverage[cn]["types_500"] & uncovered
            if len(newly) > best_marginal:
                best_marginal = len(newly)
                best_coll = r
                best_newly = newly
        if best_coll is None or best_marginal == 0:
            log(f"\n  Step {step}: No collection adds new types. STOP.")
            break
        selected.append(best_coll)
        uncovered -= best_newly
        log(f"\n  Step {step}: ADD '{best_coll['collection_name']}'")
        log(f"    Marginal: +{best_marginal} → {len(target)-len(uncovered)}/23")
        log(f"    New: {sorted(best_newly)}")
        log(f"    Assay: {best_coll['top_assay']}, Tissues: {best_coll['n_tissues']}")

    n_covered = len(target) - len(uncovered)

    log(f"\n  {'─'*60}")
    log(f"  WITHOUT HCL: {n_covered}/23 types, {len(selected)} collections")
    if uncovered:
        log(f"  UNCOVERED ({len(uncovered)}): {sorted(uncovered)}")

    # Key comparison
    log(f"\n{'─'*78}")
    log("COMPARISON")
    log(f"{'─'*78}")
    log(f"  With HCL:    20/23, 5 collections (but shared Guo lab with MCA)")
    log(f"  Without HCL: {n_covered}/23, {len(selected)} collections (lab-independent)")

    if n_covered >= 20 and len(selected) <= 4:
        log(f"\n  ✓ Lab-independent curated set is FEASIBLE")
        verdict = "CURATED_NO_HCL"
    elif n_covered >= 18 and len(selected) <= 5:
        log(f"\n  ~ Lab-independent curated set is MARGINAL ({n_covered}/23)")
        verdict = "CURATED_NO_HCL_MARGINAL"
    else:
        log(f"\n  ✗ Lab-independent curated set FAILS — need pooled approach")
        verdict = "POOLED"

    log(f"\n  VERDICT: {verdict}")

    # Describe institutions in selected set
    log(f"\n  Selected collections (no HCL):")
    for i, s in enumerate(selected, 1):
        log(f"  {i}. {s['collection_name']}")
        log(f"     Assay: {s['top_assay']} | Tissues: {s['n_tissues']} | "
            f"Cells: {s['total_cells']:,}")
        if s["doi"]:
            log(f"     DOI: {s['doi']}")

    # Check original 6
    orig_6 = {"B cell", "CD4-positive, alpha-beta T cell",
              "CD8-positive, alpha-beta T cell", "endothelial cell",
              "hepatocyte", "macrophage"}
    covered = target - uncovered
    orig_covered = orig_6 & covered
    log(f"\n  Original 6 coverage: {len(orig_covered)}/6")
    if orig_6 - covered:
        log(f"  Missing: {sorted(orig_6 - covered)}")

    # Save
    result = {
        "without_hcl_n_types": n_covered,
        "without_hcl_n_colls": len(selected),
        "without_hcl_collections": [s["collection_name"] for s in selected],
        "without_hcl_uncovered": sorted(uncovered),
        "without_hcl_verdict": verdict,
        "with_hcl_n_types": 20,
        "with_hcl_n_colls": 5,
        "shared_lab_risk": "HCL (Guo lab) shares PI and institution with MCA",
    }
    out = OUTPUT_DIR / "curated_nohcl_analysis.json"
    with open(out, "w") as f:
        json.dump(result, f, indent=2)
    log(f"\n  Saved: {out}")

    log(f"\n{'═'*78}")


if __name__ == "__main__":
    main()
