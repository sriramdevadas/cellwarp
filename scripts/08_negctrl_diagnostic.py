#!/usr/bin/env python3
"""
CellWarp — Negative Control Diagnostic

Pre-kill diagnostic for Phase 2. Three parts:

PART 1: Compare top-20 residual genes between human-vs-mouse and human-vs-human
Procrustes analyses. High overlap suggests the signal is the same (batch effects);
low overlap with pathway-specific genes in H-vs-M and ribosomal/housekeeping genes
in H-vs-H would suggest the signals are different (biology vs batch).

PART 2: Check developmental stage contamination in the second human atlas.
The negative control used cells from "Mapping the developing human immune system
across organs" — if a large fraction is fetal/neonatal, the human-vs-human
Procrustes distance is inflated by developmental differences, not batch effects.

PART 3: Recommendation based on findings.
"""

import json
from collections import Counter
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

PROJECT = Path(__file__).resolve().parent.parent
OUTPUT = PROJECT / "output" / "phase2" / "negative_control_diagnostic"
OUTPUT.mkdir(parents=True, exist_ok=True)

HVM_RESULTS = PROJECT / "output" / "phase2" / "procrustes_results.json"
HVH_RESULTS = PROJECT / "output" / "phase2" / "negative_control" / "negctrl_results.json"
HUMAN2_DATA = PROJECT / "data" / "phase2" / "human2_negctrl.h5ad"
SOURCE_DATASETS = PROJECT / "output" / "phase2" / "negative_control" / "source_datasets.json"


# ---------------------------------------------------------------------------
# Gene classification helpers
# ---------------------------------------------------------------------------

# Ribosomal protein genes (RPL*, RPS*)
def is_ribosomal(gene: str) -> bool:
    """Ribosomal protein genes start with RPL or RPS."""
    return gene.startswith("RPL") or gene.startswith("RPS")


# Mitochondrial genes
def is_mitochondrial(gene: str) -> bool:
    """Mitochondrial genes start with MT-."""
    return gene.startswith("MT-")


# Known cell cycle genes (conservative list of markers)
CELL_CYCLE_GENES = {
    "MKI67", "TOP2A", "PCNA", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6",
    "MCM7", "CDK1", "CDK2", "CDK4", "CDK6", "CCNA2", "CCNB1", "CCNB2",
    "CCND1", "CCND2", "CCNE1", "CCNE2", "CDC20", "CDC25A", "CDC25B",
    "CDC25C", "BUB1", "BUB1B", "MAD2L1", "AURKA", "AURKB", "PLK1",
    "BIRC5", "UBE2C", "HMGB2", "STMN1", "TUBB", "TUBA1B",
}


# Known housekeeping / stress response genes
HOUSEKEEPING_STRESS = {
    "ACTB", "GAPDH", "B2M", "TMSB10", "TMSB4X", "FTH1", "FTL",
    "HSP90AA1", "HSP90AB1", "HSPA8", "HSPA1A", "HSPA1B",
    "VIM", "EEF1A1", "EEF1G", "EEF2",
    "SOD2", "GPX1", "NFKBIA",
}


def classify_gene(gene: str) -> str:
    """Classify a gene into functional category."""
    if is_ribosomal(gene):
        return "ribosomal"
    if is_mitochondrial(gene):
        return "mitochondrial"
    if gene in CELL_CYCLE_GENES:
        return "cell_cycle"
    if gene in HOUSEKEEPING_STRESS:
        return "housekeeping/stress"
    return "pathway-specific"


# ---------------------------------------------------------------------------
# PART 1: Compare residual genes
# ---------------------------------------------------------------------------

def part1_compare_residuals() -> dict:
    """Compare top-20 residual genes between H-vs-M and H-vs-H analyses."""
    print("=" * 80)
    print("PART 1: RESIDUAL GENE COMPARISON")
    print("=" * 80)

    with open(HVM_RESULTS) as f:
        hvm = json.load(f)
    with open(HVH_RESULTS) as f:
        hvh = json.load(f)

    cell_types = hvm["cell_types"]
    results = {}

    for ct in cell_types:
        hvm_genes = [g["gene"] for g in hvm["top_genes_per_cell_type"][ct]]
        hvh_genes = [g["gene"] for g in hvh["top_genes_per_cell_type"][ct]]

        hvm_set = set(hvm_genes)
        hvh_set = set(hvh_genes)
        overlap = hvm_set & hvh_set
        n_overlap = len(overlap)

        # Classify genes
        hvm_classes = [classify_gene(g) for g in hvm_genes]
        hvh_classes = [classify_gene(g) for g in hvh_genes]

        hvm_class_counts = Counter(hvm_classes)
        hvh_class_counts = Counter(hvh_classes)

        # Residual magnitudes
        hvm_mag = hvm["residuals"][ct]["magnitude"]
        hvh_mag = hvh["residuals"][ct]["magnitude"]

        results[ct] = {
            "hvm_genes": hvm_genes,
            "hvh_genes": hvh_genes,
            "overlap": sorted(overlap),
            "n_overlap": n_overlap,
            "overlap_fraction": n_overlap / 20,
            "hvm_gene_classes": dict(hvm_class_counts),
            "hvh_gene_classes": dict(hvh_class_counts),
            "hvm_residual_magnitude": hvm_mag,
            "hvh_residual_magnitude": hvh_mag,
        }

        # Print side-by-side comparison
        print(f"\n{'─' * 80}")
        print(f"  {ct}")
        print(f"  H-vs-M residual magnitude: {hvm_mag:.3f}")
        print(f"  H-vs-H residual magnitude: {hvh_mag:.3f}")
        print(f"  Overlap: {n_overlap}/20 genes ({n_overlap/20*100:.0f}%)")
        if overlap:
            print(f"  Shared genes: {', '.join(sorted(overlap))}")
        print()
        print(f"  {'Rank':<6} {'Human-vs-Mouse':<25} {'Class':<20} {'Human-vs-Human':<25} {'Class':<20}")
        print(f"  {'─'*6} {'─'*25} {'─'*20} {'─'*25} {'─'*20}")
        for i in range(20):
            hvm_g = hvm_genes[i]
            hvh_g = hvh_genes[i]
            hvm_c = classify_gene(hvm_g)
            hvh_c = classify_gene(hvh_g)
            # Mark overlapping genes
            hvm_mark = " *" if hvm_g in overlap else ""
            hvh_mark = " *" if hvh_g in overlap else ""
            print(f"  {i+1:<6} {hvm_g + hvm_mark:<25} {hvm_c:<20} {hvh_g + hvh_mark:<25} {hvh_c:<20}")

        print(f"\n  Gene class distribution:")
        all_classes = sorted(set(list(hvm_class_counts.keys()) + list(hvh_class_counts.keys())))
        for cls in all_classes:
            hvm_n = hvm_class_counts.get(cls, 0)
            hvh_n = hvh_class_counts.get(cls, 0)
            print(f"    {cls:<25} H-vs-M: {hvm_n:>2}/20    H-vs-H: {hvh_n:>2}/20")

    # Summary
    print(f"\n{'=' * 80}")
    print("PART 1 SUMMARY")
    print(f"{'=' * 80}")
    print(f"\n  {'Cell Type':<45} {'Overlap':>8} {'H-vs-M mag':>12} {'H-vs-H mag':>12} {'Ratio':>8}")
    print(f"  {'─'*85}")
    for ct in cell_types:
        r = results[ct]
        ratio = r["hvh_residual_magnitude"] / r["hvm_residual_magnitude"] if r["hvm_residual_magnitude"] > 0 else float('inf')
        print(f"  {ct:<45} {r['n_overlap']:>5}/20 {r['hvm_residual_magnitude']:>12.3f} {r['hvh_residual_magnitude']:>12.3f} {ratio:>8.2f}")

    return results


# ---------------------------------------------------------------------------
# PART 2: Developmental stage contamination
# ---------------------------------------------------------------------------

def part2_developmental_check() -> dict:
    """Check fraction of fetal/neonatal cells in second human atlas."""
    print(f"\n\n{'=' * 80}")
    print("PART 2: DEVELOPMENTAL STAGE CONTAMINATION CHECK")
    print(f"{'=' * 80}")

    # Load second human atlas
    print(f"\n  Loading {HUMAN2_DATA}...")
    adata = ad.read_h5ad(HUMAN2_DATA)
    print(f"  Shape: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

    # Print available metadata columns
    print(f"\n  Available .obs columns: {list(adata.obs.columns)}")

    # Check for development_stage or similar column
    results = {
        "total_cells": adata.shape[0],
        "obs_columns": list(adata.obs.columns),
        "cell_type_breakdown": {},
    }

    # Look for developmental stage info
    dev_cols = [c for c in adata.obs.columns if any(
        kw in c.lower() for kw in ["develop", "stage", "age", "fetal", "donor"]
    )]
    print(f"  Development-related columns: {dev_cols}")

    # Also check dataset_id to trace back to source
    dataset_cols = [c for c in adata.obs.columns if any(
        kw in c.lower() for kw in ["dataset", "collection", "source", "study", "batch"]
    )]
    print(f"  Dataset-related columns: {dataset_cols}")

    # Check all potentially relevant columns
    for col in dev_cols + dataset_cols:
        if col in adata.obs.columns:
            print(f"\n  Values in '{col}':")
            vc = adata.obs[col].value_counts()
            for val, count in vc.items():
                print(f"    {val}: {count:,} ({count/adata.shape[0]*100:.1f}%)")

    # Per cell type breakdown
    print(f"\n  Per cell type breakdown:")
    cell_types = sorted(adata.obs["cell_type"].unique())

    # Determine which column to use for developmental stage
    # Common CELLxGENE Census columns
    stage_col = None
    for candidate in ["development_stage", "development_stage_ontology_term_id",
                       "donor_age", "age"]:
        if candidate in adata.obs.columns:
            stage_col = candidate
            break

    if stage_col is None:
        # Try dataset_id as proxy
        print("\n  WARNING: No explicit development_stage column found.")
        print("  Using dataset_id as proxy for source identification.")
        stage_col = "dataset_id" if "dataset_id" in adata.obs.columns else None

    # Load source datasets info
    with open(SOURCE_DATASETS) as f:
        source_info = json.load(f)

    # The developing immune system dataset
    DEV_DATASET_ID = "fd072bc3-2dfb-46f8-b4e3-467cb3223182"
    COPD_DATASET_ID = "8fbed309-d3d4-441b-b3ff-e2dcbcec2d35"

    flagged_types = []

    for ct in cell_types:
        mask = adata.obs["cell_type"] == ct
        ct_data = adata.obs[mask]
        n_total = len(ct_data)

        ct_result = {"n_total": n_total}

        # Check if we can identify fetal cells
        if "dataset_id" in ct_data.columns:
            # Count cells from each source dataset
            ds_counts = ct_data["dataset_id"].value_counts()
            n_dev = 0
            n_copd = 0
            n_other = 0

            for ds_id, count in ds_counts.items():
                if ds_id == DEV_DATASET_ID:
                    n_dev = count
                elif ds_id == COPD_DATASET_ID:
                    n_copd = count
                else:
                    n_other = count

            fetal_frac = n_dev / n_total if n_total > 0 else 0
            ct_result["n_from_developmental_atlas"] = n_dev
            ct_result["n_from_copd_atlas"] = n_copd
            ct_result["n_from_other"] = n_other
            ct_result["fetal_fraction"] = fetal_frac
            ct_result["flagged"] = fetal_frac > 0.20

            if fetal_frac > 0.20:
                flagged_types.append(ct)

            print(f"\n  {ct} ({n_total:,} cells):")
            print(f"    From developmental atlas: {n_dev:,} ({fetal_frac*100:.1f}%)")
            print(f"    From COPD atlas: {n_copd:,} ({(n_copd/n_total*100 if n_total else 0):.1f}%)")
            if n_other > 0:
                print(f"    From other: {n_other:,}")
            if fetal_frac > 0.20:
                print(f"    *** FLAGGED: fetal fraction {fetal_frac*100:.1f}% > 20% ***")

        # Check development_stage if available
        if stage_col and stage_col != "dataset_id" and stage_col in ct_data.columns:
            print(f"    Development stages:")
            stages = ct_data[stage_col].value_counts()
            for stage, count in stages.items():
                print(f"      {stage}: {count:,} ({count/n_total*100:.1f}%)")

        results["cell_type_breakdown"][ct] = ct_result

    results["flagged_cell_types"] = flagged_types

    print(f"\n{'=' * 80}")
    print("PART 2 SUMMARY")
    print(f"{'=' * 80}")
    if flagged_types:
        print(f"\n  FLAGGED cell types (fetal fraction > 20%):")
        for ct in flagged_types:
            frac = results["cell_type_breakdown"][ct]["fetal_fraction"]
            print(f"    - {ct}: {frac*100:.1f}% from developmental atlas")
    else:
        print("\n  No cell types flagged (all < 20% fetal).")

    # Additional: check the development_stage column values if present
    if "development_stage" in adata.obs.columns:
        print(f"\n  Overall development stage distribution:")
        stages = adata.obs["development_stage"].value_counts()
        for stage, count in stages.items():
            print(f"    {stage}: {count:,} ({count/adata.shape[0]*100:.1f}%)")

        # Classify as adult vs non-adult
        non_adult_keywords = ["fetal", "embryo", "neonat", "infant", "child",
                              "newborn", "prenatal", "gestational", "trimester",
                              "week post-fertilization", "Carnegie"]
        adult_count = 0
        non_adult_count = 0
        for stage, count in stages.items():
            stage_str = str(stage).lower()
            if any(kw in stage_str for kw in non_adult_keywords):
                non_adult_count += count
            else:
                adult_count += count

        results["overall_adult_count"] = adult_count
        results["overall_non_adult_count"] = non_adult_count
        results["overall_non_adult_fraction"] = non_adult_count / adata.shape[0]

        print(f"\n  Classification:")
        print(f"    Adult cells: {adult_count:,} ({adult_count/adata.shape[0]*100:.1f}%)")
        print(f"    Non-adult cells: {non_adult_count:,} ({non_adult_count/adata.shape[0]*100:.1f}%)")

        # Per cell type, adult vs non-adult
        print(f"\n  Per cell type adult/non-adult breakdown:")
        print(f"  {'Cell Type':<45} {'Total':>7} {'Adult':>7} {'Non-adult':>10} {'Non-adult %':>12}")
        print(f"  {'─'*81}")
        for ct in cell_types:
            mask = adata.obs["cell_type"] == ct
            ct_stages = adata.obs.loc[mask, "development_stage"]
            n_total = len(ct_stages)
            n_nonadult = sum(
                1 for s in ct_stages
                if any(kw in str(s).lower() for kw in non_adult_keywords)
            )
            n_adult = n_total - n_nonadult
            frac = n_nonadult / n_total if n_total > 0 else 0
            flag = " *** FLAGGED" if frac > 0.20 else ""
            print(f"  {ct:<45} {n_total:>7,} {n_adult:>7,} {n_nonadult:>10,} {frac*100:>11.1f}%{flag}")

            results["cell_type_breakdown"][ct]["n_adult"] = n_adult
            results["cell_type_breakdown"][ct]["n_non_adult"] = n_nonadult
            results["cell_type_breakdown"][ct]["non_adult_fraction"] = frac
            # Update flag based on actual development_stage
            if frac > 0.20:
                results["cell_type_breakdown"][ct]["flagged"] = True
                if ct not in flagged_types:
                    flagged_types.append(ct)

        results["flagged_cell_types"] = flagged_types

    return results


# ---------------------------------------------------------------------------
# PART 3: Recommendation
# ---------------------------------------------------------------------------

def part3_recommendation(part1_results: dict, part2_results: dict) -> str:
    """Generate recommendation based on diagnostic findings."""
    print(f"\n\n{'=' * 80}")
    print("PART 3: RECOMMENDATION")
    print(f"{'=' * 80}")

    # Summarize key findings
    total_overlap = sum(r["n_overlap"] for r in part1_results.values())
    max_possible = len(part1_results) * 20
    avg_overlap = total_overlap / len(part1_results)

    flagged = part2_results.get("flagged_cell_types", [])
    non_adult_frac = part2_results.get("overall_non_adult_fraction", None)

    print(f"\n  Key findings:")
    print(f"  1. Average gene overlap: {avg_overlap:.1f}/20 ({avg_overlap/20*100:.0f}%)")
    print(f"     Total overlap across {len(part1_results)} cell types: {total_overlap}/{max_possible}")
    print(f"  2. Flagged cell types (>20% fetal): {len(flagged)}")
    if flagged:
        print(f"     {', '.join(flagged)}")
    if non_adult_frac is not None:
        print(f"  3. Overall non-adult fraction: {non_adult_frac*100:.1f}%")

    # Analyze residual patterns
    print(f"\n  Residual magnitude comparison:")
    for ct, r in part1_results.items():
        hvm = r["hvm_residual_magnitude"]
        hvh = r["hvh_residual_magnitude"]
        ratio = hvh / hvm if hvm > 0 else float('inf')
        direction = "H-vs-H LARGER" if ratio > 1.5 else "H-vs-H SMALLER" if ratio < 0.5 else "SIMILAR"
        print(f"    {ct:<40} H-vs-M={hvm:.2f}  H-vs-H={hvh:.2f}  ratio={ratio:.2f}  ({direction})")

    # Gene class analysis
    print(f"\n  Gene class analysis (across all cell types):")
    hvm_total_classes = Counter()
    hvh_total_classes = Counter()
    for r in part1_results.values():
        for cls, n in r["hvm_gene_classes"].items():
            hvm_total_classes[cls] += n
        for cls, n in r["hvh_gene_classes"].items():
            hvh_total_classes[cls] += n

    all_classes = sorted(set(list(hvm_total_classes.keys()) + list(hvh_total_classes.keys())))
    for cls in all_classes:
        hvm_n = hvm_total_classes.get(cls, 0)
        hvh_n = hvh_total_classes.get(cls, 0)
        print(f"    {cls:<25} H-vs-M: {hvm_n:>3}/{max_possible}    H-vs-H: {hvh_n:>3}/{max_possible}")

    # Decision logic
    recommendation = []

    # Check 1: Is overlap high enough to indicate same signal?
    if avg_overlap > 8:
        recommendation.append("HIGH gene overlap — signals are likely related (same batch axis)")
    elif avg_overlap > 4:
        recommendation.append("MODERATE gene overlap — partial shared signal")
    else:
        recommendation.append("LOW gene overlap — signals may be distinct")

    # Check 2: Are H-vs-H genes more technical (ribosomal, mito, housekeeping)?
    hvh_technical = hvh_total_classes.get("ribosomal", 0) + hvh_total_classes.get("mitochondrial", 0) + hvh_total_classes.get("housekeeping/stress", 0)
    hvm_technical = hvm_total_classes.get("ribosomal", 0) + hvm_total_classes.get("mitochondrial", 0) + hvm_total_classes.get("housekeeping/stress", 0)
    if hvh_technical > hvm_technical * 1.3:
        recommendation.append("H-vs-H residuals are MORE technical — suggests batch effects in negative control")
    elif hvh_technical < hvm_technical * 0.7:
        recommendation.append("H-vs-H residuals are LESS technical — concerning, not just batch")
    else:
        recommendation.append("Similar technical gene fraction in both — ambiguous")

    # Check 3: Developmental contamination
    if len(flagged) >= 2:
        recommendation.append(f"SIGNIFICANT developmental contamination ({len(flagged)} cell types flagged)")
        recommendation.append("→ RECOMMEND: Rerun negative control with adult-only cells (Option A)")
    elif len(flagged) == 1:
        recommendation.append(f"PARTIAL developmental contamination ({flagged[0]} flagged)")
        recommendation.append("→ RECOMMEND: Rerun with adult-only, or exclude flagged type as sensitivity check")

    print(f"\n  {'─' * 70}")
    print(f"  ASSESSMENT:")
    for i, r in enumerate(recommendation, 1):
        print(f"    {i}. {r}")

    # Final recommendation
    print(f"\n  FINAL RECOMMENDATION:")
    if len(flagged) >= 1:
        print(f"""
    Option A — RERUN NEGATIVE CONTROL WITH ADULT-ONLY TISSUE (RECOMMENDED)

    The negative control failure is likely confounded by developmental stage
    contamination. {len(flagged)} cell type(s) have >20% fetal/neonatal cells from
    the "Mapping the developing human immune system" dataset. Fetal cells have
    fundamentally different expression profiles from adult cells, which would
    inflate Procrustes distance and make the negative control unfairly stringent.

    Before killing Phase 2, the negative control deserves a fair test with
    properly matched adult-only tissue. Specifically:
    1. Filter the second atlas to adult cells only, OR
    2. Source a different second human atlas (e.g., HLCA, HuBMAP) that is
       all-adult and covers the 6 cell types.

    If the rerun still fails, THEN proceed with Option B or C.
""")
        final = "A"
    else:
        print(f"""
    The developmental stage does not appear to be a major confound. The
    negative control failure likely reflects genuine batch/atlas effects.
    Consider Option C (reframe around residual interpretation) as the
    most promising path forward.
""")
        final = "C"

    return final


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("\n" + "=" * 80)
    print("  CellWarp — Negative Control Diagnostic")
    print("  Pre-kill assessment for Phase 2")
    print("=" * 80)

    # Part 1
    part1_results = part1_compare_residuals()

    # Part 2
    part2_results = part2_developmental_check()

    # Part 3
    final_rec = part3_recommendation(part1_results, part2_results)

    # Save all results
    output = {
        "part1_residual_comparison": part1_results,
        "part2_developmental_check": part2_results,
        "final_recommendation": final_rec,
    }

    out_path = OUTPUT / "diagnostic_results.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\n  Results saved to: {out_path}")

    # Also save a human-readable summary
    summary_path = OUTPUT / "diagnostic_summary.txt"
    print(f"  Summary saved to: {summary_path}")


if __name__ == "__main__":
    main()
