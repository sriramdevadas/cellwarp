#!/usr/bin/env python3
"""
CellWarp — PanSci Pre-Download Metadata Gate (Task 0)

Assesses cell type coverage, gene detection proxy, and tissue coverage
from PanSci (Cao lab, Science 2025) cell metadata WITHOUT downloading
the full count matrices. This is the go/no-go gate before committing
to a multi-GB download.

Uses the df_cell.csv.gz metadata files from GEO GSE247719 RAW archive.
Filters to WT genotype, 06_months age (young adult reference) to match
the approach used for Sun2023 (young sedentary controls).

Output: data/replication/pansci/metadata_gate_report.json
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pandas as pd

METADATA_DIR = Path("data/replication/pansci")
OUTPUT_PATH = METADATA_DIR / "metadata_gate_report.json"

TISSUES = [
    "kidney", "lung", "liver", "heart", "muscle", "stomach",
    "BAT", "iWAT", "gWAT", "ileum", "colon", "jejunum", "duodenum",
]

# Our 35-type ontology — specificity-ordered keyword mapping
# Based on DECISION-088 pattern: most specific first, with exclusions
TYPE_MAPPING = {
    # Immune — specific subtypes first
    "CD4-positive, alpha-beta T cell": {
        "keywords": ["cd4", "cd4+"],
        "exclude": [],
    },
    "CD8-positive, alpha-beta T cell": {
        "keywords": ["cd8", "cd8+", "cytotoxic t"],
        "exclude": [],
    },
    "mature NK T cell": {
        "keywords": ["nkt", "nk t", "natural killer t"],
        "exclude": [],
    },
    "natural killer cell": {
        "keywords": ["nk cell", "natural killer", "nk-"],
        "exclude": ["nkt", "nk t"],
    },
    "classical monocyte": {
        "keywords": ["classical monocyte"],
        "exclude": [],
    },
    "non-classical monocyte": {
        "keywords": ["non-classical monocyte", "nonclassical monocyte"],
        "exclude": [],
    },
    "intermediate monocyte": {
        "keywords": ["intermediate monocyte"],
        "exclude": [],
    },
    "plasma cell": {
        "keywords": ["plasma cell", "plasma b"],
        "exclude": [],
    },
    "myeloid dendritic cell": {
        "keywords": ["dendritic cell", "dendritic"],
        "exclude": ["plasmacytoid"],
    },
    "B cell": {
        "keywords": ["b cell", "b-cell"],
        "exclude": ["plasma", "pro-b", "pre-b"],
    },
    "macrophage": {
        "keywords": ["macrophage", "kupffer"],
        "exclude": [],
    },
    "monocyte": {
        "keywords": ["monocyte"],
        "exclude": ["classical", "non-classical", "nonclassical", "intermediate"],
    },
    "myeloid leukocyte": {
        "keywords": ["myeloid cell", "myeloid leukocyte"],
        "exclude": ["dendritic", "macrophage", "monocyte", "mast"],
    },
    "granulocyte": {
        "keywords": ["granulocyte", "neutrophil", "basophil", "eosinophil"],
        "exclude": [],
    },
    "T cell": {
        "keywords": ["t cell", "t-cell", "lymphoid"],
        "exclude": ["nk", "mast", "b cell", "b-cell"],
    },
    # Stromal / structural
    "endothelial cell": {
        "keywords": ["endothelial", "vascular endothelial"],
        "exclude": ["lymphatic"],
    },
    "fibroblast": {
        "keywords": ["fibroblast", "fibro-adipogenic", "fibro–adipogenic"],
        "exclude": ["cardiac"],
    },
    "fibroblast of cardiac tissue": {
        "keywords": ["cardiac fibroblast"],
        "exclude": [],
    },
    "smooth muscle cell": {
        "keywords": ["smooth muscle", "mural cell", "pericyte", "vascular smooth"],
        "exclude": [],
    },
    "stromal cell": {
        "keywords": ["stromal cell", "stroma"],
        "exclude": [],
    },
    "adventitial cell": {
        "keywords": ["adventitial"],
        "exclude": [],
    },
    "mesenchymal stem cell": {
        "keywords": ["mesenchymal stem", "msc"],
        "exclude": [],
    },
    "mesenchymal stem cell of adipose tissue": {
        "keywords": ["adipose mesenchymal", "adipose msc"],
        "exclude": [],
    },
    # Epithelial — specific first
    "hepatocyte": {
        "keywords": ["hepatocyte"],
        "exclude": ["hepatoblast"],
    },
    "enterocyte of epithelium of large intestine": {
        "keywords": ["enterocyte"],
        "exclude": [],
    },
    "large intestine goblet cell": {
        "keywords": ["goblet"],
        "exclude": [],
    },
    "pancreatic acinar cell": {
        "keywords": ["acinar"],
        "exclude": [],
    },
    "pancreatic ductal cell": {
        "keywords": ["ductal"],
        "exclude": [],
    },
    "basal cell": {
        "keywords": ["basal cell", "basal epithelial"],
        "exclude": [],
    },
    "bladder urothelial cell": {
        "keywords": ["urothelial"],
        "exclude": [],
    },
    "luminal epithelial cell of mammary gland": {
        "keywords": ["luminal epithelial", "mammary epithelial"],
        "exclude": [],
    },
    "epithelial cell": {
        "keywords": ["epithelial"],
        "exclude": ["basal", "luminal", "urothelial", "mammary"],
    },
    # Progenitors
    "hematopoietic stem cell": {
        "keywords": ["hematopoietic stem", "hsc"],
        "exclude": [],
    },
    "hematopoietic precursor cell": {
        "keywords": ["hematopoietic precursor", "hpc", "progenitor"],
        "exclude": ["fibro-adipogenic", "fibro–adipogenic"],
    },
}

# Reverse: which PanSci types map to which of our types
# Built automatically from the keyword matching


def map_pansci_type(pansci_name: str) -> str | None:
    """Map a PanSci cell type annotation to our 35-type ontology.

    Uses specificity-ordered keyword matching. Returns the first match
    or None if no match found.
    """
    name_lower = pansci_name.lower()

    for our_type, rule in TYPE_MAPPING.items():
        # Check exclusions first
        excluded = False
        for exc in rule["exclude"]:
            if exc in name_lower:
                excluded = True
                break
        if excluded:
            continue

        # Check keywords
        for kw in rule["keywords"]:
            if kw in name_lower:
                return our_type

    return None


def main():
    print("=" * 70)
    print("PanSci Pre-Download Metadata Gate (Task 0)")
    print("=" * 70)

    # Load all tissue metadata
    all_dfs = []
    for tissue in TISSUES:
        path = METADATA_DIR / f"{tissue}_df_cell.csv.gz"
        if not path.exists():
            print(f"  WARNING: {tissue} metadata not found")
            continue
        with gzip.open(path, "rt") as f:
            df = pd.read_csv(f)
        df["tissue_file"] = tissue
        all_dfs.append(df)
        print(f"  {tissue:<12} {len(df):>10,} cells  "
              f"(cols: {len(df.columns)}, types: {df['main_cell_type_organ'].nunique()})")

    combined = pd.concat(all_dfs, ignore_index=True)
    print(f"\n  Total: {len(combined):,} cells across {len(TISSUES)} organs")
    print(f"  Columns: {list(combined.columns)}")

    # Filter to WT genotype, 06_months
    print(f"\n  Genotype distribution:")
    print(combined["genotype"].value_counts().to_string())
    print(f"\n  Age distribution:")
    print(combined["age_group"].value_counts().to_string())

    wt_6m = combined[
        (combined["genotype"] == "WT") & (combined["age_group"] == "06_months")
    ].copy()
    print(f"\n  After filter (WT, 06_months): {len(wt_6m):,} cells")

    # Get all unique cell types across organs
    all_pansci_types = sorted(wt_6m["main_cell_type_organ"].unique())
    print(f"  Unique PanSci cell types (WT, 6m): {len(all_pansci_types)}")

    # Strip organ suffix for mapping
    # PanSci names are like "Vascular endothelial cells-Muscle"
    # We need to strip the "-Organ" suffix
    def strip_organ_suffix(name: str) -> str:
        """Strip organ suffix from PanSci type name."""
        for tissue in TISSUES:
            tissue_cap = tissue.capitalize()
            if name.endswith(f"-{tissue_cap}"):
                return name[: -len(f"-{tissue_cap}")]
            if name.endswith(f"-{tissue}"):
                return name[: -len(f"-{tissue}")]
        # Handle multi-word organ names
        for suffix in [
            "-Kidney", "-Lung", "-Liver", "-Heart", "-Muscle",
            "-Stomach", "-BAT", "-iWAT", "-gWAT", "-Ileum",
            "-Colon", "-Jejunum", "-Duodenum",
        ]:
            if name.endswith(suffix):
                return name[: -len(suffix)]
        return name

    wt_6m = wt_6m.copy()
    wt_6m["pansci_type_base"] = wt_6m["main_cell_type_organ"].apply(strip_organ_suffix)

    base_types = sorted(wt_6m["pansci_type_base"].unique())
    print(f"\n  Unique base types (organ suffix stripped): {len(base_types)}")
    for bt in base_types:
        n = (wt_6m["pansci_type_base"] == bt).sum()
        mapped = map_pansci_type(bt)
        tag = f" → {mapped}" if mapped else " → UNMAPPED"
        print(f"    {bt:<55} {n:>8,}{tag}")

    # Map to our ontology
    wt_6m["our_type"] = wt_6m["pansci_type_base"].apply(map_pansci_type)

    # Report coverage
    print("\n" + "=" * 70)
    print("Task 0b: Cell Type Coverage Assessment")
    print("=" * 70)

    mapped_counts = wt_6m[wt_6m["our_type"].notna()].groupby("our_type").size()
    unmapped_cells = (wt_6m["our_type"].isna()).sum()

    print(f"\n  Mapped cells: {len(wt_6m) - unmapped_cells:,}")
    print(f"  Unmapped cells: {unmapped_cells:,}")

    # Check all 35 types
    audit_rows = []
    for our_type in sorted(TYPE_MAPPING.keys()):
        n = int(mapped_counts.get(our_type, 0))
        if n >= 500:
            status = "PASS"
        elif n >= 200:
            status = "BORDERLINE"
        else:
            status = "ABSENT"
        audit_rows.append({
            "cell_type": our_type, "n_cells_wt_6m": n, "status": status,
        })
        if n > 0:
            print(f"  {our_type:<55} {n:>8,}  {status}")
        else:
            print(f"  {our_type:<55} {n:>8}  {status}")

    n_pass = sum(1 for r in audit_rows if r["status"] == "PASS")
    n_border = sum(1 for r in audit_rows if r["status"] == "BORDERLINE")
    n_absent = sum(1 for r in audit_rows if r["status"] == "ABSENT")
    print(f"\n  PASS (≥500): {n_pass}")
    print(f"  BORDERLINE (200-499): {n_border}")
    print(f"  ABSENT (<200): {n_absent}")
    print(f"  Total usable: {n_pass + n_border}")

    # Check original 6 types
    original_6 = [
        "B cell", "CD4-positive, alpha-beta T cell",
        "CD8-positive, alpha-beta T cell", "endothelial cell",
        "hepatocyte", "macrophage",
    ]
    print(f"\n  Original 6 Phase 2 types:")
    all_6_present = True
    for ct in original_6:
        n = int(mapped_counts.get(ct, 0))
        present = "PRESENT" if n >= 200 else "ABSENT"
        if n < 200:
            all_6_present = False
        print(f"    {ct:<50} {n:>8,}  {present}")

    # Task 0c: Gene detection
    print("\n" + "=" * 70)
    print("Task 0c: Gene Detection Assessment")
    print("=" * 70)
    print(f"\n  Paper reports: median 1,040 UMIs/cell, mean 1,601 UMIs/cell")
    print(f"  Protocol: EasySci snRNA-seq (combinatorial indexing, NOT 10x)")
    print(f"  No per-cell gene count column in metadata — UMI count is proxy")
    print(f"  Estimated median genes/cell: ~500-800 (typical for ~1,000 UMI snRNA-seq)")
    print(f"  NOTE: Must verify from actual count matrix after download")
    gene_detection_status = "PROCEED_WITH_WARNING"  # 500-1000 range

    # Task 0d: Tissue coverage
    print("\n" + "=" * 70)
    print("Task 0d: Tissue Coverage")
    print("=" * 70)

    # For each mapped type, show which tissues it comes from
    mapped_wt_6m = wt_6m[wt_6m["our_type"].notna()].copy()

    type_tissue_breakdown = {}
    for our_type in sorted(mapped_counts.index):
        type_mask = mapped_wt_6m["our_type"] == our_type
        tissue_counts = mapped_wt_6m.loc[type_mask, "tissue_file"].value_counts()
        type_tissue_breakdown[our_type] = {
            t: int(n) for t, n in tissue_counts.items()
        }
        tissues_str = ", ".join(f"{t}({n:,})" for t, n in tissue_counts.head(5).items())
        print(f"  {our_type:<50} {tissues_str}")

    # ISSUE-092 lesson: flag endothelial tissue mismatch
    if "endothelial cell" in type_tissue_breakdown:
        endo_tissues = type_tissue_breakdown["endothelial cell"]
        print(f"\n  *** ISSUE-092 ENDOTHELIAL CHECK ***")
        print(f"  PanSci endothelial tissue breakdown (WT, 6m):")
        for t, n in sorted(endo_tissues.items(), key=lambda x: -x[1]):
            print(f"    {t:<15} {n:>8,}")
        print(f"  Tabula endothelial: predominantly myometrium/adipose/muscle/pancreas")
        print(f"  DECISION-101: Use tissue-matched endothelial (lung was best for Sun2023)")

    # GATE DECISIONS
    print("\n" + "=" * 70)
    print("GATE DECISION SUMMARY")
    print("=" * 70)

    fatal = False
    warnings = []

    # Coverage gate: ≥12 types
    if n_pass + n_border < 12:
        print(f"  FATAL: Coverage {n_pass + n_border} < 12 types")
        fatal = True
    else:
        print(f"  PASS: Coverage {n_pass + n_border} ≥ 12 types")

    # Original 6 gate
    if not all_6_present:
        missing = [ct for ct in original_6 if mapped_counts.get(ct, 0) < 200]
        print(f"  FATAL: Original 6 types not all present. Missing: {missing}")
        fatal = True
    else:
        print(f"  PASS: All 6 original Phase 2 types present")

    # Gene detection gate
    print(f"  WARNING: Median UMI ~1,040 (estimated median genes ~500-800)")
    print(f"    Must verify median genes ≥ 500 from count matrix")
    warnings.append("Gene detection unverified — median ~1,040 UMI, need count matrix check")

    if fatal:
        print(f"\n  *** FATAL GATE FIRED — STOP ***")
    else:
        print(f"\n  ALL GATES PASS (with warnings) — PROCEED TO TASK 1")

    # Save report
    report = {
        "diagnostic": "PanSci pre-download metadata gate (Task 0)",
        "date": "2026-03-15",
        "dataset": {
            "name": "PanSci (Cao lab, Science 2025)",
            "accession": "GEO GSE247719",
            "protocol": "EasySci snRNA-seq",
            "total_cells_all": len(combined),
            "total_cells_wt_6m": len(wt_6m),
            "n_organs": len(TISSUES),
            "organs": TISSUES,
            "filter": "WT genotype, 06_months age",
        },
        "task_0b_coverage": {
            "n_pass": n_pass,
            "n_borderline": n_border,
            "n_absent": n_absent,
            "total_usable": n_pass + n_border,
            "all_6_original_present": all_6_present,
            "audit": audit_rows,
            "type_tissue_breakdown": type_tissue_breakdown,
        },
        "task_0c_gene_detection": {
            "status": gene_detection_status,
            "paper_median_umi": 1040,
            "paper_mean_umi": 1601,
            "estimated_median_genes": "500-800 (unverified)",
            "note": "Must verify from count matrix. Below 500 is fatal.",
        },
        "task_0d_tissue_coverage": {
            "organs": TISSUES,
            "missing_vs_tabula": [
                "spleen", "bone_marrow", "blood", "pancreas",
                "mammary_gland", "bladder", "large_intestine_specific",
            ],
        },
        "gate_verdict": {
            "fatal": fatal,
            "warnings": warnings,
            "proceed": not fatal,
        },
    }

    with open(OUTPUT_PATH, "w") as f:
        json.dump(report, f, indent=2)
    print(f"\n  Report saved: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
