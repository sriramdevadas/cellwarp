#!/usr/bin/env python3
"""
T1-A MCA Feasibility Check: Mouse Cell Atlas metadata-only query.

BIOLOGY: We need an independent mouse single-cell atlas (different from Tabula Muris
Senis) to replicate our cross-species Procrustes results. MCA (Han et al. 2018, Cell)
used microwell-seq at BGI/Zhejiang University — completely independent technology and
institution from the Tabula consortium (10x Genomics, CZI-funded). If MCA covers enough
of our 35 cell types at ≥500 cells each, we can use it for T1-A replication.

MATH: No heavy computation — pure metadata inventory. We map MCA cell type annotations
to our 35-type ontology via fuzzy string matching, then count cells per type.

Data source: MCA_CellAssignments.csv from figshare (dataset 5435866)
GEO accession: GSE108097
Paper: Han et al. "Mapping the Mouse Cell Atlas by Microwell-Seq" Cell 172, 1091-1107 (2018)
"""

import os
import sys
import json
import re
from collections import defaultdict

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Patch

# ── Configuration ──────────────────────────────────────────────────────────
OUTPUT_DIR = "output/validation/mca_feasibility"
os.makedirs(OUTPUT_DIR, exist_ok=True)

MCA_CSV = "/tmp/MCA_CellAssignments.csv"

# Our 35 cell types (from cell_type_inventory.csv, passes_200_gate=True)
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

# Original 6 types from Phase 2 (critical for replication)
ORIGINAL_6 = [
    "B cell",
    "CD4-positive, alpha-beta T cell",
    "CD8-positive, alpha-beta T cell",
    "endothelial cell",
    "hepatocyte",
    "macrophage",
]


def strip_tissue_suffix(annotation):
    """
    MCA annotations have format: 'CellType(Tissue)' or 'CellType_marker(Tissue)'
    Strip the tissue suffix in parentheses to get the base cell type.
    """
    # Remove (Tissue) suffix
    base = re.sub(r'\([^)]+\)\s*$', '', annotation).strip()
    return base


def map_types(mca_base_labels, our_types):
    """
    Two-pass mapping of MCA base cell type labels to our 35-type ontology.

    Pass 1: Keyword-based fuzzy matching
    Pass 2: Manual harmonization with confidence flags

    MCA uses descriptive labels like "B cell", "Hepatocyte", "Endothelial cell_Aqp1 high"
    Our types use CL ontology names like "CD4-positive, alpha-beta T cell"

    Returns: mapping dict {mca_base_label: our_label}, review_needed list, unmapped set
    """
    # Keyword mapping rules: our_type → list of keywords
    # IMPORTANT: Use word-boundary-aware matching for short keywords
    # to avoid substring false positives (e.g., "t cell" in "mast cell")
    keyword_map = {
        "B cell": ["b cell", "pre-pro b cell"],
        "CD4-positive, alpha-beta T cell": [],  # MCA has NO CD4 annotation
        "CD8-positive, alpha-beta T cell": ["t cell_cd8", "t cell_cd8b1"],
        "T cell": [],  # handled specially below with exclusion list
        "adventitial cell": ["adventitial"],
        "basal cell": ["basal cell", "basal epithelial"],
        "bladder urothelial cell": ["urothelium", "urothelial"],
        "classical monocyte": ["classical monocyte", "ly6c+ monocyte"],
        "endothelial cell": ["endothelial cell", "vascular endothelial",
                              "endothelial cells"],
        "enterocyte of epithelium of large intestine": [
            "enterocyte", "epithelium of small intestinal villi"],
        "epithelial cell": ["epithelial cell"],
        "fibroblast": ["fibroblast"],
        "fibroblast of cardiac tissue": ["cardiac fibroblast"],
        "granulocyte": ["granulocyte"],
        "hematopoietic precursor cell": [
            "hematopoietic stem progenitor",
            "monocyte progenitor cell", "monocyte progenitor",
            "megakaryocyte progenitor", "eosinophil progenitor",
            "lymphoid progenitor", "granulocyte monocyte progenitor",
            "multipotent progenitor",
        ],
        "hematopoietic stem cell": ["hematopoietic stem cell"],
        "hepatocyte": ["hepatocyte", "pericentral (pc) hepatocyte",
                        "periportal (pp) hepatocyte"],
        "intermediate monocyte": ["intermediate monocyte"],
        "large intestine goblet cell": ["goblet cell"],
        "luminal epithelial cell of mammary gland": [
            "ductal luminal cell", "luminal cell", "luminal progenitor"],
        "macrophage": ["macrophage", "alveolar macrophage",
                        "interstitial macrophage"],
        "mature NK T cell": ["nkt cell", "nk t cell"],
        "mesenchymal stem cell": ["mesenchymal stem cell", "mesenchymal stromal cell"],
        "mesenchymal stem cell of adipose tissue": ["adipose mesenchymal",
                                                     "adipose stromal"],
        "monocyte": ["monocyte"],
        "myeloid dendritic cell": ["dendritic cell", "conventional dendritic cell",
                                    "plasmacytoid dendritic cell"],
        "myeloid leukocyte": ["myeloid leukocyte"],
        "natural killer cell": ["nk cell"],
        "neutrophil": ["neutrophil"],
        "non-classical monocyte": ["non-classical monocyte", "patrolling monocyte",
                                    "ly6c- monocyte"],
        "pancreatic acinar cell": ["acinar cell"],
        "pancreatic ductal cell": ["ductal cell"],
        "plasma cell": ["plasma cell", "ig-producing b cell", "ig−producing b cell"],
        "smooth muscle cell": ["smooth muscle cell"],
        "stromal cell": ["stromal cell", "decidual stromal cell"],
    }

    # Labels that should NEVER be mapped (false positive traps)
    exclude_labels = {
        "mast cell", "mast cell_mcpt8 high",            # NOT T cells
        "pit cell_ifrd1 high", "pit cell_gm26917 high", # stomach pit cells
        "tuft cell",                                      # NOT T cells
        "spiral artery trophoblast giant cells",          # NOT T cells
        "vascular smooth muscle progenitor cell",         # NOT hematopoietic
        "neural progenitor cell",                         # NOT hematopoietic
        "muscle progenitor cell",                         # NOT hematopoietic
        "tendon stem/progenitor cell",                    # NOT hematopoietic
        "endocrine progenitor cell",                      # NOT hematopoietic
        "stem and progenitor cell",                       # too ambiguous
        "progenitor cell",                                # too ambiguous
        "progenitor cell_ptprcap high",                   # ambiguous
        "dividing dendritic cells",                       # keep under dendritic
        "dividing t cells",                               # NOT CD4-specific
    }

    # Specificity ordering: more specific types should be matched first
    # to avoid generic types consuming specific labels
    specificity_order = [
        # Most specific first
        "CD4-positive, alpha-beta T cell",
        "CD8-positive, alpha-beta T cell",
        "mature NK T cell",
        "natural killer cell",
        "classical monocyte",
        "intermediate monocyte",
        "non-classical monocyte",
        "bladder urothelial cell",
        "enterocyte of epithelium of large intestine",
        "large intestine goblet cell",
        "luminal epithelial cell of mammary gland",
        "fibroblast of cardiac tissue",
        "mesenchymal stem cell of adipose tissue",
        "hematopoietic stem cell",
        "hematopoietic precursor cell",
        "pancreatic acinar cell",
        "pancreatic ductal cell",
        "myeloid dendritic cell",
        "adventitial cell",
        "plasma cell",
        "smooth muscle cell",
        # Then broad types
        "B cell",
        "T cell",
        "endothelial cell",
        "hepatocyte",
        "macrophage",
        "fibroblast",
        "basal cell",
        "neutrophil",
        "granulocyte",
        "monocyte",
        "mesenchymal stem cell",
        "epithelial cell",
        "stromal cell",
        "myeloid leukocyte",
    ]

    mapping = {}  # mca_base_label → our_label
    review_needed = []
    already_mapped = set()  # MCA labels already claimed by a specific type

    print(f"\n{'='*70}")
    print("PASS 1 & 2: KEYWORD MATCHING + SPECIFICITY RESOLUTION")
    print(f"{'='*70}")

    for our_type in specificity_order:
        keywords = keyword_map.get(our_type, [])
        matched_labels = []

        for mca_label in mca_base_labels:
            if mca_label in already_mapped:
                continue
            mca_lower = mca_label.lower()

            # Skip globally excluded labels
            if mca_lower in exclude_labels:
                continue

            # Special handling for T cell: match labels that START with
            # "t cell" or contain " t cell" (word boundary), excluding
            # "mast cell", "pit cell", etc.
            if our_type == "T cell":
                if (mca_lower.startswith("t cell") or
                    mca_lower.startswith("abt cell") or
                    mca_lower.startswith("gdt cell") or
                    mca_lower.startswith("pre t cell") or
                    mca_lower.startswith("dpt cell") or
                    mca_lower.startswith("dividing t cell")):
                    matched_labels.append(mca_label)
                continue

            for kw in keywords:
                if kw in mca_lower:
                    matched_labels.append(mca_label)
                    break

        if matched_labels:
            print(f"\n  {our_type}:")
            for m in matched_labels:
                mapping[m] = our_type
                already_mapped.add(m)
                print(f"    → '{m}'")

    # Check for unmapped types
    unmapped_ours = set(OUR_35_TYPES) - set(mapping.values())
    print(f"\n  UNMAPPED types ({len(unmapped_ours)}):")
    for t in sorted(unmapped_ours):
        print(f"    ✗ {t}")

    # Flag ambiguous mappings
    for our_type in OUR_35_TYPES:
        mapped_labels = [k for k, v in mapping.items() if v == our_type]
        if len(mapped_labels) > 5:
            review_needed.append({
                "type": our_type,
                "n_labels": len(mapped_labels),
                "sample_labels": mapped_labels[:5],
            })

    return mapping, review_needed, unmapped_ours


def plot_coverage(coverage_df, output_path):
    """
    Bar chart showing cell count per mapped type, colored by PASS/BORDERLINE/ABSENT.
    """
    fig, ax = plt.subplots(figsize=(14, 10))

    df = coverage_df.sort_values("cell_count", ascending=True).reset_index(drop=True)

    colors = {"PASS": "#2ecc71", "BORDERLINE": "#f39c12", "ABSENT": "#e74c3c"}
    bar_colors = [colors[s] for s in df["status"]]

    ax.barh(range(len(df)), df["cell_count"], color=bar_colors, edgecolor="white", linewidth=0.5)

    # 500-cell gate line
    ax.axvline(x=500, color="black", linestyle="--", linewidth=1.5, label="500-cell gate")

    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["our_type"], fontsize=9)
    ax.set_xlabel("Cell Count in MCA", fontsize=12)
    ax.set_title("MCA Coverage of Our 35 Cell Types\n(T1-A Independent Replication Feasibility)",
                 fontsize=14, fontweight="bold")

    legend_elements = [
        Patch(facecolor="#2ecc71", label="PASS (≥500)"),
        Patch(facecolor="#f39c12", label="BORDERLINE (200-499)"),
        Patch(facecolor="#e74c3c", label="ABSENT (<200)"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=10)

    # Add count labels on bars
    max_count = df["cell_count"].max()
    for i, (count, status) in enumerate(zip(df["cell_count"], df["status"])):
        if count > 0:
            ax.text(count + max_count * 0.02, i, f"{count:,}",
                    va="center", fontsize=7)
        else:
            ax.text(5, i, "0", va="center", fontsize=7, color="#e74c3c")

    ax.set_xscale("symlog", linthresh=10)
    ax.set_xlim(-5, max_count * 3)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(
        lambda x, _: f"{int(x):,}" if x >= 1 else "0"))

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nCoverage bar chart saved to: {output_path}")


def main():
    PYTHON = sys.executable
    print("=" * 70)
    print("T1-A MCA FEASIBILITY CHECK")
    print("Mouse Cell Atlas (Han et al. 2018) — Metadata-Only Query")
    print("Data: MCA_CellAssignments.csv from figshare (dataset 5435866)")
    print("GEO: GSE108097")
    print("=" * 70)

    # ══════════════════════════════════════════════════════════════════════
    # STEP 1: Load MCA cell annotation metadata
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n{'─'*70}")
    print("STEP 1: LOAD MCA CELL ANNOTATION METADATA")
    print(f"{'─'*70}")

    if not os.path.exists(MCA_CSV):
        print(f"ERROR: {MCA_CSV} not found. Download from figshare first.")
        sys.exit(1)

    df = pd.read_csv(MCA_CSV)
    file_size_mb = os.path.getsize(MCA_CSV) / 1024 / 1024

    print(f"  File: {MCA_CSV}")
    print(f"  Size: {file_size_mb:.1f} MB")
    print(f"  Shape: {df.shape[0]:,} rows × {df.shape[1]} columns")
    print(f"  Columns: {list(df.columns)}")
    print(f"\n  First 5 rows:")
    print(df.head().to_string(index=False))

    # ══════════════════════════════════════════════════════════════════════
    # STEP 2: Inspect cell type annotation column
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n{'─'*70}")
    print("STEP 2: CELL TYPE ANNOTATION INSPECTION")
    print(f"{'─'*70}")

    annot_col = "Annotation"
    total_cells = len(df)
    unique_raw_labels = df[annot_col].nunique()

    print(f"  Annotation column: '{annot_col}'")
    print(f"  Total cells: {total_cells:,}")
    print(f"  Unique raw annotation labels: {unique_raw_labels}")

    # Extract base cell type (strip tissue suffix)
    df["base_type"] = df[annot_col].apply(strip_tissue_suffix)
    unique_base_labels = df["base_type"].nunique()

    print(f"  Unique base cell types (tissue-stripped): {unique_base_labels}")

    print(f"\n  Tissues represented: {df['Tissue'].nunique()}")
    print(f"  Tissue list: {sorted(df['Tissue'].unique())}")

    print(f"\n  Top 50 base cell types by count:")
    vc = df["base_type"].value_counts()
    for i, (ct, count) in enumerate(vc.head(50).items()):
        print(f"    {i+1:3d}. {count:>7,}  {ct}")

    if len(vc) > 50:
        print(f"    ... and {len(vc) - 50} more types")

    # ══════════════════════════════════════════════════════════════════════
    # STEP 3: Map MCA labels to our 35-type set
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n{'─'*70}")
    print("STEP 3: MAP MCA LABELS TO OUR 35-TYPE SET")
    print(f"{'─'*70}")

    mca_base_labels = df["base_type"].unique().tolist()
    mapping, review_needed, unmapped = map_types(mca_base_labels, OUR_35_TYPES)

    print(f"\n  Total MCA base labels mapped: {len(mapping)} / {len(mca_base_labels)}")
    print(f"  Our types with mappings: {len(OUR_35_TYPES) - len(unmapped)} / {len(OUR_35_TYPES)}")

    # ══════════════════════════════════════════════════════════════════════
    # STEP 4: Count cells per mapped type
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n{'─'*70}")
    print("STEP 4: CELL COUNTS PER MAPPED TYPE")
    print(f"{'─'*70}")

    # Map each cell to our type
    df["our_type"] = df["base_type"].map(mapping)

    results = []
    for our_type in OUR_35_TYPES:
        mca_labels_mapped = sorted([k for k, v in mapping.items() if v == our_type])
        count = int((df["our_type"] == our_type).sum())

        if count >= 500:
            status = "PASS"
        elif count >= 200:
            status = "BORDERLINE"
        else:
            status = "ABSENT"

        results.append({
            "our_type": our_type,
            "mca_labels": "; ".join(mca_labels_mapped) if mca_labels_mapped else "NO MAPPING",
            "n_mca_labels": len(mca_labels_mapped),
            "cell_count": count,
            "status": status,
        })

    coverage_df = pd.DataFrame(results)

    # Print detailed results
    print(f"\n{'our_type':<50} {'count':>8} {'status':<12} {'mca_labels'}")
    print("─" * 120)
    for _, row in coverage_df.sort_values("cell_count", ascending=False).iterrows():
        mca_short = row["mca_labels"]
        if len(mca_short) > 55:
            mca_short = mca_short[:52] + "..."
        print(f"{row['our_type']:<50} {row['cell_count']:>8,} {row['status']:<12} {mca_short}")

    # How many cells mapped vs unmapped
    mapped_cells = df["our_type"].notna().sum()
    unmapped_cells = df["our_type"].isna().sum()
    print(f"\n  Cells mapped to our types: {mapped_cells:,} / {total_cells:,} ({100*mapped_cells/total_cells:.1f}%)")
    print(f"  Cells unmapped: {unmapped_cells:,} ({100*unmapped_cells/total_cells:.1f}%)")

    # Show what unmapped cells are
    if unmapped_cells > 0:
        print(f"\n  Top 20 UNMAPPED MCA base types (not in our 35):")
        unmapped_df = df[df["our_type"].isna()]
        unmapped_vc = unmapped_df["base_type"].value_counts()
        for ct, count in unmapped_vc.head(20).items():
            print(f"    {count:>7,}  {ct}")

    # ══════════════════════════════════════════════════════════════════════
    # STEP 5: Coverage summary
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n{'─'*70}")
    print("STEP 5: COVERAGE SUMMARY")
    print(f"{'─'*70}")

    n_pass = int((coverage_df["status"] == "PASS").sum())
    n_borderline = int((coverage_df["status"] == "BORDERLINE").sum())
    n_absent = int((coverage_df["status"] == "ABSENT").sum())

    print(f"\n  PASS (≥500 cells):      {n_pass}/35")
    print(f"  BORDERLINE (200-499):   {n_borderline}/35")
    print(f"  ABSENT (<200):          {n_absent}/35")

    # Check original 6
    print(f"\n  Original 6 Phase 2 types coverage:")
    orig6_coverage = {}
    all_orig6_covered = True
    for ct in ORIGINAL_6:
        row = coverage_df[coverage_df["our_type"] == ct].iloc[0]
        status = row["status"]
        count = int(row["cell_count"])
        orig6_coverage[ct] = {"count": count, "status": status}
        marker = "✓" if status == "PASS" else ("~" if status == "BORDERLINE" else "✗")
        print(f"    {marker} {ct}: {count:,} cells ({status})")
        if status == "ABSENT":
            all_orig6_covered = False

    # PASS types list
    pass_types = coverage_df[coverage_df["status"] == "PASS"]["our_type"].tolist()
    borderline_types = coverage_df[coverage_df["status"] == "BORDERLINE"]["our_type"].tolist()
    absent_types = coverage_df[coverage_df["status"] == "ABSENT"]["our_type"].tolist()

    if pass_types:
        print(f"\n  PASS types ({len(pass_types)}):")
        for t in sorted(pass_types):
            c = int(coverage_df[coverage_df["our_type"] == t]["cell_count"].iloc[0])
            print(f"    ✓ {t} ({c:,})")

    if borderline_types:
        print(f"\n  BORDERLINE types ({len(borderline_types)}):")
        for t in sorted(borderline_types):
            c = int(coverage_df[coverage_df["our_type"] == t]["cell_count"].iloc[0])
            print(f"    ~ {t} ({c:,})")

    if absent_types:
        print(f"\n  ABSENT types ({len(absent_types)}):")
        for t in sorted(absent_types):
            c = int(coverage_df[coverage_df["our_type"] == t]["cell_count"].iloc[0])
            print(f"    ✗ {t} ({c:,})")

    # Feasibility verdict
    if n_pass >= 20 and all_orig6_covered:
        verdict = "FEASIBLE"
    elif n_pass >= 10 or (all_orig6_covered and n_pass + n_borderline >= 20):
        verdict = "PARTIAL"
    else:
        verdict = "INFEASIBLE"

    print(f"\n  ┌──────────────────────────────────────────────────┐")
    print(f"  │  FEASIBILITY VERDICT: {verdict:<27} │")
    print(f"  │  {n_pass} PASS / {n_borderline} BORDERLINE / {n_absent} ABSENT          │")
    print(f"  │  Original 6 covered: {'YES' if all_orig6_covered else 'NO':<27} │")
    print(f"  └──────────────────────────────────────────────────┘")

    if verdict == "PARTIAL":
        print(f"\n  RECOMMENDATION: Use reduced type set for T1-A.")
        print(f"  Suggested set: {len(pass_types)} PASS types")
        if borderline_types:
            print(f"  + {len(borderline_types)} BORDERLINE types if gate relaxed to 200")
    elif verdict == "INFEASIBLE":
        print(f"\n  RECOMMENDATION: MCA insufficient for T1-A.")
        print(f"  Alternatives: Human Cell Atlas for human-side replication,")
        print(f"  or aggregate MCA with additional mouse datasets.")
    else:
        print(f"\n  RECOMMENDATION: Proceed with MCA as T1-A independent replication dataset.")
        if absent_types:
            print(f"  Note: {len(absent_types)} types absent — use {n_pass}-type reduced set.")

    # ══════════════════════════════════════════════════════════════════════
    # STEP 6: Coverage bar chart
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n{'─'*70}")
    print("STEP 6: COVERAGE BAR CHART")
    print(f"{'─'*70}")

    plot_path = os.path.join(OUTPUT_DIR, "coverage_bar.png")
    plot_coverage(coverage_df, plot_path)

    print("Plot description: Horizontal bar chart showing cell count for each of our")
    print("35 types as available in MCA. Bars colored green (PASS, ≥500), orange")
    print("(BORDERLINE, 200-499), or red (ABSENT, <200). Dashed black vertical line")
    print("marks the 500-cell gate. X-axis is log-scaled.")

    # ══════════════════════════════════════════════════════════════════════
    # STEP 7: Technology and pipeline independence
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n{'─'*70}")
    print("STEP 7: TECHNOLOGY AND PIPELINE INDEPENDENCE")
    print(f"{'─'*70}")

    # From GEO series matrix and paper metadata (verified from GEO download):
    # Series contributors: Han X, Wang R
    # Contact: Guoji Guo, Zhejiang University School of Medicine
    # Platform: GPL17021 (Illumina HiSeq 2500, Mus musculus)
    # Protocol: Microwell-seq (custom platform developed by Guo lab)
    # Alignment: STAR → Drop-seq_tools pipeline → UMI count matrix

    tabula_key_authors = [
        "quake", "pisco", "darmanis", "krasnow", "weissman",
        "schaum", "wyss-coray", "almanzar", "tabula",
    ]

    # MCA authors from GEO: Han X, Wang R, contact: Guoji Guo
    # Full author list from paper: Han X, Wang R, Zhou Y, Fei L, Sun H, Lai S,
    # Saadatpour A, Zhou Z, Chen H, Ye F, Huang D, Xu Y, Huang W, Jiang M,
    # Jiang X, Mao J, Chen Y, Lu C, Xie J, Fang Q, Wang Y, Yue R, Li T,
    # Huang H, Orkin SH, Yuan GC, Chen M, Guo G
    mca_authors = [
        "han", "wang", "zhou", "fei", "sun", "lai", "saadatpour",
        "chen", "ye", "huang", "xu", "jiang", "mao", "lu", "xie",
        "fang", "yue", "li", "orkin", "yuan", "guo",
    ]

    overlaps = []
    for ta in tabula_key_authors:
        for ma in mca_authors:
            if ta in ma or ma in ta:
                # Only flag if it's a meaningful match (not just "wang" matching "wang")
                if ta == ma:
                    overlaps.append((ta, ma))

    print(f"\n  Technology: Microwell-seq (custom droplet-free microwell platform)")
    print(f"  Sequencing: Illumina HiSeq 2500")
    print(f"  Institution: Zhejiang University School of Medicine / BGI")
    print(f"  PI: Guoji Guo (ggjlab)")
    print(f"  Alignment pipeline: STAR + Drop-seq_tools (NOT Cell Ranger)")
    print(f"  Genome: Mus_musculus.GRCm38.88")

    print(f"\n  Tabula consortium key authors: {tabula_key_authors}")
    if overlaps:
        print(f"  ⚠ Potential author overlap: {overlaps}")
        independent = False
    else:
        print(f"  ✓ NO overlap with Tabula consortium key authors")
        independent = True

    print(f"\n  Independence assessment:")
    print(f"    Technology:  Microwell-seq vs 10x Chromium/Smart-seq2  → INDEPENDENT")
    print(f"    Institution: Zhejiang University vs CZI/Stanford       → INDEPENDENT")
    print(f"    Pipeline:    STAR+Drop-seq_tools vs Cell Ranger        → INDEPENDENT")
    print(f"    Authors:     {'NO overlap' if independent else 'OVERLAP DETECTED'}  → {'INDEPENDENT' if independent else 'CHECK REQUIRED'}")
    print(f"    Overall:     {'FULLY INDEPENDENT' if independent else 'PARTIALLY INDEPENDENT'}")

    # ══════════════════════════════════════════════════════════════════════
    # Save outputs
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n{'─'*70}")
    print("SAVING OUTPUTS")
    print(f"{'─'*70}")

    # Coverage table
    coverage_path = os.path.join(OUTPUT_DIR, "mca_coverage_table.csv")
    coverage_df.to_csv(coverage_path, index=False)
    print(f"  Coverage table: {coverage_path}")

    # Type mapping
    mapping_path = os.path.join(OUTPUT_DIR, "mca_type_mapping.json")
    with open(mapping_path, "w") as f:
        json.dump(mapping, f, indent=2)
    print(f"  Type mapping: {mapping_path}")

    # Summary report JSON
    report = {
        "total_mca_cells": int(total_cells),
        "unique_mca_raw_labels": int(unique_raw_labels),
        "unique_mca_base_labels": int(unique_base_labels),
        "tissues_represented": int(df["Tissue"].nunique()),
        "types_pass": n_pass,
        "types_borderline": n_borderline,
        "types_absent": n_absent,
        "types_pass_list": pass_types,
        "types_borderline_list": borderline_types,
        "types_absent_list": absent_types,
        "original_6_coverage": {k: v for k, v in orig6_coverage.items()},
        "all_original_6_covered": all_orig6_covered,
        "feasibility_verdict": verdict,
        "technology": "Microwell-seq (Illumina HiSeq 2500)",
        "institution": "Zhejiang University School of Medicine / BGI",
        "independent_of_tabula": independent,
        "review_needed": review_needed,
        "data_source": "MCA_CellAssignments.csv from figshare (dataset 5435866)",
        "geo_accession": "GSE108097",
        "paper": "Han et al. Cell 172, 1091-1107 (2018)",
        "doi": "10.1016/j.cell.2018.02.001",
        "cells_mapped_to_our_types": int(mapped_cells),
        "cells_unmapped": int(unmapped_cells),
    }

    report_path = os.path.join(OUTPUT_DIR, "mca_feasibility_report.json")
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)
    print(f"  Feasibility report: {report_path}")

    # Detailed mapping table (for manual review)
    detail_path = os.path.join(OUTPUT_DIR, "mca_mapping_detail.csv")
    detail_rows = []
    for mca_label, our_type in sorted(mapping.items(), key=lambda x: x[1]):
        count = int((df["base_type"] == mca_label).sum())
        detail_rows.append({
            "mca_base_label": mca_label,
            "our_type": our_type,
            "cell_count": count,
        })
    pd.DataFrame(detail_rows).to_csv(detail_path, index=False)
    print(f"  Mapping detail: {detail_path}")

    print(f"\n{'='*70}")
    print(f"MCA FEASIBILITY CHECK COMPLETE")
    print(f"Verdict: {verdict}")
    print(f"{'='*70}")

    return verdict, coverage_df, report


if __name__ == "__main__":
    verdict, coverage_df, report = main()
