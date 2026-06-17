#!/usr/bin/env python3
"""
CellWarp — Direction 1: Cancer as Geometric Deformation — Data Download

Downloads colorectal cancer (CRC) and matched normal colon tissue from CZ CELLxGENE
Census, harmonizes cell type labels to coarse categories, and saves normalized
AnnData files for downstream Procrustes analysis.

Biology
-------
We apply the same geometric framework used for cross-species comparison to the
normal-to-cancer transformation. If cancer reshapes cell type expression programs via
a coherent geometric deformation (analogous to how evolution reshapes them across
species), the Procrustes residuals should point to known cancer biology.

Pipeline
--------
1. Census inventory: discover tissue, disease, and cell type values
2. Build coarse cell type mapping (10 categories)
3. **PAUSE** — user reviews mapping before download
4. Download cells (≤3,000 per coarse type per condition)
5. Filter to shared ortholog gene space (16,959 genes)
6. Normalize (counts per 10k + log1p)
7. Save: data/cancer/colon_normal.h5ad, data/cancer/colon_tumor.h5ad

Checkpoints
-----------
    output/cancer/census_inventory.txt  — Census inventory (tissue, disease, cell types)
    data/cancer/colon_normal.h5ad       — Normalized normal colon data
    data/cancer/colon_tumor.h5ad        — Normalized CRC tumor data

Usage:
    python scripts/10_cancer_download.py

    The script will STOP after printing the coarse cell type mapping and wait for
    user confirmation (y/n) before proceeding with the download.

Expected runtime:
    5–15 minutes (network-dependent). Census queries ~2-5 min, download ~5-10 min.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# Add project root to path for src imports
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import cellxgene_census

from cellwarp.cancer_loader import (
    DATA_DIR,
    MAX_CELLS_PER_TYPE,
    MIN_CELLS_PER_TYPE,
    OBS_COLUMNS,
    OUTPUT_DIR,
    VAR_COLUMNS,
    build_coarse_mapping,
    download_condition,
    load_ortholog_gene_ids,
    normalize_cancer_data,
    print_coarse_mapping_table,
    print_final_summary,
    query_cell_type_inventory,
    query_disease_values,
    query_tissue_values,
    save_census_inventory,
    save_h5ad_atomic,
)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Keywords to search for colon/intestine tissue in Census (broad — for inventory display)
TISSUE_KEYWORDS = ["colon", "large intestine", "intestine", "colorectal"]

# Keywords to filter tissue_general values for download (restrictive — must match
# the tissue_general name itself, not just an associated tissue name).
# Only "colon" and "large intestine" are valid per the task spec.  "intestine"
# alone is too broad (includes small intestine), and spurious tissue_general values
# like "pancreas"/"placenta" appear because their tissue column cross-references
# colon anatomy.
TISSUE_GENERAL_KEYWORDS = ["colon", "large intestine"]

# Path to the shared ortholog table from Phase 1
ORTHOLOG_PATH = Path("./data/phase1/orthologs_human_mouse.csv")

# Output paths
INVENTORY_PATH = OUTPUT_DIR / "census_inventory.txt"
NORMAL_PATH = DATA_DIR / "colon_normal.h5ad"
TUMOR_PATH = DATA_DIR / "colon_tumor.h5ad"


def main() -> None:
    parser = argparse.ArgumentParser(description="Download CRC and normal colon data")
    parser.add_argument(
        "--auto-confirm", action="store_true",
        help="Skip interactive mapping confirmation (for non-interactive runs)",
    )
    args = parser.parse_args()
    auto_confirm = args.auto_confirm

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # ==================================================================
    # STEP 1: Load ortholog gene IDs (shared gene space from Phase 1)
    # ==================================================================
    print("=" * 70)
    print("STEP 1/7: Loading shared ortholog gene space")
    print("=" * 70)

    if not ORTHOLOG_PATH.exists():
        print(f"  ERROR: Ortholog table not found at {ORTHOLOG_PATH}")
        print("  Run scripts/01_download_data.py first to generate it.")
        sys.exit(1)

    ortholog_gene_ids = load_ortholog_gene_ids(ORTHOLOG_PATH)

    # ==================================================================
    # STEP 2: Census inventory — tissue, disease, cell types
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 2/7: Querying CELLxGENE Census inventory")
    print("=" * 70)

    print("  Opening Census connection (version='2025-11-08')...")
    with cellxgene_census.open_soma(census_version="2025-11-08") as census:
        # -- Tissue values --
        print("\n  Querying tissue values for colon/intestine...")
        tissue_df = query_tissue_values(census, tissue_keywords=TISSUE_KEYWORDS)

        print(f"\n  Found {len(tissue_df)} tissue entries:")
        print(f"  {'tissue_general':<30} {'tissue':<35} {'cells':>10}")
        print("  " + "-" * 78)
        for _, row in tissue_df.iterrows():
            print(
                f"  {str(row['tissue_general']):<30} "
                f"{str(row['tissue']):<35} "
                f"{row['cell_count']:>10,}"
            )

        # Determine the tissue_general values to use for download.
        # IMPORTANT: filter on tissue_general name itself, not the tissue column.
        # The inventory matches broadly (tissue OR tissue_general contains keyword),
        # which pulls in spurious tissue_general values like "pancreas" or "placenta"
        # whose tissue column cross-references colon anatomy.
        all_tg = tissue_df["tissue_general"].unique()
        tissue_general_values = sorted([
            tg for tg in all_tg
            if any(kw in tg.lower() for kw in TISSUE_GENERAL_KEYWORDS)
        ])
        if not tissue_general_values:
            print("\n  ERROR: No tissue_general values matched TISSUE_GENERAL_KEYWORDS.")
            print(f"  Keywords: {TISSUE_GENERAL_KEYWORDS}")
            print(f"  Available: {sorted(all_tg.tolist())}")
            sys.exit(1)
        print(f"\n  Using tissue_general values: {tissue_general_values}")
        print(f"  (filtered from {len(all_tg)} inventory entries using {TISSUE_GENERAL_KEYWORDS})")

        # -- Disease values --
        print(f"\n  Querying disease values for these tissues...")
        disease_df = query_disease_values(census, tissue_general_values)

        print(f"\n  Found {len(disease_df)} disease values:")
        print(f"  {'disease':<50} {'cells':>10}")
        print("  " + "-" * 62)
        for _, row in disease_df.iterrows():
            print(f"  {str(row['disease']):<50} {row['cell_count']:>10,}")

        # Disease labels — explicit rather than auto-detected.
        # "colon adenocarcinoma" has 4.5x more cells than "colorectal cancer"
        # and better coverage of the full tumor microenvironment (epithelial +
        # stromal + immune compartments). See DECISION-049.
        disease_normal = "normal"
        disease_tumor = "colon adenocarcinoma"

        # Verify both labels exist in Census
        for label in [disease_normal, disease_tumor]:
            if label not in disease_df["disease"].values:
                print(f"\n  ERROR: disease='{label}' not found in Census.")
                print("  Available values printed above.")
                sys.exit(1)

        tumor_count = disease_df.loc[
            disease_df["disease"] == disease_tumor, "cell_count"
        ].iloc[0]
        print(f"\n  Selected disease labels:")
        print(f"    Normal: '{disease_normal}'")
        print(f"    Tumor:  '{disease_tumor}' ({tumor_count:,} cells)")

        # -- Cell type inventory --
        print(f"\n  Querying cell types for normal vs tumor...")
        cell_type_df = query_cell_type_inventory(
            census,
            tissue_general_values,
            disease_normal,
            disease_tumor,
        )

        print(f"\n  Found {len(cell_type_df)} distinct cell type labels")
        print(f"  {'cell_type':<50} {'normal':>8} {'tumor':>8}")
        print("  " + "-" * 68)
        for _, row in cell_type_df.iterrows():
            if row["normal_count"] > 0 or row["tumor_count"] > 0:
                print(
                    f"  {str(row['cell_type']):<50} "
                    f"{row['normal_count']:>8,} "
                    f"{row['tumor_count']:>8,}"
                )

        # Save inventory
        save_census_inventory(tissue_df, disease_df, cell_type_df, INVENTORY_PATH)

        # ==============================================================
        # STEP 3: Build and print coarse cell type mapping
        # ==============================================================
        print("\n" + "=" * 70)
        print("STEP 3/7: Building coarse cell type mapping")
        print("=" * 70)

        mapping_df = build_coarse_mapping(cell_type_df)
        agg_df = print_coarse_mapping_table(mapping_df)

        # Determine which types pass the gate
        passing_types = agg_df.loc[agg_df["passes_gate"], "coarse_label"].tolist()
        failing_types = agg_df.loc[~agg_df["passes_gate"], "coarse_label"].tolist()

        print(f"\n  Types passing >=500 gate: {passing_types}")
        if failing_types:
            print(f"  Types FAILING >=500 gate (will be dropped): {failing_types}")

        # ==============================================================
        # STEP 4: PAUSE — wait for user confirmation
        # ==============================================================
        print("\n" + "=" * 70)
        print("REVIEW REQUIRED")
        print("=" * 70)
        print(
            "  Please review the coarse cell type mapping above.\n"
            "  The download will include only types passing the >=500 gate.\n"
        )

        if auto_confirm:
            print("  --auto-confirm: skipping interactive prompt.")
        else:
            response = input("  Proceed with download? [y/N]: ").strip().lower()
            if response != "y":
                print("  Aborted by user. No data downloaded.")
                sys.exit(0)

        print("\n  Proceeding with download...")

        # Build raw → coarse mapping dict for all raw labels
        coarse_map = dict(
            zip(mapping_df["raw_label"], mapping_df["coarse_label"])
        )

        # ==============================================================
        # STEP 5: Download normal colon data
        # ==============================================================
        print("\n" + "=" * 70)
        print("STEP 5/7: Downloading normal colon data")
        print("=" * 70)

        normal_adata = download_condition(
            census=census,
            tissue_general_values=tissue_general_values,
            disease_label=disease_normal,
            coarse_mapping=coarse_map,
            valid_coarse_types=passing_types,
            ortholog_gene_ids=ortholog_gene_ids,
            max_cells_per_type=MAX_CELLS_PER_TYPE,
        )
        print(f"  Normal data: {normal_adata.n_obs:,} cells x {normal_adata.n_vars:,} genes")

        # ==============================================================
        # STEP 6: Download tumor data
        # ==============================================================
        print("\n" + "=" * 70)
        print("STEP 6/7: Downloading CRC tumor data")
        print("=" * 70)

        tumor_adata = download_condition(
            census=census,
            tissue_general_values=tissue_general_values,
            disease_label=disease_tumor,
            coarse_mapping=coarse_map,
            valid_coarse_types=passing_types,
            ortholog_gene_ids=ortholog_gene_ids,
            max_cells_per_type=MAX_CELLS_PER_TYPE,
        )
        print(f"  Tumor data: {tumor_adata.n_obs:,} cells x {tumor_adata.n_vars:,} genes")

    # Census connection closed
    print("\n  Census connection closed.")

    # ==================================================================
    # STEP 7: Normalize and save
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 7/7: Normalizing and saving")
    print("=" * 70)

    print("\n  Normalizing normal colon data...")
    normal_adata = normalize_cancer_data(normal_adata)

    print("\n  Normalizing tumor data...")
    tumor_adata = normalize_cancer_data(tumor_adata)

    # Save
    print()
    save_h5ad_atomic(normal_adata, NORMAL_PATH)
    save_h5ad_atomic(tumor_adata, TUMOR_PATH)

    # ==================================================================
    # SUMMARY
    # ==================================================================
    summary = print_final_summary(normal_adata, tumor_adata, failing_types)

    # Save summary as JSON for programmatic access
    summary_path = OUTPUT_DIR / "download_summary.json"
    # Convert non-serializable types
    serializable = {
        k: v if not isinstance(v, list) or not any(hasattr(x, "__iter__") and not isinstance(x, str) for x in v) else [str(x) for x in v]
        for k, v in summary.items()
    }
    with open(summary_path, "w") as f:
        json.dump(serializable, f, indent=2, default=str)
    print(f"\n  Saved summary: {summary_path}")

    print("\n" + "=" * 70)
    print("CANCER DATA DOWNLOAD COMPLETE")
    print("=" * 70)
    print(f"  Normal: {NORMAL_PATH} ({normal_adata.n_obs:,} cells)")
    print(f"  Tumor:  {TUMOR_PATH} ({tumor_adata.n_obs:,} cells)")
    print(f"  Genes:  {normal_adata.n_vars:,} (shared ortholog space)")
    print(f"\n  Next step: Compute centroids and run Procrustes on normal vs tumor.")


if __name__ == "__main__":
    main()
