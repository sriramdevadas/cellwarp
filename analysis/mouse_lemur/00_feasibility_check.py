"""
Mouse Lemur (Microcebus murinus) Feasibility Check
====================================================
Before running the full CellWarp Procrustes pipeline on human-mouse lemur,
we check three prerequisites:

1. Data availability: Is Tabula Microcebus in CELLxGENE Census (build 2025-11-08)?
   If not, check alternative sources.
2. Cell type overlap: How many of the 35 primary analysis types have ≥500 cells
   in the mouse lemur data? Gate: ≥15 types.
3. Ortholog depth: How many human-mouse lemur 1:1 orthologs exist? Gate: ≥12,000.

Biology
-------
Microcebus murinus (gray mouse lemur) diverged from humans ~75 Mya — intermediate
between human-macaque (~25-30 Mya) and human-mouse (~90 Mya). Tabula Microcebus
(The Tabula Microcebus Consortium, 2022) is a multi-tissue single-cell atlas
from the same consortium lineage as Tabula Sapiens and Tabula Muris Senis.

If feasibility passes, this adds a fourth species to the CellWarp framework,
filling the evolutionary distance gap between macaque and mouse.

Math
----
No computation beyond counting. Gate criteria are binary: pass or fail.
"""

import json
import sys
from pathlib import Path

import cellxgene_census
import pandas as pd

# Add src to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
OUT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_ROOT / "data" / "phase1"

# The 35 primary analysis cell types (from residuals_ranked.csv)
PRIMARY_35_TYPES = [
    "stromal cell",
    "epithelial cell",
    "hematopoietic precursor cell",
    "hematopoietic stem cell",
    "pancreatic acinar cell",
    "basal cell",
    "T cell",
    "neutrophil",
    "fibroblast of cardiac tissue",
    "myeloid leukocyte",
    "mesenchymal stem cell of adipose tissue",
    "plasma cell",
    "mesenchymal stem cell",
    "CD4-positive, alpha-beta T cell",
    "classical monocyte",
    "macrophage",
    "B cell",
    "luminal epithelial cell of mammary gland",
    "large intestine goblet cell",
    "enterocyte of epithelium of large intestine",
    "myeloid dendritic cell",
    "monocyte",
    "natural killer cell",
    "intermediate monocyte",
    "mature NK T cell",
    "adventitial cell",
    "granulocyte",
    "fibroblast",
    "bladder urothelial cell",
    "pancreatic ductal cell",
    "smooth muscle cell",
    "hepatocyte",
    "endothelial cell",
    "non-classical monocyte",
    "CD8-positive, alpha-beta T cell",
]

MOUSE_LEMUR_ORGANISM = "Microcebus murinus"
MIN_CELLS_PER_TYPE = 500
MIN_SHARED_TYPES = 15
MIN_SHARED_GENES = 12_000

CENSUS_VERSION = "2025-11-08"


def check_census_availability(census):
    """Check if mouse lemur data exists in CELLxGENE Census."""
    print("\n--- Check 1: Census data availability ---")

    # Check organisms available
    datasets_df = census["census_info"]["datasets"].read().concat().to_pandas()

    # Look for mouse lemur datasets
    lemur_mask = datasets_df["organism"].apply(
        lambda x: "microcebus" in str(x).lower() or "lemur" in str(x).lower()
        if pd.notna(x) else False
    )
    lemur_datasets = datasets_df[lemur_mask]

    if len(lemur_datasets) == 0:
        # Try via organism column variations
        print("  No datasets with 'microcebus' or 'lemur' in organism column.")
        print("  Checking all unique organisms...")
        organisms = datasets_df["organism"].dropna().unique()
        primate_orgs = [o for o in organisms if any(
            k in o.lower() for k in ["primate", "lemur", "microcebus", "strep"]
        )]
        if primate_orgs:
            print(f"  Primate-related organisms found: {primate_orgs}")
        else:
            print("  No primate-related organisms beyond standard (human/macaque).")
            # Print all unique organisms for manual inspection
            print(f"  All organisms ({len(organisms)} total):")
            for org in sorted(organisms):
                print(f"    - {org}")
    else:
        print(f"  Found {len(lemur_datasets)} mouse lemur dataset(s):")
        for _, row in lemur_datasets.iterrows():
            print(f"    Collection: {row.get('collection_name', 'N/A')}")
            print(f"    Title: {row.get('dataset_title', 'N/A')}")
            print(f"    Cells: {row.get('cell_count', 'N/A')}")
            print(f"    Organism: {row.get('organism', 'N/A')}")
            print()

    # Also try querying obs directly for Microcebus murinus
    print("  Attempting direct obs query for 'Microcebus murinus'...")
    try:
        obs_df = cellxgene_census.get_obs(
            census,
            MOUSE_LEMUR_ORGANISM,
            value_filter="is_primary_data == True",
            column_names=["cell_type", "tissue", "tissue_general", "assay",
                          "dataset_id", "donor_id", "disease"],
        )
        print(f"  SUCCESS: {len(obs_df):,} cells found for {MOUSE_LEMUR_ORGANISM}")
        return obs_df, datasets_df
    except Exception as e:
        print(f"  Query failed: {e}")
        # Try alternative organism names
        alt_names = ["Microcebus murinus", "microcebus murinus"]
        for name in alt_names:
            try:
                obs_df = cellxgene_census.get_obs(
                    census,
                    name,
                    value_filter="is_primary_data == True",
                    column_names=["cell_type", "tissue", "assay", "dataset_id",
                                  "donor_id", "disease"],
                )
                if len(obs_df) > 0:
                    print(f"  Found {len(obs_df):,} cells with organism='{name}'")
                    return obs_df, datasets_df
            except Exception:
                continue

        print("  Mouse lemur NOT found in Census under any tried name.")
        return None, datasets_df


def analyze_cell_types(obs_df):
    """Analyze cell type overlap with the 35 primary types."""
    print("\n--- Check 2: Cell type overlap ---")

    # All cell types in mouse lemur data
    all_types = obs_df["cell_type"].value_counts()
    print(f"  Total unique cell types in mouse lemur: {len(all_types)}")

    # Check overlap with primary 35
    overlap_types = []
    missing_types = []

    print(f"\n  {'Cell Type':<55} {'Lemur Cells':>12} {'Status':>10}")
    print(f"  {'-' * 80}")

    for ct in PRIMARY_35_TYPES:
        if ct in all_types.index:
            n = all_types[ct]
            status = "PASS" if n >= MIN_CELLS_PER_TYPE else "LOW"
            if n >= MIN_CELLS_PER_TYPE:
                overlap_types.append((ct, n))
            else:
                overlap_types.append((ct, n))
            print(f"  {ct:<55} {n:>12,} {status:>10}")
        else:
            missing_types.append(ct)
            print(f"  {ct:<55} {'---':>12} {'MISSING':>10}")

    passing = [(ct, n) for ct, n in overlap_types if n >= MIN_CELLS_PER_TYPE]
    low = [(ct, n) for ct, n in overlap_types if n < MIN_CELLS_PER_TYPE]

    print(f"\n  Summary:")
    print(f"    Types found (≥500 cells): {len(passing)}")
    print(f"    Types found (<500 cells): {len(low)}")
    print(f"    Types missing entirely:   {len(missing_types)}")
    print(f"    Gate (≥{MIN_SHARED_TYPES} types):       "
          f"{'PASS' if len(passing) >= MIN_SHARED_TYPES else 'FAIL'}")

    # Also print top lemur cell types NOT in our 35
    lemur_only = [ct for ct in all_types.index if ct not in PRIMARY_35_TYPES]
    if lemur_only:
        print(f"\n  Top 20 mouse lemur cell types NOT in primary 35:")
        for ct in lemur_only[:20]:
            print(f"    {ct:<55} {all_types[ct]:>8,}")

    # Tissue breakdown
    print(f"\n  Tissue breakdown:")
    tissues = obs_df["tissue"].value_counts()
    for tissue, n in tissues.items():
        print(f"    {tissue:<40} {n:>8,}")

    # Assay breakdown
    print(f"\n  Assay breakdown:")
    assays = obs_df["assay"].value_counts()
    for assay, n in assays.items():
        print(f"    {assay:<40} {n:>8,}")

    return passing, low, missing_types, all_types


def check_orthologs():
    """Query BioMart for human-mouse lemur 1:1 orthologs."""
    print("\n--- Check 3: Ortholog depth ---")

    cache_path = OUT_DIR / "biomart_mouse_lemur_human_orthologs.csv"

    if cache_path.exists():
        print(f"  Loading cached orthologs from {cache_path}")
        orthologs = pd.read_csv(cache_path)
        print(f"  Loaded {len(orthologs):,} 1:1 ortholog pairs")
        return orthologs

    from pybiomart import Dataset

    # Mouse lemur BioMart code: mmurinus
    species_code = "mmurinus"
    print(f"  Querying Ensembl BioMart for human-mouse lemur orthologs...")
    print(f"  Species: Microcebus murinus ({species_code})")

    attributes = [
        "ensembl_gene_id",
        "external_gene_name",
        f"{species_code}_homolog_ensembl_gene",
        f"{species_code}_homolog_associated_gene_name",
        f"{species_code}_homolog_orthology_type",
    ]

    try:
        dataset = Dataset(name="hsapiens_gene_ensembl", host="http://www.ensembl.org")
        print(f"  Querying with attributes: {attributes}")
        result = dataset.query(attributes=attributes)
        print(f"  Raw query returned: {len(result):,} rows")
        print(f"  Columns: {list(result.columns)}")
    except Exception as e:
        print(f"  BioMart query failed: {e}")
        print("  Trying alternative host (useast)...")
        try:
            dataset = Dataset(
                name="hsapiens_gene_ensembl",
                host="http://useast.ensembl.org",
            )
            result = dataset.query(attributes=attributes)
            print(f"  Raw query returned: {len(result):,} rows")
        except Exception as e2:
            print(f"  Alternative host also failed: {e2}")
            return None

    # Find orthology type column
    orthology_col = None
    for col in result.columns:
        if "homology type" in col.lower() or "orthology_type" in col.lower():
            orthology_col = col
            break

    if orthology_col is None:
        print(f"  ERROR: Could not find orthology type column in: {list(result.columns)}")
        return None

    print(f"  Orthology type column: {orthology_col}")
    print(f"  Orthology type distribution:")
    print(result[orthology_col].value_counts().to_string())

    # Filter to 1:1 orthologs
    one2one = result[result[orthology_col] == "ortholog_one2one"].copy()
    print(f"\n  1:1 orthologs (raw): {len(one2one):,}")

    # Identify columns
    gene_id_col = "Gene stable ID" if "Gene stable ID" in one2one.columns else one2one.columns[0]
    gene_name_col = "Gene name" if "Gene name" in one2one.columns else one2one.columns[1]

    lemur_gene_col = None
    lemur_name_col = None
    for col in one2one.columns:
        col_lower = col.lower()
        if col_lower in ("gene stable id", "gene name"):
            continue
        if "gene stable id" in col_lower:
            lemur_gene_col = col
        elif "gene name" in col_lower and "homology" not in col_lower:
            lemur_name_col = col

    print(f"  Human gene ID column: {gene_id_col}")
    print(f"  Lemur gene ID column: {lemur_gene_col}")

    # Drop missing IDs
    before = len(one2one)
    one2one = one2one.dropna(subset=[gene_id_col])
    if lemur_gene_col:
        one2one = one2one.dropna(subset=[lemur_gene_col])
        one2one = one2one[one2one[lemur_gene_col] != ""]
    after = len(one2one)
    print(f"  After dropping missing: {after:,} (dropped {before - after:,})")

    # Deduplicate by human Ensembl ID
    before = len(one2one)
    one2one = one2one.drop_duplicates(subset=[gene_id_col])
    after = len(one2one)
    print(f"  After deduplication: {after:,} (dropped {before - after:,})")

    # Standardize column names for downstream use
    orthologs = pd.DataFrame({
        "human_ensembl_id": one2one[gene_id_col].values,
        "human_gene_name": one2one[gene_name_col].values if gene_name_col in one2one.columns else "",
        "lemur_ensembl_id": one2one[lemur_gene_col].values if lemur_gene_col else "",
        "lemur_gene_name": one2one[lemur_name_col].values if lemur_name_col and lemur_name_col in one2one.columns else "",
    })

    # Cache
    orthologs.to_csv(cache_path, index=False)
    print(f"  Cached to {cache_path}")

    # Compare with existing human-mouse space
    hm_orthologs = pd.read_csv(DATA_DIR / "orthologs_human_mouse.csv")
    hm_ids = set(hm_orthologs["human_ensembl_id"].dropna())
    lemur_ids = set(orthologs["human_ensembl_id"].dropna())

    three_way = hm_ids & lemur_ids
    print(f"\n  Human-mouse lemur 1:1 orthologs: {len(lemur_ids):,}")
    print(f"  Human-mouse 1:1 orthologs (existing): {len(hm_ids):,}")
    print(f"  Three-way intersection (H-M-ML): {len(three_way):,}")
    print(f"  Overlap %: {100 * len(three_way) / len(hm_ids):.1f}% of H-M space")

    # Also check human-macaque for four-way
    macaque_cache = PROJECT_ROOT / "data" / "macaque" / "biomart_macaque_human_orthologs.csv"
    if macaque_cache.exists():
        mac_orthologs = pd.read_csv(macaque_cache)
        mac_ids = set(mac_orthologs["human_ensembl_id"].dropna())
        four_way = hm_ids & lemur_ids & mac_ids
        print(f"  Human-macaque 1:1 orthologs: {len(mac_ids):,}")
        print(f"  Four-way intersection (H-M-Mac-ML): {len(four_way):,}")

    print(f"\n  Gene gate (≥{MIN_SHARED_GENES:,}): "
          f"{'PASS' if len(lemur_ids) >= MIN_SHARED_GENES else 'FAIL'} "
          f"(pairwise: {len(lemur_ids):,})")
    print(f"  Three-way gate (≥{MIN_SHARED_GENES:,}): "
          f"{'PASS' if len(three_way) >= MIN_SHARED_GENES else 'FAIL'} "
          f"({len(three_way):,})")

    return orthologs


def main():
    print("=" * 70)
    print("MOUSE LEMUR (Microcebus murinus) FEASIBILITY CHECK")
    print("=" * 70)
    print(f"  Census version: {CENSUS_VERSION}")
    print(f"  Primary types: {len(PRIMARY_35_TYPES)}")
    print(f"  Min cell type overlap: {MIN_SHARED_TYPES}")
    print(f"  Min shared genes: {MIN_SHARED_GENES:,}")
    print(f"  Min cells per type: {MIN_CELLS_PER_TYPE:,}")

    # ── Check 1: Census availability ──────────────────────────────────
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        obs_df, datasets_df = check_census_availability(census)

    if obs_df is None or len(obs_df) == 0:
        print("\n" + "=" * 70)
        print("FEASIBILITY RESULT: NO-GO")
        print("Reason: Mouse lemur data not found in CELLxGENE Census.")
        print("Next: Check GEO (GSE148911) and Tabula Microcebus consortium page.")
        print("=" * 70)

        # Write NO-GO report
        report = """# Mouse Lemur Feasibility Check

**Date:** 2026-04-05
**Status:** NO-GO (Census)
**Reason:** Mouse lemur (Microcebus murinus) data not available in CELLxGENE Census build 2025-11-08.

## Census Check
- Organism "Microcebus murinus" not found in Census organism list.
- No datasets with 'microcebus' or 'lemur' in organism metadata.

## Next Steps
1. Check GEO for Tabula Microcebus raw data (GSE148911 or similar)
2. Check the Tabula Microcebus consortium data portal
3. If data is available externally, download and process manually
4. Re-assess feasibility with external data source
"""
        with open(OUT_DIR / "feasibility_check.md", "w") as f:
            f.write(report)
        print(f"\n  Report saved to {OUT_DIR / 'feasibility_check.md'}")
        return False

    # ── Check 2: Cell type overlap ────────────────────────────────────
    passing, low, missing, all_types = analyze_cell_types(obs_df)

    # ── Check 3: Ortholog depth ───────────────────────────────────────
    orthologs = check_orthologs()

    # ── GO / NO-GO Decision ───────────────────────────────────────────
    n_passing = len(passing)
    n_orthologs = len(orthologs) if orthologs is not None else 0

    type_gate = n_passing >= MIN_SHARED_TYPES
    gene_gate = n_orthologs >= MIN_SHARED_GENES

    go = type_gate and gene_gate

    print("\n" + "=" * 70)
    print(f"FEASIBILITY RESULT: {'GO' if go else 'NO-GO'}")
    print("=" * 70)
    print(f"  Cell type gate (≥{MIN_SHARED_TYPES}): {n_passing} types → "
          f"{'PASS' if type_gate else 'FAIL'}")
    print(f"  Ortholog gate (≥{MIN_SHARED_GENES:,}): {n_orthologs:,} genes → "
          f"{'PASS' if gene_gate else 'FAIL'}")
    if not go:
        if not type_gate:
            print(f"  BLOCKED: Only {n_passing} cell types pass ≥500 cell gate "
                  f"(need {MIN_SHARED_TYPES})")
        if not gene_gate:
            print(f"  BLOCKED: Only {n_orthologs:,} orthologs "
                  f"(need {MIN_SHARED_GENES:,})")

    # ── Write feasibility report ──────────────────────────────────────
    # Healthy cell counts for passing types
    healthy_obs = obs_df[obs_df["disease"] == "normal"] if "disease" in obs_df.columns else obs_df
    healthy_counts = healthy_obs["cell_type"].value_counts()

    report_lines = [
        "# Mouse Lemur Feasibility Check\n",
        f"**Date:** 2026-04-05",
        f"**Status:** {'GO' if go else 'NO-GO'}",
        f"**Organism:** Microcebus murinus (gray mouse lemur)",
        f"**Divergence from human:** ~75 Mya",
        f"**Data source:** CELLxGENE Census build {CENSUS_VERSION}",
        "",
        "---",
        "",
        "## 1. Data Availability",
        "",
        f"- **Total cells in Census:** {len(obs_df):,}",
        f"- **Healthy cells:** {len(healthy_obs):,}",
        f"- **Unique cell types:** {len(all_types)}",
        f"- **Unique tissues:** {obs_df['tissue'].nunique()}",
        "",
        "### Tissue breakdown",
        "",
        "| Tissue | Cells |",
        "|--------|-------|",
    ]

    tissues = obs_df["tissue"].value_counts()
    for tissue, n in tissues.items():
        report_lines.append(f"| {tissue} | {n:,} |")

    report_lines.extend([
        "",
        "### Assay breakdown",
        "",
        "| Assay | Cells |",
        "|-------|-------|",
    ])

    assays = obs_df["assay"].value_counts()
    for assay, n in assays.items():
        report_lines.append(f"| {assay} | {n:,} |")

    report_lines.extend([
        "",
        "---",
        "",
        "## 2. Cell Type Overlap with Primary 35 Types",
        "",
        f"**Passing (≥{MIN_CELLS_PER_TYPE} cells):** {n_passing}",
        f"**Gate (≥{MIN_SHARED_TYPES}):** {'PASS' if type_gate else 'FAIL'}",
        "",
        "| # | Cell Type | Lemur Cells | Status | HM Rank (of 35) |",
        "|---|-----------|-------------|--------|------------------|",
    ])

    # Build rank lookup
    rank_lookup = {ct: i + 1 for i, ct in enumerate(PRIMARY_35_TYPES)}

    for i, (ct, n) in enumerate(sorted(passing, key=lambda x: x[1], reverse=True), 1):
        hm_rank = rank_lookup.get(ct, "?")
        report_lines.append(f"| {i} | {ct} | {n:,} | PASS | {hm_rank} |")

    if low:
        report_lines.extend([
            "",
            "### Types found but below 500-cell gate",
            "",
            "| Cell Type | Lemur Cells |",
            "|-----------|-------------|",
        ])
        for ct, n in sorted(low, key=lambda x: x[1], reverse=True):
            if n < MIN_CELLS_PER_TYPE:
                report_lines.append(f"| {ct} | {n:,} |")

    if missing:
        report_lines.extend([
            "",
            f"### Types missing ({len(missing)})",
            "",
        ])
        for ct in missing:
            report_lines.append(f"- {ct}")

    report_lines.extend([
        "",
        "---",
        "",
        "## 3. Ortholog Depth",
        "",
    ])

    if orthologs is not None:
        hm_orthologs = pd.read_csv(DATA_DIR / "orthologs_human_mouse.csv")
        hm_ids = set(hm_orthologs["human_ensembl_id"].dropna())
        lemur_ids = set(orthologs["human_ensembl_id"].dropna())
        three_way = hm_ids & lemur_ids

        report_lines.extend([
            "| Metric | Count | Gate | Status |",
            "|--------|-------|------|--------|",
            f"| Human-mouse lemur 1:1 orthologs | {len(lemur_ids):,} | ≥{MIN_SHARED_GENES:,} | "
            f"{'PASS' if len(lemur_ids) >= MIN_SHARED_GENES else 'FAIL'} |",
            f"| Human-mouse 1:1 orthologs (baseline) | {len(hm_ids):,} | — | — |",
            f"| Three-way intersection (H-M-ML) | {len(three_way):,} | ≥{MIN_SHARED_GENES:,} | "
            f"{'PASS' if len(three_way) >= MIN_SHARED_GENES else 'FAIL'} |",
            f"| Fraction of H-M space retained | {100 * len(three_way) / len(hm_ids):.1f}% | — | — |",
        ])

        macaque_cache = PROJECT_ROOT / "data" / "macaque" / "biomart_macaque_human_orthologs.csv"
        if macaque_cache.exists():
            mac_orthologs = pd.read_csv(macaque_cache)
            mac_ids = set(mac_orthologs["human_ensembl_id"].dropna())
            four_way = hm_ids & lemur_ids & mac_ids
            report_lines.extend([
                f"| Human-macaque 1:1 orthologs | {len(mac_ids):,} | — | — |",
                f"| Four-way intersection (H-M-Mac-ML) | {len(four_way):,} | ≥{MIN_SHARED_GENES:,} | "
                f"{'PASS' if len(four_way) >= MIN_SHARED_GENES else 'FAIL'} |",
            ])
    else:
        report_lines.append("**BioMart query failed.** Cannot assess ortholog depth.")

    report_lines.extend([
        "",
        "---",
        "",
        f"## 4. GO / NO-GO Decision",
        "",
        f"**Decision: {'GO' if go else 'NO-GO'}**",
        "",
        f"- Cell type gate: {'PASS' if type_gate else 'FAIL'} "
        f"({n_passing} types ≥ {MIN_SHARED_TYPES} minimum)",
        f"- Ortholog gate: {'PASS' if gene_gate else 'FAIL'} "
        f"({n_orthologs:,} genes ≥ {MIN_SHARED_GENES:,} minimum)",
        "",
    ])

    if go:
        report_lines.extend([
            "### Proceed to Step 1: Data Preparation",
            "",
            "The mouse lemur data passes both feasibility gates. Next steps:",
            "1. Download expression data from Census",
            "2. Map genes to human orthologs via BioMart",
            "3. Match cell type annotations to primary analysis types",
            "4. Run Procrustes analysis (Step 2)",
        ])
    else:
        report_lines.extend([
            "### Recommendations",
            "",
        ])
        if not type_gate:
            report_lines.append(
                f"- Only {n_passing} cell types pass the 500-cell gate. "
                f"Consider lowering the gate or broadening cell type matching."
            )
        if not gene_gate:
            report_lines.append(
                f"- Ortholog depth ({n_orthologs:,}) is below the 12,000 threshold. "
                f"Mouse lemur genome annotation may be incomplete."
            )

    report = "\n".join(report_lines) + "\n"

    with open(OUT_DIR / "feasibility_check.md", "w") as f:
        f.write(report)
    print(f"\n  Feasibility report saved to {OUT_DIR / 'feasibility_check.md'}")

    # Save cell type counts for downstream use
    type_counts = []
    for ct in all_types.index:
        in_primary = ct in PRIMARY_35_TYPES
        healthy_n = healthy_counts.get(ct, 0)
        type_counts.append({
            "cell_type": ct,
            "total_cells": all_types[ct],
            "healthy_cells": healthy_n,
            "in_primary_35": in_primary,
            "passes_gate": healthy_n >= MIN_CELLS_PER_TYPE and in_primary,
        })
    counts_df = pd.DataFrame(type_counts)
    counts_df.to_csv(OUT_DIR / "lemur_cell_type_counts.csv", index=False)
    print(f"  Cell type counts saved to {OUT_DIR / 'lemur_cell_type_counts.csv'}")

    return go


if __name__ == "__main__":
    result = main()
    sys.exit(0 if result else 1)
