"""
NHP Ortholog Feasibility Assessment — Human-Macaque 1:1 Orthologs

Biology
-------
To extend the CellWarp cross-species Procrustes framework to non-human primates
(NHP), we need to establish 1:1 ortholog mapping between human and macaque. Since
human and macaque diverged ~25 million years ago (vs ~80M for human-mouse), we
expect a LARGER shared ortholog space (more genes conserved at shorter evolutionary
distance).

We query Ensembl BioMart for human-macaque 1:1 orthologs, then assess overlap with
the existing 16,959-gene human-mouse ortholog space used in the core CellWarp pipeline.
The >=12,000 shared gene gate (Phase 1 criterion) must pass for the three-way
human-mouse-macaque intersection.

Math
----
Three-way ortholog intersection: genes that have 1:1 orthologs in BOTH mouse AND
macaque relative to human. This is the gene space needed for a three-species
Procrustes analysis. For pairwise human-macaque analysis, only the human-macaque
1:1 set is needed.

Species targets:
- Primary: Macaca fascicularis (crab-eating macaque, cynomolgus)
- Fallback: Macaca mulatta (rhesus macaque, >95% sequence identity with fascicularis)
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
from pybiomart import Dataset

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data" / "phase1"
OUTPUT_DIR = PROJECT_ROOT / "output" / "validation" / "nhp_feasibility"

HUMAN_MOUSE_ORTHOLOGS = DATA_DIR / "orthologs_human_mouse.csv"
OUTPUT_FILE = OUTPUT_DIR / "nhp_ortholog_assessment.csv"

# ---------------------------------------------------------------------------
# Step 1: Load existing human-mouse orthologs
# ---------------------------------------------------------------------------
print("=" * 70)
print("NHP ORTHOLOG FEASIBILITY ASSESSMENT")
print("=" * 70)

print("\n--- Step 1: Load existing human-mouse 1:1 orthologs ---")
hm_orthologs = pd.read_csv(HUMAN_MOUSE_ORTHOLOGS)
print(f"  Human-mouse 1:1 orthologs loaded: {len(hm_orthologs):,} genes")
print(f"  Columns: {list(hm_orthologs.columns)}")

# Extract the set of human Ensembl IDs in our existing space
human_ensembl_ids_in_space = set(hm_orthologs["human_ensembl_id"].dropna())
human_gene_names_in_space = set(hm_orthologs["human_gene_name"].dropna())
print(f"  Unique human Ensembl IDs: {len(human_ensembl_ids_in_space):,}")
print(f"  Unique human gene names: {len(human_gene_names_in_space):,}")

# ---------------------------------------------------------------------------
# Step 2: Query Ensembl BioMart for human-macaque 1:1 orthologs
# ---------------------------------------------------------------------------
print("\n--- Step 2: Query Ensembl BioMart for human-macaque orthologs ---")

# Try Macaca fascicularis first, fall back to Macaca mulatta
MACAQUE_SPECIES = [
    ("mfascicularis", "Macaca fascicularis (crab-eating macaque)"),
    ("mmulatta", "Macaca mulatta (rhesus macaque)"),
]

macaque_orthologs = None
species_used = None

for species_code, species_name in MACAQUE_SPECIES:
    print(f"\n  Trying {species_name} ({species_code})...")

    # Attributes for querying macaque orthologs from the human gene dataset
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
        print(f"  Columns returned: {list(result.columns)}")

        # Filter for 1:1 orthologs
        # pybiomart returns display names, not attribute names
        # e.g. "Crab-eating macaque homology type" or "Macaque homology type"
        orthology_col = None
        for col in result.columns:
            if "homology type" in col.lower() or "orthology_type" in col.lower():
                orthology_col = col
                break
        if orthology_col is None:
            print(f"  WARNING: Could not find orthology type column in: {list(result.columns)}")
            continue

        print(f"  Using orthology type column: {orthology_col}")
        print(f"  Orthology type distribution:")
        print(result[orthology_col].value_counts().to_string())

        one2one = result[result[orthology_col] == "ortholog_one2one"].copy()
        print(f"\n  1:1 orthologs: {len(one2one):,}")

        # Drop rows with missing gene IDs
        # pybiomart display columns: "Gene stable ID", "Gene name",
        # "Crab-eating macaque gene stable ID", "Crab-eating macaque gene name",
        # "Crab-eating macaque homology type"  (or "Macaque ..." for rhesus)
        gene_id_col = "Gene stable ID"
        macaque_gene_col = None
        macaque_name_col = None

        for col in one2one.columns:
            col_lower = col.lower()
            # Skip human columns (no species prefix)
            if col_lower in ("gene stable id", "gene name"):
                continue
            if "gene stable id" in col_lower or ("homolog" in col_lower and "ensembl" in col_lower):
                macaque_gene_col = col
            elif "gene name" in col_lower and "homology" not in col_lower:
                macaque_name_col = col

        if gene_id_col not in one2one.columns:
            # Fallback — first column is usually the human Ensembl ID
            gene_id_col = one2one.columns[0]

        print(f"  Human gene ID column: {gene_id_col}")
        print(f"  Macaque gene ID column: {macaque_gene_col}")
        print(f"  Macaque gene name column: {macaque_name_col}")

        # Drop rows with missing human or macaque gene IDs
        before = len(one2one)
        one2one = one2one.dropna(subset=[gene_id_col])
        if macaque_gene_col:
            one2one = one2one.dropna(subset=[macaque_gene_col])
            one2one = one2one[one2one[macaque_gene_col] != ""]
        after = len(one2one)
        print(f"  After dropping missing IDs: {after:,} (dropped {before - after:,})")

        # Deduplicate by human gene ID (keep first)
        before = len(one2one)
        one2one = one2one.drop_duplicates(subset=[gene_id_col])
        after = len(one2one)
        print(f"  After deduplication by human Ensembl ID: {after:,} (dropped {before - after:,})")

        # Store the detected column names alongside the data
        macaque_orthologs = one2one
        species_used = (species_code, species_name)
        detected_cols = {
            "gene_id": gene_id_col,
            "macaque_gene_id": macaque_gene_col,
            "macaque_gene_name": macaque_name_col,
        }
        # Detect human gene name column
        for col in one2one.columns:
            col_lower = col.lower()
            if col_lower == "gene name":
                detected_cols["gene_name"] = col
                break
        print(f"\n  SUCCESS: Using {species_name}")
        break

    except Exception as e:
        print(f"  FAILED for {species_name}: {e}")
        continue

if macaque_orthologs is None:
    print("\nERROR: Could not retrieve macaque orthologs from either species.")
    sys.exit(1)

# ---------------------------------------------------------------------------
# Step 3: Compute overlaps
# ---------------------------------------------------------------------------
print("\n--- Step 3: Compute overlap with existing human-mouse gene space ---")

species_code, species_name = species_used

# Use column names detected during Step 2
gene_id_col = detected_cols["gene_id"]
gene_name_col = detected_cols.get("gene_name", macaque_orthologs.columns[1])
macaque_gene_col = detected_cols["macaque_gene_id"]
macaque_name_col = detected_cols["macaque_gene_name"]

print(f"  Using columns: human_id={gene_id_col}, human_name={gene_name_col}")
print(f"                 macaque_id={macaque_gene_col}, macaque_name={macaque_name_col}")

human_macaque_ensembl_ids = set(macaque_orthologs[gene_id_col].dropna())
print(f"  Human-macaque 1:1 orthologs (human Ensembl IDs): {len(human_macaque_ensembl_ids):,}")
print(f"  Human-mouse 1:1 orthologs (human Ensembl IDs): {len(human_ensembl_ids_in_space):,}")

# Pairwise overlap: human-macaque vs human-mouse
overlap = human_macaque_ensembl_ids & human_ensembl_ids_in_space
print(f"\n  Three-way intersection (human-mouse-macaque): {len(overlap):,} genes")

# Only in human-macaque (not in human-mouse)
only_macaque = human_macaque_ensembl_ids - human_ensembl_ids_in_space
print(f"  Only in human-macaque (not in human-mouse): {len(only_macaque):,} genes")

# Only in human-mouse (not in human-macaque)
only_mouse = human_ensembl_ids_in_space - human_macaque_ensembl_ids
print(f"  Only in human-mouse (not in human-macaque): {len(only_mouse):,} genes")

# ---------------------------------------------------------------------------
# Step 4: Gate check
# ---------------------------------------------------------------------------
print("\n--- Step 4: Phase 1 gate check (>=12,000 shared genes) ---")

GATE = 12_000

pairwise_pass = len(human_macaque_ensembl_ids) >= GATE
three_way_pass = len(overlap) >= GATE

print(f"  Pairwise human-macaque: {len(human_macaque_ensembl_ids):,} genes → "
      f"{'PASS' if pairwise_pass else 'FAIL'} (gate: >= {GATE:,})")
print(f"  Three-way human-mouse-macaque: {len(overlap):,} genes → "
      f"{'PASS' if three_way_pass else 'FAIL'} (gate: >= {GATE:,})")

# ---------------------------------------------------------------------------
# Step 5: Save results
# ---------------------------------------------------------------------------
print("\n--- Step 5: Save output ---")

# Build output dataframe
records = []
for _, row in macaque_orthologs.iterrows():
    human_id = row[gene_id_col]
    human_name = row.get(gene_name_col, "")
    macaque_id = row.get(macaque_gene_col, "")
    macaque_name = row.get(macaque_name_col, "")
    in_space = human_id in human_ensembl_ids_in_space
    records.append({
        "human_ensembl_id": human_id,
        "human_gene": human_name if pd.notna(human_name) else "",
        "macaque_ensembl_id": macaque_id if pd.notna(macaque_id) else "",
        "macaque_ortholog": macaque_name if pd.notna(macaque_name) else "",
        "in_our_space": in_space,
    })

output_df = pd.DataFrame(records)
output_df.to_csv(OUTPUT_FILE, index=False)
print(f"  Saved {len(output_df):,} rows to {OUTPUT_FILE}")

# Also save a summary text file
summary_file = OUTPUT_DIR / "nhp_ortholog_summary.txt"
with open(summary_file, "w") as f:
    f.write("NHP ORTHOLOG FEASIBILITY ASSESSMENT — SUMMARY\n")
    f.write("=" * 50 + "\n")
    f.write(f"Date: 2026-03-15\n")
    f.write(f"Macaque species used: {species_name}\n")
    if species_code == "mmulatta":
        f.write("NOTE: Macaca mulatta (rhesus) used as proxy for Macaca fascicularis\n")
        f.write("      (crab-eating macaque). >95% sequence identity between species.\n")
    f.write(f"\nHuman-macaque 1:1 orthologs: {len(human_macaque_ensembl_ids):,}\n")
    f.write(f"Human-mouse 1:1 orthologs (existing): {len(human_ensembl_ids_in_space):,}\n")
    f.write(f"Three-way intersection: {len(overlap):,}\n")
    f.write(f"Only in human-macaque: {len(only_macaque):,}\n")
    f.write(f"Only in human-mouse: {len(only_mouse):,}\n")
    f.write(f"\nGate check (>= {GATE:,} genes):\n")
    f.write(f"  Pairwise human-macaque: {'PASS' if pairwise_pass else 'FAIL'}\n")
    f.write(f"  Three-way intersection: {'PASS' if three_way_pass else 'FAIL'}\n")
print(f"  Saved summary to {summary_file}")

# ---------------------------------------------------------------------------
# Final Summary
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("RESULTS SUMMARY")
print("=" * 70)
print(f"  Macaque species: {species_name}")
if species_code == "mmulatta":
    print("  NOTE: Rhesus macaque used as proxy for crab-eating macaque")
    print("        (>95% sequence identity between species)")
print(f"  Human-macaque 1:1 orthologs: {len(human_macaque_ensembl_ids):,}")
print(f"  Existing human-mouse space:  {len(human_ensembl_ids_in_space):,}")
print(f"  Three-way intersection:      {len(overlap):,}")
print(f"  Overlap percentage:          {100 * len(overlap) / len(human_ensembl_ids_in_space):.1f}% of human-mouse space")
print(f"  Pairwise gate (>={GATE:,}):   {'PASS' if pairwise_pass else 'FAIL'}")
print(f"  Three-way gate (>={GATE:,}):  {'PASS' if three_way_pass else 'FAIL'}")
print("=" * 70)
