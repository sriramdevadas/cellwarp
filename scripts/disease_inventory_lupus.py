#!/usr/bin/env python3
"""
Lupus (SLE) disease replication feasibility inventory.

Queries CELLxGENE Census for lupus / SLE data without downloading any expression
matrices.  Reports cell counts per cell_type label in lupus vs matched normal,
maps them onto the 35 cross-species types from the scaled Procrustes analysis,
and assesses feasibility for a disease-axis replication of the evolutionary
rigidity hypothesis.

Also inventories rheumatoid arthritis and NAFLD/NASH as backups.

Biology: If evolutionary rigidity (low Procrustes residual) predicts resistance
to disease-driven deformation, we should see a negative Spearman correlation
between cross-species Procrustes residual and lupus-vs-normal Procrustes
residual across matched cell types.

Math: No computation here — pure counting and feasibility assessment.
"""

import cellxgene_census
import pandas as pd
import numpy as np
import os, sys
from collections import defaultdict

OUTPUT_DIR = "output/disease_replication"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── 1. Load the 35 cross-species cell type names ────────────────────────
print("=" * 70)
print("STEP 1: Loading 35 cross-species cell types")
print("=" * 70)

human_centroids = pd.read_csv(
    "output/phase2/scaled_35types/centroids_human_35.csv", index_col=0
)
CROSS_SPECIES_TYPES = sorted(human_centroids.index.tolist())
print(f"Loaded {len(CROSS_SPECIES_TYPES)} cross-species cell types")
for ct in CROSS_SPECIES_TYPES:
    print(f"  - {ct}")

# ── 2. Query Census for lupus / SLE ─────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 2: Querying CELLxGENE Census for lupus data")
print("=" * 70)

census = cellxgene_census.open_soma(census_version="2025-11-08")

COLS = [
    "cell_type", "disease", "tissue", "tissue_general",
    "dataset_id", "assay", "is_primary_data",
    "donor_id", "sex", "development_stage",
    "suspension_type",
]

print("\nQuerying for lupus-related diseases...")
obs_lupus = cellxgene_census.get_obs(
    census, "Homo sapiens",
    value_filter=(
        "is_primary_data == True and "
        "(disease == 'systemic lupus erythematosus' or "
        " disease == 'lupus nephritis' or "
        " disease == 'lupus erythematosus')"
    ),
    column_names=COLS,
)

print(f"\nTotal lupus cells found: {len(obs_lupus):,}")

if len(obs_lupus) == 0:
    print("WARNING: No exact matches found for lupus disease labels.")
    print("The Census may use different ontology terms.")
    census.close()
    sys.exit(1)

print("\n--- Disease labels ---")
for disease, count in obs_lupus["disease"].value_counts().items():
    print(f"  {disease}: {count:,} cells")

print("\n--- Tissues ---")
for tissue, count in obs_lupus["tissue"].value_counts().items():
    print(f"  {tissue}: {count:,} cells")

print("\n--- Tissue (general) ---")
for tg, count in obs_lupus["tissue_general"].value_counts().items():
    print(f"  {tg}: {count:,} cells")

print("\n--- Assays ---")
for assay, count in obs_lupus["assay"].value_counts().items():
    print(f"  {assay}: {count:,} cells")

print("\n--- Datasets ---")
for did, count in obs_lupus["dataset_id"].value_counts().items():
    print(f"  {did}: {count:,} cells")

n_donors = obs_lupus["donor_id"].nunique()
n_datasets = obs_lupus["dataset_id"].nunique()
n_collections = 0  # collection_id removed from Census 2025-11-08
print(f"\nUnique donors: {n_donors}")
print(f"Unique datasets: {n_datasets}")

# ── 3. Cell counts: lupus vs matched normal ──────────────────────────────
print("\n" + "=" * 70)
print("STEP 3: Cell type counts — lupus vs matched normal")
print("=" * 70)

best_lupus_label = obs_lupus["disease"].value_counts().index[0]
print(f"\nBest lupus label: '{best_lupus_label}'")

lupus_df = obs_lupus[obs_lupus["disease"] == best_lupus_label].copy()
lupus_tissues = lupus_df["tissue"].unique().tolist()
print(f"Lupus tissues ({len(lupus_tissues)}): {lupus_tissues[:15]}")
if len(lupus_tissues) > 15:
    print(f"  ... and {len(lupus_tissues) - 15} more")

lupus_ct_counts = lupus_df["cell_type"].value_counts()
print(f"\n--- All lupus cell types ({len(lupus_ct_counts)} unique) ---")
for ct, count in lupus_ct_counts.items():
    marker = " ** MATCH" if ct in CROSS_SPECIES_TYPES else ""
    print(f"  {ct}: {count:,}{marker}")

# Query matched normal
print("\n--- Querying matched normal controls (same tissues) ---")
tissue_parts = [f"tissue == '{t}'" for t in lupus_tissues]
tissue_filter = " or ".join(tissue_parts)

obs_normal = cellxgene_census.get_obs(
    census, "Homo sapiens",
    value_filter=(
        f"is_primary_data == True and disease == 'normal' and ({tissue_filter})"
    ),
    column_names=COLS,
)

print(f"Total normal cells in matched tissues: {len(obs_normal):,}")
print(f"Normal donors: {obs_normal['donor_id'].nunique()}")
print(f"Normal datasets: {obs_normal['dataset_id'].nunique()}")

normal_ct_counts = obs_normal["cell_type"].value_counts()

# ── 4. Map to 35-type set ────────────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 4: Coarse mapping to 35 cross-species types")
print("=" * 70)
print("Rules: separate CD4/CD8, separate macrophage/monocyte,")
print("NO catch-all categories, only coherent biological identity.\n")

# Biologically justified coarsening map
# Census fine label → our 35-type label (or None = exclude)
COARSE_MAP = {
    # Direct matches
    "B cell": "B cell",
    "CD4-positive, alpha-beta T cell": "CD4-positive, alpha-beta T cell",
    "CD8-positive, alpha-beta T cell": "CD8-positive, alpha-beta T cell",
    "classical monocyte": "classical monocyte",
    "non-classical monocyte": "non-classical monocyte",
    "intermediate monocyte": "intermediate monocyte",
    "monocyte": "monocyte",
    "macrophage": "macrophage",
    "natural killer cell": "natural killer cell",
    "mature NK T cell": "mature NK T cell",
    "plasma cell": "plasma cell",
    "endothelial cell": "endothelial cell",
    "fibroblast": "fibroblast",
    "epithelial cell": "epithelial cell",
    "myeloid dendritic cell": "myeloid dendritic cell",
    "neutrophil": "neutrophil",
    "granulocyte": "granulocyte",
    "hematopoietic stem cell": "hematopoietic stem cell",
    "hematopoietic precursor cell": "hematopoietic precursor cell",
    "myeloid leukocyte": "myeloid leukocyte",
    "hepatocyte": "hepatocyte",
    "smooth muscle cell": "smooth muscle cell",
    "stromal cell": "stromal cell",
    "mesenchymal stem cell": "mesenchymal stem cell",
    "basal cell": "basal cell",
    # Coarsening — biologically justified
    "memory B cell": "B cell",
    "naive B cell": "B cell",
    "plasmablast": "plasma cell",
    "CD14-positive, CD16-negative classical monocyte": "classical monocyte",
    "CD16-positive, CD56-dim natural killer cell": "natural killer cell",
    "effector memory CD4-positive, alpha-beta T cell": "CD4-positive, alpha-beta T cell",
    "central memory CD4-positive, alpha-beta T cell": "CD4-positive, alpha-beta T cell",
    "naive thymus-derived CD4-positive, alpha-beta T cell": "CD4-positive, alpha-beta T cell",
    "regulatory T cell": "CD4-positive, alpha-beta T cell",
    "effector memory CD8-positive, alpha-beta T cell": "CD8-positive, alpha-beta T cell",
    "central memory CD8-positive, alpha-beta T cell": "CD8-positive, alpha-beta T cell",
    "naive thymus-derived CD8-positive, alpha-beta T cell": "CD8-positive, alpha-beta T cell",
    "conventional dendritic cell": "myeloid dendritic cell",
    "dendritic cell": "myeloid dendritic cell",
    "kidney epithelial cell": "epithelial cell",
    "kidney collecting duct epithelial cell": "epithelial cell",
    "proximal tubular epithelial cell": "epithelial cell",
    "renal intercalated cell": "epithelial cell",
    "renal principal cell": "epithelial cell",
    "renal alpha-intercalated cell": "epithelial cell",
    "kidney distal convoluted tubule epithelial cell": "epithelial cell",
    "kidney loop of Henle thick ascending limb epithelial cell": "epithelial cell",
    "glomerular endothelial cell": "endothelial cell",
    "vascular associated smooth muscle cell": "smooth muscle cell",
    # EXCLUDED — catch-all or lineage-ambiguous
    "T cell": None,
    "immune cell": None,
    "leukocyte": None,
    "lymphocyte": None,
    "mucosal invariant T cell": None,
    "gamma-delta T cell": None,
    "plasmacytoid dendritic cell": None,
    "erythrocyte": None,
    "platelet": None,
    "megakaryocyte": None,
    "podocyte": None,
    "mesangial cell": None,
    "pericyte": None,
}

coarse_lupus = defaultdict(int)
coarse_normal = defaultdict(int)
unmapped_lupus = defaultdict(int)
unmapped_normal = defaultdict(int)

for ct, count in lupus_ct_counts.items():
    if ct in COARSE_MAP:
        target = COARSE_MAP[ct]
        if target is not None:
            coarse_lupus[target] += count
    elif ct in CROSS_SPECIES_TYPES:
        coarse_lupus[ct] += count
    else:
        unmapped_lupus[ct] += count

for ct, count in normal_ct_counts.items():
    if ct in COARSE_MAP:
        target = COARSE_MAP[ct]
        if target is not None:
            coarse_normal[target] += count
    elif ct in CROSS_SPECIES_TYPES:
        coarse_normal[ct] += count
    else:
        unmapped_normal[ct] += count

# Build coarse inventory (only types with >0 in either condition)
coarse_rows = []
for ct in CROSS_SPECIES_TYPES:
    l_n = coarse_lupus.get(ct, 0)
    n_n = coarse_normal.get(ct, 0)
    if l_n > 0 or n_n > 0:
        coarse_rows.append({
            "cross_species_type": ct,
            "lupus_cells": l_n,
            "normal_cells": n_n,
            "passes_500_gate": l_n >= 500 and n_n >= 500,
        })

coarse_df = pd.DataFrame(coarse_rows).sort_values("lupus_cells", ascending=False)

print(f"{'Cell Type':<50} {'Lupus':>8} {'Normal':>8} {'Pass?':>6}")
print("-" * 75)
for _, row in coarse_df.iterrows():
    flag = "  YES" if row["passes_500_gate"] else "   no"
    print(f"{row['cross_species_type']:<50} {row['lupus_cells']:>8,} "
          f"{row['normal_cells']:>8,} {flag}")

n_pass = int(coarse_df["passes_500_gate"].sum())
print(f"\nTypes passing >=500 gate in BOTH conditions: {n_pass}")
print(f"Maximum feasible n for Spearman: {n_pass}")

# Unmapped
if unmapped_lupus:
    print(f"\n--- Unmapped lupus cell types (>50 cells) ---")
    for ct, count in sorted(unmapped_lupus.items(), key=lambda x: -x[1]):
        if count > 50:
            print(f"  {ct}: {count:,}")

# ── 5. Dataset provenance ────────────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 5: Dataset provenance and study independence")
print("=" * 70)

print(f"\n--- Lupus: {n_datasets} datasets ---")
for did in lupus_df["dataset_id"].unique():
    sub = lupus_df[lupus_df["dataset_id"] == did]
    print(f"\n  Dataset: {did}")
    print(f"    Cells: {len(sub):,}  |  Donors: {sub['donor_id'].nunique()}")
    print(f"    Tissues: {sub['tissue'].unique().tolist()[:5]}")
    print(f"    Assays: {sub['assay'].unique().tolist()}")
    # collection_id removed from Census 2025-11-08

print(f"\n--- Normal control adequacy ---")
lupus_assays = set(lupus_df["assay"].unique())
normal_assays = set(obs_normal["assay"].unique())
print(f"Lupus assays: {lupus_assays}")
print(f"Normal assays: {normal_assays}")
print(f"Shared assays: {lupus_assays & normal_assays}")
lupus_t = set(lupus_df["tissue"].unique())
normal_t = set(obs_normal["tissue"].unique())
print(f"Tissue overlap: {lupus_t & normal_t}")
print(f"Lupus-only tissues: {lupus_t - normal_t}")

# ── 6. Backup disease states ─────────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 6: Backup disease inventories")
print("=" * 70)

# --- Rheumatoid arthritis ---
print("\n--- Rheumatoid Arthritis ---")
obs_ra = cellxgene_census.get_obs(
    census, "Homo sapiens",
    value_filter="is_primary_data == True and disease == 'rheumatoid arthritis'",
    column_names=["cell_type", "disease", "tissue", "dataset_id",
                  "donor_id", "assay"],
)
print(f"Total RA cells: {len(obs_ra):,}")
print(f"Donors: {obs_ra['donor_id'].nunique()}")
print(f"Datasets: {obs_ra['dataset_id'].nunique()}")
print(f"Tissues: {obs_ra['tissue'].unique().tolist()}")
print(f"Assays: {obs_ra['assay'].unique().tolist()}")

ra_ct = obs_ra["cell_type"].value_counts()
print(f"\nRA cell types (top 30 of {len(ra_ct)}):")
for ct, count in ra_ct.head(30).items():
    marker = " **" if ct in CROSS_SPECIES_TYPES else ""
    print(f"  {ct}: {count:,}{marker}")

# Quick RA gate check
ra_pass = 0
for ct in CROSS_SPECIES_TYPES:
    if ra_ct.get(ct, 0) >= 500:
        ra_pass += 1
print(f"\nRA direct matches passing >=500: {ra_pass}")

# Also get RA normal
ra_tissues = obs_ra["tissue"].unique().tolist()
ra_tissue_filter = " or ".join([f"tissue == '{t}'" for t in ra_tissues])
obs_ra_normal = cellxgene_census.get_obs(
    census, "Homo sapiens",
    value_filter=f"is_primary_data == True and disease == 'normal' and ({ra_tissue_filter})",
    column_names=["cell_type", "tissue", "dataset_id", "donor_id"],
)
print(f"RA matched normal cells: {len(obs_ra_normal):,}")
ra_normal_ct = obs_ra_normal["cell_type"].value_counts()

ra_both_pass = 0
for ct in CROSS_SPECIES_TYPES:
    if ra_ct.get(ct, 0) >= 500 and ra_normal_ct.get(ct, 0) >= 500:
        ra_both_pass += 1
print(f"RA types passing >=500 in BOTH disease+normal: {ra_both_pass}")

# --- NAFLD/NASH ---
print("\n--- NAFLD / NASH ---")
obs_nafld = cellxgene_census.get_obs(
    census, "Homo sapiens",
    value_filter=(
        "is_primary_data == True and "
        "(disease == 'non-alcoholic fatty liver disease' or "
        " disease == 'non-alcoholic steatohepatitis')"
    ),
    column_names=["cell_type", "disease", "tissue", "dataset_id",
                  "donor_id", "assay"],
)
print(f"Total NAFLD/NASH cells: {len(obs_nafld):,}")
if len(obs_nafld) > 0:
    print(f"Disease labels: {obs_nafld['disease'].value_counts().to_dict()}")
    print(f"Donors: {obs_nafld['donor_id'].nunique()}")
    print(f"Datasets: {obs_nafld['dataset_id'].nunique()}")
    print(f"Tissues: {obs_nafld['tissue'].unique().tolist()}")
    nafld_ct = obs_nafld["cell_type"].value_counts()
    print(f"\nNAFLD/NASH cell types (top 30 of {len(nafld_ct)}):")
    for ct, count in nafld_ct.head(30).items():
        marker = " **" if ct in CROSS_SPECIES_TYPES else ""
        print(f"  {ct}: {count:,}{marker}")

    # Gate check
    nafld_tissues = obs_nafld["tissue"].unique().tolist()
    nafld_tf = " or ".join([f"tissue == '{t}'" for t in nafld_tissues])
    obs_nafld_normal = cellxgene_census.get_obs(
        census, "Homo sapiens",
        value_filter=f"is_primary_data == True and disease == 'normal' and ({nafld_tf})",
        column_names=["cell_type", "tissue", "donor_id"],
    )
    print(f"NAFLD matched normal cells: {len(obs_nafld_normal):,}")
    nafld_normal_ct = obs_nafld_normal["cell_type"].value_counts()
    nafld_both = 0
    for ct in CROSS_SPECIES_TYPES:
        if nafld_ct.get(ct, 0) >= 500 and nafld_normal_ct.get(ct, 0) >= 500:
            nafld_both += 1
    print(f"NAFLD types passing >=500 in BOTH: {nafld_both}")
else:
    print("No NAFLD/NASH data found. Trying alternate labels...")
    # Try liver disease variants
    for label in ["fatty liver disease", "steatohepatitis",
                   "liver disease", "metabolic dysfunction-associated steatotic liver disease"]:
        try:
            obs_try = cellxgene_census.get_obs(
                census, "Homo sapiens",
                value_filter=f"is_primary_data == True and disease == '{label}'",
                column_names=["cell_type"],
            )
            if len(obs_try) > 0:
                print(f"  '{label}': {len(obs_try):,} cells")
        except Exception as e:
            print(f"  '{label}': query error — {e}")

census.close()

# ── 7. Save outputs ──────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 7: Saving inventory files")
print("=" * 70)

# Save CSV
coarse_df.to_csv(f"{OUTPUT_DIR}/lupus_inventory.csv", index=False)
print(f"Saved: {OUTPUT_DIR}/lupus_inventory.csv")

# Save comprehensive text report
with open(f"{OUTPUT_DIR}/lupus_inventory.txt", "w") as f:
    f.write("LUPUS (SLE) DISEASE REPLICATION — FEASIBILITY INVENTORY\n")
    f.write(f"Generated: 2026-03-14\n")
    f.write("=" * 70 + "\n\n")

    f.write("SUMMARY\n")
    f.write("-" * 40 + "\n")
    f.write(f"Total lupus cells in Census: {len(obs_lupus):,}\n")
    f.write(f"Best disease label: '{best_lupus_label}'\n")
    f.write(f"Lupus donors: {n_donors}\n")
    f.write(f"Lupus datasets: {n_datasets}\n")
    # collection_id no longer available in Census 2025-11-08
    f.write(f"Matched normal cells: {len(obs_normal):,}\n\n")

    f.write(f"Coherent types passing >=500 in BOTH conditions: {n_pass}\n")
    f.write(f"Maximum feasible n for Spearman: {n_pass}\n\n")

    f.write("COARSE-MAPPED INVENTORY\n")
    f.write("-" * 40 + "\n")
    f.write(f"{'Cell Type':<50} {'Lupus':>8} {'Normal':>8} {'Pass':>6}\n")
    f.write("-" * 75 + "\n")
    for _, row in coarse_df.iterrows():
        flag = "YES" if row["passes_500_gate"] else "no"
        f.write(f"{row['cross_species_type']:<50} {row['lupus_cells']:>8,} "
                f"{row['normal_cells']:>8,} {flag:>6}\n")

    f.write(f"\nTypes passing >=500 gate in BOTH: {n_pass}\n")
    f.write(f"Maximum feasible n for Spearman: {n_pass}\n\n")

    # Unmapped lupus types
    if unmapped_lupus:
        f.write("UNMAPPED LUPUS CELL TYPES (>50 cells)\n")
        f.write("-" * 40 + "\n")
        for ct, count in sorted(unmapped_lupus.items(), key=lambda x: -x[1]):
            if count > 50:
                f.write(f"  {ct}: {count:,}\n")
        f.write("\n")

    # Dataset provenance
    f.write("DATASET PROVENANCE\n")
    f.write("-" * 40 + "\n")
    for did in lupus_df["dataset_id"].unique():
        sub = lupus_df[lupus_df["dataset_id"] == did]
        f.write(f"  Dataset: {did}\n")
        f.write(f"    Cells: {len(sub):,}  |  Donors: {sub['donor_id'].nunique()}\n")
        f.write(f"    Tissues: {sub['tissue'].unique().tolist()[:5]}\n")
        f.write(f"    Assays: {sub['assay'].unique().tolist()}\n")
        f.write(f"\n")

    # Normal control adequacy
    f.write("NORMAL CONTROL ADEQUACY\n")
    f.write("-" * 40 + "\n")
    f.write(f"Lupus assays: {lupus_assays}\n")
    f.write(f"Normal assays: {normal_assays}\n")
    f.write(f"Shared assays: {lupus_assays & normal_assays}\n")
    f.write(f"Tissue overlap: {lupus_t & normal_t}\n")
    f.write(f"Lupus-only tissues: {lupus_t - normal_t}\n\n")

    # Backup diseases
    f.write("BACKUP DISEASE STATES\n")
    f.write("-" * 40 + "\n")
    f.write(f"Rheumatoid Arthritis: {len(obs_ra):,} cells, "
            f"{obs_ra['donor_id'].nunique()} donors, "
            f"{obs_ra['dataset_id'].nunique()} datasets\n")
    f.write(f"  Types passing >=500 in BOTH disease+normal: {ra_both_pass}\n")
    f.write(f"NAFLD/NASH: {len(obs_nafld):,} cells")
    if len(obs_nafld) > 0:
        f.write(f", {obs_nafld['donor_id'].nunique()} donors, "
                f"{obs_nafld['dataset_id'].nunique()} datasets\n")
        f.write(f"  Types passing >=500 in BOTH disease+normal: {nafld_both}\n")
    else:
        f.write(" (no data found)\n")

print(f"Saved: {OUTPUT_DIR}/lupus_inventory.txt")

print("\n" + "=" * 70)
print("INVENTORY COMPLETE")
print("=" * 70)