#!/usr/bin/env python3
"""
COVID-19 detailed feasibility inventory for disease-axis replication.

Queries CELLxGENE Census for COVID-19 / SARS-CoV-2 data and matched normal
controls.  Reports cell counts per cell type, tissue, and dataset without
downloading any expression matrices.

Biology: Tests whether evolutionary rigidity (low cross-species Procrustes
residual) predicts resistance to disease-driven geometric deformation.
COVID-19 is the primary multi-disease replication target following lupus
rejection (n=8, single dataset).

Math: No computation — pure counting and feasibility assessment.

Key requirements:
  1. COVID-19 data: disease contains "COVID" or "SARS-CoV-2" or "COVID-19"
  2. Matched normal: same tissues as COVID data
  3. Per cell type: COVID cells, normal cells, datasets, PASS/MARGINAL/FAIL
  4. Tissue-match verification for each passing cell type
  5. Maximum coherent n (exclude catch-alls, require ≥2 datasets, tissue-matched)
  6. Top 10 COVID datasets by cell count (study dominance check)
  7. Severity metadata availability
"""

import cellxgene_census
import pandas as pd
import numpy as np
import os
import sys
from collections import defaultdict

OUTPUT_DIR = "output/disease_replication"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── 1. Load the 35 cross-species cell type names + residuals ──────────
print("=" * 70)
print("STEP 1: Loading 35 cross-species cell types and Procrustes residuals")
print("=" * 70)

human_centroids = pd.read_csv(
    "output/phase2/scaled_35types/centroids_human_35.csv", index_col=0
)
CROSS_SPECIES_TYPES = sorted(human_centroids.index.tolist())

residuals_df = pd.read_csv("output/phase2/scaled_35types/residuals_ranked.csv")
RESIDUAL_MAP = dict(zip(residuals_df["cell_type"], residuals_df["residual_magnitude"]))
RANK_MAP = dict(zip(residuals_df["cell_type"], residuals_df["rank"]))

print(f"Loaded {len(CROSS_SPECIES_TYPES)} cross-species cell types")

# ── 2. Query Census for COVID-19 data ─────────────────────────────────
print("\n" + "=" * 70)
print("STEP 2: Querying CELLxGENE Census for COVID-19 data")
print("=" * 70)

census = cellxgene_census.open_soma(census_version="2025-11-08")

COLS = [
    "cell_type", "disease", "tissue", "tissue_general",
    "dataset_id", "assay", "is_primary_data",
    "donor_id", "sex", "development_stage",
    "suspension_type",
]

# COVID disease labels — use OR to catch all variants
print("\nQuerying for COVID-19-related diseases...")
obs_covid = cellxgene_census.get_obs(
    census, "Homo sapiens",
    value_filter=(
        "is_primary_data == True and "
        "(disease == 'COVID-19' or "
        " disease == 'SARS-CoV-2')"
    ),
    column_names=COLS,
)

print(f"\nTotal COVID cells found: {len(obs_covid):,}")

if len(obs_covid) == 0:
    print("WARNING: No exact matches. Trying broader search...")
    # Try with contains-like approach by querying individual terms
    for label in ["COVID-19", "SARS-CoV-2", "coronavirus infection",
                  "severe acute respiratory syndrome",
                  "COVID-19, mild", "COVID-19, severe",
                  "SARS-CoV-2 infection"]:
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
    sys.exit(1)

# ── 2a. Summarize COVID data ──────────────────────────────────────────
print("\n--- Disease labels ---")
for disease, count in obs_covid["disease"].value_counts().items():
    if count > 0:
        print(f"  {disease}: {count:,} cells")

print("\n--- Tissues ---")
covid_tissue_counts = obs_covid["tissue"].value_counts()
for tissue, count in covid_tissue_counts.items():
    if count > 0:
        print(f"  {tissue}: {count:,} cells")

print("\n--- Tissue (general) ---")
for tg, count in obs_covid["tissue_general"].value_counts().items():
    if count > 0:
        print(f"  {tg}: {count:,} cells")

print("\n--- Assays ---")
for assay, count in obs_covid["assay"].value_counts().items():
    if count > 0:
        print(f"  {assay}: {count:,} cells")

print("\n--- Development stages ---")
for stage, count in obs_covid["development_stage"].value_counts().items():
    if count > 0:
        print(f"  {stage}: {count:,} cells")

n_covid_donors = obs_covid["donor_id"].nunique()
n_covid_datasets = obs_covid["dataset_id"].nunique()
print(f"\nUnique COVID donors: {n_covid_donors}")
print(f"Unique COVID datasets: {n_covid_datasets}")

# ── 3. Query matched normal controls ─────────────────────────────────
print("\n" + "=" * 70)
print("STEP 3: Matched normal controls (same tissues as COVID)")
print("=" * 70)

all_covid_tissues = [t for t, c in obs_covid["tissue"].value_counts().items() if c > 0]
print(f"\nCOVID tissues ({len(all_covid_tissues)}): {all_covid_tissues}")

tissue_parts = [f"tissue == '{t}'" for t in all_covid_tissues]
tissue_filter = " or ".join(tissue_parts)

obs_normal = cellxgene_census.get_obs(
    census, "Homo sapiens",
    value_filter=(
        f"is_primary_data == True and disease == 'normal' and ({tissue_filter})"
    ),
    column_names=COLS,
)

print(f"\nTotal normal cells in matched tissues: {len(obs_normal):,}")
print(f"Normal donors: {obs_normal['donor_id'].nunique()}")
print(f"Normal datasets: {obs_normal['dataset_id'].nunique()}")

print("\n--- Normal tissue breakdown ---")
for tissue, count in obs_normal["tissue"].value_counts().items():
    if count > 0:
        print(f"  {tissue}: {count:,} cells")

# ── 4. Coarse mapping to 35-type set ─────────────────────────────────
print("\n" + "=" * 70)
print("STEP 4: Coarse mapping to 35 cross-species types")
print("=" * 70)
print("Rules: separate CD4/CD8, separate macrophage/monocyte,")
print("NO catch-all categories (T cell, lymphocyte, immune cell excluded),")
print("only coherent biological identity.\n")

# Biologically justified coarsening map
# Census fine label → our 35-type label (or None = exclude)
COARSE_MAP = {
    # Direct matches (35-type set)
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
    "T cell": "T cell",
    "pancreatic acinar cell": "pancreatic acinar cell",
    "pancreatic ductal cell": "pancreatic ductal cell",
    "luminal epithelial cell of mammary gland": "luminal epithelial cell of mammary gland",
    "large intestine goblet cell": "large intestine goblet cell",
    "enterocyte of epithelium of large intestine": "enterocyte of epithelium of large intestine",
    "fibroblast of cardiac tissue": "fibroblast of cardiac tissue",
    "mesenchymal stem cell of adipose tissue": "mesenchymal stem cell of adipose tissue",
    "adventitial cell": "adventitial cell",
    "bladder urothelial cell": "bladder urothelial cell",
    # Coarsening — biologically justified
    "memory B cell": "B cell",
    "naive B cell": "B cell",
    "plasmablast": "plasma cell",
    "CD14-positive, CD16-negative classical monocyte": "classical monocyte",
    "CD16-positive, CD56-dim natural killer cell": "natural killer cell",
    "CD56-bright natural killer cell": "natural killer cell",
    "effector memory CD4-positive, alpha-beta T cell": "CD4-positive, alpha-beta T cell",
    "central memory CD4-positive, alpha-beta T cell": "CD4-positive, alpha-beta T cell",
    "naive thymus-derived CD4-positive, alpha-beta T cell": "CD4-positive, alpha-beta T cell",
    "regulatory T cell": "CD4-positive, alpha-beta T cell",
    "effector memory CD8-positive, alpha-beta T cell": "CD8-positive, alpha-beta T cell",
    "central memory CD8-positive, alpha-beta T cell": "CD8-positive, alpha-beta T cell",
    "naive thymus-derived CD8-positive, alpha-beta T cell": "CD8-positive, alpha-beta T cell",
    "conventional dendritic cell": "myeloid dendritic cell",
    "dendritic cell": "myeloid dendritic cell",
    "alveolar macrophage": "macrophage",
    "lung macrophage": "macrophage",
    "Kupffer cell": "macrophage",
    "lung endothelial cell": "endothelial cell",
    "endothelial cell of artery": "endothelial cell",
    "endothelial cell of vascular tree": "endothelial cell",
    "capillary endothelial cell": "endothelial cell",
    "vein endothelial cell": "endothelial cell",
    "blood vessel endothelial cell": "endothelial cell",
    "glomerular endothelial cell": "endothelial cell",
    "endothelial cell of lymphatic vessel": "endothelial cell",
    "epithelial cell of lung": "epithelial cell",
    "type I pneumocyte": "epithelial cell",
    "type II pneumocyte": "epithelial cell",
    "club cell": "epithelial cell",
    "ciliated cell": "epithelial cell",
    "goblet cell": "epithelial cell",
    "respiratory goblet cell": "epithelial cell",
    "bronchial epithelial cell": "epithelial cell",
    "respiratory epithelial cell": "epithelial cell",
    "nasal epithelial cell": "epithelial cell",
    "mucus secreting cell": "epithelial cell",
    "secretory cell": "epithelial cell",
    "lung ciliated cell": "epithelial cell",
    "pulmonary ionocyte": "epithelial cell",
    "vascular associated smooth muscle cell": "smooth muscle cell",
    "bronchial smooth muscle cell": "smooth muscle cell",
    "pericyte": None,  # distinct from smooth muscle
    "myofibroblast cell": "fibroblast",
    "lung fibroblast": "fibroblast",
    "fibroblast of lung": "fibroblast",
    # EXCLUDED — catch-all or lineage-ambiguous or not in 35-type set
    "immune cell": None,
    "leukocyte": None,
    "lymphocyte": None,
    "mucosal invariant T cell": None,
    "gamma-delta T cell": None,
    "plasmacytoid dendritic cell": None,
    "erythrocyte": None,
    "platelet": None,
    "megakaryocyte": None,
    "erythroid lineage cell": None,
    "mast cell": None,
    "basophil": None,
    "innate lymphoid cell": None,
    "group 2 innate lymphoid cell": None,
    "group 3 innate lymphoid cell": None,
}

# CATCH-ALL LABELS — always excluded from coherent count
CATCH_ALL_LABELS = {
    "T cell", "lymphocyte", "immune cell", "leukocyte",
    "myeloid cell", "mononuclear cell",
}


def coarsen(fine_label):
    """Map a fine Census cell type label to our 35-type label.
    Returns (target_label, is_catch_all) or (None, False) if excluded."""
    if fine_label in CATCH_ALL_LABELS:
        return None, True
    if fine_label in COARSE_MAP:
        return COARSE_MAP[fine_label], False
    if fine_label in CROSS_SPECIES_TYPES:
        return fine_label, False
    return None, False


# Build per-cell-type, per-tissue counts
print("Building per-cell-type, per-tissue inventory...")

# COVID: {coarse_type: {tissue: {"cells": n, "datasets": set, "donors": set}}}
covid_inventory = defaultdict(lambda: defaultdict(lambda: {
    "cells": 0, "datasets": set(), "donors": set()
}))
# Normal: same structure
normal_inventory = defaultdict(lambda: defaultdict(lambda: {
    "cells": 0, "datasets": set(), "donors": set()
}))

excluded_covid = defaultdict(int)
unmapped_covid = defaultdict(int)
excluded_normal = defaultdict(int)
unmapped_normal = defaultdict(int)
catchall_covid = defaultdict(int)

# Process COVID obs
for _, row in obs_covid.iterrows():
    ct = row["cell_type"]
    tissue = row["tissue"]
    did = row["dataset_id"]
    donor = row["donor_id"]
    target, is_catch_all = coarsen(ct)
    if is_catch_all:
        catchall_covid[ct] += 1
    elif target is not None:
        covid_inventory[target][tissue]["cells"] += 1
        covid_inventory[target][tissue]["datasets"].add(did)
        covid_inventory[target][tissue]["donors"].add(donor)
    elif ct in COARSE_MAP and COARSE_MAP[ct] is None:
        excluded_covid[ct] += 1
    else:
        unmapped_covid[ct] += 1

# Process normal obs
for _, row in obs_normal.iterrows():
    ct = row["cell_type"]
    tissue = row["tissue"]
    did = row["dataset_id"]
    donor = row["donor_id"]
    target, is_catch_all = coarsen(ct)
    if is_catch_all:
        pass  # ignore catch-all normals
    elif target is not None:
        normal_inventory[target][tissue]["cells"] += 1
        normal_inventory[target][tissue]["datasets"].add(did)
        normal_inventory[target][tissue]["donors"].add(donor)
    elif ct in COARSE_MAP and COARSE_MAP[ct] is None:
        excluded_normal[ct] += 1
    else:
        unmapped_normal[ct] += 1

print("Done.")

# ── 5. Per-cell-type feasibility assessment ──────────────────────────
print("\n" + "=" * 70)
print("STEP 5: Per-cell-type feasibility (PASS / MARGINAL / FAIL)")
print("=" * 70)
print("PASS: ≥500 cells in COVID AND normal (same tissue), ≥2 COVID datasets")
print("MARGINAL: ≥200 cells, 1 dataset")
print("FAIL: otherwise\n")

results = []
for ct in CROSS_SPECIES_TYPES:
    covid_total = 0
    covid_datasets_all = set()
    covid_tissues = []
    normal_total_matched = 0
    normal_datasets_all = set()
    tissue_matched = True
    tissue_details = []

    # Aggregate across tissues
    for tissue, info in covid_inventory.get(ct, {}).items():
        covid_total += info["cells"]
        covid_datasets_all.update(info["datasets"])
        covid_tissues.append(tissue)

        # Check if normal exists in SAME tissue
        normal_in_tissue = normal_inventory.get(ct, {}).get(tissue, {})
        n_normal = normal_in_tissue.get("cells", 0) if normal_in_tissue else 0
        normal_total_matched += n_normal
        if normal_in_tissue:
            normal_datasets_all.update(normal_in_tissue.get("datasets", set()))

        tissue_details.append({
            "tissue": tissue,
            "covid_cells": info["cells"],
            "covid_datasets": len(info["datasets"]),
            "covid_donors": len(info["donors"]),
            "normal_cells": n_normal,
            "tissue_match": n_normal >= 200,
        })

    # Check for tissue mismatch
    if covid_total > 0 and normal_total_matched == 0:
        tissue_matched = False

    # Determine verdict
    n_covid_ds = len(covid_datasets_all)
    if covid_total >= 500 and normal_total_matched >= 500 and n_covid_ds >= 2:
        verdict = "PASS"
    elif covid_total >= 200 and normal_total_matched >= 200:
        verdict = "MARGINAL"
    elif covid_total == 0:
        verdict = "FAIL (no COVID data)"
    else:
        verdict = "FAIL"

    residual = RESIDUAL_MAP.get(ct, None)
    rank = RANK_MAP.get(ct, None)

    results.append({
        "cell_type": ct,
        "covid_cells": covid_total,
        "normal_cells_matched": normal_total_matched,
        "covid_datasets": n_covid_ds,
        "covid_tissues": covid_tissues,
        "tissue_matched": tissue_matched,
        "tissue_details": tissue_details,
        "verdict": verdict,
        "residual": residual,
        "rank": rank,
        "is_catch_all": ct in CATCH_ALL_LABELS,
    })

results_df = pd.DataFrame(results)

# Print summary table
print(f"\n{'Cell Type':<50} {'COVID':>8} {'Normal':>8} {'DS':>4} {'Verdict':<10} "
      f"{'Resid':>6} {'Rank':>4} {'TisMatch':>8}")
print("-" * 105)
for _, row in results_df.sort_values("covid_cells", ascending=False).iterrows():
    res_str = f"{row['residual']:.1f}" if row["residual"] else "  -"
    rank_str = f"{row['rank']}" if row["rank"] else " -"
    tm = "YES" if row["tissue_matched"] else "MISMATCH" if row["covid_cells"] > 0 else "-"
    catch = " [CATCH-ALL]" if row["is_catch_all"] else ""
    print(f"{row['cell_type']:<50} {row['covid_cells']:>8,} {row['normal_cells_matched']:>8,} "
          f"{row['covid_datasets']:>4} {row['verdict']:<10} {res_str:>6} {rank_str:>4} {tm:>8}{catch}")

# Count verdicts
n_pass = len(results_df[results_df["verdict"] == "PASS"])
n_marginal = len(results_df[results_df["verdict"] == "MARGINAL"])
n_fail = len(results_df[(results_df["verdict"].str.startswith("FAIL"))])
print(f"\nPASS: {n_pass}  |  MARGINAL: {n_marginal}  |  FAIL: {n_fail}")

# ── 6. Tissue-match verification for passing types ───────────────────
print("\n" + "=" * 70)
print("STEP 6: Tissue-match verification (passing types only)")
print("=" * 70)

passing = results_df[results_df["verdict"] == "PASS"]
tissue_warnings = []
for _, row in passing.iterrows():
    ct = row["cell_type"]
    print(f"\n  {ct}:")
    for td in row["tissue_details"]:
        status = "OK" if td["tissue_match"] else "WARNING: <200 normal"
        print(f"    {td['tissue']}: COVID {td['covid_cells']:,} "
              f"({td['covid_datasets']} DS, {td['covid_donors']} donors) "
              f"| Normal {td['normal_cells']:,} | {status}")
        if not td["tissue_match"] and td["covid_cells"] >= 200:
            tissue_warnings.append((ct, td["tissue"]))

if tissue_warnings:
    print(f"\n  TISSUE WARNINGS ({len(tissue_warnings)}):")
    for ct, tissue in tissue_warnings:
        print(f"    {ct} in {tissue}: normal count < 200")
else:
    print("\n  All passing types have tissue-matched normal controls.")

# ── 7. Maximum coherent n ─────────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 7: Maximum coherent n for Spearman correlation")
print("=" * 70)
print("Exclude catch-all categories, require ≥2 independent datasets,")
print("require tissue-matched normal ≥500.\n")

# "T cell" is in the 35-type set but is a catch-all
coherent_mask = (
    (results_df["verdict"] == "PASS") &
    (~results_df["cell_type"].isin(CATCH_ALL_LABELS)) &
    (results_df["tissue_matched"])
)
coherent = results_df[coherent_mask].copy()

print(f"{'Cell Type':<50} {'COVID':>8} {'Normal':>8} {'DS':>4} {'Resid':>6} {'Rank':>4}")
print("-" * 85)
for _, row in coherent.sort_values("rank").iterrows():
    res_str = f"{row['residual']:.1f}" if row["residual"] else "  -"
    rank_str = f"{row['rank']}" if row["rank"] else " -"
    print(f"{row['cell_type']:<50} {row['covid_cells']:>8,} "
          f"{row['normal_cells_matched']:>8,} {row['covid_datasets']:>4} "
          f"{res_str:>6} {rank_str:>4}")

n_coherent = len(coherent)
print(f"\nMaximum coherent n = {n_coherent}")

if n_coherent >= 2:
    residuals = coherent["residual"].dropna().values
    print(f"Residual range: {residuals.min():.2f} – {residuals.max():.2f} "
          f"(dynamic range ratio: {residuals.max()/residuals.min():.1f}×)")

# Spearman power assessment
print(f"\nSpearman power assessment:")
if n_coherent >= 10:
    print(f"  n={n_coherent}: GOOD — can detect ρ≥0.60 at p<0.05")
elif n_coherent >= 8:
    print(f"  n={n_coherent}: MARGINAL — need ρ≥0.74 for p<0.05")
else:
    print(f"  n={n_coherent}: UNDERPOWERED — need very large effects")

# ── 8. Top 10 COVID datasets by cell count ────────────────────────────
print("\n" + "=" * 70)
print("STEP 8: Top 10 COVID datasets by cell count (study dominance check)")
print("=" * 70)

ds_counts = obs_covid.groupby("dataset_id").agg(
    cells=("cell_type", "size"),
    donors=("donor_id", "nunique"),
    tissues=("tissue", lambda x: list(x.unique())),
    assays=("assay", lambda x: list(x.unique())),
    n_tissues=("tissue", "nunique"),
).sort_values("cells", ascending=False)

total_cells = len(obs_covid)
print(f"\nTotal COVID cells: {total_cells:,} across {n_covid_datasets} datasets\n")

print(f"{'Rank':>4} {'Dataset ID':<40} {'Cells':>10} {'%Total':>7} {'Donors':>7} {'Tissues':>40}")
print("-" * 115)
cumulative = 0
for i, (did, row) in enumerate(ds_counts.head(10).iterrows()):
    pct = 100 * row["cells"] / total_cells
    cumulative += pct
    tissues_str = ", ".join(str(t) for t in row["tissues"][:3])
    if len(row["tissues"]) > 3:
        tissues_str += f"... +{len(row['tissues'])-3}"
    print(f"{i+1:>4} {did:<40} {row['cells']:>10,} {pct:>6.1f}% {row['donors']:>7} "
          f"{tissues_str:<40}")

print(f"\n  Top 1 dataset: {ds_counts.iloc[0]['cells']:,} cells "
      f"({100*ds_counts.iloc[0]['cells']/total_cells:.1f}%)")
print(f"  Top 3 datasets: {ds_counts.head(3)['cells'].sum():,} cells "
      f"({100*ds_counts.head(3)['cells'].sum()/total_cells:.1f}%)")
print(f"  Top 5 datasets: {ds_counts.head(5)['cells'].sum():,} cells "
      f"({100*ds_counts.head(5)['cells'].sum()/total_cells:.1f}%)")
print(f"  Top 10 datasets: {ds_counts.head(10)['cells'].sum():,} cells "
      f"({100*ds_counts.head(10)['cells'].sum()/total_cells:.1f}%)")

if ds_counts.iloc[0]["cells"] / total_cells > 0.5:
    print("\n  WARNING: Single dataset dominates (>50%). Multi-study coverage at risk.")
elif ds_counts.head(3)["cells"].sum() / total_cells > 0.8:
    print("\n  CAUTION: Top 3 datasets account for >80%. Moderate concentration.")
else:
    print("\n  GOOD: No single study dominates. Multi-study coverage confirmed.")

# ── 9. Severity metadata check ───────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 9: Severity metadata availability")
print("=" * 70)

# Check disease labels for severity encoding
print("\n--- All disease labels (may encode severity) ---")
all_disease_labels_raw = obs_covid["disease"].value_counts()
all_disease_labels = all_disease_labels_raw[all_disease_labels_raw > 0]
for label, count in all_disease_labels.items():
    print(f"  {label}: {count:,} cells")

# Check for severity in suspension_type or development_stage
print("\n--- Suspension types ---")
for st, count in obs_covid["suspension_type"].value_counts().items():
    if count > 0:
        print(f"  {st}: {count:,} cells")

# Try to find severity in obs metadata — Census may not have it directly
# but some datasets encode it in disease or other fields
severity_keywords = ["mild", "moderate", "severe", "critical",
                     "hospitalized", "ventilat", "asymptomatic", "ICU"]
print("\n--- Severity keyword search in disease labels ---")
has_severity = False
# Only check labels with >0 cells (disease column is categorical with all values)
nonzero_labels = {label: count for label, count in all_disease_labels.items() if count > 0}
for label, count in nonzero_labels.items():
    for kw in severity_keywords:
        if kw.lower() in label.lower():
            print(f"  FOUND: '{label}' contains '{kw}' ({count:,} cells)")
            has_severity = True

if not has_severity:
    print("  No severity information found in disease labels.")
    print("  Severity may be in dataset-level metadata (not cell-level in Census).")
    print("  Would need to check individual dataset publications for severity annotations.")

census.close()

# ── 10. Save comprehensive report ─────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 10: Saving inventory report")
print("=" * 70)

# Save CSV
csv_rows = []
for _, row in results_df.iterrows():
    csv_rows.append({
        "cell_type": row["cell_type"],
        "covid_cells": row["covid_cells"],
        "normal_cells_matched": row["normal_cells_matched"],
        "covid_datasets": row["covid_datasets"],
        "tissues": "; ".join(row["covid_tissues"]) if row["covid_tissues"] else "",
        "tissue_matched": row["tissue_matched"],
        "verdict": row["verdict"],
        "residual": row["residual"],
        "rank": row["rank"],
    })
pd.DataFrame(csv_rows).to_csv(f"{OUTPUT_DIR}/covid_inventory.csv", index=False)
print(f"Saved: {OUTPUT_DIR}/covid_inventory.csv")

# Save text report
with open(f"{OUTPUT_DIR}/covid_inventory.txt", "w") as f:
    f.write("COVID-19 DISEASE REPLICATION — DETAILED FEASIBILITY INVENTORY\n")
    f.write(f"Generated: 2026-03-14\n")
    f.write(f"Census version: stable (2025-11-08)\n")
    f.write("=" * 70 + "\n\n")

    # Executive summary
    f.write("EXECUTIVE SUMMARY\n")
    f.write("-" * 40 + "\n")

    # Determine overall feasibility
    if n_coherent >= 15:
        feasibility = "STRONG"
    elif n_coherent >= 10:
        feasibility = "GOOD"
    elif n_coherent >= 8:
        feasibility = "MARGINAL"
    else:
        feasibility = "WEAK"

    f.write(f"Feasibility: {feasibility}\n")
    f.write(f"- {n_coherent} coherent cell types pass all criteria\n")
    f.write(f"  (≥500 COVID + ≥500 normal, tissue-matched, ≥2 datasets, no catch-alls)\n")
    f.write(f"- Maximum Spearman n = {n_coherent}\n")
    f.write(f"- Total COVID cells: {len(obs_covid):,}\n")
    f.write(f"- COVID donors: {n_covid_donors}\n")
    f.write(f"- COVID datasets: {n_covid_datasets}\n")
    f.write(f"- COVID tissues: {len(all_covid_tissues)}\n")
    if n_coherent >= 2:
        residuals = coherent["residual"].dropna().values
        f.write(f"- Residual range: {residuals.min():.2f} – {residuals.max():.2f} "
                f"(dynamic range: {residuals.max()/residuals.min():.1f}×)\n")
    f.write(f"- PASS: {n_pass} | MARGINAL: {n_marginal} | FAIL: {n_fail}\n")
    f.write(f"- Study dominance: top dataset = "
            f"{100*ds_counts.iloc[0]['cells']/total_cells:.1f}% of total\n")
    sev_status = "YES — in disease labels" if has_severity else "NO — not in Census cell-level metadata"
    f.write(f"- Severity metadata: {sev_status}\n")

    comparison = f"vs lupus n=8 (single dataset, blood only)"
    f.write(f"\nComparison: {comparison}\n")
    f.write(f"Recommendation: {'GO — proceed to download' if n_coherent >= 10 else 'CONDITIONAL GO' if n_coherent >= 8 else 'NEEDS REVIEW'}\n\n")

    # COVID data overview
    f.write("COVID-19 DATA\n")
    f.write("-" * 40 + "\n")
    for disease, count in all_disease_labels.items():
        if count > 0:
            f.write(f"Disease label: '{disease}' — {count:,} cells\n")
    f.write(f"Total cells: {len(obs_covid):,}\n")
    f.write(f"Donors: {n_covid_donors}\n")
    f.write(f"Datasets: {n_covid_datasets}\n")
    f.write(f"Tissues: {', '.join(str(t) for t in all_covid_tissues)}\n")
    f.write(f"Assays: {', '.join(str(a) for a in obs_covid['assay'].unique())}\n\n")

    # Normal controls
    f.write("MATCHED NORMAL CONTROLS\n")
    f.write("-" * 40 + "\n")
    f.write(f"Total cells: {len(obs_normal):,}\n")
    f.write(f"Donors: {obs_normal['donor_id'].nunique()}\n")
    f.write(f"Datasets: {obs_normal['dataset_id'].nunique()}\n")
    f.write(f"Tissues: {', '.join(str(t) for t in obs_normal['tissue'].unique()[:10])}\n")
    f.write(f"Adequacy: {'GOOD' if len(obs_normal) > 100000 else 'ADEQUATE' if len(obs_normal) > 10000 else 'LIMITED'}\n\n")

    # Full inventory table
    f.write("COARSE-MAPPED INVENTORY (35-type set)\n")
    f.write("Mapping rules: separate CD4/CD8, separate macrophage/monocyte,\n")
    f.write("NO catch-all categories (T cell, lymphocyte, immune cell excluded).\n")
    f.write("-" * 40 + "\n")
    f.write(f"{'Cell Type':<50} {'COVID':>8} {'Normal':>8} {'DS':>4} "
            f"{'Verdict':<10} {'Resid':>6} {'Rank':>4} {'TisMatch':>8}\n")
    f.write("-" * 105 + "\n")
    for _, row in results_df.sort_values("covid_cells", ascending=False).iterrows():
        res_str = f"{row['residual']:.1f}" if row["residual"] else "  -"
        rank_str = f"{int(row['rank'])}" if row["rank"] else " -"
        tm = "YES" if row["tissue_matched"] else "MISMATCH" if row["covid_cells"] > 0 else "-"
        catch = " [CATCH-ALL]" if row["is_catch_all"] else ""
        f.write(f"{row['cell_type']:<50} {row['covid_cells']:>8,} "
                f"{row['normal_cells_matched']:>8,} {row['covid_datasets']:>4} "
                f"{row['verdict']:<10} {res_str:>6} {rank_str:>4} {tm:>8}{catch}\n")
    f.write(f"\nPASS: {n_pass}  |  MARGINAL: {n_marginal}  |  FAIL: {n_fail}\n\n")

    # Tissue-match details for passing types
    f.write("TISSUE-MATCH DETAILS (passing types)\n")
    f.write("-" * 40 + "\n")
    for _, row in passing.iterrows():
        ct = row["cell_type"]
        f.write(f"\n  {ct}:\n")
        for td in row["tissue_details"]:
            status = "OK" if td["tissue_match"] else "WARNING: <200 normal"
            f.write(f"    {td['tissue']}: COVID {td['covid_cells']:,} "
                    f"({td['covid_datasets']} DS, {td['covid_donors']} donors) "
                    f"| Normal {td['normal_cells']:,} | {status}\n")
    if tissue_warnings:
        f.write(f"\n  TISSUE WARNINGS ({len(tissue_warnings)}):\n")
        for ct, tissue in tissue_warnings:
            f.write(f"    {ct} in {tissue}: normal count < 200\n")
    else:
        f.write(f"\n  All passing types have tissue-matched normal controls.\n")
    f.write("\n")

    # Coherent set
    f.write("MAXIMUM COHERENT SET (for Spearman correlation)\n")
    f.write("-" * 40 + "\n")
    f.write(f"Criteria: PASS verdict + no catch-alls + tissue-matched\n")
    f.write(f"n = {n_coherent}\n\n")
    f.write(f"{'Cell Type':<50} {'COVID':>8} {'Normal':>8} {'DS':>4} {'Resid':>6} {'Rank':>4}\n")
    f.write("-" * 85 + "\n")
    for _, row in coherent.sort_values("rank").iterrows():
        res_str = f"{row['residual']:.1f}" if row["residual"] else "  -"
        rank_str = f"{int(row['rank'])}" if row["rank"] else " -"
        f.write(f"{row['cell_type']:<50} {row['covid_cells']:>8,} "
                f"{row['normal_cells_matched']:>8,} {row['covid_datasets']:>4} "
                f"{res_str:>6} {rank_str:>4}\n")
    if n_coherent >= 2:
        residuals = coherent["residual"].dropna().values
        f.write(f"\nResidual range: {residuals.min():.2f} – {residuals.max():.2f} "
                f"(dynamic range: {residuals.max()/residuals.min():.1f}×)\n")
    f.write("\n")

    # Top datasets
    f.write("TOP 10 COVID DATASETS BY CELL COUNT\n")
    f.write("-" * 40 + "\n")
    f.write(f"{'Rank':>4} {'Dataset ID':<40} {'Cells':>10} {'%Total':>7} {'Donors':>7} {'Tissues':<30}\n")
    f.write("-" * 100 + "\n")
    for i, (did, drow) in enumerate(ds_counts.head(10).iterrows()):
        pct = 100 * drow["cells"] / total_cells
        tissues_str = ", ".join(str(t) for t in drow["tissues"][:3])
        if len(drow["tissues"]) > 3:
            tissues_str += f"... +{len(drow['tissues'])-3}"
        f.write(f"{i+1:>4} {did:<40} {drow['cells']:>10,} {pct:>6.1f}% "
                f"{drow['donors']:>7} {tissues_str:<30}\n")
    f.write(f"\nTop 1 = {100*ds_counts.iloc[0]['cells']/total_cells:.1f}%, "
            f"Top 3 = {100*ds_counts.head(3)['cells'].sum()/total_cells:.1f}%, "
            f"Top 5 = {100*ds_counts.head(5)['cells'].sum()/total_cells:.1f}%, "
            f"Top 10 = {100*ds_counts.head(10)['cells'].sum()/total_cells:.1f}%\n\n")

    # Severity metadata
    f.write("SEVERITY METADATA\n")
    f.write("-" * 40 + "\n")
    f.write(f"Severity in disease labels: {'YES' if has_severity else 'NO'}\n")
    if has_severity:
        for label in all_disease_labels.index:
            for kw in severity_keywords:
                if kw.lower() in label.lower():
                    f.write(f"  '{label}': {all_disease_labels[label]:,} cells (keyword: {kw})\n")
    else:
        f.write("  Severity not encoded in Census cell-level metadata.\n")
        f.write("  Would need to check individual dataset publications.\n")
        f.write("  If available per-dataset, could add as stretch goal.\n")
    f.write("\n")

    # Unmapped / excluded types
    if unmapped_covid:
        f.write("UNMAPPED COVID CELL TYPES (>50 cells, not in 35-type set)\n")
        f.write("-" * 40 + "\n")
        for ct, count in sorted(unmapped_covid.items(), key=lambda x: -x[1]):
            if count > 50:
                f.write(f"  {ct}: {count:,}\n")
        f.write("\n")

    if excluded_covid:
        f.write("EXCLUDED COVID CELL TYPES (in coarse map as None)\n")
        f.write("-" * 40 + "\n")
        for ct, count in sorted(excluded_covid.items(), key=lambda x: -x[1]):
            if count > 50:
                f.write(f"  {ct}: {count:,}\n")
        f.write("\n")

    if catchall_covid:
        f.write("CATCH-ALL / LINEAGE-AMBIGUOUS COVID CELL TYPES (excluded)\n")
        f.write("-" * 40 + "\n")
        for ct, count in sorted(catchall_covid.items(), key=lambda x: -x[1]):
            f.write(f"  {ct}: {count:,}\n")
        f.write("\n")

    # Statistical assessment
    f.write("STATISTICAL ASSESSMENT\n")
    f.write("-" * 40 + "\n")
    f.write(f"COVID coherent n = {n_coherent}\n")
    f.write(f"Lupus coherent n = 8 (for comparison)\n")
    f.write(f"Cancer coherent n = 13 (for comparison)\n\n")
    if n_coherent >= 10:
        f.write(f"Spearman with n={n_coherent}: can detect ρ≥0.60 at p<0.05 (two-tailed).\n")
        f.write(f"Substantially better powered than lupus (n=8, need ρ≥0.74).\n")
    elif n_coherent >= 8:
        f.write(f"Spearman with n={n_coherent}: critical ρ for p<0.05 ≈ 0.74.\n")
        f.write(f"Similar power to lupus — marginal.\n")
    f.write(f"\nMulti-study coverage: {n_covid_datasets} datasets (vs lupus 1, cancer ~3).\n")
    f.write(f"Multi-tissue coverage: {len(all_covid_tissues)} tissues (vs lupus 1, cancer 1).\n")
    f.write(f"Multi-donor coverage: {n_covid_donors} donors (vs lupus 162, cancer ~290).\n\n")

    # Recommendation
    f.write("RECOMMENDATION\n")
    f.write("-" * 40 + "\n")
    if n_coherent >= 10:
        f.write(f"1. COVID-19 is a STRONG candidate for disease-axis replication.\n")
        f.write(f"   n={n_coherent} coherent types, multi-study, multi-tissue.\n")
        f.write(f"2. Substantially superior to lupus (n=8, single study, blood only).\n")
        f.write(f"3. Comparable or superior to cancer (n=13, single tissue).\n")
        f.write(f"4. PROCEED TO DOWNLOAD: design download script matching cancer pipeline.\n")
    elif n_coherent >= 8:
        f.write(f"1. COVID-19 is MARGINAL — n={n_coherent}, similar to lupus.\n")
        f.write(f"2. Multi-study advantage exists but cell type coverage is limited.\n")
        f.write(f"3. Consider whether marginal types can be promoted.\n")
    else:
        f.write(f"1. COVID-19 has insufficient coherent cell types (n={n_coherent}).\n")
        f.write(f"2. Consider alternative disease states.\n")

print(f"Saved: {OUTPUT_DIR}/covid_inventory.txt")

print("\n" + "=" * 70)
print("COVID-19 INVENTORY COMPLETE")
print("=" * 70)
