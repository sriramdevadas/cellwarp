#!/usr/bin/env python3
"""
CellWarp — Donor-Split Feasibility Analysis

Determines whether Tabula Sapiens has sufficient donor diversity per cell type
to support a donor-split within-species baseline control for Procrustes analysis.

Logic: Split human donors into two non-overlapping groups, compute centroids
independently, and run Procrustes between halves. Compare the within-species
obs/null ratio against the cross-species obs/null ratio. If cross-species
is less coherent (higher obs/null), there is an evolutionary component.

This script performs Step 1: feasibility check only. Reports donor counts,
cell distributions, tissue coverage, and classifies each cell type as
FEASIBLE or INFEASIBLE for donor splitting.

Inputs:
    data/phase2_scaled/human_scaled.h5ad  — normalized human data (35 types)

Outputs:
    analysis/donor_split/feasibility_table.csv — per-cell-type feasibility
"""

from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import h5py
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

HUMAN_DATA_PATH = Path("./data/phase2_scaled/human_scaled.h5ad")
MIN_DONORS_FOR_SPLIT = 4       # Need ≥4 donors to make two groups of ≥2
MIN_CELLS_PER_HALF = 100       # Minimum cells in smaller half for valid centroid
MIN_CELLS_PER_HALF_STRICT = 200  # Stricter threshold


def read_obs_from_h5ad(path: Path) -> pd.DataFrame:
    """Read obs metadata from h5ad via h5py (avoids anndata version issues)."""
    obs_data = {}
    with h5py.File(path, "r") as f:
        for key in f["obs"].keys():
            if key == "__categories":
                continue
            item = f["obs"][key]
            if isinstance(item, h5py.Dataset):
                obs_data[key] = item[:]
            elif isinstance(item, h5py.Group):
                if "codes" in item and "categories" in item:
                    codes = item["codes"][:]
                    cats = item["categories"][:]
                    if cats.dtype.kind in ("S", "O"):
                        cats = [c.decode() if isinstance(c, bytes) else c for c in cats]
                    obs_data[key] = [cats[c] if c >= 0 else None for c in codes]
    return pd.DataFrame(obs_data)


def main():
    print("=" * 80)
    print("DONOR-SPLIT FEASIBILITY ANALYSIS")
    print("=" * 80)

    # ------------------------------------------------------------------
    # Load human metadata only (skip expression matrix)
    # ------------------------------------------------------------------
    print("\n[1] Loading human metadata from existing h5ad...")
    obs = read_obs_from_h5ad(HUMAN_DATA_PATH)

    print(f"    Total cells: {len(obs):,}")
    print(f"    Columns: {list(obs.columns)}")

    cell_types_in_data = sorted(obs["cell_type"].unique())
    print(f"    Cell types in data: {len(cell_types_in_data)}")
    print(f"    Total unique donors: {obs['donor_id'].nunique()}")

    # Global donor overview
    print(f"\n    Donor IDs: {sorted(obs['donor_id'].unique())}")
    donor_summary = obs.groupby("donor_id").agg(
        n_cells=("cell_type", "size"),
        n_types=("cell_type", "nunique"),
        n_tissues=("tissue_general", "nunique"),
    ).sort_values("n_cells", ascending=False)
    print(f"\n    Per-donor overview:")
    for d, r in donor_summary.iterrows():
        print(f"      {d:<12} {r['n_cells']:>5} cells, {r['n_types']:>2} types, {r['n_tissues']:>2} tissues")

    # ------------------------------------------------------------------
    # Per-cell-type donor analysis
    # ------------------------------------------------------------------
    print("\n[2] Per-cell-type donor analysis")
    print("-" * 80)

    rows = []
    for ct in cell_types_in_data:
        ct_obs = obs[obs["cell_type"] == ct]
        total_cells = len(ct_obs)
        donors = ct_obs["donor_id"].unique()
        n_donors = len(donors)

        # Cells per donor
        cells_per_donor = ct_obs.groupby("donor_id").size()
        cpd_min = int(cells_per_donor.min())
        cpd_median = float(cells_per_donor.median())
        cpd_max = int(cells_per_donor.max())

        # Tissue coverage
        col = "tissue_general" if "tissue_general" in ct_obs.columns else "tissue"
        tissues = ct_obs[col].unique()
        n_tissues = len(tissues)
        tissue_list = sorted(tissues)

        # Tissue coverage per donor
        donor_tissues = ct_obs.groupby("donor_id")[col].nunique()
        mean_tissues_per_donor = donor_tissues.mean()

        # Simulate worst-case 50/50 donor split (greedy balanced partition)
        sorted_donors = cells_per_donor.sort_values(ascending=False)
        half1_cells = 0
        half2_cells = 0
        for d, c in sorted_donors.items():
            if half1_cells <= half2_cells:
                half1_cells += c
            else:
                half2_cells += c
        smaller_half = min(half1_cells, half2_cells)
        larger_half = max(half1_cells, half2_cells)

        # Feasibility classification
        if n_donors < MIN_DONORS_FOR_SPLIT:
            status = "INFEASIBLE"
            reason = f"only {n_donors} donor(s)"
        elif smaller_half < MIN_CELLS_PER_HALF:
            status = "INFEASIBLE"
            reason = f"smaller half only {smaller_half} cells"
        else:
            status = "FEASIBLE"
            reason = ""

        rows.append({
            "cell_type": ct,
            "total_cells": total_cells,
            "n_donors": n_donors,
            "cells_per_donor_min": cpd_min,
            "cells_per_donor_median": cpd_median,
            "cells_per_donor_max": cpd_max,
            "n_tissues": n_tissues,
            "tissues": "; ".join(tissue_list),
            "mean_tissues_per_donor": round(mean_tissues_per_donor, 1),
            "balanced_split_smaller_half": smaller_half,
            "balanced_split_larger_half": larger_half,
            "status": status,
            "reason": reason,
        })

    df = pd.DataFrame(rows).sort_values("n_donors", ascending=True).reset_index(drop=True)

    # ------------------------------------------------------------------
    # Print detailed table
    # ------------------------------------------------------------------
    print(f"\n{'Cell Type':<50} {'Cells':>6} {'Donors':>7} {'CPD min':>8} {'CPD med':>8} "
          f"{'CPD max':>8} {'Tissues':>8} {'Split½':>7} {'Status':<12}")
    print("-" * 130)
    for _, r in df.iterrows():
        print(f"{r['cell_type']:<50} {r['total_cells']:>6} {r['n_donors']:>7} "
              f"{r['cells_per_donor_min']:>8} {r['cells_per_donor_median']:>8.0f} "
              f"{r['cells_per_donor_max']:>8} {r['n_tissues']:>8} "
              f"{r['balanced_split_smaller_half']:>7} {r['status']:<12} {r['reason']}")

    # ------------------------------------------------------------------
    # Tissue detail per cell type
    # ------------------------------------------------------------------
    print("\n\n[3] Tissue coverage detail")
    print("-" * 80)
    for _, r in df.iterrows():
        print(f"  {r['cell_type']}: {r['tissues']}")

    # ------------------------------------------------------------------
    # Summary statistics
    # ------------------------------------------------------------------
    print("\n\n[4] SUMMARY")
    print("=" * 80)

    n_total = len(df)
    n_feasible = (df["status"] == "FEASIBLE").sum()
    n_infeasible = (df["status"] == "INFEASIBLE").sum()

    # Donor thresholds
    n_ge4 = (df["n_donors"] >= 4).sum()
    n_ge6 = (df["n_donors"] >= 6).sum()
    n_ge8 = (df["n_donors"] >= 8).sum()
    n_ge10 = (df["n_donors"] >= 10).sum()

    # Cell thresholds for split halves
    n_ge100 = ((df["n_donors"] >= MIN_DONORS_FOR_SPLIT) &
               (df["balanced_split_smaller_half"] >= 100)).sum()
    n_ge200 = ((df["n_donors"] >= MIN_DONORS_FOR_SPLIT) &
               (df["balanced_split_smaller_half"] >= 200)).sum()

    # Single/dual donor types
    single_donor = df[df["n_donors"] == 1]["cell_type"].tolist()
    dual_donor = df[df["n_donors"] == 2]["cell_type"].tolist()
    triple_donor = df[df["n_donors"] == 3]["cell_type"].tolist()

    print(f"  Total cell types analyzed: {n_total}")
    print(f"  FEASIBLE (>={MIN_DONORS_FOR_SPLIT} donors, >={MIN_CELLS_PER_HALF} cells/half): {n_feasible}")
    print(f"  INFEASIBLE: {n_infeasible}")
    print()
    print(f"  Donor thresholds:")
    print(f"    >=4 donors:  {n_ge4} / {n_total} cell types")
    print(f"    >=6 donors:  {n_ge6} / {n_total} cell types")
    print(f"    >=8 donors:  {n_ge8} / {n_total} cell types")
    print(f"    >=10 donors: {n_ge10} / {n_total} cell types")
    print()
    print(f"  Split-half cell thresholds (among types with >={MIN_DONORS_FOR_SPLIT} donors):")
    print(f"    >=100 cells per half: {n_ge100} / {n_total} cell types")
    print(f"    >=200 cells per half: {n_ge200} / {n_total} cell types")
    print()

    if single_donor:
        print(f"  WARNING: Single-donor types (CANNOT split): {single_donor}")
    if dual_donor:
        print(f"  WARNING: Dual-donor types (CANNOT split): {dual_donor}")
    if triple_donor:
        print(f"  WARNING: Triple-donor types (marginal): {triple_donor}")

    print()
    if n_feasible >= 15:
        print(f"  VERDICT: {n_feasible} types feasible -> PROCEED to Step 2")
    else:
        print(f"  VERDICT: Only {n_feasible} feasible -> INSUFFICIENT for Step 2")

    # ------------------------------------------------------------------
    # Save feasibility table
    # ------------------------------------------------------------------
    out_dir = Path("./analysis/donor_split")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "feasibility_table.csv"
    df.to_csv(out_path, index=False)
    print(f"\n  Feasibility table saved to: {out_path}")


if __name__ == "__main__":
    main()
