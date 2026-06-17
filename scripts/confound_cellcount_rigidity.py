#!/usr/bin/env python3
"""
CellWarp — Cell Count vs Rigidity Confound Test

Tests whether the rigidity ranking (Procrustes residual magnitude) is confounded
by per-cell-type cell count. The concern: cell types with fewer cells have noisier
centroids, which could inflate residual distances and make them appear more "flexible."

Biology: If this confound exists, the rigidity ranking reflects statistical noise
rather than genuine evolutionary divergence. If absent, centroid estimation is
sufficiently stable at the 500–2000 cell range.

Math: Spearman rank correlation between post-cap cell count and residual magnitude
(or rigidity rank). Partial correlation controlling for lineage (immune vs non-immune).

Inputs:
    output/phase2/scaled_35types/residuals_ranked.csv
    data/phase2_scaled/human_scaled.h5ad
    data/phase2_scaled/mouse_scaled.h5ad

Outputs:
    Printed report (no files saved — audit task)
"""

import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats

# ─── File paths ───────────────────────────────────────────────────────────────
RANKING_PATH = "output/phase2/scaled_35types/residuals_ranked.csv"
HUMAN_H5AD = "data/phase2_scaled/human_scaled.h5ad"
MOUSE_H5AD = "data/phase2_scaled/mouse_scaled.h5ad"

# ─── Immune vs non-immune classification ──────────────────────────────────────
# Based on standard lineage assignments for the 35 cell types
IMMUNE_TYPES = {
    "B cell", "T cell", "macrophage", "classical monocyte", "monocyte",
    "natural killer cell", "mature NK T cell", "hematopoietic precursor cell",
    "hematopoietic stem cell", "CD8-positive, alpha-beta T cell",
    "CD4-positive, alpha-beta T cell", "intermediate monocyte",
    "non-classical monocyte", "plasma cell", "myeloid dendritic cell",
    "myeloid leukocyte", "neutrophil", "granulocyte",
}

def get_cell_counts(h5ad_path, cell_type_col="cell_type"):
    """Get per-cell-type cell counts from an h5ad file."""
    adata = ad.read_h5ad(h5ad_path, backed="r")
    counts = adata.obs[cell_type_col].value_counts()
    adata.file.close()
    return counts

def partial_spearman(x, y, z):
    """
    Partial Spearman correlation between x and y, controlling for z.

    Math: Convert to ranks, then compute partial Pearson correlation on ranks.
    r_xy.z = (r_xy - r_xz * r_yz) / sqrt((1 - r_xz^2) * (1 - r_yz^2))

    For significance: use t-test with n-3 degrees of freedom.
    t = r_partial * sqrt((n - 3) / (1 - r_partial^2))
    """
    # Convert to ranks for Spearman
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    rz = stats.rankdata(z)

    # Pearson correlations on ranks
    r_xy = np.corrcoef(rx, ry)[0, 1]
    r_xz = np.corrcoef(rx, rz)[0, 1]
    r_yz = np.corrcoef(ry, rz)[0, 1]

    # Partial correlation
    numerator = r_xy - r_xz * r_yz
    denominator = np.sqrt((1 - r_xz**2) * (1 - r_yz**2))

    if denominator == 0:
        return np.nan, np.nan

    r_partial = numerator / denominator

    # t-test for significance
    n = len(x)
    df = n - 3  # controlling for 1 variable
    if abs(r_partial) >= 1.0:
        p_val = 0.0
    else:
        t_stat = r_partial * np.sqrt(df / (1 - r_partial**2))
        p_val = 2 * stats.t.sf(abs(t_stat), df)

    return r_partial, p_val

def main():
    print("=" * 72)
    print("CONFOUND TEST: Cell Count vs Rigidity Ranking")
    print("=" * 72)

    # ─── Step 1: Load rigidity ranking ────────────────────────────────────
    print(f"\n[1] Loading rigidity ranking from: {RANKING_PATH}")
    ranking = pd.read_csv(RANKING_PATH)
    print(f"    Loaded {len(ranking)} cell types")
    print(f"    Ranking direction: rank 1 = MOST FLEXIBLE (highest residual)")
    print(f"                       rank {len(ranking)} = MOST RIGID (lowest residual)")

    # ─── Step 2: Load cell counts from h5ad files ─────────────────────────
    print(f"\n[2] Loading cell counts from h5ad files")
    print(f"    Human: {HUMAN_H5AD}")
    human_counts = get_cell_counts(HUMAN_H5AD)
    print(f"    Mouse: {MOUSE_H5AD}")
    mouse_counts = get_cell_counts(MOUSE_H5AD)

    # Build count table for the 35 ranked types
    count_data = []
    for _, row in ranking.iterrows():
        ct = row["cell_type"]
        h_count = human_counts.get(ct, 0)
        m_count = mouse_counts.get(ct, 0)
        count_data.append({
            "cell_type": ct,
            "human_cells": int(h_count),
            "mouse_cells": int(m_count),
            "min_cells": int(min(h_count, m_count)),
            "total_cells": int(h_count + m_count),
        })
    counts_df = pd.DataFrame(count_data)

    # Merge with ranking
    df = ranking.merge(counts_df, on="cell_type")
    assert len(df) == 35, f"Expected 35 cell types, got {len(df)}"

    # Add lineage classification
    df["is_immune"] = df["cell_type"].isin(IMMUNE_TYPES).astype(int)
    print(f"    Immune types: {df['is_immune'].sum()}")
    print(f"    Non-immune types: {(1 - df['is_immune']).sum()}")

    # ─── Step 3: Spearman correlations ────────────────────────────────────
    print("\n" + "=" * 72)
    print("RESULTS")
    print("=" * 72)

    # Test with multiple cell count measures
    count_measures = {
        "min_cells (bottleneck species)": df["min_cells"].values,
        "human_cells": df["human_cells"].values,
        "mouse_cells": df["mouse_cells"].values,
        "total_cells (both species)": df["total_cells"].values,
    }

    # Primary analysis: residual_magnitude vs cell count
    # Note: Higher residual = more flexible = higher rank number... wait,
    # rank 1 = highest residual (most flexible). So rank and residual_magnitude
    # go in opposite directions. Let's use residual_magnitude directly.

    residual_mag = df["residual_magnitude"].values
    ranks = df["rank"].values  # rank 1 = most flexible (highest residual)

    print("\n── 3a. Spearman: residual_magnitude vs cell count ──")
    print("    (Positive ρ → more cells = higher residual = more flexible)")
    print("    (Negative ρ → more cells = lower residual = more rigid)")
    print("    (If confound exists: negative ρ expected — fewer cells → inflated residual)\n")

    primary_rho = None
    primary_p = None

    for label, count_vals in count_measures.items():
        rho, p = stats.spearmanr(count_vals, residual_mag)
        print(f"    {label}:")
        print(f"      Spearman ρ = {rho:.3f},  p = {p:.3f},  n = {len(count_vals)}")
        if label == "min_cells (bottleneck species)":
            primary_rho = rho
            primary_p = p
            if p < 0.05:
                direction = "higher" if rho > 0 else "lower"
                print(f"      *** SIGNIFICANT (p < 0.05) ***")
            else:
                print(f"      Not significant (p ≥ 0.05)")
        print()

    # ─── Step 4: Report card ──────────────────────────────────────────────
    print("\n── 4. Primary result (min_cells vs residual_magnitude) ──\n")
    print(f"    Spearman ρ  = {primary_rho:.3f}")
    print(f"    p-value     = {primary_p:.3f}")
    print(f"    n           = 35")

    if primary_rho > 0:
        print(f"    Direction   : Higher cell count is associated with HIGHER residual magnitude (more flexible)")
        print(f"                  (OPPOSITE to confound prediction)")
    elif primary_rho < 0:
        print(f"    Direction   : Higher cell count is associated with LOWER residual magnitude (more rigid)")
        print(f"                  (CONSISTENT with confound prediction)")
    else:
        print(f"    Direction   : No association")

    if primary_p < 0.05:
        print(f"\n    *** CONFOUND PRESENT — requires paper response ***")
    else:
        print(f"\n    CONFOUND ABSENT — report ρ and p in methods or supplement")

    # ─── Step 5: Extremes ─────────────────────────────────────────────────
    print("\n── 5a. 5 cell types with LOWEST cell counts ──\n")
    bottom5 = df.nsmallest(5, "min_cells")[["cell_type", "min_cells", "human_cells", "mouse_cells", "rank", "residual_magnitude"]]
    for _, row in bottom5.iterrows():
        print(f"    {row['cell_type']:50s}  min_cells={row['min_cells']:5d}  "
              f"(H={row['human_cells']:5d}, M={row['mouse_cells']:5d})  "
              f"rank={row['rank']:2d}  residual={row['residual_magnitude']:.3f}")

    print("\n── 5b. 5 cell types with HIGHEST cell counts ──\n")
    top5 = df.nlargest(5, "min_cells")[["cell_type", "min_cells", "human_cells", "mouse_cells", "rank", "residual_magnitude"]]
    for _, row in top5.iterrows():
        print(f"    {row['cell_type']:50s}  min_cells={row['min_cells']:5d}  "
              f"(H={row['human_cells']:5d}, M={row['mouse_cells']:5d})  "
              f"rank={row['rank']:2d}  residual={row['residual_magnitude']:.3f}")

    # ─── Step 6: Partial Spearman controlling for lineage ─────────────────
    print("\n── 6. Partial Spearman: residual_magnitude vs min_cells, controlling for lineage ──\n")
    print(f"    Lineage coding: immune=1, non-immune=0")

    r_partial, p_partial = partial_spearman(
        df["min_cells"].values,
        df["residual_magnitude"].values,
        df["is_immune"].values,
    )

    print(f"    Partial Spearman ρ = {r_partial:.3f}")
    print(f"    p-value            = {p_partial:.3f}")
    print(f"    n                  = 35")
    print(f"    Covariates         = 1 (immune vs non-immune)")

    if p_partial < 0.05:
        print(f"\n    *** CONFOUND PRESENT after controlling for lineage — requires paper response ***")
    else:
        print(f"\n    CONFOUND ABSENT after controlling for lineage — report in methods or supplement")

    # ─── Summary table ────────────────────────────────────────────────────
    print("\n\n" + "=" * 72)
    print("FULL TABLE: All 35 cell types sorted by rank")
    print("=" * 72)
    print(f"\n{'Rank':>4}  {'Cell Type':50s}  {'MinCells':>8}  {'Hum':>5}  {'Mus':>5}  {'Residual':>9}  {'Immune':>6}")
    print("-" * 92)
    for _, row in df.sort_values("rank").iterrows():
        print(f"{row['rank']:4d}  {row['cell_type']:50s}  {row['min_cells']:8d}  "
              f"{row['human_cells']:5d}  {row['mouse_cells']:5d}  "
              f"{row['residual_magnitude']:9.3f}  {'Y' if row['is_immune'] else 'N':>6}")

    print("\n" + "=" * 72)
    print("END OF CONFOUND TEST REPORT")
    print("=" * 72)

if __name__ == "__main__":
    main()
