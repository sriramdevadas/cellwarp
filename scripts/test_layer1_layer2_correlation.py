#!/usr/bin/env python3
"""
ANALYSIS-B: Layer 1 vs Layer 2 Per-Type Correlation

Tests whether per-type Procrustes residuals (Layer 1 — centroid conservation)
correlate with per-type ellipsoid alignment scores (Layer 2 — covariance
orientation conservation). The manuscript claims the two layers are
"geometrically distinct"; this tests whether they are also statistically
uncorrelated at the per-type level.

Biology: Layer 1 measures where each cell type sits relative to all others
(centroid position conservation). Layer 2 measures how cells spread within
each type (covariance ellipsoid orientation conservation). If these are
uncorrelated, a cell type can be rigid in one layer and flexible in the
other — the two layers capture genuinely independent aspects of identity
geometry.

Math: Spearman ρ between per-type Procrustes residual magnitude (lower =
more conserved centroid) and per-type Krzanowski S score at k=3 (higher =
more aligned ellipsoid). A non-significant correlation supports the
"distinct" claim. Note: the sign expectation is ambiguous — a negative ρ
would mean types with good centroid conservation also have good ellipsoid
conservation, while a positive ρ would mean the two trade off.

Output: output/validation/layer_correlation/
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT = Path(__file__).resolve().parent.parent
OUTPUT_SCALED = PROJECT / "output" / "phase2" / "scaled_35types"
ELLIPSOID_DIR = PROJECT / "output" / "mechanistic" / "ellipsoid_alignment"
OUTPUT_DIR = PROJECT / "output" / "validation" / "layer_correlation"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("ANALYSIS-B: Layer 1 vs Layer 2 Per-Type Correlation")
    print("=" * 70)

    # --- Load Layer 1: per-type Procrustes residual magnitudes ---
    with open(OUTPUT_SCALED / "procrustes_results_35.json") as f:
        proc_results = json.load(f)

    cell_types = proc_results["cell_types"]
    layer1 = pd.DataFrame([
        {
            "cell_type": ct,
            "residual_magnitude": proc_results["residuals"][ct]["magnitude"],
        }
        for ct in cell_types
    ])

    print(f"\nLayer 1 (centroid residuals): {len(layer1)} cell types")
    print(f"  Range: {layer1.residual_magnitude.min():.4f} – "
          f"{layer1.residual_magnitude.max():.4f}")
    print(f"  Mean:  {layer1.residual_magnitude.mean():.4f}")

    # --- Load Layer 2: per-type ellipsoid alignment scores ---
    alignment_df = pd.read_csv(ELLIPSOID_DIR / "35type_alignment_scores.csv")

    # Test at multiple k values for completeness
    k_values = [1, 3, 5]

    print(f"\nLayer 2 (ellipsoid alignment S scores):")
    for k in k_values:
        sub = alignment_df[alignment_df.k == k]
        print(f"  k={k}: S_pre range [{sub.S_pre.min():.4f}, {sub.S_pre.max():.4f}], "
              f"mean={sub.S_pre.mean():.4f}")

    # --- Compute correlations ---
    print("\n" + "=" * 70)
    print("RESULTS: Spearman correlation (Layer 1 residual vs Layer 2 S score)")
    print("=" * 70)

    results_records = []

    for k in k_values:
        for metric_name, metric_col, description in [
            ("S_pre", "S_pre", "Pre-Procrustes ellipsoid alignment"),
            ("S_post", "S_post", "Post-Procrustes ellipsoid alignment"),
        ]:
            sub = alignment_df[alignment_df.k == k][["cell_type", metric_col]].copy()
            merged = layer1.merge(sub, on="cell_type", how="inner")
            assert len(merged) == 35, f"Expected 35 types, got {len(merged)}"

            rho, p = stats.spearmanr(
                merged["residual_magnitude"],
                merged[metric_col],
            )

            results_records.append({
                "k": k,
                "metric": metric_name,
                "description": description,
                "spearman_rho": rho,
                "p_value": p,
                "n": len(merged),
                "significant_005": p < 0.05,
            })

            sig_str = " *" if p < 0.05 else ""
            print(f"\n  k={k}, {metric_name} ({description}):")
            print(f"    Spearman ρ = {rho:+.4f}")
            print(f"    p-value    = {p:.4f}{sig_str}")
            print(f"    n          = {len(merged)}")

    # --- Primary result: S_pre at k=3 ---
    primary = [r for r in results_records
               if r["k"] == 3 and r["metric"] == "S_pre"][0]

    print("\n" + "=" * 70)
    print("PRIMARY RESULT (S_pre, k=3)")
    print("=" * 70)
    print(f"  Spearman ρ = {primary['spearman_rho']:+.4f}")
    print(f"  p-value    = {primary['p_value']:.4f}")
    print(f"  n          = {primary['n']}")

    if primary["p_value"] >= 0.05:
        print("\n  → Per-type Layer 1 and Layer 2 scores are NOT significantly correlated.")
        print("    The 'geometrically distinct' language is fully supported.")
        print("    DECISION RULE: Add one sentence to Results confirming statistical")
        print("    uncorrelation at the per-type level.")
    else:
        print(f"\n  → Per-type scores ARE significantly correlated (ρ={primary['spearman_rho']:+.3f}).")
        print("    The 'geometrically distinct' language holds (different optimal rotations)")
        print("    but per-type scores show partial correlation.")
        print("    DECISION RULE: Report back to ADVISOR — do not edit text.")

    # --- Save per-type data for inspection ---
    sub_k3 = alignment_df[alignment_df.k == 3][["cell_type", "S_pre", "S_post"]].copy()
    merged_full = layer1.merge(sub_k3, on="cell_type", how="inner")
    merged_full = merged_full.sort_values("residual_magnitude")

    # Add rigidity rank (1 = most flexible / largest residual)
    merged_full["rigidity_rank"] = merged_full["residual_magnitude"].rank(
        ascending=False
    ).astype(int)

    print("\n  Per-type data (sorted by residual magnitude, ascending = more rigid):")
    print(f"  {'Rank':>4} {'Cell Type':<50} {'Resid':>8} {'S_pre':>8} {'S_post':>8}")
    print("  " + "-" * 82)
    for _, row in merged_full.iterrows():
        print(f"  {row.rigidity_rank:>4} {row.cell_type:<50} "
              f"{row.residual_magnitude:>8.4f} {row.S_pre:>8.4f} {row.S_post:>8.4f}")

    merged_full.to_csv(OUTPUT_DIR / "layer1_vs_layer2_per_type.csv", index=False)

    # --- Save results ---
    output = {
        "test": "layer1_vs_layer2_correlation",
        "layer1_metric": "Procrustes residual magnitude (lower = more conserved centroid)",
        "layer2_metric": "Krzanowski S score (higher = more aligned ellipsoid)",
        "primary_result": {
            "k": 3,
            "metric": "S_pre",
            "spearman_rho": float(primary["spearman_rho"]),
            "p_value": float(primary["p_value"]),
            "n": int(primary["n"]),
            "significant_at_005": bool(primary["significant_005"]),
        },
        "all_results": [
            {k: (float(v) if isinstance(v, (np.floating, float)) else
                 bool(v) if isinstance(v, (np.bool_,)) else
                 int(v) if isinstance(v, (np.integer,)) else v)
             for k, v in r.items()}
            for r in results_records
        ],
    }

    results_path = OUTPUT_DIR / "layer_correlation_results.json"
    with open(results_path, "w") as f:
        json.dump(output, f, indent=2)

    print(f"\nResults saved: {results_path}")
    print(f"Per-type data saved: {OUTPUT_DIR / 'layer1_vs_layer2_per_type.csv'}")


if __name__ == "__main__":
    main()
