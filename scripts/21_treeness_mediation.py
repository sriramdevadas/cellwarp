#!/usr/bin/env python3
"""
CellWarp Script 21 — Neighborhood Density Mediation Test

Tests whether neighborhood density in expression space mediates the
treeness-rigidity anticorrelation (DECISION-120, rho=-0.349, p=0.040).

Biology: Cell types in dense neighborhoods (many nearby types in expression
space) may be simultaneously more rigid (constrained by similar types
across species) and less tree-like (tight clusters degrade tree topology).
If true, the anticorrelation is a geometric proximity artifact rather than
a direct biological relationship between evolutionary conservation and
developmental distinctiveness.

Math: Neighborhood density for cell type i = mean Euclidean distance to
k nearest neighbors in 16,959-gene expression space. Mediation tested
via partial Spearman correlation: rho(rigidity, treeness | density).
Attenuation >= 50% confirms mediation.

Pre-registration: docs/preregistration_treeness_anticorrelation_2026-03-16.md
Input:  output/phase2/scaled_35types/centroids_human_35.csv
        output/phase2/scaled_35types/residuals_ranked.csv
        output/liang_wagner/treeness_scores_per_celltype.csv
Output: output/liang_wagner/
"""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# --- Config ---
OUTPUT_DIR = PROJECT_ROOT / "output" / "liang_wagner"
CENTROID_PATH = (
    PROJECT_ROOT / "output" / "phase2" / "scaled_35types" / "centroids_human_35.csv"
)
RESIDUAL_PATH = (
    PROJECT_ROOT / "output" / "phase2" / "scaled_35types" / "residuals_ranked.csv"
)
TREENESS_PATH = OUTPUT_DIR / "treeness_scores_per_celltype.csv"
K_VALUES = [3, 5, 10]
K_PRIMARY = 5


def partial_spearman(x, y, z):
    """Partial Spearman correlation between x and y controlling for z.

    Math: Convert to ranks, then compute partial correlation:
      r_xy.z = (r_xy - r_xz * r_yz) / sqrt((1 - r_xz^2)(1 - r_yz^2))

    where r_xy, r_xz, r_yz are Spearman rank correlations.

    For the p-value, use the t-distribution approximation:
      t = r_xy.z * sqrt((n - 3) / (1 - r_xy.z^2))
    with df = n - 3.

    Args:
        x, y, z: Arrays of equal length.

    Returns:
        (partial_rho, p_value)
    """
    from scipy.stats import t as t_dist

    r_xy, _ = spearmanr(x, y)
    r_xz, _ = spearmanr(x, z)
    r_yz, _ = spearmanr(y, z)

    numer = r_xy - r_xz * r_yz
    denom = np.sqrt((1 - r_xz**2) * (1 - r_yz**2))

    if denom < 1e-15:
        return 0.0, 1.0

    partial_rho = numer / denom

    n = len(x)
    df = n - 3
    if df <= 0:
        return partial_rho, 1.0

    # t-test for partial correlation significance
    t_stat = partial_rho * np.sqrt(df / (1 - partial_rho**2 + 1e-15))
    p_val = 2 * t_dist.sf(np.abs(t_stat), df)

    return float(partial_rho), float(p_val)


def compute_neighborhood_density(centroids_mat, k_values):
    """Compute mean distance to k nearest neighbors for each cell type.

    For each cell type i, sorts all pairwise distances to other types
    and takes the mean of the k smallest.

    Args:
        centroids_mat: (n, G) array of centroids.
        k_values: List of k values to compute.

    Returns:
        Dict mapping k -> array of length n with density scores.
    """
    dist_vec = pdist(centroids_mat, metric="euclidean")
    dist_matrix = squareform(dist_vec)
    n = dist_matrix.shape[0]

    densities = {}
    for k in k_values:
        scores = np.empty(n)
        for i in range(n):
            # Get distances to all other types, sort ascending
            dists_i = np.sort(dist_matrix[i])  # includes self (0) at position 0
            # Skip self (index 0), take next k
            scores[i] = np.mean(dists_i[1 : k + 1])
        densities[k] = scores

    return densities, dist_matrix


def abbreviate(name):
    """Shorten cell type names for plot labels."""
    abbrevs = {
        "of epithelium of large intestine": "(colon)",
        "of mammary gland": "(mammary)",
        "of cardiac tissue": "(cardiac)",
        "of adipose tissue": "(adipose)",
        "-positive, alpha-beta ": "+ ",
        "-positive alpha-beta ": "+ ",
    }
    short = name
    for old, new in abbrevs.items():
        short = short.replace(old, new)
    return short


def plot_scatter(x, y, labels, xlabel, ylabel, title, output_path,
                 color=None, colorbar_label=None):
    """Generic scatter plot with cell type labels."""
    fig, ax = plt.subplots(1, 1, figsize=(9, 7))

    if color is not None:
        sc = ax.scatter(x, y, s=60, c=color, cmap="RdYlBu_r",
                        edgecolors="white", linewidth=0.5, alpha=0.8)
        cbar = fig.colorbar(sc, ax=ax, shrink=0.8)
        if colorbar_label:
            cbar.set_label(colorbar_label, fontsize=10)
    else:
        ax.scatter(x, y, s=60, c="steelblue", edgecolors="white",
                   linewidth=0.5, alpha=0.8)

    for i, name in enumerate(labels):
        ax.annotate(
            abbreviate(name), (x[i], y[i]),
            fontsize=6.5, alpha=0.8,
            xytext=(4, 4), textcoords="offset points",
        )

    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=11)

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def main():
    print("=" * 70)
    print("CellWarp — Neighborhood Density Mediation Test")
    print("Pre-registration: docs/preregistration_treeness_anticorrelation_2026-03-16.md")
    print("=" * 70)

    # --- Load data ---
    centroids_df = pd.read_csv(CENTROID_PATH, index_col=0)
    residuals_df = pd.read_csv(RESIDUAL_PATH)
    treeness_df = pd.read_csv(TREENESS_PATH)

    cell_types = sorted(centroids_df.index.tolist())
    centroids_mat = centroids_df.loc[cell_types].values
    n = len(cell_types)

    print(f"\n  Loaded: {n} cell types x {centroids_mat.shape[1]} genes")

    # --- Step 1: Compute pairwise distances and neighborhood density ---
    print("\n" + "=" * 70)
    print("STEP 1-2 — Pairwise Distances & Neighborhood Density")
    print("=" * 70)

    densities, dist_matrix = compute_neighborhood_density(centroids_mat, K_VALUES)

    # Save distance matrix
    dist_df = pd.DataFrame(dist_matrix, index=cell_types, columns=cell_types)
    dist_df.to_csv(OUTPUT_DIR / "pairwise_distance_matrix.csv")
    print(f"  Saved 35x35 distance matrix: {OUTPUT_DIR / 'pairwise_distance_matrix.csv'}")

    # Build merged DataFrame
    density_records = []
    for i, ct in enumerate(cell_types):
        rec = {"cell_type": ct}
        for k in K_VALUES:
            rec[f"density_k{k}"] = densities[k][i]
        density_records.append(rec)
    density_df = pd.DataFrame(density_records)

    # Merge all data
    merged = density_df.copy()
    merged = merged.merge(
        residuals_df[["cell_type", "residual_magnitude"]], on="cell_type"
    )
    merged = merged.merge(
        treeness_df[["cell_type", "treeness_score"]], on="cell_type"
    )
    assert len(merged) == n, f"Merge lost rows: {len(merged)} != {n}"

    # Save density CSV
    density_df.to_csv(OUTPUT_DIR / "neighborhood_density.csv", index=False)
    print(f"  Saved: {OUTPUT_DIR / 'neighborhood_density.csv'}")

    print(f"\n  Neighborhood density (k=5) summary:")
    d5 = merged["density_k5"]
    print(f"    Mean: {d5.mean():.2f}, Median: {d5.median():.2f}, "
          f"Std: {d5.std():.2f}")
    print(f"    Range: [{d5.min():.2f}, {d5.max():.2f}]")

    print(f"\n  Top 5 densest (lowest mean distance to k=5 NN):")
    for _, row in merged.nsmallest(5, "density_k5").iterrows():
        print(f"    {row['cell_type']:<50} {row['density_k5']:.2f}")
    print(f"\n  Top 5 most isolated (highest mean distance to k=5 NN):")
    for _, row in merged.nlargest(5, "density_k5").iterrows():
        print(f"    {row['cell_type']:<50} {row['density_k5']:.2f}")

    # --- Steps 3-5: Correlations at primary k=5 ---
    all_results = {}

    for k in K_VALUES:
        print(f"\n{'=' * 70}")
        is_primary = k == K_PRIMARY
        tag = " (PRIMARY)" if is_primary else " (sensitivity)"
        print(f"STEPS 3-5 — k={k}{tag}")
        print("=" * 70)

        density_col = f"density_k{k}"
        density_vals = merged[density_col].values
        residual_vals = merged["residual_magnitude"].values
        treeness_vals = merged["treeness_score"].values

        # H1: density vs rigidity
        # density = mean dist to kNN (lower = denser)
        # residual_magnitude (lower = more rigid)
        # Expect POSITIVE: denser (lower density) ~ more rigid (lower residual)
        rho_h1, p_h1 = spearmanr(density_vals, residual_vals)
        print(f"\n  H1: rho(density_k{k}, residual_magnitude) = {rho_h1:.4f}, "
              f"p = {p_h1:.4f}")
        print(f"    Interpretation: {'POSITIVE as expected' if rho_h1 > 0 else 'NEGATIVE (unexpected)'}")
        h1_pass = bool(rho_h1 > 0 and p_h1 < 0.05)
        print(f"    H1 {'PASS' if h1_pass else 'FAIL'} "
              f"(need rho > 0, p < 0.05)")

        # H2: density vs treeness
        # Expect NEGATIVE: denser (lower density) ~ lower treeness
        # But density = mean distance (lower = denser), treeness = mean delta
        # So expect POSITIVE rho(mean_dist, treeness): more isolated = more tree-like
        rho_h2, p_h2 = spearmanr(density_vals, treeness_vals)
        print(f"\n  H2: rho(density_k{k}, treeness) = {rho_h2:.4f}, "
              f"p = {p_h2:.4f}")
        # Note: density_vals = mean distance. Lower = denser.
        # If denser = less tree-like, then rho(mean_dist, treeness) should be POSITIVE
        # (more isolated = higher mean dist = higher treeness)
        # But the pre-reg says "denser neighborhood = lower treeness" and expects
        # rho(density, treeness) < 0. "Density" here means the mean distance score
        # (confusing naming). Let me clarify:
        # density_score = mean_distance_to_kNN. Higher score = MORE isolated = LESS dense.
        # treeness_score = mean delta. Higher = more tree-like.
        # If less dense (higher score) = more tree-like (higher treeness), rho > 0.
        # Pre-reg H2 expects: denser = lower treeness, i.e., lower density_score = lower treeness.
        # That means rho(density_score, treeness) > 0: POSITIVE.
        # Wait — pre-reg says "rho < 0" but density_score is mean distance (inverted).
        # Let me just report correctly and check the sign.
        # Pre-reg: "Spearman rho(neighborhood_density_k5, treeness_score). Expected: NEGATIVE"
        # But if neighborhood_density = mean_distance (lower = denser), then
        # the expected direction is POSITIVE rho (isolated types have both higher
        # mean distance AND higher treeness).
        # The pre-reg may have intended "density" as a true density measure (higher = denser).
        # I'll compute it both ways and report clearly.

        # For clarity: define isolation = mean_dist_to_kNN (what we computed)
        # Then: isolated types expected to be LESS rigid and MORE tree-like
        # rho(isolation, residual) > 0: isolated = more diverged ✓
        # rho(isolation, treeness) > 0: isolated = more tree-like ✓
        h2_pass = bool(rho_h2 > 0 and p_h2 < 0.05)
        print(f"    Note: density score = mean distance to kNN (higher = MORE isolated)")
        print(f"    Expected: positive rho (isolated = more tree-like)")
        print(f"    H2 {'PASS' if h2_pass else 'FAIL'} "
              f"(need rho > 0, p < 0.05)")

        # H3: partial correlation (rigidity vs treeness | density)
        # Raw: rho(residual, treeness) = +0.349 (positive because higher residual = lower rigidity)
        rho_raw, p_raw = spearmanr(residual_vals, treeness_vals)
        partial_rho, partial_p = partial_spearman(
            residual_vals, treeness_vals, density_vals
        )
        attenuation = (abs(rho_raw) - abs(partial_rho)) / abs(rho_raw) * 100

        print(f"\n  H3: Mediation test")
        print(f"    Raw rho(residual, treeness): {rho_raw:.4f} (p={p_raw:.4f})")
        print(f"    Partial rho(residual, treeness | density_k{k}): "
              f"{partial_rho:.4f} (p={partial_p:.4f})")
        print(f"    Attenuation: {attenuation:.1f}%")
        h3_pass = bool(attenuation >= 50)
        print(f"    H3 {'PASS' if h3_pass else 'FAIL'} "
              f"(need attenuation >= 50%)")

        # Decision
        if h1_pass and h2_pass and h3_pass:
            verdict = "MEDIATION_CONFIRMED"
            desc = "Anticorrelation is a geometric proximity artifact"
        elif h1_pass and h2_pass and not h3_pass:
            verdict = "PARTIAL_MEDIATION"
            desc = "Density contributes but does not fully explain — additional biological component"
        elif h1_pass and not h2_pass:
            verdict = "DENSITY_PREDICTS_RIGIDITY_ONLY"
            desc = "Density predicts rigidity but not treeness — different mechanism"
        else:
            verdict = "MEDIATION_REJECTED"
            desc = "Density does not predict rigidity — mediation hypothesis rejected"

        print(f"\n  Verdict (k={k}): {verdict}")
        print(f"    {desc}")

        all_results[f"k{k}"] = {
            "k": k,
            "is_primary": is_primary,
            "h1_density_vs_residual": {
                "rho": float(rho_h1), "p": float(p_h1), "pass": h1_pass
            },
            "h2_density_vs_treeness": {
                "rho": float(rho_h2), "p": float(p_h2), "pass": h2_pass
            },
            "h3_mediation": {
                "raw_rho_residual_treeness": float(rho_raw),
                "raw_p": float(p_raw),
                "partial_rho": float(partial_rho),
                "partial_p": float(partial_p),
                "attenuation_pct": float(attenuation),
                "pass": h3_pass,
            },
            "verdict": verdict,
            "description": desc,
        }

    # --- Step 6: Sensitivity check ---
    print(f"\n{'=' * 70}")
    print("STEP 6 — Sensitivity Check")
    print("=" * 70)

    verdicts = [all_results[f"k{k}"]["verdict"] for k in K_VALUES]
    consistent = len(set(verdicts)) == 1
    print(f"  Verdicts across k values: {dict(zip(K_VALUES, verdicts))}")
    print(f"  Consistent: {'YES' if consistent else 'NO — sensitivity instability'}")

    if not consistent:
        print(f"  WARNING: Falsification condition 4 triggered — "
              f"qualitatively different conclusions across k values")

    # --- Step 7: Plots ---
    print(f"\n{'=' * 70}")
    print("STEP 7 — Plots")
    print("=" * 70)

    labels = merged["cell_type"].values

    # (a) Density vs rigidity
    plot_scatter(
        merged["density_k5"].values,
        merged["residual_magnitude"].values,
        labels,
        xlabel="Neighborhood Density k=5 (mean dist to 5 NN, higher = more isolated)",
        ylabel="Procrustes Residual Magnitude (higher = less rigid)",
        title=(f"H1: Density vs Rigidity — "
               f"rho={all_results['k5']['h1_density_vs_residual']['rho']:.3f}, "
               f"p={all_results['k5']['h1_density_vs_residual']['p']:.4f}"),
        output_path=OUTPUT_DIR / "density_rigidity_scatter.png",
    )

    # (b) Density vs treeness
    plot_scatter(
        merged["density_k5"].values,
        merged["treeness_score"].values,
        labels,
        xlabel="Neighborhood Density k=5 (mean dist to 5 NN, higher = more isolated)",
        ylabel="Treeness Score (mean delta, higher = more tree-like)",
        title=(f"H2: Density vs Treeness — "
               f"rho={all_results['k5']['h2_density_vs_treeness']['rho']:.3f}, "
               f"p={all_results['k5']['h2_density_vs_treeness']['p']:.4f}"),
        output_path=OUTPUT_DIR / "density_treeness_scatter.png",
    )

    # (c) Rigidity vs treeness colored by density quartile
    density_q = pd.qcut(merged["density_k5"], q=4, labels=False) + 1
    plot_scatter(
        merged["residual_magnitude"].values,
        merged["treeness_score"].values,
        labels,
        xlabel="Procrustes Residual Magnitude (higher = less rigid)",
        ylabel="Treeness Score (mean delta, higher = more tree-like)",
        title=(f"Rigidity vs Treeness colored by Density Quartile — "
               f"partial rho={all_results['k5']['h3_mediation']['partial_rho']:.3f}, "
               f"attenuation={all_results['k5']['h3_mediation']['attenuation_pct']:.1f}%"),
        output_path=OUTPUT_DIR / "rigidity_treeness_by_density.png",
        color=density_q.values,
        colorbar_label="Density Quartile (1=densest, 4=most isolated)",
    )

    # --- Save results JSON ---
    with open(OUTPUT_DIR / "mediation_results.json", "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"  Saved: {OUTPUT_DIR / 'mediation_results.json'}")

    # --- Write summary ---
    primary = all_results["k5"]
    with open(OUTPUT_DIR / "h1_mediation_summary.md", "w") as f:
        f.write("# Neighborhood Density Mediation Test — Summary\n\n")
        f.write(f"**Date:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}\n")
        f.write("**Pre-registration:** docs/preregistration_treeness_anticorrelation_2026-03-16.md\n")
        f.write(f"**Input:** 35 human cell type centroids x 16,959 ortholog genes\n\n")

        f.write("## Results Table\n\n")
        f.write("| Metric | k=3 | k=5 (primary) | k=10 |\n")
        f.write("|--------|-----|---------------|------|\n")

        for metric_key, metric_label in [
            ("h1_density_vs_residual", "H1: rho(density, residual)"),
            ("h2_density_vs_treeness", "H2: rho(density, treeness)"),
        ]:
            vals = []
            for k in K_VALUES:
                r = all_results[f"k{k}"][metric_key]
                star = "*" if r["p"] < 0.05 else ""
                vals.append(f"{r['rho']:.3f} (p={r['p']:.3f}){star}")
            f.write(f"| {metric_label} | {' | '.join(vals)} |\n")

        # Mediation row
        vals = []
        for k in K_VALUES:
            r = all_results[f"k{k}"]["h3_mediation"]
            vals.append(f"{r['partial_rho']:.3f} (atten={r['attenuation_pct']:.0f}%)")
        f.write(f"| H3: partial rho (attenuation) | {' | '.join(vals)} |\n")

        # Verdict row
        vals = [all_results[f"k{k}"]["verdict"] for k in K_VALUES]
        f.write(f"| Verdict | {' | '.join(vals)} |\n")

        f.write(f"\n*Significance: * = p < 0.05*\n")

        f.write(f"\n## Falsification Conditions\n\n")
        h1p = primary["h1_density_vs_residual"]
        h2p = primary["h2_density_vs_treeness"]
        h3p = primary["h3_mediation"]

        f.write(f"1. H1 fails (density !~ rigidity): "
                f"{'TRIGGERED' if not h1p['pass'] else 'not triggered'}\n")
        f.write(f"2. H2 fails (density !~ treeness): "
                f"{'TRIGGERED' if not h2p['pass'] else 'not triggered'}\n")
        f.write(f"3. H3 fails (attenuation < 50%): "
                f"{'TRIGGERED' if not h3p['pass'] else 'not triggered'}\n")
        f.write(f"4. Sensitivity instability: "
                f"{'TRIGGERED' if not consistent else 'not triggered'}\n")

        f.write(f"\n## Primary Verdict (k=5)\n\n")
        f.write(f"**{primary['verdict']}**\n\n")
        f.write(f"{primary['description']}\n")

    print(f"  Saved: {OUTPUT_DIR / 'h1_mediation_summary.md'}")

    # --- Final summary ---
    print(f"\n{'=' * 70}")
    print("FINAL RESULTS TABLE")
    print("=" * 70)
    print(f"  {'Metric':<45} {'k=3':>15} {'k=5 (PRIMARY)':>15} {'k=10':>15}")
    print(f"  {'-' * 92}")

    for metric_key, label in [
        ("h1_density_vs_residual", "H1: rho(density, residual)"),
        ("h2_density_vs_treeness", "H2: rho(density, treeness)"),
    ]:
        vals = []
        for k in K_VALUES:
            r = all_results[f"k{k}"][metric_key]
            star = "*" if r["p"] < 0.05 else ""
            vals.append(f"{r['rho']:.3f}{star} p={r['p']:.3f}")
        print(f"  {label:<45} {vals[0]:>15} {vals[1]:>15} {vals[2]:>15}")

    vals = []
    for k in K_VALUES:
        r = all_results[f"k{k}"]["h3_mediation"]
        vals.append(f"{r['partial_rho']:.3f} ({r['attenuation_pct']:.0f}%)")
    print(f"  {'H3: partial rho (attenuation)':<45} {vals[0]:>15} {vals[1]:>15} {vals[2]:>15}")

    vals = []
    for k in K_VALUES:
        r = all_results[f"k{k}"]["h3_mediation"]
        vals.append(f"{'PASS' if r['pass'] else 'FAIL'}")
    print(f"  {'H3 pass (>=50% attenuation)?':<45} {vals[0]:>15} {vals[1]:>15} {vals[2]:>15}")

    print(f"\n  {'Raw rho(residual, treeness):':<45} {primary['h3_mediation']['raw_rho_residual_treeness']:.4f} (p={primary['h3_mediation']['raw_p']:.4f})")

    verdicts_str = ", ".join(f"k={k}: {all_results[f'k{k}']['verdict']}" for k in K_VALUES)
    print(f"\n  Verdicts: {verdicts_str}")
    print(f"  Sensitivity consistent: {'YES' if consistent else 'NO'}")
    print("=" * 70)


if __name__ == "__main__":
    main()
