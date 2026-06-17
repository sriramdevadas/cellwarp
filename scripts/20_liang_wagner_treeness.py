#!/usr/bin/env python3
"""
CellWarp Script 20 — Liang-Wagner Treeness Analysis

Replicates the treeness test from Liang & Wagner (2015, Nature Communications)
on CellWarp single-cell centroids and correlates per-cell-type treeness with
Procrustes rigidity.

Biology: Liang & Wagner showed that cell type transcriptomes have tree-like
geometric structure — normal cells fit a phylogenetic tree topology, cancer
cells do not. CellWarp measures Procrustes rigidity — how geometrically
conserved cell type positions are across species. These are independent
geometric frameworks applied to the same question: is cell identity
geometrically structured? Convergence would be extraordinary validation
across 11 years, two technologies (bulk vs. single-cell), and two
geometric frameworks.

Math: Uses the delta statistic (four-point condition) to test tree-likeness
of each tetrad of cell types, Holland et al. (2002) analytic p-values, and
Storey's (2002) pi0 estimation for aggregate signal. Per-cell-type treeness
is mean delta across all tetrads containing that cell type.

Steps:
  1. Overall treeness test on 35 human cell type centroids (16,959 genes)
  2. Per-cell-type treeness scores
  3. Spearman correlation with Procrustes rigidity
  4. Sensitivity: repeat in PCA-reduced space

Input:  output/phase2/scaled_35types/centroids_human_35.csv
        output/phase2/scaled_35types/residuals_ranked.csv
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
from scipy.stats import spearmanr
from sklearn.decomposition import PCA

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from cellwarp.treeness import (
    compute_all_deltas,
    estimate_pi0_storey,
    holland_pvalues,
    per_celltype_treeness,
)

# --- Config ---
OUTPUT_DIR = PROJECT_ROOT / "output" / "liang_wagner"
CENTROID_PATH = (
    PROJECT_ROOT / "output" / "phase2" / "scaled_35types" / "centroids_human_35.csv"
)
RESIDUAL_PATH = (
    PROJECT_ROOT / "output" / "phase2" / "scaled_35types" / "residuals_ranked.csv"
)
PCA_VARIANCE = 0.95
RANDOM_SEED = 42


def load_data():
    """Load human centroids and Procrustes residuals."""
    centroids_df = pd.read_csv(CENTROID_PATH, index_col=0)
    residuals_df = pd.read_csv(RESIDUAL_PATH)

    print(
        f"  Loaded centroids: {centroids_df.shape[0]} cell types "
        f"x {centroids_df.shape[1]} genes"
    )
    print(f"  Loaded residuals: {residuals_df.shape[0]} cell types")

    # Verify cell type names match between sources
    centroid_types = set(centroids_df.index)
    residual_types = set(residuals_df["cell_type"])
    if centroid_types != residual_types:
        missing_in_centroids = residual_types - centroid_types
        missing_in_residuals = centroid_types - residual_types
        if missing_in_centroids:
            print(f"  WARNING: In residuals but not centroids: {missing_in_centroids}")
        if missing_in_residuals:
            print(f"  WARNING: In centroids but not residuals: {missing_in_residuals}")
    else:
        print(f"  Cell type names match between both sources (n={len(centroid_types)})")

    return centroids_df, residuals_df


def step1_overall_treeness(centroids_mat, cell_types):
    """Step 1: Overall treeness test on all 35 cell types.

    Enumerates all C(35, 4) = 52,360 tetrads, computes delta for each,
    and estimates overall tree structure signal via Storey's pi0.
    """
    print("\n" + "=" * 70)
    print("STEP 1 — Overall Treeness Test")
    print("=" * 70)

    deltas, tetrads = compute_all_deltas(centroids_mat)
    pvalues = holland_pvalues(deltas)
    pi0 = estimate_pi0_storey(pvalues)

    n_sig_005 = int(np.sum(pvalues < 0.05))
    n_sig_001 = int(np.sum(pvalues < 0.01))

    print(f"\n  Results:")
    print(f"    Tetrads analyzed: {len(tetrads):,}")
    print(
        f"    delta distribution: mean={np.mean(deltas):.4f}, "
        f"median={np.median(deltas):.4f}, std={np.std(deltas):.4f}"
    )
    print(f"    delta range: [{np.min(deltas):.4f}, {np.max(deltas):.4f}]")
    print(
        f"    delta quartiles: Q1={np.percentile(deltas, 25):.4f}, "
        f"Q3={np.percentile(deltas, 75):.4f}"
    )
    print(
        f"    Tetrads with p < 0.05: {n_sig_005:,} / {len(tetrads):,} "
        f"({100 * n_sig_005 / len(tetrads):.1f}%)"
    )
    print(
        f"    Tetrads with p < 0.01: {n_sig_001:,} / {len(tetrads):,} "
        f"({100 * n_sig_001 / len(tetrads):.1f}%)"
    )
    print(f"    Storey pi0 estimate: {pi0:.4f}")
    print(f"    Fraction with tree structure (1 - pi0): {1 - pi0:.4f}")

    return deltas, tetrads, pvalues, pi0


def step2_per_celltype_treeness(deltas, tetrads, cell_types):
    """Step 2: Per-cell-type treeness scores.

    For each cell type i, computes mean delta across all tetrads that include
    cell type i. Higher mean delta = more tree-like when this type is included.
    """
    print("\n" + "=" * 70)
    print("STEP 2 — Per-Cell-Type Treeness Scores")
    print("=" * 70)

    scores = per_celltype_treeness(deltas, tetrads, len(cell_types))

    treeness_df = pd.DataFrame(
        {"cell_type": cell_types, "treeness_score": scores}
    )
    treeness_df = treeness_df.sort_values(
        "treeness_score", ascending=False
    ).reset_index(drop=True)
    treeness_df["rank"] = treeness_df.index + 1

    print(f"\n  Per-cell-type treeness (mean delta, sorted):")
    print(f"  {'Rank':<6} {'Cell Type':<50} {'Treeness':>10}")
    print(f"  {'-' * 68}")
    for _, row in treeness_df.iterrows():
        print(
            f"  {int(row['rank']):<6} {row['cell_type']:<50} "
            f"{row['treeness_score']:>10.4f}"
        )

    return treeness_df, scores


def step3_correlation(treeness_df, residuals_df):
    """Step 3: Spearman correlation between treeness and Procrustes rigidity.

    Rigidity = inverse of Procrustes residual magnitude (smaller residual =
    more geometrically conserved across species = more rigid). We compute
    rho(treeness, rigidity) = -rho(treeness, residual_magnitude).
    """
    print("\n" + "=" * 70)
    print("STEP 3 — Spearman Correlation: Treeness vs. Procrustes Rigidity")
    print("=" * 70)

    merged = treeness_df.merge(
        residuals_df[["cell_type", "residual_magnitude"]],
        on="cell_type",
        how="inner",
    )
    n_matched = len(merged)

    rho_residual, p_residual = spearmanr(
        merged["treeness_score"], merged["residual_magnitude"]
    )
    # Rigidity = -residual_magnitude (higher = more rigid)
    rho_rigidity = -rho_residual
    p_rigidity = p_residual

    print(f"\n  Matched cell types: {n_matched}")
    print(
        f"  Spearman rho(treeness, residual_magnitude): {rho_residual:.4f}, "
        f"p={p_residual:.4f}"
    )
    print(
        f"  Spearman rho(treeness, rigidity): {rho_rigidity:.4f}, "
        f"p={p_rigidity:.4f}"
    )
    print(
        f"    (rigidity = -residual_magnitude; "
        f"more rigid = smaller Procrustes residual)"
    )

    # Decision table interpretation
    print(f"\n  Decision table:")
    if rho_rigidity >= 0.6 and p_rigidity < 0.05:
        outcome = "CONVERGENT_VALIDATION"
        interp = "Two independent frameworks agree on rigidity ranking"
        paper_use = "Include as validation result, highlight in discussion"
        print(f"    rho >= 0.6, p < 0.05 -> CONVERGENT VALIDATION")
    elif 0.3 <= rho_rigidity < 0.6:
        outcome = "PARTIAL_OVERLAP"
        interp = "Frameworks capture related but distinct properties"
        paper_use = "Note as consistent trend, do not overclaim"
        print(f"    rho = 0.3-0.59 -> PARTIAL OVERLAP")
    elif rho_rigidity < 0 and p_rigidity < 0.05:
        outcome = "ANTICORRELATED"
        interp = "Most surprising — rigid types have less tree structure"
        paper_use = "Flag immediately for advisor; potentially most interesting"
        print(f"    rho significantly negative -> ANTICORRELATED (flag for advisor)")
    else:
        outcome = "DISTINCT_PROPERTIES"
        interp = "Measuring distinct geometric properties — both real, different dimensions"
        paper_use = (
            'One sentence: "Procrustes rigidity is geometrically '
            'distinct from within-species treeness"'
        )
        print(f"    rho < 0.3, NS -> DISTINCT PROPERTIES")

    print(f"    Interpretation: {interp}")
    print(f"    Paper 1 use: {paper_use}")

    return merged, rho_rigidity, p_rigidity, n_matched, outcome, interp


def step4_pca_sensitivity(centroids_df, cell_types, residuals_df, treeness_fullspace):
    """Step 4: Sensitivity check in PCA-reduced space.

    Liang & Wagner used full expression space; our Procrustes analysis uses
    PCA-reduced space. Repeating treeness in PCA space tests whether
    dimensionality reduction substantially changes the treeness signal
    or the correlation with rigidity.
    """
    print("\n" + "=" * 70)
    print("STEP 4 — Sensitivity: PCA-Reduced Space")
    print("=" * 70)

    centroids_mat = centroids_df.loc[cell_types].values

    pca = PCA(n_components=PCA_VARIANCE, svd_solver="full", random_state=RANDOM_SEED)
    centroids_pca = pca.fit_transform(centroids_mat)
    n_components = pca.n_components_
    cumvar = np.sum(pca.explained_variance_ratio_) * 100

    print(f"\n  PCA: {n_components} components, {cumvar:.1f}% variance")

    deltas_pca, tetrads_pca = compute_all_deltas(centroids_pca)
    pvalues_pca = holland_pvalues(deltas_pca)
    pi0_pca = estimate_pi0_storey(pvalues_pca)

    scores_pca = per_celltype_treeness(deltas_pca, tetrads_pca, len(cell_types))

    # Correlation with rigidity in PCA space
    treeness_pca_df = pd.DataFrame(
        {"cell_type": cell_types, "treeness_pca": scores_pca}
    )
    merged_pca = treeness_pca_df.merge(
        residuals_df[["cell_type", "residual_magnitude"]],
        on="cell_type",
        how="inner",
    )
    rho_pca, p_pca = spearmanr(
        merged_pca["treeness_pca"], merged_pca["residual_magnitude"]
    )
    rho_pca_rigidity = -rho_pca

    # Consistency between full-space and PCA-space treeness
    treeness_full_df = pd.DataFrame(
        {"cell_type": cell_types, "treeness_full": treeness_fullspace}
    )
    merged_both = treeness_pca_df.merge(treeness_full_df, on="cell_type")
    rho_consistency, p_consistency = spearmanr(
        merged_both["treeness_pca"], merged_both["treeness_full"]
    )

    print(f"\n  PCA-space results:")
    print(
        f"    delta distribution: mean={np.mean(deltas_pca):.4f}, "
        f"median={np.median(deltas_pca):.4f}"
    )
    print(f"    pi0 estimate: {pi0_pca:.4f}")
    print(
        f"    rho(treeness_pca, rigidity): {rho_pca_rigidity:.4f}, "
        f"p={p_pca:.4f}"
    )
    print(
        f"    Consistency: rho(full-space, PCA) = {rho_consistency:.4f}, "
        f"p={p_consistency:.6f}"
    )

    pca_results = {
        "n_components": int(n_components),
        "variance_explained_pct": float(cumvar),
        "delta_mean": float(np.mean(deltas_pca)),
        "delta_median": float(np.median(deltas_pca)),
        "delta_std": float(np.std(deltas_pca)),
        "pi0": float(pi0_pca),
        "fraction_treelike": float(1 - pi0_pca),
        "rho_rigidity": float(rho_pca_rigidity),
        "p_rigidity": float(p_pca),
        "rho_consistency_with_full": float(rho_consistency),
        "p_consistency": float(p_consistency),
    }

    return pca_results, treeness_pca_df


def plot_delta_distribution(deltas, pi0, output_path):
    """Plot distribution of delta across all tetrads with H0 density overlay.

    Shows: histogram of observed delta values, null density under H0 (no tree
    constraint), and summary statistics. Rightward skew indicates tree structure.
    """
    fig, ax = plt.subplots(1, 1, figsize=(8, 5))

    ax.hist(
        deltas, bins=50, density=True, alpha=0.7, color="steelblue",
        edgecolor="white", linewidth=0.5,
    )
    ax.axvline(
        np.median(deltas), color="red", linestyle="--", linewidth=1.5,
        label=f"Median delta = {np.median(deltas):.3f}",
    )
    ax.axvline(
        np.mean(deltas), color="orange", linestyle="--", linewidth=1.5,
        label=f"Mean delta = {np.mean(deltas):.3f}",
    )

    # H0 density: f(delta) = (6 / (pi * sqrt(3))) / (1 + ((2*delta - 1) / sqrt(3))^2)
    x = np.linspace(0.01, 0.99, 200)
    null_pdf = (6 / (np.pi * np.sqrt(3))) / (
        1 + ((2 * x - 1) / np.sqrt(3)) ** 2
    )
    ax.plot(x, null_pdf, "k-", linewidth=1.5, alpha=0.5, label="H0 density (no tree)")

    ax.set_xlabel("delta (tetrad treeness statistic)", fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    ax.set_title(
        f"Tetrad delta Distribution — {len(deltas):,} tetrads, "
        f"pi0 = {pi0:.3f}",
        fontsize=13,
    )
    ax.legend(fontsize=10)
    ax.set_xlim(0, 1)

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def plot_treeness_vs_rigidity(merged, rho, p_val, output_path):
    """Scatter plot of treeness vs. Procrustes residual magnitude.

    Shows: per-cell-type treeness (y) vs residual magnitude (x), labeled.
    X-axis is residual (right = more diverged = less rigid).
    Spearman rho with rigidity (sign-flipped) annotated in title.
    """
    fig, ax = plt.subplots(1, 1, figsize=(9, 7))

    ax.scatter(
        merged["residual_magnitude"],
        merged["treeness_score"],
        s=60, alpha=0.8, c="steelblue", edgecolors="white", linewidth=0.5,
    )

    # Abbreviate long cell type names for readability
    abbreviations = {
        "of epithelium of large intestine": "(colon)",
        "of mammary gland": "(mammary)",
        "of cardiac tissue": "(cardiac)",
        "of adipose tissue": "(adipose)",
        "-positive, alpha-beta ": "+ ",
        "-positive alpha-beta ": "+ ",
    }
    for _, row in merged.iterrows():
        name = row["cell_type"]
        short = name
        for old, new in abbreviations.items():
            short = short.replace(old, new)
        ax.annotate(
            short,
            (row["residual_magnitude"], row["treeness_score"]),
            fontsize=6.5, alpha=0.8,
            xytext=(4, 4), textcoords="offset points",
        )

    ax.set_xlabel(
        "Procrustes Residual Magnitude (right = less rigid)", fontsize=12
    )
    ax.set_ylabel(
        "Treeness Score (mean delta, up = more tree-like)", fontsize=12
    )
    ax.set_title(
        f"Treeness vs. Procrustes Rigidity — "
        f"rho(treeness, rigidity) = {rho:.3f}, p = {p_val:.4f}, n = {len(merged)}",
        fontsize=11,
    )

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def write_summary(results, output_path):
    """Write human-readable markdown summary of all results."""
    r = results

    with open(output_path, "w") as f:
        f.write("# Liang-Wagner Treeness Analysis — Summary\n\n")
        f.write(f"**Date:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}\n")
        f.write(
            "**Input:** 35 human cell type centroids x 16,959 ortholog genes\n"
        )
        f.write(
            "**Method:** Liang & Wagner (2015) delta statistic, "
            "Holland et al. (2002) analytic p-values, "
            "Storey (2002) pi0 estimation\n\n"
        )

        f.write("## Step 1 — Overall Treeness\n\n")
        f.write(f"- Tetrads analyzed: {r['n_tetrads']:,}\n")
        f.write(
            f"- delta distribution: mean={r['delta_mean']:.4f}, "
            f"median={r['delta_median']:.4f}, std={r['delta_std']:.4f}\n"
        )
        f.write(f"- delta range: [{r['delta_min']:.4f}, {r['delta_max']:.4f}]\n")
        f.write(
            f"- Tetrads with p < 0.05: {r['n_sig_005']:,} "
            f"({r['pct_sig_005']:.1f}%)\n"
        )
        f.write(
            f"- Tetrads with p < 0.01: {r['n_sig_001']:,} "
            f"({r['pct_sig_001']:.1f}%)\n"
        )
        f.write(f"- Storey pi0: {r['pi0']:.4f}\n")
        f.write(f"- Tree-like fraction (1 - pi0): {r['fraction_treelike']:.4f}\n\n")

        f.write("## Step 2 — Per-Cell-Type Treeness Ranking\n\n")
        f.write("| Rank | Cell Type | Treeness (mean delta) |\n")
        f.write("|------|-----------|----------------------|\n")
        for _, row in r["treeness_df"].iterrows():
            f.write(
                f"| {int(row['rank'])} | {row['cell_type']} "
                f"| {row['treeness_score']:.4f} |\n"
            )

        f.write(f"\n## Step 3 — Correlation with Procrustes Rigidity\n\n")
        f.write(
            f"- Spearman rho(treeness, rigidity): "
            f"**{r['rho_rigidity']:.4f}** (p={r['p_rigidity']:.4f})\n"
        )
        f.write(f"- n = {r['n_matched']}\n")
        f.write(f"- Outcome: **{r['outcome']}**\n")
        f.write(f"- Interpretation: {r['interpretation']}\n\n")

        f.write("## Step 4 — PCA Sensitivity\n\n")
        pca = r["pca_results"]
        f.write(
            f"- PCA components: {pca['n_components']} "
            f"({pca['variance_explained_pct']:.1f}% variance)\n"
        )
        f.write(f"- delta mean (PCA space): {pca['delta_mean']:.4f}\n")
        f.write(f"- pi0 (PCA space): {pca['pi0']:.4f}\n")
        f.write(
            f"- rho(treeness_pca, rigidity): {pca['rho_rigidity']:.4f} "
            f"(p={pca['p_rigidity']:.4f})\n"
        )
        f.write(
            f"- Full-space vs PCA consistency: "
            f"rho={pca['rho_consistency_with_full']:.4f} "
            f"(p={pca['p_consistency']:.6f})\n\n"
        )

        f.write("## Decision\n\n")
        if r["rho_rigidity"] >= 0.6 and r["p_rigidity"] < 0.05:
            f.write(
                "**Include in Paper 1** — convergent validation across "
                "independent frameworks.\n"
            )
        elif 0.3 <= r["rho_rigidity"] < 0.6:
            f.write(
                "**Note as consistent trend** — partial overlap, "
                "do not overclaim.\n"
            )
        elif r["rho_rigidity"] < 0 and r["p_rigidity"] < 0.05:
            f.write(
                "**Flag for advisor** — anticorrelated, potentially "
                "most interesting finding.\n"
            )
        else:
            f.write(
                "**One sentence in discussion** — Procrustes rigidity is "
                "geometrically distinct from within-species treeness.\n"
            )

    print(f"  Saved: {output_path}")


def main():
    """Run Liang-Wagner treeness analysis pipeline."""
    print("=" * 70)
    print("CellWarp — Liang-Wagner Treeness Analysis")
    print("Replicating Liang & Wagner (2015) on single-cell centroids")
    print("=" * 70)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load data
    centroids_df, residuals_df = load_data()
    cell_types = sorted(centroids_df.index.tolist())
    centroids_mat = centroids_df.loc[cell_types].values

    # Step 1: Overall treeness
    deltas, tetrads, pvalues, pi0 = step1_overall_treeness(
        centroids_mat, cell_types
    )

    # Step 2: Per-cell-type treeness
    treeness_df, treeness_scores = step2_per_celltype_treeness(
        deltas, tetrads, cell_types
    )

    # Step 3: Correlation with rigidity
    merged, rho_rigidity, p_rigidity, n_matched, outcome, interp = (
        step3_correlation(treeness_df, residuals_df)
    )

    # Step 4: PCA sensitivity
    pca_results, treeness_pca_df = step4_pca_sensitivity(
        centroids_df, cell_types, residuals_df, treeness_scores
    )

    # --- Outputs ---
    print("\n" + "=" * 70)
    print("OUTPUTS")
    print("=" * 70)

    # Save per-cell-type scores CSV
    treeness_out = treeness_df[["rank", "cell_type", "treeness_score"]].copy()
    treeness_out.to_csv(
        OUTPUT_DIR / "treeness_scores_per_celltype.csv", index=False
    )
    print(f"  Saved: {OUTPUT_DIR / 'treeness_scores_per_celltype.csv'}")

    # Plot delta distribution
    plot_delta_distribution(
        deltas, pi0, OUTPUT_DIR / "tetrad_delta_distribution.png"
    )

    # Plot treeness vs rigidity
    plot_treeness_vs_rigidity(
        merged, rho_rigidity, p_rigidity,
        OUTPUT_DIR / "treeness_vs_rigidity_scatter.png",
    )

    # Save correlation JSON
    n_sig_005 = int(np.sum(pvalues < 0.05))
    n_sig_001 = int(np.sum(pvalues < 0.01))

    correlation_json = {
        "step1_overall": {
            "n_cell_types": len(cell_types),
            "n_tetrads": len(tetrads),
            "delta_mean": float(np.mean(deltas)),
            "delta_median": float(np.median(deltas)),
            "delta_std": float(np.std(deltas)),
            "delta_min": float(np.min(deltas)),
            "delta_max": float(np.max(deltas)),
            "delta_q25": float(np.percentile(deltas, 25)),
            "delta_q75": float(np.percentile(deltas, 75)),
            "n_significant_005": n_sig_005,
            "n_significant_001": n_sig_001,
            "pct_significant_005": float(100 * n_sig_005 / len(tetrads)),
            "pct_significant_001": float(100 * n_sig_001 / len(tetrads)),
            "pi0_storey": float(pi0),
            "fraction_treelike": float(1 - pi0),
        },
        "step3_correlation": {
            "rho_treeness_rigidity": float(rho_rigidity),
            "p_value": float(p_rigidity),
            "n": n_matched,
            "outcome": outcome,
            "interpretation": interp,
        },
        "step4_pca_sensitivity": pca_results,
    }

    with open(OUTPUT_DIR / "treeness_rigidity_correlation.json", "w") as f:
        json.dump(correlation_json, f, indent=2)
    print(f"  Saved: {OUTPUT_DIR / 'treeness_rigidity_correlation.json'}")

    # Write summary
    summary_results = {
        "n_tetrads": len(tetrads),
        "delta_mean": float(np.mean(deltas)),
        "delta_median": float(np.median(deltas)),
        "delta_std": float(np.std(deltas)),
        "delta_min": float(np.min(deltas)),
        "delta_max": float(np.max(deltas)),
        "n_sig_005": n_sig_005,
        "pct_sig_005": 100 * n_sig_005 / len(tetrads),
        "n_sig_001": n_sig_001,
        "pct_sig_001": 100 * n_sig_001 / len(tetrads),
        "pi0": float(pi0),
        "fraction_treelike": float(1 - pi0),
        "treeness_df": treeness_df,
        "rho_rigidity": float(rho_rigidity),
        "p_rigidity": float(p_rigidity),
        "n_matched": n_matched,
        "outcome": outcome,
        "interpretation": interp,
        "pca_results": pca_results,
    }
    write_summary(summary_results, OUTPUT_DIR / "liang_wagner_summary.md")

    # Final summary
    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    print(
        f"  Tree structure signal: delta mean={np.mean(deltas):.4f}, "
        f"pi0={pi0:.4f}"
    )
    print(
        f"  Treeness-rigidity correlation: rho={rho_rigidity:.4f}, "
        f"p={p_rigidity:.4f}"
    )
    print(f"  Outcome: {outcome}")
    print(
        f"  PCA sensitivity: rho={pca_results['rho_rigidity']:.4f}, "
        f"consistency={pca_results['rho_consistency_with_full']:.4f}"
    )
    print("=" * 70)


if __name__ == "__main__":
    main()
