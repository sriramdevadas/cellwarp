#!/usr/bin/env python3
"""
Developmental constraint analysis: does evolutionary rigidity correlate with
developmental origin?

Biology:
    D'Arcy Thompson's framework predicts that developmental constraints shape
    the geometry of evolutionary transformations. Cell types sharing a germ layer
    or developmental lineage may experience similar selective pressures, leading
    to correlated patterns of evolutionary rigidity (low Procrustes residual) or
    flexibility (high residual). Progenitor cells, which retain developmental
    plasticity, may show different evolutionary dynamics than terminally
    differentiated cells.

Math:
    For each of the 35 cell types from the scaled Procrustes analysis, we assign
    developmental annotations (germ layer, lineage, progenitor status) based on
    established developmental biology. We then test for group differences in
    Procrustes residual magnitude using non-parametric tests (Kruskal-Wallis for
    multi-group comparisons, Mann-Whitney U for binary comparisons). We also
    estimate proliferative capacity via mean MKI67 expression from centroids.

Output: ./output/phase2/developmental_constraint/
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path


# ─── Developmental annotations ───────────────────────────────────────────────
# Curated from standard developmental biology references.
# Germ layer: ectoderm, mesoderm, endoderm
# Lineage: hematopoietic, mesenchymal, epithelial, endothelial
# Progenitor: True if stem/progenitor cell, False if terminally differentiated

ANNOTATIONS = {
    # ── Mesoderm — Hematopoietic lineage ──
    "hematopoietic stem cell": {
        "germ_layer": "mesoderm", "lineage": "hematopoietic", "progenitor": True},
    "hematopoietic precursor cell": {
        "germ_layer": "mesoderm", "lineage": "hematopoietic", "progenitor": True},
    "T cell": {
        "germ_layer": "mesoderm", "lineage": "hematopoietic", "progenitor": False},
    "CD4-positive, alpha-beta T cell": {
        "germ_layer": "mesoderm", "lineage": "hematopoietic", "progenitor": False},
    "CD8-positive, alpha-beta T cell": {
        "germ_layer": "mesoderm", "lineage": "hematopoietic", "progenitor": False},
    "B cell": {
        "germ_layer": "mesoderm", "lineage": "hematopoietic", "progenitor": False},
    "plasma cell": {
        "germ_layer": "mesoderm", "lineage": "hematopoietic", "progenitor": False},
    "natural killer cell": {
        "germ_layer": "mesoderm", "lineage": "hematopoietic", "progenitor": False},
    "mature NK T cell": {
        "germ_layer": "mesoderm", "lineage": "hematopoietic", "progenitor": False},
    "macrophage": {
        "germ_layer": "mesoderm", "lineage": "hematopoietic", "progenitor": False},
    "monocyte": {
        "germ_layer": "mesoderm", "lineage": "hematopoietic", "progenitor": False},
    "classical monocyte": {
        "germ_layer": "mesoderm", "lineage": "hematopoietic", "progenitor": False},
    "non-classical monocyte": {
        "germ_layer": "mesoderm", "lineage": "hematopoietic", "progenitor": False},
    "intermediate monocyte": {
        "germ_layer": "mesoderm", "lineage": "hematopoietic", "progenitor": False},
    "myeloid leukocyte": {
        "germ_layer": "mesoderm", "lineage": "hematopoietic", "progenitor": False},
    "myeloid dendritic cell": {
        "germ_layer": "mesoderm", "lineage": "hematopoietic", "progenitor": False},
    "neutrophil": {
        "germ_layer": "mesoderm", "lineage": "hematopoietic", "progenitor": False},
    "granulocyte": {
        "germ_layer": "mesoderm", "lineage": "hematopoietic", "progenitor": False},

    # ── Mesoderm — Mesenchymal lineage ──
    "mesenchymal stem cell": {
        "germ_layer": "mesoderm", "lineage": "mesenchymal", "progenitor": True},
    "mesenchymal stem cell of adipose tissue": {
        "germ_layer": "mesoderm", "lineage": "mesenchymal", "progenitor": True},
    "fibroblast": {
        "germ_layer": "mesoderm", "lineage": "mesenchymal", "progenitor": False},
    "fibroblast of cardiac tissue": {
        "germ_layer": "mesoderm", "lineage": "mesenchymal", "progenitor": False},
    "stromal cell": {
        "germ_layer": "mesoderm", "lineage": "mesenchymal", "progenitor": False},
    "adventitial cell": {
        "germ_layer": "mesoderm", "lineage": "mesenchymal", "progenitor": False},
    "smooth muscle cell": {
        "germ_layer": "mesoderm", "lineage": "mesenchymal", "progenitor": False},

    # ── Mesoderm — Endothelial lineage ──
    "endothelial cell": {
        "germ_layer": "mesoderm", "lineage": "endothelial", "progenitor": False},

    # ── Endoderm — Epithelial lineage ──
    "hepatocyte": {
        "germ_layer": "endoderm", "lineage": "epithelial", "progenitor": False},
    "pancreatic acinar cell": {
        "germ_layer": "endoderm", "lineage": "epithelial", "progenitor": False},
    "pancreatic ductal cell": {
        "germ_layer": "endoderm", "lineage": "epithelial", "progenitor": False},
    "large intestine goblet cell": {
        "germ_layer": "endoderm", "lineage": "epithelial", "progenitor": False},
    "enterocyte of epithelium of large intestine": {
        "germ_layer": "endoderm", "lineage": "epithelial", "progenitor": False},

    # ── Ectoderm — Epithelial lineage ──
    "basal cell": {
        "germ_layer": "ectoderm", "lineage": "epithelial", "progenitor": True},
    "luminal epithelial cell of mammary gland": {
        "germ_layer": "ectoderm", "lineage": "epithelial", "progenitor": False},
    "bladder urothelial cell": {
        "germ_layer": "endoderm", "lineage": "epithelial", "progenitor": False},

    # ── Mixed / broad category ──
    "epithelial cell": {
        "germ_layer": "mixed", "lineage": "epithelial", "progenitor": False},
}


def load_residuals(residuals_path):
    """Load the 35-type ranked residuals."""
    df = pd.read_csv(residuals_path)
    return df


def annotate_cell_types(df):
    """Attach developmental annotations to each cell type."""
    records = []
    missing = []
    for _, row in df.iterrows():
        ct = row["cell_type"]
        if ct in ANNOTATIONS:
            ann = ANNOTATIONS[ct]
            records.append({
                "cell_type": ct,
                "rank": row["rank"],
                "residual_magnitude": row["residual_magnitude"],
                "pct_of_ssr": row["pct_of_ssr"],
                "germ_layer": ann["germ_layer"],
                "lineage": ann["lineage"],
                "progenitor": ann["progenitor"],
            })
        else:
            missing.append(ct)

    if missing:
        print(f"WARNING: No annotations for: {missing}")

    return pd.DataFrame(records)


def estimate_proliferation(centroids_path, mki67_ensembl="ENSG00000148773"):
    """
    Estimate proliferative capacity from mean MKI67 expression in centroids.

    MKI67 is a canonical marker of cell proliferation — it is expressed
    exclusively in actively dividing cells (all phases except G0).
    Higher mean MKI67 in the centroid suggests a larger fraction of
    cycling cells in that population.
    """
    df = pd.read_csv(centroids_path, index_col=0)
    if mki67_ensembl in df.columns:
        return df[mki67_ensembl].to_dict()
    else:
        print(f"WARNING: MKI67 ({mki67_ensembl}) not found in centroids.")
        return None


def run_statistical_tests(df):
    """
    Run non-parametric tests for group differences in residual magnitude.

    Kruskal-Wallis: tests if residual distributions differ across ≥3 groups.
    Mann-Whitney U: tests if two groups differ (progenitor vs differentiated).
    """
    results = {}

    # --- Kruskal-Wallis across germ layers (exclude 'mixed' — only 1 observation) ---
    germ_groups = df[df["germ_layer"] != "mixed"].groupby("germ_layer")["residual_magnitude"]
    germ_data = [g.values for _, g in germ_groups]
    germ_names = [name for name, _ in germ_groups]
    if len(germ_data) >= 2:
        stat, p = stats.kruskal(*germ_data)
        results["germ_layer_kw"] = {
            "test": "Kruskal-Wallis",
            "groups": germ_names,
            "group_medians": {name: float(np.median(vals)) for name, vals in zip(germ_names, germ_data)},
            "group_sizes": {name: len(vals) for name, vals in zip(germ_names, germ_data)},
            "H_statistic": float(stat),
            "p_value": float(p),
            "significant_005": p < 0.05,
        }

    # --- Kruskal-Wallis across lineages ---
    lin_groups = df.groupby("lineage")["residual_magnitude"]
    lin_data = [g.values for _, g in lin_groups]
    lin_names = [name for name, _ in lin_groups]
    if len(lin_data) >= 2:
        stat, p = stats.kruskal(*lin_data)
        results["lineage_kw"] = {
            "test": "Kruskal-Wallis",
            "groups": lin_names,
            "group_medians": {name: float(np.median(vals)) for name, vals in zip(lin_names, lin_data)},
            "group_sizes": {name: len(vals) for name, vals in zip(lin_names, lin_data)},
            "H_statistic": float(stat),
            "p_value": float(p),
            "significant_005": p < 0.05,
        }

    # --- Mann-Whitney U: progenitor vs differentiated ---
    prog = df[df["progenitor"] == True]["residual_magnitude"].values
    diff = df[df["progenitor"] == False]["residual_magnitude"].values
    if len(prog) >= 1 and len(diff) >= 1:
        stat, p = stats.mannwhitneyu(prog, diff, alternative="two-sided")
        results["progenitor_mwu"] = {
            "test": "Mann-Whitney U",
            "progenitor_n": len(prog),
            "progenitor_median": float(np.median(prog)),
            "differentiated_n": len(diff),
            "differentiated_median": float(np.median(diff)),
            "U_statistic": float(stat),
            "p_value": float(p),
            "significant_005": p < 0.05,
        }

    # --- Pairwise Mann-Whitney between germ layers ---
    from itertools import combinations
    pairwise = []
    for (n1, d1), (n2, d2) in combinations(zip(germ_names, germ_data), 2):
        if len(d1) >= 2 and len(d2) >= 2:
            stat, p = stats.mannwhitneyu(d1, d2, alternative="two-sided")
            pairwise.append({
                "group1": n1, "group2": n2,
                "median1": float(np.median(d1)), "median2": float(np.median(d2)),
                "U": float(stat), "p": float(p),
            })
    results["germ_layer_pairwise"] = pairwise

    # --- Spearman correlation: residual vs MKI67 (if available) ---
    if "mki67_human" in df.columns:
        valid = df.dropna(subset=["mki67_human"])
        if len(valid) >= 5:
            rho, p = stats.spearmanr(valid["residual_magnitude"], valid["mki67_human"])
            results["mki67_correlation"] = {
                "test": "Spearman correlation",
                "rho": float(rho),
                "p_value": float(p),
                "n": len(valid),
                "significant_005": p < 0.05,
            }
    if "mki67_mean" in df.columns:
        valid = df.dropna(subset=["mki67_mean"])
        if len(valid) >= 5:
            rho, p = stats.spearmanr(valid["residual_magnitude"], valid["mki67_mean"])
            results["mki67_mean_correlation"] = {
                "test": "Spearman correlation (mean human+mouse)",
                "rho": float(rho),
                "p_value": float(p),
                "n": len(valid),
                "significant_005": p < 0.05,
            }

    return results


def make_plots(df, outdir):
    """Generate all visualizations."""
    sns.set_style("whitegrid")
    palette_germ = {"ectoderm": "#E74C3C", "mesoderm": "#3498DB",
                    "endoderm": "#2ECC71", "mixed": "#95A5A6"}
    palette_lin = {"hematopoietic": "#3498DB", "mesenchymal": "#E67E22",
                   "epithelial": "#2ECC71", "endothelial": "#9B59B6"}

    # ── 1. Box plot by germ layer ──
    fig, ax = plt.subplots(figsize=(8, 6))
    order = ["ectoderm", "mesoderm", "endoderm", "mixed"]
    present = [g for g in order if g in df["germ_layer"].values]
    sns.boxplot(data=df, x="germ_layer", y="residual_magnitude",
                order=present, palette=palette_germ, ax=ax, width=0.5)
    sns.stripplot(data=df, x="germ_layer", y="residual_magnitude",
                  order=present, color="black", alpha=0.6, size=5, ax=ax)
    ax.set_xlabel("Germ Layer", fontsize=12)
    ax.set_ylabel("Procrustes Residual Magnitude", fontsize=12)
    ax.set_title("Evolutionary Divergence by Germ Layer of Origin", fontsize=14)
    # Add sample sizes
    for i, g in enumerate(present):
        n = (df["germ_layer"] == g).sum()
        ax.text(i, ax.get_ylim()[0] - 0.5, f"n={n}", ha="center", fontsize=10, color="gray")
    fig.tight_layout()
    fig.savefig(outdir / "boxplot_germ_layer.png", dpi=150)
    plt.close(fig)

    # ── 2. Box plot by lineage ──
    fig, ax = plt.subplots(figsize=(9, 6))
    lin_order = ["hematopoietic", "mesenchymal", "epithelial", "endothelial"]
    present_lin = [l for l in lin_order if l in df["lineage"].values]
    sns.boxplot(data=df, x="lineage", y="residual_magnitude",
                order=present_lin, palette=palette_lin, ax=ax, width=0.5)
    sns.stripplot(data=df, x="lineage", y="residual_magnitude",
                  order=present_lin, color="black", alpha=0.6, size=5, ax=ax)
    ax.set_xlabel("Developmental Lineage", fontsize=12)
    ax.set_ylabel("Procrustes Residual Magnitude", fontsize=12)
    ax.set_title("Evolutionary Divergence by Developmental Lineage", fontsize=14)
    for i, l in enumerate(present_lin):
        n = (df["lineage"] == l).sum()
        ax.text(i, ax.get_ylim()[0] - 0.5, f"n={n}", ha="center", fontsize=10, color="gray")
    fig.tight_layout()
    fig.savefig(outdir / "boxplot_lineage.png", dpi=150)
    plt.close(fig)

    # ── 3. Box plot: progenitor vs differentiated ──
    fig, ax = plt.subplots(figsize=(6, 6))
    df_prog = df.copy()
    df_prog["status"] = df_prog["progenitor"].map({True: "Progenitor/Stem", False: "Differentiated"})
    sns.boxplot(data=df_prog, x="status", y="residual_magnitude",
                order=["Progenitor/Stem", "Differentiated"],
                palette={"Progenitor/Stem": "#E74C3C", "Differentiated": "#3498DB"},
                ax=ax, width=0.4)
    sns.stripplot(data=df_prog, x="status", y="residual_magnitude",
                  order=["Progenitor/Stem", "Differentiated"],
                  color="black", alpha=0.6, size=5, ax=ax)
    ax.set_xlabel("Differentiation Status", fontsize=12)
    ax.set_ylabel("Procrustes Residual Magnitude", fontsize=12)
    ax.set_title("Evolutionary Divergence: Progenitors vs Differentiated Cells", fontsize=14)
    for i, s in enumerate(["Progenitor/Stem", "Differentiated"]):
        n = (df_prog["status"] == s).sum()
        ax.text(i, ax.get_ylim()[0] - 0.5, f"n={n}", ha="center", fontsize=10, color="gray")
    fig.tight_layout()
    fig.savefig(outdir / "boxplot_progenitor.png", dpi=150)
    plt.close(fig)

    # ── 4. Scatter: residual vs MKI67 ──
    if "mki67_mean" in df.columns:
        fig, ax = plt.subplots(figsize=(9, 7))
        colors = df["germ_layer"].map(palette_germ)
        ax.scatter(df["mki67_mean"], df["residual_magnitude"],
                   c=colors, s=60, edgecolors="black", linewidth=0.5, zorder=3)
        # Label points
        for _, row in df.iterrows():
            label = row["cell_type"]
            if len(label) > 25:
                label = label[:22] + "..."
            ax.annotate(label, (row["mki67_mean"], row["residual_magnitude"]),
                        fontsize=6, alpha=0.7,
                        xytext=(4, 2), textcoords="offset points")
        ax.set_xlabel("Mean MKI67 Expression (proxy for proliferation)", fontsize=12)
        ax.set_ylabel("Procrustes Residual Magnitude", fontsize=12)
        ax.set_title("Evolutionary Divergence vs Proliferative Capacity", fontsize=14)
        # Legend for germ layers
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor=c, label=g) for g, c in palette_germ.items()
                           if g in df["germ_layer"].values]
        ax.legend(handles=legend_elements, title="Germ Layer", loc="upper right")
        # Add Spearman rho on plot
        valid = df.dropna(subset=["mki67_mean"])
        if len(valid) >= 5:
            rho, p = stats.spearmanr(valid["residual_magnitude"], valid["mki67_mean"])
            ax.text(0.02, 0.98, f"Spearman ρ = {rho:.3f}, p = {p:.4f}",
                    transform=ax.transAxes, fontsize=10, va="top",
                    bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))
        fig.tight_layout()
        fig.savefig(outdir / "scatter_mki67_vs_residual.png", dpi=150)
        plt.close(fig)

    # ── 5. Ranked bar chart colored by germ layer ──
    fig, ax = plt.subplots(figsize=(12, 8))
    df_sorted = df.sort_values("residual_magnitude", ascending=True)
    colors = df_sorted["germ_layer"].map(palette_germ)
    bars = ax.barh(range(len(df_sorted)), df_sorted["residual_magnitude"],
                   color=colors, edgecolor="black", linewidth=0.3)
    ax.set_yticks(range(len(df_sorted)))
    ax.set_yticklabels(df_sorted["cell_type"], fontsize=8)
    ax.set_xlabel("Procrustes Residual Magnitude", fontsize=12)
    ax.set_title("Ranked Evolutionary Divergence — Colored by Germ Layer", fontsize=14)
    # Mark progenitors with a star
    for i, (_, row) in enumerate(df_sorted.iterrows()):
        if row["progenitor"]:
            ax.text(row["residual_magnitude"] + 0.2, i, "★",
                    fontsize=10, va="center", color="#E74C3C")
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=c, label=g) for g, c in palette_germ.items()
                       if g in df["germ_layer"].values]
    legend_elements.append(plt.Line2D([0], [0], marker='*', color='w',
                           markerfacecolor='#E74C3C', markersize=12,
                           label='Progenitor/Stem'))
    ax.legend(handles=legend_elements, loc="lower right", fontsize=9)
    fig.tight_layout()
    fig.savefig(outdir / "ranked_bars_germ_layer.png", dpi=150)
    plt.close(fig)

    print("  Plots saved to", outdir)


def print_summary(df, results):
    """Print human-readable summary table and statistical results."""
    print("\n" + "=" * 100)
    print("DEVELOPMENTAL CONSTRAINT ANALYSIS — SUMMARY")
    print("=" * 100)

    # Summary table
    print("\n── Per-Cell-Type Summary ──")
    print(f"{'Rank':<5} {'Cell Type':<45} {'Germ Layer':<12} {'Lineage':<15} "
          f"{'Prog?':<6} {'Residual':>10}")
    print("-" * 100)
    for _, row in df.sort_values("rank").iterrows():
        prog_str = "YES" if row["progenitor"] else ""
        print(f"{int(row['rank']):<5} {row['cell_type']:<45} {row['germ_layer']:<12} "
              f"{row['lineage']:<15} {prog_str:<6} {row['residual_magnitude']:>10.2f}")

    # Group summaries
    print("\n── Group Medians ──")
    print("\nBy Germ Layer:")
    for g, grp in df.groupby("germ_layer"):
        vals = grp["residual_magnitude"]
        print(f"  {g:<12}: median={vals.median():.2f}, mean={vals.mean():.2f}, "
              f"n={len(vals)}, range=[{vals.min():.2f}, {vals.max():.2f}]")

    print("\nBy Lineage:")
    for l, grp in df.groupby("lineage"):
        vals = grp["residual_magnitude"]
        print(f"  {l:<15}: median={vals.median():.2f}, mean={vals.mean():.2f}, "
              f"n={len(vals)}, range=[{vals.min():.2f}, {vals.max():.2f}]")

    print("\nBy Progenitor Status:")
    for p, grp in df.groupby("progenitor"):
        label = "Progenitor" if p else "Differentiated"
        vals = grp["residual_magnitude"]
        print(f"  {label:<15}: median={vals.median():.2f}, mean={vals.mean():.2f}, "
              f"n={len(vals)}, range=[{vals.min():.2f}, {vals.max():.2f}]")

    # Statistical tests
    print("\n── Statistical Tests ──")
    if "germ_layer_kw" in results:
        r = results["germ_layer_kw"]
        sig = "SIGNIFICANT" if r["significant_005"] else "not significant"
        print(f"\n  Kruskal-Wallis (germ layer): H={r['H_statistic']:.3f}, "
              f"p={r['p_value']:.4f} — {sig} at α=0.05")

    if "lineage_kw" in results:
        r = results["lineage_kw"]
        sig = "SIGNIFICANT" if r["significant_005"] else "not significant"
        print(f"  Kruskal-Wallis (lineage):    H={r['H_statistic']:.3f}, "
              f"p={r['p_value']:.4f} — {sig} at α=0.05")

    if "progenitor_mwu" in results:
        r = results["progenitor_mwu"]
        sig = "SIGNIFICANT" if r["significant_005"] else "not significant"
        print(f"  Mann-Whitney U (progenitor): U={r['U_statistic']:.1f}, "
              f"p={r['p_value']:.4f} — {sig} at α=0.05")
        print(f"    Progenitor median={r['progenitor_median']:.2f} (n={r['progenitor_n']}), "
              f"Differentiated median={r['differentiated_median']:.2f} (n={r['differentiated_n']})")

    if "germ_layer_pairwise" in results:
        print("\n  Pairwise Mann-Whitney U (germ layers):")
        for pw in results["germ_layer_pairwise"]:
            sig = "*" if pw["p"] < 0.05 else ""
            print(f"    {pw['group1']} vs {pw['group2']}: "
                  f"U={pw['U']:.1f}, p={pw['p']:.4f} "
                  f"(medians: {pw['median1']:.2f} vs {pw['median2']:.2f}) {sig}")

    if "mki67_correlation" in results:
        r = results["mki67_correlation"]
        sig = "SIGNIFICANT" if r["significant_005"] else "not significant"
        print(f"\n  Spearman (residual vs MKI67 human): ρ={r['rho']:.3f}, "
              f"p={r['p_value']:.4f} — {sig} at α=0.05")

    if "mki67_mean_correlation" in results:
        r = results["mki67_mean_correlation"]
        sig = "SIGNIFICANT" if r["significant_005"] else "not significant"
        print(f"  Spearman (residual vs MKI67 mean):  ρ={r['rho']:.3f}, "
              f"p={r['p_value']:.4f} — {sig} at α=0.05")

    # Interpretation
    print("\n── Interpretation ──")
    all_ns = True
    key_findings = []

    if "progenitor_mwu" in results and results["progenitor_mwu"]["significant_005"]:
        all_ns = False
        r = results["progenitor_mwu"]
        direction = ("MORE diverged" if r["progenitor_median"] > r["differentiated_median"]
                     else "LESS diverged")
        key_findings.append(
            f"Progenitor/stem cells are {direction} than differentiated cells "
            f"(median {r['progenitor_median']:.2f} vs {r['differentiated_median']:.2f}, "
            f"p={r['p_value']:.4f}).")
    elif "progenitor_mwu" in results:
        r = results["progenitor_mwu"]
        direction = ("more" if r["progenitor_median"] > r["differentiated_median"]
                     else "less")
        key_findings.append(
            f"Progenitor/stem cells trend {direction} diverged than differentiated "
            f"(median {r['progenitor_median']:.2f} vs {r['differentiated_median']:.2f}) "
            f"but not significant (p={r['p_value']:.4f}).")

    if "germ_layer_kw" in results and results["germ_layer_kw"]["significant_005"]:
        all_ns = False
        key_findings.append(
            f"Germ layer of origin significantly predicts divergence "
            f"(H={results['germ_layer_kw']['H_statistic']:.2f}, "
            f"p={results['germ_layer_kw']['p_value']:.4f}).")
    elif "germ_layer_kw" in results:
        key_findings.append(
            f"No significant difference across germ layers "
            f"(p={results['germ_layer_kw']['p_value']:.4f}).")

    if "lineage_kw" in results and results["lineage_kw"]["significant_005"]:
        all_ns = False
        key_findings.append(
            f"Developmental lineage significantly predicts divergence "
            f"(H={results['lineage_kw']['H_statistic']:.2f}, "
            f"p={results['lineage_kw']['p_value']:.4f}).")
    elif "lineage_kw" in results:
        key_findings.append(
            f"No significant difference across lineages "
            f"(p={results['lineage_kw']['p_value']:.4f}).")

    if "mki67_mean_correlation" in results:
        r = results["mki67_mean_correlation"]
        if r["significant_005"]:
            all_ns = False
            direction = "positively" if r["rho"] > 0 else "negatively"
            key_findings.append(
                f"Proliferative capacity (MKI67) {direction} correlates with "
                f"divergence (ρ={r['rho']:.3f}, p={r['p_value']:.4f}).")
        else:
            key_findings.append(
                f"No significant correlation between proliferation and divergence "
                f"(ρ={r['rho']:.3f}, p={r['p_value']:.4f}).")

    for f in key_findings:
        print(f"  • {f}")

    if all_ns:
        print("\n  CONCLUSION: Evolutionary divergence appears INDEPENDENT of developmental "
              "origin in this dataset. The Procrustes transformation warps cell types "
              "regardless of their germ layer, lineage, or differentiation status. This "
              "suggests that cross-species expression changes are driven by cell-type-"
              "specific regulatory evolution rather than shared developmental constraints.")
    else:
        print("\n  CONCLUSION: Developmental origin shows a statistically significant "
              "association with evolutionary divergence. See findings above for details.")

    print("=" * 100)


def main():
    base = Path(__file__).resolve().parent.parent
    residuals_path = base / "output/phase2/scaled_35types/residuals_ranked.csv"
    centroids_human = base / "output/phase2/scaled_35types/centroids_human_35.csv"
    centroids_mouse = base / "output/phase2/scaled_35types/centroids_mouse_35.csv"
    outdir = base / "output/phase2/developmental_constraint"
    outdir.mkdir(parents=True, exist_ok=True)

    # 1. Load residuals
    print("Loading 35-type Procrustes residuals...")
    residuals = load_residuals(residuals_path)
    print(f"  {len(residuals)} cell types loaded.")

    # 2. Annotate
    print("Assigning developmental annotations...")
    df = annotate_cell_types(residuals)
    print(f"  {len(df)} cell types annotated.")
    print(f"  Germ layers: {df['germ_layer'].value_counts().to_dict()}")
    print(f"  Lineages: {df['lineage'].value_counts().to_dict()}")
    print(f"  Progenitors: {df['progenitor'].sum()}, "
          f"Differentiated: {(~df['progenitor']).sum()}")

    # 3. MKI67 proliferation proxy
    print("Extracting MKI67 expression from centroids...")
    mki67_human = estimate_proliferation(centroids_human)
    mki67_mouse = estimate_proliferation(centroids_mouse)
    if mki67_human:
        df["mki67_human"] = df["cell_type"].map(mki67_human)
        print(f"  Human MKI67 range: [{df['mki67_human'].min():.4f}, {df['mki67_human'].max():.4f}]")
    if mki67_mouse:
        df["mki67_mouse"] = df["cell_type"].map(mki67_mouse)
        print(f"  Mouse MKI67 range: [{df['mki67_mouse'].min():.4f}, {df['mki67_mouse'].max():.4f}]")
    if mki67_human and mki67_mouse:
        df["mki67_mean"] = (df["mki67_human"] + df["mki67_mouse"]) / 2
        print(f"  Mean MKI67 range: [{df['mki67_mean'].min():.4f}, {df['mki67_mean'].max():.4f}]")

    # 4. Statistical tests
    print("\nRunning statistical tests...")
    results = run_statistical_tests(df)

    # 5. Generate plots
    print("Generating plots...")
    make_plots(df, outdir)

    # 6. Print summary
    print_summary(df, results)

    # 7. Save outputs
    df.to_csv(outdir / "developmental_annotations.csv", index=False)
    print(f"\nAnnotated table saved to {outdir / 'developmental_annotations.csv'}")

    # Save statistical results as text report
    import json

    def sanitize(obj):
        """Recursively convert numpy types for JSON serialization."""
        if isinstance(obj, dict):
            return {k: sanitize(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [sanitize(v) for v in obj]
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating, np.bool_)):
            return float(obj) if not isinstance(obj, np.bool_) else bool(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj

    with open(outdir / "statistical_tests.json", "w") as f:
        json.dump(sanitize(results), f, indent=2)
    print(f"Statistical results saved to {outdir / 'statistical_tests.json'}")

    print("\nDone.")


if __name__ == "__main__":
    main()
