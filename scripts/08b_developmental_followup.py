#!/usr/bin/env python3
"""
Follow-up on developmental constraint analysis: bootstrap validation and
continuous stemness correlation.

Biology:
    The binary progenitor vs differentiated comparison (p=0.0099) relies on
    only n=5 progenitors. Two follow-ups strengthen or refute this finding:

    1. Bootstrap validation: Resample with replacement 10,000 times from each
       group to estimate the robustness of the Mann-Whitney U test. Reports
       the 95% CI for the p-value and fraction of bootstraps where p < 0.05.

    2. Continuous stemness score: Instead of a binary classification, compute
       a continuous stemness score from expression of stemness-associated genes
       across all 35 cell types. This converts the n=5 vs n=30 binary test
       into a rank correlation across n=35, which is much more powerful.

    Stemness markers (adult tissue-relevant subset of pluripotency/self-renewal
    genes available in the 1:1 ortholog set):
      Core TFs: SOX2, KLF4, MYC, BMI1
      Surface:  CD44, PROM1 (CD133)
      Telomere: TERT
      RNA:      LIN28A, LIN28B
      Signaling: STAT3, LIFR
      Cytoskeleton: NES (nestin)

Math:
    Stemness score = mean of per-gene z-scores across 35 cell types.
    Z-scoring normalizes each gene to have mean=0, sd=1 across cell types,
    preventing highly-expressed genes from dominating. The average of human
    and mouse stemness scores is used to reduce species-specific noise.
    Correlation: Spearman rank correlation (non-parametric, appropriate for
    non-normal distributions).

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
import json


# ─── Stemness gene set ────────────────────────────────────────────────────────
# Curated adult stemness markers with Ensembl IDs from the ortholog table.
STEMNESS_GENES = {
    "ENSG00000181449": "SOX2",     # neural/epithelial stem cell TF
    "ENSG00000136826": "KLF4",     # intestinal/epithelial stem cell TF
    "ENSG00000136997": "MYC",      # proto-oncogene, proliferation/stemness
    "ENSG00000168283": "BMI1",     # polycomb, HSC self-renewal
    "ENSG00000026508": "CD44",     # adhesion, cancer stem cell marker
    "ENSG00000007062": "PROM1",    # CD133, stem cell surface marker
    "ENSG00000164362": "TERT",     # telomerase, stem cell maintenance
    "ENSG00000131914": "LIN28A",   # RNA-binding, stem cell maintenance
    "ENSG00000187772": "LIN28B",   # RNA-binding, stem cell maintenance
    "ENSG00000168610": "STAT3",    # JAK-STAT signaling, stemness
    "ENSG00000113594": "LIFR",     # LIF receptor, stem cell signaling
    "ENSG00000132688": "NES",      # nestin, neural/mesenchymal stem cells
}


def bootstrap_progenitor_test(df, n_boot=10000, seed=42):
    """
    Bootstrap the Mann-Whitney U test for progenitor vs differentiated.

    Resamples with replacement from each group independently, recomputes
    Mann-Whitney U each time. Reports:
    - Distribution of U statistics
    - Distribution of p-values
    - 95% CI for p-value
    - Fraction of bootstraps where p < 0.05 (robustness probability)
    """
    rng = np.random.default_rng(seed)

    prog = df[df["progenitor"] == True]["residual_magnitude"].values
    diff = df[df["progenitor"] == False]["residual_magnitude"].values
    n_prog, n_diff = len(prog), len(diff)

    # Original test
    u_orig, p_orig = stats.mannwhitneyu(prog, diff, alternative="two-sided")

    boot_u = np.zeros(n_boot)
    boot_p = np.zeros(n_boot)

    for i in range(n_boot):
        prog_boot = rng.choice(prog, size=n_prog, replace=True)
        diff_boot = rng.choice(diff, size=n_diff, replace=True)
        u, p = stats.mannwhitneyu(prog_boot, diff_boot, alternative="two-sided")
        boot_u[i] = u
        boot_p[i] = p

    # 95% CI for p-value
    p_ci_lower = np.percentile(boot_p, 2.5)
    p_ci_upper = np.percentile(boot_p, 97.5)

    # Fraction where p < 0.05
    frac_sig = np.mean(boot_p < 0.05)

    results = {
        "original_U": float(u_orig),
        "original_p": float(p_orig),
        "bootstrap_n": n_boot,
        "p_value_95CI": [float(p_ci_lower), float(p_ci_upper)],
        "p_value_median": float(np.median(boot_p)),
        "p_value_mean": float(np.mean(boot_p)),
        "fraction_significant_005": float(frac_sig),
        "U_median": float(np.median(boot_u)),
        "U_95CI": [float(np.percentile(boot_u, 2.5)),
                    float(np.percentile(boot_u, 97.5))],
    }

    return results, boot_p, boot_u


def compute_stemness_score(centroids_path, genes=STEMNESS_GENES):
    """
    Compute continuous stemness score from centroid expression of stemness markers.

    For each gene, z-score expression across cell types (mean=0, sd=1).
    Stemness score = mean z-score across all available stemness genes.
    This normalizes for different expression scales across genes.
    """
    df = pd.read_csv(centroids_path, index_col=0)

    # Find which genes are available
    available = {eid: name for eid, name in genes.items() if eid in df.columns}
    missing = {eid: name for eid, name in genes.items() if eid not in df.columns}

    if missing:
        print(f"  Missing stemness genes: {list(missing.values())}")
    print(f"  Available stemness genes: {list(available.values())} ({len(available)}/{len(genes)})")

    if not available:
        return None, {}, []

    # Extract expression for available genes
    expr = df[list(available.keys())].copy()
    expr.columns = [available[c] for c in expr.columns]

    # Z-score each gene across 35 cell types
    expr_z = (expr - expr.mean()) / expr.std()

    # Stemness score = mean z-score
    scores = expr_z.mean(axis=1).to_dict()

    # Also return per-gene z-scores for inspection
    per_gene = expr_z.to_dict(orient='index')

    return scores, per_gene, list(available.values())


def plot_bootstrap(boot_p, results, outdir):
    """Plot bootstrap p-value distribution."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Left: p-value distribution
    ax = axes[0]
    ax.hist(boot_p, bins=50, color="#3498DB", edgecolor="black", linewidth=0.3, alpha=0.8)
    ax.axvline(0.05, color="red", linestyle="--", linewidth=2, label="α = 0.05")
    ax.axvline(results["original_p"], color="green", linestyle="-", linewidth=2,
               label=f"Original p = {results['original_p']:.4f}")
    ax.set_xlabel("Bootstrap p-value", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.set_title("Bootstrap Distribution of Mann-Whitney p-value", fontsize=13)
    ax.legend(fontsize=10)
    frac = results["fraction_significant_005"]
    ax.text(0.98, 0.95, f"{frac*100:.1f}% of bootstraps\nhave p < 0.05",
            transform=ax.transAxes, ha="right", va="top", fontsize=11,
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    # Right: cumulative distribution of p-values
    ax = axes[1]
    sorted_p = np.sort(boot_p)
    cumulative = np.arange(1, len(sorted_p) + 1) / len(sorted_p)
    ax.plot(sorted_p, cumulative, color="#2ECC71", linewidth=2)
    ax.axvline(0.05, color="red", linestyle="--", linewidth=1.5, label="α = 0.05")
    ax.set_xlabel("p-value threshold", fontsize=12)
    ax.set_ylabel("Fraction of bootstraps ≤ threshold", fontsize=12)
    ax.set_title("Cumulative Distribution of Bootstrap p-values", fontsize=13)
    ax.legend(fontsize=10)
    ax.set_xlim(0, 0.5)

    fig.tight_layout()
    fig.savefig(outdir / "bootstrap_progenitor_pvalues.png", dpi=150)
    plt.close(fig)


def plot_stemness_vs_residual(df, outdir):
    """Scatter: stemness score vs Procrustes residual magnitude."""
    fig, ax = plt.subplots(figsize=(10, 8))

    palette = {"ectoderm": "#E74C3C", "mesoderm": "#3498DB",
               "endoderm": "#2ECC71", "mixed": "#95A5A6"}
    colors = df["germ_layer"].map(palette)

    # Marker style: progenitors get a different marker
    for _, row in df.iterrows():
        marker = "D" if row["progenitor"] else "o"
        size = 80 if row["progenitor"] else 50
        ax.scatter(row["stemness_score"], row["residual_magnitude"],
                   c=palette[row["germ_layer"]], s=size, marker=marker,
                   edgecolors="black", linewidth=0.5, zorder=3)

    # Label all points
    texts = []
    for _, row in df.iterrows():
        label = row["cell_type"]
        if len(label) > 30:
            label = label[:27] + "..."
        texts.append(ax.annotate(
            label, (row["stemness_score"], row["residual_magnitude"]),
            fontsize=6.5, alpha=0.8,
            xytext=(5, 3), textcoords="offset points"))

    # Spearman correlation
    rho, p = stats.spearmanr(df["stemness_score"], df["residual_magnitude"])
    ax.text(0.02, 0.98, f"Spearman ρ = {rho:.3f}\np = {p:.4f}",
            transform=ax.transAxes, fontsize=11, va="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    ax.set_xlabel("Stemness Score (mean z-scored expression of 12 stemness markers)",
                  fontsize=11)
    ax.set_ylabel("Procrustes Residual Magnitude", fontsize=11)
    ax.set_title("Evolutionary Divergence vs Continuous Stemness Score", fontsize=14)

    # Legend
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    legend_elements = [
        Patch(facecolor=c, label=g) for g, c in palette.items()
        if g in df["germ_layer"].values
    ]
    legend_elements.append(Line2D([0], [0], marker='D', color='w',
                           markerfacecolor='gray', markersize=8,
                           markeredgecolor='black', label='Progenitor/Stem'))
    legend_elements.append(Line2D([0], [0], marker='o', color='w',
                           markerfacecolor='gray', markersize=8,
                           markeredgecolor='black', label='Differentiated'))
    ax.legend(handles=legend_elements, loc="upper right", fontsize=9)

    fig.tight_layout()
    fig.savefig(outdir / "scatter_stemness_vs_residual.png", dpi=150)
    plt.close(fig)


def main():
    base = Path(__file__).resolve().parent.parent
    outdir = base / "output/phase2/developmental_constraint"
    centroids_human = base / "output/phase2/scaled_35types/centroids_human_35.csv"
    centroids_mouse = base / "output/phase2/scaled_35types/centroids_mouse_35.csv"
    annotations_path = outdir / "developmental_annotations.csv"

    # Load previous annotations
    print("Loading annotations from previous analysis...")
    df = pd.read_csv(annotations_path)
    print(f"  {len(df)} cell types loaded.")

    # ═══════════════════════════════════════════════════════════════════════════
    # 1. BOOTSTRAP VALIDATION
    # ═══════════════════════════════════════════════════════════════════════════
    print("\n" + "=" * 80)
    print("1. BOOTSTRAP: Progenitor vs Differentiated (10,000 resamples)")
    print("=" * 80)

    boot_results, boot_p, boot_u = bootstrap_progenitor_test(df, n_boot=10000)

    print(f"\n  Original test: U={boot_results['original_U']:.1f}, "
          f"p={boot_results['original_p']:.4f}")
    print(f"  Bootstrap p-value: median={boot_results['p_value_median']:.4f}, "
          f"mean={boot_results['p_value_mean']:.4f}")
    print(f"  95% CI for p-value: [{boot_results['p_value_95CI'][0]:.4f}, "
          f"{boot_results['p_value_95CI'][1]:.4f}]")
    print(f"  Fraction of bootstraps with p < 0.05: "
          f"{boot_results['fraction_significant_005']*100:.1f}%")
    print(f"  Bootstrap U: median={boot_results['U_median']:.1f}, "
          f"95% CI=[{boot_results['U_95CI'][0]:.1f}, {boot_results['U_95CI'][1]:.1f}]")

    plot_bootstrap(boot_p, boot_results, outdir)
    print("  Plot saved: bootstrap_progenitor_pvalues.png")

    # ═══════════════════════════════════════════════════════════════════════════
    # 2. CONTINUOUS STEMNESS SCORE
    # ═══════════════════════════════════════════════════════════════════════════
    print("\n" + "=" * 80)
    print("2. CONTINUOUS STEMNESS SCORE")
    print("=" * 80)

    print("\n  Computing human stemness scores...")
    scores_human, per_gene_human, genes_used = compute_stemness_score(centroids_human)
    print("  Computing mouse stemness scores...")
    scores_mouse, per_gene_mouse, _ = compute_stemness_score(centroids_mouse)

    if scores_human and scores_mouse:
        # Average human and mouse stemness scores
        df["stemness_human"] = df["cell_type"].map(scores_human)
        df["stemness_mouse"] = df["cell_type"].map(scores_mouse)
        df["stemness_score"] = (df["stemness_human"] + df["stemness_mouse"]) / 2

        print(f"\n  Stemness score range: [{df['stemness_score'].min():.3f}, "
              f"{df['stemness_score'].max():.3f}]")

        # Top/bottom 5
        df_sorted = df.sort_values("stemness_score", ascending=False)
        print("\n  Top 5 highest stemness:")
        for _, row in df_sorted.head(5).iterrows():
            prog = " (progenitor)" if row["progenitor"] else ""
            print(f"    {row['cell_type']:<45} stemness={row['stemness_score']:.3f}"
                  f"  residual={row['residual_magnitude']:.2f}{prog}")

        print("\n  Bottom 5 lowest stemness:")
        for _, row in df_sorted.tail(5).iterrows():
            prog = " (progenitor)" if row["progenitor"] else ""
            print(f"    {row['cell_type']:<45} stemness={row['stemness_score']:.3f}"
                  f"  residual={row['residual_magnitude']:.2f}{prog}")

        # Spearman correlation
        rho, p = stats.spearmanr(df["stemness_score"], df["residual_magnitude"])
        print(f"\n  Spearman correlation (stemness vs residual):")
        print(f"    ρ = {rho:.3f}, p = {p:.4f}")
        sig = "SIGNIFICANT" if p < 0.05 else "not significant"
        print(f"    {sig} at α=0.05")

        stemness_results = {
            "genes_used": genes_used,
            "n_genes": len(genes_used),
            "method": "mean z-scored expression, averaged across human and mouse",
            "spearman_rho": float(rho),
            "spearman_p": float(p),
            "significant_005": bool(p < 0.05),
            "n_cell_types": len(df),
        }

        # Also test human and mouse separately
        rho_h, p_h = stats.spearmanr(df["stemness_human"], df["residual_magnitude"])
        rho_m, p_m = stats.spearmanr(df["stemness_mouse"], df["residual_magnitude"])
        print(f"\n  Human-only:  ρ = {rho_h:.3f}, p = {p_h:.4f}")
        print(f"  Mouse-only:  ρ = {rho_m:.3f}, p = {p_m:.4f}")
        stemness_results["human_only"] = {"rho": float(rho_h), "p": float(p_h)}
        stemness_results["mouse_only"] = {"rho": float(rho_m), "p": float(p_m)}

        # Plot
        plot_stemness_vs_residual(df, outdir)
        print("\n  Plot saved: scatter_stemness_vs_residual.png")

    else:
        stemness_results = {"error": "Could not compute stemness scores"}

    # ═══════════════════════════════════════════════════════════════════════════
    # 3. COMBINED REPORT
    # ═══════════════════════════════════════════════════════════════════════════
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    rho_str = f"{stemness_results['spearman_rho']:.3f}" if 'spearman_rho' in stemness_results else "N/A"
    p_str = f"{stemness_results['spearman_p']:.4f}" if 'spearman_p' in stemness_results else "N/A"
    genes_str = ', '.join(stemness_results.get('genes_used', []))
    print(f"""
  BOOTSTRAP VALIDATION:
    Original Mann-Whitney: U={boot_results['original_U']:.1f}, p={boot_results['original_p']:.4f}
    Bootstrap 95% CI for p: [{boot_results['p_value_95CI'][0]:.4f}, {boot_results['p_value_95CI'][1]:.4f}]
    Probability effect holds (p<0.05): {boot_results['fraction_significant_005']*100:.1f}%

  CONTINUOUS STEMNESS:
    Spearman rho = {rho_str}
    p = {p_str}
    Genes: {genes_str}
""")

    if boot_results["fraction_significant_005"] >= 0.80:
        print("  BOOTSTRAP VERDICT: ROBUST — progenitor effect holds in "
              f"≥{boot_results['fraction_significant_005']*100:.0f}% of resamples.")
    elif boot_results["fraction_significant_005"] >= 0.50:
        print("  BOOTSTRAP VERDICT: MODERATE — progenitor effect holds in "
              f"{boot_results['fraction_significant_005']*100:.0f}% of resamples. "
              "Sensitive to which cells are drawn.")
    else:
        print("  BOOTSTRAP VERDICT: FRAGILE — progenitor effect holds in only "
              f"{boot_results['fraction_significant_005']*100:.0f}% of resamples.")

    if stemness_results.get("significant_005"):
        print("  STEMNESS VERDICT: CONFIRMED — continuous stemness significantly "
              "correlates with divergence across all 35 types.")
    else:
        print("  STEMNESS VERDICT: NOT CONFIRMED — no significant continuous "
              "relationship between stemness and divergence.")

    # Save all results
    all_results = {
        "bootstrap": boot_results,
        "stemness": stemness_results,
    }
    with open(outdir / "followup_results.json", "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\n  Results saved to {outdir / 'followup_results.json'}")

    # Save updated annotations with stemness scores
    df.to_csv(outdir / "developmental_annotations.csv", index=False)
    print(f"  Updated annotations saved to {outdir / 'developmental_annotations.csv'}")

    print("\nDone.")


if __name__ == "__main__":
    main()
