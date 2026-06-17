#!/usr/bin/env python3
"""
CellWarp Script 22 — H2 Gene Program Concentration Mediation Test

Tests whether gene program concentration (participation ratio from PCA on
identity genes) mediates the treeness-rigidity anticorrelation (DECISION-120).

Biology: Rigid cell types may have more concentrated identity programs —
variance dominated by few principal components — consistent with deep,
narrow Waddington attractors. Concentrated programs could simultaneously
constrain cross-species geometry (high rigidity) and degrade tree structure
(low treeness, because types cluster along shared dominant axes).

Math: Participation ratio PR = (sum(lambda_i))^2 / sum(lambda_i^2), where
lambda_i are PCA eigenvalues from within-type expression on identity genes.
Low PR = concentrated (few dominant dimensions). High PR = distributed.

Pre-registration: docs/preregistration_treeness_h2_concentration_2026-03-16.md
Label: Post-hoc exploratory, motivated by DECISION-121 result.

Steps:
  1. Reconstruct top-50 loading genes per cell type from PCA + residuals
  2. For each type: extract cells, filter to 50 identity genes, PCA, compute PR
  3. H2a: PR vs Procrustes rigidity (Spearman)
  4. H2b: PR vs treeness (Spearman)
  5. H2c: Partial correlation (rigidity vs treeness | PR), attenuation test
  6. Sensitivity: repeat with PC1 fraction instead of PR
  7. Scatter plots and outputs

Input:  data/phase2_scaled/human_scaled.h5ad
        output/phase2/scaled_35types/centroids_{human,mouse}_35.csv
        output/phase2/scaled_35types/procrustes_results_35.json
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
import anndata as ad
from scipy.stats import spearmanr, t as t_dist
from sklearn.decomposition import PCA

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# --- Config ---
OUTPUT_DIR = PROJECT_ROOT / "output" / "liang_wagner"
HUMAN_H5AD = PROJECT_ROOT / "data" / "phase2_scaled" / "human_scaled.h5ad"
CENTROID_HUMAN = (
    PROJECT_ROOT / "output" / "phase2" / "scaled_35types" / "centroids_human_35.csv"
)
CENTROID_MOUSE = (
    PROJECT_ROOT / "output" / "phase2" / "scaled_35types" / "centroids_mouse_35.csv"
)
RESULTS_JSON = (
    PROJECT_ROOT / "output" / "phase2" / "scaled_35types" / "procrustes_results_35.json"
)
RESIDUAL_PATH = (
    PROJECT_ROOT / "output" / "phase2" / "scaled_35types" / "residuals_ranked.csv"
)
TREENESS_PATH = OUTPUT_DIR / "treeness_scores_per_celltype.csv"

PCA_VARIANCE_THRESHOLD = 0.95
RANDOM_SEED = 42
N_TOP_GENES = 50


def partial_spearman(x, y, z):
    """Partial Spearman correlation between x and y controlling for z.

    Math: r_xy.z = (r_xy - r_xz * r_yz) / sqrt((1 - r_xz^2)(1 - r_yz^2))
    P-value via t-distribution with df = n - 3.
    """
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

    t_stat = partial_rho * np.sqrt(df / (1 - partial_rho**2 + 1e-15))
    p_val = 2 * t_dist.sf(np.abs(t_stat), df)

    return float(partial_rho), float(p_val)


def reconstruct_top_genes(cell_types, n_top=N_TOP_GENES):
    """Reconstruct top-N loading genes per cell type from PCA + residuals.

    Re-fits PCA on combined centroids (same as pipeline) and projects
    stored residual vectors back to gene space to get full gene rankings.

    Returns:
        Dict mapping cell_type -> list of ENSG gene IDs (top N by abs loading).
    """
    print("\n  Reconstructing top-50 loading genes per cell type...")

    # Load centroids
    human_df = pd.read_csv(CENTROID_HUMAN, index_col=0)
    mouse_df = pd.read_csv(CENTROID_MOUSE, index_col=0)
    gene_names = human_df.columns.tolist()  # ENSG IDs

    # Refit PCA on combined centroids (deterministic with same seed)
    human_mat = human_df.loc[cell_types].values
    mouse_mat = mouse_df.loc[cell_types].values
    combined = np.vstack([human_mat, mouse_mat])

    pca = PCA(
        n_components=PCA_VARIANCE_THRESHOLD,
        svd_solver="full",
        random_state=RANDOM_SEED,
    )
    pca.fit(combined)
    W = pca.components_  # (k, G)
    print(f"    PCA: {pca.n_components_} components, "
          f"{np.sum(pca.explained_variance_ratio_) * 100:.1f}% variance")

    # Load residual vectors from JSON
    with open(RESULTS_JSON, "r") as f:
        results = json.load(f)

    top_genes_per_type = {}
    for ct in cell_types:
        residual_pca = np.array(results["residuals"][ct]["vector_pca"])
        # Project to gene space
        gene_loadings = residual_pca @ W  # (G,)
        # Rank by absolute loading
        abs_loadings = np.abs(gene_loadings)
        top_indices = np.argsort(abs_loadings)[::-1][:n_top]
        top_genes_per_type[ct] = [gene_names[i] for i in top_indices]

    print(f"    Extracted top {n_top} genes for each of {len(cell_types)} cell types")
    return top_genes_per_type


def compute_participation_ratio(eigenvalues):
    """Compute participation ratio (effective dimensionality).

    PR = (sum(lambda_i))^2 / sum(lambda_i^2)
    Range: 1 (all variance on one PC) to p (uniform).
    """
    eigenvalues = eigenvalues[eigenvalues > 0]
    if len(eigenvalues) == 0:
        return 0.0
    sum_lambda = np.sum(eigenvalues)
    sum_lambda_sq = np.sum(eigenvalues**2)
    if sum_lambda_sq < 1e-15:
        return 0.0
    return float(sum_lambda**2 / sum_lambda_sq)


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
    print("CellWarp — H2 Gene Program Concentration Mediation Test")
    print("Pre-registration: docs/preregistration_treeness_h2_concentration_2026-03-16.md")
    print("=" * 70)

    # --- Load reference data ---
    residuals_df = pd.read_csv(RESIDUAL_PATH)
    treeness_df = pd.read_csv(TREENESS_PATH)

    cell_types = sorted(residuals_df["cell_type"].tolist())
    n = len(cell_types)
    print(f"\n  Cell types: {n}")

    # --- Step 1: Reconstruct top-50 loading genes ---
    top_genes_per_type = reconstruct_top_genes(cell_types, N_TOP_GENES)

    # --- Step 2: Compute PR and PC1 fraction per cell type ---
    print(f"\n{'=' * 70}")
    print("STEP 2 — Compute Participation Ratio per Cell Type")
    print("=" * 70)

    print(f"\n  Loading single-cell data: {HUMAN_H5AD}")
    adata = ad.read_h5ad(HUMAN_H5AD)
    print(f"  Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Gene name index: var_names are ENSG IDs
    gene_to_idx = {g: i for i, g in enumerate(adata.var_names)}

    records = []
    for ct in cell_types:
        # Extract cells for this type
        mask = adata.obs["cell_type"] == ct
        n_cells = int(mask.sum())

        # Get top-50 gene indices
        top_ensg = top_genes_per_type[ct]
        gene_indices = [gene_to_idx[g] for g in top_ensg if g in gene_to_idx]
        n_genes_found = len(gene_indices)

        if n_genes_found < 10:
            print(f"  WARNING: {ct} — only {n_genes_found} genes mapped, skipping")
            continue

        # Extract submatrix: cells x top genes
        X_sub = adata[mask][:, gene_indices].X
        if hasattr(X_sub, "toarray"):
            X_sub = X_sub.toarray()
        X_sub = np.asarray(X_sub, dtype=np.float64)

        # PCA on within-type expression (retain all components)
        n_components = min(n_cells - 1, n_genes_found)
        pca = PCA(n_components=n_components, random_state=RANDOM_SEED)
        pca.fit(X_sub)

        eigenvalues = pca.explained_variance_  # actual eigenvalues
        var_ratios = pca.explained_variance_ratio_

        pr = compute_participation_ratio(eigenvalues)
        pc1_frac = float(var_ratios[0]) if len(var_ratios) > 0 else 0.0

        records.append({
            "cell_type": ct,
            "participation_ratio": pr,
            "pc1_fraction": pc1_frac,
            "n_cells": n_cells,
            "n_genes": n_genes_found,
            "n_pca_components": n_components,
        })

        if ct in cell_types[:5] or ct in cell_types[-2:]:
            print(f"    {ct:<50} cells={n_cells:>5}  PR={pr:.2f}  "
                  f"PC1={pc1_frac:.3f}")

    concentration_df = pd.DataFrame(records)
    print(f"\n  Computed PR for {len(concentration_df)} / {n} cell types")

    # Print full ranking
    print(f"\n  Participation Ratio ranking (low PR = concentrated):")
    sorted_df = concentration_df.sort_values("participation_ratio")
    print(f"  {'Rank':<6} {'Cell Type':<50} {'PR':>8} {'PC1%':>8} {'Cells':>7}")
    print(f"  {'-' * 82}")
    for rank, (_, row) in enumerate(sorted_df.iterrows(), 1):
        print(f"  {rank:<6} {row['cell_type']:<50} "
              f"{row['participation_ratio']:>8.2f} "
              f"{row['pc1_fraction'] * 100:>7.1f}% "
              f"{int(row['n_cells']):>7}")

    # Save CSV
    concentration_df.to_csv(OUTPUT_DIR / "gene_concentration.csv", index=False)
    print(f"\n  Saved: {OUTPUT_DIR / 'gene_concentration.csv'}")

    # --- Merge all data ---
    merged = concentration_df.merge(
        residuals_df[["cell_type", "residual_magnitude"]], on="cell_type"
    )
    merged = merged.merge(
        treeness_df[["cell_type", "treeness_score"]], on="cell_type"
    )
    n_matched = len(merged)
    assert n_matched == n, f"Merge mismatch: {n_matched} != {n}"

    # --- Steps 3-5: Primary analysis (PR) ---
    print(f"\n{'=' * 70}")
    print("STEPS 3-5 — Primary Analysis (Participation Ratio)")
    print("=" * 70)

    pr_vals = merged["participation_ratio"].values
    residual_vals = merged["residual_magnitude"].values
    treeness_vals = merged["treeness_score"].values

    # H2a: PR vs residual_magnitude
    rho_h2a, p_h2a = spearmanr(pr_vals, residual_vals)
    h2a_pass = bool(rho_h2a > 0 and p_h2a < 0.05)
    print(f"\n  H2a: rho(PR, residual_magnitude) = {rho_h2a:.4f}, p = {p_h2a:.4f}")
    print(f"    Expected: positive (concentrated=low PR with low residual=rigid)")
    print(f"    H2a {'PASS' if h2a_pass else 'FAIL'} (need rho > 0, p < 0.05)")

    # H2b: PR vs treeness
    rho_h2b, p_h2b = spearmanr(pr_vals, treeness_vals)
    h2b_pass = bool(rho_h2b > 0 and p_h2b < 0.05)
    print(f"\n  H2b: rho(PR, treeness) = {rho_h2b:.4f}, p = {p_h2b:.4f}")
    print(f"    Expected: positive (concentrated=low PR with low treeness)")
    print(f"    H2b {'PASS' if h2b_pass else 'FAIL'} (need rho > 0, p < 0.05)")

    # H2c: Mediation
    rho_raw, p_raw = spearmanr(residual_vals, treeness_vals)
    partial_rho, partial_p = partial_spearman(residual_vals, treeness_vals, pr_vals)
    attenuation = (abs(rho_raw) - abs(partial_rho)) / abs(rho_raw) * 100
    h2c_pass = bool(attenuation >= 50)

    print(f"\n  H2c: Mediation test")
    print(f"    Raw rho(residual, treeness): {rho_raw:.4f} (p={p_raw:.4f})")
    print(f"    Partial rho(residual, treeness | PR): "
          f"{partial_rho:.4f} (p={partial_p:.4f})")
    print(f"    Attenuation: {attenuation:.1f}%")
    print(f"    H2c {'PASS' if h2c_pass else 'FAIL'} (need >= 50%)")

    # Verdict
    if h2a_pass and h2b_pass and h2c_pass:
        verdict_pr = "MEDIATION_CONFIRMED"
        desc_pr = ("Concentration mediates — attractor depth explains "
                   "the treeness-rigidity anticorrelation")
    elif h2a_pass and h2b_pass and not h2c_pass:
        verdict_pr = "PARTIAL_MEDIATION"
        desc_pr = ("Concentration contributes but does not fully explain — "
                   "additional component present")
    elif h2a_pass and not h2b_pass:
        verdict_pr = "CONCENTRATION_PREDICTS_RIGIDITY_ONLY"
        desc_pr = "Concentration predicts rigidity but not treeness"
    else:
        verdict_pr = "H2_REJECTED"
        desc_pr = "Concentration does not predict rigidity — H2 rejected"

    print(f"\n  Verdict (PR): {verdict_pr}")
    print(f"    {desc_pr}")

    # --- Step 6: Sensitivity (PC1 fraction) ---
    print(f"\n{'=' * 70}")
    print("STEP 6 — Sensitivity: PC1 Fraction")
    print("=" * 70)

    pc1_vals = merged["pc1_fraction"].values

    rho_s_h2a, p_s_h2a = spearmanr(pc1_vals, residual_vals)
    # PC1 fraction: HIGH = concentrated. Rigid types have low residual.
    # So expect NEGATIVE rho(pc1_frac, residual): concentrated (high PC1) = rigid (low residual)
    s_h2a_pass = bool(rho_s_h2a < 0 and p_s_h2a < 0.05)
    print(f"\n  H2a (PC1 frac): rho(PC1_frac, residual) = {rho_s_h2a:.4f}, "
          f"p = {p_s_h2a:.4f}")
    print(f"    Expected: negative (high PC1 frac = concentrated = low residual)")
    print(f"    {'PASS' if s_h2a_pass else 'FAIL'}")

    rho_s_h2b, p_s_h2b = spearmanr(pc1_vals, treeness_vals)
    # Concentrated (high PC1) = low treeness → expect negative rho
    s_h2b_pass = bool(rho_s_h2b < 0 and p_s_h2b < 0.05)
    print(f"\n  H2b (PC1 frac): rho(PC1_frac, treeness) = {rho_s_h2b:.4f}, "
          f"p = {p_s_h2b:.4f}")
    print(f"    Expected: negative (high PC1 frac = concentrated = low treeness)")
    print(f"    {'PASS' if s_h2b_pass else 'FAIL'}")

    partial_rho_s, partial_p_s = partial_spearman(
        residual_vals, treeness_vals, pc1_vals
    )
    attenuation_s = (abs(rho_raw) - abs(partial_rho_s)) / abs(rho_raw) * 100
    s_h2c_pass = bool(attenuation_s >= 50)
    print(f"\n  H2c (PC1 frac): partial rho = {partial_rho_s:.4f} "
          f"(p={partial_p_s:.4f}), attenuation = {attenuation_s:.1f}%")
    print(f"    {'PASS' if s_h2c_pass else 'FAIL'}")

    # Sensitivity verdict
    if s_h2a_pass and s_h2b_pass and s_h2c_pass:
        verdict_pc1 = "MEDIATION_CONFIRMED"
    elif s_h2a_pass and s_h2b_pass and not s_h2c_pass:
        verdict_pc1 = "PARTIAL_MEDIATION"
    elif s_h2a_pass and not s_h2b_pass:
        verdict_pc1 = "CONCENTRATION_PREDICTS_RIGIDITY_ONLY"
    else:
        verdict_pc1 = "H2_REJECTED"

    print(f"\n  Sensitivity verdict (PC1 frac): {verdict_pc1}")
    consistent = verdict_pr == verdict_pc1
    print(f"  Consistent with primary: {'YES' if consistent else 'NO'}")

    # --- Step 7: Plots ---
    print(f"\n{'=' * 70}")
    print("STEP 7 — Plots")
    print("=" * 70)

    labels = merged["cell_type"].values

    plot_scatter(
        pr_vals, residual_vals, labels,
        xlabel="Participation Ratio (higher = more distributed program)",
        ylabel="Procrustes Residual Magnitude (higher = less rigid)",
        title=(f"H2a: Gene Program Concentration vs Rigidity — "
               f"rho={rho_h2a:.3f}, p={p_h2a:.4f}"),
        output_path=OUTPUT_DIR / "pr_rigidity_scatter.png",
    )

    plot_scatter(
        pr_vals, treeness_vals, labels,
        xlabel="Participation Ratio (higher = more distributed program)",
        ylabel="Treeness Score (mean delta, higher = more tree-like)",
        title=(f"H2b: Gene Program Concentration vs Treeness — "
               f"rho={rho_h2b:.3f}, p={p_h2b:.4f}"),
        output_path=OUTPUT_DIR / "pr_treeness_scatter.png",
    )

    pr_quartile = pd.qcut(merged["participation_ratio"], q=4, labels=False) + 1
    plot_scatter(
        residual_vals, treeness_vals, labels,
        xlabel="Procrustes Residual Magnitude (higher = less rigid)",
        ylabel="Treeness Score (mean delta, higher = more tree-like)",
        title=(f"Rigidity vs Treeness by PR Quartile — "
               f"partial rho={partial_rho:.3f}, atten={attenuation:.1f}%"),
        output_path=OUTPUT_DIR / "rigidity_treeness_by_pr.png",
        color=pr_quartile.values,
        colorbar_label="PR Quartile (1=concentrated, 4=distributed)",
    )

    # --- Save results JSON ---
    results_json = {
        "primary_PR": {
            "h2a_pr_vs_residual": {
                "rho": float(rho_h2a), "p": float(p_h2a),
                "pass": h2a_pass,
            },
            "h2b_pr_vs_treeness": {
                "rho": float(rho_h2b), "p": float(p_h2b),
                "pass": h2b_pass,
            },
            "h2c_mediation": {
                "raw_rho": float(rho_raw), "raw_p": float(p_raw),
                "partial_rho": float(partial_rho), "partial_p": float(partial_p),
                "attenuation_pct": float(attenuation),
                "pass": h2c_pass,
            },
            "verdict": verdict_pr,
            "description": desc_pr,
        },
        "sensitivity_PC1": {
            "h2a_pc1_vs_residual": {
                "rho": float(rho_s_h2a), "p": float(p_s_h2a),
                "pass": s_h2a_pass,
            },
            "h2b_pc1_vs_treeness": {
                "rho": float(rho_s_h2b), "p": float(p_s_h2b),
                "pass": s_h2b_pass,
            },
            "h2c_mediation": {
                "partial_rho": float(partial_rho_s),
                "partial_p": float(partial_p_s),
                "attenuation_pct": float(attenuation_s),
                "pass": s_h2c_pass,
            },
            "verdict": verdict_pc1,
        },
        "consistent": consistent,
        "n": n_matched,
        "pr_summary": {
            "mean": float(np.mean(pr_vals)),
            "median": float(np.median(pr_vals)),
            "std": float(np.std(pr_vals)),
            "min": float(np.min(pr_vals)),
            "max": float(np.max(pr_vals)),
        },
    }

    with open(OUTPUT_DIR / "h2_mediation_results.json", "w") as f:
        json.dump(results_json, f, indent=2)
    print(f"  Saved: {OUTPUT_DIR / 'h2_mediation_results.json'}")

    # --- Write summary ---
    with open(OUTPUT_DIR / "h2_concentration_summary.md", "w") as f:
        f.write("# H2 Gene Program Concentration Mediation — Summary\n\n")
        f.write(f"**Date:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}\n")
        f.write("**Pre-registration:** docs/preregistration_treeness_h2_concentration_2026-03-16.md\n")
        f.write("**Label:** Post-hoc exploratory\n\n")

        f.write("## Results Table\n\n")
        f.write("| Test | PR (primary) | PC1 frac (sensitivity) |\n")
        f.write("|------|-------------|------------------------|\n")
        f.write(f"| H2a: concentration vs rigidity | "
                f"rho={rho_h2a:.3f} p={p_h2a:.3f} {'PASS' if h2a_pass else 'FAIL'} | "
                f"rho={rho_s_h2a:.3f} p={p_s_h2a:.3f} {'PASS' if s_h2a_pass else 'FAIL'} |\n")
        f.write(f"| H2b: concentration vs treeness | "
                f"rho={rho_h2b:.3f} p={p_h2b:.3f} {'PASS' if h2b_pass else 'FAIL'} | "
                f"rho={rho_s_h2b:.3f} p={p_s_h2b:.3f} {'PASS' if s_h2b_pass else 'FAIL'} |\n")
        f.write(f"| H2c: partial rho (attenuation) | "
                f"{partial_rho:.3f} ({attenuation:.0f}%) {'PASS' if h2c_pass else 'FAIL'} | "
                f"{partial_rho_s:.3f} ({attenuation_s:.0f}%) {'PASS' if s_h2c_pass else 'FAIL'} |\n")
        f.write(f"| Verdict | **{verdict_pr}** | {verdict_pc1} |\n")

        f.write(f"\n## Falsification Conditions\n\n")
        f.write(f"1. H2a NS (PR !~ rigidity): "
                f"{'TRIGGERED' if not h2a_pass else 'not triggered'}\n")
        f.write(f"2. H2b NS (PR !~ treeness): "
                f"{'TRIGGERED' if not h2b_pass else 'not triggered'}\n")
        f.write(f"3. Attenuation < 50%: "
                f"{'TRIGGERED' if not h2c_pass else 'not triggered'}\n")
        f.write(f"4. PR/PC1 inconsistent: "
                f"{'TRIGGERED' if not consistent else 'not triggered'}\n")

        f.write(f"\n## Primary Verdict\n\n")
        f.write(f"**{verdict_pr}**\n\n")
        f.write(f"{desc_pr}\n")

    print(f"  Saved: {OUTPUT_DIR / 'h2_concentration_summary.md'}")

    # --- Final summary ---
    print(f"\n{'=' * 70}")
    print("FINAL RESULTS TABLE")
    print("=" * 70)
    print(f"  {'Test':<45} {'PR (primary)':>20} {'PC1 (sensitivity)':>20}")
    print(f"  {'-' * 87}")
    print(f"  {'H2a: concentration vs rigidity':<45} "
          f"{'rho=' + f'{rho_h2a:.3f} p={p_h2a:.3f}':>20} "
          f"{'rho=' + f'{rho_s_h2a:.3f} p={p_s_h2a:.3f}':>20}")
    print(f"  {'H2b: concentration vs treeness':<45} "
          f"{'rho=' + f'{rho_h2b:.3f} p={p_h2b:.3f}':>20} "
          f"{'rho=' + f'{rho_s_h2b:.3f} p={p_s_h2b:.3f}':>20}")
    print(f"  {'H2c: partial rho (attenuation)':<45} "
          f"{f'{partial_rho:.3f} ({attenuation:.0f}%)':>20} "
          f"{f'{partial_rho_s:.3f} ({attenuation_s:.0f}%)':>20}")
    print(f"  {'Verdict':<45} {verdict_pr:>20} {verdict_pc1:>20}")
    print(f"\n  Raw rho(residual, treeness): {rho_raw:.4f} (p={p_raw:.4f})")
    print(f"  Consistent: {'YES' if consistent else 'NO'}")
    print("=" * 70)


if __name__ == "__main__":
    main()
