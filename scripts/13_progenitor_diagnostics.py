#!/usr/bin/env python3
"""
CellWarp — Progenitor Divergence Diagnostics

Two diagnostic checks on the progenitor divergence finding (p=0.0099):

Diagnostic 1 — Mouse-bias check
    Is the mouse-biased expression (log2FC=-0.67 in progenitor top genes) specific
    to progenitors, or a global artifact? Computes mean log2FC across all 16,959
    genes for each of 35 cell types. If all types show similar mouse bias, it's
    technical (likely Smart-seq2 protocol mix). If progenitors are outliers, it's
    biological.

Diagnostic 2 — Cell cycle regression
    The hematopoietic progenitors (HSC + precursor) show cytoplasmic translation
    enrichment, which could reflect cycling cells rather than genuine divergence.
    Tests whether the progenitor-divergence association (residual ~ progenitor_status)
    survives after controlling for cell cycle score via partial Spearman correlation.

Math
----
Diagnostic 1:
    log2FC_g = log2(mouse_g + 1) - log2(human_g + 1)
    mean_log2FC_ct = (1/G) * sum_g(log2FC_g)

Diagnostic 2 (partial Spearman):
    1. Rank-transform both X (residual) and Y (progenitor_status)
    2. Regress out Z (cell_cycle_score) from both: X_res = X - Z*β_X, Y_res = Y - Z*β_Y
    3. Compute Spearman(X_res, Y_res) as the partial correlation
    4. Test significance via permutation (10,000 iterations)

Inputs:
    output/phase2/scaled_35types/centroids_human_35.csv
    output/phase2/scaled_35types/centroids_mouse_35.csv
    output/phase2/developmental_constraint/developmental_annotations.csv

Outputs (all in output/phase2/progenitor_analysis/diagnostics/):
    mouse_bias_barchart.png
    mouse_bias_results.json
    cell_cycle_results.json

Usage:
    python scripts/13_progenitor_diagnostics.py
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

CENTROIDS_HUMAN = Path("output/phase2/scaled_35types/centroids_human_35.csv")
CENTROIDS_MOUSE = Path("output/phase2/scaled_35types/centroids_mouse_35.csv")
ANNOTATIONS = Path("output/phase2/developmental_constraint/developmental_annotations.csv")
OUTPUT_DIR = Path("output/phase2/progenitor_analysis/diagnostics")

# ---------------------------------------------------------------------------
# Cell cycle genes (Seurat cc.genes.updated.2019)
# ---------------------------------------------------------------------------

S_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG",
    "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP",
    "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76",
    "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51",
    "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM",
    "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8",
]

G2M_GENES = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A",
    "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF",
    "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB",
    "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP",
    "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1",
    "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR",
    "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF",
    "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
]

ALL_CC_GENES = list(set(S_GENES + G2M_GENES))


def partial_spearman(x, y, z):
    """
    Partial Spearman correlation between x and y, controlling for z.

    Math: Rank-transform x, y, z. Regress z out of both x and y via OLS.
    Compute Pearson r on the residuals (= partial Spearman).

    Returns (partial_rho, p_value_permutation).
    """
    # Rank transform
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    rz = stats.rankdata(z)

    # Regress z out of x: x_res = x - z * (z'x / z'z)
    rz_centered = rz - rz.mean()
    beta_x = np.dot(rz_centered, rx - rx.mean()) / np.dot(rz_centered, rz_centered)
    rx_res = rx - rz * beta_x

    # Regress z out of y: y_res = y - z * (z'y / z'z)
    beta_y = np.dot(rz_centered, ry - ry.mean()) / np.dot(rz_centered, rz_centered)
    ry_res = ry - rz * beta_y

    # Partial correlation = Pearson on residuals
    partial_rho, _ = stats.pearsonr(rx_res, ry_res)

    # Permutation test for significance
    n_perm = 10_000
    rng = np.random.RandomState(42)
    count = 0
    for _ in range(n_perm):
        perm_idx = rng.permutation(len(x))
        rx_perm = stats.rankdata(x[perm_idx])
        beta_x_p = np.dot(rz_centered, rx_perm - rx_perm.mean()) / np.dot(rz_centered, rz_centered)
        rx_perm_res = rx_perm - rz * beta_x_p
        r_perm, _ = stats.pearsonr(rx_perm_res, ry_res)
        if abs(r_perm) >= abs(partial_rho):
            count += 1
    p_perm = (count + 1) / (n_perm + 1)

    return partial_rho, p_perm


def main() -> None:
    t_start = time.time()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("  PROGENITOR DIVERGENCE DIAGNOSTICS")
    print("=" * 70)

    # Load data
    print("\nLoading data...")
    annot = pd.read_csv(ANNOTATIONS)
    centroids_h = pd.read_csv(CENTROIDS_HUMAN, index_col=0)
    centroids_m = pd.read_csv(CENTROIDS_MOUSE, index_col=0)
    gene_ids = centroids_h.columns.tolist()
    n_genes = len(gene_ids)
    print(f"  {len(annot)} cell types, {n_genes} genes")

    # Get gene symbol mapping from h5ad
    import anndata as ad
    h5 = ad.read_h5ad("data/phase2_scaled/human_scaled.h5ad", backed="r")
    id_to_symbol = dict(zip(h5.var_names, h5.var["feature_name"]))
    h5.file.close()
    gene_symbols = [id_to_symbol.get(gid, gid) for gid in gene_ids]
    symbol_to_idx = {s: i for i, s in enumerate(gene_symbols)}

    # ==================================================================
    # DIAGNOSTIC 1: Mouse-bias check
    # ==================================================================
    print("\n" + "=" * 70)
    print("  DIAGNOSTIC 1: Is mouse-biased expression global or progenitor-specific?")
    print("=" * 70)

    # Ensure consistent cell type order
    ct_order = sorted(centroids_h.index.tolist())
    progenitor_set = set(annot[annot["progenitor"] == True]["cell_type"].tolist())

    # Compute mean log2FC per cell type
    log2fc_per_ct = {}
    for ct in ct_order:
        h_vals = centroids_h.loc[ct].values.astype(np.float64)
        m_vals = centroids_m.loc[ct].values.astype(np.float64)
        # log2(x + 1) for pseudocount
        log2fc = np.log2(m_vals + 1) - np.log2(h_vals + 1)
        log2fc_per_ct[ct] = {
            "mean": float(np.mean(log2fc)),
            "median": float(np.median(log2fc)),
            "std": float(np.std(log2fc)),
            "is_progenitor": ct in progenitor_set,
        }

    # Build DataFrame for plotting
    fc_df = pd.DataFrame([
        {"cell_type": ct, "mean_log2fc": v["mean"], "is_progenitor": v["is_progenitor"]}
        for ct, v in log2fc_per_ct.items()
    ]).sort_values("mean_log2fc")

    # Statistics
    all_log2fc = fc_df["mean_log2fc"].values
    prog_log2fc = fc_df[fc_df["is_progenitor"]]["mean_log2fc"].values
    diff_log2fc = fc_df[~fc_df["is_progenitor"]]["mean_log2fc"].values

    overall_mean = np.mean(all_log2fc)
    prog_mean = np.mean(prog_log2fc)
    diff_mean = np.mean(diff_log2fc)
    mw_stat, mw_p = stats.mannwhitneyu(prog_log2fc, diff_log2fc, alternative="two-sided")

    print(f"\n  Overall mean log2FC (mouse - human): {overall_mean:.4f}")
    print(f"  Progenitor mean log2FC:              {prog_mean:.4f}")
    print(f"  Differentiated mean log2FC:          {diff_mean:.4f}")
    print(f"  Mann-Whitney p-value:                {mw_p:.4f}")

    # Check if progenitors are outliers
    prog_ranks = []
    for i, row in fc_df.reset_index(drop=True).iterrows():
        if row["is_progenitor"]:
            prog_ranks.append(i + 1)  # 1-indexed rank (1 = most mouse-biased)
    print(f"\n  Progenitor ranks in sorted distribution (1=most mouse-biased):")
    for ct in fc_df[fc_df["is_progenitor"]]["cell_type"]:
        rank = list(fc_df["cell_type"]).index(ct) + 1
        val = log2fc_per_ct[ct]["mean"]
        print(f"    {ct:<45} rank={rank:>2}/35  log2FC={val:+.4f}")

    if abs(overall_mean) > 0.05:
        bias_conclusion = (
            f"GLOBAL BIAS detected: overall mean log2FC = {overall_mean:+.4f} "
            f"({'mouse' if overall_mean > 0 else 'human'}-biased across ALL types)"
        )
    else:
        bias_conclusion = "No global bias: overall mean log2FC near zero"

    if mw_p < 0.05:
        prog_conclusion = (
            f"Progenitors ARE significantly different from differentiated "
            f"(p={mw_p:.4f}, prog={prog_mean:+.4f} vs diff={diff_mean:+.4f})"
        )
    else:
        prog_conclusion = (
            f"Progenitors are NOT significantly different from differentiated "
            f"(p={mw_p:.4f}, prog={prog_mean:+.4f} vs diff={diff_mean:+.4f})"
        )

    print(f"\n  CONCLUSION:")
    print(f"    {bias_conclusion}")
    print(f"    {prog_conclusion}")

    # Plot: bar chart
    fig, ax = plt.subplots(figsize=(14, 8))
    colors = ["#e74c3c" if row["is_progenitor"] else "#3498db"
              for _, row in fc_df.iterrows()]
    bars = ax.barh(range(len(fc_df)), fc_df["mean_log2fc"].values, color=colors,
                   edgecolor="white", linewidth=0.3)

    # Shorten cell type names for display
    short_names = []
    for ct in fc_df["cell_type"]:
        name = ct.replace("CD4-positive, alpha-beta ", "CD4+ ")
        name = name.replace("CD8-positive, alpha-beta ", "CD8+ ")
        name = name.replace("mesenchymal stem cell of adipose tissue", "MSC (adipose)")
        name = name.replace("mesenchymal stem cell", "MSC")
        name = name.replace("hematopoietic precursor cell", "Hemato. precursor")
        name = name.replace("hematopoietic stem cell", "HSC")
        name = name.replace("luminal epithelial cell of mammary gland", "Mammary luminal")
        name = name.replace("enterocyte of epithelium of large intestine", "Colon enterocyte")
        name = name.replace("large intestine goblet cell", "Goblet cell")
        name = name.replace("fibroblast of cardiac tissue", "Cardiac fibroblast")
        name = name.replace("bladder urothelial cell", "Urothelial cell")
        name = name.replace("pancreatic ductal cell", "Pancreatic ductal")
        name = name.replace("pancreatic acinar cell", "Pancreatic acinar")
        name = name.replace("myeloid dendritic cell", "mDC")
        name = name.replace("classical monocyte", "Classical mono.")
        name = name.replace("intermediate monocyte", "Intermediate mono.")
        name = name.replace("non-classical monocyte", "NC monocyte")
        name = name.replace("natural killer cell", "NK cell")
        name = name.replace("mature NK T cell", "NKT cell")
        name = name.replace("myeloid leukocyte", "Myeloid leukocyte")
        name = name.replace("smooth muscle cell", "Smooth muscle")
        short_names.append(name)

    ax.set_yticks(range(len(fc_df)))
    ax.set_yticklabels(short_names, fontsize=8)
    ax.axvline(0, color="black", linewidth=1.0, linestyle="-")
    ax.axvline(overall_mean, color="grey", linewidth=1.0, linestyle="--", alpha=0.7,
               label=f"Overall mean ({overall_mean:+.4f})")
    ax.set_xlabel("Mean log₂(mouse+1) − log₂(human+1)", fontsize=11)
    ax.set_title("Mouse vs Human Expression Bias Across 35 Cell Types", fontsize=13)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#e74c3c", label=f"Progenitor (n=5, mean={prog_mean:+.4f})"),
        Patch(facecolor="#3498db", label=f"Differentiated (n=30, mean={diff_mean:+.4f})"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=9,
              framealpha=0.9)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "mouse_bias_barchart.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Plot saved: {OUTPUT_DIR / 'mouse_bias_barchart.png'}")
    print(f"  Description: Horizontal bar chart, 35 cell types sorted by mean log2FC.")
    print(f"  Red bars = progenitors, blue = differentiated. Vertical line at 0 and overall mean.")

    # Save results
    diag1_results = {
        "overall_mean_log2fc": overall_mean,
        "progenitor_mean_log2fc": float(prog_mean),
        "differentiated_mean_log2fc": float(diff_mean),
        "mannwhitney_p": float(mw_p),
        "bias_conclusion": bias_conclusion,
        "progenitor_conclusion": prog_conclusion,
        "per_cell_type": log2fc_per_ct,
    }
    with open(OUTPUT_DIR / "mouse_bias_results.json", "w") as f:
        json.dump(diag1_results, f, indent=2)

    # ==================================================================
    # DIAGNOSTIC 2: Cell cycle regression
    # ==================================================================
    print("\n" + "=" * 70)
    print("  DIAGNOSTIC 2: Does progenitor divergence survive cell cycle regression?")
    print("=" * 70)

    # Find cell cycle genes in our gene space
    cc_found = [g for g in ALL_CC_GENES if g in symbol_to_idx]
    cc_not_found = [g for g in ALL_CC_GENES if g not in symbol_to_idx]
    print(f"\n  Cell cycle genes: {len(cc_found)}/{len(ALL_CC_GENES)} found in gene space")
    if cc_not_found:
        print(f"  Missing: {', '.join(cc_not_found[:10])}{'...' if len(cc_not_found) > 10 else ''}")

    cc_indices = [symbol_to_idx[g] for g in cc_found]

    # Compute cell cycle score per cell type (mean of human + mouse)
    cc_scores = {}
    cc_scores_human = {}
    cc_scores_mouse = {}
    for ct in ct_order:
        h_vals = centroids_h.loc[ct].values.astype(np.float64)
        m_vals = centroids_m.loc[ct].values.astype(np.float64)
        h_cc = np.mean(h_vals[cc_indices])
        m_cc = np.mean(m_vals[cc_indices])
        cc_scores[ct] = (h_cc + m_cc) / 2
        cc_scores_human[ct] = h_cc
        cc_scores_mouse[ct] = m_cc

    # Prepare arrays (sorted by ct_order)
    residuals = np.array([annot[annot["cell_type"] == ct]["residual_magnitude"].values[0] for ct in ct_order])
    progenitor_binary = np.array([1.0 if ct in progenitor_set else 0.0 for ct in ct_order])
    cc_array = np.array([cc_scores[ct] for ct in ct_order])

    # First: show raw correlation between cell cycle and residual
    rho_cc_resid, p_cc_resid = stats.spearmanr(cc_array, residuals)
    print(f"\n  Cell cycle score vs residual: Spearman rho={rho_cc_resid:.3f}, p={p_cc_resid:.4f}")

    # Cell cycle score: progenitor vs differentiated
    prog_cc = cc_array[progenitor_binary == 1]
    diff_cc = cc_array[progenitor_binary == 0]
    mw_cc, p_cc = stats.mannwhitneyu(prog_cc, diff_cc, alternative="two-sided")
    print(f"  Cell cycle score: progenitor mean={np.mean(prog_cc):.4f}, "
          f"differentiated mean={np.mean(diff_cc):.4f}, MW p={p_cc:.4f}")

    # Original association: residual ~ progenitor (Mann-Whitney, as reference)
    prog_resid = residuals[progenitor_binary == 1]
    diff_resid = residuals[progenitor_binary == 0]
    _, p_orig = stats.mannwhitneyu(prog_resid, diff_resid, alternative="two-sided")
    rho_orig, p_rho_orig = stats.spearmanr(progenitor_binary, residuals)
    print(f"\n  Original association (no control):")
    print(f"    Spearman rho(progenitor, residual) = {rho_orig:.3f}, p = {p_rho_orig:.4f}")
    print(f"    Mann-Whitney p = {p_orig:.4f}")

    # Partial Spearman: residual ~ progenitor, controlling for cell cycle
    print(f"\n  Computing partial Spearman (10,000 permutations)...")
    partial_rho, partial_p = partial_spearman(residuals, progenitor_binary, cc_array)
    print(f"    Partial Spearman rho = {partial_rho:.3f}")
    print(f"    Permutation p-value  = {partial_p:.4f}")

    # Also try controlling for mean expression of human cell cycle genes only
    cc_human_array = np.array([cc_scores_human[ct] for ct in ct_order])
    partial_rho_h, partial_p_h = partial_spearman(residuals, progenitor_binary, cc_human_array)
    print(f"\n  Controlling for human cell cycle only:")
    print(f"    Partial Spearman rho = {partial_rho_h:.3f}, p = {partial_p_h:.4f}")

    # And mouse only
    cc_mouse_array = np.array([cc_scores_mouse[ct] for ct in ct_order])
    partial_rho_m, partial_p_m = partial_spearman(residuals, progenitor_binary, cc_mouse_array)
    print(f"  Controlling for mouse cell cycle only:")
    print(f"    Partial Spearman rho = {partial_rho_m:.3f}, p = {partial_p_m:.4f}")

    # Conclusion
    if partial_p < 0.05:
        cc_conclusion = (
            f"SURVIVES: Progenitor divergence signal persists after cell cycle "
            f"regression (partial rho={partial_rho:.3f}, p={partial_p:.4f}). "
            f"Cell cycle is NOT the confound."
        )
    else:
        cc_conclusion = (
            f"ABSORBED: Progenitor divergence signal is reduced to non-significance "
            f"after cell cycle regression (partial rho={partial_rho:.3f}, p={partial_p:.4f}). "
            f"Cell cycle MAY be a confound."
        )

    print(f"\n  CONCLUSION:")
    print(f"    {cc_conclusion}")

    # Print cell cycle scores for progenitor types
    print(f"\n  Cell cycle scores for progenitor types:")
    for ct in ct_order:
        if ct in progenitor_set:
            print(f"    {ct:<45} cc_score={cc_scores[ct]:.4f}  "
                  f"(H={cc_scores_human[ct]:.4f}, M={cc_scores_mouse[ct]:.4f})")

    # Save results
    diag2_results = {
        "n_cc_genes_found": len(cc_found),
        "n_cc_genes_total": len(ALL_CC_GENES),
        "cc_vs_residual": {"spearman_rho": float(rho_cc_resid), "p": float(p_cc_resid)},
        "cc_progenitor_vs_differentiated": {
            "prog_mean": float(np.mean(prog_cc)),
            "diff_mean": float(np.mean(diff_cc)),
            "mannwhitney_p": float(p_cc),
        },
        "original_association": {
            "spearman_rho": float(rho_orig),
            "spearman_p": float(p_rho_orig),
            "mannwhitney_p": float(p_orig),
        },
        "partial_correlation_mean_cc": {
            "partial_rho": float(partial_rho),
            "permutation_p": float(partial_p),
        },
        "partial_correlation_human_cc": {
            "partial_rho": float(partial_rho_h),
            "permutation_p": float(partial_p_h),
        },
        "partial_correlation_mouse_cc": {
            "partial_rho": float(partial_rho_m),
            "permutation_p": float(partial_p_m),
        },
        "conclusion": cc_conclusion,
    }
    with open(OUTPUT_DIR / "cell_cycle_results.json", "w") as f:
        json.dump(diag2_results, f, indent=2)

    # ==================================================================
    # SUMMARY
    # ==================================================================
    elapsed = time.time() - t_start
    print("\n" + "=" * 70)
    print("  FINAL SUMMARY")
    print("=" * 70)
    print(f"""
  DIAGNOSTIC 1 — Mouse Expression Bias:
    Overall mean log2FC: {overall_mean:+.4f}
    Progenitor: {prog_mean:+.4f}, Differentiated: {diff_mean:+.4f}
    Mann-Whitney p = {mw_p:.4f}
    → {bias_conclusion}
    → {prog_conclusion}

  DIAGNOSTIC 2 — Cell Cycle Regression:
    CC score vs residual: rho={rho_cc_resid:.3f}, p={p_cc_resid:.4f}
    Original: rho(prog, resid) = {rho_orig:.3f}, p = {p_rho_orig:.4f}
    After CC control: partial rho = {partial_rho:.3f}, p = {partial_p:.4f}
    → {cc_conclusion}

  Runtime: {elapsed:.1f}s
  Output: {OUTPUT_DIR}/
""")


if __name__ == "__main__":
    main()
