#!/usr/bin/env python3
"""
CellWarp — Progenitor Divergence Deep Characterization

Characterizes WHY progenitor/stem cells are more evolutionarily diverged than
differentiated cells (Mann-Whitney p=0.0099, established in developmental
constraint analysis). Identifies the genes, pathways, and biological processes
that drive the elevated divergence in progenitor types specifically.

Biology
-------
The project found that 5 progenitor/stem cell types have significantly higher
Procrustes residuals than the 30 differentiated types. This script asks:
  - What genes drive progenitor-specific divergence (vs general divergence)?
  - Are these genes enriched for specific GO biological processes?
  - Is the divergence driven by human upregulation, mouse upregulation, or both?
  - Do progenitor types share a common divergence program or diverge independently?

Math
----
Progenitor specificity score per gene:
    S_g = mean_prog(|r_g|) / (mean_prog(|r_g|) + mean_diff(|r_g|))

where r_g is the gene-space residual for gene g, averaged across progenitor (prog)
or differentiated (diff) types. S_g > 0.5 means the gene contributes more to
progenitor divergence; S_g > 0.7 is "progenitor-specific."

Inputs:
    output/phase2/scaled_35types/procrustes_results_35.json
    output/phase2/scaled_35types/centroids_human_35.csv
    output/phase2/scaled_35types/centroids_mouse_35.csv
    output/phase2/developmental_constraint/developmental_annotations.csv
    data/phase2_scaled/human_scaled.h5ad  (for gene name mapping)

Outputs (all in output/phase2/progenitor_analysis/):
    progenitor_specificity_scores.csv
    top50_progenitor_genes.csv
    go_enrichment/  (GO results for shared + individual progenitor types)
    expression_heatmap.png
    pairwise_overlap_heatmap.png
    vs_differentiated_comparison.csv
    summary.json

Usage:
    python scripts/12_progenitor_deep_dive.py
"""

from __future__ import annotations

import json
import sys
import time
import warnings
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import anndata as ad
import gseapy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA

warnings.filterwarnings("ignore", category=FutureWarning)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

RESULTS_JSON = Path("output/phase2/scaled_35types/procrustes_results_35.json")
CENTROIDS_HUMAN = Path("output/phase2/scaled_35types/centroids_human_35.csv")
CENTROIDS_MOUSE = Path("output/phase2/scaled_35types/centroids_mouse_35.csv")
ANNOTATIONS = Path("output/phase2/developmental_constraint/developmental_annotations.csv")
HUMAN_H5AD = Path("data/phase2_scaled/human_scaled.h5ad")
OUTPUT_DIR = Path("output/phase2/progenitor_analysis")
GO_DIR = OUTPUT_DIR / "go_enrichment"
PCA_VARIANCE_THRESHOLD = 0.95
RANDOM_SEED = 42

# Gene function categories for heuristic grouping (Step 4)
GENE_CATEGORIES = {
    "cell_cycle": ["MKI67", "TOP2A", "CCNB1", "CCNA2", "CDK1", "CDC20",
                    "PCNA", "MCM", "AURKB", "PLK1", "BUB1", "CENP",
                    "CDKN", "CCND", "CCNE", "RB1", "E2F"],
    "ribosomal": ["RPL", "RPS", "MRPL", "MRPS"],
    "signaling": ["WNT", "FGF", "BMP", "SHH", "NOTCH", "JAK", "STAT",
                   "MAPK", "EGFR", "TGFB", "SMAD", "PDGF", "VEGF",
                   "IGF", "PI3K", "AKT", "MTOR", "ERK", "RAF", "RAS",
                   "NFKB", "TNF", "IL1", "IL6", "IL10", "CXCL", "CCL"],
    "transcription_factors": ["SOX", "PAX", "HOX", "GATA", "RUNX", "ETS",
                                "KLF", "MYC", "JUN", "FOS", "ATF", "CREB",
                                "HIF", "NF", "IRF", "POU", "LHX", "DLX",
                                "TBX", "FOXP", "FOXA", "FOXO"],
    "chromatin": ["HDAC", "HAT", "KDM", "KMT", "DNMT", "TET", "EZH",
                   "SUZ12", "ARID", "SWI", "BRD", "CHD", "SMARCA"],
    "structural": ["COL", "VIM", "FN1", "LAMA", "LAMB", "LAMC", "KRT",
                    "ACTA", "ACTB", "ACTG", "TUB", "MYH", "MYL",
                    "DES", "NES", "GFAP", "SYNE"],
    "metabolic": ["CYP", "ALD", "ATP", "NADH", "SDH", "IDH", "PDH",
                   "PKM", "ENO", "HK", "PFK", "LDH", "GLUT", "SLC"],
    "immune": ["HLA", "CD", "IL", "IFN", "TLR", "MHC", "TCR", "BCR",
               "CTLA", "PD", "LAG", "TIM", "TIGIT"],
    "growth_factors": ["EGF", "FGF", "IGF", "PDGF", "VEGF", "TGF",
                        "HGF", "NGF", "BDNF", "CSF", "EPO", "TPO"],
    "stem_cell": ["NANOG", "POU5F1", "OCT4", "SOX2", "KLF4", "LIN28",
                   "BMI1", "TERT", "CD34", "PROM1", "ALDH", "ABCG2",
                   "LGR5", "OLFM4"],
}


def categorize_gene(symbol: str) -> str:
    """Categorize gene by name prefix heuristics. Returns category or 'other'."""
    upper = symbol.upper()
    for cat, prefixes in GENE_CATEGORIES.items():
        for prefix in prefixes:
            if upper.startswith(prefix) or upper == prefix:
                return cat
    return "other"


def main() -> None:
    t_start = time.time()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    GO_DIR.mkdir(parents=True, exist_ok=True)

    # ==================================================================
    # Load data
    # ==================================================================
    print("=" * 70)
    print("  PROGENITOR DIVERGENCE DEEP CHARACTERIZATION")
    print("=" * 70)

    print("\n[1/7] Loading data...")

    # Load annotations
    annot = pd.read_csv(ANNOTATIONS)
    print(f"  Annotations: {len(annot)} cell types")

    # Load Procrustes results
    with open(RESULTS_JSON) as f:
        proc = json.load(f)
    cell_types = proc["cell_types"]
    print(f"  Procrustes results: {len(cell_types)} cell types")

    # Load centroids
    centroids_h = pd.read_csv(CENTROIDS_HUMAN, index_col=0)
    centroids_m = pd.read_csv(CENTROIDS_MOUSE, index_col=0)
    gene_ids = centroids_h.columns.tolist()
    print(f"  Centroids: {centroids_h.shape[0]} types × {len(gene_ids)} genes")

    # Get gene symbol mapping from h5ad
    h5 = ad.read_h5ad(HUMAN_H5AD, backed="r")
    id_to_symbol = dict(zip(h5.var_names, h5.var["feature_name"]))
    h5.file.close()
    gene_symbols = [id_to_symbol.get(gid, gid) for gid in gene_ids]
    print(f"  Gene symbols mapped: {sum(1 for s in gene_symbols if not s.startswith('ENSG'))}/{len(gene_symbols)}")

    # ==================================================================
    # STEP 1: Identify the 5 progenitor types
    # ==================================================================
    print("\n" + "=" * 70)
    print("  STEP 1: Identify progenitor types")
    print("=" * 70)

    progenitors = annot[annot["progenitor"] == True]["cell_type"].tolist()
    differentiated = annot[annot["progenitor"] == False]["cell_type"].tolist()
    print(f"\n  Progenitor types ({len(progenitors)}):")
    for ct in progenitors:
        row = annot[annot["cell_type"] == ct].iloc[0]
        print(f"    {row['rank']:>2}. {ct:<45} residual={row['residual_magnitude']:.3f}  "
              f"({row['pct_of_ssr']:.1f}% SSR)")

    print(f"\n  Differentiated types ({len(differentiated)})")
    prog_residuals = annot[annot["progenitor"] == True]["residual_magnitude"]
    diff_residuals = annot[annot["progenitor"] == False]["residual_magnitude"]
    print(f"  Progenitor median residual: {prog_residuals.median():.3f}")
    print(f"  Differentiated median residual: {diff_residuals.median():.3f}")
    print(f"  Ratio: {prog_residuals.median() / diff_residuals.median():.2f}x")

    # ==================================================================
    # STEP 2: Compute progenitor specificity scores
    # ==================================================================
    print("\n" + "=" * 70)
    print("  STEP 2: Compute progenitor specificity scores")
    print("=" * 70)

    # Recompute PCA to get loadings (not stored in JSON)
    print("\n  Recomputing PCA from centroids to recover loadings...")
    ct_order = sorted(centroids_h.index.tolist())
    human_mat = centroids_h.loc[ct_order].values
    mouse_mat = centroids_m.loc[ct_order].values
    combined = np.vstack([human_mat, mouse_mat])

    pca = PCA(n_components=PCA_VARIANCE_THRESHOLD, svd_solver="full",
              random_state=RANDOM_SEED)
    combined_pca = pca.fit_transform(combined)
    n_comp = pca.n_components_
    print(f"  PCA: {n_comp} components, {pca.explained_variance_ratio_.sum()*100:.1f}% variance")

    # Get PCA loadings matrix W: (n_comp, n_genes)
    W = pca.components_

    # Project each cell type's PCA residual to gene space
    print("  Projecting residuals to gene space (all 16,959 genes)...")
    gene_residuals = {}  # ct -> (n_genes,) array
    for ct in ct_order:
        r_pca = np.array(proc["residuals"][ct]["vector_pca"])
        gene_residuals[ct] = r_pca @ W  # (n_genes,)

    # Compute mean absolute residual per gene for progenitors vs differentiated
    prog_set = set(progenitors)
    diff_set = set(differentiated)

    prog_abs = np.array([np.abs(gene_residuals[ct]) for ct in ct_order if ct in prog_set])
    diff_abs = np.array([np.abs(gene_residuals[ct]) for ct in ct_order if ct in diff_set])

    mean_prog = prog_abs.mean(axis=0)  # (n_genes,)
    mean_diff = diff_abs.mean(axis=0)  # (n_genes,)

    # Progenitor specificity score
    denom = mean_prog + mean_diff
    # Avoid division by zero for genes with zero residual everywhere
    specificity = np.where(denom > 0, mean_prog / denom, 0.5)

    # Build results DataFrame
    spec_df = pd.DataFrame({
        "gene_id": gene_ids,
        "gene_symbol": gene_symbols,
        "specificity_score": specificity,
        "mean_progenitor_residual": mean_prog,
        "mean_differentiated_residual": mean_diff,
    }).sort_values("specificity_score", ascending=False).reset_index(drop=True)

    spec_df.to_csv(OUTPUT_DIR / "progenitor_specificity_scores.csv", index=False)

    n_high = (specificity > 0.7).sum()
    n_prog_dominant = (specificity > 0.5).sum()
    print(f"\n  Specificity score distribution:")
    print(f"    > 0.7 (progenitor-specific): {n_high} genes")
    print(f"    > 0.5 (progenitor-dominant): {n_prog_dominant} genes")
    print(f"    ≤ 0.5 (differentiated-dominant): {len(gene_ids) - n_prog_dominant} genes")
    print(f"    Mean specificity: {specificity.mean():.4f}")
    print(f"    Median specificity: {np.median(specificity):.4f}")

    top100_genes = spec_df.head(100)["gene_symbol"].tolist()
    top50_genes = spec_df.head(50)
    print(f"\n  Top 10 progenitor-specific genes:")
    for _, row in spec_df.head(10).iterrows():
        print(f"    {row['gene_symbol']:<15} specificity={row['specificity_score']:.4f}  "
              f"prog={row['mean_progenitor_residual']:.4f}  "
              f"diff={row['mean_differentiated_residual']:.4f}")

    # ==================================================================
    # STEP 3: GO enrichment on progenitor-specific divergence genes
    # ==================================================================
    print("\n" + "=" * 70)
    print("  STEP 3: GO enrichment on progenitor-specific divergence genes")
    print("=" * 70)

    go_library = "GO_Biological_Process_2023"

    # Shared top 100
    print(f"\n  Running GO enrichment on top 100 progenitor-specific genes...")
    try:
        enr_shared = gseapy.enrich(
            gene_list=top100_genes,
            gene_sets=go_library,
            outdir=None,
            cutoff=1.0,
            no_plot=True,
            verbose=False,
        )
        shared_results = enr_shared.results
        if shared_results is not None and not shared_results.empty:
            shared_results = shared_results.sort_values("Adjusted P-value").reset_index(drop=True)
            shared_results.to_csv(GO_DIR / "enrichment_shared_progenitor.csv", index=False)
            n_sig = (shared_results["Adjusted P-value"] < 0.05).sum()
            print(f"    → {len(shared_results)} terms tested, {n_sig} significant (p_adj < 0.05)")
            print(f"    Top 5 terms:")
            for _, row in shared_results.head(5).iterrows():
                term = row["Term"].split(" (GO:")[0][:60]
                print(f"      {term}: p_adj={row['Adjusted P-value']:.2e}, overlap={row['Overlap']}")
        else:
            shared_results = pd.DataFrame()
            print("    → No enrichment results returned")
    except Exception as e:
        print(f"    → GO enrichment failed: {e}")
        shared_results = pd.DataFrame()

    # Per-progenitor type: top 100 genes by absolute residual magnitude
    print(f"\n  Running GO enrichment per progenitor type...")
    per_prog_go = {}
    for ct in progenitors:
        # Get top 100 genes for this specific progenitor type
        ct_abs_residuals = np.abs(gene_residuals[ct])
        top_idx = np.argsort(ct_abs_residuals)[::-1][:100]
        ct_top_genes = [gene_symbols[i] for i in top_idx]

        print(f"\n    {ct} (top 100 divergence genes)...")
        try:
            enr = gseapy.enrich(
                gene_list=ct_top_genes,
                gene_sets=go_library,
                outdir=None,
                cutoff=1.0,
                no_plot=True,
                verbose=False,
            )
            df = enr.results
            if df is not None and not df.empty:
                df = df.sort_values("Adjusted P-value").reset_index(drop=True)
                safe_name = ct.replace(" ", "_").replace(",", "")
                df.to_csv(GO_DIR / f"enrichment_{safe_name}.csv", index=False)
                n_sig = (df["Adjusted P-value"] < 0.05).sum()
                per_prog_go[ct] = df
                print(f"      → {n_sig} significant terms")
                if n_sig > 0:
                    top_term = df.iloc[0]["Term"].split(" (GO:")[0][:60]
                    print(f"      Top: {top_term} (p_adj={df.iloc[0]['Adjusted P-value']:.2e})")
            else:
                per_prog_go[ct] = pd.DataFrame()
                print(f"      → No results")
        except Exception as e:
            per_prog_go[ct] = pd.DataFrame()
            print(f"      → Failed: {e}")

    # ==================================================================
    # STEP 4: Functional categorization of top 50 genes
    # ==================================================================
    print("\n" + "=" * 70)
    print("  STEP 4: Functional categorization of top 50 progenitor genes")
    print("=" * 70)

    top50_genes["category"] = top50_genes["gene_symbol"].apply(categorize_gene)
    top50_genes.to_csv(OUTPUT_DIR / "top50_progenitor_genes.csv", index=False)

    print(f"\n  {'Gene':<15} {'Specificity':>11} {'Prog Resid':>11} {'Diff Resid':>11} {'Category':<20}")
    print(f"  {'-'*70}")
    for _, row in top50_genes.iterrows():
        print(f"  {row['gene_symbol']:<15} {row['specificity_score']:>11.4f} "
              f"{row['mean_progenitor_residual']:>11.4f} "
              f"{row['mean_differentiated_residual']:>11.4f} "
              f"{row['category']:<20}")

    # Category summary
    cat_counts = top50_genes["category"].value_counts()
    print(f"\n  Category breakdown (top 50 genes):")
    for cat, count in cat_counts.items():
        pct = count / len(top50_genes) * 100
        print(f"    {cat:<25} {count:>3} ({pct:.0f}%)")

    # ==================================================================
    # STEP 5: Human vs mouse expression heatmap
    # ==================================================================
    print("\n" + "=" * 70)
    print("  STEP 5: Human vs mouse expression comparison")
    print("=" * 70)

    top50_symbols = top50_genes["gene_symbol"].tolist()
    top50_ids = top50_genes["gene_id"].tolist()

    # Extract expression values for progenitor types
    human_expr = np.zeros((len(progenitors), 50))
    mouse_expr = np.zeros((len(progenitors), 50))
    for i, ct in enumerate(progenitors):
        for j, gid in enumerate(top50_ids):
            if gid in centroids_h.columns:
                human_expr[i, j] = centroids_h.loc[ct, gid]
                mouse_expr[i, j] = centroids_m.loc[ct, gid]

    # Compute log2 fold change (human/mouse), with pseudocount
    pseudo = 0.01
    log2fc = np.log2((human_expr + pseudo) / (mouse_expr + pseudo))

    # Short cell type names for plotting
    short_names = {
        "hematopoietic precursor cell": "Hemato. precursor",
        "hematopoietic stem cell": "HSC",
        "basal cell": "Basal cell",
        "mesenchymal stem cell of adipose tissue": "MSC (adipose)",
        "mesenchymal stem cell": "MSC",
    }
    prog_short = [short_names.get(ct, ct) for ct in progenitors]

    # Plot 1: Side-by-side heatmaps (human and mouse expression)
    fig, axes = plt.subplots(1, 3, figsize=(24, 8),
                              gridspec_kw={"width_ratios": [1, 1, 0.8]})

    # Truncate gene names for display
    display_genes = [s[:12] for s in top50_symbols]

    vmin = 0
    vmax = max(human_expr.max(), mouse_expr.max())

    sns.heatmap(human_expr, ax=axes[0], cmap="YlOrRd", vmin=vmin, vmax=vmax,
                xticklabels=display_genes, yticklabels=prog_short,
                linewidths=0.3, linecolor="white",
                cbar_kws={"label": "Mean expression (log-normalized)", "shrink": 0.6})
    axes[0].set_title("Human Progenitors", fontsize=13)
    axes[0].tick_params(axis="x", labelsize=6, rotation=90)
    axes[0].tick_params(axis="y", labelsize=9)

    sns.heatmap(mouse_expr, ax=axes[1], cmap="YlOrRd", vmin=vmin, vmax=vmax,
                xticklabels=display_genes, yticklabels=prog_short,
                linewidths=0.3, linecolor="white",
                cbar_kws={"label": "Mean expression (log-normalized)", "shrink": 0.6})
    axes[1].set_title("Mouse Progenitors", fontsize=13)
    axes[1].tick_params(axis="x", labelsize=6, rotation=90)
    axes[1].tick_params(axis="y", labelsize=9)

    # Log2FC heatmap
    maxabs = np.abs(log2fc).max()
    sns.heatmap(log2fc, ax=axes[2], cmap="RdBu_r", vmin=-maxabs, vmax=maxabs,
                xticklabels=display_genes, yticklabels=prog_short,
                linewidths=0.3, linecolor="white",
                cbar_kws={"label": "log₂(human/mouse)", "shrink": 0.6})
    axes[2].set_title("log₂ Fold Change", fontsize=13)
    axes[2].tick_params(axis="x", labelsize=6, rotation=90)
    axes[2].tick_params(axis="y", labelsize=9)

    fig.suptitle("Top 50 Progenitor-Specific Divergence Genes: Expression Comparison",
                 fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "expression_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # Summary stats
    human_higher = (log2fc > 0.5).sum()
    mouse_higher = (log2fc < -0.5).sum()
    similar = log2fc.size - human_higher - mouse_higher
    print(f"\n  Expression comparison (top 50 genes × 5 progenitor types):")
    print(f"    Human higher (log2FC > 0.5):  {human_higher} entries ({human_higher/log2fc.size*100:.0f}%)")
    print(f"    Mouse higher (log2FC < -0.5): {mouse_higher} entries ({mouse_higher/log2fc.size*100:.0f}%)")
    print(f"    Similar (|log2FC| ≤ 0.5):     {similar} entries ({similar/log2fc.size*100:.0f}%)")
    print(f"    Mean log2FC: {log2fc.mean():.3f}")
    print(f"  Heatmap saved: {OUTPUT_DIR / 'expression_heatmap.png'}")

    # ==================================================================
    # STEP 6: Pairwise overlap between progenitor types
    # ==================================================================
    print("\n" + "=" * 70)
    print("  STEP 6: Pairwise gene overlap between progenitor types")
    print("=" * 70)

    # Get top 100 divergence genes per progenitor type
    prog_top100 = {}
    for ct in progenitors:
        ct_abs = np.abs(gene_residuals[ct])
        top_idx = np.argsort(ct_abs)[::-1][:100]
        prog_top100[ct] = set(gene_symbols[i] for i in top_idx)

    # Compute Jaccard index matrix
    n_prog = len(progenitors)
    jaccard_mat = np.zeros((n_prog, n_prog))
    overlap_mat = np.zeros((n_prog, n_prog))
    for i, ct1 in enumerate(progenitors):
        for j, ct2 in enumerate(progenitors):
            inter = len(prog_top100[ct1] & prog_top100[ct2])
            union = len(prog_top100[ct1] | prog_top100[ct2])
            jaccard_mat[i, j] = inter / union if union > 0 else 0
            overlap_mat[i, j] = inter

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Jaccard heatmap
    sns.heatmap(jaccard_mat, ax=axes[0], cmap="YlGnBu",
                xticklabels=prog_short, yticklabels=prog_short,
                annot=True, fmt=".2f", linewidths=0.5, linecolor="white",
                vmin=0, vmax=1,
                cbar_kws={"label": "Jaccard Index", "shrink": 0.8})
    axes[0].set_title("Jaccard Similarity of Top 100 Divergence Genes", fontsize=12)
    axes[0].tick_params(axis="x", rotation=45, labelsize=9)
    axes[0].tick_params(axis="y", rotation=0, labelsize=9)
    plt.setp(axes[0].get_xticklabels(), ha="right", rotation_mode="anchor")

    # Raw overlap count heatmap
    sns.heatmap(overlap_mat, ax=axes[1], cmap="YlGnBu",
                xticklabels=prog_short, yticklabels=prog_short,
                annot=True, fmt=".0f", linewidths=0.5, linecolor="white",
                cbar_kws={"label": "Shared Genes (out of 100)", "shrink": 0.8})
    axes[1].set_title("Raw Overlap Count (Top 100 Genes)", fontsize=12)
    axes[1].tick_params(axis="x", rotation=45, labelsize=9)
    axes[1].tick_params(axis="y", rotation=0, labelsize=9)
    plt.setp(axes[1].get_xticklabels(), ha="right", rotation_mode="anchor")

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "pairwise_overlap_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # Compute off-diagonal statistics
    off_diag_jaccard = []
    off_diag_overlap = []
    for i in range(n_prog):
        for j in range(i + 1, n_prog):
            off_diag_jaccard.append(jaccard_mat[i, j])
            off_diag_overlap.append(overlap_mat[i, j])

    mean_jaccard = np.mean(off_diag_jaccard)
    mean_overlap = np.mean(off_diag_overlap)

    print(f"\n  Pairwise overlap statistics:")
    print(f"    Mean Jaccard (off-diagonal): {mean_jaccard:.3f}")
    print(f"    Mean overlap count: {mean_overlap:.1f}/100")
    print(f"    Range: {min(off_diag_overlap):.0f} - {max(off_diag_overlap):.0f}")

    # Interpretation
    if mean_jaccard > 0.3:
        overlap_interp = "HIGH — progenitor types share a common divergence program"
    elif mean_jaccard > 0.15:
        overlap_interp = "MODERATE — partial sharing of divergence genes"
    else:
        overlap_interp = "LOW — each progenitor type diverges independently"
    print(f"    Interpretation: {overlap_interp}")
    print(f"  Heatmap saved: {OUTPUT_DIR / 'pairwise_overlap_heatmap.png'}")

    # Find genes shared by ≥3 progenitor types
    gene_prog_count = {}
    for ct in progenitors:
        for g in prog_top100[ct]:
            gene_prog_count[g] = gene_prog_count.get(g, 0) + 1
    shared_3plus = {g: c for g, c in gene_prog_count.items() if c >= 3}
    shared_4plus = {g: c for g, c in gene_prog_count.items() if c >= 4}
    shared_all5 = {g: c for g, c in gene_prog_count.items() if c == 5}

    print(f"\n  Genes in top 100 for ≥3 progenitor types: {len(shared_3plus)}")
    print(f"  Genes in top 100 for ≥4 progenitor types: {len(shared_4plus)}")
    print(f"  Genes in top 100 for all 5 progenitor types: {len(shared_all5)}")
    if shared_all5:
        print(f"    All-5 genes: {', '.join(sorted(shared_all5.keys()))}")
    if shared_4plus:
        print(f"    ≥4 genes: {', '.join(sorted(shared_4plus.keys()))}")

    # ==================================================================
    # STEP 7: Comparison with most-diverged differentiated type
    # ==================================================================
    print("\n" + "=" * 70)
    print("  STEP 7: Progenitor vs differentiated divergence gene comparison")
    print("=" * 70)

    # Find the most-diverged differentiated type
    diff_annot = annot[annot["progenitor"] == False].sort_values("residual_magnitude", ascending=False)
    top_diff_ct = diff_annot.iloc[0]["cell_type"]
    top_diff_resid = diff_annot.iloc[0]["residual_magnitude"]
    print(f"\n  Most-diverged differentiated type: {top_diff_ct} (residual={top_diff_resid:.3f})")

    # Get top 100 genes for that type
    ct_abs = np.abs(gene_residuals[top_diff_ct])
    top_idx = np.argsort(ct_abs)[::-1][:100]
    diff_top100 = set(gene_symbols[i] for i in top_idx)

    # Compare with progenitor-specific top 100
    prog_shared_top100 = set(top100_genes)
    overlap_prog_diff = prog_shared_top100 & diff_top100
    jaccard_prog_diff = len(overlap_prog_diff) / len(prog_shared_top100 | diff_top100)

    print(f"  Progenitor-specific top 100 vs {top_diff_ct} top 100:")
    print(f"    Overlap: {len(overlap_prog_diff)} genes")
    print(f"    Jaccard: {jaccard_prog_diff:.3f}")
    if overlap_prog_diff:
        print(f"    Shared genes: {', '.join(sorted(overlap_prog_diff))}")

    # Also compare with 2nd and 3rd most-diverged differentiated types
    print(f"\n  Overlap with other highly-diverged differentiated types:")
    for idx in range(min(5, len(diff_annot))):
        ct = diff_annot.iloc[idx]["cell_type"]
        ct_abs_r = np.abs(gene_residuals[ct])
        top_idx_r = np.argsort(ct_abs_r)[::-1][:100]
        diff_genes_r = set(gene_symbols[i] for i in top_idx_r)
        ov = len(prog_shared_top100 & diff_genes_r)
        jac = len(prog_shared_top100 & diff_genes_r) / len(prog_shared_top100 | diff_genes_r)
        print(f"    {ct:<45} overlap={ov:>3}  Jaccard={jac:.3f}")

    # Save comparison
    comp_df = pd.DataFrame({
        "gene": sorted(prog_shared_top100 | diff_top100),
        "in_progenitor_top100": [g in prog_shared_top100 for g in sorted(prog_shared_top100 | diff_top100)],
        "in_diff_top100": [g in diff_top100 for g in sorted(prog_shared_top100 | diff_top100)],
        "in_both": [g in overlap_prog_diff for g in sorted(prog_shared_top100 | diff_top100)],
    })
    comp_df.to_csv(OUTPUT_DIR / "vs_differentiated_comparison.csv", index=False)

    # ==================================================================
    # SUMMARY
    # ==================================================================
    print("\n" + "=" * 70)
    print("  SUMMARY")
    print("=" * 70)

    elapsed = time.time() - t_start

    summary = {
        "progenitor_types": progenitors,
        "n_progenitor": len(progenitors),
        "n_differentiated": len(differentiated),
        "progenitor_median_residual": float(prog_residuals.median()),
        "differentiated_median_residual": float(diff_residuals.median()),
        "residual_ratio": float(prog_residuals.median() / diff_residuals.median()),
        "n_genes_specificity_gt_0.7": int(n_high),
        "n_genes_specificity_gt_0.5": int(n_prog_dominant),
        "top10_progenitor_genes": spec_df.head(10)["gene_symbol"].tolist(),
        "go_enrichment_shared": {
            "n_terms_tested": len(shared_results) if not shared_results.empty else 0,
            "n_significant": int((shared_results["Adjusted P-value"] < 0.05).sum()) if not shared_results.empty else 0,
            "top_term": shared_results.iloc[0]["Term"].split(" (GO:")[0] if not shared_results.empty else None,
            "top_p_adj": float(shared_results.iloc[0]["Adjusted P-value"]) if not shared_results.empty else None,
        },
        "pairwise_overlap": {
            "mean_jaccard": float(mean_jaccard),
            "mean_overlap_count": float(mean_overlap),
            "interpretation": overlap_interp,
            "genes_in_3plus_progenitors": len(shared_3plus),
            "genes_in_all_5_progenitors": len(shared_all5),
        },
        "vs_top_differentiated": {
            "differentiated_type": top_diff_ct,
            "overlap_count": len(overlap_prog_diff),
            "jaccard": float(jaccard_prog_diff),
        },
        "expression_direction": {
            "human_higher_pct": float(human_higher / log2fc.size * 100),
            "mouse_higher_pct": float(mouse_higher / log2fc.size * 100),
            "similar_pct": float(similar / log2fc.size * 100),
            "mean_log2fc": float(log2fc.mean()),
        },
        "category_breakdown": {cat: int(count) for cat, count in cat_counts.items()},
        "runtime_seconds": round(elapsed, 1),
    }

    with open(OUTPUT_DIR / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"""
  PROGENITOR DIVERGENCE DEEP CHARACTERIZATION — RESULTS
  ======================================================

  1. PROGENITOR TYPES (n={len(progenitors)}):
     Median residual: {prog_residuals.median():.3f} vs differentiated {diff_residuals.median():.3f}
     ({prog_residuals.median()/diff_residuals.median():.2f}x higher, p=0.0099)

  2. PROGENITOR-SPECIFIC GENES:
     {n_high} genes with specificity score > 0.7
     Top 5: {', '.join(spec_df.head(5)['gene_symbol'].tolist())}

  3. GO ENRICHMENT (shared progenitor genes):
     {summary['go_enrichment_shared']['n_significant']} significant GO terms
     Top: {summary['go_enrichment_shared']['top_term']} (p_adj={summary['go_enrichment_shared']['top_p_adj']:.2e})

  4. FUNCTIONAL CATEGORIES (top 50 genes):
     {'; '.join(f'{cat}: {count}' for cat, count in cat_counts.head(5).items())}

  5. EXPRESSION DIRECTION:
     Human higher: {human_higher/log2fc.size*100:.0f}%
     Mouse higher: {mouse_higher/log2fc.size*100:.0f}%
     Similar: {similar/log2fc.size*100:.0f}%
     → Mean log2FC = {log2fc.mean():.3f}

  6. SHARED vs INDEPENDENT DIVERGENCE:
     Mean Jaccard (pairwise): {mean_jaccard:.3f}
     {overlap_interp}
     Genes shared by ≥3 progenitors: {len(shared_3plus)}
     Genes shared by all 5: {len(shared_all5)}

  7. PROGENITOR vs DIFFERENTIATED COMPARISON:
     Overlap with {top_diff_ct}: {len(overlap_prog_diff)}/100 genes
     Jaccard: {jaccard_prog_diff:.3f}
     → Progenitor divergence is {'qualitatively similar to' if jaccard_prog_diff > 0.2 else 'qualitatively distinct from'} differentiated divergence

  Runtime: {elapsed:.1f}s
  Output: {OUTPUT_DIR}/
""")


if __name__ == "__main__":
    main()
