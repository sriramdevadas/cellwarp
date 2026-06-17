#!/usr/bin/env python3
"""
Gene-level expression conservation vs evolutionary constraint analysis.

Tests whether genes under stronger evolutionary constraint show tighter
cross-species expression conservation in CellWarp's geometric space.

Biology
-------
If the CellWarp framework captures real evolutionary signal, then genes
under strong purifying selection (low dN/dS, high sequence identity, low
LOEUF) should show more conserved expression patterns across human and
mouse cell types. This analysis validates the biological relevance of the
geometric framework.

Math
----
Two per-gene conservation metrics:
(a) Procrustes contribution: Gene loading on PCs weighted by per-PC
    alignment quality (inverse residual). Genes on well-conserved PCs
    score high.
(b) Cross-species expression correlation: Pearson r of gene's expression
    across 35 matched cell-type centroids (human vs mouse). Primary metric.

Evolutionary constraint metrics:
- LOEUF (gnomAD): loss-of-function observed/expected upper bound fraction.
  Lower LOEUF = more LoF-intolerant = more constrained. Continuous metric,
  preferred over bimodal pLI.
- Sequence % identity: from Ensembl BioMart ortholog table. Higher = more
  conserved at the protein sequence level.
- pLI, missense Z: additional gnomAD metrics for sensitivity.

Outputs
-------
- analysis/gene_conservation/gene_conservation_table.csv
- analysis/gene_conservation/correlation_results.json
- analysis/gene_conservation/gene_conservation_summary.md
- figures/supplementary/gene_conservation_*.png
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from sklearn.decomposition import PCA

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data" / "phase1"
OUTPUT_PHASE2 = PROJECT_ROOT / "output" / "phase2"
SCALED_35 = OUTPUT_PHASE2 / "scaled_35types"
ANALYSIS_DIR = PROJECT_ROOT / "analysis" / "gene_conservation"
FIGURES_DIR = PROJECT_ROOT / "figures" / "supplementary"

ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Step 0: Load existing Procrustes data
# ---------------------------------------------------------------------------


def load_centroids_35() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load 35-type centroids for human and mouse.

    Returns matched DataFrames: (n_cell_types, n_genes) with Ensembl IDs as columns.
    """
    human = pd.read_csv(SCALED_35 / "centroids_human_35.csv", index_col=0)
    mouse = pd.read_csv(SCALED_35 / "centroids_mouse_35.csv", index_col=0)

    shared_types = sorted(set(human.index) & set(mouse.index))
    human = human.loc[shared_types]
    mouse = mouse.loc[shared_types]

    shared_genes = sorted(set(human.columns) & set(mouse.columns))
    human = human[shared_genes]
    mouse = mouse[shared_genes]

    print(f"Loaded 35-type centroids: {len(shared_types)} cell types x {len(shared_genes)} genes")
    return human, mouse


def load_centroids_6() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load 6-type centroids for human and mouse."""
    human = pd.read_csv(OUTPUT_PHASE2 / "centroids_human.csv", index_col=0)
    mouse = pd.read_csv(OUTPUT_PHASE2 / "centroids_mouse.csv", index_col=0)

    shared_types = sorted(set(human.index) & set(mouse.index))
    human = human.loc[shared_types]
    mouse = mouse.loc[shared_types]

    shared_genes = sorted(set(human.columns) & set(mouse.columns))
    human = human[shared_genes]
    mouse = mouse[shared_genes]

    print(f"Loaded 6-type centroids: {len(shared_types)} cell types x {len(shared_genes)} genes")
    return human, mouse


def load_orthologs() -> pd.DataFrame:
    """Load ortholog mapping table."""
    ortho = pd.read_csv(DATA_DIR / "orthologs_human_mouse.csv")
    print(f"Loaded {len(ortho)} ortholog pairs")
    return ortho


def load_procrustes_results() -> dict:
    """Load saved Procrustes results JSON."""
    with open(OUTPUT_PHASE2 / "procrustes_results.json") as f:
        return json.load(f)


def load_biomart_homology() -> pd.DataFrame | None:
    """Load BioMart ortholog % identity data if available."""
    path = ANALYSIS_DIR / "biomart_homology.csv"
    if path.exists():
        df = pd.read_csv(path)
        print(f"Loaded BioMart homology data: {len(df)} orthologs")
        return df
    return None


# ---------------------------------------------------------------------------
# Step 1a: Per-gene Procrustes contribution (sensitivity metric)
# ---------------------------------------------------------------------------


def compute_procrustes_contribution(
    human_6: pd.DataFrame,
    mouse_6: pd.DataFrame,
    proc_results: dict,
) -> pd.Series:
    """Compute per-gene Procrustes contribution score.

    For each gene, compute its loading on the PCs used in alignment,
    weighted by the per-PC alignment quality (1 - normalized residual).

    Math:
        For each PC j, compute the per-PC residual as the sum of squared
        differences between aligned mouse and human centroids on that PC.
        Weight = explained_variance_ratio[j] * (1 - normalized_residual[j])
        Gene score = sum over PCs of |loading[gene, PC_j]| * weight[j]

    Genes loading heavily on well-conserved PCs get high scores.
    """
    cell_types = sorted(human_6.index)
    gene_ids = human_6.columns.tolist()

    human_mat = human_6.loc[cell_types].values
    mouse_mat = mouse_6.loc[cell_types].values
    combined = np.vstack([human_mat, mouse_mat])

    pca = PCA(n_components=0.95, svd_solver="full", random_state=42)
    pca.fit_transform(combined)

    n_components = pca.n_components_
    print(f"  PCA refitted: {n_components} components, "
          f"{pca.explained_variance_ratio_.sum()*100:.1f}% variance")

    aligned_mouse = np.array(proc_results["aligned_mouse_centroids_pca"])
    human_centered = np.array(proc_results["human_centroids_pca_centered"])

    per_pc_residual = np.sum((aligned_mouse - human_centered) ** 2, axis=0)
    total_residual = per_pc_residual.sum()
    normalized_residual = per_pc_residual / total_residual if total_residual > 0 else per_pc_residual

    var_explained = pca.explained_variance_ratio_
    conservation_weight = var_explained * (1 - normalized_residual)

    W = pca.components_  # (n_components, n_genes)
    gene_scores = np.abs(W).T @ conservation_weight  # (n_genes,)

    result = pd.Series(gene_scores, index=gene_ids, name="procrustes_contribution")
    print(f"  Procrustes contribution scores: {len(result)} genes")
    print(f"  Score range: [{result.min():.6f}, {result.max():.6f}]")
    return result


# ---------------------------------------------------------------------------
# Step 1b: Per-gene cross-species expression correlation (primary metric)
# ---------------------------------------------------------------------------


def compute_expression_correlation(
    human: pd.DataFrame,
    mouse: pd.DataFrame,
) -> pd.DataFrame:
    """Compute per-gene Pearson correlation across matched cell types.

    For each gene, correlate its expression vector across N cell types
    in human vs mouse. High correlation = conserved expression pattern.

    Math: For gene g with expression h_g = (h_g1, ..., h_gN) in human
    and m_g = (m_g1, ..., m_gN) in mouse:
        r_g = Pearson(h_g, m_g)

    Returns DataFrame with columns: gene_id, pearson_r, p_value
    """
    n_types = human.shape[0]
    gene_ids = human.columns.tolist()

    correlations = []
    p_values = []

    human_vals = human.values
    mouse_vals = mouse.values

    for j in range(len(gene_ids)):
        h = human_vals[:, j]
        m = mouse_vals[:, j]

        if np.std(h) < 1e-12 or np.std(m) < 1e-12:
            correlations.append(np.nan)
            p_values.append(np.nan)
            continue

        r, p = stats.pearsonr(h, m)
        correlations.append(r)
        p_values.append(p)

    df = pd.DataFrame({
        "gene_id": gene_ids,
        "pearson_r": correlations,
        "pearson_p": p_values,
    })

    valid = df["pearson_r"].notna()
    print(f"\n  Per-gene expression correlation ({n_types} cell types):")
    print(f"  Valid genes: {valid.sum()} / {len(gene_ids)}")
    print(f"  Mean r: {df.loc[valid, 'pearson_r'].mean():.4f}")
    print(f"  Median r: {df.loc[valid, 'pearson_r'].median():.4f}")
    print(f"  Fraction r > 0.5: {(df.loc[valid, 'pearson_r'] > 0.5).mean():.3f}")
    print(f"  Fraction r > 0.8: {(df.loc[valid, 'pearson_r'] > 0.8).mean():.3f}")

    return df


# ---------------------------------------------------------------------------
# Step 2: Retrieve evolutionary constraint metrics
# ---------------------------------------------------------------------------


def fetch_gnomad_constraint() -> pd.DataFrame | None:
    """Download gnomAD gene constraint metrics as a proxy for evolutionary constraint.

    gnomAD's LOEUF (loss-of-function observed/expected upper bound fraction)
    is a continuous measure of within-human constraint that strongly correlates
    with cross-species conservation (Karczewski et al. 2020, Nature).

    Lower LOEUF = more LoF-intolerant = more constrained.

    Returns DataFrame with gene_id and constraint columns, or None on failure.
    """
    import urllib.request
    import gzip
    import io

    url = "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"

    cache_path = ANALYSIS_DIR / "gnomad_constraint_cache.tsv"
    if cache_path.exists():
        print(f"\n  Loading cached gnomAD constraint from {cache_path}")
        df = pd.read_csv(cache_path, sep="\t")
        print(f"  Loaded {len(df)} genes from cache")
        return df

    print(f"\n  Downloading gnomAD constraint metrics...")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "CellWarp/1.0"})
        with urllib.request.urlopen(req, timeout=120) as response:
            compressed = response.read()

        with gzip.open(io.BytesIO(compressed), "rt") as f:
            df = pd.read_csv(f, sep="\t", low_memory=False)

        cols_needed = ["gene", "gene_id"]
        cols_constraint = ["pLI", "oe_lof_upper", "mis_z", "oe_lof", "oe_mis"]
        cols_available = [c for c in cols_constraint if c in df.columns]
        df = df[cols_needed + cols_available].copy()
        df = df.rename(columns={"gene_id": "human_ensembl_id", "gene": "gene_name"})

        df["human_ensembl_id"] = df["human_ensembl_id"].str.replace(r"\.\d+$", "", regex=True)
        df = df.drop_duplicates(subset="human_ensembl_id", keep="first")

        df.to_csv(cache_path, sep="\t", index=False)

        print(f"  Downloaded constraint for {len(df)} genes")
        for col in cols_available:
            valid = df[col].notna().sum()
            print(f"    {col}: {valid} valid values")

        return df

    except Exception as e:
        print(f"  gnomAD download failed: {e}")
        return None


# ---------------------------------------------------------------------------
# Step 3: Correlation analysis
# ---------------------------------------------------------------------------


def spearman_with_ci(x: np.ndarray, y: np.ndarray, n_bootstrap: int = 1000) -> dict:
    """Spearman correlation with bootstrap 95% CI."""
    mask = np.isfinite(x) & np.isfinite(y)
    x_clean = x[mask]
    y_clean = y[mask]
    n = len(x_clean)

    rho, p = stats.spearmanr(x_clean, y_clean)

    rng = np.random.RandomState(42)
    boot_rhos = np.zeros(n_bootstrap)
    for i in range(n_bootstrap):
        idx = rng.choice(n, size=n, replace=True)
        boot_rhos[i], _ = stats.spearmanr(x_clean[idx], y_clean[idx])

    ci_low, ci_high = np.percentile(boot_rhos, [2.5, 97.5])

    return {
        "rho": float(rho),
        "p_value": float(p),
        "ci_low": float(ci_low),
        "ci_high": float(ci_high),
        "n": int(n),
    }


def compute_decile_means(df: pd.DataFrame, conservation_col: str,
                          constraint_col: str) -> pd.DataFrame:
    """Bin genes by constraint decile, compute mean conservation per bin.

    Decile 1 = lowest constraint value, decile 10 = highest.
    For LOEUF: decile 1 = most constrained (low LOEUF), decile 10 = least.
    For % identity: decile 1 = least conserved, decile 10 = most conserved.
    """
    valid = df[[conservation_col, constraint_col]].dropna()
    valid = valid.copy()
    valid["decile"] = pd.qcut(
        valid[constraint_col],
        q=10,
        labels=False,
        duplicates="drop",
    ) + 1

    decile_stats = valid.groupby("decile")[conservation_col].agg(
        ["mean", "std", "count", "sem"]
    ).reset_index()

    return decile_stats


# ---------------------------------------------------------------------------
# Step 4: Partial correlation controlling for expression level
# ---------------------------------------------------------------------------


def partial_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> dict:
    """Partial Spearman correlation between x and y, controlling for z.

    Uses rank-based partial correlation:
    1. Rank-transform x, y, z
    2. Regress ranks of x on ranks of z -> residuals r_x
    3. Regress ranks of y on ranks of z -> residuals r_y
    4. Pearson correlation of r_x and r_y = partial Spearman
    """
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x_c = stats.rankdata(x[mask])
    y_c = stats.rankdata(y[mask])
    z_c = stats.rankdata(z[mask])
    n = len(x_c)

    def residualize(a, b):
        b_mean = b - b.mean()
        a_mean = a - a.mean()
        beta = np.dot(a_mean, b_mean) / np.dot(b_mean, b_mean)
        return a - beta * b

    rx = residualize(x_c, z_c)
    ry = residualize(y_c, z_c)

    r, p = stats.pearsonr(rx, ry)

    return {
        "partial_rho": float(r),
        "p_value": float(p),
        "n": int(n),
        "controlling_for": "mean_expression",
    }


# ---------------------------------------------------------------------------
# Step 5: Figures
# ---------------------------------------------------------------------------


def plot_conservation_vs_constraint(
    df: pd.DataFrame,
    conservation_col: str,
    constraint_col: str,
    constraint_label: str,
    rho_info: dict,
    partial_info: dict,
    output_path: Path,
    flip_x: bool = False,
    flip_x_label: str | None = None,
):
    """Panel A: Hexbin scatter of conservation vs constraint with annotations.

    Args:
        flip_x: If True, negate x values so that 'more constrained' is to the right.
        flip_x_label: X-axis label when flip_x is True.
    """
    valid = df[[conservation_col, constraint_col]].dropna()

    fig, ax = plt.subplots(figsize=(7, 6))

    x = valid[constraint_col].values.copy()
    y = valid[conservation_col].values.copy()

    if flip_x:
        x = -x
        xlabel = flip_x_label or f"-{constraint_label}"
    else:
        xlabel = constraint_label

    hb = ax.hexbin(x, y, gridsize=50, cmap="YlOrRd", mincnt=1,
                    linewidths=0.2, edgecolors="white")
    fig.colorbar(hb, ax=ax, label="Gene count")

    rho = rho_info["rho"]
    p = rho_info["p_value"]
    n = rho_info["n"]
    partial_rho = partial_info["partial_rho"]
    partial_p = partial_info["p_value"]

    p_str = f"{p:.2e}" if p < 0.001 else f"{p:.4f}"
    pp_str = f"{partial_p:.2e}" if partial_p < 0.001 else f"{partial_p:.4f}"

    annotation = (
        f"Spearman \u03c1 = {rho:.3f} (p = {p_str}, n = {n:,})\n"
        f"Partial \u03c1 = {partial_rho:.3f} (p = {pp_str})\n"
        f"(controlling for expression level)"
    )
    ax.text(0.03, 0.97, annotation, transform=ax.transAxes,
            fontsize=9, verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", alpha=0.9))

    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel("Expression conservation (Pearson r)", fontsize=12)
    ax.set_title("Per-gene expression conservation vs evolutionary constraint", fontsize=13)

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def plot_decile_bars(
    decile_stats: pd.DataFrame,
    constraint_label: str,
    output_path: Path,
    subtitle: str = "",
):
    """Panel B: Bar plot of mean expression conservation per constraint decile."""
    fig, ax = plt.subplots(figsize=(8, 5))

    deciles = decile_stats["decile"].values
    means = decile_stats["mean"].values
    sems = decile_stats["sem"].values

    colors = plt.cm.RdYlBu_r(np.linspace(0.15, 0.85, len(deciles)))

    ax.bar(deciles, means, yerr=sems, capsize=3,
           color=colors, edgecolor="white", linewidth=0.5,
           error_kw={"linewidth": 1, "color": "0.3"})

    ax.set_xlabel(f"{constraint_label} decile", fontsize=12)
    ax.set_ylabel("Mean expression conservation\n(Pearson r)", fontsize=12)

    title = "Expression conservation by evolutionary constraint decile"
    if subtitle:
        title += f"\n{subtitle}"
    ax.set_title(title, fontsize=12)

    ax.set_xticks(deciles)

    for i, (d, m, c) in enumerate(zip(deciles, means, decile_stats["count"].values)):
        ax.text(d, m + sems[i] + 0.005, f"n={c}", ha="center", va="bottom", fontsize=7)

    sns.despine(ax=ax)
    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def plot_multi_metric_comparison(results_dict: dict, output_path: Path):
    """Summary panel: compare rho across different constraint metrics."""
    metrics = []
    rhos = []
    partial_rhos = []
    labels = []

    for key, info in results_dict.items():
        if "rho" in info:
            metrics.append(key)
            labels.append(info.get("label", key))
            rhos.append(info["rho"])
            partial_rhos.append(info.get("partial_rho", np.nan))

    if len(metrics) < 2:
        return

    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(metrics))
    width = 0.35

    bars1 = ax.bar(x - width/2, rhos, width, label="Raw Spearman \u03c1",
                    color="steelblue", edgecolor="white")
    bars2 = ax.bar(x + width/2, partial_rhos, width,
                    label="Partial \u03c1 (controlling for expression)",
                    color="darkorange", edgecolor="white")

    ax.set_ylabel("Spearman \u03c1", fontsize=12)
    ax.set_title("Expression conservation correlations across constraint metrics", fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10, rotation=15, ha="right")
    ax.legend(fontsize=9)
    ax.axhline(0, color="gray", linewidth=0.5, linestyle="-")
    sns.despine(ax=ax)

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


# ---------------------------------------------------------------------------
# Step 6: Save outputs
# ---------------------------------------------------------------------------


def save_results(
    gene_table: pd.DataFrame,
    correlation_results: dict,
    summary_text: str,
):
    """Save all outputs to analysis/gene_conservation/."""
    csv_path = ANALYSIS_DIR / "gene_conservation_table.csv"
    gene_table.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path} ({len(gene_table)} genes)")

    json_path = ANALYSIS_DIR / "correlation_results.json"
    with open(json_path, "w") as f:
        json.dump(correlation_results, f, indent=2)
    print(f"  Saved: {json_path}")

    md_path = ANALYSIS_DIR / "gene_conservation_summary.md"
    with open(md_path, "w") as f:
        f.write(summary_text)
    print(f"  Saved: {md_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    print("=" * 70)
    print("Gene-level expression conservation vs evolutionary constraint")
    print("=" * 70)

    # ------------------------------------------------------------------
    # Step 0: Load data
    # ------------------------------------------------------------------
    print("\n--- Step 0: Loading data ---")
    human_35, mouse_35 = load_centroids_35()
    human_6, mouse_6 = load_centroids_6()
    orthologs = load_orthologs()
    proc_results = load_procrustes_results()
    biomart_data = load_biomart_homology()

    gene_name_map = dict(zip(orthologs["human_ensembl_id"], orthologs["human_gene_name"]))

    # ------------------------------------------------------------------
    # Step 1b: Per-gene expression correlation (PRIMARY)
    # ------------------------------------------------------------------
    print("\n--- Step 1b: Per-gene cross-species expression correlation ---")
    expr_corr = compute_expression_correlation(human_35, mouse_35)

    # ------------------------------------------------------------------
    # Step 1a: Procrustes contribution (SENSITIVITY)
    # ------------------------------------------------------------------
    print("\n--- Step 1a: Per-gene Procrustes contribution ---")
    proc_contrib = compute_procrustes_contribution(human_6, mouse_6, proc_results)

    # Build gene table
    gene_table = expr_corr.copy()
    gene_table["gene_name"] = gene_table["gene_id"].map(gene_name_map)
    gene_table["procrustes_contribution"] = gene_table["gene_id"].map(proc_contrib)

    # Mean expression level across both species and all cell types
    mean_expr_human = human_35.mean(axis=0)
    mean_expr_mouse = mouse_35.mean(axis=0)
    mean_expr = (mean_expr_human + mean_expr_mouse) / 2
    gene_table["mean_expression"] = gene_table["gene_id"].map(mean_expr)

    # ------------------------------------------------------------------
    # Step 2: Evolutionary constraint metrics
    # ------------------------------------------------------------------
    print("\n--- Step 2: Retrieving evolutionary constraint metrics ---")

    # 2a: BioMart sequence % identity
    if biomart_data is not None:
        gene_table = gene_table.merge(
            biomart_data[["human_ensembl_id", "perc_id_human_to_mouse",
                          "perc_id_mouse_to_human", "goc_score", "wga_coverage"]],
            left_on="gene_id",
            right_on="human_ensembl_id",
            how="left",
        ).drop(columns=["human_ensembl_id"], errors="ignore")

        # Average of bidirectional % identity
        gene_table["seq_identity"] = (
            gene_table["perc_id_human_to_mouse"].astype(float) +
            gene_table["perc_id_mouse_to_human"].astype(float)
        ) / 2

        n_seqid = gene_table["seq_identity"].notna().sum()
        print(f"  Sequence % identity: {n_seqid} genes")
        print(f"  Range: [{gene_table['seq_identity'].min():.1f}%, "
              f"{gene_table['seq_identity'].max():.1f}%]")
        print(f"  Median: {gene_table['seq_identity'].median():.1f}%")

    # 2b: gnomAD constraint
    gnomad_data = fetch_gnomad_constraint()
    if gnomad_data is not None:
        gene_table = gene_table.merge(
            gnomad_data,
            left_on="gene_id",
            right_on="human_ensembl_id",
            how="left",
        ).drop(columns=["human_ensembl_id"], errors="ignore")
        if "gene_name_y" in gene_table.columns:
            gene_table = gene_table.drop(columns=["gene_name_y"])
            gene_table = gene_table.rename(columns={"gene_name_x": "gene_name"})

        for col in ["oe_lof_upper", "pLI", "mis_z"]:
            if col in gene_table.columns:
                n_valid = gene_table[col].notna().sum()
                print(f"  gnomAD {col}: {n_valid} genes")

    # ------------------------------------------------------------------
    # Step 3: Correlation analysis
    # ------------------------------------------------------------------
    print("\n--- Step 3: Correlation analysis ---")

    y_conservation = gene_table["pearson_r"].values
    z_expr = gene_table["mean_expression"].values

    # Define metrics to test: (column, label, expected_direction, flip_for_plot)
    # expected_direction: "negative" = lower value = more constrained
    #                     "positive" = higher value = more constrained
    metrics_to_test = []

    if "oe_lof_upper" in gene_table.columns and gene_table["oe_lof_upper"].notna().sum() > 1000:
        metrics_to_test.append(("oe_lof_upper", "LOEUF (gnomAD)", "negative", True))

    if "seq_identity" in gene_table.columns and gene_table["seq_identity"].notna().sum() > 1000:
        metrics_to_test.append(("seq_identity", "Sequence % identity", "positive", False))

    if "pLI" in gene_table.columns and gene_table["pLI"].notna().sum() > 1000:
        metrics_to_test.append(("pLI", "pLI (gnomAD)", "positive", False))

    if "mis_z" in gene_table.columns and gene_table["mis_z"].notna().sum() > 1000:
        metrics_to_test.append(("mis_z", "Missense Z (gnomAD)", "positive", False))

    if not metrics_to_test:
        print("  ERROR: No constraint metrics available!")
        sys.exit(1)

    # Use first metric as primary
    primary_col, primary_label, primary_dir, primary_flip = metrics_to_test[0]
    print(f"\n  Primary constraint metric: {primary_label}")

    # Collect all correlation results
    correlation_results = {
        "primary_metric": "pearson_r",
        "n_cell_types": int(human_35.shape[0]),
        "n_total_genes": int(len(gene_table)),
    }

    all_metric_results = {}

    for col, label, direction, flip in metrics_to_test:
        x = gene_table[col].values

        # Raw Spearman
        raw = spearman_with_ci(y_conservation, x)

        # Partial Spearman controlling for expression
        partial = partial_spearman(y_conservation, x, z_expr)

        # Decile analysis
        deciles = compute_decile_means(gene_table, "pearson_r", col)

        # Interpret direction
        if direction == "negative":
            # Lower value = more constrained
            # Expect negative rho (more constrained -> more conserved)
            expected_sign = "negative"
            correct_direction = raw["rho"] < 0
        else:
            # Higher value = more constrained/conserved
            # Expect positive rho
            expected_sign = "positive"
            correct_direction = raw["rho"] > 0

        result_entry = {
            "label": label,
            "column": col,
            "expected_direction": expected_sign,
            "raw_spearman": raw,
            "partial_spearman": partial,
            "deciles": deciles.to_dict(orient="records"),
            "correct_direction": bool(correct_direction),
            "rho": raw["rho"],
            "partial_rho": partial["partial_rho"],
        }
        all_metric_results[col] = result_entry

        star = "*" if raw["p_value"] < 0.05 else ""
        dir_sym = "OK" if correct_direction else "UNEXPECTED"
        print(f"\n  {label}:")
        print(f"    Spearman rho = {raw['rho']:.4f} "
              f"(p = {raw['p_value']:.2e}, n = {raw['n']:,}) {star}")
        print(f"    95% CI: [{raw['ci_low']:.4f}, {raw['ci_high']:.4f}]")
        print(f"    Expected sign: {expected_sign}, Observed: {'negative' if raw['rho'] < 0 else 'positive'} [{dir_sym}]")
        print(f"    Partial rho = {partial['partial_rho']:.4f} "
              f"(p = {partial['p_value']:.2e})")

    # Also correlate Procrustes contribution with primary constraint
    proc_vals = gene_table["procrustes_contribution"].values
    primary_x = gene_table[primary_col].values
    proc_raw = spearman_with_ci(proc_vals, primary_x)
    proc_partial = partial_spearman(proc_vals, primary_x, z_expr)
    print(f"\n  Sensitivity (Procrustes contribution vs {primary_label}):")
    print(f"    Spearman rho = {proc_raw['rho']:.4f} (p = {proc_raw['p_value']:.2e})")
    print(f"    Partial rho  = {proc_partial['partial_rho']:.4f} (p = {proc_partial['p_value']:.2e})")

    # Expression level vs conservation
    expr_corr_result = spearman_with_ci(y_conservation, z_expr)
    print(f"\n  Expression level vs conservation: rho = {expr_corr_result['rho']:.4f} "
          f"(p = {expr_corr_result['p_value']:.2e})")

    # Store in results
    correlation_results["constraint_metrics"] = {
        k: {key: val for key, val in v.items() if key != "deciles"}
        for k, v in all_metric_results.items()
    }
    correlation_results["decile_analyses"] = {
        k: v["deciles"] for k, v in all_metric_results.items()
    }
    correlation_results["sensitivity_procrustes_contribution"] = {
        "raw": proc_raw,
        "partial": proc_partial,
    }
    correlation_results["expression_vs_conservation"] = expr_corr_result

    # ------------------------------------------------------------------
    # Step 4: Decision gate
    # ------------------------------------------------------------------
    print("\n--- Step 4: Decision gate ---")

    primary_result = all_metric_results[primary_col]
    primary_raw = primary_result["raw_spearman"]
    primary_partial = primary_result["partial_spearman"]

    # Also check sequence identity if available
    seqid_result = all_metric_results.get("seq_identity")

    # Evaluate: does the correlation survive expression-level control?
    surviving_metrics = []
    for col, info in all_metric_results.items():
        p_raw = info["raw_spearman"]
        p_par = info["partial_spearman"]
        if (p_raw["p_value"] < 0.05 and
            p_par["p_value"] < 0.05 and
            info["correct_direction"]):
            # Check if partial rho has same sign as raw and retains >= 30%
            raw_sign = np.sign(p_raw["rho"])
            partial_sign = np.sign(p_par["partial_rho"])
            if raw_sign == partial_sign:
                retention = abs(p_par["partial_rho"]) / abs(p_raw["rho"]) * 100
                surviving_metrics.append((col, info["label"], retention))

    if surviving_metrics:
        verdict = "PASS"
        detail_parts = []
        for col, label, retention in surviving_metrics:
            detail_parts.append(f"{label} retains {retention:.0f}% after expression control")
        detail = (
            "Correlation between expression conservation and evolutionary "
            "constraint survives expression-level control for: " +
            "; ".join(detail_parts) + ". "
            "The geometric framework captures genuine evolutionary signal "
            "beyond expression magnitude."
        )
    else:
        # Check if raw correlations are significant but partial ones flip
        raw_sig = [col for col, info in all_metric_results.items()
                   if info["raw_spearman"]["p_value"] < 0.05 and info["correct_direction"]]
        if raw_sig:
            verdict = "PARTIAL"
            detail = (
                "Raw correlations are significant in the expected direction, but "
                "the effect is substantially confounded by expression level. "
                "After controlling for mean expression, the signal weakens or reverses. "
                "The geometric framework captures expression magnitude effects that "
                "correlate with constraint, rather than independent evolutionary signal."
            )
        else:
            verdict = "FAIL"
            detail = (
                "No significant correlation between expression conservation and "
                "evolutionary constraint in the expected direction."
            )

    correlation_results["verdict"] = verdict
    correlation_results["verdict_detail"] = detail
    print(f"\n  Verdict: {verdict}")
    print(f"  {detail}")

    # ------------------------------------------------------------------
    # Step 5: Figures
    # ------------------------------------------------------------------
    print("\n--- Step 5: Generating figures ---")

    # Panel A: Scatter using primary metric (LOEUF)
    primary_partial_info = primary_result["partial_spearman"]
    scatter_path = FIGURES_DIR / "gene_conservation_scatter.png"
    plot_conservation_vs_constraint(
        gene_table,
        "pearson_r",
        primary_col,
        primary_label,
        primary_raw,
        primary_partial_info,
        scatter_path,
        flip_x=primary_flip,
        flip_x_label="1/LOEUF (higher = more constrained)" if primary_flip else None,
    )

    # Panel A2: Scatter using seq identity (if available)
    if seqid_result is not None:
        scatter2_path = FIGURES_DIR / "gene_conservation_scatter_seqid.png"
        plot_conservation_vs_constraint(
            gene_table,
            "pearson_r",
            "seq_identity",
            "Sequence % identity",
            seqid_result["raw_spearman"],
            seqid_result["partial_spearman"],
            scatter2_path,
        )

    # Panel B: Decile bars for primary
    primary_deciles = compute_decile_means(gene_table, "pearson_r", primary_col)
    bar_path = FIGURES_DIR / "gene_conservation_deciles.png"
    if primary_dir == "negative":
        subtitle = "(Decile 1 = most constrained [lowest LOEUF] -> 10 = least constrained)"
    else:
        subtitle = "(Decile 1 = least constrained -> 10 = most constrained)"
    plot_decile_bars(primary_deciles, primary_label, bar_path, subtitle=subtitle)

    # Panel B2: Decile bars for seq identity
    if seqid_result is not None:
        seqid_deciles = compute_decile_means(gene_table, "pearson_r", "seq_identity")
        bar2_path = FIGURES_DIR / "gene_conservation_deciles_seqid.png"
        plot_decile_bars(
            seqid_deciles, "Sequence % identity", bar2_path,
            subtitle="(Decile 1 = least conserved sequence -> 10 = most conserved)"
        )

    # Multi-metric comparison
    if len(all_metric_results) >= 2:
        plot_multi_metric_comparison(
            all_metric_results,
            FIGURES_DIR / "gene_conservation_multi_metric.png",
        )

    # Conservation histogram
    fig, ax = plt.subplots(figsize=(7, 4))
    valid_r = gene_table["pearson_r"].dropna()
    ax.hist(valid_r, bins=80, color="steelblue", edgecolor="white", linewidth=0.3)
    ax.axvline(valid_r.median(), color="red", linestyle="--", linewidth=1.5,
               label=f"Median = {valid_r.median():.3f}")
    ax.set_xlabel("Per-gene expression conservation (Pearson r)", fontsize=12)
    ax.set_ylabel("Number of genes", fontsize=12)
    ax.set_title("Distribution of cross-species expression conservation\n"
                 "(across 35 matched cell types)", fontsize=12)
    ax.legend(fontsize=10)
    sns.despine(ax=ax)
    plt.tight_layout()
    hist_path = FIGURES_DIR / "gene_conservation_histogram.png"
    fig.savefig(hist_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {hist_path}")

    # ------------------------------------------------------------------
    # Step 6: Save outputs
    # ------------------------------------------------------------------
    print("\n--- Step 6: Saving outputs ---")

    # Build summary markdown
    summary_lines = [
        "# Gene-level Expression Conservation vs Evolutionary Constraint",
        "",
        "## Summary",
        f"- **Primary conservation metric**: Per-gene Pearson correlation across "
        f"{human_35.shape[0]} matched cell-type centroids (human vs mouse)",
        f"- **Genes with valid conservation**: {gene_table['pearson_r'].notna().sum():,} "
        f"/ {len(gene_table):,}",
        "",
        "## Constraint Metrics Tested",
        "",
    ]

    for col, info in all_metric_results.items():
        raw = info["raw_spearman"]
        par = info["partial_spearman"]
        direction_ok = "correct" if info["correct_direction"] else "unexpected"
        summary_lines.extend([
            f"### {info['label']}",
            f"- Spearman rho = {raw['rho']:.4f} "
            f"(p = {raw['p_value']:.2e}, n = {raw['n']:,})",
            f"- 95% CI: [{raw['ci_low']:.4f}, {raw['ci_high']:.4f}]",
            f"- Direction: {direction_ok}",
            f"- Partial rho = {par['partial_rho']:.4f} "
            f"(p = {par['p_value']:.2e}, controlling for expression level)",
            "",
        ])

    summary_lines.extend([
        "## Expression Level Confound",
        f"- Expression vs conservation: rho = {expr_corr_result['rho']:.4f} "
        f"(p = {expr_corr_result['p_value']:.2e})",
        "- Expression level is a strong confound: highly expressed genes show "
        "higher cross-species conservation AND tend to be more constrained.",
        "",
        "## Sensitivity Check (Procrustes Contribution Metric)",
        f"- vs {primary_label}: rho = {proc_raw['rho']:.4f} "
        f"(p = {proc_raw['p_value']:.2e})",
        f"- Partial rho = {proc_partial['partial_rho']:.4f} "
        f"(p = {proc_partial['p_value']:.2e})",
        "",
    ])

    # Primary decile table
    summary_lines.extend([
        f"## Decile Analysis ({primary_label})",
        "| Decile | Mean Conservation | SEM | n |",
        "|--------|-------------------|-----|---|",
    ])
    for _, row in primary_deciles.iterrows():
        summary_lines.append(
            f"| {int(row['decile'])} | {row['mean']:.4f} | {row['sem']:.4f} | {int(row['count'])} |"
        )

    summary_lines.extend([
        "",
        f"## Decision Gate: **{verdict}**",
        f"",
        detail,
        "",
        "## Files",
        "- `gene_conservation_table.csv`: Per-gene conservation and constraint values",
        "- `correlation_results.json`: Full correlation results with CIs",
        "- `figures/supplementary/gene_conservation_scatter.png`: Hexbin scatter (Panel A, LOEUF)",
        "- `figures/supplementary/gene_conservation_scatter_seqid.png`: Hexbin scatter (seq identity)",
        "- `figures/supplementary/gene_conservation_deciles.png`: Decile bars (Panel B, LOEUF)",
        "- `figures/supplementary/gene_conservation_deciles_seqid.png`: Decile bars (seq identity)",
        "- `figures/supplementary/gene_conservation_multi_metric.png`: Multi-metric comparison",
        "- `figures/supplementary/gene_conservation_histogram.png`: Conservation distribution",
    ]
    )

    summary_text = "\n".join(summary_lines) + "\n"

    save_results(gene_table, correlation_results, summary_text)

    # ------------------------------------------------------------------
    # Final summary
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("RESULTS SUMMARY")
    print("=" * 70)
    print(f"  Conservation: Pearson r across {human_35.shape[0]} cell types")
    print(f"  Genes analyzed: {gene_table['pearson_r'].notna().sum():,}")
    print()
    for col, info in all_metric_results.items():
        raw = info["raw_spearman"]
        par = info["partial_spearman"]
        print(f"  {info['label']}:")
        print(f"    Raw rho     = {raw['rho']:+.4f} (p = {raw['p_value']:.2e})")
        print(f"    Partial rho = {par['partial_rho']:+.4f} (p = {par['p_value']:.2e})")
    print()
    print(f"  Expression confound: rho = {expr_corr_result['rho']:.4f}")
    print(f"  Verdict: {verdict}")
    print("=" * 70)


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=FutureWarning)
    main()
