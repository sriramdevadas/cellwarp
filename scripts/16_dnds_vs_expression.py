#!/usr/bin/env python3
"""
CellWarp — Sequence Conservation vs Expression Divergence

Adversarial validation: tests whether Procrustes gene loadings (expression-level
conservation) correlate with protein sequence conservation (percent identity) for
human-mouse gene pairs.

Biology
-------
Protein percent identity measures amino acid conservation between orthologous
proteins. Higher identity = more conserved sequence. If Procrustes expression
divergence tracks sequence conservation, the geometric signal is downstream of
sequence evolution. If they are decoupled, Procrustes geometry captures
regulatory evolution that sequence conservation misses — a stronger result.

Note: dN/dS ratios were originally planned but are no longer available in
Ensembl Compara (confirmed NULL in releases 109–115). Protein percent identity
is the standard alternative, directly measuring amino acid conservation. It
correlates strongly with dN/dS where both are available (high identity ≈ low
dN/dS = purifying selection). We use (100 - %identity) as "sequence divergence"
for intuitive comparison with expression divergence.

Math
----
Procrustes gene loading magnitude per gene = |g_j| = |r_i · W_j| where r_i
is the residual vector in PCA space and W_j is the j-th column of the PCA
loadings matrix. Averaged across all 35 cell types, this gives each gene's
mean contribution to cross-species expression divergence.

Sequence divergence = 100 - protein percent identity (higher = more diverged).

Spearman correlation tests monotonic association between sequence divergence
and expression divergence without assuming linearity.

Pipeline
--------
1. Load protein percent identity from Ensembl Compara 115 FTP download
2. Reconstruct full gene-space Procrustes loadings from centroids + residuals
3. Genome-wide Spearman correlation
4. Per-cell-type Spearman correlations
5. Extremes analysis: decoupled genes + GO enrichment
6. Summary interpretation

Inputs:
    output/phase2/scaled_35types/procrustes_results_35.json
    output/phase2/scaled_35types/centroids_human_35.csv
    output/phase2/scaled_35types/centroids_mouse_35.csv
    output/phase2/dnds/ensembl_homology_mouse_human.tsv  (Ensembl 115 FTP)
    data/phase1/orthologs_human_mouse.csv
    output/phase2/developmental_constraint/developmental_annotations.csv

Outputs (all in output/phase2/dnds/):
    dnds_raw.csv                      — processed sequence conservation data
    gene_loadings_all.csv             — per-gene Procrustes loadings (all genes × 35 types)
    genome_wide_scatter.png           — seq divergence vs mean expression divergence
    per_celltype_correlations.png     — bar chart of per-cell-type ρ values
    per_celltype_correlations.csv     — per-cell-type ρ and p-values
    extremes_high_seq_high_expr.csv   — conserved sequence, divergent expression
    extremes_low_seq_low_expr.csv     — divergent sequence, conserved expression
    go_enrichment_decoupled.csv       — GO enrichment on decoupled gene sets
    dnds_summary.json                 — all key statistics

Usage:
    python scripts/16_dnds_vs_expression.py
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
import seaborn as sns
from scipy import stats
from sklearn.decomposition import PCA

import gseapy

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

OUTPUT_DIR = Path("./output/phase2/dnds")
RESULTS_JSON = Path("./output/phase2/scaled_35types/procrustes_results_35.json")
CENTROIDS_HUMAN = Path("./output/phase2/scaled_35types/centroids_human_35.csv")
CENTROIDS_MOUSE = Path("./output/phase2/scaled_35types/centroids_mouse_35.csv")
ORTHOLOGS_CSV = Path("./data/phase1/orthologs_human_mouse.csv")
ANNOTATIONS_CSV = Path("./output/phase2/developmental_constraint/developmental_annotations.csv")
HOMOLOGY_TSV = Path("./output/phase2/dnds/ensembl_homology_mouse_human.tsv")

PCA_VARIANCE_THRESHOLD = 0.95
RANDOM_SEED = 42


# ---------------------------------------------------------------------------
# Step 1: Load sequence conservation data
# ---------------------------------------------------------------------------

def step1_load_sequence_conservation(gene_ids: list[str]) -> pd.DataFrame:
    """
    Load protein percent identity from Ensembl Compara 115 FTP download.

    The file contains mouse → human 1:1 orthologs with:
    - gene_stable_id: mouse Ensembl gene ID
    - homology_gene_stable_id: human Ensembl gene ID
    - identity: % aa identity of mouse protein to human
    - homology_identity: % aa identity of human protein to mouse

    We compute sequence divergence = 100 - mean(identity, homology_identity)
    as a measure of protein sequence divergence (higher = more diverged).

    Args:
        gene_ids: List of human Ensembl gene IDs from our gene space.

    Returns:
        DataFrame with columns: human_ensembl_id, mouse_ensembl_id,
        pct_identity, seq_divergence.
    """
    print("\n" + "=" * 70)
    print("STEP 1: Load sequence conservation from Ensembl Compara 115")
    print("=" * 70)

    df = pd.read_csv(HOMOLOGY_TSV, sep="\t")
    print(f"  Loaded {len(df):,} mouse-human 1:1 orthologs from Ensembl 115")

    # Rename for clarity: mouse IDs are gene_stable_id, human are homology_gene_stable_id
    df = df.rename(columns={
        "gene_stable_id": "mouse_ensembl_id",
        "homology_gene_stable_id": "human_ensembl_id",
        "identity": "mouse_to_human_identity",
        "homology_identity": "human_to_mouse_identity",
    })

    # Convert to numeric
    df["mouse_to_human_identity"] = pd.to_numeric(df["mouse_to_human_identity"], errors="coerce")
    df["human_to_mouse_identity"] = pd.to_numeric(df["human_to_mouse_identity"], errors="coerce")

    # Mean percent identity (average of bidirectional)
    df["pct_identity"] = (df["mouse_to_human_identity"] + df["human_to_mouse_identity"]) / 2

    # Sequence divergence = 100 - pct_identity (higher = more diverged)
    df["seq_divergence"] = 100 - df["pct_identity"]

    # Filter to our gene space
    gene_set = set(gene_ids)
    df = df[df["human_ensembl_id"].isin(gene_set)].copy()

    # Drop duplicates (keep first)
    df = df.drop_duplicates(subset="human_ensembl_id", keep="first")

    # Handle edge cases
    n_null = df["pct_identity"].isna().sum()
    n_perfect = (df["pct_identity"] == 100).sum()

    print(f"  Matched to our gene space: {len(df):,} / {len(gene_ids):,}")
    print(f"  Null identity: {n_null:,}")
    print(f"  Perfect identity (100%): {n_perfect:,}")
    print(f"  Mean identity: {df['pct_identity'].mean():.1f}%")
    print(f"  Median identity: {df['pct_identity'].median():.1f}%")
    print(f"  Range: [{df['pct_identity'].min():.1f}%, {df['pct_identity'].max():.1f}%]")
    print(f"  Mean seq divergence: {df['seq_divergence'].mean():.1f}%")

    # Drop NaN
    df = df.dropna(subset=["pct_identity"])
    print(f"  After dropping NaN: {len(df):,} genes")

    # Save
    df[["human_ensembl_id", "mouse_ensembl_id", "pct_identity", "seq_divergence"]].to_csv(
        OUTPUT_DIR / "dnds_raw.csv", index=False
    )
    print(f"  Saved: {OUTPUT_DIR / 'dnds_raw.csv'}")

    return df


# ---------------------------------------------------------------------------
# Step 2: Compute per-gene Procrustes loading magnitudes
# ---------------------------------------------------------------------------

def step2_gene_loadings() -> tuple[pd.DataFrame, list[str], list[str]]:
    """
    Reconstruct full gene-space Procrustes loadings for all 16,959 genes
    across all 35 cell types.

    Re-fits PCA on centroids (same parameters as original pipeline), then
    projects each cell type's residual vector from PCA space back to gene
    space using the PCA components matrix.

    Returns:
        Tuple of (loadings_df, gene_ensembl_ids, cell_types) where:
        - loadings_df: DataFrame (n_genes × 35+1) with per-cell-type absolute
          loadings and mean_expression_divergence column.
        - gene_ensembl_ids: List of Ensembl IDs in gene order.
        - cell_types: List of cell type names.
    """
    print("\n" + "=" * 70)
    print("STEP 2: Compute per-gene Procrustes loading magnitudes")
    print("=" * 70)

    # Load centroids
    h_cent = pd.read_csv(CENTROIDS_HUMAN, index_col=0)
    m_cent = pd.read_csv(CENTROIDS_MOUSE, index_col=0)
    print(f"  Human centroids: {h_cent.shape}")
    print(f"  Mouse centroids: {m_cent.shape}")

    gene_ensembl_ids = h_cent.columns.tolist()

    # Load Procrustes results for residual vectors
    with open(RESULTS_JSON) as f:
        results = json.load(f)

    cell_types = results["cell_types"]

    # Re-fit PCA (same as original pipeline)
    cell_types_sorted = sorted(cell_types)
    human_mat = h_cent.loc[cell_types_sorted].values
    mouse_mat = m_cent.loc[cell_types_sorted].values
    combined = np.vstack([human_mat, mouse_mat])

    pca = PCA(
        n_components=PCA_VARIANCE_THRESHOLD,
        svd_solver="full",
        random_state=RANDOM_SEED,
    )
    pca.fit(combined)
    print(f"  PCA: {pca.n_components_} components, {sum(pca.explained_variance_ratio_)*100:.1f}% variance")

    W = pca.components_  # (k, n_genes)

    # Get gene names from orthologs for display
    orthologs = pd.read_csv(ORTHOLOGS_CSV)
    id_to_name = dict(zip(orthologs["human_ensembl_id"], orthologs["human_gene_name"]))

    # Project each cell type's residual to gene space
    loading_data = {"ensembl_gene_id": gene_ensembl_ids}
    loading_data["gene_name"] = [id_to_name.get(g, g) for g in gene_ensembl_ids]

    for ct in cell_types:
        residual_pca = np.array(results["residuals"][ct]["vector_pca"])
        gene_loadings = residual_pca @ W  # (n_genes,)
        loading_data[ct] = np.abs(gene_loadings)

    loadings_df = pd.DataFrame(loading_data)

    # Compute mean expression divergence across all 35 cell types
    ct_cols = [c for c in loadings_df.columns if c not in ("ensembl_gene_id", "gene_name")]
    loadings_df["mean_expression_divergence"] = loadings_df[ct_cols].mean(axis=1)

    print(f"  Gene loadings computed: {len(gene_ensembl_ids):,} genes x {len(cell_types)} cell types")
    print(f"  Mean expression divergence range: [{loadings_df['mean_expression_divergence'].min():.6f}, "
          f"{loadings_df['mean_expression_divergence'].max():.6f}]")

    # Save
    loadings_df.to_csv(OUTPUT_DIR / "gene_loadings_all.csv", index=False)
    print(f"  Saved: {OUTPUT_DIR / 'gene_loadings_all.csv'}")

    return loadings_df, gene_ensembl_ids, cell_types


# ---------------------------------------------------------------------------
# Step 3: Genome-wide correlation
# ---------------------------------------------------------------------------

def step3_genome_wide(seq_df: pd.DataFrame, loadings_df: pd.DataFrame) -> dict:
    """
    Merge sequence divergence with mean expression divergence and compute
    Spearman correlation.

    Args:
        seq_df: DataFrame with sequence divergence per gene.
        loadings_df: DataFrame with per-gene loadings and mean_expression_divergence.

    Returns:
        Dict with correlation statistics.
    """
    print("\n" + "=" * 70)
    print("STEP 3: Genome-wide correlation (seq divergence vs expression divergence)")
    print("=" * 70)

    # Merge on human Ensembl ID
    merged = pd.merge(
        seq_df[["human_ensembl_id", "pct_identity", "seq_divergence"]],
        loadings_df[["ensembl_gene_id", "mean_expression_divergence"]],
        left_on="human_ensembl_id",
        right_on="ensembl_gene_id",
        how="inner",
    )
    print(f"  Merged genes: {len(merged):,}")

    # Spearman correlation: seq_divergence vs expression divergence
    rho, p_value = stats.spearmanr(merged["seq_divergence"], merged["mean_expression_divergence"])
    print(f"\n  Spearman correlation (seq_divergence vs expr_divergence):")
    print(f"    rho = {rho:.4f}")
    print(f"    p-value = {p_value:.2e}")
    print(f"    n = {len(merged):,}")

    # Also try pct_identity directly (negative expected if seq conservation → expr conservation)
    rho_inv, p_inv = stats.spearmanr(merged["pct_identity"], merged["mean_expression_divergence"])
    print(f"\n  Spearman (pct_identity vs expr_divergence) — expect negative:")
    print(f"    rho = {rho_inv:.4f}")
    print(f"    p-value = {p_inv:.2e}")

    # Pearson on seq_divergence (already linear scale)
    r_pearson, p_pearson = stats.pearsonr(
        merged["seq_divergence"], merged["mean_expression_divergence"]
    )
    print(f"\n  Pearson correlation:")
    print(f"    r = {r_pearson:.4f}")
    print(f"    p-value = {p_pearson:.2e}")

    # Plot: density scatter
    fig, ax = plt.subplots(figsize=(8, 6))

    x = merged["seq_divergence"].values
    y = merged["mean_expression_divergence"].values

    hb = ax.hexbin(x, y, gridsize=60, cmap="YlOrRd", mincnt=1)
    fig.colorbar(hb, ax=ax, label="Gene count")

    # Regression line
    slope, intercept, _, _, _ = stats.linregress(x, y)
    x_line = np.linspace(x.min(), x.max(), 100)
    ax.plot(x_line, slope * x_line + intercept, "k--", linewidth=1.5, alpha=0.7)

    ax.set_xlabel("Sequence divergence (100 - % protein identity)", fontsize=12)
    ax.set_ylabel("Mean expression divergence\n(Procrustes loading magnitude)", fontsize=11)
    ax.set_title(
        f"Sequence vs Expression Conservation\n"
        f"Spearman rho = {rho:.3f}, p = {p_value:.2e}, n = {len(merged):,}",
        fontsize=12,
    )

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "genome_wide_scatter.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Plot saved: {OUTPUT_DIR / 'genome_wide_scatter.png'}")
    print(f"  Description: Hexbin density scatter of sequence divergence (100-%identity) "
          f"vs mean Procrustes loading magnitude across {len(merged):,} genes. "
          f"Dashed line = linear regression.")

    result = {
        "spearman_rho": float(rho),
        "spearman_p": float(p_value),
        "spearman_rho_identity": float(rho_inv),
        "spearman_p_identity": float(p_inv),
        "pearson_r": float(r_pearson),
        "pearson_p": float(p_pearson),
        "n_genes": len(merged),
    }

    # Save merged data for downstream steps
    merged.to_csv(OUTPUT_DIR / "merged_seq_expr_divergence.csv", index=False)

    return result


# ---------------------------------------------------------------------------
# Step 4: Per-cell-type correlation
# ---------------------------------------------------------------------------

def step4_per_celltype(
    seq_df: pd.DataFrame,
    loadings_df: pd.DataFrame,
    cell_types: list[str],
) -> pd.DataFrame:
    """
    Compute Spearman correlation between sequence divergence and each cell
    type's gene loadings.

    Args:
        seq_df: DataFrame with seq_divergence per gene.
        loadings_df: DataFrame with per-gene, per-cell-type loadings.
        cell_types: List of cell type names.

    Returns:
        DataFrame with per-cell-type correlation statistics.
    """
    print("\n" + "=" * 70)
    print("STEP 4: Per-cell-type correlation")
    print("=" * 70)

    seq_usable = seq_df[["human_ensembl_id", "seq_divergence"]].copy()

    rows = []
    for ct in sorted(cell_types):
        ct_data = loadings_df[["ensembl_gene_id", ct]].copy()
        merged = pd.merge(
            seq_usable, ct_data,
            left_on="human_ensembl_id", right_on="ensembl_gene_id",
            how="inner",
        )
        rho, p_val = stats.spearmanr(merged["seq_divergence"], merged[ct])
        rows.append({
            "cell_type": ct,
            "spearman_rho": rho,
            "p_value": p_val,
            "n_genes": len(merged),
            "significant": p_val < 0.05,
        })

    corr_df = pd.DataFrame(rows).sort_values("spearman_rho", ascending=False)
    corr_df = corr_df.reset_index(drop=True)

    print(f"\n  {'Cell Type':<50} {'rho':>8} {'p-value':>12} {'Sig':>5}")
    print(f"  {'-' * 77}")
    for _, row in corr_df.iterrows():
        sig = "*" if row["significant"] else ""
        print(f"  {row['cell_type']:<50} {row['spearman_rho']:>8.4f} {row['p_value']:>12.2e} {sig:>5}")

    n_sig = corr_df["significant"].sum()
    print(f"\n  Significant (p < 0.05): {n_sig}/{len(cell_types)}")
    print(f"  Median rho: {corr_df['spearman_rho'].median():.4f}")
    print(f"  Range: [{corr_df['spearman_rho'].min():.4f}, {corr_df['spearman_rho'].max():.4f}]")

    # Load lineage annotations for coloring
    annotations = pd.read_csv(ANNOTATIONS_CSV)
    lineage_map = dict(zip(annotations["cell_type"], annotations["lineage"]))
    corr_df["lineage"] = corr_df["cell_type"].map(lineage_map)

    lineage_colors = {
        "hematopoietic": "#e74c3c",
        "epithelial": "#2ecc71",
        "mesenchymal": "#3498db",
        "endothelial": "#9b59b6",
    }

    # Plot: bar chart sorted by rho
    corr_sorted = corr_df.sort_values("spearman_rho")

    fig, ax = plt.subplots(figsize=(10, 10))
    colors = [lineage_colors.get(lin, "#95a5a6") for lin in corr_sorted["lineage"]]
    ax.barh(range(len(corr_sorted)), corr_sorted["spearman_rho"], color=colors, edgecolor="white")

    ax.set_yticks(range(len(corr_sorted)))
    ct_labels = [ct.replace("-positive, alpha-beta ", "+ ") for ct in corr_sorted["cell_type"]]
    ax.set_yticklabels(ct_labels, fontsize=8)
    ax.set_xlabel("Spearman rho (seq divergence vs gene loading)", fontsize=11)
    ax.set_title("Per-Cell-Type Sequence-Expression Conservation Coupling", fontsize=12)
    ax.axvline(0, color="black", linewidth=0.5)

    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=c, label=l) for l, c in lineage_colors.items()]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=9)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "per_celltype_correlations.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    corr_df.to_csv(OUTPUT_DIR / "per_celltype_correlations.csv", index=False)
    print(f"\n  Plot saved: {OUTPUT_DIR / 'per_celltype_correlations.png'}")
    print(f"  Description: Horizontal bar chart of Spearman rho (seq divergence vs loading) "
          f"for 35 cell types, colored by lineage category.")

    return corr_df


# ---------------------------------------------------------------------------
# Step 5: Extremes analysis
# ---------------------------------------------------------------------------

def step5_extremes(seq_df: pd.DataFrame, loadings_df: pd.DataFrame) -> dict:
    """
    Identify genes where sequence and expression evolution are decoupled.

    Two categories:
    - HIGH sequence conservation (identity > 95%) + HIGH expression divergence
    - LOW sequence conservation (identity < 70%) + LOW expression divergence

    Runs GO enrichment on each set.

    Returns:
        Dict with extremes statistics and GO results.
    """
    print("\n" + "=" * 70)
    print("STEP 5: Extremes analysis — decoupled genes")
    print("=" * 70)

    merged = pd.merge(
        seq_df[["human_ensembl_id", "pct_identity", "seq_divergence"]],
        loadings_df[["ensembl_gene_id", "gene_name", "mean_expression_divergence"]],
        left_on="human_ensembl_id",
        right_on="ensembl_gene_id",
        how="inner",
    )

    # Compute percentile thresholds for expression divergence
    expr_75 = merged["mean_expression_divergence"].quantile(0.75)
    expr_25 = merged["mean_expression_divergence"].quantile(0.25)

    print(f"  Expression divergence 75th percentile: {expr_75:.6f}")
    print(f"  Expression divergence 25th percentile: {expr_25:.6f}")
    print(f"  Identity distribution: median={merged['pct_identity'].median():.1f}%")

    # Category 1: HIGH sequence conservation + HIGH expression divergence
    # (protein is locked down but expression pattern has diverged)
    cat1 = merged[
        (merged["pct_identity"] > 95) &
        (merged["mean_expression_divergence"] > expr_75)
    ].sort_values("mean_expression_divergence", ascending=False)

    # Category 2: LOW sequence conservation + LOW expression divergence
    # (protein is evolving fast but expression pattern is conserved)
    cat2 = merged[
        (merged["pct_identity"] < 70) &
        (merged["mean_expression_divergence"] < expr_25)
    ].sort_values("pct_identity")

    print(f"\n  Category 1 — Conserved sequence, divergent expression "
          f"(identity > 95%, expr > 75th pctile):")
    print(f"    {len(cat1):,} genes")
    print(f"\n    Top 20:")
    for i, (_, row) in enumerate(cat1.head(20).iterrows()):
        print(f"      {i+1:>3}. {row['gene_name']:<15} identity={row['pct_identity']:.1f}%  "
              f"expr_div={row['mean_expression_divergence']:.6f}")

    print(f"\n  Category 2 — Divergent sequence, conserved expression "
          f"(identity < 70%, expr < 25th pctile):")
    print(f"    {len(cat2):,} genes")
    print(f"\n    Top 20:")
    for i, (_, row) in enumerate(cat2.head(20).iterrows()):
        print(f"      {i+1:>3}. {row['gene_name']:<15} identity={row['pct_identity']:.1f}%  "
              f"expr_div={row['mean_expression_divergence']:.6f}")

    # Save
    cat1.to_csv(OUTPUT_DIR / "extremes_high_seq_high_expr.csv", index=False)
    cat2.to_csv(OUTPUT_DIR / "extremes_low_seq_low_expr.csv", index=False)

    # GO enrichment on each set
    print("\n  Running GO enrichment on decoupled gene sets...")

    go_results = {}
    for label, df_sub, desc in [
        ("conserved_seq_divergent_expr", cat1,
         "conserved sequence + divergent expression"),
        ("divergent_seq_conserved_expr", cat2,
         "divergent sequence + conserved expression"),
    ]:
        genes = df_sub["gene_name"].dropna().tolist()
        if len(genes) < 5:
            print(f"\n    {desc}: too few genes ({len(genes)}) for GO enrichment")
            go_results[label] = {"n_genes": len(genes), "n_terms": 0, "top_terms": []}
            continue

        print(f"\n    {desc}: {len(genes)} genes")
        try:
            enr = gseapy.enrich(
                gene_list=genes,
                gene_sets="GO_Biological_Process_2023",
                outdir=None,
                cutoff=1.0,
                no_plot=True,
                verbose=False,
            )
            enr_df = enr.results
            if enr_df is not None and not enr_df.empty:
                enr_df = enr_df.sort_values("Adjusted P-value")
                n_sig = (enr_df["Adjusted P-value"] < 0.05).sum()
                print(f"      {n_sig} significant GO terms (p_adj < 0.05)")

                top_terms = []
                for _, row in enr_df.head(10).iterrows():
                    term_short = row["Term"].split(" (GO:")[0]
                    sig = "*" if row["Adjusted P-value"] < 0.05 else ""
                    print(f"        {term_short[:60]:<62} p_adj={row['Adjusted P-value']:.2e} {sig}")
                    top_terms.append({
                        "term": row["Term"],
                        "p_adj": float(row["Adjusted P-value"]),
                        "overlap": row["Overlap"],
                        "genes": row.get("Genes", ""),
                    })

                go_results[label] = {
                    "n_genes": len(genes),
                    "n_terms": int(n_sig),
                    "top_terms": top_terms,
                }
                enr_df.to_csv(OUTPUT_DIR / f"go_enrichment_{label}.csv", index=False)
            else:
                print(f"      No GO enrichment results")
                go_results[label] = {"n_genes": len(genes), "n_terms": 0, "top_terms": []}
        except Exception as e:
            print(f"      GO enrichment failed: {e}")
            go_results[label] = {"n_genes": len(genes), "n_terms": 0, "top_terms": [], "error": str(e)}

    # Plot extremes on the scatter
    fig, ax = plt.subplots(figsize=(9, 7))

    x = merged["seq_divergence"].values
    y = merged["mean_expression_divergence"].values
    ax.hexbin(x, y, gridsize=60, cmap="Greys", mincnt=1, alpha=0.5)

    if not cat1.empty:
        ax.scatter(
            cat1["seq_divergence"], cat1["mean_expression_divergence"],
            c="red", s=8, alpha=0.6,
            label=f"Conserved seq, divergent expr (n={len(cat1)})",
            zorder=5,
        )
    if not cat2.empty:
        ax.scatter(
            cat2["seq_divergence"], cat2["mean_expression_divergence"],
            c="blue", s=8, alpha=0.6,
            label=f"Divergent seq, conserved expr (n={len(cat2)})",
            zorder=5,
        )

    ax.set_xlabel("Sequence divergence (100 - % identity)", fontsize=12)
    ax.set_ylabel("Mean expression divergence", fontsize=11)
    ax.set_title("Decoupled Genes: Sequence vs Expression Conservation", fontsize=12)
    ax.legend(fontsize=9, loc="upper left")

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "extremes_scatter.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Extremes scatter saved: {OUTPUT_DIR / 'extremes_scatter.png'}")

    return {
        "cat1_conserved_seq_divergent_expr": {
            "n_genes": len(cat1),
            "threshold_identity": 95.0,
            "threshold_expr_pctile": 75,
            "threshold_expr_value": float(expr_75),
            "top_20_genes": cat1.head(20)["gene_name"].tolist(),
            "go": go_results.get("conserved_seq_divergent_expr", {}),
        },
        "cat2_divergent_seq_conserved_expr": {
            "n_genes": len(cat2),
            "threshold_identity": 70.0,
            "threshold_expr_pctile": 25,
            "threshold_expr_value": float(expr_25),
            "top_20_genes": cat2.head(20)["gene_name"].tolist(),
            "go": go_results.get("divergent_seq_conserved_expr", {}),
        },
    }


# ---------------------------------------------------------------------------
# Step 6: Summary
# ---------------------------------------------------------------------------

def step6_summary(
    genome_wide: dict,
    per_celltype: pd.DataFrame,
    extremes: dict,
    n_genes_total: int,
    n_seq_usable: int,
) -> None:
    """
    Print summary interpretation and save all results to JSON.
    """
    print("\n" + "=" * 70)
    print("STEP 6: Summary and Interpretation")
    print("=" * 70)

    rho = genome_wide["spearman_rho"]
    p = genome_wide["spearman_p"]

    # Determine outcome
    if abs(rho) < 0.1:
        outcome = "B"
        interpretation = (
            f"WEAK/NULL CORRELATION (rho = {rho:.3f}, p = {p:.2e}). "
            f"Expression divergence and sequence conservation are largely decoupled. "
            f"The Procrustes geometry captures regulatory evolution that sequence "
            f"conservation alone misses. This strengthens the case for geometric "
            f"morphometrics as a complementary lens on cross-species divergence."
        )
    elif rho > 0.1 and p < 0.05:
        outcome = "A"
        strength = "moderate" if rho < 0.3 else "strong"
        complement = "substantial" if rho < 0.3 else "some"
        interpretation = (
            f"POSITIVE CORRELATION (rho = {rho:.3f}, p = {p:.2e}). "
            f"Expression divergence partially tracks sequence divergence — a {strength} "
            f"coupling. The geometric signal is partly downstream of sequence evolution. "
            f"However, {complement} expression divergence is independent of sequence "
            f"conservation, indicating regulatory evolution contributes beyond protein coding."
        )
    elif rho < -0.1 and p < 0.05:
        outcome = "C"
        interpretation = (
            f"NEGATIVE CORRELATION (rho = {rho:.3f}, p = {p:.2e}). "
            f"Genes with more conserved protein sequence show MORE expression divergence. "
            f"This may reflect compensatory regulatory evolution — genes constrained at "
            f"the sequence level may diverge more in expression patterns across species."
        )
    else:
        outcome = "B"
        interpretation = (
            f"NON-SIGNIFICANT CORRELATION (rho = {rho:.3f}, p = {p:.2e}). "
            f"No reliable relationship between sequence and expression conservation. "
            f"Procrustes geometry is independent of sequence-level conservation."
        )

    print(f"\n  OUTCOME: {outcome}")
    print(f"\n  {interpretation}")

    # Per-cell-type summary
    n_sig_ct = int(per_celltype["significant"].sum())
    median_rho = float(per_celltype["spearman_rho"].median())
    print(f"\n  Per-cell-type correlations:")
    print(f"    Significant (p < 0.05): {n_sig_ct}/35")
    print(f"    Median rho: {median_rho:.4f}")

    if "lineage" in per_celltype.columns:
        print(f"\n    By lineage:")
        for lin in sorted(per_celltype["lineage"].dropna().unique()):
            subset = per_celltype[per_celltype["lineage"] == lin]
            print(f"      {lin:<20} median rho = {subset['spearman_rho'].median():.4f} "
                  f"(n={len(subset)})")

    # Extremes
    cat1 = extremes["cat1_conserved_seq_divergent_expr"]
    cat2 = extremes["cat2_divergent_seq_conserved_expr"]
    print(f"\n  Decoupled genes:")
    print(f"    Conserved seq + divergent expr: {cat1['n_genes']} genes")
    print(f"    Divergent seq + conserved expr: {cat2['n_genes']} genes")

    if cat1.get("go", {}).get("n_terms", 0) > 0:
        top = cat1["go"]["top_terms"][0]
        print(f"    Cat1 top GO: {top['term'].split(' (GO:')[0]} (p_adj={top['p_adj']:.2e})")
    if cat2.get("go", {}).get("n_terms", 0) > 0:
        top = cat2["go"]["top_terms"][0]
        print(f"    Cat2 top GO: {top['term'].split(' (GO:')[0]} (p_adj={top['p_adj']:.2e})")

    # Save summary JSON
    summary = {
        "outcome": outcome,
        "interpretation": interpretation,
        "metric": "protein_percent_identity",
        "metric_note": (
            "dN/dS ratios are no longer available in Ensembl Compara (releases 109-115). "
            "Protein percent identity is used instead — it directly measures amino acid "
            "conservation and correlates strongly with dN/dS where both are available."
        ),
        "genome_wide": genome_wide,
        "per_celltype_summary": {
            "n_significant": n_sig_ct,
            "n_total": len(per_celltype),
            "median_rho": median_rho,
            "range_rho": [float(per_celltype["spearman_rho"].min()),
                          float(per_celltype["spearman_rho"].max())],
        },
        "extremes": extremes,
        "data_summary": {
            "n_genes_total": n_genes_total,
            "n_seq_usable": n_seq_usable,
        },
    }

    with open(OUTPUT_DIR / "dnds_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\n  Summary saved: {OUTPUT_DIR / 'dnds_summary.json'}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    t_start = time.time()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("SEQUENCE CONSERVATION vs EXPRESSION DIVERGENCE")
    print("Adversarial validation of Procrustes geometric signal")
    print("=" * 70)

    # Step 2 first (doesn't need Ensembl)
    loadings_df, gene_ids, cell_types = step2_gene_loadings()

    # Step 1: Load sequence conservation
    seq_df = step1_load_sequence_conservation(gene_ids)
    n_seq_usable = len(seq_df)

    # Step 3: Genome-wide correlation
    genome_wide = step3_genome_wide(seq_df, loadings_df)

    # Step 4: Per-cell-type correlations
    per_celltype = step4_per_celltype(seq_df, loadings_df, cell_types)

    # Step 5: Extremes analysis
    extremes = step5_extremes(seq_df, loadings_df)

    # Step 6: Summary
    step6_summary(genome_wide, per_celltype, extremes, len(gene_ids), n_seq_usable)

    t_total = time.time() - t_start
    print(f"\n  Total runtime: {t_total:.1f}s ({t_total/60:.1f} min)")
    print(f"  All outputs: {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
