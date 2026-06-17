#!/usr/bin/env python3
"""
PPI Network Centrality vs Evolutionary Rigidity (T3-A)
======================================================
Tests whether protein-protein interaction network centrality predicts
evolutionary rigidity across 35 cell types. First protein-level test
after six expression-level nulls.

Biology: If rigid cell types use more network-central proteins (positive ρ),
rigidity reflects deep embeddedness in PPI architecture. If rigid types
use more peripheral proteins (negative ρ), rigidity reflects dedicated
specialized function over network flexibility.

Math: For each gene, compute degree centrality, PageRank, and betweenness
centrality in the STRING PPI network at three confidence thresholds
(400, 700, 900). Aggregate gene-level centrality to cell-type-level
scores via three methods (all expressed genes, top deformation genes,
identity genes). Correlate 27 centrality × rigidity combinations via
Spearman, with FDR correction.
"""

import gzip
import json
import os
import sys
import urllib.request

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from sklearn.decomposition import PCA

# ─── Paths ───────────────────────────────────────────────────────────
OUT_DIR = "output/mechanistic/ppi_centrality"
CENTROIDS_CSV = "output/phase2/scaled_35types/centroids_human_35.csv"
RESIDUALS_CSV = "output/phase2/scaled_35types/residuals_ranked.csv"
PROCRUSTES_JSON = "output/phase2/scaled_35types/procrustes_results_35.json"
PCA_NPZ = "output/phase2/scaled_35types/pca_centroids_35.npz"
ORTHOLOGS_CSV = "data/phase1/orthologs_human_mouse.csv"

STRING_LINKS_URL = (
    "https://stringdb-downloads.org/download/"
    "protein.links.v12.0/9606.protein.links.v12.0.txt.gz"
)
STRING_INFO_URL = (
    "https://stringdb-downloads.org/download/"
    "protein.info.v12.0/9606.protein.info.v12.0.txt.gz"
)
STRING_LINKS_FILE = "data/string/9606.protein.links.v12.0.txt.gz"
STRING_INFO_FILE = "data/string/9606.protein.info.v12.0.txt.gz"

RANDOM_SEED = 42
CONFIDENCE_THRESHOLDS = [400, 700, 900]
BETWEENNESS_K = 500  # approximate betweenness pivots
TOP_DEFORMATION_GENES = 200
TOP_IDENTITY_GENES = 500
N_BOOTSTRAP = 1000

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs("data/string", exist_ok=True)

np.random.seed(RANDOM_SEED)


# ═══════════════════════════════════════════════════════════════════
# STEP 1: Download STRING data
# ═══════════════════════════════════════════════════════════════════

def download_string():
    """Download STRING v12.0 human protein links and info files."""
    print("=" * 70)
    print("STEP 1: Download STRING v12.0 data")
    print("=" * 70)

    for url, path in [(STRING_LINKS_URL, STRING_LINKS_FILE),
                      (STRING_INFO_URL, STRING_INFO_FILE)]:
        if os.path.exists(path):
            size_mb = os.path.getsize(path) / 1e6
            print(f"  Already exists: {path} ({size_mb:.1f} MB)")
        else:
            print(f"  Downloading: {os.path.basename(path)}...")
            urllib.request.urlretrieve(url, path)
            size_mb = os.path.getsize(path) / 1e6
            print(f"  Downloaded: {size_mb:.1f} MB")


# ═══════════════════════════════════════════════════════════════════
# STEP 2: Parse STRING and map to our gene space
# ═══════════════════════════════════════════════════════════════════

def parse_string():
    """Parse STRING protein info to get ENSP→gene symbol mapping,
    then filter to our 16,959 ortholog gene space.

    Returns:
        protein_to_gene: dict mapping STRING protein ID to gene symbol
        gene_to_protein: dict mapping gene symbol to STRING protein ID
        our_genes_in_string: set of our genes present in STRING
        links_df: DataFrame with columns [protein1, protein2, combined_score]
    """
    print("\n" + "=" * 70)
    print("STEP 2: Parse STRING and map to ortholog gene space")
    print("=" * 70)

    # Load our gene space: ENSG → gene symbol
    orthologs = pd.read_csv(ORTHOLOGS_CSV)
    ensg_to_symbol = dict(zip(orthologs["human_ensembl_id"],
                              orthologs["human_gene_name"]))

    # Our centroids use ENSG IDs; get the gene symbol set
    centroids = pd.read_csv(CENTROIDS_CSV, nrows=0)
    our_ensg_ids = [c for c in centroids.columns if c.startswith("ENSG")]
    our_symbols = set()
    ensg_to_sym_our = {}
    for ensg in our_ensg_ids:
        sym = ensg_to_symbol.get(ensg)
        if sym:
            our_symbols.add(sym)
            ensg_to_sym_our[ensg] = sym

    print(f"  Our ortholog space: {len(our_ensg_ids):,} ENSG IDs")
    print(f"  Mapped to symbols: {len(our_symbols):,} gene symbols")

    # Parse STRING protein.info to map ENSP → gene symbol
    print("  Parsing STRING protein info...")
    ensp_to_gene = {}
    with gzip.open(STRING_INFO_FILE, "rt") as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                ensp_id = parts[0]  # e.g. 9606.ENSP00000000233
                gene_name = parts[1]  # preferred name = gene symbol
                ensp_to_gene[ensp_id] = gene_name

    print(f"  STRING proteins: {len(ensp_to_gene):,}")

    # Filter to our genes
    protein_to_gene = {}
    gene_to_protein = {}
    for ensp, gene in ensp_to_gene.items():
        if gene in our_symbols:
            protein_to_gene[ensp] = gene
            gene_to_protein[gene] = ensp  # last one wins if duplicates

    our_genes_in_string = set(protein_to_gene.values())
    print(f"  Our genes in STRING: {len(our_genes_in_string):,} / {len(our_symbols):,}")

    # Report coverage at each threshold
    print("\n  Parsing STRING links...")
    links = []
    with gzip.open(STRING_LINKS_FILE, "rt") as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split()
            if len(parts) == 3:
                p1, p2, score = parts[0], parts[1], int(parts[2])
                # Only keep edges where both proteins are in our gene set
                if p1 in protein_to_gene and p2 in protein_to_gene:
                    g1 = protein_to_gene[p1]
                    g2 = protein_to_gene[p2]
                    if g1 != g2:  # no self-loops
                        links.append((g1, g2, score))

    links_df = pd.DataFrame(links, columns=["gene1", "gene2", "combined_score"])
    print(f"  Edges in our gene space: {len(links_df):,}")

    # Coverage at each threshold
    for thresh in CONFIDENCE_THRESHOLDS:
        subset = links_df[links_df["combined_score"] >= thresh]
        genes_at_thresh = set(subset["gene1"]) | set(subset["gene2"])
        n_in = len(genes_at_thresh & our_symbols)
        print(f"  Genes at threshold {thresh}: {n_in:,} / {len(our_symbols):,} "
              f"({100*n_in/len(our_symbols):.1f}%)")

    return protein_to_gene, gene_to_protein, our_genes_in_string, links_df, our_symbols, ensg_to_sym_our


# ═══════════════════════════════════════════════════════════════════
# STEP 3: Build networks and compute centrality
# ═══════════════════════════════════════════════════════════════════

def build_networks_and_centrality(links_df, our_symbols):
    """Build networks at 3 confidence thresholds, compute centrality metrics.

    Returns:
        centrality_data: dict mapping (threshold, metric) → dict of gene → score
        network_stats: list of dicts with network statistics
    """
    print("\n" + "=" * 70)
    print("STEP 3: Build networks and compute centrality")
    print("=" * 70)

    centrality_data = {}
    network_stats = []

    for thresh in CONFIDENCE_THRESHOLDS:
        print(f"\n  --- Threshold {thresh} ---")
        subset = links_df[links_df["combined_score"] >= thresh]

        G = nx.Graph()
        # Add all our genes as nodes (even if no edges at this threshold)
        G.add_nodes_from(our_symbols)
        for _, row in subset.iterrows():
            G.add_edge(row["gene1"], row["gene2"])

        n_nodes = G.number_of_nodes()
        n_edges = G.number_of_edges()
        degrees = [d for _, d in G.degree()]
        mean_deg = np.mean(degrees)
        max_possible_edges = n_nodes * (n_nodes - 1) / 2
        density = n_edges / max_possible_edges if max_possible_edges > 0 else 0

        stats_dict = {
            "threshold": thresh,
            "n_nodes": n_nodes,
            "n_edges": n_edges,
            "mean_degree": round(mean_deg, 2),
            "density": round(density, 6),
        }
        network_stats.append(stats_dict)
        print(f"    Nodes: {n_nodes:,}, Edges: {n_edges:,}")
        print(f"    Mean degree: {mean_deg:.2f}, Density: {density:.6f}")

        # Degree centrality (normalized)
        print("    Computing degree centrality...")
        dc = nx.degree_centrality(G)
        centrality_data[(thresh, "degree")] = dc

        # PageRank
        print("    Computing PageRank...")
        pr = nx.pagerank(G, alpha=0.85)
        centrality_data[(thresh, "pagerank")] = pr

        # Betweenness centrality (approximate, k=500 pivots)
        print(f"    Computing betweenness centrality (approximate, k={BETWEENNESS_K})...")
        bc = nx.betweenness_centrality(G, k=BETWEENNESS_K, seed=RANDOM_SEED)
        centrality_data[(thresh, "betweenness")] = bc

    return centrality_data, network_stats


# ═══════════════════════════════════════════════════════════════════
# STEP 4: Prepare cell-type gene sets for aggregation
# ═══════════════════════════════════════════════════════════════════

def prepare_gene_sets(ensg_to_sym):
    """Prepare the three gene sets for cell-type aggregation:
    A) All expressed genes per cell type
    B) Top 200 deformation genes per cell type
    C) Identity genes expressed per cell type

    Returns:
        expressed_genes: dict of cell_type → set of gene symbols
        deformation_genes: dict of cell_type → list of top 200 gene symbols
        identity_genes: set of top 500 identity gene symbols
    """
    print("\n" + "=" * 70)
    print("STEP 4: Prepare gene sets for aggregation")
    print("=" * 70)

    # Load centroids
    centroids = pd.read_csv(CENTROIDS_CSV, index_col=0)
    gene_cols = [c for c in centroids.columns if c.startswith("ENSG")]
    cell_types = list(centroids.index)

    # Map ENSG columns to gene symbols
    gene_symbols = [ensg_to_sym.get(g, g) for g in gene_cols]
    centroids_sym = centroids[gene_cols].copy()
    centroids_sym.columns = gene_symbols

    # --- Method A: All expressed genes (mean expression > 0) ---
    print("  Method A: Expressed genes per cell type")
    expressed_genes = {}
    for ct in cell_types:
        expr = centroids_sym.loc[ct]
        expressed = set(expr[expr > 0].index)
        expressed_genes[ct] = expressed

    n_expressed = [len(v) for v in expressed_genes.values()]
    print(f"    Range: {min(n_expressed):,} – {max(n_expressed):,} expressed genes")

    # --- Method B: Top 200 deformation genes per cell type ---
    print("  Method B: Top 200 deformation genes per cell type")

    # Refit PCA on combined human+mouse centroids (same as original pipeline)
    # Original pipeline stacked 70 centroids → PCA → 33 components
    mouse_centroids = pd.read_csv(
        "output/phase2/scaled_35types/centroids_mouse_35.csv", index_col=0)
    X_human = centroids[gene_cols].values.astype(np.float64)
    X_mouse = mouse_centroids[gene_cols].values.astype(np.float64)
    X_combined = np.vstack([X_human, X_mouse])
    pca = PCA(n_components=0.95, svd_solver="full", random_state=RANDOM_SEED)
    pca.fit(X_combined)
    W = pca.components_  # (n_comp, n_genes)
    print(f"    PCA: {pca.n_components_} components (combined human+mouse), "
          f"{sum(pca.explained_variance_ratio_)*100:.1f}% variance")

    # Load residual vectors from Procrustes JSON
    with open(PROCRUSTES_JSON) as f:
        proc = json.load(f)

    deformation_genes = {}
    for ct in cell_types:
        residual_pca = np.array(proc["residuals"][ct]["vector_pca"])
        # Project from PCA space to gene space
        gene_loadings = residual_pca @ W  # (n_genes,)
        abs_loadings = np.abs(gene_loadings)
        top_idx = np.argsort(abs_loadings)[::-1][:TOP_DEFORMATION_GENES]
        top_syms = [gene_symbols[i] for i in top_idx]
        deformation_genes[ct] = top_syms

    print(f"    {TOP_DEFORMATION_GENES} genes per cell type")

    # --- Method C: Identity genes (top 500 by variance across centroids) ---
    print("  Method C: Top 500 identity genes")
    gene_variances = centroids_sym.var(axis=0)
    top_identity = gene_variances.nlargest(TOP_IDENTITY_GENES).index.tolist()
    identity_genes = set(top_identity)
    print(f"    {len(identity_genes)} identity genes")
    print(f"    Examples: {top_identity[:10]}")

    return expressed_genes, deformation_genes, identity_genes, cell_types, centroids_sym


# ═══════════════════════════════════════════════════════════════════
# STEP 5: Aggregate centrality to cell-type level
# ═══════════════════════════════════════════════════════════════════

def aggregate_centrality(centrality_data, expressed_genes, deformation_genes,
                         identity_genes, cell_types, centroids_sym, our_symbols):
    """Aggregate gene-level centrality to cell-type-level scores.

    For each cell type, compute mean centrality across:
      A) All expressed genes
      B) Top 200 deformation genes
      C) Identity genes expressed in that cell type

    Returns:
        scores: DataFrame with columns [cell_type, threshold, metric, method, score]
    """
    print("\n" + "=" * 70)
    print("STEP 5: Aggregate centrality to cell-type level")
    print("=" * 70)

    rows = []
    for thresh in CONFIDENCE_THRESHOLDS:
        for metric in ["degree", "pagerank", "betweenness"]:
            cent = centrality_data[(thresh, metric)]

            for ct in cell_types:
                # Method A: all expressed genes
                genes_a = expressed_genes[ct] & our_symbols
                vals_a = [cent.get(g, 0) for g in genes_a]
                score_a = np.mean(vals_a) if vals_a else 0

                # Method B: top 200 deformation genes
                genes_b = [g for g in deformation_genes[ct] if g in our_symbols]
                vals_b = [cent.get(g, 0) for g in genes_b]
                score_b = np.mean(vals_b) if vals_b else 0

                # Method C: identity genes expressed in this cell type
                expr = centroids_sym.loc[ct]
                genes_c = [g for g in identity_genes
                           if g in our_symbols and expr.get(g, 0) > 0]
                vals_c = [cent.get(g, 0) for g in genes_c]
                score_c = np.mean(vals_c) if vals_c else 0

                rows.append({
                    "cell_type": ct, "threshold": thresh,
                    "metric": metric, "method": "A_expressed",
                    "score": score_a
                })
                rows.append({
                    "cell_type": ct, "threshold": thresh,
                    "metric": metric, "method": "B_deformation",
                    "score": score_b
                })
                rows.append({
                    "cell_type": ct, "threshold": thresh,
                    "metric": metric, "method": "C_identity",
                    "score": score_c
                })

    scores_df = pd.DataFrame(rows)
    print(f"  Aggregated {len(scores_df)} scores "
          f"({len(cell_types)} types × {len(CONFIDENCE_THRESHOLDS)} thresholds "
          f"× 3 metrics × 3 methods)")
    return scores_df


# ═══════════════════════════════════════════════════════════════════
# STEP 6: Correlation analysis
# ═══════════════════════════════════════════════════════════════════

def correlation_analysis(scores_df):
    """Compute Spearman correlation between cell-type centrality and
    Procrustes residual magnitude for all 27 combinations.

    Returns:
        results_df: DataFrame with ρ, p-value, FDR q-value for each combo
    """
    print("\n" + "=" * 70)
    print("STEP 6: Spearman correlation analysis (27 combinations)")
    print("=" * 70)

    # Load residual magnitudes
    residuals = pd.read_csv(RESIDUALS_CSV)
    resid_map = dict(zip(residuals["cell_type"], residuals["residual_magnitude"]))

    results = []
    for thresh in CONFIDENCE_THRESHOLDS:
        for metric in ["degree", "pagerank", "betweenness"]:
            for method in ["A_expressed", "B_deformation", "C_identity"]:
                mask = ((scores_df["threshold"] == thresh) &
                        (scores_df["metric"] == metric) &
                        (scores_df["method"] == method))
                subset = scores_df[mask].copy()
                subset["residual"] = subset["cell_type"].map(resid_map)
                subset = subset.dropna(subset=["residual", "score"])

                rho, pval = stats.spearmanr(subset["score"], subset["residual"])

                results.append({
                    "threshold": thresh,
                    "metric": metric,
                    "method": method,
                    "combo": f"{metric}_{method}_t{thresh}",
                    "spearman_rho": round(rho, 4),
                    "p_value": pval,
                    "n_cell_types": len(subset),
                })

    results_df = pd.DataFrame(results)

    # FDR correction (Benjamini-Hochberg)
    from statsmodels.stats.multitest import multipletests
    reject, qvals, _, _ = multipletests(results_df["p_value"],
                                         method="fdr_bh")
    results_df["fdr_q"] = qvals
    results_df["significant"] = reject

    # Print summary
    print(f"\n  {'Combination':<35} {'ρ':>8} {'p':>10} {'q_FDR':>10} {'Sig':>5}")
    print("  " + "-" * 70)
    for _, row in results_df.sort_values("p_value").iterrows():
        sig = " ***" if row["significant"] else ""
        print(f"  {row['combo']:<35} {row['spearman_rho']:>8.4f} "
              f"{row['p_value']:>10.4f} {row['fdr_q']:>10.4f}{sig}")

    n_sig = results_df["significant"].sum()
    print(f"\n  Significant after FDR correction: {n_sig} / {len(results_df)}")

    return results_df


# ═══════════════════════════════════════════════════════════════════
# STEP 7: Sensitivity analysis on best combination
# ═══════════════════════════════════════════════════════════════════

def sensitivity_analysis(scores_df, results_df):
    """For the best-performing combination:
    - Top 5 / bottom 5 cell types by centrality
    - Exclude 5 progenitor types and re-test
    - Bootstrap 1,000 resamples for 95% CI on ρ
    """
    print("\n" + "=" * 70)
    print("STEP 7: Sensitivity analysis")
    print("=" * 70)

    residuals = pd.read_csv(RESIDUALS_CSV)
    resid_map = dict(zip(residuals["cell_type"], residuals["residual_magnitude"]))

    # Find best combination by |ρ|
    best_idx = results_df["spearman_rho"].abs().idxmax()
    best = results_df.loc[best_idx]
    print(f"  Best combination: {best['combo']}")
    print(f"    ρ = {best['spearman_rho']:.4f}, p = {best['p_value']:.4f}, "
          f"q_FDR = {best['fdr_q']:.4f}")

    # Get scores for best combination
    mask = ((scores_df["threshold"] == best["threshold"]) &
            (scores_df["metric"] == best["metric"]) &
            (scores_df["method"] == best["method"]))
    subset = scores_df[mask].copy()
    subset["residual"] = subset["cell_type"].map(resid_map)
    subset = subset.dropna(subset=["residual", "score"])
    subset = subset.sort_values("score", ascending=False)

    # Top 5 and bottom 5 by centrality
    print("\n  Top 5 by centrality score:")
    for _, row in subset.head(5).iterrows():
        rank_idx = residuals[residuals["cell_type"] == row["cell_type"]]["rank"].values
        r = rank_idx[0] if len(rank_idx) > 0 else "?"
        print(f"    {row['cell_type']:<45} centrality={row['score']:.6f}  "
              f"rigidity_rank={r}")

    print("\n  Bottom 5 by centrality score:")
    for _, row in subset.tail(5).iterrows():
        rank_idx = residuals[residuals["cell_type"] == row["cell_type"]]["rank"].values
        r = rank_idx[0] if len(rank_idx) > 0 else "?"
        print(f"    {row['cell_type']:<45} centrality={row['score']:.6f}  "
              f"rigidity_rank={r}")

    # Progenitor exclusion
    progenitor_types = [
        "hematopoietic stem cell",
        "hematopoietic precursor cell",
        "mesenchymal stem cell",
        "mesenchymal stem cell of adipose tissue",
        "basal cell",  # often progenitor-like
    ]
    non_prog = subset[~subset["cell_type"].isin(progenitor_types)]
    rho_excl, p_excl = stats.spearmanr(non_prog["score"], non_prog["residual"])
    print(f"\n  Excluding {len(progenitor_types)} progenitor types "
          f"(n={len(non_prog)}):")
    print(f"    ρ = {rho_excl:.4f}, p = {p_excl:.4f}")

    # Bootstrap
    print(f"\n  Bootstrap ({N_BOOTSTRAP} resamples)...")
    rng = np.random.RandomState(RANDOM_SEED)
    boot_rhos = []
    scores_arr = subset["score"].values
    resid_arr = subset["residual"].values
    n = len(scores_arr)

    for _ in range(N_BOOTSTRAP):
        idx = rng.choice(n, size=n, replace=True)
        r, _ = stats.spearmanr(scores_arr[idx], resid_arr[idx])
        boot_rhos.append(r)

    boot_rhos = np.array(boot_rhos)
    ci_lo = np.percentile(boot_rhos, 2.5)
    ci_hi = np.percentile(boot_rhos, 97.5)
    print(f"    95% CI on ρ: [{ci_lo:.4f}, {ci_hi:.4f}]")
    print(f"    Bootstrap mean ρ: {np.mean(boot_rhos):.4f}")
    print(f"    Bootstrap SD: {np.std(boot_rhos):.4f}")

    sensitivity = {
        "best_combo": best["combo"],
        "best_rho": float(best["spearman_rho"]),
        "best_p": float(best["p_value"]),
        "best_fdr_q": float(best["fdr_q"]),
        "progenitor_excluded_rho": round(rho_excl, 4),
        "progenitor_excluded_p": round(p_excl, 4),
        "progenitor_excluded_n": len(non_prog),
        "bootstrap_ci_95": [round(ci_lo, 4), round(ci_hi, 4)],
        "bootstrap_mean_rho": round(float(np.mean(boot_rhos)), 4),
        "bootstrap_sd": round(float(np.std(boot_rhos)), 4),
    }

    return sensitivity, subset, best


# ═══════════════════════════════════════════════════════════════════
# STEP 8: Visualizations and output
# ═══════════════════════════════════════════════════════════════════

def generate_outputs(results_df, sensitivity, best_subset, best_row,
                     network_stats, our_symbols_count, string_gene_count):
    """Generate heatmap, scatter plot, summary table, and text summary."""
    print("\n" + "=" * 70)
    print("STEP 8: Generate outputs")
    print("=" * 70)

    # --- 1. Heatmap: 27 combinations × Spearman ρ ---
    pivot = results_df.pivot_table(
        index=["metric", "method"],
        columns="threshold",
        values="spearman_rho"
    )
    # Create readable labels
    label_map = {
        "A_expressed": "All expressed",
        "B_deformation": "Top-200 deformation",
        "C_identity": "Top-500 identity",
    }
    metric_map = {"degree": "Degree", "pagerank": "PageRank",
                  "betweenness": "Betweenness"}

    new_index = []
    for metric, method in pivot.index:
        new_index.append(f"{metric_map[metric]} × {label_map[method]}")
    pivot.index = new_index

    fig, ax = plt.subplots(figsize=(8, 7))
    vmax = max(abs(pivot.values.min()), abs(pivot.values.max()))
    sns.heatmap(pivot, annot=True, fmt=".3f", center=0,
                cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                linewidths=0.5, ax=ax)
    ax.set_title("PPI Centrality vs Evolutionary Rigidity\nSpearman ρ across 35 cell types",
                 fontsize=13, fontweight="bold")
    ax.set_xlabel("STRING confidence threshold")
    ax.set_ylabel("Centrality metric × Gene aggregation method")
    plt.tight_layout()
    heatmap_path = os.path.join(OUT_DIR, "centrality_rigidity_heatmap.png")
    plt.savefig(heatmap_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {heatmap_path}")

    # --- 2. Scatter plot for best combination ---
    residuals = pd.read_csv(RESIDUALS_CSV)
    resid_map = dict(zip(residuals["cell_type"], residuals["rank"]))

    fig, ax = plt.subplots(figsize=(10, 8))
    x = best_subset["score"].values
    y = best_subset["residual"].values
    ax.scatter(x, y, s=60, alpha=0.7, edgecolors="k", linewidths=0.5)

    # Label all points
    for _, row in best_subset.iterrows():
        name = row["cell_type"]
        # Abbreviate long names
        if len(name) > 25:
            name = name[:22] + "..."
        ax.annotate(name, (row["score"], row["residual"]),
                    fontsize=6, alpha=0.8,
                    xytext=(4, 4), textcoords="offset points")

    # Fit line
    slope, intercept = np.polyfit(x, y, 1)
    x_fit = np.linspace(x.min(), x.max(), 100)
    ax.plot(x_fit, slope * x_fit + intercept, "r--", alpha=0.5)

    rho = sensitivity["best_rho"]
    p = sensitivity["best_p"]
    ax.set_title(f"Best combination: {sensitivity['best_combo']}\n"
                 f"Spearman ρ = {rho:.4f}, p = {p:.4f}",
                 fontsize=12, fontweight="bold")
    ax.set_xlabel("Mean PPI centrality score", fontsize=11)
    ax.set_ylabel("Procrustes residual magnitude (rigidity)", fontsize=11)
    plt.tight_layout()
    scatter_path = os.path.join(OUT_DIR, "best_combo_scatter.png")
    plt.savefig(scatter_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {scatter_path}")

    # --- 3. Summary table ---
    summary = results_df[["combo", "threshold", "metric", "method",
                           "spearman_rho", "p_value", "fdr_q",
                           "significant", "n_cell_types"]].copy()
    summary = summary.sort_values("p_value")
    table_path = os.path.join(OUT_DIR, "correlation_summary.csv")
    summary.to_csv(table_path, index=False)
    print(f"  Saved: {table_path}")

    # --- 4. Full results JSON ---
    full_results = {
        "network_stats": network_stats,
        "gene_coverage": {
            "our_ortholog_genes": our_symbols_count,
            "genes_in_string": string_gene_count,
            "pct_coverage": round(100 * string_gene_count / our_symbols_count, 1),
        },
        "correlation_results": summary.to_dict(orient="records"),
        "n_significant_fdr": int(results_df["significant"].sum()),
        "sensitivity": sensitivity,
        "methods_note": (
            "Betweenness centrality computed using approximate method "
            f"with k={BETWEENNESS_K} pivots (NetworkX). "
            "Exact computation infeasible at this network size."
        ),
    }
    json_path = os.path.join(OUT_DIR, "ppi_centrality_results.json")
    with open(json_path, "w") as f:
        json.dump(full_results, f, indent=2, default=str)
    print(f"  Saved: {json_path}")

    # --- 5. Text summary ---
    print("\n" + "=" * 70)
    print("SUMMARY: PPI Centrality vs Evolutionary Rigidity (T3-A)")
    print("=" * 70)

    n_sig = int(results_df["significant"].sum())
    print(f"\n  Network: STRING v12.0 (human, {our_symbols_count:,} ortholog genes)")
    for ns in network_stats:
        print(f"    Threshold {ns['threshold']}: {ns['n_nodes']:,} nodes, "
              f"{ns['n_edges']:,} edges, mean degree {ns['mean_degree']}")

    print(f"\n  27 Spearman correlations (3 thresholds × 3 metrics × 3 methods)")
    print(f"  Significant after FDR correction: {n_sig} / 27")

    if n_sig > 0:
        sig_rows = results_df[results_df["significant"]]
        print(f"\n  Significant combinations:")
        for _, row in sig_rows.iterrows():
            direction = "positive" if row["spearman_rho"] > 0 else "negative"
            print(f"    {row['combo']}: ρ={row['spearman_rho']:.4f}, "
                  f"p={row['p_value']:.4f}, q={row['fdr_q']:.4f} ({direction})")

    # Direction interpretation
    best_rho = sensitivity["best_rho"]
    if abs(best_rho) < 0.2:
        interp = ("WEAK/NULL relationship — PPI network centrality does not "
                  "meaningfully predict evolutionary rigidity.")
    elif best_rho > 0:
        interp = ("POSITIVE correlation — rigid cell types tend to use more "
                  "network-central proteins. Rigidity may reflect deep embeddedness "
                  "in protein interaction architecture.")
    else:
        interp = ("NEGATIVE correlation — rigid cell types tend to use more "
                  "peripheral/specialized proteins. Rigidity may reflect dedicated "
                  "function over network flexibility.")

    print(f"\n  Best combination: {sensitivity['best_combo']}")
    print(f"    ρ = {best_rho:.4f}, p = {sensitivity['best_p']:.4f}, "
          f"q_FDR = {sensitivity['best_fdr_q']:.4f}")
    print(f"    95% CI (bootstrap): [{sensitivity['bootstrap_ci_95'][0]:.4f}, "
          f"{sensitivity['bootstrap_ci_95'][1]:.4f}]")
    print(f"    Without progenitors: ρ = {sensitivity['progenitor_excluded_rho']:.4f}, "
          f"p = {sensitivity['progenitor_excluded_p']:.4f}")
    print(f"\n  Interpretation: {interp}")

    # Plot descriptions
    print("\n  PLOT DESCRIPTIONS:")
    print(f"  - {heatmap_path}: Heatmap showing Spearman ρ for all 27 combinations.")
    print(f"    Rows are centrality metric × aggregation method, columns are STRING")
    print(f"    confidence thresholds. Red = positive ρ, blue = negative ρ.")
    print(f"  - {scatter_path}: Scatter plot of mean PPI centrality vs Procrustes")
    print(f"    residual magnitude for the best combination ({sensitivity['best_combo']}).")
    print(f"    Each point is a cell type, labeled. Red dashed line = linear fit.")


# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════

def main():
    # Step 1: Download STRING
    download_string()

    # Step 2: Parse STRING and map to our gene space
    (protein_to_gene, gene_to_protein, our_genes_in_string,
     links_df, our_symbols, ensg_to_sym) = parse_string()

    # Step 3: Build networks and compute centrality
    centrality_data, network_stats = build_networks_and_centrality(
        links_df, our_symbols)

    # Step 4: Prepare gene sets
    (expressed_genes, deformation_genes, identity_genes,
     cell_types, centroids_sym) = prepare_gene_sets(ensg_to_sym)

    # Step 5: Aggregate centrality to cell-type level
    scores_df = aggregate_centrality(
        centrality_data, expressed_genes, deformation_genes,
        identity_genes, cell_types, centroids_sym, our_symbols)

    # Step 6: Correlation analysis
    results_df = correlation_analysis(scores_df)

    # Step 7: Sensitivity analysis
    sensitivity, best_subset, best_row = sensitivity_analysis(
        scores_df, results_df)

    # Step 8: Generate outputs
    generate_outputs(results_df, sensitivity, best_subset, best_row,
                     network_stats,
                     our_symbols_count=len(our_symbols),
                     string_gene_count=len(our_genes_in_string))

    print("\n  Done. All outputs saved to:", OUT_DIR)


if __name__ == "__main__":
    main()
