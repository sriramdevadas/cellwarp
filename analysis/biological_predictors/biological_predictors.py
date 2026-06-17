#!/usr/bin/env python3
"""
CellWarp — Biological Predictors of Cell Type Rigidity

Asks: what biological properties of a cell type predict its Procrustes
rigidity (residual magnitude)? This addresses a reviewer concern:
"What do I now know about biology?"

Biology: Cell types that are more "rigid" (low Procrustes residuals) maintain
their relative transcriptomic geometry across species. We test whether
standard biological covariates — transcriptomic complexity, cell cycle
activity, tissue breadth, marker gene counts, TF diversity, evolutionary
conservation, inter-donor variability — explain the rigidity pattern.

Math: Spearman rank correlations (univariate), elastic net regression and
random forest (multivariate), with leave-one-out cross-validation to
guard against overfitting with n=35.

Inputs:
    output/phase2/scaled_35types/residuals_ranked.csv
    output/phase2/mechanistic/housekeeping/hk_ratio_vs_residual.csv
    output/phase2/mechanistic/tf_complexity/tf_complexity_vs_residual.csv
    output/phase2/diagnostics/expression_level_vs_rigidity/normalized_rankings.csv
    output/phase2/diagnostics/interdonor_variance/interdonor_variance_by_celltype.csv
    output/phase2/cell_type_inventory.csv
    output/validation/t3e_chromatin/rigidity_conservation_merged.csv
    output/phase2/dnds/merged_seq_expr_divergence.csv
    data/phase2_scaled/human_raw_aligned.h5ad

Outputs:
    analysis/biological_predictors/feature_table.csv
    analysis/biological_predictors/univariate_correlations.csv
    analysis/biological_predictors/multivariate_model_results.json
    analysis/biological_predictors/biological_predictors_summary.md
    figures/supplementary/biological_predictors_*.pdf/.png
"""

import sys
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats
from scipy.spatial.distance import pdist
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ElasticNetCV
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import r2_score

warnings.filterwarnings("ignore", category=FutureWarning)

# ─── Paths ───────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from cellwarp.figure_style import (apply_style, save_figure, add_panel_label,
                          short_name, lineage_color, format_p,
                          LINEAGE_MAP, SHORT_NAMES, LINEAGE_COLORS,
                          COL1, COL15, COL2, C_BLUE, C_ORANGE, C_GRAY,
                          C_DARKGRAY, C_LIGHTGRAY, FONT_SIZE_ANNOT,
                          FONT_SIZE_LABEL, FONT_SIZE_LEGEND, DPI,
                          add_lineage_legend, clean_spine)

RESIDUALS_CSV = PROJECT_ROOT / "output/phase2/scaled_35types/residuals_ranked.csv"
HK_CSV = PROJECT_ROOT / "output/phase2/mechanistic/housekeeping/hk_ratio_vs_residual.csv"
TF_CSV = PROJECT_ROOT / "output/phase2/mechanistic/tf_complexity/tf_complexity_vs_residual.csv"
EXPR_CSV = PROJECT_ROOT / "output/phase2/diagnostics/expression_level_vs_rigidity/normalized_rankings.csv"
INTERDONOR_CSV = PROJECT_ROOT / "output/phase2/diagnostics/interdonor_variance/interdonor_variance_by_celltype.csv"
INVENTORY_CSV = PROJECT_ROOT / "output/phase2/cell_type_inventory.csv"
CHROMATIN_CSV = PROJECT_ROOT / "output/validation/t3e_chromatin/rigidity_conservation_merged.csv"
DNDS_CSV = PROJECT_ROOT / "output/phase2/dnds/merged_seq_expr_divergence.csv"
HUMAN_H5AD = PROJECT_ROOT / "data/phase2_scaled/human_raw_aligned.h5ad"

OUT_DIR = PROJECT_ROOT / "analysis/biological_predictors"
FIG_DIR = PROJECT_ROOT / "figures/supplementary"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

SEED = 42
np.random.seed(SEED)

# ─── Germ layer mapping (developmental origin) ──────────────────────────────
# Based on standard developmental biology: hematopoietic and mesenchymal from
# mesoderm, gut/liver/pancreas from endoderm, skin/mammary from ectoderm.
GERM_LAYER = {
    # Mesoderm — hematopoietic lineage
    'B cell': 'mesoderm', 'T cell': 'mesoderm', 'macrophage': 'mesoderm',
    'classical monocyte': 'mesoderm', 'monocyte': 'mesoderm',
    'natural killer cell': 'mesoderm', 'mature NK T cell': 'mesoderm',
    'hematopoietic precursor cell': 'mesoderm', 'hematopoietic stem cell': 'mesoderm',
    'CD8-positive, alpha-beta T cell': 'mesoderm',
    'CD4-positive, alpha-beta T cell': 'mesoderm',
    'intermediate monocyte': 'mesoderm', 'non-classical monocyte': 'mesoderm',
    'plasma cell': 'mesoderm', 'myeloid dendritic cell': 'mesoderm',
    'myeloid leukocyte': 'mesoderm', 'neutrophil': 'mesoderm',
    'granulocyte': 'mesoderm',
    # Mesoderm — stromal / endothelial
    'adventitial cell': 'mesoderm', 'fibroblast': 'mesoderm',
    'fibroblast of cardiac tissue': 'mesoderm',
    'mesenchymal stem cell': 'mesoderm',
    'mesenchymal stem cell of adipose tissue': 'mesoderm',
    'smooth muscle cell': 'mesoderm', 'stromal cell': 'mesoderm',
    'endothelial cell': 'mesoderm',
    # Endoderm — gut, liver, pancreas, bladder
    'hepatocyte': 'endoderm', 'pancreatic acinar cell': 'endoderm',
    'pancreatic ductal cell': 'endoderm',
    'enterocyte of epithelium of large intestine': 'endoderm',
    'large intestine goblet cell': 'endoderm',
    'bladder urothelial cell': 'endoderm', 'epithelial cell': 'endoderm',
    # Ectoderm — skin, mammary
    'basal cell': 'ectoderm',
    'luminal epithelial cell of mammary gland': 'ectoderm',
}

# ─── Tirosh et al. 2016 cell cycle marker genes ─────────────────────────────
# S phase genes
S_GENES = [
    'MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG',
    'GINS2', 'MCM6', 'CDCA7', 'DTL', 'PRIM1', 'UHRF1', 'MLF1IP',
    'HELLS', 'RFC2', 'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 'WDR76',
    'SLBP', 'CCNE2', 'UBR7', 'POLD3', 'MSH2', 'ATAD2', 'RAD51',
    'RRM2', 'CDC45', 'CDC6', 'EXO1', 'TIPIN', 'DSCC1', 'BLM',
    'CASP8AP2', 'USP1', 'CLSPN', 'POLA1', 'CHAF1B', 'BRIP1', 'E2F8',
]
# G2/M phase genes
G2M_GENES = [
    'HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A',
    'NDC80', 'CKS2', 'NUF2', 'CKS1B', 'MKI67', 'TMPO', 'CENPF',
    'TACC3', 'FAM64A', 'SMC4', 'CCNB2', 'CKAP2L', 'CKAP2', 'AURKB',
    'BUB1', 'KIF11', 'ANP32E', 'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP',
    'CDCA3', 'HN1', 'CDC20', 'TTK', 'CDC25C', 'KIF2C', 'RANGAP1',
    'NCAPD2', 'DLGAP5', 'CDCA2', 'CDCA8', 'ECT2', 'KIF23', 'HMMR',
    'AURKA', 'PSRC1', 'ANLN', 'LBR', 'CKAP5', 'CENPE', 'CTCF',
    'NEK2', 'G2E3', 'GAS2L3', 'CBX5', 'CENPA',
]


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 0: Load residuals (rigidity scores)
# ═══════════════════════════════════════════════════════════════════════════════
def load_residuals():
    """Load per-cell-type Procrustes residual magnitudes (rigidity scores).

    Higher residual = more divergent = more "flexible" across species.
    Lower residual = more conserved = more "rigid".
    """
    df = pd.read_csv(RESIDUALS_CSV)
    return df.set_index("cell_type")[["residual_magnitude"]].rename(
        columns={"residual_magnitude": "rigidity_score"}
    )


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1: Compute / retrieve biological features
# ═══════════════════════════════════════════════════════════════════════════════
def load_precomputed_features():
    """Load biological features already computed in prior analyses.

    Returns a DataFrame indexed by cell_type with columns:
    - avg_mean_expr: mean expression level per cell
    - hk_ratio_human/mouse: housekeeping gene expression ratio
    - n_active_tfs_human/mouse: number of active TFs
    - entropy_human/mouse: TF activity entropy (diversity)
    - mean_interdonor_var: mean inter-donor expression variance
    - human_count, mouse_count, min_count: cell counts in dataset
    - mean_phastCons: chromatin conservation score (phastCons)
    - lineage: broad lineage category
    - is_progenitor: whether cell type is a progenitor/stem cell
    """
    frames = {}

    # Expression level
    expr = pd.read_csv(EXPR_CSV)
    frames["expr"] = expr.set_index("cell_type")[["avg_mean_expr"]]

    # Housekeeping ratio
    hk = pd.read_csv(HK_CSV)
    frames["hk"] = hk.set_index("cell_type")[
        ["hk_ratio_human", "hk_ratio_mouse", "lineage", "progenitor"]
    ].rename(columns={"progenitor": "is_progenitor"})

    # TF complexity
    tf = pd.read_csv(TF_CSV)
    frames["tf"] = tf.set_index("cell_type")[
        ["n_active_tfs_human", "n_active_tfs_mouse",
         "entropy_human", "entropy_mouse"]
    ]

    # Inter-donor variance
    idv = pd.read_csv(INTERDONOR_CSV)
    frames["idv"] = idv.set_index("cell_type")[["mean_interdonor_var"]]

    # Cell counts (filter to passing 35 types)
    inv = pd.read_csv(INVENTORY_CSV)
    inv = inv[inv["passes_200_gate"]]
    frames["inv"] = inv.set_index("cell_type")[
        ["human_count", "mouse_count", "min_count"]
    ]

    # Chromatin conservation (phastCons, 2kb window, placental 20-way)
    chrom = pd.read_csv(CHROMATIN_CSV)
    frames["chrom"] = chrom.set_index("cell_type")[["mean_phastCons"]]

    # Merge all
    merged = frames["expr"]
    for key in ["hk", "tf", "idv", "inv", "chrom"]:
        merged = merged.join(frames[key], how="outer")

    print(f"  Loaded precomputed features for {len(merged)} cell types")
    return merged


def compute_h5ad_features(h5ad_path):
    """Compute biological features from the raw AnnData expression matrix.

    Computes per cell type:
    (a) Transcriptomic complexity: genes detected in >10% of cells
    (b) Cell cycle fraction: proportion of cells in S or G2/M phase
    (c) Within-type heterogeneity: mean variance across top 10 PCs
    (d) Tissue breadth: number of distinct tissues
    (e) Number of marker genes: genes with FDR<0.05 and log2FC>1

    Biology: These features capture fundamental properties of cell type
    identity — how many genes a type expresses (complexity), how
    proliferative it is (cell cycle), how variable individual cells are
    (heterogeneity), and how tissue-restricted it is (breadth).

    Math: Complexity = count(genes with >0 expression in >10% cells).
    Cell cycle = Tirosh scoring via mean expression of S/G2M gene sets.
    Heterogeneity = trace of covariance in PCA space (top 10 PCs).
    Markers = Wilcoxon rank-sum, Benjamini-Hochberg FDR correction.
    """
    import scanpy as sc

    print(f"  Loading {h5ad_path.name} ...")
    adata = ad.read_h5ad(h5ad_path)
    print(f"    Shape: {adata.shape[0]} cells × {adata.shape[1]} genes")
    print(f"    Cell types: {adata.obs['cell_type'].nunique()}")

    # Get gene names (feature_name column holds symbols)
    if "feature_name" in adata.var.columns:
        gene_symbols = adata.var["feature_name"].values
    else:
        gene_symbols = adata.var_names.values

    cell_types = sorted(adata.obs["cell_type"].unique())
    results = {ct: {} for ct in cell_types}

    # --- (a) Transcriptomic complexity ---
    print("  Computing transcriptomic complexity ...")
    for ct in cell_types:
        mask = adata.obs["cell_type"] == ct
        X_ct = adata[mask].X
        if hasattr(X_ct, "toarray"):
            X_ct = X_ct.toarray()
        n_cells = X_ct.shape[0]
        # Fraction of cells expressing each gene (> 0)
        frac_expressing = (X_ct > 0).sum(axis=0) / n_cells
        n_genes_detected = int((frac_expressing > 0.10).sum())
        results[ct]["n_genes_detected"] = n_genes_detected

    # --- (b) Cell cycle scoring ---
    print("  Scoring cell cycle phases (Tirosh et al. markers) ...")
    # Normalize for scoring (log1p of library-size-normalized counts)
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)

    # Map gene symbols for cell cycle scoring
    if "feature_name" in adata_norm.var.columns:
        adata_norm.var_names = adata_norm.var["feature_name"].values
        adata_norm.var_names_make_unique()

    s_genes_found = [g for g in S_GENES if g in adata_norm.var_names]
    g2m_genes_found = [g for g in G2M_GENES if g in adata_norm.var_names]
    print(f"    S-phase genes found: {len(s_genes_found)}/{len(S_GENES)}")
    print(f"    G2/M genes found: {len(g2m_genes_found)}/{len(G2M_GENES)}")

    sc.tl.score_genes_cell_cycle(adata_norm, s_genes=s_genes_found,
                                  g2m_genes=g2m_genes_found)

    for ct in cell_types:
        mask = adata_norm.obs["cell_type"] == ct
        phases = adata_norm.obs.loc[mask, "phase"]
        n_total = len(phases)
        n_cycling = int(((phases == "S") | (phases == "G2M")).sum())
        results[ct]["cell_cycle_fraction"] = n_cycling / n_total

    # --- (c) Within-type heterogeneity (PCA variance) ---
    print("  Computing within-type heterogeneity (PCA variance) ...")
    # Use the log-normalized data, compute PCA per type
    # Use top-2000 genes by variance (avoids skmisc dependency)
    all_var = np.array(adata_norm.X.toarray().var(axis=0) if hasattr(adata_norm.X, "toarray")
                       else adata_norm.X.var(axis=0)).flatten()
    top_gene_idx = np.argsort(all_var)[-2000:]
    for ct in cell_types:
        mask = adata_norm.obs["cell_type"] == ct
        X_ct = adata_norm[mask].X[:, top_gene_idx]
        if hasattr(X_ct, "toarray"):
            X_ct = X_ct.toarray()
        # Standardize
        X_ct = X_ct - X_ct.mean(axis=0)
        # SVD for top 10 PCs
        n_pcs = min(10, X_ct.shape[0] - 1, X_ct.shape[1])
        if n_pcs < 2:
            results[ct]["within_type_heterogeneity"] = np.nan
            continue
        from sklearn.decomposition import PCA as skPCA
        pca = skPCA(n_components=n_pcs, random_state=SEED)
        X_pca = pca.fit_transform(X_ct)
        # Total variance in PCA space = sum of variances along each PC
        total_var = X_pca.var(axis=0).sum()
        results[ct]["within_type_heterogeneity"] = float(total_var)

    # --- (d) Tissue breadth ---
    print("  Computing tissue breadth ...")
    tissue_col = "tissue_general" if "tissue_general" in adata.obs.columns else "tissue"
    for ct in cell_types:
        mask = adata.obs["cell_type"] == ct
        n_tissues = adata.obs.loc[mask, tissue_col].nunique()
        results[ct]["tissue_breadth"] = n_tissues

    # --- (e) Number of marker genes ---
    print("  Computing marker genes (Wilcoxon rank-sum, FDR<0.05, log2FC>1) ...")
    # Use log-normalized data for marker detection
    sc.tl.rank_genes_groups(adata_norm, groupby="cell_type", method="wilcoxon",
                            use_raw=False, key_added="markers")

    for ct in cell_types:
        result_df = sc.get.rank_genes_groups_df(adata_norm, group=ct,
                                                 key="markers")
        # Filter: FDR < 0.05 and log2FC > 1
        sig = result_df[(result_df["pvals_adj"] < 0.05) &
                        (result_df["logfoldchanges"] > 1.0)]
        results[ct]["n_marker_genes"] = len(sig)

    # Clean up
    del adata_norm

    # Convert to DataFrame
    feat_df = pd.DataFrame(results).T
    feat_df.index.name = "cell_type"
    print(f"  Computed h5ad features for {len(feat_df)} cell types")
    return feat_df, adata, gene_symbols


def compute_marker_conservation(adata, gene_symbols, cell_types, dnds_path):
    """Compute mean evolutionary conservation of each type's marker genes.

    For each cell type, identifies its top 100 marker genes and looks up
    their protein sequence divergence (from Ensembl Compara orthologs).
    Lower mean divergence = more conserved markers.

    Biology: If a cell type's defining genes are highly conserved between
    human and mouse, we predict it should also be more rigid in
    Procrustes space (conserved expression geometry).

    Math: mean(seq_divergence) over top-100 marker genes per type,
    where seq_divergence = 100 - percent_identity from Ensembl.
    """
    import scanpy as sc

    print("  Computing marker gene conservation (seq_divergence) ...")

    # Load dN/dS (sequence divergence) data
    dnds = pd.read_csv(dnds_path)
    # Map Ensembl IDs to seq_divergence
    dnds_map = dict(zip(dnds["human_ensembl_id"], dnds["seq_divergence"]))

    # Get Ensembl IDs from adata
    if "feature_name" in adata.var.columns:
        ensembl_ids = adata.var_names.values  # var_names are Ensembl IDs
    else:
        ensembl_ids = adata.var_names.values

    # Normalize for marker detection
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)
    if "feature_name" in adata_norm.var.columns:
        # Keep Ensembl IDs as var_names for lookup
        pass

    sc.tl.rank_genes_groups(adata_norm, groupby="cell_type", method="wilcoxon",
                            use_raw=False, key_added="markers_cons")

    results = {}
    for ct in cell_types:
        result_df = sc.get.rank_genes_groups_df(adata_norm, group=ct,
                                                 key="markers_cons")
        # Top 100 marker genes by score
        top_genes = result_df.head(100)["names"].values
        # Look up divergence
        divs = [dnds_map[g] for g in top_genes if g in dnds_map]
        if len(divs) >= 10:
            results[ct] = {"mean_marker_seq_divergence": np.mean(divs),
                           "n_markers_with_dnds": len(divs)}
        else:
            results[ct] = {"mean_marker_seq_divergence": np.nan,
                           "n_markers_with_dnds": len(divs)}

    del adata_norm
    df = pd.DataFrame(results).T
    df.index.name = "cell_type"
    n_valid = df["mean_marker_seq_divergence"].notna().sum()
    print(f"    Conservation scores for {n_valid}/{len(cell_types)} types")
    return df


def add_germ_layer(df):
    """Add developmental origin (germ layer) as a categorical feature.

    Biology: Mesoderm gives rise to blood/immune, stromal, endothelial cells.
    Endoderm gives rise to gut, liver, pancreas epithelial cells.
    Ectoderm gives rise to skin, neural, mammary epithelial cells.
    """
    df["germ_layer"] = df.index.map(GERM_LAYER)
    # Encode as dummy variables for regression
    df["is_endoderm"] = (df["germ_layer"] == "endoderm").astype(int)
    df["is_ectoderm"] = (df["germ_layer"] == "ectoderm").astype(int)
    return df


def build_feature_table():
    """Build the complete feature table: 35 cell types × all features.

    Merges precomputed features with features computed from the h5ad.
    """
    print("\n" + "=" * 70)
    print("STEP 0: Loading rigidity scores")
    print("=" * 70)
    residuals = load_residuals()
    print(f"  Loaded residuals for {len(residuals)} cell types")
    print(f"  Range: {residuals['rigidity_score'].min():.2f} – "
          f"{residuals['rigidity_score'].max():.2f}")

    print("\n" + "=" * 70)
    print("STEP 1: Loading / computing biological features")
    print("=" * 70)

    # Precomputed features
    print("\n--- Precomputed features ---")
    precomp = load_precomputed_features()

    # h5ad-derived features
    print("\n--- Computing features from expression data ---")
    h5ad_feats, adata, gene_symbols = compute_h5ad_features(HUMAN_H5AD)

    # Marker gene conservation
    print("\n--- Computing marker gene conservation ---")
    cell_types = sorted(adata.obs["cell_type"].unique())
    conservation = compute_marker_conservation(adata, gene_symbols,
                                                cell_types, DNDS_CSV)
    del adata  # Free memory

    # Merge everything
    print("\n--- Merging all features ---")
    feat = residuals.join(precomp, how="left")
    feat = feat.join(h5ad_feats, how="left")
    feat = feat.join(conservation, how="left")
    feat = add_germ_layer(feat)

    # Add derived features
    feat["mean_hk_ratio"] = (feat["hk_ratio_human"] + feat["hk_ratio_mouse"]) / 2
    feat["mean_n_active_tfs"] = (feat["n_active_tfs_human"] + feat["n_active_tfs_mouse"]) / 2
    feat["mean_entropy"] = (feat["entropy_human"] + feat["entropy_mouse"]) / 2
    feat["log_min_count"] = np.log10(feat["min_count"])

    # Add lineage from figure_style
    feat["lineage"] = feat.index.map(LINEAGE_MAP)

    print(f"\n  Final feature table: {feat.shape[0]} cell types × {feat.shape[1]} columns")
    print(f"  Features with missing values:")
    for col in feat.columns:
        n_miss = feat[col].isna().sum()
        if n_miss > 0:
            print(f"    {col}: {n_miss} missing")

    return feat


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2: Univariate correlations
# ═══════════════════════════════════════════════════════════════════════════════

# Features to test (display name → column name)
UNIVARIATE_FEATURES = {
    "Transcriptomic complexity": "n_genes_detected",
    "Mean expression level": "avg_mean_expr",
    "Cell cycle fraction": "cell_cycle_fraction",
    "Within-type heterogeneity": "within_type_heterogeneity",
    "Tissue breadth": "tissue_breadth",
    "N marker genes": "n_marker_genes",
    "Mean HK ratio": "mean_hk_ratio",
    "Mean active TFs": "mean_n_active_tfs",
    "TF entropy": "mean_entropy",
    "Inter-donor variance": "mean_interdonor_var",
    "Chromatin conservation": "mean_phastCons",
    "Marker seq divergence": "mean_marker_seq_divergence",
    "Log min cell count": "log_min_count",
    "Is progenitor": "is_progenitor",
    "Is endoderm": "is_endoderm",
}


def run_univariate(feat):
    """Compute Spearman ρ between each feature and rigidity score.

    Math: Spearman rank correlation — nonparametric, robust to outliers.
    Two-sided test for H0: ρ = 0.
    """
    print("\n" + "=" * 70)
    print("STEP 2: Univariate correlations")
    print("=" * 70)

    y = feat["rigidity_score"]
    rows = []

    for name, col in UNIVARIATE_FEATURES.items():
        if col not in feat.columns:
            print(f"  SKIP: {name} ({col}) — not in feature table")
            continue
        x = feat[col].dropna()
        common = y.index.intersection(x.index)
        if len(common) < 10:
            print(f"  SKIP: {name} — only {len(common)} non-NA values")
            continue

        rho, p = stats.spearmanr(x.loc[common], y.loc[common])
        rows.append({
            "feature": name,
            "column": col,
            "spearman_rho": rho,
            "p_value": p,
            "n": len(common),
            "abs_rho": abs(rho),
        })

    result = pd.DataFrame(rows).sort_values("abs_rho", ascending=False)
    result = result.reset_index(drop=True)

    # Print table
    print(f"\n{'Feature':<30} {'ρ':>8} {'p':>12} {'n':>4}  Interpretation")
    print("-" * 80)
    for _, row in result.iterrows():
        interp = ""
        if abs(row["spearman_rho"]) > 0.4:
            interp = "*** STRONG"
        elif abs(row["spearman_rho"]) > 0.3:
            interp = "** MODERATE"
        elif abs(row["spearman_rho"]) > 0.2:
            interp = "* WEAK"

        # Flag confounds
        if row["column"] in ("log_min_count", "tissue_breadth") and abs(row["spearman_rho"]) > 0.3:
            interp += " ⚠️ POTENTIAL CONFOUND"

        print(f"  {row['feature']:<28} {row['spearman_rho']:>8.3f} "
              f"{row['p_value']:>12.4e} {row['n']:>4}  {interp}")

    # Decision gate
    max_rho = result["abs_rho"].max()
    print(f"\n  Max |ρ| = {max_rho:.3f}")
    if max_rho < 0.3:
        print("  DECISION: No feature correlates strongly (|ρ| < 0.3).")
        print("  → Rigidity is NOT driven by standard biological covariates.")
        print("  → This rules out confounds — report in Supplementary.")
    elif max_rho >= 0.4:
        print("  DECISION: Strong correlation(s) found (|ρ| ≥ 0.4).")
        top = result.iloc[0]
        print(f"  → Top predictor: {top['feature']} (ρ = {top['spearman_rho']:.3f})")
        if top["column"] in ("log_min_count", "tissue_breadth"):
            print("  ⚠️ WARNING: Top predictor is a POTENTIAL CONFOUND!")
            print("  → This may indicate sampling artifacts, not biology.")
    else:
        print("  DECISION: Moderate correlation(s) found (0.3 ≤ |ρ| < 0.4).")

    return result


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3: Multivariate analysis
# ═══════════════════════════════════════════════════════════════════════════════
def run_multivariate(feat, univariate_result):
    """Fit multivariate models: elastic net, random forest, LOO-CV.

    With n=35 and many predictors, we use regularization (elastic net)
    and tree-based methods (random forest) to avoid overfitting.
    LOO-CV gives an honest estimate of predictive R².

    Math:
    - Elastic net: min ||y - Xβ||² + α(λ||β||₁ + (1-λ)||β||²₂)
      with α and λ chosen by 5-fold CV.
    - Random forest: 100 trees, permutation importance.
    - LOO-CV R²: for each observation, fit on remaining 34, predict held-out.
    """
    print("\n" + "=" * 70)
    print("STEP 3: Multivariate analysis")
    print("=" * 70)

    y = feat["rigidity_score"]

    # Select numeric features (exclude identifiers, categoricals)
    feature_cols = [col for _, col in UNIVARIATE_FEATURES.items()
                    if col in feat.columns]
    # Drop features with too many missing values
    feature_cols = [c for c in feature_cols if feat[c].notna().sum() >= 30]

    X = feat[feature_cols].copy()
    # Impute remaining NAs with median
    for c in X.columns:
        if X[c].isna().any():
            X[c] = X[c].fillna(X[c].median())

    print(f"  Predictors: {len(feature_cols)}")
    for c in feature_cols:
        print(f"    {c}")

    # Standardize
    scaler = StandardScaler()
    X_scaled = pd.DataFrame(scaler.fit_transform(X), index=X.index,
                             columns=X.columns)

    results = {}

    # --- Elastic Net with CV ---
    print("\n  --- Elastic Net (5-fold CV for alpha) ---")
    enet = ElasticNetCV(l1_ratio=[0.1, 0.3, 0.5, 0.7, 0.9],
                         cv=5, random_state=SEED, max_iter=10000)
    enet.fit(X_scaled, y)
    y_pred_enet = enet.predict(X_scaled)
    r2_train = r2_score(y, y_pred_enet)
    print(f"    Training R² = {r2_train:.3f}")
    print(f"    Best alpha = {enet.alpha_:.4f}, l1_ratio = {enet.l1_ratio_:.2f}")

    # Coefficients
    coef_df = pd.DataFrame({
        "feature": feature_cols,
        "coefficient": enet.coef_,
        "abs_coefficient": np.abs(enet.coef_),
    }).sort_values("abs_coefficient", ascending=False)
    print("\n    Elastic Net coefficients (standardized):")
    for _, row in coef_df.iterrows():
        marker = "●" if row["abs_coefficient"] > 0.01 else "○"
        print(f"      {marker} {row['feature']:<35} {row['coefficient']:>8.4f}")

    results["elastic_net"] = {
        "r2_train": r2_train,
        "alpha": float(enet.alpha_),
        "l1_ratio": float(enet.l1_ratio_),
        "coefficients": dict(zip(feature_cols,
                                  [float(c) for c in enet.coef_])),
        "intercept": float(enet.intercept_),
    }

    # --- Random Forest ---
    print("\n  --- Random Forest (100 trees) ---")
    rf = RandomForestRegressor(n_estimators=100, random_state=SEED,
                                max_features="sqrt", min_samples_leaf=3)
    rf.fit(X_scaled, y)
    y_pred_rf = rf.predict(X_scaled)
    r2_rf_train = r2_score(y, y_pred_rf)
    print(f"    Training R² = {r2_rf_train:.3f}")

    importance = pd.DataFrame({
        "feature": feature_cols,
        "importance": rf.feature_importances_,
    }).sort_values("importance", ascending=False)
    print("\n    Feature importance (MDI):")
    for _, row in importance.iterrows():
        bar = "█" * int(row["importance"] * 50)
        print(f"      {row['feature']:<35} {row['importance']:.4f} {bar}")

    results["random_forest"] = {
        "r2_train": r2_rf_train,
        "feature_importance": dict(zip(importance["feature"],
                                        [float(v) for v in importance["importance"]])),
    }

    # --- Leave-One-Out Cross-Validation ---
    print("\n  --- Leave-One-Out Cross-Validation ---")
    loo = LeaveOneOut()
    y_pred_loo_enet = np.zeros(len(y))
    y_pred_loo_rf = np.zeros(len(y))

    for train_idx, test_idx in loo.split(X_scaled):
        X_tr, X_te = X_scaled.iloc[train_idx], X_scaled.iloc[test_idx]
        y_tr = y.iloc[train_idx]

        # Elastic net
        enet_loo = ElasticNetCV(l1_ratio=[0.1, 0.3, 0.5, 0.7, 0.9],
                                 cv=5, random_state=SEED, max_iter=10000)
        enet_loo.fit(X_tr, y_tr)
        y_pred_loo_enet[test_idx] = enet_loo.predict(X_te)

        # Random forest
        rf_loo = RandomForestRegressor(n_estimators=100, random_state=SEED,
                                        max_features="sqrt", min_samples_leaf=3)
        rf_loo.fit(X_tr, y_tr)
        y_pred_loo_rf[test_idx] = rf_loo.predict(X_te)

    r2_loo_enet = r2_score(y, y_pred_loo_enet)
    r2_loo_rf = r2_score(y, y_pred_loo_rf)
    print(f"    Elastic Net LOO-CV R² = {r2_loo_enet:.3f}")
    print(f"    Random Forest LOO-CV R² = {r2_loo_rf:.3f}")

    results["loo_cv"] = {
        "elastic_net_r2": float(r2_loo_enet),
        "random_forest_r2": float(r2_loo_rf),
    }
    results["loo_predictions"] = {
        "cell_types": list(y.index),
        "observed": [float(v) for v in y.values],
        "predicted_enet": [float(v) for v in y_pred_loo_enet],
        "predicted_rf": [float(v) for v in y_pred_loo_rf],
    }

    # Determine best model
    best_model = "elastic_net" if r2_loo_enet > r2_loo_rf else "random_forest"
    best_r2 = max(r2_loo_enet, r2_loo_rf)
    best_preds = y_pred_loo_enet if best_model == "elastic_net" else y_pred_loo_rf
    print(f"\n    Best model (LOO-CV): {best_model} (R² = {best_r2:.3f})")

    if best_r2 < 0.05:
        print("    → Multivariate model has minimal predictive power.")
        print("    → Biological features do not jointly predict rigidity well.")
    elif best_r2 > 0.3:
        print("    → Moderate predictive power — biological features explain some variance.")

    return results, y, y_pred_loo_enet, y_pred_loo_rf, importance, coef_df


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 4: Publication-quality figures
# ═══════════════════════════════════════════════════════════════════════════════
def make_figures(feat, univariate_result, multivariate_results,
                 y_obs, y_pred_enet, y_pred_rf, importance, coef_df):
    """Generate three-panel supplementary figure.

    Panel A: Correlation matrix heatmap — all features × rigidity.
    Panel B: Top 2–3 scatter plots (strongest univariate predictors).
    Panel C: Observed vs predicted rigidity (LOO cross-validated).
    """
    print("\n" + "=" * 70)
    print("STEP 4: Generating figures")
    print("=" * 70)
    apply_style()

    # ─── Panel A: Correlation matrix ─────────────────────────────────────
    print("  Panel A: Correlation matrix heatmap ...")
    cols_for_corr = ["rigidity_score"] + [
        col for _, col in UNIVARIATE_FEATURES.items()
        if col in feat.columns and feat[col].notna().sum() >= 30
    ]
    corr_data = feat[cols_for_corr].dropna(axis=0, how="any")
    corr_matrix = corr_data.corr(method="spearman")

    # Short display names for heatmap
    heatmap_names = {
        "rigidity_score": "Rigidity\nscore",
        "n_genes_detected": "Genes\ndetected",
        "avg_mean_expr": "Mean\nexpression",
        "cell_cycle_fraction": "Cell cycle\nfraction",
        "within_type_heterogeneity": "Within-type\nheterogeneity",
        "tissue_breadth": "Tissue\nbreadth",
        "n_marker_genes": "N marker\ngenes",
        "mean_hk_ratio": "HK gene\nratio",
        "mean_n_active_tfs": "Active\nTFs",
        "mean_entropy": "TF\nentropy",
        "mean_interdonor_var": "Inter-donor\nvariance",
        "mean_phastCons": "Chromatin\nconservation",
        "mean_marker_seq_divergence": "Marker seq\ndivergence",
        "log_min_count": "Log cell\ncount",
        "is_progenitor": "Is\nprogenitor",
        "is_endoderm": "Is\nendoderm",
    }
    display_names = [heatmap_names.get(c, c) for c in corr_matrix.columns]
    corr_display = corr_matrix.copy()
    corr_display.index = display_names
    corr_display.columns = display_names

    fig_a, ax_a = plt.subplots(figsize=(COL2 * 0.65, COL2 * 0.55))
    mask = np.zeros_like(corr_display, dtype=bool)
    # Show full matrix (not just lower triangle) for clarity
    sns.heatmap(corr_display, annot=True, fmt=".2f", cmap="RdBu_r",
                center=0, vmin=-1, vmax=1, square=True,
                linewidths=0.5, linecolor="white",
                annot_kws={"size": 5.5},
                cbar_kws={"shrink": 0.7, "label": "Spearman ρ"},
                ax=ax_a)
    ax_a.tick_params(axis="both", labelsize=6)
    ax_a.set_title("Inter-feature and rigidity correlations", fontsize=8, pad=10)
    add_panel_label(ax_a, "A", x=-0.02, y=1.02)
    save_figure(fig_a, str(FIG_DIR / "biological_predictors_panelA"))

    # ─── Panel B: Top scatter plots ──────────────────────────────────────
    print("  Panel B: Top univariate scatter plots ...")
    # Pick top 3 by |ρ|
    top3 = univariate_result.head(3)
    n_plots = len(top3)

    fig_b, axes_b = plt.subplots(1, n_plots, figsize=(COL2, COL2 * 0.3))
    if n_plots == 1:
        axes_b = [axes_b]

    for i, (_, row) in enumerate(top3.iterrows()):
        ax = axes_b[i]
        col = row["column"]
        x = feat[col]
        y_vals = feat["rigidity_score"]

        # Color by lineage
        colors = [lineage_color(ct) for ct in feat.index]

        ax.scatter(x, y_vals, c=colors, s=25, alpha=0.8, edgecolors="white",
                   linewidth=0.3, zorder=3)

        # Label points
        for ct in feat.index:
            if pd.notna(x.loc[ct]):
                ax.annotate(short_name(ct), (x.loc[ct], y_vals.loc[ct]),
                            fontsize=4, alpha=0.7, ha="center", va="bottom",
                            xytext=(0, 2), textcoords="offset points")

        # Regression line (skip for binary features)
        valid = x.dropna().index.intersection(y_vals.dropna().index)
        if len(valid) > 5 and x.loc[valid].nunique() > 2:
            slope, intercept, _, _, _ = stats.linregress(
                np.asarray(x.loc[valid], dtype=float),
                np.asarray(y_vals.loc[valid], dtype=float))
            x_line = np.linspace(x.loc[valid].min(), x.loc[valid].max(), 100)
            ax.plot(x_line, slope * x_line + intercept, "--",
                    color=C_GRAY, linewidth=0.8, zorder=2)

        ax.set_xlabel(row["feature"], fontsize=FONT_SIZE_LABEL)
        if i == 0:
            ax.set_ylabel("Rigidity score\n(Procrustes residual)", fontsize=FONT_SIZE_LABEL)
        else:
            ax.set_ylabel("")

        # Annotate ρ and p
        rho_str = f"ρ = {row['spearman_rho']:.3f}"
        p_str = format_p(row["p_value"])
        ax.text(0.05, 0.95, f"{rho_str}\n{p_str}",
                transform=ax.transAxes, fontsize=FONT_SIZE_ANNOT,
                va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          edgecolor=C_LIGHTGRAY, alpha=0.9))

        clean_spine(ax)
        if i == 0:
            add_panel_label(ax, "B", x=-0.15, y=1.08)

    fig_b.tight_layout(w_pad=2)
    save_figure(fig_b, str(FIG_DIR / "biological_predictors_panelB"))

    # ─── Panel C: Observed vs predicted (LOO-CV) ────────────────────────
    print("  Panel C: Observed vs predicted (LOO-CV) ...")

    # Use best model
    r2_enet = multivariate_results["loo_cv"]["elastic_net_r2"]
    r2_rf = multivariate_results["loo_cv"]["random_forest_r2"]
    if r2_enet >= r2_rf:
        y_pred = y_pred_enet
        model_name = "Elastic Net"
        r2_cv = r2_enet
    else:
        y_pred = y_pred_rf
        model_name = "Random Forest"
        r2_cv = r2_rf

    fig_c, ax_c = plt.subplots(figsize=(COL1, COL1))
    colors = [lineage_color(ct) for ct in y_obs.index]
    ax_c.scatter(y_obs, y_pred, c=colors, s=30, alpha=0.8,
                 edgecolors="white", linewidth=0.3, zorder=3)

    # Label points
    for ct in y_obs.index:
        ax_c.annotate(short_name(ct), (y_obs.loc[ct], y_pred[y_obs.index.get_loc(ct)]),
                       fontsize=4, alpha=0.7, ha="center", va="bottom",
                       xytext=(0, 2), textcoords="offset points")

    # Diagonal
    lo = min(y_obs.min(), min(y_pred)) - 0.5
    hi = max(y_obs.max(), max(y_pred)) + 0.5
    ax_c.plot([lo, hi], [lo, hi], "--", color=C_LIGHTGRAY, linewidth=0.8, zorder=1)
    ax_c.set_xlim(lo, hi)
    ax_c.set_ylim(lo, hi)

    ax_c.set_xlabel("Observed rigidity score", fontsize=FONT_SIZE_LABEL)
    ax_c.set_ylabel(f"Predicted ({model_name}, LOO-CV)", fontsize=FONT_SIZE_LABEL)
    ax_c.text(0.05, 0.95, f"LOO-CV R² = {r2_cv:.3f}",
              transform=ax_c.transAxes, fontsize=FONT_SIZE_ANNOT,
              va="top", ha="left",
              bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                        edgecolor=C_LIGHTGRAY, alpha=0.9))
    add_panel_label(ax_c, "C", x=-0.12, y=1.08)
    add_lineage_legend(ax_c, loc="lower right")
    clean_spine(ax_c)
    save_figure(fig_c, str(FIG_DIR / "biological_predictors_panelC"))

    # ─── Bonus: Feature importance bar chart ─────────────────────────────
    print("  Bonus: Feature importance bar chart ...")
    fig_d, ax_d = plt.subplots(figsize=(COL1, COL1 * 0.9))
    imp = importance.sort_values("importance", ascending=True)
    display = [heatmap_names.get(f, f).replace("\n", " ") for f in imp["feature"]]
    ax_d.barh(range(len(imp)), imp["importance"], color=C_BLUE, height=0.7)
    ax_d.set_yticks(range(len(imp)))
    ax_d.set_yticklabels(display, fontsize=6)
    ax_d.set_xlabel("Feature importance (MDI)", fontsize=FONT_SIZE_LABEL)
    ax_d.set_title("Random Forest feature importance", fontsize=8)
    clean_spine(ax_d)
    add_panel_label(ax_d, "D", x=-0.25, y=1.02)
    save_figure(fig_d, str(FIG_DIR / "biological_predictors_panelD"))

    print("\n  Figure descriptions:")
    print("    Panel A: Heatmap of Spearman correlations between all features")
    print("             and rigidity score. Blue = positive, red = negative.")
    print("    Panel B: Scatter plots of top 3 univariate predictors vs rigidity.")
    print("             Points colored by lineage. Dashed line = linear fit.")
    print("    Panel C: Observed vs LOO-CV predicted rigidity from best model.")
    print("             Diagonal = perfect prediction.")
    print("    Panel D: Random forest feature importance (mean decrease impurity).")


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 5: Save outputs
# ═══════════════════════════════════════════════════════════════════════════════
def save_outputs(feat, univariate_result, multivariate_results):
    """Save all analysis outputs to analysis/biological_predictors/."""
    print("\n" + "=" * 70)
    print("STEP 5: Saving outputs")
    print("=" * 70)

    # Feature table
    feat.to_csv(OUT_DIR / "feature_table.csv")
    print(f"  Saved: feature_table.csv ({feat.shape[0]} rows × {feat.shape[1]} cols)")

    # Univariate correlations
    univariate_result.to_csv(OUT_DIR / "univariate_correlations.csv", index=False)
    print(f"  Saved: univariate_correlations.csv ({len(univariate_result)} features)")

    # Multivariate results
    with open(OUT_DIR / "multivariate_model_results.json", "w") as f:
        json.dump(multivariate_results, f, indent=2)
    print("  Saved: multivariate_model_results.json")

    # Summary markdown
    write_summary(feat, univariate_result, multivariate_results)


def write_summary(feat, univariate_result, multivariate_results):
    """Write human-readable summary of the biological predictors analysis."""
    max_rho = univariate_result["abs_rho"].max()
    top = univariate_result.iloc[0]
    best_loo = max(multivariate_results["loo_cv"]["elastic_net_r2"],
                   multivariate_results["loo_cv"]["random_forest_r2"])

    # Check for confounds
    confound_cols = {"log_min_count", "tissue_breadth"}
    top_is_confound = top["column"] in confound_cols

    lines = [
        "# Biological Predictors of Cell Type Rigidity",
        "",
        "## Summary",
        "",
        f"Tested {len(UNIVARIATE_FEATURES)} biological features as predictors of",
        f"Procrustes rigidity across {len(feat)} cell types.",
        "",
        "## Univariate Results (Spearman ρ with rigidity score)",
        "",
        "| Feature | ρ | p-value | Interpretation |",
        "|---------|---|---------|----------------|",
    ]
    for _, row in univariate_result.iterrows():
        interp = ""
        if row["abs_rho"] > 0.4:
            interp = "Strong"
        elif row["abs_rho"] > 0.3:
            interp = "Moderate"
        elif row["abs_rho"] > 0.2:
            interp = "Weak"
        else:
            interp = "Negligible"
        if row["column"] in confound_cols and row["abs_rho"] > 0.3:
            interp += " ⚠️ CONFOUND"
        lines.append(f"| {row['feature']} | {row['spearman_rho']:.3f} | "
                     f"{row['p_value']:.4e} | {interp} |")

    lines += [
        "",
        f"**Strongest univariate predictor:** {top['feature']} "
        f"(ρ = {top['spearman_rho']:.3f}, p = {top['p_value']:.4e})",
        "",
    ]

    if top_is_confound:
        lines += [
            "⚠️ **WARNING:** The strongest predictor is a potential CONFOUND.",
            f"{top['feature']} may reflect sampling artifacts rather than biology.",
            "",
        ]

    lines += [
        "## Multivariate Results",
        "",
        f"- Elastic Net training R² = "
        f"{multivariate_results['elastic_net']['r2_train']:.3f}",
        f"- Random Forest training R² = "
        f"{multivariate_results['random_forest']['r2_train']:.3f}",
        f"- Elastic Net LOO-CV R² = "
        f"{multivariate_results['loo_cv']['elastic_net_r2']:.3f}",
        f"- Random Forest LOO-CV R² = "
        f"{multivariate_results['loo_cv']['random_forest_r2']:.3f}",
        "",
    ]

    if best_loo < 0.05:
        lines += [
            "**Conclusion:** Multivariate models have minimal predictive power.",
            "Biological features do not jointly predict rigidity well.",
            "This is informative: rigidity rankings are NOT explained by",
            "standard biological covariates (complexity, cell cycle, tissue",
            "breadth, etc.), ruling them out as confounds.",
        ]
    elif best_loo > 0.3:
        lines += [
            "**Conclusion:** Moderate predictive power from biological features.",
            "Some variance in rigidity is explained by measurable biological",
            "properties of cell types.",
        ]
    else:
        lines += [
            "**Conclusion:** Weak predictive power from biological features.",
            "Rigidity rankings are only partially explained by standard",
            "biological covariates.",
        ]

    lines += [
        "",
        "## Decision Gate",
        "",
    ]
    if max_rho < 0.3:
        lines += [
            "No feature correlates with rigidity at |ρ| > 0.3.",
            "This analysis shows that rigidity is not easily explained by",
            "standard biological covariates. Report in Supplementary as:",
            "\"rigidity is not driven by transcriptomic complexity, cell cycle,",
            "or tissue breadth.\"",
        ]
    elif max_rho >= 0.4:
        lines += [
            f"Strong correlation found: {top['feature']} (|ρ| = {top['abs_rho']:.3f}).",
            "This is a meaningful biological finding worth reporting.",
        ]
    else:
        lines += [
            f"Moderate correlation found: {top['feature']} (|ρ| = {top['abs_rho']:.3f}).",
            "Marginal finding — consider for Supplementary.",
        ]

    lines += [
        "",
        "## Germ Layer Distribution",
        "",
        f"- Mesoderm: {(feat['germ_layer'] == 'mesoderm').sum()} types",
        f"- Endoderm: {(feat['germ_layer'] == 'endoderm').sum()} types",
        f"- Ectoderm: {(feat['germ_layer'] == 'ectoderm').sum()} types",
        "",
        "Note: Strong mesoderm bias (25/35 types) limits germ layer analysis.",
        "",
        "## Files Generated",
        "",
        "- `feature_table.csv` — 35 rows × all features + rigidity score",
        "- `univariate_correlations.csv` — all Spearman correlations",
        "- `multivariate_model_results.json` — model coefficients, R², predictions",
        "- `figures/supplementary/biological_predictors_panel{A,B,C,D}.{pdf,png}`",
    ]

    text = "\n".join(lines) + "\n"
    (OUT_DIR / "biological_predictors_summary.md").write_text(text)
    print("  Saved: biological_predictors_summary.md")


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════
def main():
    print("=" * 70)
    print("CellWarp — Biological Predictors of Cell Type Rigidity")
    print("=" * 70)

    # Build feature table (Steps 0–1)
    feat = build_feature_table()

    # Save intermediate feature table
    feat.to_csv(OUT_DIR / "feature_table.csv")

    # Univariate (Step 2)
    univariate_result = run_univariate(feat)

    # Multivariate (Step 3)
    (multivariate_results, y_obs, y_pred_enet, y_pred_rf,
     importance, coef_df) = run_multivariate(feat, univariate_result)

    # Figures (Step 4)
    make_figures(feat, univariate_result, multivariate_results,
                 y_obs, y_pred_enet, y_pred_rf, importance, coef_df)

    # Save outputs (Step 5)
    save_outputs(feat, univariate_result, multivariate_results)

    print("\n" + "=" * 70)
    print("DONE — Biological Predictors Analysis Complete")
    print("=" * 70)


if __name__ == "__main__":
    main()
