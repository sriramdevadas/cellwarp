#!/usr/bin/env python3
"""
Cancer CNV Diagnostic (T2-A Validation)
========================================

Tests whether cancer Procrustes deformation loadings cluster on chromosomal
arms with known recurrent copy number variations (CNVs) in colon adenocarcinoma.

Biology
-------
Colorectal cancer (CRC) has well-characterized recurrent CNVs:
  - GAINS: chr7 (whole), chr8q, chr13 (whole), chr20q
  - LOSSES: chr17p, chr18q

If Procrustes deformation loadings are enriched on these arms, the geometric
signal may partially reflect CNV dosage effects rather than genuine cell
identity transformation. This is especially critical for epithelial cells
(the tumor cell population), which carry the CNVs.

Math
----
For each cell type i:
  1. Residual vector in PCA space: r_i (dim = n_components)
  2. Project to gene space: g_i = r_i @ W, where W = PCA components (n_components x n_genes)
  3. Absolute loading per gene: |g_ij|
  4. Group genes by chromosomal arm
  5. Mann-Whitney U test: |g| on CNV arms vs non-CNV arms
  6. Effect size: median(|g| on CNV arms) / median(|g| on non-CNV arms)

Verdict thresholds:
  - CLEAN: p > 0.05 AND ratio < 1.5
  - POSSIBLY CONFOUNDED: p in [0.05, 0.15) OR ratio in [1.5, 2.0)
  - CONFOUNDED: p < 0.05 AND ratio > 2.0

Inputs:
    output/cancer/scaled/cancer_scaled_results.json
    output/cancer/scaled/pca_cancer_scaled.npz
    output/cancer/scaled/centroids_normal_scaled.csv
    data/phase1/orthologs_human_mouse.csv
    pybiomart (Ensembl BioMart) for chromosomal positions

Outputs (all in output/cancer/cnv_diagnostic/):
    gene_chr_annotations.csv          — gene-to-chromosome arm mapping
    cnv_diagnostic_results.json       — full statistical results
    manhattan_{cell_type}.png         — per-cell-type Manhattan plots
    manhattan_combined.png            — mean loading across all cell types
    cnv_diagnostic_summary.txt        — human-readable verdict report
"""

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from scipy import stats

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
RESULTS_JSON = PROJECT_ROOT / "output" / "cancer" / "scaled" / "cancer_scaled_results.json"
PCA_NPZ = PROJECT_ROOT / "output" / "cancer" / "scaled" / "pca_cancer_scaled.npz"
CENTROIDS_CSV = PROJECT_ROOT / "output" / "cancer" / "scaled" / "centroids_normal_scaled.csv"
ORTHOLOGS_CSV = PROJECT_ROOT / "data" / "phase1" / "orthologs_human_mouse.csv"
OUTPUT_DIR = PROJECT_ROOT / "output" / "cancer" / "cnv_diagnostic"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Known CRC CNV regions
# ---------------------------------------------------------------------------
# Arms with recurrent gains or losses in colon adenocarcinoma
CRC_CNV_ARMS = {"7p", "7q", "8q", "13q", "17p", "18q", "20q"}
# chr7 and chr13 have whole-chromosome gains, so both arms included

# Human centromere positions (GRCh38, approximate midpoints in bp)
# Sources: UCSC Genome Browser centromere annotations
CENTROMERE_POS = {
    "1": 123_400_000,
    "2": 93_900_000,
    "3": 90_900_000,
    "4": 50_000_000,
    "5": 48_800_000,
    "6": 58_100_000,
    "7": 60_100_000,
    "8": 45_200_000,
    "9": 43_000_000,
    "10": 39_800_000,
    "11": 53_400_000,
    "12": 35_500_000,
    "13": 17_700_000,
    "14": 17_200_000,
    "15": 19_000_000,
    "16": 36_800_000,
    "17": 25_100_000,
    "18": 18_500_000,
    "19": 26_200_000,
    "20": 28_100_000,
    "21": 12_000_000,
    "22": 15_000_000,
    "X": 61_000_000,
    "Y": 10_400_000,
}


# ---------------------------------------------------------------------------
# Step 1: Reconstruct full gene-space deformation loadings
# ---------------------------------------------------------------------------
def step1_gene_loadings():
    """
    Load cancer Procrustes results and project residual vectors from PCA
    space back to full gene space (17,180 genes x 14 cell types).

    Returns:
        loadings_df: DataFrame with columns [ensembl_gene_id, gene_name, <cell_types>...]
        cell_types: list of cell type names
    """
    print("=" * 70)
    print("STEP 1: Reconstruct gene-space deformation loadings")
    print("=" * 70)

    with open(RESULTS_JSON) as f:
        results = json.load(f)

    pca_data = np.load(PCA_NPZ)
    W = pca_data["components"]  # (n_components, n_genes)

    # Get gene order from centroids CSV header
    centroids = pd.read_csv(CENTROIDS_CSV, nrows=0)
    gene_ensembl_ids = [c for c in centroids.columns if c.startswith("ENSG")]

    cell_types = results["cell_types"]
    print(f"  {len(cell_types)} cell types, {len(gene_ensembl_ids):,} genes, "
          f"{W.shape[0]} PCA components")

    # Map Ensembl IDs to gene names
    orthologs = pd.read_csv(ORTHOLOGS_CSV)
    id_to_name = dict(zip(orthologs["human_ensembl_id"], orthologs["human_gene_name"]))

    loading_data = {
        "ensembl_gene_id": gene_ensembl_ids,
        "gene_name": [id_to_name.get(g, g) for g in gene_ensembl_ids],
    }

    for ct in cell_types:
        residual_pca = np.array(results["residuals"][ct]["vector_pca"])
        gene_loadings = residual_pca @ W  # (n_genes,)
        loading_data[ct] = np.abs(gene_loadings)

    loadings_df = pd.DataFrame(loading_data)
    print(f"  Loadings shape: {loadings_df.shape}")

    return loadings_df, cell_types


# ---------------------------------------------------------------------------
# Step 2: Get chromosomal positions via pybiomart
# ---------------------------------------------------------------------------
def step2_chr_annotations(gene_ids: list[str]) -> pd.DataFrame:
    """
    Query Ensembl BioMart for chromosome and start position of each human gene.
    Map to chromosomal arm (p or q) using centromere positions.

    Args:
        gene_ids: List of Ensembl gene IDs (ENSG...).

    Returns:
        DataFrame with columns [ensembl_gene_id, chromosome, start_position, arm, chr_arm]
    """
    print("\n" + "=" * 70)
    print("STEP 2: Fetch chromosomal positions from Ensembl BioMart")
    print("=" * 70)

    cache_file = OUTPUT_DIR / "gene_chr_annotations.csv"
    if cache_file.exists():
        print(f"  Loading cached annotations from {cache_file}")
        chr_df = pd.read_csv(cache_file, dtype={"chromosome": str})
        matched = chr_df["ensembl_gene_id"].isin(gene_ids).sum()
        print(f"  Cached: {len(chr_df):,} genes, {matched:,} match our gene set")
        if matched >= len(gene_ids) * 0.9:
            return chr_df
        print("  Cache insufficient, re-querying...")

    from pybiomart import Dataset

    dataset = Dataset(
        name="hsapiens_gene_ensembl",
        host="http://www.ensembl.org",
    )

    print("  Querying Ensembl BioMart for gene coordinates (all human genes)...")
    print("  (This may take 1-2 minutes)")
    result = dataset.query(
        attributes=[
            "ensembl_gene_id",
            "chromosome_name",
            "start_position",
        ],
    )
    print(f"  BioMart returned {len(result):,} rows, columns: {list(result.columns)}")
    # pybiomart returns display-name columns; rename to standard names
    result.columns = ["ensembl_gene_id", "chromosome", "start_position"]
    # Filter to our gene set locally
    result = result[result["ensembl_gene_id"].isin(gene_ids)].copy()
    result["chromosome"] = result["chromosome"].astype(str)

    # Keep only standard chromosomes (1-22, X, Y)
    valid_chr = {str(i) for i in range(1, 23)} | {"X", "Y"}
    result = result[result["chromosome"].isin(valid_chr)].copy()

    # Assign arm (p = before centromere, q = after)
    def assign_arm(row):
        chrom = str(row["chromosome"])
        pos = row["start_position"]
        centro = CENTROMERE_POS.get(chrom)
        if centro is None:
            return "?"
        return "p" if pos < centro else "q"

    result["arm"] = result.apply(assign_arm, axis=1)
    result["chr_arm"] = result["chromosome"] + result["arm"]

    # De-duplicate (some Ensembl IDs map to multiple entries; keep first)
    result = result.drop_duplicates(subset="ensembl_gene_id", keep="first")

    result.to_csv(cache_file, index=False)
    print(f"  Fetched {len(result):,} genes on standard chromosomes")
    print(f"  Coverage: {len(result):,}/{len(gene_ids):,} "
          f"({100 * len(result) / len(gene_ids):.1f}%)")

    # Summary per chromosome
    chr_counts = result["chromosome"].value_counts()
    print(f"  Chromosomes: {len(chr_counts)} (1-22 + X)")

    return result


# ---------------------------------------------------------------------------
# Step 3: CNV enrichment test per cell type
# ---------------------------------------------------------------------------
def step3_cnv_enrichment(loadings_df: pd.DataFrame, chr_df: pd.DataFrame,
                         cell_types: list[str]) -> dict:
    """
    For each cell type, test whether deformation loadings are enriched on
    known CRC CNV chromosomal arms vs non-CNV arms.

    Uses Mann-Whitney U test (non-parametric, no distributional assumptions).

    Args:
        loadings_df: Gene-level absolute deformation loadings.
        chr_df: Gene-to-chromosome arm mapping.
        cell_types: Cell type names.

    Returns:
        Dict with per-cell-type and overall results.
    """
    print("\n" + "=" * 70)
    print("STEP 3: CNV enrichment test — CRC CNV arms vs non-CNV arms")
    print("=" * 70)

    # Merge loadings with chromosomal annotations
    merged = pd.merge(
        loadings_df, chr_df[["ensembl_gene_id", "chromosome", "start_position", "chr_arm"]],
        on="ensembl_gene_id", how="inner",
    )
    print(f"  Merged: {len(merged):,} genes with both loadings and chr position")

    merged["is_cnv"] = merged["chr_arm"].isin(CRC_CNV_ARMS)
    n_cnv = merged["is_cnv"].sum()
    n_non = (~merged["is_cnv"]).sum()
    print(f"  CNV arm genes: {n_cnv:,}, non-CNV arm genes: {n_non:,}")

    results = {}
    verdicts = {}

    print(f"\n  {'Cell Type':<45} {'U-stat':>10} {'p-value':>10} "
          f"{'Ratio':>8} {'Verdict'}")
    print(f"  {'-' * 45} {'-' * 10} {'-' * 10} {'-' * 8} {'-' * 20}")

    for ct in cell_types:
        cnv_loadings = merged.loc[merged["is_cnv"], ct].values
        non_loadings = merged.loc[~merged["is_cnv"], ct].values

        # Mann-Whitney U: test if CNV arm loadings are greater
        stat, p_two = stats.mannwhitneyu(cnv_loadings, non_loadings, alternative="greater")
        p_value = p_two  # one-sided (greater)

        median_cnv = np.median(cnv_loadings)
        median_non = np.median(non_loadings)
        ratio = median_cnv / median_non if median_non > 0 else float("inf")

        mean_cnv = np.mean(cnv_loadings)
        mean_non = np.mean(non_loadings)
        mean_ratio = mean_cnv / mean_non if mean_non > 0 else float("inf")

        # Verdict
        if p_value < 0.05 and ratio > 2.0:
            verdict = "CONFOUNDED"
        elif p_value < 0.15 or ratio > 1.5:
            verdict = "POSSIBLY CONFOUNDED"
        else:
            verdict = "CLEAN"

        verdicts[ct] = verdict
        results[ct] = {
            "u_statistic": float(stat),
            "p_value": float(p_value),
            "median_cnv": float(median_cnv),
            "median_non_cnv": float(median_non),
            "median_ratio": float(ratio),
            "mean_cnv": float(mean_cnv),
            "mean_non_cnv": float(mean_non),
            "mean_ratio": float(mean_ratio),
            "n_cnv_genes": int(len(cnv_loadings)),
            "n_non_cnv_genes": int(len(non_loadings)),
            "verdict": verdict,
        }

        print(f"  {ct:<45} {stat:>10.0f} {p_value:>10.4f} {ratio:>8.3f} {verdict}")

    # Combined (mean loading across all cell types)
    ct_cols = cell_types
    merged["mean_loading"] = merged[ct_cols].mean(axis=1)
    cnv_mean = merged.loc[merged["is_cnv"], "mean_loading"].values
    non_mean = merged.loc[~merged["is_cnv"], "mean_loading"].values
    stat_comb, p_comb = stats.mannwhitneyu(cnv_mean, non_mean, alternative="greater")
    ratio_comb = np.median(cnv_mean) / np.median(non_mean)

    if p_comb < 0.05 and ratio_comb > 2.0:
        verdict_comb = "CONFOUNDED"
    elif p_comb < 0.15 or ratio_comb > 1.5:
        verdict_comb = "POSSIBLY CONFOUNDED"
    else:
        verdict_comb = "CLEAN"

    print(f"\n  {'COMBINED (all cell types)':<45} {stat_comb:>10.0f} {p_comb:>10.4f} "
          f"{ratio_comb:>8.3f} {verdict_comb}")

    results["_combined"] = {
        "u_statistic": float(stat_comb),
        "p_value": float(p_comb),
        "median_cnv": float(np.median(cnv_mean)),
        "median_non_cnv": float(np.median(non_mean)),
        "median_ratio": float(ratio_comb),
        "mean_cnv": float(np.mean(cnv_mean)),
        "mean_non_cnv": float(np.mean(non_mean)),
        "mean_ratio": float(np.mean(cnv_mean) / np.mean(non_mean)),
        "verdict": verdict_comb,
    }

    # Per-arm breakdown
    print("\n  Per-arm mean absolute loading (averaged across all cell types):")
    arm_means = merged.groupby("chr_arm")["mean_loading"].agg(["mean", "count"])
    arm_means = arm_means.sort_values("mean", ascending=False)
    print(f"  {'Arm':<8} {'Mean |loading|':>15} {'N genes':>10} {'CNV?':>6}")
    print(f"  {'-' * 8} {'-' * 15} {'-' * 10} {'-' * 6}")
    for arm_name, row in arm_means.iterrows():
        is_cnv_flag = "*" if arm_name in CRC_CNV_ARMS else ""
        print(f"  {arm_name:<8} {row['mean']:>15.6f} {int(row['count']):>10} {is_cnv_flag:>6}")

    results["_per_arm"] = {
        arm: {"mean_loading": float(row["mean"]), "n_genes": int(row["count"]),
              "is_cnv": arm in CRC_CNV_ARMS}
        for arm, row in arm_means.iterrows()
    }

    # Overall verdict
    n_confounded = sum(1 for v in verdicts.values() if v == "CONFOUNDED")
    n_possibly = sum(1 for v in verdicts.values() if v == "POSSIBLY CONFOUNDED")
    n_clean = sum(1 for v in verdicts.values() if v == "CLEAN")

    if n_confounded >= 3:
        overall = ("SYSTEMATIC CNV CONFOUND — ≥3 cell types confounded. "
                   "Cancer deformation scores need CNV-corrected reanalysis.")
    elif n_confounded <= 2:
        overall = (f"LARGELY CLEAN — {n_confounded} confounded, {n_possibly} possibly "
                   f"confounded, {n_clean} clean. CNV is not the dominant signal.")
    else:
        overall = "INCONCLUSIVE"

    results["_overall"] = {
        "n_confounded": n_confounded,
        "n_possibly_confounded": n_possibly,
        "n_clean": n_clean,
        "verdict": overall,
    }

    # Epithelial-specific highlight
    epi_key = "epithelial cell"
    epi_res = results.get(epi_key, {})
    print(f"\n  *** EPITHELIAL CELL (critical test — highest deformation, 35.5% SSR) ***")
    print(f"      p-value: {epi_res.get('p_value', 'N/A'):.4f}")
    print(f"      Median ratio (CNV/non-CNV): {epi_res.get('median_ratio', 'N/A'):.3f}")
    print(f"      Verdict: {epi_res.get('verdict', 'N/A')}")

    print(f"\n  OVERALL VERDICT: {overall}")

    return results, merged


# ---------------------------------------------------------------------------
# Step 4: Manhattan-style plots
# ---------------------------------------------------------------------------
def step4_manhattan_plots(merged: pd.DataFrame, cell_types: list[str],
                          results: dict):
    """
    Generate Manhattan-style plots showing mean absolute deformation loading
    by genomic position, with CRC CNV arms highlighted in red.

    One plot per cell type + one combined plot.
    """
    print("\n" + "=" * 70)
    print("STEP 4: Manhattan-style visualization")
    print("=" * 70)

    # Sort by chromosome (numeric order) then position
    chr_order = [str(i) for i in range(1, 23)] + ["X"]
    merged = merged[merged["chromosome"].isin(chr_order)].copy()
    merged["chr_num"] = merged["chromosome"].map(
        {c: i for i, c in enumerate(chr_order)}
    )
    merged = merged.sort_values(["chr_num", "start_position"]).copy()

    # Compute cumulative genome position for x-axis
    chr_offsets = {}
    offset = 0
    for chrom in chr_order:
        chr_offsets[chrom] = offset
        chr_data = merged[merged["chromosome"] == chrom]
        if len(chr_data) > 0:
            offset += chr_data["start_position"].max() + 5_000_000  # gap

    merged["genome_pos"] = merged.apply(
        lambda r: chr_offsets.get(r["chromosome"], 0) + r["start_position"], axis=1
    )

    # Chromosome midpoints for x-axis labels
    chr_midpoints = {}
    for chrom in chr_order:
        chr_data = merged[merged["chromosome"] == chrom]
        if len(chr_data) > 0:
            chr_midpoints[chrom] = chr_data["genome_pos"].median()

    def _make_manhattan(ax, x, y, chr_arms, title, ylabel, results_entry=None):
        """Plot a single Manhattan panel."""
        colors = []
        for arm in chr_arms:
            if arm in CRC_CNV_ARMS:
                colors.append("#D62728")  # red for CNV arms
            else:
                colors.append("#1F77B4")  # blue for non-CNV

        ax.scatter(x, y, c=colors, s=1.5, alpha=0.3, rasterized=True)

        # Smoothed line (rolling median in windows)
        window = max(len(x) // 200, 50)
        df_temp = pd.DataFrame({"x": x, "y": y}).sort_values("x")
        rolling_median = df_temp["y"].rolling(window=window, center=True).median()
        ax.plot(df_temp["x"].values, rolling_median.values, color="black",
                linewidth=0.8, alpha=0.7)

        ax.set_xlabel("Chromosome", fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.set_title(title, fontsize=10, fontweight="bold")

        # Chromosome labels
        ax.set_xticks(list(chr_midpoints.values()))
        ax.set_xticklabels(list(chr_midpoints.keys()), fontsize=6, rotation=0)
        ax.tick_params(axis="y", labelsize=7)

        # Add vertical lines between chromosomes
        for chrom in chr_order:
            if chrom in chr_offsets:
                ax.axvline(chr_offsets[chrom], color="gray", linewidth=0.3, alpha=0.3)

        # Annotate stats
        if results_entry:
            p = results_entry.get("p_value", None)
            ratio = results_entry.get("median_ratio", None)
            verdict = results_entry.get("verdict", "")
            color_v = {"CLEAN": "green", "POSSIBLY CONFOUNDED": "orange",
                       "CONFOUNDED": "red"}.get(verdict, "black")
            ax.text(0.02, 0.95, f"p={p:.4f}, ratio={ratio:.2f}\n{verdict}",
                    transform=ax.transAxes, fontsize=7, verticalalignment="top",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                              edgecolor=color_v, alpha=0.9),
                    color=color_v)

        # Legend
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], marker="o", color="w", markerfacecolor="#D62728",
                   markersize=5, label="CRC CNV arm"),
            Line2D([0], [0], marker="o", color="w", markerfacecolor="#1F77B4",
                   markersize=5, label="Non-CNV arm"),
        ]
        ax.legend(handles=legend_elements, loc="upper right", fontsize=6)

    # --- Per cell type plots ---
    for ct in cell_types:
        fig, ax = plt.subplots(figsize=(14, 3.5))
        _make_manhattan(
            ax,
            merged["genome_pos"].values,
            merged[ct].values,
            merged["chr_arm"].values,
            f"Cancer deformation loading — {ct}",
            "|Deformation loading|",
            results.get(ct),
        )
        fig.tight_layout()
        safe_name = ct.replace(" ", "_").replace(",", "").replace("-", "_")
        fname = OUTPUT_DIR / f"manhattan_{safe_name}.png"
        fig.savefig(fname, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {fname.name}")

    # --- Combined plot ---
    fig, ax = plt.subplots(figsize=(14, 3.5))
    _make_manhattan(
        ax,
        merged["genome_pos"].values,
        merged["mean_loading"].values,
        merged["chr_arm"].values,
        "Cancer deformation loading — COMBINED (mean across all cell types)",
        "Mean |deformation loading|",
        results.get("_combined"),
    )
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "manhattan_combined.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: manhattan_combined.png")

    # --- Summary multi-panel (epithelial + combined) ---
    fig, axes = plt.subplots(2, 1, figsize=(14, 7), sharex=True)
    _make_manhattan(
        axes[0],
        merged["genome_pos"].values,
        merged["epithelial cell"].values,
        merged["chr_arm"].values,
        "Epithelial cell (highest deformation — 35.5% SSR)",
        "|Deformation loading|",
        results.get("epithelial cell"),
    )
    _make_manhattan(
        axes[1],
        merged["genome_pos"].values,
        merged["mean_loading"].values,
        merged["chr_arm"].values,
        "Combined (mean across all 14 cell types)",
        "Mean |deformation loading|",
        results.get("_combined"),
    )
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "manhattan_epithelial_vs_combined.png", dpi=150,
                bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: manhattan_epithelial_vs_combined.png")

    print(f"\n  Description: Manhattan-style plots with genomic position on x-axis "
          f"and absolute deformation loading on y-axis. Red dots = genes on known CRC "
          f"CNV arms (7p, 7q, 8q, 13q, 17p, 18q, 20q). Black rolling median line "
          f"shows local trend. If CNV confound is present, red clusters should show "
          f"elevated loading compared to blue (non-CNV) regions.")


# ---------------------------------------------------------------------------
# Step 5: Write summary report
# ---------------------------------------------------------------------------
def step5_summary_report(results: dict, cell_types: list[str]):
    """Write human-readable summary report and JSON results."""
    print("\n" + "=" * 70)
    print("STEP 5: Summary report")
    print("=" * 70)

    # JSON
    json_path = OUTPUT_DIR / "cnv_diagnostic_results.json"
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"  Saved: {json_path.name}")

    # Text report
    lines = []
    lines.append("=" * 70)
    lines.append("CANCER CNV DIAGNOSTIC — T2-A VALIDATION")
    lines.append("=" * 70)
    lines.append("")
    lines.append("Question: Do cancer Procrustes deformation loadings cluster on")
    lines.append("chromosomal arms with known CRC copy number variations?")
    lines.append("")
    lines.append(f"CRC CNV arms tested: {', '.join(sorted(CRC_CNV_ARMS))}")
    lines.append("")
    lines.append("-" * 70)
    lines.append("PER-CELL-TYPE RESULTS")
    lines.append("-" * 70)
    lines.append(f"{'Cell Type':<45} {'p-value':>8} {'Ratio':>8} {'Verdict'}")
    lines.append(f"{'-' * 45} {'-' * 8} {'-' * 8} {'-' * 20}")

    for ct in cell_types:
        r = results[ct]
        lines.append(f"{ct:<45} {r['p_value']:>8.4f} {r['median_ratio']:>8.3f} "
                      f"{r['verdict']}")

    r_comb = results["_combined"]
    lines.append("")
    lines.append(f"{'COMBINED':<45} {r_comb['p_value']:>8.4f} "
                 f"{r_comb['median_ratio']:>8.3f} {r_comb['verdict']}")

    lines.append("")
    lines.append("-" * 70)
    lines.append("EPITHELIAL CELL — CRITICAL TEST")
    lines.append("-" * 70)
    epi = results.get("epithelial cell", {})
    lines.append(f"Epithelial cell has the highest cancer deformation (10.07, 35.5% SSR).")
    lines.append(f"Epithelial tumor cells carry the CNVs.")
    lines.append(f"  p-value (Mann-Whitney, one-sided): {epi.get('p_value', 'N/A'):.4f}")
    lines.append(f"  Median loading ratio (CNV/non-CNV): {epi.get('median_ratio', 'N/A'):.3f}")
    lines.append(f"  Mean loading ratio (CNV/non-CNV):   {epi.get('mean_ratio', 'N/A'):.3f}")
    lines.append(f"  Verdict: {epi.get('verdict', 'N/A')}")

    lines.append("")
    lines.append("-" * 70)
    lines.append("OVERALL VERDICT")
    lines.append("-" * 70)
    ov = results["_overall"]
    lines.append(f"  Confounded:          {ov['n_confounded']}/{len(cell_types)}")
    lines.append(f"  Possibly confounded: {ov['n_possibly_confounded']}/{len(cell_types)}")
    lines.append(f"  Clean:               {ov['n_clean']}/{len(cell_types)}")
    lines.append(f"")
    lines.append(f"  {ov['verdict']}")
    lines.append("")
    lines.append("=" * 70)

    report = "\n".join(lines)

    report_path = OUTPUT_DIR / "cnv_diagnostic_summary.txt"
    with open(report_path, "w") as f:
        f.write(report)
    print(f"  Saved: {report_path.name}")

    # Print to console
    print("\n" + report)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("CANCER CNV DIAGNOSTIC (T2-A VALIDATION)")
    print("Testing CRC copy number variation confound in Procrustes loadings")
    print("=" * 70)

    # Step 1: Gene-space loadings
    loadings_df, cell_types = step1_gene_loadings()

    # Step 2: Chromosomal annotations
    chr_df = step2_chr_annotations(loadings_df["ensembl_gene_id"].tolist())

    # Step 3: CNV enrichment test
    results, merged = step3_cnv_enrichment(loadings_df, chr_df, cell_types)

    # Step 4: Manhattan plots
    step4_manhattan_plots(merged, cell_types, results)

    # Step 5: Summary report
    step5_summary_report(results, cell_types)

    print("\nDone. All outputs in:", OUTPUT_DIR)


if __name__ == "__main__":
    main()
