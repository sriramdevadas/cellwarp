"""
CellWarp — Niche Adaptation Hypothesis Test for Progenitor Divergence

Hypothesis: Progenitor-specific divergence genes are enriched in niche-response
categories (cytokine signaling, growth factor receptors, ECM/adhesion, metabolic
sensors) because progenitors diverge via their species-specific niche interfaces,
not their intrinsic identity programs.

Biology
-------
Stem/progenitor cells reside in specialized microenvironments (niches) that
provide signals controlling self-renewal, quiescence, and differentiation.
If cross-species progenitor divergence is driven by niche adaptation, the genes
most specifically diverged in progenitors (high progenitor_specificity_score)
should be enriched for niche-interface pathways: cytokine receptors, ECM
components, integrins, growth factor receptors, and inflammatory mediators.

Negative controls: intrinsic proliferation (MYC targets) and differentiated
metabolic output (oxidative phosphorylation, bile acid metabolism) should NOT
be enriched in progenitor-specific divergence.

Math
----
Fisher's exact test on 2×2 contingency tables:
    |              | In gene set | Not in gene set |
    | Foreground   |      a      |        b        |
    | Background   |      c      |        d        |

Odds ratio = (a*d) / (b*c). FDR correction via Benjamini-Hochberg.
"""

from __future__ import annotations

import json
from pathlib import Path

import gseapy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

OUTPUT_DIR = Path("output/phase2/progenitor_analysis/niche_hypothesis")
SPECIFICITY_PATH = Path("output/phase2/progenitor_analysis/progenitor_specificity_scores.csv")

# Thresholds
TOP_N = 200
PROGENITOR_THRESHOLD = 0.7  # top specificity scores
# Bottom threshold: computed symmetrically as genes with lowest scores

# Gene set definitions — map short names to (library, full_name)
NICHE_GENE_SETS = {
    "ECM_Organization": ("Reactome_Pathways_2024", "Extracellular Matrix Organization"),
    "Interleukin_Signaling": ("Reactome_Pathways_2024", "Signaling by Interleukins"),
    "Growth_Factor_Receptors": ("Reactome_Pathways_2024",
                                 "Diseases of Signal Transduction by Growth Factor Receptors and Second Messengers"),
    "Inflammatory_Response": ("MSigDB_Hallmark_2020", "Inflammatory Response"),
    "TGF_beta_Signaling": ("MSigDB_Hallmark_2020", "TGF-beta Signaling"),
    "Integrin_Interactions": ("Reactome_Pathways_2024", "Integrin Cell Surface Interactions"),
}

NEGATIVE_CONTROL_SETS = {
    "Myc_Targets_V1": ("MSigDB_Hallmark_2020", "Myc Targets V1"),
}

DIFFERENTIATION_OUTPUT_SETS = {
    "Oxidative_Phosphorylation": ("MSigDB_Hallmark_2020", "Oxidative Phosphorylation"),
    "Bile_Acid_Metabolism": ("MSigDB_Hallmark_2020", "Bile Acid Metabolism"),
}


# ---------------------------------------------------------------------------
# Step 1: Load gene sets from MSigDB via gseapy
# ---------------------------------------------------------------------------

def load_gene_sets() -> dict[str, list[str]]:
    """Fetch all required gene sets from Enrichr libraries."""
    all_sets = {}
    all_defs = {**NICHE_GENE_SETS, **NEGATIVE_CONTROL_SETS, **DIFFERENTIATION_OUTPUT_SETS}

    # Cache libraries to avoid redundant downloads
    lib_cache: dict[str, dict] = {}
    for short_name, (lib_name, full_name) in all_defs.items():
        if lib_name not in lib_cache:
            print(f"  Fetching library: {lib_name}")
            lib_cache[lib_name] = gseapy.get_library(name=lib_name)
        library = lib_cache[lib_name]
        if full_name not in library:
            raise KeyError(f"Gene set '{full_name}' not found in {lib_name}")
        genes = library[full_name]
        all_sets[short_name] = genes
        print(f"    {short_name}: {len(genes)} genes")

    return all_sets


# ---------------------------------------------------------------------------
# Step 2 & 3: Fisher's exact test for enrichment
# ---------------------------------------------------------------------------

def fishers_test(
    foreground: set[str],
    background_all: set[str],
    gene_set: set[str],
) -> dict:
    """
    2x2 Fisher's exact test for enrichment of gene_set in foreground vs background.

    Args:
        foreground: Set of foreground gene symbols.
        background_all: Set of ALL gene symbols in the universe.
        gene_set: Set of gene symbols in the pathway.

    Returns:
        Dict with odds_ratio, p_value, a, b, c, d (contingency table values).
    """
    fg = foreground & background_all  # ensure foreground is subset
    gs = gene_set & background_all    # genes in set that are in our universe

    a = len(fg & gs)       # foreground AND in gene set
    b = len(fg - gs)       # foreground NOT in gene set
    c = len(gs - fg)       # background (not foreground) AND in gene set
    d = len(background_all - fg - gs)  # background NOT in gene set

    oddsratio, p_value = stats.fisher_exact([[a, b], [c, d]], alternative="greater")

    return {
        "overlap": a,
        "fg_not_in_set": b,
        "bg_in_set": c,
        "bg_not_in_set": d,
        "odds_ratio": oddsratio,
        "p_value": p_value,
        "fg_size": len(fg),
        "gene_set_size": len(gs),
    }


def run_enrichment_tests(
    foreground_genes: set[str],
    all_genes: set[str],
    gene_sets: dict[str, list[str]],
    label: str,
) -> pd.DataFrame:
    """
    Run Fisher's exact test for each gene set against foreground.

    Returns DataFrame with results, FDR-corrected q-values.
    """
    results = []
    for set_name, genes in gene_sets.items():
        gene_set = set(genes)
        res = fishers_test(foreground_genes, all_genes, gene_set)
        res["gene_set"] = set_name
        res["category"] = _categorize(set_name)
        results.append(res)

    df = pd.DataFrame(results)
    # BH FDR correction
    from statsmodels.stats.multitest import multipletests
    _, q_values, _, _ = multipletests(df["p_value"], method="fdr_bh")
    df["q_value"] = q_values
    df = df.sort_values("p_value").reset_index(drop=True)

    # Print results
    print(f"\n{'='*75}")
    print(f"  Fisher's exact test: {label}")
    print(f"  Foreground: {len(foreground_genes)} genes")
    print(f"{'='*75}")
    print(f"  {'Gene Set':<35s} {'Overlap':>7s} {'OR':>8s} {'p-value':>12s} {'q-value':>12s} {'Category':<15s}")
    print(f"  {'-'*35} {'-'*7} {'-'*8} {'-'*12} {'-'*12} {'-'*15}")
    for _, row in df.iterrows():
        sig = "***" if row["q_value"] < 0.001 else "**" if row["q_value"] < 0.01 else "*" if row["q_value"] < 0.05 else ""
        print(f"  {row['gene_set']:<35s} {row['overlap']:>4d}/{row['gene_set_size']:<3d} "
              f"{row['odds_ratio']:>7.2f} {row['p_value']:>12.2e} {row['q_value']:>12.2e} {row['category']:<12s} {sig}")

    return df


def _categorize(set_name: str) -> str:
    """Categorize gene sets for summary."""
    if set_name in NICHE_GENE_SETS:
        return "niche"
    elif set_name in NEGATIVE_CONTROL_SETS:
        return "neg_control"
    elif set_name in DIFFERENTIATION_OUTPUT_SETS:
        return "diff_output"
    return "unknown"


# ---------------------------------------------------------------------------
# Step 4: Specificity score distribution within gene sets
# ---------------------------------------------------------------------------

def compute_set_specificity(
    spec_df: pd.DataFrame,
    gene_sets: dict[str, list[str]],
) -> pd.DataFrame:
    """
    For each gene set, compute mean and SEM of progenitor_specificity_score
    for genes in that set (restricted to our 16,959-gene universe).
    """
    genome_mean = spec_df["specificity_score"].mean()
    genome_std = spec_df["specificity_score"].std()

    results = []
    all_symbols = set(spec_df["gene_symbol"].values)

    for set_name, genes in gene_sets.items():
        overlap_genes = set(genes) & all_symbols
        if len(overlap_genes) == 0:
            continue
        scores = spec_df[spec_df["gene_symbol"].isin(overlap_genes)]["specificity_score"]
        mean_score = scores.mean()
        sem = scores.std() / np.sqrt(len(scores))
        # One-sample t-test vs genome background
        t_stat, t_pval = stats.ttest_1samp(scores, genome_mean)

        results.append({
            "gene_set": set_name,
            "n_genes_in_universe": len(overlap_genes),
            "mean_specificity": mean_score,
            "sem": sem,
            "std": scores.std(),
            "genome_mean": genome_mean,
            "delta": mean_score - genome_mean,
            "t_stat": t_stat,
            "t_pval": t_pval,
            "category": _categorize(set_name),
        })

    df = pd.DataFrame(results).sort_values("mean_specificity", ascending=False).reset_index(drop=True)
    return df


def plot_specificity_bars(
    set_spec_df: pd.DataFrame,
    genome_mean: float,
    output_path: Path,
) -> str:
    """Bar chart of mean specificity score per gene set with error bars."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df = set_spec_df.sort_values("mean_specificity", ascending=True).copy()

    # Color by category
    colors = []
    for _, row in df.iterrows():
        if row["category"] == "niche":
            colors.append("#2196F3")  # blue
        elif row["category"] == "neg_control":
            colors.append("#F44336")  # red
        else:
            colors.append("#FF9800")  # orange

    fig, ax = plt.subplots(figsize=(10, 6))
    y_pos = range(len(df))
    bars = ax.barh(y_pos, df["mean_specificity"], xerr=df["sem"],
                   color=colors, edgecolor="black", linewidth=0.5,
                   capsize=3, alpha=0.85)

    ax.set_yticks(y_pos)
    labels = []
    for _, row in df.iterrows():
        sig = ""
        if row["t_pval"] < 0.001:
            sig = " ***"
        elif row["t_pval"] < 0.01:
            sig = " **"
        elif row["t_pval"] < 0.05:
            sig = " *"
        labels.append(f"{row['gene_set']} (n={row['n_genes_in_universe']}){sig}")
    ax.set_yticklabels(labels, fontsize=9)

    # Genome background line
    ax.axvline(genome_mean, color="black", linestyle="--", linewidth=1.5, label=f"Genome mean ({genome_mean:.3f})")
    # Thresholds
    ax.axvline(0.6, color="green", linestyle=":", alpha=0.5, label="Progenitor-enriched (>0.6)")
    ax.axvline(0.4, color="purple", linestyle=":", alpha=0.5, label="Differentiated-enriched (<0.4)")

    ax.set_xlabel("Mean Progenitor Specificity Score", fontsize=11)
    ax.set_title("Niche-Response Gene Set Specificity for Progenitor Divergence", fontsize=12)
    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)

    # Add category legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#2196F3", edgecolor="black", label="Niche-response"),
        Patch(facecolor="#F44336", edgecolor="black", label="Negative control (intrinsic)"),
        Patch(facecolor="#FF9800", edgecolor="black", label="Differentiation output"),
    ]
    ax2 = ax.twinx()
    ax2.set_yticks([])
    ax2.legend(handles=legend_elements, loc="upper right", fontsize=8, framealpha=0.9)

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    # Text description
    niche_sets = df[df["category"] == "niche"]
    above_06 = niche_sets[niche_sets["mean_specificity"] > 0.6]
    desc = (
        f"Bar chart: {len(df)} gene sets. Genome background mean = {genome_mean:.3f}. "
        f"{len(above_06)}/{len(niche_sets)} niche-response sets have mean specificity > 0.6 "
        f"(progenitor-divergence-enriched). "
        f"* p<0.05, ** p<0.01, *** p<0.001 (one-sample t-test vs genome mean)."
    )
    return desc


# ---------------------------------------------------------------------------
# Step 5: Summary
# ---------------------------------------------------------------------------

def write_summary(
    progenitor_results: pd.DataFrame,
    differentiated_results: pd.DataFrame,
    set_spec_df: pd.DataFrame,
    output_dir: Path,
) -> dict:
    """Generate summary statistics and write report."""

    # Niche sets enriched in progenitor divergence
    prog_niche = progenitor_results[progenitor_results["category"] == "niche"]
    prog_niche_sig = prog_niche[prog_niche["q_value"] < 0.05]

    diff_niche = differentiated_results[differentiated_results["category"] == "niche"]
    diff_niche_sig = diff_niche[diff_niche["q_value"] < 0.05]

    # Negative controls
    prog_neg = progenitor_results[progenitor_results["category"] == "neg_control"]
    diff_neg = differentiated_results[differentiated_results["category"] == "neg_control"]

    prog_diff_output = progenitor_results[progenitor_results["category"] == "diff_output"]
    diff_diff_output = differentiated_results[differentiated_results["category"] == "diff_output"]

    # Specificity analysis
    niche_spec = set_spec_df[set_spec_df["category"] == "niche"]
    niche_above_06 = niche_spec[niche_spec["mean_specificity"] > 0.6]

    summary = {
        "hypothesis": "Progenitor divergence is driven by niche adaptation (cytokine/ECM/growth factor signaling)",
        "progenitor_divergence": {
            "n_niche_sets_tested": len(prog_niche),
            "n_niche_sets_significant_q05": len(prog_niche_sig),
            "significant_niche_sets": prog_niche_sig["gene_set"].tolist() if len(prog_niche_sig) > 0 else [],
            "strongest_set": prog_niche.iloc[0]["gene_set"] if len(prog_niche) > 0 else None,
            "strongest_odds_ratio": float(prog_niche.iloc[0]["odds_ratio"]) if len(prog_niche) > 0 else None,
            "strongest_q_value": float(prog_niche.iloc[0]["q_value"]) if len(prog_niche) > 0 else None,
        },
        "differentiated_divergence": {
            "n_niche_sets_significant_q05": len(diff_niche_sig),
            "significant_niche_sets": diff_niche_sig["gene_set"].tolist() if len(diff_niche_sig) > 0 else [],
        },
        "negative_controls": {
            "Myc_Targets_V1_progenitor_q": float(prog_neg.iloc[0]["q_value"]) if len(prog_neg) > 0 else None,
            "Myc_Targets_V1_differentiated_q": float(diff_neg.iloc[0]["q_value"]) if len(diff_neg) > 0 else None,
            "OxPhos_progenitor_q": float(prog_diff_output[prog_diff_output["gene_set"] == "Oxidative_Phosphorylation"]["q_value"].values[0]) if len(prog_diff_output[prog_diff_output["gene_set"] == "Oxidative_Phosphorylation"]) > 0 else None,
            "Bile_Acid_progenitor_q": float(prog_diff_output[prog_diff_output["gene_set"] == "Bile_Acid_Metabolism"]["q_value"].values[0]) if len(prog_diff_output[prog_diff_output["gene_set"] == "Bile_Acid_Metabolism"]) > 0 else None,
        },
        "specificity_analysis": {
            "genome_mean": float(niche_spec.iloc[0]["genome_mean"]) if len(niche_spec) > 0 else None,
            "n_niche_sets_above_06": len(niche_above_06),
            "sets_above_06": niche_above_06["gene_set"].tolist() if len(niche_above_06) > 0 else [],
        },
    }

    # Determine verdict
    niche_enriched_in_prog = len(prog_niche_sig) > 0
    niche_not_in_diff = len(diff_niche_sig) == 0 or len(diff_niche_sig) < len(prog_niche_sig)
    neg_control_clean = all(prog_neg["q_value"] > 0.05) if len(prog_neg) > 0 else True

    if niche_enriched_in_prog and neg_control_clean:
        verdict = "SUPPORTED"
        explanation = (
            f"Niche-response gene sets ARE enriched in progenitor-specific divergence "
            f"({len(prog_niche_sig)}/{len(prog_niche)} significant). "
            f"MYC targets (intrinsic proliferation) are NOT enriched (negative control PASSES). "
        )
    elif niche_enriched_in_prog and not neg_control_clean:
        verdict = "PARTIALLY_SUPPORTED"
        explanation = (
            f"Niche-response gene sets are enriched in progenitor divergence, "
            f"but negative controls also show enrichment — signal may not be niche-specific."
        )
    else:
        verdict = "NOT_SUPPORTED"
        explanation = (
            f"No significant enrichment of niche-response gene sets in progenitor-specific divergence "
            f"({len(prog_niche_sig)}/{len(prog_niche)} significant). "
            f"Evidence does not support the niche adaptation hypothesis."
        )

    # Compare progenitor vs differentiated
    if len(diff_niche_sig) > 0:
        explanation += (
            f" Note: {len(diff_niche_sig)} niche sets also enriched in differentiated divergence — "
            f"niche response may not be progenitor-specific."
        )
    else:
        explanation += (
            f" Niche-response sets NOT enriched in differentiated-cell divergence — "
            f"enrichment is progenitor-specific."
        )

    summary["verdict"] = verdict
    summary["explanation"] = explanation

    # Save
    json_path = output_dir / "niche_hypothesis_results.json"
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)

    return summary


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 75)
    print("  NICHE ADAPTATION HYPOTHESIS TEST")
    print("  Testing whether progenitor divergence genes are enriched")
    print("  in niche-response pathways (ECM, cytokines, growth factors)")
    print("=" * 75)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # -----------------------------------------------------------------------
    # Load data
    # -----------------------------------------------------------------------
    print("\n[1/5] Loading data...")
    spec_df = pd.read_csv(SPECIFICITY_PATH)
    print(f"  Loaded {len(spec_df)} genes with specificity scores")
    print(f"  Score range: [{spec_df['specificity_score'].min():.4f}, {spec_df['specificity_score'].max():.4f}]")
    print(f"  Mean: {spec_df['specificity_score'].mean():.4f}, Median: {spec_df['specificity_score'].median():.4f}")

    all_gene_symbols = set(spec_df["gene_symbol"].values)
    print(f"  Universe size: {len(all_gene_symbols)} genes")

    # Top 200 progenitor-specific divergence genes
    prog_top = spec_df.nlargest(TOP_N, "specificity_score")
    prog_genes = set(prog_top["gene_symbol"].values)
    prog_min_score = prog_top["specificity_score"].min()
    print(f"\n  Progenitor-specific foreground: {len(prog_genes)} genes (score >= {prog_min_score:.4f})")

    # Bottom 200 differentiated-specific divergence genes
    diff_top = spec_df.nsmallest(TOP_N, "specificity_score")
    diff_genes = set(diff_top["gene_symbol"].values)
    diff_max_score = diff_top["specificity_score"].max()
    print(f"  Differentiated-specific foreground: {len(diff_genes)} genes (score <= {diff_max_score:.4f})")

    # -----------------------------------------------------------------------
    # Fetch gene sets
    # -----------------------------------------------------------------------
    print("\n[2/5] Fetching gene sets from MSigDB/Reactome...")
    gene_sets = load_gene_sets()

    # Report overlap with our gene universe
    print("\n  Gene set overlap with our 16,959-gene universe:")
    for name, genes in gene_sets.items():
        overlap = len(set(genes) & all_gene_symbols)
        print(f"    {name}: {overlap}/{len(genes)} genes in universe ({100*overlap/len(genes):.0f}%)")

    # -----------------------------------------------------------------------
    # Step 2: Fisher's tests — progenitor divergence
    # -----------------------------------------------------------------------
    print("\n[3/5] Fisher's exact tests — progenitor-specific divergence genes...")
    prog_results = run_enrichment_tests(prog_genes, all_gene_symbols, gene_sets, "PROGENITOR-specific divergence (top 200)")

    # -----------------------------------------------------------------------
    # Step 3: Fisher's tests — differentiated divergence
    # -----------------------------------------------------------------------
    print("\n[4/5] Fisher's exact tests — differentiated-specific divergence genes...")
    diff_results = run_enrichment_tests(diff_genes, all_gene_symbols, gene_sets, "DIFFERENTIATED-specific divergence (bottom 200)")

    # -----------------------------------------------------------------------
    # Step 4: Specificity distribution within gene sets
    # -----------------------------------------------------------------------
    print("\n[5/5] Specificity score distribution within gene sets...")
    set_spec = compute_set_specificity(spec_df, gene_sets)

    print(f"\n  {'Gene Set':<35s} {'N':>5s} {'Mean':>8s} {'SEM':>8s} {'Delta':>8s} {'t-pval':>12s}")
    print(f"  {'-'*35} {'-'*5} {'-'*8} {'-'*8} {'-'*8} {'-'*12}")
    for _, row in set_spec.iterrows():
        sig = "***" if row["t_pval"] < 0.001 else "**" if row["t_pval"] < 0.01 else "*" if row["t_pval"] < 0.05 else ""
        print(f"  {row['gene_set']:<35s} {row['n_genes_in_universe']:>5d} {row['mean_specificity']:>8.4f} "
              f"{row['sem']:>8.4f} {row['delta']:>+8.4f} {row['t_pval']:>12.2e} {sig}")

    genome_mean = spec_df["specificity_score"].mean()

    # Plot
    desc = plot_specificity_bars(set_spec, genome_mean, OUTPUT_DIR / "specificity_by_gene_set.png")
    print(f"\n  Plot: {desc}")

    # -----------------------------------------------------------------------
    # Save all results
    # -----------------------------------------------------------------------
    prog_results.to_csv(OUTPUT_DIR / "fisher_progenitor_divergence.csv", index=False)
    diff_results.to_csv(OUTPUT_DIR / "fisher_differentiated_divergence.csv", index=False)
    set_spec.to_csv(OUTPUT_DIR / "specificity_distribution.csv", index=False)

    # Overlapping genes detail
    overlap_records = []
    for set_name, genes in gene_sets.items():
        gs = set(genes)
        prog_overlap = prog_genes & gs & all_gene_symbols
        diff_overlap = diff_genes & gs & all_gene_symbols
        for g in prog_overlap:
            s = spec_df[spec_df["gene_symbol"] == g]["specificity_score"].values[0]
            overlap_records.append({"gene": g, "gene_set": set_name, "group": "progenitor", "specificity_score": s})
        for g in diff_overlap:
            s = spec_df[spec_df["gene_symbol"] == g]["specificity_score"].values[0]
            overlap_records.append({"gene": g, "gene_set": set_name, "group": "differentiated", "specificity_score": s})
    if overlap_records:
        pd.DataFrame(overlap_records).to_csv(OUTPUT_DIR / "overlap_gene_details.csv", index=False)

    # -----------------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------------
    summary = write_summary(prog_results, diff_results, set_spec, OUTPUT_DIR)

    print("\n" + "=" * 75)
    print("  VERDICT:", summary["verdict"])
    print("=" * 75)
    print(f"  {summary['explanation']}")

    # Strongest supporting sets
    prog_niche_sig = prog_results[(prog_results["category"] == "niche") & (prog_results["q_value"] < 0.05)]
    if len(prog_niche_sig) > 0:
        print(f"\n  Strongest niche-response sets in progenitor divergence:")
        for _, row in prog_niche_sig.iterrows():
            print(f"    {row['gene_set']}: OR={row['odds_ratio']:.2f}, q={row['q_value']:.2e}, "
                  f"overlap={row['overlap']}/{row['gene_set_size']}")

    # Negative control
    prog_neg = prog_results[prog_results["category"] == "neg_control"]
    if len(prog_neg) > 0:
        row = prog_neg.iloc[0]
        neg_status = "NOT enriched (PASS)" if row["q_value"] > 0.05 else "ENRICHED (FAIL)"
        print(f"\n  Negative control (MYC Targets V1): {neg_status}")
        print(f"    OR={row['odds_ratio']:.2f}, q={row['q_value']:.2e}")

    # Differentiation output controls
    prog_do = prog_results[prog_results["category"] == "diff_output"]
    for _, row in prog_do.iterrows():
        status = "NOT enriched (expected)" if row["q_value"] > 0.05 else "ENRICHED (unexpected)"
        print(f"  Differentiation output ({row['gene_set']}): {status}, q={row['q_value']:.2e}")

    print(f"\n  All results saved to: {OUTPUT_DIR}/")
    print(f"  Files: fisher_progenitor_divergence.csv, fisher_differentiated_divergence.csv,")
    print(f"         specificity_distribution.csv, specificity_by_gene_set.png,")
    print(f"         niche_hypothesis_results.json, overlap_gene_details.csv")

    return summary


if __name__ == "__main__":
    main()
