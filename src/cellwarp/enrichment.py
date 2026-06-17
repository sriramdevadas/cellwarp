"""
CellWarp — GO Enrichment Module

Performs Gene Ontology (GO) over-representation analysis on Procrustes residual
gene sets using gseapy's Enrichr interface. Tests whether genes driving
cell-type-specific cross-species divergence are enriched for known biological
processes.

Biology
-------
After Procrustes alignment, each cell type has a residual vector in gene space
representing what the global geometric transformation could NOT explain. The
top genes in this residual capture cell-type-specific evolutionary divergence.
If these genes are enriched for coherent biological processes (e.g., liver
metabolism genes for hepatocytes, complement genes for macrophages), it
validates that the Procrustes residuals reflect real biology — not noise.

Math
----
GO enrichment uses the hypergeometric test (Fisher's exact test):

    P(X ≥ k) = Σ_{i=k}^{min(K,n)} C(K,i) C(N-K,n-i) / C(N,n)

where N = background gene count, K = genes in the GO term, n = query gene
list size, k = overlap. The p-values are corrected for multiple testing
using Benjamini-Hochberg FDR.

Phase 3 gate criterion: ≥3 of 6 cell types must have GO enrichment with
adjusted p-value < 0.05.
"""

from __future__ import annotations

import json
import random
import time
from pathlib import Path
from typing import Any

import gseapy
from gseapy.enrichr import EnrichrAPIError
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

GO_LIBRARY = "GO_Biological_Process_2023"
P_ADJ_CUTOFF = 0.05
N_TOP_TERMS = 10

# gseapy 1.1.12's built-in urllib3 retry whitelist excludes 429, so wrap
# the enrich() call ourselves. Numerical results are unaffected: Enrichr's
# hypergeometric test is deterministic for a fixed gene list + library.
_ENRICHR_MAX_RETRIES = 5
_ENRICHR_BASE_DELAY_S = 2.0
_ENRICHR_MAX_SLEEP_S = 60.0


def _is_retryable_enrichr_error(exc: Exception) -> bool:
    msg = str(exc)
    for code in ("429", "500", "502", "503", "504"):
        if code in msg:
            return True
    return False


def _enrichr_retry_sleep(attempt: int) -> float:
    base = _ENRICHR_BASE_DELAY_S * (2 ** attempt)
    jitter = random.uniform(0, base * 0.25)
    return min(base + jitter, _ENRICHR_MAX_SLEEP_S)


# ---------------------------------------------------------------------------
# Core enrichment
# ---------------------------------------------------------------------------


def run_enrichment(
    gene_list: list[str],
    gene_sets: str = GO_LIBRARY,
    background: list[str] | int | None = None,
    cutoff: float = P_ADJ_CUTOFF,
) -> pd.DataFrame:
    """
    Run GO over-representation analysis on a gene list.

    Uses gseapy.enrich() which performs a hypergeometric test against the
    specified gene set library (default: GO_Biological_Process_2023 from
    Enrichr).

    Args:
        gene_list: List of gene symbols to test for enrichment.
        gene_sets: Enrichr library name or custom gene sets dict/gmt.
        background: Background genes for the hypergeometric test. If None
            and using an Enrichr library, the library's own background is used.
        cutoff: Adjusted p-value cutoff for reporting (affects plots, not
            the returned DataFrame which includes all terms).

    Returns:
        DataFrame with enrichment results. Key columns: Term, Adjusted P-value,
        Overlap, Genes, Combined Score. Empty DataFrame if no enrichment found.
    """
    last_exc: Exception | None = None
    for attempt in range(_ENRICHR_MAX_RETRIES):
        try:
            enr = gseapy.enrich(
                gene_list=gene_list,
                gene_sets=gene_sets,
                background=background,
                outdir=None,
                cutoff=1.0,  # return all terms; we filter ourselves
                no_plot=True,
                verbose=False,
            )
            results = enr.results
            if results is None or results.empty:
                return pd.DataFrame()
            return results.sort_values("Adjusted P-value").reset_index(drop=True)
        except EnrichrAPIError as e:
            last_exc = e
            if _is_retryable_enrichr_error(e) and attempt < _ENRICHR_MAX_RETRIES - 1:
                sleep_s = _enrichr_retry_sleep(attempt)
                print(
                    f"    [retry] Enrichr error ({e}); sleeping {sleep_s:.1f}s "
                    f"before attempt {attempt + 2}/{_ENRICHR_MAX_RETRIES}"
                )
                time.sleep(sleep_s)
                continue
            raise
    assert last_exc is not None
    raise last_exc


def run_enrichment_all_cell_types(
    top_genes_per_ct: dict[str, list[str]],
    gene_sets: str = GO_LIBRARY,
    background: list[str] | int | None = None,
    checkpoint_dir: str | Path | None = None,
) -> dict[str, pd.DataFrame]:
    """
    Run GO enrichment for every cell type's residual gene list.

    Args:
        top_genes_per_ct: Dict mapping cell type name → list of gene symbols.
        gene_sets: Enrichr library name.
        background: Background genes (see run_enrichment).
        checkpoint_dir: If provided, each cell type's results are written to
            {checkpoint_dir}/_partial_{ct}.csv as the loop runs, so a mid-loop
            failure preserves completed work. Final atomic consolidation is
            the caller's responsibility.

    Returns:
        Dict mapping cell type → enrichment results DataFrame.
    """
    results = {}
    if checkpoint_dir is not None:
        checkpoint_dir = Path(checkpoint_dir)
        checkpoint_dir.mkdir(parents=True, exist_ok=True)
    for ct, genes in sorted(top_genes_per_ct.items()):
        print(f"\n  Running GO enrichment for {ct} ({len(genes)} genes)...")
        df = run_enrichment(genes, gene_sets=gene_sets, background=background)
        n_sig = (df["Adjusted P-value"] < P_ADJ_CUTOFF).sum() if not df.empty else 0
        print(f"    → {len(df)} terms tested, {n_sig} significant (p_adj < {P_ADJ_CUTOFF})")
        results[ct] = df
        if checkpoint_dir is not None:
            safe_ct = ct.replace("/", "_").replace(" ", "_")
            df.to_csv(checkpoint_dir / f"_partial_{safe_ct}.csv", index=False)
    return results


# ---------------------------------------------------------------------------
# Evaluation against Phase 3 gate
# ---------------------------------------------------------------------------


def evaluate_gate(
    enrichment_results: dict[str, pd.DataFrame],
    p_adj_cutoff: float = P_ADJ_CUTOFF,
    min_cell_types: int = 3,
) -> dict[str, Any]:
    """
    Evaluate Phase 3 GO enrichment gate criterion.

    Gate: ≥min_cell_types cell types must have at least one GO term with
    adjusted p-value < p_adj_cutoff.

    Args:
        enrichment_results: Dict mapping cell type → enrichment DataFrame.
        p_adj_cutoff: Significance threshold.
        min_cell_types: Minimum number of cell types required to pass.

    Returns:
        Dict with gate evaluation results.
    """
    ct_sig = {}
    for ct, df in sorted(enrichment_results.items()):
        if df.empty:
            ct_sig[ct] = {"n_significant": 0, "top_term": None, "top_p_adj": None}
        else:
            n_sig = int((df["Adjusted P-value"] < p_adj_cutoff).sum())
            top = df.iloc[0] if not df.empty else None
            ct_sig[ct] = {
                "n_significant": n_sig,
                "top_term": top["Term"] if top is not None else None,
                "top_p_adj": float(top["Adjusted P-value"]) if top is not None else None,
            }

    n_passing = sum(1 for v in ct_sig.values() if v["n_significant"] > 0)
    passed = n_passing >= min_cell_types

    return {
        "gate_passed": passed,
        "n_cell_types_with_enrichment": n_passing,
        "min_required": min_cell_types,
        "p_adj_cutoff": p_adj_cutoff,
        "per_cell_type": ct_sig,
    }


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def plot_dotplot(
    df: pd.DataFrame,
    cell_type: str,
    output_path: Path,
    n_top: int = N_TOP_TERMS,
) -> str:
    """
    Generate a dot plot for GO enrichment results of a single cell type.

    X-axis: -log10(adjusted p-value). Dot size: gene ratio (overlap / term size).
    Color: combined score.

    Args:
        df: Enrichment results DataFrame from run_enrichment.
        cell_type: Cell type name (for title).
        output_path: Path to save figure.
        n_top: Number of top terms to show.

    Returns:
        Text description of the plot.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if df.empty or (df["Adjusted P-value"] >= 1.0).all():
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, f"No enriched GO terms for {cell_type}",
                ha="center", va="center", fontsize=12, transform=ax.transAxes)
        ax.set_axis_off()
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return f"No enriched GO terms for {cell_type}."

    top = df.head(n_top).copy()
    # Parse overlap: "3/100" → gene_ratio = 0.03
    top["overlap_n"] = top["Overlap"].apply(lambda x: int(str(x).split("/")[0]))
    top["overlap_total"] = top["Overlap"].apply(lambda x: int(str(x).split("/")[1]))
    top["gene_ratio"] = top["overlap_n"] / top["overlap_total"]
    top["neg_log10_padj"] = -np.log10(top["Adjusted P-value"].clip(lower=1e-50))

    # Shorten term names for display
    top["Term_short"] = top["Term"].apply(
        lambda t: t.split(" (GO:")[0] if " (GO:" in t else t
    )
    top["Term_short"] = top["Term_short"].apply(
        lambda t: t[:60] + "..." if len(t) > 60 else t
    )

    fig, ax = plt.subplots(figsize=(10, max(4, 0.5 * len(top))))
    scatter = ax.scatter(
        top["neg_log10_padj"],
        range(len(top)),
        s=top["gene_ratio"] * 1000,
        c=top["Combined Score"],
        cmap="YlOrRd",
        edgecolors="black",
        linewidths=0.5,
        alpha=0.8,
    )
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top["Term_short"], fontsize=9)
    ax.set_xlabel("-log₁₀(Adjusted P-value)", fontsize=11)
    ax.set_title(f"GO Biological Process Enrichment: {cell_type}", fontsize=12)
    ax.invert_yaxis()
    ax.axvline(-np.log10(0.05), color="grey", linestyle="--", alpha=0.5, label="p_adj=0.05")

    # Color bar
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label("Combined Score", fontsize=9)

    # Size legend
    for ratio, label in [(0.05, "5%"), (0.1, "10%"), (0.2, "20%")]:
        ax.scatter([], [], s=ratio * 1000, c="grey", alpha=0.5, label=f"Gene ratio: {label}")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    # Text description
    n_sig = (top["Adjusted P-value"] < 0.05).sum()
    desc = (
        f"Dot plot for {cell_type}: {n_sig}/{len(top)} shown terms significant "
        f"(p_adj < 0.05). Top term: '{top.iloc[0]['Term_short']}' "
        f"(p_adj={top.iloc[0]['Adjusted P-value']:.2e})."
    )
    return desc


def plot_summary_heatmap(
    enrichment_results: dict[str, pd.DataFrame],
    output_path: Path,
    n_top: int = 5,
) -> str:
    """
    Generate a summary heatmap of top GO terms across all cell types.

    Each row is a GO term, each column is a cell type. Color intensity is
    -log10(adjusted p-value). Only includes terms that are in the top N
    for at least one cell type.

    Args:
        enrichment_results: Dict mapping cell type → enrichment DataFrame.
        output_path: Path to save figure.
        n_top: Number of top terms per cell type to include.

    Returns:
        Text description of the plot.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Collect top terms from each cell type
    all_terms = set()
    for ct, df in enrichment_results.items():
        if not df.empty:
            top_terms = df.head(n_top)["Term"].tolist()
            all_terms.update(top_terms)

    if not all_terms:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, "No enriched GO terms across any cell type",
                ha="center", va="center", fontsize=12, transform=ax.transAxes)
        ax.set_axis_off()
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return "No enriched GO terms across any cell type."

    all_terms = sorted(all_terms)
    cell_types = sorted(enrichment_results.keys())

    # Build matrix: -log10(p_adj) for each term × cell type
    matrix = pd.DataFrame(0.0, index=all_terms, columns=cell_types)
    for ct, df in enrichment_results.items():
        if df.empty:
            continue
        for _, row in df.iterrows():
            if row["Term"] in all_terms:
                matrix.loc[row["Term"], ct] = -np.log10(
                    max(row["Adjusted P-value"], 1e-50)
                )

    # Shorten term names
    short_names = [
        (t.split(" (GO:")[0][:55] + "..." if len(t.split(" (GO:")[0]) > 55
         else t.split(" (GO:")[0]) if " (GO:" in t else t[:55]
        for t in matrix.index
    ]
    # Shorten cell type names for columns
    short_ct = [ct.replace("-positive, alpha-beta ", "+ ") for ct in matrix.columns]

    fig, ax = plt.subplots(figsize=(max(8, len(cell_types) * 1.8), max(6, len(all_terms) * 0.35)))
    sns.heatmap(
        matrix.values,
        xticklabels=short_ct,
        yticklabels=short_names,
        cmap="YlOrRd",
        linewidths=0.5,
        linecolor="white",
        ax=ax,
        cbar_kws={"label": "-log₁₀(Adjusted P-value)", "shrink": 0.6},
    )
    ax.set_title("GO Enrichment Summary Across Cell Types", fontsize=13)
    ax.tick_params(axis="y", labelsize=8)
    ax.tick_params(axis="x", labelsize=9, rotation=45)
    plt.setp(ax.get_xticklabels(), ha="right", rotation_mode="anchor")

    # Significance line in colorbar
    sig_line = -np.log10(0.05)
    cbar = ax.collections[0].colorbar
    cbar.ax.axhline(sig_line, color="black", linestyle="--", linewidth=1)
    cbar.ax.text(1.5, sig_line, "p=0.05", va="center", fontsize=8, transform=cbar.ax.get_yaxis_transform())

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    n_nonzero = (matrix > 0).sum().sum()
    desc = (
        f"Heatmap: {len(all_terms)} GO terms × {len(cell_types)} cell types. "
        f"{n_nonzero} non-zero entries. "
        f"Color intensity = -log10(adjusted p-value). "
        f"Black dashed line on colorbar = p_adj=0.05 threshold."
    )
    return desc


# ---------------------------------------------------------------------------
# Interpretation
# ---------------------------------------------------------------------------


# Curated biological context for interpreting GO results per cell type.
# Genes are human symbols from the Procrustes residual analysis.
CELL_TYPE_BIOLOGY = {
    "B cell": {
        "expected_processes": [
            "immune response", "antigen presentation", "complement",
            "phagocytosis", "innate immune response",
        ],
        "context": (
            "B cell residual genes include complement (C1QB, C1QA, C1QC), "
            "MHC class II (CD74), and monocyte/macrophage markers (TYROBP, "
            "FCER1G, CD163). This may reflect cross-contamination from the "
            "macrophage-B cell axis in PCA space (ISSUE-011)."
        ),
    },
    "CD4-positive, alpha-beta T cell": {
        "expected_processes": [
            "B cell activation", "antigen presentation", "ribosome",
            "translation", "immune response",
        ],
        "context": (
            "CD4+ T cell residual is dominated by B cell markers (CD79A, "
            "MS4A1, BANK1) and ribosomal proteins (RPL39, RPS15A). The B "
            "cell markers likely reflect PCA projection artifacts (ISSUE-012)."
        ),
    },
    "CD8-positive, alpha-beta T cell": {
        "expected_processes": [
            "cytotoxic T cell activity", "immune response", "cell killing",
        ],
        "context": (
            "CD8+ T cell has the smallest residual magnitude (0.43, only "
            "0.1% of SSR). Gene loadings are very small (<0.06). GO "
            "enrichment may be weak due to low signal."
        ),
    },
    "endothelial cell": {
        "expected_processes": [
            "angiogenesis", "cell adhesion", "blood vessel development",
            "extracellular matrix organization",
        ],
        "context": (
            "Endothelial residual genes include VIM, ANXA2, SPARCL1, "
            "S100A6 — cytoskeletal/extracellular matrix genes consistent "
            "with vascular biology."
        ),
    },
    "hepatocyte": {
        "expected_processes": [
            "lipid metabolism", "cholesterol biosynthesis", "coagulation",
            "complement activation", "drug metabolism",
        ],
        "context": (
            "Hepatocyte residual genes include ALB (albumin), FGB/FGG "
            "(fibrinogen), CPS1 (urea cycle), APOB/APOC3/APOA1 "
            "(apolipoproteins), CYP2E1 (cytochrome P450). These are core "
            "liver functions — strong biological signal."
        ),
    },
    "macrophage": {
        "expected_processes": [
            "innate immune response", "complement activation", "phagocytosis",
            "inflammatory response", "antigen presentation",
        ],
        "context": (
            "Macrophage has the largest residual (60.7% of SSR). Top genes "
            "are complement (C1QB, C1QA), innate immune receptors (TYROBP, "
            "FCER1G), and lysosomal proteases (CTSB, CTSD). Strong "
            "biological signal for innate immunity divergence."
        ),
    },
}


def interpret_enrichment(
    ct: str,
    df: pd.DataFrame,
    n_top: int = N_TOP_TERMS,
) -> str:
    """
    Generate a human-readable interpretation of GO enrichment for one cell type.

    Args:
        ct: Cell type name.
        df: Enrichment results DataFrame.
        n_top: Number of top terms to report.

    Returns:
        Multi-line string with interpretation.
    """
    lines = [f"\n{'='*70}", f"  {ct}", f"{'='*70}"]

    # Add biological context
    bio = CELL_TYPE_BIOLOGY.get(ct, {})
    if bio.get("context"):
        lines.append(f"  Context: {bio['context']}")

    if df.empty:
        lines.append("  Result: NO enriched GO terms found.")
        lines.append("  Interpretation: Residual genes for this cell type do not")
        lines.append("  converge on identifiable biological processes.")
        return "\n".join(lines)

    n_sig = (df["Adjusted P-value"] < P_ADJ_CUTOFF).sum()
    lines.append(f"  Significant terms (p_adj < {P_ADJ_CUTOFF}): {n_sig}")

    if n_sig == 0:
        lines.append("  Result: No terms reach significance after FDR correction.")
        top_term = df.iloc[0]
        lines.append(
            f"  Best term: '{top_term['Term'].split(' (GO:')[0]}' "
            f"(p_adj={top_term['Adjusted P-value']:.3e})"
        )
        return "\n".join(lines)

    lines.append(f"\n  Top {min(n_top, n_sig)} significant terms:")
    for i, (_, row) in enumerate(df.head(min(n_top, len(df))).iterrows()):
        if row["Adjusted P-value"] >= P_ADJ_CUTOFF:
            break
        term_short = row["Term"].split(" (GO:")[0]
        genes = row.get("Genes", "")
        lines.append(
            f"    {i+1}. {term_short}"
        )
        lines.append(
            f"       p_adj={row['Adjusted P-value']:.2e}, "
            f"overlap={row['Overlap']}, "
            f"genes: {genes}"
        )

    # Check if expected processes are represented
    expected = bio.get("expected_processes", [])
    if expected:
        sig_terms = df[df["Adjusted P-value"] < P_ADJ_CUTOFF]["Term"].str.lower().tolist()
        matched = []
        for exp in expected:
            if any(exp.lower() in t for t in sig_terms):
                matched.append(exp)
        if matched:
            lines.append(f"\n  Expected biology confirmed: {', '.join(matched)}")
        else:
            lines.append(f"\n  None of expected processes ({', '.join(expected[:3])}...) found in top terms.")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Save results
# ---------------------------------------------------------------------------


def save_enrichment_results(
    enrichment_results: dict[str, pd.DataFrame],
    gate_eval: dict[str, Any],
    output_dir: Path,
) -> None:
    """
    Save enrichment results as JSON and per-cell-type CSVs.

    Args:
        enrichment_results: Dict mapping cell type → enrichment DataFrame.
        gate_eval: Gate evaluation dict from evaluate_gate.
        output_dir: Directory to save results.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save per-cell-type CSVs
    for ct, df in enrichment_results.items():
        safe_name = ct.replace(" ", "_").replace(",", "").replace("+", "plus")
        csv_path = output_dir / f"enrichment_{safe_name}.csv"
        df.to_csv(csv_path, index=False)

    # Save gate evaluation as JSON
    gate_path = output_dir / "go_enrichment_summary.json"
    with open(gate_path, "w") as f:
        json.dump(gate_eval, f, indent=2)
    print(f"\n  Gate evaluation saved: {gate_path}")

    # Save top terms summary as CSV
    summary_rows = []
    for ct, df in sorted(enrichment_results.items()):
        if df.empty:
            continue
        for _, row in df.head(N_TOP_TERMS).iterrows():
            summary_rows.append({
                "cell_type": ct,
                "term": row["Term"],
                "adjusted_p_value": row["Adjusted P-value"],
                "overlap": row["Overlap"],
                "genes": row.get("Genes", ""),
                "combined_score": row.get("Combined Score", None),
            })
    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        summary_path = output_dir / "top_terms_all_cell_types.csv"
        summary_df.to_csv(summary_path, index=False)
        print(f"  Top terms summary saved: {summary_path}")
