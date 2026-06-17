#!/usr/bin/env python3
"""
CellWarp — GO Enrichment Stability Check (T2-C Validation)

Tests whether GO enrichment results are stable across different gene list
size thresholds (top-20, top-50, top-100). For each of the 6 original cell
types, we compare enrichment results at the three cutoffs to check that
the same biological themes emerge regardless of how many genes are included.

Biology
-------
If the Procrustes residual signal is real, the top biological themes should
be consistent: the most significant GO terms at top-20 should also appear
at top-50 and top-100, possibly with more granular terms emerging at larger
gene sets. Instability would suggest the enrichment signal is driven by a
few genes at the boundary of the cutoff.

Pass criterion: same top biological theme per cell type at all three
thresholds (e.g., hepatocyte always shows liver metabolism terms).

Inputs
------
- output/phase2/scaled_35types/centroids_{human,mouse}_35.csv
- output/phase2/scaled_35types/procrustes_results_35.json (residual vectors)
- data/phase1/human_aligned.h5ad (for gene symbol mapping)

Outputs (→ output/phase3/go_enrichment/stability/)
-------
- enrichment_{cell_type}_top{N}.csv  — Full results per threshold
- stability_report.txt               — Human-readable comparison
- stability_summary.json             — Machine-readable summary
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from cellwarp.enrichment import run_enrichment, P_ADJ_CUTOFF
from cellwarp.procrustes import RANDOM_SEED, PCA_VARIANCE_THRESHOLD

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ORIGINAL_6 = [
    "B cell",
    "CD4-positive, alpha-beta T cell",
    "CD8-positive, alpha-beta T cell",
    "endothelial cell",
    "hepatocyte",
    "macrophage",
]
THRESHOLDS = [20, 50, 100]
OUTPUT_DIR = PROJECT_ROOT / "output" / "phase3" / "go_enrichment" / "stability"


def safe_name(ct: str) -> str:
    """Convert cell type name to filesystem-safe string."""
    return ct.replace(" ", "_").replace(",", "").replace("+", "plus")


# ---------------------------------------------------------------------------
# Stage 1: Reconstruct full gene-space loadings
# ---------------------------------------------------------------------------


def reconstruct_gene_loadings() -> tuple[dict[str, np.ndarray], list[str]]:
    """
    Refit PCA on 35-type centroids and project residual vectors back to gene
    space to get full 16,959-dimensional loading vectors for each cell type.

    Returns:
        Tuple of (gene_loadings_dict, gene_symbols) where:
        - gene_loadings_dict: maps cell type → (16959,) array of gene loadings
        - gene_symbols: ordered list of gene symbols
    """
    print("[1/4] Reconstructing gene-space loadings from PCA + residuals...")

    # Load centroids
    h_centroids = pd.read_csv(
        PROJECT_ROOT / "output/phase2/scaled_35types/centroids_human_35.csv",
        index_col=0,
    )
    m_centroids = pd.read_csv(
        PROJECT_ROOT / "output/phase2/scaled_35types/centroids_mouse_35.csv",
        index_col=0,
    )
    gene_ids = list(h_centroids.columns)
    print(f"  Centroids: {h_centroids.shape[0]} cell types × {len(gene_ids)} genes")

    # Get gene symbols from h5ad
    import anndata as ad

    human_h5 = ad.read_h5ad(
        PROJECT_ROOT / "data/phase1/human_aligned.h5ad", backed="r"
    )
    id_to_symbol = dict(zip(human_h5.var_names, human_h5.var["feature_name"]))
    human_h5.file.close()
    gene_symbols = [id_to_symbol.get(gid, gid) for gid in gene_ids]
    n_mapped = sum(1 for s in gene_symbols if not s.startswith("ENSG"))
    print(f"  Gene symbols mapped: {n_mapped}/{len(gene_symbols)}")

    # Refit PCA (same parameters as original pipeline)
    ct_order = sorted(h_centroids.index.tolist())
    h_mat = h_centroids.loc[ct_order].values
    m_mat = m_centroids.loc[ct_order].values
    combined = np.vstack([h_mat, m_mat])

    pca = PCA(
        n_components=PCA_VARIANCE_THRESHOLD,
        svd_solver="full",
        random_state=RANDOM_SEED,
    )
    combined_pca = pca.fit_transform(combined)
    n_types = len(ct_order)
    human_pca = combined_pca[:n_types]
    mouse_pca = combined_pca[n_types:]
    print(
        f"  PCA: {pca.n_components_} components "
        f"({np.sum(pca.explained_variance_ratio_) * 100:.1f}% variance)"
    )

    # Load residual vectors from JSON
    with open(
        PROJECT_ROOT / "output/phase2/scaled_35types/procrustes_results_35.json"
    ) as f:
        proc_data = json.load(f)

    # Project residuals to gene space for each of the 6 cell types
    W = pca.components_  # (k, G)
    gene_loadings = {}
    for ct in ORIGINAL_6:
        r_pca = np.array(proc_data["residuals"][ct]["vector_pca"])
        g = r_pca @ W  # (G,)
        gene_loadings[ct] = g
        print(f"  {ct}: residual mag={np.linalg.norm(r_pca):.4f}, "
              f"max gene loading={np.max(np.abs(g)):.4f}")

    return gene_loadings, gene_symbols


# ---------------------------------------------------------------------------
# Stage 2: Extract top-N genes and run enrichment
# ---------------------------------------------------------------------------


def run_enrichment_at_thresholds(
    gene_loadings: dict[str, np.ndarray],
    gene_symbols: list[str],
) -> dict[str, dict[int, pd.DataFrame]]:
    """
    For each cell type and threshold, extract top-N genes and run GO enrichment.

    Returns:
        Nested dict: results[cell_type][threshold] = enrichment DataFrame
    """
    print("\n[2/4] Running GO enrichment at 3 thresholds...")

    results = {}
    for ct in ORIGINAL_6:
        results[ct] = {}
        loadings = gene_loadings[ct]

        # Sort genes by absolute loading
        order = np.argsort(np.abs(loadings))[::-1]

        for n_top in THRESHOLDS:
            top_idx = order[:n_top]
            top_genes = [gene_symbols[i] for i in top_idx]
            print(f"  {ct} top-{n_top}: running enrichment ({len(top_genes)} genes)...")

            df = run_enrichment(top_genes, gene_sets="GO_Biological_Process_2023")
            n_sig = (df["Adjusted P-value"] < P_ADJ_CUTOFF).sum() if not df.empty else 0
            print(f"    → {n_sig} significant terms")
            results[ct][n_top] = df

    return results


# ---------------------------------------------------------------------------
# Stage 3: Compare and report
# ---------------------------------------------------------------------------


def shorten_term(term: str) -> str:
    """Remove GO ID from term name."""
    return term.split(" (GO:")[0] if " (GO:" in term else term


def compare_stability(
    results: dict[str, dict[int, pd.DataFrame]],
) -> dict[str, dict]:
    """
    For each cell type, compare GO enrichment across thresholds.

    Returns:
        Summary dict with stability metrics per cell type.
    """
    print("\n[3/4] Comparing enrichment stability across thresholds...")

    summary = {}
    report_lines = []
    report_lines.append("=" * 75)
    report_lines.append("  GO ENRICHMENT STABILITY CHECK — T2-C Validation")
    report_lines.append("=" * 75)
    report_lines.append(
        f"  Thresholds tested: {THRESHOLDS}"
    )
    report_lines.append(
        f"  Significance: adjusted p-value < {P_ADJ_CUTOFF}"
    )
    report_lines.append("")

    all_stable = True

    for ct in ORIGINAL_6:
        report_lines.append(f"\n{'─' * 75}")
        report_lines.append(f"  {ct}")
        report_lines.append(f"{'─' * 75}")

        ct_summary = {"thresholds": {}}

        # Get significant terms at each threshold
        sig_terms = {}
        top5_terms = {}
        for n_top in THRESHOLDS:
            df = results[ct][n_top]
            if df.empty:
                sig_terms[n_top] = set()
                top5_terms[n_top] = []
            else:
                sig_df = df[df["Adjusted P-value"] < P_ADJ_CUTOFF]
                sig_terms[n_top] = set(sig_df["Term"].tolist())
                top5_terms[n_top] = sig_df.head(5)["Term"].tolist()

            n_sig = len(sig_terms[n_top])
            ct_summary["thresholds"][str(n_top)] = {
                "n_significant": n_sig,
                "top5_terms": [shorten_term(t) for t in top5_terms[n_top]],
            }

            report_lines.append(f"\n  Top-{n_top}: {n_sig} significant terms")
            if top5_terms[n_top]:
                for i, t in enumerate(top5_terms[n_top], 1):
                    p = df[df["Term"] == t]["Adjusted P-value"].values[0]
                    report_lines.append(
                        f"    {i}. {shorten_term(t)} (p_adj={p:.2e})"
                    )
            else:
                report_lines.append("    (no significant terms)")

        # Overlap analysis: what fraction of top-20 sig terms also in top-50/100?
        report_lines.append(f"\n  Overlap analysis:")

        base_terms = sig_terms[20]
        if base_terms:
            overlap_50 = base_terms & sig_terms[50]
            overlap_100 = base_terms & sig_terms[100]
            frac_50 = len(overlap_50) / len(base_terms) if base_terms else 0
            frac_100 = len(overlap_100) / len(base_terms) if base_terms else 0

            ct_summary["overlap_top20_in_top50"] = frac_50
            ct_summary["overlap_top20_in_top100"] = frac_100
            ct_summary["n_top20_sig"] = len(base_terms)

            report_lines.append(
                f"    Top-20 sig terms retained in top-50: "
                f"{len(overlap_50)}/{len(base_terms)} ({frac_50:.0%})"
            )
            report_lines.append(
                f"    Top-20 sig terms retained in top-100: "
                f"{len(overlap_100)}/{len(base_terms)} ({frac_100:.0%})"
            )

            # Check for terms lost at larger gene sets
            lost_50 = base_terms - sig_terms[50]
            lost_100 = base_terms - sig_terms[100]
            if lost_50:
                report_lines.append(
                    f"    Lost at top-50: {[shorten_term(t) for t in list(lost_50)[:3]]}"
                )
            if lost_100:
                report_lines.append(
                    f"    Lost at top-100: {[shorten_term(t) for t in list(lost_100)[:3]]}"
                )
        else:
            ct_summary["overlap_top20_in_top50"] = None
            ct_summary["overlap_top20_in_top100"] = None
            ct_summary["n_top20_sig"] = 0
            report_lines.append("    (no significant terms at top-20 to compare)")

        # New themes at larger gene sets
        new_50 = sig_terms[50] - sig_terms[20]
        new_100 = sig_terms[100] - sig_terms[20]
        ct_summary["n_new_at_top50"] = len(new_50)
        ct_summary["n_new_at_top100"] = len(new_100)

        if new_50:
            report_lines.append(
                f"\n  New themes at top-50 ({len(new_50)} new terms):"
            )
            # Show top 3 new terms by p-value
            df50 = results[ct][50]
            new_df = df50[df50["Term"].isin(new_50)].head(3)
            for _, row in new_df.iterrows():
                report_lines.append(
                    f"    + {shorten_term(row['Term'])} (p_adj={row['Adjusted P-value']:.2e})"
                )

        if new_100 - new_50:
            exclusively_new_100 = new_100 - new_50
            report_lines.append(
                f"\n  New themes only at top-100 ({len(exclusively_new_100)} additional):"
            )
            df100 = results[ct][100]
            new_df = df100[df100["Term"].isin(exclusively_new_100)].head(3)
            for _, row in new_df.iterrows():
                report_lines.append(
                    f"    + {shorten_term(row['Term'])} (p_adj={row['Adjusted P-value']:.2e})"
                )

        # Biological theme stability: is the #1 term at top-20 still in top-5
        # at top-50 and top-100?
        report_lines.append(f"\n  Theme stability:")
        if top5_terms[20]:
            top1_term_20 = top5_terms[20][0]
            top1_short = shorten_term(top1_term_20)

            in_top5_50 = top1_term_20 in top5_terms[50]
            in_top5_100 = top1_term_20 in top5_terms[100]

            # Also check if it's significant at larger thresholds (even if not top-5)
            in_sig_50 = top1_term_20 in sig_terms[50]
            in_sig_100 = top1_term_20 in sig_terms[100]

            ct_summary["top1_term"] = top1_short
            ct_summary["top1_in_top5_at_50"] = in_top5_50
            ct_summary["top1_in_top5_at_100"] = in_top5_100
            ct_summary["top1_sig_at_50"] = in_sig_50
            ct_summary["top1_sig_at_100"] = in_sig_100

            stable = in_sig_50 and in_sig_100
            ct_summary["theme_stable"] = stable
            if not stable:
                all_stable = False

            report_lines.append(f"    Top-1 theme at top-20: {top1_short}")
            report_lines.append(
                f"    In top-5 at top-50: {'YES' if in_top5_50 else 'NO'} | "
                f"Significant: {'YES' if in_sig_50 else 'NO'}"
            )
            report_lines.append(
                f"    In top-5 at top-100: {'YES' if in_top5_100 else 'NO'} | "
                f"Significant: {'YES' if in_sig_100 else 'NO'}"
            )
            report_lines.append(
                f"    VERDICT: {'STABLE' if stable else 'UNSTABLE'}"
            )
        else:
            ct_summary["theme_stable"] = None
            report_lines.append("    (no significant terms at top-20)")

        summary[ct] = ct_summary

    # Overall verdict
    n_stable = sum(
        1 for v in summary.values()
        if v.get("theme_stable") is True
    )
    n_tested = sum(
        1 for v in summary.values()
        if v.get("theme_stable") is not None
    )

    report_lines.append(f"\n{'=' * 75}")
    report_lines.append(f"  OVERALL STABILITY VERDICT")
    report_lines.append(f"{'=' * 75}")
    report_lines.append(
        f"  Cell types with stable top theme: {n_stable}/{n_tested}"
    )
    report_lines.append(
        f"  Pass criterion: same top biological theme at all three thresholds"
    )
    report_lines.append(
        f"  Result: {'PASS' if all_stable else 'PARTIAL — see details'}"
    )

    # Print overlap summary table
    report_lines.append(f"\n  {'Cell Type':<45} {'20→50':>7} {'20→100':>7} {'Stable':>7}")
    report_lines.append(f"  {'─' * 70}")
    for ct in ORIGINAL_6:
        s = summary[ct]
        o50 = f"{s['overlap_top20_in_top50']:.0%}" if s.get("overlap_top20_in_top50") is not None else "N/A"
        o100 = f"{s['overlap_top20_in_top100']:.0%}" if s.get("overlap_top20_in_top100") is not None else "N/A"
        stable = "YES" if s.get("theme_stable") is True else ("NO" if s.get("theme_stable") is False else "N/A")
        report_lines.append(f"  {ct:<45} {o50:>7} {o100:>7} {stable:>7}")

    report_lines.append(f"{'=' * 75}")

    report = "\n".join(report_lines)
    print(report)

    return summary, report


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    """Run GO enrichment stability check."""
    print("=" * 75)
    print("  CellWarp — GO Enrichment Stability Check (T2-C)")
    print("=" * 75)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Stage 1: Reconstruct gene-space loadings
    gene_loadings, gene_symbols = reconstruct_gene_loadings()

    # Stage 2: Run enrichment at 3 thresholds
    results = run_enrichment_at_thresholds(gene_loadings, gene_symbols)

    # Save per-threshold CSVs
    print("\n[3/4] Saving per-threshold enrichment CSVs...")
    for ct in ORIGINAL_6:
        for n_top in THRESHOLDS:
            df = results[ct][n_top]
            fname = f"enrichment_{safe_name(ct)}_top{n_top}.csv"
            df.to_csv(OUTPUT_DIR / fname, index=False)

    # Stage 3: Compare stability
    summary, report = compare_stability(results)

    # Save report and summary
    print("\n[4/4] Saving outputs...")
    with open(OUTPUT_DIR / "stability_report.txt", "w") as f:
        f.write(report)
    print(f"  Report: {OUTPUT_DIR / 'stability_report.txt'}")

    # JSON summary
    with open(OUTPUT_DIR / "stability_summary.json", "w") as f:
        json.dump(
            {
                "thresholds": THRESHOLDS,
                "p_adj_cutoff": P_ADJ_CUTOFF,
                "per_cell_type": summary,
            },
            f,
            indent=2,
        )
    print(f"  Summary: {OUTPUT_DIR / 'stability_summary.json'}")
    print("\nDone.")


if __name__ == "__main__":
    main()
