#!/usr/bin/env python3
"""
V3 — CellMarker Hypergeometric Test + Background Matching: Independent Validation
==================================================================================
Verifies the CellMarker enrichment statistics for the global 500-gene identity set.

NO imports from src/ — everything reimplemented from scratch.

Canonical targets:
  Foreground: fold = 4.49×, p = 2.10×10⁻¹³, overlaps = 34
  Background: fold = 3.22×, p = 1.16×10⁻¹²

Pass criteria:
  Foreground fold within ±0.05, p within 1 OOM, overlaps exactly 34
  Background fold within ±0.05, p within 1 OOM

Additional confirmations:
  - Global test uses 500-gene centroid-variance set
  - Per-type test (6/6) uses 50-gene centroid-deviation set (out of scope for Tier 1)
  - These are distinct gene sets — confirm no conflation

Author: V3 software validation (independent of src/)
Date: 2026-03-21
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# ── File paths ─────────────────────────────────────────────────────────
BASE = Path(__file__).resolve().parent.parent

CENTROID_35 = BASE / "output/phase2/scaled_35types/centroids_human_35.csv"
ORTHOLOGS = BASE / "data/phase1/orthologs_human_mouse.csv"
CELLMARKER_HUMAN = BASE / "data/validation/cellmarker/cellmarker_human_filtered.csv"
OUTPUT_DIR = BASE / "output/validation/v3_cellmarker"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ======================================================================
# Independent hypergeometric test
# ======================================================================

def hypergeom_enrichment(
    test_genes: set, marker_genes: set, background_size: int
) -> dict:
    """
    One-sided hypergeometric test for enrichment.

    H0: The overlap between test_genes and marker_genes is no greater
        than expected by chance.
    H1: The overlap is greater than expected.

    Parameters:
        N = background_size (total gene universe)
        K = |marker_genes| (CellMarker genes in universe)
        n = |test_genes| (identity genes)
        k = |test_genes ∩ marker_genes| (observed overlap)

    P(X ≥ k) = 1 - P(X ≤ k-1) = sf(k-1, N, K, n)
    fold = k / expected = k / (n × K / N)
    """
    K = len(marker_genes)
    n = len(test_genes)
    overlap = test_genes & marker_genes
    k = len(overlap)
    expected = n * K / background_size
    fold = k / expected if expected > 0 else 0.0
    p_val = float(stats.hypergeom.sf(k - 1, background_size, K, n))

    return {
        "k": k,
        "K": K,
        "n": n,
        "N": background_size,
        "expected": expected,
        "fold": fold,
        "p_value": p_val,
        "overlap_genes": sorted(overlap),
    }


# ======================================================================
# Main
# ======================================================================

def main():
    print("=" * 70)
    print("  V3 — CellMarker Hypergeometric Test: Independent Validation")
    print("=" * 70)

    # ── Load centroids ────────────────────────────────────────────────
    print("\nStep 0: Load data")
    centroids = pd.read_csv(CENTROID_35, index_col=0)
    assert centroids.shape[0] == 35, f"Expected 35 types, got {centroids.shape[0]}"
    print(f"  Centroids: {centroids.shape[0]} types × {centroids.shape[1]} genes")

    gene_ids = list(centroids.columns)
    n_genes = len(gene_ids)
    print(f"  Gene universe N = {n_genes}")

    # ── Load ortholog mapping ─────────────────────────────────────────
    orthologs = pd.read_csv(ORTHOLOGS)
    ens_to_symbol = dict(zip(orthologs["human_ensembl_id"], orthologs["human_gene_name"]))
    gene_symbols = [ens_to_symbol.get(g, g) for g in gene_ids]
    ens_to_sym_map = dict(zip(gene_ids, gene_symbols))

    # ── Load CellMarker ───────────────────────────────────────────────
    cellmarker = pd.read_csv(CELLMARKER_HUMAN)
    print(f"\n  CellMarker file: {CELLMARKER_HUMAN}")
    print(f"  CellMarker rows: {len(cellmarker)}")
    print(f"  CellMarker columns: {list(cellmarker.columns)}")

    # Record database metadata
    cm_unique_genes = set(cellmarker["gene_symbol"].dropna().unique())
    print(f"  CellMarker unique gene symbols (all): {len(cm_unique_genes)}")

    # Filter to our gene space
    cellmarker_in_bg = cm_unique_genes & set(gene_symbols)
    K = len(cellmarker_in_bg)
    print(f"  CellMarker genes in 16,959-gene background (K): {K}")

    # Back-calculate expected K from canonical fold:
    # fold = k / (n * K / N) → K = k * N / (n * fold)
    k_expected = 34
    fold_expected = 4.49
    K_backcalc = k_expected * n_genes / (500 * fold_expected)
    print(f"  K back-calculated from canonical fold: {K_backcalc:.1f}")
    print(f"  K actual: {K}")
    if abs(K - round(K_backcalc)) > 2:
        print(f"  *** FLAG: K mismatch — back-calculated {K_backcalc:.0f} vs actual {K} ***")

    # ══════════════════════════════════════════════════════════════════
    # Step 1: Compute global 500-gene identity set
    # ══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("  Step 1: Global 500-gene identity set (centroid variance)")
    print(f"{'=' * 70}")

    centroid_matrix = centroids.values  # (35, n_genes)
    gene_variances = np.var(centroid_matrix, axis=0)
    top500_idx = np.argsort(gene_variances)[::-1][:500]
    top500_symbols = set(ens_to_sym_map[gene_ids[i]] for i in top500_idx)

    print(f"  Top 500 genes by centroid variance across 35 types")
    print(f"  Variance range: {gene_variances[top500_idx[0]]:.6f} (max) "
          f"to {gene_variances[top500_idx[499]]:.6f} (500th)")
    print(f"  Unique symbols: {len(top500_symbols)}")

    # ══════════════════════════════════════════════════════════════════
    # Step 2: Confirm distinct gene sets (global vs per-type)
    # ══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("  Step 2: Confirm global vs per-type gene sets are distinct")
    print(f"{'=' * 70}")

    # Per-type: top 50 genes by absolute deviation from global mean centroid
    global_mean = np.mean(centroid_matrix, axis=0)
    validated_types = [
        "B cell", "CD4-positive, alpha-beta T cell",
        "CD8-positive, alpha-beta T cell", "endothelial cell",
        "hepatocyte", "macrophage",
    ]
    per_type_top50 = {}
    for ct in validated_types:
        ct_centroid = centroids.loc[ct].values
        deviation = np.abs(ct_centroid - global_mean)
        top50_idx_ct = np.argsort(deviation)[::-1][:50]
        per_type_top50[ct] = set(ens_to_sym_map[gene_ids[i]] for i in top50_idx_ct)

    # Check: are these distinct from the global 500?
    all_per_type_genes = set()
    for genes in per_type_top50.values():
        all_per_type_genes |= genes

    overlap_global_pertype = top500_symbols & all_per_type_genes
    print(f"  Global 500-gene set: {len(top500_symbols)} genes (centroid variance)")
    print(f"  Per-type union (6 × 50): {len(all_per_type_genes)} unique genes (centroid deviation)")
    print(f"  Overlap: {len(overlap_global_pertype)} genes")
    print(f"  CONFIRMATION: These are distinct gene sets derived by different methods:")
    print(f"    - Global: top 500 by VARIANCE across all 35 type centroids")
    print(f"    - Per-type: top 50 by absolute DEVIATION of each type from global mean")
    print(f"  Per-type 6/6 test is OUT OF SCOPE for Tier 1 V3 validation.")

    # ══════════════════════════════════════════════════════════════════
    # Step 3a: FOREGROUND — Global enrichment test
    # ══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("  Step 3a: FOREGROUND — Global 500-gene enrichment")
    print(f"{'=' * 70}")

    fg = hypergeom_enrichment(top500_symbols, cellmarker_in_bg, n_genes)

    print(f"  k (observed overlap):    {fg['k']}")
    print(f"  K (CellMarker in bg):    {fg['K']}")
    print(f"  n (identity genes):      {fg['n']}")
    print(f"  N (gene universe):       {fg['N']}")
    print(f"  Expected overlap:        {fg['expected']:.2f}")
    print(f"  Fold enrichment:         {fg['fold']:.3f}")
    print(f"  p-value:                 {fg['p_value']:.4e}")
    print(f"  Overlapping genes ({fg['k']}):")
    for i, g in enumerate(fg["overlap_genes"]):
        print(f"    {i+1:>2}. {g}")

    # Pass/fail foreground
    fg_overlap_pass = fg["k"] == 34
    fg_fold_pass = abs(fg["fold"] - 4.49) <= 0.05
    fg_p_oom = abs(np.log10(fg["p_value"]) - np.log10(2.10e-13))
    fg_p_pass = fg_p_oom <= 1.0

    print(f"\n  FOREGROUND ASSESSMENT:")
    print(f"    Overlap = {fg['k']}: {'PASS' if fg_overlap_pass else '*** FAIL ***'} (expected 34)")
    print(f"    Fold = {fg['fold']:.3f}: {'PASS' if fg_fold_pass else '*** FAIL ***'} "
          f"(expected 4.49 ±0.05, diff={abs(fg['fold'] - 4.49):.3f})")
    print(f"    p = {fg['p_value']:.4e}: {'PASS' if fg_p_pass else '*** FAIL ***'} "
          f"(expected 2.10e-13, OOM diff={fg_p_oom:.2f})")

    # ══════════════════════════════════════════════════════════════════
    # Step 3b: BACKGROUND — Expression-matched control
    # ══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("  Step 3b: BACKGROUND — Expression-matched control")
    print(f"{'=' * 70}")

    # Mean expression per gene across 35 type centroids
    mean_expr = np.mean(centroid_matrix, axis=0)
    identity_idx_set = set(top500_idx)
    matched_bg_idx = set()

    for idx in top500_idx:
        target_expr = mean_expr[idx]
        lo = target_expr * 0.9
        hi = target_expr * 1.1
        candidates = [
            (j, abs(mean_expr[j] - target_expr))
            for j in range(n_genes)
            if j not in identity_idx_set and lo <= mean_expr[j] <= hi
        ]
        candidates.sort(key=lambda x: (x[1], gene_ids[x[0]]))
        matched_bg_idx.update(c[0] for c in candidates[:10])

    matched_bg_symbols = set(ens_to_sym_map[gene_ids[i]] for i in matched_bg_idx)
    universe = top500_symbols | matched_bg_symbols
    universe_size = len(universe)

    print(f"  Expression-matched background genes: {len(matched_bg_symbols)}")
    print(f"  Restricted universe: {universe_size} genes")

    cellmarker_in_universe = cellmarker_in_bg & universe
    print(f"  CellMarker genes in restricted universe: {len(cellmarker_in_universe)}")

    bg = hypergeom_enrichment(
        top500_symbols & universe, cellmarker_in_universe, universe_size
    )

    print(f"\n  k (observed overlap):    {bg['k']}")
    print(f"  K (CellMarker in univ):  {bg['K']}")
    print(f"  n (identity in univ):    {bg['n']}")
    print(f"  N (restricted universe): {bg['N']}")
    print(f"  Expected overlap:        {bg['expected']:.2f}")
    print(f"  Fold enrichment:         {bg['fold']:.3f}")
    print(f"  p-value:                 {bg['p_value']:.4e}")

    # Pass/fail background
    bg_fold_pass = abs(bg["fold"] - 3.22) <= 0.05
    bg_p_oom = abs(np.log10(bg["p_value"]) - np.log10(1.16e-12))
    bg_p_pass = bg_p_oom <= 1.0

    print(f"\n  BACKGROUND ASSESSMENT:")
    print(f"    Fold = {bg['fold']:.3f}: {'PASS' if bg_fold_pass else '*** FAIL ***'} "
          f"(expected 3.22 ±0.05, diff={abs(bg['fold'] - 3.22):.3f})")
    print(f"    p = {bg['p_value']:.4e}: {'PASS' if bg_p_pass else '*** FAIL ***'} "
          f"(expected 1.16e-12, OOM diff={bg_p_oom:.2f})")

    # ══════════════════════════════════════════════════════════════════
    # Summary
    # ══════════════════════════════════════════════════════════════════
    all_pass = fg_overlap_pass and fg_fold_pass and fg_p_pass and bg_fold_pass and bg_p_pass
    verdict = "PASS" if all_pass else "FAIL"

    print(f"\n{'=' * 70}")
    print(f"  V3 VERDICT: {verdict}")
    print(f"{'=' * 70}")

    if not all_pass:
        failures = []
        if not fg_overlap_pass:
            failures.append(f"Foreground overlap: {fg['k']} (expected 34)")
        if not fg_fold_pass:
            failures.append(f"Foreground fold: {fg['fold']:.3f} (expected 4.49)")
        if not fg_p_pass:
            failures.append(f"Foreground p-value OOM diff: {fg_p_oom:.2f}")
        if not bg_fold_pass:
            failures.append(f"Background fold: {bg['fold']:.3f} (expected 3.22)")
        if not bg_p_pass:
            failures.append(f"Background p-value OOM diff: {bg_p_oom:.2f}")
        for f in failures:
            print(f"  ! {f}")

    print(f"\n  CellMarker database: CellMarker 2.0")
    print(f"  Filter: wet-lab validated (Experiment source)")
    print(f"  Download date: 2026-03-16")
    print(f"  File: {CELLMARKER_HUMAN}")

    # ── Save results ──────────────────────────────────────────────────
    output = {
        "validation": "V3 — CellMarker Hypergeometric Test",
        "date": "2026-03-21",
        "verdict": verdict,
        "foreground": {
            "k": fg["k"],
            "K": fg["K"],
            "n": fg["n"],
            "N": fg["N"],
            "expected": fg["expected"],
            "fold": fg["fold"],
            "p_value": fg["p_value"],
            "overlap_genes": fg["overlap_genes"],
            "pass_overlap": fg_overlap_pass,
            "pass_fold": fg_fold_pass,
            "pass_p": fg_p_pass,
        },
        "background": {
            "k": bg["k"],
            "K": bg["K"],
            "n": bg["n"],
            "N": bg["N"],
            "expected": bg["expected"],
            "fold": bg["fold"],
            "p_value": bg["p_value"],
            "universe_size": universe_size,
            "matched_bg_genes": len(matched_bg_symbols),
            "pass_fold": bg_fold_pass,
            "pass_p": bg_p_pass,
        },
        "gene_set_confirmation": {
            "global_method": "top 500 by centroid variance across 35 types",
            "per_type_method": "top 50 by absolute deviation from global mean centroid",
            "are_distinct_methods": True,
            "overlap_global_pertype": len(overlap_global_pertype),
            "per_type_validation_in_scope": False,
        },
        "cellmarker_metadata": {
            "version": "CellMarker 2.0",
            "filter": "Experiment (wet-lab validated only)",
            "download_date": "2026-03-16",
            "download_commit": "0f2c22d",
            "source_url": "http://xteam.xbio.top/CellMarker/",
            "file": str(CELLMARKER_HUMAN.relative_to(BASE)),
        },
    }

    # Add Methods text discrepancy flag
    output["methods_text_flags"] = [
        "Methods line 93 states '1,794 unique expression-matched background genes' and "
        "'2,294 genes total' — these are 6-TYPE numbers (DECISION-114). "
        "The fold (3.22) and p (1.16e-12) are from the 35-type stored result "
        "(DECISION-141), which claims 1,964 bg genes / 2,464 universe. "
        "Methods is mixing numbers from two different runs.",
        "The stored 35-type background result (1,964 bg genes) CANNOT BE "
        "REPRODUCED from the current cellmarker_35type_rerun.py code, which "
        "produces 208 bg genes. The script was committed 2 days after results "
        "were generated (b3e08b0, 2026-03-20 vs DECISION-141 dated 2026-03-18). "
        "The intermediate code version is lost.",
    ]

    out_path = OUTPUT_DIR / "v3_cellmarker_results.json"

    def json_default(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, (np.bool_, bool)):
            return bool(obj)
        return str(obj)

    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, default=json_default)
    print(f"\n  Results saved to: {out_path}")


if __name__ == "__main__":
    main()
