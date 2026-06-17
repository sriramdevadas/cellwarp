"""
CellMarker validation rerun using 35-type centroids.

Corrects the original issue073 validation which used 6-type centroids
(output/phase2/centroids_human.csv) to match the locked Methods description
claiming the identity gene set was defined across all 35 cell types.

Method:
  Step 1: Global 500-gene identity set — top 500 genes by centroid variance
          across all 35 cell type centroids in output/phase2/scaled_35types/centroids_human_35.csv.
  Step 2: Per-cell-type top 50 genes — absolute deviation of each type's
          centroid from the 35-type global mean centroid.
  Step 3: CellMarker 2.0 enrichment (hypergeometric) with expression-matched
          background control.

Output: output/validation/cellmarker_35type_rerun/
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# ── Paths ──────────────────────────────────────────────────────────────
CENTROID_35 = Path("output/phase2/scaled_35types/centroids_human_35.csv")
ORTHOLOGS = Path("data/phase1/orthologs_human_mouse.csv")
CELLMARKER_HUMAN = Path("data/validation/cellmarker/cellmarker_human_filtered.csv")
OUTPUT_DIR = Path("output/validation/cellmarker_35type_rerun")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 6 validated cell types (same as original)
VALIDATED_TYPES = [
    "B cell",
    "CD4-positive, alpha-beta T cell",
    "CD8-positive, alpha-beta T cell",
    "endothelial cell",
    "hepatocyte",
    "macrophage",
]

# ── Step 0: Load data ─────────────────────────────────────────────────
print("=" * 70)
print("CellMarker Validation Rerun — 35-Type Centroids")
print("=" * 70)

# Load 35-type centroids
centroids = pd.read_csv(CENTROID_35, index_col=0)
print(f"\nCentroid matrix: {centroids.shape[0]} cell types × {centroids.shape[1]} genes")
assert centroids.shape[0] == 35, f"Expected 35 types, got {centroids.shape[0]}"

# Load ortholog mapping (ensembl → gene symbol)
orthologs = pd.read_csv(ORTHOLOGS)
ens_to_symbol = dict(zip(orthologs["human_ensembl_id"], orthologs["human_gene_name"]))

# Gene IDs in centroid matrix
gene_ids = list(centroids.columns)
gene_symbols = [ens_to_symbol.get(g, g) for g in gene_ids]
ens_to_sym_map = dict(zip(gene_ids, gene_symbols))
n_genes = len(gene_ids)
print(f"Gene space: {n_genes} genes")

# Load CellMarker
cellmarker = pd.read_csv(CELLMARKER_HUMAN)
cellmarker_all_genes = set(cellmarker["gene_symbol"].dropna().unique())
# Intersect with our gene space
cellmarker_in_bg = cellmarker_all_genes & set(gene_symbols)
print(f"CellMarker genes in background: {len(cellmarker_in_bg)}")

# Per-cell-type CellMarker markers
cellmarker_per_type = {}
type_name_map = {
    "B cell": ["B cell"],
    "CD4-positive, alpha-beta T cell": ["CD4+ T cell", "CD4-positive, alpha-beta T cell",
                                         "Helper T cell", "CD4+ memory T cell"],
    "CD8-positive, alpha-beta T cell": ["CD8+ T cell", "CD8-positive, alpha-beta T cell",
                                         "Cytotoxic T cell", "CD8+ cytotoxic T cell"],
    "endothelial cell": ["Endothelial cell", "Vascular endothelial cell",
                         "Endothelial cell of lymphatic vessel"],
    "hepatocyte": ["Hepatocyte"],
    "macrophage": ["Macrophage", "Tissue-resident macrophage"],
}
for cw_type, cm_names in type_name_map.items():
    markers = set()
    for cm_name in cm_names:
        matches = cellmarker[cellmarker["cell_type"].str.lower() == cm_name.lower()]
        markers.update(matches["gene_symbol"].dropna().unique())
    cellmarker_per_type[cw_type] = markers & set(gene_symbols)

print("\nCellMarker markers per validated type:")
for t, m in cellmarker_per_type.items():
    print(f"  {t}: {len(m)} markers")

# ── Step 1: Global 500-gene identity set ──────────────────────────────
print("\n" + "=" * 70)
print("STEP 1: Global 500-gene identity set (variance across 35 centroids)")
print("=" * 70)

centroid_matrix = centroids.values  # (35, n_genes)
gene_variances = np.var(centroid_matrix, axis=0)  # variance across 35 types per gene
top500_idx = np.argsort(gene_variances)[::-1][:500]
top500_ensembl = [gene_ids[i] for i in top500_idx]
top500_symbols = set(ens_to_sym_map[g] for g in top500_ensembl)

print(f"Top 500 genes selected by centroid variance across 35 types")
print(f"Variance range: {gene_variances[top500_idx[0]]:.6f} (max) to {gene_variances[top500_idx[499]]:.6f} (500th)")

# Compare with original 6-type top 500
# Load 6-type centroids for comparison
centroids_6 = pd.read_csv("output/phase2/centroids_human.csv", index_col=0)
gene_ids_6 = list(centroids_6.columns)
gene_var_6 = np.var(centroids_6.values, axis=0)
top500_idx_6 = np.argsort(gene_var_6)[::-1][:500]
top500_symbols_6 = set(ens_to_sym_map.get(gene_ids_6[i], gene_ids_6[i]) for i in top500_idx_6)

overlap_500 = top500_symbols & top500_symbols_6
print(f"\nComparison with 6-type top 500:")
print(f"  6-type top 500: {len(top500_symbols_6)} unique symbols")
print(f"  35-type top 500: {len(top500_symbols)} unique symbols")
print(f"  Overlap: {len(overlap_500)} genes ({len(overlap_500)/500*100:.1f}%)")

# ── Step 2: Per-cell-type top 50 genes ────────────────────────────────
print("\n" + "=" * 70)
print("STEP 2: Per-cell-type top 50 genes (deviation from 35-type mean)")
print("=" * 70)

global_mean = np.mean(centroid_matrix, axis=0)  # mean across 35 types
print(f"Global mean centroid computed from {centroid_matrix.shape[0]} types")

per_type_top50 = {}
for ct in VALIDATED_TYPES:
    ct_centroid = centroids.loc[ct].values
    deviation = np.abs(ct_centroid - global_mean)
    top50_idx_ct = np.argsort(deviation)[::-1][:50]
    top50_syms = [ens_to_sym_map[gene_ids[i]] for i in top50_idx_ct]
    per_type_top50[ct] = set(top50_syms)
    print(f"  {ct}: top 50 deviation genes extracted")

# ── Step 3: CellMarker enrichment ─────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 3: CellMarker enrichment tests")
print("=" * 70)


def hypergeom_enrichment(test_genes, marker_genes, background_size):
    """One-sided hypergeometric test for enrichment."""
    # K = markers in background, n = test set size, k = observed overlap
    K = len(marker_genes)
    n = len(test_genes)
    k = len(test_genes & marker_genes)
    expected = n * K / background_size
    enrichment = k / expected if expected > 0 else 0.0
    # P(X >= k) = 1 - P(X <= k-1) = sf(k-1)
    p_val = stats.hypergeom.sf(k - 1, background_size, K, n)
    return k, expected, enrichment, p_val


# 3a. Global enrichment (500 genes vs all CellMarker)
print("\n--- 3a. Global 500-gene enrichment ---")
k_global, exp_global, enr_global, p_global = hypergeom_enrichment(
    top500_symbols, cellmarker_in_bg, n_genes
)
print(f"  Observed overlap: {k_global} / 500")
print(f"  Expected overlap: {exp_global:.2f}")
print(f"  Enrichment ratio: {enr_global:.3f}")
print(f"  p-value: {p_global:.4e}")
print(f"  ORIGINAL (6-type): 33/500, expected 7.58, ratio 4.355, p=1.12e-12")

# Identify overlapping genes
overlap_genes_global = sorted(top500_symbols & cellmarker_in_bg)
print(f"  Overlapping genes: {', '.join(overlap_genes_global[:20])}{'...' if len(overlap_genes_global) > 20 else ''}")

# 3b. Per-cell-type enrichment (top 50 deviation genes)
print("\n--- 3b. Per-cell-type enrichment (n=50 deviation genes) ---")
per_type_results = []
for ct in VALIDATED_TYPES:
    ct_genes = per_type_top50[ct]
    ct_markers = cellmarker_per_type[ct]
    k_ct, exp_ct, enr_ct, p_ct = hypergeom_enrichment(ct_genes, ct_markers, n_genes)
    passed = p_ct < 0.05 and enr_ct > 1.5
    overlap_ct = sorted(ct_genes & ct_markers)
    per_type_results.append({
        "cell_type": ct,
        "n_loading_genes": 50,
        "n_cellmarker_markers": len(ct_markers),
        "overlap": k_ct,
        "enrichment_ratio": round(enr_ct, 2),
        "p_value": p_ct,
        "pass": "PASS" if passed else "FAIL",
        "overlapping_genes": overlap_ct,
    })
    status = "PASS" if passed else "FAIL"
    print(f"  {ct}: overlap={k_ct}, enrichment={enr_ct:.2f}, p={p_ct:.4e} — {status}")
    if overlap_ct:
        print(f"    Genes: {', '.join(overlap_ct)}")

n_pass = sum(1 for r in per_type_results if r["pass"] == "PASS")
print(f"\n  Result: {n_pass}/6 cell types pass")
print(f"  ORIGINAL (6-type, n=50): 6/6 pass")

# 3c. Expression-matched background control
print("\n--- 3c. Expression-matched background control ---")
# For each of the 500 identity genes, find up to 10 non-identity genes
# within ±10% of mean expression level
mean_expr = np.mean(centroid_matrix, axis=0)  # mean expression per gene across 35 types
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
cellmarker_in_universe = cellmarker_in_bg & universe

print(f"  Expression-matched background: {len(matched_bg_symbols)} unique genes")
print(f"  Total universe: {universe_size} genes")
print(f"  CellMarker genes in universe: {len(cellmarker_in_universe)}")

k_matched, exp_matched, enr_matched, p_matched = hypergeom_enrichment(
    top500_symbols & universe, cellmarker_in_universe, universe_size
)
print(f"  Observed overlap: {k_matched}")
print(f"  Expected overlap: {exp_matched:.2f}")
print(f"  Enrichment ratio: {enr_matched:.3f}")
print(f"  p-value: {p_matched:.4e}")
print(f"  ORIGINAL (6-type): universe=2294, ratio=2.857, p=1.44e-10")

# ── Step 4: Save results ──────────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 4: Save results")
print("=" * 70)

results = {
    "metadata": {
        "centroid_source": str(CENTROID_35),
        "n_cell_types": 35,
        "n_genes": n_genes,
        "cellmarker_source": str(CELLMARKER_HUMAN),
        "comparison_note": "Rerun of issue073 CellMarker validation using 35-type centroids",
    },
    "gene_set_comparison": {
        "overlap_6type_vs_35type_top500": len(overlap_500),
        "overlap_pct": round(len(overlap_500) / 500 * 100, 1),
    },
    "human_primary_35type": {
        "observed_overlap": k_global,
        "expected_overlap": round(exp_global, 2),
        "enrichment_ratio": round(enr_global, 3),
        "p_value": p_global,
        "n_identity_genes": 500,
        "n_cellmarker_genes": len(cellmarker_in_bg),
        "n_background_genes": n_genes,
        "overlapping_genes": overlap_genes_global,
    },
    "human_expression_matched_35type": {
        "universe_size": universe_size,
        "n_cellmarker_in_universe": len(cellmarker_in_universe),
        "n_identity_genes": 500,
        "observed_overlap": k_matched,
        "expected_overlap": round(exp_matched, 2),
        "enrichment_ratio": round(enr_matched, 3),
        "p_value": p_matched,
        "n_matched_bg_genes": len(matched_bg_symbols),
    },
    "per_cell_type_35type": per_type_results,
    "comparison_with_original": {
        "original_global_enrichment": 4.355,
        "original_global_p": 1.12e-12,
        "original_per_type_pass": "6/6",
        "original_expression_matched_enrichment": 2.857,
        "original_expression_matched_p": 1.44e-10,
        "rerun_global_enrichment": round(enr_global, 3),
        "rerun_global_p": p_global,
        "rerun_per_type_pass": f"{n_pass}/6",
        "rerun_expression_matched_enrichment": round(enr_matched, 3),
        "rerun_expression_matched_p": p_matched,
    },
}

# Add aliases for downstream consumers (validate.py, generate_phase1_figures.py).
# These flatten selected fields without removing the original nested keys.
results["global_enrichment"] = {
    "enrichment": results["human_primary_35type"]["enrichment_ratio"],
    "p_value":    results["human_primary_35type"]["p_value"],
    "observed":   results["human_primary_35type"]["observed_overlap"],
    "expected":   results["human_primary_35type"]["expected_overlap"],
}
results["expression_matched"] = {
    "enrichment": results["human_expression_matched_35type"]["enrichment_ratio"],
    "p_value":    results["human_expression_matched_35type"]["p_value"],
    "observed":   results["human_expression_matched_35type"]["observed_overlap"],
}
results["per_type_pass"] = results["comparison_with_original"]["rerun_per_type_pass"]

with open(OUTPUT_DIR / "cellmarker_35type_results.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print(f"Saved: {OUTPUT_DIR / 'cellmarker_35type_results.json'}")

# Save top 500 gene list for reference
top500_df = pd.DataFrame({
    "ensembl_id": top500_ensembl,
    "gene_symbol": [ens_to_sym_map[g] for g in top500_ensembl],
    "centroid_variance_35type": [gene_variances[i] for i in top500_idx],
})
top500_df.to_csv(OUTPUT_DIR / "identity_genes_35type_top500.csv", index=False)
print(f"Saved: {OUTPUT_DIR / 'identity_genes_35type_top500.csv'}")

# ── Summary ───────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("SUMMARY — 35-type vs 6-type CellMarker Validation")
print("=" * 70)
print(f"\n{'Metric':<45} {'6-type (original)':>20} {'35-type (rerun)':>20}")
print("-" * 85)
print(f"{'Global enrichment fold':<45} {'4.355':>20} {enr_global:>20.3f}")
print(f"{'Global p-value':<45} {'1.12e-12':>20} {p_global:>20.4e}")
print(f"{'Observed / expected overlap':<45} {'33 / 7.58':>20} {f'{k_global} / {exp_global:.2f}':>20}")
print(f"{'Per-type pass rate':<45} {'6/6':>20} {f'{n_pass}/6':>20}")
print(f"{'Expression-matched enrichment':<45} {'2.857':>20} {enr_matched:>20.3f}")
print(f"{'Expression-matched p-value':<45} {'1.44e-10':>20} {p_matched:>20.4e}")
print(f"{'Gene set overlap (6 vs 35)':<45} {f'{len(overlap_500)}/500 ({len(overlap_500)/500*100:.1f}%)':>20}")

if enr_global > 4.355 * 1.1:
    verdict = "STRONGER"
elif enr_global > 4.355 * 0.9:
    verdict = "COMPARABLE"
else:
    verdict = "WEAKER"
print(f"\nOverall verdict: {verdict}")
