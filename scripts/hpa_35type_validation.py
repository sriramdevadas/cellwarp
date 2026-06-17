"""
Human Protein Atlas (HPA) held-out marker validation — held-out replication of
the CellMarker enrichment analysis (scripts/cellmarker_35type_rerun.py).

Methodology mirrors the CellMarker analysis exactly:
  - Same 35-type centroids → same global 500-gene identity set (top variance
    across 35 centroids) and same per-cell-type top-50 deviation genes.
  - Same one-sided hypergeometric test against the same gene-space background.
  - Same expression-matched negative control.

The only thing swapped is the marker database. HPA marker assignments are
derived from:
  - "RNA single cell type specific nCPM"   — single-cell tissue atlases
  - "RNA blood cell specific nTPM"         — blood-cell atlas (B / T cells)
  - Gated on "RNA single cell type specificity" ∈ {Cell type enriched,
    Group enriched, Cell type enhanced} or "RNA blood cell specificity" ∈
    {Immune cell enriched, Immune cell enhanced, Group enriched} — these are
    HPA's analogue of CellMarker's experimentally-validated marker filter.

HPA cell-type label → CellWarp validated type mapping (frozen below): all 6
CellWarp validated types are covered. Aliases include the natural tissue-resident
or memory/naive subsets for T/B cells and Kupffer cells for macrophages.

Output schema matches output/validation/cellmarker_35type_rerun/cellmarker_35type_results.json
so the two databases are directly comparable.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

CENTROID_35 = Path("output/phase2/scaled_35types/centroids_human_35.csv")
ORTHOLOGS = Path("data/phase1/orthologs_human_mouse.csv")
HPA_TSV = Path("data/validation/hpa/proteinatlas.tsv")
DOWNLOAD_DATE_FILE = Path("data/validation/hpa/download_date.txt")
OUTPUT_DIR = Path("analysis/validation/hpa")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Same 6 validated types as the CellMarker analysis
VALIDATED_TYPES = [
    "B cell",
    "CD4-positive, alpha-beta T cell",
    "CD8-positive, alpha-beta T cell",
    "endothelial cell",
    "hepatocyte",
    "macrophage",
]

# CellWarp validated-type → HPA cell-type labels.
# Source columns:
#   * "RNA single cell type specific nCPM"  (sc_) for non-blood cells
#   * "RNA blood cell specific nTPM"        (bc_) for B / T cell subtypes
HPA_MAP = {
    "B cell": {
        "bc_": {"memory B-cell", "naive B-cell"},
    },
    "CD4-positive, alpha-beta T cell": {
        "bc_": {"naive CD4 T-cell", "memory CD4 T-cell"},
    },
    "CD8-positive, alpha-beta T cell": {
        "bc_": {"naive CD8 T-cell", "memory CD8 T-cell"},
    },
    "endothelial cell": {
        "sc_": {"Vascular endothelial cells", "Lymphatic endothelial cells"},
    },
    "hepatocyte": {
        "sc_": {"Hepatocytes"},
    },
    "macrophage": {
        # Kupffer = liver-resident macrophage (matches CellMarker's
        # "Tissue-resident macrophage" expansion in cellmarker_35type_rerun.py)
        "sc_": {"Macrophages", "Kupffer cells"},
    },
}

SC_ENRICHED_CLASSES = {"Cell type enhanced", "Cell type enriched", "Group enriched"}
BC_ENRICHED_CLASSES = {"Immune cell enhanced", "Immune cell enriched", "Group enriched"}


def parse_specific_field(value: str) -> dict[str, float]:
    """Parse 'Hepatocytes: 220.7;B-cells: 12.5' → {'Hepatocytes': 220.7, 'B-cells': 12.5}."""
    out: dict[str, float] = {}
    if pd.isna(value):
        return out
    for entry in str(value).split(";"):
        if ":" in entry:
            label, val = entry.split(":", 1)
            try:
                out[label.strip()] = float(val.strip())
            except ValueError:
                out[label.strip()] = float("nan")
    return out


def build_hpa_markers(hpa_df: pd.DataFrame, gene_symbols_in_bg: set[str]) -> tuple[set[str], dict[str, set[str]]]:
    """Construct HPA marker sets.

    Returns:
        hpa_marker_pool : all HPA-marker genes (gene is enriched/enhanced in
                          any single cell type or blood cell type) intersected
                          with the centroid-matrix gene background.
        per_type        : CellWarp validated type → set of HPA-marker symbols
                          mapped from HPA labels via HPA_MAP, intersected with bg.
    """
    pool_global: set[str] = set()
    per_type: dict[str, set[str]] = {ct: set() for ct in VALIDATED_TYPES}

    sc_spec_col = "RNA single cell type specificity"
    sc_ncpm_col = "RNA single cell type specific nCPM"
    bc_spec_col = "RNA blood cell specificity"
    bc_ntpm_col = "RNA blood cell specific nTPM"

    for _, row in hpa_df.iterrows():
        gene = row["Gene"]
        if pd.isna(gene) or gene not in gene_symbols_in_bg:
            continue

        sc_spec = row.get(sc_spec_col)
        if pd.notna(sc_spec) and sc_spec in SC_ENRICHED_CLASSES:
            sc_labels = set(parse_specific_field(row.get(sc_ncpm_col)).keys())
            if sc_labels:
                pool_global.add(gene)
                for cw_type, src_map in HPA_MAP.items():
                    if "sc_" in src_map and sc_labels & src_map["sc_"]:
                        per_type[cw_type].add(gene)

        bc_spec = row.get(bc_spec_col)
        if pd.notna(bc_spec) and bc_spec in BC_ENRICHED_CLASSES:
            bc_labels = set(parse_specific_field(row.get(bc_ntpm_col)).keys())
            if bc_labels:
                pool_global.add(gene)
                for cw_type, src_map in HPA_MAP.items():
                    if "bc_" in src_map and bc_labels & src_map["bc_"]:
                        per_type[cw_type].add(gene)

    return pool_global, per_type


def hypergeom_enrichment(test_genes: set[str], marker_genes: set[str], background_size: int):
    """One-sided hypergeometric test for enrichment (matches CellMarker script)."""
    K = len(marker_genes)
    n = len(test_genes)
    k = len(test_genes & marker_genes)
    expected = n * K / background_size if background_size > 0 else 0.0
    enrichment = k / expected if expected > 0 else 0.0
    p_val = stats.hypergeom.sf(k - 1, background_size, K, n)
    return k, expected, enrichment, p_val


def main() -> None:
    download_date = (
        DOWNLOAD_DATE_FILE.read_text().strip()
        if DOWNLOAD_DATE_FILE.exists()
        else "unknown"
    )

    print("=" * 70)
    print("HPA Held-Out Marker Validation — 35-Type Centroids")
    print(f"HPA download date: {download_date}")
    print("=" * 70)

    centroids = pd.read_csv(CENTROID_35, index_col=0)
    assert centroids.shape[0] == 35, f"Expected 35 types, got {centroids.shape[0]}"
    print(f"\nCentroid matrix: {centroids.shape[0]} types × {centroids.shape[1]} genes")

    orthologs = pd.read_csv(ORTHOLOGS)
    ens_to_symbol = dict(zip(orthologs["human_ensembl_id"], orthologs["human_gene_name"]))

    gene_ids = list(centroids.columns)
    gene_symbols = [ens_to_symbol.get(g, g) for g in gene_ids]
    ens_to_sym_map = dict(zip(gene_ids, gene_symbols))
    gene_symbols_in_bg = set(gene_symbols)
    n_genes = len(gene_ids)
    print(f"Gene space: {n_genes} genes (ortholog background)")

    hpa = pd.read_csv(HPA_TSV, sep="\t", low_memory=False)
    print(f"HPA TSV: {hpa.shape[0]} genes × {hpa.shape[1]} columns")

    hpa_pool, hpa_per_type = build_hpa_markers(hpa, gene_symbols_in_bg)
    print(f"\nHPA marker pool (any cell-type enriched/enhanced): {len(hpa_pool)} genes in bg")
    print("HPA markers per validated type:")
    coverage = {}
    for ct, markers in hpa_per_type.items():
        coverage[ct] = len(markers)
        print(f"  {ct}: {len(markers)} markers")

    # ── Step 1: same global 500-gene identity set as CellMarker analysis ────
    centroid_matrix = centroids.values
    gene_variances = np.var(centroid_matrix, axis=0)
    top500_idx = np.argsort(gene_variances)[::-1][:500]
    top500_ensembl = [gene_ids[i] for i in top500_idx]
    top500_symbols = {ens_to_sym_map[g] for g in top500_ensembl}
    print(f"\nGlobal identity set: top 500 genes by centroid variance across 35 types")

    # ── Step 2: same per-type top-50 deviation genes ────────────────────────
    global_mean = np.mean(centroid_matrix, axis=0)
    per_type_top50: dict[str, set[str]] = {}
    for ct in VALIDATED_TYPES:
        ct_centroid = centroids.loc[ct].values
        deviation = np.abs(ct_centroid - global_mean)
        top50_idx_ct = np.argsort(deviation)[::-1][:50]
        per_type_top50[ct] = {ens_to_sym_map[gene_ids[i]] for i in top50_idx_ct}

    # ── Step 3: enrichment tests ────────────────────────────────────────────
    print("\n--- 3a. Global 500-gene enrichment vs HPA marker pool ---")
    k_global, exp_global, enr_global, p_global = hypergeom_enrichment(
        top500_symbols, hpa_pool, n_genes
    )
    print(f"  K (HPA markers in background): {len(hpa_pool)}")
    print(f"  Observed overlap: {k_global} / 500")
    print(f"  Expected overlap: {exp_global:.2f}")
    print(f"  Enrichment ratio: {enr_global:.3f}")
    print(f"  p-value: {p_global:.4e}")
    print(f"  CellMarker baseline (35-type): 257 K, 33 obs, ratio 4.486, p 2.10e-13")
    overlap_genes_global = sorted(top500_symbols & hpa_pool)

    print("\n--- 3b. Per-cell-type enrichment ---")
    per_type_results = []
    for ct in VALIDATED_TYPES:
        ct_genes = per_type_top50[ct]
        ct_markers = hpa_per_type[ct]
        k_ct, exp_ct, enr_ct, p_ct = hypergeom_enrichment(ct_genes, ct_markers, n_genes)
        passed = bool(p_ct < 0.05 and enr_ct > 1.5)
        overlap_ct = sorted(ct_genes & ct_markers)
        per_type_results.append({
            "cell_type": ct,
            "n_loading_genes": 50,
            "n_hpa_markers": len(ct_markers),
            "observed_overlap": int(k_ct),
            "expected_overlap": round(float(exp_ct), 4),
            "enrichment_ratio": round(float(enr_ct), 3),
            "p_value": float(p_ct),
            "pass": "PASS" if passed else "FAIL",
            "overlapping_genes": overlap_ct,
        })
        status = "PASS" if passed else "FAIL"
        print(f"  {ct}: overlap={k_ct}, enrichment={enr_ct:.2f}, p={p_ct:.4e} — {status}")
        if overlap_ct:
            print(f"    Genes: {', '.join(overlap_ct)}")

    n_pass = sum(1 for r in per_type_results if r["pass"] == "PASS")
    print(f"\n  Per-type result: {n_pass}/{len(VALIDATED_TYPES)} cell types pass")

    print("\n--- 3c. Expression-matched background control ---")
    mean_expr = np.mean(centroid_matrix, axis=0)
    identity_idx_set = set(top500_idx)
    matched_bg_idx: set[int] = set()
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

    matched_bg_symbols = {ens_to_sym_map[gene_ids[i]] for i in matched_bg_idx}
    universe = top500_symbols | matched_bg_symbols
    universe_size = len(universe)
    hpa_in_universe = hpa_pool & universe

    k_matched, exp_matched, enr_matched, p_matched = hypergeom_enrichment(
        top500_symbols & universe, hpa_in_universe, universe_size
    )
    print(f"  Universe size (identity ∪ matched-bg): {universe_size}")
    print(f"  HPA markers in universe: {len(hpa_in_universe)}")
    print(f"  Observed overlap: {k_matched}")
    print(f"  Expected overlap: {exp_matched:.2f}")
    print(f"  Enrichment ratio: {enr_matched:.3f}")
    print(f"  p-value: {p_matched:.4e}")
    print(f"  CellMarker baseline (35-type): universe ~2294, ratio 3.32, p 1.15e-12")

    # ── Save ────────────────────────────────────────────────────────────────
    results = {
        "metadata": {
            "centroid_source": str(CENTROID_35),
            "n_cell_types": 35,
            "n_genes": n_genes,
            "marker_source": "Human Protein Atlas",
            "marker_file": str(HPA_TSV),
            "hpa_download_date": download_date,
            "hpa_columns_used": [
                "RNA single cell type specificity",
                "RNA single cell type specific nCPM",
                "RNA blood cell specificity",
                "RNA blood cell specific nTPM",
            ],
            "hpa_enriched_classes_single_cell": sorted(SC_ENRICHED_CLASSES),
            "hpa_enriched_classes_blood_cell": sorted(BC_ENRICHED_CLASSES),
            "hpa_cellwarp_mapping": {k: {src: sorted(v) for src, v in m.items()}
                                     for k, m in HPA_MAP.items()},
            "test": "one-sided hypergeometric (scipy.stats.hypergeom.sf)",
            "comparison_note": "Held-out replication of CellMarker 35-type "
                               "validation with the same 50 deviation genes "
                               "per type and the same ortholog background.",
        },
        "coverage": {
            "validated_types_total": len(VALIDATED_TYPES),
            "validated_types_covered": sum(1 for ct in VALIDATED_TYPES
                                           if coverage[ct] > 0),
            "per_type_marker_counts": coverage,
        },
        "global_enrichment": {
            "K_markers_in_background": len(hpa_pool),
            "observed_overlap": int(k_global),
            "expected_overlap": round(float(exp_global), 4),
            "enrichment_ratio": round(float(enr_global), 3),
            "p_value": float(p_global),
            "n_identity_genes": 500,
            "n_background_genes": n_genes,
            "overlapping_genes": overlap_genes_global,
        },
        "expression_matched_control": {
            "universe_size": universe_size,
            "n_markers_in_universe": len(hpa_in_universe),
            "n_identity_genes": 500,
            "observed_overlap": int(k_matched),
            "expected_overlap": round(float(exp_matched), 4),
            "enrichment_ratio": round(float(enr_matched), 3),
            "p_value": float(p_matched),
            "n_matched_bg_genes": len(matched_bg_symbols),
        },
        "per_cell_type": per_type_results,
        "cellmarker_baseline_for_comparison": {
            "global_K": 257,
            "global_enrichment_ratio": 4.486,
            "global_p_value": 2.10e-13,
            "per_type_pass": "5/6",
            "expression_matched_ratio": 3.32,
            "expression_matched_p": 1.15e-12,
            "note": "Numbers copied from output/validation/cellmarker_35type_rerun/cellmarker_35type_results.json",
        },
    }

    out_path = OUTPUT_DIR / "hpa_35type_results.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nSaved: {out_path}")

    summary = OUTPUT_DIR / "hpa_35type_summary.md"
    lines = [
        "# HPA Held-Out Marker Validation — 35-Type Summary",
        "",
        f"HPA download date: {download_date}",
        "",
        "## Pooled (global) enrichment",
        "",
        f"- K (HPA markers in background): **{len(hpa_pool)}**",
        f"- Observed overlap: **{k_global}** / 500",
        f"- Expected overlap: {exp_global:.2f}",
        f"- Fold-enrichment: **{enr_global:.3f}**",
        f"- p-value: **{p_global:.4e}**",
        "",
        "## Expression-matched control",
        "",
        f"- Universe size: {universe_size}",
        f"- HPA markers in universe: {len(hpa_in_universe)}",
        f"- Observed overlap: **{k_matched}** / {len(top500_symbols & universe)}",
        f"- Expected overlap: {exp_matched:.2f}",
        f"- Fold-enrichment: **{enr_matched:.3f}**",
        f"- p-value: **{p_matched:.4e}**",
        "",
        "## Per-cell-type",
        "",
        "| Cell type | HPA markers | Overlap | Fold | p | Result |",
        "|---|---|---|---|---|---|",
    ]
    for r in per_type_results:
        lines.append(
            f"| {r['cell_type']} | {r['n_hpa_markers']} | {r['observed_overlap']}/50 | "
            f"{r['enrichment_ratio']:.2f} | {r['p_value']:.2e} | {r['pass']} |"
        )
    lines.append("")
    lines.append("## Coverage")
    lines.append("")
    lines.append("HPA covers all 6 validated cell types via the following labels:")
    lines.append("")
    for cw_ct, src_map in HPA_MAP.items():
        labels = []
        for src, lbls in src_map.items():
            labels.extend(sorted(lbls))
        lines.append(f"- **{cw_ct}** ← {', '.join(repr(l) for l in labels)}")
    summary.write_text("\n".join(lines))
    print(f"Saved: {summary}")


if __name__ == "__main__":
    main()
