#!/usr/bin/env python3
"""
CellWarp — Phase 3: GO Enrichment Validation (Script 06)

Performs Gene Ontology over-representation analysis on the Procrustes residual
gene lists from Phase 2. For each of 6 cell types, the top genes driving
cell-type-specific cross-species divergence are tested against
GO_Biological_Process_2023 via gseapy's Enrichr interface.

Biology
-------
Phase 2 Procrustes analysis identified the genes driving cell-type-specific
divergence between human and mouse — what a global rotation + scaling could
NOT explain. If these residual gene sets are enriched for coherent biological
processes, it validates that the geometric signal reflects real evolutionary
biology (metabolic adaptation, immune system divergence, etc.) rather than
technical noise.

Phase 3 Gate Criterion
----------------------
≥3 of 6 cell types must have at least one GO term with adjusted p-value < 0.05.

Inputs
------
- ./output/phase2/procrustes_results.json (top_genes_per_cell_type)

Outputs (→ ./output/phase3/go_enrichment/)
-------
- enrichment_{cell_type}.csv        — Full enrichment results per cell type
- top_terms_all_cell_types.csv      — Summary table of top 10 terms per type
- go_enrichment_summary.json        — Gate evaluation + summary statistics
- dotplot_{cell_type}.png           — Dot plot per cell type
- summary_heatmap.png               — Heatmap of top terms across all types
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from cellwarp.enrichment import (
    evaluate_gate,
    interpret_enrichment,
    plot_dotplot,
    plot_summary_heatmap,
    run_enrichment_all_cell_types,
    save_enrichment_results,
)


def main() -> None:
    """Run GO enrichment analysis on Procrustes residual gene sets."""

    print("=" * 70)
    print("  CellWarp — Phase 3: GO Enrichment Validation")
    print("=" * 70)

    # ------------------------------------------------------------------
    # 1. Load residual top genes from Phase 2
    # ------------------------------------------------------------------
    print("\n[1/5] Loading Procrustes residual genes...")

    results_path = PROJECT_ROOT / "output" / "phase2" / "procrustes_results.json"
    if not results_path.exists():
        print(f"  ERROR: {results_path} not found. Run script 04 first.")
        sys.exit(1)

    with open(results_path) as f:
        procrustes_data = json.load(f)

    top_genes_data = procrustes_data["top_genes_per_cell_type"]
    cell_types = procrustes_data["cell_types"]

    # Extract gene lists per cell type
    top_genes_per_ct: dict[str, list[str]] = {}
    for ct in cell_types:
        genes = [g["gene"] for g in top_genes_data[ct]]
        top_genes_per_ct[ct] = genes
        print(f"  {ct:<45} {len(genes)} genes")

    print(f"\n  Total cell types: {len(cell_types)}")
    print(f"  Genes per cell type: {len(next(iter(top_genes_per_ct.values())))}")

    # ------------------------------------------------------------------
    # 2. Run GO enrichment per cell type
    # ------------------------------------------------------------------
    print("\n[2/5] Running GO enrichment (GO_Biological_Process_2023)...")

    enrichment_results = run_enrichment_all_cell_types(
        top_genes_per_ct,
        gene_sets="GO_Biological_Process_2023",
    )

    # ------------------------------------------------------------------
    # 3. Print human-readable interpretations
    # ------------------------------------------------------------------
    print("\n[3/5] Interpreting results...")

    for ct in cell_types:
        interpretation = interpret_enrichment(ct, enrichment_results[ct])
        print(interpretation)

    # ------------------------------------------------------------------
    # 4. Generate plots
    # ------------------------------------------------------------------
    print("\n[4/5] Generating plots...")

    output_dir = PROJECT_ROOT / "output" / "phase3" / "go_enrichment"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Dot plots per cell type
    for ct in cell_types:
        safe_name = ct.replace(" ", "_").replace(",", "").replace("+", "plus")
        plot_path = output_dir / f"dotplot_{safe_name}.png"
        desc = plot_dotplot(enrichment_results[ct], ct, plot_path)
        print(f"  {ct}: {desc}")

    # Summary heatmap
    heatmap_path = output_dir / "summary_heatmap.png"
    desc = plot_summary_heatmap(enrichment_results, heatmap_path)
    print(f"  Summary heatmap: {desc}")

    # ------------------------------------------------------------------
    # 5. Evaluate gate and save
    # ------------------------------------------------------------------
    print("\n[5/5] Evaluating Phase 3 GO enrichment gate...")

    gate_eval = evaluate_gate(enrichment_results)
    save_enrichment_results(enrichment_results, gate_eval, output_dir)

    # Print gate summary
    print("\n" + "=" * 70)
    print("  PHASE 3 GO ENRICHMENT GATE EVALUATION")
    print("=" * 70)
    print(f"  Criterion: ≥3 of {len(cell_types)} cell types with GO enrichment p_adj < 0.05")
    print(f"  Result: {gate_eval['n_cell_types_with_enrichment']}/{len(cell_types)} cell types pass")
    print()

    for ct, info in sorted(gate_eval["per_cell_type"].items()):
        status = "PASS" if info["n_significant"] > 0 else "FAIL"
        if info["top_term"]:
            term_short = info["top_term"].split(" (GO:")[0][:50]
            print(
                f"  [{status}] {ct:<45} "
                f"{info['n_significant']:>3} sig terms  "
                f"top: {term_short} (p={info['top_p_adj']:.2e})"
            )
        else:
            print(f"  [{status}] {ct:<45} {info['n_significant']:>3} sig terms")

    print()
    if gate_eval["gate_passed"]:
        print("  >>> GATE: PASSED — GO enrichment validates Procrustes residuals <<<")
    else:
        print("  >>> GATE: FAILED — insufficient GO enrichment across cell types <<<")
    print("=" * 70)


if __name__ == "__main__":
    main()
