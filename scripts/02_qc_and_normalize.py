#!/usr/bin/env python3
"""
CellWarp — Phase 1 QC and Normalization Script

Loads the aligned .h5ad files from Phase 1 data download, applies quality control
filters, normalizes expression, identifies highly variable genes, and computes
PCA + UMAP for visualization. Each species is processed independently.

Biology
-------
Before comparing gene expression across species, we must ensure data quality.
Low-quality cells (few detected genes) are removed. Expression counts are
normalized to a common scale so that differences between cells reflect biology,
not sequencing depth. UMAP plots confirm that cell types form distinct clusters.

Pipeline per species:
    1. Filter cells with <200 detected genes
    2. Compute QC metrics (genes/cell, counts/cell)
    3. Normalize: counts per 10,000 + log1p
    4. Select top 2,000 highly variable genes (HVGs)
    5. PCA: 50 components on HVG subset
    6. UMAP: 2D visualization
    7. Save UMAP plot and print statistics

Inputs:
    data/phase1/human_aligned.h5ad
    data/phase1/mouse_aligned.h5ad

Outputs:
    data/phase1/human_qc.h5ad              — QC'd + normalized human data
    data/phase1/mouse_qc.h5ad              — QC'd + normalized mouse data
    output/phase1_qc/umap_human.png        — UMAP plot (human)
    output/phase1_qc/umap_mouse.png        — UMAP plot (mouse)

Usage:
    python scripts/02_qc_and_normalize.py
"""

from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import anndata as ad

from cellwarp.data_loader import MIN_CELLS_PER_TYPE, save_h5ad_atomic
from cellwarp.qc import run_full_qc_pipeline


def main() -> None:
    data_dir = Path("./data/phase1")
    output_dir = Path("./output/phase1_qc")
    output_dir.mkdir(parents=True, exist_ok=True)

    human_aligned_path = data_dir / "human_aligned.h5ad"
    mouse_aligned_path = data_dir / "mouse_aligned.h5ad"
    human_qc_path = data_dir / "human_qc.h5ad"
    mouse_qc_path = data_dir / "mouse_qc.h5ad"

    # ==================================================================
    # Human QC
    # ==================================================================
    print("=" * 70)
    print("HUMAN — QC and Normalization")
    print("=" * 70)

    if human_qc_path.exists():
        print(f"  Checkpoint found: {human_qc_path}")
        human = ad.read_h5ad(human_qc_path)
        print(f"  Loaded: {human.n_obs:,} cells × {human.n_vars:,} genes")
    else:
        print(f"  Loading: {human_aligned_path}")
        human = ad.read_h5ad(human_aligned_path)
        human = run_full_qc_pipeline(human, "Human", output_dir)
        save_h5ad_atomic(human, human_qc_path)

    # ==================================================================
    # Mouse QC
    # ==================================================================
    print("\n" + "=" * 70)
    print("MOUSE — QC and Normalization")
    print("=" * 70)

    if mouse_qc_path.exists():
        print(f"  Checkpoint found: {mouse_qc_path}")
        mouse = ad.read_h5ad(mouse_qc_path)
        print(f"  Loaded: {mouse.n_obs:,} cells × {mouse.n_vars:,} genes")
    else:
        print(f"  Loading: {mouse_aligned_path}")
        mouse = ad.read_h5ad(mouse_aligned_path)
        mouse = run_full_qc_pipeline(mouse, "Mouse", output_dir)
        save_h5ad_atomic(mouse, mouse_qc_path)

    # ==================================================================
    # Summary
    # ==================================================================
    print("\n" + "=" * 70)
    print("QC COMPLETE — FINAL SUMMARY")
    print("=" * 70)

    print(f"\n  Human: {human.n_obs:,} cells × {human.n_vars:,} genes")
    print(f"  Mouse: {mouse.n_obs:,} cells × {mouse.n_vars:,} genes")

    # Gate check: ≥500 cells per type per species
    all_pass = True
    print(f"\n  {'Cell Type':<40} {'Human':>8} {'Mouse':>8} {'Gate':>6}")
    print(f"  {'-' * 64}")

    human_cts = human.obs["cell_type"].value_counts().to_dict()
    mouse_cts = mouse.obs["cell_type"].value_counts().to_dict()
    all_types = sorted(set(list(human_cts.keys()) + list(mouse_cts.keys())))

    for ct in all_types:
        h = human_cts.get(ct, 0)
        m = mouse_cts.get(ct, 0)
        gate = "PASS" if h >= MIN_CELLS_PER_TYPE and m >= MIN_CELLS_PER_TYPE else "FAIL"
        if gate == "FAIL":
            all_pass = False
        print(f"  {ct:<40} {h:>8,} {m:>8,} {gate:>6}")

    print(f"\n  Cells ≥{MIN_CELLS_PER_TYPE}/type/species after QC: "
          f"{'PASS' if all_pass else 'FAIL'}")

    # UMAP text description
    print(f"\n  UMAP plots saved to: {output_dir}/")
    print("  Text description of UMAP plots:")
    print("  - Check that each cell type forms a distinct cluster")
    print("  - Check that clusters do not overlap significantly")
    print("  - Immune types (CD4+ T, CD8+ T, B, macrophage) may cluster")
    print("    closer to each other than to hepatocytes/endothelial")
    print("  - Inspect plots visually to confirm: output/phase1_qc/umap_*.png")

    if all_pass:
        print("\n  QC criteria MET. Proceed to script 03 (SAMap validation).")
    else:
        print("\n  WARNING: Some cell types below 500-cell threshold after QC.")
        print("  Review the warnings above and consider adjusting QC parameters.")


if __name__ == "__main__":
    main()
