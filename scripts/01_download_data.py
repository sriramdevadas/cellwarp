#!/usr/bin/env python3
"""
CellWarp — Phase 1 Data Download Script

Downloads single-cell RNA-seq data from CZ CELLxGENE Census for human (Tabula
Sapiens) and mouse (Tabula Muris Senis), maps genes to shared 1:1 ortholog space
via Ensembl BioMart, and saves aligned .h5ad files for downstream QC and analysis.

Biology
-------
We pull five homologous cell types from two species to test whether cross-species
gene expression differences follow coherent geometric transformations (Procrustes
analysis in Phase 2). Raw counts are saved here; normalization happens in script 02.

Checkpoints
-----------
Intermediate results are saved after each major step. If the script crashes or is
interrupted, re-running will skip completed steps automatically.

    data/phase1/orthologs_human_mouse.csv   — BioMart 1:1 ortholog table
    data/phase1/human_raw.h5ad              — Raw human data (all genes)
    data/phase1/mouse_raw.h5ad              — Raw mouse data (all genes)
    data/phase1/human_aligned.h5ad          — Human data, shared ortholog space
    data/phase1/mouse_aligned.h5ad          — Mouse data, shared ortholog space

Usage:
    python scripts/01_download_data.py

Expected runtime:
    10–30 minutes (network-dependent). BioMart ~1-2 min, Census ~5-15 min/species.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Add project root to path for src imports
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import anndata as ad
import cellxgene_census

from cellwarp.data_loader import (
    CELL_TYPE_MAP,
    DATA_DIR,
    HUMAN_COLLECTION,
    HUMAN_ORGANISM,
    MOUSE_COLLECTION,
    MOUSE_ORGANISM,
    build_obs_value_filter,
    discover_cell_type_names,
    download_species_data,
    fetch_orthologs,
    filter_to_shared_orthologs,
    get_dataset_ids_for_collection,
    print_download_summary,
    save_h5ad_atomic,
    subsample_per_cell_type,
)


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    # ── Checkpoint paths ──
    ortholog_path = DATA_DIR / "orthologs_human_mouse.csv"
    human_raw_path = DATA_DIR / "human_raw.h5ad"
    mouse_raw_path = DATA_DIR / "mouse_raw.h5ad"
    human_aligned_path = DATA_DIR / "human_aligned.h5ad"
    mouse_aligned_path = DATA_DIR / "mouse_aligned.h5ad"

    # ==================================================================
    # STEP 1: Fetch human-mouse 1:1 orthologs
    # ==================================================================
    print("=" * 70)
    print("STEP 1/5: Fetching human-mouse 1:1 orthologs from Ensembl BioMart")
    print("=" * 70)

    orthologs = fetch_orthologs(cache_path=ortholog_path)
    print(f"  Ortholog table: {len(orthologs):,} gene pairs\n")

    # ==================================================================
    # STEP 2: Connect to Census, discover datasets and cell types
    # ==================================================================
    print("=" * 70)
    print("STEP 2/5: Connecting to CZ CELLxGENE Census")
    print("=" * 70)

    # Check if raw data is already cached — if so, we can skip Census entirely
    # for Steps 2-4 and go straight to ortholog filtering.
    human_adata: ad.AnnData | None = None
    mouse_adata: ad.AnnData | None = None

    if human_raw_path.exists() and mouse_raw_path.exists():
        print(f"  Both raw checkpoints found — loading from cache")
        print(f"  Loading {human_raw_path}...")
        human_adata = ad.read_h5ad(human_raw_path)
        print(f"  Human: {human_adata.n_obs:,} cells × {human_adata.n_vars:,} genes")
        print(f"  Loading {mouse_raw_path}...")
        mouse_adata = ad.read_h5ad(mouse_raw_path)
        print(f"  Mouse: {mouse_adata.n_obs:,} cells × {mouse_adata.n_vars:,} genes")
        print("  Skipping Census connection (Steps 2-4).\n")
    else:
        print("  Opening Census connection (version='2025-11-08')...")
        with cellxgene_census.open_soma(census_version="2025-11-08") as census:
            # -- Discover dataset IDs --
            print(f"\n  Looking up '{HUMAN_COLLECTION}' datasets...")
            human_dataset_ids = get_dataset_ids_for_collection(
                census, HUMAN_COLLECTION
            )
            print(f"\n  Looking up '{MOUSE_COLLECTION}' datasets...")
            mouse_dataset_ids = get_dataset_ids_for_collection(
                census, MOUSE_COLLECTION
            )

            # -- Discover cell type names --
            print(f"\n  Verifying cell type names in {HUMAN_ORGANISM}...")
            human_ct_map = discover_cell_type_names(
                census, HUMAN_ORGANISM, human_dataset_ids, CELL_TYPE_MAP
            )
            print(f"\n  Verifying cell type names in {MOUSE_ORGANISM}...")
            mouse_ct_map = discover_cell_type_names(
                census, MOUSE_ORGANISM, mouse_dataset_ids, CELL_TYPE_MAP
            )

            # Determine shared cell types (present in both species)
            shared_project_names = sorted(
                set(human_ct_map.keys()) & set(mouse_ct_map.keys())
            )
            if len(shared_project_names) < len(CELL_TYPE_MAP):
                missing = set(CELL_TYPE_MAP.keys()) - set(shared_project_names)
                print(
                    f"\n  WARNING: {len(shared_project_names)}/{len(CELL_TYPE_MAP)} "
                    f"cell types found in both species."
                )
                for m in sorted(missing):
                    print(f"    MISSING in at least one species: {m}")

            if len(shared_project_names) == 0:
                print("\n  FATAL: No cell types found in both species. Cannot proceed.")
                sys.exit(1)

            # Get Census names for confirmed types
            human_cell_names = [human_ct_map[p] for p in shared_project_names]
            mouse_cell_names = [mouse_ct_map[p] for p in shared_project_names]

            # ==============================================================
            # STEP 3: Download human data
            # ==============================================================
            print("\n" + "=" * 70)
            print("STEP 3/5: Downloading human data from Tabula Sapiens")
            print("=" * 70)

            if human_raw_path.exists():
                print(f"  Checkpoint found: {human_raw_path}")
                human_adata = ad.read_h5ad(human_raw_path)
                print(
                    f"  Loaded: {human_adata.n_obs:,} cells × "
                    f"{human_adata.n_vars:,} genes"
                )
            else:
                human_filter = build_obs_value_filter(
                    human_cell_names, human_dataset_ids
                )
                human_adata = download_species_data(
                    census, HUMAN_ORGANISM, human_filter
                )
                human_adata = subsample_per_cell_type(human_adata)
                save_h5ad_atomic(human_adata, human_raw_path)

            # ==============================================================
            # STEP 4: Download mouse data
            # ==============================================================
            print("\n" + "=" * 70)
            print("STEP 4/5: Downloading mouse data from Tabula Muris Senis")
            print("=" * 70)

            if mouse_raw_path.exists():
                print(f"  Checkpoint found: {mouse_raw_path}")
                mouse_adata = ad.read_h5ad(mouse_raw_path)
                print(
                    f"  Loaded: {mouse_adata.n_obs:,} cells × "
                    f"{mouse_adata.n_vars:,} genes"
                )
            else:
                mouse_filter = build_obs_value_filter(
                    mouse_cell_names, mouse_dataset_ids
                )
                mouse_adata = download_species_data(
                    census, MOUSE_ORGANISM, mouse_filter
                )
                mouse_adata = subsample_per_cell_type(mouse_adata)
                save_h5ad_atomic(mouse_adata, mouse_raw_path)

        # Census connection closed here
        print("\n  Census connection closed.")

    assert human_adata is not None and mouse_adata is not None

    # ==================================================================
    # STEP 5: Filter to shared 1:1 orthologs
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 5/5: Aligning to shared 1:1 ortholog gene space")
    print("=" * 70)

    if human_aligned_path.exists() and mouse_aligned_path.exists():
        print(f"  Checkpoints found — loading from cache")
        human_aligned = ad.read_h5ad(human_aligned_path)
        mouse_aligned = ad.read_h5ad(mouse_aligned_path)
        print(
            f"  Human aligned: {human_aligned.n_obs:,} cells × "
            f"{human_aligned.n_vars:,} genes"
        )
        print(
            f"  Mouse aligned: {mouse_aligned.n_obs:,} cells × "
            f"{mouse_aligned.n_vars:,} genes"
        )
    else:
        human_aligned, mouse_aligned = filter_to_shared_orthologs(
            human_adata, mouse_adata, orthologs
        )
        save_h5ad_atomic(human_aligned, human_aligned_path)
        save_h5ad_atomic(mouse_aligned, mouse_aligned_path)

    # ==================================================================
    # SUMMARY
    # ==================================================================
    print("\n" + "=" * 70)
    print("DOWNLOAD COMPLETE — SUMMARY")
    print("=" * 70)

    summary = print_download_summary(human_aligned, mouse_aligned)

    print("\n" + "-" * 70)
    print("PHASE 1 GATE CRITERIA (data download portion):")
    print(
        f"  Cells ≥500/type/species: "
        f"{'PASS' if summary['gate_cells_pass'] else 'FAIL'}"
    )
    print(
        f"  Shared genes ≥12,000:    "
        f"{'PASS' if summary['gate_genes_pass'] else 'FAIL'}"
    )
    print("-" * 70)

    if summary["gate_cells_pass"] and summary["gate_genes_pass"]:
        print("\nData download criteria MET. Proceed to script 02 (QC & normalize).")
    else:
        print("\nWARNING: Data download criteria NOT fully met.")
        print("Review the warnings above and consider adjusting cell type selections.")
        sys.exit(1)


if __name__ == "__main__":
    main()
