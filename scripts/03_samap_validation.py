#!/usr/bin/env python3
"""
Script 03: SAMap Validation of Cell Type Pairings

Runs SAMap (Tarashansky et al. 2021) to independently identify cross-species
cell type correspondences between human (Tabula Sapiens) and mouse (Tabula Muris
Senis), then compares SAMap's top pairings to our 6 manually curated pairings.

Biology: If SAMap independently identifies the same cell type correspondences
we assigned manually, it confirms that our pairings reflect genuine biological
homology rather than arbitrary labeling. This is a sanity check — SAMap is NOT
part of our core Procrustes pipeline.

Usage:
    python scripts/03_samap_validation.py

Inputs:
    data/phase1/human_aligned.h5ad  — raw counts, ortholog-aligned
    data/phase1/mouse_aligned.h5ad  — raw counts, ortholog-aligned

Outputs:
    output/phase1_samap/samap_mapping_scores.csv   — full pairwise score table
    output/phase1_samap/samap_comparison_report.txt — human-readable report
    output/phase1_samap/samap_pairing_details.csv   — per-pairing comparison
    output/phase1_samap/samap_heatmap.png            — annotated heatmap
"""

import os
import sys
import time

# Add project root to path
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

from cellwarp.samap_validate import (
    MANUAL_PAIRINGS,
    compare_pairings,
    get_cell_type_scores,
    run_samap,
    save_results,
)

# --- Configuration ---
HUMAN_PATH = "data/phase1/human_aligned.h5ad"
MOUSE_PATH = "data/phase1/mouse_aligned.h5ad"
OUTPUT_DIR = "output/phase1_samap"
HU_ID = "hu"
MO_ID = "mo"
N_ITERS = 3


def main():
    print("=" * 60)
    print("CellWarp Phase 1 — SAMap Cell Type Validation")
    print("=" * 60)

    # Check inputs exist
    for path in [HUMAN_PATH, MOUSE_PATH]:
        if not os.path.exists(path):
            print(f"ERROR: Input file not found: {path}")
            sys.exit(1)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # --- Step 1: Run SAMap ---
    print("\n[Step 1/4] Running SAMap alignment...")
    t0 = time.time()
    sm = run_samap(
        HUMAN_PATH,
        MOUSE_PATH,
        hu_id=HU_ID,
        mo_id=MO_ID,
        n_iters=N_ITERS,
    )
    elapsed = time.time() - t0
    print(f"  SAMap completed in {elapsed:.0f}s")

    # --- Step 2: Extract mapping scores ---
    print("\n[Step 2/4] Extracting cell type mapping scores...")
    scores_df = get_cell_type_scores(sm, hu_id=HU_ID, mo_id=MO_ID)
    print("\nPairwise mapping scores (human x mouse):")
    print(scores_df.to_string(float_format=lambda x: f"{x:.4f}"))

    # --- Step 3: Compare to manual pairings ---
    print("\n[Step 3/4] Comparing SAMap pairings to manual pairings...")
    comparison = compare_pairings(scores_df, MANUAL_PAIRINGS)

    print(f"\n  Confirmed: {comparison['n_confirmed']}/{comparison['n_total']} "
          f"({comparison['fraction_confirmed']:.0%})")
    for d in comparison["details"]:
        status = "OK" if d["confirmed"] else "MISMATCH"
        print(f"    [{status}] {d['human_type']}: "
              f"manual={d['mouse_type_manual']}, "
              f"samap={d['mouse_type_samap']} "
              f"(score={d['score_manual']:.4f})")

    # --- Step 4: Save results ---
    print(f"\n[Step 4/4] Saving outputs to {OUTPUT_DIR}/")
    save_results(scores_df, comparison, OUTPUT_DIR, MANUAL_PAIRINGS)

    # --- Summary ---
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    if comparison["all_confirmed"]:
        print(f"  ALL {comparison['n_total']}/{comparison['n_total']} manual pairings "
              "confirmed by SAMap.")
        print("  Phase 1 SAMap validation gate: PASSED")
    else:
        mismatches = [d for d in comparison["details"] if not d["confirmed"]]
        print(f"  {comparison['n_confirmed']}/{comparison['n_total']} pairings confirmed.")
        print(f"  {len(mismatches)} mismatches detected:")
        for m in mismatches:
            print(f"    - {m['human_type']}: expected {m['mouse_type_manual']}, "
                  f"got {m['mouse_type_samap']}")
        print("  Review mismatches before proceeding to Phase 2.")

    print("\nOutputs saved to:")
    for f in sorted(os.listdir(OUTPUT_DIR)):
        if f != ".gitkeep":
            fpath = os.path.join(OUTPUT_DIR, f)
            size = os.path.getsize(fpath)
            print(f"  {fpath} ({size:,} bytes)")

    # Text description of the heatmap for non-visual review
    print("\nHeatmap description:")
    print("  The heatmap shows a 6x6 matrix of SAMap alignment scores.")
    print("  Rows = human cell types, columns = mouse cell types.")
    print("  Manual pairings are outlined with black borders.")
    if comparison["all_confirmed"]:
        print("  The highest score in each row falls on the diagonal (manual pairing),")
        print("  confirming that SAMap independently agrees with our assignments.")
    print("=" * 60)


if __name__ == "__main__":
    main()
