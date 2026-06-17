#!/usr/bin/env python3
"""Option K scoping — how many Qu author labels map cleanly to TS-35.

No Procrustes yet. Produces:
  - Enumeration of all 17 Qu finalcluster values + 7 majorcluster values
  - For each, map to Tabula Sapiens 35-type set (exists / granularity mismatch / no match)
  - Raw Qu cell counts per candidate target
  - Projected n for expanded macaque primary (conservative / moderate / aggressive)
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parent.parent.parent
QU_META = PROJECT / "data/macaque/qu_2022/qu_metadata.csv"
HUMAN_CSV = PROJECT / "output/phase2/scaled_35types/centroids_human_35.csv"

# D1 targets already mapped
D1_TYPES = {
    "bladder urothelial cell", "endothelial cell", "epithelial cell",
    "fibroblast", "hepatocyte", "smooth muscle cell", "stromal cell",
}

# Proposed mappings for additional Qu types. "status" ∈
# {"clean", "granularity", "organ-stratified", "no match"}.
PROPOSED = [
    # (qu_final, qu_major_or_None, target_ts35, status, notes)
    ("B cell",                "B cell",          "B cell",              "clean",
     "Qu fine = TS fine, both are 'B cell'"),
    ("Macrophage",            "Myeloid cell",    "macrophage",           "clean",
     "Qu fine = TS fine"),
    ("Monocyte",              "Myeloid cell",    "monocyte",             "clean-generic",
     "Qu is generic; TS has generic 'monocyte' AND subtypes (classical/intermediate/"
     "non-classical). Use TS 'monocyte' for parity; avoid subtype mismatch."),
    ("Plasma cell",           "Plasma cell",     "plasma cell",          "clean",
     "Both fine labels match"),
    ("Neutrophil",            "Myeloid cell",    "neutrophil",           "clean",
     "Qu fine = TS fine (TS also has 'granulocyte' broader; use 'neutrophil' tight)"),
    ("Dendritic cell",        "Myeloid cell",    "myeloid dendritic cell","granularity",
     "Qu doesn't split mDC/pDC; TS is 'myeloid dendritic cell' specifically. Only 689 cells"
     " — below most per-type gates."),
    ("Effector T cell",       "T cell",          "T cell",                "granularity",
     "Qu doesn't distinguish CD4/CD8 effector. TS has CD4+ T, CD8+ T, and generic 'T cell'. "
     "Map to 'T cell' (generic) with caveat, or skip."),
    ("NKT cell",              "T cell",          "mature NK T cell",      "granularity",
     "Qu 'NKT cell' at 38,139 is 22% of cells — far beyond invariant NKT biology. Likely a "
     "broad T/NK bucket. TS 'mature NK T cell' is narrower. High risk of mismatch; skip."),
    # Organ-stratified — possible but requires per-organ rules, not clean 1:1
    ("Secretory cell",        "Epithelial cell", "(organ-dependent)",    "organ-stratified",
     "Colon secretory → goblet; other organs → broader epithelial. Breaks one-to-one clean map."),
    ("Ciliated cell",         "Epithelial cell", "(no clean match)",      "no match",
     "Airway epithelial subtype; TS doesn't have 'ciliated cell' as a 35-type"),
    ("FibSmo cell",           "Stromal cell",    "(no clean match)",      "no match",
     "Hybrid label not in TS 35; already absorbed into 'stromal cell' via majorcluster"),
    ("Mast cell",              "Myeloid cell",    "(no clean match)",      "no match",
     "Not in TS 35 (DECISIONS.md notes mast cell absent from primary). Count 341 below gate."),
    ("Spermatid",             "Spermatid",       "(no clean match)",      "no match",
     "Testis germ cell, not in TS 35"),
]


def main():
    df = pd.read_csv(QU_META, low_memory=False)
    n_total = len(df)
    print("=" * 80)
    print(f"Option K — Qu expansion scoping ({n_total:,} cells total)")
    print("=" * 80)

    print("\n[1] Tabula Sapiens 35-type inventory")
    human = pd.read_csv(HUMAN_CSV, index_col=0)
    ts35 = sorted(human.index.tolist())
    print(f"  TS 35 types ({len(ts35)}):")
    for t in ts35:
        print(f"    - {t}")

    print("\n[2] Qu finalcluster inventory + counts")
    fc = df["finalcluster"].value_counts(dropna=False)
    print(f"  17 finalcluster values, {fc.sum():,} cells")
    for k, v in fc.items():
        in_d1 = "D1" if str(k) in {"Endothelial cell", "Epithelial cell", "Fibroblasts",
                                   "Smooth muscle cell"} else ""
        print(f"    {str(k):<25s} {v:>7,}  {in_d1}")

    print("\n[3] Qu majorcluster inventory + counts")
    mc = df["majorcluster"].value_counts(dropna=False)
    for k, v in mc.items():
        print(f"    {str(k):<25s} {v:>7,}")

    print("\n[4] D1 targets already mapped (recap)")
    for t in sorted(D1_TYPES):
        print(f"    - {t}")

    print("\n[5] Proposed additional mappings (Qu → TS 35)")
    print(f"  {'Qu label':<20s}  {'→':<3}  {'TS 35 target':<30s}  "
          f"{'Qu count':>10s}  {'Status':<20s}  Notes")
    # Qu count includes all Qu cells with that finalcluster (pre-organ filter)
    counts_fc = df["finalcluster"].value_counts().to_dict()
    rows_accum = []
    for qu_fc, _, tgt, status, notes in PROPOSED:
        n = counts_fc.get(qu_fc, 0)
        in_ts = tgt in ts35
        print(f"  {qu_fc:<20s}  {'→':<3}  {tgt:<30s}  {n:>10,}  {status:<20s}  {notes}")
        rows_accum.append({
            "qu_label": qu_fc, "ts35_target": tgt, "in_ts35": in_ts,
            "qu_count": n, "status": status, "notes": notes,
        })

    # Tally projected n under each stringency
    clean = [r for r in rows_accum if r["status"] == "clean"]
    clean_generic = [r for r in rows_accum if r["status"] == "clean-generic"]
    granularity = [r for r in rows_accum if r["status"] == "granularity"]

    print("\n[6] Projected n for expanded macaque primary")
    n_d1 = len(D1_TYPES)
    print(f"  D1 baseline: n = {n_d1} (Qu already mapped)")
    n_conservative = n_d1 + len(clean)
    print(f"\n  Conservative (D1 + clean fine-matches): n = {n_conservative}")
    for r in clean:
        print(f"    + {r['qu_label']:<20s} → {r['ts35_target']:<25s} ({r['qu_count']:,} cells)")

    n_moderate = n_conservative + len(clean_generic)
    print(f"\n  Moderate (+ clean-generic): n = {n_moderate}")
    for r in clean_generic:
        print(f"    + {r['qu_label']:<20s} → {r['ts35_target']:<25s} ({r['qu_count']:,} cells; note generic-to-generic)")

    # Aggressive = moderate + keep granularity-mismatch "Effector T cell → T cell"
    n_aggressive = n_moderate + len([r for r in granularity
                                     if r["ts35_target"] == "T cell"])
    print(f"\n  Aggressive (+ T cell generic from Effector T): n = {n_aggressive}")
    for r in granularity:
        if r["ts35_target"] == "T cell":
            print(f"    + {r['qu_label']:<20s} → {r['ts35_target']:<25s} ({r['qu_count']:,} cells; granularity caveat)")

    print("\n[7] Types flagged for exclusion")
    for r in rows_accum:
        if r["status"] in {"no match", "organ-stratified"} or (
            r["status"] == "granularity" and r["ts35_target"] != "T cell"
        ):
            print(f"    - {r['qu_label']:<20s}  (count {r['qu_count']:,}, {r['status']}, {r['ts35_target']})")

    # Per-organ check for clean types to identify if any are concentrated
    print("\n[8] Cross-organ distribution for clean fine-matches (sanity)")
    for r in clean + clean_generic:
        qu_fc = r["qu_label"]
        sub = df[df["finalcluster"] == qu_fc]
        top = sub["Organ"].value_counts().head(5)
        print(f"\n  {qu_fc} ({len(sub):,} total):")
        for org, n in top.items():
            print(f"    {org:<15s} {n:>6,}")

    # Save scoping CSV
    out_csv = PROJECT / "output/macaque_pipeline/option_k_scoping.csv"
    pd.DataFrame(rows_accum).to_csv(out_csv, index=False)
    print(f"\n  Scoping table: {out_csv}")


if __name__ == "__main__":
    main()
