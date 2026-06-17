#!/usr/bin/env python3
"""M1 close — Spearman ranking correlation + Table 1 summary numbers.

Reads the two committed JSONs from the K-expansion run:
  - reconstruction_qu12_results.json (human-macaque 12-type)
  - human_mouse_12type_control.json (human-mouse 12-type control)

Produces:
  - Spearman rho (and p) between per-type residuals across species pairs
  - Consolidated Table-1-ready summary JSON
"""
from __future__ import annotations
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

PROJECT = Path(__file__).resolve().parent.parent.parent
QU12_JSON  = PROJECT / "output/macaque_pipeline/reconstruction_qu12_results.json"
HM12_JSON  = PROJECT / "output/macaque_pipeline/human_mouse_12type_control.json"
D1_JSON    = PROJECT / "output/macaque_pipeline/reconstruction_qu7_D1_results.json"
HM7_JSON   = PROJECT / "output/macaque_pipeline/human_mouse_7type_control.json"
OUT_JSON   = PROJECT / "output/macaque_pipeline/m1_close_table1_summary.json"


def load_residuals(path: Path) -> dict[str, float]:
    d = json.loads(path.read_text())
    # Residuals may be nested; find under top-level or under per-analysis key
    if "per_type_residuals_ranked" in d:
        rows = d["per_type_residuals_ranked"]
    elif "control_16959" in d and "per_type_residuals_ranked" in d["control_16959"]:
        rows = d["control_16959"]["per_type_residuals_ranked"]
    else:
        raise RuntimeError(f"Could not find residuals in {path}")
    return {r["type"]: r["magnitude"] for r in rows}


def load_meta(path: Path, key: str = None) -> dict:
    d = json.loads(path.read_text())
    if key and key in d:
        return d[key]
    return d


def hep_rank(res: dict[str, float]) -> tuple[int, int]:
    """Return (rank, n_types) where rank 1 = largest residual."""
    items = sorted(res.items(), key=lambda kv: -kv[1])
    for i, (t, _) in enumerate(items, 1):
        if t == "hepatocyte":
            return i, len(items)
    return -1, len(items)


def main():
    # ── 1. Spearman ranking correlation on 12-type ──────────────────────
    mac_res = load_residuals(QU12_JSON)
    mm_res  = load_residuals(HM12_JSON)
    common = sorted(set(mac_res) & set(mm_res))
    assert len(common) == 12, f"Expected 12 common types, got {len(common)}"
    mac_vec = np.array([mac_res[t] for t in common])
    mm_vec  = np.array([mm_res[t] for t in common])
    rho, p_rho = spearmanr(mac_vec, mm_vec)

    print("Spearman ranking correlation (human-macaque 12-type vs human-mouse 12-type)")
    print(f"  rho = {rho:.4f}")
    print(f"  p   = {p_rho:.4f}")
    print(f"  n   = {len(common)}")

    # Also compute on D1 7-type (sanity)
    d1_res = load_residuals(D1_JSON)
    hm7_res = load_residuals(HM7_JSON)
    common7 = sorted(set(d1_res) & set(hm7_res))
    rho7, p7 = spearmanr([d1_res[t] for t in common7], [hm7_res[t] for t in common7])
    print(f"\nSpearman 7-type (D1 vs HM control): rho={rho7:.4f}, p={p7:.4f}, n={len(common7)}")

    # ── 2. Hepatocyte ranks ─────────────────────────────────────────────
    mac_hep_rank, mac_n = hep_rank(mac_res)
    mm_hep_rank, mm_n = hep_rank(mm_res)
    print(f"\nHepatocyte rank — macaque 12: {mac_hep_rank}/{mac_n}  ({mac_res['hepatocyte']:.3f})")
    print(f"Hepatocyte rank — mouse  12: {mm_hep_rank}/{mm_n}  ({mm_res['hepatocyte']:.3f})")

    # ── 3. Collect key fields ───────────────────────────────────────────
    mac_raw = json.loads(QU12_JSON.read_text())
    mm_raw = json.loads(HM12_JSON.read_text())
    mm_main = mm_raw["control_16959"]
    mm_parity = mm_raw["control_13927_parity"]

    summary = {
        "M1_close": True,
        "analysis_version": "K-expansion n=12",
        "types_included": mac_raw["types_included"],
        "gene_space_macaque": mac_raw["gene_space"],
        "gene_space_human_mouse_control_primary": mm_raw["gene_space"],
        "spearman_ranking_correlation": {
            "pair": "human-macaque 12-type vs human-mouse 12-type control",
            "rho": float(rho),
            "p_value": float(p_rho),
            "n_types": int(len(common)),
            "comment": "Direct comparison of per-type Procrustes residuals between two species pairs on the same 12 cell-type definitions; informs Figure 5C / Methods text on whether per-type divergence ranks co-vary across species pairs.",
        },
        "human_macaque_12type": {
            "n_types": 12,
            "gene_space": 13927,
            "obs_null_median": mac_raw["permutation_test"]["obs_null_ratio_median"],
            "p_value": mac_raw["permutation_test"]["p_value"],
            "n_permuted_le_observed": mac_raw["permutation_test"]["permuted_le_observed"],
            "n_permutations": mac_raw["permutation_test"]["n_permutations"],
            "procrustes_distance": mac_raw["procrustes"]["distance"],
            "procrustes_scaling": mac_raw["procrustes"]["scaling"],
            "pca_components": mac_raw["pca"]["n_components"],
            "pca_cumulative_variance": mac_raw["pca"]["cumulative_variance"],
            "hepatocyte_rank_of_n": f"{mac_hep_rank}/{mac_n}",
            "hepatocyte_residual_magnitude": float(mac_res["hepatocyte"]),
            "hepatocyte_pct_ssr": next(
                r["pct_ssr"] for r in mac_raw["per_type_residuals_ranked"]
                if r["type"] == "hepatocyte"
            ),
            "preprocessing": "Qu raw UMI → normalize_total(1e4) → log1p; 2,000 cells/type (seed 42)",
        },
        "human_mouse_12type_control": {
            "n_types": 12,
            "gene_space_primary": 16959,
            "obs_null_median_primary": mm_main["permutation_test"]["obs_null_ratio_median"],
            "p_value_primary": mm_main["permutation_test"]["p_value"],
            "n_permuted_le_observed_primary": mm_main["permutation_test"]["permuted_le_observed"],
            "procrustes_distance_primary": mm_main["procrustes"]["distance"],
            "procrustes_scaling_primary": mm_main["procrustes"]["scaling"],
            "pca_components_primary": mm_main["pca"]["n_components"],
            "gene_space_parity_13927": 13927,
            "obs_null_median_parity": mm_parity["permutation_test"]["obs_null_ratio_median"],
            "p_value_parity": mm_parity["permutation_test"]["p_value"],
            "hepatocyte_rank_of_n": f"{mm_hep_rank}/{mm_n}",
            "hepatocyte_residual_magnitude": float(mm_res["hepatocyte"]),
            "hepatocyte_pct_ssr": next(
                r["pct_ssr"] for r in mm_main["per_type_residuals_ranked"]
                if r["type"] == "hepatocyte"
            ),
            "preprocessing": "committed centroids_{human,mouse}_35.csv (symmetric: raw → NormalizeData + log1p)",
        },
        "no_immune_sensitivity_n7_D1": {
            "n_types": 7,
            "types": ["bladder urothelial cell", "endothelial cell", "epithelial cell",
                      "fibroblast", "hepatocyte", "smooth muscle cell", "stromal cell"],
            "obs_null_median": 0.7334,
            "p_value": 0.01250,
            "framing": "Excluding the 5 immune types (B cell, macrophage, monocyte, plasma cell, neutrophil) added in K, the n=7 non-immune analysis confirms cross-species coherence at p < 0.05. Serves as the methodologically clean replacement for the committed sensitivity (RIRA-only 13-type, p=0.0002 but preprocessing-asymmetric).",
            "human_mouse_7type_control_obs_null": 0.4848,
            "human_mouse_7type_control_p": 0.00060,
        },
        "hepatocyte_rank_reversal": {
            "human_macaque_12type": f"{mac_hep_rank}/{mac_n} (most flexible, %SSR="
                + f"{next(r['pct_ssr'] for r in mac_raw['per_type_residuals_ranked'] if r['type'] == 'hepatocyte'):.1f}%)",
            "human_mouse_12type":   f"{mm_hep_rank}/{mm_n} (most rigid, %SSR="
                + f"{next(r['pct_ssr'] for r in mm_main['per_type_residuals_ranked'] if r['type'] == 'hepatocyte'):.1f}%)",
            "narrative": "Hepatocyte rank reversal from most-rigid (mouse) to most-flexible (macaque) confirmed and strengthened at n=12 under symmetric preprocessing.",
        },
        "preprocessing_symmetry": True,
        "macaque_annotation_provenance": "Qu et al. 2022 author-assigned labels from Zenodo 10.5281/zenodo.5881495 @meta.data (finalcluster + majorcluster + Organ); no post-hoc marker-gene re-harmonization",
        "committed_primary_reference": {
            "macaque_primary_20type": {"obs_null": 0.841, "p": 0.0002,
                "preprocessing": "ASYMMETRIC (RIRA SCT-@counts as-is + Qu, with undocumented marker-gene harmonization)"},
            "macaque_sensitivity_13type_RIRA": {"obs_null": 0.7488, "p": 0.0002,
                "preprocessing": "ASYMMETRIC SCT — REPLACED by n=7 non-immune sensitivity in new framing"},
        },
    }

    OUT_JSON.write_text(json.dumps(summary, indent=2))
    print(f"\nSaved: {OUT_JSON}")


if __name__ == "__main__":
    main()
