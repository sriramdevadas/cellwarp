#!/usr/bin/env python3
"""
CellWarp — HVG-only Procrustes Robustness Check (ANALYSIS-D)

Tests whether the Procrustes coherence signal survives restriction to highly
variable genes (HVGs), addressing the reviewer concern that the full 16,959-gene
space may be padded with uninformative housekeeping genes.

Biology
-------
If geometric coherence is driven by cell-type-discriminating gene programs rather
than thousands of stably expressed housekeeping genes, the signal should survive
or strengthen when restricted to HVGs.  Conversely, if the signal depends on the
stabilising effect of many low-variance genes, HVG restriction would weaken it.

Math
----
HVGs are selected from the joint human+mouse centroid variance (mean of per-species
variance per gene), ranked descending.  Three thresholds are tested: top 2,000,
3,000, and 5,000 genes.  For each threshold the full Procrustes pipeline (joint PCA,
10,000-iteration permutation test) is run and the rigidity ranking is compared to
the full-space ranking by Spearman ρ.

Decision rule:
    All three p<0.05 and ranking ρ>0.7  →  signal robust to HVG filtering
    Signal weakens / disappears          →  flagged, not patched
"""

from __future__ import annotations

import contextlib
import io
import json
import sys
import time
import warnings
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA

from cellwarp.procrustes import (
    PCA_VARIANCE_THRESHOLD,
    RANDOM_SEED,
    N_PERMUTATIONS,
    pca_reduce_centroids,
    procrustes_align,
    permutation_test,
    compute_residual_vectors,
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

CENTROID_DIR = PROJECT_ROOT / "output/phase2/scaled_35types"
OUTPUT_DIR = PROJECT_ROOT / "output/validation/hvg_robustness"
FULL_SPACE_RESIDUALS = CENTROID_DIR / "residuals_ranked.csv"

HVG_THRESHOLDS = [2_000, 3_000, 5_000]


def select_hvgs_from_centroids(
    h_cent: pd.DataFrame,
    m_cent: pd.DataFrame,
    n_top: int,
) -> list[str]:
    """Select top-n HVGs by mean per-species variance across cell types.

    Uses variance of centroid expression across cell types (between-type variance)
    which captures genes that discriminate cell types — the same property HVG
    methods capture at the cell level.
    """
    h_var = h_cent.var(axis=0)
    m_var = m_cent.var(axis=0)
    mean_var = (h_var + m_var) / 2
    top_genes = mean_var.nlargest(n_top).index.tolist()
    return top_genes


def compute_rigidity_ranking(
    h_cent: pd.DataFrame,
    m_cent: pd.DataFrame,
) -> pd.Series:
    """Run Procrustes silently and return per-type residual magnitudes."""
    with contextlib.redirect_stdout(io.StringIO()):
        h_pca, m_pca, pca_model, cell_types = pca_reduce_centroids(
            h_cent, m_cent, variance_threshold=PCA_VARIANCE_THRESHOLD
        )
        proc = procrustes_align(h_pca, m_pca)
        residuals = compute_residual_vectors(proc, cell_types)

    magnitudes = {ct: float(np.linalg.norm(v)) for ct, v in residuals.items()}
    return pd.Series(magnitudes).sort_values()


def main() -> None:
    t_start = time.time()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("ANALYSIS-D — HVG-ONLY PROCRUSTES ROBUSTNESS CHECK")
    print("=" * 70)

    # ---- Load full-space centroids ----
    print("\n  Loading 35-type centroids (full 16,959-gene space)...")
    h_cent = pd.read_csv(CENTROID_DIR / "centroids_human_35.csv", index_col=0)
    m_cent = pd.read_csv(CENTROID_DIR / "centroids_mouse_35.csv", index_col=0)
    n_types = h_cent.shape[0]
    n_genes_full = h_cent.shape[1]
    print(f"  {n_types} types × {n_genes_full:,} genes")

    # ---- Load full-space rigidity ranking ----
    print("\n  Loading full-space rigidity ranking...")
    full_residuals = pd.read_csv(FULL_SPACE_RESIDUALS)
    full_ranking = full_residuals.set_index("cell_type")["residual_magnitude"]
    full_ranking = full_ranking.rank()
    print(f"  Full-space ranking loaded ({len(full_ranking)} types)")

    # ---- Run for each HVG threshold ----
    all_results = {}

    for n_hvg in HVG_THRESHOLDS:
        print(f"\n{'─' * 70}")
        print(f"  HVG threshold: top {n_hvg:,} genes")
        print(f"{'─' * 70}")

        # Select HVGs
        hvg_genes = select_hvgs_from_centroids(h_cent, m_cent, n_hvg)
        print(f"  Selected {len(hvg_genes)} HVGs by between-type variance")

        # Subset centroids
        h_hvg = h_cent[hvg_genes]
        m_hvg = m_cent[hvg_genes]
        print(f"  Centroid shape: {h_hvg.shape}")

        # Joint PCA
        print("  Running joint PCA...")
        h_pca, m_pca, pca_model, cell_types = pca_reduce_centroids(
            h_hvg, m_hvg, variance_threshold=PCA_VARIANCE_THRESHOLD
        )
        n_components = pca_model.n_components_
        cumvar = float(np.cumsum(pca_model.explained_variance_ratio_)[-1])

        # Procrustes
        print("  Running Procrustes alignment...")
        proc = procrustes_align(h_pca, m_pca)

        # Permutation test
        print(f"  Running permutation test ({N_PERMUTATIONS:,} iterations)...")
        p_value, null_dist = permutation_test(
            h_pca, m_pca, n_permutations=N_PERMUTATIONS, seed=RANDOM_SEED
        )

        obs_null_ratio = proc.distance / float(np.median(null_dist))

        # Rigidity ranking
        print("  Computing HVG rigidity ranking...")
        hvg_ranking_raw = compute_rigidity_ranking(h_hvg, m_hvg)
        hvg_ranks = hvg_ranking_raw.rank()

        # Spearman ρ against full-space ranking
        shared_types = sorted(set(full_ranking.index) & set(hvg_ranks.index))
        rho, rho_p = stats.spearmanr(
            full_ranking.loc[shared_types],
            hvg_ranks.loc[shared_types],
        )

        # Summary
        print(f"\n  --- Results for top {n_hvg:,} HVGs ---")
        print(f"  PCA components:    {n_components} ({cumvar*100:.1f}% variance)")
        print(f"  Procrustes dist:   {proc.distance:.6f}")
        print(f"  Scaling:           {proc.scaling:.6f}")
        print(f"  Obs/null ratio:    {obs_null_ratio:.4f}")
        print(f"  p-value:           {p_value:.6f}")
        print(f"  Ranking ρ vs full: {rho:.4f} (p={rho_p:.6f})")

        result = {
            "n_hvg": n_hvg,
            "n_genes_selected": len(hvg_genes),
            "pca_components": int(n_components),
            "pca_cumvar": float(cumvar),
            "procrustes_distance": float(proc.distance),
            "scaling": float(proc.scaling),
            "obs_null_ratio": float(obs_null_ratio),
            "p_value": float(p_value),
            "null_median": float(np.median(null_dist)),
            "null_mean": float(np.mean(null_dist)),
            "null_2_5_pct": float(np.percentile(null_dist, 2.5)),
            "null_97_5_pct": float(np.percentile(null_dist, 97.5)),
            "ranking_rho_vs_full": float(rho),
            "ranking_rho_p": float(rho_p),
            "per_type_residuals": {
                ct: float(hvg_ranking_raw[ct]) for ct in sorted(hvg_ranking_raw.index)
            },
        }
        all_results[str(n_hvg)] = result

        # Save null distribution
        np.save(
            OUTPUT_DIR / f"null_distribution_hvg{n_hvg}.npy",
            null_dist,
        )

    # ---- Overall decision ----
    print("\n" + "=" * 70)
    print("OVERALL DECISION")
    print("=" * 70)

    all_sig = all(all_results[str(n)]["p_value"] < 0.05 for n in HVG_THRESHOLDS)
    all_rho_high = all(
        all_results[str(n)]["ranking_rho_vs_full"] > 0.7 for n in HVG_THRESHOLDS
    )

    print(f"  Full-space obs/null:         0.522")
    for n_hvg in HVG_THRESHOLDS:
        r = all_results[str(n_hvg)]
        print(
            f"  HVG {n_hvg:>5,} obs/null: {r['obs_null_ratio']:.4f}  "
            f"p={r['p_value']:.6f}  ρ={r['ranking_rho_vs_full']:.4f}"
        )

    if all_sig and all_rho_high:
        decision = (
            "PASS — Signal robust to HVG filtering at all three thresholds. "
            "All p<0.05, all ranking ρ>0.7."
        )
    elif all_sig:
        decision = (
            "PARTIAL — All tests significant but some ranking ρ≤0.7. "
            "Signal survives but ranking is altered by HVG restriction."
        )
    else:
        decision = (
            "FAIL — Signal weakened or lost at one or more thresholds. "
            "Flagged for review; not patched."
        )

    print(f"\n  DECISION: {decision}")

    # ---- Save combined results ----
    output = {
        "analysis": "ANALYSIS-D: HVG-only Procrustes robustness check",
        "date": time.strftime("%Y-%m-%d"),
        "full_space_reference": {
            "n_genes": n_genes_full,
            "obs_null_ratio": 0.522,
            "p_value": 0.0001,
        },
        "hvg_results": all_results,
        "decision": decision,
        "all_significant": all_sig,
        "all_rho_above_0_7": all_rho_high,
        "parameters": {
            "hvg_thresholds": HVG_THRESHOLDS,
            "hvg_selection_method": "mean between-type variance across species",
            "n_permutations": N_PERMUTATIONS,
            "random_seed": RANDOM_SEED,
            "pca_variance_threshold": PCA_VARIANCE_THRESHOLD,
        },
        "runtime_seconds": time.time() - t_start,
    }

    with open(OUTPUT_DIR / "hvg_robustness.json", "w") as f:
        json.dump(output, f, indent=2)

    # Save comparison table as CSV
    rows = []
    for n_hvg in HVG_THRESHOLDS:
        r = all_results[str(n_hvg)]
        rows.append({
            "gene_space": f"HVG_{n_hvg}",
            "n_genes": r["n_genes_selected"],
            "obs_null_ratio": r["obs_null_ratio"],
            "p_value": r["p_value"],
            "scaling": r["scaling"],
            "ranking_rho_vs_full": r["ranking_rho_vs_full"],
            "ranking_rho_p": r["ranking_rho_p"],
        })
    rows.append({
        "gene_space": "Full_16959",
        "n_genes": n_genes_full,
        "obs_null_ratio": 0.522,
        "p_value": 0.0001,
        "scaling": 0.9023,
        "ranking_rho_vs_full": 1.0,
        "ranking_rho_p": 0.0,
    })
    pd.DataFrame(rows).to_csv(OUTPUT_DIR / "hvg_comparison.csv", index=False)

    print(f"\n  Results saved to {OUTPUT_DIR}/")
    print(f"  Runtime: {time.time() - t_start:.1f}s")


if __name__ == "__main__":
    main()
