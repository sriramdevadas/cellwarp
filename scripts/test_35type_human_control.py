#!/usr/bin/env python3
"""
CellWarp — 35-type Human-versus-Human Negative Control (ANALYSIS-C)

Compares Tabula Sapiens centroids against CellHint human centroids to measure
within-species geometric coherence as a baseline for the cross-species result.

Biology
-------
If cross-species Procrustes coherence (obs/null=0.522) merely reflects shared
biological organisation (cell types form a conserved constellation regardless of
species), then two independent human atlases should show comparable coherence.
The negative control quantifies this baseline.

Math
----
Same pipeline as primary analysis: joint PCA → Procrustes → permutation test
(10,000 iterations) → LOOCV.  Reference = Tabula Sapiens, target = CellHint Human.

Decision rule:
    obs/null >> 0.522  →  cross-species signal stronger than within-species (win)
    obs/null ≈ 0.522   →  cannot separate from within-species (problem)
    overlap < 20 types →  stop, underpowered
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
from sklearn.decomposition import PCA

from cellwarp.procrustes import (
    PCA_VARIANCE_THRESHOLD,
    RANDOM_SEED,
    N_PERMUTATIONS,
    pca_reduce_centroids,
    procrustes_align,
    permutation_test,
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

TS_CENTROIDS = PROJECT_ROOT / "output/phase2/scaled_35types/centroids_human_35.csv"
CH_CENTROIDS = PROJECT_ROOT / "output/validation/cellhint_replication/centroids_cellhint.csv"
OUTPUT_DIR = PROJECT_ROOT / "output/validation/human_control_35"


def run_loocv(h_cent: pd.DataFrame, m_cent: pd.DataFrame) -> list[dict]:
    """LOOCV: train on n-1 types, predict held-out position."""
    cell_types = sorted(h_cent.index.tolist())
    n_types = len(cell_types)
    results = []

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)

        for i, held_out in enumerate(cell_types):
            train_types = [ct for ct in cell_types if ct != held_out]
            h_train = h_cent.loc[train_types].values
            m_train = m_cent.loc[train_types].values
            combined = np.vstack([h_train, m_train])

            pca = PCA(
                n_components=PCA_VARIANCE_THRESHOLD,
                svd_solver="full",
                random_state=RANDOM_SEED,
            )
            combined_pca = pca.fit_transform(combined)
            n_train = len(train_types)
            h_pca = combined_pca[:n_train]
            m_pca = combined_pca[n_train:]

            with contextlib.redirect_stdout(io.StringIO()):
                proc = procrustes_align(h_pca, m_pca)

            held_h_pca = pca.transform(
                h_cent.loc[held_out].values.reshape(1, -1)
            ).flatten()
            held_m_pca = pca.transform(
                m_cent.loc[held_out].values.reshape(1, -1)
            ).flatten()

            held_m_centered = held_m_pca - proc.translation_target
            predicted = (
                proc.scaling * (held_m_centered @ proc.rotation)
                + proc.translation_ref
            )

            error = float(np.linalg.norm(predicted - held_h_pca))
            null_distances = np.array([
                np.linalg.norm(predicted - h_pca[j]) for j in range(n_train)
            ])
            null_dist = float(np.mean(null_distances))
            ratio = error / null_dist if null_dist > 0 else float("inf")

            results.append({
                "cell_type": held_out,
                "error": error,
                "null_mean": null_dist,
                "ratio": ratio,
                "correct": ratio < 1.0,
            })

    return results


def main() -> None:
    t_start = time.time()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("ANALYSIS-C — 35-TYPE HUMAN-VERSUS-HUMAN NEGATIVE CONTROL")
    print("=" * 70)

    # ---- Load centroids ----
    print("\n  Loading centroids...")
    ts = pd.read_csv(TS_CENTROIDS, index_col=0)
    ch = pd.read_csv(CH_CENTROIDS, index_col=0)
    print(f"  Tabula Sapiens: {ts.shape[0]} types × {ts.shape[1]:,} genes")
    print(f"  CellHint Human: {ch.shape[0]} types × {ch.shape[1]:,} genes")

    # ---- Find overlapping types ----
    shared_types = sorted(set(ts.index) & set(ch.index))
    ts_only = sorted(set(ts.index) - set(ch.index))
    ch_only = sorted(set(ch.index) - set(ts.index))

    print(f"\n  Overlapping types: {len(shared_types)}")
    print(f"  TS-only types ({len(ts_only)}): {ts_only}")
    print(f"  CellHint-only types ({len(ch_only)}): {ch_only}")

    if len(shared_types) < 20:
        msg = f"STOP: Only {len(shared_types)} overlapping types — below 20-type threshold."
        print(f"\n  {msg}")
        result = {
            "status": "STOPPED_UNDERPOWERED",
            "n_overlap": len(shared_types),
            "shared_types": shared_types,
            "ts_only": ts_only,
            "ch_only": ch_only,
            "reason": msg,
        }
        with open(OUTPUT_DIR / "human_control_35.json", "w") as f:
            json.dump(result, f, indent=2)
        print(f"\n  Saved to {OUTPUT_DIR / 'human_control_35.json'}")
        return

    # ---- Subset to shared types and genes ----
    shared_genes = sorted(set(ts.columns) & set(ch.columns))
    print(f"  Shared genes: {len(shared_genes):,}")

    ts_sub = ts.loc[shared_types, shared_genes]
    ch_sub = ch.loc[shared_types, shared_genes]

    n_types = len(shared_types)

    # ---- Joint PCA ----
    print(f"\n  Running joint PCA on {2 * n_types} centroids...")
    ts_pca, ch_pca, pca_model, cell_types = pca_reduce_centroids(
        ts_sub, ch_sub, variance_threshold=PCA_VARIANCE_THRESHOLD
    )
    n_components = pca_model.n_components_
    cumvar = np.cumsum(pca_model.explained_variance_ratio_)[-1]

    # ---- Procrustes alignment ----
    print("\n  Running Procrustes alignment (TS = reference, CellHint = target)...")
    proc = procrustes_align(ts_pca, ch_pca)

    # ---- Permutation test ----
    print(f"\n  Running permutation test ({N_PERMUTATIONS:,} iterations)...")
    p_value, null_dist = permutation_test(
        ts_pca, ch_pca, n_permutations=N_PERMUTATIONS, seed=RANDOM_SEED
    )

    obs_null_ratio = proc.distance / float(np.median(null_dist))
    null_2_5 = float(np.percentile(null_dist, 2.5))
    null_97_5 = float(np.percentile(null_dist, 97.5))

    # ---- LOOCV ----
    print(f"\n  Running LOOCV ({n_types} folds)...")
    loocv_results = run_loocv(ts_sub, ch_sub)
    n_correct = sum(1 for r in loocv_results if r["correct"])
    mean_ratio = np.mean([r["ratio"] for r in loocv_results])

    # ---- Summary ----
    print("\n" + "=" * 70)
    print("RESULTS SUMMARY")
    print("=" * 70)
    print(f"  N types:           {n_types}")
    print(f"  N shared genes:    {len(shared_genes):,}")
    print(f"  PCA components:    {n_components} ({cumvar*100:.1f}% variance)")
    print(f"  Procrustes dist:   {proc.distance:.6f}")
    print(f"  Scaling:           {proc.scaling:.6f}")
    print(f"  Null median:       {np.median(null_dist):.6f}")
    print(f"  Obs/null ratio:    {obs_null_ratio:.4f}")
    print(f"  p-value:           {p_value:.6f}")
    print(f"  Null 2.5-97.5%:    [{null_2_5:.4f}, {null_97_5:.4f}]")
    print(f"  LOOCV:             {n_correct}/{n_types} correct (mean ratio={mean_ratio:.4f})")
    print()

    # ---- Decision rule ----
    primary_obs_null = 0.522
    sixtype_obs_null = 0.607
    if obs_null_ratio > primary_obs_null * 1.2:
        decision = "WIN — H-vs-H obs/null substantially exceeds cross-species"
    elif obs_null_ratio <= primary_obs_null:
        decision = "PROBLEM — H-vs-H obs/null at or below cross-species"
    else:
        decision = "MODERATE — H-vs-H obs/null between cross-species and 1.2× cross-species"

    print(f"  Primary cross-species obs/null:  {primary_obs_null}")
    print(f"  6-type H-vs-H obs/null:          {sixtype_obs_null}")
    print(f"  This H-vs-H obs/null:            {obs_null_ratio:.4f}")
    print(f"  DECISION: {decision}")

    # ---- Save results ----
    result = {
        "analysis": "ANALYSIS-C: 35-type human-versus-human negative control",
        "date": time.strftime("%Y-%m-%d"),
        "atlases": {
            "reference": "Tabula Sapiens (primary human atlas)",
            "target": "CellHint Human (Xu et al., Cell 2023)",
        },
        "overlap": {
            "n_shared_types": n_types,
            "shared_types": cell_types,
            "ts_only_types": ts_only,
            "ch_only_types": ch_only,
            "n_shared_genes": len(shared_genes),
        },
        "procrustes": {
            "distance": float(proc.distance),
            "distance_squared": float(proc.distance_squared),
            "scaling": float(proc.scaling),
            "obs_null_ratio": float(obs_null_ratio),
            "p_value": float(p_value),
            "null_median": float(np.median(null_dist)),
            "null_mean": float(np.mean(null_dist)),
            "null_2_5_percentile": null_2_5,
            "null_97_5_percentile": null_97_5,
            "pca_components": int(n_components),
            "pca_cumvar": float(cumvar),
        },
        "loocv": {
            "n_correct": n_correct,
            "n_total": n_types,
            "fraction_correct": n_correct / n_types,
            "mean_ratio": float(mean_ratio),
            "per_type": loocv_results,
        },
        "comparison": {
            "primary_cross_species_obs_null": primary_obs_null,
            "sixtype_hvh_obs_null": sixtype_obs_null,
            "this_obs_null": float(obs_null_ratio),
            "decision": decision,
        },
        "parameters": {
            "n_permutations": N_PERMUTATIONS,
            "random_seed": RANDOM_SEED,
            "pca_variance_threshold": PCA_VARIANCE_THRESHOLD,
        },
        "runtime_seconds": time.time() - t_start,
    }

    np.save(OUTPUT_DIR / "null_distribution_hvh.npy", null_dist)

    with open(OUTPUT_DIR / "human_control_35.json", "w") as f:
        json.dump(result, f, indent=2)

    loocv_df = pd.DataFrame(loocv_results)
    loocv_df.to_csv(OUTPUT_DIR / "loocv_results_hvh.csv", index=False)

    print(f"\n  Results saved to {OUTPUT_DIR}/")
    print(f"  Runtime: {time.time() - t_start:.1f}s")


if __name__ == "__main__":
    main()
