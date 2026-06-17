#!/usr/bin/env python3
"""
CellWarp — HCA-vs-Tabula Within-Human Centroid Stability Diagnostic

Tests whether the T1-A replication failure (p=0.542) reflects:
  (A) MCA microwell-seq technology too sparse for centroid geometry, OR
  (B) Tabula Sapiens batch effect drives the primary result.

Method: Compare HCA human centroids to Tabula Sapiens human centroids using the
identical Procrustes pipeline (same gene space, same PCA, same permutation test).
This is NOT a cross-species test — it is a within-human centroid stability check.

If HCA-vs-Tabula centroids are geometrically similar (obs/null ≈ 1.0), the Tabula
atlas is geometrically idiosyncratic even within humans — Hypothesis B supported.
If HCA-vs-Tabula centroids align significantly (obs/null < 0.8), Tabula geometry
is stable within humans and MCA is the problem — Hypothesis A supported.

Also runs Tabula human-vs-mouse on the same restricted type set for direct comparison.

Output: output/validation/hca_centroid_comparison/
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Add project root to path for src imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import anndata as ad
from cellwarp.procrustes import (
    compute_centroids,
    pca_reduce_centroids,
    procrustes_align,
    permutation_test,
    compute_residual_vectors,
    RANDOM_SEED,
)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

OUTPUT_DIR = Path("output/validation/hca_centroid_comparison")
HCA_PATH = Path("data/replication/hca_t1a.h5ad")
TABULA_HUMAN_PATH = Path("data/phase1/human_qc.h5ad")
TABULA_MOUSE_PATH = Path("data/phase1/mouse_qc.h5ad")

N_PERMUTATIONS = 10_000
VARIANCE_THRESHOLD = 0.95


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # ==================================================================
    # STEP 1: Load data
    # ==================================================================
    print("=" * 70)
    print("STEP 1: Load datasets")
    print("=" * 70)

    print("\nLoading HCA human (non-Tabula, pooled Census)...")
    hca = ad.read_h5ad(HCA_PATH)
    print(f"  Shape: {hca.shape}")
    hca_ct_col = "our_cell_type_label"
    hca_types = set(hca.obs[hca_ct_col].unique())
    print(f"  Cell types ({len(hca_types)}): {sorted(hca_types)}")

    print("\nLoading Tabula Sapiens human...")
    tabula_h = ad.read_h5ad(TABULA_HUMAN_PATH)
    print(f"  Shape: {tabula_h.shape}")
    tabula_ct_col = "cell_type"
    tabula_h_types = set(tabula_h.obs[tabula_ct_col].unique())
    print(f"  Cell types ({len(tabula_h_types)}): {sorted(tabula_h_types)}")

    print("\nLoading Tabula Muris Senis mouse...")
    tabula_m = ad.read_h5ad(TABULA_MOUSE_PATH)
    print(f"  Shape: {tabula_m.shape}")
    tabula_m_types = set(tabula_m.obs[tabula_ct_col].unique())
    print(f"  Cell types ({len(tabula_m_types)}): {sorted(tabula_m_types)}")

    # ==================================================================
    # STEP 2: Find shared cell types across all three datasets
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 2: Find shared cell types")
    print("=" * 70)

    # Shared between HCA and Tabula human (for the HCA-vs-Tabula comparison)
    shared_hca_tabula = sorted(hca_types & tabula_h_types)
    print(f"\nHCA ∩ Tabula human: {len(shared_hca_tabula)} types")
    for ct in shared_hca_tabula:
        n_hca = (hca.obs[hca_ct_col] == ct).sum()
        n_tab = (tabula_h.obs[tabula_ct_col] == ct).sum()
        print(f"  {ct}: HCA={n_hca}, Tabula={n_tab}")

    # Shared across all three (for comparable Tabula H-vs-M)
    shared_all = sorted(hca_types & tabula_h_types & tabula_m_types)
    print(f"\nHCA ∩ Tabula human ∩ Tabula mouse: {len(shared_all)} types")
    for ct in shared_all:
        n_hca = (hca.obs[hca_ct_col] == ct).sum()
        n_th = (tabula_h.obs[tabula_ct_col] == ct).sum()
        n_tm = (tabula_m.obs[tabula_ct_col] == ct).sum()
        print(f"  {ct}: HCA={n_hca}, Tabula_H={n_th}, Tabula_M={n_tm}")

    # Use the shared_all set for both comparisons (direct comparability)
    n_types = len(shared_all)
    print(f"\nUsing {n_types} types for both comparisons (direct comparability)")

    # Verify gene spaces match
    assert list(hca.var_names) == list(tabula_h.var_names), "Gene space mismatch HCA vs Tabula human"
    assert list(tabula_h.var_names) == list(tabula_m.var_names), "Gene space mismatch Tabula human vs mouse"
    gene_names = list(hca.var_names)
    n_genes = len(gene_names)
    print(f"Gene space: {n_genes} shared ortholog genes (verified identical)")

    # Filter to shared types
    hca_sub = hca[hca.obs[hca_ct_col].isin(shared_all)].copy()
    tabula_h_sub = tabula_h[tabula_h.obs[tabula_ct_col].isin(shared_all)].copy()
    tabula_m_sub = tabula_m[tabula_m.obs[tabula_ct_col].isin(shared_all)].copy()

    # Standardize column name for compute_centroids
    hca_sub.obs["cell_type"] = hca_sub.obs[hca_ct_col]

    # ==================================================================
    # STEP 3: Compute centroids
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 3: Compute centroids (mean expression over 16,959 genes)")
    print("=" * 70)

    print("\n--- HCA human centroids ---")
    hca_centroids = compute_centroids(hca_sub, "cell_type")

    print("\n--- Tabula Sapiens human centroids ---")
    tabula_h_centroids = compute_centroids(tabula_h_sub, "cell_type")

    print("\n--- Tabula Muris Senis mouse centroids ---")
    tabula_m_centroids = compute_centroids(tabula_m_sub, "cell_type")

    # ==================================================================
    # STEP 4: Comparison A — HCA human vs Tabula human
    # ==================================================================
    print("\n" + "=" * 70)
    print("COMPARISON A: HCA human vs Tabula Sapiens human")
    print("  (within-human centroid stability check)")
    print("=" * 70)

    print("\n--- PCA on combined HCA + Tabula human centroids ---")
    hca_pca, tab_h_pca, pca_a, types_a = pca_reduce_centroids(
        hca_centroids, tabula_h_centroids, variance_threshold=VARIANCE_THRESHOLD
    )

    print("\n--- Procrustes: Tabula → HCA (align Tabula onto HCA reference) ---")
    # procrustes_align(X=reference, Y=target) → aligns Y onto X
    result_a = procrustes_align(hca_pca, tab_h_pca)

    print("\n--- Permutation test (10,000 iterations) ---")
    p_a, null_a = permutation_test(hca_pca, tab_h_pca, N_PERMUTATIONS, RANDOM_SEED)

    obs_null_a = result_a.distance / np.median(null_a)

    print("\n--- Per-type residuals ---")
    residuals_a = compute_residual_vectors(result_a, types_a)
    residual_mags_a = {ct: float(np.linalg.norm(residuals_a[ct])) for ct in types_a}

    # ==================================================================
    # STEP 5: Comparison B — Tabula human vs Tabula mouse (same N types)
    # ==================================================================
    print("\n" + "=" * 70)
    print("COMPARISON B: Tabula human vs Tabula mouse (same type set)")
    print("  (primary result restricted to same N types for comparability)")
    print("=" * 70)

    print("\n--- PCA on combined Tabula human + mouse centroids ---")
    tab_h_pca_b, tab_m_pca_b, pca_b, types_b = pca_reduce_centroids(
        tabula_h_centroids, tabula_m_centroids, variance_threshold=VARIANCE_THRESHOLD
    )

    print("\n--- Procrustes: Tabula mouse → Tabula human ---")
    result_b = procrustes_align(tab_h_pca_b, tab_m_pca_b)

    print("\n--- Permutation test (10,000 iterations) ---")
    p_b, null_b = permutation_test(tab_h_pca_b, tab_m_pca_b, N_PERMUTATIONS, RANDOM_SEED)

    obs_null_b = result_b.distance / np.median(null_b)

    print("\n--- Per-type residuals ---")
    residuals_b = compute_residual_vectors(result_b, types_b)
    residual_mags_b = {ct: float(np.linalg.norm(residuals_b[ct])) for ct in types_b}

    # ==================================================================
    # STEP 6: Summary comparison table
    # ==================================================================
    print("\n" + "=" * 70)
    print("SUMMARY COMPARISON TABLE")
    print("=" * 70)

    header = f"{'Comparison':<45} {'n types':>7} {'p':>8} {'obs/null':>9} {'scaling':>8} {'distance':>10}"
    print(header)
    print("-" * len(header))
    print(
        f"{'HCA human vs Tabula human':<45} {n_types:>7} {p_a:>8.4f} {obs_null_a:>9.3f} "
        f"{result_a.scaling:>8.3f} {result_a.distance:>10.2f}"
    )
    print(
        f"{'Tabula human vs Tabula mouse (same N types)':<45} {n_types:>7} {p_b:>8.4f} {obs_null_b:>9.3f} "
        f"{result_b.scaling:>8.3f} {result_b.distance:>10.2f}"
    )

    # Per-type residuals comparison
    print(f"\n{'Cell Type':<45} {'HCA-vs-Tab resid':>17} {'Tab H-vs-M resid':>17}")
    print("-" * 80)
    for ct in types_a:
        print(f"{ct:<45} {residual_mags_a[ct]:>17.4f} {residual_mags_b[ct]:>17.4f}")

    # ==================================================================
    # STEP 7: Interpretation
    # ==================================================================
    print("\n" + "=" * 70)
    print("INTERPRETATION")
    print("=" * 70)

    if obs_null_a >= 0.95:
        interpretation = (
            "HYPOTHESIS B SUPPORTED: HCA-vs-Tabula is near-null (obs/null >= 0.95). "
            "Tabula centroids are geometrically idiosyncratic even within humans. "
            "Primary cross-species result is in serious trouble."
        )
        verdict = "B_SUPPORTED"
    elif obs_null_a < 0.80 and p_a < 0.05:
        interpretation = (
            "HYPOTHESIS A SUPPORTED: HCA-vs-Tabula is significant (obs/null < 0.80, p < 0.05). "
            "Human centroids are geometrically stable across independent atlases. "
            "MCA failure is technology-specific. Primary result survives — need 10x mouse replication."
        )
        verdict = "A_SUPPORTED"
    else:
        interpretation = (
            f"AMBIGUOUS: obs/null={obs_null_a:.3f}, p={p_a:.4f}. "
            "Does not clearly support either hypothesis. "
            "Document carefully, do not over-interpret."
        )
        verdict = "AMBIGUOUS"

    print(f"\n  {interpretation}")
    print(f"\n  Verdict: {verdict}")

    # ==================================================================
    # STEP 8: Save outputs
    # ==================================================================
    print("\n" + "=" * 70)
    print("STEP 8: Save outputs")
    print("=" * 70)

    results = {
        "diagnostic": "HCA-vs-Tabula within-human centroid stability check",
        "date": "2026-03-15",
        "purpose": "Distinguish Hypothesis A (MCA sparse) vs B (Tabula batch effect) after T1-A failure (p=0.542)",
        "comparison_a": {
            "name": "HCA human vs Tabula Sapiens human",
            "type": "within-human centroid stability",
            "n_types": n_types,
            "cell_types": types_a,
            "n_genes": n_genes,
            "p_value": float(p_a),
            "distance": float(result_a.distance),
            "obs_null_ratio": float(obs_null_a),
            "scaling": float(result_a.scaling),
            "null_median": float(np.median(null_a)),
            "null_mean": float(np.mean(null_a)),
            "n_permutations": N_PERMUTATIONS,
            "pca_components": int(pca_a.n_components_),
            "pca_variance": float(np.sum(pca_a.explained_variance_ratio_)),
            "per_type_residuals": {
                ct: {
                    "magnitude": residual_mags_a[ct],
                    "pct_ssr": float(residual_mags_a[ct]**2 / result_a.distance_squared * 100),
                }
                for ct in types_a
            },
        },
        "comparison_b": {
            "name": "Tabula human vs Tabula mouse (same N types)",
            "type": "cross-species (primary result restricted)",
            "n_types": n_types,
            "cell_types": types_b,
            "n_genes": n_genes,
            "p_value": float(p_b),
            "distance": float(result_b.distance),
            "obs_null_ratio": float(obs_null_b),
            "scaling": float(result_b.scaling),
            "null_median": float(np.median(null_b)),
            "null_mean": float(np.mean(null_b)),
            "n_permutations": N_PERMUTATIONS,
            "pca_components": int(pca_b.n_components_),
            "pca_variance": float(np.sum(pca_b.explained_variance_ratio_)),
            "per_type_residuals": {
                ct: {
                    "magnitude": residual_mags_b[ct],
                    "pct_ssr": float(residual_mags_b[ct]**2 / result_b.distance_squared * 100),
                }
                for ct in types_b
            },
        },
        "interpretation": {
            "verdict": verdict,
            "detail": interpretation,
            "decision_logic": {
                "near_null": "obs/null >= 0.95 → Hypothesis B (Tabula idiosyncratic)",
                "significant": "obs/null < 0.80 AND p < 0.05 → Hypothesis A (MCA sparse)",
                "ambiguous": "otherwise → ambiguous",
            },
        },
        "random_seed": RANDOM_SEED,
    }

    # Save JSON
    results_path = OUTPUT_DIR / "hca_centroid_comparison.json"
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"  Results: {results_path}")

    # Save null distributions
    np.save(OUTPUT_DIR / "null_a_hca_vs_tabula.npy", null_a)
    # NOTE: When the cell-type triple-intersection reduces to the same 6 types
    # as the sibling sun2023_replication workflow's comparison_b, the saved
    # null_b_tabula_hvm.npy will be byte-identical with the sibling's
    # null_b_tabula.npy. Both are diagnostic preserves; neither feeds Table 1
    # or manuscript claims. See R21 MD5-16 investigation.
    np.save(OUTPUT_DIR / "null_b_tabula_hvm.npy", null_b)
    print(f"  Null distributions saved")

    # ==================================================================
    # STEP 9: Generate plots
    # ==================================================================
    print("\n--- Generating plots ---")

    # Plot 1: Side-by-side null distribution histograms
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    ax.hist(null_a, bins=50, alpha=0.7, color="steelblue", edgecolor="white", label="Null distribution")
    ax.axvline(result_a.distance, color="red", linewidth=2, label=f"Observed (d={result_a.distance:.2f})")
    ax.set_title(f"HCA human vs Tabula human\np={p_a:.4f}, obs/null={obs_null_a:.3f}", fontsize=11)
    ax.set_xlabel("Procrustes distance")
    ax.set_ylabel("Count")
    ax.legend(fontsize=9)

    ax = axes[1]
    ax.hist(null_b, bins=50, alpha=0.7, color="darkorange", edgecolor="white", label="Null distribution")
    ax.axvline(result_b.distance, color="red", linewidth=2, label=f"Observed (d={result_b.distance:.2f})")
    ax.set_title(f"Tabula human vs Tabula mouse ({n_types} types)\np={p_b:.4f}, obs/null={obs_null_b:.3f}", fontsize=11)
    ax.set_xlabel("Procrustes distance")
    ax.set_ylabel("Count")
    ax.legend(fontsize=9)

    plt.suptitle("HCA Centroid Comparison Diagnostic", fontsize=13, fontweight="bold")
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "null_distributions.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Plot: null_distributions.png")

    # Plot 2: Per-type residual bar comparison
    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(n_types)
    width = 0.35
    bars1 = ax.bar(x - width / 2, [residual_mags_a[ct] for ct in types_a], width,
                   label="HCA-vs-Tabula (within-human)", color="steelblue", alpha=0.8)
    bars2 = ax.bar(x + width / 2, [residual_mags_b[ct] for ct in types_b], width,
                   label="Tabula H-vs-M (cross-species)", color="darkorange", alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([ct.replace("CD4-positive, alpha-beta ", "CD4+ ").replace(
        "CD8-positive, alpha-beta ", "CD8+ ") for ct in types_a], rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("Residual magnitude (PCA space)")
    ax.set_title(f"Per-Cell-Type Residuals: Within-Human vs Cross-Species ({n_types} types)")
    ax.legend()
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "residual_comparison.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Plot: residual_comparison.png")

    # Summary table as CSV
    summary_df = pd.DataFrame({
        "Comparison": ["HCA human vs Tabula human",
                       f"Tabula human vs Tabula mouse (same {n_types} types)"],
        "n_types": [n_types, n_types],
        "p_value": [p_a, p_b],
        "obs_null_ratio": [obs_null_a, obs_null_b],
        "scaling": [result_a.scaling, result_b.scaling],
        "distance": [result_a.distance, result_b.distance],
    })
    summary_df.to_csv(OUTPUT_DIR / "comparison_table.csv", index=False)
    print(f"  Table: comparison_table.csv")

    # ==================================================================
    # Final summary
    # ==================================================================
    print("\n" + "=" * 70)
    print("FINAL RESULT")
    print("=" * 70)
    print(f"""
  ┌────────────────────────────────────────────────────────────────────┐
  │ Comparison                                n  p       obs/null  s  │
  ├────────────────────────────────────────────────────────────────────┤
  │ HCA human vs Tabula human             {n_types:>4}  {p_a:<7.4f} {obs_null_a:<9.3f} {result_a.scaling:.3f}│
  │ Tabula human vs Tabula mouse          {n_types:>4}  {p_b:<7.4f} {obs_null_b:<9.3f} {result_b.scaling:.3f}│
  └────────────────────────────────────────────────────────────────────┘

  VERDICT: {verdict}
  {interpretation}
""")


if __name__ == "__main__":
    main()
