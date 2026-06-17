#!/usr/bin/env python3
"""
ANALYSIS-A: Lineage-Stratified Permutation Test for Procrustes Coherence

Tests whether the cross-species geometric coherence survives a stricter null
hypothesis: within-lineage label permutation. The standard null permutes all
35 cell type labels freely. The lineage-stratified null permutes labels ONLY
within lineage blocks (immune with immune, epithelial with epithelial, etc.),
testing whether positional correspondence is conserved beyond what lineage
membership alone predicts.

Biology: If the global permutation test is significant only because immune types
cluster together and non-immune types cluster together, the within-lineage test
will fail — the signal is lineage-level, not per-type. If the within-lineage
test passes, the geometric coherence reflects per-type positional precision that
goes beyond lineage structure.

Math: For each of 10,000 iterations, within each lineage block of size n_b,
randomly permute which mouse centroid maps to which human centroid (n_b!
configurations per block). Compute Procrustes distance for the block-permuted
configuration. The p-value is the fraction of permuted distances ≤ observed.

Output: output/validation/lineage_stratified/
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT = Path(__file__).resolve().parent.parent
OUTPUT_SCALED = PROJECT / "output" / "phase2" / "scaled_35types"
OUTPUT_DIR = PROJECT / "output" / "validation" / "lineage_stratified"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

SEED = 42
N_PERM = 10_000


# ---------------------------------------------------------------------------
# Lineage block definitions
# ---------------------------------------------------------------------------
# Assigned using Cell Ontology lineage + standard biological grouping.
# Ambiguities documented inline.

LINEAGE_BLOCKS = {
    "immune_hematopoietic": [
        "B cell",
        "CD4-positive, alpha-beta T cell",
        "CD8-positive, alpha-beta T cell",
        "T cell",                         # broad T cell category
        "classical monocyte",
        "granulocyte",
        "hematopoietic precursor cell",   # Note: hematopoietic lineage progenitor
        "hematopoietic stem cell",        # Note: hematopoietic lineage stem cell
        "intermediate monocyte",
        "macrophage",
        "mature NK T cell",
        "monocyte",
        "myeloid dendritic cell",
        "myeloid leukocyte",
        "natural killer cell",
        "neutrophil",
        "non-classical monocyte",
        "plasma cell",
    ],
    "epithelial": [
        "basal cell",
        "bladder urothelial cell",
        "enterocyte of epithelium of large intestine",
        "epithelial cell",                # broad epithelial category
        "large intestine goblet cell",
        "luminal epithelial cell of mammary gland",
        "pancreatic acinar cell",
        "pancreatic ductal cell",
    ],
    "stromal_mesenchymal": [
        "adventitial cell",
        "fibroblast",
        "fibroblast of cardiac tissue",
        "mesenchymal stem cell",          # Note: mesenchymal lineage, grouped with stromal
        "mesenchymal stem cell of adipose tissue",  # Same rationale
        "smooth muscle cell",
        "stromal cell",
    ],
    "endothelial": [
        "endothelial cell",               # Singleton — cannot be permuted within block
    ],
    "metabolic_parenchymal": [
        "hepatocyte",                     # Singleton — cannot be permuted within block
    ],
}

# Document ambiguities
AMBIGUITIES = {
    "hematopoietic precursor cell": (
        "Could be classified as stem/progenitor. Assigned to immune/hematopoietic "
        "because Cell Ontology places it under hematopoietic cell (CL:0000988) and "
        "it differentiates exclusively into blood/immune lineages."
    ),
    "hematopoietic stem cell": (
        "Could be classified as stem/progenitor. Assigned to immune/hematopoietic "
        "for the same reason as hematopoietic precursor cell."
    ),
    "mesenchymal stem cell": (
        "Could be classified as stem/progenitor. Assigned to stromal/mesenchymal "
        "because its differentiation products (fibroblasts, adipocytes, osteoblasts) "
        "are stromal, and its transcriptomic profile clusters with stromal types."
    ),
    "mesenchymal stem cell of adipose tissue": (
        "Same rationale as mesenchymal stem cell."
    ),
    "endothelial cell": (
        "Singleton block. Only one endothelial type in the 35-type ontology. "
        "Cannot be permuted within block — stays correctly matched in all "
        "iterations, making the null more conservative."
    ),
    "hepatocyte": (
        "Singleton block. Only parenchymal/metabolic type in the ontology. "
        "Same singleton treatment as endothelial cell."
    ),
}


# ---------------------------------------------------------------------------
# Procrustes distance (silent, fast)
# ---------------------------------------------------------------------------
def _procrustes_distance(X: np.ndarray, Y: np.ndarray) -> float:
    """Compute Procrustes distance (same math as src/procrustes.py)."""
    X_c = X - X.mean(axis=0)
    Y_c = Y - Y.mean(axis=0)

    M = X_c.T @ Y_c
    U, sigma, Vt = np.linalg.svd(M)
    V = Vt.T

    d = np.linalg.det(V @ U.T)
    D_diag = np.ones(X.shape[1])
    D_diag[-1] = np.sign(d)

    ss_Y = np.sum(Y_c ** 2)
    trace_sigma_D = np.sum(sigma * D_diag)
    s = trace_sigma_D / ss_Y

    Y_aligned = s * (Y_c @ (V * D_diag) @ U.T)
    return np.sqrt(np.sum((X_c - Y_aligned) ** 2))


# ---------------------------------------------------------------------------
# Stratified permutation
# ---------------------------------------------------------------------------
def stratified_permutation(
    X: np.ndarray,
    Y: np.ndarray,
    block_indices: list[list[int]],
    n_perm: int = N_PERM,
    seed: int = SEED,
) -> tuple[float, np.ndarray]:
    """
    Within-lineage permutation test.

    For each iteration, permute cell type labels only within each lineage
    block. Blocks of size 1 are never permuted (identity permutation).

    Args:
        X: Human centroids in PCA space (n × k).
        Y: Mouse centroids in PCA space (n × k).
        block_indices: List of lists, each containing row indices belonging
                       to one lineage block.
        n_perm: Number of permutations.
        seed: Random seed.

    Returns:
        (p_value, null_distribution)
    """
    n = X.shape[0]
    rng = np.random.RandomState(seed)

    observed = _procrustes_distance(X, Y)

    null_distances = np.zeros(n_perm)
    for i in range(n_perm):
        # Build permuted index: within each block, permute independently
        perm = np.arange(n)
        for block in block_indices:
            if len(block) > 1:
                perm[block] = rng.permutation(block)
        null_distances[i] = _procrustes_distance(X, Y[perm])

    n_leq = int(np.sum(null_distances <= observed))
    p_value = (n_leq + 1) / (n_perm + 1)

    return p_value, null_distances


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("ANALYSIS-A: Lineage-Stratified Permutation Test")
    print("=" * 70)

    # --- Load PCA centroids ---
    saved = np.load(OUTPUT_SCALED / "pca_centroids_35.npz", allow_pickle=True)
    cell_types = list(saved["cell_types"])
    X = saved["human"]  # (35, 33)
    Y = saved["mouse"]  # (35, 33)
    n_types = len(cell_types)

    print(f"\nLoaded {n_types} cell types, {X.shape[1]} PCA components")

    # --- Validate lineage assignments ---
    all_assigned = []
    for block_name, members in LINEAGE_BLOCKS.items():
        all_assigned.extend(members)

    assert sorted(all_assigned) == sorted(cell_types), (
        f"Lineage block assignments don't match cell type list!\n"
        f"Missing: {set(cell_types) - set(all_assigned)}\n"
        f"Extra: {set(all_assigned) - set(cell_types)}"
    )

    # --- Build block index arrays ---
    ct_to_idx = {ct: i for i, ct in enumerate(cell_types)}
    block_indices = []
    print("\nLineage blocks:")
    for block_name, members in LINEAGE_BLOCKS.items():
        indices = sorted([ct_to_idx[m] for m in members])
        block_indices.append(indices)
        n_perms_block = 1
        for j in range(1, len(indices) + 1):
            n_perms_block *= j
        print(f"  {block_name:<30} {len(members):>2} types  "
              f"({n_perms_block:>20,} possible permutations)")

    # Count total permutations for the stratified null
    import math
    total_stratified = 1
    total_global = math.factorial(n_types)
    for block_name, members in LINEAGE_BLOCKS.items():
        total_stratified *= math.factorial(len(members))
    print(f"\n  Total stratified permutations: {total_stratified:,.0f}")
    print(f"  Total global permutations:     {total_global:,.0f}")
    print(f"  Ratio (stratified/global):     {total_stratified/total_global:.2e}")

    # --- Document ambiguities ---
    print("\nAmbiguous assignments:")
    for ct, reason in AMBIGUITIES.items():
        print(f"  {ct}: {reason}")

    # --- Observed Procrustes distance ---
    obs_dist = _procrustes_distance(X, Y)
    print(f"\nObserved Procrustes distance: {obs_dist:.6f}")

    # --- Run stratified permutation test ---
    print(f"\nRunning within-lineage permutation test ({N_PERM:,} iterations)...")
    t0 = time.time()
    p_strat, null_strat = stratified_permutation(X, Y, block_indices, N_PERM, SEED)
    elapsed_strat = time.time() - t0
    print(f"  Runtime: {elapsed_strat:.1f}s")

    obs_null_strat = obs_dist / np.median(null_strat)
    print(f"\n  WITHIN-LINEAGE NULL:")
    print(f"    Observed distance:  {obs_dist:.6f}")
    print(f"    Null mean:          {np.mean(null_strat):.6f}")
    print(f"    Null median:        {np.median(null_strat):.6f}")
    print(f"    Null std:           {np.std(null_strat):.6f}")
    print(f"    Null 2.5th-97.5th:  [{np.percentile(null_strat, 2.5):.6f}, "
          f"{np.percentile(null_strat, 97.5):.6f}]")
    print(f"    obs/null ratio:     {obs_null_strat:.4f}")
    print(f"    p-value:            {p_strat:.6f}")
    print(f"    Significant (α=0.05): {'YES' if p_strat < 0.05 else 'NO'}")
    print(f"    Significant (α=0.01): {'YES' if p_strat < 0.01 else 'NO'}")

    # --- Load global null for comparison ---
    null_global = np.load(
        OUTPUT_SCALED / "null_distribution_35.npy"
    )
    obs_null_global = obs_dist / np.median(null_global)

    print(f"\n  GLOBAL (UNRESTRICTED) NULL (loaded from primary analysis):")
    print(f"    Null mean:          {np.mean(null_global):.6f}")
    print(f"    Null std:           {np.std(null_global):.6f}")
    print(f"    obs/null ratio:     {obs_null_global:.4f}")

    # --- Compare distributions ---
    ks_stat, ks_p = stats.ks_2samp(null_strat, null_global)
    mean_diff = np.mean(null_strat) - np.mean(null_global)
    std_ratio = np.std(null_strat) / np.std(null_global)

    print(f"\n  DISTRIBUTION COMPARISON (stratified vs global):")
    print(f"    Stratified null mean:  {np.mean(null_strat):.6f}")
    print(f"    Global null mean:      {np.mean(null_global):.6f}")
    print(f"    Mean difference:       {mean_diff:+.6f} "
          f"({'tighter' if mean_diff < 0 else 'wider'})")
    print(f"    Stratified null std:   {np.std(null_strat):.6f}")
    print(f"    Global null std:       {np.std(null_global):.6f}")
    print(f"    Std ratio (strat/glob): {std_ratio:.4f}")
    print(f"    KS test:               D={ks_stat:.4f}, p={ks_p:.2e}")

    # How much tighter is the within-lineage null
    pct_tighter = (1 - np.mean(null_strat) / np.mean(null_global)) * 100
    print(f"    Within-lineage null is {abs(pct_tighter):.1f}% "
          f"{'lower' if pct_tighter > 0 else 'higher'} than global null")

    # --- Summary ---
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Global permutation:        obs/null = {obs_null_global:.4f}, p = 0.0001")
    print(f"  Within-lineage permutation: obs/null = {obs_null_strat:.4f}, "
          f"p = {p_strat:.4f}")
    if p_strat < 0.05:
        print("  → Geometric coherence SURVIVES the within-lineage test.")
        print("    The signal is per-type, not merely lineage-level.")
    else:
        print("  → Geometric coherence does NOT survive the within-lineage test.")
        print("    The signal is lineage-level, not per-type.")

    # --- Save results ---
    results = {
        "test": "lineage_stratified_permutation",
        "n_cell_types": n_types,
        "n_pca_components": int(X.shape[1]),
        "n_permutations": N_PERM,
        "random_seed": SEED,
        "observed_procrustes_distance": float(obs_dist),
        "lineage_blocks": {
            name: {
                "members": members,
                "n_types": len(members),
                "indices": sorted([ct_to_idx[m] for m in members]),
            }
            for name, members in LINEAGE_BLOCKS.items()
        },
        "ambiguities": AMBIGUITIES,
        "stratified_null": {
            "p_value": float(p_strat),
            "obs_null_ratio": float(obs_null_strat),
            "null_mean": float(np.mean(null_strat)),
            "null_median": float(np.median(null_strat)),
            "null_std": float(np.std(null_strat)),
            "null_min": float(np.min(null_strat)),
            "null_max": float(np.max(null_strat)),
            "null_percentile_2_5": float(np.percentile(null_strat, 2.5)),
            "null_percentile_97_5": float(np.percentile(null_strat, 97.5)),
            "significant_at_005": bool(p_strat < 0.05),
            "significant_at_001": bool(p_strat < 0.01),
        },
        "global_null": {
            "obs_null_ratio": float(obs_null_global),
            "null_mean": float(np.mean(null_global)),
            "null_std": float(np.std(null_global)),
        },
        "distribution_comparison": {
            "mean_difference": float(mean_diff),
            "std_ratio_strat_over_global": float(std_ratio),
            "pct_tighter": float(pct_tighter),
            "ks_statistic": float(ks_stat),
            "ks_p_value": float(ks_p),
        },
    }

    results_path = OUTPUT_DIR / "lineage_stratified_results.json"
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved: {results_path}")

    # Save null distribution
    np.save(OUTPUT_DIR / "null_distribution_stratified.npy", null_strat)
    print(f"Null distribution saved: {OUTPUT_DIR / 'null_distribution_stratified.npy'}")

    # Save lineage assignments as CSV for documentation
    rows = []
    for block_name, members in LINEAGE_BLOCKS.items():
        for ct in members:
            rows.append({
                "cell_type": ct,
                "lineage_block": block_name,
                "ambiguous": ct in AMBIGUITIES,
                "ambiguity_note": AMBIGUITIES.get(ct, ""),
            })
    pd.DataFrame(rows).sort_values(
        ["lineage_block", "cell_type"]
    ).to_csv(OUTPUT_DIR / "lineage_assignments.csv", index=False)
    print(f"Lineage assignments saved: {OUTPUT_DIR / 'lineage_assignments.csv'}")


if __name__ == "__main__":
    main()
