#!/usr/bin/env python3
"""
Parent-child CL subset sensitivity.

Identifies parent-child cell-type pairs in Table S5 per the Cell Ontology
(CL) hierarchy. Constructs two sensitivity variants:
  Variant A — drop parents (broad terms), keep children + independents
  Variant B — drop children (specific terms), keep parents + independents

For each variant, refits the joint PCA at 95% variance on the subset
centroids, runs full Procrustes alignment with 10,000-iteration label
permutation, computes per-type residuals, and reports the Spearman ρ
of the resulting rigidity ranking against the corresponding subset of
the primary 35-type ranking.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA

PROJECT = Path(__file__).resolve().parents[3]
OUT_DIR = Path(__file__).resolve().parent
PHASE2_DIR = PROJECT / "output" / "phase2" / "scaled_35types"

SEED = 42
N_PERM = 10_000
PCA_VARIANCE = 0.95

# Parent-child pairs identified from Table S5 (CL IDs verified against the file)
# Parents map to the list of their child cell types (CL IDs in comments)
PARENT_CHILD = {
    # T cell (CL:0000084) — broad term, with children:
    "T cell": [
        "CD4-positive, alpha-beta T cell",  # CL:0000624
        "CD8-positive, alpha-beta T cell",  # CL:0000625
        "mature NK T cell",                  # CL:0000814
    ],
    # monocyte (CL:0000576) — broad term, with children:
    "monocyte": [
        "classical monocyte",      # CL:0000860
        "intermediate monocyte",   # CL:0002393
        "non-classical monocyte",  # CL:0000875
    ],
    # mesenchymal stem cell (CL:0000134) — with subtype:
    "mesenchymal stem cell": [
        "mesenchymal stem cell of adipose tissue",  # CL:0002570
    ],
    # fibroblast (CL:0000057) — with subtype:
    "fibroblast": [
        "fibroblast of cardiac tissue",  # CL:0002548
    ],
}

# myeloid leukocyte (CL:0000766) ↔ macrophage (CL:0000235): flagged as
# broad/broader relationship rather than strict parent/child. macrophage
# is more specific, myeloid leukocyte is the broader umbrella that
# contains macrophages plus other myeloid types (neutrophil, monocyte,
# myeloid dendritic cell, ...).
BROAD_BROADER = {
    "broader_term": "myeloid leukocyte",
    "more_specific_term": "macrophage",
}


def procrustes_distance(X, Y) -> float:
    X_c = X - X.mean(axis=0)
    Y_c = Y - Y.mean(axis=0)
    M = X_c.T @ Y_c
    U, sigma, Vt = np.linalg.svd(M)
    V = Vt.T
    sign = 1.0 if np.linalg.det(V @ U.T) > 0 else -1.0
    k = X_c.shape[1]
    D = np.ones(k); D[-1] = sign
    ss_Y = float(np.sum(Y_c ** 2))
    s = float(np.sum(sigma * D) / ss_Y) if ss_Y > 0 else 0.0
    Y_aligned = s * (Y_c @ (V * D) @ U.T)
    return float(np.sqrt(np.sum((X_c - Y_aligned) ** 2)))


def procrustes_full(X, Y):
    """Return (distance, rotation R, scaling s, Y_aligned)."""
    X_c = X - X.mean(axis=0)
    Y_c = Y - Y.mean(axis=0)
    M = X_c.T @ Y_c
    U, sigma, Vt = np.linalg.svd(M)
    V = Vt.T
    k = X_c.shape[1]
    D = np.eye(k)
    if np.linalg.det(V @ U.T) < 0:
        D[-1, -1] = -1.0
    R = V @ D @ U.T
    ss_Y = float(np.sum(Y_c ** 2))
    s = float(np.sum(sigma * np.diag(D)) / ss_Y) if ss_Y > 0 else 0.0
    Y_aligned = s * (Y_c @ R)
    d = float(np.sqrt(np.sum((X_c - Y_aligned) ** 2)))
    return d, R, s, Y_aligned, X_c


def run_variant(variant_name: str, types_to_keep: list[str],
                hc: pd.DataFrame, mc: pd.DataFrame,
                primary_residuals: pd.DataFrame) -> dict:
    print(f"\n{'=' * 70}")
    print(f"VARIANT {variant_name}: {len(types_to_keep)} types")
    print('=' * 70)

    cell_types = sorted(types_to_keep)
    hc_s = hc.loc[cell_types]
    mc_s = mc.loc[cell_types]

    # Joint PCA at 95% variance
    combined = np.vstack([hc_s.values, mc_s.values])
    pca = PCA(n_components=PCA_VARIANCE, svd_solver="full", random_state=SEED)
    combined_pca = pca.fit_transform(combined)
    n_pc = pca.n_components_
    H_pca = combined_pca[: len(cell_types)]
    M_pca = combined_pca[len(cell_types) :]
    print(f"  PCA: {n_pc} components, "
          f"{np.sum(pca.explained_variance_ratio_)*100:.1f}% variance")

    # Procrustes + rotation
    obs_d, R, scaling, M_aligned, H_c = procrustes_full(H_pca, M_pca)
    print(f"  Observed Procrustes distance: {obs_d:.4f}, scaling={scaling:.4f}")

    # Permutation test
    rng = np.random.RandomState(SEED)
    null = np.zeros(N_PERM)
    n_types = len(cell_types)
    t0 = time.time()
    for i in range(N_PERM):
        perm = rng.permutation(n_types)
        null[i] = procrustes_distance(H_pca, M_pca[perm])
    elapsed = time.time() - t0
    null_mean = float(null.mean())
    obs_null_ratio = obs_d / null_mean
    p_value = float((np.sum(null <= obs_d) + 1) / (N_PERM + 1))
    print(f"  Null mean: {null_mean:.4f},  obs/null = {obs_null_ratio:.4f},  "
          f"p = {p_value:.2e}  ({elapsed:.0f}s)")

    # Per-type residuals
    residual_mags = []
    for i, ct in enumerate(cell_types):
        residual_mags.append({
            "cell_type": ct,
            "residual_magnitude": float(np.linalg.norm(M_aligned[i] - H_c[i])),
        })
    res_df = pd.DataFrame(residual_mags)
    res_df = res_df.sort_values("residual_magnitude").reset_index(drop=True)
    res_df["rank"] = res_df.index + 1

    # Spearman vs primary ranking (subset)
    primary_subset = primary_residuals[primary_residuals["cell_type"].isin(cell_types)].copy()
    primary_subset = primary_subset.sort_values("residual_magnitude").reset_index(drop=True)
    primary_subset["primary_rank_subset"] = primary_subset.index + 1
    merged = res_df.merge(primary_subset[["cell_type", "primary_rank_subset"]],
                          on="cell_type", how="inner")
    rho, p_rho = stats.spearmanr(merged["rank"], merged["primary_rank_subset"])

    print(f"  Ranking Spearman vs primary subset ({len(merged)} types): "
          f"ρ = {rho:+.4f}, p = {p_rho:.4f}")

    return {
        "variant": variant_name,
        "n_types": n_types,
        "cell_types": cell_types,
        "n_pca_components": int(n_pc),
        "cumulative_variance": float(np.sum(pca.explained_variance_ratio_)),
        "obs_distance": obs_d,
        "null_mean": null_mean,
        "obs_null_ratio": float(obs_null_ratio),
        "p_value_permutation": p_value,
        "scaling": scaling,
        "residuals": res_df.to_dict("records"),
        "ranking_spearman_vs_primary_subset": {
            "rho": float(rho),
            "p_value": float(p_rho),
            "n_overlap": int(len(merged)),
        },
        "runtime_sec": elapsed,
    }


def main() -> None:
    print("=" * 70)
    print("Parent-child CL subset sensitivity")
    print("=" * 70)

    hc = pd.read_csv(PHASE2_DIR / "centroids_human_35.csv", index_col=0)
    mc = pd.read_csv(PHASE2_DIR / "centroids_mouse_35.csv", index_col=0)
    all_types = sorted(hc.index.tolist())
    print(f"All 35 types: {len(all_types)}")

    primary_residuals = pd.read_csv(PHASE2_DIR / "residuals_ranked.csv")

    # Verify parent and child names are present
    missing_parents = [p for p in PARENT_CHILD if p not in all_types]
    missing_children = [c for cs in PARENT_CHILD.values() for c in cs
                        if c not in all_types]
    assert not missing_parents, f"Parents missing: {missing_parents}"
    assert not missing_children, f"Children missing: {missing_children}"
    assert BROAD_BROADER["broader_term"] in all_types
    assert BROAD_BROADER["more_specific_term"] in all_types

    parents = set(PARENT_CHILD.keys())
    children = set(c for cs in PARENT_CHILD.values() for c in cs)

    # Variant A — drop parents (and the broader of macrophage/myeloid leukocyte)
    drop_A = parents | {BROAD_BROADER["broader_term"]}
    variant_A_types = [t for t in all_types if t not in drop_A]

    # Variant B — drop children (and the more specific of macrophage/myeloid leukocyte)
    drop_B = children | {BROAD_BROADER["more_specific_term"]}
    variant_B_types = [t for t in all_types if t not in drop_B]

    print(f"\nDropped in Variant A (parents + broader): {sorted(drop_A)}")
    print(f"Dropped in Variant B (children + more specific): {sorted(drop_B)}")
    print(f"Variant A kept: {len(variant_A_types)} types")
    print(f"Variant B kept: {len(variant_B_types)} types")

    result_A = run_variant("A_drop_parents", variant_A_types, hc, mc, primary_residuals)
    result_B = run_variant("B_drop_children", variant_B_types, hc, mc, primary_residuals)

    summary = {
        "metadata": {
            "script": str(Path(__file__).resolve().relative_to(PROJECT)),
            "seed": SEED,
            "n_permutations": N_PERM,
            "pca_variance_threshold": PCA_VARIANCE,
            "primary_reference": {
                "n_types": 35,
                "obs_null_ratio": 0.522,
                "p_value": "<1e-6",
            },
        },
        "parent_child_pairs": PARENT_CHILD,
        "broad_broader_pair": BROAD_BROADER,
        "variant_A_drop_parents": result_A,
        "variant_B_drop_children": result_B,
        "comparison": {
            "primary_35_obs_null": 0.522,
            "variant_A_obs_null": result_A["obs_null_ratio"],
            "variant_A_p_perm": result_A["p_value_permutation"],
            "variant_A_ranking_spearman": result_A["ranking_spearman_vs_primary_subset"]["rho"],
            "variant_B_obs_null": result_B["obs_null_ratio"],
            "variant_B_p_perm": result_B["p_value_permutation"],
            "variant_B_ranking_spearman": result_B["ranking_spearman_vs_primary_subset"]["rho"],
        },
    }

    out_path = OUT_DIR / "results.json"
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\n✓ Saved {out_path}")

    # Compact CSV summary
    rows = []
    for v, r in [("A_drop_parents", result_A), ("B_drop_children", result_B)]:
        rows.append({
            "variant": v,
            "n_types": r["n_types"],
            "obs_null_ratio": round(r["obs_null_ratio"], 4),
            "p_perm": r["p_value_permutation"],
            "ranking_spearman_vs_primary": round(
                r["ranking_spearman_vs_primary_subset"]["rho"], 4),
            "spearman_p": r["ranking_spearman_vs_primary_subset"]["p_value"],
        })
    pd.DataFrame(rows).to_csv(OUT_DIR / "summary.csv", index=False)
    print(f"✓ Saved {OUT_DIR/'summary.csv'}")


if __name__ == "__main__":
    main()
