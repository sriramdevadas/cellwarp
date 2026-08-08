#!/usr/bin/env python3
"""
CellWarp — Leave-One-Out Cross-Validation on Cell Types (Phase 3)

For each of 35 cell types, trains Procrustes on the remaining 34 and predicts
the held-out type's position. Tests whether the learned geometric transformation
generalizes beyond the training set.

Biology
-------
If the Procrustes transformation captures a genuine, coherent geometric
relationship between species' cell type landscapes, then a transformation
learned from 34 cell types should predict where the 35th cell type will fall
in the other species. This is the strongest test that the geometric signal
is not overfitting — it demonstrates predictive structure.

Math
----
For each held-out cell type i ∈ {1, ..., 35}:
  1. TRAIN: Remove type i from both species. PCA on 68 remaining centroids.
     Procrustes on 34 training types → get R, s, μ_X, μ_Y.
  2. PREDICT: Project held-out mouse centroid into training PCA space.
     Apply learned transformation: ŷ_i = s·(y_i − μ_Y)·R + μ_X
  3. ERROR: e_i = ‖ŷ_i − x_i‖  (Euclidean distance in PCA space to actual
     held-out human centroid, also projected into training PCA space).
  4. NULL: Mean distance from ŷ_i to all 34 OTHER human centroids in PCA
     space. This represents the expected error if the held-out mouse centroid
     were randomly paired with a different human cell type.
  5. RATIO: r_i = e_i / null_i.  r_i < 1.0 means the model predicts the
     correct position better than random pairing.

Gate criterion: >50% of cell types have ratio < 1.0 (predicted better than chance).

Inputs:
    output/phase2/scaled_35types/centroids_human_35.csv
    output/phase2/scaled_35types/centroids_mouse_35.csv

Outputs (all in output/phase3/loocv/):
    loocv_results.csv          — per-cell-type error, null, ratio
    loocv_summary.json         — gate evaluation
    loocv_barchart.png         — error/null ratio per cell type

Usage:
    python scripts/08_loocv.py
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

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from cellwarp.procrustes import (
    PCA_VARIANCE_THRESHOLD,
    procrustes_align,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CENTROID_DIR = Path("output/phase2/scaled_35types")
OUTPUT_DIR = Path("output/phase3/loocv")
LOOCV_GATE_FRACTION = 0.50  # >50% of types must have ratio < 1.0

# Cell category assignments for coloring the bar chart.
# Categories: Immune, Stem/Progenitor, Structural, Epithelial.
CELL_CATEGORIES = {
    "B cell": "Immune",
    "CD4-positive, alpha-beta T cell": "Immune",
    "CD8-positive, alpha-beta T cell": "Immune",
    "T cell": "Immune",
    "classical monocyte": "Immune",
    "granulocyte": "Immune",
    "intermediate monocyte": "Immune",
    "macrophage": "Immune",
    "mature NK T cell": "Immune",
    "monocyte": "Immune",
    "myeloid dendritic cell": "Immune",
    "myeloid leukocyte": "Immune",
    "natural killer cell": "Immune",
    "neutrophil": "Immune",
    "non-classical monocyte": "Immune",
    "plasma cell": "Immune",
    "hematopoietic precursor cell": "Stem/Progenitor",
    "hematopoietic stem cell": "Stem/Progenitor",
    "mesenchymal stem cell": "Stem/Progenitor",
    "mesenchymal stem cell of adipose tissue": "Stem/Progenitor",
    "basal cell": "Stem/Progenitor",
    "endothelial cell": "Structural",
    "fibroblast": "Structural",
    "fibroblast of cardiac tissue": "Structural",
    "adventitial cell": "Structural",
    "smooth muscle cell": "Structural",
    "stromal cell": "Structural",
    "hepatocyte": "Epithelial",
    "enterocyte of epithelium of large intestine": "Epithelial",
    "large intestine goblet cell": "Epithelial",
    "pancreatic acinar cell": "Epithelial",
    "pancreatic ductal cell": "Epithelial",
    "bladder urothelial cell": "Epithelial",
    "luminal epithelial cell of mammary gland": "Epithelial",
    "epithelial cell": "Epithelial",
}

CATEGORY_COLORS = {
    "Epithelial": "#9467bd",
    "Immune": "#3574b0",
    "Stem/Progenitor": "#e2783c",
    "Structural": "#4ea34e",
}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    t_start = time.time()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("PHASE 3 — LEAVE-ONE-OUT CROSS-VALIDATION (35 CELL TYPES)")
    print("=" * 70)

    # ---- Load pre-computed centroids ----
    print("\n  Loading pre-computed centroids...")
    h_cent = pd.read_csv(CENTROID_DIR / "centroids_human_35.csv", index_col=0)
    m_cent = pd.read_csv(CENTROID_DIR / "centroids_mouse_35.csv", index_col=0)

    cell_types = sorted(h_cent.index.tolist())
    n_types = len(cell_types)
    n_genes = h_cent.shape[1]

    print(f"  Human centroids: {n_types} types × {n_genes:,} genes")
    print(f"  Mouse centroids: {m_cent.shape[0]} types × {m_cent.shape[1]:,} genes")
    assert sorted(m_cent.index.tolist()) == cell_types, "Cell type mismatch!"

    # ---- LOOCV loop ----
    print(f"\n  Running LOOCV ({n_types} folds)...")

    results = []

    # Suppress spurious RuntimeWarnings from np.linalg.det inside procrustes_align.
    # (ISSUE-022) Nothing overflows. det() receives V @ U.T from an SVD: a 33x33
    # orthogonal matrix, condition number 1, determinant +/-1 -- measured
    # 1.0000000000000058, some 308 orders below the float64 ceiling. The warnings
    # are FPU status flags raised inside the LAPACK routine and are backend
    # specific: three fire under an Accelerate-backed numpy (divide by zero,
    # overflow, invalid) and none at all under an OpenBLAS-backed one, on the same
    # input with the same result. The value and its sign are correct either way,
    # which is what makes suppressing them safe.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)

        for i, held_out in enumerate(cell_types):
            # 1. Remove held-out type from training set
            train_types = [ct for ct in cell_types if ct != held_out]
            h_train = h_cent.loc[train_types]
            m_train = m_cent.loc[train_types]

            # 2. PCA on 68 training centroids (34 human + 34 mouse)
            h_mat = h_train.values  # (34, G)
            m_mat = m_train.values  # (34, G)
            combined = np.vstack([h_mat, m_mat])  # (68, G)

            pca = PCA(
                n_components=PCA_VARIANCE_THRESHOLD,
                svd_solver="full",
                random_state=42,
            )
            combined_pca = pca.fit_transform(combined)
            n_train = len(train_types)
            h_pca = combined_pca[:n_train]  # (34, k)
            m_pca = combined_pca[n_train:]  # (34, k)

            # 3. Procrustes on 34 training types (silent)
            with contextlib.redirect_stdout(io.StringIO()):
                proc = procrustes_align(h_pca, m_pca)

            # 4. Project held-out centroids into training PCA space
            held_h_gene = h_cent.loc[held_out].values.reshape(1, -1)
            held_m_gene = m_cent.loc[held_out].values.reshape(1, -1)
            held_h_pca = pca.transform(held_h_gene).flatten()  # (k,)
            held_m_pca = pca.transform(held_m_gene).flatten()  # (k,)

            # 5. Apply Procrustes transform to held-out mouse centroid
            # Transform: predicted = s * (y - μ_Y) @ R + μ_X
            held_m_centered = held_m_pca - proc.translation_target
            predicted = (
                proc.scaling * (held_m_centered @ proc.rotation)
                + proc.translation_ref
            )

            # 6. Error = Euclidean distance from predicted to actual human
            error = float(np.linalg.norm(predicted - held_h_pca))

            # 7. Null = mean distance from predicted position to all 34
            #    training human centroids (world coords = h_pca directly)
            null_distances = np.array([
                np.linalg.norm(predicted - h_pca[j])
                for j in range(n_train)
            ])
            null_dist = float(np.mean(null_distances))

            ratio = error / null_dist if null_dist > 0 else float("inf")
            category = CELL_CATEGORIES.get(held_out, "Other")

            results.append({
                "cell_type": held_out,
                "category": category,
                "loocv_error": error,
                "null_distance": null_dist,
                "ratio": ratio,
                "n_pca_components": int(pca.n_components_),
                "better_than_chance": ratio < 1.0,
            })

            status = "BETTER" if ratio < 1.0 else "worse"
            print(
                f"    [{i + 1:>2}/{n_types}] {held_out:<45} "
                f"error={error:.4f}  null={null_dist:.4f}  "
                f"ratio={ratio:.4f}  ({status})"
            )

    # ---- Compute summary statistics ----
    df = pd.DataFrame(results)
    n_better = int(df["better_than_chance"].sum())
    frac_better = n_better / n_types
    loocv_pass = frac_better > LOOCV_GATE_FRACTION

    mean_ratio = float(df["ratio"].mean())
    median_ratio = float(df["ratio"].median())
    best_idx = int(df["ratio"].idxmin())
    worst_idx = int(df["ratio"].idxmax())

    # ---- Save results ----
    df.to_csv(OUTPUT_DIR / "loocv_results.csv", index=False)

    summary = {
        "n_cell_types": n_types,
        "n_better_than_chance": n_better,
        "fraction_better_than_chance": float(frac_better),
        "gate_criterion": {
            "threshold": LOOCV_GATE_FRACTION,
            "fraction": float(frac_better),
            "pass": bool(loocv_pass),
        },
        "per_type_summary": {
            "mean_ratio": mean_ratio,
            "median_ratio": median_ratio,
            "min_ratio": float(df["ratio"].min()),
            "max_ratio": float(df["ratio"].max()),
            "best_predicted": df.loc[best_idx, "cell_type"],
            "worst_predicted": df.loc[worst_idx, "cell_type"],
        },
    }

    with open(OUTPUT_DIR / "loocv_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # ---- Plot bar chart ----
    df_sorted = df.sort_values("ratio", ascending=True).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(12, 9))
    ax.barh(
        range(len(df_sorted)),
        df_sorted["ratio"].values,
        color=[
            CATEGORY_COLORS.get(c, "gray") for c in df_sorted["category"]
        ],
        edgecolor="white",
        linewidth=0.5,
    )
    ax.axvline(
        1.0, color="red", linestyle="--", linewidth=2, label="Null (ratio = 1.0)"
    )
    ax.set_yticks(range(len(df_sorted)))
    ax.set_yticklabels(df_sorted["cell_type"].values, fontsize=8)
    ax.set_xlabel("LOOCV Error / Null Distance", fontsize=12)
    ax.set_title(
        f"LOOCV Error Ratio per Cell Type (35 types)\n"
        f"{n_better}/{n_types} ({frac_better * 100:.0f}%) predicted better "
        f"than chance — {'PASS' if loocv_pass else 'FAIL'} "
        f"(gate: >{LOOCV_GATE_FRACTION * 100:.0f}%)",
        fontsize=12,
    )

    # Legend for categories
    legend_patches = [
        Patch(facecolor=CATEGORY_COLORS[cat], label=cat)
        for cat in sorted(CATEGORY_COLORS.keys())
    ]
    legend_patches.append(
        plt.Line2D(
            [0], [0], color="red", linestyle="--", linewidth=2, label="Null (1.0)"
        )
    )
    ax.legend(handles=legend_patches, loc="lower right", fontsize=9)
    ax.invert_yaxis()
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "loocv_barchart.png", dpi=150)
    plt.close(fig)

    # ---- Summary ----
    total_time = time.time() - t_start

    print("\n" + "=" * 70)
    print("LOOCV — RESULTS")
    print("=" * 70)
    print(f"\n  Cell types: {n_types}")
    print(
        f"  Types predicted better than chance: "
        f"{n_better}/{n_types} ({frac_better * 100:.0f}%)"
    )
    print(f"  Mean ratio: {mean_ratio:.4f}")
    print(f"  Median ratio: {median_ratio:.4f}")
    print(
        f"  Best predicted: {df.loc[best_idx, 'cell_type']} "
        f"(ratio={df['ratio'].min():.4f})"
    )
    print(
        f"  Worst predicted: {df.loc[worst_idx, 'cell_type']} "
        f"(ratio={df['ratio'].max():.4f})"
    )
    print(
        f"\n  GATE CRITERION: >{LOOCV_GATE_FRACTION * 100:.0f}% of types "
        f"predicted better than chance"
    )
    print(f"    Fraction: {frac_better * 100:.0f}%")
    print(f"    Result: {'PASS' if loocv_pass else 'FAIL'}")

    print(f"\n  By category:")
    for cat in sorted(CATEGORY_COLORS.keys()):
        cat_df = df[df["category"] == cat]
        if len(cat_df) > 0:
            cat_better = int(cat_df["better_than_chance"].sum())
            print(
                f"    {cat:<20} {cat_better}/{len(cat_df)} better than chance "
                f"(mean ratio={cat_df['ratio'].mean():.4f})"
            )

    print(f"\n  Total runtime: {total_time:.1f}s ({total_time / 60:.1f} min)")
    print(f"  Output: {OUTPUT_DIR}/")

    print(
        f"\n  Plot description: Horizontal bar chart showing LOOCV error/null "
        f"ratio per cell type, sorted by ratio (best predicted at top)."
    )
    print(
        f"  Bars colored by cell category (Immune=blue, Stem/Progenitor=orange, "
        f"Structural=green, Epithelial=purple)."
    )
    print(
        f"  Red dashed line at 1.0 marks the null threshold. "
        f"Types left of line are predicted better than random pairing."
    )


if __name__ == "__main__":
    main()
