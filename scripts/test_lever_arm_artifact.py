# DIAGNOSTIC ONLY — not part of canonical validation pipeline
#
# Test whether Procrustes residual magnitude correlates with distance from
# the centroid-of-centroids in PCA space ("lever-arm artifact" check).
#
# Adversarial claim: peripheral cell types in PCA space mechanically produce
# larger Procrustes residuals because they sit further from the rotation
# centre, amplifying any rotational misfit. If true, the rigidity ranking
# partly reflects geometry of the PCA embedding, not biology.
#
# Method:
#   1. Load 35-type PCA centroids and re-run Procrustes alignment.
#   2. Compute the centroid-of-centroids (global mean of the 35 HUMAN
#      centered positions — the Procrustes reference frame).
#   3. For each cell type, compute Euclidean distance of the human centroid
#      from the global mean (human frame chosen because Procrustes aligns
#      mouse ONTO human; the human positions define the reference geometry).
#   4. Spearman-correlate distance with residual magnitude.
#
# Decision rule:
#   |ρ| < 0.3 AND p > 0.05 → no lever-arm artifact
#   |ρ| ≥ 0.3 OR  p ≤ 0.05 → confound detected, escalate to ADVISOR

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from cellwarp.procrustes import procrustes_align, compute_residual_vectors  # noqa: E402

PCA_PATH = PROJECT_ROOT / "output" / "phase2" / "scaled_35types" / "pca_centroids_35.npz"
RESIDUALS_PATH = PROJECT_ROOT / "output" / "phase2" / "scaled_35types" / "residuals_ranked.csv"
OUTPUT_DIR = PROJECT_ROOT / "output" / "diagnostics" / "lever_arm"


def main():
    print("=" * 70)
    print("LEVER-ARM ARTIFACT DIAGNOSTIC")
    print("=" * 70)

    # ------------------------------------------------------------------
    # 1. Load PCA centroids (pre-Procrustes)
    # ------------------------------------------------------------------
    data = np.load(PCA_PATH, allow_pickle=True)
    human_pca = data["human"]          # (35, 33)
    mouse_pca = data["mouse"]          # (35, 33)
    cell_types = data["cell_types"].tolist()

    n_types, n_dims = human_pca.shape
    print(f"\nLoaded {n_types} cell types × {n_dims} PCA dimensions")

    # ------------------------------------------------------------------
    # 2. Re-run Procrustes alignment (mouse → human)
    # ------------------------------------------------------------------
    result = procrustes_align(human_pca, mouse_pca)

    # Compute residuals (aligned_mouse - centered_human)
    residuals = compute_residual_vectors(result, cell_types)
    residual_mags = {ct: np.linalg.norm(residuals[ct]) for ct in cell_types}

    # ------------------------------------------------------------------
    # 3. Compute distance from centroid-of-centroids (human reference frame)
    # ------------------------------------------------------------------
    # Procrustes centers the reference (human) by subtracting its mean.
    # After centering, the centroid-of-centroids is the origin (0 vector).
    # Distance from origin = norm of each centered human position.
    #
    # Choice: use human (reference) positions because Procrustes rotation
    # is fitted to align mouse onto human. The human frame defines the
    # geometry that the rotation operates on.
    human_centered = result.centered_reference  # (35, 33), mean-subtracted
    dist_from_origin = np.array([
        np.linalg.norm(human_centered[i]) for i in range(n_types)
    ])

    # ------------------------------------------------------------------
    # 4. Build comparison table and correlate
    # ------------------------------------------------------------------
    df = pd.DataFrame({
        "cell_type": cell_types,
        "residual_magnitude": [residual_mags[ct] for ct in cell_types],
        "distance_from_centroid_of_centroids": dist_from_origin,
    }).sort_values("residual_magnitude", ascending=False).reset_index(drop=True)

    rho, p_value = stats.spearmanr(
        df["distance_from_centroid_of_centroids"],
        df["residual_magnitude"],
    )

    # ------------------------------------------------------------------
    # 5. Also check mouse-side distances as sensitivity
    # ------------------------------------------------------------------
    mouse_aligned = result.aligned_target  # (35, 33)
    mouse_dist = np.array([
        np.linalg.norm(mouse_aligned[i]) for i in range(n_types)
    ])
    rho_mouse, p_mouse = stats.spearmanr(mouse_dist, df["residual_magnitude"].values)

    # Mean of human + mouse distances
    mean_dist = (dist_from_origin + mouse_dist) / 2.0
    rho_mean, p_mean = stats.spearmanr(mean_dist, df["residual_magnitude"].values)

    # ------------------------------------------------------------------
    # 6. Report
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("RESULTS")
    print("=" * 70)

    print(f"\n  Distance metric: Euclidean distance from centroid-of-centroids")
    print(f"  Reference frame: human (Procrustes reference), centered (origin = global mean)")
    print(f"  n = {n_types} cell types")

    print(f"\n  --- Primary (human reference distances) ---")
    print(f"  Spearman ρ = {rho:.3f}")
    print(f"  p-value   = {p_value:.4f}")

    print(f"\n  --- Sensitivity: mouse (aligned) distances ---")
    print(f"  Spearman ρ = {rho_mouse:.3f}")
    print(f"  p-value   = {p_mouse:.4f}")

    print(f"\n  --- Sensitivity: mean(human, mouse) distances ---")
    print(f"  Spearman ρ = {rho_mean:.3f}")
    print(f"  p-value   = {p_mean:.4f}")

    # Decision
    print("\n" + "=" * 70)
    if abs(rho) < 0.3 and p_value > 0.05:
        verdict = "PASS"
        interpretation = (
            "No evidence for lever-arm artifact. Distance from the "
            "centroid-of-centroids does not predict Procrustes residual "
            f"magnitude (Spearman ρ={rho:.3f}, p={p_value:.3f}, n={n_types}). "
            "No text change needed."
        )
    else:
        verdict = "FLAG"
        interpretation = (
            "Lever-arm confound detected. Distance from centroid-of-centroids "
            f"correlates with residual magnitude (Spearman ρ={rho:.3f}, "
            f"p={p_value:.3f}, n={n_types}). ADVISOR must decide whether to "
            "partial out distance from treeness result or add limitation sentence."
        )

    print(f"  VERDICT: {verdict}")
    print(f"  {interpretation}")
    print("=" * 70)

    # ------------------------------------------------------------------
    # 7. Save outputs
    # ------------------------------------------------------------------
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    df["distance_mouse_aligned"] = mouse_dist
    df["distance_mean"] = mean_dist
    df.to_csv(OUTPUT_DIR / "lever_arm_distances.csv", index=False)

    import json
    results = {
        "test": "lever_arm_artifact",
        "n_types": n_types,
        "n_pca_dims": n_dims,
        "distance_frame": "human_centered_reference",
        "primary": {
            "rho": float(rho),
            "p_value": float(p_value),
        },
        "sensitivity_mouse_aligned": {
            "rho": float(rho_mouse),
            "p_value": float(p_mouse),
        },
        "sensitivity_mean_distance": {
            "rho": float(rho_mean),
            "p_value": float(p_mean),
        },
        "verdict": verdict,
        "interpretation": interpretation,
    }
    with open(OUTPUT_DIR / "lever_arm_results.json", "w") as f:
        json.dump(results, f, indent=2)

    print(f"\n  Saved: {OUTPUT_DIR / 'lever_arm_distances.csv'}")
    print(f"  Saved: {OUTPUT_DIR / 'lever_arm_results.json'}")


if __name__ == "__main__":
    main()
