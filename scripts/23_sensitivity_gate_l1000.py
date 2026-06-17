"""
CellWarp — Sensitivity Gate: Fractal Geometry Test (L1000 Landmark Genes)

Pre-registration Step 1 from:
  docs/preregistration_dilirank_hepatocyte_2026-03-16.md

Biology
-------
The L1000 platform measures 978 "landmark" genes chosen to capture most
transcriptomic variance via inference. This test asks: does the CellWarp
rigidity ranking — computed over ~17k ortholog genes — survive when
restricted to only these 978 genes (5.8% of gene space)?

If yes, cell identity geometry is "fractal" — preserved at dramatically
reduced dimensionality. This confirms the L1000 platform can capture
the geometric signal, unblocking DILIrank validation.

Math
----
1. Subset human and mouse 35-type centroid matrices to landmark genes
   present in ortholog space.
2. PCA-reduce the subsetted centroids (same pipeline as full analysis).
3. Run Procrustes alignment + permutation test.
4. Extract per-cell-type residual magnitudes → rigidity ranking.
5. Spearman ρ between landmark-space ranking and full-space ranking.

Threshold: ρ ≥ 0.6 → PASS (fractal geometry confirmed).
           ρ < 0.6 → FAIL (L1000 platform limitation).

NO DILI DATA IS ACCESSED IN THIS SCRIPT.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# Add src/ to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from cellwarp.procrustes import pca_reduce_centroids, procrustes_align, permutation_test

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT = Path(__file__).resolve().parent.parent
CENTROIDS_HUMAN = PROJECT / "output/phase2/scaled_35types/centroids_human_35.csv"
CENTROIDS_MOUSE = PROJECT / "output/phase2/scaled_35types/centroids_mouse_35.csv"
ORTHOLOGS = PROJECT / "data/phase1/orthologs_human_mouse.csv"
LANDMARK_GENES = Path("/tmp/l1000_landmark_genes.txt")
FULL_RANKING = PROJECT / "output/phase2/scaled_35types/residuals_ranked.csv"
OUTPUT_DIR = PROJECT / "output/landmark_sensitivity"

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
print("=" * 70)
print("SENSITIVITY GATE: Fractal Geometry Test (L1000 978 Landmark Genes)")
print("=" * 70)

# Load landmark gene symbols
landmark_symbols = set(
    pd.read_csv(LANDMARK_GENES, header=None, names=["gene_symbol"])["gene_symbol"]
)
print(f"\nL1000 landmark genes loaded: {len(landmark_symbols)}")

# Load ortholog mapping (human gene symbol → Ensembl ID)
orthologs = pd.read_csv(ORTHOLOGS)
print(f"Ortholog pairs loaded: {len(orthologs)}")

# Map landmark symbols to human Ensembl IDs present in ortholog space
landmark_orthologs = orthologs[
    orthologs["human_gene_name"].isin(landmark_symbols)
]
landmark_ensembl_ids = set(landmark_orthologs["human_ensembl_id"])
print(f"Landmark genes in ortholog space: {len(landmark_ensembl_ids)} of {len(landmark_symbols)}")
print(f"  Coverage: {len(landmark_ensembl_ids) / len(landmark_symbols) * 100:.1f}%")

# Report which landmark genes are NOT in ortholog space
missing = landmark_symbols - set(landmark_orthologs["human_gene_name"])
if missing:
    print(f"  Missing from ortholog space: {len(missing)} genes")
    if len(missing) <= 20:
        print(f"    {sorted(missing)}")

# Load full-space centroids
print("\nLoading centroid matrices...")
human_full = pd.read_csv(CENTROIDS_HUMAN, index_col=0)
mouse_full = pd.read_csv(CENTROIDS_MOUSE, index_col=0)
print(f"  Human centroids: {human_full.shape}")
print(f"  Mouse centroids: {mouse_full.shape}")

# Verify landmark Ensembl IDs are in centroid columns
available_landmark = landmark_ensembl_ids & set(human_full.columns)
print(f"\nLandmark genes in centroid matrices: {len(available_landmark)}")

# ---------------------------------------------------------------------------
# Subset centroids to landmark genes
# ---------------------------------------------------------------------------
landmark_cols = sorted(available_landmark)
human_lm = human_full[landmark_cols]
mouse_lm = mouse_full[landmark_cols]
print(f"\nSubsetted centroids: {human_lm.shape[0]} types × {human_lm.shape[1]} landmark genes")
print(f"  Gene space reduction: {human_full.shape[1]} → {human_lm.shape[1]} "
      f"({human_lm.shape[1] / human_full.shape[1] * 100:.1f}%)")

# ---------------------------------------------------------------------------
# PCA + Procrustes on landmark gene subset
# ---------------------------------------------------------------------------
print("\n" + "-" * 70)
print("Running PCA + Procrustes on landmark gene subset")
print("-" * 70)

human_pca, mouse_pca, pca_model, cell_types = pca_reduce_centroids(
    human_lm, mouse_lm, variance_threshold=0.95
)

result = procrustes_align(human_pca, mouse_pca)
p_value, null_dist = permutation_test(human_pca, mouse_pca)

# ---------------------------------------------------------------------------
# Extract residual magnitudes → rigidity ranking
# ---------------------------------------------------------------------------
print("\n" + "-" * 70)
print("Per-cell-type residuals in landmark gene space")
print("-" * 70)

residual_data = []
total_ssr = result.distance_squared
for i, ct in enumerate(cell_types):
    r = result.aligned_target[i] - result.centered_reference[i]
    mag = np.linalg.norm(r)
    pct = (mag ** 2 / total_ssr * 100) if total_ssr > 0 else 0.0
    residual_data.append({
        "cell_type": ct,
        "residual_magnitude_landmark": mag,
        "pct_of_ssr_landmark": pct,
    })

lm_ranking = pd.DataFrame(residual_data)
lm_ranking = lm_ranking.sort_values("residual_magnitude_landmark").reset_index(drop=True)
lm_ranking["rank_landmark"] = range(1, len(lm_ranking) + 1)
# Rank 1 = smallest residual = most rigid
# (Matching the convention in residuals_ranked.csv where rank 35 = most rigid)
# Actually, residuals_ranked.csv ranks 1 = largest residual. Let me check.
# From the file: rank 1 = stromal cell (16.98, largest), rank 35 = CD8+ T (5.38, smallest)
# So rank 1 = most diverged, rank 35 = most conserved/rigid.
# For Spearman ρ we just need the residual magnitudes, not ranks.

# ---------------------------------------------------------------------------
# Load full-space ranking and compute Spearman ρ
# ---------------------------------------------------------------------------
print("\n" + "-" * 70)
print("Spearman ρ: landmark-space vs full-space rigidity")
print("-" * 70)

full_ranking = pd.read_csv(FULL_RANKING)

# Merge on cell_type
merged = pd.merge(
    full_ranking[["cell_type", "residual_magnitude"]],
    lm_ranking[["cell_type", "residual_magnitude_landmark"]],
    on="cell_type",
)
print(f"\nCell types matched: {len(merged)} of {len(full_ranking)}")

rho, p_spearman = stats.spearmanr(
    merged["residual_magnitude"],
    merged["residual_magnitude_landmark"],
)

print(f"\n  Spearman ρ = {rho:.4f}")
print(f"  p-value    = {p_spearman:.2e}")
print(f"  n          = {len(merged)}")

# ---------------------------------------------------------------------------
# Threshold check
# ---------------------------------------------------------------------------
THRESHOLD = 0.6

print("\n" + "=" * 70)
if rho >= THRESHOLD:
    print(f"  SENSITIVITY GATE: PASS (ρ = {rho:.4f} ≥ {THRESHOLD})")
    print(f"  FRACTAL GEOMETRY CONFIRMED")
    print(f"  Identity geometry preserved at {human_lm.shape[1] / human_full.shape[1] * 100:.1f}% of gene space")
    print(f"  DILIrank analysis UNBLOCKED")
else:
    print(f"  SENSITIVITY GATE: FAIL (ρ = {rho:.4f} < {THRESHOLD})")
    print(f"  L1000 PLATFORM LIMITATION")
    print(f"  Landmark genes do not preserve identity geometry")
    print(f"  DILIrank analysis ABORTED")
print("=" * 70)

# ---------------------------------------------------------------------------
# Save output
# ---------------------------------------------------------------------------
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Save landmark-space ranking
lm_ranking.to_csv(OUTPUT_DIR / "residuals_ranked_landmark.csv", index=False)

# Save merged comparison
merged["rank_full"] = merged["residual_magnitude"].rank(ascending=False).astype(int)
merged["rank_landmark"] = merged["residual_magnitude_landmark"].rank(ascending=False).astype(int)
merged["rank_diff"] = abs(merged["rank_full"] - merged["rank_landmark"])
merged = merged.sort_values("rank_full")
merged.to_csv(OUTPUT_DIR / "rigidity_ranking_comparison.csv", index=False)

# Save summary
import json
summary = {
    "test": "Sensitivity Gate — Fractal Geometry Test",
    "date": "2026-03-16",
    "preregistration": "docs/preregistration_dilirank_hepatocyte_2026-03-16.md",
    "l1000_landmark_genes_total": int(len(landmark_symbols)),
    "landmark_genes_in_ortholog_space": int(len(landmark_ensembl_ids)),
    "landmark_genes_in_centroids": int(len(available_landmark)),
    "full_space_genes": int(human_full.shape[1]),
    "gene_space_reduction_pct": round(float(len(available_landmark)) / float(human_full.shape[1]) * 100, 1),
    "procrustes_distance_landmark": round(float(result.distance), 6),
    "procrustes_p_value_landmark": round(float(p_value), 6),
    "spearman_rho": round(float(rho), 4),
    "spearman_p_value": float(f"{p_spearman:.2e}"),
    "threshold": float(THRESHOLD),
    "result": "PASS" if rho >= THRESHOLD else "FAIL",
    "n_cell_types": int(len(merged)),
    "pca_components_landmark": int(pca_model.n_components_),
}
with open(OUTPUT_DIR / "sensitivity_gate_result.json", "w") as f:
    json.dump(summary, f, indent=2)

# Print top/bottom of ranking comparison
print(f"\nRanking comparison (top 5 most diverged, bottom 5 most rigid):")
print(f"  {'Cell Type':<45} {'Full Rank':>9} {'LM Rank':>8} {'Δ':>4}")
print(f"  {'-' * 68}")
for _, row in merged.head(5).iterrows():
    print(f"  {row['cell_type']:<45} {int(row['rank_full']):>9} {int(row['rank_landmark']):>8} {int(row['rank_diff']):>4}")
print(f"  {'...'}")
for _, row in merged.tail(5).iterrows():
    print(f"  {row['cell_type']:<45} {int(row['rank_full']):>9} {int(row['rank_landmark']):>8} {int(row['rank_diff']):>4}")

print(f"\nMean absolute rank difference: {merged['rank_diff'].mean():.1f}")
print(f"Max rank difference: {int(merged['rank_diff'].max())}")
print(f"\nOutput saved to: {OUTPUT_DIR}")
print("\nDone.")
