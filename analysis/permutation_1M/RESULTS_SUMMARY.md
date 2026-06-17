# Procrustes Permutation Test — 1,000,000 Permutations

**Motivation:** increase from 10⁴ to 10⁶ permutations for a tighter p-value bound.
**Analysis:** PRIMARY (35 cell types × 33 PCs, human vs mouse)

## Result

| Metric | Value |
|---|---|
| Observed Procrustes distance | 61.153 |
| obs/null ratio | 0.523 |
| Permutations | 1,000,000 |
| Null values ≤ observed | **0** |
| **p-value** | **p < 10⁻⁶** |
| Runtime | 116.4 s (1.94 min) |

The observed distance (61.15) is **26 standard deviations** below the null mean (116.94). The closest any of 1,000,000 null permutations came was 98.88 — still 38 units above observed. The p-value remains floored even at 10⁶ resolution.

## Null distribution

| Statistic | Value |
|---|---|
| Mean | 116.935 |
| Median | 117.091 |
| Std | 2.150 |
| Min | 98.884 |
| Max | 124.591 |
| 0.01th percentile | 106.463 |
| 0.1th percentile | 108.749 |
| 1st percentile | 111.221 |
| 5th percentile | 113.160 |

## Comparison with 10K result

| | 10K permutations | 1M permutations |
|---|---|---|
| Null ≤ observed | 0 / 10,000 | 0 / 1,000,000 |
| p-value | < 1 × 10⁻⁴ | < 1 × 10⁻⁶ |
| Null minimum | — | 98.884 |
| Gap (null min − observed) | — | 37.731 |

The 10K floor (p = 0.0001) was genuine — the true p-value is far below it. At 10⁶ permutations, the null distribution does not even approach the observed distance. A normal approximation on the null (mean = 116.94, σ = 2.15) gives a z-score of −25.9, corresponding to p ≈ 10⁻¹⁴⁸, though extreme tail extrapolation should be interpreted cautiously.

## Method

Same pipeline as `scripts/08_scaled_procrustes.py`:
- Input: PCA-reduced centroids (35 cell types × 33 PCs) from `output/phase2/scaled_35types/pca_centroids_35.npz`
- Null model: shuffle mouse↔human cell type label pairings (row permutation of mouse centroid matrix)
- p-value = (#{null ≤ observed} + 1) / (B + 1) — conservative formulation
- Seed: 42 (first 10K permutations reproduce the original 10K test exactly)
- Procrustes: OPA without reflection, implemented from scratch (same as `src/procrustes.py`)

## Output files

- `results_1M.json` — structured results
- `null_distribution_1M.npy` — 1M null distances (numpy array, 7.6 MB)
- `RESULTS_SUMMARY.md` — this file
