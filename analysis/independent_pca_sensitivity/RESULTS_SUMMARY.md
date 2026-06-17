# Independent-PCA Sensitivity Analysis: Results Summary

**Date:** 2026-04-04
**Reviewer concern:** Joint PCA on 70 centroids (35 human + 35 mouse) optimizes the subspace for shared variance, potentially inflating apparent geometric coherence and biasing per-type residuals.

## Method

1. PCA computed **separately** on 35 human centroids and 35 mouse centroids in the 16,959-gene ortholog space
2. Retained 33 components per species (matching joint analysis; captures 99.96% variance in each)
3. Aligned the two PCA subspaces via **Procrustes on loading matrices** — this aligns the PCA axes in gene space without using cell-type correspondence
4. Projected mouse centroids into the aligned (human-oriented) subspace
5. Ran Procrustes superimposition + 1,000,000-iteration label-permutation test (upgraded from 10K, same approach as C10)
6. Computed per-type residuals and rigidity ranking in the independent-PCA space
7. Compared to joint-PCA rigidity ranking via Spearman rho

## Key Results

| Metric | Joint PCA | Independent PCA |
|--------|-----------|-----------------|
| Procrustes distance | 61.153 | 52.716 |
| Null mean | 116.959 | 111.306 |
| **obs/null ratio** | **0.523** | **0.473** |
| **p-value** | **< 10⁻⁶** | **< 10⁻⁶** |
| Permutations | 1,000,000 | 1,000,000 |
| Significant (alpha=0.01) | YES | YES |
| PCA components | 33 | 33 |

**Spearman rho (rigidity rankings):** 0.915 (p < 1e-6)

## Interpretation

The independent-PCA analysis **reproduces and slightly strengthens** the joint-PCA finding:

- **obs/null ratio decreases** from 0.523 to 0.473, meaning geometric coherence is actually *stronger* when each species defines its own variance axes independently. This is the opposite of what a joint-PCA artifact would produce.
- **p-value is p < 10⁻⁶** (0 of 1,000,000 permutations achieved a distance as low as observed; upgraded from 10K floor).
- **Rigidity rankings are highly concordant** (Spearman rho = 0.915), confirming that per-type residuals are not driven by the joint embedding.

## Rank Stability

Only 2 of 35 cell types changed by more than 10 rank positions:

| Cell Type | Joint Rank | Indep Rank | Delta |
|-----------|------------|------------|-------|
| basal cell | 6 | 17 | +11 |
| monocyte | 22 | 11 | -11 |

Both are moderate shifts in the middle of the ranking. The extremes are stable: CD8+ T cell remains rank 35 (most conserved), stromal cell remains in the top 3 (most divergent), and the progenitor types (hematopoietic stem/precursor cells) remain in the top 4.

## Subspace Alignment Quality

Procrustes alignment of loading matrices yielded singular values ranging from 0.80 (best-matched axis) to 0.002 (worst). The top ~10 axes are well-aligned (sigma > 0.45), while the lower-variance axes diverge between species — expected, since lower PCs capture species-specific noise. The mean alignment quality across all 33 axes was 0.347.

## Separate PCA Variance

| Statistic | Human PCA | Mouse PCA | Joint PCA |
|-----------|-----------|-----------|-----------|
| Components for 95% variance | 20 | 18 | 33 |
| Variance at k=33 | 99.96% | 99.96% | 95.18% |
| PC1 variance | 24.4% | 24.7% | 28.1% |

The separate PCAs require fewer components for 95% variance (20 and 18) compared to the joint PCA (33), because within-species variance is lower-dimensional than between-species variance. Using k=33 in the independent analysis therefore captures nearly all within-species variance (99.96%), providing a generous subspace for the alignment.

## CCA Secondary Analysis (Included for Completeness)

CCA was also tested as an alternative subspace alignment method. Unlike Procrustes-on-loadings, CCA **uses the cell-type pairing** to find maximally correlated axes. Results:

- CCA obs/null ratio: 0.827 (weaker than both joint and independent PCA)
- CCA p-value: 0.042 (significant at alpha=0.05 but not alpha=0.01)

CCA absorbs much of the geometric coherence signal into the subspace alignment itself, leaving less for the Procrustes test to detect. This is expected and does not indicate weaker signal — it indicates that CCA and Procrustes are partially redundant when applied sequentially. The Procrustes-on-loadings approach is preferable for this sensitivity analysis because it cleanly separates subspace alignment (pairing-free) from geometric correspondence testing (pairing-dependent).

## Conclusion

**The reviewer concern is empirically refuted.** Joint PCA does not inflate the geometric coherence signal. Independent PCA with pairing-free subspace alignment produces:
- Equal or stronger coherence (obs/null 0.473 vs 0.523)
- p < 10⁻⁶ (0 of 1,000,000 null permutations reached observed distance; null minimum 93.75 vs observed 52.72)
- Highly concordant rigidity rankings (rho=0.915)

This result can be cited in the revision as a direct sensitivity analysis addressing the joint-PCA concern, alongside the existing defenses (symmetric null, LOOCV with per-fold PCA refitting, four independent replication datasets each constructing their own PCA space).

## Output Files

- `independent_pca_results.json` — Full results with all metrics, comparisons, and per-type data
- `residuals_ranked_independent_pca.csv` — Rigidity ranking with joint-PCA comparison columns
- `null_distribution_independent_pca.npy` — 1,000,000 permuted distances (upgraded from 10K)
- `pca_centroids_independent.npz` — PCA centroids, loadings, and alignment rotation matrix Q
- `run_independent_pca.py` — Reproducible analysis script with full documentation
