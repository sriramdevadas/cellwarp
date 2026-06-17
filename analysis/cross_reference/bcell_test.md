# B Cell Cross-Atlas Ranking Consistency: Permutation Test

## Question
Is B cell's mean rank shift of 3.14 across 7 replications significantly better
than expected by chance?

## Method
- 10,000 permutations
- For each permutation: within each replication, shuffle rank assignments among
  the types present (preserving subset composition)
- Rank shift = |primary_subset_rank - replication_rank| (primary re-ranked
  within each replication's subset)
- For each type with 7+ replications, compute mean rank shift
- Test statistic: minimum mean shift achieved by ANY type with 7+ reps
- P-value: fraction of permutations where any type achieves mean shift ≤ 3.14

## Result: NOT SIGNIFICANT

| Metric | Value |
|--------|-------|
| B cell observed mean shift | 3.14 |
| B cell replications | 7 |
| Null mean (B cell) | 4.25 ± 0.99 |
| B cell z-score | -1.12 |
| **B cell type-specific p** | **0.150** |
| **Global p (any type ≤ 3.14)** | **0.407** |

B cell's observed mean shift falls at the 41st percentile of the null
distribution of minimum-across-types. After correcting for the fact that we
have 6 types with 7+ replications (any one of which could be the "winner"),
the result is entirely consistent with chance.

## All Types with 7+ Replications

| Cell Type | Observed | Null Mean ± SD | Type p | z |
|-----------|:--------:|:--------------:|:------:|:-:|
| Macrophage | 3.00 | 4.47 ± 1.09 | 0.097 | -1.35 |
| B cell | 3.14 | 4.25 ± 0.99 | 0.150 | -1.12 |
| CD4+ T cell | 3.57 | 4.86 ± 1.27 | 0.171 | -1.02 |
| **CD8+ T cell** | **4.43** | **7.71 ± 1.83** | **0.040** | **-1.79** |
| Endothelial cell | 5.71 | 6.73 ± 1.74 | 0.298 | -0.58 |
| Fibroblast | 6.00 | 4.95 ± 1.30 | 0.804 | +0.81 |

## Expected Mean Shift Under Null
- Per-type null average: 5.50
- Heuristic (n_types/3): varies 4.0–7.3 per replication (replication sizes 12–22)
- B cell's null expectation (4.25) is below average because it sits near the
  middle of the primary ranking (rank 17/35), making small shifts easier to
  achieve by chance in smaller subsets

## Surprise Finding
CD8+ T cell is the only type approaching significance (p = 0.040, z = -1.79).
Its observed shift of 4.43 vs null 7.71 reflects the fact that CD8+ T cell sits
at rank 35 (most rigid), giving it a very high expected shift under random
permutation — yet it consistently lands near the rigid end across all 7
replications. This is NOT significant after Bonferroni correction (6 tests,
threshold 0.0083), but it is the strongest signal.

## Interpretation
B cell's cross-atlas consistency (mean shift 3.14) is a modest result, not a
statistically exceptional one. The low shift is partly explained by B cell's
mid-table primary rank (17/35), which mechanically limits possible deviations
in smaller replication subsets. No individual cell type achieves significance
after multiple-testing correction.

---
*Generated: 2026-04-06 | Script: bcell_permutation_test.py | 10,000 permutations, seed=42*
