# Within-Species Controls: Type-Count-Matched Reanalysis

**Date:** 2026-04-06
**Goal:** Test whether cross-species Procrustes coherence (obs/null = 0.522) is
distinguishable from the within-species baseline when controlling for the number
of cell types in each comparison.

## Key Finding

**No within-species pair has n ≥ 15 cell types.** The within-species pairs range
from n = 6 to n = 12 (median 7), while the cross-species comparisons use n = 15–35.
A direct type-count-matched comparison is therefore **not possible** with the
current data.

## What We Can Report

### Type-count distribution (24 within-species pairs)

| n_types | Count |
|---------|-------|
| 6       | 11    |
| 7       | 3     |
| 8       | 5     |
| 9       | 1     |
| 10      | 1     |
| 11      | 1     |
| 12      | 2     |

### Analysis at available thresholds

| Threshold  | Pairs | Mean obs/null | Median obs/null | 95% Bootstrap CI     | Cross-species (0.522) in CI? | Frac < 0.522 |
|------------|-------|---------------|-----------------|----------------------|------------------------------|--------------|
| n ≥ 12     | 2     | 0.349         | 0.349           | [0.269, 0.428]       | No                           | 2/2 (100%)   |
| n ≥ 10     | 4     | 0.415         | 0.427           | [0.309, 0.510]       | No                           | 3/4 (75%)    |
| n ≥ 8      | 10    | 0.449         | 0.458           | [0.378, 0.518]       | No                           | 8/10 (80%)   |
| n ≥ 6      | 24    | 0.466         | 0.457           | [0.407, 0.528]       | **Yes** (barely)             | 17/24 (71%)  |

### Effect sizes (all 24 pairs)

- **Cohen's d** (cross-species − within-species mean) / sd = 0.367 (small–medium)
- **Spearman correlation** (n_types vs obs/null): r = −0.183, p = 0.39 (not significant)
  - No evidence that higher type counts systematically change the obs/null ratio

### Largest within-species pairs (n ≥ 10)

| Pair | n_types | obs/null | p-value |
|------|---------|----------|---------|
| blood vs spleen | 12 | 0.269 | 0.0001 |
| lymph node vs spleen | 12 | 0.428 | 0.0001 |
| adipose tissue vs spleen | 11 | 0.426 | 0.0001 |
| blood vs lymph node | 10 | 0.538 | 0.0003 |

## Interpretation

1. **At n ≥ 8 (10 pairs), cross-species is marginally distinguishable.** The 95% CI
   upper bound (0.518) falls just below 0.522. Within-species pairs at higher type
   counts show stronger coherence (lower obs/null) than cross-species, consistent
   with the expected hierarchy.

2. **At all 24 pairs (n ≥ 6), the CI barely includes 0.522.** The cross-species
   ratio sits at the 71st percentile of the within-species distribution — higher
   than most within-species pairs, but not outside the range.

3. **The type-count mismatch is a genuine limitation.** Within-species pairs max
   out at n = 12 while cross-species uses n = 15–35. Although the Spearman test
   shows no significant n_types → obs/null correlation (p = 0.39), the sample is
   small and the power to detect such an effect is low.

4. **The hierarchy still holds directionally.** Self-comparison (0.033) << within-
   species (0.466) < cross-species (0.522) << null (1.0). Evolution degrades
   geometric coherence relative to within-species, even though the difference is
   modest (Cohen's d = 0.37).

## Honest Assessment

The cross-species Procrustes coherence is **real** (far below permutation null,
p = 0.0001) and **directionally weaker** than within-species coherence (mean 0.522
vs 0.466), consistent with evolutionary divergence. However, the distinction
between cross-species and within-species is **modest** — the cross-species value
falls within or near the upper edge of the within-species distribution. The lack
of type-count-matched comparisons (no within-species pairs at n ≥ 15) prevents a
fully controlled comparison.

**Bottom line:** The data support the hierarchy self < within-species < cross-species
< null, but the within-vs-cross gap is small (Cohen's d = 0.37) and the cross-
species ratio is not dramatically separated from the within-species distribution.

## Files

- `matched_control_results.json` — Full results with per-threshold statistics
- `matched_control_analysis.py` — Analysis script
- `matched_control_summary.md` — This file
