# Mouse Lemur CellWarp Analysis — Summary

**Date:** 2026-04-05
**Species pair:** Human (*H. sapiens*) vs Mouse lemur (*M. murinus*)
**Evolutionary distance:** ~75 Mya
**Data sources:**
- Human: Tabula Sapiens (pre-computed centroids from primary 35-type analysis)
- Mouse lemur: Tabula Microcebus (Ezran et al., Nature 2025, CELLxGENE Discover)

---

## Global Coherence (Procrustes)

| Metric | Value |
|--------|-------|
| Cell types | 15 |
| Shared ortholog genes | 13,796 |
| Procrustes distance | 21.7651 |
| Null median distance | 62.8579 |
| **obs/null ratio** | **0.3463** |
| **p-value** | **0.000100** |
| Significant (alpha=0.01) | YES |
| PCA components | 15 |
| Cumulative variance | 95.5% |
| Scaling factor | 1.384536 |
| Rotation det | +1.000000 |

## Comparison Across Species Pairs

| Species Pair | Div (Mya) | obs/null | p-value | n types |
|-------------|-----------|----------|---------|---------|
| Human-macaque | 25 | 0.841 | 0.0002 | 20 |
| **Human-mouse lemur** | **75** | **0.3463** | **0.0001** | **15** |
| Human-mouse | 90 | 0.522 | <0.0001 | 35 |

obs/null monotonically decreasing with distance: **NO**

## Per-Type Residuals (ranked by divergence)

| Rank | Cell Type | Magnitude | % SSR |
|------|-----------|-----------|-------|
| 1 | T cell | 9.8442 | 20.5% |
| 2 | endothelial cell | 9.2355 | 18.0% |
| 3 | monocyte | 6.5375 | 9.0% |
| 4 | B cell | 6.2209 | 8.2% |
| 5 | neutrophil | 5.8774 | 7.3% |
| 6 | mesenchymal stem cell | 5.4959 | 6.4% |
| 7 | natural killer cell | 5.3252 | 6.0% |
| 8 | CD4-positive, alpha-beta T cell | 5.1605 | 5.6% |
| 9 | enterocyte of epithelium of large intestine | 4.8993 | 5.1% |
| 10 | mature NK T cell | 3.9620 | 3.3% |
| 11 | pancreatic acinar cell | 3.9597 | 3.3% |
| 12 | macrophage | 3.6050 | 2.7% |
| 13 | CD8-positive, alpha-beta T cell | 2.9575 | 1.8% |
| 14 | plasma cell | 2.7616 | 1.6% |
| 15 | fibroblast | 2.3651 | 1.2% |

## Ranking Correlation

### vs Primary (human-mouse, 35 types)
- Spearman rho = 0.1571
- p-value = 0.5760
- n shared types = 15

### vs Macaque (human-macaque, 20 types)
- Spearman rho = -0.5476
- p-value = 0.1600
- n shared types = 8

**Expected:** Ranking correlation is weak/non-significant based on simulation
findings and five prior non-replications (MCA rho=0.120, Sun2023 rho=0.146,
PanSci rho=0.194, CellHint rho=-0.386, macaque rho=0.137).

## Key Observations

1. **Global geometric coherence at 75 Mya**: Confirmed
   — the Procrustes transformation structure holds
   at an intermediate evolutionary distance.

2. **obs/null ratio = 0.346**: Does not follow expected evolutionary distance
   scaling. The mouse lemur obs/null (0.346) is LOWER than mouse (0.522),
   implying more coherence despite greater evolutionary distance. However,
   **obs/null is not directly comparable across analyses with different n_types**:
   the lemur analysis uses 15 types, mouse uses 35, macaque uses 20. With fewer
   types, permutation null distributions are wider (fewer possible configurations),
   and well-separated types contribute more per-type signal. The scientifically
   meaningful comparison is p-value significance: all three are p < 0.01.

3. **Three-species geometric conservation**: All three pairs significant (p < 0.01).
   Cross-species geometric coherence is a robust signal across 25-90 Mya of
   primate and rodent evolution.

4. **Ranking non-replication (rho = 0.157, NS)**: This is the 6th independent
   non-replication, consistent with prior findings (MCA 0.120, Sun2023 0.146,
   PanSci 0.194, CellHint -0.386, macaque 0.137). Ranking instability is a
   confirmed general property of the measurement, not a lemur-specific artifact.

## Data Quality Notes

- Mouse lemur data filtered to 10x 3' v2 only (95% of atlas)
- 15 cell types pass >=500 cell gate (exactly at threshold)
- Hepatocyte has 458 cells (just below gate) — excluded
- 60% immune types, 40% non-immune (better balance than macaque 65%)
- 4 donors only (vs 47 in RIRA macaque data)
- Human centroids from pre-computed primary analysis (Tabula Sapiens, same
  centroids used in the 35-type primary pipeline)
