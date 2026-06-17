# Bootstrap Ranking Confidence Intervals

**Date:** 2026-04-05
**Bootstrap iterations:** 1000
**Random seed:** 42
**Resampling:** With replacement, same n per type
**PCA:** Fresh per iteration (not projecting onto original PC space)
**Runtime:** 1206s

## Key Finding

Of 35 cell types:
- **35** have stable rankings (CI width <= 10)
- **0** have unstable rankings (CI width > 15)

## Stability Categories

### STABLE_FLEXIBLE (12 types)

- **stromal cell**: median rank 1, 95% CI [1, 1], width 0
- **epithelial cell**: median rank 2, 95% CI [2, 2], width 0
- **hematopoietic precursor cell**: median rank 3, 95% CI [3, 3], width 0
- **hematopoietic stem cell**: median rank 4, 95% CI [4, 4], width 0
- **pancreatic acinar cell**: median rank 5, 95% CI [5, 6], width 1
- **basal cell**: median rank 6, 95% CI [5, 6], width 1
- **T cell**: median rank 8, 95% CI [7, 11], width 4
- **neutrophil**: median rank 9, 95% CI [7, 11], width 4
- **myeloid leukocyte**: median rank 9, 95% CI [7, 12], width 5
- **fibroblast of cardiac tissue**: median rank 9, 95% CI [7, 12], width 5
- **mesenchymal stem cell of adipose tissue**: median rank 10, 95% CI [7, 12], width 5
- **plasma cell**: median rank 12, 95% CI [9, 12], width 3

### STABLE_RIGID (8 types)

- **fibroblast**: median rank 28, 95% CI [27, 29], width 2
- **bladder urothelial cell**: median rank 29, 95% CI [27, 29], width 2
- **pancreatic ductal cell**: median rank 30, 95% CI [30, 31], width 1
- **smooth muscle cell**: median rank 31, 95% CI [30, 31], width 1
- **hepatocyte**: median rank 32, 95% CI [32, 34], width 2
- **endothelial cell**: median rank 33, 95% CI [32, 34], width 2
- **non-classical monocyte**: median rank 34, 95% CI [32, 34], width 2
- **CD8-positive, alpha-beta T cell**: median rank 35, 95% CI [35, 35], width 0

### STABLE_MIDDLE (15 types)

- **mesenchymal stem cell**: median rank 13, 95% CI [13, 15], width 2
- **CD4-positive, alpha-beta T cell**: median rank 14, 95% CI [13, 16], width 3
- **macrophage**: median rank 16, 95% CI [14, 19], width 5
- **B cell**: median rank 16, 95% CI [15, 18], width 3
- **classical monocyte**: median rank 16, 95% CI [14, 18], width 4
- **luminal epithelial cell of mammary gland**: median rank 18, 95% CI [16, 19], width 3
- **large intestine goblet cell**: median rank 19, 95% CI [17, 21], width 4
- **enterocyte of epithelium of large intestine**: median rank 20, 95% CI [19, 21], width 2
- **myeloid dendritic cell**: median rank 21, 95% CI [18, 25], width 7
- **natural killer cell**: median rank 22, 95% CI [21, 24], width 3
- **monocyte**: median rank 22, 95% CI [21, 25], width 4
- **intermediate monocyte**: median rank 24, 95% CI [23, 26], width 3
- **mature NK T cell**: median rank 25, 95% CI [23, 27], width 4
- **adventitial cell**: median rank 26, 95% CI [24, 27], width 3
- **granulocyte**: median rank 27, 95% CI [22, 29], width 7

### MODERATE (0 types)

(none)

### UNSTABLE (0 types)

(none)

## Design Decisions

1. **Fresh PCA per iteration** (not projecting onto original PC space): The resampled centroids shift positions, so the principal axes of variation also shift. Refitting PCA tests the stability of the FULL pipeline end-to-end. Projecting onto the original PC space would understate uncertainty by holding the coordinate system fixed.

2. **Resampling with replacement** (classic nonparametric bootstrap): This tests sensitivity to which specific cells are sampled. The expected fraction of unique cells per type per iteration is ~63.2%.

3. **Same n per type**: Each bootstrap iteration preserves the original cell count per type, avoiding confounding sample-size effects with biological variation.

## Procrustes Distance Stability

- Mean distance: 61.2274
- SD: 0.2259
- CV: 0.0037
- Range: [60.4516, 61.9486]

## PCA Components

- Mean components: 33.0
- Range: [33, 33]

## Full Ranking Table

| Cell Type | Orig Rank | Median | 95% CI | Width | SD | Category |
|---|---|---|---|---|---|---|
| stromal cell | 1 | 1 | [1, 1] | 0 | 0.1 | STABLE_FLEXIBLE |
| epithelial cell | 2 | 2 | [2, 2] | 0 | 0.1 | STABLE_FLEXIBLE |
| hematopoietic precursor cell | 3 | 3 | [3, 3] | 0 | 0.1 | STABLE_FLEXIBLE |
| hematopoietic stem cell | 4 | 4 | [4, 4] | 0 | 0.1 | STABLE_FLEXIBLE |
| pancreatic acinar cell | 5 | 5 | [5, 6] | 1 | 0.5 | STABLE_FLEXIBLE |
| basal cell | 6 | 6 | [5, 6] | 1 | 0.5 | STABLE_FLEXIBLE |
| T cell | 7 | 8 | [7, 11] | 4 | 1.2 | STABLE_FLEXIBLE |
| neutrophil | 8 | 9 | [7, 11] | 4 | 1.3 | STABLE_FLEXIBLE |
| myeloid leukocyte | 10 | 9 | [7, 12] | 5 | 1.7 | STABLE_FLEXIBLE |
| fibroblast of cardiac tissue | 9 | 9 | [7, 12] | 5 | 1.4 | STABLE_FLEXIBLE |
| mesenchymal stem cell of adipose tissue | 11 | 10 | [7, 12] | 5 | 1.4 | STABLE_FLEXIBLE |
| plasma cell | 12 | 12 | [9, 12] | 3 | 0.9 | STABLE_FLEXIBLE |
| mesenchymal stem cell | 13 | 13 | [13, 15] | 2 | 0.6 | STABLE_MIDDLE |
| CD4-positive, alpha-beta T cell | 14 | 14 | [13, 16] | 3 | 0.7 | STABLE_MIDDLE |
| macrophage | 16 | 16 | [14, 19] | 5 | 1.3 | STABLE_MIDDLE |
| B cell | 17 | 16 | [15, 18] | 3 | 1.0 | STABLE_MIDDLE |
| classical monocyte | 15 | 16 | [14, 18] | 4 | 1.1 | STABLE_MIDDLE |
| luminal epithelial cell of mammary gland | 18 | 18 | [16, 19] | 3 | 0.8 | STABLE_MIDDLE |
| large intestine goblet cell | 19 | 19 | [17, 21] | 4 | 0.8 | STABLE_MIDDLE |
| enterocyte of epithelium of large intestine | 20 | 20 | [19, 21] | 2 | 0.6 | STABLE_MIDDLE |
| myeloid dendritic cell | 21 | 21 | [18, 25] | 7 | 1.6 | STABLE_MIDDLE |
| natural killer cell | 23 | 22 | [21, 24] | 3 | 0.9 | STABLE_MIDDLE |
| monocyte | 22 | 22 | [21, 25] | 4 | 1.3 | STABLE_MIDDLE |
| intermediate monocyte | 24 | 24 | [23, 26] | 3 | 0.9 | STABLE_MIDDLE |
| mature NK T cell | 25 | 25 | [23, 27] | 4 | 0.9 | STABLE_MIDDLE |
| adventitial cell | 26 | 26 | [24, 27] | 3 | 1.0 | STABLE_MIDDLE |
| granulocyte | 27 | 27 | [22, 29] | 7 | 1.7 | STABLE_MIDDLE |
| fibroblast | 28 | 28 | [27, 29] | 2 | 0.7 | STABLE_RIGID |
| bladder urothelial cell | 29 | 29 | [27, 29] | 2 | 0.6 | STABLE_RIGID |
| pancreatic ductal cell | 30 | 30 | [30, 31] | 1 | 0.5 | STABLE_RIGID |
| smooth muscle cell | 31 | 31 | [30, 31] | 1 | 0.5 | STABLE_RIGID |
| hepatocyte | 32 | 32 | [32, 34] | 2 | 0.7 | STABLE_RIGID |
| endothelial cell | 33 | 33 | [32, 34] | 2 | 0.8 | STABLE_RIGID |
| non-classical monocyte | 34 | 34 | [32, 34] | 2 | 0.7 | STABLE_RIGID |
| CD8-positive, alpha-beta T cell | 35 | 35 | [35, 35] | 0 | 0.1 | STABLE_RIGID |

## Pairwise Ordering Reliability

- Total pairs: 595
- Reliable ordering (P>0.9 or P<0.1): 556 (93%)
- Coin flip (0.4<P<0.6): 7 (1%)
