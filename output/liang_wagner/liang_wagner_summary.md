# Liang-Wagner Treeness Analysis — Summary

**Date:** 2026-05-08 20:29
**Input:** 35 human cell type centroids x 16,959 ortholog genes
**Method:** Liang & Wagner (2015) delta statistic, Holland et al. (2002) analytic p-values, Storey (2002) pi0 estimation

## Step 1 — Overall Treeness

- Tetrads analyzed: 52,360
- delta distribution: mean=0.2478, median=0.1569, std=0.2489
- delta range: [0.0000, 1.0000]
- Tetrads with p < 0.05: 703 (1.3%)
- Tetrads with p < 0.01: 140 (0.3%)
- Storey pi0: 1.0000
- Tree-like fraction (1 - pi0): 0.0000

## Step 2 — Per-Cell-Type Treeness Ranking

| Rank | Cell Type | Treeness (mean delta) |
|------|-----------|----------------------|
| 1 | hematopoietic stem cell | 0.3298 |
| 2 | granulocyte | 0.3213 |
| 3 | plasma cell | 0.3168 |
| 4 | hepatocyte | 0.3122 |
| 5 | neutrophil | 0.3075 |
| 6 | pancreatic acinar cell | 0.2924 |
| 7 | enterocyte of epithelium of large intestine | 0.2871 |
| 8 | large intestine goblet cell | 0.2744 |
| 9 | hematopoietic precursor cell | 0.2740 |
| 10 | macrophage | 0.2650 |
| 11 | luminal epithelial cell of mammary gland | 0.2626 |
| 12 | myeloid dendritic cell | 0.2616 |
| 13 | myeloid leukocyte | 0.2601 |
| 14 | basal cell | 0.2529 |
| 15 | stromal cell | 0.2520 |
| 16 | epithelial cell | 0.2498 |
| 17 | pancreatic ductal cell | 0.2496 |
| 18 | endothelial cell | 0.2484 |
| 19 | fibroblast of cardiac tissue | 0.2421 |
| 20 | bladder urothelial cell | 0.2343 |
| 21 | B cell | 0.2342 |
| 22 | natural killer cell | 0.2316 |
| 23 | T cell | 0.2284 |
| 24 | CD4-positive, alpha-beta T cell | 0.2192 |
| 25 | mesenchymal stem cell of adipose tissue | 0.2184 |
| 26 | smooth muscle cell | 0.2174 |
| 27 | mature NK T cell | 0.2090 |
| 28 | non-classical monocyte | 0.2089 |
| 29 | intermediate monocyte | 0.2080 |
| 30 | CD8-positive, alpha-beta T cell | 0.2069 |
| 31 | adventitial cell | 0.2043 |
| 32 | classical monocyte | 0.2034 |
| 33 | monocyte | 0.1987 |
| 34 | mesenchymal stem cell | 0.1961 |
| 35 | fibroblast | 0.1951 |

## Step 3 — Correlation with Procrustes Rigidity

- Spearman rho(treeness, rigidity): **-0.3487** (p=0.0401)
- n = 35
- Outcome: **ANTICORRELATED**
- Interpretation: Most surprising — rigid types have less tree structure

## Step 4 — PCA Sensitivity

- PCA components: 20 (95.6% variance)
- delta mean (PCA space): 0.2490
- pi0 (PCA space): 1.0000
- rho(treeness_pca, rigidity): -0.3594 (p=0.0340)
- Full-space vs PCA consistency: rho=0.9908 (p=0.000000)

## Decision

**Flag for advisor** — anticorrelated, potentially most interesting finding.
