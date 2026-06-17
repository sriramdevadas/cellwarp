# Cross-Reference: Bootstrap Stability vs Cross-Atlas Consistency

**Date:** 2026-04-06
**Purpose:** Check whether types that are bootstrap-stable (narrow CI)
are also consistent when measured in independent atlases.

---

## Overall Correlation

| Metric | Value |
|--------|-------|
| Types in ≥2 replications | 20 |
| Spearman ρ (CI width vs mean rank shift) | **-0.410** |
| Spearman p-value | 0.0727 |
| Pearson r | -0.331 |
| Pearson p-value | 0.1537 |

**Interpretation:** Weak negative trend but not
statistically significant (p = 0.073). The two stability measures are
at best loosely coupled.

---

## (a) Stable in Bootstrap AND Consistent Across Atlases

Criteria: bootstrap CI width ≤ 3, mean rank shift ≤ Q25 (3.0), present in ≥2 replications.

| Cell Type | Primary Rank | Bootstrap CI | Category | Mean Shift | N Replications |
|-----------|:------------|:------------|:---------|:----------|:--------------|
| mesenchymal stem cell | 13 | [13, 15] | STABLE_MIDDLE | 3.0 | 2 |
| B cell | 17 | [15, 18] | STABLE_MIDDLE | 2.9 | 7 |

## (b) Stable in Bootstrap BUT Inconsistent Across Atlases

Criteria: bootstrap CI width ≤ 3, mean rank shift ≥ Q75 (6.3), present in ≥2 replications.

| Cell Type | Primary Rank | Bootstrap CI | Category | Mean Shift | N Replications |
|-----------|:------------|:------------|:---------|:----------|:--------------|
| stromal cell | 1 | [1, 1] | STABLE_FLEXIBLE | 9.5 | 2 |
| epithelial cell | 2 | [2, 2] | STABLE_FLEXIBLE | 7.0 | 5 |
| hepatocyte | 32 | [32, 34] | STABLE_RIGID | 9.8 | 6 |

## (c) Pan-Census Zero-Shift Types

These types had exactly zero rank shift in the Pan-Census (22-type) replication.

| Cell Type | Primary Rank | Bootstrap CI | Width | Category | Census Rank Shift |
|-----------|:------------|:------------|:------|:---------|:-----------------|
| myeloid leukocyte | 10 | [7.0, 12.0] | 5.0 | STABLE_FLEXIBLE | 0 |
| plasma cell | 12 | [9.0, 12.0] | 3.0 | STABLE_FLEXIBLE | 0 |
| B cell | 17 | [15.0, 18.0] | 3.0 | STABLE_MIDDLE | 0 |

## (d) CD8+ T Cell — Detailed Cross-Atlas Profile

- **Primary rank:** 35 of 35 (most rigid)
- **Bootstrap:** median rank = 35.0, CI = [35.0, 35.0], width = 0.0
- **Category:** STABLE_RIGID

| Replication | Rank | N Types | Primary Subset Rank | Rank Shift |
|-------------|:-----|:--------|:-------------------|:----------|
| Sun2023 | 14 | 15 | 15 | 1 |
| PanSci | 14 | 16 | 16 | 2 |
| CellHint | 8 | 15 | 15 | 7 |
| CellHint_harmonized | 11 | 12 | 12 | 1 |
| Pan_Census | 10 | 22 | 22 | 12 |
| Macaque | — | — | — | — |
| Mouse_lemur | 13 | 15 | 15 | 2 |

**Mean rank shift:** 4.2 (SD = 4.4, n = 6)

**Key finding:** CD8+ T cell has bootstrap CI width = 0 (perfectly locked at rank 35 under resampling) but shows substantial cross-atlas variability. It is always among the most rigid types in the primary analysis bootstrap, but its relative rank shifts considerably when measured in different atlases.

---

## Full Per-Type Summary (≥2 replications)

| Cell Type | Primary | Boot CI | Boot Cat | Mean Shift | SD Shift | N Repl |
|-----------|:--------|:--------|:---------|:----------|:---------|:-------|
| myeloid leukocyte | 10 | 5 | STABLE_FLEXIBLE | 0.0 | 0.0 | 2 |
| myeloid dendritic cell | 21 | 7 | STABLE_MIDDLE | 2.8 | 1.7 | 4 |
| B cell | 17 | 3 | STABLE_MIDDLE | 2.9 | 2.5 | 7 |
| macrophage | 16 | 5 | STABLE_MIDDLE | 2.9 | 3.1 | 7 |
| monocyte | 22 | 4 | STABLE_MIDDLE | 3.0 | 2.5 | 7 |
| mesenchymal stem cell | 13 | 2 | STABLE_MIDDLE | 3.0 | 2.8 | 2 |
| smooth muscle cell | 31 | 1 | STABLE_RIGID | 3.2 | 2.4 | 6 |
| natural killer cell | 23 | 3 | STABLE_MIDDLE | 3.4 | 1.8 | 5 |
| neutrophil | 8 | 4 | STABLE_FLEXIBLE | 3.8 | 4.3 | 4 |
| CD4-positive, alpha-beta T cell | 14 | 3 | STABLE_MIDDLE | 4.0 | 4.1 | 6 |
| CD8-positive, alpha-beta T cell | 35 | 0 | STABLE_RIGID | 4.2 | 4.4 | 6 |
| fibroblast | 28 | 2 | STABLE_RIGID | 4.9 | 3.6 | 7 |
| plasma cell | 12 | 3 | STABLE_FLEXIBLE | 5.4 | 3.0 | 7 |
| endothelial cell | 33 | 2 | STABLE_RIGID | 6.0 | 4.2 | 7 |
| T cell | 7 | 4 | STABLE_FLEXIBLE | 6.2 | 4.1 | 5 |
| mature NK T cell | 25 | 4 | STABLE_MIDDLE | 6.5 | 6.4 | 2 |
| epithelial cell | 2 | 0 | STABLE_FLEXIBLE | 7.0 | 4.8 | 5 |
| granulocyte | 27 | 7 | STABLE_MIDDLE | 7.3 | 4.6 | 3 |
| stromal cell | 1 | 0 | STABLE_FLEXIBLE | 9.5 | 12.0 | 2 |
| hepatocyte | 32 | 2 | STABLE_RIGID | 9.8 | 4.8 | 6 |

---
*Generated: 2026-04-06*