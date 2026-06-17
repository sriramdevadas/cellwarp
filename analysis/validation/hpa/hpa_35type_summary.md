# HPA Held-Out Marker Validation — 35-Type Summary

HPA download date: 2026-05-15T22:35:23Z

## Pooled (global) enrichment

- K (HPA markers in background): **14559**
- Observed overlap: **495** / 500
- Expected overlap: 429.24
- Fold-enrichment: **1.153**
- p-value: **8.3114e-27**

## Expression-matched control

- Universe size: 2738
- HPA markers in universe: 2326
- Observed overlap: **495** / 500
- Expected overlap: 424.76
- Fold-enrichment: **1.165**
- p-value: **5.2041e-32**

## Per-cell-type

| Cell type | HPA markers | Overlap | Fold | p | Result |
|---|---|---|---|---|---|
| B cell | 368 | 21/50 | 19.36 | 2.50e-22 | PASS |
| CD4-positive, alpha-beta T cell | 211 | 9/50 | 14.47 | 9.70e-09 | PASS |
| CD8-positive, alpha-beta T cell | 180 | 15/50 | 28.27 | 2.20e-18 | PASS |
| endothelial cell | 549 | 41/50 | 25.33 | 3.59e-53 | PASS |
| hepatocyte | 759 | 30/50 | 13.41 | 3.87e-28 | PASS |
| macrophage | 656 | 44/50 | 22.75 | 2.19e-56 | PASS |

## Coverage

HPA covers all 6 validated cell types via the following labels:

- **B cell** ← 'memory B-cell', 'naive B-cell'
- **CD4-positive, alpha-beta T cell** ← 'memory CD4 T-cell', 'naive CD4 T-cell'
- **CD8-positive, alpha-beta T cell** ← 'memory CD8 T-cell', 'naive CD8 T-cell'
- **endothelial cell** ← 'Lymphatic endothelial cells', 'Vascular endothelial cells'
- **hepatocyte** ← 'Hepatocytes'
- **macrophage** ← 'Kupffer cells', 'Macrophages'