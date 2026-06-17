# Step 0B — RIRA Unknown Cell Inspection

**Date:** 2026-03-16

## Summary

- Total RIRA cells: 596,848
- Unknown: 167,180 (28.0%)
- Annotated: 429,668 (72.0%)

## 1. Distribution Across Tissues

| Tissue | Total | Unknown | Unknown % |
|---|---|---|---|
| PBMC | 98,968 | 37,259 | 37.6% |
| Spleen | 93,564 | 33,700 | 36.0% |
| Bone marrow | 76,123 | 32,890 | 43.2% |
| Mes. LN | 115,717 | 24,273 | 21.0% |
| PLN | 131,386 | 21,358 | 16.3% |
| Lung | 35,037 | 10,486 | 29.9% |
| Liver | 45,796 | 7,158 | 15.6% |

Unknown cells are distributed across all 7 tissues (diffuse, not concentrated).

## 2. Distribution Across Donors

- Donors with Unknown cells: 47/47
- Unknown % per donor: min 6.7%, max 91.6%, median 19.7%
- Spread across all donors — not a batch-specific artifact.

## 3. Gene Detection: Unknown vs Annotated

| Metric | Unknown | Annotated | Ratio |
|---|---|---|---|
| Median genes/cell | 711 | 905 | 0.786 |
| Mean genes/cell | 818 | 998 | 0.820 |
| Median UMI/cell | 1149 | 2146 | 0.535 |

## 4. scGateConsensus for Unknown Cells

| Label | Count |
|---|---|
| Ambiguous | 69,103 |
| T_NK | 22,539 |
| Erythrocyte | 17,395 |
| Uncalled | 16,344 |
| Bcell | 16,101 |
| Platelet | 14,782 |
| Myeloid | 8,335 |
| Stromal | 1,565 |
| Epithelial | 1,016 |

## 5. UMAP Coordinates

No UMAP coordinates in metadata. Spatial overlap analysis skipped.

## 6. Exclusion Recommendation

**EXCLUDE** all Unknown cells from centroid computation.

**Justification:** Unknown cells failed the atlas authors' scGate-based cell type classification. They have lower gene detection (median 711 vs 905, ratio 0.786) and lower UMI counts, consistent with lower-quality cells. They are diffuse across all tissues and donors (not batch-specific). The largest scGateConsensus category is 'Ambiguous' (69K cells) — cells where multiple markers gave conflicting signals. Including these in type-specific centroids would contaminate the signal with unclassified or mis-assigned cells.
