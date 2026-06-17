# Step 0C — Qu et al. Per-Cell-Type Donor Counts

**Date:** 2026-03-16

## Samples

Qu et al. has 20 samples across 16 tissues (4 with replicates).
No explicit donor/animal IDs in 10x metadata files.
Paper states: adult M. fascicularis, 1-2 animals (DECISION-099).

| GSM | Tissue | Cells |
|---|---|---|
| GSM5901076 | Adipose | 7,706 |
| GSM5901077 | Aorta | 11,538 |
| GSM5901078 | Bladder | 8,823 |
| GSM5901079 | Breast | 11,003 |
| GSM5901080 | Colon1 | 11,015 |
| GSM5901081 | Colon2 | 11,398 |
| GSM5901082 | Heart | 10,977 |
| GSM5901083 | Kidney | 13,459 |
| GSM5901084 | Liver1 | 8,901 |
| GSM5901085 | Liver2 | 12,105 |
| GSM5901086 | Lung | 12,612 |
| GSM5901087 | Muscle | 14,939 |
| GSM5901088 | Spleen | 10,089 |
| GSM5901089 | Stomach | 8,881 |
| GSM5901090 | Testis1 | 11,745 |
| GSM5901091 | Testis2 | 10,536 |
| GSM5901092 | Tongue | 5,736 |
| GSM5901093 | Trachea | 8,485 |
| GSM5901094 | Uterus1 | 16,544 |
| GSM5901095 | Uterus2 | 24,390 |

**Total: 230,882 cells**

## Donor Count Assessment

**Critical limitation:** Qu et al. metadata contains NO explicit donor/animal IDs.
Per DECISION-099, the paper reports 1-2 animals. Without metadata confirmation,
the true biological donor count cannot be determined from the data files alone.

**Proxy used:** Number of tissue samples contributing cells of each type.
This is an UPPER BOUND — if all samples come from 1 animal, the true donor count is 1.

## Per-Type Assessment

| Cell Type | Source Tissues | Total Cells (all tissues) | Sample Count | LOW-CONFIDENCE |
|---|---|---|---|---|
| hepatocyte | Liver tissue only | 21,006 | 2 | YES |
| endothelial cell | Endothelial cells present in most tissues | 208,601 | 18 | NO |
| fibroblast | Mesenchymal cells in multiple tissues | 168,683 | 14 | NO |
| smooth muscle cell | Vascular/visceral smooth muscle | 112,051 | 9 | NO |
| epithelial cell | Epithelial cells in multiple tissues | 132,346 | 11 | NO |
| bladder urothelial cell | Bladder tissue only | 8,823 | 1 | YES |
| stromal cell | Stromal cells across multiple tissues | 164,551 | 13 | NO |
| pancreatic acinar cell | NHPCA-sourced (no pancreas in Qu) | 0 | 0 | N/A (not in Qu) |
| pancreatic ductal cell | NHPCA-sourced (no pancreas in Qu) | 0 | 0 | N/A (not in Qu) |

## Summary

- Types with ≥3 tissue sources (optimistic case): 5
- Types with <3 tissue sources: 2
- Types not in Qu: 2

**Under worst-case assumption (1 animal):** ALL Qu-sourced types are LOW-CONFIDENCE.
**Under best-case assumption (different animal per replicate):** Types with ≥3 tissue sources
could have ≥3 donors, but this requires tissue replicates to come from different animals.

**Recommendation:** Flag ALL Qu-sourced types as LOW-CONFIDENCE until donor count is
confirmed from the paper's methods section or supplementary materials.
