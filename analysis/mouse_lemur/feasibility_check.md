# Mouse Lemur Feasibility Check

**Date:** 2026-04-05
**Status:** GO
**Organism:** Microcebus murinus (gray mouse lemur)
**Divergence from human:** ~75 Mya
**Data source:** CELLxGENE Discover (NOT Census — mouse lemur not in Census build 2025-11-08)
**Dataset:** Tabula Microcebus LCA complete (Ezran et al., Nature 2025)
**DOI:** 10.1038/s41586-025-09113-9

---

## 1. Data Availability

Mouse lemur is **not** in CELLxGENE Census (build 2025-11-08). Available Census
organisms: Homo sapiens, Mus musculus, Macaca mulatta, Callithrix jacchus,
Pan troglodytes. Marmoset and chimp are brain-only — not suitable for
multi-tissue analysis.

Tabula Microcebus data downloaded from CELLxGENE Discover:
- **URL:** https://cellxgene.cziscience.com/collections/a137437b-d284-4a27-b1e9-36958a8f92c1
- **File:** LCA_complete (dataset_id: a392ab34-9016-4f48-b45d-5b3a9cfa39fe)
- **Size:** 1.54 GB

| Metric | Value |
|--------|-------|
| Total cells | 244,081 |
| Genes | 16,400 |
| Unique cell types | 143 |
| Tissues | 27 (40 tissue terms) |
| Donors | 4 (L1, L2, L3, L4) |
| Assay (10x 3' v2) | 231,752 cells (95.0%) |
| Assay (Smart-seq2) | 12,329 cells (5.0%) |

### Top tissues by cell count

| Tissue | Cells |
|--------|-------|
| bone marrow | 35,959 |
| lung | 34,084 |
| skin of body | 25,333 |
| blood | 22,067 |
| kidney | 14,855 |
| bone element | 14,851 |
| spleen | 11,152 |
| liver | 8,365 |
| pancreas | 8,172 |

---

## 2. Cell Type Overlap with Primary 35 Types

**Passing (≥500 cells): 15**
**Gate (≥15): PASS**

| # | Cell Type | Lemur Cells | HM Rank (of 35) |
|---|-----------|-------------|------------------|
| 1 | neutrophil | 58,814 | 8 |
| 2 | CD4-positive, alpha-beta T cell | 17,076 | 14 |
| 3 | macrophage | 9,492 | 16 |
| 4 | monocyte | 9,260 | 22 |
| 5 | CD8-positive, alpha-beta T cell | 8,326 | 35 (most rigid) |
| 6 | B cell | 7,650 | 17 |
| 7 | fibroblast | 7,476 | 28 |
| 8 | natural killer cell | 5,036 | 23 |
| 9 | plasma cell | 4,531 | 12 |
| 10 | T cell | 3,478 | 7 |
| 11 | mature NK T cell | 1,965 | 25 |
| 12 | pancreatic acinar cell | 1,053 | 5 |
| 13 | endothelial cell | 902 | 33 |
| 14 | mesenchymal stem cell | 901 | 13 |
| 15 | enterocyte of epithelium of large intestine | 802 | 20 |

### Types found but below 500-cell gate (8 types)

| Cell Type | Lemur Cells | HM Rank |
|-----------|-------------|---------|
| hepatocyte | 458 | 32 |
| hematopoietic precursor cell | 444 | 3 |
| bladder urothelial cell | 404 | 29 |
| epithelial cell | 259 | 2 |
| pancreatic ductal cell | 204 | 30 |
| stromal cell | 178 | 1 (most flexible) |
| smooth muscle cell | 97 | 31 |
| large intestine goblet cell | 55 | 19 |

**NOTE:** Hepatocyte has 458 cells — just 42 below the 500 gate. This is a
near-miss. Hepatocyte is the most divergent type in human-macaque Procrustes
(21.5% of SSR). Its absence from the passing set limits comparison with the
macaque result.

### Types missing entirely (12 types)

hematopoietic stem cell, basal cell, fibroblast of cardiac tissue,
myeloid leukocyte, mesenchymal stem cell of adipose tissue,
luminal epithelial cell of mammary gland, classical monocyte,
myeloid dendritic cell, intermediate monocyte, non-classical monocyte,
adventitial cell, granulocyte.

### Composition analysis

- **Immune types**: 9/15 (60%) — B cell, CD4+ T, CD8+ T, T cell, macrophage,
  monocyte, NK, NKT, plasma cell
- **Non-immune types**: 6/15 (40%) — neutrophil, fibroblast, endothelial,
  mesenchymal stem cell, pancreatic acinar, enterocyte

Better immune/non-immune balance than macaque (65% immune). However, monocyte
subtypes (classical, intermediate, non-classical) are collapsed into generic
"monocyte" — reduces resolution for myeloid lineage.

---

## 3. Ortholog Depth

| Metric | Count | Gate | Status |
|--------|-------|------|--------|
| Human-mouse lemur 1:1 orthologs (BioMart) | 16,655 | ≥12,000 | **PASS** |
| Human-mouse 1:1 orthologs (baseline) | 17,187 | — | — |
| Three-way intersection (H-M-ML, BioMart) | 14,983 | ≥12,000 | **PASS** |
| Lemur orthologs present in h5ad | 14,441 | — | — |
| Usable genes (in lemur data & human data & 1:1) | 13,796 | ≥12,000 | **PASS** |
| Human-macaque 1:1 orthologs | 19,123 | — | — |

Ortholog depth is slightly lower than human-macaque (16,655 vs 19,123), as
expected for greater evolutionary distance (~75 Mya vs ~25-30 Mya). Still
comfortably above the 12,000 gene gate.

Gene space mapping: ENSMICG-prefixed Ensembl IDs in h5ad map directly to
BioMart lemur_ensembl_id. No gene symbol concordance issues (mapping by
Ensembl ID, not gene name).

---

## 4. GO / NO-GO Decision

**Decision: GO**

| Gate | Criterion | Actual | Status |
|------|-----------|--------|--------|
| Cell type overlap | ≥15 types with ≥500 cells | 15 | **PASS** (at threshold) |
| Ortholog depth | ≥12,000 usable genes | 13,796 | **PASS** |
| Technology | 10x compatible | 10x 3' v2 (95%) | **PASS** |

### Caveats

1. **Cell type gate is exactly at threshold** (15/15). One fewer passing type
   would trigger NO-GO. Hepatocyte at 458 is a near-miss.
2. **Immune-heavy** (60% immune types), though better than macaque (65%).
3. **Only 4 donors** — less power than RIRA macaque (47 donors for immune types).
4. **Mixed assays** — 95% 10x, 5% Smart-seq2. Should filter to 10x only for
   technology consistency, though this may drop some cells below the 500 gate.
5. **Monocyte subtypes collapsed** — classical/intermediate/non-classical all map
   to generic "monocyte". Reduces myeloid resolution.

### Proceed to Step 1: Data Preparation

Next steps:
1. Filter to 10x 3' v2 cells only (avoid technology mixing)
2. Map lemur Ensembl IDs to human orthologs
3. Normalize (library-size + log1p)
4. Run Procrustes with the 15 passing types (or adjust gate if 10x-only filtering drops types)

---

## Appendix: Data File Locations

| File | Path | Size |
|------|------|------|
| Tabula Microcebus h5ad | data/mouse_lemur/tabula_microcebus_LCA_complete.h5ad | 1.54 GB |
| Ortholog mapping | analysis/mouse_lemur/biomart_mouse_lemur_human_orthologs.csv | ~1 MB |
| This report | analysis/mouse_lemur/feasibility_check.md | — |
