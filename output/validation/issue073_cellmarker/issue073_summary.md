# ISSUE-073: CellMarker Validation — Identity Gene Circularity Test

## PRIMARY RESULT (human):
- Observed overlap: 33 / 500 identity genes in CellMarker
- Expected overlap: 7.58
- Enrichment ratio: 4.355
- p-value: 1.12e-12
- **Verdict: PASS**

## EXPRESSION-CONTROLLED (human):
- Matched background universe: 2294 genes (500 identity + 1794 expression-matched controls)
- Matched background enrichment ratio: 2.857
- p-value: 1.44e-10
- Unmatched ratio: 4.355 → matched ratio: 2.857
- **Mean-variance trap: RULED OUT** (enrichment holds after expression-level control)

## MOUSE REPLICATION:
- Observed overlap: 17 / 500 identity genes in CellMarker (mouse)
- Expected overlap: 4.66
- Enrichment ratio: 3.648
- p-value: 4.37e-06
- **Verdict: PASS**

## PER CELL TYPE:
- 3 / 6 cell types pass individual enrichment test (p<0.05 AND ratio>1.5)
- Failing types: ['CD4-positive, alpha-beta T cell', 'CD8-positive, alpha-beta T cell', 'endothelial cell']
- No CellMarker data: []

### Per cell type table
| cell_type | n_loading | n_markers | overlap | enrichment | p_value | pass |
|-----------|-----------|-----------|---------|------------|---------|------|
| B cell | 20 | 34 | 1 | 94.21 | 1.06e-02 | PASS |
| CD4-positive, alpha-beta T cell | 20 | 40 | 0 | 0.0 | 1.00e+00 | FAIL |
| CD8-positive, alpha-beta T cell | 20 | 45 | 0 | 0.0 | 1.00e+00 | FAIL |
| endothelial cell | 20 | 49 | 0 | 0.0 | 1.00e+00 | FAIL |
| hepatocyte | 20 | 16 | 2 | 339.14 | 1.32e-05 | PASS |
| macrophage | 20 | 26 | 1 | 141.31 | 7.06e-03 | PASS |

## Overlapping identity genes (33)
AIF1, ALB, CD163, CD2, CD44, CD63, CD74, CD79A, CD83, CLU, CPS1, CXCR4, EMCN, ENG, EZR, FLT1, GLUL, ICAM1, IL1RAP, LGALS3, LY9, LYZ, MS4A1, NRP1, PECAM1, PTPRC, S100A6, SLC8A1, TEAD1, TF, TM4SF1, VIM, VWF

## Method
- CellMarker 2.0 (http://xteam.xbio.top/CellMarker/)
- Filtered to marker_source = "Experiment" only (wet-lab validated: flow cytometry, IHC, ISH)
- Human: 1101 unique marker genes across 327 cell types
- Mouse: 477 unique marker genes across 218 cell types
- Background: 16957 1:1 ortholog genes (full CellWarp gene space)
- Identity genes: top 500 by centroid variance across 6 human cell type centroids (Tabula Sapiens)
- Expression-matched control: 10 background genes matched within ±10% mean expression per identity gene
- Test: one-sided hypergeometric (scipy.stats.hypergeom)
- Note: top_genes_per_cell_type stored only 20 loading genes per type (not 50 as spec requested)


## PER CELL TYPE RERUN (n=50 loading genes)

Method: top 50 genes by absolute deviation of cell type centroid from mean centroid
across all 6 types (fallback method, computed directly from centroids_human.csv).
Previous run used only 20 loading genes from procrustes_results.json.

| cell_type | n_loading | n_markers | overlap | enrichment | p_value | pass |
|-----------|-----------|-----------|---------|------------|---------|------|
| CD4-positive, alpha-beta T cell | 50 | 40 | 1 | 30.72 | 3.21e-02 | PASS |
| CD8-positive, alpha-beta T cell | 50 | 45 | 3 | 84.48 | 5.26e-06 | PASS |
| endothelial cell | 50 | 49 | 3 | 77.98 | 6.83e-06 | PASS |

### Overlapping genes
- **CD4-positive, alpha-beta T cell**: PTPRC
- **CD8-positive, alpha-beta T cell**: CD2, CD8A, PTPRC
- **endothelial cell**: EMCN, PECAM1, VWF

### Comparison to original (n=20)
All three types had 0 overlap at n=20. At n=50:
- CD4-positive, alpha-beta T cell: 1 overlap (was 0)
- CD8-positive, alpha-beta T cell: 3 overlap (was 0)
- endothelial cell: 3 overlap (was 0)
