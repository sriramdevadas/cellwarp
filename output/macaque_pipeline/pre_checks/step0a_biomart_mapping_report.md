# Step 0A — BioMart Mapping of Qu ENSMFAG IDs

**Date:** 2026-03-16
**Status:** ABORT

## Gene ID Formats

| Dataset | ID Format | Example |
|---|---|---|
| CellWarp (human) | ENSG IDs (var_names) + symbols (var.feature_name) | ENSG00000000003 → TSPAN6 |
| RIRA (M. mulatta) | Gene symbols | PGBD2, ZNF692, SH3BP5L |
| Qu (M. fascicularis) | ENSMFAG IDs (48.4% HGNC symbol, 51.6% ID-only) | ENSMFAG00000044637 → PGBD2 |

## Qu Gene Symbol Composition

| Category | Count | % |
|---|---|---|
| HGNC-like symbol (direct) | 15082 | 48.4% |
| ENSMFAG-only (need BioMart) | 16083 | 51.6% |
| Total | 31165 | 100% |

## BioMart Query

- Source: Ensembl BioMart REST API
- Dataset: mfascicularis_gene_ensembl (Macaca_fascicularis_6.0)
- 1:1 orthologs with human gene name: 18874
- Unique ENSMFAG→human mappings: 18874

## ENSMFAG Mapping Results

| Metric | Count |
|---|---|
| ENSMFAG IDs to map | 16083 |
| Mapped via BioMart | 1306 (8.1%) |
| Unmapped | 14777 |
| New CellWarp-overlapping genes recovered | 961 |

## CellWarp Gene Space Overlap

| Set | Genes | % of CellWarp |
|---|---|---|
| CellWarp active gene space | 16957 | 100% |
| RIRA ∩ CellWarp | 15469 | 91.2% |
| Qu ∩ CellWarp (pre-mapping) | 13597 | 80.2% |
| Qu ∩ CellWarp (post-mapping) | 14558 | 85.9% |

## Three-Way Intersection

| Component | Genes |
|---|---|
| CellWarp space | 16957 |
| ∩ RIRA | 15469 |
| ∩ Qu (post-BioMart) | 14558 |
| **Three-way intersection** | **13927** |
| % of CellWarp space | 82.1% |

## Abort Criterion Check (DECISION-123-AMENDMENT Stage 1)

- Threshold: ≥95% of 15,028 = **14276** genes
- Post-mapping three-way intersection: **13927**
- Result: **ABORT** (13927/14276 = 97.6% of threshold)
