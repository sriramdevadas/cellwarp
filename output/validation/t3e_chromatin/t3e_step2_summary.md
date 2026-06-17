# T3-E Step 2: Regulatory Sequence Conservation vs Procrustes Rigidity

Generated: 2026-05-08 20:34

## 1. UCSC Tracks

**Primary:** phastCons20way — 20 species (17 primates + treeshrew + mouse + dog).
Includes Mus musculus (mm10). Primate-dominated alignment; conservation scores
primarily reflect primate constraint with mouse as an outgroup anchor.

**Sensitivity:** phastCons100way — 100 vertebrates (mammals + birds + reptiles +
amphibians + fish). Broadest taxonomic scope.

Both tracks verified via md5 checksum.

## 2. Primary Spearman Result

| Metric | Value |
|--------|-------|
| Spearman ρ | -0.0583 |
| p-value | 0.739546 |
| n (cell types) | 35 |
| 95% CI | [-0.3841, 0.2804] |
| Track | placental_20way (phastCons20way) |
| Option | A (cell-type-specific top-50 loading genes) |
| Window | ±2kb |

**Pre-registered threshold: 8TH NULL TRIGGERED** — ρ = -0.058 < 0.35. Close chromatin/regulatory sequence as proximate mechanism.

## 3. Sensitivity Analysis

| Track | Option | Window (kb) | ρ | p | Conclusion |
|-------|--------|-------------|------|-------|-----------|
| placental_20way | A | 1 | -0.1090 | 0.533228 | NULL_TRIGGERED |
| placental_20way | A | 2 | -0.0583 | 0.739546 | NULL_TRIGGERED |
| placental_20way | A | 5 | -0.0095 | 0.956697 | NULL_TRIGGERED |
| placental_20way | B | 1 | 0.1930 | 0.266653 | NULL_TRIGGERED |
| placental_20way | B | 2 | 0.2017 | 0.245312 | NULL_TRIGGERED |
| placental_20way | B | 5 | 0.1829 | 0.292921 | NULL_TRIGGERED |
| 100way_vertebrate | A | 1 | -0.0798 | 0.648485 | NULL_TRIGGERED |
| 100way_vertebrate | A | 2 | -0.0356 | 0.839228 | NULL_TRIGGERED |
| 100way_vertebrate | A | 5 | -0.0616 | 0.725085 | NULL_TRIGGERED |
| 100way_vertebrate | B | 1 | 0.1832 | 0.292170 | NULL_TRIGGERED |
| 100way_vertebrate | B | 2 | 0.1854 | 0.286204 | NULL_TRIGGERED |
| 100way_vertebrate | B | 5 | 0.0675 | 0.700006 | NULL_TRIGGERED |

**All sensitivity results are consistent with primary conclusion.**

## 4. Data Quality Flags

- Total conservation scores computed: 420
- Genes resolved via Ensembl REST: 658/670 (12 MT-excluded, 0 failed)
- Option A entries with <80% gene coverage: 0

## 5. Partial Correlation

Partial ρ (controlling for mean expression level) = -0.1582, p = 0.371490
Partial correlation does not substantially change the conclusion.

## Files Generated

- `conservation_scores.csv` — All conservation scores
- `rigidity_conservation_merged.csv` — Primary analysis merged data (35 rows)
- `spearman_primary_result.json` — Primary Spearman result
- `sensitivity_table.csv` — Sensitivity grid (12 rows)
- `identity_gene_tss_hg38.bed` — TSS coordinates (BED6)
- `scatter_primary.png` — Primary result scatter plot
- `sensitivity_heatmap.png` — Sensitivity heatmap (placental)
- `sensitivity_heatmap_100way.png` — Sensitivity heatmap (100way)
- `download_log.txt` — BigWig download sizes and checksums
- `ucsc_track_availability.txt` — UCSC track availability report
