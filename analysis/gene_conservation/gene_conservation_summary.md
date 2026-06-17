# Gene-level Expression Conservation vs Evolutionary Constraint

## Summary
- **Primary conservation metric**: Per-gene Pearson correlation across 35 matched cell-type centroids (human vs mouse)
- **Genes with valid conservation**: 15,940 / 16,959

## Constraint Metrics Tested

### LOEUF (gnomAD)
- Spearman rho = -0.0877 (p = 8.06e-28, n = 15,499)
- 95% CI: [-0.1037, -0.0714]
- Direction: correct
- Partial rho = -0.0101 (p = 2.11e-01, controlling for expression level)

### Sequence % identity
- Spearman rho = 0.0493 (p = 4.77e-10, n = 15,940)
- 95% CI: [0.0340, 0.0650]
- Direction: correct
- Partial rho = -0.0322 (p = 4.66e-05, controlling for expression level)

### pLI (gnomAD)
- Spearman rho = 0.0243 (p = 2.50e-03, n = 15,499)
- 95% CI: [0.0072, 0.0401]
- Direction: correct
- Partial rho = -0.0450 (p = 2.16e-08, controlling for expression level)

### Missense Z (gnomAD)
- Spearman rho = 0.0460 (p = 9.49e-09, n = 15,548)
- 95% CI: [0.0317, 0.0612]
- Direction: correct
- Partial rho = -0.0211 (p = 8.40e-03, controlling for expression level)

## Expression Level Confound
- Expression vs conservation: rho = 0.2612 (p = 6.95e-247)
- Expression level is a strong confound: highly expressed genes show higher cross-species conservation AND tend to be more constrained.

## Sensitivity Check (Procrustes Contribution Metric)
- vs LOEUF (gnomAD): rho = -0.2938 (p = 0.00e+00)
- Partial rho = -0.0269 (p = 6.97e-04)

## Decile Analysis (LOEUF (gnomAD))
| Decile | Mean Conservation | SEM | n |
|--------|-------------------|-----|---|
| 1 | 0.3675 | 0.0076 | 1551 |
| 2 | 0.3847 | 0.0075 | 1561 |
| 3 | 0.3852 | 0.0077 | 1541 |
| 4 | 0.3701 | 0.0076 | 1549 |
| 5 | 0.3779 | 0.0078 | 1551 |
| 6 | 0.3701 | 0.0077 | 1556 |
| 7 | 0.3597 | 0.0079 | 1541 |
| 8 | 0.3551 | 0.0081 | 1551 |
| 9 | 0.3141 | 0.0081 | 1553 |
| 10 | 0.2926 | 0.0084 | 1545 |

## Decision Gate: **PARTIAL**

Raw correlations are significant in the expected direction, but the effect is substantially confounded by expression level. After controlling for mean expression, the signal weakens or reverses. The geometric framework captures expression magnitude effects that correlate with constraint, rather than independent evolutionary signal.

## Files
- `gene_conservation_table.csv`: Per-gene conservation and constraint values
- `correlation_results.json`: Full correlation results with CIs
- `figures/supplementary/gene_conservation_scatter.png`: Hexbin scatter (Panel A, LOEUF)
- `figures/supplementary/gene_conservation_scatter_seqid.png`: Hexbin scatter (seq identity)
- `figures/supplementary/gene_conservation_deciles.png`: Decile bars (Panel B, LOEUF)
- `figures/supplementary/gene_conservation_deciles_seqid.png`: Decile bars (seq identity)
- `figures/supplementary/gene_conservation_multi_metric.png`: Multi-metric comparison
- `figures/supplementary/gene_conservation_histogram.png`: Conservation distribution
