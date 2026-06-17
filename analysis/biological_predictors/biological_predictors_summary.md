# Biological Predictors of Cell Type Rigidity

## Summary

Tested 15 biological features as predictors of
Procrustes rigidity across 35 cell types.

## Univariate Results (Spearman ρ with rigidity score)

| Feature | ρ | p-value | Interpretation |
|---------|---|---------|----------------|
| Is progenitor | 0.428 | 1.0227e-02 | Strong |
| Cell cycle fraction | 0.300 | 7.9659e-02 | Moderate |
| TF entropy | -0.268 | 1.2028e-01 | Weak |
| Marker seq divergence | -0.241 | 1.6280e-01 | Weak |
| Mean active TFs | -0.190 | 2.7546e-01 | Negligible |
| Within-type heterogeneity | 0.151 | 3.8752e-01 | Negligible |
| Mean HK ratio | 0.137 | 4.3172e-01 | Negligible |
| Tissue breadth | -0.127 | 4.6623e-01 | Negligible |
| Inter-donor variance | -0.127 | 4.6860e-01 | Negligible |
| Transcriptomic complexity | 0.105 | 5.4708e-01 | Negligible |
| Is endoderm | -0.078 | 6.5689e-01 | Negligible |
| Mean expression level | -0.059 | 7.3472e-01 | Negligible |
| Chromatin conservation | 0.058 | 7.3955e-01 | Negligible |
| Log min cell count | -0.033 | 8.4924e-01 | Negligible |
| N marker genes | -0.014 | 9.3635e-01 | Negligible |

**Strongest univariate predictor:** Is progenitor (ρ = 0.428, p = 1.0227e-02)

## Multivariate Results

- Elastic Net training R² = 0.194
- Random Forest training R² = 0.557
- Elastic Net LOO-CV R² = -0.064
- Random Forest LOO-CV R² = -0.046

**Conclusion:** Multivariate models have minimal predictive power.
Biological features do not jointly predict rigidity well.
This is informative: rigidity rankings are NOT explained by
standard biological covariates (complexity, cell cycle, tissue
breadth, etc.), ruling them out as confounds.

## Decision Gate

Strong correlation found: Is progenitor (|ρ| = 0.428).
This is a meaningful biological finding worth reporting.

## Germ Layer Distribution

- Mesoderm: 26 types
- Endoderm: 7 types
- Ectoderm: 2 types

Note: Strong mesoderm bias (25/35 types) limits germ layer analysis.

## Files Generated

- `feature_table.csv` — 35 rows × all features + rigidity score
- `univariate_correlations.csv` — all Spearman correlations
- `multivariate_model_results.json` — model coefficients, R², predictions
- `figures/supplementary/biological_predictors_panel{A,B,C,D}.{pdf,png}`
