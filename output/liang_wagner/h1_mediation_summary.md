# Neighborhood Density Mediation Test — Summary

**Date:** 2026-03-16 14:50
**Pre-registration:** docs/preregistration_treeness_anticorrelation_2026-03-16.md
**Input:** 35 human cell type centroids x 16,959 ortholog genes

## Results Table

| Metric | k=3 | k=5 (primary) | k=10 |
|--------|-----|---------------|------|
| H1: rho(density, residual) | 0.069 (p=0.692) | 0.095 (p=0.586) | 0.095 (p=0.586) |
| H2: rho(density, treeness) | 0.753 (p=0.000)* | 0.722 (p=0.000)* | 0.628 (p=0.000)* |
| H3: partial rho (attenuation) | 0.452 (atten=-30%) | 0.406 (atten=-17%) | 0.373 (atten=-7%) |
| Verdict | MEDIATION_REJECTED | MEDIATION_REJECTED | MEDIATION_REJECTED |

*Significance: * = p < 0.05*

## Falsification Conditions

1. H1 fails (density !~ rigidity): TRIGGERED
2. H2 fails (density !~ treeness): not triggered
3. H3 fails (attenuation < 50%): TRIGGERED
4. Sensitivity instability: not triggered

## Primary Verdict (k=5)

**MEDIATION_REJECTED**

Density does not predict rigidity — mediation hypothesis rejected
