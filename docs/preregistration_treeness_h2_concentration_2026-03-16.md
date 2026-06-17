CellWarp Treeness-Rigidity Anticorrelation — H2 Gene Program Concentration
Date: 2026-03-16
Author: [you]
Status: Committed before analysis begins
Label: Post-hoc exploratory, motivated by DECISION-121 result

Background
DECISION-120 established that Liang-Wagner treeness is anticorrelated with
Procrustes rigidity (rho=-0.349, p=0.040, n=35). DECISION-121 rejected the
neighborhood density mediation hypothesis — the anticorrelation is not a
geometric proximity artifact. This pre-registration tests whether gene
program concentration explains the anticorrelation.

Hypothesis
Rigid cell types have more concentrated identity gene programs — a larger
fraction of within-type variance explained by few principal components,
lower effective dimensionality — while flexible types have more distributed
programs. In the Waddington attractor framing: rigid types occupy deep,
narrow attractors (concentrated programs); flexible types occupy shallow,
broad attractors (distributed programs).

Concentrated programs may simultaneously:
(a) constrain cross-species geometry (high rigidity — the narrow attractor
    is conserved because it is simple),
(b) degrade tree structure (low treeness — types with concentrated programs
    cluster tightly along shared dominant axes, violating the four-point
    condition).

If both hold, concentration mediates the treeness-rigidity anticorrelation.

Operationalization
Concentration is measured from the per-cell single-cell expression matrix
(not centroids) for each cell type in human Tabula Sapiens, filtered to
the top-50 loading genes per cell type (the CellWarp identity gene set).
We ask about concentration of the IDENTITY program specifically, not the
full transcriptome.

Primary measure: Participation ratio (PR)
  PR = (sum(lambda_i))^2 / sum(lambda_i^2)
  where lambda_i are PCA eigenvalues from within-type expression.
  PR is a standard effective dimensionality measure.
  Low PR = concentrated program (variance dominated by few PCs).
  High PR = distributed program (variance spread across many PCs).
  Range: 1 (all variance on one PC) to p (uniform across p PCs).

Secondary measure: Fraction of variance explained by PC1 alone.
  Simpler, less sensitive, reported as sensitivity check.

Both computed on cells x 50-gene matrices after PCA retaining all
components.

Tests (Spearman, n=35)

H2a: PR vs Procrustes rigidity (residual magnitude).
  Prediction: POSITIVE rho. Rigid types (low residual) have low PR
  (concentrated). Since we correlate PR with residual_magnitude:
  low PR pairs with low residual = positive rho.

H2b: PR vs treeness score.
  Prediction: POSITIVE rho. Concentrated types (low PR) have lower
  treeness. Since low PR = low treeness: positive rho.
  (Note: this prediction arises because concentrated types cluster
  along shared dominant axes, degrading tree topology.)

H2c (mediation): Partial Spearman correlation of residual_magnitude
  vs treeness controlling for PR. Compare to raw rho=0.349.
  Attenuation = (|raw rho| - |partial rho|) / |raw rho| x 100.
  Threshold: attenuation >= 50% for mediation confirmed.

Input Data
- Human Tabula Sapiens single-cell expression: data/phase1/human_aligned.h5ad
  (or the scaled 35-type download file)
- Top-50 loading genes per cell type: from existing CellWarp Procrustes
  analysis outputs (output/phase2/scaled_35types/)
- Procrustes residuals: output/phase2/scaled_35types/residuals_ranked.csv
- Treeness scores: output/liang_wagner/treeness_scores_per_celltype.csv

Falsification Conditions
Any ONE of the following closes H2:

1. H2a NS: PR does NOT predict rigidity (p >= 0.05 or wrong sign).
   Concentration is unrelated to rigidity — cannot mediate.

2. H2b NS: PR does NOT predict treeness (p >= 0.05 or wrong sign).
   Concentration is unrelated to treeness — cannot mediate.

3. Attenuation < 50%: partial mediation only. Report as contributing
   factor, defer full explanation to Paper 2.

4. Attenuation >= 50%: H2 confirmed as primary explanation.

Decision Table
| H2a | H2b | H2c (>=50%) | Interpretation |
|-----|-----|-------------|---------------|
| PASS | PASS | PASS | Concentration mediates — attractor depth explains anticorrelation |
| PASS | PASS | FAIL | Concentration contributes but does not fully explain |
| PASS | FAIL | — | Concentration predicts rigidity but not treeness |
| FAIL | — | — | Concentration does not predict rigidity — H2 rejected |

No post-hoc exclusions, alternative concentration metrics beyond
{PR, PC1 fraction}, or repeated analysis with different gene sets.
