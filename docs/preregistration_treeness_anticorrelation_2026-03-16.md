CellWarp Treeness-Rigidity Anticorrelation — Mediation Pre-Registration
Date: 2026-03-16
Author: [you]
Status: Committed before mediation analysis begins

Background
DECISION-120 established that Liang-Wagner treeness (within-species
developmental distinctiveness) is anticorrelated with Procrustes rigidity
(cross-species functional conservation): rho=-0.349, p=0.040, n=35. This
pre-registration tests whether neighborhood density in expression space
mediates this anticorrelation — i.e., whether cell types in dense
neighborhoods (many nearby types) are simultaneously more rigid (constrained
by similar types) and less tree-like (tight clustering degrades tree
topology).

Primary Hypothesis (H1)
Neighborhood density predicts Procrustes rigidity. Cell types with denser
local neighborhoods in expression space (lower mean distance to k nearest
neighbors) will tend to have smaller Procrustes residuals (higher rigidity).
Rationale: types in dense neighborhoods share expression programs with
nearby types, constraining their cross-species geometry.

Outcome measure: Spearman rho(neighborhood_density_k5, residual_magnitude).
Expected direction: POSITIVE (denser = lower mean distance; lower distance
correlates with lower residual = more rigid). Primary k=5.

Confirmation: rho > 0, p < 0.05, n=35.

Secondary Hypothesis (H2)
Neighborhood density predicts treeness. Cell types with denser
neighborhoods will have lower treeness scores (mean delta). Rationale:
tight clusters of similar types violate the four-point tree condition
because pairwise distances within the cluster are similar, driving
partition sums toward equality (delta toward 0).

Outcome measure: Spearman rho(neighborhood_density_k5, treeness_score).
Expected direction: NEGATIVE (denser neighborhood = lower treeness).

Confirmation: rho < 0, p < 0.05, n=35.

Mediation Hypothesis (H3)
Neighborhood density mediates the treeness-rigidity anticorrelation.
After controlling for neighborhood density via partial Spearman
correlation, the rigidity-treeness correlation will be substantially
attenuated. This would demonstrate that the anticorrelation is explained
by geometric proximity structure rather than a direct biological
relationship.

Outcome measure: Partial Spearman rho(residual_magnitude, treeness_score |
neighborhood_density_k5). Compare to raw rho=-0.349.

Mediation threshold: attenuation >= 50%. Attenuation defined as
(|raw rho| - |partial rho|) / |raw rho| * 100.

Confirmation: attenuation >= 50%.

Neighborhood Density Definition
For each cell type i, compute mean Euclidean distance to its k nearest
neighbors among the other 34 cell types in full 16,959-gene expression
space. Primary k=5. Sensitivity: k=3, k=10.

Input Data
- Human Tabula Sapiens centroids: output/phase2/scaled_35types/centroids_human_35.csv
  (35 cell types x 16,959 ortholog genes — same file used for treeness analysis)
- Procrustes residuals: output/phase2/scaled_35types/residuals_ranked.csv
- Treeness scores: output/liang_wagner/treeness_scores_per_celltype.csv

Falsification Conditions
Any ONE of the following kills the mediation hypothesis:

1. H1 fails: neighborhood density does NOT predict rigidity (p >= 0.05 or
   wrong sign). If density is unrelated to rigidity, it cannot mediate.

2. H2 fails: neighborhood density does NOT predict treeness (p >= 0.05 or
   wrong sign). If density is unrelated to treeness, it cannot mediate.

3. H3 fails: attenuation < 50%. Controlling for density does not
   substantially reduce the anticorrelation, meaning density is not the
   primary driver.

4. Sensitivity instability: k=3 and k=10 produce qualitatively different
   conclusions from k=5 (mediation confirmed at one k but not others).

Decision Table
| H1 | H2 | H3 (>=50%) | Interpretation |
|----|----|----|-------------|
| PASS | PASS | PASS | Mediation confirmed — anticorrelation is a geometric proximity artifact |
| PASS | PASS | FAIL | Density contributes but does not explain — anticorrelation has additional biological component |
| PASS | FAIL | — | Density predicts rigidity but not treeness — different mechanism |
| FAIL | — | — | Density does not predict rigidity — mediation hypothesis rejected |

Analysis will NOT be repeated with alternative density metrics, k values
beyond {3, 5, 10}, or post-hoc exclusions. The analysis is run once on
the pre-specified data and metrics.
