# Analysis Plan: DILIrank Hepatocyte Rigidity Validation
## LINCS/L1000 Geometric Deformation Analysis

*This is an internal analysis plan documenting the intended approach; it was not registered with an external preregistration repository (consistent with the manuscript's disclosure).*

---

### Section 1 — Hypothesis Statement

Evolutionarily rigid cell types occupy geometrically constrained
regions of transcriptomic space. This constraint implies a lack of
geometric slack: the identity geometry of rigid cell types resists
transcriptional deformation under pharmacological perturbation until
a threshold is reached, at which point identity geometry collapses
catastrophically — shattering rather than bending.

The specific prediction is a non-linear threshold effect, not a
linear dose-response correlation. Most-Concern DILI drugs (DILIrank
v2) will be enriched in the UPPER TAIL of the relative Procrustes
deformation distribution computed from LINCS L1000 HepG2 signatures.
That is: drugs that cause severe liver injury do not gradually
increase hepatocyte geometric displacement across DILI severity
classes. Instead, most drugs produce modest displacement (the rigid
geometry absorbs the perturbation), but drugs classified as
Most-Concern are disproportionately represented among those that
produce catastrophic geometric displacement — the identity geometry
shatters.

This shattering claim is falsified if Test 2 (Fisher's exact test,
upper-tail enrichment) fails, regardless of Test 1 outcome. A
passing Test 1 (Mann-Whitney shift) without a passing Test 2
indicates that Most-Concern drugs produce greater average
displacement, but the displacement is distributed across the
severity range rather than concentrated in the tail. This would
support a toxicity enrichment claim but not the shattering mechanism.

---

### Section 2 — Biological Rationale

**Why hepatocyte rigidity predicts DILI.** Hepatocyte is among the
most evolutionarily rigid cell types in the CellWarp framework,
contributing only 1.9% of the sum of squared residuals in the
original 6-type Procrustes analysis. This rigidity is confirmed
across three independent datasets (Tabula, Sun2023, PanSci) and
is robust to dimensionality reduction, subsampling, and all ten
mechanistic confound tests. Drug-induced liver injury (DILI) is the
single most documented drug toxicity endpoint, with standardized
severity classifications (DILIrank). The cell-type-to-organ mapping
is clean and unambiguous: hepatocytes are the primary parenchymal
cell of the liver, and DILI is defined by hepatocyte damage. No
intermediate tissue-level inference is required.

**Why Procrustes deformation is orthogonal to standard differential
expression.** A drug that displaces hepatocyte identity geometry —
even with modest individual gene expression changes — is doing
something categorically different from a drug that activates a
stress response pathway without displacing identity. Procrustes
distance captures geometric displacement of the entire identity
configuration; differential expression captures per-gene activation
state. These are orthogonal signals, as demonstrated by the
CellWarp three-way framework: cancer (identity transformation)
and COVID-19 (activation state change) produce different geometric
signatures despite comparable numbers of differentially expressed
genes (DECISION-084). Procrustes distance is computed over the
entire 978-gene L1000 landmark set, not a subset. This is a
deliberate design choice to measure whole-identity geometry
displacement rather than pathway-specific or cytoskeletal changes.
Cytotoxicity-driven structural collapse (e.g. cytoskeletal gene
changes from cell death — Lmod1 class) would be distributed across
the landmark set and not produce a coherent directional Procrustes
displacement. Coherent directional displacement is the signal;
diffuse collapse is noise.

**HepG2 limitation.** HepG2 cells are an immortalized hepatocellular
carcinoma line that lacks full CYP450 metabolic enzyme activity
compared to primary hepatocytes. Drugs requiring metabolic activation
(pro-drugs, reactive metabolite formation) to become hepatotoxic may
show zero geometric deformation in HepG2 despite carrying a
Most-Concern DILI classification — a false negative, not a false
positive. Since the analysis is optimized against false positives
(the shattering claim requires upper-tail enrichment, not universal
detection), this limitation strengthens the integrity of any signal
observed: drugs that DO produce geometric displacement in
metabolically limited HepG2 cells are doing so through direct
mechanisms, not metabolic byproducts. Pre-specified mitigation:
CYP450-dependent drug stratification using DrugBank curated
substrate annotations (Section 4, Step 4). Primary hepatocyte
validation would be expected to strengthen, not weaken, the result
— an explicit invitation for wet-lab follow-up.

---

### Section 3 — Data Sources

1. **DILIrank v2** — FDA drug-induced liver injury severity
   classification. ~192 Most-Concern drugs, plus Less-Concern,
   Ambiguous-DILI-Concern, and No-DILI-Concern categories. Chen M,
   Suzuki A, Thakkar S, Yu K, Hu C, Tong W. "DILIrank: the largest
   reference drug list ranked by the risk for developing drug-induced
   liver injury in humans." Drug Discovery Today. 2016;21(4):648-653.
   v2 update via FDA NCTR. Download source: FDA National Center for
   Toxicological Research (NCTR) DILIrank dataset.

2. **LINCS L1000 Phase II** — Library of Integrated Network-Based
   Cellular Signatures. GEO accession GSE70138 and/or CLUE.io Level
   5 signatures. HepG2 cell line signatures only. 978 landmark genes.
   Access date: [DATE TO BE FILLED ON ACCESS].

3. **CellWarp hepatocyte rigidity score** — Procrustes residual from
   current pipeline output. Three-dataset confirmed (Tabula, Sun2023,
   PanSci). Commit hash: [SPECIFY COMMIT HASH ON EXECUTION].

4. **DrugBank CYP450 substrate annotations** — Curated substrate and
   inhibitor classifications for cytochrome P450 enzymes (CYP1A2,
   CYP2C9, CYP2C19, CYP2D6, CYP3A4). Version: current at time of
   access. DrugBank chosen over ChEMBL for curation quality and
   minimization of researcher degrees of freedom during annotation.
   DrugBank provides binary curated classifications; ChEMBL would
   require threshold decisions on binding affinity cutoffs, introducing
   post-hoc flexibility.

---

### Section 4 — Analysis Pipeline

Steps are executed in strict order. No steps may be added after data
access.

**Step 1 — SENSITIVITY GATE (fractal geometry test).** Execute before
any DILI data is accessed. Compute hepatocyte rigidity ranking using
only the 978 L1000 landmark genes (intersected with ortholog space).
Compute Spearman rho against full-space rigidity ranking across all
cell types with confirmed rigidity scores.
- THRESHOLD: rho >= 0.6.
- PASS: Fractal geometry confirmed — identity geometry is preserved
  at 5.8% of gene space (978/16,896 shared orthologs). Report rho
  value in paper as a named finding ("fractal geometry of cell
  identity").
- FAIL: L1000 platform has landmark gene bias that obscures rigid
  cell type geometry. Abort DILI analysis entirely. Document as a
  platform characterization finding (the L1000 landmark set does not
  capture identity geometry). Do not proceed to Step 2.

**Step 2 — DATA LINKAGE.** Match DILIrank v2 drugs to LINCS L1000
compounds by InChIKey (preferred) or canonical name with manual
verification of ambiguous matches. Report: total DILIrank drugs,
total L1000 compounds in HepG2, matched count by DILI category,
match rate.
- THRESHOLD: >= 25% of Most-Concern drugs matched AND minimum n=50
  Most-Concern drugs in absolute terms.
- Below threshold: Insufficient coverage for landmark claim. Abort
  analysis and document coverage gap.

**Step 3 — RELATIVE DEFORMATION COMPUTATION.** For each matched drug
with HepG2 L1000 signatures, compute Procrustes distance from
treated HepG2 signature to untreated HepG2 centroid in 978-gene
landmark space. Use Level 5 (replicate-collapsed) signatures where
available.
- PRIMARY concentration: 10 uM (standard L1000 screening dose).
- CONSISTENCY CHECK: 1 uM.
- DISCORDANCE RULE (pre-specified): If 10 uM and 1 uM results are
  discordant (defined as: drug falls in opposite halves of the
  deformation distribution at the two concentrations), use 10 uM
  result, flag drug in supplementary data. No post-hoc concentration
  selection permitted.

**Step 4 — CYP450 STRATIFICATION.** Using DrugBank curated substrate
annotations, classify each Most-Concern drug as CYP450-dependent
(substrate of any major CYP450 enzyme: CYP1A2, CYP2C9, CYP2C19,
CYP2D6, CYP3A4) or CYP450-independent. All subsequent tests run
twice:
- (A) Full Most-Concern set.
- (B) CYP450-excluded Most-Concern set (removing CYP450-dependent
  drugs, which may show false-negative deformation in HepG2).
Both are pre-specified analyses. Neither is primary over the other —
both reported side by side.

**Step 5 — FOUR STATISTICAL TESTS IN ORDER.** No additional tests
permitted after data access. All tests run on both (A) full and (B)
CYP450-excluded sets.

**Test 1 (Primary hypothesis): Mann-Whitney U, one-tailed.**
H0: Most-Concern and No-DILI-Concern drugs have the same
distribution of relative Procrustes deformation distances.
H1: Most-Concern drugs have larger deformation distances.
Exclude Ambiguous-DILI-Concern drugs from this test (they are
neither clean positives nor clean negatives). alpha = 0.05.

**Test 2 (Shattering model): Fisher's exact test, one-tailed.**
Contingency table: Most-Concern vs No-DILI-Concern drugs x top
quartile vs bottom three quartiles of deformation distribution.
H0: Most-Concern drugs are not enriched in the top quartile.
H1: Most-Concern drugs are enriched in the top quartile
(upper-tail concentration). alpha = 0.05.
This test is required to support the shattering mechanism claim.
Test 1 passing without Test 2 = enrichment real, shattering
mechanism not supported.

**Test 3 (Bimodality): Hartigan's dip test.**
Applied within Most-Concern drugs only. Tests whether the
deformation distribution among Most-Concern drugs is bimodal.
- Bimodal (p < 0.05): Strongest shattering evidence — two
  populations of Most-Concern drugs, those absorbed by rigid
  geometry and those that shatter it.
- Unimodal (p >= 0.05): Threshold framing weakened but not fatal
  to Test 1/2 results.

**Test 4 (Covariate falsification): Partial correlation.**
Partial Spearman correlation between Procrustes deformation distance
and DILI severity class (ordinal: No-Concern=0, Less-Concern=1,
Most-Concern=2; Ambiguous excluded), controlling for mean baseline
expression level of top 500 most variable genes in untreated HepG2.
Signal must survive this control — i.e., partial correlation must
remain significant (p < 0.05) and in the same direction as the
zero-order correlation.
- PASS: Geometric displacement predicts DILI severity independent
  of baseline expression magnitude.
- FAIL: Signal is driven by parts list (baseline expression level),
  not geometric architecture. Hard abort (see Section 5).

---

### Section 5 — Falsification Conditions

Three hard aborts, one conditional weakening.

**HARD ABORT 1:** Step 1 sensitivity gate fails (rho < 0.6).
L1000 landmark gene space does not preserve identity geometry.
Analysis is invalid at the platform level. Do not access DILI data.
Document as a platform characterization finding.

**HARD ABORT 2:** Step 2 match rate insufficient (< 25% of
Most-Concern drugs matched OR n < 50 Most-Concern drugs in absolute
terms). Analysis is underpowered for a landmark claim. Document
coverage gap. Do not proceed to statistical tests.

**HARD ABORT 3:** Test 4 covariate falsification fails. Procrustes
deformation—DILI association is driven by baseline expression
magnitude (parts list), not geometric architecture. The Toxicity Map
claim is not viable because the signal reduces to differential
expression by another name. Landmark version is not viable. Paper
reverts to genuinely satisfying version without DILI claim. Advisor
review required before any further Landmark track work.

**CONDITIONAL WEAKENING:** Test 1 passes, Test 2 fails. Most-Concern
drugs produce greater average geometric displacement (enrichment is
real), but the displacement is distributed across the severity range
rather than concentrated in the upper tail. Shattering mechanism
claim is retired. Toxicity Map claim survives in weakened linear
form: "rigid cell types show graded geometric displacement under
pharmacological perturbation that correlates with clinical toxicity
severity." Advisor reviews framing before paper writing proceeds.

---

### Section 6 — Reporting Commitment

All four statistical test results are reported in the paper
regardless of direction or significance. The sensitivity gate result
(Step 1) is reported regardless of outcome — fractal geometry
confirmation or platform limitation, both are findings worth
documenting. CYP450 stratification results (both full and excluded
sets) are reported regardless of direction or whether stratification
changes the outcome. No selective reporting under any circumstances.

If any Hard Abort condition is triggered, the abort and its specific
reason are reported in the methods section of whatever paper version
follows. Pre-registration aborts are scientific results, not failures
to suppress.

---

### Section 7 — Replication Design (Cardiotoxicity)

**PLACEHOLDER ONLY.** Do not activate until cardiomyocyte rigidity is
confirmed in the 100-type ontology expansion (DECISION-124, shelf
lifted after macaque pipeline completes).

At that point a separate pre-registration document will be filed with
the following planned design:
- **Endpoint:** CredibleMeds QTc prolongation annotations.
- **Cell line:** iPSC-derived cardiomyocyte L1000 signatures (if
  available in LINCS).
- **Statistical design:** Identical structure to hepatocyte analysis
  (sensitivity gate, data linkage, deformation computation, CYP450
  stratification, four statistical tests, same falsification
  conditions).
- **Cell-type-to-organ mapping:** Cardiomyocyte -> cardiotoxicity
  (QTc prolongation, cardiac arrest). Clean mapping, no intermediate
  tissue inference.

---

### Section 8 — Registration Metadata

- **Date:** 2026-03-16
- **Registered by:** [AUTHOR]
- **Analysis to be executed by:** [AUTHOR]
- **Advisor review:** Completed 2026-03-16
- **Data access begins:** AFTER this document is committed to the
  repository
- **Repository commit:** This document MUST be committed before any
  LINCS or DILIrank data is downloaded or examined. The commit hash
  serves as the timestamp of pre-registration.
