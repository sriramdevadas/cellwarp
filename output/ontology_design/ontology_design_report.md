# 100-Type Ontology Design Report

**Date:** 2026-03-16
**Status:** Design exercise, for review. No expression data downloaded.
**Purpose:** Expand from 35 to 100 cell types with principled selection criteria.
**Scope:** Paper 3 / resource paper foundation; informs macaque scout compatibility.

---

## 1. Baseline Coverage Gaps

The current 35-type ontology covers **4 of 12 major organ systems** with meaningful
depth. The immune system alone accounts for 46% of types (16/35), driven by blood
being the most accessible tissue in both Tabula atlases.

| Organ System | Current Types | Gap Severity |
|---|---|---|
| Immune | 16 | OVERREPRESENTED |
| Epithelial | 7 | Moderate coverage |
| Structural/Mesenchymal | 5 | Moderate coverage |
| Hematopoietic/Progenitor | 4 | Adequate |
| Metabolic/Hepatic | 2 | Thin |
| Vascular/Endothelial | 1 | Thin |
| **Neural** | **0** | **CRITICAL GAP** |
| **Endocrine** | **0** | **CRITICAL GAP** |
| **Cardiac/Muscle** | **0** | **CRITICAL GAP** |
| **Respiratory** | **0** | **SEVERE GAP** |
| **Renal** | **0** | **SEVERE GAP** |
| **Dermal** | **0** | **MODERATE GAP** |

Seven types are mammal-specific (absent in zebrafish) and confirmed for human-mouse-
macaque ontology. Four progenitor/stem types cluster in the flexible end of the
rigidity spectrum, consistent with the attractor hypothesis.

Full assessment: `output/ontology_design/baseline_coverage_assessment.md`

---

## 2. Candidate Pool

### Source A — CellMarker 2.0
Estimated ~180 human cell types with ≥5 wet-lab validated markers. This is the quality
filter ensuring each type has sufficient literature backing for identity gene validation
(same standard as DECISION-114 CellMarker validation).

### Source B — Tabula Sapiens (≥500 cells)
36 types confirmed from cell_type_inventory.csv. Tabula Sapiens covers 24 tissues but
**lacks brain, thyroid, adrenal, stomach, and parathyroid tissue**. This is the primary
bottleneck for organ system expansion.

### Source C — Tabula Muris Senis (≥500 cells)
35 types confirmed. TMS covers 23 tissues including brain (advantage over TS for neural
types) but has lower cell counts for many types.

### Strict Intersection (TS ≥500 ∩ TMS ≥500)
**35 types** — exactly the current ontology. This is the ceiling with single-atlas
Tabula queries.

### Census-Pooled Estimated Intersection
**~80-120 types** when using full CELLxGENE Census (pools hundreds of datasets per
species). Requires empirical Census metadata queries to confirm. Neural types require
Allen Brain Cell Atlas as supplementary source.

### Candidate Pool Enumerated
105 types enumerated across all sources (35 current + 12 near-miss + 58 Census/
supplementary). Full list: `output/ontology_design/candidate_pool.csv`

---

## 3. Proposed 100-Type Ontology

### Organ System Balance

| Organ System | Current | Proposed | Change | % of 100 |
|---|---|---|---|---|
| Immune | 16 | 20 | +4 | 20% |
| Epithelial | 7 | 15 | +8 | 15% |
| Structural/Mesenchymal | 5 | 12 | +7 | 12% |
| Neural | 0 | 12 | +12 | 12% |
| Vascular/Endothelial | 1 | 8 | +7 | 8% |
| Metabolic/Hepatic | 2 | 8 | +6 | 8% |
| Hematopoietic/Progenitor | 4 | 10 | +6 | 10% |
| Endocrine | 0 | 8 | +8 | 8% |
| Other | 0 | 7 | +7 | 7% |
| **Total** | **35** | **100** | **+65** | **100%** |

The immune system drops from 46% to 20% of the ontology. Neural, endocrine, and
vascular systems gain meaningful representation. All major organ systems are covered.

### Selection Criteria Applied (in priority order)

1. **Current 35 retained** — all validated, foundational (Criterion 1)
2. **Organ system balance** — 12 neural, 8 epithelial, 8 endocrine, 7 vascular, 7
   structural, 6 metabolic, 6 hematopoietic, 4 immune, 7 other new types (Criterion 2)
3. **Biological diversity** — new types span predicted rigidity range: neural types
   (predicted rigid: highly conserved across mammals), endocrine types (predicted
   flexible: rapidly evolving hormone systems), structural types (mixed) (Criterion 3)
4. **Macaque compatibility** — flagged per type (see breakdown below) (Criterion 4)
5. **Disease relevance** — prioritized types with known disease roles (Criterion 5)
6. **Annotation quality** — avoided catch-all labels; flagged existing catch-alls
   (T cell, epithelial cell, myeloid leukocyte, granulocyte) (Criterion 6)

### Full Proposed Ontology

See `output/ontology_design/proposed_100_type_ontology.csv` for the complete table with:
- CellWarp standardized name
- CellMarker 2.0 top 3 markers
- Organ system assignment
- Tabula Sapiens cell count (confirmed or estimated)
- Tabula Muris Senis cell count (confirmed or estimated)
- Macaque compatibility flag
- Data source notes
- Priority tier
- Annotation challenge notes

### Highlighted New Types by Disease Relevance

| Disease Area | New Types Added | Count |
|---|---|---|
| **Neurodegenerative** | astrocyte, microglial cell, OPC, dopaminergic neuron | 4 |
| **Cancer** | keratinocyte (melanoma), cholangiocyte (CCA), cardiac myocyte (cardiac), melanocyte (melanoma) | 4 |
| **Autoimmune** | regulatory T cell, follicular helper T cell, plasmacytoid DC | 3 |
| **Metabolic** | pancreatic beta cell (diabetes), adipocyte (obesity), hepatic stellate cell (fibrosis) | 3 |
| **Respiratory** | alveolar type 2 (COVID/ARDS), alveolar macrophage (COVID), ciliated epithelial | 3 |
| **Renal** | proximal tubule, podocyte | 2 |
| **Cardiovascular** | pericyte, arterial endothelial (atherosclerosis) | 2 |

---

## 4. Macaque Compatibility Breakdown

| Flag | Count | Description |
|---|---|---|
| **CONFIRMED** | 38 | Known to exist in NHPCA or Allen Brain macaque data. Includes all neural types via Allen Brain Cell Atlas macaque coverage, plus 27/35 current types confirmed in NHPCA. |
| **LIKELY** | 46 | Mammalian-universal type expected to exist in macaque. No annotation confirmation yet. Includes all common immune subtypes, epithelial types, structural types. |
| **UNCERTAIN** | 8 | May be absent, differently annotated, or too rare. Includes: monocyte subtypes (classical/intermediate/non-classical — NHPCA annotation gap), hepatic stellate cell, parathyroid cell, pancreatic delta cell, chromaffin cell, Langerhans cell. |
| **ABSENT_NHPCA** | 8 | Confirmed absent from NHPCA atlas: classical monocyte, intermediate monocyte, non-classical monocyte, myeloid leukocyte, cardiac fibroblast, HSC, large intestine goblet, mammary luminal. These are annotation-resolution gaps, not biological absences — the cell types exist in macaque but NHPCA uses coarser labels. |

### Macaque Atlas Strategy

The NHPCA (Han et al. 2022) provides 1.14M cells across 45 tissues but has two
limitations:
1. **Gene detection:** 93% snRNA-seq, median 1,324 genes/cell (below 2,000 gate)
2. **Annotation resolution:** Coarse labels miss monocyte subtypes, some tissue-
   specific types

For the 100-type ontology, the recommended macaque strategy is:
- **NHPCA** for cell type coverage verification (~65 types likely PASS)
- **Allen Brain Cell Atlas (macaque)** for all 12 neural types + brain endothelial
- **Qu et al. 2022 atlas** (10x Chromium, 174K cells) as technology-quality supplement
  for types where NHPCA gene detection is insufficient
- **Gene detection empirical test:** Download a subset (e.g., hepatocyte + CD8+ T)
  from NHPCA and run Procrustes against human/mouse. If p<0.01, gene detection is
  sufficient despite being below nominal gate. If p>0.05, NHPCA is not usable and
  Qu et al. becomes primary.

---

## 5. Feasibility Assessment

| Dimension | Status | Detail |
|---|---|---|
| Gene space | **PASS** | Ortholog space is species-level; adding types doesn't change it |
| Compute time | **PASS** | ~25 min for full pipeline at n=100 (~2.5× current) |
| Memory | **PASS** | ~16 GB peak, within M-series MacBook capacity |
| Data availability | **CONDITIONAL** | 35 types confirmed; 60 more require Census queries + Allen Brain |
| High-risk types | **5 flagged** | AT1 pneumocyte, hepatic stellate, delta cell, parathyroid, enteroendocrine |
| Realistic ceiling | **93-95 types** | After removing high-risk types |

### External Dependencies

| Dependency | Types Resolved | Status |
|---|---|---|
| Allen Brain Cell Atlas (human) | 13 (all neural + brain endothelial) | Available in Census |
| CELLxGENE Census pooling | ~30 (near-miss + cross-dataset) | Requires metadata query |
| Human Pancreas Analysis Program | 2 (beta + alpha cell) | Available via GEO/Census |
| Census gastric datasets | 2 (chief + parietal cell) | Limited availability |
| Census adrenal/thyroid datasets | 3 (adrenocortical, chromaffin, thyroid follicular) | Limited availability |

Full feasibility analysis: `output/ontology_design/feasibility_check.md`

---

## 6. Recommended Sequencing — Implementation Priority

### Phase A: High-Priority Expansion (35 → 70 types, +35 new)

Add types with highest confidence of data availability AND macaque compatibility AND
disease relevance. These can proceed with Census queries immediately.

| Priority | Types | Rationale |
|---|---|---|
| **A1: Neural core (6)** | oligodendrocyte, astrocyte, microglial cell, OPC, excitatory neuron, inhibitory neuron | Fills largest gap. Allen Brain confirmed for human + macaque. Neurodegenerative disease axis. |
| **A2: Immune completion (4)** | regulatory T cell, mast cell, basophil, plasmacytoid DC | Autoimmune disease axis. Census pooling straightforward. |
| **A3: Vascular expansion (4)** | vein endothelial, arterial endothelial, capillary endothelial, lymphatic endothelial | Endothelial subtypes test tissue-context rigidity. Near-miss data available. |
| **A4: Epithelial organs (6)** | alveolar type 2, keratinocyte, ciliated epithelial, proximal tubule, club cell, podocyte | Respiratory, renal, dermal coverage. COVID-relevant. |
| **A5: Metabolic expansion (5)** | cholangiocyte, Kupffer cell, enterocyte (small intestine), pancreatic beta cell, pancreatic alpha cell | GI + hepatic + endocrine coverage. Diabetes axis. |
| **A6: Structural (4)** | pericyte, myofibroblast, adipocyte, mesothelial cell | Fibrosis, obesity axes. Predicted flexible types. |
| **A7: Hematopoietic (4)** | erythrocyte, megakaryocyte, thymocyte, erythroid progenitor | Blood lineage completion. Near-miss data. |
| **A8: Neural subtypes (2)** | Schwann cell, ependymal cell | Peripheral + ventricular glia. |

### Phase B: Lower-Priority Expansion (70 → 90 types, +20 new)

Types requiring more extensive Census searches or with lower disease relevance.

| Priority | Types | Rationale |
|---|---|---|
| **B1: Neural subtypes (4)** | Purkinje, medium spiny neuron, dopaminergic, cerebellar granule | Specific neuronal classes. Allen Brain data. |
| **B2: Vascular subtypes (3)** | brain endothelial, sinusoidal endothelial, endocardial | Tissue-specific endothelium. |
| **B3: Structural (3)** | osteoblast, satellite cell, chondrocyte | Musculoskeletal coverage. May be rare. |
| **B4: Metabolic/GI (3)** | gastric chief cell, parietal cell, hepatic stellate cell | Gastric + hepatic coverage. Limited data. |
| **B5: Endocrine (3)** | thyroid follicular, adrenocortical, chromaffin | Endocrine coverage. Census-dependent. |
| **B6: Hematopoietic (2)** | common myeloid progenitor, pro-B cell | Progenitor lineage. Annotation-dependent. |
| **B7: Other (2)** | melanocyte, goblet cell (airway) | Skin + respiratory. |

### Phase C: Provisional Types (90 → 100, +10 aspirational)

Types at highest risk of data unavailability. Include in design but prepared to replace.

| Priority | Types | Risk |
|---|---|---|
| **C1** | type I pneumocyte | HIGH — fragile, poor capture in scRNA-seq |
| **C2** | Paneth cell | MEDIUM — may be borderline in mouse |
| **C3** | pancreatic delta cell | HIGH — rare in both species |
| **C4** | enteroendocrine cell | HIGH — rare in both species |
| **C5** | parathyroid cell | HIGH — no standard atlas tissue |
| **C6** | Langerhans cell | HIGH — both species very low counts |
| **C7** | cardiac myocyte | MEDIUM — mouse near-miss |
| **C8** | alveolar macrophage | LOW — annotation overlap with macrophage |
| **C9** | follicular helper T cell | LOW — annotation overlap with CD4+ T |
| **C10** | dendritic cell type 2 | LOW — annotation overlap with myeloid DC |

**Replacement candidates** from the 105-type candidate pool if provisional types fail:
tuft cell, intercalated cell (kidney), collecting duct cell, bronchial smooth muscle
cell, tendon cell.

---

## 7. Key Risks and Mitigations

### Risk 1: Annotation Inconsistency
**Problem:** Cell type labels vary across atlases (e.g., "macrophage" vs "Kupffer cell",
"CD4 T cell" vs "helper T cell"). Census pooling amplifies this.
**Mitigation:** Standardize to Cell Ontology (CL) identifiers. Build a mapping table
from atlas-specific labels to CellWarp ontology names before data download. Reject
types where mapping is ambiguous.

### Risk 2: Endocrine Ceiling
**Problem:** 5/8 endocrine types have very low cell counts (<500 in Tabula). Even Census
pooling may not resolve delta cell, enteroendocrine, or parathyroid.
**Mitigation:** Designate these as provisional. If endocrine representation drops below
5 types, redirect slots to better-powered organ systems (e.g., additional neural
subtypes or tissue-specific immune types).

### Risk 3: Batch Effects from Census Pooling
**Problem:** Pooling cells from multiple labs introduces batch effects that Procrustes
may interpret as biological signal (same failure mode as MCA microwell-seq).
**Mitigation:** Apply the same negative control framework used for replication datasets.
For each Census-pooled cell type, compute within-dataset centroids and verify
convergence before pooling. Reject types where inter-dataset variance exceeds
inter-species variance.

### Risk 4: NHPCA Gene Detection for Macaque
**Problem:** NHPCA median 1,324 genes/cell may be insufficient for Procrustes (MCA
failed at ~700 genes/cell).
**Mitigation:** Empirical test on 2-3 cell types before full expansion. If NHPCA fails,
pivot to Qu et al. 2022 atlas (10x Chromium technology) as primary macaque source,
accepting reduced tissue coverage (16 tissues vs 45).

### Risk 5: Neural Atlas Integration
**Problem:** Allen Brain Cell Atlas uses different preprocessing than Tabula series.
Brain cells may have different QC characteristics.
**Mitigation:** Apply CellWarp standard QC pipeline (src/qc.py) uniformly to all
sources. Verify that Allen Brain centroids are not systematically offset from Tabula
centroids for shared types (e.g., endothelial cells appear in both brain and
non-brain atlases — use as internal consistency check).

---

## 8. Open Decision Points

This document is for review. The following decisions are requested:

1. **Census pooling approval:** The strict Tabula-only intersection yields 35 types
   (current ontology). Reaching 100 requires pooling across Census datasets. Is this
   methodologically acceptable, or does it introduce too much batch effect risk?

2. **Allen Brain dependency:** 12 neural types require Allen Brain Cell Atlas as a
   supplementary human source. This adds a major new data dependency. Is this
   acceptable, or should the neural count be reduced?

3. **Endocrine realistic targets:** Only 3/8 proposed endocrine types have high
   confidence of reaching ≥500 cells (beta cell, alpha cell, thyroid follicular).
   Should we aim for 5 endocrine types and redistribute 3 slots, or keep 8 as
   aspirational?

4. **Provisional types:** 5-7 types are flagged as HIGH RISK. Should these be included
   in the design (with replacement plan) or excluded from the initial target?

5. **Macaque compatibility threshold:** 8 types are ABSENT from NHPCA and 8 more are
   UNCERTAIN. Should the ontology require all 100 types to have macaque compatibility,
   or is ~85% coverage acceptable for Paper 2?

6. **Catch-all labels:** The current 35 includes 4 generic catch-all types (T cell,
   epithelial cell, myeloid leukocyte, granulocyte) with known annotation ambiguity.
   Should these be retained for continuity with the 35-type results, or replaced with
   more specific types?

---

## Appendix: Outputs Generated

| File | Description |
|---|---|
| `output/ontology_design/baseline_coverage_assessment.md` | Step 0: Current 35-type organ system mapping, gaps, zebrafish-absent types, progenitor status |
| `output/ontology_design/candidate_pool.csv` | Step 1: 105 candidates from CellMarker × TS × TMS intersection |
| `output/ontology_design/proposed_100_type_ontology.csv` | Step 2: Full 100-type table with markers, counts, macaque flags, notes |
| `output/ontology_design/feasibility_check.md` | Step 3: Gene space, compute scaling, atlas gap analysis |
| `output/ontology_design/ontology_design_report.md` | Step 4: This summary document |
