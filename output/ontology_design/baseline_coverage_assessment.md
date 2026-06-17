# Baseline Coverage Assessment — Current 35-Type Ontology

**Date:** 2026-03-16
**Purpose:** Document what the current 35 types cover, what they miss, and which are
mammal-specific candidates for expansion.

---

## 1. Current 35 Types by Organ System

### Immune (16 types — 46% of ontology)
| # | Cell Type | Rigidity Rank | Residual (% SSR) |
|---|-----------|---------------|-------------------|
| 1 | B cell | 17 | 2.59% |
| 2 | CD4-positive, alpha-beta T cell | 14 | 2.79% |
| 3 | CD8-positive, alpha-beta T cell | 35 (most rigid) | 0.77% |
| 4 | T cell | 7 | 3.64% |
| 5 | natural killer cell | 23 | 2.02% |
| 6 | mature NK T cell | 25 | 1.83% |
| 7 | macrophage | 16 | 2.61% |
| 8 | classical monocyte | 15 | 2.64% |
| 9 | intermediate monocyte | 24 | 1.90% |
| 10 | non-classical monocyte | 34 | 0.93% |
| 11 | monocyte | 22 | 2.03% |
| 12 | myeloid dendritic cell | 21 | 2.12% |
| 13 | myeloid leukocyte | 10 | 3.54% |
| 14 | neutrophil | 8 | 3.60% |
| 15 | granulocyte | 27 | 1.68% |
| 16 | plasma cell | 12 | 3.28% |

### Hematopoietic/Progenitor (4 types — 11%)
| # | Cell Type | Rigidity Rank | Residual (% SSR) |
|---|-----------|---------------|-------------------|
| 17 | hematopoietic precursor cell | 3 | 6.30% |
| 18 | hematopoietic stem cell | 4 | 5.98% |
| 19 | mesenchymal stem cell | 13 | 2.86% |
| 20 | mesenchymal stem cell of adipose tissue | 11 | 3.49% |

### Structural/Mesenchymal (5 types — 14%)
| # | Cell Type | Rigidity Rank | Residual (% SSR) |
|---|-----------|---------------|-------------------|
| 21 | stromal cell | 1 (most flexible) | 7.71% |
| 22 | fibroblast | 28 | 1.51% |
| 23 | fibroblast of cardiac tissue | 9 | 3.55% |
| 24 | adventitial cell | 26 | 1.75% |
| 25 | smooth muscle cell | 31 | 1.20% |

### Epithelial (7 types — 20%)
| # | Cell Type | Rigidity Rank | Residual (% SSR) |
|---|-----------|---------------|-------------------|
| 26 | epithelial cell | 2 | 7.08% |
| 27 | basal cell | 6 | 4.36% |
| 28 | luminal epithelial cell of mammary gland | 18 | 2.49% |
| 29 | large intestine goblet cell | 19 | 2.40% |
| 30 | enterocyte of epithelium of large intestine | 20 | 2.31% |
| 31 | bladder urothelial cell | 29 | 1.46% |
| 32 | pancreatic ductal cell | 30 | 1.23% |

### Metabolic/Hepatic (2 types — 6%)
| # | Cell Type | Rigidity Rank | Residual (% SSR) |
|---|-----------|---------------|-------------------|
| 33 | hepatocyte | 32 | 0.98% |
| 34 | pancreatic acinar cell | 5 | 4.41% |

### Vascular/Endothelial (1 type — 3%)
| # | Cell Type | Rigidity Rank | Residual (% SSR) |
|---|-----------|---------------|-------------------|
| 35 | endothelial cell | 33 | 0.95% |

---

## 2. Coverage Gaps

| Organ System | Current Count | Gap Severity | Notes |
|---|---|---|---|
| **Neural** | **0** | **CRITICAL** | Zero brain/CNS types. No neurons, glia, or neural support cells. Major organ system entirely unrepresented. Tabula Sapiens lacks brain tissue — requires supplementary human atlas (e.g., Allen Brain Cell Atlas). |
| **Endocrine** | **0** | **CRITICAL** | Zero endocrine types. Pancreatic beta/alpha/delta cells all below 500-cell gate in Tabula Sapiens (102/48/9 cells respectively). No thyroid, adrenal, or parathyroid representation. |
| **Cardiac/Muscle** | **0 dedicated** | **SEVERE** | Smooth muscle cell (structural) and cardiac fibroblast are present but no cardiomyocytes (dropped — absent from Tabula Sapiens). No skeletal muscle satellite cells. |
| **Respiratory** | **0** | **SEVERE** | No alveolar cells (type 1 or type 2), no airway epithelium (ciliated, club, goblet). AT2 cells abundant in TS (11,594) but only 292 in TMS. |
| **Renal** | **0** | **SEVERE** | No kidney-specific types. Proximal tubule, podocyte, collecting duct — all absent. |
| **Reproductive** | **0** | **MODERATE** | No reproductive-specific types. Sex-specific cells (Sertoli, granulosa) create asymmetric coverage across donors. Lower priority for cross-species framework. |
| **Dermal** | **0** | **MODERATE** | No keratinocytes, melanocytes, or skin-specific types despite skin being a major organ. |
| **Gastric** | **0** | **MODERATE** | No stomach-specific types (chief cells, parietal cells). |

### Summary
The current ontology covers **4 of 12 major organ systems** with meaningful depth
(immune, epithelial/GI, structural, metabolic/hepatic). The immune system alone
accounts for 46% of types — a massive overrepresentation driven by blood being the
most accessible tissue in both Tabula atlases. Neural, endocrine, cardiac, respiratory,
renal, reproductive, dermal, and gastric systems have zero dedicated representation.

---

## 3. Zebrafish-Absent Types (Mammal-Specific)

Seven of the current 35 types have **no zebrafish homolog** regardless of atlas
availability (confirmed in T3-F feasibility, DECISION-117):

| Cell Type | Reason for Absence | Expansion Candidate? |
|---|---|---|
| Luminal epithelial cell of mammary gland | Zebrafish lack mammary glands | Yes — mammalian universal |
| Bladder urothelial cell | Zebrafish lack urinary bladder | Yes — mammalian universal |
| Large intestine goblet cell | No large/small intestine distinction in teleosts | Yes — mammalian universal |
| Enterocyte of epithelium of large intestine | Same as above | Yes — mammalian universal |
| Mature NK T cell | NKT cells are a mammalian specialization | Yes — mammalian |
| Classical monocyte | Monocyte subtyping not applicable to teleosts | Yes — mammalian |
| Intermediate monocyte | Same as above | Yes — mammalian |
| Non-classical monocyte | Same as above | Yes — mammalian |

**Note:** The prompt specifies 7 types; actual count from T3-F report is 8 (CD4+ T and
CD8+ T cells are also listed as annotation-level mismatches, and 3 monocyte subtypes
are listed). The 7 most structurally absent (no homologous tissue/cell lineage) are
the top 7 rows. All 8 are retained in the 100-type ontology because they exist in
human, mouse, and macaque.

---

## 4. Progenitor/Stem Cell Types — Qualification Status

| Cell Type | Human Count | Mouse Count | Rigidity Rank | Qualification Status |
|---|---|---|---|---|
| Hematopoietic stem cell | 858 | 3,342 | 4 (flexible) | **PARTIALLY QUALIFIED.** DECISION-099: absent from NHPCA macaque atlas. T2-D progenitor analysis: p=0.038 but embryonic confound (rank jump 4→1 in TMS). Adult mouse HSC source still needed. |
| Hematopoietic precursor cell | 2,140 | 2,144 | 3 (flexible) | **QUALIFIED.** Passes all gates. Second-most flexible type. |
| Mesenchymal stem cell | 23,499 | 16,142 | 13 | **QUALIFIED.** Large populations in both species. Mid-range rigidity. |
| Mesenchymal stem cell of adipose tissue | 3,777 | 6,626 | 11 | **QUALIFIED.** Tissue-specific MSC subtype. NHPCA macaque: 22,259 cells (PASS). |
| Basal cell* | 33,526 | 3,726 | 6 (flexible) | **QUALIFIED as epithelial progenitor.** Basal cells serve as tissue-resident stem cells in skin, airway, and prostate. High flexibility consistent with progenitor behavior. |

*Basal cell has progenitor capacity but is primarily classified as epithelial in the
current ontology. If the user counts 5 progenitor types, basal cell is the likely
fifth candidate.

### Pattern
All 4 dedicated progenitor/stem types rank in the top 13 (flexible half) of the
rigidity spectrum. This is consistent with the attractor hypothesis: progenitor types
occupy shallow, broad attractors — high evolutionary flexibility.

---

## 5. Macaque Compatibility (NHPCA Atlas, DECISION-099)

| Status | Count | Types |
|---|---|---|
| PASS (≥200 cells in NHPCA) | 27 | B cell, CD4+ T, CD8+ T, T cell, NK, NKT, macrophage, monocyte, myeloid DC, neutrophil, granulocyte, plasma cell, endothelial, hepatocyte, epithelial, basal, stromal, fibroblast, smooth muscle, adventitial, MSC, MSC adipose, bladder urothelial, enterocyte, pancreatic acinar, pancreatic ductal, hematopoietic precursor |
| ABSENT from NHPCA | 8 | Classical monocyte, intermediate monocyte, non-classical monocyte, myeloid leukocyte, fibroblast of cardiac tissue, HSC, large intestine goblet, mammary luminal |

**Key concern:** NHPCA uses 93% snRNA-seq with median 1,324 genes/cell — below the
2,000-gene quality gate. Gene detection may be insufficient for Procrustes analysis
(cf. MCA failure at ~700 genes). Qu et al. atlas uses 10x Chromium (better technology)
but requires 14.3 GB download and metadata parsing.
