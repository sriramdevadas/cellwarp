# Sun et al. 2023 — Pre-Download Metadata Validation

**Date:** 2026-03-15
**Paper:** Sun S et al. "A single-cell transcriptomic atlas of exercise-induced
anti-inflammatory and geroprotective effects across the body"
**Journal:** Innovation (Cambridge) 4(1):100380, 2023
**DOI:** 10.1016/j.xinn.2023.100380 | PMID: 36747595 | PMC: PMC9898793
**Data:** OMIX002605 (processed), GSA CRA007207 (raw)
**GitHub:** github.com/wxb1998/Mouse-exercise-Project

---

## Question 1: Sedentary Control Separability

### VERDICT: CLEAN

Sedentary control cells are clearly labeled and separable at both file and
metadata levels.

**Experimental design (2×2 factorial):**
- YC = Young Control (sedentary, 2-month-old male C57BL/6J)
- YE = Young Exercise (2-month-old, 12 months voluntary wheel running)
- OC = Old Control (sedentary, 16-month-old)
- OE = Old Exercise (16-month-old, 12 months voluntary wheel running)
- LPS subgroups: separate files with "LPS" in filename

**OMIX002605 file naming:** `{Condition}-{Tissue}` (e.g., `YC-Liver`, `OC-Kidney`).
84 total files. No replicate numbers in processed files (3 biological replicates
per condition-tissue merged during processing).

**Metadata columns in processed data:**
- `group`: YC / YE / OC / OE (condition only)
- `tissue`: Liver / Kidney / Lung / etc. (tissue only)
- `batch`: YC_Liver / OE_Kidney / etc. (combined)
- `leiden_anno`: Cell type annotation (101 types)

**Extraction method:** Download only files with `YC-` or `OC-` prefix, excluding
any with `LPS`. Or filter integrated object: `adata.obs['group'].isin(['YC', 'OC'])`.

**Caveats:**
- Per-mouse resolution lost in processed files (3 replicates merged per condition-tissue)
- LPS samples are separate files — clean exclusion possible
- scRNA-seq (9 tissues) and snRNA-seq (5 tissues) run through separate pipelines

---

## Question 2: Cell Type Coverage Against Our 35-Type Ontology

### VERDICT: 14 PASS / 8 BORDERLINE / 13 ABSENT — PASSES ≥12 THRESHOLD (MARGINAL)

**Tissues profiled:**
- scRNA-seq (9): lung, aorta, kidney, liver, small intestine, testis, spleen, bone marrow, peripheral blood
- snRNA-seq (5): brain, cerebellum, spinal cord, heart, skeletal muscle

**Cell types identified:** 101 main types (63 scRNA + 38 snRNA), 305 clusters

**Sedentary control cell estimate:** ~127K cells (YC+OC combined) or ~63K per
age-matched arm, across 14 tissues.

### Mapping Table

| # | Our Type | Sun Match | Label(s) | Tissue | Proto | ≥500? | Status |
|---|----------|----------|----------|--------|-------|-------|--------|
| 1 | B cell | YES | BC, ProBC, LProBC | spleen, BM, blood, liver | 10x | Yes | **PASS** |
| 2 | CD4+ T cell | YES | CD4_Naive, CD4_Mem | spleen, blood, BM | 10x | Yes | **PASS** |
| 3 | CD8+ T cell | YES | CD8_Naive, CD8_Mem, CD8_CTL | spleen, blood, BM | 10x | Yes | **PASS** |
| 4 | T cell | YES | TC, ProTC, CD4+CD8+TC, CD4-CD8-TC | multiple | both | Yes | **PASS** |
| 5 | adventitial cell | NO | — | — | — | — | **ABSENT** |
| 6 | basal cell | NO | — | — | — | — | **ABSENT** |
| 7 | bladder urothelial | NO | — | no bladder | — | — | **ABSENT** |
| 8 | classical monocyte | PARTIAL | Mono, Mono_BM | blood, BM | 10x | Monocytes not subtyped | **BORDERLINE** |
| 9 | endothelial cell | YES | EC_Liver, EC_Lung, EC_Aorta, EC_Kidney, EC (sn) | all | both | Yes (7 variants) | **PASS** |
| 10 | enterocyte (large intestine) | PARTIAL | Epi_Intestine | small intestine only | 10x | Wrong segment | **BORDERLINE** |
| 11 | epithelial cell | YES | Epi_Lung, AT1, AT2 | lung, intestine | 10x | Yes | **PASS** |
| 12 | fibroblast | YES | Fib_Lung, Fib_Aorta, Fib_Intestine, Fib_Testis, Fib_Muscle, Fib_Heart | multiple | both | Yes (6 variants) | **PASS** |
| 13 | cardiac fibroblast | YES | Fib_Heart | heart | snRNA | Count uncertain | **BORDERLINE** |
| 14 | granulocyte | YES | Neu, Bas, Mast | blood, BM, multiple | 10x | Yes (poolable) | **PASS** |
| 15 | hematopoietic precursor | PARTIAL | Progenitor | BM | 10x | Count uncertain | **BORDERLINE** |
| 16 | HSC | NO | — | Progenitor may include but not separated | — | — | **ABSENT** |
| 17 | hepatocyte | YES | Hep | liver | 10x | Yes (dominant) | **PASS** |
| 18 | intermediate monocyte | NO | — | not subtyped | — | — | **ABSENT** |
| 19 | large intestine goblet | NO | — | no large intestine | — | — | **ABSENT** |
| 20 | luminal epithelial mammary | NO | — | no mammary | — | — | **ABSENT** |
| 21 | macrophage | YES | Kup, AMac, Mac1, Mac2 | liver, lung, multiple | both | Yes (4+ variants) | **PASS** |
| 22 | mature NKT | YES | NKT | spleen, blood | 10x | Rare, count uncertain | **BORDERLINE** |
| 23 | MSC | NO | — | — | — | — | **ABSENT** |
| 24 | MSC of adipose | NO | — | — | — | — | **ABSENT** |
| 25 | monocyte | YES | Mono, Mono_BM | blood, BM | 10x | Yes | **PASS** |
| 26 | myeloid dendritic | YES | mDC | spleen, blood | 10x | Count uncertain | **BORDERLINE** |
| 27 | myeloid leukocyte | PARTIAL | — | poolable from Mono+Mac+Neu+DC | — | — | **BORDERLINE** |
| 28 | NK cell | YES | NK | spleen, blood | 10x | Yes | **PASS** |
| 29 | neutrophil | YES | Neu, ProNeu | blood, BM | 10x | Yes | **PASS** |
| 30 | non-classical monocyte | NO | — | not subtyped | — | — | **ABSENT** |
| 31 | pancreatic acinar | NO | — | no pancreas | — | — | **ABSENT** |
| 32 | pancreatic ductal | NO | — | no pancreas | — | — | **ABSENT** |
| 33 | plasma cell | YES | Pla | spleen, BM | 10x | Rare, count uncertain | **BORDERLINE** |
| 34 | smooth muscle | YES | SMC | aorta, (sn: muscle) | both | Yes | **PASS** |
| 35 | stromal cell | NO | — | — | — | — | **ABSENT** |

### Summary

| Status | Count | Percentage |
|--------|-------|------------|
| PASS | 14 | 40% |
| BORDERLINE | 8 | 23% |
| ABSENT | 13 | 37% |

**Missing tissues that kill types:** pancreas (2 types), mammary (1), bladder (1),
large intestine (2), skin/airway (1 — basal). Monocyte subtyping gaps kill 3 more
(classical, intermediate, non-classical as separate categories).

**Note:** Coverage is LOWER than MCA (19/35 PASS) despite more tissues, because
MCA has broader tissue sampling (51 tissues vs 14) and our 35-type ontology includes
many specialized types from tissues Sun et al. did not profile. However, 14 types
includes ALL 6 original Phase 2 types (B cell, CD4+ T, CD8+ T, endothelial,
hepatocyte, macrophage), which is the most critical subset for replication.

---

## Question 3: Author Overlap with Tabula Consortium

### VERDICT: CLEAN — Zero overlap

**Sun et al. 2023 authors (18):**
Shuhui Sun, Shuai Ma, Yusheng Cai, Si Wang, Jie Ren, Yuanhan Yang, Jiale Ping,
Xuebao Wang, Yiyuan Zhang, Haoteng Yan, Wei Li, Concepcion Rodriguez Esteban,
Yan Yu, Feifei Liu, Juan Carlos Izpisua Belmonte, Weiqi Zhang, Jing Qu, Guang-Hui Liu

**Affiliations:** Chinese Academy of Sciences (Beijing), Xuanwu Hospital (Beijing),
Altos Labs (San Diego — Rodriguez Esteban and Izpisua Belmonte only)

**Tabula Muris Senis (2020 Nature) authors (134+ named):**
Full consortium at Stanford University and Chan Zuckerberg Biohub (San Francisco).
Key names: Stephen Quake, Angela Oliveira Pisco, Spyros Darmanis, Norma Neff,
Jim Karkanias, Tony Wyss-Coray, Irving Weissman, Mark Krasnow, Nicholas Schaum,
et al.

**Comparison result:**
- Exact name matches: **ZERO**
- Same-last-name + same-first-initial matches: **ZERO** (Wang, Zhang, Li, Liu, Yu
  appear in both but are completely different people with different first names)
- Institutional overlap: **ZERO** (CAS/Beijing vs Stanford/CZ Biohub)
- Izpisua Belmonte on any Tabula paper: **NO** (PubMed confirmed)
- Rodriguez Esteban on any Tabula paper: **NO** (PubMed confirmed)

---

## FINAL VERDICT: PROCEED

| Question | Result | Gate |
|----------|--------|------|
| Q1: Sedentary control separability | **CLEAN** | PASS |
| Q2: Type coverage (≥12 required) | **14 PASS** | PASS (marginal) |
| Q3: Author overlap | **CLEAN** | PASS |

**Decision: PROCEED with Sun et al. 2023 download.**

All three gates pass. Q2 is marginal (14 vs threshold 12) but includes all 6
original Phase 2 types. The 8 BORDERLINE types may push the count to 17-20 after
downloading and checking actual cell counts.

### Recommended download strategy
1. Download YC (young control) files ONLY from OMIX002605 for the 9 scRNA-seq tissues
2. Skip snRNA-seq tissues (brain, cerebellum, spinal cord, heart, skeletal muscle) —
   protocol mismatch concern, and these tissues contribute few of our target types
3. Also download OC (old control) as backup — doubles cell counts if YC alone is thin
4. Verify file format inside tar archives before committing to full download
5. Map cell type annotations to our 35-type ontology
6. Check per-type cell counts against ≥500 gate

### Caveats to document
- Exercise study, not baseline atlas — using control arm only
- Only male mice (C57BL/6J males)
- 3 biological replicates merged in processed files — no per-mouse resolution
- Small intestine profiled (not large intestine) — enterocyte/goblet mapping is approximate
- 10x Chromium 3' v3 for 9 scRNA-seq tissues — protocol-matched to TMS 10x component
- CellRanger + CellBender pipeline (not identical to our CPM + log1p normalization —
  will need to re-normalize from counts)
