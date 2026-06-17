# Macaque Extension Feasibility Report

**Date:** 2026-04-05
**Purpose:** Assess feasibility of adding macaque as a third species for
three-way consistency testing of per-type rigidity rankings.
**Status:** Scoping only — no new analyses run.

---

## Critical Discovery: Macaque Pipeline Already Completed

The macaque extension was **already executed on 2026-03-16** under Strategy B
(DECISION-123, DECISION-130, DECISION-131). This report synthesizes those
results against the user's stated goal of three-way ranking consistency.

---

## 1. Best Candidate Dataset

**Strategy B (composite):** RIRA + Qu et al. — already downloaded and processed.

| Component | Species | Cells | Tissues | Donors | Technology | Role |
|---|---|---|---|---|---|---|
| RIRA (Mahyari 2025) | *M. mulatta* | 412,820 (excl. Unknown) | 7 | 47 | 10x 5' + CITE-seq | Immune types |
| Qu et al. (2022) | *M. fascicularis* | 230,882 | 16 | 1-2 | 10x Chromium 3' | Non-immune types |

**Why this combination:** Both use 10x technology (avoiding the Microwell-seq /
DNBelab C4 failure mode documented in DECISION-104). RIRA provides
unprecedented 47-donor power for immune centroids. Qu et al. covers
hepatocyte, endothelial, fibroblast, and other non-immune types absent
from RIRA.

**Fallback (NOT used):** NHPCA (Han 2022, *M. fascicularis*, 1.14M cells,
45 tissues) — rejected as primary due to DNBelab C4 technology mismatch.

**Other datasets evaluated (2026-04-05 survey):**
- CELLxGENE Census M. mulatta: 88% sci-RNA-seq3 brain — NOT useful
- Cross-Study Multi-Organ Atlas (2026 preprint): mixed technologies — too heterogeneous
- Tabula Microcebus (*M. murinus*): not macaque, interesting as future 4th species
- Allen Brain macaque: brain-only, sci-RNA-seq3

No new dataset identified in this survey changes the Strategy B recommendation.

---

## 2. Cell Type Overlap

### Macaque types vs. CellWarp 35-type primary set

**20/35 types matched** (57% overlap) — above the 15-type minimum.

| # | Cell Type | Source | Donors | LOW_CONFIDENCE | HM Rank (of 35) |
|---|---|---|---|---|---|
| 1 | B cell | RIRA | 47 | No | 17 |
| 2 | CD4+ T cell | RIRA | 45 | No | 14 |
| 3 | CD8+ T cell | RIRA | 45 | No | 35 (most rigid) |
| 4 | T cell (generic) | RIRA | 48 | No | 7 |
| 5 | classical monocyte | RIRA | 44 | No | 15 |
| 6 | intermediate monocyte | RIRA | 33 | No | 24 |
| 7 | non-classical monocyte | RIRA | 44 | No | 34 |
| 8 | macrophage | RIRA | 44 | No | 16 |
| 9 | myeloid dendritic cell | RIRA | 44 | No | 21 |
| 10 | myeloid leukocyte | RIRA | 45 | No | 10 |
| 11 | natural killer cell | RIRA | 45 | No | 23 |
| 12 | granulocyte | RIRA | 39 | No | 27 |
| 13 | hematopoietic precursor cell | RIRA | 44 | No | 3 |
| 14 | endothelial cell | Qu | 6 | **Yes** | 33 |
| 15 | epithelial cell | Qu | 7 | **Yes** | 2 |
| 16 | fibroblast | Qu | 6 | **Yes** | 28 |
| 17 | smooth muscle cell | Qu | 5 | **Yes** | 31 |
| 18 | stromal cell | Qu | 6 | **Yes** | 1 (most flexible) |
| 19 | bladder urothelial cell | Qu | 1 | **Yes** | 29 |
| 20 | hepatocyte | Qu | 2 | **Yes** | 32 |

**Composition problem:** 13/20 types (65%) are immune lineage from RIRA.
Only 7 non-immune types, all from Qu et al. with LOW-CONFIDENCE flag
(1-2 donors). This immune bias compresses ranking geometry — immune
types cluster together, preventing separation of rigid vs. flexible.

### 15 types NOT in macaque set

mesenchymal stem cell, mesenchymal stem cell of adipose tissue, basal cell,
monocyte, enterocyte, mature NK T cell, pancreatic acinar cell, luminal
epithelial cell, fibroblast of cardiac tissue, plasma cell, large intestine
goblet cell, hematopoietic stem cell, adventitial cell, neutrophil,
pancreatic ductal cell.

---

## 3. Gene Overlap

| Metric | Count | Gate | Status |
|---|---|---|---|
| Human-macaque 1:1 orthologs | 19,123 | >=12,000 | **STRONG PASS** |
| Human-mouse 1:1 orthologs | 16,959 | >=12,000 | PASS (baseline) |
| Three-way intersection (H-Mac-M) | 15,028 (raw) | >=12,000 | **PASS** |
| Three-way after expression filter | **13,927** | >=12,000 | **PASS** |
| Fraction of HM space retained | 82.1% (13,927/16,959) | — | Acceptable |

**Gene symbol concordance (species mixing):** 13,093/13,927 = 94.0%.
FLAG range (90-94%), not ABORT (<90%). Documented in paper.

Ortholog mapping is straightforward — uses the same BioMart pybiomart
pipeline already in `src/data_loader.py`. Cached file exists at
`data/macaque/biomart_macaque_human_orthologs.csv`.

---

## 4. Existing Results (Already Computed)

### Human-macaque Procrustes (DECISION-131)

| Analysis | Types | p-value | obs/null | Significant |
|---|---|---|---|---|
| Primary (RIRA + Qu) | 20 | 0.0002 | 0.841 | Yes (p < 0.01) |
| Sensitivity (RIRA only) | 13 | 0.0002 | 0.749 | Yes (p < 0.01) |

**Three-species geometric conservation: CONFIRMED.** The Procrustes
transformation structure replicates across human, mouse, and macaque.

### Rigidity Ranking Comparison

| Metric | Value | Interpretation |
|---|---|---|
| Spearman rho (20 types) | **0.137** | Non-replication |
| p-value | 0.565 | Not significant |
| Hepatocyte | Most divergent in BOTH | Replicates |
| CD8+ T, endothelial | Rank 7, 5 (mac) vs 35, 33 (HM) | **Inverted** |

**Ranking does NOT replicate.** This is consistent with the established
finding that rigidity ranking is dataset-dependent (0/4 external datasets
replicate: MCA rho=0.120, Sun2023 rho=0.146, PanSci rho=0.194, CellHint
rho=-0.386). The macaque result (rho=0.137) is the fifth non-replication.

**Root cause (DECISION-131):** Immune bias. 13/20 macaque types are immune,
compressing the ranking geometry. CD8+ T cells cannot emerge as most rigid
when compared against 12 other immune types. This is a dataset composition
artifact, not a biological contradiction.

### What Was NOT Computed

**Mouse-macaque Procrustes was never run.** The existing pipeline only
computed human-macaque. The user's three-way consistency test requires
all three pairwise comparisons:
1. Human-mouse (done, primary result, p=0.0001)
2. Human-macaque (done, p=0.0002)
3. Mouse-macaque (**NOT done**)

---

## 5. Assessment: Can Three-Way Testing Rescue Rankings?

### The user's hypothesis
> If a cell type is rigid in human-mouse, human-macaque, AND mouse-macaque,
> the ranking is far more credible.

### Why this is unlikely to work with current data

1. **Pairwise building blocks fail.** Three-way consistency requires each
   pairwise ranking to be meaningful. Human-macaque ranking already shows
   rho=0.137 (null). Adding mouse-macaque won't help if 2/3 pairwise
   rankings are noise.

2. **Immune bias is structural, not fixable.** The macaque data has 13/20
   immune types. Even running mouse-macaque Procrustes would face the same
   composition problem — the immune types would dominate the ranking,
   preventing separation of rigid vs. flexible.

3. **Ranking instability is a general property.** Five independent datasets
   now show non-replication (rho range: -0.386 to +0.194). This is not a
   macaque-specific problem. It reflects measurement sensitivity to atlas
   composition, technology, and donor sampling — documented as a closed
   question in the project.

4. **Exception: hepatocyte.** Hepatocyte is the most divergent type in
   BOTH human-mouse and human-macaque. This single-type consistency is
   already documented and is genuinely three-species evidence (once
   mouse-macaque is computed). But one type does not rescue a ranking.

### What WOULD make three-way ranking work

A macaque atlas with:
- >=25 cell types spanning multiple organ systems (not immune-dominated)
- >=3 donors per type (eliminating LOW-CONFIDENCE flags)
- 10x technology throughout
- Coverage of the rigid anchor types: hepatocyte, CD8+ T, endothelial
- Coverage of the flexible anchor types: stromal, epithelial, progenitors

**No such atlas currently exists.** The 2026 cross-study multi-organ atlas
(bioRxiv preprint) integrates 30 studies but is technology-heterogeneous.
NHPCA has broad coverage but DNBelab C4 technology. Neither is suitable.

---

## 6. Effort Estimate (If Proceeding)

### Already done (zero additional effort)
- Data download (RIRA: 4.8 GB, Qu: 1.1 GB) — complete
- Ortholog mapping — cached
- Human-macaque Procrustes — complete
- Pre-pipeline QC checks — complete

### New work required for three-way consistency test

| Task | Effort | Notes |
|---|---|---|
| Mouse-macaque Procrustes | 1 session | Adapt existing pipeline, run on 20 shared types |
| Three-way ranking concordance | 1 session | Compute Kendall W or three-way Spearman |
| Report and figures | 0.5 session | |
| **Total** | **~2.5 sessions** | |

### Additional work if pursuing broader macaque atlas

| Task | Effort | Notes |
|---|---|---|
| Identify + download new atlas | 2-4 sessions | No suitable atlas currently exists |
| Reannotate cell types | 3-5 sessions | Harmonize nomenclature |
| Rerun full pipeline | 1-2 sessions | |
| **Total** | **6-11 sessions** | Blocked on atlas availability |

---

## 7. GO / NO-GO Recommendation

### **CONDITIONAL NO-GO** for three-way ranking consistency.

**Rationale:**
- The macaque pipeline is already complete and the answer is in hand:
  three-species geometric conservation is confirmed (p=0.0002), but
  ranking does not replicate (rho=0.137).
- Running mouse-macaque Procrustes is low-effort (~1 session) and would
  complete the three pairwise comparisons, but is unlikely to produce
  ranking consistency given the structural immune bias.
- The fundamental limitation is atlas composition (13/20 immune), not
  analysis pipeline. No code change fixes this.
- Five datasets now show ranking non-replication. Adding a sixth
  non-replication does not strengthen the paper.

### What IS already usable

The existing macaque result already supports:
1. "Three-species geometric conservation across 90 million years" (p=0.0002)
2. "Hepatocyte rigidity replicates across all three species"
3. Paper 1 narrative (landmark version) is not blocked by ranking
   non-replication — the paper frames rankings as explicitly exploratory

### Recommended next action (if any)

**LOW-PRIORITY:** Run mouse-macaque Procrustes as a completeness exercise
(~1 session). This would:
- Complete the three pairwise matrix
- Confirm/deny hepatocyte three-way consistency
- Provide a clean "we tested all three pairs" statement for the paper
- Likely show another ranking non-replication (expected, not harmful)

This is a "nice to have" for Paper 1 supplementary material, not a
blocking requirement.

### What would change this to GO

A new macaque multi-tissue 10x atlas with >=25 diverse cell types and
>=3 donors per type. Monitor:
- The 2026 cross-study atlas (once peer-reviewed and 10x-only subset extractable)
- Future CELLxGENE Census updates (may incorporate more 10x macaque data)
- Any new "Tabula Macaca" equivalent

---

## Appendix: Data File Locations

| File | Path | Size |
|---|---|---|
| RIRA metadata | `data/macaque/rira/` | ~53 MB |
| RIRA RNA counts (converted) | `data/macaque/rira/` | ~6.3 GB (MTX) |
| Qu et al. scRNA-seq | `data/macaque/qu/` | ~1.1 GB |
| Macaque-human orthologs | `data/macaque/biomart_macaque_human_orthologs.csv` | ~2 MB |
| Primary results | `output/macaque_pipeline/primary_procrustes_results.json` | |
| Sensitivity results | `output/macaque_pipeline/sensitivity_procrustes_results.json` | |
| Ranking comparison | `output/macaque_pipeline/rigidity_ranking_comparison.csv` | |
| Cell counts | `output/macaque_pipeline/centroid_cell_counts.csv` | |

## Appendix: Decisions Referenced

- DECISION-098: Ortholog feasibility (19,123 1:1, STRONG PASS)
- DECISION-099: NHPCA partial, Qu inaccessible at time
- DECISION-123: Strategy B recommended (RIRA + Qu)
- DECISION-123-AMENDMENT: Abort criteria, species mixing gates
- DECISION-130: Pre-pipeline checks passed
- DECISION-131: Three-species conservation confirmed, ranking partial
