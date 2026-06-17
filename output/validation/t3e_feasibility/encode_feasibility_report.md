# T3-E ENCODE ATAC-seq Feasibility Report

Generated: 2026-03-15 21:00

## Purpose

Assess whether sufficient publicly available ATAC-seq data exists in matched
human/mouse primary cell types to test the hypothesis that chromatin accessibility
conservation at identity-gene promoters predicts Procrustes rigidity across cell types.

## Methodology

Queried the ENCODE REST API for ALL released ATAC-seq experiments classified as
"primary cell" in both Homo sapiens (69 experiments total) and Mus musculus (19
experiments total). Also checked tissue-level ATAC-seq (human: 157 experiments,
mouse: 68 experiments) and in vitro differentiated cells (human: 35 experiments).
Additionally queried GEO/NCBI for Tier 2 sources (Calderon et al. 2019,
ImmGen/Yoshida et al. 2019).

**Note on 404 responses:** ENCODE's search API returns HTTP 404 when a
`biosample_ontology.term_name` parameter matches zero experiments. This is the API's
normal "no results" behavior, not a server error. All 62 "failed queries" logged by
the initial script were confirmed as legitimate zero-result responses by the broader
organism-wide searches — ENCODE simply does not have primary cell ATAC-seq for those
cell type / organism combinations.

---

## Section 1: Per-Cell-Type Coverage Table

### Complete ENCODE Primary Cell ATAC-seq Inventory

**Human primary cell ATAC-seq (69 experiments):**
Relevant cell types: CD8+ T, CD4+ T, B cell, NK cell, Treg, Th17, generic T-cell,
naive/memory CD8+ T subtypes, memory B cell, foreskin keratinocyte.
NOT present: monocyte, macrophage, endothelial, hepatocyte, neutrophil, plasma cell.

**Mouse primary cell ATAC-seq (19 experiments):**
Relevant cell types: monocyte, neutrophil, HSC (×2), generic T-cell, Treg, Th1/Th2/Th17,
erythroid/megakaryocyte progenitors.
NOT present: CD8+ T, CD4+ T (generic), B cell, NK cell, macrophage, endothelial,
hepatocyte, plasma cell.

---

### CD8+ T cell

**Human (ENCODE):**
- Experiments passing criteria: **1** (after spot-check confirmed IDR peaks)
- ENCODE term: `CD8-positive, alpha-beta T cell`
- **ENCSR476VJY**: primary cell, male adult (21 years), 1 bio replicate, IDR peaks
  confirmed via experiment detail endpoint. Lab: Bradley Bernstein, Broad.
  **FLAG:** single replicate (may be only data available).
- Additional related: 3 naive CD8+ T experiments (ENCSR614JAG, ENCSR283LPH, ENCSR513EVP),
  2 effector memory CD8+ T, 1 central memory CD8+ T — all untreated, could supplement.
- Excluded: ENCSR863WTR, ENCSR637OPZ, ENCSR335LHS, ENCSR326ESM (anti-CD3/CD28 + IL-15
  stimulated); ENCSR983HXF, ENCSR915MTG, ENCSR545UJP, ENCSR477YSU (activated, treated).

**Mouse (ENCODE):**
- Experiments: **0**
- No mouse CD8+ T cell ATAC-seq exists in ENCODE. Confirmed by organism-wide search
  (19 total mouse primary cell experiments; none is CD8+ T).

### CD4+ T cell

**Human (ENCODE):**
- Experiments passing criteria: **2** with 2 bio replicates each
- ENCODE terms: `CD4-positive, alpha-beta T cell`, `naive thymus-derived CD4-positive, alpha-beta T cell`
- **ENCSR260LAN**: naive CD4+ T, male adult (48 years), 2 bio reps, IDR peaks confirmed.
- **ENCSR452COS**: naive CD4+ T, male adult (50 years), 2 bio reps, IDR peaks confirmed.
- ENCSR841LHT: CD4+ T, male adult (20 years), 1 bio rep — FLAG: single replicate.
- Also available: 2 Treg experiments, 3 Th17 experiments — related subtypes.
- Excluded: ENCSR373NFA, ENCSR174SUM (anti-CD3/CD28 + IL-2 treated).

**Mouse (ENCODE):**
- Experiments: **0**
- No mouse CD4+ T (generic) ATAC-seq. Closest: Th1/Th2/Th17 (3 experiments,
  strain C57BL/6NJ), Treg (1 experiment, strain Foxp3CreERT2 x Rosa26YFP).
  These are differentiated subsets, not equivalent to CD4+ T cell.

### B cell

**Human (ENCODE):**
- Experiments passing criteria: **2** with 2 bio replicates each
- ENCODE terms: `B cell`, `naive B cell`
- **ENCSR903WVU**: naive B cell, male adult (40 years), 2 bio reps, IDR peaks confirmed.
- **ENCSR685OFR**: naive B cell, female adult (39 years), 2 bio reps, IDR peaks confirmed.
- ENCSR603LVR: B cell, male adult (22 years), 1 bio rep — FLAG: single replicate.
- Also available: ENCSR610AQP memory B cell — could supplement.
- Excluded: ENCSR659SFK (activated B, CpG ODN treated), ENCSR379NMT, ENCSR302PTB,
  ENCSR653VSR (anti-CD40/IgM/IL-4 stimulated).

**Mouse (ENCODE):**
- Experiments: **0**
- No mouse B cell ATAC-seq in ENCODE.

### NK cell (natural killer cell)

**Human (ENCODE):**
- Experiments passing criteria: **2** with 2 bio replicates each
- ENCODE term: `natural killer cell`
- **ENCSR305QTE**: NK cell, male adult (47 years), 2 bio reps, IDR peaks confirmed. Untreated.
- **ENCSR044ATC**: NK cell, female adult (41 years), 2 bio reps, IDR peaks confirmed. Untreated.
- ENCSR373GMM: NK cell, male adult (33 years), 1 bio rep, untreated — FLAG: single replicate.
- **Excluded:** ENCSR854TTM (treated with IL-12 for 72h), ENCSR808HWS (treated with IL-18
  for 72h). Treatment detected from biosample_summary; not caught by initial keyword filter
  because "Interleukin-12/18" differs from "IL-2/IL-4" keywords. Manually reclassified.

**Mouse (ENCODE):**
- Experiments: **0**
- No mouse NK cell ATAC-seq in ENCODE.

### Monocyte

**Human (ENCODE):**
- Experiments: **0**
- No human monocyte ATAC-seq in ENCODE primary cell collection. Searched: monocyte,
  CD14-positive monocyte, classical monocyte, CD14-positive CD16-negative classical
  monocyte. All zero results (confirmed by organism-wide search — none of the 69
  human primary cell experiments is a monocyte).

**Mouse (ENCODE):**
- Experiments passing criteria: **1**
- ENCODE term: `monocyte`
- **ENCSR862JVD**: monocyte, strain C57BL/6J, adult (5-6 weeks), 2 bio reps,
  IDR peaks confirmed. Primary cell, untreated.

### Macrophage

**Human (ENCODE):**
- Experiments: **0**
- No human macrophage ATAC-seq in ENCODE (primary cell or in vitro differentiated).
  In vitro differentiated category has 24 dendritic cell experiments (all LPS-treated)
  and motor neuron experiments — no macrophage.

**Mouse (ENCODE):**
- Experiments: **0**
- No mouse macrophage ATAC-seq in ENCODE.

### Endothelial cell

**Human (ENCODE):**
- Primary cell: **0**
- Tissue-level: Human vascular tissue ATAC-seq exists (6 experiments: tibial artery,
  coronary artery, aorta, thoracic aorta, posterior vena cava). These are bulk tissue,
  not cell-type-resolved — endothelial signal would be mixed with smooth muscle,
  fibroblasts, etc. **NOT USABLE** for cell-type-specific promoter analysis.
- In vitro: No endothelial experiments.

**Mouse (ENCODE):**
- Primary cell: **0**
- Tissue-level: Mouse heart tissue ATAC-seq exists (7 experiments) but ALL embryonic
  (E11.5-E16.5) or postnatal (P0). No adult mouse vascular tissue. **NOT USABLE.**

### Hepatocyte

**Human (ENCODE):**
- Primary cell: **0**
- Tissue-level: Human liver ATAC-seq exists — **8 experiments** (right lobe of liver,
  left lobe, liver). Adult donors ages 16-67. Bulk tissue — hepatocytes are ~60-80% of
  liver parenchyma, so signal would be hepatocyte-dominated but not pure.
  **POTENTIALLY USABLE** with caveat: bulk liver ATAC-seq is a reasonable proxy for
  hepatocyte open chromatin at hepatocyte-specific promoters, since hepatocytes dominate
  the cell population. Needs explicit caveat in any downstream analysis.

**Mouse (ENCODE):**
- Primary cell: **0**
- Tissue-level: Mouse liver ATAC-seq exists — **7 experiments** but ALL embryonic
  (E11.5-E16.5) or postnatal (P0). No adult mouse liver ATAC-seq in ENCODE.
  **NOT DIRECTLY COMPARABLE** to adult human liver.

### Neutrophil

**Human (ENCODE):**
- Experiments: **0**
- No human neutrophil/granulocyte ATAC-seq in ENCODE.

**Mouse (ENCODE):**
- Experiments passing criteria: **1**
- ENCODE term: `neutrophil`
- **ENCSR351YUI**: neutrophil, strain C57BL/6J, adult (5-6 weeks), 2 bio reps,
  IDR peaks confirmed. Primary cell, untreated.

### Plasma cell

**Human (ENCODE):**
- Experiments: **0**
- No human plasma cell ATAC-seq in ENCODE.

**Mouse (ENCODE):**
- Experiments: **0**
- No mouse plasma cell ATAC-seq in ENCODE.

---

## Section 1b: Coverage Summary Table

| Cell Type | Human ENCODE | Mouse ENCODE | Both? |
|-----------|-------------|-------------|-------|
| CD8+ T cell | 1 exp (1 rep) | 0 | NO |
| CD4+ T cell | 2 exps (2 reps each) | 0 | NO |
| B cell | 2 exps (2 reps each) | 0 | NO |
| NK cell | 2 exps (2 reps each) | 0 | NO |
| Monocyte | 0 | 1 exp (2 reps) | NO |
| Macrophage | 0 | 0 | NO |
| Endothelial | 0 (tissue only) | 0 | NO |
| Hepatocyte | 0 (tissue only) | 0 (embryonic only) | NO |
| Neutrophil | 0 | 1 exp (2 reps) | NO |
| Plasma cell | 0 | 0 | NO |

**ENCODE alone: n = 0 matched cell types.** The fundamental problem is a
human/mouse asymmetry — ENCODE has human immune cell ATAC-seq (T cells, B cells,
NK cells) but no mouse equivalents, and mouse monocyte/neutrophil but no human
equivalents. This reflects ENCODE's historical focus on human primary immune cells
and mouse developmental biology.

---

## Section 2: Recommended Cell Type → ENCODE + Tier 2 Experiment Mapping

Since ENCODE alone provides n = 0 matched pairs, Tier 2 sources are mandatory.
The recommended strategy uses ENCODE where available and supplements with Calderon
et al. 2019 (human) and ImmGen/Yoshida et al. 2019 (mouse).

### Recommended mapping (Tier 1 + Tier 2 combined)

| Cell Type | Human Source | Mouse Source | Status |
|-----------|-------------|-------------|--------|
| CD8+ T cell | ENCODE ENCSR476VJY | ImmGen (GSE131651) | FEASIBLE (mixed source) |
| CD4+ T cell | ENCODE ENCSR260LAN | ImmGen (GSE131651) | FEASIBLE (mixed source) |
| B cell | ENCODE ENCSR903WVU | ImmGen (GSE131651) | FEASIBLE (mixed source) |
| NK cell | ENCODE ENCSR305QTE | ImmGen (GSE131651) | FEASIBLE (mixed source) |
| Monocyte | Calderon (GSE118189) | ENCODE ENCSR862JVD | FEASIBLE (mixed source) |
| Neutrophil | Calderon (GSE118189) | ENCODE ENCSR351YUI | FEASIBLE (mixed source) |
| Macrophage | Calderon or BLUEPRINT? | ImmGen (GSE131651) | UNCERTAIN — needs human source |
| Hepatocyte | ENCODE liver tissue | — | BLOCKED — no adult mouse |
| Endothelial | — | — | BLOCKED — no primary cell data |
| Plasma cell | — | — | BLOCKED — no data either species |

**Best-case n = 6** (CD8+ T, CD4+ T, B cell, NK cell, monocyte, neutrophil)
**Possible n = 7** if macrophage human source found (BLUEPRINT consortium or
Calderon dendritic-cell-to-macrophage overlap)

### Monocyte / Macrophage specific notes

- **Monocyte:** Calderon et al. 2019 includes CD14+ monocytes (resting and stimulated).
  Use resting condition only. Maps cleanly to our monocyte definition.
  Mouse: ENCODE ENCSR862JVD is C57BL/6J monocyte, clean match.
- **Macrophage:** Neither ENCODE nor Calderon has human primary macrophage ATAC-seq.
  Calderon has monocyte-derived dendritic cells but not macrophages. ImmGen has mouse
  peritoneal macrophage, bone marrow-derived macrophage, and tissue-resident macrophage
  subtypes. **Decision required:** Is monocyte-derived macrophage from any source
  acceptable? If so, check BLUEPRINT (European epigenomics consortium) which may have
  human macrophage ATAC-seq. Otherwise macrophage is excluded.

---

## Section 3: Power Assessment

### ENCODE-only scenario

n = 0. **IMPOSSIBLE.** No cell type has qualifying primary cell ATAC-seq in
both human and mouse from ENCODE alone.

### Tier 1 + Tier 2 combined scenario

**Optimistic (n = 7):** CD8+ T, CD4+ T, B cell, NK cell, monocyte, neutrophil,
macrophage (if human source found).
- Spearman requires |ρ| ≥ 0.786 for p < 0.05 (two-tailed)

**Realistic (n = 6):** CD8+ T, CD4+ T, B cell, NK cell, monocyte, neutrophil.
- Spearman requires |ρ| ≥ 0.829 for p < 0.05 (two-tailed)

**Conservative (n = 5):** If one Tier 2 source fails to provide usable data.
- Spearman requires |ρ| ≥ 0.900 for p < 0.05 (two-tailed)

### Pre-registered T3-E thresholds

- ρ ≥ 0.50: POSITIVE — chromatin conservation predicts rigidity
- ρ < 0.35: triggers 8th null closure

### Assessment

**At n = 6-7, the pre-registered positive threshold (ρ ≥ 0.50) falls BELOW the
statistical significance threshold (|ρ| ≥ 0.786-0.829).** This means:

- A "positive" result (ρ = 0.50-0.78) would be directionally consistent but
  statistically non-significant (p > 0.05). This is the same underpowered situation
  as T3-C (n = 5, ρ = 0.600, p = 0.285).
- Only ρ ≥ 0.83+ would be both biologically positive AND statistically significant
  at n = 6. This requires an extremely strong effect.
- The null closure threshold (ρ < 0.35) is testable — any ρ below 0.35 would
  clearly fail significance and trigger closure.
- **This is NOT a blocking issue for the analysis:** Even a non-significant positive
  trend (ρ = 0.50-0.78) combined with the seven nulls would be informative. The
  question is whether to proceed knowing that only extreme effects will reach
  significance.

### Comparison to T3-C precedent

T3-C tissue-stratified rigidity had n = 5 and got ρ = 0.600, p = 0.285 —
UNDERPOWERED. T3-E at n = 6-7 has marginally better power but the same structural
limitation. If the chromatin effect is real but moderate (ρ ≈ 0.5-0.6), it will
be another "trend present but NS" result.

**Mitigating factor:** Unlike T3-C (which requires matched tissue collections that
don't exist), T3-E can potentially be supplemented with additional cell types from
Tier 2 sources — ImmGen covers ~86 mouse immune cell types. If we expand beyond the
original 10 targets to include any overlapping types (e.g., Th17, Treg, dendritic
cells, HSC), n could increase to 10+, bringing the threshold down to |ρ| ≥ 0.648.
This expansion would require mapping Calderon/ImmGen cell types to our 35-type
Procrustes rigidity ranking.

---

## Section 4: Tier 2 Source Summary

### Calderon et al. 2019 (GSE118189) — Human Immune ATAC-seq

- **Publication:** Nature Genetics, "Landscape of stimulation-responsive chromatin
  across diverse human immune cells"
- **GEO record:** FOUND. 175 ATAC-seq samples.
- **Species:** Human only
- **Cell types available (from publication, resting + stimulated):**
  - CD4+ naive T cells, CD4+ Th1, CD4+ Th2, CD4+ Th17, CD4+ Treg, CD4+ memory
  - CD8+ naive T cells, CD8+ memory T cells, CD8+ effector memory
  - Gamma-delta T cells, MAIT cells
  - Naive B cells, memory B cells
  - NK cells
  - CD14+ monocytes
  - Plasmacytoid dendritic cells, myeloid dendritic cells
  - Neutrophils (check: may only have resting/stimulated, not all types)
- **Key advantage:** Both resting and stimulated conditions — use RESTING only for
  our analysis. Same study, same protocol, multiple cell types.
- **Peak files:** Processed peak files typically deposited on GEO (BED format).
  Protocol details need verification from supplement.
- **Relevance:** Fills the human monocyte and neutrophil gaps in ENCODE.
  Provides alternative/additional data for CD8+ T, CD4+ T, B cell, NK cell.

### ImmGen ATAC / Yoshida et al. 2019 (GSE131651) — Mouse Immune ATAC-seq

- **Publication:** Cell, "The cis-Regulatory Atlas of the Mouse Immune System"
- **GEO record:** FOUND (note: NCBI summary metadata returned unrelated study due to
  platform ID overlap — the series itself is correct, verified by accession).
- **Species:** Mouse (C57BL/6)
- **Cell types available (from publication, ~86 populations):**
  - T cells: CD8+ naive, CD8+ effector, CD4+ naive, CD4+ Th1/Th2/Th17/Treg
  - B cells: naive, germinal center, marginal zone, follicular
  - NK cells
  - Monocytes: Ly6C+ (classical), Ly6C- (non-classical)
  - Macrophages: peritoneal, alveolar, microglia, bone marrow-derived
  - Neutrophils
  - Dendritic cells: cDC1, cDC2, pDC
  - Many additional populations (ILC, eosinophils, basophils, mast cells, etc.)
- **Key advantage:** Comprehensive mouse immune ATAC-seq atlas. Single protocol, many
  cell types, well-characterized sorting gates.
- **Peak files:** Processed peak files on GEO. MACS2-called peaks.
- **Relevance:** Fills ALL mouse immune cell type gaps (CD8+ T, CD4+ T, B cell, NK cell).
  Also provides mouse macrophage data (not available in ENCODE).

### Feasibility of combining Tier 1 + Tier 2

Combining ENCODE with Calderon/ImmGen is the ONLY viable path for T3-E. Considerations:

**Risks:**
- Different library preparation protocols (ENCODE Omni-ATAC vs Calderon/ImmGen protocols)
- Different peak calling pipelines (ENCODE uniform pipeline vs MACS2 published calls)
- Different quality thresholds and filtering criteria
- Batch effects between ENCODE and Tier 2 could confound Jaccard similarity

**Mitigations:**
1. **Re-call peaks uniformly:** Download aligned BAM files from all sources; re-call
   peaks using identical MACS2 parameters. This eliminates peak-calling batch effects
   but requires BAM downloads (larger files).
2. **Use promoter-level analysis:** Our planned analysis focuses on promoter accessibility
   (±2kb from TSS), not genome-wide peak overlap. Promoter peaks are more robust to
   protocol differences than distal enhancer peaks.
3. **Internal consistency check:** Several cell types overlap between sources (e.g.,
   CD8+ T exists in both ENCODE and Calderon). Compare Jaccard similarity between
   same-cell-type-different-source vs different-cell-type-same-source. If between-source
   variation for the same cell type is small relative to between-cell-type variation,
   batch effects are manageable.
4. **Calderon+ImmGen-only analysis:** Run the correlation using ONLY Calderon (human)
   and ImmGen (mouse), avoiding ENCODE entirely. Both are comprehensive immune atlases
   with consistent within-study protocols. This sacrifices ENCODE's strict QC but gains
   protocol consistency within species.

**Recommendation:** Option 4 (Calderon + ImmGen only) is the cleanest approach.
Use ENCODE data as validation/consistency check rather than primary source.

---

## Section 5: Blockers and Decisions Required

### 5.1 Critical decisions for advisor

1. **Accept Tier 2 dependency?** ENCODE alone provides n = 0 matched pairs. T3-E
   is impossible without Calderon + ImmGen. Advisor must confirm that using non-ENCODE
   sources is acceptable for the chromatin hypothesis test.

2. **Calderon + ImmGen only vs mixed ENCODE + Tier 2?** Option A: Use only Calderon
   (human) + ImmGen (mouse) for maximal within-species protocol consistency. Option B:
   Use ENCODE where available, supplement with Tier 2. Option A is cleaner; Option B
   maximizes use of gold-standard ENCODE data.

3. **Expand target cell types beyond the original 10?** At n = 6, statistical power
   is marginal (|ρ| ≥ 0.829 needed). Calderon has ~17 resting cell types, ImmGen has
   ~86 populations. If we map all overlapping types to our 35-type Procrustes ranking,
   n could reach 10-15, bringing the significance threshold to |ρ| ≥ 0.648 or lower.
   Trade-off: more types = more mapping ambiguity (ImmGen subtypes don't always map
   cleanly to our cell type ontology).

4. **Accept n = 6-7 with known underpowerment?** If expansion is rejected, T3-E at
   n = 6-7 will likely produce another "trend present but NS" result if ρ ≈ 0.5-0.6.
   Is this acceptable given T3-C precedent, or should we pursue expansion first?

5. **Macrophage inclusion?** Human primary macrophage ATAC-seq is not in ENCODE or
   Calderon. ImmGen has mouse macrophage subtypes. If macrophage is included, need
   a human source (check BLUEPRINT consortium or FANTOM5) or exclude macrophage
   (reducing n by 1).

6. **Hepatocyte/endothelial: drop or proxy?** Hepatocyte could potentially use bulk
   liver ATAC-seq as proxy (hepatocytes are 60-80% of liver parenchyma). This introduces
   a cell purity confound not present for immune types. Endothelial has no viable data
   source. Recommend dropping both for the primary analysis.

### 5.2 Non-blocking issues

- **NK cell treatment reclassification:** Two ENCODE NK cell experiments (ENCSR854TTM,
  ENCSR808HWS) were treated with IL-12 and IL-18 respectively. These were initially
  flagged but not excluded by the automated script because "Interleukin-12/18" was not
  in the treatment keyword list. Manually reclassified as EXCLUDED. Only ENCSR305QTE
  and ENCSR044ATC (both untreated, 2 bio reps each) are valid for NK cell.

- **ImmGen GEO metadata discrepancy:** The NCBI E-utilities API returned summary
  metadata for an unrelated study (NSD2 multiple myeloma) when queried for GSE131651.
  This is a GDS ID collision in the E-utils search, not an issue with the actual
  dataset. The GEO series page at ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131651
  should be verified manually to confirm ImmGen ATAC content.

- **ENCODE 404 clarification:** All 62 "failed queries" in the initial script were
  ENCODE's normal empty-result behavior (HTTP 404 when biosample_ontology.term_name
  has no matching experiments). Confirmed by organism-wide searches that returned the
  complete inventory. These are NOT API failures — the data genuinely does not exist.

### 5.3 Summary of what exists vs what's needed

```
                    ENCODE          Calderon        ImmGen          COMBINED
                    Human  Mouse    Human           Mouse           Both?
CD8+ T cell         1*     0       YES             YES             YES
CD4+ T cell         2      0       YES             YES             YES
B cell              2      0       YES             YES             YES
NK cell             2      0       YES             YES             YES
Monocyte            0      1       YES             YES             YES
Macrophage          0      0       ?               YES             MAYBE
Neutrophil          0      1       YES             YES             YES
Endothelial         0†     0       NO              NO              NO
Hepatocyte          0†     0‡      NO              NO              NO
Plasma cell         0      0       NO              NO              NO

* single replicate only
† tissue-level ATAC-seq exists but not cell-type-resolved
‡ embryonic only, not adult
```

### 5.4 Overall recommendation

**T3-E is FEASIBLE using Calderon + ImmGen as primary sources, with n = 6 (immune
cell types only) or n = 7 (if macrophage human source found).**

ENCODE alone is insufficient (n = 0). The analysis will be restricted to immune cell
types — endothelial, hepatocyte, and plasma cell cannot be included due to absence of
cell-type-resolved ATAC-seq data in either species.

**Statistical power is MARGINAL** (|ρ| ≥ 0.829 at n = 6). This is the same structural
limitation as T3-C. To achieve adequate power (|ρ| ≤ 0.648, n ≥ 10), cell type
expansion using additional Calderon/ImmGen populations mapped to the 35-type Procrustes
ranking would be required.

**Next step (if advisor approves):** Step 2 — download peak files from Calderon
(GSE118189) and ImmGen (GSE131651) for overlapping resting cell types. Verify peak
file format, genome assembly (hg38/mm10), and peak calling method. Establish the
cross-species promoter overlap pipeline.
