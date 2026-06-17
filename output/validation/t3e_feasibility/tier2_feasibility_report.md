# T3-E Step 1b: Calderon + ImmGen Coverage Assessment

Generated: 2026-03-15

## Purpose

Assess cell type coverage of Calderon et al. 2019 (human) and ImmGen ATAC / Yoshida
et al. 2019 (mouse) for the T3-E chromatin accessibility hypothesis test. Map available
cell types against our 35-type Tabula Procrustes set to determine realistic sample size
for the rigidity–chromatin Spearman correlation.

---

## Correction to Step 1 Report

**GSE131651 is WRONG.** The task specified GSE131651 as the ImmGen ATAC accession. This
is actually "NSD2 overexpression drives clustered chromatin and transcriptional changes
in a subset of insulated domains" — a human multiple myeloma study (48 samples, Hi-C/
ChIP-seq/RNA-seq). The correct ImmGen ATAC accession is **GSE100738**.

**Calderon does NOT have neutrophils.** The Step 1 ENCODE feasibility report (Section 2)
listed neutrophils as available from Calderon (GSE118189). This is incorrect. Calderon
et al. sorted cells from PBMCs via Ficoll density gradient, which excludes granulocytes
including neutrophils. The 25 Calderon cell types (confirmed from GEO sample titles) do
not include neutrophils. This reduces the Step 1 "best-case n=6" estimate by 1.

---

## Section 1: Calderon et al. 2019 — Human Immune ATAC-seq

**GEO accession:** GSE118189
**Paper:** Calderon D et al. "Landscape of stimulation-responsive chromatin across
diverse human immune cells." Nature Genetics 51:1494-1505 (2019).
**Total samples:** 175 ATAC-seq
**Donors:** 6 (IDs: 1001, 1002, 1003, 1004, 1008, 1010) — not all donors have all types
**Source:** Peripheral blood (PBMCs via Ficoll gradient) from healthy donors
**Naming convention:** `{donor}-{celltype}-{U|S}` where U=unstimulated, S=stimulated
**Supplementary data:** `GSE118189_ATAC_counts.txt.gz` — consensus peak count matrix
**Platforms:** GPL20301 (HiSeq 2500), GPL24676

### Unstimulated cell types (25 unique, for our analysis)

| # | Cell Type | Unstim. Donors | Notes |
|---|-----------|----------------|-------|
| 1 | Bulk_B | 4 | All B cells |
| 2 | Mem_B | 4 | Memory B |
| 3 | Naive_B | 4 | Naive B |
| 4 | Plasmablasts | 2 | Early antibody-secreting; NOT plasma cell |
| 5 | CD8pos_T | 4 | Bulk CD8+ T |
| 6 | Central_memory_CD8pos_T | 4 | |
| 7 | Effector_memory_CD8pos_T | 4 | |
| 8 | Naive_CD8_T | 4 | |
| 9 | Gamma_delta_T | 4 | Not in our 35-type set |
| 10 | Effector_CD4pos_T | 4 | |
| 11 | Follicular_T_Helper | 5 | Tfh |
| 12 | Memory_Teffs | 4 | Memory CD4+ effector T |
| 13 | Memory_Tregs | 4 | |
| 14 | Naive_Teffs | 4 | Naive CD4+ effector T |
| 15 | Naive_Tregs | 3 | |
| 16 | Regulatory_T | 4 | |
| 17 | Th1_precursors | 4 | |
| 18 | Th2_precursors | 4 | |
| 19 | Th17_precursors | 3 | |
| 20 | Monocytes | 3 | CD14+ classical monocytes |
| 21 | Myeloid_DCs | 3 | Conventional myeloid DC (BDCA1+/CD1c+) |
| 22 | pDCs | 3 | Plasmacytoid DC; not in our 35-type set |
| 23 | Immature_NK | 5 | CD56bright |
| 24 | Mature_NK | 5 | CD56dim |
| 25 | Memory_NK | 6 | |

**Key limitations:**
- Immune cells only (PBMCs) — no hepatocytes, endothelial, epithelial, stromal, etc.
- No neutrophils (excluded by Ficoll gradient)
- No macrophages (tissue-resident, not circulating in PBMCs)
- Multiple subtypes for each major lineage → mapping to our broad Tabula labels is ambiguous

---

## Section 2: ImmGen ATAC — Yoshida et al. 2019 — Mouse Immune ATAC-seq

**GEO accession:** GSE100738 (NOT GSE131651)
**Paper:** Yoshida H, Lareau CA, Ramirez RN et al. (Immunological Genome Project).
"The cis-Regulatory Atlas of the Mouse Immune System." Cell 176(4):897-912 (2019).
**PubMed:** 30686579
**Total samples:** 154 ATAC-seq
**Populations:** 86 unique immunocyte types
**Strain:** C57BL/6
**Protocol:** Fast-ATAC, 10,000 cells per library (flow-sorted to high purity)
**Supplementary data:** BigWig files (.bw), RAW.tar — **no processed peak BED/narrowPeak
files deposited** (peaks must be called from raw data or BigWig)
**Platform:** GPL19057

### Resting cell types (excluding experimentally stimulated)

**Excluded (stimulated/activated):**
- NKT.Sp.LPS.3hr/18hr/3d (LPS-stimulated NKT)
- T.4.Sp.aCD3+CD40.18hr (anti-CD3/CD40 stimulated CD4+ T)
- MF.pIC.Alv.Lu (polyIC-treated alveolar macrophage)
- T8.IEL/MP/TE/Tcm/Tem.LCMV.* (all LCMV infection-model CD8 T effector/memory)
- GN.Thio.PC (thioglycollate-elicited granulocytes)

**Resting populations by lineage:**

| Lineage | Populations | Reps each | Tissue |
|---------|-------------|-----------|--------|
| **B cells** | B.Fo.Sp, B.Sp, B.MZ.Sp, B.mem.Sp, B.Fem.Sp, B.T1/T2/T3.Sp, B.FrE.BM, B.GC.CB/CC.Sp, B.PB.Sp, B1b.PC | 2-3 | Spleen, BM, PerC |
| **Plasma cells** | B.PC.Sp, B.PC.BM | 2 | Spleen, BM |
| **B cell development** | proB.CLP/FrA/FrBC.BM, preB.FrD.BM | 2 | BM |
| **CD4+ T (naive)** | T.4.Nve.Sp, T.4.Nve.Fem.Sp | 2 | Spleen |
| **CD8+ T (naive)** | T.8.Nve.Sp, T8.TN.P14.Sp | 2 | Spleen |
| **Thymic development** | preT.DN1/2a/2b/3.Th, T.DN4/ISP/DP/4/8.Th | 2 | Thymus |
| **Gamma-delta T** | Tgd.Sp + 6 tissue-specific subtypes | 2 | Spleen, Thymus, LN |
| **Treg** | Treg.4.25hi.Sp, Treg.4.FP3+.Nrplo.Co | 2 | Spleen, Colon |
| **NK cells** | NK.27-11b+, NK.27+11b-, NK.27+11b+ (BM + Sp) | 2 | BM, Spleen |
| **NKT** | NKT.Sp (resting only) | 2 | Spleen |
| **Classical monocyte** | Mo.6C+II-.Bl | 3 | Blood |
| **Non-classical monocyte** | Mo.6C-II-.Bl | 3 | Blood |
| **Macrophages** | MF.PC, MF.Fem.PC, MF.226+/ICAM+.PC, MF.Alv.Lu, MF.LP.SI, MF.microglia.CNS, MF.RP.Sp | 2-4 | PerC, Lung, SI, CNS, Spleen |
| **Granulocyte/Neutrophil** | GN.BM, GN.Sp | 2 | BM, Spleen |
| **Dendritic cells** | DC.4+.Sp (cDC2), DC.8+.Sp (cDC1), DC.pDC.Sp, DC.103+11b-/+.SI | 2 | Spleen, SI |
| **ILCs** | ILC2.SI, ILC3 subtypes (3 types) | 2 | SI |
| **HSC/progenitors** | LTHSC.34-/34+.BM, STHSC.150-.BM, MPP3/4.BM | 1 | BM |
| **Stromal/endothelial** | BEC.SLN, LEC.SLN, FRC.SLN, IAP.SLN | 2 | LN |
| **Thymic epithelium** | Ep.MEChi.Th | 2 | Thymus |

**Key features:**
- Comprehensive mouse immune atlas — broadest available ATAC-seq resource
- Includes non-immune stromal/endothelial from lymph node (limited tissue context)
- HSC/progenitor populations have only 1 replicate each (rare cell limitation)
- No hepatocytes, no general endothelial (only LN-specific), no epithelial (only thymic)
- Multiple tissue-resident macrophage subtypes available

---

## Section 3: Three-Column Mapping Table — Our 35 Types vs Calderon vs ImmGen

### Key for Status column:
- **CLEAN PAIR** — both datasets have a clear, unambiguous match
- **AMBIGUOUS** — data exists in both species but mapping requires judgment call
- **MISSING (human)** — no Calderon match; ImmGen has mouse data
- **MISSING (mouse)** — Calderon has human data; no ImmGen match
- **MISSING (both)** — neither dataset covers this type

| # | Our 35-Type Label | Rigidity Rank | Calderon (human, unstim.) | ImmGen (mouse, resting) | Status |
|---|-------------------|---------------|--------------------------|------------------------|--------|
| 1 | stromal cell | 1 (least rigid) | — | FRC.SLN (LN-specific) | MISSING (human); mouse is too tissue-specific |
| 2 | epithelial cell | 2 | — | Ep.MEChi.Th (thymic only) | MISSING (both practical) |
| 3 | hematopoietic precursor cell | 3 | — | MPP3.48+.BM (1 rep), MPP4.135+.BM (1 rep) | MISSING (human) |
| 4 | hematopoietic stem cell | 4 | — | LTHSC.34-/34+.BM, STHSC.150-.BM (1 rep each) | MISSING (human) |
| 5 | pancreatic acinar cell | 5 | — | — | MISSING (both) |
| 6 | basal cell | 6 | — | — | MISSING (both) |
| 7 | T cell | 7 | Multiple subtypes | Multiple subtypes | **AMBIGUOUS** — generic label; double-count with CD4/CD8 entries |
| 8 | neutrophil | 8 | — | GN.BM (2), GN.Sp (2) | MISSING (human) |
| 9 | fibroblast of cardiac tissue | 9 | — | — | MISSING (both) |
| 10 | myeloid leukocyte | 10 | — | — | MISSING (both) — too broad |
| 11 | MSC of adipose tissue | 11 | — | — | MISSING (both) |
| 12 | plasma cell | 12 | Plasmablasts-U (2 donors) | B.PC.Sp (2), B.PC.BM (2) | **AMBIGUOUS** — plasmablast ≠ plasma cell |
| 13 | mesenchymal stem cell | 13 | — | — | MISSING (both) |
| 14 | CD4+ alpha-beta T cell | 14 | Naive_Teffs-U (4 donors) + 5 other subtypes | T.4.Nve.Sp (2 reps) | **AMBIGUOUS** — Calderon has 6 CD4 subtypes |
| 15 | classical monocyte | 15 | Monocytes-U (3 donors) | Mo.6C+II-.Bl (3 reps) | **CLEAN PAIR** |
| 16 | macrophage | 16 | — | MF.PC (2), MF.Alv.Lu (2), MF.RP.Sp (2), +5 subtypes | MISSING (human) |
| 17 | B cell | 17 | Bulk_B-U (4) / Naive_B-U (4) | B.Fo.Sp (3) / B.Sp (2) | **CLEAN PAIR** |
| 18 | luminal epithelial cell (mammary) | 18 | — | — | MISSING (both) |
| 19 | large intestine goblet cell | 19 | — | — | MISSING (both) |
| 20 | enterocyte (large intestine) | 20 | — | — | MISSING (both) |
| 21 | myeloid dendritic cell | 21 | Myeloid_DCs-U (3 donors) | DC.4+.Sp (2) / DC.8+.Sp (2) | **AMBIGUOUS** — ImmGen splits cDC1/cDC2 |
| 22 | monocyte | 22 | Monocytes-U (3 donors) | Mo.6C+II-.Bl (3) | **AMBIGUOUS** — same data as classical monocyte (#15) |
| 23 | natural killer cell | 23 | Mature_NK-U (5 donors) | NK.27+11b+.Sp (2) + subtypes | **CLEAN PAIR** |
| 24 | intermediate monocyte | 24 | — | — | MISSING (both) |
| 25 | mature NK T cell | 25 | — | NKT.Sp (2 reps) | MISSING (human) |
| 26 | adventitial cell | 26 | — | — | MISSING (both) |
| 27 | granulocyte | 27 | — | GN.BM (2), GN.Sp (2) | MISSING (human) |
| 28 | fibroblast | 28 | — | FRC.SLN (LN-specific, 2 reps) | MISSING (human); mouse too specific |
| 29 | bladder urothelial cell | 29 | — | — | MISSING (both) |
| 30 | pancreatic ductal cell | 30 | — | — | MISSING (both) |
| 31 | smooth muscle cell | 31 | — | — | MISSING (both) |
| 32 | hepatocyte | 32 | — | — | MISSING (both) |
| 33 | endothelial cell | 33 | — | BEC.SLN (2) / LEC.SLN (2) | MISSING (human); mouse is LN-specific |
| 34 | non-classical monocyte | 34 | — | Mo.6C-II-.Bl (3 reps) | MISSING (human) |
| 35 | CD8+ alpha-beta T cell | 35 (most rigid) | CD8pos_T-U (4 donors) | T.8.Nve.Sp (2 reps) | **CLEAN PAIR** |

---

## Section 4: Status Summary

### CLEAN PAIRS (n = 4)

| Our Type | Rank | Residual | Calderon Match | ImmGen Match |
|----------|------|----------|----------------|--------------|
| CD8+ alpha-beta T cell | 35 | 5.376 | CD8pos_T-U (4 donors) | T.8.Nve.Sp (2 reps) |
| B cell | 17 | 9.846 | Bulk_B-U (4 donors) | B.Fo.Sp (3 reps) |
| natural killer cell | 23 | 8.687 | Mature_NK-U (5 donors) | NK.27+11b+.Sp (2 reps) |
| classical monocyte | 15 | 9.929 | Monocytes-U (3 donors) | Mo.6C+II-.Bl (3 reps) |

### AMBIGUOUS PAIRS — Requiring Advisor Decision (n = 3 independent)

**1. CD4-positive, alpha-beta T cell (rank 14, residual 10.212)**
- **Calderon:** 6 CD4+ subtypes available unstimulated: Naive_Teffs (4 donors),
  Effector_CD4pos_T (4), Memory_Teffs (4), Follicular_T_Helper (5),
  Regulatory_T (4), Memory_Tregs (4)
- **ImmGen:** T.4.Nve.Sp (2 reps) — naive CD4+ only
- **Issue:** Our Tabula "CD4-positive, alpha-beta T cell" is a broad label including all
  activation states. Mapping options: (a) use Naive_Teffs ↔ T.4.Nve.Sp (both naive,
  cleanest biological match); (b) aggregate/average across Calderon subtypes.
- **Recommendation for advisor:** Option (a) using naive-only is the cleanest. This
  slightly misrepresents our Tabula type (which pools naive + memory + effector) but
  ensures comparable sorting between species.

**2. Myeloid dendritic cell (rank 21, residual 8.899)**
- **Calderon:** Myeloid_DCs-U (3 donors) — likely BDCA1+ cDC2
- **ImmGen:** DC.4+.Sp (2 reps, = cDC2) and DC.8+.Sp (2 reps, = cDC1)
- **Issue:** Calderon's "Myeloid_DCs" maps to cDC2. ImmGen has both cDC1 and cDC2. Our
  Tabula "myeloid dendritic cell" likely includes both subtypes.
- **Recommendation for advisor:** Use DC.4+.Sp (cDC2) ↔ Myeloid_DCs-U (cDC2) for
  matched subtype. Or average DC.4+ and DC.8+ for the mouse side to better represent
  our broad "myeloid DC" label. Either way, the mapping is defensible.

**3. Plasma cell (rank 12, residual 11.084)**
- **Calderon:** Plasmablasts-U (2 donors only)
- **ImmGen:** B.PC.Sp (2 reps), B.PC.BM (2 reps) — genuine plasma cells
- **Issue:** Plasmablasts are biologically distinct from plasma cells. Plasmablasts
  are proliferating, recently differentiated, short-lived antibody-secreting cells.
  Plasma cells are long-lived, non-dividing, terminally differentiated. Chromatin
  accessibility profiles may differ substantially.
- **Recommendation for advisor:** Exclude unless willing to accept the biological
  mismatch caveat. Only 2 donors on the human side is also a data quality concern.

### Non-independent ambiguous entries (do not count as additional n)

- **monocyte (rank 22):** Maps to same Calderon/ImmGen data as classical monocyte
  (rank 15). Cannot count both — choose one label. Classical monocyte is cleaner.
- **T cell (rank 7):** Generic label that overlaps with CD4+ T (#14) and CD8+ T (#35).
  No independent data — subsumed by the subtype mappings.

---

## Section 5: Power Projection

### Scenario analysis

| Scenario | n | Types Included | |ρ| threshold (p<0.05) |
|----------|---|----------------|------------------------|
| Clean only | 4 | CD8+ T, B, NK, class. mono. | 1.000 — **IMPOSSIBLE** |
| + CD4+ T resolved | 5 | + CD4+ T (naive) | 0.900 |
| + myeloid DC resolved | 6 | + myeloid DC (cDC2) | 0.829 |
| + plasma cell resolved | 7 | + plasma cell | 0.786 |

### Pre-registered T3-E thresholds (from Step 1 report)

- ρ ≥ 0.50: POSITIVE — chromatin conservation predicts rigidity
- ρ < 0.35: triggers 8th null closure

### Reference power thresholds (from task specification)

| n | |ρ| threshold |
|---|-------------|
| 8 | 0.738 |
| 9 | 0.683 |
| 10 | 0.648 |
| 12 | 0.591 |
| 15 | 0.521 |

### Assessment

**At n = 4 (clean pairs only): UNTESTABLE.** Spearman with n=4 cannot achieve p<0.05
at any ρ value. The pre-registered positive threshold (ρ ≥ 0.50) requires n ≥ 10 to
be statistically significant.

**At n = 6 (best case with 2 ambiguous resolved): MARGINALLY TESTABLE.** Requires
|ρ| ≥ 0.829, which demands a near-perfect monotonic relationship. The pre-registered
ρ ≥ 0.50 threshold falls below the significance boundary.

**At n = 7 (all 3 ambiguous resolved): MARGINALLY TESTABLE.** Requires |ρ| ≥ 0.786.
Still above the pre-registered positive threshold.

**Comparison to Step 1 estimate:** Step 1 projected "best-case n = 6" and "realistic
n = 6." This assessment finds **realistic n = 4-6**, with 4 clean and 2-3 ambiguous.
The Step 1 estimate included neutrophils (not in Calderon) and did not fully assess
mapping ambiguity for CD4+ T and myeloid DC.

**Pre-registered ρ ≥ 0.50 is NOT statistically testable at any achievable n (4-7).**
The minimum n for ρ ≥ 0.50 to reach p < 0.05 is approximately n = 10. At n = 6-7,
only extreme effects (ρ ≥ 0.83+) would be significant.

---

## Section 6: Structural Problem — Immune-Only Coverage

Both Calderon and ImmGen are **immune cell atlases**. Our 35-type Procrustes set includes
17 non-immune types (stromal, epithelial, hepatocyte, endothelial, fibroblast, smooth
muscle, etc.) that are absent from both datasets. Even among the 18 immune/blood types
in our set, only 4-7 have matched human-mouse data.

This is not a gap that can be filled by additional samples — it is a structural limitation
of using immune-focused ATAC-seq atlases for a test that spans the full cell type spectrum.

The rigidity ranking spans from stromal (most flexible, rank 1, residual 16.98) to CD8+ T
(most rigid, rank 35, residual 5.38). Our 4 clean pairs span ranks 15-35 (classical
monocyte through CD8+ T) — the rigid half of the ranking only. The flexible half
(ranks 1-14) has zero representation. This range restriction severely limits the
Spearman correlation's ability to detect a relationship.

---

## Section 7: ImmGen Cell Types NOT in Our 35-Type Set

The following ImmGen populations have no match in our Procrustes set but could be
relevant if the 35-type set is expanded:

| ImmGen Type | Description | Potential Calderon Match? |
|-------------|-------------|-------------------------|
| Tgd.Sp | Gamma-delta T cell | Yes: Gamma_delta_T-U (4 donors) |
| Treg.4.25hi.Sp | CD4+ Treg | Yes: Regulatory_T-U (4 donors) |
| DC.pDC.Sp | Plasmacytoid DC | Yes: pDCs-U (3 donors) |
| ILC2.SI | ILC2 | No |
| ILC3.*.SI | ILC3 subtypes | No |
| B.MZ.Sp | Marginal zone B | No (Calderon = PBMCs, no MZ B) |
| B1b.PC | B1b cells | No |

**Gamma-delta T, Treg, and pDC** are cross-species pairs with both Calderon and ImmGen
data, but they are not in our current 35-type Procrustes ranking. Including them would
require recomputing Procrustes with an expanded type set — a substantial reanalysis.

---

## Section 8: Failed Retrievals and Issues

1. **GSE131651 accession error:** Task specified GSE131651 as ImmGen ATAC. This is
   NSD2/multiple myeloma (PMID 31649247). Correct accession found via NCBI E-utilities
   search for "ImmGen ATAC mouse" → GDS ID 200100738 → **GSE100738**.

2. **Calderon paper full-text:** PMC article (PMC6426066) returned HTTP 410 (Gone).
   Nature article page returned HTTP 303. Cell type list confirmed from GEO sample
   titles instead (25 unique types, all accounted for).

3. **ImmGen peak files:** GEO supplementary files include only BigWig (.bw) and RAW.tar
   archives. No pre-called peak BED/narrowPeak files are deposited. Peak calling from
   raw data (BAM/FastQ) will be required for the actual analysis. This is a non-trivial
   preprocessing step not anticipated in Step 1.

4. **Calderon peak format:** A consensus peak count matrix (`GSE118189_ATAC_counts.txt.gz`)
   is available. This provides pre-computed accessibility counts at consensus peaks across
   all samples/conditions, potentially eliminating the need to download individual BAM
   files for the Calderon side. Format/genome assembly needs verification at download.

5. **ImmGen publication:** The correct citation is Yoshida et al. 2019 **Cell** 176:897-912
   (PMID 30686579), not Yoshida et al. 2021 Science as stated in the task. The second
   linked PubMed ID (31541153) is the ImmGen sexual dimorphism paper (Gal-Oz et al. 2019
   Nature Communications).

---

## Section 9: Summary and Recommendation

| Metric | Value |
|--------|-------|
| Clean cross-species pairs | **4** (CD8+ T, B cell, NK cell, classical monocyte) |
| Ambiguous pairs (advisor decision) | **3** (CD4+ T, myeloid DC, plasma cell) |
| Best-case n | **7** (if all ambiguous resolved) |
| Realistic n | **5-6** (CD4+ T and myeloid DC resolvable; plasma cell questionable) |
| |ρ| threshold at n=6 | 0.829 |
| Pre-registered ρ≥0.50 testable? | **NO** — requires n≥10 |
| Rigidity rank coverage | Ranks 12-35 only (rigid half); no flexible types |
| Non-immune types covered | **0 of 17** |

**The Calderon + ImmGen pathway produces n = 4-7, concentrated in the rigid half of
the ranking, with the pre-registered positive threshold statistically untestable.**

This is a weaker position than Step 1 projected. The two key downgrades are:
1. Calderon does not have neutrophils (removing one assumed clean pair)
2. Mapping ambiguity for CD4+ T and myeloid DC was not previously assessed

**Advisor decisions needed:**
1. Accept CD4+ T mapping (naive subset in both species)?
2. Accept myeloid DC mapping (cDC2 in both species)?
3. Accept or reject plasma cell (plasmablast ≠ plasma cell)?
4. Given n = 5-7 and untestable positive threshold: proceed with T3-E as
   directional/descriptive analysis, or seek additional data sources to reach n ≥ 10?
5. Consider expanding the 35-type Procrustes set to include gamma-delta T, Treg, and
   pDC — which have matched Calderon/ImmGen data but no current rigidity scores?

---

## Raw Metadata Checkpoints

Saved to `output/validation/t3e_feasibility/tier2_raw/`:
- `calderon_geo_metadata.txt` — full Calderon cell type list with donor counts
- `immgen_geo_metadata.txt` — full ImmGen cell type list with replicate counts and
  activation status annotations
