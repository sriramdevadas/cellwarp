# Independent 10x Mouse Atlas Feasibility Search — Candidate List

**Date:** 2026-03-15
**Purpose:** Identify viable mouse scRNA-seq atlases for T1-A cross-species Procrustes
replication after MCA microwell-seq failure (p=0.542, obs/null=1.003).
**Requirement:** 10x Chromium (or comparable high-sensitivity protocol), adult mouse,
multi-tissue, >=10/35 cell types, independent of Tabula consortium.

---

## CRITICAL FINDING: NO DROP-IN REPLACEMENT EXISTS

No single independent multi-tissue adult mouse 10x Chromium scRNA-seq atlas exists
outside the Tabula consortium. This is a genuine gap in the field. The only comprehensive
multi-tissue adult mouse atlases are:
- Tabula Muris / Tabula Muris Senis (10x + Smart-seq2) — our reference, excluded
- Mouse Cell Atlas (microwell-seq) — already failed T1-A

All other candidates use non-10x protocols or cover too few tissues. The search covered:
GEO, PubMed, Single Cell Portal (Broad), CELLxGENE Census, ArrayExpress, 10x Genomics
datasets, Parse Biosciences, ENCODE, ImmGen, CZI-funded projects, DISCO database,
Jackson Laboratory, IMPC, commercial tissue banks, and spatial transcriptomics companions.

---

## RANKED CANDIDATE LIST

### RANK 1: PanSci (Cao Lab, Science 2025) — STRONGEST OVERALL

| Field | Value |
|---|---|
| Paper | Zhang Z et al. "A panoramic view of cell population dynamics in mammalian aging" |
| Journal | Science 387(6731), Jan 2025 |
| DOI | 10.1126/science.adn3949 |
| PMID | 39607904 |
| Protocol | **EasySci (snRNA-seq combinatorial indexing) — NOT 10x** |
| Tissues (13) | Kidney, lung, heart, **liver**, skeletal muscle, stomach, brown adipose, inguinal WAT, perigonadal WAT, ileum, colon, jejunum, duodenum |
| Missing tissues | **Spleen, bone marrow, blood, pancreas, bladder, mammary** |
| Mouse strain/age | C57BL/6 wild-type, 6/12/23 months (sex-balanced) + Rag1-KO and SCID |
| Cell count | **21,786,931 nuclei** (all conditions) |
| Median UMI/cell | **~1,040** (low — 2-5x less than 10x) |
| Cell types | 239 organ-specific main types, 3,925 sub-clusters |
| GEO | **GSE247719** |
| Data format | h5ad (per-tissue files, 1.8-13.8 GB each; all-cells 29.8 GB) |
| Zenodo | Code only (10.5281/zenodo.13755852) |
| UCSC Browser | mouse-pansci.cells.ucsc.edu |
| Azimuth ref | Available (Mouse-PanSci) |
| Author overlap | **NONE** (Rockefeller, UCSC, NYU — zero Tabula/CZ Biohub) |
| Access method | **GEO free** |
| Est. cost | $0 |
| Est. types | **~20-25/35** |

**Cell type coverage assessment:**

| Our Type | PanSci Match | Tissue |
|---|---|---|
| B cell | Lymphoid cells_B cells | Multi-tissue immune |
| CD4+ T cell | CD4+ naive T cells (sub-cluster) | Multi-tissue immune |
| CD8+ T cell | CD8+ cytotoxic T cells (sub-cluster) | Multi-tissue immune |
| Endothelial | Endothelial cells, Vascular EC, Lymphatic EC | All tissues |
| Hepatocyte | Hepatocytes | Liver |
| Macrophage | Myeloid cells_Alveolar mac, Interstitial mac | Multi-tissue |
| Fibroblast | Present in multiple tissues | Heart, lung, intestine |
| Smooth muscle | Present | Multiple |
| Epithelial | Present | Lung, intestine, stomach |
| NK cell | Present in immune compartment | Likely from tissue-resident |
| Plasma cell | Lymphoid cells_Plasma cells | Likely multi-tissue |
| Enterocyte | Present | Colon, ileum, jejunum, duodenum |
| Goblet cell | Present | Intestinal tissues |
| Stromal | Present | Multiple |
| MSC | Possible in adipose | WAT tissues |
| Neutrophil | Present in myeloid | Multi-tissue |
| Monocyte | Present in myeloid | Multi-tissue |
| Dendritic cell | Present in myeloid | Multi-tissue |
| Basal cell | Possible in stomach/intestine | — |
| HSC/precursor | **ABSENT** (no bone marrow) | — |
| Pancreatic types | **ABSENT** (no pancreas) | — |
| Bladder urothelial | **ABSENT** (no bladder) | — |
| Mammary luminal | **ABSENT** (no mammary) | — |

**Strengths:**
- Massive scale (21.8M cells) ensures huge per-type counts for centroid estimation
- Liver included (hepatocyte — most biologically interpretable residual)
- Intestinal tissues (enterocyte, goblet — hard to find elsewhere)
- h5ad format directly compatible with our pipeline
- No batch correction applied — transparent
- Multiple ages allow age-matching to TMS if desired
- Immune cell subset h5ad files pre-extracted (B_cell, T_cell, Myeloid — 2-3 GB each)

**Concerns:**
- **Protocol mismatch:** EasySci is snRNA-seq (nuclei only, not whole cells). Loses
  cytoplasmic transcripts. Different technical noise profile from 10x scRNA-seq.
  Same class of concern as MCA microwell-seq.
- **Low UMI depth:** Median 1,040 UMI/cell vs 2,000-5,000 for 10x. But massive cell
  numbers may compensate via centroid averaging.
- **Missing spleen/BM/blood:** Immune cell types come from tissue-resident populations
  only, not from dedicated immune organs. T/B/NK cell counts may be lower.
- **CD4/CD8 T cell split:** Available only at sub-cluster level, not main type.
  Need T_cell h5ad file and sub-cluster annotation column.
- **Custom annotations:** NOT Cell Ontology terms. Mapping required but straightforward.

**KEY QUESTION:** Does the EasySci sensitivity limitation cause the same centroid
compression problem as MCA microwell-seq? The MCA scaling factor was 0.267 (3.7x
compression). If PanSci shows similar compression, the Procrustes signal may again
be below detection threshold. The empirical test is: compute PanSci centroids and
check the scaling factor. If scaling ~0.5-1.5, proceed; if scaling <0.3, same problem.

---

### RANK 2: Sun et al. 2023 (CAS, Innovation) — BEST 10x PROTOCOL MATCH

| Field | Value |
|---|---|
| Paper | Sun S et al. "A single-cell transcriptomic atlas of exercise-induced anti-inflammatory and geroprotective effects across the body" |
| Journal | Innovation (Cambridge) 4(1):100380, 2023 |
| DOI | 10.1016/j.xinn.2023.100380 |
| PMID | 36747595 |
| Protocol | **10x Chromium 3' v3 (9 tissues) + snRNA-seq (5 tissues)** |
| 10x tissues (9) | Lung, aorta, kidney, **liver**, small intestine, testis, **spleen**, **bone marrow**, peripheral blood |
| snRNA-seq tissues (5) | Brain, cerebellum, spinal cord, heart, skeletal muscle |
| Mouse strain/age | C57BL/6J, 2 months (young) and 16 months (old) |
| Cell count | **507,636 total** (~63K per condition, ~4,500/tissue/condition estimated) |
| Cell types | 101 main types, 305 clusters |
| Data (raw) | GSA CRA007207 (6.7 TB FASTQ) |
| Data (processed) | **OMIX002605** (84 tar files, ~4.5 GB total) |
| GitHub | github.com/wxb1998/Mouse-exercise-Project |
| Author overlap | **NONE** (Chinese Academy of Sciences, Beijing + Altos Labs) |
| Access method | **OMIX open access** (Chinese National Genomics Data Center) |
| Est. cost | $0 |
| Est. types | **~17-20/35** |

**Cell type coverage (from code analysis):**

| Our Type | Sun et al. Match | Tissue | Protocol |
|---|---|---|---|
| B cell | BC | Spleen, BM, liver | 10x |
| CD4+ T cell | CD4_Mem, CD4_Naive | Spleen, BM, blood | 10x |
| CD8+ T cell | CD8_Naive, CD8_Mem, CD8_CTL | Spleen, BM, blood | 10x |
| Endothelial | EC_Liver, EC_Lung, EC_Kidney, EC_Aorta | All | 10x |
| Hepatocyte | Hep | Liver | 10x |
| Macrophage | Mac1, Mac2, Kup, AMac | Multiple | 10x |
| Monocyte | Mono, Mono_BM | BM, blood | 10x |
| NK cell | NK | Multiple | 10x |
| Neutrophil | Neu, ProNeu | BM, spleen | 10x |
| Dendritic cell | mDC, pDC | Multiple | 10x |
| Plasma cell | Pla | Spleen, BM | 10x |
| Fibroblast | Fib_Lung, Fib_Intestine, Fib_Aorta | Multiple | 10x |
| Smooth muscle | SMC | Aorta | 10x |
| Epithelial | Epi_Lung, AT1, AT2 | Lung | 10x |
| HSC/precursor | Progenitor | BM | 10x |
| Cardiac fibroblast | Fib_Heart | Heart | snRNA-seq |
| Granulocyte | Bas, Mast, Neu | Multiple | 10x |
| Enterocyte | Intestine profiled but not named | Small intestine | 10x |
| Basal | NOT annotated | — | — |
| Goblet | NOT annotated | — | — |
| Pancreatic types | **ABSENT** | — | — |
| Bladder urothelial | **ABSENT** | — | — |
| Mammary luminal | **ABSENT** | — | — |
| MSC | NOT explicit | — | — |

**Strengths:**
- **10x Chromium 3' v3 for 9 key tissues** — protocol-matched to TMS 10x component
- **Spleen + bone marrow + blood** — dedicated immune organs (PanSci lacks these)
- All original 6 cell types covered with explicit annotations
- CD4/CD8 T cell split at main annotation level
- Liver, lung, kidney, intestine all profiled with 10x
- Exercise study design means sedentary controls exist as proper baseline
- No author overlap with Tabula

**Concerns:**
- **Data access:** OMIX002605 is a Chinese national repository. Download speeds may be
  slow outside China. File format inside tar archives is unknown without downloading.
- **Cell counts in control arm:** With ~507K total across ~8 conditions × 14 tissues,
  the YC (young sedentary control) arm may have only ~4,500 cells/tissue. Some cell
  types may fall below 500-cell gate.
- **Exercise study, not atlas:** Designed around exercise intervention, not baseline
  characterization. Controls should be adequate but annotation depth may focus on
  exercise-relevant biology.
- **5 tissues use snRNA-seq:** Heart and skeletal muscle are snRNA-seq only. If we
  need cardiac fibroblast, this introduces a protocol mix.
- **Custom annotations:** "Hep", "BC", "CD8_CTL" — not CL ontology. Mapping required.
- **Only male mice** (C57BL/6J males).

**KEY ADVANTAGE OVER PANSCI:** 10x Chromium for the tissues we care about most (liver,
spleen, BM, lung, kidney). This directly addresses the concern that MCA failed due
to protocol mismatch. If Sun et al. 10x data also fails Procrustes, the argument
shifts from "protocol problem" to "result may be Tabula-specific."

---

### RANK 3: Calico/Kimmel 2019 — ONLY TRUE 10x CANDIDATE (BUT TOO NARROW)

| Field | Value |
|---|---|
| Paper | Kimmel JC et al. "Murine single-cell RNA-seq reveals cell-identity- and tissue-specific trajectories of aging" |
| Journal | Genome Research 29(12):2088-2103, 2019 |
| PMID | 31754020 |
| Protocol | **10x Chromium 3' v2** |
| Tissues (3) | Kidney, lung, spleen |
| Mouse strain/age | C57BL/6J, 7-8 months (young) and 22-23 months (old) |
| Cell count | ~55,293 cells |
| GEO | **GSE132901** (21 samples, open access) |
| Website | mca.research.calicolabs.com |
| GitHub | github.com/calico/2019_murine_cell_aging |
| Author overlap | **NONE** (Calico Life Sciences, South San Francisco) |
| Access method | **GEO free** |
| Est. cost | $0 |
| Est. types | **~10-12/35** |

**Strengths:** True 10x Chromium, fully independent, easy GEO access, well-characterized.
**Fatal weakness:** Only 3 tissues. **No liver** (= no hepatocyte, our most interpretable
residual). Missing heart, intestine, bone marrow, blood, pancreas, mammary, bladder.
Cannot reach the minimum 10 types needed for a meaningful Procrustes replication.

**Potential use:** Supplement a primary replication dataset. Combine Calico spleen (immune
types) with another source's liver, etc. But then it's a composite approach (see Rank 5).

**NOTE:** Cell type annotations were trained using a Tabula Muris neural network classifier.
This creates a methodological dependency: annotations are TM-derived even though data is
independent. Would need to re-annotate for clean independence.

---

### RANK 4: Parse 5M Mouse Atlas — LARGE BUT PROTOCOL AND COVERAGE GAPS

| Field | Value |
|---|---|
| Dataset | 5 Million Mouse Single Cell Atlas |
| Source | Parse Biosciences Trailmaker repository |
| URL | parsebiosciences.com/datasets/5-million-mouse-single-cell-atlas-from-7-tissues/ |
| Protocol | **Evercode WT Penta (SPLiT-seq derivative, snRNA-seq) — NOT 10x** |
| Tissues (7) | Brain, colon, eye, heart, kidney, liver, quadriceps muscle |
| Mouse strain/age | Adult (specific strain/age not stated) |
| Cell count | **5,011,382 nuclei**, 211 cell types |
| Author overlap | **NONE** |
| Access method | Free download (Trailmaker account, CC BY-NC 4.0) |
| Est. cost | $0 |
| Est. types | **~12-15/35** |

**Missing tissues:** Spleen, bone marrow, blood, lung, pancreas, mammary, bladder.
**Protocol concern:** SPLiT-seq variant (combinatorial barcoding, snRNA-seq). Different
sensitivity profile from both 10x and microwell-seq.
**License concern:** CC BY-NC 4.0 — may restrict publication in some contexts.
**Not peer-reviewed:** Commercial technology demonstration, not a published study.

---

### RANK 5: Composite 10x Approach — TECHNICALLY FEASIBLE BUT BATCH-CONFOUNDED

Assemble a virtual multi-tissue atlas from 5-6 independent tissue-specific 10x studies:

| Tissue | Best Candidate | GEO | Cells | Protocol |
|---|---|---|---|---|
| Liver | Guilliams et al. 2022 Cell | GSE192742 | ~120K | 10x CITE-seq |
| Spleen | ImmGen (Broad/HMS) | SCP306 | ~10K | 10x 3' |
| Bone marrow | ImmGen | SCP978 | ~varies | 10x 3' |
| Lung | LungMAP CellRef | Multiple | ~40K | 10x (mixed versions) |
| Kidney | Park et al. 2018 Science | GSE107585 | ~58K | 10x |
| Colon | Sirvinskas et al. 2022 | GSE168448 | ~24K | 10x |

**Estimated coverage:** ~18-20/35 types at high confidence + ~10 more at "maybe" level.

**CRITICAL RISK:** Each tissue from a different lab. Tissue-of-origin is perfectly
confounded with batch-of-origin. A reviewer would rightly ask whether any Procrustes
signal reflects genuine biology or the batch structure. This is the same critique the
composite approach is supposed to ADDRESS, not reproduce. **Scientifically weaker than
a single-source atlas.** Only pursue as last resort.

---

## EXCLUDED CANDIDATES

| Dataset | Reason for Exclusion |
|---|---|
| Tabula Muris (original, 2018) | Same consortium as TMS (Quake PI, CZ Biohub) |
| Tabula Microcebus | Wrong species (mouse lemur, a primate) |
| MCA 2.0/3.0 (Han et al.) | Still microwell-seq — same protocol problem as MCA 1.0 |
| Cao et al. 2019 (MOCA) | Embryonic (E9.5-E13.5), sci-RNA-seq (not 10x) |
| Allen Brain Cell Atlas | Brain only — no immune types, no hepatocyte |
| BICCN | Brain only (NIH Brain Initiative) |
| ENCODE4 mouse postnatal | 5 tissues (3 brain), Parse snRNA-seq, too few types |
| Mouse gastrulation atlas | Embryonic (E6.5-E8.5, Pijuan-Sala 2019) |
| ImmGen (standalone) | Immune cells only — no hepatocyte, endothelial, epithelial |
| LifeTime Initiative (EU) | Policy framework, not a dataset |
| IMPC | Embryonic + knockout backgrounds |
| Jackson Laboratory | No multi-tissue atlas published |
| DISCO database | Human-focused aggregator, not independent data |
| Mouse GTEx equivalent | Does not exist |
| BioIVT / commercial tissue banks | No scRNA-seq data produced |
| 10x Genomics reference datasets | Single-tissue demos only (brain, PBMC) |
| Nanostring CosMx | Brain-only mouse data; no companion scRNA-seq |
| Vizgen MERFISH | Used TMS as reference, not independent 10x |
| CZ Biohub (non-Tabula) | All mouse work IS the Tabula consortium |

---

## TABULA AUTHOR OVERLAP CHECK

**Tabula consortium key personnel to exclude:**
- Stephen Quake (CZ Biohub, Stanford) — PI of TM and TMS
- Angela Oliveira Pisco (CZ Biohub) — lead computational biologist
- Spyros Darmanis (CZ Biohub) — experimental lead
- Norma Neff (CZ Biohub) — sequencing lead

**Verification status:**

| Candidate | Overlap? | Verified |
|---|---|---|
| PanSci (Cao lab) | NONE — Rockefeller/UCSC/NYU | Full author list checked |
| Sun et al. 2023 (CAS) | NONE — Chinese Academy of Sciences | Full author list checked |
| Calico/Kimmel 2019 | NONE — Calico Life Sciences | Full author list checked |
| Parse 5M Atlas | NONE — Parse Biosciences | Commercial dataset |
| ImmGen | NONE — Broad/HMS consortium | Known independent consortium |
| Guilliams 2022 (liver) | NONE likely — VIB-UGent, Belgium | **NEEDS MANUAL CHECK** |
| Park 2018 (kidney) | NONE likely — Washington University | **NEEDS MANUAL CHECK** |
| Sirvinskas 2022 (colon) | NONE likely — Fritz Lipmann Inst, Germany | **NEEDS MANUAL CHECK** |
| LungMAP CellRef | NONE likely — NIH/NHLBI consortium | **NEEDS MANUAL CHECK** |

---

## RANKED RECOMMENDATION

### PURSUE FIRST: Sun et al. 2023 (Rank 2)

**Rationale:** The ONLY candidate with genuine 10x Chromium 3' v3 data for the tissues
we need most (liver, spleen, bone marrow, lung, kidney). This directly addresses the
protocol mismatch hypothesis. If 10x data from an independent lab ALSO fails Procrustes,
the argument fundamentally shifts from "MCA is too sparse" to "primary result may be
Tabula-specific." Either outcome is informative.

**Action items:**
1. Download a sample tar file from OMIX002605 to verify data format
2. Check cell counts in YC (young sedentary control) arm per tissue per type
3. Verify >=500 cells/type for >=10 of our 35 types in control arm
4. If feasible, run identical Procrustes pipeline (Sun 10x control → HCA pooled human)

**Risk:** Data access from Chinese repository may be slow/difficult. Cell counts in
control arm may be too thin (~4,500 cells/tissue estimated).

### PURSUE SECOND: PanSci (Rank 1)

**Rationale:** Far more cells (21.8M vs 507K) and better tissue coverage (13 vs 14
organs), but EasySci protocol has the same class of concern as MCA microwell-seq
(low per-cell sensitivity, snRNA-seq only). PanSci is the stronger dataset IF the
Procrustes signal is robust to protocol differences. Since MCA already failed,
PanSci would need to be argued as "higher sensitivity than MCA despite being non-10x."

**Action items:**
1. Download liver h5ad (6 GB) + immune subset h5ad files (~7 GB total)
2. Compute centroids for 6-month wild-type C57BL/6 (age-matched to TMS)
3. Check scaling factor — if ~0.5-1.5, proceed; if <0.3, same problem as MCA
4. If scaling OK, run full Procrustes pipeline

**Risk:** Same sensitivity concern as MCA. Median 1,040 UMI may produce centroid
compression. But 21.8M cells provides massive averaging power.

### PURSUE THIRD: Calico + Composite (Ranks 3+5)

**Rationale:** Calico is the only genuine 10x dataset but covers too few tissues.
Supplementing with tissue-specific 10x studies (liver, intestine, bone marrow) could
reach sufficient coverage. However, this introduces batch confounding.

**Action items:** Only pursue if Ranks 1 and 2 both fail or are infeasible.

---

## STRATEGIC NOTE

The ideal outcome is Sun et al. 10x data producing a clear result (significant or null).
If significant: primary result replicated with independent 10x data (strongest possible
validation). If null: combined with MCA null, suggests primary result IS Tabula-specific
or protocol-specific to the exact TMS processing pipeline (not just "10x" in general).
Either outcome resolves the ambiguity.

PanSci is valuable as a SECOND replication — if Sun et al. 10x succeeds, PanSci tests
whether the signal generalizes to yet another protocol. If Sun et al. 10x fails, PanSci
failure would confirm technology dependence; PanSci success would suggest the signal
exists but requires massive cell counts to detect.
