# Macaque Single-Cell RNA-seq Atlas Survey

**Date:** 2026-04-05
**Purpose:** Research survey of macaque scRNA-seq atlas datasets for potential
cross-species transcriptomic comparison in CellWarp.
**Status:** Research only -- no code written, no data downloaded.

---

## 1. NHPCA -- Non-Human Primate Cell Atlas (Han et al. 2022)

| Field | Value |
|---|---|
| **Species** | *Macaca fascicularis* (cynomolgus macaque) |
| **Cells** | 1,144,706 (after QC: >=500 genes, <=10% mito) |
| **Tissues** | 45 organs/tissues across all major physiological systems |
| **Cell types** | 113 major cell types |
| **Sequencing technology** | **DNBelab C4** (BGI droplet-based platform); predominantly snRNA-seq, scRNA-seq for select tissues (colon, duodenum, spleen, stomach, lymph node, bone marrow, PBMC) |
| **Data format** | Count matrices downloadable from NHPCA portal |
| **Download** | Raw: CNGB accession CNP0001469; Count matrices: https://db.cngb.org/nhpca/download; Interactive: https://db.cngb.org/nhpca/ |
| **Cell Ontology terms** | Not explicitly stated; uses custom annotation nomenclature based on marker genes |
| **Ortholog mapping** | Not pre-computed; would need BioMart mapping from M. fascicularis gene IDs |
| **Publication** | Han L, Wei X, Liu C, et al. "Cell transcriptomic atlas of the non-human primate Macaca fascicularis." *Nature* 604, 723-731 (2022). |
| **Donors** | 8 animals (5 male, 3 female) |

**CellWarp compatibility notes:**
- **TECHNOLOGY RISK: HIGH.** DNBelab C4 is same class of concern as MCA/Microwell-seq
  failure (DECISION-104, p=0.542). Median ~1,324 genes/cell vs ~2,000-5,000 for 10x.
  Already flagged in DECISION-099 and DECISION-123 as technology mismatch risk.
- Best used as VALIDATION/FALLBACK, not primary atlas (per existing project strategy).
- Broadest tissue coverage of any macaque atlas (45 tissues).
- Code: https://github.com/single-cell-BGI/NHPCA

---

## 2. Monkey Atlas (Qu et al. 2022)

| Field | Value |
|---|---|
| **Species** | *Macaca fascicularis* (cynomolgus macaque) |
| **Cells** | 174,233 (scRNA-seq) + 66,566 (scATAC-seq) = ~240,799 total |
| **Tissues** | 16 (trachea, spleen, stomach, kidney, uterus, tongue, testis, muscle, lung, liver, heart, colon, breast, bladder, adipose, aorta) |
| **Cell types** | 40 transcriptionally distinct clusters (17 major types with subtypes) |
| **Sequencing technology** | **10x Genomics Chromium** Single Cell 3' Solution (both scRNA-seq and scATAC-seq) |
| **Data format** | Count matrices; standard 10x format (barcodes/features/matrix) |
| **Download** | GEO: GSE196792 (scRNA-seq), GSE196791 (scATAC-seq); Zenodo: 10.5281/zenodo.5881495; Interactive: https://biobigdata.nju.edu.cn/MonkeyAtlas/ |
| **Cell Ontology terms** | Not explicitly stated; uses standard cell type nomenclature |
| **Ortholog mapping** | Yes -- paper reports **12,971 one-to-one orthologous genes** between macaque and human (also used for cross-species comparison with mouse) |
| **Publication** | Qu J, et al. "A reference single-cell regulomic and transcriptomic map of cynomolgus monkeys." *Nature Communications* 13, 4069 (2022). |
| **Donors** | 1 male + 1 female, 4 years old |

**CellWarp compatibility notes:**
- **TECHNOLOGY: COMPATIBLE.** 10x Chromium -- same platform as Tabula Sapiens/TMS.
- Already downloaded in CellWarp (GSE196792, 230,882 cells across 20 samples).
- LOW donor count (1-2 animals) -- non-immune centroids less robust (flagged in DECISION-123).
- Covers hepatocyte, endothelial, and other non-immune types needed for CellWarp.
- Already integrated into CellWarp macaque pipeline as part of Strategy B.
- Includes scATAC-seq data (potential future use for chromatin analysis).

---

## 3. RIRA -- Rhesus Immune Reference Atlas (Mahyari et al. 2025)

| Field | Value |
|---|---|
| **Species** | *Macaca mulatta* (rhesus macaque) |
| **Cells** | 426,664 (after QC) |
| **Tissues** | 7 (PBMC, spleen, lymph nodes, liver, lung, bone marrow, mesenteric lymph node) |
| **Cell types** | Multiple immune subsets: T_NK (287,807), B cell (84,412), Myeloid (40,601), Non-Immune (10,195), Erythrocyte (6,653), Unknown (167,180) |
| **Sequencing technology** | **10x Genomics 5' capture** (v2/HT chemistry) + CITE-seq + TCR sequencing |
| **Data format** | RDS (R sparse matrix), convertible to MTX |
| **Download** | GEO: GSE277821; BioProject: PRJNA1163395; GitHub: BimberLab/RIRA |
| **Cell Ontology terms** | Uses established immunologic definitions anchored to surface protein markers (CITE-seq) and canonical lineage genes |
| **Ortholog mapping** | Gene reference: MMul_10 (GCF_003339765.1), NCBI gene build 103. Uses HGNC-like symbols. |
| **Publication** | Mahyari E, et al. "Enhanced interpretation of immune cell phenotype and function through a rhesus macaque single-cell atlas." *Cell Genomics* 5(5):100849 (2025). |
| **Donors** | **47 animals** (ages 2-7 years) -- strongest donor diversity of any macaque atlas |

**CellWarp compatibility notes:**
- **TECHNOLOGY: COMPATIBLE.** 10x Genomics 5' capture.
- Already downloaded in CellWarp (596,848 cells, 34,606 genes).
- Excellent for immune cell types (CD8+ T, CD4+ T, B cells, macrophages).
- 47-donor diversity is outstanding for centroid robustness.
- 167,180 "Unknown" cells (28%) -- recommended EXCLUDE from centroid computation.
- Immune-only; needs Qu et al. for non-immune types (hepatocyte, endothelial).
- Already integrated into CellWarp macaque pipeline as part of Strategy B.

---

## 4. Tabula Microcebus (Tabula Microcebus Consortium 2025)

| Field | Value |
|---|---|
| **Species** | *Microcebus murinus* (gray mouse lemur) -- PRIMATE, not macaque |
| **Cells** | ~226,000 |
| **Tissues** | 27 organs |
| **Cell types** | >750 molecularly distinct cell types (768 with tissue-of-origin separation) |
| **Sequencing technology** | Mixed: **10x Chromium** (droplet-based) + **Smart-seq2** (plate-based) |
| **Data format** | **h5ad** (primary); also .mat (MATLAB); conversion scripts for R/Seurat |
| **Download** | Figshare: https://figshare.com/articles/dataset/Tabula_Microcebus_v1_0/14468196; Main file: LCA_complete_wRaw_toPublish.h5ad; Project: https://figshare.com/projects/Tabula_Microcebus/112227 |
| **Cell Ontology terms** | Expert-curated annotations with computational clustering |
| **Ortholog mapping** | Cross-species comparison with human and mouse performed in companion paper. Specific 1:1 ortholog count not published; Ensembl has M. murinus genome but ortholog depth likely lower than macaque due to ~75Myr divergence. |
| **Publication** | Tabula Microcebus Consortium. "A molecular cell atlas of mouse lemur, an emerging model primate." *Nature* (2025). Companion: "Mouse lemur cell atlas informs primate genes, physiology and disease." *Nature* (2025). |
| **Donors** | 4 animals, clinically and histologically characterized |

**CellWarp compatibility notes:**
- PRIMATE but NOT macaque -- mouse lemur diverged ~75Myr from humans (vs ~25Myr for macaque).
- 10x data subset is technology-compatible; Smart-seq2 subset needs separation.
- h5ad format is directly compatible with CellWarp pipeline.
- 750+ cell types is extraordinary annotation depth.
- Ortholog depth unknown but likely intermediate between macaque (~15-19K) and zebrafish (~9.7K).
- Interesting as a FOURTH species (after macaque) for evolutionary distance curve, but NOT a substitute for macaque.
- Data use restrictions: "not publish analyses of genes, cell types or transcriptomic data on a whole atlas or tissue scale prior to initial publication" -- **check if still active** (papers published July 2025).

---

## 5. CELLxGENE Census -- Macaca mulatta

| Field | Value |
|---|---|
| **Species** | *Macaca mulatta* (rhesus macaque) |
| **Cells** | 1,525,089 total (1,372,898 unique) |
| **Tissues** | 29 tissues, 2 tissue general categories |
| **Cell types** | 54 annotated cell types |
| **Assays** | 2 assay types |
| **Technology composition** | **~88% sci-RNA-seq3 brain data (BICCN)** -- dominated by brain tissue |
| **Data format** | TileDB-SOMA (accessible via cellxgene_census Python API); sliceable to AnnData/h5ad |
| **Download** | Via cellxgene_census Python API: `cellxgene_census.open_soma()` |
| **Cell Ontology terms** | Yes -- Census schema requires Cell Ontology annotations |
| **Ortholog mapping** | Census provides standardized gene annotations |
| **Availability** | Census LTS release 2025-11-08 |

**CellWarp compatibility notes:**
- **CRITICALLY LIMITED:** 88% of data is sci-RNA-seq3 brain tissue from BICCN.
- sci-RNA-seq3 is a low-sensitivity combinatorial indexing protocol -- same class of
  concern as Microwell-seq and DNBelab C4.
- Non-brain tissue coverage is minimal.
- Already assessed: "RIRA NOT in Census (Nov 2025 snapshot). Census cannot substitute
  for direct GEO download."
- Census is useful for QUERYING metadata, not as a primary macaque data source.
- Other primates in Census: *Callithrix jacchus* (2.28M cells), *Pan troglodytes* (158K cells, dorsolateral prefrontal cortex only).

---

## 6. Cross-Study Multi-Organ Atlas (2026 preprint)

| Field | Value |
|---|---|
| **Species** | *Macaca fascicularis* (cynomolgus macaque) |
| **Studies integrated** | 30 publicly available studies |
| **Tissues/Organs** | 57 anatomical regions, 43 organs, 14 physiological systems |
| **Cell types** | Not specified in search results (uses UCE foundation model annotation) |
| **Sequencing technology** | Mixed (integrates multiple technologies across 30 studies) |
| **Data format** | Not yet assessed (preprint) |
| **Download** | bioRxiv preprint; data availability TBD |
| **Cell Ontology terms** | Uses Universal Cell Embeddings (UCE) for cross-species harmonization against Tabula Sapiens V2 |
| **Publication** | "A Cross-Study Multi-Organ Cell Atlas of Macaca fascicularis Informed by Human Foundation Model Annotation: A Resource for Translational Target Assessment." *bioRxiv* (2026). |

**CellWarp compatibility notes:**
- VERY recent preprint (March 2026) -- not yet peer-reviewed.
- Integrative atlas across 30 studies means MIXED TECHNOLOGIES (includes DNBelab C4, 10x, possibly others).
- UCE annotation harmonized with Tabula Sapiens V2 is interesting for CellWarp compatibility.
- Could be useful as a validation resource once published.
- Technology heterogeneity is a concern for CellWarp's sensitivity to protocol differences.

---

## 7. Other Individual Macaque Tissue Atlases

### Allen Brain / BICCN Macaque Brain Atlas
- *Macaca mulatta*, brain only, 4.2M cells, multi-omic (scRNA-seq + snATAC-seq)
- Science Advances (2023): Chiou et al.
- Brain-only -- not useful for CellWarp multi-tissue comparison
- Uses sci-RNA-seq3 (low sensitivity)

### Macaque Cardiopulmonary Aging Atlas
- *Macaca fascicularis*, heart + lung only
- Cell Research (2020)
- Single-tissue pair -- limited utility

---

## Ortholog Analysis

### Human-Macaque One-to-One Orthologs

| Source | Count | Species Pair | Method |
|---|---|---|---|
| CellWarp DECISION-098 | **19,123** | Human-Macaque (combined M. fascicularis + M. mulatta) | BioMart (cached in project) |
| Qu et al. 2022 | **12,971** | Human-M. fascicularis | BioMart |
| PMC11017088 (2024) | **~15,000** | Human-Macaque | BioMart, Ensembl v.110 |
| Yan et al. 2011 Nat Biotech | **14,978** | Human-Macaque (3 species alignment) | Genome alignment |
| CellWarp 3-way intersection | **15,028** | Human-Macaque-Mouse | BioMart (87.4% of human-mouse space) |

**Key finding:** The 19,123 count from CellWarp's own BioMart query is the highest,
likely reflecting the most recent Ensembl release. The ~15,000 figure from Qu et al.
and the 2024 paper used earlier Ensembl versions. Variation is expected across Ensembl
releases as gene models are updated.

### Comparison with Human-Mouse Orthologs

| Species Pair | 1:1 Orthologs | Source |
|---|---|---|
| Human-Mouse | **16,959** (shared in CellWarp) / **17,187** (BioMart raw) | CellWarp DECISION-012 |
| Human-Mouse (Ensembl Compara) | ~15,893 | Literature (80% of human protein-coding genes) |
| Human-Macaque | **19,123** | CellWarp DECISION-098 |
| Human-Macaque-Mouse (3-way) | **15,028** | CellWarp PROGRESS.md |

### Does Macaque Pass the CellWarp >12,000 Minimum?

**YES -- STRONGLY.** At 19,123 human-macaque 1:1 orthologs, macaque exceeds the 12,000
minimum by 59%. Even the three-way intersection (human-macaque-mouse) at 15,028 exceeds
the threshold. This is the deepest ortholog space of any species evaluated for CellWarp,
as expected given the ~25Myr divergence time (vs ~80Myr for human-mouse, ~450Myr for
human-zebrafish).

### Macaque-Mouse Orthologs

No specific published count was found for M. mulatta - M. musculus 1:1 orthologs.
However, given that human-macaque is ~19K and human-mouse is ~17K, the macaque-mouse
count should be in the ~15-17K range (estimated from the three-way intersection of 15,028
plus additional genes in the macaque-mouse pair not shared with human).

### Ortholog Mapping Tools

- **Ensembl BioMart** is the standard tool. CellWarp already uses pybiomart for queries.
  Filter: `homology_type = "ortholog_one2one"`.
- Dataset names: `mmulatta_gene_ensembl` (M. mulatta), `mfascicularis_gene_ensembl`
  (M. fascicularis), `hsapiens_gene_ensembl` (human).
- Alternative: `orthologsBioMART` Python/R wrapper for simpler queries.
- CellWarp already has cached ortholog table at `data/phase1/orthologs_human_mouse.csv`
  with 17,187 pairs (16,959 after intersection with expressed genes).

---

## Summary Assessment for CellWarp

### Already in the pipeline (Strategy B, DECISION-123):
1. **RIRA** (M. mulatta, 10x, immune types, 47 donors) -- PRIMARY for immune
2. **Qu et al.** (M. fascicularis, 10x, non-immune types, 1-2 donors) -- PRIMARY for non-immune
3. **NHPCA** (M. fascicularis, DNBelab C4, 45 tissues) -- FALLBACK/VALIDATION

### New resources identified in this survey:
4. **CELLxGENE Census M. mulatta** -- NOT useful (88% brain sci-RNA-seq3)
5. **Cross-Study Multi-Organ Atlas (2026 preprint)** -- potential future validation resource
6. **Tabula Microcebus** -- interesting as 4th species, not macaque substitute
7. **RIRA is now published** (Cell Genomics May 2025) -- was preprint at time of DECISION-123

### Key risks already documented:
- Species mixing (M. mulatta + M. fascicularis) -- DECISION-123-AMENDMENT abort criteria
- Qu et al. low donor count (1-2) -- LOW-CONFIDENCE flag
- DNBelab C4 technology mismatch for NHPCA -- same class as MCA failure
- 28% Unknown cells in RIRA -- EXCLUDE from centroids

### Conclusion:
The existing CellWarp Strategy B (RIRA + Qu et al.) remains the best approach.
No new dataset discovered in this survey changes the recommendation. The 2026
cross-study atlas is interesting but too heterogeneous in technology. Tabula
Microcebus is a candidate for future 4th-species work but is not a macaque
substitute and has unknown ortholog depth.
