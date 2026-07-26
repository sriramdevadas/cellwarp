# Large Data Sources — Not in Git
These directories are excluded from git due to size. Some are downloaded
programmatically by the scripts that use them; others must be fetched by hand
from the sources below. Inputs obtained through CZ CELLxGENE Census are not
listed here; the Census release they are pinned to is recorded in
docs/declarations.txt and in README.md.

## data/replication/sun2023/ (1.6GB)
Source: OMIX002605 (CNCB; https://ngdc.cncb.ac.cn/omix/release/OMIX002605)
GSA mirror: CRA007207
Format: .tar.gz Cell Ranger filtered matrices
Protocol: 10x Chromium 3' v3, 8 young (YC) tissues
Required path: data/replication/sun2023/<tissue>/filtered_feature_bc_matrix/
Used by: scripts/02_replication_sun2023.py

## data/replication/pansci/ (13GB)
Source: GEO GSE247719 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247719)
Format: h5ad
Required path: data/replication/pansci/pansci.h5ad
Used by: scripts/05_replication_pansci.py

## data/mouse_lemur/

Source: Tabula Microcebus consortium (Ezran et al., 2025; Nature 644, 173-184)
DOI: 10.1038/s41586-025-09113-9
Access: CZ CELLxGENE Discover
  collection a137437b-d284-4a27-b1e9-36958a8f92c1
  dataset a392ab34-9016-4f48-b45d-5b3a9cfa39fe
Download date: 2026-04-05 (CELLxGENE Discover does not version-pin)
Assay filter: 10x 3' v2 only (95% of atlas; manuscript MC:389)
Required path: data/mouse_lemur/tabula_microcebus_LCA_complete.h5ad (244,081 cells; gitignored, downloaded at first run)
Tracked metadata anchor: data/replication/tabula_microcebus_metadata.csv
Cross-consistency test: tests/test_tabula_microcebus_metadata.py (8 assertions linking CSV ↔ manuscript ↔ feasibility doc)
Used by: analysis/mouse_lemur/ pipeline

## data/replication/raw_downloads/ (2.9GB)
Contents: Raw Cell Ranger outputs from replication runs
Reproducible from accession numbers above

## data/ucsc/ (12GB)
Source: UCSC Genome Browser hg38 phastCons tracks
Required files (place at data/ucsc/<filename>):
  - phastCons_placental.bw — placental-mammal 20way (PRIMARY t3e track):
    https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons20way/hg38.phastCons20way.bw
    (rename to phastCons_placental.bw on disk)
  - phastCons100way.bw — 100-way vertebrate (sensitivity track):
    https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.phastCons100way.bw
Sentinel: scripts/t3e_step2_compute.py reads BASE/data/ucsc/{filename}

## data/macaque/ (13GB)
Source: Qu et al. 2022 (Macaca fascicularis), single-source: GEO GSE196792
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196792
The published human-macaque comparison (Figure 5; 12 cell types matched to Tabula Sapiens)
uses Qu et al. as the sole macaque source; macaque cells are processed from raw 10x UMI
counts, and cell-type assignments come from the author-provided finalcluster metadata in
Qu's Zenodo deposition (10.5281/zenodo.5881495).
Deprecated: an earlier design combined Qu et al. with the Rhesus Immune-Cell Reference
Atlas (RIRA, Macaca mulatta immune subset). That combination used asymmetric cross-atlas
preprocessing (RIRA in deposition-form SCTransform-corrected counts; Qu reprocessed from
raw counts) and was deprecated in favor of the single-source Qu design with symmetric
raw-count preprocessing; the RIRA code and intermediates remain in the repository for
provenance but are not the published result.
Required paths:
  data/macaque/qu_2022/extracted/<GSM>_<tissue>_{barcodes,features,matrix}.{tsv.gz,mtx.gz}  (published source; from the GSE196792 tar)
  data/macaque/biomart_macaque_human_orthologs.csv  (produced by scripts/fetch_macaque_orthologs.py)
  data/macaque/rira/converted/{barcodes,genes,matrix}.{tsv.gz,mtx.gz}  (deprecated combined analysis only)
  data/macaque/rira/rira_metadata.csv  (produced by analysis/macaque/extract_rira_metadata.R; deprecated combined analysis only)
Used by: analysis/macaque/reconstruct_macaque_pipeline.py (canonical loader; single-source Qu is the published path)

## Ortholog tables (BioMart-derived)

Three BioMart-derived ortholog tables underpin the cross-species analyses; all three are archived in the repository so reproduction does not depend on live BioMart queries (Ensembl BioMart returns can shift across releases).

**Human-mouse 1:1 orthologs**: `data/phase1/orthologs_human_mouse.csv`
- Source: Ensembl BioMart release 115, accessed 2026-03-12
- Query: `hsapiens_gene_ensembl` dataset, `mmusculus_homolog_*` cross-reference attributes
- Filter: `orthology_type == "ortholog_one2one"` applied at query time
- 17,187 rows (the post-1:1-filter table cited in manuscript MC:229)
- Downstream atlas-intersection filtering yields the 16,959-gene working space

**Human-macaque orthologs**: `data/macaque/biomart_macaque_human_orthologs.csv`
- Source: Ensembl BioMart release 115, accessed 2026-03-15
- Query: `mfascicularis_gene_ensembl` dataset, `hsapiens_homolog_*` cross-reference attributes (primary key is M. fascicularis Ensembl gene ID)
- 40,305 rows (raw BioMart return, unfiltered)
- Producer script: `scripts/fetch_macaque_orthologs.py` (reproduces the query against current BioMart; archived CSV is canonical)
- Downstream 1:1 filtering and human-name mapping by `scripts/nhp_ortholog_assessment.py` yields the 13,927-gene operating space cited in manuscript MC:273

**Human-mouse lemur 1:1 orthologs**: `analysis/mouse_lemur/biomart_mouse_lemur_human_orthologs.csv`
- Source: Ensembl BioMart, live `pybiomart` query against ensembl.org; release NOT pinned in code or commit; queried on or about commit date 2026-04-05 (Ensembl 116 was expected April 2026, so the served release - 115 or early 116 - is not determinable from records)
- Query: `hsapiens_gene_ensembl` dataset, `mmurinus_homolog_*` cross-reference attributes (primary key is human Ensembl gene ID; mouse lemur = Microcebus murinus, BioMart code `mmurinus`)
- Filter: 1:1 selection via `mmurinus_homolog_orthology_type` at query time, deduplicated by human Ensembl ID
- 16,655 rows (1:1 orthologs; cached columns: human_ensembl_id, human_gene_name, lemur_ensembl_id - schema differs from the other two tables)
- Producer script: `analysis/mouse_lemur/00_feasibility_check.py` (`check_orthologs()`, live BioMart query with CSV cache); consumed by `analysis/mouse_lemur/01_run_pipeline.py`. Archived CSV is canonical.

The human-mouse and human-macaque access dates correspond to manuscript MC:228 and MC:273 respectively; the human-mouse lemur table records no access date (commit-dated 2026-04-05). Re-querying BioMart at a different release may return materially different ortholog assignments; archived CSVs are the canonical source for reproduction.

## data/h3k27ac/ (601MB)
Source: H3K27ac ChIP-seq for 6 immune cell types (T3-E null #9 sensitivity test).
Human (ENCODE narrowPeak, GRCh38): downloaded programmatically by
scripts/t3e_step3b_enhancer.py via the ENCODE REST API for experiments
ENCSR835OJV (CD8+ T), ENCSR391EQV (NK), ENCSR000ASJ (Monocyte),
ENCSR000AUP (B cell), ENCSR120WKZ (CD4+ T), ENCSR267YXV (Neutrophil).
Mouse (Lara-Astiaso 2014 bigWig, mm9): downloaded programmatically from
GEO supplementary files for samples GSM1441277–GSM1441283.
Required path: data/h3k27ac/{human,mouse}/<celltype>_*.{bed.gz,bigWig}
Fetch by running scripts/t3e_step3b_enhancer.py directly. Note that
reproduce/run_all.sh gates that step on a sentinel file
(data/h3k27ac/SENTINEL_FETCHED) that no script in this repository writes, so
the step always skips under run_all.sh; running the script directly is
unaffected.

## data/annotations/ (UCSC refGene)
Source: UCSC Genome Browser refGene tables, used for transcription-start-site
coordinates in the T3-E enhancer analysis (S1 Text mechanistic null 9).
Downloaded programmatically by scripts/t3e_step3b_enhancer.py from
https://hgdownload.soe.ucsc.edu/goldenPath/{hg38,mm10}/database/refGene.txt.gz
Required path: data/annotations/{hg38,mm10}_refGene.txt.gz
Release note: that URL is a rolling UCSC path serving the current table, and
the download is not md5-verified or date-stamped, so no refGene release is
recorded. The reported statistic is re-derivable in kind but not
byte-reproducibly; a re-run may pick up newer TSS coordinates. By contrast the
phastCons download under data/ucsc/ is logged with date and verified md5 in
output/validation/t3e_chromatin/download_log.txt.

## data/string/ (STRING v12.0)
Source: STRING v12.0 human protein-protein interaction network (taxon 9606),
used for the PPI network-centrality analysis (S1 Text mechanistic null 7).
Downloaded programmatically by scripts/19_ppi_centrality.py, which skips the
download if the files are already present:
  https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz
  https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz
Required path: data/string/9606.protein.{links,info}.v12.0.txt.gz
The v12.0 release is pinned in the download URLs and in the on-disk filenames.

## data/dilirank/ (9.9GB)
Required files (place at data/dilirank/<filename>):
  - dilirank_v2.xlsx — DILIrank v2, FDA NCTR (request via
    https://www.fda.gov/science-research/liver-toxicity-knowledge-base-ltkb/dilirank-dataset)
  - lincs_l2_epsilon.gctx — LINCS L1000 Phase I Level 2 GEX, epsilon plate
  - lincs_l2_delta.gctx   — LINCS L1000 Phase I Level 2 GEX, delta plate
    Both .gctx from GEO GSE92742 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742),
    978-landmark-gene Level 2 files.
Used by: scripts/24_dilirank_analysis.py (Paper 1 / Supplementary Text S1 DILI validation)
Note: This dataset is not exercised by reproduce/run_all.sh; analysis is
optional and the pipeline skips gracefully if files are absent.

## CellMarker (retained derived subset: data/validation/cellmarker/cellmarker_human_filtered.csv)
The repository retains only a small derived subset, cellmarker_human_filtered.csv (columns: gene_symbol, cell_type). It contains the experimentally-validated human marker entries (flow cytometry, immunohistochemistry, in-situ hybridization) drawn from CellMarker 2.0; see Materials and Methods for the exact filtering criteria.

This subset is an input to the conserved-contribution gate reported in Results section 5 (Fig 5), where it supplies the marker reference set for the enrichment test. It is additionally read by retained validation tooling for analyses the current paper does not report; SCOPE.md records that boundary.

The raw CellMarker database is not redistributed here; obtain it from CellMarker 2.0 (http://bio-bigdata.hrbmu.edu.cn/CellMarker/, downloaded 2026-03-16) and apply the same filter to regenerate. Reference: Hu et al., CellMarker 2.0, Nucleic Acids Research, 2023.

## ChEMBL / DrugBank drug-target annotations (not retained, not reconstructible)
Mechanistic null 10 in S1 Text — drug-target conservation against per-type divergence, reported as an exploratory and non-significant test — used drug-target and mechanism tables derived from ChEMBL and DrugBank. Those derived tables were removed from the tree and are not redistributed, and no database release, version, accession or download date was recorded anywhere in the repository, so this input cannot be reconstructed from the deposit. The reported statistic is available as output/t3g/primary_correlation_results.json, which retains the derived counts (446 drug-target genes, 420 of them inside the 16,959-gene ortholog set); those counts are not sufficient to identify a release. Re-running the test requires a fresh derivation against a newly recorded ChEMBL and DrugBank release. Literature references: Mendez et al., ChEMBL, Nucleic Acids Research, 2019; Wishart et al., DrugBank 5.0, Nucleic Acids Research, 2018 — these are the database papers, not the releases used.

The DILI analysis under data/dilirank/ separately uses DrugBank-sourced CYP450 substrate annotations obtained from a curated Figshare dataset; that is a different input, it is recorded in the analysis output, and it is unaffected by the gap above.

## Mouse Cell Atlas (MCA) cell assignments
Source: figshare 5435866 (MCA_CellAssignments.csv). Used only by the T1a Mouse-Cell-Atlas
feasibility check (scripts/t1a_mca_feasibility.py), which is not a manuscript figure.
Reconstitute by running scripts/12_t1a_mca_download.py, which downloads the file into
data/replication/.

## MSigDB gene-set catalogs (data/phase3/catalogs/*.gmt)
Source: MSigDB via the R package msigdbr 26.1.0 (MSigDB release 2026.1.Hs). The .gmt catalogs (c2_reactome, c5_gobp, hallmark) are not redistributed here; regenerate them by running phase3_extract_catalogs.R with msigdbr 26.1.0 installed, which produces deterministic output for that pinned version. The fgsea step these catalogs were written for is not implemented in this repository and no GSEA result is deposited; the figure build in reproduce/run_all.sh reads the deposited phase3 outputs and requires no .gmt.

## Macaque atlas QC tables (NHPCA; feasibility scouting; not redistributed)
During early non-human-primate feasibility scouting, sample- and cluster-level QC
summary tables from the Non-Human Primate Cell Atlas (NHPCA) - the adult Macaca
fascicularis whole-body atlas of Han et al. (2022), "Cell transcriptomic atlas of the
non-human primate Macaca fascicularis," Nature 604, 723-731
(doi:10.1038/s41586-022-04587-3) - were examined. The atlas spans ~1.14 million cells
across 45 tissues; raw data are deposited at the CNGB Nucleotide Sequence Archive under
accession CNP0001469 (portal: https://db.cngb.org/nhpca/). These supplementary QC
tables are third-party material and are not redistributed here; obtain them from that
record. They are not used by any analysis in this repository; the macaque comparison
uses single-source Qu et al. 2022 (GEO GSE196792) expression data and the BioMart
macaque-human orthologs under data/macaque/ (the RIRA combination was deprecated).
