# DEPOSIT_MANIFEST

Status: FINAL - as-built deposit manifest. Counts and scope reflect the public repository root that the two Zenodo archives are built from. Both records are reserved drafts (DOIs minted, not yet published); publish at acceptance.

## Purpose

Persist the Class-O / Class-D / Class-R deposit scheme so the Zenodo upload is a checklist rather than a re-derivation. This document is the single source of truth for what goes in each Zenodo record, what is deliberately dropped, and how dropped material is regenerated. It reflects the as-built archives: the CODE record is the full tracked repository snapshot; the DATA record is the curated original-output selection.

Two Zenodo records, both reserved as drafts (DOIs minted, not published). Do NOT delete either draft - the reserved DOI is lost if the draft is deleted. Do NOT publish until acceptance. At upload, set the CODE-to-DATA and DATA-to-CODE cross-links (each record's "related identifiers" pointing at the other's DOI).

## CODE record

- DOI: 10.5281/zenodo.20735612
- Upload type: Software
- License: MIT
- Upload mode: MANUAL (not a GitHub-linked release). Zenodo only auto-applies .zenodo.json on GitHub-linked releases, so .zenodo.json will NOT auto-populate this manual upload - enter the fields manually, using .zenodo.json as the reference content.
- Cross-link: relates to the DATA record (DOI 10.5281/zenodo.20735640) and the GitHub repository.

Contents:
- The full git-tracked repository snapshot at the public root: every tracked file (1,100 files). This includes src/, scripts/, analysis/, reproduce/, docs/, tests/, output/, figures/, assets/, and tracked data/ (the three BioMart ortholog tables, the retained CellMarker derived subset cellmarker_human_filtered.csv, the centroid inputs, the replication manifests and download logs), plus all root config and metadata files (CITATION.cff, .zenodo.json, DATA_SOURCES.md, this manifest, etc.).
- The CODE record is a superset of the DATA record; the DATA record re-publishes the original outputs as a standalone CC0 dataset.

Metadata (manual entry; reference .zenodo.json for full content):
- Creators: Sriram Devadas; ORCID 0009-0002-9180-1390; affiliation Independent Researcher.
- Employer is never named anywhere - Independent Researcher only.
- Title / description: as in .zenodo.json (confirm .zenodo.json is current at freeze).

## DATA record

- DOI: 10.5281/zenodo.20735640
- Upload type: Dataset
- License: CC0 1.0 Universal (public domain dedication)
- Upload mode: MANUAL. No CC0 metadata scaffold exists in the repo, so the full metadata is captured here.
- Cross-link: relates to the CODE record (DOI 10.5281/zenodo.20735612) and the GitHub repository.

Metadata (manual entry at upload):
- Title (FINAL): "CellWarp: original analysis outputs (Class-O)"
- Upload type: Dataset
- License: CC0 1.0 Universal
- Creators: Sriram Devadas; ORCID 0009-0002-9180-1390; affiliation Independent Researcher.
- Description (FINAL): Original analysis outputs from the CellWarp cross-species transcriptomic-geometry study. Contains computed outputs only: all tracked files under output/, all figures under figures/, and the data-type result files under analysis/ (.csv/.json/.npy/.png/.tsv/.npz). Raw third-party atlases (Class-R) and derived intermediates (Class-D centroids) are not included and are regenerable - see the dropped-classes regeneration map below and DATA_SOURCES.md. The code that produced these outputs, plus the input reference tables, is deposited separately under DOI 10.5281/zenodo.20735612 (MIT).

Contents (Class-O original analysis outputs), as built (761 files):
- output/ : all 496 tracked files. The 5 empty .gitkeep placeholders are dropped from the deposit; ~424 further files under output/ are gitignored scratch/large artifacts and are intentionally excluded.
- figures/ : all 165 tracked files (main and supplementary figure PDFs and PNGs).
- analysis/ : the 105 data-type result files (.csv/.json/.npy/.png/.tsv/.npz), EXCLUDING analysis/mouse_lemur/biomart_mouse_lemur_human_orthologs.csv (a BioMart reference table, which travels with the CODE record). analysis/ code (.py/.R/.md) is not in the DATA record.
- EXCLUDED from DATA: everything under data/ (input reference tables and the 2 Class-D centroid files in data/centroids/), and all .py/.R/.md/.gitkeep files.

## Dropped classes + regeneration map

Class-R (raw third-party) - NOT deposited:
- Gitignored .h5ad atlas files under data/ and analysis/ (including the 100G+ and 37G atlases).
- Regenerate by re-download from the source accessions below.

Class-D (derived intermediates) - NOT deposited in the DATA record:
- The 2 tracked centroid files under data/centroids/ (pansci_16type_centroids.csv, sun2023_15type_centroids.csv) ship in the CODE record only; additional centroid intermediates are gitignored.
- Regenerate by re-running the released pipeline against the re-downloaded raw inputs.

Removed third-party material (referenced by accession, not redistributed):
- d53 removals: ChEMBL target/mechanism tables, tms_facs metadata, raw CellMarker human/mouse tables, Mouse Cell Atlas (MCA) cell assignments, and MSigDB .gmt catalogs - removed from the tree; see DATA_SOURCES.md for accessions. (The small derived subset data/validation/cellmarker/cellmarker_human_filtered.csv is retained in the CODE record; see DATA_SOURCES.md.)
- d58 removals: macaque whole-body atlas QC tables (NHPCA, Han et al. 2022, CNGB CNP0001469) - examined during NHP feasibility scouting, never used by any analysis; removed from the tree and referenced by accession in DATA_SOURCES.md. Also removed: data/phase1/orthologs_human_zebrafish.csv, an unused BioMart orphan (regenerable from Ensembl BioMart).

Per-source regeneration (accessions from DATA_SOURCES.md):
- Sun et al.:        OMIX002605 / CRA007207
- PanSci:            GSE247719
- Macaque (analysis): GSE196792 (Qu et al. 2022) plus RIRA
- Macaque (scouting QC tables, not used): CNP0001469 (NHPCA, Han et al. 2022)
- Mouse lemur:       DOI 10.1038/s41586-025-09113-9 (Tabula Microcebus)
- DILI:              GSE92742
- Ortholog tables: Ensembl BioMart. Human-mouse and human-macaque are pinned to release 115 (accessed 2026-03-12 / 2026-03-15; see DATA_SOURCES.md). Human-mouse lemur is a live pybiomart query (generator analysis/mouse_lemur/00_feasibility_check.py), release NOT pinned, commit-dated 2026-04-05; archived CSV is canonical (see DATA_SOURCES.md note).

## Authoritative regeneration sources

Primary (canonical): DATA_SOURCES.md (repo root).

Corroborating:
- analysis/census_replication/census_datasets_full.csv   (1,845 rows)
- data/replication/pan_census_manifest.csv               (15 rows)
- data/replication/tabula_microcebus_metadata.csv        (1 row)
- docs/submission/manuscript_combined.txt

## Snapshot procedure (at acceptance)

Source: the single parentless root commit on origin/main (the public repository root; third-party-free history). The working tree is clean at that root; the two archives below are built directly from `git ls-files` at that commit.

Steps:
1. CODE archive: a zip of exactly `git ls-files` at the root (1,100 files). DATA archive: the selection above (761 files: all output/ minus the 5 .gitkeep, all figures/, and the 105 analysis/ data-type files minus the ortholog reference table).
2. Upload each set to its reserved Zenodo draft (CODE -> 20735612, DATA -> 20735640).
3. Enter metadata manually (CODE: from .zenodo.json; DATA: from this manifest) and set the CODE-to-DATA / DATA-to-CODE related-identifier cross-links.
4. Publish both records at acceptance. Until then: keep drafts reserved, unpublished, undeleted.
