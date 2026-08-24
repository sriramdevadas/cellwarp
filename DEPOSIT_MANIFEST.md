# DEPOSIT_MANIFEST

Status: CURRENT-STATE manifest. It records what each Zenodo record holds now, and carries its version-1 history in a labelled section rather than in place of it. Both records are published at **version 2**, built from commit 5e77414 (1,094 tracked files) and published 2026-08-24: CODE 10.5281/zenodo.22083132 and DATA 10.5281/zenodo.22083465. Version 1 remains published, built from 4fca942 (1,100 tracked files) on 2026-06-17, and remains the deposit for the earlier PLOS Computational Biology submission. The manuscript cites the CONCEPT DOIs 10.5281/zenodo.20735611 and 10.5281/zenodo.20735639, which resolve to the newest version, which is why publishing v2 required no change to submitted text. Counts below that are labelled v1 are v1's; derive any new version's counts at build time.

## Purpose

Persist the Class-O / Class-D / Class-R deposit scheme so the Zenodo upload is a checklist rather than a re-derivation. This document is the single source of truth for what goes in each Zenodo record, what is deliberately dropped, and how dropped material is regenerated. It reflects the as-built archives: the CODE record is the full tracked repository snapshot; the DATA record is the curated original-output selection.

Two Zenodo records, each now published at version 2. Published records are permanent, so the next deposit is a new version of the same record rather than a fresh one; that preserves the concept DOIs the manuscript cites. At the v2 upload, set the CODE-to-DATA and DATA-to-CODE cross-links using the concept DOIs, and give the DATA record the GitHub related identifier it currently lacks.

## CODE record

- DOI (cite this): **concept 10.5281/zenodo.20735611** - always resolves to the newest version.
- Current version: **v2 10.5281/zenodo.22083132** (commit 5e77414, published 2026-08-24).
- Superseded version: v1 10.5281/zenodo.20735612 (commit 4fca942, published 2026-06-17). Still live and findable; it carries the superseded PLOS Computational Biology paper title. Named here as history - do not quote it as the record's DOI.
- Upload type: Software
- License: MIT
- Upload mode: MANUAL (not a GitHub-linked release). Zenodo only auto-applies .zenodo.json on GitHub-linked releases, so .zenodo.json will NOT auto-populate this manual upload - enter the fields manually, using .zenodo.json as the reference content.
- Cross-link, as actually set on v2: `References -> 10.5281/zenodo.20735639` (the DATA **concept** DOI) and `IsIdenticalTo -> https://github.com/sriramdevadas/cellwarp`. Cross-links point at concept DOIs, not version DOIs, so they keep resolving to the newest version of the other record.

Contents:
- The full git-tracked repository snapshot at the public root: every tracked file (1,100 files). This includes src/, scripts/, analysis/, reproduce/, docs/, tests/, output/, figures/, assets/, and tracked data/ (the three BioMart ortholog tables, the retained CellMarker derived subset cellmarker_human_filtered.csv, the centroid inputs, the replication manifests and download logs), plus all root config and metadata files (CITATION.cff, .zenodo.json, DATA_SOURCES.md, this manifest, etc.).
- The CODE record is a superset of the DATA record; the DATA record re-publishes the original outputs as a standalone CC0 dataset.

Metadata (manual entry; reference .zenodo.json for full content):
- Creators: Sriram Devadas; ORCID 0009-0002-9180-1390; affiliation Independent Researcher, Acton, MA, USA.
- Affiliation is "Independent Researcher" everywhere. This work was conducted independently, without institutional funding, resources or involvement, and that affiliation should be used in every deposit and metadata record.
- Title (FINAL): "CellWarp: analysis code for the cross-species transcriptomic-geometry study"
- Description: the `description` field of `.zenodo.json`. Take the description from there and
  nothing else - in particular **do not take the title from `.zenodo.json`**, whose `title`
  field is `"CellWarp"`, the software name, not a deposit title. That is the instruction this
  line used to carry, and following it literally would title v2 `CellWarp`, with nothing in the
  repository to compare the published title against.
- Why the title changes between versions: v1 carries the superseded PLOS Computational Biology
  paper title, which appears nowhere in this repository - it was typed at the upload form. The
  FINAL title above contains no paper title at all, so it cannot go stale again if the
  manuscript title moves in review, and it reuses the phrasing already used in the DATA record's
  description below. The link to the article is a related identifier set at acceptance, when
  there is an article DOI to point at; it does not belong in the title.

## DATA record

- DOI (cite this): **concept 10.5281/zenodo.20735639** - always resolves to the newest version.
- Current version: **v2 10.5281/zenodo.22083465** (commit 5e77414, published 2026-08-24).
- Superseded version: v1 10.5281/zenodo.20735640 (commit 4fca942, published 2026-06-17). Still live and findable; it carries the superseded PLOS Computational Biology paper title. Named here as history - do not quote it as the record's DOI.
- Upload type: Dataset
- License: CC0 1.0 Universal (public domain dedication)
- Upload mode: MANUAL. No CC0 metadata scaffold exists in the repo, so the full metadata is captured here.
- Cross-link, as actually set on v2: `References -> 10.5281/zenodo.20735611` (the CODE **concept** DOI) and `IsIdenticalTo -> https://github.com/sriramdevadas/cellwarp`. Cross-links point at concept DOIs, not version DOIs, so they keep resolving to the newest version of the other record.

Metadata (manual entry at upload):
- Title (FINAL): "CellWarp: original analysis outputs (Class-O)"
- Upload type: Dataset
- License: CC0 1.0 Universal
- Creators: Sriram Devadas; ORCID 0009-0002-9180-1390; affiliation Independent Researcher, Acton, MA, USA.
- Description (FINAL): Original analysis outputs from the CellWarp cross-species transcriptomic-geometry study. Contains computed outputs only: all tracked files under output/, all figures under figures/, and the data-type result files under analysis/ (.csv/.json/.npy/.png/.tsv/.npz). Raw third-party atlases (Class-R) and derived intermediates (Class-D centroids) are not included and are regenerable - see the dropped-classes regeneration map below and DATA_SOURCES.md. The code that produced these outputs, plus the input reference tables, is deposited separately under DOI 10.5281/zenodo.20735611 (MIT).

Contents (Class-O original analysis outputs), as built (761 files):
- output/ : all 496 tracked files. The 5 empty .gitkeep placeholders are dropped from the deposit; ~424 further files under output/ are gitignored scratch/large artifacts and are intentionally excluded.
- figures/ : all 165 tracked files (main and supplementary figure PDFs and PNGs).
- analysis/ : the 105 data-type result files (.csv/.json/.npy/.png/.tsv/.npz), EXCLUDING analysis/mouse_lemur/biomart_mouse_lemur_human_orthologs.csv (a BioMart reference table, which travels with the CODE record). analysis/ code (.py/.R/.md) is not in the DATA record.
- EXCLUDED from DATA: everything under data/ (input reference tables and the 2 Class-D centroid files in data/centroids/); the .gitkeep placeholders wherever they occur; and, under analysis/ only, the code files (.py/.R/.md) that the data-type filter in the previous bullet already drops. This last exclusion is scoped to analysis/ and does NOT reach output/. An earlier version of this line claimed .md was dropped everywhere, which contradicted both the bullet above it and step 1 of the snapshot procedure, and is settled against the published v1 archive rather than by preference: cellwarp-data-4fca942.zip contains 20 output/*.md files, its output/ entry count is 491 = 496 tracked minus 5 .gitkeep, and 491 + 165 + 105 = 761, the archive's own total. 'All output/ minus the .gitkeep placeholders' is therefore the operative rule.

## Dropped classes + regeneration map

Class-R (raw third-party) - NOT deposited:
- Gitignored .h5ad atlas files under data/ and analysis/ (including the 100G+ and 37G atlases).
- Regenerate by re-download from the source accessions below.

Class-D (derived intermediates) - NOT deposited in the DATA record:
- The 2 tracked centroid files under data/centroids/ (pansci_16type_centroids.csv, sun2023_15type_centroids.csv) ship in the CODE record only; additional centroid intermediates are gitignored.
- Regenerate by re-running the released pipeline against the re-downloaded raw inputs.

Removed third-party material (referenced by accession, not redistributed):
- d53 removals: ChEMBL target/mechanism tables, tms_facs metadata, raw CellMarker human/mouse tables, Mouse Cell Atlas (MCA) cell assignments, and MSigDB .gmt catalogs - removed from the tree. Where a release or accession was recorded, DATA_SOURCES.md carries it: CellMarker, MCA and MSigDB each have one. ChEMBL and DrugBank do not: DATA_SOURCES.md records that no release, version, accession or download date was captured, so that input is not reconstructible from the deposit. (The small derived subset data/validation/cellmarker/cellmarker_human_filtered.csv is retained in the CODE record; see DATA_SOURCES.md.)
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
- docs/submission/manuscript_combined.txt                (in the v1 CODE archive; since
  retired from the repository. Retrievable via `git show`. Its Tabula Microcebus deposit
  anchors are carried by DATA_SOURCES.md, which is the primary source above.)

## Snapshot procedure -- executed for v2 at submission, and the template for v3

**This procedure has been run.** v2 of both records was built and published on 2026-08-24 from commit 5e77414 and the steps below are the record of what was done, in past tense where they describe it and generic where they are the template for a v3.

As built, both versions side by side:

```
                     v1                              v2
  commit             4fca942                         5e77414
  published          2026-06-17                      2026-08-24
  CODE entries       1,100  (= git ls-tree -r)       1,094  (= git ls-files)
  DATA entries         761                             717
    output/            491  (496 tracked - 5)          435  (440 tracked - 5)
    figures/           165                             162
    analysis/          105                             120
  CODE version DOI   10.5281/zenodo.20735612         10.5281/zenodo.22083132
  DATA version DOI   10.5281/zenodo.20735640         10.5281/zenodo.22083465
```

Both v2 counts were derived at build time, not copied from this file, and the archives were named for the release commit: `cellwarp-code-5e77414.zip` (md5 `2db39f94...`, 1,094 entries) and `cellwarp-data-5e77414.zip` (md5 `6a1c646b...`, 717 entries), neither carrying a wrapper directory. **Derive any v3 counts the same way rather than reading either column above.**

Steps:
1. CODE archive: a zip of exactly `git ls-files` at the release commit. DATA archive: the selection above - all output/ minus the .gitkeep placeholders, all figures/, and the analysis/ data-type files minus the ortholog reference table. Record both counts at build time; the table above is the v1 and v2 record, not a target to copy.
   - `docs/submission/plosone/coverletter.txt` stays out of the CODE archive: a cover letter is
     correspondence with the editor, not part of the public archived record. Building the archive
     from `git ls-files` drops it, and it is also ignored - `.gitignore:196` is `cover*letter*.txt`,
     which does match that filename (`git check-ignore -v` names that rule), so the file never
     appears in `git status` and a bare `git add -A` will not stage it. Do not `git add -f` it.
2. Upload each set as a NEW VERSION of its existing record, identified by its **concept** DOI (CODE 20735611, DATA 20735639) rather than by any version DOI, so the concept DOIs resolve to the new version. Naming a version DOI here would name whichever version was current when the line was written, which is how this document previously came to point at v1.
3. Enter metadata manually (CODE: from .zenodo.json; DATA: from this manifest) and set the CODE-to-DATA / DATA-to-CODE related-identifier cross-links.
4. Publish both new versions **at submission, not at acceptance** -- decided, and done on
   2026-08-24. This reversed the earlier
   instruction, and the reason is that v1 does not describe this paper. Measured at the
   submission tree: 113 of Gate 1's 232 checks read artifacts that do not exist at 4fca942
   (15 distinct files); v1's own `reproduce/validate.py` defines 30 checks against this paper's
   232; and `docs/submission/plosone/` is 39 tracked files here and 0 at 4fca942. A reviewer
   following the concept DOI during review would land on an archive that cannot substantiate
   the submitted text. Version 1 of each record stays published and remains the deposit for the
   earlier PLOS Computational Biology submission.
