# DEPOSIT_MANIFEST

Status: CURRENT-STATE manifest. It records what each Zenodo record holds now, and carries its earlier history in labelled sections rather than in place of it.

**The two records are no longer at the same version, and nothing should be read as though they were.**

- **CODE is at version 3**, built from commit `98e67a7` (1,097 tracked files) and published 2026-08-26: **10.5281/zenodo.22104908**. Superseded but permanently resolvable: v2 `10.5281/zenodo.22083132` (5e77414, 1,094 files, 2026-08-24) and v1 `10.5281/zenodo.20735612` (4fca942, 1,100 files, 2026-06-17). Every prior citation of either keeps working.
- **DATA is at version 2** and was deliberately not re-deposited: `10.5281/zenodo.22083465` (5e77414, 717 files, 2026-08-24). v1 `10.5281/zenodo.20735640` remains published and remains the deposit for the earlier PLOS Computational Biology submission.

The manuscript cites the CONCEPT DOIs 10.5281/zenodo.20735611 and 10.5281/zenodo.20735639, which resolve to the newest version of each record, which is why publishing a new version requires no change to submitted text. Counts below that are labelled with a version are that version's; derive any new version's counts at build time.

## Purpose

Persist the Class-O / Class-D / Class-R deposit scheme so the Zenodo upload is a checklist rather than a re-derivation. This document is the single source of truth for what goes in each Zenodo record, what is deliberately dropped, and how dropped material is regenerated. It reflects the as-built archives: the CODE record is the full tracked repository snapshot; the DATA record is the curated original-output selection.

Two Zenodo records: CODE at version 3, DATA at version 2. Published records are permanent, so each deposit is a new version of the same record rather than a fresh one; that preserves the concept DOIs the manuscript cites. The CODE-to-DATA and DATA-to-CODE cross-links are set with the concept DOIs and were carried through the v3 upload.

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

## Snapshot procedure -- executed for v2 and for CODE v3, and the template for the next one

**This procedure has been run twice.** v2 of both records was built and published on 2026-08-24 from commit 5e77414; **CODE v3 alone** was built and published on 2026-08-26 from commit `98e67a7`. The steps below are the record of what was done, in past tense where they describe it and generic where they are the template for a next version.

### Why CODE went to v3 and DATA did not

**The records are asymmetric by measurement, not by oversight.** Between `5e77414` and `98e67a7`
exactly one file under `output/`, `figures/` or `analysis/` changed:
`analysis/biological_predictors/biological_predictors.py`, a one-line change to a `print` statement.
It is a `.py` under `analysis/`, which the DATA selection rule excludes.

Re-running that rule at both commits:

```
at 5e77414 (v2 as built)   717 files      <- reproduces the published v2 count exactly
at 98e67a7 (v3 CODE)       717 files
file-list differences        0
files whose BYTES changed    0
```

The 717 at `5e77414` is the check that makes the second row trustworthy: the rule was validated
against a published archive before being used to answer the question at the newer commit. A DATA v3
would have been byte-identical to DATA v2, so it was not cut. **A reader must not infer that the two
records move together, because from 2026-08-26 they do not.**

### What v3 changed, and why the archive looks different

CODE v3 exists because every reader-path repair between dispatches 76 and 84 landed after `5e77414`,
so the published archive told a reviewer to run a command that fails on a stock image. Fourteen files
differ from v2 and **no pipeline script is among them** except the `print` statement above; the
manuscript differs by one line, with `reproduce/MANUSCRIPT_MD5` moved to match, so Gate 4 stays
internally consistent inside the archive.

Two deliberate divergences from the v1 and v2 archives, recorded here so neither is later read as a
packaging error:

- **v3 unpacks into a single directory**, `cellwarp-code-98e67a7/`, built with
  `git archive --prefix=`. v1 and v2 carry no wrapper and scatter their contents into whatever
  directory they are opened in. Nothing keys on the layout: `reproduce/validate.py:19` and
  `scripts/build_submission_packet.py:58` both resolve from `Path(__file__).resolve().parent.parent`
  rather than from the working directory.
- **v3 is smaller than v2 while carrying three more files** -- 101,528,695 bytes against
  101,906,462. That is `-9` compression, not truncation. **The entry count is the identity check;
  the byte count is not.**

Route C -- a reviewer arriving from the concept DOI and running the gates -- was measured against the
v3 archive in a clean container before it was uploaded: 232/232, 195 passed, 30/30 pairs, 3 of 3.

As built, both versions side by side:

```
                     v1                        v2                        v3  (CODE only)
  commit             4fca942                   5e77414                   98e67a7
  published          2026-06-17                2026-08-24                2026-08-26
  CODE entries       1,100 (git ls-tree -r)    1,094 (git ls-files)      1,097 (git ls-files)
  CODE version DOI   10.5281/zenodo.20735612   10.5281/zenodo.22083132   10.5281/zenodo.22104908
  CODE wrapper dir   none                      none                      cellwarp-code-98e67a7/
  CODE zip bytes     not recorded              101,906,462               101,528,695
  CODE zip md5       not recorded              2db39f94...               98efed5492ace31e776670f74f1a66bd
  DATA entries         761                       717                     not re-deposited
    output/            491 (496 tracked - 5)     435 (440 tracked - 5)     "
    figures/           165                       162                       "
    analysis/          105                       120                       "
  DATA version DOI   10.5281/zenodo.20735640   10.5281/zenodo.22083465   stays at v2
```

Every count in that table was derived at build time, not copied from this file, and each archive was named for the release commit: `cellwarp-code-5e77414.zip` (md5 `2db39f94...`, 1,094 entries), `cellwarp-data-5e77414.zip` (md5 `6a1c646b...`, 717 entries) and `cellwarp-code-98e67a7.zip` (md5 `98efed5492ace31e776670f74f1a66bd`, sha256 `d9e23f41...`, 1,097 entries). v1 and v2 carry no wrapper directory; v3 does. **Derive any new counts the same way rather than reading any column above.**

Steps:
1. CODE archive: a zip of exactly `git ls-files` at the release commit. **v3 was built with**
   `git archive --format=zip -9 --prefix=cellwarp-code-<commit>/ -o cellwarp-code-<commit>.zip HEAD`.
   That is not an approximation of `git ls-files`: there is no `.gitattributes`, so `git archive` and
   `git ls-files` were verified to yield the same 1,097 paths with zero differences before the build,
   and the zip's path list with the prefix stripped `diff`s clean against `git ls-files`. DATA archive: the selection above - all output/ minus the .gitkeep placeholders, all figures/, and the analysis/ data-type files minus the ortholog reference table. Record both counts at build time; the table above is the v1 and v2 record, not a target to copy.
   - `docs/submission/plosone/coverletter.txt` stays out of the CODE archive: a cover letter is
     correspondence with the editor, not part of the public archived record. Building the archive
     from `git ls-files` drops it, and it is also ignored - `.gitignore:196` is `cover*letter*.txt`,
     which does match that filename (`git check-ignore -v` names that rule), so the file never
     appears in `git status` and a bare `git add -A` will not stage it. Do not `git add -f` it.
2. Upload each set as a NEW VERSION of its existing record, identified by its **concept** DOI (CODE 20735611, DATA 20735639) rather than by any version DOI, so the concept DOIs resolve to the new version. Naming a version DOI here would name whichever version was current when the line was written, which is how this document previously came to point at v1.
3. Enter metadata manually (CODE: from .zenodo.json; DATA: from this manifest) and set the CODE-to-DATA / DATA-to-CODE related-identifier cross-links. **`.zenodo.json` is a partial source and following this step literally loses information twice.** Take the `description` from it and nothing else:
   - its `title` is `"CellWarp"`, the software name, not a deposit title -- the deposit title is the FINAL one recorded above;
   - it carries **only one** of the two related identifiers (`isSupplementTo` -> the GitHub URL) and **not** `references` -> `10.5281/zenodo.20735639`, so entering metadata from that file alone silently drops the CODE-to-DATA cross-link.
   At the v3 upload both survived, verified at DataCite after publication: `IsSupplementTo`, `References`, three `HasVersion` entries, no `IsIdenticalTo`, `version` `98e67a7`.
4. Publish both new versions **at submission, not at acceptance** -- decided, and done on
   2026-08-24. This reversed the earlier
   instruction, and the reason is that v1 does not describe this paper. Measured at the
   submission tree: 113 of Gate 1's 232 checks read artifacts that do not exist at 4fca942
   (15 distinct files); v1's own `reproduce/validate.py` defines 30 checks against this paper's
   232; and `docs/submission/plosone/` is 39 tracked files here and 0 at 4fca942. A reviewer
   following the concept DOI during review would land on an archive that cannot substantiate
   the submitted text. Version 1 of each record stays published and remains the deposit for the
   earlier PLOS Computational Biology submission.
