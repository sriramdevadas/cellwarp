# Supplementary Materials Manifest

**Target journal:** PLOS ONE
**Updated:** 2026-07-23

Authoritative legends and "related to" cross-references for every supplementary
figure and table are in the Supporting Information section of the current
manuscript (`docs/submission/plosone/manuscript_combined.txt`). This file is a
directory manifest only; legends and figure references are not duplicated here,
so nothing in it re-stales when the main figures are renumbered.

## Supplementary items
- 5 supplementary figures: S1-S5 Fig, deposited in `figures/submission/supplementary/`
  (`figS1`-`figS5`).
- 12 supplementary tables: S1-S7 Table and S9-S13 Table (no S8). These 12 items
  span 13 files, because S9 Table is provided as two CSVs
  (`table_S9_genestd_standardization.csv` + `table_S9_schemeB_CPC1_markers.csv`).
  Item count 12, file count 13.
  S13 Table is the inventory of statistical tests and the multiple-comparison
  family (`table_S13_test_inventory.xlsx`). It moved here from
  `docs/submission/figures_for_review/`, the retired packet's output directory,
  where it had been the only item not designated by this manifest. It is the one
  supplementary item with no from-scratch producer; `scripts/table1_formatting.py`
  applies its edits in place and reproduces it byte-for-byte, and
  `tests/test_submission_packet_consistency.py::test_table_1_lock_md5` pins it.

## Files in `docs/supplementary_materials/`
- `table_S*.csv` and `table_S*.xlsx`: per-table data exports (CSV or XLSX per
  table); each table's legend is in the manuscript Supporting Information.
- `figure_S8_markernull.pdf` / `.png`: legacy file, not a current SI item -- a
  byte-identical copy of the deposited
  `figures/submission/supplementary/figS5_markernull.pdf` (S5 Fig) under its old
  "S8" name; see Legacy artifacts.
- `supplementary_legends.md`: retired stub; authoritative legends are in the manuscript.
- `MANIFEST.md`: this file.

No combined Supporting Information PDF is built or deposited; Supporting
Information travels as the separate PLOS-named items listed above, and
`scripts/assemble_supplementary_pdf.py` is retired.

## Legacy artifacts (tracked, not current SI items)
- `table_S8_marker_ortholog_retention.csv`: the S8 Table was cut from the current
  paper.
- `figure_S8_markernull.pdf` / `.png`: duplicate of the deposited `figS5_markernull`
  under its old "S8" name.

Both remain tracked because they are pinned by
`scripts/build_submission_packet.py` and the packet consistency test.
