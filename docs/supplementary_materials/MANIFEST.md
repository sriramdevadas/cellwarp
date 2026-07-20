# Supplementary Materials Manifest

**Target journal:** PLOS ONE
**Updated:** 2026-06-11

Authoritative legends and "related to Figure N" cross-references for every supplementary
figure (S1-S8) and table (S1-S12) are in the manuscript Supporting Information section
(`docs/submission/manuscript_combined.txt`). This file is a directory manifest only;
legends and figure references are not duplicated here, so nothing in it re-stales when the
main figures are renumbered.

## Supplementary items
- 5 supplementary figures (S1-S5)
- 12 supplementary tables (S1-S7, S9-S12; no S8)

## Files in `docs/supplementary_materials/`
- `table_S*.csv` and `table_S*.xlsx`: per-table data exports (CSV or XLSX per table); each table's legend is in the manuscript Supporting Information.
- `figure_S8_markernull.pdf` / `.png`: standalone Figure S8 rendering (marker-similarity null).
- `supplementary_legends.md`: retired stub; authoritative legends are in the manuscript.
- `MANIFEST.md`: this file.

## Rebuilding the combined PDF
`scripts/assemble_supplementary_pdf.py` regenerates `supplementary_materials.pdf` from the
manuscript Supporting Information legends, the pre-rendered `figS1`-`figS8` PDFs in
`figures/submission/supplementary/`, and vendored fonts. Figures S1-S8 are pre-rendered
composites, not built at assemble time.
