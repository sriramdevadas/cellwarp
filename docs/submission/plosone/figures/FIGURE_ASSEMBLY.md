> **Status: historical record, not a description of the built state.** This was the working
> plan for the figure-assembly pass, and it is retained as the record of that pass.
> `reproduce/figure_script_map.md` is the current authority for how each figure and table is
> produced. Two statements below no longer match the producers: the panel tables cite `.pdf`
> sources where `build_main_figures.py` reads the `.png` siblings, and the Fig 1D row says the
> panel is built from `analysis/mouse_lemur/procrustes_results.json` where the producer reads
> `null_distribution.npy` and hard-codes the statistics. The "Remaining to finish figures"
> list is the plan as it stood, not a live checklist. Nothing below has been corrected or
> ticked off: a record edited to look current stops being a record.

# Figure assembly — PLOS ONE (five consolidated main figures)

Maps each of the five new main figures to its source panels, flags panels that
must be built fresh, and records PLOS spec + the two known legibility fixes.
Fig 2C (the only genuinely new figure) is **built** here
(`Fig2C_bg_replication.{tiff,pdf,png}`, `build_fig2c_bg.py`). The other four
consolidate existing panels in `figures/panels/` (old numbering) plus a few
panels promoted from SI or built fresh, as noted.

## PLOS ONE figure spec (verified at journals.plos.org/plosone/s/figures)

- Format **TIFF or EPS** (TIFF preferred); **300–600 dpi**; **RGB (8-bit) or grayscale**.
- Width **6.68–19.05 cm** (789–2250 px @300dpi); height **≤ 22.23 cm** (≤ 2625 px).
- **≤ 10 MB per file**; flattened; **LZW compression**; single page; no alpha; no embedded titles/captions.
- Fonts **Arial, Times, or Symbol only, 8–12 pt**, legible at print size.
- Known legibility fixes to apply (CC-WRITE §4): the LOOCV panel (old Fig 1E) and
  the old Fig 6 axis/tick fonts were too small — rebuild any panel derived from them
  at ≥ 8 pt Arial before promoting/demoting. Avoid Unicode superscripts in Arial
  (missing glyphs); write "p ≤ 0.0001", not "p ≤ 10⁻⁴", inside figures.

## Fig 1 — The configuration of cell types is conserved

| Sub-panel | Source | Action |
|---|---|---|
| 1A pipeline schematic | `panels/fig1a_pipeline_schematic.pdf` | reuse |
| 1B 1M-permutation null | `panels/fig1b_null_1M.pdf` | reuse |
| 1C lineage-stratified null | `panels/fig1c_lineage_stratified.pdf` | reuse |
| 1D mouse-lemur null | **BUILD** from `analysis/mouse_lemur/procrustes_results.json` (obs/null 0.346, p<1e-4, n=15) | new panel; label **~75 Mya** (not 90 — see NUMBER_DIFF C1) |

Demote to SI: `fig1d_bootstrap`, `fig1e_loocv` → S1/S3 Fig (rebuild fig1e at ≥8 pt).

## Fig 2 — Both layers conserved + primate replication

| Sub-panel | Source | Action |
|---|---|---|
| 2A ellipsoid-alignment heatmap | `panels/fig3a_ellipsoid_heatmap.pdf` | reuse |
| 2B pre/post-rotation compression | `panels/fig3b_pre_post.pdf` | reuse |
| 2C **BG three-pair replication** | `Fig2C_bg_replication.tiff` (built here) | **done** |

Demote to SI: `fig3c_layer_nulls`, `fig3d_layer_scatter` → S Fig.

## Fig 3 — Configuration robust across atlases

| Sub-panel | Source | Action |
|---|---|---|
| replication obs/null summary | `panels/fig4d_replication_summary.pdf` (or `fig2d_replication_summary.pdf`) | reuse as the single consolidated panel (obs/null 0.45–0.55 + pan-Census 0.811 + two non-replications) |

Source nulls `fig4a_sun2023`, `fig4b_pansci`, `fig4c_cellhint` → SI (S Fig); the
within-human microwell diagnostic (`fig4e_human_control`) is the S2 Fig referenced in Results §3.

## Fig 4 — Per-type divergence not resolvable (climax figure)

| Sub-panel | Source | Action |
|---|---|---|
| 4A bootstrap CI forest | `analysis/bootstrap_rankings/` (currently S3A via `scripts/composite_figS3.py`) | **PROMOTE** from SI; rebuild at ≥8 pt |
| 4B within/cross-atlas inversion scatter | S3B (ρ = −0.410) | **PROMOTE** from SI |
| 4C simulation recovery ceiling | `analysis/simulation_study/simulation_figures.py` (ρ ≈ 0.42; currently S1D) | **PROMOTE** from SI |

Make legible and self-explanatory (this is the paper's boundary result).

## Fig 5 — Conserved identity genes

| Sub-panel | Source | Action |
|---|---|---|
| conservation-quartile geometry | `analysis/conserved_contribution/make_figure7.py` (old Fig 2A/2D; obs/null 0.384 vs matched-random 0.525) | reuse |
| master-TF enrichment | `analysis/conserved_contribution/make_figure7.py` Fig 2C (0.94 vs 0.54 / 0.76) and/or `panels/fig2b_cellmarker.pdf` | reuse |

## Not used (cut / other project)

`fig7a_treeness`, `fig7b_density` (treeness line, cut); `fig5*/fig6a_macaque*`
(macaque → SI per Results §3 / S1 Text); `fig6b_mechanistic_nulls` → S1 Table/S Fig;
`figs4d_dili`/`figS6_dili` (DILI, banked — SCOPE.md).

## Legend text

Findings-first legends for all five figures live inline in the manuscript
(`manuscript_combined.txt`), each immediately after the paragraph that first cites
it, with the caption title on its own line above the legend body. There is no
separate FIGURE LEGENDS section. The first sentence of each states the finding,
not "Panel A shows…".

## Remaining to finish figures

1. Build Fig 1D (mouse-lemur null) fresh.
2. Promote Fig 4A/B/C from SI, rebuild at ≥8 pt Arial.
3. Compose Fig 1, 2, 3, 5 from the panels above into single flattened TIFFs (≤19.05 cm, LZW).
4. Apply the fig1e / old-fig6 font fix to any promoted/demoted panel.
5. Re-export all five as TIFF (RGB, 300–600 dpi, LZW) — Fig 2C already conforms.
