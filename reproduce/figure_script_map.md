# Figure and Table to Script Mapping

Every display item in the current CellWarp paper (PLOS ONE submission: 5 main
figures, S1-S5 Fig, and 11 supplementary tables S1-S7 and S9-S12), mapped to the
script that produces the deposited file.

The manuscript is `docs/submission/plosone/manuscript_combined.txt`; its figure
and table legends are the authority for panel content. Old-to-new figure/table
renumbering history (this paper descends from a 7-figure PCOMPBIOL draft) is
recorded separately in `docs/submission/plosone/NUMBER_DIFF.md` and is out of
scope here.

## Main figures

The five main composites are assembled by
`docs/submission/plosone/figures/build_main_figures.py` (PDF + PNG, 300 dpi) into
`docs/submission/plosone/figures/`, which holds the deposited figures; the TIFFs
uploaded to PLOS are rasterised from those PDFs (see Submission TIFFs below). The
manuscript is the authority for their legends. Panels marked "embedded" are
pre-rendered PNGs pulled in by the assembler; panels marked "built in
build_main_figures.py" are drawn fresh from data at assembly time.

`docs/submission/plosone/build_review_pdf.py` concatenates the manuscript text and
these five figure PDFs into an internal-review PDF. That output is deliberately
not tracked, so any copy on disk is stale unless it has just been rebuilt.

### Fig 1 (configuration conserved) -> `Fig1_configuration_conserved.{pdf,png}`

| Panel | Content | Producing script | Depends on |
|---|---|---|---|
| 1A | Pipeline schematic | `scripts/generate_phase2_figures.py` (embedded `figures/panels/fig1a_pipeline_schematic.png`) | none |
| 1B | 1M-permutation null (obs/null 0.52) | `scripts/generate_phase1_figures.py` (embedded `figures/panels/fig1b_null_1M.png`) | `scripts/permutation_1M.py` -> `analysis/permutation_1M/results_1M.json` + `null_distribution_1M.npy` |
| 1C | Lineage-stratified null (0.67) | `scripts/generate_phase1_figures.py` (embedded `figures/panels/fig1c_lineage_stratified.png`) | `scripts/test_lineage_stratified_permutation.py` |
| 1D | Human-mouse-lemur null (0.35, n=15) | built in `build_main_figures.py` (fresh; no panel file) | `analysis/mouse_lemur/null_distribution.npy` (see Known gaps: obs value hard-coded) |

### Fig 2 (two layers + primate replication) -> `Fig2_two_layers_bg.{pdf,png}`

| Panel | Content | Producing script | Depends on |
|---|---|---|---|
| 2A | Per-cell-type Krzanowski S heatmap (k=1,3,5) | `scripts/generate_phase3_figures.py` (embedded `figures/panels/fig3a_ellipsoid_heatmap.png`) | `scripts/t3b_ellipsoid_alignment.py` -> `output/mechanistic/ellipsoid_alignment/` |
| 2B | Aggregate S before/after centroid-optimal rotation | `scripts/generate_phase3_figures.py` (embedded `figures/panels/fig3b_pre_post.png`) | `scripts/t3b_ellipsoid_alignment.py` |
| 2C | Basal-ganglia three-pair Layer-2 replication | `docs/submission/plosone/figures/build_fig2c_bg.py` (embedded `Fig2C_bg_replication.png`) | vendored `docs/submission/plosone/figures/bg_results/layer2_results_{pair}.json` + `layer2_cpc1_drivers_{pair}_W2_schemeB.csv` |

`build_fig2c_bg.py` also writes the standalone `Fig2C_bg_replication.{pdf,tiff}`;
the `.tiff` is the PLOS-spec deposit (RGB, LZW, 300 dpi). `reproduce/validate.py`
carries twelve basal-ganglia self-consistency checks over the same `bg_results/`
JSONs.

### Fig 3 (configuration robust) -> `Fig3_configuration_robust.{pdf,png}`

Single panel.

| Panel | Content | Producing script | Depends on |
|---|---|---|---|
| (whole) | Replication obs/null bar chart | `scripts/generate_phase1_figures.py` (embedded `figures/panels/fig4d_replication_summary.png`) | all replication scripts: `scripts/16_sun2023_replication.py`, `scripts/pansci_replication.py`, `scripts/33_cellhint_replication.py`, `analysis/census_replication/`, plus the Andrews and MCA x HCA non-replications |

### Fig 4 (per-type not resolvable) -> `Fig4_pertype_not_resolvable.{pdf,png}`

Built fresh; no embedded panels.

| Panel | Content | Producing script | Depends on |
|---|---|---|---|
| 4A | Within-atlas precision (bootstrap CI forest) | built in `build_main_figures.py` | `analysis/bootstrap_rankings/bootstrap_summary.csv` <- `analysis/bootstrap_rankings/bootstrap_ranking_analysis.py` |
| 4B | Within-atlas precision vs cross-atlas rank shift (rho = -0.41, n=20) | built in `build_main_figures.py` | `analysis/cross_reference/master_ranking_table.csv` <- `analysis/cross_reference/cross_reference_analysis.py` |
| 4C | Simulation recovery ceiling (rho ~ 0.42) | built in `build_main_figures.py` | `analysis/simulation_study/simulation_results.json` <- `analysis/simulation_study/simulation_study.py` |

### Fig 5 (conserved identity genes) -> `Fig5_conserved_identity_genes.{pdf,png}`

A four-panel figure produced whole by `make_figure7.py`; `build_main_figures.py`
wraps its PNG.

| Panel | Content | Producing script | Depends on |
|---|---|---|---|
| 5A-5D | C distribution / Hartigan dip (A); C vs expression and specificity (B); master-TF enrichment vs matched backgrounds (C); donor-split reproducibility (D) | `analysis/conserved_contribution/make_figure7.py` -> `figures/main/fig7_conserved_contribution.{pdf,png}`; `build_main_figures.py` embeds the `.png` | `analysis/conserved_contribution/run_gate.py` (`gate_results.json`), `run_robustness.py` (`robustness_results.json`), `donor_stability/donor_stability_results.json`, `gene_conservation_core.csv`, `donor_stability/agg_*.npz` (see Known gaps) |

## Submission TIFFs

`docs/submission/plosone/figures/build_submission_tiffs.py` rasterises the five
main-figure PDFs above to `Fig1.tif` through `Fig5.tif` in the same directory, at
300 dpi, RGB, LZW, alpha composited on white. It reads only the `.pdf` and its
sibling `.png` (the PNG is an independently rasterised copy, used to check the
flatten) and regenerates nothing: an explicit table of expected pixel dimensions
aborts the build if a figure has changed size. The `.tif` files are the artifacts
uploaded to PLOS, whose file names must match the in-text citations; the `.pdf`
and `.png` remain the deposited figures.

## Supplementary figures (S1-S5 Fig)

Deposited in `figures/submission/supplementary/`.

| S Fig | Content | Producing script | Depends on |
|---|---|---|---|
| S1 Fig | Pipeline validation: independent PCA + simulation study | `scripts/build_submission_figures.py` (composites `figures/supplementary/figS1_independent_pca.pdf` + `figS7_simulation_study_polished.pdf`) | `analysis/independent_pca_sensitivity/run_independent_pca.py`; `analysis/simulation_study/simulation_figures.py` |
| S2 Fig | Parameter, protocol, and negative-control sensitivity | `scripts/build_submission_figures.py` (panel labels) + `scripts/56_add_figs2_panel_f.py` (panel F) | `scripts/17_pca_sensitivity.py`, `scripts/18_pca_sensitivity_v2.py`, `scripts/14_smartseq2_sensitivity.py`, `analysis/expanded_negative_controls/expanded_negative_controls.py`, replication outputs |
| S3 Fig | Bootstrap stability of per-type divergence rankings | `scripts/composite_figS3.py` (invoked by `build_submission_figures.py`) | `analysis/bootstrap_rankings/bootstrap_ranking_analysis.py` |
| S4 Fig | Matched-scale 6-type human-vs-human negative control | `scripts/49_build_figS7_matched_scale.py` (producer filename retains the old `figS7` stem) | `scripts/test_35type_human_control.py` -> `output/phase2/negative_control_v2/` |
| S5 Fig | Marker-similarity-stratified null | `analysis/sensitivity_analyses/markernull.py` (producer writes `figure_S8_markernull.{pdf,png}` to `docs/supplementary_materials/`, not the deposited `figS5_markernull.pdf`; see Known gaps) | primary centroids; species-averaged gene-space |

## Supplementary tables (S1-S7, S9-S12)

The current paper has no S8 Table (the former ortholog-retention table was cut).
Canonical sources live in `docs/supplementary_materials/`. For several tables an
analysis script computes the values and the deposited CSV/XLSX is then written or
edited in place by the materializer
`scripts/46_synthesis_pass_supplementary_table_edits.py`; both are listed.

| Table | Content | Content producer | Canonical writer -> file |
|---|---|---|---|
| S1 Table | Biological predictors, cross-atlas consistency, three-species summary | `analysis/biological_predictors/biological_predictors.py` + `analysis/ranking_replication/ranking_replication_analysis.py` (+ `scripts/create_table_S1.py`) | `scripts/46_synthesis_pass_supplementary_table_edits.py` (`edit_table_s1`) -> `table_S1.xlsx` |
| S2 Table | Simulation parameters + bootstrap ranking CIs | `analysis/simulation_study/simulation_study.py` + `analysis/bootstrap_rankings/bootstrap_ranking_analysis.py` (+ `scripts/create_table_S2.py`) | `scripts/46_synthesis_pass_supplementary_table_edits.py` (`edit_table_s2`) -> `table_S2.xlsx` |
| S3 Table | CellHint rank-reversal analysis | `analysis/cellhint_investigation/investigate_rank_reversal.py` | `scripts/46_synthesis_pass_supplementary_table_edits.py` (`edit_table_s3`) -> `table_S3.csv` |
| S4 Table | Progressive CellHint harmonization | `analysis/harmonized_replication/harmonized_replication.py` | `scripts/46_synthesis_pass_supplementary_table_edits.py` (`edit_table_s4`) -> `table_S4.csv` |
| S5 Table | 35-type cell-type matching | `scripts/08_cell_type_inventory.py` -> `output/phase2/cell_type_inventory*.csv` | `scripts/46_synthesis_pass_supplementary_table_edits.py` (`edit_table_s5`) -> `table_S5.csv` |
| S6 Table | CPC1 driver genes | `scripts/generate_table_S6.py` | self -> `Table_S6_CPC1_driver_genes.xlsx` |
| S7 Table | Layer-1 coherence under ribosomal/housekeeping exclusion | `analysis/sensitivity_analyses/layer1_exclusion.py` (writes `layer1_exclusion_ranking_*.csv` in its analysis dir) | NO SCRIPTED WRITER of the canonical `table_S7_layer1_housekeeping_exclusion.csv` (see Known gaps) |
| S9 Table | Per-gene standardization (2 CSVs) | `analysis/sensitivity_analyses/genestd_standardization.py` | self -> `table_S9_genestd_standardization.csv` + `table_S9_schemeB_CPC1_markers.csv` |
| S10 Table | Marker-similarity-stratified permutation null | `analysis/sensitivity_analyses/markernull.py` | self -> `table_S10_markernull.csv` |
| S11 Table | Per-gene cross-species conservation score C | `analysis/conserved_contribution/make_table_s11.py` (dep `run_gate.py`) | self -> `table_S11_gene_conservation.csv` |
| S12 Table | Software environment and version-pinned dependencies | NO IN-REPO PRODUCER; hand-authored, source `requirements.txt` (see Known gaps) | hand-authored -> `table_S12_software_environment.csv` |

## Supporting-information controls with no display item

S1 Text §10 reports a control that produces neither a figure nor a table, so it has no row
above; its numbers are quoted in the text alone.

| Item | Content | Producing script | Depends on |
|---|---|---|---|
| S1 Text §10 | Selection/derangement circularity control: real conserved obs/null 0.384 against a derangement sigma-null 0.991 +- 0.021 (z = -29.5) and a label-shuffle cross-check 0.983 +- 0.024 (z = -25.3), N = 1,000 draws each | `analysis/selection_null/selection_null.py` (baseline lock: `repro_baseline.py`) | `output/phase2/scaled_35types/centroids_{human,mouse}_35.csv` via `analysis/conserved_contribution/gate_lib.py`; unmodified `src/cellwarp/procrustes.py` -> `analysis/selection_null/outputs/selection_null_summary_{derangement,labelshuffle}.json`, `sigma_null_draws_*.csv`, `sigma_perms_*.npy` |

The wrapper imports the published pipeline and only re-selects its inputs. It is deterministic:
re-run against this tree it returns both draw CSVs and both sigma arrays byte-identical, and both
summary JSONs identical apart from `runtime_sec`. See `analysis/selection_null/README.md`.

`reproduce/validate.py` gates six of these values by reading the two deposited summaries: both
sigma-null means, the derangement 1st percentile, both z values, and the derangement count of
draws at or below the real value. The two pre-registered conditions are stored as booleans while
the harness compares numbers, so what is checked is the numeric backing of each: the threshold
condition 1 compares against, and the z condition 2 bounds. The real 0.384 that condition 1
compares is already gated from `gate_results.json`.

## Known gaps

- The three current-figure producers `docs/submission/plosone/figures/build_main_figures.py`, `docs/submission/plosone/figures/build_fig2c_bg.py`, and `analysis/conserved_contribution/make_figure7.py` are not invoked by `reproduce/run_all.sh`; the deposited main figures are assembled outside the full-reproduction pipeline.
- Fig 5's assembly path (`build_main_figures.py`) embeds `figures/main/fig7_conserved_contribution.png`, i.e. it depends on an artifact under `figures/main/` (outside the `docs/submission/plosone/figures/` tree) that `make_figure7.py` writes.
- `build_main_figures.py` hard-codes Fig 1D's observed distance (the `obs` constant in its `axD` block) and the `obs/null 0.35 / p < 0.0001 / n = 15` label text and `~75 Mya` title, rather than reading them from `analysis/mouse_lemur/procrustes_results.json`. Every one of those six values has been checked against that JSON and they agree, so this is a maintainability gap and not a correctness one: the panel is not showing a wrong number, but an edit to the JSON would not reach the figure.
- `Fig5_conserved_identity_genes.pdf` is a raster wrap of a PNG, even though `make_figure7.py` also emits a native vector `figures/main/fig7_conserved_contribution.pdf`.
- Fig 5D is not bit-reproducible from tracked data alone: `make_figure7.py` reads gitignored Census aggregates (`analysis/conserved_contribution/donor_stability/agg_*.npz`) for the donor-split recompute.
- The canonical `table_S7_layer1_housekeeping_exclusion.csv` has no scripted writer (the analysis script writes only the per-variant ranking CSVs in its own directory); `table_S12_software_environment.csv` is hand-authored, its authoritative source being `requirements.txt`.
- `analysis/sensitivity_analyses/markernull.py` writes its figure to `docs/supplementary_materials/figure_S8_markernull.{pdf,png}` (old stem and directory), not to the deposited `figures/submission/supplementary/figS5_markernull.pdf`.
- `scripts/build_submission_packet.py` and `tests/` still pin the old 7-figure packet (Figure_1..7, Table_S8), which does not correspond to the current display items.
- `reproduce/validate.py` has a `.csv` branch in `load_value`, but it hands `{"df": DataFrame}` to `_resolve_key`, which can only navigate dicts and lists; any key raises. No check uses it, so it is unreachable code. This is why the S1 Text §10 mechanism figure (Q75 of C collapsing to ~0.078 under derangement) is not gated: that value lives only in the per-draw CSV, while the deposited summaries carry the real Q75 alone.
