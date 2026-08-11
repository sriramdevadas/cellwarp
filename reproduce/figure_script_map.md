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
| 4C | Simulation recovery ceiling (rho ~ 0.45 at the calibrated signal) | built in `build_main_figures.py` | `analysis/simulation_study/simulation_results.json` for the deposited grid, and `analysis/simulation_study/sweep_spread_results.json` for the calibrated curve and the ceiling line the panel draws <- `simulation_study.py` and `sweep_spread.py` |

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

Three analyses produce neither a figure nor a table, so they have no row above; their numbers are
quoted in the text alone.

| Item | Content | Producing script | Depends on |
|---|---|---|---|
| S1 Text §10 | Selection/derangement circularity control: real conserved obs/null 0.384 against a derangement sigma-null 0.991 +- 0.021 (z = -29.5) and a label-shuffle cross-check 0.983 +- 0.024 (z = -25.3), N = 1,000 draws each | `analysis/selection_null/selection_null.py` (baseline lock: `repro_baseline.py`) | `output/phase2/scaled_35types/centroids_{human,mouse}_35.csv` via `analysis/conserved_contribution/gate_lib.py`; unmodified `src/cellwarp/procrustes.py` -> `analysis/selection_null/outputs/selection_null_summary_{derangement,labelshuffle}.json`, `sigma_null_draws_*.csv`, `sigma_perms_*.npy` |
| Results §8, Methods (Simulation study) | Planted-spread sweep: rank-recovery ceiling at the calibrated signal 3.683416443528231 (median rho 0.4494 at 200 cells, deposited spread sigma = 1.0), and the ceiling across a swept planted spread (0.4494 at the deposited ~65x spread up to 0.4955 at 25x), each spread point re-calibrated so obs/null stays at 0.522; zero-spread negative control 0.0059 | `analysis/simulation_study/sweep_spread.py` | `analysis/simulation_study/simulation_study.py`, exec'd unmodified for `run_single()` and `calibrate()` (its `main()` is `__name__`-guarded) -> `analysis/simulation_study/sweep_spread_results.json` |
| S14 Table (Fig 2C support) | The primate Layer-2 statistic against retained dimension: 18 rows, 3 pairs × 2 weightings × k ∈ {1,3,5}, carrying the retained joint-PCA dimension p (37/48, 40/49, 39/47), the chance level k/p, and for both rotation arms the statistic, permutation-null mean, observed-minus-null margin and p. The margin columns are the point: they are not uniform in k, and k = 1 carries the smallest margin in all twelve pair-weighting-arm combinations (smallest of all, Human-Marmoset W0 post-rotation, 0.0137). **Regenerates nothing**: the retained dimension is already the `k` field of each weighting block in the vendored results, so no basal-ganglia npz is read and no `layer2_results_*.json` is touched; every S and null is copied through and asserted against its source before the CSV is written. Gated on the table's own content rather than its inputs — the summary JSON is built by reading the written CSV back off disk, because validate.py's `.csv` branch is unreachable | `analysis/layer2_dimension_table/build_layer2_dimension_table.py` | `docs/submission/plosone/figures/bg_results/layer2_results_{Human_Macaque,Human_Marmoset,Macaque_Marmoset}.json` -> `docs/supplementary_materials/table_S14_layer2_dimension.csv` and `analysis/layer2_dimension_table/table_S14_summary.json` |
| Methods (primate replication), S1 Text | Ontology restriction of the primate Layer-2 covariance replication. Cross-species type correspondence in the basal-ganglia atlas is inherited from a consensus taxonomy built by cross-species integration, so the matching *between* types is not independent of an alignment step the way the centroids are — the safeguard argument in Methods covers the centroids only. Each pair is restricted to the types an independent Cell Ontology lookup corroborates and the permutation null is rebuilt inside the restricted subset. Two readings are carried, both flagged in the artifact: **permissive** (any CL-bearing cell_type abbreviation) leaves 50/48/48 types, **strict** (excluding the five transmitter-class terms Cholinergic, Dopa, GABA, Glut, Sero, whose CL ids name a transmitter class rather than a basal-ganglia subclass) leaves 24/20/20 against 31/32/32. Strict is used throughout because permissive discards only four or five types per pair. Findings: all 36 cells (3 pairs × 2 weightings × k ∈ {1,3,5} × pre/post) stay above the restricted null, each at the 10,000-permutation floor; the restricted margin is *wider* than the full-set margin in 34 of 36, the two exceptions narrowing by 0.003 and 0.004; and the direction reverses with k — non-corroborated types carry higher S at k = 1 in 11 of 12 combinations, corroborated types higher at k = 5 in 11 of 12. **`analysis/bg/layer2_analyze.py` is not edited and no `layer2_results_*.json` is regenerated**: the per-type diagonal that producer collapses at its lines 119-120 is reproduced in the new script, and all 36 deposited means are asserted against `np.diag()` before any restricted value is computed | `analysis/bg/ontology_restriction.py` (basal-ganglia deposit) | `analysis/hmba_data/layer2_stats{,_Human_Marmoset,_Macaque_Marmoset}.npz` plus `analysis/hmba_meta/{cluster_annotation_term,cluster_annotation_to_abbreviation_map,abbreviation_term}.csv` -> `analysis/bg/results/ontology_restriction_results.json`, vendored byte-identically to `docs/submission/plosone/figures/bg_results/ontology_restriction_results.json`, which is what `reproduce/validate.py` reads |
| Results §5, Fig 5A/5C | Detection-breadth sensitivity of the conservation score C. **The pipeline has no detection-breadth criterion for C** — `gate_lib.per_gene_corr` filters on `np.std > 0` in each species alone, which admits a gene detected in three of the 35 centroids, and `analysis/biological_predictors/biological_predictors.py:268`'s `frac_expressing > 0.10` is a per-cell rate requiring raw cell matrices that are not deposited. This analysis **supplies one for sensitivity only and does not filter any deposited result**: a gene counts as detected in a type when its centroid value there is > 0, the most generous reading, so every shortfall is a lower bound. Findings: breadth correlates *positively* with C (Spearman +0.2501), so sparsely detected genes sit at negative C rather than in the high-C tail; Fig 5A's terminal bar is the 60-bin histogram's last bin [0.9738, 1.0000] holding 221 genes at median breadth 26, not the two-or-three-centroid genes the shape suggests; the conserved set (C ≥ 0.5919, 3,985 genes) has median breadth 35 with 0.9% below breadth 5; and the master-TF enrichment survives the strictest filter — restricted to the 8,435 genes detected in all 35 types in both species, 24 of the 73 TFs, observed 0.9213 against nulls 0.5049 and 0.7149. The unfiltered row is computed first and asserted against `gate_results.json` `check3a.median_Crank` before any filtered row, so a filtered difference cannot be a reimplementation artifact | `analysis/conserved_contribution/breadth_sensitivity.py` | `output/phase2/scaled_35types/centroids_{human,mouse}_35.csv` via `analysis/conserved_contribution/gate_lib.py` (frozen `matched_draws`, `expr_bins` and `POSITIVE_CONTROL_TFS`, not re-derived) -> `analysis/conserved_contribution/breadth_sensitivity_results.json` and `breadth_per_gene.csv` |
| Results §4 | 95% confidence intervals on the four cross-atlas rank correlations, Bonett-Wright Fisher-z with `SE = 1.06/sqrt(n-3)`: Sun2023 rho +0.1464 n 15 -> [-0.424, +0.633]; PanSci +0.1941 n 16 -> [-0.362, +0.649]; pan-Census -0.0525 n 22 -> [-0.485, +0.400]; CellHint (matched 15-type baseline) -0.1393 n 15 -> [-0.629, +0.430]. Every interval spans both zero and 0.20, so the four do not separate "no cross-atlas ranking signal" from "a moderate one". The 1.06 is the Spearman constant and differs deliberately from the Pearson `1/sqrt(n-3)` in `scripts/generate_phase2_figures.py`, `scripts/t3e_step2_compute.py` and `scripts/t3e_step3b_enhancer.py` | `analysis/ranking_replication/cross_atlas_ci.py` | Four separate tracked artifacts, none of which holds more than one of the four correlations: `output/validation/sun2023_replication_expanded/sun2023_expanded.json` and `output/validation/pansci_replication/pansci_replication.json` (`rigidity_ranking.rho`, `.n_matched_types`), `analysis/census_replication/replication_results.json` (`ranking_correlation.spearman_rho`, `.n_types`), and `analysis/harmonized_replication/sensitivity_analysis.csv` (row `0_unharmonized`) -> `analysis/ranking_replication/cross_atlas_ci_results.json`. **The CellHint arm reads the CSV, not `cellhint_replication.json`:** that file's `rigidity_ranking.rho` is -0.386, the pre-PCA-matching artifact S1 Text line 48 explains away, not the -0.139 the paper reports. The CSV labels the -0.139 row `0_unharmonized`, where "unharmonized" means no ontology or tissue restriction -- a different sense of the word from S1 Text's. The producer asserts both values at the read site |

| S1 Text §2, Results §1 | Parent-and-child landmark sensitivity: re-running the primary with the broader term of each parent/child pair dropped (variant A, 30 landmarks, obs/null 0.5210) and with the more specific term dropped (variant B, 26 landmarks, obs/null 0.5441), both at the permutation floor, with rankings agreeing with the primary subset at rho 0.9359 and 0.9111 | `analysis/sensitivity/parent_child/run.py` | `output/phase2/scaled_35types/` centroids; writes `results.json` and `summary.csv`, neither of which is a deposited display item |

The wrapper imports the published pipeline and only re-selects its inputs. It is deterministic:
re-run against this tree it returns both draw CSVs and both sigma arrays byte-identical, and both
summary JSONs identical apart from `runtime_sec`. See `analysis/selection_null/README.md`.

`reproduce/validate.py` gates six of these values by reading the two deposited summaries: both
sigma-null means, the derangement 1st percentile, both z values, and the derangement count of
draws at or below the real value. The two pre-registered conditions are stored as booleans while
the harness compares numbers, so what is checked is the numeric backing of each: the threshold
condition 1 compares against, and the z condition 2 bounds. The real 0.384 that condition 1
compares is already gated from `gate_results.json`.

The spread sweep exists because the deposited `RECOVERY_SIGNALS` grid does not contain the
calibrated signal: `simulation_study.py` calibrates to 3.683416443528231, records it rounded as
`calibration.estimated_real_signal` and runs the stability experiment there, but evaluates
ranking recovery only at 3.0 and 5.0. The sweep seeds identically to that grid
(`rep + 30_000 + n_cells * 100`), so its draws are paired with it rather than independent. It is
deterministic: re-run from the repository it reproduces the earlier run on all 246 leaf values,
differing only in the wall-clock `sec` fields, which are kept as a record. Wall time about 480 s.

`reproduce/validate.py` gates six simulation values: the calibrated signal and the two grid
medians Fig 4C plots, from `simulation_results.json`; and from `sweep_spread_results.json` the
calibrated-signal median, the upper endpoint of the swept-spread range, and the zero-spread
negative control, which must stay near zero and would catch a change that broke the recovery
metric.

## Figure rasterisation and reproducibility

Nothing in this section touches a gated value. Every number `reproduce/validate.py`
checks, and every JSON, CSV, `.npy` and table workbook, is unaffected; all four gates
pass under either interpreter described below. What follows concerns only the bytes of
the deposited figure files.

Figure rendering depends on which FreeType matplotlib links, because measured text
extents set both `tight_layout` geometry and glyph rasterisation:

- The documented install path (`reproduce/README.md`: `python3.12 -m venv .venv`, then
  `.venv/bin/python`) and the Docker image use a matplotlib wheel linking its
  **vendored FreeType 2.6.1**. That is the environment
  `reproduce/environment_ground_truth.txt` was recorded from, and it reproduces the
  deposited figures' rendering. Rebuilding S4 Fig there returns an identical MediaBox,
  a PDF differing from the deposited file in 10 bytes and a PNG differing in 5 -- the
  embedded creation timestamp, and the matplotlib version string carried in `/Creator`,
  `/Producer` and the PNG `Software` chunk, which has since moved 3.10.8 to 3.10.9.
  Rebuilding S5 Fig there differs in those same 5 PNG bytes.
- A conda environment built from `environment.yml` may resolve matplotlib against a
  **system FreeType** instead (2.14.3 as measured here). The version pins above are all
  satisfied and text still measures differently. Rebuilding S4 Fig there moves the tight
  bounding box -- MediaBox `301.3628047521` becomes `301.0692252066` -- and the PDF
  differs from the deposited file in 56,694 bytes; `figure_S8_markernull.png` differs in
  30,705 of its 1,520,000 pixels (2.02%), spread over the whole panel, because every text
  element is re-rasterised. Nothing plotted moves and no gate fails.

So the deposited figures regenerate in kind under either path and byte-for-byte only
under the first -- the same distinction `DATA_SOURCES.md` draws for the UCSC refGene
tables and the DoRothEA regulon. A reader who rebuilds a figure and gets a different
file should read `matplotlib.ft2font.__freetype_version__` before looking for a
substantive cause.

Timestamps are suppressed only where a producer has been changed to do it:
`scripts/49_build_figS7_matched_scale.py` passes `metadata={"CreationDate": None}` to its
PDF `savefig`, so repeated runs in one environment are byte-identical. The other figure
producers still embed a creation timestamp and differ run to run for that reason alone.

Which interpreter built what: commit `d9a3183` rebuilt S4 Fig under the
vendored-FreeType interpreter deliberately, not under the conda environment its gates
were run in, because the conda environment would have re-rendered the whole panel.

## Known gaps

- **`reproduce/run_all.sh` rebuilds no deposited figure at all.** All eight producers that write one are absent from it: `docs/submission/plosone/figures/build_main_figures.py`, `build_fig2c_bg.py`, `analysis/conserved_contribution/make_figure7.py`, `scripts/build_submission_figures.py`, `scripts/56_add_figs2_panel_f.py`, `scripts/49_build_figS7_matched_scale.py`, `analysis/sensitivity_analyses/markernull.py` and `docs/submission/plosone/figures/build_submission_tiffs.py`. The one figure script the pipeline does invoke, `scripts/composite_figS3.py`, writes an intermediate rather than a deposited file. So the deposited figures — main and supplementary alike — are assembled outside the full-reproduction pipeline, and a reader who runs `run_all.sh` end to end regenerates none of them.
- Fig 5's assembly path (`build_main_figures.py`) embeds `figures/main/fig7_conserved_contribution.png`, i.e. it depends on an artifact under `figures/main/` (outside the `docs/submission/plosone/figures/` tree) that `make_figure7.py` writes.
- `build_main_figures.py` hard-codes Fig 1D's observed distance (the `obs` constant in its `axD` block) and the `obs/null 0.35 / p < 0.0001 / n = 15` label text and `~75 Mya` title, rather than reading them from `analysis/mouse_lemur/procrustes_results.json`. Every one of those six values has been checked against that JSON and they agree, so this is a maintainability gap and not a correctness one: the panel is not showing a wrong number, but an edit to the JSON would not reach the figure.
- `Fig5_conserved_identity_genes.pdf` is a raster wrap of a PNG, even though `make_figure7.py` also emits a native vector `figures/main/fig7_conserved_contribution.pdf`.
- Fig 5D is not bit-reproducible from tracked data alone: `make_figure7.py` reads gitignored Census aggregates (`analysis/conserved_contribution/donor_stability/agg_*.npz`) for the donor-split recompute.
- The canonical `table_S7_layer1_housekeeping_exclusion.csv` has no scripted writer (the analysis script writes only the per-variant ranking CSVs in its own directory); `table_S12_software_environment.csv` is hand-authored, its authoritative source being `requirements.txt`.
- **S2 Fig is a three-stage chain and the stages must run in order.**
  (1) `analysis/expanded_negative_controls/negative_control_figure.py` writes
  `figures/supplementary/negative_control_distributions.pdf`, which is panel E;
  (2) `scripts/build_submission_figures.py` composites panels A-E from that file plus
  the panel PNGs; (3) `scripts/56_add_figs2_panel_f.py` appends panel F, reading and
  rewriting the same deposited path. Stage 3 is now invoked by stage 2, and stage 3
  asserts on entry that its input carries exactly A-E, so it cannot append panel F
  twice. **Stage 1 is not chained and must be current before stage 2 runs.** It is in
  `reproduce/run_all.sh` at line 107, but its output is tracked, so a stale copy is
  what actually gets embedded: until D66 that tracked copy still carried a
  Cell Press house-style phrase that its own producer had already corrected, and
  rebuilding the chain would have put the phrase back into a submitted figure. If
  stage 1's wording ever changes, re-run it and commit its output before rebuilding
  S2 Fig.
- `analysis/sensitivity_analyses/markernull.py` writes its figure to `docs/supplementary_materials/figure_S8_markernull.{pdf,png}` (old stem and directory), not to the deposited `figures/submission/supplementary/figS5_markernull.pdf`. **That hop is copied by hand and enforced by nothing.** `figure_S8_markernull` occurs zero times in `scripts/build_submission_packet.py` and zero times in `tests/test_submission_packet_consistency.py`; `figS5_markernull.pdf` appears in the packet manifest only as a Group E *canonical*, whose mirror `Figure_S5.pdf` is then checked. So `--verify` confirms the deposited figure matches its packet copy while never checking that either matches what the producer wrote, and a producer re-run that is not hand-copied forward leaves the deposited S5 Fig stale with all four gates green. There is no deposited `.png` for S5 Fig; the PDF is the only deposited artifact.
- S5 Fig's panel geometry is a function of the length of the crossover annotation in panel A. `fig.tight_layout()` runs after that text is drawn and counts it in the tight bounding box, and the text overhangs panel A's right edge, so any edit changing the width of its widest line re-lays out the figure. Measured when the parenthetical was removed: panel A's axes box went `x1 = 0.464464` to `0.467754` (`x0` unchanged, so it is purely the overhang), `tight_layout` re-solved the whole 1x2 grid, and 70,064 of 1,520,000 PNG pixels moved (4.61%) -- 30,917 in panel A and 39,147 in panel B, which was not edited at all. No plotted value changes and no mark-producing operator appears or disappears (507 before and after); only coordinates transform. Expect this from any future edit to that annotation.
- `scripts/build_submission_packet.py` and `tests/` still pin the old 7-figure packet (Figure_1..7, Table_S8), which does not correspond to the current display items.
- `reproduce/validate.py` has a `.csv` branch in `load_value`, but it hands `{"df": DataFrame}` to `_resolve_key`, which can only navigate dicts and lists; any key raises. No check uses it, so it is unreachable code. This is why the S1 Text §10 mechanism figure (Q75 of C collapsing to ~0.078 under derangement) is not gated: that value lives only in the per-draw CSV, while the deposited summaries carry the real Q75 alone. (`analysis/ranking_replication/cross_atlas_ci.py` needs one CSV value and reads it in the producer, emitting JSON for the gate, which is the way round this branch until it is fixed.)
- **Two tracked artifacts disagree on the Layer-1/Layer-2 per-type correlation at k = 5, post-rotation.** `output/mechanistic/ellipsoid_alignment/35type_rigidity_correlation.csv` gives `-0.2162464985994398`; `output/validation/layer_correlation/layer_correlation_results.json` gives `-0.2207282913165266` for the same quantity. The two producers — `scripts/t3b_ellipsoid_alignment.py` (via `correlate_with_rigidity()`, writing the CSV) and `scripts/test_layer1_layer2_correlation.py` (ANALYSIS-B, writing the JSON) — agree bit-for-bit on the other five (k, metric) cells, including k = 3 pre-rotation, which is the cell S13 T31 cites (rho -0.266, p 0.123, n 35), so no submitted number depends on the disagreement. **Not diagnosed, deliberately deferred.** It matters because `t3b_ellipsoid_alignment.py` runs unguarded at `reproduce/run_all.sh:165`, so a full-pipeline run rewrites the CSV side of the pair while nothing rewrites or checks the JSON side: neither artifact is read by any of the 150 `validate.py` checks, so whichever value is wrong, no gate is watching it.
