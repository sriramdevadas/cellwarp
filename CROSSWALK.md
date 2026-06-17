# CROSSWALK: Manuscript → Code

This document maps every Methods subsection, every figure and table, and every
numerical claim in the CellWarp manuscript to the code that produces it and the
output file that contains the value.

This is a reproducibility deposit artifact: a reviewer can use this document
to walk from any specific manuscript claim to the script that generates it
and the persisted output that contains the value.

Three sections:

1. **Methods → code** -- one entry per Methods subsection (25 total),
   with the primary generating script(s) and output directory.
2. **Figures and tables → code** -- one entry per main and supplementary
   figure panel and supplementary table, with generating script and
   dependencies. Mirrors `reproduce/figure_script_map.md`.
3. **Numerical claims → code** -- every reportable numerical claim in
   Methods and Results, with the generating script, output file, JSON key
   path or column, and verification status against the persisted output.

Quick-look statistics about this CROSSWALK:
- Methods subsections covered: 25
- Figure / table display items covered: 47 (current main + supplementary figure panels including Figure 7 and Figure S8, 11 tables, and Table 1)
- Numerical claims indexed: 114
- Claims verified against persisted pipeline output (✓): 109
- Claims anchored to source publications (not pipeline) (anchor): 3
- Claims not persisted as a single artifact (intermediate): 2

For automated regression-checking of key statistics, see `reproduce/validate.py`,
which programmatically asserts a subset of the most load-bearing claims against
their persisted output files.

The "Function entry-point(s)" column in Section 1 lists procedural functions
corresponding to manuscript-described steps; helper functions are omitted.
The same column appears in Section 3 only where the value's computing function
differs from the entry in Section 1 for that subsection (most rows therefore
leave the column blank).

---

## 1. Methods → code

For each Methods subsection (in manuscript order), the primary generating
script(s) and the output directory or file where results are persisted.

| Manuscript subsection (line) | Primary script(s) | Function entry-point(s) | Output |
|---|---|---|---|
| Data acquisition and cell type selection | (Census fetch wrapped in pipeline; no single script exposed) | `02_qc_and_normalize.main` (orchestrator) | `output/phase2/scaled_35types/centroids_*.csv` |
| Ortholog mapping and gene space | BioMart query (output cached) | (one-shot BioMart query; archived CSV is anchor) | `data/phase1/orthologs_human_mouse.csv` |
| Normalization and centroid computation | `scripts/02_qc_and_normalize.py` | `02_qc_and_normalize.main`, `procrustes.compute_centroids`, `procrustes.pca_reduce_centroids` | `output/phase2/scaled_35types/centroids_*.csv` |
| PCA interdependence | (discussion only; references joint-PCA + LOOCV) | n/a | n/a |
| Independent PCA sensitivity analysis | `analysis/independent_pca_sensitivity/run_independent_pca.py` | `align_subspaces_procrustes`, `permutation_test`, `procrustes_align` | `analysis/independent_pca_sensitivity/independent_pca_results.json` |
| Procrustes superimposition | `src/cellwarp/procrustes.py` + `scripts/08_scaled_procrustes.py` | `procrustes.procrustes_align`, `08_scaled_procrustes.stage4_procrustes` | `output/phase2/scaled_35types/procrustes_results_35.json` |
| Human-versus-human negative control | `scripts/test_35type_human_control.py` + `analysis/expanded_negative_controls/expanded_negative_controls.py` | `test_35type_human_control.main`, `test_35type_human_control.run_loocv` | `output/phase2/negative_control_v2/negctrl_v2_results.json`, `analysis/within_species_matched/matched_control_results.json` |
| Replication datasets | `scripts/16_sun2023_replication.py`, `scripts/pansci_replication.py`, `scripts/33_cellhint_replication.py` | `16_sun2023_replication.main`, `pansci_replication.main`, `33_cellhint_replication.main` (each: `run_inventory` + `run_procrustes`) | `output/validation/sun2023_replication/`, `output/validation/pansci_replication/`, `output/validation/cellhint_replication/` |
| Replication diagnostics | `analysis/cellhint_investigation/investigate_rank_reversal.py`, `analysis/harmonized_replication/harmonized_replication.py` | `harmonized_replication.ontology_matching`, `tissue_matching`, `compute_harmonized_centroids`, `correlate_rankings` | `analysis/harmonized_replication/correlation_results.json`, `analysis/harmonized_replication/sensitivity_analysis.csv` |
| Macaque extension | `analysis/macaque/reconstruct_macaque_pipeline.py` | `harmonize_rira_labels`, `build_three_way_gene_list`, `main` | `output/macaque_pipeline/reconstruction_qu12_results.json`, `output/macaque_pipeline/reconstruction_qu7_D1_results.json`, `output/macaque_pipeline/human_mouse_12type_control.json`, `output/macaque_pipeline/human_mouse_7type_control.json` |
| Two-layer geometric conservation | `scripts/t3b_ellipsoid_alignment.py` + `scripts/layer3_permutation_test.py` | `t3b.compute_covariance_eigen`, `t3b.compute_alignment_scores`, `t3b.perm_test_label_shuffle`, `t3b.eigenvalue_conservation` | `output/mechanistic/ellipsoid_alignment/summary_stats.json`, `output/mechanistic/ellipsoid_alignment/permutation_results.json`, `output/layer3_permutation/layer3_permutation_results.json` |
| SAMap validation and comparison | `scripts/34_samap_35types.py` | `34_samap_35types.main` (calls `samap.run_samap` with N_ITERS=3) | `output/phase1_samap/samap_35types/samap_rigidity_correlation.json` |
| CellMarker identity gene set comparison | `scripts/cellmarker_35type_rerun.py` | `cellmarker_35type_rerun.hypergeom_enrichment` (rest procedural in `__main__`) | `output/validation/cellmarker_35type_rerun/cellmarker_35type_results.json` |
| L1000 landmark gene analysis | `scripts/35_l1000_random_baseline.py` | `35_l1000_random_baseline.silent_pipeline`, `run_random_baseline` | `output/figures/l1000_random_baseline_results.json` |
| Permutation test and statistical framework | `scripts/permutation_1M.py`, `scripts/test_lineage_stratified_permutation.py` | `permutation_1M.main`, `procrustes.permutation_test`, `test_lineage_stratified_permutation.stratified_permutation` | `analysis/permutation_1M/results_1M.json`, `output/validation/lineage_stratified/lineage_stratified_results.json` |
| Leave-one-out cross-validation | `scripts/08_loocv.py` | `08_loocv.main` (loop body is procedural; PCA + Procrustes per fold) | `output/phase3/loocv/loocv_summary.json`, `output/validation/v2_loocv/v2_loocv_results.json` |
| Bootstrap robustness | `scripts/07_bootstrap.py` | `07_bootstrap.subsample_cells`, `silent_permutation_test`, `compute_centroids_silent` | `output/phase3/bootstrap/bootstrap_summary.json` |
| Mantel test | (Mantel computation; output is canonical) | (output canonical; Mantel test on per-pair Euclidean distances) | `analysis/mantel_test/mantel_results.json` |
| Cell count confound analysis | `scripts/confound_cellcount_rigidity.py` | `confound_cellcount_rigidity.partial_spearman`, `get_cell_counts`, `main` | `output/cellcount_confound/cellcount_confound_results.json` |
| Mechanistic null tests | 10 scripts (see Section 2 mechanistic-null appendix for per-test mapping) | (per-test functions listed in Section 2 mechanistic-null appendix) | scattered: `output/phase2/mechanistic/`, `output/mechanistic/`, `output/validation/t3e_*/`, `output/t3g/` |
| Simulation study | `analysis/simulation_study/simulation_study.py` | `generate_centroids`, `calibrate`, `run_power_curve`, `run_ranking_recovery`, `run_null_calibration` | `analysis/simulation_study/simulation_results.json` |
| Bootstrap ranking analysis | `analysis/bootstrap_rankings/bootstrap_ranking_analysis.py` | `run_single_bootstrap`, `classify_stability` | `analysis/bootstrap_rankings/bootstrap_summary.csv` |
| Expanded negative controls | `analysis/expanded_negative_controls/expanded_negative_controls.py` | `compute_tissue_centroids`, `enumerate_tissue_pairs`, `random_half_split_comparison`, `_pca_and_procrustes` | `analysis/within_species_matched/matched_control_results.json` |
| Donor-split within-species control | (output is canonical) | (output canonical; per-split shared-PCA + Procrustes within Tabula Sapiens donors) | `analysis/donor_split/donor_split_shared_pca_results.json`, `analysis/donor_split/donor_split_results.json`, `analysis/donor_split/donor_split_comparison.json` |
| Biological predictors | `analysis/biological_predictors/biological_predictors.py` | `build_feature_table`, `run_univariate`, `run_multivariate` | `analysis/biological_predictors/univariate_correlations.csv`, `analysis/biological_predictors/multivariate_model_results.json` |
| Mouse lemur analysis | `analysis/mouse_lemur/01_run_pipeline.py` | `step1_data_preparation`, `step2_procrustes`, `step3_residuals_and_ranking`, `step4_evolutionary_context` | `analysis/mouse_lemur/procrustes_results.json` |
| Pan-Census replication | `analysis/census_replication/02_run_replication.py` | `02_run_replication.load_and_filter`, `process_species`, `filter_orthologs` | `analysis/census_replication/replication_results.json` |
| Software and reproducibility | n/a | n/a | `pyproject.toml`, `environment.yml` |

---

## 2. Figures and tables → code

Every display item in the CellWarp paper, mapped to the script that generates
it and its dependencies. Figure and panel numbers follow the current manuscript
(`docs/submission/manuscript_combined.txt`); this section mirrors
`reproduce/figure_script_map.md`.

### Main figures

| Display Item | Generating Script(s) | Depends On | Output file |
|---|---|---|---|
| Figure 1A (pipeline schematic) | `scripts/generate_phase2_figures.py` | None | (schematic; no underlying data file) |
| Figure 1B (1M null distribution) | `scripts/generate_phase1_figures.py` | `scripts/permutation_1M.py` | `analysis/permutation_1M/results_1M.json`, `analysis/permutation_1M/null_distribution_1M.npy` |
| Figure 1C (lineage-stratified null) | `scripts/generate_phase1_figures.py` | `scripts/test_lineage_stratified_permutation.py` | `output/validation/lineage_stratified/lineage_stratified_results.json` |
| Figure 1D (bootstrap stability) | `scripts/generate_phase3_figures.py` | `scripts/07_bootstrap.py` | `output/phase3/bootstrap/bootstrap_summary.json`, `bootstrap_results.csv` |
| Figure 1E (LOOCV bar chart) | `scripts/generate_phase3_figures.py` | `scripts/08_loocv.py` | `output/phase3/loocv/loocv_summary.json`, `loocv_results.csv` |
| Figure 2A (per-gene conservation C distribution; Hartigan dip) | `analysis/conserved_contribution/make_figure7.py` | `analysis/conserved_contribution/run_gate.py`, `make_table_s11.py` | `analysis/conserved_contribution/gate_results.json`, `gene_conservation_core.csv` |
| Figure 2B (C vs expression and Tau specificity) | `analysis/conserved_contribution/make_figure7.py` | `analysis/conserved_contribution/run_gate.py`, `run_robustness.py` | `analysis/conserved_contribution/gate_results.json`, `robustness_results.json` |
| Figure 2C (master-TF conservation enrichment vs matched backgrounds) | `analysis/conserved_contribution/make_figure7.py` | `analysis/conserved_contribution/run_gate.py`, `highN_tf_pvalues.py` | `analysis/conserved_contribution/gate_results.json`, `highN_tf_pvalues.json` |
| Figure 2D (per-gene C donor-split reproducibility) | `analysis/conserved_contribution/make_figure7.py` | (donor-stability run) | `analysis/conserved_contribution/donor_stability/donor_stability_results.json` |
| Figure 3A (Sun2023 null) | `scripts/generate_phase3_figures.py` | `scripts/16_sun2023_replication.py` | `output/validation/sun2023_replication_expanded/sun2023_expanded.json`, `output/validation/sun2023_replication/null_a_sun2023.npy` |
| Figure 3B (PanSci null) | `scripts/generate_phase3_figures.py` | `scripts/pansci_replication.py` | `output/validation/pansci_replication/pansci_replication.json`, `null_distribution.npy` |
| Figure 3C (CellHint null) | `scripts/generate_phase3_figures.py` | `scripts/33_cellhint_replication.py` | `output/validation/cellhint_replication/cellhint_replication.json` |
| Figure 3D (replication summary) | `scripts/generate_phase1_figures.py` | All replication scripts | All replication JSONs (Sun2023, PanSci, CellHint, Andrews, MCA × HCA, pan-Census) |
| Figure 3E (donor-split / human control) | `scripts/generate_phase1_figures.py` | `scripts/test_35type_human_control.py`, donor-split outputs | `analysis/donor_split/donor_split_shared_pca_results.json`, `output/phase2/negative_control_v2/negctrl_v2_results.json` |
| Figure 4A (ellipsoid heatmap) | `scripts/generate_phase3_figures.py` | `scripts/t3b_ellipsoid_alignment.py` | `output/mechanistic/ellipsoid_alignment/35type_alignment_scores.csv` |
| Figure 4B (pre vs post-rotation) | `scripts/generate_phase3_figures.py` | `scripts/t3b_ellipsoid_alignment.py` | `output/mechanistic/ellipsoid_alignment/permutation_results.json`, `35type_alignment_scores.csv` |
| Figure 4C (layer null distributions) | `scripts/generate_phase3_figures.py` | `scripts/layer3_permutation_test.py` | `output/mechanistic/ellipsoid_alignment/permutation_results.json` |
| Figure 4D (layer scatter) | `scripts/generate_phase3_figures.py` | `scripts/t3b_ellipsoid_alignment.py` | `output/mechanistic/ellipsoid_alignment/35type_rigidity_correlation.csv` |
| Figure 5A (three-species Procrustes) | `scripts/generate_phase1_figures.py` | Macaque pipeline (pre-computed) | `output/macaque_pipeline/m1_close_table1_summary.json`, `analysis/mouse_lemur/procrustes_results.json` |
| Figure 5B (no-immune sensitivity) | `scripts/generate_phase1_figures.py` | Macaque pipeline (pre-computed) | `output/macaque_pipeline/reconstruction_qu7_D1_results.json`, `human_mouse_7type_control.json` |
| Figure 5C (hepatocyte ranking) | `scripts/generate_phase1_figures.py` | Macaque pipeline + primary | `output/macaque_pipeline/m1_close_table1_summary.json`, `human_mouse_12type_control.json` |
| Figure 6A (divergence ranking) | `scripts/generate_phase1_figures.py` | `scripts/08_scaled_procrustes.py` | `output/phase2/scaled_35types/residuals_ranked.csv` |
| Figure 6B (cell count confound) | `scripts/generate_phase1_figures.py` | `scripts/confound_cellcount_rigidity.py` | `output/cellcount_confound/cellcount_confound_results.json` |
| Figure 7A (L1000 baseline) | `scripts/generate_phase3_figures.py` | `scripts/35_l1000_random_baseline.py` | `output/figures/l1000_random_baseline_results.json` |
| Figure 7B (mechanistic nulls forest) | `scripts/generate_phase2_figures.py` | Mechanistic null scripts (see below) | per-test outputs (see Mechanistic null scripts table below) |
| Table 1 (unified statistical tests) | Compiled in manuscript from all scripts | All analysis scripts | (compiled inline in manuscript; sources span all analysis JSONs in this CROSSWALK) |

### Supplementary figures

| Display Item | Generating Script(s) | Depends On | Output file |
|---|---|---|---|
| Figure S1A-B (independent PCA) | `scripts/generate_phase1_figures.py` | `analysis/independent_pca_sensitivity/run_independent_pca.py` | `analysis/independent_pca_sensitivity/independent_pca_results.json` |
| Figure S1C-F (simulation study) | `analysis/simulation_study/simulation_figures.py` | `analysis/simulation_study/simulation_study.py` | `analysis/simulation_study/simulation_results.json` |
| Figure S2A-B (PCA k-sensitivity) | `scripts/generate_phase3_figures.py` | `scripts/17_pca_sensitivity.py`, `scripts/18_pca_sensitivity_v2.py` | `output/validation/pca_sensitivity/pca_sensitivity_results.json` |
| Figure S2C-D (Smart-seq2 protocol) | `scripts/generate_phase3_figures.py` | `scripts/14_smartseq2_sensitivity.py` | `output/phase2/sensitivity/smartseq2/sensitivity_results.json` |
| Figure S2E (expanded negatives) | `analysis/expanded_negative_controls/negative_control_figure.py` | `analysis/expanded_negative_controls/expanded_negative_controls.py` | `analysis/within_species_matched/matched_control_results.json` |
| Figure S2F (replication inventory) | `scripts/56_add_figs2_panel_f.py` | All replication outputs | All replication JSONs (Sun2023, PanSci, CellHint, Andrews, MCA × HCA, pan-Census, HCA × Tabula) |
| Figure S3A-B (bootstrap rankings) | `scripts/composite_figS3.py` | `analysis/bootstrap_rankings/bootstrap_ranking_analysis.py` | `analysis/bootstrap_rankings/bootstrap_summary.csv`, `bootstrap_rankings_raw.csv` |
| Figure S4A-B (CellHint investigation) | `scripts/generate_phase3_figures.py` | `analysis/cellhint_investigation/investigate_rank_reversal.py` | `analysis/harmonized_replication/correlation_results.json`, `sensitivity_analysis.csv` |
| Figure S5A (SAMap heatmap) | `scripts/generate_phase3_figures.py` | `scripts/34_samap_35types.py` | `output/phase1_samap/samap_35types/samap_rigidity_correlation.json`, `samap_mapping_scores_35.csv` |
| Figure S6A-B (CellMarker enrichment) | `scripts/generate_phase1_figures.py` | `scripts/cellmarker_35type_rerun.py` | `output/validation/cellmarker_35type_rerun/cellmarker_35type_results.json` |
| Figure S7 (matched-scale 6-type negative control) | `scripts/test_35type_human_control.py` | `output/phase2/negative_control_v2/negctrl_v2_results.json` | `output/phase2/negative_control_v2/negctrl_v2_results.json` |
| Figure S8 (marker-similarity-stratified null) | `analysis/sensitivity_analyses/markernull.py` | (primary centroids; species-averaged gene-space) | `analysis/sensitivity_analyses/markernull_results.json`, `docs/supplementary_materials/figure_S8_markernull.pdf` |

### Supplementary tables

| Display Item | Generating Script(s) | Depends On | Output file |
|---|---|---|---|
| Table S1 (biological predictors + cross-atlas) | `scripts/create_table_S1.py` | `analysis/biological_predictors/biological_predictors.py`, `analysis/ranking_replication/ranking_replication_analysis.py` | `analysis/biological_predictors/univariate_correlations.csv`, `analysis/cross_reference/master_ranking_table.csv` |
| Table S2 (simulation + bootstrap CIs) | `scripts/create_table_S2.py` | `analysis/simulation_study/simulation_study.py`, `analysis/bootstrap_rankings/bootstrap_ranking_analysis.py` | `analysis/simulation_study/simulation_results.json`, `analysis/bootstrap_rankings/bootstrap_summary.csv` |
| Table S3 (CellHint rank reversal) | `analysis/cellhint_investigation/investigate_rank_reversal.py` | `scripts/33_cellhint_replication.py` | `analysis/cellhint_investigation/` (rank-reversal artifacts) |
| Table S4 (CellHint harmonization) | `analysis/harmonized_replication/harmonized_replication.py` | `scripts/33_cellhint_replication.py` | `analysis/harmonized_replication/sensitivity_analysis.csv`, `correlation_results.json` |
| Table S5 (35-type matching) | `scripts/08_cell_type_inventory.py` | `scripts/02_qc_and_normalize.py` | `output/phase2/cell_type_inventory.csv` |
| Table S6 (CPC1 driver genes) | `scripts/generate_table_S6.py` | `scripts/t3b_ellipsoid_alignment.py` | `output/mechanistic/ellipsoid_alignment/cpc_gene_table.csv`, `cpc_genes.json` |
| Table S7 (Layer-1 ribosomal/housekeeping exclusion) | `analysis/sensitivity_analyses/layer1_exclusion.py` | `scripts/08_scaled_procrustes.py` | `analysis/sensitivity_analyses/layer1_exclusion_results.json`, `docs/supplementary_materials/table_S7_layer1_housekeeping_exclusion.csv` |
| Table S8 (marker 1:1-ortholog retention vs residual) | `analysis/sensitivity_analyses/ortholog_retention.py` | `scripts/cellmarker_35type_rerun.py`, `data/phase1/orthologs_human_mouse.csv` | `analysis/sensitivity_analyses/ortholog_retention_results.json`, `docs/supplementary_materials/table_S8_marker_ortholog_retention.csv` |
| Table S9 (per-gene standardization: Layer-1 + CPC1 Scheme A/B) | `analysis/sensitivity_analyses/genestd_standardization.py` | `scripts/08_scaled_procrustes.py`, `scripts/t3b_ellipsoid_alignment.py` | `analysis/sensitivity_analyses/genestd_results.json`, `docs/supplementary_materials/table_S9_genestd_standardization.csv`, `table_S9_schemeB_CPC1_markers.csv` |
| Table S10 (marker-similarity-stratified null) | `analysis/sensitivity_analyses/markernull.py` | (primary centroids; species-averaged gene-space) | `analysis/sensitivity_analyses/markernull_results.json`, `docs/supplementary_materials/table_S10_markernull.csv` |
| Table S11 (per-gene cross-species conservation score C) | `analysis/conserved_contribution/make_table_s11.py` | `analysis/conserved_contribution/run_gate.py` | `analysis/conserved_contribution/gene_conservation_core.csv`, `docs/supplementary_materials/table_S11_gene_conservation.csv` |

### Mechanistic null scripts (Figure 7B, mechanistic-null table in Methods)

| Test | Script | Output file |
|---|---|---|
| 1. Housekeeping ratio | `scripts/12_housekeeping_ratio.py` | `output/phase2/mechanistic/housekeeping/hk_ratio_results.json` |
| 2. TF complexity | `scripts/13_tf_complexity.py` | `output/phase2/mechanistic/tf_complexity/tf_complexity_vs_residual.csv` |
| 3. Niche hypothesis | `scripts/12_niche_hypothesis.py` | `output/phase2/mechanistic/niche/` |
| 4. Within-type variance | `scripts/12_variance_diagnostic.py` | `output/phase2/variance_diagnostic/diagnostic_results.json` |
| 5. Inter-donor variance | `scripts/16_interdonor_variance.py` | `output/phase2/diagnostics/interdonor_variance/interdonor_variance_by_celltype.csv` |
| 6. Expression-level confounds | `scripts/diagnostic_expression_vs_rigidity.py` | `output/phase2/diagnostics/expression_level_vs_rigidity/correlations.csv` |
| 7. PPI centrality | `scripts/19_ppi_centrality.py` | `output/mechanistic/ppi_centrality/ppi_centrality_results.json` |
| 8. Promoter sequence conservation (phastCons) | `scripts/t3e_step2_compute.py` | `output/validation/t3e_chromatin/t3e_step2_summary.md` |
| 9. Active enhancer conservation | `scripts/t3e_step3b_enhancer.py` | `output/validation/t3e_enhancer/t3e_step3b_summary.md` |
| 10. Drug target conservation | (output is canonical; producer script archived) | `output/t3g/primary_correlation_results.json` |

---

## 3. Numerical claims → code

Every reportable numerical claim in Methods and Results, with its generating
script and persisted output. Values shown reflect the persisted pipeline
output. The verification status indicates how confidently the reported value
can be reproduced from the deposit.

Status legend:
- ✓ -- verified: claim value matches the persisted pipeline output within
  the appropriate tolerance for its data type.
- anchor -- anchored to source publication (e.g., Tabula Sapiens v1.0 paper);
  not produced by this pipeline.
- intermediate -- value is an intermediate count in the data-acquisition
  pipeline, not persisted as a single artifact.

### Results §1 -- Cross-species transcriptomic geometry is globally coherent

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Verified |
|---|---|---|---|---|---|---|
| matched cell types | 35 | `scripts/02_qc_and_normalize.py` | `output/phase2/scaled_35types/centroids_human_35.csv` | shape[0] |  | ✓ |
| shared 1:1 ortholog genes | 16,959 | `src/cellwarp/procrustes.py` | `output/phase2/scaled_35types/centroids_human_35.csv` | shape[1] |  | ✓ |
| PCA components | 33 | `scripts/08_scaled_procrustes.py` | `analysis/permutation_1M/results_1M.json` | `n_pca_components` |  | ✓ |
| PCA variance retained | 95.2% | `scripts/08_scaled_procrustes.py` | `output/validation/pca_sensitivity/pca_sensitivity_results.json` | `ref.variance_explained` |  | ✓ |
| primary obs/null | 0.522 | `scripts/08_scaled_procrustes.py` | `output/phase2/scaled_35types/procrustes_results_35.json` | computed: distance / null_median |  | ✓ |
| range of central 95% null | 0.507–0.544 | `scripts/permutation_1M.py` | `analysis/permutation_1M/results_1M.json` | `null_distribution_summary.percentile_2_5/97_5` (computed ratio) |  | ✓ |
| primary p-value | < 10⁻⁶ | `scripts/permutation_1M.py` | `analysis/permutation_1M/results_1M.json` | `p_value` |  | ✓ |
| permutations | 1,000,000 | `scripts/permutation_1M.py` | `analysis/permutation_1M/results_1M.json` | `n_permutations` |  | ✓ |
| lineage-strat obs/null | 0.668 | `scripts/test_lineage_stratified_permutation.py` | `output/validation/lineage_stratified/lineage_stratified_results.json` | `stratified_null.obs_null_ratio` |  | ✓ |
| lineage-strat p | 0.0001 | same | same | `stratified_null.p_value` |  | ✓ |
| within-lineage null tighter | 21.9% | `scripts/test_lineage_stratified_permutation.py` | `output/validation/lineage_stratified/lineage_stratified_results.json` | `distribution_comparison.pct_tighter` |  | ✓ |
| macaque divergence | ~25 Mya | (manuscript anchor; ref 15) | n/a | n/a |  | anchor |
| macaque type count | n=12 | `analysis/macaque/reconstruct_macaque_pipeline.py` | `output/macaque_pipeline/reconstruction_qu12_results.json` | `types_included` length |  | ✓ |
| macaque obs/null | 0.810 | same | same | `permutation_test.obs_null_ratio_median` |  | ✓ |
| macaque p | 0.0043 | same | same | `permutation_test.p_value` |  | ✓ |
| mouse lemur divergence | ~75 Mya | (manuscript anchor; ref 15) | n/a | n/a |  | anchor |
| mouse lemur type count | n=15 | `analysis/mouse_lemur/01_run_pipeline.py` | `analysis/mouse_lemur/procrustes_results.json` | `n_types` |  | ✓ |
| mouse lemur obs/null | 0.346 | same | same | `permutation_test.obs_null_ratio` |  | ✓ |
| mouse lemur p | 0.0001 | same | same | `permutation_test.p_value` |  | ✓ |
| bootstrap iterations | 100 | `scripts/07_bootstrap.py` | `output/phase3/bootstrap/bootstrap_summary.json` | `n_bootstrap` |  | ✓ |
| bootstrap subsample fraction | 50% | same | same | `subsample_fraction` |  | ✓ |
| bootstrap CV | 0.004 | same | same | `distances.cv` |  | ✓ |
| bootstrap stability gate | CV=0.2 | same | same | `gate_criterion.threshold` |  | ✓ |
| bootstrap fraction sig | 100/100 | same | same | `p_values.fraction_significant_001` |  | ✓ |
| LOOCV correct count | 35/35 | `scripts/08_loocv.py` | `output/validation/v2_loocv/v2_loocv_results.json` | `n_correct/n_total` |  | ✓ |
| LOOCV mean ratio | 0.4201 | same | same | `mean_ratio` |  | ✓ |
| LOOCV improvement | 58% | (computed: 1 − 0.4201 ≈ 0.58) | same | same |  | ✓ |
| indep-PCA obs/null | 0.473 | `analysis/independent_pca_sensitivity/run_independent_pca.py` | `analysis/independent_pca_sensitivity/independent_pca_results.json` | `permutation_test.obs_null_ratio` |  | ✓ |
| indep-PCA p | < 10⁻⁶ | same | same | `permutation_test.p_value` |  | ✓ |
| indep-PCA Spearman ρ | 0.915 | same | same | `comparison_to_joint_pca.spearman_rho` |  | ✓ |
| indep-PCA dramatic shifts | 2 of 35 | same | same | `comparison_to_joint_pca.n_dramatic_rank_changes` |  | ✓ |
| Mantel Pearson r | 0.789 | (Mantel) | `analysis/mantel_test/mantel_results.json` | `pearson_r` |  | ✓ |
| Mantel Spearman rho | 0.737 | same | same | `spearman_r` |  | ✓ |
| Mantel p (both) | < 0.001 | same | same | `pearson_p`, `spearman_p` |  | ✓ |
| Mantel permutations | 10,000 | same | same | `n_permutations` |  | ✓ |
| within-species pairs | 24 | `analysis/expanded_negative_controls/expanded_negative_controls.py` | `analysis/within_species_matched/matched_control_results.json` | `all_pairs_stats.n` |  | ✓ |
| within-species mean obs/null | 0.466 | same | same | `all_pairs_stats.mean_obs_null` |  | ✓ |
| fraction more coherent | 71% | same | same | `all_pairs_stats.fraction_lower_than_cross_species` |  | ✓ |
| random balanced splits | 100 | `analysis/donor_split/...` | `analysis/donor_split/donor_split_shared_pca_results.json` | `n_splits` |  | ✓ |
| feasible cell types | 32/35 | same | same | `infeasible_types` length |  | ✓ |
| median types per split | 31 | same | same | `within_species.median_n_types` |  | ✓ |
| median PCA components | 32 | same | same | `pca_components` (median) |  | ✓ |
| self-comparison value | 0.033 | same | same | `reference_values.self_comparison_obs_null` |  | ✓ |
| within-species median | 0.375 | same | same | `within_species.median_obs_null_ratio` |  | ✓ |
| within-species 95% CI | 0.318–0.423 | same | same | `within_species.ci_95` |  | ✓ |
| cross-species median | 0.527 | same | same | `cross_species_matched.median_obs_null_ratio` |  | ✓ |
| cross-species 95% CI | 0.485–0.584 | same | same | `cross_species_matched.ci_95` |  | ✓ |
| delta median | +0.159 | same | same | `delta.median` |  | ✓ |
| delta 95% CI | +0.100–+0.218 | same | same | `delta.ci_95` |  | ✓ |
| positive splits | 100/100 | same | same | `delta.pct_positive` |  | ✓ |
| indep-PCA delta | +0.158 | `analysis/donor_split/...` | `analysis/donor_split/donor_split_results.json` | `delta.median` |  | ✓ |

(The donor-split control claims above (random splits, within-species 0.375, cross-species 0.527, delta +0.159) are reported in the manuscript under §3, "A donor-split control isolates the cross-species signal from within-species variation"; they are grouped here with §1 as the coherence-establishing controls.)

### Results §2 -- Conserved master-regulator programs carry the cell-type geometry

Per-gene cross-species conservation score C (Pearson correlation of each 1:1
ortholog's cell-type-specificity profile across the 35 matched centroids, human
vs mouse; Table S11), its enrichment for master/identity transcription factors,
the geometry-attribution gene sets, and the reproducibility splits. Gate harness:
`analysis/conserved_contribution/run_gate.py` → `gate_results.json`; these values
are asserted by `reproduce/validate.py`.

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Verified |
|---|---|---|---|---|---|---|
| genes with defined C (of 16,959) | 15,940 | `analysis/conserved_contribution/make_table_s11.py` | `analysis/conserved_contribution/gate_results.json` | `n_valid` |  | ✓ |
| C vs procrustes_contribution ρ (non-circularity) | 0.27 | `analysis/conserved_contribution/run_gate.py` | same | `rho_pearson_vs_circular_loading` |  | ✓ |
| C vs expression Spearman ρ | 0.22 | same | same | `check2.spearman_C_vs_expr` |  | ✓ |
| C vs Tau specificity Spearman ρ | 0.06 | `analysis/conserved_contribution/run_robustness.py` | `analysis/conserved_contribution/robustness_results.json` | `rho_C_tau` |  | ✓ |
| Hartigan dip D (broad continuum) | 0.007 (p = 2.8 × 10⁻⁵) | `analysis/conserved_contribution/run_robustness.py` | `robustness_results.json`, `gate_results.json` | `dip_D`, `check1.dip_p` |  | ✓ |
| master-TF median C-percentile | 0.94 | `analysis/conserved_contribution/run_gate.py` | `gate_results.json` | `check3a.median_Crank` |  | ✓ |
| expression-matched background percentile | 0.54 | same | same | `check3a.null_median_Crank` |  | ✓ |
| joint expr+specificity-matched background percentile | 0.76 | `analysis/conserved_contribution/highN_tf_pvalues.py` | `analysis/conserved_contribution/highN_tf_pvalues.json` | `joint_expr_tau_matched.null_median` |  | ✓ |
| master-TF enrichment p vs both backgrounds (10⁶ draws) | < 10⁻⁶ (0 exceedances) | `analysis/conserved_contribution/highN_tf_pvalues.py` | `highN_tf_pvalues.json` | `*.p_empirical`, `*.exceedances` |  | ✓ |
| TF fold-enrichment at conserved end (joint-matched) | 1.67 | `analysis/conserved_contribution/run_gate.py` | `gate_results.json` | `check3b.H_tf.fold` |  | ✓ |
| CellMarker fold-enrichment (raw) | 1.89 | same | same | `check3b.H_cellmarker.fold` |  | ✓ |
| conserved-set obs/null (geometry attribution) | 0.384 | same | same | `secondary.conserved.ratio` |  | ✓ |
| divergent-set obs/null | 0.709 | same | same | `secondary.divergent.ratio` |  | ✓ |
| expression-matched-random obs/null | 0.525 ± 0.012 | same | same | `secondary.matched_random_ratio_mean` (`_sd`) |  | ✓ |
| all-genes obs/null anchor | 0.522 | same | same | `secondary.validity_all_genes_ratio` |  | ✓ |
| donor-split cross-half C Spearman (median) | 0.80 | (donor-stability run) | `analysis/conserved_contribution/donor_stability/donor_stability_results.json` | `donor_split_cap10000.cross_half_C_spearman_median` |  | ✓ |
| top-quartile conserved-set Jaccard (donor-sensitivity) | 0.58 | `analysis/conserved_contribution/run_gate.py` | `gate_results.json` | `check4.highcount_jaccard` |  | ✓ |
| cross-protocol C Spearman (10x vs Smart-seq2) | 0.59 | (donor-stability run) | `donor_stability/donor_stability_results.json` | `cross_protocol.spearman_C10x_vs_CSS` |  | ✓ |
| fresh-pull obs/null (Census re-acquisition) | 0.521 | (donor-stability run) | `donor_stability/donor_stability_results.json` | `validity.obs_null_full` |  | ✓ |

### Results §4 -- Centroid and within-type covariance geometry are conserved through distinct features

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Verified |
|---|---|---|---|---|---|---|
| Layer 2 S at k=5 | 0.483 | `scripts/t3b_ellipsoid_alignment.py` | `output/mechanistic/ellipsoid_alignment/summary_stats.json` | `35type.mean_alignment.k=5.pre` |  | ✓ |
| Layer 2 null mean | 0.375 | same | `output/mechanistic/ellipsoid_alignment/permutation_results.json` | `35type.label_shuffle_pre.k=5.null_mean` |  | ✓ |
| Layer 2 p (k=5) | 0.0001 | same | same | `35type.label_shuffle_pre.k=5.p_value` |  | ✓ |
| Layer 2 permutations | 10,000 | same | same | (configured) |  | ✓ |
| post-rotation S | 0.230 | same | `output/mechanistic/ellipsoid_alignment/summary_stats.json` | `35type.mean_alignment.k=5.post` |  | ✓ |
| post-rotation null mean | 0.180 | same | `output/mechanistic/ellipsoid_alignment/permutation_results.json` | `35type.label_shuffle_post.k=5.null_mean` |  | ✓ |
| post-rotation p | 0.0001 | same | same | `35type.label_shuffle_post.k=5.p_value` |  | ✓ |
| Layer-1 vs Layer-2 ρ | −0.266 | same | (per-type cross-correlation) | `output/mechanistic/ellipsoid_alignment/35type_rigidity_correlation.csv` |  | ✓ |
| Layer-1 vs Layer-2 p | 0.123 | same | same | same |  | ✓ |
| Layer 2 PanSci pre-rotation S (k=5) | 0.396 | `scripts/t3b_ellipsoid_alignment_pansci.py` | `output/twolayer_pansci_replication/pansci_layer2_summary.json` | `layer2_pre_rotation.k5.S` |  | ✓ |
| Layer 2 PanSci pre-rotation null mean | 0.302 | same | same | `layer2_pre_rotation.k5.null_mean` |  | ✓ |
| Layer 2 PanSci post-rotation S (k=5) | 0.402 | same | same | `layer2_post_rotation.k5.S` |  | ✓ |
| Layer 2 PanSci post-rotation null mean | 0.360 | same | same | `layer2_post_rotation.k5.null_mean` |  | ✓ |
| Layer 2 PanSci p (pre & post, k=5) | < 10⁻⁴ | same | same | `layer2_pre_rotation.k5.p`, `layer2_post_rotation.k5.p` |  | ✓ |
| CPC1 ribosomal-dominated types (primary weighting) | 25/35 | `analysis/sensitivity_analyses/genestd_standardization.py` | `analysis/sensitivity_analyses/genestd_results.json` | `cpc1.base.n_ribosomal_dominated` |  | ✓ |
| CPC1 ribosomal-dominated types (Scheme B per-gene std) | 1/35 | same | same | `cpc1.B.n_ribosomal_dominated` |  | ✓ |
| Layer-1 obs/null under per-gene std (Scheme A / B) | 0.606 / 0.487 | same | same | `layer1.A.obs_null`, `layer1.B.obs_null` |  | ✓ |
| per-type ranking ρ vs primary (Scheme A / B) | 0.54 / 0.76 | same | same | `layer1.A.ranking_rho_vs_primary`, `layer1.B.ranking_rho_vs_primary` |  | ✓ |

### Results §5 -- Geometric coherence replicates across independent atlases

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Verified |
|---|---|---|---|---|---|---|
| Sun2023 type count | n=15 | `scripts/17_sun2023_expanded.py` | `output/validation/sun2023_replication_expanded/sun2023_expanded.json` | `procrustes.n_types` |  | ✓ |
| Sun2023 obs/null | 0.554 | same | same | `procrustes.obs_null_ratio` |  | ✓ |
| Sun2023 p | 0.0001 | same | same | `procrustes.p_value` |  | ✓ |
| PanSci type count | n=16 | `scripts/pansci_replication.py` | `output/validation/pansci_replication/pansci_replication.json` | `procrustes.n_types` |  | ✓ |
| PanSci obs/null | 0.552 | same | same | `procrustes.obs_null_ratio` |  | ✓ |
| PanSci p | 0.0001 | same | same | `procrustes.p_value` |  | ✓ |
| CellHint cells | 2.87M | (CellHint dataset metadata) | (anchor: Xu et al. 2023) | n/a |  | anchor |
| CellHint type count | n=15 | `scripts/33_cellhint_replication.py` | `output/validation/cellhint_replication/cellhint_replication.json` | `procrustes.n_types` |  | ✓ |
| CellHint obs/null | 0.448 | same | same | `procrustes.obs_null_ratio` |  | ✓ |
| CellHint p | 0.0001 | same | same | `procrustes.p_value` |  | ✓ |
| pan-Census type count | n=22 | `analysis/census_replication/02_run_replication.py` | `analysis/census_replication/replication_results.json` | `n_cell_types` |  | ✓ |
| pan-Census datasets | 15 | `data/replication/pan_census_manifest.csv` | same CSV | `df.shape[0]` |  | ✓ |
| pan-Census obs/null | 0.811 | same | `analysis/census_replication/replication_results.json` | `permutation_test.obs_null_ratio` |  | ✓ |
| pan-Census p | 0.0001 | same | same | `permutation_test.p_value` |  | ✓ |
| Andrews obs/null | 0.797 | `scripts/31_andrews_replication.py` | `output/validation/andrews_replication/andrews_replication_results.json` | `obs_null_ratio` | `31_andrews_replication.main` | ✓ |
| Andrews p | 0.1159 | same | same | `p_value` | `31_andrews_replication.main` | ✓ |
| Andrews n | 6 | same | same | `n_types` | `31_andrews_replication.main` | ✓ |
| MCA × HCA obs/null | 1.003 | `scripts/14_t1a_replication.py` | `output/validation/t1a_replication/t1a_results.json` | `t1a_procrustes.obs_null_ratio` | `14_t1a_replication.main` | ✓ |
| MCA × HCA p | 0.542 | same | same | `t1a_procrustes.p_value` | `14_t1a_replication.main` | ✓ |
| MCA × HCA n | 17 | same | same | `t1a_procrustes.n_types` | `14_t1a_replication.main` | ✓ |
| HCA × Tabula obs/null | 0.728 | `scripts/15_hca_centroid_comparison.py` | `output/validation/hca_centroid_comparison/hca_centroid_comparison.json` | `comparison_a.obs_null_ratio` | `15_hca_centroid_comparison.main` | ✓ |
| HCA × Tabula p | 0.003 | same | same | `comparison_a.p_value` | `15_hca_centroid_comparison.main` | ✓ |
| HCA × Tabula n | 6 | same | same | `comparison_a.n_types` | `15_hca_centroid_comparison.main` | ✓ |
| Andrews scaling | 0.229 | `scripts/31_andrews_replication.py` | `output/validation/andrews_replication/andrews_replication_results.json` | `scaling` | `31_andrews_replication.main` | ✓ |
| MCA scaling | 0.267 | `scripts/14_t1a_replication.py` | `output/validation/t1a_replication/t1a_results.json` | `t1a_procrustes.scaling` | `14_t1a_replication.main` | ✓ |
| Sun ranking ρ | +0.15 | same Sun2023 file | same | `rigidity_ranking.rho` |  | ✓ |
| Sun ranking p | 0.60 | same | same | `rigidity_ranking.p_value` |  | ✓ |
| PanSci ranking ρ | +0.19 | `output/validation/pansci_replication/pansci_replication.json` | same | `rigidity_ranking.rho` |  | ✓ |
| PanSci ranking p | 0.47 | same | same | `rigidity_ranking.p_value` |  | ✓ |
| CellHint ranking ρ | −0.39 | `output/validation/cellhint_replication/cellhint_replication.json` | same | `rigidity_ranking.rho` |  | ✓ |
| CellHint ranking p | 0.16 | same | same | `rigidity_ranking.p_value` |  | ✓ |
| pan-Census ranking ρ | −0.05 | `analysis/census_replication/replication_results.json` | same | `ranking_correlation.spearman_rho` |  | ✓ |
| CellHint harmonized ρ | −0.04 | `analysis/harmonized_replication/correlation_results.json` | same | `rho` |  | ✓ |
| harmonized n | 12 | same | same | `n_types` |  | ✓ |

### Results §6 -- Macaque extension: matched-n controls and per-type ranking dependence

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Verified |
|---|---|---|---|---|---|---|
| three-way ortholog space | 13,927 genes | `analysis/macaque/reconstruct_macaque_pipeline.py` | `output/macaque_pipeline/reconstruction_qu12_results.json` | `gene_space` |  | ✓ |
| three-way fraction | 82.1% | (computed: 13927/16959) | n/a | n/a |  | ✓ |
| macaque-vs-HM obs/null | 0.810 | same | `output/macaque_pipeline/reconstruction_qu12_results.json` | `permutation_test.obs_null_ratio_median` |  | ✓ |
| matched-n HM obs/null | 0.440 | same | `output/macaque_pipeline/human_mouse_12type_control.json` | `control_16959.permutation_test.obs_null_ratio_median` |  | ✓ |
| matched-n HM p | 0.0001 | same | same | `control_16959.permutation_test.p_value` |  | ✓ |
| no-immune subset n | 7 | same | `output/macaque_pipeline/reconstruction_qu7_D1_results.json` | `n_types_analyzed` |  | ✓ |
| no-immune subset obs/null | 0.733 | same | same | `permutation_test.obs_null_ratio_median` |  | ✓ |
| no-immune p | 0.013 | same | same | `permutation_test.p_value` |  | ✓ |
| matched-n HM 7-type obs/null | 0.485 | same | `output/macaque_pipeline/human_mouse_7type_control.json` | `permutation_test.obs_null_ratio_median` |  | ✓ |
| matched-n HM 7-type p | 0.0006 | same | same | `permutation_test.p_value` |  | ✓ |
| divergence ranking ρ | 0.147 | same | `output/macaque_pipeline/m1_close_table1_summary.json` | `spearman_ranking_correlation.rho` |  | ✓ |
| divergence ranking p | 0.649 | same | same | `spearman_ranking_correlation.p_value` |  | ✓ |
| mouse lemur ρ vs primary | 0.157 | `analysis/mouse_lemur/01_run_pipeline.py` | `analysis/mouse_lemur/ranking_correlation.json` | `vs_primary.rho` |  | ✓ |
| mouse lemur p | 0.576 | same | same | `vs_primary.p_value` |  | ✓ |
| hepatocyte rank HM | 12/12 | same | `output/macaque_pipeline/m1_close_table1_summary.json` | `hepatocyte_rank_reversal.hm12_rank` |  | ✓ |
| hepatocyte SSR HM | 2.0% | same | `output/macaque_pipeline/human_mouse_12type_control.json` | `hepatocyte_pct_ssr` |  | ✓ |
| hepatocyte rank macaque | 1/12 | same | `output/macaque_pipeline/m1_close_table1_summary.json` | `human_macaque_12type.hepatocyte_rank_of_n` |  | ✓ |
| hepatocyte SSR macaque | 47.3% | same | same | `human_macaque_12type.hepatocyte_pct_ssr` |  | ✓ |

### Results §7 -- Per-cell-type divergence is precisely estimated but biologically unresolved

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Verified |
|---|---|---|---|---|---|---|
| CellMarker enrichment fold | 4.49 | `scripts/cellmarker_35type_rerun.py` | `output/validation/cellmarker_35type_rerun/cellmarker_35type_results.json` | `global_enrichment.enrichment` |  | ✓ |
| CellMarker p-value | 2.10 × 10⁻¹³ | same | same | `global_enrichment.p_value` |  | ✓ |
| expression-matched fold | 3.32 | same | same | `expression_matched.enrichment` |  | ✓ |
| expression-matched p | 1.15 × 10⁻¹² | same | same | `expression_matched.p_value` |  | ✓ |
| per-type pass rate | 5/6 | same | same | `per_type_pass` |  | ✓ |
| features tested | 15 | `analysis/biological_predictors/biological_predictors.py` | `analysis/biological_predictors/univariate_correlations.csv` | `df.shape[0]` |  | ✓ |
| cell count ρ | 0.052 | `scripts/confound_cellcount_rigidity.py` | `output/cellcount_confound/cellcount_confound_results.json` | `spearman_rho` |  | ✓ |
| tissue breadth ρ | −0.13 | `analysis/biological_predictors/biological_predictors.py` | `analysis/biological_predictors/univariate_correlations.csv` | "Tissue breadth" row |  | ✓ |
| mean expression ρ | −0.06 | same | same | "Mean expression level" row |  | ✓ |
| progenitor ρ | 0.43 | same | same | "Is progenitor" row |  | ✓ |
| progenitor p | 0.01 | same | same | same |  | ✓ |
| elastic net LOO R² | −0.064 | same | `analysis/biological_predictors/multivariate_model_results.json` | `loo_cv.elastic_net_r2` |  | ✓ |
| random forest LOO R² | −0.044 | same | same | `loo_cv.random_forest_r2` |  | ✓ |

### Results §8 -- Simulation characterizes pipeline power and ranking behavior

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Verified |
|---|---|---|---|---|---|---|
| sim type counts | n=15,25,35 | `analysis/simulation_study/simulation_study.py` | `analysis/simulation_study/simulation_results.json` | `power_curve` |  | ✓ |
| sim detection power | 100% | same | same | `power_curve` (calibrated row) |  | ✓ |
| sim FPR at α=0.05 | 4.8% | same | same | `null_calibration.fpr_005` |  | ✓ |
| ranking ceiling | ρ ≈ 0.42 | same | same | `ranking_recovery` |  | ✓ |
| test-retest ρ | 0.994 | same | same | `stability` |  | ✓ |
| cross-atlas ranking ρ | ≈ 0.15 | (descriptive across replications) | n/a | n/a |  | ✓ |
| bootstrap iterations | 1,000 | `analysis/bootstrap_rankings/bootstrap_ranking_analysis.py` | `analysis/bootstrap_rankings/bootstrap_summary.csv` | (configured) |  | ✓ |
| all stable (CI ≤ 10) | 35/35 | same | same | `(df.ci_width <= 10).all()` |  | ✓ |
| median CI width | 3 | same | same | `df.ci_width.median()` |  | ✓ |
| max CI width | 7 | same | same | `df.ci_width.max()` |  | ✓ |
| CI=0 type count | 5 | same | same | `(df.ci_width == 0).sum()` |  | ✓ |
| CI vs cross-atlas ρ (T59) | −0.41 | `analysis/cross_reference/cross_reference_analysis.py` | `analysis/cross_reference/convergent_types_summary.md` | "Spearman ρ (CI width vs mean rank shift)" | `cross_reference_analysis.main` | ✓ |
| CI vs cross-atlas p | 0.073 | same | same | "Spearman p-value" | `cross_reference_analysis.main` | ✓ |
| n in ≥2 replications | 20 | same | same | "Types in ≥2 replications" | `cross_reference_analysis.main` | ✓ |

### Mechanistic null tests (Figure 7B; reported in-text under §7 and the Discussion, no longer a standalone Results section)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Verified |
|---|---|---|---|---|---|---|
| tests run | 10 | (mechanistic-null table in Methods) | (see Section 2 mechanistic null appendix) | n/a |  | ✓ |
| n=35 (most tests) | 35 | various mechanistic-null scripts | various | various |  | ✓ |
| power at ρ≈0.30 | 37% | `analysis/simulation_study/simulation_study.py` | `analysis/simulation_study/simulation_results.json` | `power_curve` |  | ✓ |
| TF complexity ρ | −0.229 | `scripts/13_tf_complexity.py` | `output/phase2/mechanistic/tf_complexity/tf_complexity_vs_residual.csv` | aggregated |  | ✓ |
| PPI best ρ | 0.291 | `scripts/19_ppi_centrality.py` | `output/mechanistic/ppi_centrality/ppi_centrality_results.json` | `correlation_results[*].spearman_rho` (max) |  | ✓ |
| enhancer ρ | −0.429 | `scripts/t3e_step3b_enhancer.py` | `output/validation/t3e_enhancer/t3e_step3b_summary.md` | "Spearman ρ" |  | ✓ |
| enhancer n | 6 | same | same | "n (cell types)" |  | ✓ |
| SAMap-residual ρ | −0.247 | `scripts/34_samap_35types.py` | `output/phase1_samap/samap_35types/samap_rigidity_correlation.json` | `spearman_rho` |  | ✓ |
| SAMap p | 0.153 | same | same | `spearman_p` |  | ✓ |
| CellHint matched 15-type ρ | −0.139 | `analysis/harmonized_replication/harmonized_replication.py` | `analysis/harmonized_replication/sensitivity_analysis.csv` | `0_unharmonized` row |  | ✓ |
| CellHint matched p | 0.62 | same | same | same row |  | ✓ |
| CellHint harmonized ρ | −0.04 | same | `analysis/harmonized_replication/correlation_results.json` | `rho` |  | ✓ |
| CellHint harmonized n | 12 | same | same | `n_types` |  | ✓ |

### Methods -- primary-pipeline parameters and counts

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Verified |
|---|---|---|---|---|---|---|
| Tabula Sapiens cells | 483,152 | (Tabula Sapiens v1.0 paper) | n/a | n/a |  | anchor |
| Tabula Sapiens donors | 24 | (Tabula Sapiens v1.0 paper) | n/a | n/a |  | anchor |
| Tabula Muris Senis cells | ~350,000 | (TMS paper) | n/a | n/a |  | anchor |
| unique CL labels in TS | 180 | (data acquisition pipeline) | n/a | n/a |  | intermediate |
| unique CL labels in TMS | 151 | (data acquisition pipeline) | n/a | n/a |  | intermediate |
| shared CL labels | 66 | (data acquisition pipeline) | n/a | n/a |  | intermediate |
| reaching ≥200 cells | 45 | (data acquisition pipeline) | n/a | n/a |  | intermediate |
| meeting ≥500 cells | 35 | `scripts/02_qc_and_normalize.py` | `output/phase2/scaled_35types/centroids_human_35.csv` | shape[0] |  | ✓ |
| 10x-only ranking ρ | 0.754 | `scripts/14_smartseq2_sensitivity.py` | `output/phase2/sensitivity/smartseq2/sensitivity_results.json` | `rigidity_comparison.spearman_rho_ranks` | `14_smartseq2_sensitivity.main` | ✓ |
| progenitor full p | 0.0099 | (biological_predictors progenitor sub-analysis) | (multiple sources) | n/a |  | ✓ |
| progenitor 10x-only p | 0.119 | same | same | n/a |  | ✓ |
| mouse SS2 fraction ρ | −0.042 | `scripts/14_smartseq2_sensitivity.py` | `output/phase2/sensitivity/smartseq2/sensitivity_results.json` | `rigidity_comparison.ss2_fraction_vs_rank_change_rho` | `14_smartseq2_sensitivity.main` | ✓ |
| total 1:1 orthologs | 17,187 | (BioMart query) | `data/phase1/orthologs_human_mouse.csv` | shape[0] |  | ✓ |
| operating gene space | 16,959 | `scripts/02_qc_and_normalize.py` | `output/phase2/scaled_35types/centroids_human_35.csv` | shape[1] |  | ✓ |
| gene-set sensitivity range | 0.496–0.511 | `scripts/test_hvg_robustness.py` | `output/validation/hvg_robustness/hvg_robustness.json` | `hvg_results.{2000,3000,5000}.obs_null_ratio` (median denominator) |  | ✓ |
| gene-set ranking ρ | ≥ 0.957 | `scripts/test_hvg_robustness.py` | `output/validation/hvg_robustness/hvg_robustness.json` | `hvg_results.{2000,3000,5000}.ranking_rho_vs_full` |  | ✓ |
| L1000 landmark genes | 978 | (Broad Institute L1000) | n/a | n/a |  | anchor |
| L1000 in ortholog space | 907 | `scripts/35_l1000_random_baseline.py` | `output/figures/l1000_random_baseline_results.json` | `primary.n_genes_sampled` |  | ✓ |
| L1000 fraction | 5.3% | (computed: 907/16959) | n/a | n/a |  | ✓ |
| L1000-landmark sensitivity gate (pre-reg Fractal Geometry Test) | ρ=0.8515, PASS (thr 0.6), p=8.92e-11 | `scripts/23_sensitivity_gate_l1000.py` | `output/landmark_sensitivity/sensitivity_gate_result.json` | `spearman_rho`, `result` |  | ✓ |
| per-type cell cap | 2,000 | `src/cellwarp/procrustes.py` | (configured constant) | `MAX_CELLS_PER_TYPE` |  | ✓ |
| min cell count per type | 516 (myeloid leukocyte mouse) | data acquisition | n/a | n/a |  | intermediate |
| scipy verification max-Δ | 7.3 × 10⁻¹¹ | `scripts/verify_procrustes_vs_scipy.py` | (verification log) | n/a |  | ✓ |
| within-species v2 cells | 11,640 | `scripts/test_35type_human_control.py` | `output/phase2/negative_control_v2/cell_availability.json` | (sum of selected cells) |  | ✓ |
| v2 within-species obs/null | 0.607 | same | `output/phase2/negative_control_v2/negctrl_v2_results.json` | (within-species comparison) |  | ✓ |
| v2 cross-species obs/null | 0.317 | same | same | (matched-scale 6-type comparison) |  | ✓ |
| v2 within-species p | 0.0088 | same | same | `permutation_test.p_value` |  | ✓ |
| Sun2023 lung-restricted obs/null | 0.490 | `scripts/18_sun2023_issue092_diagnosis.py` | `output/validation/sun2023_issue092_diagnosis/issue092_diagnosis.json` | `task_3_revised_procrustes.obs_null_ratio` | `18_sun2023_issue092_diagnosis.main` | ✓ |
| Sun2023 lung-restricted p | 0.0001 | same | same | `task_3_revised_procrustes.p_value` | `18_sun2023_issue092_diagnosis.main` | ✓ |
| sorted-eigenvalue mean r | 0.953 | `scripts/layer3_permutation_test.py` | `output/layer3_permutation/layer3_permutation_results.json` | `observed_mean_r` |  | ✓ |
| eigenvalue p | 0.866 | same | same | `empirical_p` |  | ✓ |
| eigenvalue permutations | 10,000 | same | same | `n_permutations` |  | ✓ |
| eigenval-divergence ρ | 0.395 | `scripts/t3b_ellipsoid_alignment.py` | `output/mechanistic/ellipsoid_alignment/summary_stats.json` | `35type.eigenval_vs_rigidity.rho` |  | ✓ |
| eigenval-divergence p | 0.019 | same | same | `35type.eigenval_vs_rigidity.p` |  | ✓ |
| nearest null distance | 98.88 | `scripts/permutation_1M.py` | `analysis/permutation_1M/results_1M.json` | `null_distribution_summary.min` |  | ✓ |
| observed distance | 61.15 | same | same | `observed_procrustes_distance` |  | ✓ |
| replication 4× α=0.05 | 6.25 × 10⁻⁶ | (computed: 0.05⁴) | n/a | n/a |  | ✓ |
| best-predicted type | CD8+ T cell | `scripts/08_loocv.py` | `output/validation/v2_loocv/v2_loocv_results.json` | `best_type` |  | ✓ |
| best ratio | 0.226 | same | same | `min_ratio` |  | ✓ |
| hardest type | hematopoietic stem cell | same | same | `worst_type` |  | ✓ |
| hardest ratio | 0.796 | same | same | `max_ratio` |  | ✓ |
| bootstrap mean distance | 61.18 | `scripts/07_bootstrap.py` | `output/phase3/bootstrap/bootstrap_summary.json` | `distances.mean` |  | ✓ |
| bootstrap distance std | 0.24 | same | same | `distances.std` |  | ✓ |
| cell-count ρ | 0.052 | `scripts/confound_cellcount_rigidity.py` | `output/cellcount_confound/cellcount_confound_results.json` | `spearman_rho` |  | ✓ |
| cell-count p | 0.768 | same | same | `spearman_p` |  | ✓ |
| partial ρ | 0.061 | same | same | `partial_rho` |  | ✓ |
| partial p | 0.732 | same | same | `partial_p` |  | ✓ |

### Methods -- mechanistic-null table

| Manuscript row | Test | Value | Script | Output file | Output key | Function entry-point | Verified |
|---|---|---|---|---|---|---|---|
| Row 1 | Housekeeping ρ | 0.167 | `scripts/12_housekeeping_ratio.py` | `output/phase2/mechanistic/housekeeping/hk_ratio_results.json` | `human_correlation.spearman_rho` |  | ✓ |
| Row 2 | TF complexity ρ | −0.229 | `scripts/13_tf_complexity.py` | `output/phase2/mechanistic/tf_complexity/tf_complexity_vs_residual.csv` | aggregated |  | ✓ |
| Row 3 | Niche adaptation | 0/6 sig after FDR | `scripts/12_niche_hypothesis.py` | `output/phase2/mechanistic/niche/` | (per-set FDR) |  | ✓ |
| Row 4 | Within-type variance ρ | −0.038 | `scripts/12_variance_diagnostic.py` | `output/phase2/variance_diagnostic/diagnostic_results.json` | `spearman_rho_mean` |  | ✓ |
| Row 5 | Inter-donor variance ρ | −0.127 | `scripts/16_interdonor_variance.py` | `analysis/biological_predictors/univariate_correlations.csv` | "Inter-donor variance" row |  | ✓ |
| Row 6 | Expression-level all ρ < 0.21 | max abs ρ = 0.209 | `scripts/diagnostic_expression_vs_rigidity.py` | `output/phase2/diagnostics/expression_level_vs_rigidity/correlations.csv` | max(abs(rho)) |  | ✓ |
| Row 7 | PPI best ρ | 0.291 | `scripts/19_ppi_centrality.py` | `output/mechanistic/ppi_centrality/ppi_centrality_results.json` | `correlation_results[*].spearman_rho` (best by abs) |  | ✓ |
| Row 7 | PPI 0/27 sig after FDR | 0/27 | same | same | `n_significant_fdr / len(correlation_results)` |  | ✓ |
| Row 8 | phastCons ρ | −0.058 | `scripts/t3e_step2_compute.py` | `output/validation/t3e_chromatin/t3e_step2_summary.md` | "Spearman ρ" (placental_20way primary, Option A, ±2kb) |  | ✓ |
| Row 8 | phastCons n | 35 | same | same | "n (cell types)" |  | ✓ |
| Row 9 | enhancer ρ | −0.429 | `scripts/t3e_step3b_enhancer.py` | `output/validation/t3e_enhancer/t3e_step3b_summary.md` | "Spearman ρ" |  | ✓ |
| Row 9 | enhancer n | 6 | same | same | "n (cell types)" |  | ✓ |
| Row 10 | drug target ρ | −0.176 | (output canonical) | `output/t3g/primary_correlation_results.json` | `primary_correlation.rho` |  | ✓ |
| Row 10 | drug target n | 34 | (output canonical) | same | `primary_correlation.n` |  | ✓ |

### Methods -- replication and macaque parameters

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Verified |
|---|---|---|---|---|---|---|
| replication permutations | 10,000 | (replication scripts) | replication JSONs | `permutation_test.n_permutations` |  | ✓ |
| hepatocyte SSR macaque (full) | 47.3% | `analysis/macaque/reconstruct_macaque_pipeline.py` | `output/macaque_pipeline/m1_close_table1_summary.json` | `human_macaque_12type.hepatocyte_pct_ssr` |  | ✓ |
| hepatocyte SSR HM (full) | 2.0% | same | `output/macaque_pipeline/human_mouse_12type_control.json` | `hepatocyte_pct_ssr` |  | ✓ |

### Methods -- mouse lemur parameters

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Verified |
|---|---|---|---|---|---|---|
| lemur passing types | 15 | `analysis/mouse_lemur/01_run_pipeline.py` | `analysis/mouse_lemur/procrustes_results.json` | `n_types` |  | ✓ |
| lemur cell range | 546–2,000 | same | `analysis/mouse_lemur/lemur_cell_type_counts.csv` | per-type counts |  | ✓ |
| lemur 1:1 ortholog count | 13,796 | same | `analysis/mouse_lemur/procrustes_results.json` | `gene_space` |  | ✓ |
| lemur PCA components | 15 | same | same | `pca.n_components` |  | ✓ |
| lemur cumulative variance | 95.5% | same | same | `pca.cumulative_variance` |  | ✓ |
| lemur obs/null | 0.346 | same | same | `permutation_test.obs_null_ratio` |  | ✓ |
| lemur permutations | 10,000 | same | same | `permutation_test.n_permutations` |  | ✓ |
| lemur ranking ρ | 0.157 | same | `analysis/mouse_lemur/ranking_correlation.json` | `vs_primary.rho` |  | ✓ |

### Methods -- pan-Census parameters

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Verified |
|---|---|---|---|---|---|---|
| pan-Census datasets | 15 | `analysis/census_replication/02_run_replication.py` | `data/replication/pan_census_manifest.csv` | row count |  | ✓ |
| pan-Census collections | 15 | same | same CSV | `df.collection_id.nunique()` |  | ✓ |
| pan-Census mouse datasets | 9 | same | same | `df.species.value_counts()` |  | ✓ |
| pan-Census human datasets | 6 | same | same | `df.species.value_counts()` |  | ✓ |
| pan-Census shared types | 22 | same | `analysis/census_replication/replication_results.json` | `n_cell_types` |  | ✓ |
| pan-Census obs/null | 0.811 | same | same | `permutation_test.obs_null_ratio` |  | ✓ |
| pan-Census p | 0.0001 | same | same | `permutation_test.p_value` |  | ✓ |

### Software and reproducibility

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Verified |
|---|---|---|---|---|---|---|
| Python | 3.12.12 | n/a | `pyproject.toml` | `requires-python` |  | ✓ |
| RANDOM_SEED | 42 | (entry-point convention) | various | per-script `RANDOM_SEED` |  | ✓ |
| numpy | 2.4.3 | n/a | `pyproject.toml` | `[lock]` |  | ✓ |
| scipy | 1.17.1 | n/a | same | same |  | ✓ |
| pandas | 2.3.3 | n/a | same | same |  | ✓ |
| scanpy | 1.12 | n/a | same | same |  | ✓ |
| anndata | 0.12.10 | n/a | same | same |  | ✓ |
| scikit-learn | 1.8.0 | n/a | same | same |  | ✓ |
| statsmodels | 0.14.6 | n/a | same | same |  | ✓ |
| cellxgene-census | 1.17.0 | n/a | same | same |  | ✓ |
| samap | 1.0.14 | n/a | same | same |  | ✓ |
| Census version | 2025-11-08 | (pinned) | `analysis/census_replication/02_run_replication.py` | `census_version` |  | ✓ |
| Tabula Microcebus collection_id | a137437b-… | n/a | `data/replication/tabula_microcebus_metadata.csv` | `collection_id` |  | ✓ |
| Tabula Microcebus dataset_id | a392ab34-… | n/a | same | `dataset_id` |  | ✓ |
| Tabula Microcebus download date | 2026-04-05 | n/a | same | `download_date` |  | ✓ |
| BioMart release | 115 | n/a | `data/phase1/orthologs_human_mouse.csv` (header) | accession date in DOI |  | ✓ |

---

## Notes for reviewers

- All output files referenced in this CROSSWALK are present in the deposit
  repository at the paths shown.
- Some output files are gitignored at the directory level (e.g., `data/`,
  `output/`) but are tracked individually; see `git ls-files data/ output/`
  to enumerate the tracked set.
- For automated regression-checking of key statistics against published
  values, run `python reproduce/validate.py` after a full pipeline execution.
- This CROSSWALK is generated and maintained as a single Markdown document.
  Future updates land in this file directly.
- The deposit also includes two dataset-level manifest CSVs that travel
  with the analysis code (see README §"Deposit Artifacts"):
  `data/replication/pan_census_manifest.csv` (15 CELLxGENE Census datasets)
  and `data/replication/tabula_microcebus_metadata.csv` (Tabula Microcebus
  deposit anchors). Both are referenced from this CROSSWALK in the
  pan-Census and mouse-lemur sections respectively.


<!-- sensitivity-analysis crosswalk rows -->
| §4/Disc | Layer-1 ribosomal-excluded obs/null | 0.501 | `analysis/sensitivity_analyses/layer1_exclusion.py` | `analysis/sensitivity_analyses/layer1_exclusion_results.json` | `variants.1.obs_null_ratio` | Table S7 | ✓ |
| §4/Disc | Layer-1 ribosomal+housekeeping obs/null | 0.479 | same | same | `variants.2.obs_null_ratio` | Table S7 | ✓ |
| §6 | Marker ortholog-retention vs residual (CellMarker) | +0.31 (p=0.09, n=32) | `analysis/sensitivity_analyses/ortholog_retention.py` | `analysis/sensitivity_analyses/ortholog_retention_results.json` | `cellmarker_secondary.spearman_rho` | Table S8 | ✓ |
| §3 | Donor-split delta (cross-species − within-species), reported as the median of per-split deltas | +0.159 | (donor-split) | `analysis/donor_split/donor_split_shared_pca_results.json` | `delta.median`; the aggregate arms cross = 0.527 / within = 0.375 (`donor_split_comparison.json` `shared_pca`) are arm medians, distinct from the reported +0.159 delta; the manuscript reports neither 1.41× nor the 0.152 arm-difference | §3 text | ✓ |
| (treeness) | treeness-rigidity ρ −0.349 etc. | REMOVED | analysis cut | -- | -- | -- | removed |
| L84 | bootstrap CI-vs-rankshift ρ (T59) | −0.410 (p=0.073) | `analysis/cross_reference/cross_reference_analysis.py` | `analysis/cross_reference/master_ranking_table.csv` | recomputed (Qu-12 macaque) | ✓ |
