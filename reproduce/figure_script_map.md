# Figure and Table to Script Mapping

Every display item in the CellWarp paper, mapped to the script that generates it.
Figure and panel numbers follow the current manuscript after the S5-CW1b rebuild.

## Main Figures

| Display Item | Generating Script(s) | Depends On |
|---|---|---|
| Figure 1A (pipeline schematic) | scripts/generate_phase2_figures.py | None |
| Figure 1B (1M null distribution) | scripts/generate_phase1_figures.py | scripts/permutation_1M.py |
| Figure 1C (lineage-stratified null) | scripts/generate_phase1_figures.py | scripts/test_lineage_stratified_permutation.py |
| Figure 1D (bootstrap stability) | scripts/generate_phase3_figures.py | scripts/07_bootstrap.py |
| Figure 1E (LOOCV bar chart) | scripts/generate_phase3_figures.py | scripts/08_loocv.py |
| Figure 2A (per-gene conservation C distribution; Hartigan dip) | analysis/conserved_contribution/make_figure7.py | analysis/conserved_contribution/run_gate.py, make_table_s11.py |
| Figure 2B (C vs expression and Tau specificity) | analysis/conserved_contribution/make_figure7.py | analysis/conserved_contribution/run_gate.py, run_robustness.py |
| Figure 2C (master-TF conservation enrichment vs matched backgrounds) | analysis/conserved_contribution/make_figure7.py | analysis/conserved_contribution/run_gate.py, highN_tf_pvalues.py |
| Figure 2D (per-gene C donor-split reproducibility) | analysis/conserved_contribution/make_figure7.py | analysis/conserved_contribution/donor_stability run |
| Figure 3A (Sun2023 null) | scripts/generate_phase3_figures.py | scripts/16_sun2023_replication.py |
| Figure 3B (PanSci null) | scripts/generate_phase3_figures.py | scripts/pansci_replication.py |
| Figure 3C (CellHint null) | scripts/generate_phase3_figures.py | scripts/33_cellhint_replication.py |
| Figure 3D (replication summary) | scripts/generate_phase1_figures.py | All replication scripts |
| Figure 3E (donor-split / human control) | scripts/generate_phase1_figures.py | scripts/test_35type_human_control.py |
| Figure 4A (ellipsoid heatmap) | scripts/generate_phase3_figures.py | scripts/t3b_ellipsoid_alignment.py |
| Figure 4B (pre vs post-rotation) | scripts/generate_phase3_figures.py | scripts/t3b_ellipsoid_alignment.py |
| Figure 4C (layer null distributions) | scripts/generate_phase3_figures.py | scripts/layer3_permutation_test.py |
| Figure 4D (layer scatter) | scripts/generate_phase3_figures.py | scripts/t3b_ellipsoid_alignment.py |
| Figure 5A (three-species Procrustes) | scripts/generate_phase1_figures.py | Macaque pipeline (pre-computed) |
| Figure 5B (no-immune sensitivity) | scripts/generate_phase1_figures.py | Macaque pipeline (pre-computed) |
| Figure 5C (hepatocyte ranking) | scripts/generate_phase1_figures.py | Macaque pipeline + primary |
| Figure 6A (divergence ranking) | scripts/generate_phase1_figures.py | scripts/08_scaled_procrustes.py |
| Figure 6B (cell count confound) | scripts/generate_phase1_figures.py | scripts/confound_cellcount_rigidity.py |
| Figure 7A (L1000 baseline) | scripts/generate_phase3_figures.py | scripts/35_l1000_random_baseline.py |
| Figure 7B (mechanistic nulls forest) | scripts/generate_phase2_figures.py | Mechanistic null scripts (see below) |
| Table 1 (unified statistical tests) | Compiled in manuscript from all scripts | All analysis scripts |

## Supplementary Figures

| Display Item | Generating Script(s) | Depends On |
|---|---|---|
| Figure S1A-B (independent PCA) | scripts/generate_phase1_figures.py | analysis/independent_pca_sensitivity/run_independent_pca.py |
| Figure S1C-F (simulation study) | analysis/simulation_study/simulation_figures.py | analysis/simulation_study/simulation_study.py |
| Figure S2A-B (PCA k-sensitivity) | scripts/generate_phase3_figures.py | scripts/17_pca_sensitivity.py, scripts/18_pca_sensitivity_v2.py |
| Figure S2C-D (Smart-seq2 protocol) | scripts/generate_phase3_figures.py | scripts/14_smartseq2_sensitivity.py |
| Figure S2E (expanded negatives) | analysis/expanded_negative_controls/negative_control_figure.py | analysis/expanded_negative_controls/expanded_negative_controls.py |
| Figure S2F (replication inventory) | scripts/56_add_figs2_panel_f.py | All replication outputs |
| Figure S3A-B (bootstrap rankings) | scripts/composite_figS3.py | analysis/bootstrap_rankings/bootstrap_ranking_analysis.py |
| Figure S4A-B (CellHint investigation) | scripts/generate_phase3_figures.py | analysis/cellhint_investigation/investigate_rank_reversal.py |
| Figure S5A (SAMap heatmap) | scripts/generate_phase3_figures.py | scripts/34_samap_35types.py |
| Figure S6A-B (CellMarker enrichment) | scripts/generate_phase1_figures.py | scripts/cellmarker_35type_rerun.py |
| Figure S7 (matched-scale 6-type negative control) | scripts/test_35type_human_control.py | output/phase2/negative_control_v2/negctrl_v2_results.json |
| Figure S8 (marker-similarity-stratified null) | analysis/sensitivity_analyses/markernull.py | (primary centroids; species-averaged gene-space) |

## Supplementary Tables

| Display Item | Generating Script(s) | Depends On |
|---|---|---|
| Table S1 (biological predictors + cross-atlas) | scripts/create_table_S1.py | analysis/biological_predictors/biological_predictors.py, analysis/ranking_replication/ranking_replication_analysis.py |
| Table S2 (simulation + bootstrap CIs) | scripts/create_table_S2.py | analysis/simulation_study/simulation_study.py, analysis/bootstrap_rankings/bootstrap_ranking_analysis.py |
| Table S3 (CellHint rank reversal) | analysis/cellhint_investigation/investigate_rank_reversal.py | scripts/33_cellhint_replication.py |
| Table S4 (CellHint harmonization) | analysis/harmonized_replication/harmonized_replication.py | scripts/33_cellhint_replication.py |
| Table S5 (35-type matching) | scripts/08_cell_type_inventory.py | scripts/02_qc_and_normalize.py |
| Table S6 (CPC1 driver genes) | scripts/generate_table_S6.py | scripts/t3b_ellipsoid_alignment.py |
| Table S7 (Layer-1 ribosomal/housekeeping exclusion) | analysis/sensitivity_analyses/layer1_exclusion.py | scripts/08_scaled_procrustes.py |
| Table S8 (marker 1:1-ortholog retention vs residual) | analysis/sensitivity_analyses/ortholog_retention.py | scripts/cellmarker_35type_rerun.py, data/phase1/orthologs_human_mouse.csv |
| Table S9 (per-gene standardization: Layer-1 + CPC1 Scheme A/B) | analysis/sensitivity_analyses/genestd_standardization.py | scripts/08_scaled_procrustes.py, scripts/t3b_ellipsoid_alignment.py |
| Table S10 (marker-similarity-stratified null) | analysis/sensitivity_analyses/markernull.py | (primary centroids; species-averaged gene-space) |
| Table S11 (per-gene cross-species conservation score C) | analysis/conserved_contribution/make_table_s11.py | analysis/conserved_contribution/run_gate.py |

## Mechanistic Null Scripts (Figure 7B)

| Test | Script |
|---|---|
| Housekeeping ratio | scripts/12_housekeeping_ratio.py |
| TF complexity | scripts/13_tf_complexity.py |
| Niche hypothesis | scripts/12_niche_hypothesis.py |
| Variance diagnostic | scripts/12_variance_diagnostic.py |
| Interdonor variance | scripts/16_interdonor_variance.py |
| PPI centrality | scripts/19_ppi_centrality.py |
| Chromatin conservation | scripts/t3e_step2_compute.py |
| Enhancer conservation | scripts/t3e_step3b_enhancer.py |
| Expression vs rigidity | scripts/diagnostic_expression_vs_rigidity.py |
| Ribosomal confound | scripts/16_ribosomal_confound_test.py |
