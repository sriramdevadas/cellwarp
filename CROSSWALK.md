# CROSSWALK: Manuscript → Code

This document maps every Methods subsection and every figure and table in the
current CellWarp manuscript (PLOS ONE), together with the numerical values behind
its reportable claims, to the code that produces them and the output file that
contains each value.

This is a reproducibility deposit artifact: a reviewer can use this document to
walk from any specific manuscript claim to the script that generates it and the
persisted output that contains the value.

Three sections:

1. **Methods → code** -- one entry per current Methods subsection (9 total),
   with the primary generating script(s) and output.
2. **Figures and tables → code** -- one entry per current main- and
   supplementary-figure panel and supplementary table, with generating script
   and dependencies.
3. **Numerical claims → code** -- the numerical values behind the reportable
   claims in the manuscript and Supporting Information, with the generating
   script, output file, key or column, and a Status showing how far each is
   established.

Quick-look statistics about this CROSSWALK (re-derived from the current sections):
- Methods subsections covered: 9
- Figure / table display items covered: 5 main figures + 5 supporting figures + 11 supplementary tables (35 panel/table rows in Section 2: 15 main-figure panel rows, 9 supplementary-figure panel rows, 11 tables)
- Numerical claims indexed: 291
- Status breakdown: validate.py 39, mapped 232, computed 9, anchor 6, intermediate 5

Verification. `reproduce/validate.py` programmatically asserts its checks
against their persisted output files; most correspond to claim rows in Section
3 (Status = validate.py), and the rest assert values outside Section 3's
current scope -- regression guards retained for the cut SAMap / CellMarker /
L1000 analyses, plus a few with no Section-3 counterpart. Every
remaining Section-3 value carries a Status recording how far it is established:
mapped to its generating script and output file with the path resolved, computed
arithmetically from other cited values, anchored to a source publication, or an
unpersisted intermediate. None of these is independently value-re-verified here.

All output files referenced in this CROSSWALK are present in the deposit
repository at the paths shown.

For automated regression-checking of the load-bearing statistics, run
`python reproduce/validate.py` after a full pipeline execution.

The "Function entry-point(s)" column in Section 1 lists procedural functions
corresponding to manuscript-described steps; helper functions are omitted. The
same column appears in Section 3 only where the value's computing function differs
from the entry in Section 1 for that subsection (most rows therefore leave the
column blank).

---

## 1. Methods → code

For each Methods subsection of the current manuscript
(`docs/submission/plosone/manuscript_combined.txt`, in manuscript order), the
primary generating script(s) and the output where results are persisted. One row
per current subsection.

Controls the manuscript describes inside a Methods subsection (independent-PCA and
Mantel artifact-ruling controls, bootstrap, leave-one-out cross-validation,
per-gene-standardized CPC1) are folded into that subsection's script list. Content
that survives only in Supporting Information rather than as a standalone Methods
subsection (macaque extension, cross-atlas replications, donor-split and expanded
negative controls, the ten mechanistic nulls, bootstrap-ranking stability,
biological predictors) is not given its own Methods row here; it is indexed in
Section 3. The former SAMap, CellMarker, and L1000 analyses are cut from the
current paper and are omitted.

| Manuscript subsection (line) | Primary script(s) | Function entry-point(s) | Output |
|---|---|---|---|
| Ethics statement (L129) | n/a (computational reanalysis of public data; no code) | n/a | n/a |
| Data and cell-type matching (L135) | `scripts/02_qc_and_normalize.py` (Census fetch, QC, CP10K normalization, centroids); `scripts/08_cell_type_inventory.py` (35-type Cell-Ontology matching, S5 Table); Ensembl BioMart query (1:1 orthologs, archived CSV); `analysis/mouse_lemur/01_run_pipeline.py` (Tabula Microcebus data) | `02_qc_and_normalize.main`, `procrustes.compute_centroids`; `08_cell_type_inventory.main` | `output/phase2/scaled_35types/centroids_*.csv`; `data/phase1/orthologs_human_mouse.csv`; `output/phase2/cell_type_inventory.csv` |
| The Procrustes framework (L143) | `src/cellwarp/procrustes.py` + `scripts/08_scaled_procrustes.py` (joint PCA, superimposition, per-type residuals); `scripts/permutation_1M.py`, `scripts/test_lineage_stratified_permutation.py` (permutation + lineage-stratified nulls); `scripts/07_bootstrap.py`, `scripts/08_loocv.py` (bootstrap, LOOCV robustness); `analysis/independent_pca_sensitivity/run_independent_pca.py`, `analysis/mantel_test/` (artifact-ruling controls, S1 Text §1) | `procrustes.procrustes_align`, `procrustes.pca_reduce_centroids`, `procrustes.permutation_test`, `08_scaled_procrustes.stage4_procrustes` | `output/phase2/scaled_35types/procrustes_results_35.json`; `analysis/permutation_1M/results_1M.json`; `output/validation/lineage_stratified/lineage_stratified_results.json`; `output/phase3/bootstrap/bootstrap_summary.json`, `output/validation/v2_loocv/v2_loocv_results.json`; `analysis/independent_pca_sensitivity/independent_pca_results.json`; `analysis/mantel_test/mantel_results.json` |
| The two-layer decomposition (L151) | `scripts/t3b_ellipsoid_alignment.py` (Layer-2 Krzanowski S, CPC1 drivers) + `scripts/layer3_permutation_test.py` (eigenvalue-conservation null); `analysis/sensitivity_analyses/genestd_standardization.py` (per-gene-standardized CPC1, S1 Text §3) | `t3b.compute_covariance_eigen`, `t3b.compute_alignment_scores`, `t3b.perm_test_label_shuffle`, `t3b.eigenvalue_conservation` | `output/mechanistic/ellipsoid_alignment/summary_stats.json`, `permutation_results.json`; `output/layer3_permutation/layer3_permutation_results.json`; `analysis/sensitivity_analyses/genestd_results.json` |
| Primate replication (basal ganglia) (L157) | `docs/submission/plosone/figures/build_fig2c_bg.py` reproduces the Fig 2C panel here from the vendored `docs/submission/plosone/figures/bg_results/`; the upstream two-layer statistics come from the self-contained basal-ganglia deposit declared in the manuscript's Data and code availability statement (L181), archived at Zenodo, its DOI not yet minted | n/a (vendored result JSONs; the upstream BG deposit -- length-residualization and low-N robustness -- is an external repository) | `docs/submission/plosone/figures/bg_results/layer2_results_{pair}.json`, `layer2_cpc1_drivers_{pair}_W2_schemeB.csv`; `docs/submission/plosone/figures/Fig2C_bg_replication.{tiff,pdf,png}` |
| Simulation (L165) | `analysis/simulation_study/simulation_study.py` | `generate_centroids`, `calibrate`, `run_power_curve`, `run_ranking_recovery`, `run_null_calibration` | `analysis/simulation_study/simulation_results.json` |
| Conserved-contribution and identity-gene analysis (L169) | `analysis/conserved_contribution/run_gate.py` (conservation score C, geometry attribution, master-TF enrichment gate); `run_robustness.py` (Hartigan dip, Tau); `highN_tf_pvalues.py` (matched-null empirical p-values); `make_table_s11.py` (Table S11); `make_figure7.py` (Fig 5) | `run_gate.main`, `gate_lib.build_gene_table`, `run_robustness.main` | `analysis/conserved_contribution/gate_results.json`, `robustness_results.json`, `highN_tf_pvalues.json`, `gene_conservation_core.csv`, `donor_stability/donor_stability_results.json` |
| Data and code availability (L179) | n/a (deposit metadata + version-pinned manifest) | n/a | `pyproject.toml`, `requirements.txt`; `data/replication/pan_census_manifest.csv`, `data/replication/tabula_microcebus_metadata.csv`; `docs/supplementary_materials/table_S12_software_environment.csv` |
| Use of generative AI (L183) | n/a (no code) | n/a | n/a |

---

## 2. Figures and tables → code

Every display item in the current CellWarp paper, mapped to the script that
generates it and its dependencies. Figure and panel numbers follow the current
manuscript (`docs/submission/plosone/manuscript_combined.txt`).
`reproduce/figure_script_map.md` is the authority for how each display item is
produced; this section covers the same ground from the manuscript's side and does
not undertake to track it line for line. Old-to-new renumbering (this paper
descends from a 7-figure PCOMPBIOL draft) is recorded in
`docs/submission/plosone/NUMBER_DIFF.md` and is out of scope here.

### Main figures

The five main composites are assembled by
`docs/submission/plosone/figures/build_main_figures.py` (PDF + PNG) into
`docs/submission/plosone/figures/`: Fig 1 → `Fig1_configuration_conserved`,
Fig 2 → `Fig2_two_layers_bg`, Fig 3 → `Fig3_configuration_robust`,
Fig 4 → `Fig4_pertype_not_resolvable`, Fig 5 → `Fig5_conserved_identity_genes`.
The "Generating Script(s)" column names the panel or data producer; the "Output
file" column names the backing data artifact, as elsewhere in this section.

| Display Item | Generating Script(s) | Depends On | Output file |
|---|---|---|---|
| Fig 1A (pipeline schematic) | `scripts/generate_phase2_figures.py` (panel) | None | `figures/panels/fig1a_pipeline_schematic.png` (schematic; no data file) |
| Fig 1B (1M null distribution, obs/null 0.52) | `scripts/generate_phase1_figures.py` (panel) | `scripts/permutation_1M.py` | `analysis/permutation_1M/results_1M.json`, `analysis/permutation_1M/null_distribution_1M.npy` |
| Fig 1C (lineage-stratified null) | `scripts/generate_phase1_figures.py` (panel) | `scripts/test_lineage_stratified_permutation.py` | `output/validation/lineage_stratified/lineage_stratified_results.json` |
| Fig 1D (human-mouse-lemur null, n=15) | built fresh in `docs/submission/plosone/figures/build_main_figures.py` | `analysis/mouse_lemur/` pipeline | `analysis/mouse_lemur/null_distribution.npy`, `analysis/mouse_lemur/procrustes_results.json` (obs value hard-coded in assembler; see Known gaps) |
| Fig 2A (per-type Krzanowski S heatmap, k=1,3,5) | `scripts/generate_phase3_figures.py` (panel) | `scripts/t3b_ellipsoid_alignment.py` | `output/mechanistic/ellipsoid_alignment/35type_alignment_scores.csv` |
| Fig 2B (aggregate S pre vs post rotation) | `scripts/generate_phase3_figures.py` (panel) | `scripts/t3b_ellipsoid_alignment.py` | `output/mechanistic/ellipsoid_alignment/permutation_results.json`, `35type_alignment_scores.csv` |
| Fig 2C (basal-ganglia three-pair Layer-2 replication) | `docs/submission/plosone/figures/build_fig2c_bg.py` | vendored basal-ganglia results | `docs/submission/plosone/figures/bg_results/layer2_results_{pair}.json`, `layer2_cpc1_drivers_{pair}_W2_schemeB.csv` (standalone deposit `Fig2C_bg_replication.tiff`) |
| Fig 3 (replication summary bar chart) | `scripts/generate_phase1_figures.py` (panel) | All replication scripts | All replication JSONs (Sun2023, PanSci, CellHint, Andrews, MCA × HCA, pan-Census) |
| Fig 4A (within-atlas precision, bootstrap CI forest) | built fresh in `docs/submission/plosone/figures/build_main_figures.py` | `analysis/bootstrap_rankings/bootstrap_ranking_analysis.py` | `analysis/bootstrap_rankings/bootstrap_summary.csv` |
| Fig 4B (within-atlas precision vs cross-atlas rank shift, rho = -0.41) | built fresh in `docs/submission/plosone/figures/build_main_figures.py` | `analysis/cross_reference/cross_reference_analysis.py` | `analysis/cross_reference/master_ranking_table.csv` |
| Fig 4C (simulation recovery ceiling, rho ~ 0.42) | built fresh in `docs/submission/plosone/figures/build_main_figures.py` | `analysis/simulation_study/simulation_study.py` | `analysis/simulation_study/simulation_results.json` |
| Fig 5A (per-gene conservation C distribution; Hartigan dip) | `analysis/conserved_contribution/make_figure7.py` | `analysis/conserved_contribution/run_gate.py`, `make_table_s11.py` | `analysis/conserved_contribution/gate_results.json`, `gene_conservation_core.csv` |
| Fig 5B (C vs expression and Tau specificity) | `analysis/conserved_contribution/make_figure7.py` | `analysis/conserved_contribution/run_gate.py`, `run_robustness.py` | `analysis/conserved_contribution/gate_results.json`, `robustness_results.json` |
| Fig 5C (master-TF conservation enrichment vs matched backgrounds) | `analysis/conserved_contribution/make_figure7.py` | `analysis/conserved_contribution/run_gate.py`, `highN_tf_pvalues.py` | `analysis/conserved_contribution/gate_results.json`, `highN_tf_pvalues.json` |
| Fig 5D (per-gene C donor-split reproducibility) | `analysis/conserved_contribution/make_figure7.py` | (donor-stability run) | `analysis/conserved_contribution/donor_stability/donor_stability_results.json` |

### Supplementary figures

Deposited in `figures/submission/supplementary/`.

| Display Item | Generating Script(s) | Depends On | Output file |
|---|---|---|---|
| S1 Fig A-B (independent PCA) | `scripts/build_submission_figures.py` | `analysis/independent_pca_sensitivity/run_independent_pca.py` | `analysis/independent_pca_sensitivity/independent_pca_results.json` |
| S1 Fig C-F (simulation study) | `scripts/build_submission_figures.py` (composites `figures/supplementary/figS7_simulation_study_polished.pdf`) | `analysis/simulation_study/simulation_figures.py`, `simulation_study.py` | `analysis/simulation_study/simulation_results.json` |
| S2 Fig A-B (PCA k-sensitivity) | `scripts/build_submission_figures.py` | `scripts/17_pca_sensitivity.py`, `scripts/18_pca_sensitivity_v2.py` | `output/validation/pca_sensitivity/pca_sensitivity_results.json` |
| S2 Fig C-D (Smart-seq2 protocol) | `scripts/build_submission_figures.py` | `scripts/14_smartseq2_sensitivity.py` | `output/phase2/sensitivity/smartseq2/sensitivity_results.json` |
| S2 Fig E (expanded negatives) | `scripts/build_submission_figures.py` | `analysis/expanded_negative_controls/expanded_negative_controls.py` | `analysis/within_species_matched/matched_control_results.json` |
| S2 Fig F (replication inventory) | `scripts/56_add_figs2_panel_f.py` | All replication outputs | All replication JSONs (Sun2023, PanSci, CellHint, Andrews, MCA × HCA, pan-Census, HCA × Tabula) |
| S3 Fig A-B (bootstrap rankings) | `scripts/composite_figS3.py` (invoked by `build_submission_figures.py`) | `analysis/bootstrap_rankings/bootstrap_ranking_analysis.py` | `analysis/bootstrap_rankings/bootstrap_summary.csv`, `bootstrap_rankings_raw.csv` |
| S4 Fig (matched-scale 6-type negative control) | `scripts/49_build_figS7_matched_scale.py` (producer retains the old `figS7` stem) | `scripts/test_35type_human_control.py` | `output/phase2/negative_control_v2/negctrl_v2_results.json` |
| S5 Fig (marker-similarity-stratified null) | `analysis/sensitivity_analyses/markernull.py` (producer writes `figure_S8_markernull.{pdf,png}` to `docs/supplementary_materials/`, not the deposited `figS5_markernull.pdf`; see Known gaps) | (primary centroids; species-averaged gene-space) | `analysis/sensitivity_analyses/markernull_results.json` |

### Supplementary tables

The current paper has no S8 Table (the former ortholog-retention table was cut).
For S1-S5 the deposited canonical file is written or edited in place by the
materializer `scripts/46_synthesis_pass_supplementary_table_edits.py`; both the
content producer and the canonical file are listed.

| Display Item | Generating Script(s) | Depends On | Output file |
|---|---|---|---|
| Table S1 (biological predictors + cross-atlas + three-species) | `analysis/biological_predictors/biological_predictors.py`, `analysis/ranking_replication/ranking_replication_analysis.py`; canonical via `scripts/46_synthesis_pass_supplementary_table_edits.py` (`edit_table_s1`) | `scripts/create_table_S1.py` | `analysis/biological_predictors/univariate_correlations.csv`, `analysis/cross_reference/master_ranking_table.csv` → `docs/supplementary_materials/table_S1.xlsx` |
| Table S2 (simulation + bootstrap CIs) | `analysis/simulation_study/simulation_study.py`, `analysis/bootstrap_rankings/bootstrap_ranking_analysis.py`; canonical via `46_synthesis…` (`edit_table_s2`) | `scripts/create_table_S2.py` | `analysis/simulation_study/simulation_results.json`, `analysis/bootstrap_rankings/bootstrap_summary.csv` → `docs/supplementary_materials/table_S2.xlsx` |
| Table S3 (CellHint rank reversal) | `analysis/cellhint_investigation/investigate_rank_reversal.py`; canonical via `46_synthesis…` (`edit_table_s3`) | `scripts/33_cellhint_replication.py` | `analysis/cellhint_investigation/` (rank-reversal artifacts) → `docs/supplementary_materials/table_S3.csv` |
| Table S4 (CellHint harmonization) | `analysis/harmonized_replication/harmonized_replication.py`; canonical via `46_synthesis…` (`edit_table_s4`) | `scripts/33_cellhint_replication.py` | `analysis/harmonized_replication/sensitivity_analysis.csv`, `correlation_results.json` → `docs/supplementary_materials/table_S4.csv` |
| Table S5 (35-type matching) | `scripts/08_cell_type_inventory.py`; canonical via `46_synthesis…` (`edit_table_s5`) | `scripts/02_qc_and_normalize.py` | `output/phase2/cell_type_inventory.csv` → `docs/supplementary_materials/table_S5.csv` |
| Table S6 (CPC1 driver genes) | `scripts/generate_table_S6.py` | `scripts/t3b_ellipsoid_alignment.py` | `output/mechanistic/ellipsoid_alignment/cpc_gene_table.csv`, `cpc_genes.json` → `docs/supplementary_materials/Table_S6_CPC1_driver_genes.xlsx` |
| Table S7 (Layer-1 ribosomal/housekeeping exclusion) | `analysis/sensitivity_analyses/layer1_exclusion.py` | `scripts/08_scaled_procrustes.py` | `analysis/sensitivity_analyses/layer1_exclusion_results.json` + per-variant `layer1_exclusion_ranking_*.csv`; canonical `docs/supplementary_materials/table_S7_layer1_housekeeping_exclusion.csv` has NO scripted writer (see Known gaps) |
| Table S9 (per-gene standardization: Layer-1 + CPC1 Scheme A/B) | `analysis/sensitivity_analyses/genestd_standardization.py` | `scripts/08_scaled_procrustes.py`, `scripts/t3b_ellipsoid_alignment.py` | `analysis/sensitivity_analyses/genestd_results.json` → `docs/supplementary_materials/table_S9_genestd_standardization.csv`, `table_S9_schemeB_CPC1_markers.csv` |
| Table S10 (marker-similarity-stratified null) | `analysis/sensitivity_analyses/markernull.py` | (primary centroids; species-averaged gene-space) | `analysis/sensitivity_analyses/markernull_results.json` → `docs/supplementary_materials/table_S10_markernull.csv` |
| Table S11 (per-gene cross-species conservation score C) | `analysis/conserved_contribution/make_table_s11.py` | `analysis/conserved_contribution/run_gate.py` | `analysis/conserved_contribution/gene_conservation_core.csv` → `docs/supplementary_materials/table_S11_gene_conservation.csv` |
| Table S12 (software environment) | NO IN-REPO PRODUCER; hand-authored, source `requirements.txt` (see Known gaps) | (none) | `docs/supplementary_materials/table_S12_software_environment.csv` |

### Known gaps

The gaps in how figures and tables are produced are listed once, in the "Known
gaps" section of `reproduce/figure_script_map.md`, which is the authority for that
ground. The "(see Known gaps)" notes in the tables above refer to it.

They used to be restated here as well. The restatement has been removed rather
than resynchronised: it had already drifted from the original by the time anyone
noticed, which is what a second copy eventually does. One list is easier to keep
true than two that promise to agree.

---

## 3. Numerical claims → code

The numerical values behind the reportable claims in the manuscript and Supporting
Information, each mapped to its generating script and persisted output, and grouped
by the location it now occupies in the
current five-figure paper. Subsection headers name that current location; `####`
sub-headers mark rows that spill to a different location than their block's primary
destination. Old-to-new figure/section renumbering history is in
`docs/submission/plosone/NUMBER_DIFF.md`.

Status legend:
- **validate.py (<check name>)** -- value asserted by `reproduce/validate.py` (the named check).
- **mapped** -- value present in the deposit and reproducible from the cited output; not one of validate.py's automated checks.
- **computed** -- derived arithmetically from other rows or cited values; not separately persisted.
- **anchor** -- anchored to a source publication, not produced by this pipeline.
- **intermediate** -- an intermediate count in data acquisition, not persisted as a single artifact.

### Global cross-species coherence and its controls -- Results: configuration conserved (Fig 1); controls in Methods and S1 Text §1/§7/§9/§11


#### Configuration conserved -- primary coherence and lineage-stratified null (Results §1, Fig 1B/1C)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| matched cell types | 35 | `scripts/02_qc_and_normalize.py` | `output/phase2/scaled_35types/centroids_human_35.csv` | shape[0] |  | mapped |
| shared 1:1 ortholog genes | 16,959 | `src/cellwarp/procrustes.py` | `output/phase2/scaled_35types/centroids_human_35.csv` | shape[1] |  | mapped |
| PCA components | 33 | `scripts/08_scaled_procrustes.py` | `analysis/permutation_1M/results_1M.json` | `n_pca_components` |  | mapped |
| PCA variance retained | 95.2% | `scripts/08_scaled_procrustes.py` | `output/validation/pca_sensitivity/pca_sensitivity_results.json` | `ref.variance_explained` |  | mapped |
| primary obs/null | 0.522 | `scripts/08_scaled_procrustes.py` | `output/phase2/scaled_35types/procrustes_results_35.json` | computed: distance / null_median |  | validate.py (Global coherence obs/null (35 types)) |
| range of central 95% null | 0.507–0.544 | `scripts/permutation_1M.py` | `analysis/permutation_1M/results_1M.json` | `null_distribution_summary.percentile_2_5/97_5` (computed ratio) |  | mapped |
| primary p-value | < 10⁻⁶ | `scripts/permutation_1M.py` | `analysis/permutation_1M/results_1M.json` | `p_value` |  | validate.py (Global coherence p-value (1M primary)) |
| permutations | 1,000,000 | `scripts/permutation_1M.py` | `analysis/permutation_1M/results_1M.json` | `n_permutations` |  | mapped |
| lineage-strat obs/null | 0.668 | `scripts/test_lineage_stratified_permutation.py` | `output/validation/lineage_stratified/lineage_stratified_results.json` | `stratified_null.obs_null_ratio` |  | mapped |
| lineage-strat p | 0.0001 | same | same | `stratified_null.p_value` |  | mapped |
| within-lineage null tighter | 21.9% | `scripts/test_lineage_stratified_permutation.py` | `output/validation/lineage_stratified/lineage_stratified_results.json` | `distribution_comparison.pct_tighter` |  | mapped |

#### Macaque extension (→ S1 Text §7)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| macaque divergence | ~25 Mya | (manuscript anchor; ref 15) | n/a | n/a |  | anchor |
| macaque type count | n=12 | `analysis/macaque/reconstruct_macaque_pipeline.py` | `output/macaque_pipeline/reconstruction_qu12_results.json` | `types_included` length |  | mapped |
| macaque obs/null | 0.810 | same | same | `permutation_test.obs_null_ratio_median` |  | mapped |
| macaque p | 0.0043 | same | same | `permutation_test.p_value` |  | mapped |

#### Mouse-lemur replication (→ Results §1, Fig 1D)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| mouse lemur divergence | ~75 Mya | (manuscript anchor; ref 15) | n/a | n/a |  | anchor |
| mouse lemur type count | n=15 | `analysis/mouse_lemur/01_run_pipeline.py` | `analysis/mouse_lemur/procrustes_results.json` | `n_types` |  | mapped |
| mouse lemur obs/null | 0.346 | same | same | `permutation_test.obs_null_ratio` |  | mapped |
| mouse lemur p | 0.0001 | same | same | `permutation_test.p_value` |  | mapped |

#### Bootstrap and leave-one-out cross-validation (→ S1 Text §11)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| bootstrap iterations | 100 | `scripts/07_bootstrap.py` | `output/phase3/bootstrap/bootstrap_summary.json` | `n_bootstrap` |  | mapped |
| bootstrap subsample fraction | 50% | same | same | `subsample_fraction` |  | mapped |
| bootstrap CV | 0.004 | same | same | `distances.cv` |  | mapped |
| bootstrap stability gate | CV=0.2 | same | same | `gate_criterion.threshold` |  | mapped |
| bootstrap fraction sig | 100/100 | same | same | `p_values.fraction_significant_001` |  | mapped |
| LOOCV correct count | 35/35 | `scripts/08_loocv.py` | `output/validation/v2_loocv/v2_loocv_results.json` | `n_correct/n_total` |  | mapped |
| LOOCV mean ratio | 0.4201 | same | same | `mean_ratio` |  | mapped |
| LOOCV improvement | 58% | (computed: 1 − 0.4201 ≈ 0.58) | same | same |  | computed |

#### Independent-PCA sensitivity (→ S1 Text §1)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| indep-PCA obs/null | 0.473 | `analysis/independent_pca_sensitivity/run_independent_pca.py` | `analysis/independent_pca_sensitivity/independent_pca_results.json` | `permutation_test.obs_null_ratio` |  | validate.py (Independent PCA obs/null) |
| indep-PCA p | < 10⁻⁶ | same | same | `permutation_test.p_value` |  | mapped |
| indep-PCA Spearman ρ | 0.915 | same | same | `comparison_to_joint_pca.spearman_rho` |  | mapped |
| indep-PCA dramatic shifts | 2 of 35 | same | same | `comparison_to_joint_pca.n_dramatic_rank_changes` |  | mapped |

#### Mantel test (→ Methods: Procrustes framework / S1 Text §1)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| Mantel Pearson r | 0.789 | (Mantel) | `analysis/mantel_test/mantel_results.json` | `pearson_r` |  | mapped |
| Mantel Spearman rho | 0.737 | same | same | `spearman_r` |  | mapped |
| Mantel p (both) | < 0.001 | same | same | `pearson_p`, `spearman_p` |  | mapped |
| Mantel permutations | 10,000 | same | same | `n_permutations` |  | mapped |

#### Within-species and donor-split controls (→ S1 Text §9)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| within-species pairs | 24 | `analysis/expanded_negative_controls/expanded_negative_controls.py` | `analysis/within_species_matched/matched_control_results.json` | `all_pairs_stats.n` |  | mapped |
| within-species mean obs/null | 0.466 | same | same | `all_pairs_stats.mean_obs_null` |  | mapped |
| fraction more coherent | 71% | same | same | `all_pairs_stats.fraction_lower_than_cross_species` |  | mapped |
| random balanced splits | 100 | `analysis/donor_split/...` | `analysis/donor_split/donor_split_shared_pca_results.json` | `n_splits` |  | mapped |
| feasible cell types | 32/35 | same | same | `infeasible_types` length |  | mapped |
| median types per split | 31 | same | same | `within_species.median_n_types` |  | mapped |
| median PCA components | 32 | same | same | `pca_components` (median) |  | mapped |
| self-comparison value | 0.033 | same | same | `reference_values.self_comparison_obs_null` |  | mapped |
| within-species median | 0.375 | same | same | `within_species.median_obs_null_ratio` |  | mapped |
| within-species 95% CI | 0.318–0.423 | same | same | `within_species.ci_95` |  | mapped |
| cross-species median | 0.527 | same | same | `cross_species_matched.median_obs_null_ratio` |  | mapped |
| cross-species 95% CI | 0.485–0.584 | same | same | `cross_species_matched.ci_95` |  | mapped |
| delta median | +0.159 | same | same | `delta.median` |  | mapped |
| delta 95% CI | +0.100–+0.218 | same | same | `delta.ci_95` |  | mapped |
| positive splits | 100/100 | same | same | `delta.pct_positive` |  | mapped |
| indep-PCA delta | +0.158 | `analysis/donor_split/...` | `analysis/donor_split/donor_split_results.json` | `delta.median` |  | mapped |

### Conserved identity genes: per-gene conservation score C and master-TF enrichment -- Results §5 (Fig 5); extended in S1 Text §12

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| genes with defined C (of 16,959) | 15,940 | `analysis/conserved_contribution/make_table_s11.py` | `analysis/conserved_contribution/gate_results.json` | `n_valid` |  | mapped |
| C vs procrustes_contribution ρ (non-circularity) | 0.27 | `analysis/conserved_contribution/run_gate.py` | same | `rho_pearson_vs_circular_loading` |  | mapped |
| C vs expression Spearman ρ | 0.22 | same | same | `check2.spearman_C_vs_expr` |  | validate.py (Conserved-contribution: C vs expression Spearman) |
| C vs Tau specificity Spearman ρ | 0.06 | `analysis/conserved_contribution/run_robustness.py` | `analysis/conserved_contribution/robustness_results.json` | `rho_C_tau` |  | validate.py (Conserved-contribution: C vs specificity (Tau) Spearman) |
| Hartigan dip D (broad continuum) | 0.007 (p = 2.8 × 10⁻⁵) | `analysis/conserved_contribution/run_robustness.py` | `robustness_results.json`, `gate_results.json` | `dip_D`, `check1.dip_p` |  | validate.py (Conserved-contribution: Hartigan dip statistic D (broad continuum)) |
| master-TF median C-percentile | 0.94 | `analysis/conserved_contribution/run_gate.py` | `gate_results.json` | `check3a.median_Crank` |  | validate.py (Conserved-contribution: master-TF median C-percentile) |
| expression-matched background percentile | 0.54 | same | same | `check3a.null_median_Crank` |  | mapped |
| joint expr+specificity-matched background percentile | 0.76 | `analysis/conserved_contribution/highN_tf_pvalues.py` | `analysis/conserved_contribution/highN_tf_pvalues.json` | `joint_expr_tau_matched.null_median` |  | mapped |
| master-TF enrichment p vs both backgrounds (10⁶ draws) | < 10⁻⁶ (0 exceedances) | `analysis/conserved_contribution/highN_tf_pvalues.py` | `highN_tf_pvalues.json` | `*.p_empirical`, `*.exceedances` |  | mapped |
| TF fold-enrichment at conserved end (joint-matched) | 1.67 | `analysis/conserved_contribution/run_gate.py` | `gate_results.json` | `check3b.H_tf.fold` |  | mapped |
| CellMarker fold-enrichment (raw) | 1.89 | same | same | `check3b.H_cellmarker.fold` |  | mapped |
| conserved-set obs/null (geometry attribution) | 0.384 | same | same | `secondary.conserved.ratio` |  | validate.py (Conserved-contribution: conserved-set obs/null (geometry attribution)) |
| divergent-set obs/null | 0.709 | same | same | `secondary.divergent.ratio` |  | validate.py (Conserved-contribution: divergent-set obs/null (geometry attribution)) |
| expression-matched-random obs/null | 0.525 ± 0.012 | same | same | `secondary.matched_random_ratio_mean` (`_sd`) |  | validate.py (Conserved-contribution: expr-matched-random obs/null (geometry attribution)) |
| all-genes obs/null anchor | 0.522 | same | same | `secondary.validity_all_genes_ratio` |  | validate.py (Conserved-contribution: all-genes obs/null anchor (geometry attribution)) |
| donor-split cross-half C Spearman (median) | 0.80 | (donor-stability run) | `analysis/conserved_contribution/donor_stability/donor_stability_results.json` | `donor_split_cap10000.cross_half_C_spearman_median` |  | validate.py (Conserved-contribution: donor-split cross-half C Spearman) |
| top-quartile conserved-set Jaccard (donor-sensitivity) | 0.58 | `analysis/conserved_contribution/run_gate.py` | `gate_results.json` | `check4.highcount_jaccard` |  | mapped |
| cross-protocol C Spearman (10x vs Smart-seq2) | 0.59 | (donor-stability run) | `donor_stability/donor_stability_results.json` | `cross_protocol.spearman_C10x_vs_CSS` |  | validate.py (Conserved-contribution: cross-protocol C Spearman) |
| fresh-pull obs/null (Census re-acquisition) | 0.521 | (donor-stability run) | `donor_stability/donor_stability_results.json` | `validity.obs_null_full` |  | validate.py (Conserved-contribution: fresh-pull obs/null) |

#### Selection/derangement circularity control (→ S1 Text §10)

Tests whether selecting genes on the conservation score C manufactures the conserved-gene
geometry. The wrapper imports the published pipeline unmodified and only re-selects its inputs;
per draw it deranges the 35 mouse centroid rows, recomputes C, re-selects the top quartile, and
re-runs the obs/null. Pre-registered PASS conditions are stated in the script and echoed into
each summary under `preregistered_criteria_SURFACED_not_declared`. Six of the values below are
gated by `reproduce/validate.py`, which reads the deposited summaries. The conditions themselves
are stored as booleans and the harness compares numbers, so what is gated is their numeric
backing: the 1st percentile that condition 1 compares against (the real 0.384 it compares is
already gated from `gate_results.json`, above), and the z that condition 2 bounds.

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| real conserved obs/null (reproduced from scratch) | 0.384 | `analysis/selection_null/selection_null.py` | `analysis/selection_null/outputs/selection_null_summary_derangement.json` | `real.conserved_obs_null` | `selection_null.obs_null_ratio` | mapped (same quantity gated from `gate_results.json`) |
| full-space obs/null (reproduced from scratch) | 0.522 | same | same | `real.full_space_obs_null` | same | mapped |
| derangement sigma-null mean (± sd) | 0.991 ± 0.021 | same | same | `sigma_null.mean`, `sigma_null.sd` | `selection_null.make_draws` (mode=derangement) | validate.py (Selection null: derangement sigma-null mean) |
| derangement sigma-null 1st percentile | 0.927 | same | same | `sigma_null.p01` |  | validate.py (Selection null: derangement sigma-null 1st percentile) |
| derangement z | −29.5 | same | same | `real_position.z` |  | validate.py (Selection null: derangement z (real vs sigma-null)) |
| derangement draws at or below real | 0 of 1000 | same | same | `real_position.n_draws_at_or_below_real` |  | validate.py (Selection null: derangement draws at or below real) |
| label-shuffle sigma-null mean (± sd) | 0.983 ± 0.024 | same | `selection_null_summary_labelshuffle.json` | `sigma_null.mean`, `sigma_null.sd` | `make_draws` (mode=labelshuffle) | validate.py (Selection null: label-shuffle sigma-null mean) |
| label-shuffle z | −25.3 | same | same | `real_position.z` |  | validate.py (Selection null: label-shuffle z (real vs sigma-null)) |
| label-shuffle draws at or below real | 0 of 1000 | same | same | `real_position.n_draws_at_or_below_real` |  | mapped |
| conserved quartile size (every draw) | 3,985 of 15,940 valid | same | both summaries | `substrate.n_conserved_quartile`, `sigma_null.n_conserved_per_draw` |  | mapped |
| Q75 of C collapses under derangement | 0.59 → 0.08 | same | `sigma_null_draws_derangement.csv` | `substrate.Q75_real`; per-draw `q75` column |  | mapped (collapsed value is CSV-only; see Known gaps) |
| both pre-registered PASS conditions met | true, both modes | same | both summaries | `preregistered_criteria_SURFACED_not_declared.real_below_1st_percentile`, `.z_le_minus3` |  | mapped (numeric backing gated above) |

### Two-layer decomposition: centroid position and within-type covariance -- Results §2 (Fig 2)


#### Human–mouse two layers (Results §2, Fig 2A/2B)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| Layer 2 S at k=5 | 0.483 | `scripts/t3b_ellipsoid_alignment.py` | `output/mechanistic/ellipsoid_alignment/summary_stats.json` | `35type.mean_alignment.k=5.pre` |  | mapped |
| Layer 2 null mean | 0.375 | same | `output/mechanistic/ellipsoid_alignment/permutation_results.json` | `35type.label_shuffle_pre.k=5.null_mean` |  | mapped |
| Layer 2 p (k=5) | 0.0001 | same | same | `35type.label_shuffle_pre.k=5.p_value` |  | mapped |
| Layer 2 permutations | 10,000 | same | same | (configured) |  | mapped |
| post-rotation S | 0.230 | same | `output/mechanistic/ellipsoid_alignment/summary_stats.json` | `35type.mean_alignment.k=5.post` |  | mapped |
| post-rotation null mean | 0.180 | same | `output/mechanistic/ellipsoid_alignment/permutation_results.json` | `35type.label_shuffle_post.k=5.null_mean` |  | mapped |
| post-rotation p | 0.0001 | same | same | `35type.label_shuffle_post.k=5.p_value` |  | mapped |
| Layer-1 vs Layer-2 ρ | −0.266 | same | `output/mechanistic/ellipsoid_alignment/35type_rigidity_correlation.csv` | `rho` (k=3, metric=pre) |  | mapped |
| Layer-1 vs Layer-2 p | 0.123 | same | same | same |  | mapped |

#### PanSci Layer-2 replication (→ Results §2 / S1 Text §6)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| Layer 2 PanSci pre-rotation S (k=5) | 0.396 | `scripts/t3b_ellipsoid_alignment_pansci.py` | `output/twolayer_pansci_replication/pansci_layer2_summary.json` | `layer2_pre_rotation.k5.S` |  | mapped |
| Layer 2 PanSci pre-rotation null mean | 0.302 | same | same | `layer2_pre_rotation.k5.null_mean` |  | mapped |
| Layer 2 PanSci post-rotation S (k=5) | 0.402 | same | same | `layer2_post_rotation.k5.S` |  | mapped |
| Layer 2 PanSci post-rotation null mean | 0.360 | same | same | `layer2_post_rotation.k5.null_mean` |  | mapped |
| Layer 2 PanSci p (pre & post, k=5) | < 10⁻⁴ | same | same | `layer2_pre_rotation.k5.p`, `layer2_post_rotation.k5.p` |  | mapped |

#### CPC1 drivers under per-gene standardization (→ S1 Text §3)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| CPC1 ribosomal-dominated types (primary weighting) | 25/35 | `analysis/sensitivity_analyses/genestd_standardization.py` | `analysis/sensitivity_analyses/genestd_results.json` | `cpc1.base.n_ribosomal_dominated` |  | mapped |
| CPC1 ribosomal-dominated types (Scheme B per-gene std) | 1/35 | same | same | `cpc1.B.n_ribosomal_dominated` |  | mapped |
| Layer-1 obs/null under per-gene std (Scheme A / B) | 0.606 / 0.487 | same | same | `layer1.A.obs_null`, `layer1.B.obs_null` |  | mapped |
| per-type ranking ρ vs primary (Scheme A / B) | 0.54 / 0.76 | same | same | `layer1.A.ranking_rho_vs_primary`, `layer1.B.ranking_rho_vs_primary` |  | mapped |

### Primate covariance replication -- basal ganglia, three pairs -- Results §2 (Fig 2C)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| BG Human–Macaque: post-rotation compression (k=5, W2) | 0.628 (< 1, compresses) | (vendored: basal-ganglia deposit) | `docs/submission/plosone/figures/bg_results/layer2_results_Human_Macaque.json` | `W2_schemeB.layer2.k5.compression_ratio_post_over_pre` |  | validate.py (BG Human-Macaque: post-rotation compression (k=5, W2)) |
| BG Human–Macaque: post-rotation permutation p (k=5, W2) | 9.999 × 10⁻⁵ | (vendored: basal-ganglia deposit) | `docs/submission/plosone/figures/bg_results/layer2_results_Human_Macaque.json` | `W2_schemeB.layer2.k5.p_post` |  | validate.py (BG Human-Macaque: post-rotation permutation p (k=5, W2)) |
| BG Human–Macaque: identity-marker rank-1 drivers (W2) | 18 (of 55) | (vendored: basal-ganglia deposit) | `docs/submission/plosone/figures/bg_results/layer2_results_Human_Macaque.json` | `W2_schemeB.rank1_class_counts.canonical identity marker` |  | validate.py (BG Human-Macaque: identity-marker driver count (W2)) |
| BG Human–Macaque: matched cell types | 55 | (vendored: basal-ganglia deposit) | `docs/submission/plosone/figures/bg_results/layer2_results_Human_Macaque.json` | `n_types` |  | validate.py (BG Human-Macaque: n_types) |
| BG Human–Marmoset: post-rotation compression (k=5, W2) | 0.594 (< 1, compresses) | (vendored: basal-ganglia deposit) | `docs/submission/plosone/figures/bg_results/layer2_results_Human_Marmoset.json` | `W2_schemeB.layer2.k5.compression_ratio_post_over_pre` |  | validate.py (BG Human-Marmoset: post-rotation compression (k=5, W2)) |
| BG Human–Marmoset: post-rotation permutation p (k=5, W2) | 9.999 × 10⁻⁵ | (vendored: basal-ganglia deposit) | `docs/submission/plosone/figures/bg_results/layer2_results_Human_Marmoset.json` | `W2_schemeB.layer2.k5.p_post` |  | validate.py (BG Human-Marmoset: post-rotation permutation p (k=5, W2)) |
| BG Human–Marmoset: identity-marker rank-1 drivers (W2) | 7 (of 52) | (vendored: basal-ganglia deposit) | `docs/submission/plosone/figures/bg_results/layer2_results_Human_Marmoset.json` | `W2_schemeB.rank1_class_counts.canonical identity marker` |  | validate.py (BG Human-Marmoset: identity-marker driver count (W2)) |
| BG Human–Marmoset: matched cell types | 52 | (vendored: basal-ganglia deposit) | `docs/submission/plosone/figures/bg_results/layer2_results_Human_Marmoset.json` | `n_types` |  | validate.py (BG Human-Marmoset: n_types) |
| BG Macaque–Marmoset: post-rotation compression (k=5, W2) | 0.704 (< 1, compresses) | (vendored: basal-ganglia deposit) | `docs/submission/plosone/figures/bg_results/layer2_results_Macaque_Marmoset.json` | `W2_schemeB.layer2.k5.compression_ratio_post_over_pre` |  | validate.py (BG Macaque-Marmoset: post-rotation compression (k=5, W2)) |
| BG Macaque–Marmoset: post-rotation permutation p (k=5, W2) | 9.999 × 10⁻⁵ | (vendored: basal-ganglia deposit) | `docs/submission/plosone/figures/bg_results/layer2_results_Macaque_Marmoset.json` | `W2_schemeB.layer2.k5.p_post` |  | validate.py (BG Macaque-Marmoset: post-rotation permutation p (k=5, W2)) |
| BG Macaque–Marmoset: identity-marker rank-1 drivers (W2) | 5 (of 52) | (vendored: basal-ganglia deposit) | `docs/submission/plosone/figures/bg_results/layer2_results_Macaque_Marmoset.json` | `W2_schemeB.rank1_class_counts.canonical identity marker` |  | validate.py (BG Macaque-Marmoset: identity-marker driver count (W2)) |
| BG Macaque–Marmoset: matched cell types | 52 | (vendored: basal-ganglia deposit) | `docs/submission/plosone/figures/bg_results/layer2_results_Macaque_Marmoset.json` | `n_types` |  | validate.py (BG Macaque-Marmoset: n_types) |

### Cross-atlas replications of the configuration -- Results §3 (Fig 3); inventory in S1 Text §6 / S2 Fig(F)


#### Replication observed-to-null (Results §3, Fig 3)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| Sun2023 type count | n=15 | `scripts/17_sun2023_expanded.py` | `output/validation/sun2023_replication_expanded/sun2023_expanded.json` | `procrustes.n_types` |  | mapped |
| Sun2023 obs/null | 0.554 | same | same | `procrustes.obs_null_ratio` |  | mapped |
| Sun2023 p | 0.0001 | same | same | `procrustes.p_value` |  | mapped |
| PanSci type count | n=16 | `scripts/pansci_replication.py` | `output/validation/pansci_replication/pansci_replication.json` | `procrustes.n_types` |  | mapped |
| PanSci obs/null | 0.552 | same | same | `procrustes.obs_null_ratio` |  | mapped |
| PanSci p | 0.0001 | same | same | `procrustes.p_value` |  | mapped |
| CellHint cells | 2.87M | (CellHint dataset metadata) | (anchor: Xu et al. 2023) | n/a |  | anchor |
| CellHint type count | n=15 | `scripts/33_cellhint_replication.py` | `output/validation/cellhint_replication/cellhint_replication.json` | `procrustes.n_types` |  | mapped |
| CellHint obs/null | 0.448 | same | same | `procrustes.obs_null_ratio` |  | validate.py (CellHint obs/null) |
| CellHint p | 0.0001 | same | same | `procrustes.p_value` |  | validate.py (CellHint p-value) |
| pan-Census type count | n=22 | `analysis/census_replication/02_run_replication.py` | `analysis/census_replication/replication_results.json` | `n_cell_types` |  | mapped |
| pan-Census datasets | 15 | `data/replication/pan_census_manifest.csv` | same CSV | `df.shape[0]` |  | mapped |
| pan-Census obs/null | 0.811 | same | `analysis/census_replication/replication_results.json` | `permutation_test.obs_null_ratio` |  | mapped |
| pan-Census p | 0.0001 | same | same | `permutation_test.p_value` |  | mapped |
| Andrews obs/null | 0.797 | `scripts/31_andrews_replication.py` | `output/validation/andrews_replication/andrews_replication_results.json` | `obs_null_ratio` | `31_andrews_replication.main` | mapped |
| Andrews p | 0.1159 | same | same | `p_value` | `31_andrews_replication.main` | mapped |
| Andrews n | 6 | same | same | `n_types` | `31_andrews_replication.main` | mapped |
| MCA × HCA obs/null | 1.003 | `scripts/14_t1a_replication.py` | `output/validation/t1a_replication/t1a_results.json` | `t1a_procrustes.obs_null_ratio` | `14_t1a_replication.main` | mapped |
| MCA × HCA p | 0.542 | same | same | `t1a_procrustes.p_value` | `14_t1a_replication.main` | mapped |
| MCA × HCA n | 17 | same | same | `t1a_procrustes.n_types` | `14_t1a_replication.main` | mapped |
| HCA × Tabula obs/null | 0.728 | `scripts/15_hca_centroid_comparison.py` | `output/validation/hca_centroid_comparison/hca_centroid_comparison.json` | `comparison_a.obs_null_ratio` | `15_hca_centroid_comparison.main` | mapped |
| HCA × Tabula p | 0.003 | same | same | `comparison_a.p_value` | `15_hca_centroid_comparison.main` | mapped |
| HCA × Tabula n | 6 | same | same | `comparison_a.n_types` | `15_hca_centroid_comparison.main` | mapped |
| Andrews scaling | 0.229 | `scripts/31_andrews_replication.py` | `output/validation/andrews_replication/andrews_replication_results.json` | `scaling` | `31_andrews_replication.main` | mapped |
| MCA scaling | 0.267 | `scripts/14_t1a_replication.py` | `output/validation/t1a_replication/t1a_results.json` | `t1a_procrustes.scaling` | `14_t1a_replication.main` | mapped |

#### Per-replication cross-atlas ranking (→ Results §4 / S3 Fig B)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| Sun ranking ρ | +0.15 | same Sun2023 file | same | `rigidity_ranking.rho` |  | mapped |
| Sun ranking p | 0.60 | same | same | `rigidity_ranking.p_value` |  | mapped |
| PanSci ranking ρ | +0.19 | `output/validation/pansci_replication/pansci_replication.json` | same | `rigidity_ranking.rho` |  | mapped |
| PanSci ranking p | 0.47 | same | same | `rigidity_ranking.p_value` |  | mapped |
| CellHint ranking ρ | −0.39 | `output/validation/cellhint_replication/cellhint_replication.json` | same | `rigidity_ranking.rho` |  | mapped |
| CellHint ranking p | 0.16 | same | same | `rigidity_ranking.p_value` |  | mapped |
| pan-Census ranking ρ | −0.05 | `analysis/census_replication/replication_results.json` | same | `ranking_correlation.spearman_rho` |  | mapped |
| CellHint harmonized ρ | −0.04 | `analysis/harmonized_replication/correlation_results.json` | same | `rho` |  | mapped |
| harmonized n | 12 | same | same | `n_types` |  | mapped |

### Macaque extension: matched-n controls and per-type ranking dependence -- S1 Text §7

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| three-way ortholog space | 13,927 genes | `analysis/macaque/reconstruct_macaque_pipeline.py` | `output/macaque_pipeline/reconstruction_qu12_results.json` | `gene_space` |  | mapped |
| three-way fraction | 82.1% | (computed: 13927/16959) | n/a | n/a |  | computed |
| macaque-vs-HM obs/null | 0.810 | same | `output/macaque_pipeline/reconstruction_qu12_results.json` | `permutation_test.obs_null_ratio_median` |  | mapped |
| matched-n HM obs/null | 0.440 | same | `output/macaque_pipeline/human_mouse_12type_control.json` | `control_16959.permutation_test.obs_null_ratio_median` |  | mapped |
| matched-n HM p | 0.0001 | same | same | `control_16959.permutation_test.p_value` |  | mapped |
| no-immune subset n | 7 | same | `output/macaque_pipeline/reconstruction_qu7_D1_results.json` | `n_types_analyzed` |  | mapped |
| no-immune subset obs/null | 0.733 | same | same | `permutation_test.obs_null_ratio_median` |  | mapped |
| no-immune p | 0.013 | same | same | `permutation_test.p_value` |  | mapped |
| matched-n HM 7-type obs/null | 0.485 | same | `output/macaque_pipeline/human_mouse_7type_control.json` | `permutation_test.obs_null_ratio_median` |  | mapped |
| matched-n HM 7-type p | 0.0006 | same | same | `permutation_test.p_value` |  | mapped |
| divergence ranking ρ | 0.147 | same | `output/macaque_pipeline/m1_close_table1_summary.json` | `spearman_ranking_correlation.rho` |  | mapped |
| divergence ranking p | 0.649 | same | same | `spearman_ranking_correlation.p_value` |  | mapped |

#### Mouse-lemur ranking (→ Methods: mouse-lemur parameters / Results §1)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| mouse lemur ρ vs primary | 0.157 | `analysis/mouse_lemur/01_run_pipeline.py` | `analysis/mouse_lemur/ranking_correlation.json` | `vs_primary.rho` |  | mapped |
| mouse lemur p | 0.576 | same | same | `vs_primary.p_value` |  | mapped |

#### Hepatocyte rank/SSR reversal (S1 Text §7)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| hepatocyte rank HM | 12/12 | same | `output/macaque_pipeline/m1_close_table1_summary.json` | `hepatocyte_rank_reversal.hm12_rank` |  | mapped |
| hepatocyte SSR HM | 2.0% | same | `output/macaque_pipeline/human_mouse_12type_control.json` | `hepatocyte_pct_ssr` |  | mapped |
| hepatocyte rank macaque | 1/12 | same | `output/macaque_pipeline/m1_close_table1_summary.json` | `human_macaque_12type.hepatocyte_rank_of_n` |  | mapped |
| hepatocyte SSR macaque | 47.3% | same | same | `human_macaque_12type.hepatocyte_pct_ssr` |  | mapped |

### Biological predictors of per-type divergence -- S1 Table (per-type unresolvability = Results §4)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| features tested | 15 | `analysis/biological_predictors/biological_predictors.py` | `analysis/biological_predictors/univariate_correlations.csv` | `df.shape[0]` |  | mapped |
| cell count ρ | 0.052 | `scripts/confound_cellcount_rigidity.py` | `output/cellcount_confound/cellcount_confound_results.json` | `spearman_rho` |  | mapped |
| tissue breadth ρ | −0.13 | `analysis/biological_predictors/biological_predictors.py` | `analysis/biological_predictors/univariate_correlations.csv` | "Tissue breadth" row |  | mapped |
| mean expression ρ | −0.06 | same | same | "Mean expression level" row |  | mapped |
| progenitor ρ | 0.43 | same | same | "Is progenitor" row |  | mapped |
| progenitor p | 0.01 | same | same | same |  | mapped |
| elastic net LOO R² | −0.064 | same | `analysis/biological_predictors/multivariate_model_results.json` | `loo_cv.elastic_net_r2` |  | mapped |
| random forest LOO R² | −0.044 | same | same | `loo_cv.random_forest_r2` |  | mapped |

### Per-type divergence precise within-atlas, not resolvable across -- Results §4 (Fig 4); simulation in Methods; bootstrap in S3 Fig


#### Simulation power and ranking ceiling (→ Methods: Simulation; Fig 4C)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| sim type counts | n=15,25,35 | `analysis/simulation_study/simulation_study.py` | `analysis/simulation_study/simulation_results.json` | `power_curve` |  | mapped |
| sim detection power | 100% | same | same | `power_curve` (calibrated row) |  | mapped |
| sim FPR at α=0.05 | 4.8% | same | same | `null_calibration.fpr_005` |  | mapped |
| ranking ceiling | ρ ≈ 0.42 | same | same | `ranking_recovery` |  | mapped |
| test-retest ρ | 0.994 | same | same | `stability` |  | mapped |
| cross-atlas ranking ρ | ≈ 0.15 | (descriptive across replications) | n/a | n/a |  | computed |

#### Bootstrap ranking stability (→ Results §4 / S3 Fig A)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| bootstrap iterations | 1,000 | `analysis/bootstrap_rankings/bootstrap_ranking_analysis.py` | `analysis/bootstrap_rankings/bootstrap_summary.csv` | (configured) |  | mapped |
| all stable (CI ≤ 10) | 35/35 | same | same | `(df.ci_width <= 10).all()` |  | mapped |
| median CI width | 3 | same | same | `df.ci_width.median()` |  | mapped |
| max CI width | 7 | same | same | `df.ci_width.max()` |  | mapped |
| CI=0 type count | 5 | same | same | `(df.ci_width == 0).sum()` |  | mapped |

#### Within- vs cross-atlas inversion (→ Results §4, Fig 4B / S3 Fig B)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| CI-width vs cross-atlas rank-shift ρ (reported in Fig 4B / S3 Fig B) | −0.410 (p = 0.073, n = 20) | `analysis/cross_reference/cross_reference_analysis.py` | `analysis/cross_reference/convergent_types_summary.md` | "Spearman ρ (CI width vs mean rank shift)" = −0.410; "Spearman p-value" = 0.0727 | `cross_reference_analysis.main` | mapped |
| CI vs cross-atlas p | 0.073 | same | same | "Spearman p-value" | `cross_reference_analysis.main` | mapped |
| n in ≥2 replications | 20 | same | same | "Types in ≥2 replications" | `cross_reference_analysis.main` | mapped |

### Ten mechanistic null tests for the per-type ranking -- summary claims (S1 Text §8)


#### Mechanistic nulls (S1 Text §8)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| tests run | 10 | (this section: per-test table below) | (see per-test mechanistic-null table below) | n/a |  | computed |
| n=35 (most tests) | 35 | various mechanistic-null scripts | various | various |  | computed |
| power at ρ≈0.30 | 37% | `analysis/simulation_study/simulation_study.py` | `analysis/simulation_study/simulation_results.json` | `power_curve` |  | mapped |
| TF complexity ρ | −0.229 | `scripts/13_tf_complexity.py` | `output/phase2/mechanistic/tf_complexity/tf_complexity_vs_residual.csv` | aggregated |  | mapped |
| PPI best ρ | 0.291 | `scripts/19_ppi_centrality.py` | `output/mechanistic/ppi_centrality/ppi_centrality_results.json` | `correlation_results[*].spearman_rho` (max) |  | mapped |
| enhancer ρ | −0.429 | `scripts/t3e_step3b_enhancer.py` | `output/validation/t3e_enhancer/t3e_step3b_summary.md` | "Spearman ρ" |  | mapped |
| enhancer n | 6 | same | same | "n (cell types)" |  | mapped |

#### CellHint replication diagnostics (→ S1 Text §6 / S3–S4 Table)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| CellHint matched 15-type ρ | −0.139 | `analysis/harmonized_replication/harmonized_replication.py` | `analysis/harmonized_replication/sensitivity_analysis.csv` | `0_unharmonized` row |  | mapped |
| CellHint matched p | 0.62 | same | same | same row |  | mapped |
| CellHint harmonized ρ | −0.04 | same | `analysis/harmonized_replication/correlation_results.json` | `rho` |  | mapped |
| CellHint harmonized n | 12 | same | same | `n_types` |  | mapped |

### Primary-pipeline parameters and counts -- Methods


#### Data, gene-space and pipeline parameters (Methods)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| Tabula Sapiens cells | 483,152 | (Tabula Sapiens v1.0 paper) | n/a | n/a |  | anchor |
| Tabula Sapiens donors | 24 | (Tabula Sapiens v1.0 paper) | n/a | n/a |  | anchor |
| Tabula Muris Senis cells | ~350,000 | (TMS paper) | n/a | n/a |  | anchor |
| unique CL labels in TS | 180 | (data acquisition pipeline) | n/a | n/a |  | intermediate |
| unique CL labels in TMS | 151 | (data acquisition pipeline) | n/a | n/a |  | intermediate |
| shared CL labels | 66 | (data acquisition pipeline) | n/a | n/a |  | intermediate |
| reaching ≥200 cells | 45 | (data acquisition pipeline) | n/a | n/a |  | intermediate |
| meeting ≥500 cells | 35 | `scripts/02_qc_and_normalize.py` | `output/phase2/scaled_35types/centroids_human_35.csv` | shape[0] |  | mapped |
| 10x-only ranking ρ | 0.754 | `scripts/14_smartseq2_sensitivity.py` | `output/phase2/sensitivity/smartseq2/sensitivity_results.json` | `rigidity_comparison.spearman_rho_ranks` | `14_smartseq2_sensitivity.main` | mapped |
| progenitor full p | 0.0099 | (biological_predictors progenitor sub-analysis) | (multiple sources) | n/a |  | computed |
| progenitor 10x-only p | 0.119 | same | same | n/a |  | computed |
| mouse SS2 fraction ρ | −0.042 | `scripts/14_smartseq2_sensitivity.py` | `output/phase2/sensitivity/smartseq2/sensitivity_results.json` | `rigidity_comparison.ss2_fraction_vs_rank_change_rho` | `14_smartseq2_sensitivity.main` | mapped |
| total 1:1 orthologs | 17,187 | (BioMart query) | `data/phase1/orthologs_human_mouse.csv` | shape[0] |  | mapped |
| operating gene space | 16,959 | `scripts/02_qc_and_normalize.py` | `output/phase2/scaled_35types/centroids_human_35.csv` | shape[1] |  | mapped |
| gene-set sensitivity range | 0.496–0.511 | `scripts/test_hvg_robustness.py` | `output/validation/hvg_robustness/hvg_robustness.json` | `hvg_results.{2000,3000,5000}.obs_null_ratio` (median denominator) |  | mapped |
| gene-set ranking ρ | ≥ 0.957 | `scripts/test_hvg_robustness.py` | `output/validation/hvg_robustness/hvg_robustness.json` | `hvg_results.{2000,3000,5000}.ranking_rho_vs_full` |  | mapped |

#### Cap, verification and negative-control parameters (Methods / S1 Text §9 / S1 Text §6)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| per-type cell cap | 2,000 | `src/cellwarp/procrustes.py` | (configured constant) | `MAX_CELLS_PER_TYPE` |  | mapped |
| min cell count per type | 516 (myeloid leukocyte mouse) | data acquisition | n/a | n/a |  | intermediate |
| scipy verification max-Δ | 7.3 × 10⁻¹¹ | `scripts/verify_procrustes_vs_scipy.py` | (verification log) | n/a |  | computed |
| within-species v2 cells | 11,640 | `scripts/test_35type_human_control.py` | `output/phase2/negative_control_v2/cell_availability.json` | (sum of selected cells) |  | mapped |
| v2 within-species obs/null | 0.607 | same | `output/phase2/negative_control_v2/negctrl_v2_results.json` | (within-species comparison) |  | mapped |
| v2 cross-species obs/null | 0.317 | same | same | (matched-scale 6-type comparison) |  | mapped |
| v2 within-species p | 0.0088 | same | same | `permutation_test.p_value` |  | mapped |
| Sun2023 lung-restricted obs/null | 0.490 | `scripts/18_sun2023_issue092_diagnosis.py` | `output/validation/sun2023_issue092_diagnosis/issue092_diagnosis.json` | `task_3_revised_procrustes.obs_null_ratio` | `18_sun2023_issue092_diagnosis.main` | mapped |
| Sun2023 lung-restricted p | 0.0001 | same | same | `task_3_revised_procrustes.p_value` | `18_sun2023_issue092_diagnosis.main` | mapped |

#### Two-layer eigenvalue statistics (→ Results §2 / Methods)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| sorted-eigenvalue mean r | 0.953 | `scripts/layer3_permutation_test.py` | `output/layer3_permutation/layer3_permutation_results.json` | `observed_mean_r` |  | mapped |
| eigenvalue p | 0.866 | same | same | `empirical_p` |  | mapped |
| eigenvalue permutations | 10,000 | same | same | `n_permutations` |  | mapped |
| eigenval-divergence ρ | 0.395 | `scripts/t3b_ellipsoid_alignment.py` | `output/mechanistic/ellipsoid_alignment/summary_stats.json` | `35type.eigenval_vs_rigidity.rho` |  | validate.py (Ellipsoid eigenval-residual rho (35 types)) |
| eigenval-divergence p | 0.019 | same | same | `35type.eigenval_vs_rigidity.p` |  | validate.py (Ellipsoid eigenval-residual p-value (35 types)) |

#### Primary distance, LOOCV, bootstrap, cell-count (Methods)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| nearest null distance | 98.88 | `scripts/permutation_1M.py` | `analysis/permutation_1M/results_1M.json` | `null_distribution_summary.min` |  | mapped |
| observed distance | 61.15 | same | same | `observed_procrustes_distance` |  | mapped |
| replication 4× α=0.05 | 6.25 × 10⁻⁶ | (computed: 0.05⁴) | n/a | n/a |  | computed |
| best-predicted type | CD8+ T cell | `scripts/08_loocv.py` | `output/validation/v2_loocv/v2_loocv_results.json` | `best_type` |  | mapped |
| best ratio | 0.226 | same | same | `min_ratio` |  | mapped |
| hardest type | hematopoietic stem cell | same | same | `worst_type` |  | mapped |
| hardest ratio | 0.796 | same | same | `max_ratio` |  | mapped |
| bootstrap mean distance | 61.18 | `scripts/07_bootstrap.py` | `output/phase3/bootstrap/bootstrap_summary.json` | `distances.mean` |  | mapped |
| bootstrap distance std | 0.24 | same | same | `distances.std` |  | mapped |
| cell-count ρ | 0.052 | `scripts/confound_cellcount_rigidity.py` | `output/cellcount_confound/cellcount_confound_results.json` | `spearman_rho` |  | mapped |
| cell-count p | 0.768 | same | same | `spearman_p` |  | mapped |
| partial ρ | 0.061 | same | same | `partial_rho` |  | mapped |
| partial p | 0.732 | same | same | `partial_p` |  | mapped |

### Ten mechanistic null tests -- per-test table (S1 Text §8)

| Manuscript row | Test | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|---|
| Row 1 | Housekeeping ρ | 0.167 | `scripts/12_housekeeping_ratio.py` | `output/phase2/mechanistic/housekeeping/hk_ratio_results.json` | `human_correlation.spearman_rho` |  | validate.py (Housekeeping ratio vs residual Spearman rho (human)) |
| Row 2 | TF complexity ρ | −0.229 | `scripts/13_tf_complexity.py` | `output/phase2/mechanistic/tf_complexity/tf_complexity_vs_residual.csv` | aggregated |  | mapped |
| Row 3 | Niche adaptation | 0/6 sig after FDR | `scripts/12_niche_hypothesis.py` | `output/phase2/progenitor_analysis/niche_hypothesis/niche_hypothesis_results.json` | `progenitor_divergence.n_niche_sets_significant_q05` (=0) of `n_niche_sets_tested` (=6) |  | mapped |
| Row 4 | Within-type variance ρ | −0.038 | `scripts/12_variance_diagnostic.py` | `output/phase2/variance_diagnostic/diagnostic_results.json` | `spearman_rho_mean` |  | mapped |
| Row 5 | Inter-donor variance ρ | −0.127 | `scripts/16_interdonor_variance.py` | `analysis/biological_predictors/univariate_correlations.csv` | "Inter-donor variance" row |  | mapped |
| Row 6 | Expression-level all ρ < 0.21 | max abs ρ = 0.209 | `scripts/diagnostic_expression_vs_rigidity.py` | `output/phase2/diagnostics/expression_level_vs_rigidity/correlations.csv` | max(abs(rho)) |  | mapped |
| Row 7 | PPI best ρ | 0.291 | `scripts/19_ppi_centrality.py` | `output/mechanistic/ppi_centrality/ppi_centrality_results.json` | `correlation_results[*].spearman_rho` (best by abs) |  | mapped |
| Row 7 | PPI 0/27 sig after FDR | 0/27 | same | same | `n_significant_fdr / len(correlation_results)` |  | mapped |
| Row 8 | phastCons ρ | −0.058 | `scripts/t3e_step2_compute.py` | `output/validation/t3e_chromatin/t3e_step2_summary.md` | "Spearman ρ" (placental_20way primary, Option A, ±2kb) |  | mapped |
| Row 8 | phastCons n | 35 | same | same | "n (cell types)" |  | mapped |
| Row 9 | enhancer ρ | −0.429 | `scripts/t3e_step3b_enhancer.py` | `output/validation/t3e_enhancer/t3e_step3b_summary.md` | "Spearman ρ" |  | mapped |
| Row 9 | enhancer n | 6 | same | same | "n (cell types)" |  | mapped |
| Row 10 | drug target ρ | −0.176 | (output canonical) | `output/t3g/primary_correlation_results.json` | `primary_correlation.rho` |  | mapped |
| Row 10 | drug target n | 34 | (output canonical) | same | `primary_correlation.n` |  | mapped |

### Replication and macaque parameters -- Methods (macaque → S1 Text §7)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| replication permutations | 10,000 | (replication scripts) | replication JSONs | `permutation_test.n_permutations` |  | mapped |
| hepatocyte SSR macaque (full) | 47.3% | `analysis/macaque/reconstruct_macaque_pipeline.py` | `output/macaque_pipeline/m1_close_table1_summary.json` | `human_macaque_12type.hepatocyte_pct_ssr` |  | mapped |
| hepatocyte SSR HM (full) | 2.0% | same | `output/macaque_pipeline/human_mouse_12type_control.json` | `hepatocyte_pct_ssr` |  | mapped |

### Mouse-lemur pipeline parameters -- Methods (result → Results §1, Fig 1D)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| lemur passing types | 15 | `analysis/mouse_lemur/01_run_pipeline.py` | `analysis/mouse_lemur/procrustes_results.json` | `n_types` |  | mapped |
| lemur cell range | 546–2,000 | same | `analysis/mouse_lemur/lemur_cell_type_counts.csv` | per-type counts |  | mapped |
| lemur 1:1 ortholog count | 13,796 | same | `analysis/mouse_lemur/procrustes_results.json` | `gene_space` |  | mapped |
| lemur PCA components | 15 | same | same | `pca.n_components` |  | mapped |
| lemur cumulative variance | 95.5% | same | same | `pca.cumulative_variance` |  | mapped |
| lemur obs/null | 0.346 | same | same | `permutation_test.obs_null_ratio` |  | mapped |
| lemur permutations | 10,000 | same | same | `permutation_test.n_permutations` |  | mapped |
| lemur ranking ρ | 0.157 | same | `analysis/mouse_lemur/ranking_correlation.json` | `vs_primary.rho` |  | mapped |

### Pan-Census replication parameters -- Methods (result → Results §3)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| pan-Census datasets | 15 | `analysis/census_replication/02_run_replication.py` | `data/replication/pan_census_manifest.csv` | row count |  | mapped |
| pan-Census collections | 15 | same | same CSV | `df.collection_id.nunique()` |  | mapped |
| pan-Census mouse datasets | 9 | same | same | `df.species.value_counts()` |  | mapped |
| pan-Census human datasets | 6 | same | same | `df.species.value_counts()` |  | mapped |
| pan-Census shared types | 22 | same | `analysis/census_replication/replication_results.json` | `n_cell_types` |  | mapped |
| pan-Census obs/null | 0.811 | same | same | `permutation_test.obs_null_ratio` |  | mapped |
| pan-Census p | 0.0001 | same | same | `permutation_test.p_value` |  | mapped |

### Software environment and reproducibility -- Methods / S12 Table (S1 Text §13)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| Python | 3.12.12 | n/a | `pyproject.toml` | `requires-python` |  | mapped |
| RANDOM_SEED | 42 | (entry-point convention) | various | per-script `RANDOM_SEED` |  | mapped |
| numpy | 2.4.3 | n/a | `pyproject.toml` | `[lock]` |  | mapped |
| scipy | 1.17.1 | n/a | same | same |  | mapped |
| pandas | 2.3.3 | n/a | same | same |  | mapped |
| scanpy | 1.12 | n/a | same | same |  | mapped |
| anndata | 0.12.10 | n/a | same | same |  | mapped |
| scikit-learn | 1.8.0 | n/a | same | same |  | mapped |
| statsmodels | 0.14.6 | n/a | same | same |  | mapped |
| cellxgene-census | 1.17.0 | n/a | same | same |  | mapped |
| samap | 1.0.14 | n/a | same | same |  | mapped |
| Census version | 2025-11-08 | (pinned) | `analysis/census_replication/02_run_replication.py` | `census_version` |  | mapped |
| Tabula Microcebus collection_id | a137437b-… | n/a | `data/replication/tabula_microcebus_metadata.csv` | `collection_id` |  | mapped |
| Tabula Microcebus dataset_id | a392ab34-… | n/a | same | `dataset_id` |  | mapped |
| Tabula Microcebus download date | 2026-04-05 | n/a | same | `download_date` |  | mapped |
| BioMart release | 115 | n/a | `data/phase1/orthologs_human_mouse.csv` (header) | accession date in DOI |  | mapped |

### Sensitivity-analysis rows -- Layer-1 exclusion (S7 Table / S1 Text §4) and donor-split delta (S1 Text §9)


#### Layer-1 ribosomal / housekeeping exclusion (→ S7 Table / S1 Text §4)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| Layer-1 ribosomal-excluded obs/null | 0.501 | `analysis/sensitivity_analyses/layer1_exclusion.py` | `analysis/sensitivity_analyses/layer1_exclusion_results.json` | `variants.1.obs_null_ratio` |  | validate.py (Layer-1 ribosomal-excluded obs/null) |
| Layer-1 ribosomal+housekeeping obs/null | 0.479 | same | same | `variants.2.obs_null_ratio` |  | validate.py (Layer-1 ribosomal+housekeeping-excluded obs/null) |

#### Donor-split delta (→ S1 Text §9)

| Claim (excerpt) | Value | Script | Output file | Output key | Function entry-point | Status |
|---|---|---|---|---|---|---|
| Donor-split delta (cross-species − within-species), reported as the median of per-split deltas | +0.159 | (donor-split) | `analysis/donor_split/donor_split_shared_pca_results.json` | `delta.median`; the aggregate arms cross = 0.527 / within = 0.375 (`donor_split_comparison.json` `shared_pca`) are arm medians, distinct from the reported +0.159 delta; the manuscript reports neither 1.41× nor the 0.152 arm-difference |  | mapped |

---

## Notes for reviewers

- Some output files are gitignored at the directory level (e.g., `data/`,
  `output/`) but are tracked individually; see `git ls-files data/ output/`
  to enumerate the tracked set.
- This CROSSWALK is generated and maintained as a single Markdown document.
  Future updates land in this file directly.
- The deposit also includes two dataset-level manifest CSVs that travel
  with the analysis code (see README §"Deposit Artifacts"):
  `data/replication/pan_census_manifest.csv` (15 CELLxGENE Census datasets)
  and `data/replication/tabula_microcebus_metadata.csv` (Tabula Microcebus
  deposit anchors). Both are referenced from this CROSSWALK in the
  pan-Census and mouse-lemur sections respectively.
