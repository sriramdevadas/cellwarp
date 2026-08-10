# SCOPE — what is in the paper, and what is not

This repository holds more analysis than the PLOS ONE manuscript reports. That is
by design: the paper is a principled subset of a longer research program, and
several dated pre-specifications on disk cover analyses that are **reported
elsewhere or held in reserve, not in this paper**. This file classifies every
`scripts/*.py`, `analysis/*/`, `output/*/`, and preregistration so that the
boundary between reported and unreported work is explicit and auditable. Its
purpose is to make clear that unreported analyses are scoped-out on stated
grounds, not selectively suppressed.

## Categories

- **PAPER** — produces a display item, a reported statistic, or is required
  pipeline / build / validation tooling for the PLOS ONE manuscript.
- **BANKED** — a completed, pre-specified analysis held in reserve. Not reported
  in this paper, but executed and deployable if a reviewer raises the
  corresponding point.
- **OTHER PROJECT** — belongs to a distinct research line with its own separate
  (future) write-up: cancer geometry, the disease/activation axis, the
  treeness–rigidity line, the aging axis, the evolutionary-ratio (dN/dS) line,
  developmental constraint, and the 100-type ontology expansion.
- **EXPLORATORY** — scouting, feasibility, and diagnostic scripts that informed
  the paper's design but produce no reported number and make no standalone claim.

---

## Dated pre-specifications on disk (the transparency-critical list)

Six dated analysis plans are committed with the code. Two are reported in the
paper; **four are not**, and are listed here so their absence from the manuscript
is on the record rather than inferred.

| Preregistration file | Date | In paper? | Category | Disposition |
|---|---|---|---|---|
| `docs/submission/figures_for_review/Supplementary_Preregistration.md` | — | **Yes** | PAPER | The deposited pre-specification of record, mirrored into the review-packet directory from `docs/preregistration_conserved_contribution_2026-06-05.md` by `scripts/build_submission_packet.py`. |
| `docs/preregistration_conserved_contribution_2026-06-05.md` | 2026-06-05 | **Yes** | PAPER | Pre-specifies the conserved-contribution / master-TF analysis (Fig 5, Results §2, §5). |
| `docs/preregistration_aging_axis_2026-03-16.md` | 2026-03-16 | **No** | OTHER PROJECT | Aging-axis project (within-mouse aging as a directional centroid shift). **Not executed** — no output or script on disk. See "Aging/DILI check" below. |
| `docs/preregistration_dilirank_hepatocyte_2026-03-16.md` | 2026-03-16 | **No** | BANKED | Hepatocyte-rigidity → drug-induced-liver-injury landmark test. **Executed** (`output/dilirank/`), returned a marginal/null result that did not survive the pre-specified covariate-falsification gate; correctly not reported as a positive finding. See below. |
| `docs/preregistration_treeness_anticorrelation_2026-03-16.md` | 2026-03-16 | **No** | OTHER PROJECT | Treeness–rigidity mediation (neighborhood density). Treeness line; panels cut from this paper. |
| `docs/preregistration_treeness_h2_concentration_2026-03-16.md` | 2026-03-16 | **No** | OTHER PROJECT | Treeness–rigidity mediation (gene-program concentration). Post-hoc exploratory; treeness line. |

The manuscript's disclosure sentence (Methods, Statistical analysis)
is written to reflect this: no analysis was registered with an *external* public
registry (e.g. OSF), but dated internal pre-specifications were deposited with the
code, and some of them cover analyses reported elsewhere — pointing here.

### Aging / DILI check (resolves the §4 gate)

**Question:** do the `aging_axis` or `dilirank_hepatocyte` preregistrations
overlap the ten mechanistic-null tests reported for the per-type divergence
ranking (Results §4; the ten nulls are reported in S1 Text §8), and does the aging plan test
age-as-a-confound-on-the-cross-species-residual?

**Answer — no eleventh null is needed:**

- **No overlap.** The ten mechanistic nulls each correlate a biological covariate
  against the human–mouse per-type Procrustes residual (housekeeping ratio, TF
  network complexity, niche gene sets, within-type variance, inter-donor
  variance, expression-level confounds, PPI centrality, promoter phastCons,
  H3K27ac enhancer conservation, drug-target conservation). Neither
  preregistration is among them or tests a covariate against the residual in that
  form.
- **Aging does not test the named confound.** The aging plan tests (a) whether
  mouse somatic cells show a directional centroid shift with age (a *within-mouse*
  aging axis), and (b, exploratory) whether that aging displacement vector
  *aligns in direction* with the cross-species transform ("does a cell age in the
  direction it evolved"). Neither is the confound the manuscript names — whether
  *age-composition differences between the human (adult) and mouse (aged) atlases
  reorder the per-type cross-species residual*. That compositional confound is
  disclosed in the manuscript as acknowledged-but-untested (Methods, atlas
  description; Discussion, per-type bounds); the aging plan would not resolve it
  even if run.
- **Aging was not executed.** Only the plan is committed; there is no aging output
  or script. There is therefore no result to fold into §4.
- **DILI was executed and is null under its own gate.** Test 1 (Mann-Whitney)
  p ≈ 0.054 (both arms, not significant at α = 0.05); Test 2 (Fisher upper-tail)
  p ≈ 0.076 full / 0.020 CYP450-excluded; Test 3 (dip) not significant; Test 4
  (covariate-falsification partial correlation) p ≈ 0.087 / 0.102 — **fails** the
  pre-specified requirement that the partial correlation remain significant,
  triggering the plan's Hard Abort 3. Per its own reporting commitment it is a
  banked negative/marginal result, not reported as a finding. (Its Step-1
  sensitivity gate — L1000-landmark preservation of the ranking, ρ = 0.852 — is
  checked by `reproduce/validate.py` and so is retained as validation tooling,
  but its result is not reported anywhere in the manuscript body or SI; the
  toxicity claim is not reported either.)

Conclusion: §4's ten-null framing stands; aging, DILI, and the two treeness plans
are scoped out here rather than in the manuscript body.

---

## `scripts/` classification

### PAPER — pipeline, analysis, figures, tables, validation, build

| Script | Role |
|---|---|
| `00_verify_env.py`, `01_download_data.py`, `02_qc_and_normalize.py` | Environment / data acquisition / QC + normalization pipeline |
| `04_procrustes_analysis.py`, `04b_procrustes_sensitivity.py`, `04c_negative_control.py` | Core Procrustes analysis + sensitivity + early negative control |
| `08_scaled_procrustes.py` | Scaled Procrustes, 35 types (per-type ranking; Fig 4-family) |
| `permutation_1M.py` | 1M-permutation headline null (Fig 1B) |
| `permutation_1M_independent_pca.py` | Independent-PCA 1M null (Fig S1A) |
| `test_lineage_stratified_permutation.py` | Lineage-stratified null (Fig 1C) |
| `07_bootstrap.py` | Configuration bootstrap stability (Methods; S1 Text §11, text-only — no figure) |
| `08_loocv.py` | Leave-one-out CV (Methods; S1 Text §11, text-only — no figure) |
| `09_mantel_test.py` | Mantel pairwise-distance conservation |
| `09_negative_control_v2.py`, `test_35type_human_control.py`, `49_build_figS7_matched_scale.py`, `polish_figS7.py` | Human-vs-human negative controls / matched-scale (S4 Fig, Fig 3-family) |
| `16_sun2023_replication.py`, `pansci_replication.py`, `pansci_metadata_gate.py`, `33_cellhint_replication.py`, `14_t1a_replication.py`, `13_replication_inventory.py` | Direct replications (Fig 3A–C, replication inventory) |
| `12_t1a_mca_download.py`, `13_t1a_hca_download.py`, `31_andrews_replication.py`, `15_hca_centroid_comparison.py` | Replication data acquisition + Andrews / MCA×HCA non-replications + within-human diagnostic (Fig S2F, §"config robust") |
| `56_add_figs2_panel_f.py`, `patch_figs2_panel_f_values.py` | Replication-inventory panel (Fig S2F) |
| `14_smartseq2_sensitivity.py` | Smart-seq2 protocol sensitivity (Fig S2C–D) |
| `17_pca_sensitivity.py`, `18_pca_sensitivity_v2.py` | PCA-dimension sensitivity (Fig S2A–B) |
| `t3b_ellipsoid_alignment.py`, `layer3_permutation_test.py`, `test_layer1_layer2_correlation.py` | Two-layer decomposition (Fig 2-family: heatmap, compression, null, scatter) |
| `t3b_ellipsoid_alignment_pansci.py` | Layer-2 PanSci cross-protocol test |
| `generate_table_S6.py` | CPC1 driver genes (Table S6) |
| `12_housekeeping_ratio.py`, `13_tf_complexity.py`, `12_niche_hypothesis.py`, `12_variance_diagnostic.py`, `16_interdonor_variance.py`, `19_ppi_centrality.py`, `t3e_step2_phastcons.py`, `t3e_step2_compute.py`, `t3e_step3b_enhancer.py`, `diagnostic_expression_vs_rigidity.py` | The ten mechanistic-null tests + expression diagnostic (S1 Text §8) |
| `16_ribosomal_confound_test.py` | Ribosomal-contribution sensitivity (also cited by evolutionary-ratio line) |
| `35_l1000_random_baseline.py` | L1000 landmark-preservation; validation tooling gated by `reproduce/validate.py`, not reported in the manuscript |
| `confound_cellcount_rigidity.py` | Cell-count confound (Fig 4-family / old Fig 6B) |
| `34_samap_35types.py`, `03_samap_validation.py` | SAMap-vs-residual check: `34_samap_35types` is validation tooling gated by `reproduce/validate.py` (its figure, old Fig S5, is cut from the current paper; result not reported); `03_samap_validation` is its superseded 6-type predecessor |
| `cellmarker_35type_rerun.py`, `cellmarker_background_validated.py`, `hpa_35type_validation.py`, `52_overlay_figS6_pvalues.py`, `54_rebuild_figS6_matplotlib.py` | `cellmarker_35type_rerun` is validation tooling gated by `reproduce/validate.py` (its figure, old Fig S6, is cut from the current paper; result not reported); the rest — CellMarker background variant, HPA identity check, and the cut-Fig-S6 rebuild scripts — are not gated and not reported |
| `fetch_macaque_orthologs.py`, `nhp_ortholog_assessment.py`, `47_rerun_macaque_permutations_save_null.py`, `48_build_fig6_K12.py` | Macaque extension (S1 Text §7) |
| `41_donor_split_analysis.py`, `42_donor_split_shared_pca.py`, `43_generate_fig2e_donor_split.py` | Donor-split within-species control (Fig 3-family) |
| `08_cell_type_inventory.py` | 35-type matching inventory (Table S5) |
| `create_table_S1.py`, `create_table_S2.py`, `46_synthesis_pass_supplementary_table_edits.py`, `task_a_fix_s2_labels.py` | Table S1 / S2 build and supplementary-table edits |
| `table1_formatting.py` | Edits S13 Table (`docs/supplementary_materials/table_S13_test_inventory.xlsx`) in place and idempotently. Not a builder: S13 has no from-scratch producer in the tree, so this pass is what makes the tracked artifact reproducible from tracked code. Its lock is `TABLE_1_LOCK_MD5`, named for the table's former number |
| `test_hvg_robustness.py` | HVG-only robustness check |
| `generate_phase1_figures.py`, `generate_phase2_figures.py`, `generate_phase3_figures.py`, `composite_figS3.py`, `57_build_main_composites.py`, `build_submission_figures.py` | Figure generation / composition |
| `v1_procrustes_validation.py`, `v2_loocv_validation.py`, `v3_cellmarker_validation.py`, `verify_procrustes_vs_scipy.py` | Independent re-implementations used as reproduce cross-checks |
| `citation_audit.py`, `citation_renumber.py`, `fix_citations.py`, `convert_manuscript_to_docx.py`, `build_manuscript_pdf.py`, `assemble_supplementary_pdf.py`, `build_submission_packet.py`, `build_krt_docx.py`, `scripts/office/` (`validate.py`, `soffice.py`) | Manuscript / packet / PDF build tooling. Six entries are Cell-Systems-era, retired but retained, and none produces a tracked artifact: `build_krt_docx.py` (key-resources-table builder); `assemble_supplementary_pdf.py` (combined-SI-PDF assembler); `citation_audit.py` and `citation_renumber.py` (superscript-citation audit and renumbering, which the PLOS ONE bracketed style does not use -- the audit already parsed zero citations before the parent manuscript was retired and exited 0 while doing so); `convert_manuscript_to_docx.py` and `build_manuscript_pdf.py` (the parent render chain, parent -> .docx -> .pdf). All six read `docs/submission/manuscript_combined.txt` or its render, both retired from the repository; the live PLOS ONE DOCX is built by `docs/submission/plosone/build_manuscript_docx.py` instead |
| `06_go_enrichment.py` | GO enrichment of Procrustes residuals; TIER 1 core-pipeline step, invoked unguarded at `reproduce/run_all.sh:73`, writing the 15 tracked files in `output/phase3/go_enrichment/` via `src/cellwarp/enrichment.py` helpers (the GO panels themselves are not in the paper) |
| `12_progenitor_deep_dive.py` | Sole producer of `output/phase2/progenitor_analysis/progenitor_specificity_scores.csv`, read by `12_niche_hypothesis.py` (mechanistic null 3, S1 Text §8); not invoked by `reproduce/run_all.sh`, so the tracked CSV is the deposited copy of that input |

### OTHER PROJECT

| Script | Line |
|---|---|
| `10_cancer_download.py`, `11_cancer_procrustes.py`, `12_cancer_scaled.py`, `18_cancer_cnv_diagnostic.py`, `identity_vs_state_analysis.py` | Cancer geometry / identity-vs-state axis |
| `13_covid_procrustes.py`, `disease_inventory_covid.py`, `disease_inventory_lupus.py` | Disease / activation-state axis (COVID, lupus) |
| `20_liang_wagner_treeness.py`, `21_treeness_mediation.py`, `22_treeness_h2_concentration.py` | Treeness–rigidity line (panels cut) |
| `08_developmental_constraint.py`, `08b_developmental_followup.py` | Developmental-constraint line |
| `16_dnds_vs_expression.py`, `50_evolutionary_ratio_analysis.py` | Evolutionary-ratio (dN/dS) line |

### BANKED

| Script | Held-in-reserve analysis |
|---|---|
| `23_sensitivity_gate_l1000.py`, `24_dilirank_analysis.py` | DILI landmark preregistration (Step 1 gate + Steps 2–5). The Step-1 fractal-geometry ρ is checked by `reproduce/validate.py` but not reported in the manuscript; the toxicity claim is banked (null under its own gate). |
| `v7_tost_equivalence.py` | TOST equivalence for the SAMap-vs-residual correlation; deployable if a reviewer asks whether the null correlation is a true null. |

### EXPLORATORY — scouting, feasibility, diagnostics (no reported number)

| Script | Note |
|---|---|
| `t1a_mca_feasibility.py`, `17_sun2023_expanded.py`, `18_sun2023_issue092_diagnosis.py` | Atlas feasibility / verification checks |
| `phase3_extract_catalogs.R` | MSigDB catalog extraction for the phase3 GSEA step; the `.gmt` catalogs it writes are not redistributed (see `DATA_SOURCES.md`) |

---

## `analysis/` classification

| Directory | Category | Role |
|---|---|---|
| `conserved_contribution/` | PAPER | Master-TF / conserved-contribution (Fig 5, §2, §5, Table S11). Includes four out-of-sample controls behind S1 Text's answer to the selection-circularity objection: `block1_form_a.py` (C retention across donor halves), `block3_form_b.py` (the same split tested on the geometry rather than on C), `item1_retention.py` (the paired retention ratio formed per case), and `block2_w3_tfcensus.py` (master-TF enrichment repeated against a full TF census). `block3_form_b.py` and `block2_w3_tfcensus.py` each read one input that is not deposited; their headers say which |
| `independent_pca_sensitivity/` | PAPER | Independent-PCA sensitivity (Fig S1A–B) |
| `simulation_study/` | PAPER | Plant-and-recover simulation (Fig 4-family, §"simulation", Table S2). Includes `sweep_spread.py`, which evaluates rank recovery at the calibrated signal — absent from the deposited `RECOVERY_SIGNALS` grid — and sweeps the planted spread with re-calibration at each point; gated by `reproduce/validate.py`. Also includes `paired_signal.py` and `paired_spread.py`, which report the same recovery per replicate rather than pooled, so the signal and spread points are paired across the same drawn configurations |
| `bootstrap_rankings/` | PAPER | Bootstrap ranking CIs (Fig S3, Table S2) |
| `permutation_1M/` | PAPER | 1M headline null (Fig 1B) |
| `mantel_test/` | PAPER | Mantel test |
| `donor_split/` | PAPER | Donor-split control (Fig 3-family) |
| `expanded_negative_controls/` | PAPER | Within-species tissue-pair controls (Fig S2E) |
| `within_species_matched/` | PAPER | Matched-scale human-vs-human control (S4 Fig) |
| `census_replication/` | PAPER | Pan-Census replication (Fig S2F). Includes `item2_assay_composition.py`, which settles the assay composition of the primary human atlas that the deposited mouse-only protocol breakdown cannot; it requires the optional Census extra and the network |
| `cellhint_investigation/`, `harmonized_replication/`, `ranking_replication/` | PAPER | Cross-atlas ranking / harmonization (Tables S1, S3, S4; Fig S4). `ranking_replication/` also holds `block2_matched_n.py`, which builds the matched-n primary baseline each replication is read against; `task2_residual_mechanism.py`, which tests whether the replicated types are the more conserved ones — the mechanism Results asserts but does not quantify; and `cross_atlas_ci.py`, which puts a Bonett–Wright Fisher-z interval on each of the four cross-atlas correlations, reading them from the four separate artifacts that hold them (no upstream producer holds all four). It uses `SE = 1.06/√(n−3)`, the Spearman constant, deliberately unlike the Pearson `1/√(n−3)` in `generate_phase2_figures.py`, `t3e_step2_compute.py` and `t3e_step3b_enhancer.py`; gated by `reproduce/validate.py`, not invoked by `run_all.sh` |
| `biological_predictors/` | PAPER | Biological-predictor correlates (Table S1, §4) |
| `sensitivity_analyses/` | PAPER | Ribosomal/housekeeping exclusion, per-gene standardization, marker-null (Tables S7/S9/S10, S5 Fig) |
| `sensitivity/layer2_no_ribosomal/` | PAPER | Layer 2 under ribosomal-protein exclusion (S1 Text §4; Results §2; Methods), gated by `reproduce/validate.py` |
| `sensitivity/parent_child/` | PAPER | Parent-and-child landmark sensitivity (S1 Text §2; Results §1), gated by `reproduce/validate.py` |
| `sensitivity/mt_cellcycle/`, `sensitivity/mt_cellcycle_bilateral/` | BANKED | MT / cell-cycle sensitivity; completed but not referenced by the manuscript or `reproduce/validate.py` |
| `macaque/` | PAPER | Macaque extension (underpowered 12-type human–macaque comparison, S1 Text §7) |
| `matched_three_species/` | BANKED | Three-way ortholog intersection; NO-GO per its own `intersection_report.md` (8 types < 10-type threshold) |
| `mouse_lemur/` | PAPER | Mouse-lemur extension (§1) |
| `gene_conservation/` | PAPER | Gene-conservation-score construction (S1 Text supplementary methods) |
| `validation/` (i.e. `analysis/validation/`) | EXPLORATORY | HPA 35-type identity check (sole content `hpa/`); not read by `reproduce/validate.py` and not reported in the manuscript |
| `evolutionary_ratio/` | OTHER PROJECT | Evolutionary-ratio (dN/dS) line |
| `macaque_atlas_research/`, `third_species/` | EXPLORATORY | Macaque / third-species atlas surveys and feasibility reports |
| `cross_reference/` | PAPER | Builds `master_ranking_table.csv`, the source behind Fig 4B (CI-width vs cross-atlas rank-shift, ρ = −0.41) and Table S1 |

---

## `output/` classification

| Directory | Category | Role |
|---|---|---|
| `phase1_qc/`, `phase2/`, `phase3/`, `figures/`, `supplementary/`, `validation/` | PAPER | Core pipeline outputs, figures, and validation (this `validation/` is `output/validation/`, read by `reproduce/validate.py`) |
| `paper_audit/` | EXPLORATORY | Internal numbers audit (`master_numbers.csv`); not consumed by any tool and not reported |
| `phase1_samap/` | PAPER | SAMap-vs-residual check output; gated by `reproduce/validate.py` (its figure, old Fig S5, is cut from the current paper; not reported) |
| `mechanistic/` | PAPER | Mechanistic nulls + ellipsoid alignment (Fig 4-family) |
| `layer3_permutation/`, `twolayer_pansci_replication/` | PAPER | Two-layer decomposition + PanSci cross-protocol |
| `cellcount_confound/` | PAPER | Cell-count confound |
| `landmark_sensitivity/` | PAPER | L1000 landmark-preservation check; validation tooling gated by `reproduce/validate.py`, not reported in the manuscript |
| `macaque_pipeline/` | PAPER | Macaque extension |
| `t3g/` | PAPER | Drug-target-conservation null (mechanistic Null 10) |
| `disease_replication/` | OTHER PROJECT | Disease / COVID axis |
| `cancer/` | OTHER PROJECT | Cancer geometry |
| `liang_wagner/` | OTHER PROJECT | Treeness line |
| `ontology_design/` | OTHER PROJECT | 100-type ontology expansion (future work) |
| `dilirank/` | BANKED | DILI landmark test (null under its own gate) |

---

*Every statistic checked by `reproduce/validate.py` resolves either to a PAPER
artifact above or to the vendored basal-ganglia inputs in
`docs/submission/plosone/figures/bg_results/`, which sit outside the three trees
classified here. BANKED, OTHER PROJECT, and EXPLORATORY items make no claim in the
PLOS ONE manuscript; they are catalogued here so the reported set is legible as a
principled subset rather than a selection.*
