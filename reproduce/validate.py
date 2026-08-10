#!/usr/bin/env python3
"""
Validate CellWarp reproduction outputs against paper-reported values.

Reads output files and compares key statistics to manuscript values.
Prints a pass/fail summary with exact numbers.

Reproducibility tolerance: results reproduce to ~6 significant figures on the
same platform. Stochastic / environment-sensitive steps (Leiden cell-type
assignment; SAMap graph mapping) can vary slightly across platforms and
library versions; cell-type membership, headline results, and all scientific
conclusions are unaffected. Tolerances below are sized accordingly.
"""
import csv
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# ──────────────────────────────────────────────────────────────────
# EXPECTED VALUES FROM MANUSCRIPT
#
# Every number reported in the paper should have a corresponding check.
# Values come from the output files that produced the paper's results.
# ──────────────────────────────────────────────────────────────────

CHECKS = [
    # ── Main result (Results section 1) ──
    {
        "name": "Global coherence obs/null (35 types)",
        "file": "output/phase2/scaled_35types/procrustes_results_35.json",
        "compute": "obs_null",  # distance / null_median
        "key_distance": "procrustes.distance",
        "key_null_median": "permutation_test.null_distribution_summary.median",
        "expected": 0.522,
        "tolerance": 0.005,
        "paper_ref": "Results section 1, Figure 1",
    },
    {
        "name": "Global coherence p-value (10K run)",
        "file": "output/phase2/scaled_35types/procrustes_results_35.json",
        "key": "permutation_test.p_value",
        "expected_below": 0.001,
        "paper_ref": "Methods: replications and secondary tests use 10,000 permutations (floor p < 1e-4)",
    },
    {
        "name": "Global coherence p-value (1M primary)",
        "file": "analysis/permutation_1M/results_1M.json",
        "key": "p_value",
        "expected_below": 1e-6,
        "paper_ref": "Results section 1, Figure 1B; Methods: primary test, 1,000,000 permutations",
    },

    # ── Independent PCA (Figure S1A-B) ──
    {
        "name": "Independent PCA obs/null",
        "file": "analysis/independent_pca_sensitivity/independent_pca_results.json",
        "key": "permutation_test.obs_null_ratio",
        "expected": 0.473,
        "tolerance": 0.005,
        "paper_ref": "S1 Fig A",
    },

    # ── CellHint replication (Figure 3C) ──
    {
        "name": "CellHint obs/null",
        "file": "output/validation/cellhint_replication/cellhint_replication.json",
        "key": "procrustes.obs_null_ratio",
        "expected": 0.448,
        "tolerance": 0.005,
        "paper_ref": "Figure 3, S2 Fig(F)",
    },
    {
        "name": "CellHint p-value",
        "file": "output/validation/cellhint_replication/cellhint_replication.json",
        "key": "procrustes.p_value",
        "expected_below": 0.001,
        "paper_ref": "Figure 3",
    },

    # ── Smart-seq2 protocol sensitivity (Figure S2C-D) ──
    {
        "name": "10x-only obs/null ratio",
        "file": "output/phase2/sensitivity/smartseq2/sensitivity_results.json",
        "key": "procrustes_10x_only.obs_null_ratio",
        "expected": 0.508,
        "tolerance": 0.01,
        "paper_ref": "S2 Fig C-D",
    },

    # ── SAMap validation (validation tooling; old Fig S5, cut) ──
    {
        "name": "SAMap-residual correlation rho",
        "file": "output/phase1_samap/samap_35types/samap_rigidity_correlation.json",
        "key": "spearman_rho",
        "expected": -0.247,
        "tolerance": 0.02,
        "paper_ref": "not in current paper (SAMap validation tooling; old Fig S5, cut)",
    },

    # ── CellMarker (validation tooling; old Fig S6, cut) ──
    {
        "name": "CellMarker enrichment fold",
        "file": "output/validation/cellmarker_35type_rerun/cellmarker_35type_results.json",
        "key": "global_enrichment.enrichment",
        "expected": 4.49,
        "tolerance": 0.1,
        "paper_ref": "not in current paper (CellMarker validation tooling; old Fig S6, cut)",
    },
    {
        "name": "CellMarker enrichment p-value",
        "file": "output/validation/cellmarker_35type_rerun/cellmarker_35type_results.json",
        "key": "global_enrichment.p_value",
        # Tightened to a two-sided range that catches order-of-magnitude
        # drift in either direction. Pipeline value is ~2.10e-13.
        "expected_in_range": [1e-14, 1e-12],
        "paper_ref": "not in current paper (CellMarker validation tooling; old Fig S6, cut)",
    },
    {
        "name": "CellMarker expression-matched p-value",
        "file": "output/validation/cellmarker_35type_rerun/cellmarker_35type_results.json",
        "key": "expression_matched.p_value",
        # Pipeline value is ~1.15e-12 (V2 selection).
        "expected_in_range": [1e-13, 1e-11],
        "paper_ref": "not in current paper (CellMarker validation tooling; old Fig S6, cut)",
    },

    # ── Mechanistic null tests — housekeeping (Methods §Mechanistic null tests, row 1) ──
    {
        "name": "Housekeeping ratio vs residual Spearman rho (human)",
        "file": "output/phase2/mechanistic/housekeeping/hk_ratio_results.json",
        "key": "human_correlation.spearman_rho",
        "expected": 0.167,
        "tolerance": 0.005,
        "paper_ref": "S1 Text §8 (mechanistic null 1: housekeeping ratio)",
    },

    # ── Treeness checks removed (analysis cut) ──

    # ── PCA k-sensitivity (Figure S2A) ──
    {
        "name": "PCA sensitivity: all k values significant",
        "file": "output/validation/pca_sensitivity/pca_sensitivity_results.json",
        "key": "10.p_value",
        "expected_below": 0.001,
        "paper_ref": "S2 Fig A-B",
    },

    # ── L1000 random baseline (validation tooling; old Fig 7A, cut) ──
    {
        "name": "L1000 random baseline empirical p-value",
        "file": "output/figures/l1000_random_baseline_results.json",
        "key": "primary.empirical_p_value",
        "expected_below": 0.05,
        "paper_ref": "not in current paper (L1000 validation tooling; old Fig 7A, cut)",
    },
    {
        "name": "L1000 observed rho",
        "file": "output/figures/l1000_random_baseline_results.json",
        "key": "primary.observed_l1000_rho",
        "expected": 0.852,
        "tolerance": 0.02,
        "paper_ref": "not in current paper (L1000 validation tooling; old Fig 7A, cut)",
    },

    # ── Ellipsoid alignment (two-layer covariance; Results section 2, Figure 2) ──
    {
        "name": "Ellipsoid eigenval-residual rho (35 types)",
        "file": "output/mechanistic/ellipsoid_alignment/summary_stats.json",
        "key": "35type.eigenval_vs_rigidity.rho",
        "expected": 0.395,
        "tolerance": 0.02,
        "paper_ref": "Results section 2, Figure 2 (two-layer covariance)",
    },
    {
        "name": "Ellipsoid eigenval-residual p-value (35 types)",
        "file": "output/mechanistic/ellipsoid_alignment/summary_stats.json",
        "key": "35type.eigenval_vs_rigidity.p",
        "expected_below": 0.05,
        "paper_ref": "Results section 2, Figure 2 (two-layer covariance)",
    },

    # ── Layer-1 housekeeping/ribosomal exclusion robustness ──
    {
        "name": "Layer-1 ribosomal-excluded obs/null",
        "file": "analysis/sensitivity_analyses/layer1_exclusion_results.json",
        "key": "variants.1.obs_null_ratio",
        "expected": 0.5012,
        "tolerance": 0.01,
        "paper_ref": "Methods: ribosomal/housekeeping exclusion robustness",
    },
    {
        "name": "Layer-1 ribosomal+housekeeping-excluded obs/null",
        "file": "analysis/sensitivity_analyses/layer1_exclusion_results.json",
        "key": "variants.2.obs_null_ratio",
        "expected": 0.4788,
        "tolerance": 0.01,
        "paper_ref": "Methods: ribosomal/housekeeping exclusion robustness",
    },
    # ── Conserved-contribution gene set (Results section 5, Figure 5) ──
    {
        "name": "Conserved-contribution: master-TF median C-percentile",
        "file": "analysis/conserved_contribution/gate_results.json",
        "key": "check3a.median_Crank",
        "expected": 0.94,
        "tolerance": 0.02,
        "paper_ref": "Results section 5, Figure 5C",
    },
    {
        "name": "Conserved-contribution: master-TF positive-control p-value",
        "file": "analysis/conserved_contribution/gate_results.json",
        "key": "check3a.p_value",
        "expected_below": 0.01,
        "paper_ref": "Results section 5, Figure 5C. Note this is the 1,000-draw "
                     "run, whose floor is 1/1001; the 10^-6 bound Figure 5C "
                     "states comes from the high-N re-run gated just below",
    },

    # ── Fig 5C's stated bound: the 1,000,000-draw re-run ──
    # make_figure7.py recomputes the medians and the rug from source and states
    # "p < 10^-6 vs both" as a literal; this file is the only record of the
    # computation behind that literal, and it was ungated. The bound is checked
    # rather than the value because the p is the (0+1)/(n+1) draw floor.
    {
        "name": "Fig 5C: master TFs vs expression-matched null, p below the stated bound",
        "file": "analysis/conserved_contribution/highN_tf_pvalues.json",
        "key": "expression_matched.p_empirical",
        "expected_below": 1e-6,
        "paper_ref": "Figure 5C title (p < 10^-6, first of the two arms)",
    },
    {
        "name": "Fig 5C: master TFs vs expression+Tau-matched null, p below the stated bound",
        "file": "analysis/conserved_contribution/highN_tf_pvalues.json",
        "key": "joint_expr_tau_matched.p_empirical",
        "expected_below": 1e-6,
        "paper_ref": "Figure 5C title. This is the second arm, and it is what "
                     "makes the drawn claim 'vs both' rather than 'vs one'",
    },
    {
        "name": "Fig 5C: testable master TFs",
        "file": "analysis/conserved_contribution/highN_tf_pvalues.json",
        "key": "n_testable_tfs",
        "expected": 73,
        "tolerance": 0,
        "paper_ref": "Figure 5C title (73 master TFs). The panel counts these "
                     "itself from the ortholog table, so this is a second "
                     "record of the same number rather than its only source",
    },
    {
        "name": "Conserved-contribution: C vs expression Spearman",
        "file": "analysis/conserved_contribution/gate_results.json",
        "key": "check2.spearman_C_vs_expr",
        "expected": 0.22,
        "tolerance": 0.03,
        "paper_ref": "Results section 5, Figure 5B",
    },
    {
        "name": "Conserved-contribution: C vs specificity (Tau) Spearman",
        "file": "analysis/conserved_contribution/robustness_results.json",
        "key": "rho_C_tau",
        "expected": 0.06,
        "tolerance": 0.03,
        "paper_ref": "Results section 5, Figure 5B",
    },
    {
        "name": "Conserved-contribution: donor-split cross-half C Spearman",
        "file": "analysis/conserved_contribution/donor_stability/donor_stability_results.json",
        "key": "donor_split_cap10000.cross_half_C_spearman_median",
        "expected": 0.80,
        "tolerance": 0.03,
        "paper_ref": "Results section 5, Figure 5D",
    },
    {
        "name": "Conserved-contribution: cross-protocol C Spearman",
        "file": "analysis/conserved_contribution/donor_stability/donor_stability_results.json",
        "key": "cross_protocol.spearman_C10x_vs_CSS",
        "expected": 0.59,
        "tolerance": 0.03,
        "paper_ref": "Results section 5",
    },
    {
        "name": "Conserved-contribution: fresh-pull obs/null",
        "file": "analysis/conserved_contribution/donor_stability/donor_stability_results.json",
        "key": "validity.obs_null_full",
        "expected": 0.521,
        "tolerance": 0.01,
        "paper_ref": "Results section 5 (fresh Census re-acquisition)",
    },

    # ── Conserved-contribution geometry attribution (Results section 2, §10) ──
    {
        "name": "Conserved-contribution: conserved-set obs/null (geometry attribution)",
        "file": "analysis/conserved_contribution/gate_results.json",
        "key": "secondary.conserved.ratio",
        "expected": 0.384,
        "tolerance": 0.005,
        "paper_ref": "Results section 5 (geometry attribution); gate_results.json[secondary]",
    },
    {
        "name": "Conserved-contribution: divergent-set obs/null (geometry attribution)",
        "file": "analysis/conserved_contribution/gate_results.json",
        "key": "secondary.divergent.ratio",
        "expected": 0.709,
        "tolerance": 0.01,
        "paper_ref": "Results section 5 (geometry attribution); gate_results.json[secondary]",
    },
    {
        "name": "Conserved-contribution: expr-matched-random obs/null (geometry attribution)",
        "file": "analysis/conserved_contribution/gate_results.json",
        "key": "secondary.matched_random_ratio_mean",
        "expected": 0.525,
        "tolerance": 0.02,
        "paper_ref": "Results section 5 (geometry attribution; expr-matched random, 0.525 +/- 0.012)",
    },
    {
        "name": "Conserved-contribution: all-genes obs/null anchor (geometry attribution)",
        "file": "analysis/conserved_contribution/gate_results.json",
        "key": "secondary.validity_all_genes_ratio",
        "expected": 0.522,
        "tolerance": 0.005,
        "paper_ref": "Results section 5 (geometry attribution; all-genes anchor)",
    },
    {
        "name": "Conserved-contribution: Hartigan dip statistic D (broad continuum)",
        "file": "analysis/conserved_contribution/robustness_results.json",
        "key": "dip_D",
        "expected": 0.007,
        "tolerance": 0.001,
        "paper_ref": "Results section 5, Figure 5A (dip D = 0.007, p = 2.8e-5)",
    },

    # ── S1 Text §10: selection/derangement circularity control ──
    # Gates the two pre-specified PASS conditions by their numeric backing, since the
    # conditions themselves are stored as booleans: condition 1 is "real below the
    # sigma-null 1st percentile" (the real value 0.384 is already gated above from
    # gate_results.json; the threshold is gated here), condition 2 is "z <= -3".
    # Tolerances are set to catch a broken control, not to pin this host: a control that
    # stopped destroying the correspondence would collapse the sigma-null from ~0.99
    # toward the full-space 0.52, which is ~47x the 0.01 band used here.
    {
        "name": "Selection null: derangement sigma-null mean",
        "file": "analysis/selection_null/outputs/selection_null_summary_derangement.json",
        "key": "sigma_null.mean",
        "expected": 0.991,
        "tolerance": 0.01,
        "paper_ref": "S1 Text §10 (derangement sigma-null 0.991 ± 0.021)",
    },
    {
        "name": "Selection null: derangement sigma-null 1st percentile",
        "file": "analysis/selection_null/outputs/selection_null_summary_derangement.json",
        "key": "sigma_null.p01",
        "expected": 0.927,
        "tolerance": 0.01,
        "paper_ref": "S1 Text §10 (pre-specified PASS condition 1: real below the 1st percentile)",
    },
    {
        "name": "Selection null: derangement z (real vs sigma-null)",
        "file": "analysis/selection_null/outputs/selection_null_summary_derangement.json",
        "key": "real_position.z",
        "expected": -29.5,
        "tolerance": 0.5,
        "paper_ref": "S1 Text §10 (pre-specified PASS condition 2: z <= -3; reported z = -29.5)",
    },
    {
        "name": "Selection null: derangement draws at or below real",
        "file": "analysis/selection_null/outputs/selection_null_summary_derangement.json",
        "key": "real_position.n_draws_at_or_below_real",
        "expected": 0,
        "tolerance": 0,
        "paper_ref": "S1 Text §10 (0 of 1,000 draws at or below the real value)",
    },
    {
        "name": "Selection null: label-shuffle sigma-null mean",
        "file": "analysis/selection_null/outputs/selection_null_summary_labelshuffle.json",
        "key": "sigma_null.mean",
        "expected": 0.983,
        "tolerance": 0.01,
        "paper_ref": "S1 Text §10 (label-shuffle cross-check sigma-null 0.983 ± 0.024)",
    },
    {
        "name": "Selection null: label-shuffle z (real vs sigma-null)",
        "file": "analysis/selection_null/outputs/selection_null_summary_labelshuffle.json",
        "key": "real_position.z",
        "expected": -25.3,
        "tolerance": 0.5,
        "paper_ref": "S1 Text §10 (label-shuffle cross-check z = -25.3)",
    },

    # ── Fig 2C: basal-ganglia two-layer replication, self-consistency (current paper) ──
    # Structural checks over the vendored BG results; no manuscript values transcribed.
    # (The conserved-contribution checks above are Fig 5C, not this; "Figure 2C" here is the current basal-ganglia panel.)
    {
        "name": "BG Human-Macaque: post-rotation compression (k=5, W2)",
        "file": "docs/submission/plosone/figures/bg_results/layer2_results_Human_Macaque.json",
        "key": "W2_schemeB.layer2.k5.compression_ratio_post_over_pre",
        "expected_below": 1.0,
        "paper_ref": "Results section 2, Figure 2C",
    },
    {
        "name": "BG Human-Macaque: post-rotation permutation p (k=5, W2)",
        "file": "docs/submission/plosone/figures/bg_results/layer2_results_Human_Macaque.json",
        "key": "W2_schemeB.layer2.k5.p_post",
        "expected_below": 0.00011,
        "paper_ref": "Results section 2, Figure 2C",
    },
    {
        "name": "BG Human-Macaque: identity-marker driver count (W2)",
        "file": "docs/submission/plosone/figures/bg_results/layer2_results_Human_Macaque.json",
        "key": "W2_schemeB.rank1_class_counts.canonical identity marker",
        "expected": 18,
        "tolerance": 0,
        "paper_ref": "Results section 2, Figure 2C (18/55)",
    },
    {
        "name": "BG Human-Macaque: n_types",
        "file": "docs/submission/plosone/figures/bg_results/layer2_results_Human_Macaque.json",
        "key": "n_types",
        "expected": 55,
        "tolerance": 0,
        "paper_ref": "Results section 2, Figure 2C (18/55)",
    },
    # ── The Macaca-restricted arms of the same replication (Methods; S1 Text) ──
    # The macaque arm pools two Macaca species at genus level. Restricting it to the
    # six M. mulatta donors and recomputing raises these fractions; both are stated in
    # Methods and S1 Text, so both are gated as the pooled ones are, numerator and
    # denominator. The producers and the full outputs live in the basal-ganglia record;
    # these two results files are vendored here exactly as the pooled three are.
    {
        "name": "BG Human-Macaque (M. mulatta restricted): identity-marker driver count (W2)",
        "file": "docs/submission/plosone/figures/bg_results/layer2_results_Human_Macaque_mulatta.json",
        "key": "W2_schemeB.rank1_class_counts.canonical identity marker",
        "expected": 25,
        "tolerance": 0,
        "paper_ref": "Methods and S1 Text (25 of 54)",
    },
    {
        "name": "BG Human-Macaque (M. mulatta restricted): n_types",
        "file": "docs/submission/plosone/figures/bg_results/layer2_results_Human_Macaque_mulatta.json",
        "key": "n_types",
        "expected": 54,
        "tolerance": 0,
        "paper_ref": "Methods and S1 Text (25 of 54; one sparsely sampled type falls "
                     "below the 100-cell inclusion threshold under restriction, so the "
                     "denominator is 54 against the pooled 55)",
    },
    {
        "name": "BG Macaque-Marmoset (M. mulatta restricted): identity-marker driver count (W2)",
        "file": "docs/submission/plosone/figures/bg_results/layer2_results_Macaque_Marmoset_mulatta.json",
        "key": "W2_schemeB.rank1_class_counts.canonical identity marker",
        "expected": 8,
        "tolerance": 0,
        "paper_ref": "Methods and S1 Text (8 of 51)",
    },
    {
        "name": "BG Macaque-Marmoset (M. mulatta restricted): n_types",
        "file": "docs/submission/plosone/figures/bg_results/layer2_results_Macaque_Marmoset_mulatta.json",
        "key": "n_types",
        "expected": 51,
        "tolerance": 0,
        "paper_ref": "Methods and S1 Text (8 of 51; same threshold effect, 51 against "
                     "the pooled 52)",
    },

    {
        "name": "BG Human-Marmoset: post-rotation compression (k=5, W2)",
        "file": "docs/submission/plosone/figures/bg_results/layer2_results_Human_Marmoset.json",
        "key": "W2_schemeB.layer2.k5.compression_ratio_post_over_pre",
        "expected_below": 1.0,
        "paper_ref": "Results section 2, Figure 2C",
    },
    {
        "name": "BG Human-Marmoset: post-rotation permutation p (k=5, W2)",
        "file": "docs/submission/plosone/figures/bg_results/layer2_results_Human_Marmoset.json",
        "key": "W2_schemeB.layer2.k5.p_post",
        "expected_below": 0.00011,
        "paper_ref": "Results section 2, Figure 2C",
    },
    {
        "name": "BG Human-Marmoset: identity-marker driver count (W2)",
        "file": "docs/submission/plosone/figures/bg_results/layer2_results_Human_Marmoset.json",
        "key": "W2_schemeB.rank1_class_counts.canonical identity marker",
        "expected": 7,
        "tolerance": 0,
        "paper_ref": "Results section 2, Figure 2C (7/52)",
    },
    {
        "name": "BG Human-Marmoset: n_types",
        "file": "docs/submission/plosone/figures/bg_results/layer2_results_Human_Marmoset.json",
        "key": "n_types",
        "expected": 52,
        "tolerance": 0,
        "paper_ref": "Results section 2, Figure 2C (7/52)",
    },
    {
        "name": "BG Macaque-Marmoset: post-rotation compression (k=5, W2)",
        "file": "docs/submission/plosone/figures/bg_results/layer2_results_Macaque_Marmoset.json",
        "key": "W2_schemeB.layer2.k5.compression_ratio_post_over_pre",
        "expected_below": 1.0,
        "paper_ref": "Results section 2, Figure 2C",
    },
    {
        "name": "BG Macaque-Marmoset: post-rotation permutation p (k=5, W2)",
        "file": "docs/submission/plosone/figures/bg_results/layer2_results_Macaque_Marmoset.json",
        "key": "W2_schemeB.layer2.k5.p_post",
        "expected_below": 0.00011,
        "paper_ref": "Results section 2, Figure 2C",
    },
    {
        "name": "BG Macaque-Marmoset: identity-marker driver count (W2)",
        "file": "docs/submission/plosone/figures/bg_results/layer2_results_Macaque_Marmoset.json",
        "key": "W2_schemeB.rank1_class_counts.canonical identity marker",
        "expected": 5,
        "tolerance": 0,
        "paper_ref": "Results section 2, Figure 2C (5/52)",
    },
    {
        "name": "BG Macaque-Marmoset: n_types",
        "file": "docs/submission/plosone/figures/bg_results/layer2_results_Macaque_Marmoset.json",
        "key": "n_types",
        "expected": 52,
        "tolerance": 0,
        "paper_ref": "Results section 2, Figure 2C (5/52)",
    },

    # ── Simulation study: rank-recovery ceiling (Fig 4C, Methods) ──
    # The deposited RECOVERY_SIGNALS grid does not contain the calibrated signal,
    # so the grid values and the calibrated-signal value are gated separately.
    # Both are deterministic, hence the tight tolerances.
    {
        "name": "Simulation: calibrated signal strength",
        "file": "analysis/simulation_study/simulation_results.json",
        "key": "calibration.estimated_real_signal",
        "expected": 3.683,
        "tolerance": 0.001,
        "paper_ref": "Methods, Simulation study (calibrated signal ~ 3.68)",
    },
    {
        "name": "Simulation: grid median recovery rho, signal 3.0, 200 cells",
        "file": "analysis/simulation_study/simulation_results.json",
        "key": "ranking_recovery.17.median_rho",
        "expected": 0.4224,
        "tolerance": 0.001,
        "paper_ref": "Fig 4C (plotted curve, signal 3.0); Methods, Simulation study",
    },
    {
        "name": "Simulation: grid median recovery rho, signal 5.0, 200 cells",
        "file": "analysis/simulation_study/simulation_results.json",
        "key": "ranking_recovery.21.median_rho",
        "expected": 0.4293,
        "tolerance": 0.001,
        "paper_ref": "Fig 4C (plotted curve, signal 5.0); Methods, Simulation study",
    },
    {
        "name": "Simulation sweep: median recovery rho at the calibrated signal, 200 cells",
        "file": "analysis/simulation_study/sweep_spread_results.json",
        "key": "sweep.3.recovery.1.median_rho",
        "expected": 0.4494,
        "tolerance": 0.001,
        "paper_ref": "Results section 8 and Methods, Simulation study: rank-recovery "
                     "ceiling at the calibrated signal (deposited spread, sigma = 1.0)",
    },
    {
        "name": "Simulation sweep: spread-range upper endpoint (25x spread, 200 cells)",
        "file": "analysis/simulation_study/sweep_spread_results.json",
        "key": "sweep.2.recovery.1.median_rho",
        "expected": 0.4955,
        "tolerance": 0.001,
        "paper_ref": "Results section 8 and Methods, Simulation study: upper endpoint of the "
                     "planted-spread ceiling range (lower endpoint is the calibrated-signal check above)",
    },
    {
        "name": "Simulation sweep: zero-spread negative control near zero",
        "file": "analysis/simulation_study/sweep_spread_results.json",
        "key": "sweep.0.recovery.1.median_rho",
        "expected_below": 0.05,
        "paper_ref": "No manuscript value; guards the recovery metric. With no planted "
                     "divergence there is no ranking to recover (observed 0.0059)",
    },

    # ── Parent-and-child landmark sensitivity (S1 Text §2; Results §1) ──
    # Ontology parent labels sit alongside their own child types among the 35
    # landmarks. Both reduced sets clear the permutation floor. n_types is gated
    # beside each ratio because obs/null is not comparable across landmark counts,
    # so the ratio alone would gate half a claim.
    {
        "name": "Parent-child variant A (drop parents): obs/null at 30 types",
        "file": "analysis/sensitivity/parent_child/results.json",
        "key": "variant_A_drop_parents.obs_null_ratio",
        "expected": 0.5210,
        "tolerance": 0.001,
        "paper_ref": "S1 Text §2 (obs/null = 0.521 on the remaining 30 types); "
                     "Results §1 (0.52 at 30 types)",
    },
    {
        "name": "Parent-child variant A: landmark count",
        "file": "analysis/sensitivity/parent_child/results.json",
        "key": "variant_A_drop_parents.n_types",
        "expected": 30,
        "tolerance": 0,
        "paper_ref": "S1 Text §2 (dropping the five parent labels whose children "
                     "are also present leaves 30 types)",
    },
    {
        "name": "Parent-child variant A: permutation p at the floor",
        "file": "analysis/sensitivity/parent_child/results.json",
        "key": "variant_A_drop_parents.p_value_permutation",
        "expected_below": 0.00011,
        "paper_ref": "S1 Text §2 and Results §1 (each p < 1e-4; 10,000 permutations)",
    },
    {
        "name": "Parent-child variant A: ranking rho vs the primary subset",
        "file": "analysis/sensitivity/parent_child/results.json",
        "key": "variant_A_drop_parents.ranking_spearman_vs_primary_subset.rho",
        "expected": 0.9359,
        "tolerance": 0.001,
        "paper_ref": "S1 Text §2 (ranking Spearman 0.936 against the full set)",
    },
    {
        "name": "Parent-child variant B (drop children): obs/null at 26 types",
        "file": "analysis/sensitivity/parent_child/results.json",
        "key": "variant_B_drop_children.obs_null_ratio",
        "expected": 0.5441,
        "tolerance": 0.001,
        "paper_ref": "S1 Text §2 (obs/null = 0.544 on 26 types); Results §1 (0.54 at 26)",
    },
    {
        "name": "Parent-child variant B: landmark count",
        "file": "analysis/sensitivity/parent_child/results.json",
        "key": "variant_B_drop_children.n_types",
        "expected": 26,
        "tolerance": 0,
        "paper_ref": "S1 Text §2 (dropping the nine children instead leaves 26 types)",
    },
    {
        "name": "Parent-child variant B: permutation p at the floor",
        "file": "analysis/sensitivity/parent_child/results.json",
        "key": "variant_B_drop_children.p_value_permutation",
        "expected_below": 0.00011,
        "paper_ref": "S1 Text §2 and Results §1 (each p < 1e-4; 10,000 permutations)",
    },
    {
        "name": "Parent-child variant B: ranking rho vs the primary subset",
        "file": "analysis/sensitivity/parent_child/results.json",
        "key": "variant_B_drop_children.ranking_spearman_vs_primary_subset.rho",
        "expected": 0.9111,
        "tolerance": 0.001,
        "paper_ref": "S1 Text §2 (ranking Spearman 0.911 against the full set)",
    },

    # ── Marker-similarity-stratified null (S1 Text §2; Results §1; S10 Table, S5 Fig) ──
    # S10 Table publishes the 100,000-permutation columns, so the sweep and its
    # p-values are gated on obs_null_100k / p_100k. The K=1 anchor is the one
    # exception and reads the 10k field deliberately: see its own paper_ref.
    {
        "name": "Marker-null K=1 degenerate anchor recovers the primary (negative control)",
        "file": "analysis/sensitivity_analyses/markernull_results.json",
        "key": "anchors.k1.obs_null_10k",
        "expected": 0.5222043226858066,
        "tolerance": 1e-9,
        "paper_ref": "Negative control, not a manuscript value. At K=1 the single group "
                     "holds all 35 types, so the restricted null IS the unrestricted null "
                     "and this must reproduce the primary bit-for-bit. Reads the 10k field "
                     "because the primary was run at 10,000 permutations; the 100k value "
                     "(0.5222677) is correctly different. The tight tolerance is the point: "
                     "a loose one would gate nothing here",
    },
    {
        "name": "Marker-null Ward K=5: obs/null",
        "file": "analysis/sensitivity_analyses/markernull_results.json",
        "key": "ward_sweep.5.obs_null_100k",
        "expected": 0.7208,
        "tolerance": 0.001,
        "paper_ref": "Results §1 (0.72 at five groups); S1 Text §2 (0.72 at K = 5); S10 Table",
    },
    {
        "name": "Marker-null Ward K=8: obs/null",
        "file": "analysis/sensitivity_analyses/markernull_results.json",
        "key": "ward_sweep.8.obs_null_100k",
        "expected": 0.8074,
        "tolerance": 0.001,
        "paper_ref": "S1 Text §2 (0.81 at K = 8); S10 Table",
    },
    {
        "name": "Marker-null Ward K=10: obs/null",
        "file": "analysis/sensitivity_analyses/markernull_results.json",
        "key": "ward_sweep.10.obs_null_100k",
        "expected": 0.9033,
        "tolerance": 0.001,
        "paper_ref": "Results §1 (0.90 at ten); S1 Text §2 (0.90 at K = 10); S10 Table",
    },
    {
        "name": "Marker-null Ward K=15: obs/null",
        "file": "analysis/sensitivity_analyses/markernull_results.json",
        "key": "ward_sweep.15.obs_null_100k",
        "expected": 0.9631,
        "tolerance": 0.001,
        "paper_ref": "S1 Text §2 (0.96 at K = 15); S10 Table",
    },
    {
        "name": "Marker-null Ward K=14: still clears p = 0.05",
        "file": "analysis/sensitivity_analyses/markernull_results.json",
        "key": "ward_sweep.14.p_100k",
        "expected_below": 0.05,
        "paper_ref": "S1 Text §2 (the observed configuration outperforms the finer null "
                     "through K = 14). Paired with the K=15 check below: together they "
                     "gate the crossover from both sides, so neither can pass vacuously",
    },
    {
        "name": "Marker-null Ward K=15: no longer clears p = 0.05",
        "file": "analysis/sensitivity_analyses/markernull_results.json",
        "key": "ward_sweep.15.p_100k",
        "expected_in_range": [0.05, 1.0],
        "paper_ref": "Results §1 (not significant at fifteen); S1 Text §2 (no longer clears "
                     "p = 0.05 at K = 15). The one check in this file that must FAIL to pass",
    },
    {
        "name": "Marker-null K=15: non-singleton group count",
        "file": "analysis/sensitivity_analyses/markernull_results.json",
        "key": "ward_sweep.15.n_nonsingleton",
        "expected": 9,
        "tolerance": 0,
        "paper_ref": "S1 Text §2 (the nine non-singleton groups at K = 15; the six "
                     "singletons S1 also cites are 15 minus this count); S10 Table",
    },
    {
        "name": "Marker-null K=15: log10 within-group permutation space",
        "file": "analysis/sensitivity_analyses/markernull_results.json",
        "key": "ward_sweep.15.log10_permspace",
        "expected": 9.0771,
        "tolerance": 0.001,
        "paper_ref": "S1 Text §2 (about 1.2 billion within-group permutations against the "
                     "100,000 drawn; 10^9.0771 = 1.194e9). This is the evidence that the "
                     "K=15 null is not exhausted, so the loss of signal there is real",
    },
    {
        "name": "Marker-null: per-type residual vs marker-distinctness rho",
        "file": "analysis/sensitivity_analyses/markernull_results.json",
        "key": "monotonicity.spearman_rho",
        "expected": -0.1361,
        "tolerance": 0.001,
        "paper_ref": "S1 Text §2 (the per-type residual is independent of "
                     "marker-distinctness, Spearman rho = -0.14). NOTE the key name: this "
                     "block is the n = 35 per-type correlation, not a monotonicity test "
                     "over K. Its stored p implies df = 33, which is the check on that",
    },
    {
        "name": "Marker-null: per-type residual vs marker-distinctness p",
        "file": "analysis/sensitivity_analyses/markernull_results.json",
        "key": "monotonicity.spearman_p",
        "expected": 0.4355,
        "tolerance": 0.001,
        "paper_ref": "S1 Text §2 (p = 0.44). See the note on the preceding check",
    },

    # ── Layer 2: aggregate Krzanowski S, observed (Results §2, Fig 2B, Methods) ──
    # These six are the panel's own numbers and the paper's second headline result.
    # Nothing read them before this group: validate.py opened summary_stats.json
    # only for eigenval_vs_rigidity, and never opened permutation_results.json.
    {
        "name": "Layer 2 (35 types): observed S pre-rotation, k=1",
        "file": "output/mechanistic/ellipsoid_alignment/summary_stats.json",
        "key": "35type.mean_alignment.k=1.pre",
        "expected": 0.1291,
        "tolerance": 0.001,
        "paper_ref": "Methods, two-layer decomposition (0.129 before rotation at k = 1); Fig 2B",
    },
    {
        "name": "Layer 2 (35 types): observed S post-rotation, k=1",
        "file": "output/mechanistic/ellipsoid_alignment/summary_stats.json",
        "key": "35type.mean_alignment.k=1.post",
        "expected": 0.0256,
        "tolerance": 0.001,
        "paper_ref": "Results §2 and Methods (0.026 after rotation at k = 1); Fig 2B",
    },
    {
        "name": "Layer 2 (35 types): observed S pre-rotation, k=3",
        "file": "output/mechanistic/ellipsoid_alignment/summary_stats.json",
        "key": "35type.mean_alignment.k=3.pre",
        "expected": 0.3850,
        "tolerance": 0.001,
        "paper_ref": "Fig 2B bar label (0.385)",
    },
    {
        "name": "Layer 2 (35 types): observed S post-rotation, k=3",
        "file": "output/mechanistic/ellipsoid_alignment/summary_stats.json",
        "key": "35type.mean_alignment.k=3.post",
        "expected": 0.1671,
        "tolerance": 0.001,
        "paper_ref": "Fig 2B bar label (0.167)",
    },
    {
        "name": "Layer 2 (35 types): observed S pre-rotation, k=5",
        "file": "output/mechanistic/ellipsoid_alignment/summary_stats.json",
        "key": "35type.mean_alignment.k=5.pre",
        "expected": 0.4826,
        "tolerance": 0.001,
        "paper_ref": "Abstract, Results §2, Fig 2 caption, Methods (S = 0.483 at k = 5)",
    },
    {
        "name": "Layer 2 (35 types): observed S post-rotation, k=5",
        "file": "output/mechanistic/ellipsoid_alignment/summary_stats.json",
        "key": "35type.mean_alignment.k=5.post",
        "expected": 0.2295,
        "tolerance": 0.001,
        "paper_ref": "Abstract, Results §2, Fig 2 caption, Methods (S drops to 0.230 at k = 5)",
    },

    # ── Layer 2: permutation nulls and p-values (Results §2, Fig 2 caption, Methods) ──
    # The k=1 pair is a two-sided control in the same sense as the marker-null
    # K=14/K=15 pair: the paper claims Layer 2 is significant at k=3 and k=5 and
    # NOT at k=1, so one check requires significance and its partner requires the
    # absence of it. The 0.119246 null is the value that was written into Methods
    # as 0.115 from memory and corrected in cf089eb5; nothing has gated it until now.
    {
        "name": "Layer 2: permutation null mean, pre-rotation k=1",
        "file": "output/mechanistic/ellipsoid_alignment/permutation_results.json",
        "key": "35type.label_shuffle_pre.k=1.null_mean",
        "expected": 0.1192,
        "tolerance": 0.001,
        "paper_ref": "Methods (0.129 against a null of 0.119 before rotation)",
    },
    {
        "name": "Layer 2: pre-rotation k=1 is NOT significant",
        "file": "output/mechanistic/ellipsoid_alignment/permutation_results.json",
        "key": "35type.label_shuffle_pre.k=1.p_value",
        "expected_in_range": [0.05, 1.0],
        "paper_ref": "Methods (significant at k = 5 and k = 3, but not k = 1). Must FAIL to "
                     "pass; partner to the k=3 and k=5 floor checks below",
    },
    {
        "name": "Layer 2: pre-rotation k=3 permutation p at the floor",
        "file": "output/mechanistic/ellipsoid_alignment/permutation_results.json",
        "key": "35type.label_shuffle_pre.k=3.p_value",
        "expected_below": 0.00011,
        "paper_ref": "Results §2 (the agreement is significant at k = 3 and k = 5)",
    },
    {
        "name": "Layer 2: permutation null mean, pre-rotation k=5",
        "file": "output/mechanistic/ellipsoid_alignment/permutation_results.json",
        "key": "35type.label_shuffle_pre.k=5.null_mean",
        "expected": 0.3749,
        "tolerance": 0.001,
        "paper_ref": "Fig 2 caption and Methods (S = 0.483 vs null 0.375)",
    },
    {
        "name": "Layer 2: pre-rotation k=5 permutation p at the floor",
        "file": "output/mechanistic/ellipsoid_alignment/permutation_results.json",
        "key": "35type.label_shuffle_pre.k=5.p_value",
        "expected_below": 0.00011,
        "paper_ref": "Methods (Layer 2 significant at k = 5, p < 1e-4)",
    },
    {
        "name": "Layer 2: permutation null mean, post-rotation k=1",
        "file": "output/mechanistic/ellipsoid_alignment/permutation_results.json",
        "key": "35type.label_shuffle_post.k=1.null_mean",
        "expected": 0.0296,
        "tolerance": 0.001,
        "paper_ref": "Results §2 and Methods (0.026 against 0.030 after rotation)",
    },
    {
        "name": "Layer 2: post-rotation k=1 p-value",
        "file": "output/mechanistic/ellipsoid_alignment/permutation_results.json",
        "key": "35type.label_shuffle_post.k=1.p_value",
        "expected": 0.7651,
        "tolerance": 0.001,
        "paper_ref": "Results §2, Fig 2 caption and Methods (p = 0.77; the post-rotation "
                     "point estimate falls just below its null)",
    },
    {
        "name": "Layer 2: post-rotation k=3 permutation p at the floor",
        "file": "output/mechanistic/ellipsoid_alignment/permutation_results.json",
        "key": "35type.label_shuffle_post.k=3.p_value",
        "expected_below": 0.00011,
        "paper_ref": "Results §2 (significant at k = 3 and k = 5 after rotation)",
    },
    {
        "name": "Layer 2: permutation null mean, post-rotation k=5",
        "file": "output/mechanistic/ellipsoid_alignment/permutation_results.json",
        "key": "35type.label_shuffle_post.k=5.null_mean",
        "expected": 0.1801,
        "tolerance": 0.001,
        "paper_ref": "Results §2, Fig 2 caption and Methods (well above its null, 0.180)",
    },
    {
        "name": "Layer 2: post-rotation k=5 permutation p at the floor",
        "file": "output/mechanistic/ellipsoid_alignment/permutation_results.json",
        "key": "35type.label_shuffle_post.k=5.p_value",
        "expected_below": 0.00011,
        "paper_ref": "Results §2, Fig 2 caption and Methods (p < 1e-4 after rotation)",
    },
    {
        "name": "Layer 2: ratio to null, pre-rotation k=5",
        "file": "output/mechanistic/ellipsoid_alignment/permutation_results.json",
        "compute": "obs_null",
        "key_distance": "35type.label_shuffle_pre.k=5.observed",
        "key_null_median": "35type.label_shuffle_pre.k=5.null_mean",
        "expected": 1.2873,
        "tolerance": 0.001,
        "paper_ref": "Results §2 and Fig 2 caption (1.29 before rotation). NOTE the harness "
                     "field is named key_null_median but the paper's ratio is to the null "
                     "MEAN, which is what is resolved here. This is the quantity that was "
                     "written as 1.28 by dividing the rounded display values",
    },
    {
        "name": "Layer 2: ratio to null, post-rotation k=5",
        "file": "output/mechanistic/ellipsoid_alignment/permutation_results.json",
        "compute": "obs_null",
        "key_distance": "35type.label_shuffle_post.k=5.observed",
        "key_null_median": "35type.label_shuffle_post.k=5.null_mean",
        "expected": 1.2746,
        "tolerance": 0.001,
        "paper_ref": "Results §2 and Fig 2 caption (1.27 after rotation). See the note on "
                     "the preceding check",
    },

    # ── Lineage-stratified null (Results §1, Fig 1B/1C, Methods, S1 Text §2) ──
    {
        "name": "Lineage-stratified null: obs/null",
        "file": "output/validation/lineage_stratified/lineage_stratified_results.json",
        "key": "stratified_null.obs_null_ratio",
        "expected": 0.6683,
        "tolerance": 0.001,
        "paper_ref": "Results §1 and Methods (obs/null 0.668); Fig 1C caption (0.67); "
                     "S1 Text §2; Table 1 T02",
    },
    {
        "name": "Lineage-stratified null: permutation p at the floor",
        "file": "output/validation/lineage_stratified/lineage_stratified_results.json",
        "key": "stratified_null.p_value",
        "expected_below": 0.00011,
        "paper_ref": "Results §1, Fig 1C caption, Methods, S1 Text §2 (p < 1e-4)",
    },
    {
        "name": "Lineage-stratified null: how much tighter than the global null",
        "file": "output/validation/lineage_stratified/lineage_stratified_results.json",
        "key": "distribution_comparison.pct_tighter",
        "expected": 21.8614,
        "tolerance": 0.001,
        "paper_ref": "S1 Text §2 (the within-lineage null is 21.9 percent tighter than the "
                     "global null)",
    },
    {
        "name": "Observed human-mouse Procrustes distance",
        "file": "output/validation/lineage_stratified/lineage_stratified_results.json",
        "key": "observed_procrustes_distance",
        "expected": 61.1530,
        "tolerance": 0.001,
        "paper_ref": "Fig 1B caption (observed distance 61.15) and Methods (observed "
                     "distance 61.15; nearest null 98.88)",
    },
    {
        "name": "Canonical headline obs/null, full precision (cross-artifact anchor)",
        "file": "output/validation/lineage_stratified/lineage_stratified_results.json",
        "key": "global_null.obs_null_ratio",
        "expected": 0.5222043226858066,
        "tolerance": 1e-9,
        "paper_ref": "Results §1 (0.52). This file, markernull_results.json's K=1 anchor and "
                     "the matched-n baselines producer all carry this identical value, so it "
                     "is the canonical form. t1a_results.json's tabula_35_obs_null is "
                     "0.5222043484, 2.6e-8 away, and must NOT be used for the 35-type value",
    },

    # ── Mouse-lemur replication (Results §1, Fig 1D; Table 1 T13) ──
    # Fig 1D is drawn by docs/submission/plosone/figures/build_main_figures.py,
    # which reads only null_distribution.npy from this directory and states the
    # rest as text. These four checks put the stated values under Gate 1 so the
    # panel can read them instead.
    {
        "name": "Mouse lemur: obs/null",
        "file": "analysis/mouse_lemur/procrustes_results.json",
        "key": "permutation_test.obs_null_ratio",
        "expected": 0.3463,
        "tolerance": 0.001,
        "paper_ref": "Results §1 and Fig 1D caption (obs/null 0.35); Table 1 T13",
    },
    {
        "name": "Mouse lemur: permutation p below the stated bound",
        "file": "analysis/mouse_lemur/procrustes_results.json",
        "key": "permutation_test.p_value",
        "expected_below": 0.0001,
        "paper_ref": "Results §1 and Fig 1D caption (p < 10^-4). Phrased as a bound "
                     "rather than a value because the stored p is the (0+1)/(n+1) "
                     "permutation floor, which is strictly below 1e-4, and the "
                     "claim the paper makes is the bound",
    },
    {
        "name": "Mouse lemur: matched cell types",
        "file": "analysis/mouse_lemur/procrustes_results.json",
        "key": "n_types",
        "expected": 15,
        "tolerance": 0,
        "paper_ref": "Fig 1D caption (n = 15)",
    },
    {
        "name": "Mouse lemur: divergence time",
        "file": "analysis/mouse_lemur/procrustes_results.json",
        "key": "divergence_mya",
        "expected": 75,
        "tolerance": 0,
        "paper_ref": "Fig 1D panel title (~75 Mya), Fig 1 caption (~75 Myr) and "
                     "Results §1 (roughly 75 million years). NUMBER_DIFF.md C1 "
                     "records a draft that said 90 and was corrected against this "
                     "file, so the number is gated rather than left to drift back",
    },

    # ── Planted-spread description (Methods, Simulation study) ──
    # Deterministic: measure_range() draws under RandomState(0) over 20,000 draws.
    {
        "name": "Simulation: planted spread, median max/min ratio",
        "file": "analysis/simulation_study/sweep_spread_results.json",
        "key": "deposited.max_min_median",
        "expected": 65.0130,
        "tolerance": 0.001,
        "paper_ref": "Methods, Simulation study (a 65-fold ratio between the largest and "
                     "smallest planted divergence)",
    },
    {
        "name": "Simulation: planted spread, median 95th/5th percentile ratio",
        "file": "analysis/simulation_study/sweep_spread_results.json",
        "key": "deposited.p95_p5_median",
        "expected": 20.5899,
        "tolerance": 0.001,
        "paper_ref": "Methods, Simulation study (20.6-fold between the 5th and 95th "
                     "percentiles)",
    },

    # ── Matched-scale 6-type controls (Results §4, S4 Fig; Table 1 T16/T17) ──
    # The two arms live in separate files with identical structure. The 35-type
    # value in t1a_results.json is deliberately not used; see the anchor above.
    {
        "name": "Matched-scale: human-mouse 6-type obs/null",
        "file": "output/phase2/procrustes_results.json",
        "compute": "obs_null",
        "key_distance": "procrustes.distance",
        "key_null_median": "permutation_test.null_distribution_summary.median",
        "expected": 0.3166,
        "tolerance": 0.001,
        "paper_ref": "Results §4 and S4 Fig caption (obs/null = 0.317); S13 Table T17",
    },
    {
        "name": "Matched-scale: human-mouse 6-type permutation p",
        "file": "output/phase2/procrustes_results.json",
        "key": "permutation_test.p_value",
        "expected": 0.0035,
        "tolerance": 0.001,
        "paper_ref": "S4 Fig caption (p = 0.0035); S13 Table T17",
    },
    {
        "name": "Matched-scale: human-human 6-type obs/null",
        "file": "output/phase2/negative_control_v2/negctrl_v2_results.json",
        "compute": "obs_null",
        "key_distance": "procrustes.distance",
        "key_null_median": "permutation_test.null_distribution_summary.median",
        "expected": 0.6066,
        "tolerance": 0.001,
        "paper_ref": "Results §4, S4 Fig caption and S1 Text §9 (obs/null = 0.607); "
                     "Table 1 T16",
    },
    {
        "name": "Matched-scale: human-human 6-type permutation p",
        "file": "output/phase2/negative_control_v2/negctrl_v2_results.json",
        "key": "permutation_test.p_value",
        "expected": 0.0088,
        "tolerance": 0.001,
        "paper_ref": "S4 Fig caption (p = 0.0088); S13 Table T16",
    },

    # ── Lineage blocks: the strata the lineage-stratified null permutes within ──
    # S1 Text §2 states these five sizes as facts and rests the conservativeness
    # argument on the two singletons. They are gated on their own account rather
    # than through the ratio above: a re-run that changed a block would move
    # obs_null_ratio, but metadata edited WITHOUT a re-run - a rename, a recount -
    # would leave S1 wrong while every gated statistic sat still.
    {
        "name": "Lineage block size: immune/hematopoietic",
        "file": "output/validation/lineage_stratified/lineage_stratified_results.json",
        "key": "lineage_blocks.immune_hematopoietic.n_types",
        "expected": 18,
        "tolerance": 0,
        "paper_ref": "S1 Text §2 (immune/hematopoietic 18 types); S1 Text §2 also rests the "
                     "discriminative-power claim on this being the largest block",
    },
    {
        "name": "Lineage block size: epithelial",
        "file": "output/validation/lineage_stratified/lineage_stratified_results.json",
        "key": "lineage_blocks.epithelial.n_types",
        "expected": 8,
        "tolerance": 0,
        "paper_ref": "S1 Text §2 (epithelial 8)",
    },
    {
        "name": "Lineage block size: stromal/mesenchymal",
        "file": "output/validation/lineage_stratified/lineage_stratified_results.json",
        "key": "lineage_blocks.stromal_mesenchymal.n_types",
        "expected": 7,
        "tolerance": 0,
        "paper_ref": "S1 Text §2 (stromal/mesenchymal 7)",
    },
    {
        "name": "Lineage block size: endothelial (singleton)",
        "file": "output/validation/lineage_stratified/lineage_stratified_results.json",
        "key": "lineage_blocks.endothelial.n_types",
        "expected": 1,
        "tolerance": 0,
        "paper_ref": "S1 Text §2 (endothelial 1). The singleton count is load-bearing: "
                     "singletons stay matched under all permutations, which is why S1 "
                     "calls the stratified null more conservative",
    },
    {
        "name": "Lineage block size: metabolic/parenchymal (singleton)",
        "file": "output/validation/lineage_stratified/lineage_stratified_results.json",
        "key": "lineage_blocks.metabolic_parenchymal.n_types",
        "expected": 1,
        "tolerance": 0,
        "paper_ref": "S1 Text §2 (metabolic/parenchymal 1). See the note on the preceding "
                     "check",
    },
    {
        "name": "Lineage-stratified null: landmark count",
        "file": "output/validation/lineage_stratified/lineage_stratified_results.json",
        "key": "n_cell_types",
        "expected": 35,
        "tolerance": 0,
        "paper_ref": "Results §1 and throughout (35 matched cell types). Gated here beside "
                     "the block sizes, which must sum to it",
    },
    {
        "name": "Lineage-stratified null: permutation count",
        "file": "output/validation/lineage_stratified/lineage_stratified_results.json",
        "key": "n_permutations",
        "expected": 10000,
        "tolerance": 0,
        "paper_ref": "S1 Text §2 (10,000 iterations); Methods (replications, extensions and "
                     "secondary tests used 10,000 permutations, floor p < 1e-4). The floor "
                     "the p-value check above compares against is a function of this",
    },

    # ── Layer 2 under ribosomal-protein exclusion (S1 Text §4; Results §2; Methods) ──
    # The producer is analysis/sensitivity/layer2_no_ribosomal/run.py. SCOPE.md still
    # classifies this directory as BANKED and "not referenced by the manuscript" --
    # that row is false and its rewrite is queued separately.
    #
    # DO NOT gate anything from this file's primary_comparison block: it holds
    # ROUNDED copies of the primary Layer-2 values (0.483 / 0.375 / 0.23 / 0.18) and
    # two p-values stored as the string '<1e-4'. A check against a rounded copy passes
    # at any tolerance the rounding permits and detects nothing.
    {
        "name": "Ribosomal-excluded Layer 2: pre-rotation S, k=5",
        "file": "analysis/sensitivity/layer2_no_ribosomal/results.json",
        "key": "krzanowski_pre_rotation.5.observed_S",
        "expected": 0.4455,
        "tolerance": 0.001,
        "paper_ref": "S1 Text §4 (after ribosomal exclusion, pre-rotation S = 0.446)",
    },
    {
        "name": "Ribosomal-excluded Layer 2: pre-rotation null mean, k=5",
        "file": "analysis/sensitivity/layer2_no_ribosomal/results.json",
        "key": "krzanowski_pre_rotation.5.null_mean",
        "expected": 0.3287,
        "tolerance": 0.001,
        "paper_ref": "S1 Text §4 (null 0.329)",
    },
    {
        "name": "Ribosomal-excluded Layer 2: pre-rotation p at the floor, k=5",
        "file": "analysis/sensitivity/layer2_no_ribosomal/results.json",
        "key": "krzanowski_pre_rotation.5.p_value",
        "expected_below": 0.00011,
        "paper_ref": "S1 Text §4 (both p < 1e-4)",
    },
    {
        "name": "Ribosomal-excluded Layer 2: post-rotation S, k=5",
        "file": "analysis/sensitivity/layer2_no_ribosomal/results.json",
        "key": "krzanowski_post_rotation.5.observed_S",
        "expected": 0.2345,
        "tolerance": 0.001,
        "paper_ref": "S1 Text §4 and Methods (excluding ribosomal-protein genes entirely "
                     "leaves the compression significant, post-rotation S = 0.234)",
    },
    {
        "name": "Ribosomal-excluded Layer 2: post-rotation null mean, k=5",
        "file": "analysis/sensitivity/layer2_no_ribosomal/results.json",
        "key": "krzanowski_post_rotation.5.null_mean",
        "expected": 0.1800,
        "tolerance": 0.001,
        "paper_ref": "S1 Text §4 and Methods (vs null 0.180)",
    },
    {
        "name": "Ribosomal-excluded Layer 2: post-rotation p at the floor, k=5",
        "file": "analysis/sensitivity/layer2_no_ribosomal/results.json",
        "key": "krzanowski_post_rotation.5.p_value",
        "expected_below": 0.00011,
        "paper_ref": "S1 Text §4 and Methods (p < 1e-4)",
    },

    # ── The k=1 claim: exclusion moves it ABOVE its null (Results §2) ──
    # Results §2 says the post-rotation k=1 value sits slightly BELOW its null in the
    # primary (0.026 against 0.030), and that excluding ribosomal-protein genes moves
    # k=1 above its null. Four checks gate that reversal by its two point estimates
    # and their two nulls, so the claim is covered whichever rotation it refers to.
    # The two p-value checks that follow are the two-sided partner: the reversal is a
    # point-estimate claim, NOT a significance claim, and the paper does not make one.
    {
        "name": "Ribosomal-excluded Layer 2: pre-rotation S, k=1",
        "file": "analysis/sensitivity/layer2_no_ribosomal/results.json",
        "key": "krzanowski_pre_rotation.1.observed_S",
        "expected": 0.0636,
        "tolerance": 0.001,
        "paper_ref": "Results §2 (excluding ribosomal-protein genes moves k = 1 above its "
                     "null); this is the observed side of that comparison, pre-rotation",
    },
    {
        "name": "Ribosomal-excluded Layer 2: pre-rotation null mean, k=1",
        "file": "analysis/sensitivity/layer2_no_ribosomal/results.json",
        "key": "krzanowski_pre_rotation.1.null_mean",
        "expected": 0.0490,
        "tolerance": 0.001,
        "paper_ref": "Results §2; the null side of the same comparison. Observed 0.0636 "
                     "exceeds this, which is the claim",
    },
    {
        "name": "Ribosomal-excluded Layer 2: post-rotation S, k=1",
        "file": "analysis/sensitivity/layer2_no_ribosomal/results.json",
        "key": "krzanowski_post_rotation.1.observed_S",
        "expected": 0.0487,
        "tolerance": 0.001,
        "paper_ref": "Results §2; post-rotation side. In the primary this value sits BELOW "
                     "its null (0.026 against 0.030); here it sits above",
    },
    {
        "name": "Ribosomal-excluded Layer 2: post-rotation null mean, k=1",
        "file": "analysis/sensitivity/layer2_no_ribosomal/results.json",
        "key": "krzanowski_post_rotation.1.null_mean",
        "expected": 0.0393,
        "tolerance": 0.001,
        "paper_ref": "Results §2; the null side. Observed 0.0487 exceeds this, completing "
                     "the reversal the sentence asserts",
    },
    {
        "name": "Ribosomal-excluded Layer 2: pre-rotation k=1 is NOT significant",
        "file": "analysis/sensitivity/layer2_no_ribosomal/results.json",
        "key": "krzanowski_pre_rotation.1.p_value",
        "expected_in_range": [0.05, 1.0],
        "paper_ref": "No manuscript value. Results §2 claims only that k = 1 moves above its "
                     "null, not that it becomes significant. This check exists so that a "
                     "future run which made it significant would fail rather than quietly "
                     "strengthen a claim the paper does not make. Must FAIL to pass",
    },
    {
        "name": "Ribosomal-excluded Layer 2: post-rotation k=1 is NOT significant",
        "file": "analysis/sensitivity/layer2_no_ribosomal/results.json",
        "key": "krzanowski_post_rotation.1.p_value",
        "expected_in_range": [0.05, 1.0],
        "paper_ref": "See the preceding check. Must FAIL to pass",
    },

    # ── Structural control on the exclusion itself ──
    {
        "name": "Ribosomal exclusion actually removed the ribosomal drivers",
        "file": "analysis/sensitivity/layer2_no_ribosomal/results.json",
        "key": "cpc1_summary.n_new_top1_still_ribosomal",
        "expected": 0,
        "tolerance": 0,
        "paper_ref": "No manuscript value; a structural control. After excluding the "
                     "ribosomal-protein genes, no type can retain a ribosomal rank-1 CPC1 "
                     "driver. Trivially true when the filter ran, and false when it did "
                     "not, which is exactly what makes it worth checking",
    },

    # ── Cross-atlas rank-correlation confidence intervals (Results §4) ──
    # Producer: analysis/ranking_replication/cross_atlas_ci.py. Bonett-Wright
    # Fisher-z intervals, SE = 1.06/sqrt(n-3), on the four correlations Results §4
    # reports. The claim being gated is that every interval spans both zero and
    # 0.20, i.e. these samples do not separate "no cross-atlas ranking signal"
    # from "a moderate one". Each pair carries four checks: the two booleans that
    # ARE the claim, plus rho and n, so an upstream artifact cannot drift beneath
    # an interval that still reads as correct.
    #
    # The booleans are stored as JSON booleans and load_value() coerces with
    # float(), so true arrives as 1.0 and false as 0.0; expected 1 / tolerance 0
    # therefore fails loudly the moment an interval stops spanning its target.
    # This follows the S1 Text §10 precedent, where pre-registered conditions are
    # likewise stored as booleans and gated numerically.
    {
        "name": "Cross-atlas CI (Sun2023): rho as read from the source artifact",
        "file": "analysis/ranking_replication/cross_atlas_ci_results.json",
        "key": "pairs.sun2023.rho",
        "expected": 0.1464,
        "tolerance": 0.0001,
        "paper_ref": "Results §4 and S1 Text §6 (+0.15 Sun2023); input to the interval below",
    },
    {
        "name": "Cross-atlas CI (Sun2023): n as read from the source artifact",
        "file": "analysis/ranking_replication/cross_atlas_ci_results.json",
        "key": "pairs.sun2023.n",
        "expected": 15,
        "tolerance": 0,
        "paper_ref": "Results §4 (n = 15 to 22 across the four pairs); the SE is a "
                     "function of n, so n is gated alongside rho",
    },
    {
        "name": "Cross-atlas CI (Sun2023): interval contains zero",
        "file": "analysis/ranking_replication/cross_atlas_ci_results.json",
        "key": "pairs.sun2023.ci_contains_zero",
        "expected": 1,
        "tolerance": 0,
        "paper_ref": "Results §4: the 95% interval spans zero (computed [-0.424, +0.633])",
    },
    {
        "name": "Cross-atlas CI (Sun2023): interval contains 0.20",
        "file": "analysis/ranking_replication/cross_atlas_ci_results.json",
        "key": "pairs.sun2023.ci_contains_0_20",
        "expected": 1,
        "tolerance": 0,
        "paper_ref": "Results §4: the same interval also spans a moderate correlation "
                     "of 0.20, which is what makes the pair uninformative rather than negative",
    },
    {
        "name": "Cross-atlas CI (PanSci): rho as read from the source artifact",
        "file": "analysis/ranking_replication/cross_atlas_ci_results.json",
        "key": "pairs.pansci.rho",
        "expected": 0.1941,
        "tolerance": 0.0001,
        "paper_ref": "Results §4 and S1 Text §6 (+0.19 PanSci); input to the interval below",
    },
    {
        "name": "Cross-atlas CI (PanSci): n as read from the source artifact",
        "file": "analysis/ranking_replication/cross_atlas_ci_results.json",
        "key": "pairs.pansci.n",
        "expected": 16,
        "tolerance": 0,
        "paper_ref": "Results §4. PanSci's matched-type count is 16, not 15; it is the "
                     "one pair whose n differs from the 15 the other two 15-type arms carry",
    },
    {
        "name": "Cross-atlas CI (PanSci): interval contains zero",
        "file": "analysis/ranking_replication/cross_atlas_ci_results.json",
        "key": "pairs.pansci.ci_contains_zero",
        "expected": 1,
        "tolerance": 0,
        "paper_ref": "Results §4: the 95% interval spans zero (computed [-0.362, +0.649])",
    },
    {
        "name": "Cross-atlas CI (PanSci): interval contains 0.20",
        "file": "analysis/ranking_replication/cross_atlas_ci_results.json",
        "key": "pairs.pansci.ci_contains_0_20",
        "expected": 1,
        "tolerance": 0,
        "paper_ref": "Results §4: the same interval also spans 0.20. This is the largest "
                     "of the four correlations, so it is the pair most likely to lose the "
                     "zero end if the estimate moved",
    },
    {
        "name": "Cross-atlas CI (pan-Census): rho as read from the source artifact",
        "file": "analysis/ranking_replication/cross_atlas_ci_results.json",
        "key": "pairs.pan_census.rho",
        "expected": -0.0525,
        "tolerance": 0.0001,
        "paper_ref": "Results §4 and S1 Text §6 (-0.05 pan-Census); input to the interval below",
    },
    {
        "name": "Cross-atlas CI (pan-Census): n as read from the source artifact",
        "file": "analysis/ranking_replication/cross_atlas_ci_results.json",
        "key": "pairs.pan_census.n",
        "expected": 22,
        "tolerance": 0,
        "paper_ref": "Results §4: the upper end of the n = 15 to 22 range, so the "
                     "narrowest of the four intervals",
    },
    {
        "name": "Cross-atlas CI (pan-Census): interval contains zero",
        "file": "analysis/ranking_replication/cross_atlas_ci_results.json",
        "key": "pairs.pan_census.ci_contains_zero",
        "expected": 1,
        "tolerance": 0,
        "paper_ref": "Results §4: the 95% interval spans zero (computed [-0.485, +0.400])",
    },
    {
        "name": "Cross-atlas CI (pan-Census): interval contains 0.20",
        "file": "analysis/ranking_replication/cross_atlas_ci_results.json",
        "key": "pairs.pan_census.ci_contains_0_20",
        "expected": 1,
        "tolerance": 0,
        "paper_ref": "Results §4: the same interval also spans 0.20. This is the "
                     "narrowest interval of the four, so it is the binding case for the "
                     "sentence: if any pair fails to reach 0.20 it is this one",
    },
    {
        "name": "Cross-atlas CI (CellHint): rho is the matched 15-type baseline, not the artifact",
        "file": "analysis/ranking_replication/cross_atlas_ci_results.json",
        "key": "pairs.cellhint.rho",
        "expected": -0.1393,
        "tolerance": 0.0001,
        "paper_ref": "Results §4 and S1 Text §6 (-0.14 CellHint at matched 15-type "
                     "baseline). Deliberately NOT -0.386: that is the pre-PCA-matching "
                     "artifact S1 Text explains away, and it is what "
                     "cellhint_replication.json carries as rigidity_ranking.rho. This "
                     "check is what stops the wrong row being substituted",
    },
    {
        "name": "Cross-atlas CI (CellHint): n as read from the source artifact",
        "file": "analysis/ranking_replication/cross_atlas_ci_results.json",
        "key": "pairs.cellhint.n",
        "expected": 15,
        "tolerance": 0,
        "paper_ref": "Results §4; S4 Table Level 0 (matched 15-type PCA baseline)",
    },
    {
        "name": "Cross-atlas CI (CellHint): interval contains zero",
        "file": "analysis/ranking_replication/cross_atlas_ci_results.json",
        "key": "pairs.cellhint.ci_contains_zero",
        "expected": 1,
        "tolerance": 0,
        "paper_ref": "Results §4: the 95% interval spans zero (computed [-0.629, +0.430])",
    },
    {
        "name": "Cross-atlas CI (CellHint): interval contains 0.20",
        "file": "analysis/ranking_replication/cross_atlas_ci_results.json",
        "key": "pairs.cellhint.ci_contains_0_20",
        "expected": 1,
        "tolerance": 0,
        "paper_ref": "Results §4: the same interval also spans 0.20, despite the point "
                     "estimate being negative",
    },

    # ── Donor-split within-species control (S1 Text §9) ──
    # Producer: scripts/42_donor_split_shared_pca.py, the shared-PCA variant.
    # S1 Text §9 states these four numbers in submitted text and nothing read them
    # until now.
    #
    # The artifact is donor_split_shared_pca_results.json and NOT
    # analysis/conserved_contribution/donor_stability/donor_stability_results.json.
    # Both are called "the donor split" and both say "100 splits", and they are
    # different objects: §9 partitions the 24 Tabula Sapiens donors into halves and
    # compares human-half-1 vs human-half-2 against human-half-1 vs mouse, while
    # §12's donor_stability run partitions both species' donors for the per-gene
    # conservation score C. donor_stability_results.json carries no delta at all --
    # its donor-split key is cross_half_C_spearman_median. S1 Text §9 warns about
    # exactly this collision in its own text.
    #
    # Script 41 is the superseded independent-PCA predecessor whose delta is the
    # +0.158 sensitivity value §9 also quotes; 42 is the reported one.
    {
        "name": "Donor split (S1 Text §9): median delta, cross-species minus within-species",
        "file": "analysis/donor_split/donor_split_shared_pca_results.json",
        "key": "delta.median",
        "expected": 0.1588,
        "tolerance": 0.0001,
        "paper_ref": "S1 Text §9 (median delta = +0.159); shared-PCA donor split, 100 splits, seed 42",
    },
    {
        "name": "Donor split (S1 Text §9): delta 95% CI lower bound",
        "file": "analysis/donor_split/donor_split_shared_pca_results.json",
        "key": "delta.ci_95.0",
        "expected": 0.1003,
        "tolerance": 0.0001,
        "paper_ref": "S1 Text §9 (95% CI +0.100 to +0.218). The lower bound is the "
                     "load-bearing end: the hierarchy claim needs it above zero",
    },
    {
        "name": "Donor split (S1 Text §9): delta 95% CI upper bound",
        "file": "analysis/donor_split/donor_split_shared_pca_results.json",
        "key": "delta.ci_95.1",
        "expected": 0.2180,
        "tolerance": 0.0001,
        "paper_ref": "S1 Text §9 (95% CI +0.100 to +0.218)",
    },
    {
        "name": "Donor split (S1 Text §9): cross-species exceeded within-species in every split",
        "file": "analysis/donor_split/donor_split_shared_pca_results.json",
        "key": "delta.pct_positive",
        "expected": 1,
        "tolerance": 0,
        "paper_ref": "S1 Text §9 (100/100 splits). §9 reads this as direction-robustness "
                     "across non-independent resamples of one 24-donor atlas, not as "
                     "independent replication",
    },
    {
        "name": "Donor split (S1 Text §9): number of random balanced splits",
        "file": "analysis/donor_split/donor_split_shared_pca_results.json",
        "key": "n_splits",
        "expected": 100,
        "tolerance": 0,
        "paper_ref": "S1 Text §9 (100 random balanced splits, seed 42). Gated alongside "
                     "pct_positive because 100/100 is a ratio and both halves must hold",
    },

    # ── Detection-breadth sensitivity of the master-TF enrichment ──
    # Producer: analysis/conserved_contribution/breadth_sensitivity.py.
    # C is a Pearson correlation across 35 centroids, and the pipeline imposes no
    # detection-breadth requirement on the genes entering it: gate_lib.per_gene_corr
    # filters on np.std > 0 per species alone, so a gene detected in three types and
    # absent from thirty-two is admitted. The producer supplies a breadth criterion for
    # sensitivity only (centroid > 0, the most generous reading) and re-runs the Fig 5C
    # enrichment on genes detected in all 35 types in both species.
    #
    # These five gate the strictest filter, which is the row the manuscript will state.
    # The producer additionally asserts, before computing any filtered row, that its
    # unfiltered row reproduces gate_results.json check3a.median_Crank exactly and both
    # deposited null medians within sampler tolerance -- without that the filtered
    # numbers would be uninterpretable, because a difference could be the filter or
    # could be a reimplementation of the statistic.
    {
        "name": "Breadth sensitivity: master-TF median C-percentile, all-35 filter",
        "file": "analysis/conserved_contribution/breadth_sensitivity_results.json",
        "key": "enrichment.strictest.obs_median_C_percentile",
        "expected": 0.9213,
        "tolerance": 0.0001,
        "paper_ref": "Results §5 / Fig 5C sensitivity: the enrichment survives restricting "
                     "to genes detected in all 35 types in both species (unfiltered 0.9377)",
    },
    {
        "name": "Breadth sensitivity: expression-matched null median, all-35 filter",
        "file": "analysis/conserved_contribution/breadth_sensitivity_results.json",
        "key": "enrichment.strictest.expression_matched.null_median",
        "expected": 0.5049,
        "tolerance": 0.0005,
        "paper_ref": "Results §5 / Fig 5C sensitivity (unfiltered null median 0.54). "
                     "Tolerance is looser than the observed value's because this is a "
                     "sampled null, not a deterministic statistic",
    },
    {
        "name": "Breadth sensitivity: expr x Tau-matched null median, all-35 filter",
        "file": "analysis/conserved_contribution/breadth_sensitivity_results.json",
        "key": "enrichment.strictest.joint_expr_tau_matched.null_median",
        "expected": 0.7149,
        "tolerance": 0.0005,
        "paper_ref": "Results §5 / Fig 5C sensitivity (unfiltered null median 0.76). The "
                     "stricter of the two nulls, so the binding comparison",
    },
    {
        "name": "Breadth sensitivity: genes retained by the all-35 filter",
        "file": "analysis/conserved_contribution/breadth_sensitivity_results.json",
        "key": "enrichment.strictest.n_genes",
        "expected": 8435,
        "tolerance": 0,
        "paper_ref": "Results §5 / Fig 5C sensitivity: 8,435 of 15,940 genes, i.e. the "
                     "filter discards 47 percent of the pool",
    },
    {
        "name": "Breadth sensitivity: master TFs retained by the all-35 filter",
        "file": "analysis/conserved_contribution/breadth_sensitivity_results.json",
        "key": "enrichment.strictest.n_tfs",
        "expected": 24,
        "tolerance": 0,
        "paper_ref": "Results §5 / Fig 5C sensitivity: 24 of the 73 master TFs survive. "
                     "Gated because the enrichment holding on 24 TFs is a weaker claim "
                     "than on 73, and the count is what makes that legible",
    },

    # ── Ontology restriction of the primate covariance replication ──
    # Producer: analysis/bg/ontology_restriction.py in the basal-ganglia deposit;
    # the results JSON is vendored here byte-identically, as the five
    # layer2_results_*.json already are.
    #
    # Cross-species type correspondence in that atlas comes from a consensus taxonomy
    # built by cross-species integration, so the matching BETWEEN types is not
    # independent of an alignment step the way the centroids are. The producer restricts
    # each pair to the types an independent Cell Ontology lookup corroborates and re-runs
    # the permutation null inside the restricted subset. It asserts all 36 deposited
    # means against np.diag() before computing any restricted value, and it neither edits
    # layer2_analyze.py nor regenerates any layer2_results_*.json.
    #
    # Numbers gated here are the ones Methods and S1 Text state. All are counts, so the
    # tolerance is 0 throughout: each is a statement about how many cells behaved a
    # certain way, and "approximately 34 of 36" would mean nothing.
    {
        "name": "Ontology restriction: restricted S above its own null in every cell",
        "file": "docs/submission/plosone/figures/bg_results/ontology_restriction_results.json",
        "key": "summary.n_restricted_above_own_null",
        "expected": 36,
        "tolerance": 0,
        "paper_ref": "Methods (primate replication) and S1 Text: all 36 cells -- 3 pairs x "
                     "2 weightings x k in {1,3,5} x pre/post -- stay above a null built on "
                     "the restricted subset. This is the claim the Methods sentence makes",
    },
    {
        "name": "Ontology restriction: restricted margin wider than full-set margin",
        "file": "docs/submission/plosone/figures/bg_results/ontology_restriction_results.json",
        "key": "summary.n_restricted_margin_wider_than_full",
        "expected": 34,
        "tolerance": 0,
        "paper_ref": "S1 Text: 34 of 36. The two exceptions narrow by 0.003 and 0.004, so "
                     "the count and not an average is what carries the claim",
    },
    {
        "name": "Ontology restriction: k=1 direction, non-corroborated higher",
        "file": "docs/submission/plosone/figures/bg_results/ontology_restriction_results.json",
        "key": "summary.direction_counts_by_k.1.non_corroborated_higher",
        "expected": 11,
        "tolerance": 0,
        "paper_ref": "S1 Text: at k = 1 the non-corroborated stratum carries higher S in 11 "
                     "of 12 pair-weighting-arm combinations",
    },
    {
        "name": "Ontology restriction: k=5 direction, corroborated higher",
        "file": "docs/submission/plosone/figures/bg_results/ontology_restriction_results.json",
        "key": "summary.direction_counts_by_k.5.corroborated_higher",
        "expected": 11,
        "tolerance": 0,
        "paper_ref": "S1 Text: at k = 5 the corroborated stratum carries higher S in 11 of "
                     "12. Gated alongside the k=1 count because the reversal between them "
                     "is the finding, and one count alone does not show it",
    },
    {
        "name": "Ontology restriction: Human-Macaque strict strata (corroborated)",
        "file": "docs/submission/plosone/figures/bg_results/ontology_restriction_results.json",
        "key": "strata.Human_Macaque.strict_corroborated",
        "expected": 24,
        "tolerance": 0,
        "paper_ref": "S1 Text: 24 of 55 under the strict reading",
    },
    {
        "name": "Ontology restriction: Human-Macaque strict strata (not corroborated)",
        "file": "docs/submission/plosone/figures/bg_results/ontology_restriction_results.json",
        "key": "strata.Human_Macaque.strict_not",
        "expected": 31,
        "tolerance": 0,
        "paper_ref": "S1 Text: the complement, 31 of 55. Both halves are gated because the "
                     "Methods sentence says the restriction discards more than half",
    },
    {
        "name": "Ontology restriction: Human-Marmoset strict strata (corroborated)",
        "file": "docs/submission/plosone/figures/bg_results/ontology_restriction_results.json",
        "key": "strata.Human_Marmoset.strict_corroborated",
        "expected": 20,
        "tolerance": 0,
        "paper_ref": "S1 Text: 20 of 52 under the strict reading",
    },
    {
        "name": "Ontology restriction: Human-Marmoset strict strata (not corroborated)",
        "file": "docs/submission/plosone/figures/bg_results/ontology_restriction_results.json",
        "key": "strata.Human_Marmoset.strict_not",
        "expected": 32,
        "tolerance": 0,
        "paper_ref": "S1 Text: the complement, 32 of 52",
    },
    {
        "name": "Ontology restriction: Macaque-Marmoset strict strata (corroborated)",
        "file": "docs/submission/plosone/figures/bg_results/ontology_restriction_results.json",
        "key": "strata.Macaque_Marmoset.strict_corroborated",
        "expected": 20,
        "tolerance": 0,
        "paper_ref": "S1 Text: 20 of 52 under the strict reading",
    },
    {
        "name": "Ontology restriction: Macaque-Marmoset strict strata (not corroborated)",
        "file": "docs/submission/plosone/figures/bg_results/ontology_restriction_results.json",
        "key": "strata.Macaque_Marmoset.strict_not",
        "expected": 32,
        "tolerance": 0,
        "paper_ref": "S1 Text: the complement, 32 of 52",
    },

    # ── Per-gene standardization: the four numbers the abstract's claim rests on ──
    # Producer: analysis/sensitivity_analyses/genestd_standardization.py. Its JSON
    # already carried all four; only the checks were missing. Scheme A is a z-score
    # across the 70 centroids, Scheme B a z-score across cells per species
    # (pp.scale-style) -- the S9 Table caption defines both, and the attribution
    # matters because the abstract's "moves off housekeeping structure toward
    # cell-identity markers" is a Scheme B claim, the 1-of-35 arm, not Scheme A's 0.
    {
        "name": "Per-gene standardization Scheme A: Layer-1 obs/null",
        "file": "analysis/sensitivity_analyses/genestd_results.json",
        "key": "layer1.A.obs_null",
        "expected": 0.606,
        "tolerance": 0.001,
        "paper_ref": "S9 Table caption (obs/null 0.522 -> 0.606 [A]); Scheme A is the "
                     "z-score across the 70 centroids",
    },
    {
        "name": "Per-gene standardization Scheme B: Layer-1 obs/null",
        "file": "analysis/sensitivity_analyses/genestd_results.json",
        "key": "layer1.B.obs_null",
        "expected": 0.487,
        "tolerance": 0.001,
        "paper_ref": "S9 Table caption (0.487 [B]); Scheme B is the per-species z-score "
                     "across cells",
    },
    {
        "name": "Per-gene standardization Scheme A: ribosomal-dominated CPC1 types",
        "file": "analysis/sensitivity_analyses/genestd_results.json",
        "key": "cpc1.A.n_ribosomal_dominated",
        "expected": 0,
        "tolerance": 0,
        "paper_ref": "S9 Table caption (25/35 -> 0/35 [A]). Gated with the Scheme B count "
                     "because the pair is the claim; either alone is uninterpretable",
    },
    {
        "name": "Per-gene standardization Scheme B: ribosomal-dominated CPC1 types",
        "file": "analysis/sensitivity_analyses/genestd_results.json",
        "key": "cpc1.B.n_ribosomal_dominated",
        "expected": 1,
        "tolerance": 0,
        "paper_ref": "S9 Table caption (1/35 [B]). This is the arm the abstract's "
                     "identity-marker claim rests on",
    },

    # ── Reported parameters that no artifact carried as a gateable scalar ──
    # Producer: analysis/reported_parameters/reported_parameters.py. Nothing new is
    # computed there: each field is read or reduced from a tracked file, because these
    # values lived in a CSV column (validate.py's .csv branch is unreachable), a Python
    # literal, or a per-type list that has to be reduced before it is a number.
    {
        "name": "Primate replication: minimum per-type cells, Human-Macaque",
        "file": "analysis/reported_parameters/reported_parameters.json",
        "key": "primate_replication.pairs.Human_Macaque.min_cells_either_arm",
        "expected": 108,
        "tolerance": 0,
        "paper_ref": "Results §2 says '52 to 55 matched cell types, a few hundred to over "
                     "100,000 cells each'. This is the actual floor, and the producer "
                     "additionally asserts no type falls below the declared threshold of 100",
    },
    {
        "name": "Primate replication: minimum per-type cells, Human-Marmoset",
        "file": "analysis/reported_parameters/reported_parameters.json",
        "key": "primate_replication.pairs.Human_Marmoset.min_cells_either_arm",
        "expected": 108,
        "tolerance": 0,
        "paper_ref": "Results §2; same floor as Human-Macaque, the same cell type",
    },
    {
        "name": "Primate replication: minimum per-type cells, Macaque-Marmoset",
        "file": "analysis/reported_parameters/reported_parameters.json",
        "key": "primate_replication.pairs.Macaque_Marmoset.min_cells_either_arm",
        "expected": 110,
        "tolerance": 0,
        "paper_ref": "Results §2",
    },
    {
        "name": "Primate replication: maximum per-type cells, Human-Macaque",
        "file": "analysis/reported_parameters/reported_parameters.json",
        "key": "primate_replication.pairs.Human_Macaque.max_cells_either_arm",
        "expected": 145491,
        "tolerance": 0,
        "paper_ref": "Results §2, the 'over 100,000' end of the range",
    },
    {
        "name": "Primate replication: maximum per-type cells, Macaque-Marmoset",
        "file": "analysis/reported_parameters/reported_parameters.json",
        "key": "primate_replication.pairs.Macaque_Marmoset.max_cells_either_arm",
        "expected": 83269,
        "tolerance": 0,
        "paper_ref": "Results §2. This pair's maximum is BELOW 100,000, so the range "
                     "sentence does not hold pair by pair",
    },
    {
        "name": "Mouse-lemur: ortholog pairs from BioMart",
        "file": "analysis/reported_parameters/reported_parameters.json",
        "key": "mouse_lemur.n_ortholog_pairs",
        "expected": 16655,
        "tolerance": 0,
        "paper_ref": "Methods (mouse-lemur extension), about to be expanded. 16,655 pairs "
                     "map to the 13,796-gene space gated separately",
    },
    {
        "name": "Mouse-lemur: per-type cell threshold",
        "file": "analysis/reported_parameters/reported_parameters.json",
        "key": "mouse_lemur.min_cells_per_type",
        "expected": 500,
        "tolerance": 0,
        "paper_ref": "Methods (mouse-lemur extension). Read from the producer's module "
                     "constant, not transcribed",
    },
    {
        "name": "Mouse-lemur: per-type cell cap",
        "file": "analysis/reported_parameters/reported_parameters.json",
        "key": "mouse_lemur.max_cells_per_type",
        "expected": 2000,
        "tolerance": 0,
        "paper_ref": "Methods (mouse-lemur extension); the same 2,000 cap the primary uses",
    },
    {
        "name": "Mouse-lemur: retained PCA components",
        "file": "analysis/mouse_lemur/procrustes_results.json",
        "key": "pca.n_components",
        "expected": 15,
        "tolerance": 0,
        "paper_ref": "Methods (mouse-lemur extension). Equal to the matched-type count, "
                     "which is why it is gated: 15 centred points span at most 14 "
                     "dimensions, so this is the saturated end of the design",
    },
    {
        "name": "Mouse-lemur: shared ortholog gene space",
        "file": "analysis/mouse_lemur/procrustes_results.json",
        "key": "gene_space",
        "expected": 13796,
        "tolerance": 0,
        "paper_ref": "Methods (mouse-lemur extension); the 16,655 BioMart pairs reduce to "
                     "this many genes present in both matrices",
    },
    {
        "name": "Bootstrap rank stability: types meeting the CI-width criterion",
        "file": "analysis/reported_parameters/reported_parameters.json",
        "key": "bootstrap_rank_stability.n_stable",
        "expected": 35,
        "tolerance": 0,
        "paper_ref": "S1 Text and Fig 4A: 'All 35 types classified stable (95% CI width "
                     "<= 10)'. The criterion is 10 of 35 ranks",
    },
    {
        "name": "Bootstrap rank stability: maximum CI width",
        "file": "analysis/reported_parameters/reported_parameters.json",
        "key": "bootstrap_rank_stability.ci_width_max",
        "expected": 7,
        "tolerance": 0,
        "paper_ref": "S1 Text: 'median width 3, maximum 7'. The maximum is the load-bearing "
                     "one: it is what puts every type inside the criterion",
    },
    {
        "name": "Bootstrap rank stability: median CI width",
        "file": "analysis/reported_parameters/reported_parameters.json",
        "key": "bootstrap_rank_stability.ci_width_median",
        "expected": 3,
        "tolerance": 0,
        "paper_ref": "S1 Text ('median width 3') and Fig 4A caption ('median 95% confidence "
                     "interval = 3 ranks')",
    },
]


def _resolve_key(data, key_path):
    """Navigate nested dict/list by dot-separated key path."""
    value = data
    for part in key_path.split("."):
        if isinstance(value, dict):
            value = value[part]
        elif isinstance(value, list):
            value = value[int(part)]
        else:
            raise KeyError(f"Cannot navigate into {type(value)} with key '{part}'")
    return value


def load_value(check):
    """Load a value from an output file."""
    filepath = REPO_ROOT / check["file"]
    if not filepath.exists():
        return None, f"File not found: {filepath}"

    if filepath.suffix == ".json":
        with open(filepath) as f:
            data = json.load(f)

        # Handle computed values (e.g., obs/null = distance / null_median)
        if check.get("compute") == "obs_null":
            try:
                dist = float(_resolve_key(data, check["key_distance"]))
                null_median = float(_resolve_key(data, check["key_null_median"]))
                return dist / null_median, None
            except (KeyError, IndexError, TypeError) as e:
                return None, f"Compute obs_null failed: {e}"

        try:
            value = _resolve_key(data, check["key"])
        except (KeyError, IndexError, TypeError) as e:
            return None, f"Key '{check['key']}' not found: {e}"
        return float(value), None
    elif filepath.suffix == ".csv":
        import pandas as pd
        df = pd.read_csv(filepath)
        try:
            value = _resolve_key({"df": df}, check["key"])
        except Exception as e:
            return None, f"CSV extraction failed: {e}"
        return float(value), None
    else:
        return None, f"Unknown file type: {filepath.suffix}"


def run_checks():
    passed = 0
    failed = 0
    errors = 0

    print("=" * 72)
    print("  CellWarp Reproduction Validation")
    print("=" * 72)
    print()

    for check in CHECKS:
        name = check["name"]
        value, err = load_value(check)

        if err:
            print(f"  ERROR   {name}")
            print(f"          {err}")
            print()
            errors += 1
            continue

        if "expected_below" in check:
            ok = value < check["expected_below"]
            detail = f"got {value:.5e}, need < {check['expected_below']:.5e}"
        elif "expected_in_range" in check:
            low, high = check["expected_in_range"]
            ok = low <= value <= high
            detail = f"got {value:.2e}, need in [{low:.2e}, {high:.2e}]"
        else:
            diff = abs(value - check["expected"])
            ok = diff <= check["tolerance"]
            detail = (f"got {value:.4f}, expected {check['expected']:.4f}, "
                      f"diff {diff:.4f}, tol {check['tolerance']}")

        status = "PASS" if ok else "FAIL"
        if ok:
            passed += 1
        else:
            failed += 1

        symbol = "+" if ok else "X"
        print(f"  {symbol} {status:5s}  {name}")
        print(f"          {detail}")
        print(f"          ref: {check['paper_ref']}")
        print()

    print("=" * 72)
    total = passed + failed + errors
    print(f"  {passed}/{total} passed | {failed} failed | {errors} errors")
    print("=" * 72)

    if failed > 0 or errors > 0:
        print("\nValidation FAILED. See above for details.")
        sys.exit(1)
    else:
        print("\nAll checks passed. Reproduction successful.")
        sys.exit(0)


if __name__ == "__main__":
    run_checks()
