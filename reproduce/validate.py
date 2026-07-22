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
        "name": "Global coherence p-value",
        "file": "output/phase2/scaled_35types/procrustes_results_35.json",
        "key": "permutation_test.p_value",
        "expected_below": 0.001,
        "paper_ref": "Results section 1",
    },

    # ── Independent PCA (Figure S1A-B) ──
    {
        "name": "Independent PCA obs/null",
        "file": "analysis/independent_pca_sensitivity/independent_pca_results.json",
        "key": "permutation_test.obs_null_ratio",
        "expected": 0.473,
        "tolerance": 0.005,
        "paper_ref": "Figure S1A",
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
        "paper_ref": "Figure S2C-D",
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
        "paper_ref": "Figure S2A-B",
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
        "paper_ref": "Results section 5, Figure 5C",
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
            detail = f"got {value:.2e}, need < {check['expected_below']:.2e}"
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
