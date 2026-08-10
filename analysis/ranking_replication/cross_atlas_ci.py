#!/usr/bin/env python3
"""
95% confidence intervals on the four cross-atlas rank correlations.

Question: Results section 4 reports four cross-atlas Spearman correlations
(+0.15 Sun2023, +0.19 PanSci, -0.05 pan-Census, -0.14 CellHint) and reads them
together as a non-replication. At n = 15 to 22 that reading needs an interval
and not only a point estimate, because the claim being made is about what the
data rule out: each interval spans both zero and 0.20, so these samples are
equally consistent with no cross-atlas ranking signal and with a moderate one.
This script computes those intervals from the deposited correlations, so the
number that reaches the submitted text is produced by tracked code rather than
by hand.

Method: Fisher z-transform. z = artanh(rho); interval z +/- 1.959964 * SE;
back-transformed with tanh.

    SE = 1.06 / sqrt(n - 3)        <- Bonett-Wright, for a SPEARMAN rho

DO NOT "CORRECT" THE 1.06 TO 1.0. Three scripts in this tree use the Pearson
standard error, SE = 1 / sqrt(n - 3):

    scripts/generate_phase2_figures.py:193
    scripts/t3e_step2_compute.py:233
    scripts/t3e_step3b_enhancer.py:553      (spearman_ci)

Those are correct where they stand and wrong here. 1/sqrt(n-3) is the
asymptotic standard error of the Fisher transform of a PEARSON r; the rank
transform behind a Spearman rho costs precision that it does not price in.
Bonett and Wright (2000), Psychometrika 65:23-28, give 1.06 as the Spearman
constant. Using 1.0 here would report intervals about 6% narrower on the z
scale than the data support -- the direction that makes a non-replication look
better resolved than it is, which is exactly the error this analysis exists to
avoid. The constant differs from its neighbours deliberately.

Inputs (all tracked). Four separate artifacts, three JSON and one CSV; no
single upstream producer holds all four correlations:

    output/validation/sun2023_replication_expanded/sun2023_expanded.json
    output/validation/pansci_replication/pansci_replication.json
    analysis/census_replication/replication_results.json
    analysis/harmonized_replication/sensitivity_analysis.csv

Outputs (tracked):

    analysis/ranking_replication/cross_atlas_ci_results.json

Deterministic: arithmetic over four deposited numbers. No resampling, no seed,
no network, runtime under a second. Not invoked by reproduce/run_all.sh, on the
same footing as analysis/simulation_study/sweep_spread.py. Gated by
reproduce/validate.py.
"""
import csv
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
OUT = HERE / "cross_atlas_ci_results.json"

# Phi^-1(0.975): the two-sided 95% normal deviate, spelled out so the interval
# does not silently depend on which scipy is installed.
Z_CRIT = 1.959963984540054

# Bonett-Wright Spearman constant. Read the header before touching this.
BW_SE_NUMERATOR = 1.06

# The "moderate correlation" the manuscript's sentence names. An interval that
# contains this as well as zero is what makes the four uninformative rather
# than negative.
TARGET = 0.20

CELLHINT_CSV = "analysis/harmonized_replication/sensitivity_analysis.csv"
CELLHINT_LEVEL = "0_unharmonized"
CELLHINT_EXPECTED = -0.139        # the matched 15-type baseline the paper reports
CELLHINT_WRONG_ROW = -0.386       # the pre-matching artifact; see read_cellhint()


def fisher_z_ci(rho, n):
    """Bonett-Wright Fisher-z 95% interval for a Spearman rho. -> (se, lo, hi)."""
    if n <= 3:
        raise SystemExit(f"n = {n} leaves no residual degrees of freedom for SE")
    se = BW_SE_NUMERATOR / math.sqrt(n - 3)
    z = math.atanh(rho)
    return se, math.tanh(z - Z_CRIT * se), math.tanh(z + Z_CRIT * se)


def _dig(data, dotted):
    value = data
    for part in dotted.split("."):
        value = value[part]
    return value


def read_json_pair(relpath, rho_key, n_key):
    data = json.loads((REPO / relpath).read_text())
    return float(_dig(data, rho_key)), int(_dig(data, n_key)), relpath, rho_key, n_key


def read_cellhint():
    """
    The trap in this file. Two different numbers are both called the
    "unharmonized" CellHint correlation, and the wrong one is a plausible read
    that nothing else here would catch.

      -0.386  The pre-PCA-matching artifact. S1 Text line 48 calls it exactly
              that: "The apparent unharmonized CellHint rho = -0.39 is an
              artifact of comparing residuals across incompatible PCA
              subspaces, resolved by same-dimension recomputation (S3 Table)."
              It is also what
              output/validation/cellhint_replication/cellhint_replication.json
              carries as rigidity_ranking.rho -- the obvious place to look and
              the wrong number -- and it is S4 Table's "Original" row.

      -0.139  The matched 15-type baseline. This is the CellHint arm the
              manuscript and S1 Text actually report, and S4 Table's "Level 0"
              row. sensitivity_analysis.csv labels it 0_unharmonized, where
              "unharmonized" means no ontology or tissue restriction applied --
              NOT the pre-matching artifact S1 Text means by the same word.

    So the CSV and S1 Text use "unharmonized" for different objects. The row is
    selected by its own label and both values are asserted below: the one we
    want, and the one we must not have picked up.
    """
    path = REPO / CELLHINT_CSV
    with path.open() as handle:
        rows = [r for r in csv.DictReader(handle) if r["level"] == CELLHINT_LEVEL]
    if len(rows) != 1:
        raise SystemExit(
            f"{CELLHINT_CSV}: expected exactly one {CELLHINT_LEVEL!r} row, found {len(rows)}"
        )
    rho = float(rows[0]["rho"])
    n = int(rows[0]["n_types"])
    if abs(rho - CELLHINT_EXPECTED) > 5e-4:
        raise SystemExit(
            f"{CELLHINT_CSV} row {CELLHINT_LEVEL!r}: rho = {rho:.6f}, expected the "
            f"matched 15-type baseline {CELLHINT_EXPECTED}. Refusing to build an "
            f"interval on an unrecognised CellHint value."
        )
    if abs(rho - CELLHINT_WRONG_ROW) < 5e-3:
        raise SystemExit(
            f"{CELLHINT_CSV} row {CELLHINT_LEVEL!r}: rho = {rho:.6f} is the "
            f"pre-PCA-matching artifact {CELLHINT_WRONG_ROW}, not the matched "
            f"15-type baseline the manuscript reports. See read_cellhint()."
        )
    return rho, n, CELLHINT_CSV, f"level=={CELLHINT_LEVEL} -> rho", f"level=={CELLHINT_LEVEL} -> n_types"


SOURCES = [
    ("sun2023", "Sun2023",
     lambda: read_json_pair(
         "output/validation/sun2023_replication_expanded/sun2023_expanded.json",
         "rigidity_ranking.rho", "rigidity_ranking.n_matched_types")),
    ("pansci", "PanSci",
     lambda: read_json_pair(
         "output/validation/pansci_replication/pansci_replication.json",
         "rigidity_ranking.rho", "rigidity_ranking.n_matched_types")),
    ("pan_census", "pan-Census",
     lambda: read_json_pair(
         "analysis/census_replication/replication_results.json",
         "ranking_correlation.spearman_rho", "ranking_correlation.n_types")),
    ("cellhint", "CellHint (matched 15-type baseline)", read_cellhint),
]


def main():
    print("=" * 78)
    print("CROSS-ATLAS RANK CORRELATIONS: 95% CONFIDENCE INTERVALS")
    print("=" * 78)
    print(f"  SE = {BW_SE_NUMERATOR} / sqrt(n - 3)   (Bonett-Wright, Spearman)")
    print(f"  z_crit = {Z_CRIT}")
    print()

    pairs = {}
    for key, label, reader in SOURCES:
        rho, n, source, rho_key, n_key = reader()
        se, lo, hi = fisher_z_ci(rho, n)
        contains_zero = lo < 0.0 < hi
        contains_target = lo < TARGET < hi
        pairs[key] = {
            "label": label,
            "source": source,
            "rho_key": rho_key,
            "n_key": n_key,
            "rho": rho,
            "n": n,
            "se": se,
            "ci_lower": lo,
            "ci_upper": hi,
            "ci_contains_zero": contains_zero,
            "ci_contains_0_20": contains_target,
        }
        print(f"  {label:38s} rho = {rho:+.4f}  n = {n:2d}  SE = {se:.6f}")
        print(f"  {'':38s} 95% CI [{lo:+.6f}, {hi:+.6f}]   "
              f"contains 0: {contains_zero}   contains {TARGET}: {contains_target}")
        print(f"  {'':38s} source: {source}")
        print()

    result = {
        "method": {
            "transform": "Fisher z (artanh), back-transformed with tanh",
            "se_formula": "1.06 / sqrt(n - 3)",
            "se_rationale": (
                "Bonett-Wright Spearman standard error. NOT the Pearson "
                "1/sqrt(n-3) used by scripts/generate_phase2_figures.py, "
                "scripts/t3e_step2_compute.py and scripts/t3e_step3b_enhancer.py; "
                "see this script's module docstring before changing it."
            ),
            "reference": "Bonett DG, Wright TA (2000). Psychometrika 65:23-28",
            "z_crit": Z_CRIT,
            "confidence": 0.95,
            "target": TARGET,
        },
        "pairs": pairs,
        "summary": {
            "n_pairs": len(pairs),
            "all_contain_zero": all(p["ci_contains_zero"] for p in pairs.values()),
            "all_contain_0_20": all(p["ci_contains_0_20"] for p in pairs.values()),
            "n_range": [min(p["n"] for p in pairs.values()),
                        max(p["n"] for p in pairs.values())],
        },
    }

    OUT.write_text(json.dumps(result, indent=2) + "\n")
    print("=" * 78)
    print(f"  all four intervals contain 0:    {result['summary']['all_contain_zero']}")
    print(f"  all four intervals contain {TARGET}: {result['summary']['all_contain_0_20']}")
    print(f"  n range: {result['summary']['n_range']}")
    print("=" * 78)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
