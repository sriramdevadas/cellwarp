#!/usr/bin/env python3
"""
CellWarp — TOST Equivalence Test for SAMap vs Rigidity Correlation

Tests whether the observed Spearman ρ between SAMap correspondence scores
and CellWarp rigidity scores is bounded within |ρ| < 0.50, using the
Two One-Sided Tests (TOST) procedure on Fisher z-transformed correlations.

Biology
-------
SAMap (Tarashansky et al. 2021) and CellWarp measure different aspects of
cross-species cell type relationships. A non-significant correlation
(p=0.153) does not prove orthogonality — it only fails to prove dependence.
TOST provides a positive statistical claim: the correlation is bounded
within a pre-specified equivalence region (|ρ| < 0.50), meaning the two
methods measure sufficiently different properties to be considered
independent for the purposes of the paper's argument.

Math
----
1. Convert observed ρ to Fisher z:
       z_obs = 0.5 * ln((1 + ρ) / (1 - ρ))

2. Convert equivalence bound Δ = 0.50 to Fisher z:
       z_bound = 0.5 * ln((1 + Δ) / (1 - Δ))

3. Standard error of Fisher z for Spearman (same as Pearson asymptotic):
       SE = 1 / sqrt(n - 3)

4. Two one-sided t-tests against ±z_bound:
       t1 = (z_obs - (-z_bound)) / SE   (tests H1: ρ > -Δ)
       t2 = (z_bound - z_obs) / SE      (tests H2: ρ < +Δ)

5. TOST p-value = max(p1, p2), where p1 and p2 are one-tailed p-values
   from t-distribution with df = n - 3.

6. Equivalence declared if TOST p < 0.05.

Inputs:
    output/phase1_samap/samap_35types/samap_rigidity_correlation.json

Outputs:
    output/phase1_samap/samap_35types/tost_equivalence_result.json

Usage:
    python scripts/v7_tost_equivalence.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

RANDOM_SEED = 42  # Unused here but included for project consistency

PROJECT_ROOT = Path(__file__).resolve().parent.parent
INPUT_FILE = PROJECT_ROOT / "output/phase1_samap/samap_35types/samap_rigidity_correlation.json"
OUTPUT_FILE = PROJECT_ROOT / "output/phase1_samap/samap_35types/tost_equivalence_result.json"

EQUIVALENCE_BOUND = 0.50  # |ρ| < 0.50
ALPHA = 0.05


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("=" * 70)
    print("V7 — TOST EQUIVALENCE TEST: SAMap vs CellWarp Rigidity")
    print("=" * 70)

    # ---- Load observed result ----
    print(f"\n  Loading: {INPUT_FILE}")
    with open(INPUT_FILE) as f:
        data = json.load(f)

    rho = data["spearman_rho"]
    n = data["n"]

    print(f"  Observed Spearman ρ: {rho:.6f}")
    print(f"  n (cell types):      {n}")
    print(f"  Equivalence bound:   |ρ| < {EQUIVALENCE_BOUND}")

    # ---- Fisher z-transform ----
    z_obs = 0.5 * np.log((1 + rho) / (1 - rho))
    z_bound = 0.5 * np.log((1 + EQUIVALENCE_BOUND) / (1 - EQUIVALENCE_BOUND))
    se = 1.0 / np.sqrt(n - 3)
    df = n - 3

    print(f"\n  Fisher z-transform:")
    print(f"    z_obs   = {z_obs:.6f}")
    print(f"    z_bound = {z_bound:.6f}")
    print(f"    SE      = {se:.6f}")
    print(f"    df      = {df}")

    # ---- Two one-sided tests ----
    # Test 1: H0: ρ ≤ -Δ  vs  H1: ρ > -Δ
    t1 = (z_obs - (-z_bound)) / se
    p1 = stats.t.sf(t1, df)  # one-tailed: P(T > t1)

    # Test 2: H0: ρ ≥ +Δ  vs  H1: ρ < +Δ
    t2 = (z_bound - z_obs) / se
    p2 = stats.t.sf(t2, df)  # one-tailed: P(T > t2)

    tost_p = max(p1, p2)
    equivalence = tost_p < ALPHA

    print(f"\n  One-sided test 1 (H1: ρ > -{EQUIVALENCE_BOUND}):")
    print(f"    t1 = {t1:.6f}")
    print(f"    p1 = {p1:.6f}")

    print(f"\n  One-sided test 2 (H1: ρ < +{EQUIVALENCE_BOUND}):")
    print(f"    t2 = {t2:.6f}")
    print(f"    p2 = {p2:.6f}")

    print(f"\n  TOST p-value: {tost_p:.6f}")
    print(f"  α threshold:  {ALPHA}")

    # ---- Conclusion ----
    print("\n" + "=" * 70)
    if equivalence:
        print(f"  CONCLUSION: Equivalence DECLARED (TOST p={tost_p:.4f} < {ALPHA})")
        print(f"  The SAMap-CellWarp correlation is bounded within |ρ| < {EQUIVALENCE_BOUND}.")
        print(f"  The two methods measure different biological properties.")
    else:
        print(f"  CONCLUSION: Equivalence NOT declared (TOST p={tost_p:.4f} ≥ {ALPHA})")
        print(f"  Cannot confirm the correlation is bounded within |ρ| < {EQUIVALENCE_BOUND}.")
    print("=" * 70)

    # ---- Save results ----
    result = {
        "analysis": "TOST equivalence test — SAMap vs CellWarp rigidity",
        "input_file": str(INPUT_FILE.relative_to(PROJECT_ROOT)),
        "observed_spearman_rho": float(rho),
        "n": int(n),
        "equivalence_bound": float(EQUIVALENCE_BOUND),
        "fisher_z_obs": float(z_obs),
        "fisher_z_bound": float(z_bound),
        "se": float(se),
        "df": int(df),
        "t1": float(t1),
        "p1": float(p1),
        "t2": float(t2),
        "p2": float(p2),
        "tost_p": float(tost_p),
        "alpha": float(ALPHA),
        "equivalence_declared": bool(equivalence),
        "random_seed": RANDOM_SEED,
    }

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_FILE, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\n  Results saved: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
