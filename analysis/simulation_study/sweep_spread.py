#!/usr/bin/env python3
"""
Planted-spread sweep for the rank-recovery ceiling, re-calibrating at each spread.

Question: the deposited simulation reports rank recovery over RECOVERY_SIGNALS =
[0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]. That grid does not contain the calibrated
signal. calibrate() interpolates between the 3.0 and 4.0 calibration points to
reproduce the observed obs/null of 0.522, giving 3.683416443528231, and the
deposited run records it (calibration.estimated_real_signal, rounded to 3.683)
and runs the stability experiment there -- but not the recovery experiment. So
recovery at the signal the model was calibrated to is not among the deposited
values. This script evaluates it, and sweeps the one free assumption the ceiling
depends on: RIGIDITY_SPREAD, the log-normal sigma of the planted per-type
divergence.

Design: simulation_study.py is exec'd, not edited. Its main() is __name__-guarded,
so the module's definitions load without running anything, and run_single() and
calibrate() are used unmodified. RIGIDITY_SPREAD is rebound in that namespace per
sweep point, and calibrate() is re-run there, so the signal tracks the spread and
obs/null stays at 0.522 rather than the sweep moving two quantities at once.

Seeding is seed = rep + 30_000 + n_cells * 100, identical to simulation_study.py's
own recovery loop. Every sweep point therefore draws the same 100 configurations as
the deposited grid, and the two are paired rather than independent.

Spread points: a zero-spread negative control; three sigmas bisected to median
max/min planted ratios of 10x, 25x and 50x; and the deposited sigma = 1.0, whose
measured median max/min is 65.0x. "25x" is not a literal anywhere in the deposited
script, so each sigma is reported with its measured range under three readings.

Outputs (tracked):
  analysis/simulation_study/sweep_spread_results.json

Wall time about 480 s: 5 spread points x 4 cell counts x 100 replicates, plus
calibrate() over 11 signals x 30 replicates at each of the 4 non-zero spreads.
Deterministic -- every seed is a closed-form function of replicate and cell count.
Gated by reproduce/validate.py.
"""
import json
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SRC = HERE / "simulation_study.py"
OUT = HERE

ns = {"__name__": "_sim_defs", "__file__": str(SRC)}
exec(compile(SRC.read_text(), str(SRC), "exec"), ns)      # main() is guarded, so nothing runs
print("deposited module exec'd (main() guarded, not run)")
print(f"  deposited RIGIDITY_SPREAD = {ns['RIGIDITY_SPREAD']}   "
      f"REAL_OBS_NULL = {ns['REAL_OBS_NULL']}   N_REPS = {ns['N_REPS']}")


def measure_range(sigma, n_types=35, n_draws=20000, seed=0):
    """What multiplicative range does this sigma actually produce?"""
    rng = np.random.RandomState(seed)
    mm, pp, iq = [], [], []
    for _ in range(n_draws):
        v = np.exp(sigma * rng.randn(n_types))
        v /= v.mean()
        mm.append(v.max() / v.min() if v.min() > 0 else np.inf)
        pp.append(np.percentile(v, 95) / np.percentile(v, 5))
        iq.append(np.percentile(v, 75) / np.percentile(v, 25))
    return dict(sigma=sigma,
                max_min_median=float(np.median(mm)),
                p95_p5_median=float(np.median(pp)),
                iqr_ratio_median=float(np.median(iq)))


def sigma_for_target(target, key="max_min_median", lo=0.01, hi=3.0):
    """Bisect for the sigma whose median max/min equals `target`."""
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if measure_range(mid, n_draws=4000)[key] < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


print("\n" + "=" * 78)
print("WHAT THE DEPOSITED SIGMA ACTUALLY PRODUCES")
print("=" * 78)
base = measure_range(1.0)
print(f"  sigma=1.0  ->  max/min {base['max_min_median']:.1f}x   "
      f"p95/p5 {base['p95_p5_median']:.1f}x   IQR-ratio {base['iqr_ratio_median']:.2f}x")
print("  the script's own comment says '~10x range'; the manuscript says '~25x spread'")

targets = {}
for t in (10.0, 25.0, 50.0):
    s = sigma_for_target(t)
    targets[t] = s
    m = measure_range(s)
    print(f"  target {t:5.0f}x max/min -> sigma {s:.4f}  "
          f"(check: max/min {m['max_min_median']:.1f}x, p95/p5 {m['p95_p5_median']:.1f}x)")

# sweep points: the three targets, the deposited sigma, and the 1x negative control
SWEEP = [("1x NEGATIVE CONTROL (no spread)", 0.0),
         ("10x", targets[10.0]),
         ("25x", targets[25.0]),
         ("deposited sigma=1.0 (~65x)", 1.0),
         ("50x", targets[50.0])]

R = {"deposited": dict(base), "targets": {str(k): v for k, v in targets.items()},
     "n_reps": ns["N_REPS"], "recovery_n_cells": ns["RECOVERY_N_CELLS"], "sweep": []}

for label, sigma in SWEEP:
    t0 = time.time()
    print("\n" + "=" * 78)
    print(f"SPREAD POINT: {label}   RIGIDITY_SPREAD = {sigma:.4f}")
    print("=" * 78)
    ns["RIGIDITY_SPREAD"] = sigma
    rng_meta = measure_range(sigma) if sigma > 0 else dict(
        sigma=0.0, max_min_median=1.0, p95_p5_median=1.0, iqr_ratio_median=1.0)

    if sigma == 0.0:
        # calibration is meaningless with no planted spread; use the deposited signal
        real_signal = 3.683
        cal = None
        print(f"  (no spread: calibration skipped, deposited signal {real_signal} used)")
    else:
        real_signal, cal = ns["calibrate"]()
        print(f"  RE-CALIBRATED signal for this spread: {real_signal:.3f}")

    rec = []
    for nc in ns["RECOVERY_N_CELLS"]:
        rhos = []
        for rep in range(ns["N_REPS"]):
            seed = rep + 30_000 + nc * 100
            rhos.append(ns["run_single"](35, nc, real_signal, seed)["ranking_recovery_rho"])
        rhos = np.array(rhos)
        rec.append(dict(n_cells=nc, mean_rho=float(rhos.mean()), std_rho=float(rhos.std()),
                        median_rho=float(np.median(rhos)),
                        se=float(rhos.std() / np.sqrt(len(rhos)))))
        print(f"    n_cells={nc:5d}  rho = {rhos.mean():+.4f} +/- {rhos.std():.4f}  "
              f"median {np.median(rhos):+.4f}  SE {rhos.std()/np.sqrt(len(rhos)):.4f}")

    meds = [r["median_rho"] for r in rec]
    monotone = all(meds[i] >= meds[i + 1] for i in range(len(meds) - 1))
    print(f"  ceiling (median rho at n_cells=200): {rec[1]['median_rho']:+.4f}")
    print(f"  ordering monotonically decreasing with n: {monotone}")
    R["sweep"].append(dict(label=label, sigma=sigma, range=rng_meta,
                           calibrated_signal=float(real_signal),
                           calibration=cal, recovery=rec,
                           monotone_decreasing=bool(monotone),
                           sec=round(time.time() - t0, 1)))
    (OUT / "sweep_spread_results.json").write_text(json.dumps(R, indent=2, default=float))

print("\n" + "=" * 78)
print("SUMMARY: ceiling vs planted spread (all signals RE-CALIBRATED to obs/null 0.522)")
print("=" * 78)
print(f"  {'spread point':32s} {'sigma':>7s} {'max/min':>9s} {'signal':>8s} "
      f"{'rho@50':>8s} {'rho@200':>8s} {'rho@500':>8s} {'rho@2000':>9s}")
for s in R["sweep"]:
    m = {r["n_cells"]: r["median_rho"] for r in s["recovery"]}
    print(f"  {s['label']:32s} {s['sigma']:7.4f} {s['range']['max_min_median']:9.1f} "
          f"{s['calibrated_signal']:8.3f} {m[50]:8.4f} {m[200]:8.4f} {m[500]:8.4f} {m[2000]:9.4f}")
(OUT / "sweep_spread_results.json").write_text(json.dumps(R, indent=2, default=float))
print(f"\nwrote {OUT/'sweep_spread_results.json'}")
