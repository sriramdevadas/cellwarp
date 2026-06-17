#!/usr/bin/env python3
"""
Recovery-only slice at signal = 7.0.

Re-uses the same constants, RNG seeding scheme, and run_single() from
simulation_study.py so the resulting mean_rho values are directly comparable
to the existing ranking_recovery entries in simulation_results.json.

Runs 4 n_cells values × 100 reps = 400 pipeline calls.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from simulation_study import (  # noqa: E402
    RECOVERY_N_CELLS,
    N_REPS,
    run_single,
)

SIGNAL = 7.0


def main() -> None:
    t0 = time.time()
    out = {
        "signal_strength": SIGNAL,
        "n_types": 35,
        "metric": "mean_rho (mean across N_REPS replicates of Spearman ρ "
                  "between planted and recovered per-type rankings)",
        "n_reps": N_REPS,
        "results": [],
    }
    print(f"Recovery sweep at signal = {SIGNAL}  (n_types=35, n_reps={N_REPS})")
    for nc in RECOVERY_N_CELLS:
        rhos = []
        for rep in range(N_REPS):
            seed = rep + 30_000 + nc * 100  # matches simulation_study.run_ranking_recovery
            r = run_single(35, nc, SIGNAL, seed)
            rhos.append(r["ranking_recovery_rho"])
        arr = np.array(rhos)
        row = {
            "signal_strength": SIGNAL,
            "n_cells_per_type": nc,
            "mean_rho": float(np.mean(arr)),
            "std_rho": float(np.std(arr)),
            "median_rho": float(np.median(arr)),
        }
        out["results"].append(row)
        print(f"  sig={SIGNAL:.1f}  n_cells={nc:5d}  "
              f"ρ = {row['mean_rho']:.3f} ± {row['std_rho']:.3f}")
    out["runtime_sec"] = time.time() - t0

    out_path = HERE / "ranking_recovery_signal_7.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved → {out_path}")
    print(f"Runtime: {out['runtime_sec']:.1f}s")


if __name__ == "__main__":
    main()
