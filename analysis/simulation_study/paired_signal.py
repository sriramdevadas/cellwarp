#!/usr/bin/env python3
"""
Per-replicate rank recovery at n = 200 under three planted signals.

Execs the deposited simulation_study.py rather than editing it; its main() is
__name__-guarded, so the definitions load without running anything. Reports recovery
per replicate rather than pooled, so the three signal points are paired across the
same drawn configurations.

Backs the rank-recovery ceiling described in the Methods simulation section and the
recovery values S1 Text reports alongside it.

All inputs are tracked.

Only the two path constants differ from the version that produced those values; the
code is unchanged.
"""
import json, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
SRC = HERE / "simulation_study.py"
OUT = HERE
ns = {"__name__": "_sim_defs", "__file__": str(SRC)}
exec(compile(SRC.read_text(), str(SRC), "exec"), ns)
print("RIGIDITY_SPREAD =", ns["RIGIDITY_SPREAD"], " N_REPS =", ns["N_REPS"])

NC = 200
SIGNALS = {"3.0": 3.0, "3.683416": 3.683416443528231, "5.0": 5.0}
res = {}
for lab, sig in SIGNALS.items():
    t0 = time.time(); rhos = []
    for rep in range(ns["N_REPS"]):
        seed = rep + 30_000 + NC * 100          # identical to simulation_study.py:424
        rhos.append(ns["run_single"](35, NC, sig, seed)["ranking_recovery_rho"])
    a = np.array(rhos)
    res[lab] = dict(signal=sig, rhos=[float(x) for x in a],
                    mean=float(a.mean()), std=float(a.std()), median=float(np.median(a)))
    print(f"  signal {lab:10s} mean={a.mean():.8f} median={np.median(a):.8f} "
          f"sd={a.std():.6f}  [{time.time()-t0:.1f}s]")
(OUT/"paired_signal.json").write_text(json.dumps(res, indent=1))
print("wrote", OUT/"paired_signal.json")
