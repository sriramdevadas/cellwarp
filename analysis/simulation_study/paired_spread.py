#!/usr/bin/env python3
"""
Per-replicate rank recovery at n = 200 for the four planted-spread points.

Each spread point is evaluated at its own recalibrated signal, taken from
sweep_spread_results.json in this directory, so the spread moves without dragging
the signal with it. Execs the deposited simulation_study.py unmodified, as
sweep_spread.py does.

Backs the spread sensitivity of the rank-recovery ceiling described alongside the
simulation in Methods and S1 Text.

All inputs are tracked.

Only the three path constants differ from the version that produced those values;
the code is unchanged.
"""
import json, time
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
SRC = HERE / "simulation_study.py"
OUT = HERE
ns={"__name__":"_sim_defs","__file__":str(SRC)}
exec(compile(SRC.read_text(),str(SRC),"exec"),ns)
sw=json.load(open(HERE / "sweep_spread_results.json"))["sweep"]
pts=[(x["label"],x["sigma"],x["calibrated_signal"]) for x in sw if x["sigma"]>0]
NC=200; res={}
for lab,sig,cal in pts:
    ns["RIGIDITY_SPREAD"]=sig
    t0=time.time(); rhos=[]
    for rep in range(ns["N_REPS"]):
        rhos.append(ns["run_single"](35,NC,cal,rep+30_000+NC*100)["ranking_recovery_rho"])
    a=np.array(rhos)
    res[lab]=dict(sigma=sig,cal_signal=cal,rhos=[float(x) for x in a],
                  mean=float(a.mean()),std=float(a.std()),median=float(np.median(a)))
    print(f"  {lab:<34} median={np.median(a):.8f} mean={a.mean():.8f}  [{time.time()-t0:.1f}s]")
(OUT/"paired_spread.json").write_text(json.dumps(res,indent=1))
