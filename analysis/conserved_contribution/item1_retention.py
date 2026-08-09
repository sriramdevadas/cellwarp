#!/usr/bin/env python3
"""
Paired retention ratio, formed per case from that case's own arms.

    advantage_out_of_sample = F - T
    advantage_in_sample     = F - I
    retention               = (F - T) / (F - I)

Lower obs/null is stronger and F is the weakest arm, so advantages are measured
downward from F, and each ratio is formed per split per direction from that case's
own F rather than from a pooled value.

Produces the three values S1 Text reports in one sentence: "the out-of-sample set
retains a median 90% of the advantage (interquartile range 79% to 99%)".

NOT TRACKED: reads block3_form_b_results.json, written by block3_form_b.py in this
directory. That file is not deposited; run block3_form_b.py first.

Only the output path constant differs from the version that produced those values;
the code is unchanged.
"""
import json
from pathlib import Path

import numpy as np

OUT = HERE
d = json.load(open(OUT / "block3_form_b_results.json"))


def arm(rec, prefix):
    k = [x for x in rec if x.startswith(prefix)][0]
    return rec[k]["ratio"]


cases = []
for rec in d["splits"]:
    for dirn, tag in (("select_A_eval_B", "A->B"), ("select_B_eval_A", "B->A")):
        T, I, F = arm(rec[dirn], "T "), arm(rec[dirn], "I "), arm(rec[dirn], "F ")
        cases.append(dict(split=rec["split"], direction=tag, T=T, I=I, F=F,
                          adv_out=F - T, adv_in=F - I,
                          retention=(F - T) / (F - I)))

n = len(cases)
ret = np.array([c["retention"] for c in cases])
aout = np.array([c["adv_out"] for c in cases])
ain = np.array([c["adv_in"] for c in cases])

print("=" * 104)
print("(1) PER-CASE TABLE -- retention formed from each case's own three arms")
print("=" * 104)
print(f"{'split':>5s} {'dir':>5s} {'T':>8s} {'I':>8s} {'F':>8s} {'F-T':>9s} {'F-I':>9s} {'retention':>10s}")
for c in cases:
    print(f"{c['split']:5d} {c['direction']:>5s} {c['T']:8.4f} {c['I']:8.4f} {c['F']:8.4f} "
          f"{c['adv_out']:+9.4f} {c['adv_in']:+9.4f} {c['retention']:10.4f}")


def q(x, p):
    return float(np.percentile(x, p))


print("\n" + "=" * 104)
print(f"(2) RETENTION across {n} cases")
print("=" * 104)
print(f"  median        : {np.median(ret):.4f}")
print(f"  IQR           : [{q(ret,25):.4f}, {q(ret,75):.4f}]   (width {q(ret,75)-q(ret,25):.4f})")
print(f"  full range    : [{ret.min():.4f}, {ret.max():.4f}]")
print(f"  mean +/- sd   : {ret.mean():.4f} +/- {ret.std(ddof=1):.4f}")

print("\n" + "=" * 104)
print("(3) CASES WITH retention > 1  (out-of-sample beat in-sample -- noise, not signal)")
print("=" * 104)
gt1 = [c for c in cases if c["retention"] > 1]
print(f"  count: {len(gt1)} / {n}   ({100*len(gt1)/n:.1f}%)")
print(f"  (Dispatch 32 reported I < T in 30/40, so 10/40 exceeding 1 is the same fact)")
for c in gt1:
    print(f"    split {c['split']:2d} {c['direction']}  retention {c['retention']:.4f}  "
          f"(F-T {c['adv_out']:+.4f} > F-I {c['adv_in']:+.4f})")
print("  These are cases where the out-of-sample selection happened to land on a better")
print("  gene set than the in-sample one. A retention above 1 is not evidence of anything;")
print("  it is sampling noise in a quantity whose expectation is at most 1.")

print("\n" + "=" * 104)
print("(4) DENOMINATOR DIAGNOSTICS -- is the quotient well-conditioned?")
print("=" * 104)
print(f"  F-I  min    : {ain.min():+.4f}")
print(f"  F-I  median : {np.median(ain):+.4f}")
print(f"  F-I  max    : {ain.max():+.4f}")
print(f"  F-I  any <= 0 ? {int((ain <= 0).sum())} cases")
print(f"  smallest |F-I| as a fraction of the median: {abs(ain).min()/abs(np.median(ain)):.3f}")
order = np.argsort(ain)
print(f"\n  five smallest denominators and their retention:")
for i in order[:5]:
    c = cases[i]
    print(f"    split {c['split']:2d} {c['direction']}  F-I {c['adv_in']:+.4f}  "
          f"F-T {c['adv_out']:+.4f}  retention {c['retention']:.4f}")
UNSTABLE = 0.02
uns = [c for c in cases if abs(c["adv_in"]) < UNSTABLE]
print(f"\n  cases with |F-I| < {UNSTABLE} (would make the ratio unstable): {len(uns)}")
if uns:
    for c in uns:
        print(f"    REPORTED SEPARATELY, NOT TRIMMED: split {c['split']} {c['direction']} "
              f"F-I {c['adv_in']:+.4f} retention {c['retention']:.4f}")
else:
    print("    none -- every denominator is far from zero, so the quotient is well-conditioned")

print("\n" + "=" * 104)
print("(5) THE TWO COMPONENTS SEPARATELY, in obs/null units")
print("=" * 104)
for nm, v in (("F-T  (out-of-sample advantage)", aout), ("F-I  (in-sample advantage)", ain)):
    print(f"  {nm:34s} median {np.median(v):+.4f}   IQR [{q(v,25):+.4f}, {q(v,75):+.4f}]   "
          f"range [{v.min():+.4f}, {v.max():+.4f}]")

print("\n" + "=" * 104)
print("NEGATIVE CONTROL (note 7): substitute F for I, i.e. (F-T)/(F-F)")
print("=" * 104)
c0 = cases[0]
print(f"  case: split {c0['split']} {c0['direction']}   F = {c0['F']:.6f}")
num = c0["F"] - c0["T"]
den = c0["F"] - c0["F"]
print(f"  numerator (F-T) = {num:.6f}   denominator (F-F) = {den:.6f}")
try:
    with np.errstate(divide="raise", invalid="raise"):
        bad = np.float64(num) / np.float64(den)
    print(f"  RESULT: {bad}")
    print(f"  -> NEGCTL {'ok (non-finite)' if not np.isfinite(bad) else 'WRONG: produced a finite, printable number'}")
except FloatingPointError as e:
    print(f"  -> NEGCTL ok: raised FloatingPointError({e}) rather than returning a number")
# and the plain-Python path, which is what a careless script would use
try:
    bad2 = num / den
    print(f"  plain python float division: {bad2}  -> WOULD HAVE PRINTED A NUMBER")
except ZeroDivisionError as e:
    print(f"  plain python float division: raised ZeroDivisionError({e}) -- safe")

R = dict(
    n_cases=n,
    cases=cases,
    retention=dict(median=float(np.median(ret)), iqr=[q(ret, 25), q(ret, 75)],
                   min=float(ret.min()), max=float(ret.max()),
                   mean=float(ret.mean()), sd=float(ret.std(ddof=1)),
                   n_gt_1=len(gt1)),
    adv_out=dict(median=float(np.median(aout)), iqr=[q(aout, 25), q(aout, 75)],
                 min=float(aout.min()), max=float(aout.max())),
    adv_in=dict(median=float(np.median(ain)), iqr=[q(ain, 25), q(ain, 75)],
                min=float(ain.min()), max=float(ain.max()),
                n_le_zero=int((ain <= 0).sum()), n_unstable=len(uns),
                unstable_threshold=UNSTABLE),
)
(OUT / "item1_retention_results.json").write_text(json.dumps(R, indent=2))
print(f"\nwrote {OUT/'item1_retention_results.json'}")
