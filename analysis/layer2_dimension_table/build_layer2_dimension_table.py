#!/usr/bin/env python3
"""
S14 Table: the primate Layer-2 statistic against retained dimension, per pair, weighting
and subspace dimension k.

Eighteen rows: three primate pairs x two weightings x k in {1, 3, 5}. For each, the
retained joint-PCA dimension p, the chance level k/p that follows from it, and both
rotation arms' statistic, permutation-null mean, observed-minus-null margin and p-value.

WHY THE MARGIN COLUMNS ARE HERE

The paper's Fig 2C claim is that post-rotation sits below pre-rotation in every pair at
every k, each above its own null. That is true, and the margin columns are what let a
reader see how much room each cell actually has. They are not uniform: k = 1 carries the
smallest observed-minus-null margin in all twelve pair-weighting-arm combinations, and in
the primary human-mouse pair the k = 1 statistic does not clear its null at all. A reader
who can see margin against k sees that without being told it.

WHAT THIS READS, AND WHAT IT DOES NOT TOUCH

Input is the three vendored basal-ganglia results files under
docs/submission/plosone/figures/bg_results/. Those five files are read by Gate 1 and by
Fig 2C's producer; nothing here regenerates them. The retained dimension is already in
them, as the `k` field of each weighting block, so this needs no basal-ganglia npz and no
recomputation of any statistic -- every S and null value in the table is copied through
unchanged, and the reproduction is asserted against the source before the file is written.

Outputs (tracked):
  docs/supplementary_materials/table_S14_layer2_dimension.csv
  analysis/layer2_dimension_table/table_S14_summary.json

The JSON exists only so reproduce/validate.py can gate the table's OWN content. Its
.csv branch is unreachable -- it hands {"df": DataFrame} to a resolver that can only
navigate dicts and lists, so any key raises (see figure_script_map.md, Known gaps). The
summary is therefore built by reading the CSV back off disk after writing it, not from
the in-memory rows, so what the gate checks is the file a reader receives.

Deterministic, no network, runtime under a second. Not invoked by reproduce/run_all.sh.
Gated by reproduce/validate.py.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
BG = REPO / "docs/submission/plosone/figures/bg_results"
OUT = REPO / "docs/supplementary_materials/table_S14_layer2_dimension.csv"

PAIRS = ("Human_Macaque", "Human_Marmoset", "Macaque_Marmoset")
SCHEMES = ("W0_unscaled", "W2_schemeB")
KS = (1, 3, 5)

COLUMNS = [
    "pair", "weighting", "k",
    "retained_dimension_p", "chance_level_k_over_p",
    "S_pre", "null_mean_pre", "margin_pre", "p_pre",
    "S_post", "null_mean_post", "margin_post", "p_post",
]


def main() -> None:
    rows = []
    for pair in PAIRS:
        src = json.loads((BG / f"layer2_results_{pair}.json").read_text())
        for scheme in SCHEMES:
            block = src[scheme]
            p_dim = int(block["k"])
            if block["scheme"] != scheme:
                raise SystemExit(f"{pair} {scheme}: block self-label is {block['scheme']!r}")
            for k in KS:
                d = block["layer2"][f"k{k}"]
                rows.append({
                    "pair": pair.replace("_", "-"),
                    "weighting": scheme,
                    "k": k,
                    "retained_dimension_p": p_dim,
                    "chance_level_k_over_p": round(k / p_dim, 4),
                    "S_pre": d["S_pre"],
                    "null_mean_pre": d["null_mean_pre"],
                    "margin_pre": d["S_pre"] - d["null_mean_pre"],
                    "p_pre": d["p_pre"],
                    "S_post": d["S_post"],
                    "null_mean_post": d["null_mean_post"],
                    "margin_post": d["S_post"] - d["null_mean_post"],
                    "p_post": d["p_post"],
                })

    if len(rows) != 18:
        raise SystemExit(f"expected 18 rows, built {len(rows)}")

    # Reproduce-before-writing: every S and null in the table must be the deposited
    # value unchanged. A table that quietly differs from its source is worse than no
    # table, and this is the only place the two can be compared.
    for r in rows:
        src = json.loads((BG / f"layer2_results_{r['pair'].replace('-', '_')}.json").read_text())
        d = src[r["weighting"]]["layer2"][f"k{r['k']}"]
        for col, key in (("S_pre", "S_pre"), ("null_mean_pre", "null_mean_pre"),
                         ("p_pre", "p_pre"), ("S_post", "S_post"),
                         ("null_mean_post", "null_mean_post"), ("p_post", "p_post")):
            if r[col] != d[key]:
                raise SystemExit(
                    f"{r['pair']} {r['weighting']} k={r['k']} {col}: table {r[col]!r} "
                    f"differs from the deposited {d[key]!r}; refusing to write")
        if r["retained_dimension_p"] != int(src[r["weighting"]]["k"]):
            raise SystemExit(f"{r['pair']} {r['weighting']}: dimension mismatch")

    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(rows)

    # Summary for the gate, built by reading the CSV back rather than from `rows`, so
    # what validate.py checks is the file on disk and not this script's intent.
    with OUT.open(newline="") as fh:
        back = list(csv.DictReader(fh))
    if [c for c in back[0]] != COLUMNS:
        raise SystemExit(f"CSV header round-trip mismatch: {list(back[0])}")
    summary = {
        "table": str(OUT.relative_to(REPO)),
        "note": "read back from the written CSV; validate.py cannot read a .csv directly",
        "n_rows": len(back),
        "columns": COLUMNS,
        "retained_dimension": {f"{r['pair']}|{r['weighting']}": int(r["retained_dimension_p"])
                               for r in back},
        "chance_level": {f"{r['pair']}|{r['weighting']}|k{r['k']}":
                         float(r["chance_level_k_over_p"]) for r in back},
        "min_margin_pre": min(float(r["margin_pre"]) for r in back),
        "min_margin_post": min(float(r["margin_post"]) for r in back),
        "n_rows_k1": sum(1 for r in back if r["k"] == "1"),
    }
    (HERE / "table_S14_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    print(f"reproduction: all {len(rows)} rows match the deposited values")
    print(f"{'pair':<18}{'weighting':<13}{'k':>2}{'p':>4}{'k/p':>8}"
          f"{'margin_pre':>12}{'margin_post':>13}")
    for r in rows:
        print(f"{r['pair']:<18}{r['weighting']:<13}{r['k']:>2}{r['retained_dimension_p']:>4}"
              f"{r['chance_level_k_over_p']:>8.4f}{r['margin_pre']:>12.6f}{r['margin_post']:>13.6f}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
