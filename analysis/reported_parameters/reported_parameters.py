#!/usr/bin/env python3
"""
Reported parameters and summary statistics that the submitted texts state but no
artifact carried as a gateable scalar.

Every value here already appears, or is about to appear, in the manuscript or S1 Text.
The problem this solves is not that they were unknown but that they lived in places
reproduce/validate.py cannot read: a CSV column (its .csv branch is unreachable -- see
figure_script_map.md, Known gaps), a Python literal in a producer, or a per-type list
that has to be reduced before it is a number.

NOTHING HERE IS A NEW RESULT. Every field is read or derived from a tracked file in
this repository, and the source of each is recorded alongside it:

  * per-type cell counts for the three primate pairs, reduced to min and max, from the
    vendored bg_results/layer2_results_*.json `per_type_n` blocks;
  * mouse-lemur parameters, read from analysis/mouse_lemur/01_run_pipeline.py by
    parsing its module constants rather than transcribing them, plus the ortholog row
    count measured from its tracked BioMart CSV;
  * the bootstrap rank-stability criterion and its observed distribution, from
    analysis/cross_reference/master_ranking_table.csv.

WHAT IS DELIBERATELY NOT HERE

The primate per-type threshold itself (min_cells=100) is a literal in the
basal-ganglia deposit's analysis/bg/layer2_covariance.py:60-62, outside this
repository, so this script does not restate it as though this repository measured it.
It asserts the consequence instead -- that no matched type in any pair falls below 100
cells in either arm -- and records the measured minima, which is what the text states.

Values already carried by a gateable JSON scalar are NOT duplicated here. The four
per-gene-standardization numbers are in analysis/sensitivity_analyses/genestd_results.json
(layer1.A/B.obs_null, cpc1.A/B.n_ribosomal_dominated) and the mouse-lemur type count,
PCA dimension, gene space and permutation count are in
analysis/mouse_lemur/procrustes_results.json; both are gated directly.

Outputs (tracked):
  analysis/reported_parameters/reported_parameters.json

Deterministic, no network, runtime under a second. Not invoked by reproduce/run_all.sh.
Gated by reproduce/validate.py.
"""
from __future__ import annotations

import ast
import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
OUT = HERE / "reported_parameters.json"

BG_RESULTS = REPO / "docs/submission/plosone/figures/bg_results"
ML_PRODUCER = REPO / "analysis/mouse_lemur/01_run_pipeline.py"
ML_ORTHOLOGS = REPO / "analysis/mouse_lemur/biomart_mouse_lemur_human_orthologs.csv"
RANKING_TABLE = REPO / "analysis/cross_reference/master_ranking_table.csv"

PRIMATE_PAIRS = ("Human_Macaque", "Human_Marmoset", "Macaque_Marmoset")
PRIMATE_THRESHOLD_DECLARED = 100      # analysis/bg/layer2_covariance.py:60-62, other repo
STABLE_CRITERION_MAX_WIDTH = 10       # S1 Text: "All 35 types classified stable (95% CI width <= 10)"


def module_constants(path: Path, names: set[str]) -> dict:
    """Read module-level literal assignments by parsing the source, not importing it."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    out = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 \
                and isinstance(node.targets[0], ast.Name) and node.targets[0].id in names:
            try:
                out[node.targets[0].id] = ast.literal_eval(node.value)
            except ValueError:
                out[node.targets[0].id] = ast.unparse(node.value)
    missing = names - set(out)
    if missing:
        raise SystemExit(f"{path.name}: expected constants absent: {sorted(missing)}")
    return out


def main() -> None:
    result = {
        "analysis": "reported parameters and summary statistics not otherwise gateable",
        "note": "no new results; every field is read or derived from a tracked file in "
                "this repository, and its source is recorded beside it",
        "primate_replication": {},
        "mouse_lemur": {},
        "bootstrap_rank_stability": {},
    }

    # ── primate per-type cell counts ────────────────────────────────────────
    pr = {"source": "docs/submission/plosone/figures/bg_results/layer2_results_{pair}.json, "
                    "per_type_n",
          "threshold_declared_by_producer": PRIMATE_THRESHOLD_DECLARED,
          "threshold_source": "analysis/bg/layer2_covariance.py:60-62 in the basal-ganglia "
                              "deposit (another repository); this repository asserts its "
                              "consequence rather than restating the literal",
          "pairs": {}}
    for pair in PRIMATE_PAIRS:
        d = json.loads((BG_RESULTS / f"layer2_results_{pair}.json").read_text())
        per = d["per_type_n"]
        mins = [min(r["n_ref"], r["n_tgt"]) for r in per]
        maxs = [max(r["n_ref"], r["n_tgt"]) for r in per]
        lo = min(mins)
        if lo < PRIMATE_THRESHOLD_DECLARED:
            raise SystemExit(
                f"{pair}: a matched type has {lo} cells, below the producer's declared "
                f"threshold of {PRIMATE_THRESHOLD_DECLARED}")
        pr["pairs"][pair] = {
            "n_types": d["n_types"],
            "min_cells_either_arm": int(lo),
            "max_cells_either_arm": int(max(maxs)),
            "min_cells_type": min(per, key=lambda r: min(r["n_ref"], r["n_tgt"]))["cell_type"],
            "all_types_at_or_above_threshold": True,
        }
    result["primate_replication"] = pr

    # ── mouse-lemur parameters ──────────────────────────────────────────────
    c = module_constants(ML_PRODUCER, {"MIN_CELLS_PER_TYPE", "MAX_CELLS_PER_TYPE",
                                       "N_PERMUTATIONS", "RANDOM_SEED"})
    orth = pd.read_csv(ML_ORTHOLOGS)
    src = ML_PRODUCER.read_text(encoding="utf-8")
    human_arm = "output/phase2/scaled_35types/centroids_human_35.csv"
    if human_arm not in src:
        raise SystemExit("mouse-lemur producer no longer references the primary human centroids")
    result["mouse_lemur"] = {
        "source": "analysis/mouse_lemur/01_run_pipeline.py module constants (parsed, not "
                  "transcribed) and its tracked BioMart CSV",
        "human_arm": human_arm,
        "human_arm_note": "the primary 35-type Tabula Sapiens centroids, restricted to the "
                          "matched types; the lemur arm is built from cells",
        "ortholog_source": "analysis/mouse_lemur/biomart_mouse_lemur_human_orthologs.csv",
        "n_ortholog_pairs": int(len(orth)),
        "min_cells_per_type": c["MIN_CELLS_PER_TYPE"],
        "max_cells_per_type": c["MAX_CELLS_PER_TYPE"],
        "n_permutations": c["N_PERMUTATIONS"],
        "random_seed": c["RANDOM_SEED"],
    }

    # ── bootstrap rank stability ────────────────────────────────────────────
    t = pd.read_csv(RANKING_TABLE)
    w = t["bootstrap_CI_width"]
    result["bootstrap_rank_stability"] = {
        "source": "analysis/cross_reference/master_ranking_table.csv, bootstrap_CI_width "
                  "and bootstrap_category",
        "criterion_max_ci_width": STABLE_CRITERION_MAX_WIDTH,
        "criterion_source": "S1 Text: 'All 35 types classified stable (95% CI width <= 10; "
                            "median width 3, maximum 7)'",
        "n_types": int(len(t)),
        "n_stable": int((w <= STABLE_CRITERION_MAX_WIDTH).sum()),
        "ci_width_median": float(w.median()),
        "ci_width_max": float(w.max()),
        "ci_width_min": float(w.min()),
        "categories": {k: int(v) for k, v in t["bootstrap_category"].value_counts().items()},
    }

    OUT.write_text(json.dumps(result, indent=2) + "\n")
    p = result["primate_replication"]["pairs"]
    for k, v in p.items():
        print(f"  {k:<18} n={v['n_types']}  min={v['min_cells_either_arm']} "
              f"max={v['max_cells_either_arm']}  ({v['min_cells_type']})")
    m = result["mouse_lemur"]
    print(f"  mouse-lemur: orthologs={m['n_ortholog_pairs']} min_cells={m['min_cells_per_type']} "
          f"cap={m['max_cells_per_type']} nperm={m['n_permutations']}")
    b = result["bootstrap_rank_stability"]
    print(f"  bootstrap: {b['n_stable']}/{b['n_types']} stable at width <= "
          f"{b['criterion_max_ci_width']}; median {b['ci_width_median']} max {b['ci_width_max']}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
