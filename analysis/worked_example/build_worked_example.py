#!/usr/bin/env python3
"""Worked example: the one cell type the main text carries through all three layers.

Emits the values Results section 4 states about hepatocyte -- its Layer-1 residual and
divergence rank, its within-atlas 95% CI width on that rank, its Layer-2 subspace
similarity at k = 5 on both sides of the centroid-optimal rotation, and its rank in each
cross-atlas replication -- into a JSON that reproduce/validate.py can read.

WHY A PRODUCER AND NOT A GATE ENTRY PER CSV

The four sources are CSVs, and validate.py cannot read one: its load_value has a .csv
branch, but that branch hands {"df": DataFrame} to a resolver that can only navigate
dicts and lists, so any key raises (reproduce/figure_script_map.md, Known gaps). Nothing
here is a new result. Every field is a cell of a tracked CSV, or a count derived from
the declared table below, with its source recorded beside it.

WHY THE CROSS-FILE ASSERTIONS ARE THE POINT

A single extracted row is only as trustworthy as the join that produced it, and the four
files were written by four different producers that never see each other. So before any
value is emitted this asserts, over all 35 types and not just the worked example, that:

  1. the four CSVs carry the same 35 cell_type strings;
  2. residuals_ranked.rank == bootstrap_summary.original_rank == master.primary_rank,
     row for row;
  3. bootstrap_summary.ci_width == master.bootstrap_CI_width, row for row;
  4. the two S columns' 35-type means reproduce the deposited
     summary_stats.json 35type.mean_alignment."k=5".{pre,post}, which is what Gate 1
     already checks and what Fig 2B draws.

Assertion 4 ties the per-type Layer-2 column to a value that is already gated, so the
extracted per-type S sits under the same statistic the paper reports rather than beside
it. If any of the four fails, nothing is written.

WHAT "INDEPENDENT COMPARISON" MEANS HERE

master_ranking_table.csv has seven replication columns, but they are not seven
comparisons. CellHint_rank and CellHint_harmonized_rank are the same source atlas at two
harmonization levels, so they are one comparison counted twice; and five of the seven
hold one arm of the primary fixed and substitute the other. Which is which cannot be
derived from the CSV -- it is read from the producers -- so it is declared in
REPLICATION_COLUMNS with each source path, and every count below is derived from that
table rather than written in. The declared column set is asserted against the columns the
CSV actually carries, so a replication added upstream cannot go silently uncounted.

Output (tracked):
  analysis/worked_example/worked_example_hepatocyte.json

Deterministic, no network, runtime well under a second. Not invoked by
reproduce/run_all.sh. Gated by reproduce/validate.py.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent.parent.parent
OUT = Path(__file__).resolve().parent / "worked_example_hepatocyte.json"

WORKED_EXAMPLE = "hepatocyte"
N_TYPES = 35
K = 5
TOL = 1e-12

RESIDUALS = "output/phase2/scaled_35types/residuals_ranked.csv"
BOOTSTRAP = "analysis/bootstrap_rankings/bootstrap_summary.csv"
ALIGNMENT = "output/mechanistic/ellipsoid_alignment/35type_alignment_scores.csv"
MASTER = "analysis/cross_reference/master_ranking_table.csv"
SUMMARY = "output/mechanistic/ellipsoid_alignment/summary_stats.json"

# The two rank-named columns of master_ranking_table.csv that are not replications.
NOT_A_REPLICATION = ("primary_rank", "bootstrap_median_rank")

# One entry per replication column. `comparison` is the identity two columns share when
# they are the same source data: CellHint appears twice, at two harmonization levels, and
# counts once. `substitutes` names the arm of the primary human-mouse pair that this
# comparison replaces; `holds_fixed` names the arm it keeps, and is None only when the
# comparison shares neither arm with the primary. Read from the producers named in
# `source`, not inferred from the table.
REPLICATION_COLUMNS = (
    {
        "column": "Sun2023_rank",
        "comparison": "Sun2023",
        "substitutes": "mouse",
        "holds_fixed": "human",
        "source": "output/validation/sun2023_replication_expanded/ranking_comparison.csv",
        "source_column": "sun2023_residual",
        "note": "Sun et al. 2023 young sedentary control mouse, 10x Chromium 3' v3, "
                "against the primary's Tabula Sapiens human arm "
                "(scripts/16_sun2023_replication.py)",
    },
    {
        "column": "PanSci_rank",
        "comparison": "PanSci",
        "substitutes": "mouse",
        "holds_fixed": "human",
        "source": "output/validation/pansci_replication/ranking_comparison.csv",
        "source_column": "pansci_residual",
        "note": "PanSci 6-month wild-type mouse, EasySci combinatorial-indexing snRNA-seq "
                "and not 10x, against the primary's human arm "
                "(scripts/pansci_replication.py)",
    },
    {
        "column": "CellHint_rank",
        "comparison": "CellHint",
        "substitutes": "human",
        "holds_fixed": "mouse",
        "source": "output/validation/cellhint_replication/ranking_comparison.csv",
        "source_column": "cellhint_residual",
        "note": "CellHint human meta-atlas, 9 tissues independent of Tabula Sapiens, "
                "against the primary's Tabula Muris mouse arm "
                "(scripts/33_cellhint_replication.py)",
    },
    {
        "column": "CellHint_harmonized_rank",
        "comparison": "CellHint",
        "substitutes": "human",
        "holds_fixed": "mouse",
        "source": "analysis/harmonized_replication/harmonized_residuals_cellhint.csv",
        "source_column": "residual_magnitude",
        "note": "THE SAME CellHint ATLAS as the row above, after progressive "
                "harmonization (ontology matching, tissue restriction, cell-count "
                "capping) and re-ranked on the restricted type set "
                "(analysis/harmonized_replication/harmonized_replication.py). A second "
                "reading of one comparison, not a second comparison",
    },
    {
        "column": "Pan_Census_rank",
        "comparison": "Pan_Census",
        "substitutes": "both",
        "holds_fixed": None,
        "source": "analysis/census_replication/ranking_comparison.csv",
        "source_column": "replication_residual",
        "note": "Both arms rebuilt from CELLxGENE Census collections that exclude Tabula "
                "Sapiens, Tabula Muris Senis, CellHint, PanSci and Sun2023 "
                "(analysis/census_replication/02_run_replication.py). The only column "
                "here that shares neither arm with the primary",
    },
    {
        "column": "Macaque_rank",
        "comparison": "Macaque",
        "substitutes": "mouse",
        "holds_fixed": "human",
        "source": "output/macaque_pipeline/reconstruction_qu12_results.json",
        "source_column": "per_type_residuals_ranked",
        "note": "Human-macaque over 12 types (Qu 2022), on the primary's human centroids "
                "(analysis/macaque/reconstruction_rira13_report.md)",
    },
    {
        "column": "Mouse_lemur_rank",
        "comparison": "Mouse_lemur",
        "substitutes": "mouse",
        "holds_fixed": "human",
        "source": "analysis/mouse_lemur/per_type_residuals.csv",
        "source_column": "residual_magnitude",
        "note": "Human-mouse lemur over 15 types, reading "
                "output/phase2/scaled_35types/centroids_human_35.csv for the human arm "
                "(analysis/mouse_lemur/01_run_pipeline.py:72)",
    },
)


def require(condition, message):
    if not condition:
        raise SystemExit("ASSERTION FAILED: " + message)


def main():
    resid = pd.read_csv(BASE / RESIDUALS)
    boot = pd.read_csv(BASE / BOOTSTRAP)
    align = pd.read_csv(BASE / ALIGNMENT)
    master = pd.read_csv(BASE / MASTER)
    summary = json.loads((BASE / SUMMARY).read_text())

    # ── 1. the same 35 cell_type strings in all four ──────────────────────────
    names = {
        RESIDUALS: set(resid.cell_type),
        BOOTSTRAP: set(boot.cell_type),
        ALIGNMENT: set(align.cell_type),
        MASTER: set(master.cell_type),
    }
    for path, got in names.items():
        require(len(got) == N_TYPES,
                "%s carries %d distinct cell types, expected %d" % (path, len(got), N_TYPES))
    reference = names[RESIDUALS]
    for path, got in names.items():
        require(got == reference,
                "%s does not carry the same cell_type strings as %s: only there %s, "
                "missing %s" % (path, RESIDUALS, sorted(got - reference),
                                sorted(reference - got)))
    require(WORKED_EXAMPLE in reference,
            "the worked example %r is not one of the %d types" % (WORKED_EXAMPLE, N_TYPES))

    # ── 2. the rank column is one column under three names ────────────────────
    ranks = (resid[["cell_type", "rank"]]
             .merge(boot[["cell_type", "original_rank"]], on="cell_type")
             .merge(master[["cell_type", "primary_rank"]], on="cell_type"))
    require(len(ranks) == N_TYPES, "rank join lost rows: %d of %d" % (len(ranks), N_TYPES))
    require((ranks["rank"] == ranks["original_rank"]).all(),
            "residuals_ranked.rank and bootstrap_summary.original_rank disagree at %s"
            % ranks.loc[ranks["rank"] != ranks["original_rank"], "cell_type"].tolist())
    require((ranks["rank"] == ranks["primary_rank"]).all(),
            "residuals_ranked.rank and master_ranking_table.primary_rank disagree at %s"
            % ranks.loc[ranks["rank"] != ranks["primary_rank"], "cell_type"].tolist())

    # ── 3. the CI-width column is one column under two names ──────────────────
    widths = boot[["cell_type", "ci_width"]].merge(
        master[["cell_type", "bootstrap_CI_width"]], on="cell_type")
    require(len(widths) == N_TYPES,
            "CI-width join lost rows: %d of %d" % (len(widths), N_TYPES))
    require((widths.ci_width == widths.bootstrap_CI_width).all(),
            "bootstrap_summary.ci_width and master_ranking_table.bootstrap_CI_width "
            "disagree at %s"
            % widths.loc[widths.ci_width != widths.bootstrap_CI_width, "cell_type"].tolist())

    # ── 4. the per-type S columns roll up to the gated aggregate ──────────────
    at_k = align[align.k == K]
    require(len(at_k) == N_TYPES,
            "alignment scores carry %d rows at k=%d, expected %d" % (len(at_k), K, N_TYPES))
    mean_pre = float(at_k.S_pre.mean())
    mean_post = float(at_k.S_post.mean())
    deposited = summary["35type"]["mean_alignment"]["k=%d" % K]
    require(abs(mean_pre - deposited["pre"]) < TOL,
            "mean S_pre at k=%d is %.17g, deposited %.17g" % (K, mean_pre, deposited["pre"]))
    require(abs(mean_post - deposited["post"]) < TOL,
            "mean S_post at k=%d is %.17g, deposited %.17g"
            % (K, mean_post, deposited["post"]))

    # ── the declared replication table matches the file ───────────────────────
    in_file = [c for c in master.columns
               if c.endswith("_rank") and c not in NOT_A_REPLICATION]
    declared = [d["column"] for d in REPLICATION_COLUMNS]
    require(sorted(in_file) == sorted(declared),
            "master_ranking_table.csv carries replication columns %s but this producer "
            "declares %s; a comparison added upstream would go uncounted"
            % (sorted(in_file), sorted(declared)))

    # ── extract the worked example ────────────────────────────────────────────
    r = resid[resid.cell_type == WORKED_EXAMPLE].iloc[0]
    b = boot[boot.cell_type == WORKED_EXAMPLE].iloc[0]
    a = at_k[at_k.cell_type == WORKED_EXAMPLE].iloc[0]
    m = master[master.cell_type == WORKED_EXAMPLE].iloc[0]

    present = []
    for d in REPLICATION_COLUMNS:
        value = m[d["column"]]
        if pd.isna(value):
            continue
        present.append({
            "column": d["column"],
            "comparison": d["comparison"],
            "rank": int(value),
            "substitutes": d["substitutes"],
            "holds_fixed": d["holds_fixed"],
            "source": d["source"],
            "source_column": d["source_column"],
            "note": d["note"],
        })

    require(len(present) == int(m.n_replications_present),
            "%s has %d non-null replication ranks but n_replications_present says %d"
            % (WORKED_EXAMPLE, len(present), int(m.n_replications_present)))

    comparisons = {p["comparison"] for p in present}
    first = [p for p in present if p["rank"] == 1]
    no_shared_arm = {p["comparison"] for p in present if p["holds_fixed"] is None}

    payload = {
        "analysis": "worked example carried through all three layers by Results section 4",
        "note": "no new result; every field is a cell of a tracked CSV or a count "
                "derived from this producer's declared replication table, with the "
                "source recorded beside it",
        "cell_type": WORKED_EXAMPLE,
        "n_types": N_TYPES,
        "layer1": {
            "residual_magnitude": float(r.residual_magnitude),
            "divergence_rank": int(r["rank"]),
            "rank_convention": "1 = most divergent of %d" % N_TYPES,
            "source": RESIDUALS,
        },
        "within_atlas": {
            "ci_width_ranks": float(b.ci_width),
            "ci_lower": float(b.ci_lower),
            "ci_upper": float(b.ci_upper),
            "category": str(b.category),
            "source": BOOTSTRAP,
            "source_note": "equal row for row to master_ranking_table.bootstrap_CI_width, "
                           "asserted above",
        },
        "layer2": {
            "k": K,
            "S_pre": float(a.S_pre),
            "S_post": float(a.S_post),
            "column_mean_S_pre": mean_pre,
            "column_mean_S_post": mean_post,
            "source": ALIGNMENT,
            "aggregate_source": SUMMARY,
            "aggregate_note": "the two column means reproduce "
                              "35type.mean_alignment.k=%d.{pre,post} exactly; they are "
                              "the two bars Fig 2B draws at k=%d" % (K, K),
        },
        "cross_atlas": {
            "mean_rank_shift": float(m.mean_rank_shift),
            "sd_rank_shift": float(m.sd_rank_shift),
            "shift_definition": "mean over replications of |primary_subset_rank - "
                                "replication_rank|, the primary re-ranked within each "
                                "replication's own matched subset "
                                "(analysis/cross_reference/cross_reference_analysis.py:120)",
            "per_replication_ranks": present,
            # The same ranks keyed by column name, so a gate can address one without a
            # list index that would silently follow a reordering of the table above.
            "ranks_by_column": {p["column"]: p["rank"] for p in present},
            "n_replication_columns_present": len(present),
            "n_distinct_comparisons": len(comparisons),
            "n_columns_ranking_first": len(first),
            "n_distinct_comparisons_ranking_first": len({p["comparison"] for p in first}),
            "n_comparisons_sharing_no_arm_with_primary": len(no_shared_arm),
            "counting_note": "columns are not comparisons: CellHint_rank and "
                             "CellHint_harmonized_rank are one source atlas at two "
                             "harmonization levels. Every count here is derived from "
                             "REPLICATION_COLUMNS in the producer, whose declared column "
                             "set is asserted against the file",
        },
    }

    OUT.write_text(json.dumps(payload, indent=2) + "\n")
    print("wrote", OUT.relative_to(BASE))
    print("  %s: rank %d of %d, residual %.4f, CI width %g ranks"
          % (WORKED_EXAMPLE, payload["layer1"]["divergence_rank"], N_TYPES,
             payload["layer1"]["residual_magnitude"],
             payload["within_atlas"]["ci_width_ranks"]))
    print("  Layer 2 at k=%d: S_pre %.4f, S_post %.4f (35-type means %.4f / %.4f)"
          % (K, payload["layer2"]["S_pre"], payload["layer2"]["S_post"],
             mean_pre, mean_post))
    print("  cross-atlas: %d columns present -> %d distinct comparisons; "
          "%d columns rank it first -> %d distinct comparisons; "
          "%d comparison shares neither arm with the primary"
          % (payload["cross_atlas"]["n_replication_columns_present"],
             payload["cross_atlas"]["n_distinct_comparisons"],
             payload["cross_atlas"]["n_columns_ranking_first"],
             payload["cross_atlas"]["n_distinct_comparisons_ranking_first"],
             payload["cross_atlas"]["n_comparisons_sharing_no_arm_with_primary"]))
    print("  four cross-file assertions passed over all %d types" % N_TYPES)


if __name__ == "__main__":
    main()
