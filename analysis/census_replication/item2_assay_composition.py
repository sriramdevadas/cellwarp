#!/usr/bin/env python3
"""
Assay composition of the primary human population.

The deposited protocol_breakdown.csv is mouse-only, so the repository alone cannot
settle what fraction of the primary human atlas is 10x Chromium whole-cell data.
This queries Census for the same population the pipeline downloads, building the
filter with the published cellwarp.data_loader.build_obs_value_filter rather than a
reimplementation, so the population queried is the population the pipeline uses.

Produces the assay denominators reported in the Methods data description and in
S2 Text's opening characterisation of the primary human atlas.

NOT TRACKED: requires the cellxgene_census package, which is a declared optional
extra rather than part of the documented gate environment, and reaches the network.
The two deposited sensitivity artifacts it cross-checks against are tracked.

Only the two path constants differ from the version that produced those values; the
code is unchanged.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import cellxgene_census

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
OUT = HERE
sys.path.insert(0, str(REPO / "src"))
from cellwarp.data_loader import (  # noqa: E402
    HUMAN_ORGANISM, MOUSE_ORGANISM, HUMAN_COLLECTION, MOUSE_COLLECTION,
    CELL_TYPE_MAP, MAX_CELLS_PER_TYPE,
    get_dataset_ids_for_collection, build_obs_value_filter,
)

CENSUS_VERSION = "2025-11-08"          # literal, as pull_aggregate.py:23 pins
R = {"census_version": CENSUS_VERSION, "max_cells_per_type": MAX_CELLS_PER_TYPE}

TYPES35 = list(pd.read_csv(REPO / "output/phase2/scaled_35types/centroids_human_35.csv",
                           index_col=0).index)
TYPES6 = sorted({c for v in CELL_TYPE_MAP.values() for c in v})
print(f"primary 35-type list : {len(TYPES35)} types")
print(f"CELL_TYPE_MAP list   : {len(TYPES6)} types -> {TYPES6}")


def proto_of(assay: str) -> str:
    """pull_aggregate.py:40-41, verbatim."""
    return "10x" if "10x" in assay else "SS"


with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
    print(f"\nopened Census {CENSUS_VERSION}")

    for tag, organism, collection in (("human", HUMAN_ORGANISM, HUMAN_COLLECTION),
                                      ("mouse", MOUSE_ORGANISM, MOUSE_COLLECTION)):
        print(f"\n{'='*78}\n{tag.upper()}  ({collection})\n{'='*78}")
        ds = get_dataset_ids_for_collection(census, collection)

        # published filter string, printed so the population is auditable
        vf_pub = build_obs_value_filter(TYPES35, ds)
        print(f"\nPUBLISHED FILTER (build_obs_value_filter, 35 types), first 220 chars:")
        print(f"  {vf_pub[:220]}...")

        # ---- one query on the dataset restriction only, so attrition is measurable ----
        obs = cellxgene_census.get_obs(
            census, organism,
            value_filter=" or ".join([f"dataset_id == '{d}'" for d in ds]),
            column_names=["cell_type", "assay", "is_primary_data", "disease"],
        )
        n0 = len(obs)
        prim = obs["is_primary_data"].astype(bool)
        norm = obs["disease"].astype(str) == "normal"
        in35 = obs["cell_type"].astype(str).isin(set(TYPES35))
        in6 = obs["cell_type"].astype(str).isin(set(TYPES6))

        print(f"\n-- FILTER ATTRITION (note 21: what each clause removes) --")
        print(f"  cells in the {len(ds)} collection datasets          : {n0:>10,}")
        s1 = int(prim.sum())
        print(f"  after is_primary_data == True                : {s1:>10,}   (-{n0-s1:,})")
        s2 = int((prim & norm).sum())
        print(f"  after disease == 'normal'                    : {s2:>10,}   (-{s1-s2:,})")
        s3 = int((prim & norm & in35).sum())
        print(f"  after cell_type in [35 primary types]        : {s3:>10,}   (-{s2-s3:,})")
        s6 = int((prim & norm & in6).sum())
        print(f"  after cell_type in [6 CELL_TYPE_MAP types]   : {s6:>10,}")

        # POSITIVE CONTROL on the reconstruction: run the published filter itself
        obs_pub = cellxgene_census.get_obs(census, organism, value_filter=vf_pub,
                                           column_names=["cell_type"])
        print(f"\n  POSITIVE CONTROL -- published filter executed directly: {len(obs_pub):,}")
        print(f"  python-side reconstruction                            : {s3:,}")
        print(f"  -> {'MATCH' if len(obs_pub) == s3 else 'MISMATCH -- reconstruction is wrong'}")

        sel = obs[prim & norm & in35]

        # ---- assay distribution ----
        print(f"\n-- ASSAY DISTRIBUTION over the 35-type primary population --")
        vc = sel["assay"].astype(str).value_counts()
        for a, n in vc.items():
            print(f"  {a:28s} {n:>10,}   {100*n/len(sel):6.2f}%   -> {proto_of(a)}")
        p10 = int(sum(n for a, n in vc.items() if proto_of(a) == "10x"))
        pss = int(len(sel) - p10)
        print(f"  {'TOTAL':28s} {len(sel):>10,}")
        print(f"  collapsed: 10x {p10:,} ({100*p10/len(sel):.2f}%) | "
              f"non-10x {pss:,} ({100*pss/len(sel):.2f}%)")

        sel6 = obs[prim & norm & in6]
        vc6 = sel6["assay"].astype(str).value_counts()
        p10_6 = int(sum(n for a, n in vc6.items() if proto_of(a) == "10x"))
        print(f"\n-- same, restricted to the 6 CELL_TYPE_MAP types --")
        for a, n in vc6.items():
            print(f"  {a:28s} {n:>10,}   {100*n/len(sel6):6.2f}%")
        print(f"  collapsed: 10x {p10_6:,} ({100*p10_6/len(sel6):.2f}%) | "
              f"non-10x {len(sel6)-p10_6:,} ({100*(len(sel6)-p10_6)/len(sel6):.2f}%)")

        # ---- per-cell-type assay counts ----
        tab = (sel.assign(proto=sel["assay"].astype(str).map(proto_of))
                  .groupby(["cell_type", "proto"]).size().unstack(fill_value=0))
        for c in ("10x", "SS"):
            if c not in tab.columns:
                tab[c] = 0
        tab = tab[["10x", "SS"]]
        tab["n_total"] = tab.sum(1)
        tab["frac_SS"] = tab["SS"] / tab["n_total"]
        tab = tab.reindex(TYPES35)
        print(f"\n-- PER-CELL-TYPE assay counts (35 types, pre-subsample) --")
        print(f"  {'cell_type':45s} {'10x':>9s} {'SS':>9s} {'total':>9s} {'frac_SS':>8s}")
        for t, row in tab.iterrows():
            print(f"  {str(t)[:45]:45s} {int(row['10x']):>9,} {int(row['SS']):>9,} "
                  f"{int(row['n_total']):>9,} {row['frac_SS']:>8.4f}")

        R[tag] = dict(
            n_datasets=len(ds), n_cells_in_datasets=n0,
            after_primary=s1, after_disease=s2, after_types35=s3, after_types6=s6,
            published_filter_count=int(len(obs_pub)),
            reconstruction_matches=bool(len(obs_pub) == s3),
            assay_counts={str(a): int(n) for a, n in vc.items()},
            n_10x=p10, n_nonpaper10x=pss,
            frac_10x=float(p10 / len(sel)), frac_non10x=float(pss / len(sel)),
            assay_counts_6type={str(a): int(n) for a, n in vc6.items()},
            frac_10x_6type=float(p10_6 / len(sel6)),
            per_type={str(t): dict(n_10x=int(r["10x"]), n_SS=int(r["SS"]),
                                   n_total=int(r["n_total"]), frac_SS=float(r["frac_SS"]))
                      for t, r in tab.iterrows()},
        )

# ---- MOUSE POSITIVE CONTROL vs the deposited post-subsample breakdown ----
print(f"\n{'='*78}\nMOUSE POSITIVE CONTROL vs deposited protocol_breakdown.csv\n{'='*78}")
dep = pd.read_csv(REPO / "output/phase2/sensitivity/smartseq2/protocol_breakdown.csv")
sens = json.load(open(REPO / "output/phase2/sensitivity/smartseq2/sensitivity_results.json"))
fresh = pd.DataFrame(R["mouse"]["per_type"]).T.reset_index().rename(columns={"index": "cell_type"})
m = dep.merge(fresh, on="cell_type", how="inner", suffixes=("_dep", "_fresh"))
print(f"types matched: {len(m)} / {len(dep)}")
from scipy import stats as st
rho, p = st.spearmanr(m["fraction_smartseq2"], m["frac_SS"])
pear = float(np.corrcoef(m["fraction_smartseq2"], m["frac_SS"])[0, 1])
print(f"\ndeposited total {sens['mouse_cells']['total']:,}, SS2 {sens['mouse_cells']['smartseq2']:,} "
      f"= {sens['mouse_cells']['smartseq2_fraction']*100:.2f}%   (POST-subsample, cap {MAX_CELLS_PER_TYPE})")
print(f"fresh     total {int(fresh['n_total'].sum()):,}, SS2 {int(fresh['n_SS'].sum()):,} "
      f"= {100*fresh['n_SS'].sum()/fresh['n_total'].sum():.2f}%   (PRE-subsample)")
print(f"\nSHAPE CHECK -- per-type Smart-seq2 fraction, deposited vs fresh:")
print(f"  Spearman rho = {rho:.4f}  p = {p:.3e}   Pearson r = {pear:.4f}   n = {len(m)}")
m["abs_diff"] = (m["fraction_smartseq2"] - m["frac_SS"]).abs()
print(f"  |diff|: median {m['abs_diff'].median():.4f}  max {m['abs_diff'].max():.4f}")
print(f"\n  largest disagreements:")
for _, r in m.nlargest(6, "abs_diff").iterrows():
    print(f"    {r['cell_type'][:44]:44s} dep {r['fraction_smartseq2']:.4f}  "
          f"fresh {r['frac_SS']:.4f}  diff {r['fraction_smartseq2']-r['frac_SS']:+.4f}")
R["mouse_positive_control"] = dict(
    n_types_matched=int(len(m)), spearman_rho=float(rho), spearman_p=float(p),
    pearson_r=pear, median_abs_diff=float(m["abs_diff"].median()),
    max_abs_diff=float(m["abs_diff"].max()),
    deposited_total=int(sens["mouse_cells"]["total"]),
    deposited_ss2=int(sens["mouse_cells"]["smartseq2"]),
    deposited_ss2_frac=float(sens["mouse_cells"]["smartseq2_fraction"]),
    fresh_total=int(fresh["n_total"].sum()), fresh_ss2=int(fresh["n_SS"].sum()),
    fresh_ss2_frac=float(fresh["n_SS"].sum() / fresh["n_total"].sum()),
)

(OUT / "item2_assay_results.json").write_text(json.dumps(R, indent=2))
print(f"\nwrote {OUT/'item2_assay_results.json'}")
