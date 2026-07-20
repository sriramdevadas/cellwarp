#!/usr/bin/env python3
"""Build the submission packet (R21).

Materializes mirrored packet copies from their canonical sources. Subsumes
the earlier scripts/55_materialize_symlinks_for_submission.py (now archived)
and lifts canonical-to-mirror copy operations out of the per-figure
producer scripts (regenerate_figure_s4a.py, generate_phase3_figures.py,
46_synthesis_pass_supplementary_table_edits.py, 49_build_figS7_matched_scale.py).

Manifest: 30 rules covering 30 pairs.
  - Group A: 6 main-figure packet mirrors (figures/main/*.pdf → docs/...Figure_*.pdf)
  - Group B: 4 supp-figure packet mirrors (figures/submission/supplementary/figS*.pdf
             → docs/...Figure_S*.pdf)
  - Group C: 6 table packet mirrors (docs/supplementary_materials/table_S*.{csv,xlsx}
             → docs/submission/figures_for_review/Table_S*.{csv,xlsx})
  - Group D: 4 producer-tree / panel-promotion mirrors:
             - figS4_matched_scale_control.png canonical (submission/supp) →
               legacy mirror
             - figS4_matched_scale_control.pdf canonical (submission/supp) →
               legacy mirror
             - output/cancer/scaled/cross_analysis_scaled.png → figures/panels/suppl_text_s1_cancer.png
             - output/disease_replication/covid/covid_cross_analysis.png →
               figures/panels/suppl_text_s1_covid.png
  - Group E: gene-std / marker-null consolidation (figS5, table_S9/S10)
  - Group F: conserved-contribution integration (Figure 7, Table S11, pre-registration)

Idempotency: rerunning produces no diff after a clean run.

Consistency test: tests/test_submission_packet_consistency.py parametrizes
over the same manifest; this script and that test must stay in sync.

CLI:
    python scripts/build_submission_packet.py            # default --verify
    python scripts/build_submission_packet.py --verify   # md5-check all pairs
    python scripts/build_submission_packet.py --dry-run  # enumerate mirrors that would change
    python scripts/build_submission_packet.py --rebuild  # cp canonical → mirror; verify post

Exit codes: 0 success, 1 any drift or missing canonical, 2 verification mismatch
after rebuild.
"""
from __future__ import annotations

import argparse
import hashlib
import shutil
import sys
from pathlib import Path
from typing import List, Optional

PROJECT = Path(__file__).resolve().parent.parent


# Manifest: list of (canonical_relpath, [mirror_relpath, ...]).
# Keep in sync with tests/test_submission_packet_consistency.py.
MATERIALIZATION_RULES: List[dict] = [
    # ─── Group A: Main figures ───
    {"canonical": "figures/main/fig1_global_coherence.pdf",
     "mirrors": ["docs/submission/figures_for_review/Figure_1.pdf"]},
    {"canonical": "figures/main/fig2_two_layer.pdf",
     "mirrors": ["docs/submission/figures_for_review/Figure_4.pdf"]},
    {"canonical": "figures/main/fig3_replication.pdf",
     "mirrors": ["docs/submission/figures_for_review/Figure_3.pdf"]},
    {"canonical": "figures/main/fig4_human_macaque.pdf",
     "mirrors": ["docs/submission/figures_for_review/Figure_5.pdf"]},
    {"canonical": "figures/main/fig5_rigidity_ranking.pdf",
     "mirrors": ["docs/submission/figures_for_review/Figure_6.pdf"]},
    {"canonical": "figures/main/fig6_l1000_nulls.pdf",
     "mirrors": ["docs/submission/figures_for_review/Figure_7.pdf"]},
    # ─── Group B: Supplementary figures ───
    {"canonical": "figures/submission/supplementary/figS1_pipeline_validation.pdf",
     "mirrors": ["docs/submission/figures_for_review/Figure_S1.pdf"]},
    {"canonical": "figures/submission/supplementary/figS2_parameter_protocol_sensitivity.pdf",
     "mirrors": ["docs/submission/figures_for_review/Figure_S2.pdf"]},
    {"canonical": "figures/submission/supplementary/figS3_bootstrap_rankings.pdf",
     "mirrors": ["docs/submission/figures_for_review/Figure_S3.pdf"]},
    {"canonical": "figures/submission/supplementary/figS4_matched_scale_control.pdf",
     "mirrors": ["docs/submission/figures_for_review/Figure_S4.pdf"]},
    # ─── Group C: Tables ───
    # Table_1.xlsx is canonical-self (no in-repo producer post-R21); pinned by
    # tests/test_submission_packet_consistency.py::test_table_1_lock_md5().
    {"canonical": "docs/supplementary_materials/table_S1.xlsx",
     "mirrors": ["docs/submission/figures_for_review/Table_S1.xlsx"]},
    {"canonical": "docs/supplementary_materials/table_S2.xlsx",
     "mirrors": ["docs/submission/figures_for_review/Table_S2.xlsx"]},
    {"canonical": "docs/supplementary_materials/table_S3.csv",
     "mirrors": ["docs/submission/figures_for_review/Table_S3.csv"]},
    {"canonical": "docs/supplementary_materials/table_S4.csv",
     "mirrors": ["docs/submission/figures_for_review/Table_S4.csv"]},
    {"canonical": "docs/supplementary_materials/table_S5.csv",
     "mirrors": ["docs/submission/figures_for_review/Table_S5.csv"]},
    {"canonical": "docs/supplementary_materials/Table_S6_CPC1_driver_genes.xlsx",
     "mirrors": ["docs/submission/figures_for_review/Table_S6.xlsx"]},
    {"canonical": "docs/supplementary_materials/table_S7_layer1_housekeeping_exclusion.csv",
     "mirrors": ["docs/submission/figures_for_review/Table_S7.csv"]},
    {"canonical": "docs/supplementary_materials/table_S8_marker_ortholog_retention.csv",
     "mirrors": ["docs/submission/figures_for_review/Table_S8.csv"]},
    {"canonical": "docs/supplementary_materials/table_S12_software_environment.csv",
     "mirrors": ["docs/submission/figures_for_review/Table_S12.csv"]},
    # ─── Group D: Producer-tree / panel-promotion mirrors ───
    # Fig S4 .png pair (A.5 Task 6: missing from R20 pin)
    {"canonical": "figures/submission/supplementary/figS4_matched_scale_control.png",
     "mirrors": ["figures/supplementary/figS4_matched_scale_control.png"]},
    # Fig S4 .pdf legacy mirror — moved from scripts/49 dual-tree write to here
    {"canonical": "figures/submission/supplementary/figS4_matched_scale_control.pdf",
     "mirrors": ["figures/supplementary/figS4_matched_scale_control.pdf"]},
    # Panel promotions — moved from scripts/generate_phase3_figures.py to here
    {"canonical": "output/cancer/scaled/cross_analysis_scaled.png",
     "mirrors": ["figures/panels/suppl_text_s1_cancer.png"]},
    {"canonical": "output/disease_replication/covid/covid_cross_analysis.png",
     "mirrors": ["figures/panels/suppl_text_s1_covid.png"]},
    # ─── Group E: Stage-5 consolidation artifacts (gene-std / marker-null) ───
    {"canonical": "figures/submission/supplementary/figS5_markernull.pdf",
     "mirrors": ["docs/submission/figures_for_review/Figure_S5.pdf"]},
    {"canonical": "docs/supplementary_materials/table_S9_genestd_standardization.csv",
     "mirrors": ["docs/submission/figures_for_review/Table_S9.csv"]},
    {"canonical": "docs/supplementary_materials/table_S9_schemeB_CPC1_markers.csv",
     "mirrors": ["docs/submission/figures_for_review/Table_S9_schemeB_CPC1_markers.csv"]},
    {"canonical": "docs/supplementary_materials/table_S10_markernull.csv",
     "mirrors": ["docs/submission/figures_for_review/Table_S10.csv"]},
    # ─── Group F: Conserved-contribution (Figure 7, Table S11, pre-registration) ───
    {"canonical": "figures/main/fig7_conserved_contribution.pdf",
     "mirrors": ["docs/submission/figures_for_review/Figure_2.pdf"]},
    {"canonical": "docs/supplementary_materials/table_S11_gene_conservation.csv",
     "mirrors": ["docs/submission/figures_for_review/Table_S11.csv"]},
    {"canonical": "docs/preregistration_conserved_contribution_2026-06-05.md",
     "mirrors": ["docs/submission/figures_for_review/Supplementary_Preregistration.md"]},
]


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def resolve(rel: str) -> Path:
    return PROJECT / rel


def verify(rules: List[dict]) -> int:
    """Verify all (canonical, mirror) pairs are byte-equal. Return nonzero on any mismatch."""
    failures = []
    pair_count = 0
    for rule in rules:
        canonical = resolve(rule["canonical"])
        if not canonical.exists():
            failures.append(f"CANONICAL MISSING: {rule['canonical']}")
            continue
        c_md5 = md5(canonical)
        for mirror_rel in rule["mirrors"]:
            pair_count += 1
            mirror = resolve(mirror_rel)
            if not mirror.exists():
                failures.append(f"MIRROR MISSING: {mirror_rel}")
                continue
            m_md5 = md5(mirror)
            if c_md5 != m_md5:
                failures.append(
                    f"DRIFT: canonical={rule['canonical']} ({c_md5}) "
                    f"!= mirror={mirror_rel} ({m_md5})"
                )
    if failures:
        print(f"VERIFY FAIL: {len(failures)} / {pair_count} pairs", file=sys.stderr)
        for f in failures:
            print(f"  {f}", file=sys.stderr)
        return 1
    print(f"VERIFY OK: {pair_count} / {pair_count} pairs match")
    return 0


def dry_run(rules: List[dict]) -> int:
    """Enumerate mirrors that would change on rebuild. No writes."""
    would_change = []
    would_create = []
    pair_count = 0
    for rule in rules:
        canonical = resolve(rule["canonical"])
        if not canonical.exists():
            print(f"  SKIP (canonical missing): {rule['canonical']}", file=sys.stderr)
            continue
        c_md5 = md5(canonical)
        for mirror_rel in rule["mirrors"]:
            pair_count += 1
            mirror = resolve(mirror_rel)
            if not mirror.exists():
                would_create.append((rule["canonical"], mirror_rel))
            else:
                if md5(mirror) != c_md5:
                    would_change.append((rule["canonical"], mirror_rel))
    print(f"DRY-RUN: {pair_count} pairs evaluated")
    print(f"  would create:  {len(would_create)}")
    print(f"  would change:  {len(would_change)}")
    print(f"  in sync:       {pair_count - len(would_create) - len(would_change)}")
    for c, m in would_create:
        print(f"  CREATE  {m}   ← {c}")
    for c, m in would_change:
        print(f"  UPDATE  {m}   ← {c}")
    return 0


def rebuild(rules: List[dict]) -> int:
    """For each rule: ensure mirrors are byte-equal to canonical via shutil.copy2.
    Re-verifies post-rebuild. Idempotent."""
    n_writes = 0
    pair_count = 0
    for rule in rules:
        canonical = resolve(rule["canonical"])
        if not canonical.exists():
            print(f"FAIL: canonical missing: {rule['canonical']}", file=sys.stderr)
            return 1
        c_md5 = md5(canonical)
        for mirror_rel in rule["mirrors"]:
            pair_count += 1
            mirror = resolve(mirror_rel)
            mirror.parent.mkdir(parents=True, exist_ok=True)
            need_write = (not mirror.exists()) or md5(mirror) != c_md5
            if need_write:
                shutil.copy2(canonical, mirror)
                n_writes += 1
                print(f"  WROTE   {mirror_rel}   ← {rule['canonical']}")
            # post-write md5 verify
            if md5(mirror) != c_md5:
                print(f"FAIL: post-write md5 mismatch: {mirror_rel}", file=sys.stderr)
                return 2
    print(f"REBUILD OK: {n_writes} writes / {pair_count} pairs (rest already in sync)")
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Materialize and/or verify submission-packet mirrors."
    )
    g = parser.add_mutually_exclusive_group()
    g.add_argument("--dry-run", action="store_true", help="enumerate mirrors that would change; no writes")
    g.add_argument("--verify", action="store_true", help="verify all (canonical, mirror) pairs are byte-equal (default)")
    g.add_argument("--rebuild", action="store_true", help="cp canonical → mirror for any drifted/missing mirror; verify after")
    args = parser.parse_args(argv)

    rules = MATERIALIZATION_RULES
    if args.dry_run:
        return dry_run(rules)
    if args.rebuild:
        return rebuild(rules)
    # default = verify
    return verify(rules)


if __name__ == "__main__":
    sys.exit(main())
