"""R20 architectural alpha: packet consistency.

For every file in docs/submission/figures_for_review/, assert that its
md5 matches the canonical counterpart it was copied from. Pins prevent
stale-packet drift of the kind that shipped a pre-R15-cascade Table 1
to Cell Systems in R20 Part A.

When canonical content updates, update the pinned canonical path here
and re-run the packet replication (typically in the round that updates
the content).

R21 extension: pin set updated to match the 22-rule manifest in
scripts/build_submission_packet.py. Figure_S4 canonical re-pointed to
the producer's true canonical (figS5_polished). Table_1.xlsx pin
removed (no in-repo producer post-R21 cascade revert; protected by
test_table_1_lock_md5 below). Added Group D pins for the figS7 .png
pair and the figS7 .pdf legacy mirror.
"""

import hashlib
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent

# Each entry: packet_relpath -> canonical_relpath
# Canonical is the source of truth; packet is a copy that must match.
# NOT in sync with scripts/build_submission_packet.py MATERIALIZATION_RULES, and
# deliberately not resynchronised. That manifest holds 30 (canonical, mirror)
# pairs; this map holds 26. The four omitted are its Group E rules:
#   Figure_S5.pdf                     <- figures/submission/supplementary/figS5_markernull.pdf
#   Table_S9.csv                      <- docs/supplementary_materials/table_S9_genestd_standardization.csv
#   Table_S9_schemeB_CPC1_markers.csv <- docs/supplementary_materials/table_S9_schemeB_CPC1_markers.csv
#   Table_S10.csv                     <- docs/supplementary_materials/table_S10_markernull.csv
# This map is a strict subset: nothing here contradicts the manifest. The packet
# pins the superseded 7-figure layout and is retained, not maintained, so the
# stricter check is build_submission_packet.py --verify, which covers all 30.
PACKET_CANONICAL_MAP = {
    "docs/submission/figures_for_review/Figure_1.pdf":
        "figures/main/fig1_global_coherence.pdf",
    "docs/submission/figures_for_review/Figure_4.pdf":
        "figures/main/fig2_two_layer.pdf",
    "docs/submission/figures_for_review/Figure_3.pdf":
        "figures/main/fig3_replication.pdf",
    "docs/submission/figures_for_review/Figure_5.pdf":
        "figures/main/fig4_human_macaque.pdf",
    "docs/submission/figures_for_review/Figure_6.pdf":
        "figures/main/fig5_rigidity_ranking.pdf",
    "docs/submission/figures_for_review/Figure_7.pdf":
        "figures/main/fig6_l1000_nulls.pdf",
    "docs/submission/figures_for_review/Figure_S1.pdf":
        "figures/submission/supplementary/figS1_pipeline_validation.pdf",
    "docs/submission/figures_for_review/Figure_S2.pdf":
        "figures/submission/supplementary/figS2_parameter_protocol_sensitivity.pdf",
    "docs/submission/figures_for_review/Figure_S3.pdf":
        "figures/submission/supplementary/figS3_bootstrap_rankings.pdf",
    "docs/submission/figures_for_review/Figure_S4.pdf":
        "figures/submission/supplementary/figS4_matched_scale_control.pdf",
    # Fig S4 .png pair (R21 A.5 Task 6): canonical at submission/supp side.
    "figures/supplementary/figS4_matched_scale_control.png":
        "figures/submission/supplementary/figS4_matched_scale_control.png",
    # Fig S4 .pdf legacy mirror (R21 B.2: moved from scripts/49 dual-tree).
    "figures/supplementary/figS4_matched_scale_control.pdf":
        "figures/submission/supplementary/figS4_matched_scale_control.pdf",
    # Table_1.xlsx pin REMOVED (R21 B.2): no in-repo producer post-R21 cascade
    # revert; protected by test_table_1_lock_md5() below.
    "docs/submission/figures_for_review/Table_S1.xlsx":
        "docs/supplementary_materials/table_S1.xlsx",
    "docs/submission/figures_for_review/Table_S2.xlsx":
        "docs/supplementary_materials/table_S2.xlsx",
    "docs/submission/figures_for_review/Table_S3.csv":
        "docs/supplementary_materials/table_S3.csv",
    "docs/submission/figures_for_review/Table_S4.csv":
        "docs/supplementary_materials/table_S4.csv",
    "docs/submission/figures_for_review/Table_S5.csv":
        "docs/supplementary_materials/table_S5.csv",
    "docs/submission/figures_for_review/Table_S6.xlsx":
        "docs/supplementary_materials/Table_S6_CPC1_driver_genes.xlsx",
    "docs/submission/figures_for_review/Table_S7.csv":
        "docs/supplementary_materials/table_S7_layer1_housekeeping_exclusion.csv",
    "docs/submission/figures_for_review/Table_S8.csv":
        "docs/supplementary_materials/table_S8_marker_ortholog_retention.csv",
    # Panel-promotion mirrors (R21 B.2: moved from scripts/generate_phase3_figures.py).
    "figures/panels/suppl_text_s1_cancer.png":
        "output/cancer/scaled/cross_analysis_scaled.png",
    "figures/panels/suppl_text_s1_covid.png":
        "output/disease_replication/covid/covid_cross_analysis.png",
    # Conserved-contribution (Stage-2 integration: Figure 7, Table S11, pre-registration).
    "docs/submission/figures_for_review/Figure_2.pdf":
        "figures/main/fig7_conserved_contribution.pdf",
    "docs/submission/figures_for_review/Table_S11.csv":
        "docs/supplementary_materials/table_S11_gene_conservation.csv",
    "docs/submission/figures_for_review/Supplementary_Preregistration.md":
        "docs/preregistration_conserved_contribution_2026-06-05.md",
    "docs/submission/figures_for_review/Table_S12.csv":
        "docs/supplementary_materials/table_S12_software_environment.csv",
}


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


@pytest.mark.parametrize(
    "packet_rel,canonical_rel",
    sorted(PACKET_CANONICAL_MAP.items()),
    ids=lambda p: Path(p).name if isinstance(p, str) else str(p),
)
def test_packet_matches_canonical(packet_rel: str, canonical_rel: str) -> None:
    packet = REPO_ROOT / packet_rel
    canonical = REPO_ROOT / canonical_rel
    assert packet.exists(), f"Packet file missing: {packet_rel}"
    assert canonical.exists(), f"Canonical file missing: {canonical_rel}"
    packet_md5 = _md5(packet)
    canonical_md5 = _md5(canonical)
    assert packet_md5 == canonical_md5, (
        f"Packet drift: {packet_rel} md5={packet_md5} "
        f"!= canonical {canonical_rel} md5={canonical_md5}. "
        f"Re-run packet replication (scripts/build_submission_packet.py --rebuild)."
    )


# R20 close lock for Table_1.xlsx. No in-repo script reproduces this canonical
# (Table-1 producer scripts archived in R21; see scripts/archive/README.md).
# 5b follow-up: scripts/44b_k53_table1_update.py will close the reproducibility
# gap and re-establish a canonical→mirror chain at that point.
TABLE_1_LOCK_MD5 = "61e6988f058987dbb0118b20dde16bfd"  # D28 re-lock: T11 corrected-p cell "significant" -> "direction-robust (100/100; resampling)" (Review 6 Minor 1); prior D22 lock d375352f (reconcile T34/T35 to manuscript Table 1: 500-gene identity-set labels + n + caption note); reproducible byte-for-byte via scripts/table1_formatting.py


def test_table_1_lock_md5() -> None:
    """Lock md5 of Table_1.xlsx pending 5b authoring of 44b script.
    Pinned R20 close value; protects against silent drift while no in-repo
    producer materializes this canonical."""
    path = REPO_ROOT / "docs/submission/figures_for_review/Table_1.xlsx"
    assert path.exists(), f"Table_1.xlsx missing at {path}"
    digest = _md5(path)
    assert digest == TABLE_1_LOCK_MD5, (
        f"Table_1.xlsx md5 drift: got {digest}, expected R20 lock {TABLE_1_LOCK_MD5}"
    )


def test_md5_16_byte_identical() -> None:
    """MD5-16 LEGITIMATE_SHARED_FIXTURE byte-identity guard. If both null
    files exist, they should be byte-identical (deterministic shared
    infrastructure: scripts 15 and 16 compute the same Tabula-h-vs-m
    Procrustes null on the same 6-type restricted subset with the same
    RANDOM_SEED=42). See R21 MD5-16 investigation report."""
    a = REPO_ROOT / "output/validation/hca_centroid_comparison/null_b_tabula_hvm.npy"
    b = REPO_ROOT / "output/validation/sun2023_replication/null_b_tabula.npy"
    if not (a.exists() and b.exists()):
        pytest.skip(f"MD5-16 source files not present (a={a.exists()}, b={b.exists()})")
    assert _md5(a) == _md5(b), (
        f"MD5-16 byte-identity drift: {_md5(a)} vs {_md5(b)}"
    )
