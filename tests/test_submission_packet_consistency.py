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
TABLE_1_LOCK_MD5 = "9e00ccf02e6102f9bc0058defb79b18a"  # D68 re-lock (membership rule): the footnote asserted that the nine excluded IDs "are descriptive, aggregate, or diagnostic quantities that make no inferential significance claim", which is true of four of them. T03, T04, T09, T10 and T18 are excluded because they report a pass rate over repeated sub-tests rather than one p-value for one hypothesis -- their Raw p cells read "100/100 sig", "35/35 < 1.0", "all p < 0.01 (100/100)" twice and "20/24 sig" -- and six of the nine (T03, T04, T09, T10, T18, T62) are Confirmatory, so the old sentence also implied a status the Status column contradicts. Replaced with the rule itself, so the nine are derivable rather than listed, plus a second sentence naming the four rows that are exceptions in the other direction (T11, T44, T47, T48). Diff against the prior lock is cell A66 and nothing else: 0 other cell values, 0 fills, 0 styles, 0 row heights, 0 column widths, merged ranges identical (A66:K69 pair unchanged), dimensions A1:K69 both sides, one sheet both sides, 64 numbered rows both sides. Verified unchanged across the re-lock: column K (which tests/test_table1_callouts.py asserts), every Raw p, every Bonferroni corrected p (all 55 in-family values), every Status, and the k=55 header. Two writer changes were needed together and neither works alone. (a) The new block's sentinel is "Membership rule:", a string that exists only in the text it writes; keyed on anything already in the tracked footnote it would have been a silent no-op, reporting success with the bytes untouched, this lock still matching and all four gates green -- the exact failure the D61 note at the T34/T35 block records as having nearly shipped. (b) The D57 block's guard was keyed on "T44, T47 and T48 also report" and the new F2 text drops the "also"; left alone it would have stopped matching, fired on the NEXT run, and re-appended the superseded sentence after the "given in that column instead." anchor F2 keeps, so the script would have ceased to be a fixed point. Re-keyed on the substring both wordings share. Confirmed by running the writer three times: runs 2 and 3 are byte-identical to run 1, and "Membership rule:", "T44, T47 and T48", "given in that column instead." and "Four rows are exceptions" each occur exactly once. Prior D67 lock aaa8760c (sheet name): the worksheet tab was "Table 1" while the file, the manuscript caption (manuscript_combined.txt:330 "S13 Table."), the legend and this workbook's own footnote all say S13 Table -- and the tab is the first thing a reader opening the workbook sees. Renamed to "S13 Table" by table1_formatting.py, keyed on the old name so it runs against an unmigrated copy and is a fixed point once renamed. Diff against the prior lock is the sheet name and nothing else: 0 cell values, 0 fills, 0 styles, 0 row heights, 0 column widths, merged ranges identical (A66:K66, A69:K69), dimensions A1:K69 both sides, one sheet both sides. Dependents updated in the same commit: table1_formatting.py's two selectors and tests/test_table1_callouts.py's SHEET constant, which were the only tracked references that select the sheet by name. Prior D61 lock fbd40534 (statuses): T42-T47, the six mechanistic nulls S1_Text.txt:60 records as "pre-specified as confound diagnostics", move Exploratory -> Confirmatory, which is what the footnote's own definition of Confirmatory (pre-specified or direct replication) makes them. T48-T51 stay Exploratory: a closure threshold locked before invocation pre-commits the decision rule, not the hypothesis. k does not move -- family membership is the em-dash in the Bonferroni column and none of T42-T47 is among the nine excluded. Fills are now DERIVED from Status rather than written beside it, so no row can contradict the legend again; the derivation reproduces the colour of all 63 rows that already agreed and fixes T68, which was green on an Exploratory row because D57 copied a template style when appending T67-T69. Diff vs the prior lock: 6 cell values (J42-J47), 77 fgColor cells (7 rows x 11), 33 bgColor cells (T09/T10/T11 normalised 00E2EFDA -> 00000000, ignored under a solid pattern); merges, row heights and column widths unchanged. T68's Status is deliberately left as it stands: no submitted text assigns a design status to the marker-distinctness correlation, and the D57 literal that made T67 Confirmatory and T68 Exploratory recorded no ground for the difference. Prior D61 lock 962d83ae: footnote only, one cell (A66), no row, value, fill, merge, height or width changes. (a) the T34/T35 sentence claimed the per-type top-50 pass-rate "is reported as" 5 of 6 / 6 of 6 and cited Figure S6; the result appears in none of the three submitted texts and old Fig S6 was cut (so annotated on reproduce/validate.py's three CellMarker checks), which made S13 the only place in the submission carrying those numbers -- restated as a pass-rate that is not reported in this paper. (b) "Table S4" -> "S4 Table", the last inverted supplementary reference in the workbook; D56 normalised column K but LEGACY_SUPP_RE is whole-cell-anchored and could not reach footnote prose, and no block of table1_formatting.py had ever authored that sentence. Note also that the D22 block was guarded by `if "T34/T35:" not in _fv`, which the tracked footnote already satisfied: correcting the literal alone would have been a silent no-op with all four gates green, so the block is now keyed on the superseded sentence. This is the first lock whose value IS the script's own output -- the trailer below was false at D57 and earlier, when the tracked bytes had never come from a run. Prior D57 lock 7909c020: append T67-T69 for three tests the paper reports but the sheet omitted (marker-similarity null at K=15, per-type residual vs marker-distinctness, PanSci Layer-2 pre-rotation); inferential family 52 -> 55 with every corrected p recomputed as the exact product and capped at 1; header k=52 -> k=55; footnote restated for 64 tests. No row changes verdict. Prior D56 lock 319d8c3c: normalise 18 supplementary references to the number-first order the submitted texts use (S1 Fig A, S4 Table); the inverted Fig S1 / Table S1 form occurs nowhere in the three submitted files. Column K only; no value, row or callout target changes. Prior D55 lock b446569a: repair the Figure/Table column against the five-figure paper (38 of 61 rows), T64 ceiling 0.42 -> 0.45 with the calibrated-signal qualifier, T52 n 34 -> 35 from output/t3g/primary_correlation_results.json, and one footnote sentence stating that the nine excluded IDs are the em-dash rows in the Bonferroni column; the column is now also validated by tests/test_table1_callouts.py, which derives the valid display items from the manuscript rather than freezing a list. Prior D28 lock 61e6988f (T11 corrected-p cell "significant" -> "direction-robust (100/100; resampling)", Review 6 Minor 1); prior D22 lock d375352f (reconcile T34/T35 to manuscript Table 1). Reproducible byte-for-byte via scripts/table1_formatting.py


def test_table_1_lock_md5() -> None:
    """Lock md5 of Table_1.xlsx pending 5b authoring of 44b script.
    Pinned R20 close value; protects against silent drift while no in-repo
    producer materializes this canonical."""
    path = REPO_ROOT / "docs/supplementary_materials/table_S13_test_inventory.xlsx"
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
