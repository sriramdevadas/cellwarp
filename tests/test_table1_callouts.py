"""
Table 1's Figure/Table column must name display items the paper actually has.

WHAT THIS CHECKS, AND WHAT IT DOES NOT
--------------------------------------
This is a DANGLING-REFERENCE check, not a correctness check. It catches a cell
that names a figure or table which does not exist -- the failure that let Table 1
keep pointing at Fig 6B, Fig 7A, Fig 7B, Fig S6 and Fig S7 through four correction
commits after the paper was re-assembled from seven main figures to five.

It CANNOT catch a cell that names a display item which exists but shows something
else. Table 1 pointed T30 at Fig 4C for the Layer-2 result; Fig 4C is a real panel,
it is simply the recovery ceiling now. Nine rows were wrong that way, and no
mechanical check sees them, because the mapping from a row's content to the right
panel is semantic. A green run here does not mean the callouts are right.

WHY THE VALID SET IS DERIVED AND NOT LISTED
-------------------------------------------
The valid names are parsed from the manuscript's caption lines at test time. A
hardcoded list of names would freeze on the day it was written and then outlive
every renumbering made around it, which is precisely how the errors this file
guards against survived: TABLE_1_LOCK_MD5 pinned Table 1's bytes without
validating its content. Deriving the set means the test tracks the paper.

Coupling tests/ to the manuscript follows the suite's own precedent:
test_manifest_sync.py couples to environment.yml and pyproject.toml, and
test_pan_census_manifest.py couples to a script's UUID lists.
"""
from __future__ import annotations

import re
from pathlib import Path

import openpyxl
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
TABLE_1 = REPO_ROOT / "docs/supplementary_materials/table_S13_test_inventory.xlsx"
MANUSCRIPT = REPO_ROOT / "docs/submission/plosone/manuscript_combined.txt"

SHEET = "S13 Table"   # D67: tab renamed from "Table 1" to match the caption
FIG_COL = 11  # "Figure/Table"
ID_COL = 1

# A caption line opens with the display item's name and a full stop, e.g.
# "Fig 1." or "S4 Table." or "S1 Text.". Anchored, so a mid-sentence mention
# cannot masquerade as a caption.
CAPTION_RE = re.compile(r"^(Fig \d+|S\d+ (?:Fig|Table|Text))\.")

# One difference between how the sheet writes a reference and how the manuscript
# captions it, and it is not an error:
#
#   panel letters   the sheet cites a panel, the manuscript captions the whole
#                   figure -- "Fig 2B" and "S2 Fig C" against "Fig 2" and "S2 Fig"
#
# The sheet now writes supplementary items number first, as the manuscript does.
# The inverted form ("Fig S4", "Table S3") is still recognised, because a workbook
# is hand-maintained and the older form may reappear in an edit; accepting both
# costs nothing and rejecting a correct reference on word order alone would be a
# false failure.
MAIN_FIG_RE = re.compile(r"^(Fig \d+)[A-Z]?$")               # Fig 2, Fig 2B
SUPP_RE = re.compile(r"^S(\d+) (Fig|Table|Text)(?: [A-Z])?$")  # S2 Fig, S2 Fig C, S4 Table
LEGACY_SUPP_RE = re.compile(r"^(Fig|Table|Text) S(\d+)[A-Z]?$")  # Fig S4, Fig S2C

EM_DASH = "—"


def _caption_names() -> set[str]:
    """Every display item the manuscript captions, in the manuscript's own form."""
    names = set()
    for line in MANUSCRIPT.read_text(encoding="utf-8").splitlines():
        m = CAPTION_RE.match(line)
        if m:
            names.add(m.group(1))
    return names


def _sheet_rows() -> list[tuple[str, str]]:
    """(id, figure_cell) for every numbered test row."""
    ws = openpyxl.load_workbook(TABLE_1, data_only=True)[SHEET]
    rows = []
    for row in ws.iter_rows(min_row=2, values_only=True):
        tid = row[ID_COL - 1]
        if tid and re.fullmatch(r"T\d+", str(tid).strip()):
            cell = row[FIG_COL - 1]
            rows.append((str(tid).strip(), "" if cell is None else str(cell).strip()))
    return rows


def _normalise(cell: str) -> str:
    """Map a sheet cell onto the manuscript's caption form, or return it unchanged."""
    main = MAIN_FIG_RE.match(cell)
    if main:                    # "Fig 2" -> "Fig 2";  "Fig 2B" -> "Fig 2"
        return main.group(1)
    supp = SUPP_RE.match(cell)
    if supp:                    # "S2 Fig" -> "S2 Fig";  "S2 Fig C" -> "S2 Fig"
        return "S%s %s" % (supp.group(1), supp.group(2))
    legacy = LEGACY_SUPP_RE.match(cell)
    if legacy:                  # "Fig S4" -> "S4 Fig";  "Fig S2C" -> "S2 Fig"
        return "S%s %s" % (legacy.group(2), legacy.group(1))
    return cell


def test_manuscript_captions_parse() -> None:
    """
    Guard the guard. If the caption regex ever stops matching, every downstream
    assertion would compare against an empty valid set and pass vacuously.
    """
    names = _caption_names()
    assert names, "no caption lines parsed from the manuscript"
    figs = {n for n in names if n.startswith("Fig ")}
    sis = names - figs
    assert len(figs) == 5, f"expected 5 main-figure captions, parsed {len(figs)}: {sorted(figs)}"
    assert len(sis) == 19, f"expected 19 supporting-information captions, parsed {len(sis)}"


def test_table_1_has_rows() -> None:
    """Guard the guard: an empty sheet would pass every per-row check."""
    rows = _sheet_rows()
    assert len(rows) == 64, f"expected 64 numbered test rows, found {len(rows)}"


@pytest.mark.parametrize("tid,cell", _sheet_rows(), ids=lambda v: v if isinstance(v, str) else str(v))
def test_figure_cell_is_not_blank(tid: str, cell: str) -> None:
    """
    Every row states a display item or states that it has none. A blank is neither,
    and is invisible to a check that only validates non-empty cells.
    """
    assert cell, (
        f"{tid} has an empty Figure/Table cell. Use the em-dash to say a row has no "
        f"display item; a blank cannot be distinguished from an oversight."
    )


@pytest.mark.parametrize("tid,cell", _sheet_rows(), ids=lambda v: v if isinstance(v, str) else str(v))
def test_figure_cell_names_a_real_display_item(tid: str, cell: str) -> None:
    """
    Each cell is the em-dash, or names a display item the manuscript captions.
    Panel letters resolve to their parent figure; the sheet's "Fig S4" resolves to
    the manuscript's "S4 Fig".
    """
    if cell == EM_DASH:
        return
    resolved = _normalise(cell)
    valid = _caption_names()
    assert resolved in valid, (
        f"{tid} names {cell!r} (resolved to {resolved!r}), which the manuscript does "
        f"not caption. Captioned items are: {sorted(valid)}"
    )
