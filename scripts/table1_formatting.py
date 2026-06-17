#!/usr/bin/env python3
"""
Stage-3a Table 1 integrity + de-rigidify (reproducible synthesis pass).

Table_1.xlsx (docs/submission/figures_for_review/Table_1.xlsx) is the tracked
artifact; it has no from-scratch repro-pipeline producer (the historical builders
live in the gitignored scripts/archive/). This tracked pass applies the Stage-3a
edits in place, idempotently, so the result is reproducible from tracked code.

Changes (k = 52 unchanged; no rows added/removed):
  B  de-rigidify: "rigidity" -> "Procrustes residual" in the 7 descriptions that
     contain it (T33, T36, T37, T40, T53, T60, T61).
  C1 T59 value -0.393 -> -0.410, raw p 0.087 -> 0.073 (Stage-0 V6 recompute).
  C2 T59 Figure/Table "Fig S3D" -> "Fig S3B" (Figure S3 has only panels A, B).
  C3 T53 value +0.247 -> -0.247 (Stage-2b SAMap recompute vs Procrustes residual).
  C5 T44 Test Type "Spearman + FDR" -> "Fisher's exact + FDR" (the test is Fisher's).
  C6 footnote: state the exclusion principle for the k = 52 family.
  D22 reconcile T34/T35 to the manuscript Table 1: description -> overall 500-gene
     identity-set CellMarker enrichment (primary / expression-matched background),
     n -> "500 genes" (was the pre-d19 "top 50 genes x 6 types" / "6 types" labels),
     and add the T34/T35 caption note. Values (4.49, 3.32) and p-values unchanged.
  D28 Minor 1: T11 corrected-p (Bonferroni, col I) cell "significant" ->
     "direction-robust (100/100; resampling)", matching the body / Methods /
     Fig 3E framing of the +0.159 donor-split delta as direction-robustness, and
     the manuscript Table 1 source row. (Minor 3 / tiny-p sci-notation display was
     already correct in-repo: T34/T35/T38 raw + Bonferroni cells are 0.00E+00.)
Bonferroni columns for T59/T53 stay at 1.0 (p*52 > 1, capped) — no change.
"""
from pathlib import Path
import re
import openpyxl

XLSX = Path(__file__).resolve().parent.parent / "docs/submission/figures_for_review/Table_1.xlsx"
COL = {"id": 1, "desc": 3, "testtype": 4, "n": 5, "value": 7, "rawp": 8, "fig": 11}


def main():
    wb = openpyxl.load_workbook(XLSX)
    ws = wb["Table 1"]
    rowmap = {}
    footnote_row = None
    for r in range(1, ws.max_row + 1):
        cid = ws.cell(r, COL["id"]).value
        if cid is None:
            continue
        cid = str(cid).strip()
        if re.match(r"T\d+$", cid):
            rowmap[cid] = r
        elif cid.startswith("Bonferroni correction applied"):
            footnote_row = r

    # ---- B: de-rigidify descriptions ("rigidity" -> "Procrustes residual") ----
    derig = 0
    for cid, r in rowmap.items():
        d = ws.cell(r, COL["desc"]).value
        if d and "rigidity" in d.lower():
            ws.cell(r, COL["desc"]).value = re.sub(r"rigidity", "Procrustes residual", d, flags=re.I)
            derig += 1
    assert derig in (0, 7), f"expected 0 (idempotent re-run) or 7 rigidity descriptions, found {derig}"

    # ---- C1: T59 value + raw p ----
    ws.cell(rowmap["T59"], COL["value"]).value = -0.410
    ws.cell(rowmap["T59"], COL["rawp"]).value = 0.073
    # ---- C2: T59 panel ref ----
    ws.cell(rowmap["T59"], COL["fig"]).value = "Fig S3B"
    # ---- C3: T53 SAMap sign ----
    ws.cell(rowmap["T53"], COL["value"]).value = -0.247
    # ---- C5: T44 test method ----
    tt = ws.cell(rowmap["T44"], COL["testtype"]).value
    if tt and "spearman" in str(tt).lower():
        ws.cell(rowmap["T44"], COL["testtype"]).value = "Fisher's exact + FDR"

    # ---- C6: exclusion-principle footnote wording ----
    if footnote_row:
        f = ws.cell(footnote_row, COL["id"]).value
        f = f.replace(
            "Bonferroni correction applied over k = 52 inferential tests "
            "(excluding T03, T04, T09, T10, T18, T62, T63, T64, T65 — descriptive/exploratory, aggregate, or diagnostic).",
            "Bonferroni correction applied over the k = 52 principal inferential tests "
            "(k = 52 = 61 − 9). The 9 excluded IDs (T03, T04, T09, T10, T18, T62, T63, T64, T65) are descriptive, "
            "aggregate, or diagnostic quantities that make no inferential significance claim and are therefore not "
            "part of the multiple-testing family.")
        ws.cell(footnote_row, COL["id"]).value = f

    # ---- D22: reconcile T34/T35 to the manuscript Table 1 (idempotent) ----
    # The manuscript (d19) describes the inventoried CellMarker enrichments as the
    # overall 500-gene identity set (n = 500 genes); the spreadsheet retained the
    # pre-d19 "top 50 genes x 6 types" / "6 types" labels. Reconcile description + n
    # and add the T34/T35 caption note. Values (4.49, 3.32) + p-values unchanged.
    ws.cell(rowmap["T34"], COL["desc"]).value = (
        "Overall 500-gene identity-set CellMarker enrichment (primary 16,959-gene background)")
    ws.cell(rowmap["T35"], COL["desc"]).value = (
        "Overall 500-gene identity-set CellMarker enrichment (expression-matched background)")
    for _tid in ("T34", "T35"):
        ws.cell(rowmap[_tid], COL["n"]).value = "500 genes"
    if footnote_row:
        _fv = ws.cell(footnote_row, COL["id"]).value
        if "T34/T35:" not in _fv:
            _t3435 = (
                "T34/T35: the inventoried fold-enrichments are for the overall 500-gene "
                "identity set; the per-type top-50 centroid-deviation test on the 6 validation "
                "types is reported as a pass-rate (5 of 6 under CellMarker, 6 of 6 under the "
                "held-out HPA reference; Figure S6, Methods), not as a single inventory p-value. ")
            _anchor = "T29, T54, and the T55-T57 treeness tests"
            _fv = _fv.replace(_anchor, _t3435 + _anchor, 1)
            ws.cell(footnote_row, COL["id"]).value = _fv

    # ---- D28: Minor 1 — T11 corrected-p (Bonferroni, col I) label (idempotent) ----
    # The donor-split delta (+0.159) is framed in the body / Methods / Fig 3E caption
    # / Limitations as direction-robustness (100/100 resampling splits), not a
    # calibrated effect. Align the table's corrected-p cell with that framing and with
    # the manuscript Table 1 source row. Column I (9) = "Bonferroni p (k=52)".
    ws.cell(rowmap["T11"], 9).value = "direction-robust (100/100; resampling)"

    # ---- E.2: figure-column renumber (D3 main-figure map; idempotent) ----
    FIG_SET = {"T12": "Fig 5A", "T14": "Fig 5B", "T15": "Fig 5C", "T30": "Fig 4C",
               "T31": "Fig 4D", "T36": "Fig 6B", "T37": "Fig 6B", "T38": "Fig 7A",
               "T39": "Fig 7A", "T42": "Fig 7B", "T43": "Fig 7B", "T44": "Fig 7B",
               "T45": "Fig 7B", "T46": "Fig 7B", "T47": "Fig 7B", "T48": "Fig 7B",
               "T49": "Fig 7B", "T50": "Fig 7B", "T51": "Fig 7B"}
    for tid, val in FIG_SET.items():
        ws.cell(rowmap[tid], COL["fig"]).value = val

    # ---- Stage-5.5 C: presentation formatting (reproducible cell sizing) ----
    # Column widths so the Description / Status / Category cells are not clipped.
    from openpyxl.styles import Alignment
    widths = {"A": 6, "B": 18, "C": 72, "D": 24, "E": 7, "F": 11,
              "G": 11, "H": 13, "I": 14, "J": 27, "K": 14}
    for col, w in widths.items():
        ws.column_dimensions[col].width = w
    # Wrap text + top-align the prose columns; size each data row to its lines.
    wrap_cols = {3: 72, 10: 27, 2: 18, 4: 24}   # Description, Status, Category, Test Type
    for r in rowmap.values():
        n_lines = 1
        for ci, w in wrap_cols.items():
            ws.cell(r, ci).alignment = Alignment(wrap_text=True, vertical="top")
            txt = str(ws.cell(r, ci).value or "")
            n_lines = max(n_lines, -(-len(txt) // max(1, w - 2)))
        ws.row_dimensions[r].height = 15 * n_lines + 4
    # Footnote: span the full table width, wrap, and a tall row so it shows fully.
    if footnote_row:
        ws.merge_cells(start_row=footnote_row, start_column=1,
                       end_row=footnote_row, end_column=11)
        fc = ws.cell(footnote_row, 1)
        fc.alignment = Alignment(wrap_text=True, vertical="top", horizontal="left")
        ws.row_dimensions[footnote_row].height = 150

    # ---- E.2 rev: deterministic save (eliminate wall-clock timestamps) ----
    from datetime import datetime as _dt
    import zipfile as _zf
    import os as _os
    _FIXED = _dt(2026, 1, 1, 0, 0, 0)
    _ISO = "2026-01-01T00:00:00Z"
    wb.properties.created = _FIXED
    wb.properties.modified = _FIXED

    wb.save(XLSX)

    # Rewrite the .xlsx zip: fixed entry date_time for every member + fixed
    # docProps/core.xml dcterms timestamps (same content/names/order/compression).
    with _zf.ZipFile(XLSX, "r") as _zin:
        _members = [(i, _zin.read(i.filename)) for i in _zin.infolist()]
    _tmp = str(XLSX) + ".tmp"
    with _zf.ZipFile(_tmp, "w") as _zout:
        for _info, _data in _members:
            if _info.filename == "docProps/core.xml":
                _s = _data.decode("utf-8")
                _s = re.sub(r"(<dcterms:created[^>]*>)[^<]*(</dcterms:created>)",
                            r"\g<1>" + _ISO + r"\g<2>", _s)
                _s = re.sub(r"(<dcterms:modified[^>]*>)[^<]*(</dcterms:modified>)",
                            r"\g<1>" + _ISO + r"\g<2>", _s)
                _data = _s.encode("utf-8")
            _ni = _zf.ZipInfo(_info.filename, date_time=(2026, 1, 1, 0, 0, 0))
            _ni.compress_type = _info.compress_type
            _ni.external_attr = _info.external_attr
            _ni.create_system = _info.create_system
            _zout.writestr(_ni, _data)
    _os.replace(_tmp, XLSX)

    # ---- verify ----
    wb2 = openpyxl.load_workbook(XLSX)
    w = wb2["Table 1"]
    rig = sum(1 for r in range(1, w.max_row + 1)
              if w.cell(r, COL["desc"]).value and "rigidity" in str(w.cell(r, COL["desc"]).value).lower())
    print(f"  de-rigidified {derig} descriptions; remaining 'rigidity' in descriptions: {rig}")
    print(f"  T59 value={w.cell(rowmap['T59'],COL['value']).value} rawp={w.cell(rowmap['T59'],COL['rawp']).value} fig={w.cell(rowmap['T59'],COL['fig']).value}")
    print(f"  T53 value={w.cell(rowmap['T53'],COL['value']).value}")
    print(f"  T44 testtype={w.cell(rowmap['T44'],COL['testtype']).value}")
    for _tid in ("T34", "T35"):
        print(f"  {_tid} desc={w.cell(rowmap[_tid],COL['desc']).value!r} n={w.cell(rowmap[_tid],COL['n']).value!r} value={w.cell(rowmap[_tid],COL['value']).value}")
    print(f"  T11 corrected-p (I{rowmap['T11']})={w.cell(rowmap['T11'],9).value!r}")
    print(f"  T34/T35 note present: {'T34/T35:' in str(w.cell([r for r in range(1,w.max_row+1) if str(w.cell(r,1).value or '').startswith('Bonferroni correction applied')][0],1).value)}")
    print(f"  saved {XLSX}")


if __name__ == "__main__":
    main()
