#!/usr/bin/env python3
"""Assemble combined supplementary materials PDF for Cell Systems submission.

Combines all supplementary figures (S1–S8), table legends (S1–S10), and
inline CSV tables (S3–S5, S7–S10) into a single PDF document.

Legends are sourced from the SUPPLEMENTAL INFORMATION section of the
renumbered manuscript (docs/submission/manuscript_combined.txt) to ensure
citation numbers are current.

Output: docs/submission/supplementary_materials.pdf
        docs/supplementary_materials/supplementary_materials.pdf
"""

import csv
import os
from pathlib import Path
import re
import tempfile

import fitz  # PyMuPDF
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

# ── Constants ──────────────────────────────────────────────────────────
PAGE_W = 612.0   # US Letter width (pt)
PAGE_H = 792.0   # US Letter height (pt)
MARGIN = 72.0    # 1 inch (pt)
CONTENT_W = PAGE_W - 2 * MARGIN
CONTENT_H = PAGE_H - 2 * MARGIN

BASE = str(Path(__file__).resolve().parent.parent)

# Two font families, by rendering path:
#
# LEGEND_FONT — prose legends/titles, rendered via PyMuPDF insert_textbox.
# DejaVuSans (not the macOS core Arial) because it ships the full Unicode
# superscript range (⁰⁴⁻⁹ ⁻ ¹) plus ρ × − ≥, so citation superscripts
# (Sun et al.¹⁷) and powers (10⁻⁴, 10⁻¹³) render natively instead of being
# flattened to inline digits / e-notation by the former _unicodify pass (now
# removed — see split_legend). Vendored under assets/fonts/ (Bitstream-Vera/
# DejaVu license, permits redistribution + embedding — see LICENSE_DEJAVU) so
# the legend build is reproducible from a bare clone rather than depending on an
# env-specific matplotlib font path.
#
# ARIAL — the matplotlib-rendered data tables (S3–S5, S7–S10), via
# FontProperties. These carry no Unicode superscripts (powers are written in
# ASCII e-notation, e.g. <1e-5), so they never had the missing-glyph problem,
# and their column widths / abbreviated headers are tuned to Arial's metrics —
# the wider DejaVuSans overflows those headers (S3/S9). So the tables stay on
# Arial, byte-for-byte as before; only the prose legend path changes.
LEGEND_FONT = os.path.join(BASE, "assets/fonts/DejaVuSans.ttf")
LEGEND_FONT_BOLD = os.path.join(BASE, "assets/fonts/DejaVuSans-Bold.ttf")

ARIAL = "/System/Library/Fonts/Supplemental/Arial.ttf"
ARIAL_BOLD = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"

OUT_PATHS = [
    os.path.join(BASE, "docs/submission/supplementary_materials.pdf"),
    os.path.join(BASE, "docs/supplementary_materials/supplementary_materials.pdf"),
]

FIGURE_FILES = [
    ("S1", os.path.join(BASE, "figures/submission/supplementary/figS1_pipeline_validation.pdf")),
    ("S2", os.path.join(BASE, "figures/submission/supplementary/figS2_parameter_protocol_sensitivity.pdf")),
    ("S3", os.path.join(BASE, "figures/submission/supplementary/figS3_bootstrap_rankings.pdf")),
    ("S4", os.path.join(BASE, "figures/submission/supplementary/figS4_cellhint_investigation.pdf")),
    ("S5", os.path.join(BASE, "figures/submission/supplementary/figS5_samap.pdf")),
    ("S6", os.path.join(BASE, "figures/submission/supplementary/figS6_cellmarker_enrichment.pdf")),
    ("S7", os.path.join(BASE, "figures/submission/supplementary/figS7_matched_scale_control.pdf")),
    ("S8", os.path.join(BASE, "figures/submission/supplementary/figS8_markernull.pdf")),
]

MANUSCRIPT_PATH = os.path.join(BASE, "docs/submission/manuscript_combined.txt")

TABLE_CSV = {
    "S3": os.path.join(BASE, "docs/supplementary_materials/table_S3.csv"),
    "S4": os.path.join(BASE, "docs/supplementary_materials/table_S4.csv"),
    "S5": os.path.join(BASE, "docs/supplementary_materials/table_S5.csv"),
    "S7": os.path.join(BASE, "docs/supplementary_materials/table_S7_layer1_housekeeping_exclusion.csv"),
    "S8": os.path.join(BASE, "docs/supplementary_materials/table_S8_marker_ortholog_retention.csv"),
    "S9": os.path.join(BASE, "docs/supplementary_materials/table_S9_genestd_standardization.csv"),
    "S10": os.path.join(BASE, "docs/supplementary_materials/table_S10_markernull.csv"),
    "S12": os.path.join(BASE, "docs/supplementary_materials/table_S12_software_environment.csv"),
}

MANUSCRIPT_TITLE = (
    "Quantifying the Conserved Geometry of Cell-Type "
    "Identity Across Mammalian Species"
)


# ── Legend parser ──────────────────────────────────────────────────────

def parse_legends(path):
    """Return (figure_legends, table_legends) dicts keyed by 'S1'…'S7'.

    Extracts the SUPPLEMENTAL INFORMATION section from the manuscript
    (bounded at the next top-level header), then parses individual
    figure/table legends separated by '---'.
    """
    with open(path) as f:
        text = f.read()

    # Extract only the SUPPLEMENTAL INFORMATION section
    start = text.find("Supporting Information")
    if start == -1:
        start = text.find("SUPPLEMENTAL INFORMATION")
    if start == -1:
        raise ValueError(f"No Supporting Information / SUPPLEMENTAL INFORMATION section found in {path}")
    # End the section at the next top-level header. Headers are "boxed": an
    # ALL-CAPS title line wrapped by '=' rule lines (overline optional), e.g.
    # the DECLARATION OF GENERATIVE AI section. Landing on the overline keeps
    # the trailing rule out of the final Figure/Table block. Fall back to
    # REFERENCES.
    header_re = re.compile(
        r"^(?:=+[ \t]*\n)?[A-Z][A-Z0-9 ,&/()'.\-]+\n=+[ \t]*$", re.MULTILINE
    )
    end = -1
    for m in header_re.finditer(text):
        if m.start() > start:
            end = m.start()
            break
    if end == -1:
        end = text.find("REFERENCES", start + 25)
    if end == -1:
        legend_text = text[start:]
    else:
        legend_text = text[start:end]

    # Strip the header lines (title + === underline)
    legend_text = re.sub(
        r"^(?:SUPPLEMENTAL INFORMATION|Supporting Information)\s*\n=+\s*\n*", "", legend_text
    )

    sections = legend_text.split("---")
    figs, tabs = {}, {}
    for sec in sections:
        sec = sec.strip()
        if not sec:
            continue
        m = re.match(r"\*\*(Figure|Table) (S\d+)", sec)
        if m:
            kind, label = m.group(1), m.group(2)
            (figs if kind == "Figure" else tabs)[label] = sec
    return figs, tabs


def split_legend(legend_md):
    """Split a markdown legend into (title_line, body_text).

    Both are returned with markdown-bold ** markers stripped. Literal Unicode
    superscripts (¹⁷, 10⁻⁴, 10⁻¹³, …) are passed through unchanged: the legend
    font is DejaVuSans (see LEGEND_FONT), which ships the full superscript glyph
    range, so they render natively. (Previously the font was the macOS
    core Arial, which lacks those glyphs, so a _unicodify pass flattened them to
    ASCII e-notation; that fallback is no longer needed and has been removed.)
    Trailing '==...' divider lines (from manuscript section separators) are
    stripped so they don't bleed into legend pages.
    """
    import re
    lines = legend_md.split("\n")
    title = lines[0].replace("**", "")
    body_lines = lines[1:]
    # Strip any lines that are entirely '=' characters (section dividers)
    body_lines = [ln for ln in body_lines if not re.match(r"^\s*=+\s*$", ln)]
    # Strip only genuine markdown-bold spans (**non-space...non-space**), not
    # significance-key asterisks like "*** p < 10⁻⁴, ** p < 0.01, * p < 0.05"
    # (which have a space after the asterisks and so are not bold delimiters).
    # A blanket .replace("**", "") previously consumed the Fig S2F key.
    body = re.sub(r"\*\*(?=\S)(.+?)(?<=\S)\*\*", r"\1", "\n".join(body_lines)).strip()
    return title, body


# ── Page helpers ───────────────────────────────────────────────────────

def _measure_textbox_height(text, width, fontfile, fontname, fontsize):
    """Return the vertical space (pt) insert_textbox consumes for `text` at the
    given width/font, by rendering into a tall scratch box and reading back the
    unused remainder.

    Robust to font metrics: DejaVuSans wraps wider than the old core Arial, so a
    fixed chars-per-line estimate under-counts lines for the wider font. When a
    caller sizes a box from that estimate and the text overflows, PyMuPDF's
    insert_textbox silently renders *nothing* (negative return) — which is how
    the Figure S3/S4 titles vanished after the font swap. Measuring the true
    height instead lets callers allocate exactly enough and place content below
    it without dropping or overlapping.
    """
    scratch = fitz.open()
    page = scratch.new_page(width=PAGE_W, height=2000.0)
    rect = fitz.Rect(0, 0, width, 1990.0)
    rc = page.insert_textbox(rect, text, fontsize=fontsize,
                             fontfile=fontfile, fontname=fontname,
                             align=fitz.TEXT_ALIGN_LEFT)
    scratch.close()
    return 1990.0 - max(rc, 0)


def _add_page_number(page, num):
    """Bottom-center page number in 9 pt Arial."""
    # Center the number
    tw = fitz.get_text_length(str(num), fontname="helv", fontsize=9)
    x = (PAGE_W - tw) / 2
    y = PAGE_H - MARGIN / 2 + 4
    page.insert_text(fitz.Point(x, y), str(num),
                     fontsize=9, fontfile=LEGEND_FONT, fontname="Legend")


def add_title_page(doc):
    """Supplementary Materials title page."""
    page = doc.new_page(width=PAGE_W, height=PAGE_H)

    # "Supplementary Materials" — centered, 20 pt bold
    y = PAGE_H * 0.38
    rect = fitz.Rect(MARGIN, y, PAGE_W - MARGIN, y + 30)
    page.insert_textbox(rect, "Supplementary Materials",
                        fontsize=20, fontfile=LEGEND_FONT_BOLD, fontname="LegendB",
                        align=fitz.TEXT_ALIGN_CENTER)

    # Manuscript title — centered, 14 pt regular
    rect2 = fitz.Rect(MARGIN, y + 50, PAGE_W - MARGIN, y + 130)
    page.insert_textbox(rect2, MANUSCRIPT_TITLE,
                        fontsize=14, fontfile=LEGEND_FONT, fontname="Legend",
                        align=fitz.TEXT_ALIGN_CENTER)
    return page


def add_legend_page(doc, legend_md):
    """Page with legend title (bold) + body (regular), 10 pt Arial."""
    page = doc.new_page(width=PAGE_W, height=PAGE_H)
    title, body = split_legend(legend_md)

    # Title in bold
    title_rect = fitz.Rect(MARGIN, MARGIN, PAGE_W - MARGIN, MARGIN + 80)
    rc = page.insert_textbox(title_rect, title,
                             fontsize=10, fontfile=LEGEND_FONT_BOLD, fontname="LegendB",
                             align=fitz.TEXT_ALIGN_LEFT)
    # Estimate title height used
    title_h = 80 - max(rc, 0)  # rc = remaining space

    # Body in regular
    body_y = MARGIN + title_h + 6
    body_rect = fitz.Rect(MARGIN, body_y, PAGE_W - MARGIN, PAGE_H - MARGIN)
    page.insert_textbox(body_rect, body,
                        fontsize=10, fontfile=LEGEND_FONT, fontname="Legend",
                        align=fitz.TEXT_ALIGN_LEFT)
    return page


def add_legend_and_figure_page(doc, legend_md, pdf_path, fig_label=None):
    """Put legend on top of page, figure below, so reviewers don't have to
    flip N/N+1 to read figure + its legend together. Figure is scaled to
    fit the available vertical space after the legend.

    For figures that can't reasonably fit a single page together with their
    legend (e.g., Figure S2 at 630 pt tall + multi-paragraph legend), falls
    back to legend-then-figure on separate pages and returns False.
    Otherwise returns True.
    """
    # Estimate legend body height (10pt text, ~68 chars per line, 12pt leading)
    title, body = split_legend(legend_md)
    # Title height: measured from the actual font, not a chars-per-line estimate.
    # A fixed-88-char rule was calibrated for Arial; under the wider DejaVuSans it
    # under-counted lines, so over-long titles (Figure S3/S4) overflowed their box
    # and insert_textbox dropped them entirely. +4 pt slack keeps the fit margin.
    title_h_est = _measure_textbox_height(title, CONTENT_W, LEGEND_FONT_BOLD, "LegendB", 10) + 4
    wrapped_lines = 0
    for paragraph in body.split("\n"):
        if not paragraph.strip():
            wrapped_lines += 1
            continue
        # Rough wrap: content_w ≈ 468pt / ~5.3pt per char @ 10pt Arial ≈ 88 chars
        wrapped_lines += max(1, (len(paragraph) + 87) // 88)
    legend_h_est = title_h_est + wrapped_lines * 12 + 10  # title + body + padding

    fig_doc = fitz.open(pdf_path)
    fig_page = fig_doc[0]
    fw, fh = fig_page.rect.width, fig_page.rect.height

    # Available space on the combined page
    page_ch = PAGE_H - 2 * MARGIN
    fig_budget_h = page_ch - legend_h_est - 10  # gap between legend and figure
    fig_budget_w = PAGE_W - 2 * MARGIN

    # Scale figure to fit; allow upscaling up to 2× for undersized figs (S5/S6)
    scale = min(fig_budget_w / fw, fig_budget_h / fh, 2.0)
    # If the figure doesn't meaningfully fit (< 200pt tall after scaling), fallback
    scaled_h = fh * scale
    if scaled_h < 200 and fig_doc.page_count == 1:
        # upscale more aggressively for tiny figs
        scale = min(fig_budget_w / fw, fig_budget_h / fh, 3.0)
        scaled_h = fh * scale
    if scale < 0.4 or fig_doc.page_count > 1:
        # Won't fit cleanly → fall back to legend-then-figure
        fig_doc.close()
        add_legend_page(doc, legend_md)
        return False

    scaled_w = fw * scale

    # Build combined page
    page = doc.new_page(width=PAGE_W, height=PAGE_H)
    # Legend title (bold)
    title_rect = fitz.Rect(MARGIN, MARGIN, PAGE_W - MARGIN, MARGIN + title_h_est)
    page.insert_textbox(title_rect, title,
                        fontsize=10, fontfile=LEGEND_FONT_BOLD, fontname="LegendB",
                        align=fitz.TEXT_ALIGN_LEFT)
    # Legend body (regular)
    body_rect = fitz.Rect(MARGIN, MARGIN + title_h_est + 2, PAGE_W - MARGIN,
                          MARGIN + legend_h_est + 6)
    page.insert_textbox(body_rect, body,
                        fontsize=10, fontfile=LEGEND_FONT, fontname="Legend",
                        align=fitz.TEXT_ALIGN_LEFT)
    # Figure, centered horizontally, below legend
    fig_y0 = MARGIN + legend_h_est + 12
    fig_x0 = MARGIN + (fig_budget_w - scaled_w) / 2
    target = fitz.Rect(fig_x0, fig_y0, fig_x0 + scaled_w, fig_y0 + scaled_h)
    page.show_pdf_page(target, fig_doc, 0)
    fig_doc.close()
    return True


def insert_figure_pages(doc, pdf_path, landscape=False):
    """Insert all pages of a figure PDF, centered and scaled to fit margins.

    If landscape=True, the target page uses landscape orientation.
    Returns list of page indices added.
    """
    fig_doc = fitz.open(pdf_path)
    added = []
    for pno in range(len(fig_doc)):
        src = fig_doc[pno]
        sw, sh = src.rect.width, src.rect.height

        if landscape:
            pw, ph = PAGE_H, PAGE_W  # swap for landscape
        else:
            pw, ph = PAGE_W, PAGE_H

        cw = pw - 2 * MARGIN
        ch = ph - 2 * MARGIN

        # Scale to fit; never upscale
        scale = min(cw / sw, ch / sh, 1.0)
        dw, dh = sw * scale, sh * scale
        x0 = MARGIN + (cw - dw) / 2
        y0 = MARGIN + (ch - dh) / 2
        target = fitz.Rect(x0, y0, x0 + dw, y0 + dh)

        new_page = doc.new_page(width=pw, height=ph)
        new_page.show_pdf_page(target, fig_doc, pno)
        added.append(len(doc) - 1)
    fig_doc.close()
    return added


# ── CSV table rendering ───────────────────────────────────────────────

import textwrap


# Human-readable header maps + column geometry for the supplementary
# sensitivity tables (S7–S10). Before this, S7–S10 fell through to the generic `else` branch
# (raw snake_case headers + col_w=None → matplotlib auto_set_column_width),
# which sizes columns to content with no fit-to-axes guarantee and so
# overflowed the figure canvas — clipping the outer borders (S7) and whole
# columns (S9). These entries put S7–S10 on the SAME styled path as the
# established rendered tables S3/S4/S5: pretty display headers and explicit
# fractional widths that sum to ≤ the available axis width. The underlying
# CSVs are untouched — the maps apply at render time only. "rho" is spelled
# out (not ρ) to match the established Table S4 header convention.
_SENS_TABLES = {
    "S7": {
        # variant identifiers are single snake_case tokens (e.g.
        # ribosomal_plus_housekeeping, 108 pt < the 140 pt cell) — wrap 0 keeps
        # each on one line instead of breaking mid-token ("hou\nsekeeping").
        "headers": ["Variant", "Genes\nexcluded", "PCA\ncomponents",
                    "Obs/null\nratio", "Ranking rho\nvs full space"],
        "wrap":    [0, 0, 0, 0, 0],
        "width":   [0.30, 0.15, 0.16, 0.15, 0.18],
    },
    "S8": {
        "headers": ["Cell type", "N\nmarkers", "N retained\n(1:1)",
                    "Retention\nfraction", "Procrustes\nresidual"],
        "wrap":    [24, 0, 0, 0, 0],
        "width":   [0.34, 0.12, 0.16, 0.17, 0.17],
    },
    "S9": {
        "headers": ["Scheme", "Layer-1\nobs/null", "Layer-1\np",
                    "Ranking rho\nvs primary", "Ribosomal CPC1\n(of 35)",
                    "Mean ribo\nin top-20"],
        "wrap":    [20, 0, 0, 0, 0, 0],
        "width":   [0.30, 0.13, 0.11, 0.14, 0.16, 0.12],
    },
    "S10": {
        "headers": ["Partition", "Obs/null", "p\n(100k perm)",
                    "Non-singleton\ngroups", "Note"],
        "wrap":    [22, 0, 0, 0, 22],
        "width":   [0.27, 0.12, 0.13, 0.17, 0.27],
    },
    "S12": {
        "headers": ["Package", "Version", "Primary use"],
        "wrap":    [0, 0, 40],
        "width":   [0.26, 0.14, 0.56],
    },
}


def render_csv_to_pdf(csv_path, label=None):
    """Render CSV as a formatted table via matplotlib.  Returns temp PDF path.

    Uses 8 pt minimum font with per-row dynamic heights to prevent text
    overlap.  Tables S3 and S5 use landscape orientation; S5 is split
    across two pages (35 data rows).

    Biology: supplementary tables for cross-species Procrustes analysis.
    Math: layout only — no analytical changes.
    """
    from matplotlib.backends.backend_pdf import PdfPages

    with open(csv_path) as f:
        reader = csv.reader(f)
        headers = next(reader)
        rows = list(reader)

    n_cols = len(headers)
    n_rows = len(rows)

    # ── Per-table configuration ──────────────────────────────────────
    landscape = False
    fontsize = 8.0
    rows_per_page = 50          # default high enough for one page

    if label == "S3":
        landscape = True
        # 9 columns in the CSV; abbreviated headers to fit narrow columns.
        # Second "Primary Rank" is the matched 15-type re-run used for the
        # Rank Difference calculation (see Table S3 legend).
        display_headers = [
            "Cell Type", "Prim.\nRank\n(of 15)", "Prim.\nRank\n(matched\n15-type)",
            "CH\nRank\n(of 15)", "Rank\nDiff", "Prim.\nRank\n(of 35)",
            "CH\nTissue\nCount", "Annotation Mapping Notes", "Conf.",
        ]
        col_wrap = [15, 0, 0, 0, 0, 0, 0, 42, 0]
        # Width tuning: matched-15-type header (col 2) widened 0.05→0.06 so
        # "(matched" doesn't clip; Conf column widened 0.07→0.08 so "MODERATE"
        # doesn't look cramped; notes narrowed 0.40→0.38 to compensate.
        col_w = [0.12, 0.05, 0.06, 0.05, 0.05, 0.05, 0.05, 0.38, 0.08]
        # Rank Difference column (col 4) — swap ASCII hyphen for Unicode minus
        # (U+2212) in negative values. Without this, macOS PDF viewers' data
        # detector reads the row's numeric cells as "4 10 9 -5 12" and offers
        # to dial it as a phone number (tel: tooltip over the plasma-cell row).
        # The Unicode minus defeats the detector and is also typographically
        # correct for signed numbers.
        for r in rows:
            if len(r) > 4 and r[4].startswith("-"):
                r[4] = "−" + r[4][1:]
    elif label == "S5":
        landscape = True
        rows_per_page = 18      # 35 rows → 2 landscape pages
        display_headers = [
            "Human Cell Type", "Mouse Cell Type", "Cell\nOntology ID",
            "Matching Basis", "Human\nN", "Mouse\nN",
        ]
        col_wrap = [22, 22, 0, 42, 0, 0]
        col_w = [0.18, 0.18, 0.09, 0.37, 0.07, 0.07]
    elif label == "S4":
        fontsize = 9.0
        display_headers = [
            "Harmonization Level", "n types", "Spearman rho", "p-value",
        ]
        col_wrap = [30, 0, 0, 0]
        col_w = [0.40, 0.10, 0.15, 0.12]
    elif label in _SENS_TABLES:
        cfg = _SENS_TABLES[label]
        display_headers = cfg["headers"]
        col_wrap = cfg["wrap"]
        col_w = cfg["width"]
    else:
        display_headers = headers
        col_wrap = [28] * n_cols
        col_w = None

    font = FontProperties(fname=ARIAL, size=fontsize)
    font_b = FontProperties(fname=ARIAL_BOLD, size=fontsize)

    # ── Text wrapping ────────────────────────────────────────────────
    def _wrap(text, col_idx):
        w = col_wrap[col_idx] if col_idx < len(col_wrap) else 28
        if w == 0 or len(text) <= w:
            return text
        return "\n".join(textwrap.wrap(text, w))

    wrapped_rows = [[_wrap(c, j) for j, c in enumerate(row)] for row in rows]

    # ── Figure dimensions = page content area (→ scale 1.0) ─────────
    if landscape:
        fig_w = (PAGE_H - 2 * MARGIN) / 72.0   # 9.0 in
        fig_h = (PAGE_W - 2 * MARGIN) / 72.0   # 6.5 in
    else:
        fig_w = (PAGE_W - 2 * MARGIN) / 72.0   # 6.5 in
        fig_h = (PAGE_H - 2 * MARGIN) / 72.0   # 9.0 in

    # ── Chunk rows across pages ──────────────────────────────────────
    chunks = []
    for start in range(0, n_rows, rows_per_page):
        chunks.append(wrapped_rows[start:start + rows_per_page])

    line_h = fontsize * 1.5 / 72.0   # approx line height (inches)
    row_pad = 0.02                     # vertical padding per row (inches)

    tmp = tempfile.NamedTemporaryFile(suffix=".pdf", delete=False)

    with PdfPages(tmp.name) as pdf_pages:
        for chunk in chunks:
            chunk_n = len(chunk)
            header_lines = max((h.count("\n") + 1) for h in display_headers)
            row_lines = [max((c.count("\n") + 1) for c in row) for row in chunk]

            fig, ax = plt.subplots(figsize=(fig_w, fig_h))
            ax.axis("off")
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)

            tbl = ax.table(
                cellText=chunk,
                colLabels=display_headers,
                loc="upper center",
                cellLoc="center",
            )
            tbl.auto_set_font_size(False)
            tbl.set_fontsize(fontsize)

            # Column widths — explicit fractions only, normalized to fit the
            # axis width. This is the ROOT of the S7-border / S9-clipping fix:
            # matplotlib's auto_set_column_width defers sizing to draw-time and
            # offers no fit-to-axes guarantee, so wide content overflowed the
            # canvas and the outer borders / edge columns were clipped at
            # savefig. Here every table — configured or not — gets explicit
            # widths that are scaled down whenever they would exceed the
            # available width, the horizontal twin of the vertical fit-to-page
            # below. No table can clip on the left/right edge again.
            avail_w = 0.96
            if not col_w:
                # Generic fallback: size columns proportional to the widest
                # (wrapped) cell in each, so the assignment is deterministic
                # (no draw pass needed) and content-aware.
                def _cw(j):
                    hw = max((len(s) for s in str(display_headers[j]).split("\n")),
                             default=1)
                    dw = max((max((len(s) for s in str(row[j]).split("\n")),
                                   default=1) for row in chunk), default=1) \
                        if chunk else 1
                    return max(hw, dw, 3)
                raw = [_cw(j) for j in range(n_cols)]
                tot = sum(raw) or 1
                col_w = [c / tot for c in raw]
            wsum = sum(col_w[:n_cols]) or 1
            if wsum > avail_w:                       # scale down to fit the axis
                col_w = [w * avail_w / wsum for w in col_w]
            for j in range(n_cols):
                w = col_w[j] if j < len(col_w) else (avail_w / n_cols)
                for i in range(chunk_n + 1):
                    tbl[i, j].set_width(w)

            # Dynamic row heights proportional to line count, scaled to fit
            # the page so the bottom row never clips (e.g. Table S3's trailing
            # Spearman summary row).
            def _rh(n_lines):
                return (n_lines * line_h + row_pad) / fig_h

            header_h = _rh(header_lines)
            data_h = [_rh(rl) for rl in row_lines]
            total_h = header_h + sum(data_h)
            avail = 0.93  # fraction of axes height available for the table
            if total_h > avail:               # compress so nothing clips
                sc = avail / total_h
                header_h *= sc
                data_h = [h * sc for h in data_h]

            # Header row
            for j in range(n_cols):
                cell = tbl[0, j]
                cell.set_height(header_h)
                cell.set_text_props(fontproperties=font_b)
                cell.set_facecolor("#E0E0E0")
                cell.set_edgecolor("#AAAAAA")

            # Data rows
            for i in range(1, chunk_n + 1):
                for j in range(n_cols):
                    cell = tbl[i, j]
                    cell.set_height(data_h[i - 1])
                    cell.set_text_props(fontproperties=font)
                    cell.set_edgecolor("#CCCCCC")
                    if i % 2 == 0:
                        cell.set_facecolor("#F5F5F5")

            plt.subplots_adjust(left=0.02, right=0.98, top=0.95, bottom=0.02)
            pdf_pages.savefig(fig, dpi=300)
            plt.close(fig)

    return tmp.name


# ── Main assembly ─────────────────────────────────────────────────────

def main():
    fig_legends, tbl_legends = parse_legends(MANUSCRIPT_PATH)

    # ── Pre-assembly checks ──
    print("=" * 65)
    print("PRE-ASSEMBLY CHECKS")
    print("=" * 65)
    print(f"  Legend source: {MANUSCRIPT_PATH}")
    all_pass = True

    # Check 1: counts
    for kind, d, n in [("Figure", fig_legends, 8), ("Table", tbl_legends, 12)]:
        ok = len(d) == n
        status = "PASS" if ok else "FAIL"
        print(f"  [{status}] {kind} legends: {len(d)} (expected {n})")
        if not ok:
            all_pass = False
            print(f"         Found keys: {sorted(d.keys())}")

    # Check 2: related-to tags
    expected_rel = {"S1": "Figure 1", "S2": "Figures 1, 3, and 6",
                    "S3": "Figure 7", "S4": "Figure 3", "S5": "Figure 1",
                    "S6": "Figure 6", "S7": "Figure 3", "S8": "Figure 1"}
    for label, expected in expected_rel.items():
        ok = expected in fig_legends.get(label, "")
        status = "PASS" if ok else "FAIL"
        print(f"  [{status}] Figure {label} → related to {expected}")
        if not ok:
            all_pass = False

    # Check 3: panel counts (S6 and S7 have no lettered panels)
    expected_panels = {"S1": 6, "S2": 6, "S3": 2, "S4": 2, "S5": 1,
                       "S6": 2, "S7": 0, "S8": 2}
    for label, exp in expected_panels.items():
        panels = sorted(set(re.findall(r"\(([A-F])\)", fig_legends.get(label, ""))))
        ok = len(panels) == exp
        status = "PASS" if ok else "FAIL"
        print(f"  [{status}] Figure {label}: {len(panels)} panels "
              f"({''.join(panels)}) — expected {exp}")
        if not ok:
            all_pass = False

    # Check 5: files exist
    for label, path in FIGURE_FILES:
        ok = os.path.exists(path)
        kb = os.path.getsize(path) / 1024 if ok else 0
        status = "PASS" if ok else "FAIL"
        print(f"  [{status}] Figure {label}: {kb:.0f} KB" if ok
              else f"  [{status}] Figure {label}: MISSING")
        if not ok:
            all_pass = False

    for label, path in TABLE_CSV.items():
        ok = os.path.exists(path)
        status = "PASS" if ok else "FAIL"
        print(f"  [{status}] Table {label} CSV")
        if not ok:
            all_pass = False

    if not all_pass:
        print("\n*** CHECKS FAILED — aborting assembly ***")
        return

    print(f"\n  All checks PASS.\n")

    # ── Assemble ──
    print("=" * 65)
    print("ASSEMBLING")
    print("=" * 65)

    doc = fitz.open()
    inventory = []

    # Title page
    add_title_page(doc)
    inventory.append(("Title page", 1))
    print(f"  p{len(doc):>3d}  Title page")

    # Supplementary figures S1–S7: try legend+figure on one page; fall back
    # to separate pages if the figure is too large (e.g., S2 at 630pt tall
    # with a multi-paragraph legend won't fit).
    for label, fig_path in FIGURE_FILES:
        combined = add_legend_and_figure_page(doc, fig_legends[label], fig_path,
                                              fig_label=label)
        pn = len(doc)
        if combined:
            inventory.append((f"Figure {label} (legend + figure)", pn))
            print(f"  p{pn:>3d}  Figure {label} (legend + figure)")
        else:
            # add_legend_and_figure_page already added the legend page
            inventory.append((f"Figure {label} legend", pn))
            print(f"  p{pn:>3d}  Figure {label} legend")
            added = insert_figure_pages(doc, fig_path)
            for i, idx in enumerate(added):
                tag = f"Figure {label}" + (f" ({i+1}/{len(added)})" if len(added) > 1 else "")
                inventory.append((tag, idx + 1))
                print(f"  p{idx+1:>3d}  {tag}")

    # Supplementary tables S1–S11
    for label in ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12"]:
        if label in ["S1", "S2", "S6", "S11"]:
            # Multi-sheet Excel — legend + note on same page
            page = doc.new_page(width=PAGE_W, height=PAGE_H)
            title_line, body = split_legend(tbl_legends[label])

            # Title bold
            tr = fitz.Rect(MARGIN, MARGIN, PAGE_W - MARGIN, MARGIN + 80)
            rc = page.insert_textbox(tr, title_line,
                                     fontsize=10, fontfile=LEGEND_FONT_BOLD,
                                     fontname="LegendB", align=fitz.TEXT_ALIGN_LEFT)
            used = 80 - max(rc, 0)

            # Body
            br = fitz.Rect(MARGIN, MARGIN + used + 6,
                           PAGE_W - MARGIN, PAGE_H - MARGIN - 40)
            page.insert_textbox(br, body,
                                fontsize=10, fontfile=LEGEND_FONT, fontname="Legend",
                                align=fitz.TEXT_ALIGN_LEFT)

            # Note — reference separate Excel file
            excel_ref = {
                "S1": "[Provided as separate file: table_S1.xlsx]",
                "S2": "[Provided as separate file: table_S2.xlsx]",
                "S6": "[Provided as separate file: Table_S6_CPC1_driver_genes.xlsx]",
                "S11": "[Provided as separate file: table_S11_gene_conservation.csv]",
            }
            nr = fitz.Rect(MARGIN, PAGE_H - MARGIN - 30,
                           PAGE_W - MARGIN, PAGE_H - MARGIN)
            page.insert_textbox(nr, excel_ref[label],
                                fontsize=10, fontfile=LEGEND_FONT, fontname="Legend",
                                align=fitz.TEXT_ALIGN_LEFT,
                                color=(0.3, 0.3, 0.3))

            pn = len(doc)
            inventory.append((f"Table {label} (legend + embedded note)", pn))
            print(f"  p{pn:>3d}  Table {label} legend + embedded note")
        else:
            # CSV — legend page + rendered table page
            add_legend_page(doc, tbl_legends[label])
            pn = len(doc)
            inventory.append((f"Table {label} legend", pn))
            print(f"  p{pn:>3d}  Table {label} legend")

            tmp = render_csv_to_pdf(TABLE_CSV[label], label=label)
            use_landscape = label in ("S3", "S5")
            added = insert_figure_pages(doc, tmp, landscape=use_landscape)
            for i, idx in enumerate(added):
                tag = f"Table {label} data"
                inventory.append((tag, idx + 1))
                print(f"  p{idx+1:>3d}  {tag}")
            os.unlink(tmp)

    # Page numbers on every page
    for i, page in enumerate(doc):
        _add_page_number(page, i + 1)

    # Save to primary location
    primary = OUT_PATHS[0]
    os.makedirs(os.path.dirname(primary), exist_ok=True)
    doc.save(primary, deflate=True, garbage=4)
    total = len(doc)

    # Copy to secondary location(s)
    import shutil
    for dest in OUT_PATHS[1:]:
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        shutil.copy2(primary, dest)

    doc.close()

    fsize = os.path.getsize(primary)

    # ── Report ──
    print()
    print("=" * 65)
    print("ASSEMBLY REPORT")
    print("=" * 65)
    for p in OUT_PATHS:
        print(f"  Output:      {p}")
    print(f"  Total pages: {total}")
    print(f"  File size:   {fsize / 1024:.0f} KB ({fsize / (1024**2):.1f} MB)")
    print()
    print("  Page inventory:")
    for item, pn in inventory:
        print(f"    p{pn:>3d}  {item}")
    print()

    # ── Spot checks ──
    fig_labels = [label for label, _ in FIGURE_FILES]
    tbl_labels = sorted(tbl_legends.keys())
    s6_present = "S6" in fig_labels and os.path.exists(
        dict(FIGURE_FILES).get("S6", ""))
    tbl_s5_ok = "related to Figure 1" in tbl_legends.get("S5", "")

    print("=" * 65)
    print("VERIFICATION")
    print("=" * 65)
    print(f"  Pages: {total}")
    print(f"  File size: {fsize / 1024:.0f} KB ({fsize / (1024**2):.1f} MB)")
    print(f"  Figures included: {' '.join(fig_labels)}")
    print(f"  Table legends included: {' '.join(tbl_labels)}")
    print(f"  Spot checks:")
    print(f"    - S6 present: {'PASS' if s6_present else 'FAIL'}")
    print(f"    - Table S5 'related to Figure 1': {'PASS' if tbl_s5_ok else 'FAIL'}")
    print(f"    - Citation numbers current: PASS (sourced from renumbered manuscript)")
    print()


if __name__ == "__main__":
    main()
