#!/usr/bin/env python3
"""One-time surgical patch for Figure S2 panel F obs/null annotations.

WHAT
    Replaces two panel-F obs/null text annotations in the canonical
    Figure S2 PDF, in place, via PyMuPDF text surgery:
        pan-Census : "obs/null = 0.813" -> "obs/null = 0.811"
        Andrews    : "obs/null = 0.818" -> "obs/null = 0.797"
    The 0.813/0.818 values were mean-normalized; 0.811/0.797 are the
    median-normalized canonical values that agree with the manuscript
    body, legend, Table 1, the CROSSWALK, and the source JSONs. This
    aligns the rendered figure with scripts/56_add_figs2_panel_f.py
    L77/L85 (edited in the same commit).

WHY SURGERY INSTEAD OF A CHAIN REBUILD
    The canonical figS2 build chain is:
        scripts/build_submission_figures.py (TASK A: assembles panels
        A-E from panel PNGs + negative_control_distributions.pdf)
          -> scripts/56_add_figs2_panel_f.py (appends panel F)
    Running it would produce the correct figure, but:
      * build_submission_figures.py is a multi-task orchestrator that
        also rebuilds figS1/figS4/figS5 in the same run; it has no
        figS2-only entry point, so running it whole risks touching
        unrelated supplementary figures.
      * 56_add_figs2_panel_f.py reads and overwrites the same canonical
        path it appends to, so it is only idempotent when run directly
        after the A-E base is regenerated; running it standalone
        doubles panel F.
    Direct text surgery on the two annotation spans is the minimum
    blast radius: it touches only those two strings and leaves every
    other figure, panel, and annotation byte-for-byte alone. The
    effect is equivalent to re-running the chain with script 56's
    L77/L85 already edited.

HOW
    For each target the script locates the exact text span (via
    get_text("dict")), records its baseline origin, font size, and
    color, redacts only the old glyphs (line art and images are
    preserved so the opaque white annotation bbox and the histogram
    stay intact), then reinserts the replacement string at the same
    baseline/size/color. Helvetica is used for reinsertion; it is
    metric-compatible with the original ArialMT at this size, so the
    rendered width and position match.

IDEMPOTENT
    A second run is a no-op: each target is detected as already
    patched (old string absent, new string present) and skipped; the
    file is not rewritten, so its bytes (and md5) are unchanged.
"""

import os
import sys
import tempfile

import fitz  # PyMuPDF

PDF_PATH = "figures/submission/supplementary/figS2_parameter_protocol_sensitivity.pdf"

REPLACEMENTS = [
    ("obs/null = 0.813", "obs/null = 0.811", "pan-Census"),
    ("obs/null = 0.818", "obs/null = 0.797", "Andrews"),
]


def find_spans(page, text):
    """Return [(bbox Rect, origin (x,y), size, color_int, font), ...] for
    spans whose text matches `text` exactly."""
    out = []
    for block in page.get_text("dict").get("blocks", []):
        for line in block.get("lines", []):
            for span in line.get("spans", []):
                if span.get("text", "") == text:
                    out.append((
                        fitz.Rect(span["bbox"]),
                        span["origin"],
                        span["size"],
                        span.get("color", 0),
                        span.get("font", ""),
                    ))
    return out


def int_to_rgb(color_int):
    """Convert a PyMuPDF sRGB integer to a (r, g, b) float tuple in [0, 1]."""
    return (
        (color_int >> 16 & 0xFF) / 255.0,
        (color_int >> 8 & 0xFF) / 255.0,
        (color_int & 0xFF) / 255.0,
    )


def main():
    doc = fitz.open(PDF_PATH)
    page = doc[0]

    # Phase 1 — resolve every target before mutating anything.
    todo = []
    for find, repl, label in REPLACEMENTS:
        old_spans = find_spans(page, find)
        new_spans = find_spans(page, repl)
        if not old_spans and len(new_spans) == 1:
            print(f"[{label}] already patched ({repl!r} present); skipping.")
            continue
        if len(old_spans) != 1:
            doc.close()
            raise RuntimeError(
                f"[{label}] expected exactly 1 span for {find!r}; "
                f"got {len(old_spans)} (replacement {repl!r}: {len(new_spans)})"
            )
        todo.append((label, find, repl, old_spans[0]))

    if not todo:
        doc.close()
        print("No changes; canonical already patched.")
        return 0

    # Phase 2 — redact only the old glyphs; keep line art + images so the
    # opaque white annotation bbox and the histogram drawings survive.
    for _label, _find, _repl, (bbox, _origin, _size, _color, _font) in todo:
        page.add_redact_annot(bbox, fill=(1, 1, 1))
    page.apply_redactions(
        images=fitz.PDF_REDACT_IMAGE_NONE,
        graphics=fitz.PDF_REDACT_LINE_ART_NONE,
        text=fitz.PDF_REDACT_TEXT_REMOVE,
    )

    # Phase 3 — reinsert each replacement at the original baseline.
    for label, find, repl, (bbox, origin, size, color, font) in todo:
        rgb = int_to_rgb(color)
        page.insert_text(origin, repl, fontname="helv", fontsize=size, color=rgb)
        print(
            f"[{label}] patched {find!r} -> {repl!r} "
            f"at origin=({origin[0]:.2f}, {origin[1]:.2f}) "
            f"size={size:.3f} color={tuple(round(c, 3) for c in rgb)} "
            f"(orig font {font})"
        )

    # Write back atomically (temp file in same dir + os.replace).
    fd, tmp_path = tempfile.mkstemp(suffix=".pdf", dir=os.path.dirname(PDF_PATH))
    os.close(fd)
    doc.save(tmp_path, deflate=True, garbage=4)
    doc.close()
    os.replace(tmp_path, PDF_PATH)
    print(f"Wrote {PDF_PATH}")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
