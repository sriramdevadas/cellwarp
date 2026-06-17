#!/usr/bin/env python3
"""Overlay corrected p-value annotations on figS6_cellmarker_enrichment.pdf.

Replaces the two 'p < 1e-6' labels with the actual computed p-values:
  Primary (4.49-fold):            p = 2.1e-13  (2.1×10⁻¹³)
  Expression-matched (3.22-fold): p = 1.2e-12  (1.2×10⁻¹²)

Strategy: white-rectangle over existing text, then redraw using Helvetica
at matching size. Leaves the underlying plot intact.
"""
import shutil
from pathlib import Path

import fitz  # PyMuPDF

PROJECT = Path(__file__).resolve().parent.parent
TARGETS = [
    PROJECT / "figures/submission/supplementary/figS6_cellmarker_enrichment.pdf",
    PROJECT / "figures/supplementary/figS6_cellmarker_enrichment.pdf",
    PROJECT / "docs/submission/figures_for_review/Figure_S6.pdf",
]

# Located via page.get_text("dict"):
#   Primary "p < 1e-6"            at bbox (48.90, 19.12, 74.78, 26.94)
#   Expression-matched "p < 1e-6" at bbox (152.52, 42.72, 178.40, 50.54)
# ASCII e-notation matches the existing panel rendering convention
# (the panel already used '1e-6' / 'p <' form). Unicode superscripts
# don't render in the embedded core font.
REPLACEMENTS = [
    {
        "bbox": (48.90, 19.11, 74.78, 26.95),
        "new_text": "p = 2.1e-13",
    },
    {
        "bbox": (152.52, 42.71, 178.40, 50.55),
        "new_text": "p = 1.2e-12",
    },
]
FONT_SIZE = 6.5


def patch_pdf(path: Path):
    if not path.exists():
        print(f"  skip (missing): {path}")
        return
    tmp = path.with_suffix(".tmp.pdf")
    shutil.copy(path, tmp)
    doc = fitz.open(str(tmp))
    page = doc[0]
    # First pass: add redaction annots for the old text
    for rep in REPLACEMENTS:
        x0, y0, x1, y1 = rep["bbox"]
        pad = 0.6
        rect = fitz.Rect(x0 - pad, y0 - pad, x1 + pad, y1 + pad)
        page.add_redact_annot(rect, fill=(1, 1, 1))
    # Apply redactions — this removes the underlying text content
    page.apply_redactions()
    # Second pass: insert new text at each location
    for rep in REPLACEMENTS:
        x0, y0, x1, y1 = rep["bbox"]
        page.insert_text(
            fitz.Point(x0 - 2, y1 - 1),
            rep["new_text"],
            fontname="helv",
            fontsize=FONT_SIZE,
            color=(0, 0, 0),
        )
    doc.save(str(path), incremental=False, deflate=True)
    doc.close()
    tmp.unlink()
    print(f"  patched: {path}")


def main():
    for t in TARGETS:
        patch_pdf(t)


if __name__ == "__main__":
    main()
