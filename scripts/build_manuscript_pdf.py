#!/usr/bin/env python3
"""Build the submission/preprint manuscript PDF from the .docx, with neutral
document metadata.

The .pdf is a derived render (gitignored, like the supplementary PDF) but it is
uploaded to the journal portal and posted to bioRxiv, so its embedded metadata
ships publicly. docx2pdf drives MS Word, whose macOS PDF export can stamp the
machine's Office user / company / app-version into the PDF Info dict and XMP
packet. This script regenerates the PDF and then overwrites that metadata with a
single neutral author and title, clears the producer/creator fingerprint, and
strips the XMP packet — so no machine user, employer, or institutional identity
is ever embedded. Reproducible: re-run after any docx rebuild.

    .venv/bin/python scripts/build_manuscript_pdf.py
"""
from pathlib import Path
import os
import tempfile

import fitz  # PyMuPDF
from docx2pdf import convert

REPO = Path(__file__).resolve().parent.parent
DOCX = REPO / "docs/submission/cellwarp_manuscript.docx"
PDF = REPO / "docs/submission/cellwarp_manuscript.pdf"

TITLE = "Quantifying the Conserved Geometry of Cell-Type Identity Across Mammalian Species"
AUTHOR = "Sriram Devadas"

# Neutral Info dict: a single human author + the manuscript title; every other
# field blank (no producer/creator app-version fingerprint, no company).
CLEAN_META = {
    "title": TITLE,
    "author": AUTHOR,
    "subject": "",
    "keywords": "",
    "creator": "",
    "producer": "",
    "trapped": "",
    "format": "",
    "encryption": "",
}


def main():
    assert DOCX.exists(), f"docx not found: {DOCX}"
    convert(str(DOCX), str(PDF))
    assert PDF.exists(), f"docx2pdf did not produce {PDF}"

    doc = fitz.open(str(PDF))
    doc.set_metadata(CLEAN_META)          # overwrite the Info dictionary
    doc.del_xml_metadata()                # drop the XMP packet (Word author/company)
    tmp = tempfile.NamedTemporaryFile(suffix=".pdf", delete=False, dir=str(PDF.parent))
    tmp.close()
    doc.save(tmp.name, garbage=4, deflate=True)   # full rewrite, no incremental trailer
    doc.close()
    os.replace(tmp.name, str(PDF))

    # Verify
    v = fitz.open(str(PDF))
    m = v.metadata
    print(f"PDF rebuilt: {PDF}  ({PDF.stat().st_size} bytes, {v.page_count} pages)")
    print(f"  author   : {m.get('author')!r}")
    print(f"  title    : {m.get('title')!r}")
    print(f"  producer : {m.get('producer')!r}")
    print(f"  creator  : {m.get('creator')!r}")
    v.close()


if __name__ == "__main__":
    main()
