#!/usr/bin/env python3
"""Assemble the internal-review PDF: manuscript text (paginated) + the five main
figures. No LaTeX/LibreOffice needed — text pages via matplotlib, figures appended
with pdfunite. Review-artifact quality (monospace text); not the typeset submission."""
import subprocess, textwrap, pathlib, shutil
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

HERE = pathlib.Path(__file__).resolve().parent
TXT = HERE/"manuscript_combined.txt"
FIGDIR = HERE/"figures"
FIGS = [FIGDIR/f"{n}.pdf" for n in
        ["Fig1_configuration_conserved","Fig2_two_layers_bg","Fig3_configuration_robust",
         "Fig4_pertype_not_resolvable","Fig5_conserved_identity_genes"]]
OUT = HERE/"CellWarp_PLOSONE_review.pdf"
TEXT_PDF = HERE/"_manuscript_text.pdf"

WRAP = 96          # chars per line
LINES_PER_PAGE = 60
FONT = "monospace"; SIZE = 7.2

def paginate(txt):
    out = []
    for para in txt.split("\n"):
        if para.strip() == "":
            out.append(""); continue
        wrapped = textwrap.wrap(para, WRAP, break_long_words=False, break_on_hyphens=False) or [""]
        out.extend(wrapped)
    # split into pages
    pages = [out[i:i+LINES_PER_PAGE] for i in range(0, len(out), LINES_PER_PAGE)]
    return pages

def render_text_pdf(pages):
    with PdfPages(TEXT_PDF) as pdf:
        for pg in pages:
            fig = plt.figure(figsize=(8.5, 11)); fig.patch.set_facecolor("white")
            fig.text(0.08, 0.955, "\n".join(pg), family=FONT, fontsize=SIZE, va="top", ha="left")
            pdf.savefig(fig); plt.close(fig)

def main():
    pages = paginate(TXT.read_text())
    render_text_pdf(pages)
    print(f"text pages: {len(pages)}")
    parts = [str(TEXT_PDF)] + [str(f) for f in FIGS if f.exists()]
    missing = [f.name for f in FIGS if not f.exists()]
    if missing: print("WARNING missing figure PDFs:", missing)
    if shutil.which("pdfunite"):
        subprocess.run(["pdfunite"] + parts + [str(OUT)], check=True)
    else:
        subprocess.run(["gs","-dBATCH","-dNOPAUSE","-q","-sDEVICE=pdfwrite",
                        f"-sOutputFile={OUT}"] + parts, check=True)
    TEXT_PDF.unlink(missing_ok=True)
    print("wrote", OUT, "(", len(parts), "parts: text +", len(parts)-1, "figures )")

if __name__ == "__main__":
    main()
