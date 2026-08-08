#!/usr/bin/env python3
"""
Add panel F to Figure S2: replication inventory.

Panel F shows mini permutation null distributions for all 6 attempted replications
(Sun 2023, PanSci, CellHint, pan-Census; Andrews, MCA × HCA), plus a paired
within-human diagnostic sub-panel (HCA × Tabula Sapiens) attached visually to the
MCA × HCA cell.

Output: regenerates figures/submission/supplementary/figS2_parameter_protocol_sensitivity.pdf
with panels A-E unchanged and a new panel F appended below.
"""

import json
import os
from pathlib import Path

import fitz
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

# ── Style ──────────────────────────────────────────────────────────
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["pdf.fonttype"] = 42  # embed TrueType (vector text)
mpl.rcParams["ps.fonttype"] = 42

PROJECT = Path(__file__).resolve().parent.parent
SUPP_DIR = PROJECT / "figures/submission/supplementary"
EXISTING_PDF = SUPP_DIR / "figS2_parameter_protocol_sensitivity.pdf"  # current A-E PDF
OUT_PDF = SUPP_DIR / "figS2_parameter_protocol_sensitivity.pdf"  # rewrite same path
TMP_PANEL_F = Path("/tmp/figS2_panel_f.pdf")

# Cell Systems §8b: avoid red-green pairs.
# Using muted teal for successes, muted amber for failures, neutral gray for diagnostic.
COLOR_SUCCESS = "#2C7A7B"   # muted teal
COLOR_FAILURE = "#C05621"   # muted amber/burnt orange
COLOR_DIAG    = "#4A5568"   # neutral slate gray
COLOR_OBS     = "#1A202C"   # near-black for observed line
COLOR_HIST    = "#A0AEC0"   # light gray histogram fill
COLOR_HIST_F  = "#FBD38D"   # light amber for failure histograms
COLOR_HIST_D  = "#CBD5E0"   # very light gray for diagnostic histogram

# 6 main + 1 diagnostic (HCA × TS)
DATASETS = [
    {
        "key": "sun",
        "title": "Sun et al. 2023",
        "type": "success",
        "json": PROJECT / "output/validation/sun2023_replication_expanded/sun2023_expanded.json",
        "npy":  PROJECT / "output/validation/sun2023_replication_expanded/null_distribution.npy",
        "obs":  34.81, "obs_null": 0.554, "p": 0.0001, "n": 15,
    },
    {
        "key": "pansci",
        "title": "PanSci",
        "type": "success",
        "json": PROJECT / "output/validation/pansci_replication/pansci_replication.json",
        "npy":  PROJECT / "output/validation/pansci_replication/null_distribution.npy",
        "obs":  38.32, "obs_null": 0.552, "p": 0.0001, "n": 16,
    },
    {
        "key": "cellhint",
        "title": "CellHint",
        "type": "success",
        "json": PROJECT / "output/validation/cellhint_replication/cellhint_replication.json",
        "npy":  PROJECT / "output/validation/cellhint_replication/null_distribution.npy",
        "obs":  28.15, "obs_null": 0.448, "p": 0.0001, "n": 15,
    },
    {
        "key": "pancensus",
        "title": "pan-Census",
        "type": "success",
        "json": PROJECT / "analysis/census_replication/replication_results.json",
        "npy":  PROJECT / "analysis/census_replication/null_distribution.npy",
        "obs":  60.14, "obs_null": 0.811, "p": 0.0001, "n": 22,
    },
    {
        "key": "andrews",
        "title": "Andrews et al.",
        "type": "failure",
        "json": PROJECT / "output/validation/andrews_replication/andrews_replication_results.json",
        "npy":  PROJECT / "output/validation/andrews_replication/null_distribution.npy",
        "obs":  10.84, "obs_null": 0.797, "p": 0.1159, "n": 6,
    },
    {
        "key": "mca_hca",
        "title": "MCA × HCA",
        "type": "failure",
        "json": PROJECT / "output/validation/t1a_replication/t1a_results.json",
        "npy":  PROJECT / "output/validation/t1a_replication/null_distribution.npy",
        "obs":  32.82, "obs_null": 1.003, "p": 0.542, "n": 17,
    },
]

DIAGNOSTIC = {
    "key": "hca_ts",
    "title": "HCA × Tabula Sapiens",
    "subtitle": "(within-human diagnostic)",
    "json": PROJECT / "output/validation/hca_centroid_comparison/hca_centroid_comparison.json",
    "npy":  PROJECT / "output/validation/hca_centroid_comparison/null_a_hca_vs_tabula.npy",
    "obs":  14.32, "obs_null": 0.728, "p": 0.003, "n": 6,
}


def load_null(npy_path: Path) -> np.ndarray:
    """Load 10k-element null distribution array; sanity-check shape."""
    arr = np.load(npy_path)
    if arr.ndim != 1:
        raise ValueError(f"{npy_path}: expected 1D array, got {arr.shape}")
    if arr.size != 10000:
        print(f"  WARN: {npy_path.name} has {arr.size} elements (expected 10000)")
    return arr


def format_p(p: float) -> str:
    # Permutation-test floor (p ≤ 1/(10000+1) ≈ 0.0001) rendered as mathtext
    # "$p < 10^{-4}$" (honest to floor, not pseudo-exact).
    if p <= 0.0001:
        return r"$p < 10^{-4}$"
    if p < 0.01:
        return f"p = {p:.4f}"
    return f"p = {p:.3f}"


def sig_marker(p: float) -> str:
    if p <= 0.0001:
        return "***"
    if p <= 0.01:
        return "**"
    if p <= 0.05:
        return "*"
    return "n.s."


def plot_mini_panel(ax, ds: dict, *, color_hist: str, edge_color: str, is_diagnostic: bool = False):
    null = load_null(ds["npy"])
    obs = ds["obs"]

    # KDE-like histogram: normalized, 40 bins
    counts, _bin_edges, _patches = ax.hist(
        null, bins=40, density=True,
        color=color_hist, edgecolor="#4A5568", linewidth=0.4,
        alpha=0.85,
    )
    # Vertical headroom: bottom ~2/3 holds histogram, top ~1/3 reserved for the
    # obs/null + p annotation bbox at y=0.94. Critical for Andrews / MCA × HCA
    # where histograms span the full x-range and upper-left corner still has
    # bars at non-trivial density. 1.5× peak chosen so label bbox (which extends
    # roughly y∈[0.83, 0.96] in axes coords given ~7pt 2-line text) sits in
    # clear space above peak (peak occupies y∈[0, 0.667]).
    max_density = float(counts.max()) if len(counts) > 0 else 1.0
    ax.set_ylim(0, max_density * 1.5)
    # Observed vertical line
    ax.axvline(obs, color=COLOR_OBS, linewidth=1.4, linestyle="-")

    # Title (dataset name + n)
    title = f"{ds['title']} (n = {ds['n']})"
    if is_diagnostic:
        title = f"{ds['title']}\n{ds['subtitle']} (n = {ds['n']})"
    ax.set_title(title, fontsize=8, pad=2.5,
                 fontweight="normal" if is_diagnostic else "bold")

    # Bottom annotation: obs/null and p — upper-LEFT (was upper-right where
    # histograms peak). For the success panels (Sun/PanSci/CellHint/pan-Census)
    # the obs line sits far left and histograms sit far right; upper-left is
    # mostly empty. For Andrews / MCA × HCA the obs sits inside the histogram
    # but upper-left is still less dense than upper-right.
    pmark = sig_marker(ds["p"])
    ann = f"obs/null = {ds['obs_null']:.3f}\n{format_p(ds['p'])} {pmark}"
    ax.text(
        0.03, 0.94, ann,
        transform=ax.transAxes, fontsize=7,
        ha="left", va="top", zorder=10,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                  edgecolor="#A0AEC0", linewidth=0.4, alpha=1.0),
    )

    # Style
    ax.set_xlabel("Procrustes distance", fontsize=7, labelpad=1.5)
    ax.set_ylabel("Density", fontsize=7, labelpad=1.5)
    ax.tick_params(axis="both", labelsize=6, length=2.5, pad=1.5)
    for s in ax.spines.values():
        s.set_linewidth(0.5)

    # Edge color frame to distinguish success/failure/diagnostic
    for spine_name in ("top", "right", "bottom", "left"):
        ax.spines[spine_name].set_edgecolor(edge_color)
        ax.spines[spine_name].set_linewidth(1.4 if is_diagnostic else 1.6)
    if is_diagnostic:
        # Dashed border to subordinate the diagnostic
        for spine_name in ("top", "right", "bottom", "left"):
            ax.spines[spine_name].set_linestyle("--")


def build_panel_f_pdf():
    """Generate panel F as a standalone vector PDF."""
    # Layout: title bar at top + 3×2 grid (Sun/PanSci/CellHint on row 1; Pan-Census/Andrews/MCA×HCA on row 2)
    # + paired diagnostic sub-panel attached visually below MCA×HCA cell.
    # Using gridspec for fine control.

    # Panel F target dimensions: width = ~136mm (matches existing figS2 page width)
    # Width inches = 136 / 25.4 ≈ 5.35"
    # Height: 2 rows × ~55mm + subordinate row ~30mm ≈ 5.5" (incl. caption space)
    fig = plt.figure(figsize=(5.35, 4.55), dpi=300)

    # Create gridspec: 3 rows, 3 cols. Last row spans only col 3 for diagnostic.
    # hspace=0.85 + row-3 ratio reduced 0.75 → 0.50: prior 0.60 still left row 3
    # encroaching row 2 axes. The diagnostic is a single small panel and
    # doesn't need to match the height of the main rows; shrinking row 3
    # plus more hspace gives it a clear band beneath row 2.
    gs = fig.add_gridspec(
        nrows=3, ncols=3,
        left=0.075, right=0.97,
        top=0.94, bottom=0.07,
        hspace=0.85, wspace=0.42,
        height_ratios=[1.0, 1.0, 0.50],
    )

    # Row 1: 3 successes (Sun, PanSci, CellHint)
    # Row 2: pan-Census (success), Andrews (failure), MCA×HCA (failure)
    grid_assignments = [
        # (dataset, row, col)
        (DATASETS[0], 0, 0),  # Sun
        (DATASETS[1], 0, 1),  # PanSci
        (DATASETS[2], 0, 2),  # CellHint
        (DATASETS[3], 1, 0),  # pan-Census
        (DATASETS[4], 1, 1),  # Andrews
        (DATASETS[5], 1, 2),  # MCA × HCA
    ]

    edge_colors = {"success": COLOR_SUCCESS, "failure": COLOR_FAILURE}
    hist_colors = {"success": COLOR_HIST,    "failure": COLOR_HIST_F}

    for ds, r, c in grid_assignments:
        ax = fig.add_subplot(gs[r, c])
        plot_mini_panel(
            ax, ds,
            color_hist=hist_colors[ds["type"]],
            edge_color=edge_colors[ds["type"]],
        )

    # Diagnostic sub-panel: row 2 (under MCA×HCA), col 2, with explicit "paired" arrow
    ax_diag = fig.add_subplot(gs[2, 2])
    plot_mini_panel(
        ax_diag, DIAGNOSTIC,
        color_hist=COLOR_HIST_D,
        edge_color=COLOR_DIAG,
        is_diagnostic=True,
    )

    # Connector from MCA × HCA (gs[1,2]) to diagnostic (gs[2,2]).
    # ConnectionPatch in figure coords: dashed vertical line + arrow tip,
    # anchored at the bottom-center of MCA × HCA and the top-center of the
    # diagnostic, making the pairing visually unambiguous.
    from matplotlib.patches import ConnectionPatch
    fig.canvas.draw()  # populate axes positions before reading them
    mca_ax = fig.axes[5]  # 6th subplot = gs[1,2] = MCA × HCA
    mca_pos = mca_ax.get_position()
    diag_pos = ax_diag.get_position()
    conn_x = (mca_pos.x0 + mca_pos.x1) / 2
    conn = ConnectionPatch(
        xyA=(conn_x, mca_pos.y0),
        xyB=(conn_x, diag_pos.y1),
        coordsA="figure fraction", coordsB="figure fraction",
        arrowstyle="->", linestyle="--", linewidth=0.8,
        color=COLOR_DIAG, mutation_scale=8,
    )
    fig.add_artist(conn)
    # Brief italic label to right of the connector.
    label_y = (mca_pos.y0 + diag_pos.y1) / 2
    fig.text(
        conn_x + 0.02, label_y,
        "paired",
        ha="left", va="center",
        fontsize=6.5, style="italic", color=COLOR_DIAG,
    )

    # Legend strip at top: success / failure / diagnostic
    leg_ax = fig.add_axes([0.075, 0.965, 0.895, 0.025])
    leg_ax.axis("off")
    legend_specs = [
        ("Success (replicates)", COLOR_SUCCESS, "-"),
        ("Failure (does not replicate)", COLOR_FAILURE, "-"),
        ("Within-human diagnostic", COLOR_DIAG, "--"),
    ]
    for i, (label, color, ls) in enumerate(legend_specs):
        x_anchor = 0.0 + i * 0.34
        leg_ax.plot([x_anchor + 0.005, x_anchor + 0.04], [0.5, 0.5],
                    color=color, linewidth=1.6, linestyle=ls,
                    transform=leg_ax.transAxes, clip_on=False)
        leg_ax.text(x_anchor + 0.05, 0.5, label,
                    transform=leg_ax.transAxes, fontsize=7,
                    va="center", ha="left")

    fig.savefig(TMP_PANEL_F, format="pdf", bbox_inches=None, pad_inches=0.05)
    plt.close(fig)
    print(f"  Wrote panel F PDF: {TMP_PANEL_F} ({os.path.getsize(TMP_PANEL_F):,} bytes)")


def assemble_combined_pdf():
    """Combine existing figS2 (panels A-E) + new panel F into a single multi-page or
    single-page PDF with panel F appended below the existing content."""
    # Strategy: read existing FigS2 page dimensions; create a new wider/taller page
    # that embeds the existing PDF as the top section + panel F PDF below + 'F' label.

    src = fitz.open(str(EXISTING_PDF))
    src_page = src[0]

    # Entry assertion: this script reads and rewrites the same path it appends to,
    # so it is only correct on an A-E base. Run standalone against its own output it
    # would append panel F a second time -- the hazard this script's sibling
    # scripts/patch_figs2_panel_f_values.py documents in its docstring, and the
    # reason that patch used text surgery rather than a chain rebuild. Asserting the
    # input shape here is also the idempotency guard the chain previously lacked.
    _labels = sorted({
        s["text"].strip()
        for b in src_page.get_text("dict")["blocks"] if b["type"] == 0
        for l in b["lines"] for s in l["spans"]
        if len(s["text"].strip()) == 1 and s["text"].strip().isalpha()
        and s["text"].strip().isupper() and s["size"] >= 8.0
        and "bold" in s["font"].lower()
    })
    assert _labels == list("ABCDE"), (
        f"{EXISTING_PDF.name} carries panel labels {_labels}, expected "
        f"['A','B','C','D','E']. This script appends panel F to an A-E base; "
        f"running it against its own output would duplicate panel F. Rebuild the "
        f"A-E base with scripts/build_submission_figures.py first."
    )

    src_w, src_h = src_page.rect.width, src_page.rect.height
    print(f"  Existing FigS2 page: {src_w:.1f} × {src_h:.1f} pts (panels {''.join(_labels)})")

    pf_doc = fitz.open(str(TMP_PANEL_F))
    pf_page = pf_doc[0]
    pf_w, pf_h = pf_page.rect.width, pf_page.rect.height
    print(f"  Panel F PDF page:    {pf_w:.1f} × {pf_h:.1f} pts")

    # Scale panel F to match src content width
    LEFT = 3.6
    CONTENT_W = src_w - 2 * LEFT
    pf_scale = CONTENT_W / pf_w
    pf_display_h = pf_h * pf_scale
    LABEL_FS = 9.0
    LABEL_OFFSET = 3.4
    ROW_GAP = 6.0

    # Compose new page
    new_top_y = 0.0
    src_y0 = 0.0
    src_y1 = src_h
    label_f_y = src_y1 + ROW_GAP + LABEL_FS - 1
    pf_y0 = src_y1 + ROW_GAP + LABEL_FS + LABEL_OFFSET
    pf_y1 = pf_y0 + pf_display_h

    new_page_w = src_w
    new_page_h = pf_y1 + 4.0

    print(f"  New page: {new_page_w:.1f} × {new_page_h:.1f} pts")
    print(f"  Panel F: y={pf_y0:.1f}–{pf_y1:.1f}; label F at y={label_f_y:.1f}")

    new_doc = fitz.open()
    new_page = new_doc.new_page(width=new_page_w, height=new_page_h)

    # Embed existing figS2 (panels A-E) at top
    src_rect = fitz.Rect(0, src_y0, src_w, src_y1)
    new_page.show_pdf_page(src_rect, src, 0)

    # Embed panel F
    pf_rect = fitz.Rect(LEFT, pf_y0, LEFT + CONTENT_W, pf_y1)
    new_page.show_pdf_page(pf_rect, pf_doc, 0)

    # Add bold "F" label
    new_page.insert_text(
        fitz.Point(LEFT, label_f_y),
        "F",
        fontname="hebo",  # Helvetica-Bold
        fontsize=LABEL_FS,
        color=(0, 0, 0),
    )

    # Save to a new tmp file then move into place (avoid in-place corruption)
    tmp_out = Path("/tmp/figS2_with_panel_f.pdf")
    new_doc.save(str(tmp_out), garbage=4, deflate=True)
    new_doc.close()
    src.close()
    pf_doc.close()

    # Move into place
    os.replace(tmp_out, OUT_PDF)
    print(f"  Wrote combined PDF: {OUT_PDF} ({os.path.getsize(OUT_PDF):,} bytes)")


def main():
    print("=" * 70)
    print("Generating Figure S2 panel F + combined figS2 PDF")
    print("=" * 70)

    # Sanity-check all source files exist
    for ds in DATASETS:
        assert ds["json"].exists(), f"Missing JSON: {ds['json']}"
        assert ds["npy"].exists(), f"Missing npy: {ds['npy']}"
    assert DIAGNOSTIC["json"].exists(), f"Missing JSON: {DIAGNOSTIC['json']}"
    assert DIAGNOSTIC["npy"].exists(), f"Missing npy: {DIAGNOSTIC['npy']}"

    print("\n[1] Generating panel F vector PDF...")
    build_panel_f_pdf()

    print("\n[2] Assembling combined figS2 PDF (panels A-E + F)...")
    assemble_combined_pdf()

    print("\n[3] Final verification...")
    final = fitz.open(str(OUT_PDF))
    page = final[0]
    print(f"  Page size: {page.rect.width:.1f} × {page.rect.height:.1f} pts "
          f"({page.rect.width / 72 * 25.4:.1f} × {page.rect.height / 72 * 25.4:.1f} mm)")
    print(f"  File size: {os.path.getsize(OUT_PDF):,} bytes "
          f"({os.path.getsize(OUT_PDF) / 1024 / 1024:.2f} MB)")
    images = page.get_images(full=True)
    print(f"  Embedded images: {len(images)}")
    fonts = set()
    for b in page.get_text("dict")["blocks"]:
        if b.get("type") == 0:
            for line in b["lines"]:
                for span in line["spans"]:
                    fonts.add(span["font"])
    print(f"  Fonts used: {sorted(fonts)}")
    final.close()
    print("\nDone.")


if __name__ == "__main__":
    main()
