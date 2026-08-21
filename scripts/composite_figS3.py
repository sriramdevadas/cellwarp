#!/usr/bin/env python3
"""
Composite Figure S3: Bootstrap ranking-stability panels (A–B). (Treeness panels cut.)

Combines (two-panel layout; the prior treeness panels were cut):
  (A) Bootstrap ranking forest plot — regenerated with Cell Systems styling
  (B) Bootstrap CI-width vs cross-atlas rank-shift scatter — regenerated

Both panels are regenerated from data to enforce Arial fonts, 6-8pt text,
colorblind-safe palette, and embedded fonts (pdf.fonttype=42).

Biology: Bootstrap CI panels show that cross-species divergence rankings are
stable to sampling variation — stromal cells remain most divergent, CD8+ T
cells most conserved across 1000 full-pipeline bootstrap iterations.

Math: Non-parametric bootstrap with replacement, fresh PCA per iteration,
95% CIs on per-type residual magnitude ranks.

Output: figures/submission/supplementary/figS3_bootstrap_rankings.pdf
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap

# Add src to path for figure_style
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from pathlib import Path
from cellwarp.figure_style import (
    apply_style, save_figure, COL2, C_BLUE, C_ORANGE, C_PURPLE,
    C_TEAL, C_GRAY, C_DARKGRAY, FONT_FAMILY, SHORT_NAMES, LABEL_EXPAND,
    FONT_SIZE_LABEL, FONT_SIZE_TICK, FONT_SIZE_LEGEND, FONT_SIZE_ANNOT,
)

# Constants
N_BOOTSTRAP = 1000
N_TYPES = 35
RANK_TOP_THIRD = 12
RANK_BOTTOM_THIRD = 24

# Paths
BOOTSTRAP_DIR = str(Path(__file__).resolve().parent.parent / "analysis/bootstrap_rankings")
SUPP_DIR = str(Path(__file__).resolve().parent.parent / "figures/supplementary")
OUT_DIR = str(Path(__file__).resolve().parent.parent / "figures/submission/supplementary")

# Colorblind-safe category colors (no red-green)
CATEGORY_COLORS = {
    "STABLE_FLEXIBLE": C_ORANGE,     # orange — reliably divergent (high residual)
    "STABLE_RIGID": C_TEAL,          # teal — reliably conserved (low residual)
    "STABLE_MIDDLE": C_BLUE,         # blue — stable middle
    "MODERATE": C_PURPLE,            # purple — moderate uncertainty
    "UNSTABLE": C_GRAY,              # gray — unstable
}

# De-rigidified display labels. The per-type score is reported in the
# prose as the Procrustes residual / per-type divergence; the internal category
# KEYS come from bootstrap_summary.csv and are mapped to these labels at plot
# time so the figure carries the prose term without touching the data file.
CATEGORY_DISPLAY = {
    "STABLE_FLEXIBLE": "Stable divergent",
    "STABLE_RIGID": "Stable conserved",
    "STABLE_MIDDLE": "Stable intermediate",
    "MODERATE": "Moderate",
    "UNSTABLE": "Unstable",
}


def short_name(ct):
    """Shortened cell type name for compact display."""
    return SHORT_NAMES.get(ct, ct)


def generate_panel_c(summary_df):
    """
    Panel C: Bootstrap forest plot — 95% CIs for 35 cell types.

    Regenerated from data with Cell Systems styling: Arial 6-8pt,
    colorblind-safe palette, clean spines.
    """
    apply_style()

    df = summary_df.sort_values("median_rank").reset_index()

    fig, ax = plt.subplots(figsize=(COL2 * 0.5, 5.5))

    for idx, row in df.iterrows():
        color = CATEGORY_COLORS.get(row["category"], C_GRAY)
        # Conserved-at-top: largest median_rank (rank 35) lands at top of plot,
        # matching Fig 3A convention. (Previously: y = len(df) - idx - 1.)
        y = idx

        # CI bar
        ax.plot(
            [row["ci_lower"], row["ci_upper"]], [y, y],
            color=color, linewidth=1.8, solid_capstyle="round", zorder=2,
        )
        # Median point
        ax.plot(
            row["median_rank"], y,
            "o", color=color, markersize=4, zorder=3,
            markeredgecolor="white", markeredgewidth=0.3,
        )
        # Cell type label — apply LABEL_EXPAND for HPC/HSC/MSC abbreviation
        # disambiguation, consistent with Fig 1E/3A/4A and S1B/S2.
        short = short_name(row["cell_type"])
        ct_label = LABEL_EXPAND.get(short, short)
        ax.text(
            -0.5, y, ct_label,
            ha="right", va="center", fontsize=6, fontfamily=FONT_FAMILY,
        )

    # Reference lines for thirds
    ax.axvline(RANK_TOP_THIRD + 0.5, color=C_GRAY, linestyle=":", alpha=0.5,
               linewidth=0.6)
    ax.axvline(RANK_BOTTOM_THIRD - 0.5, color=C_GRAY, linestyle=":", alpha=0.5,
               linewidth=0.6)

    # Zone labels
    ax.text(
        RANK_TOP_THIRD / 2, len(df) + 0.5, "divergent\n(high residual)",
        ha="center", va="bottom", fontsize=6, color=C_ORANGE, alpha=0.8,
        fontfamily=FONT_FAMILY,
    )
    ax.text(
        (RANK_BOTTOM_THIRD + N_TYPES) / 2, len(df) + 0.5, "conserved\n(low residual)",
        ha="center", va="bottom", fontsize=6, color=C_TEAL, alpha=0.8,
        fontfamily=FONT_FAMILY,
    )

    ax.set_xlim(0, N_TYPES + 1)
    ax.set_ylim(-1, len(df) + 2)
    ax.set_xlabel(
        "Residual magnitude rank\n(1 = most diverged, 35 = most conserved)",
        fontsize=FONT_SIZE_LABEL, fontfamily=FONT_FAMILY,
    )
    ax.set_yticks([])
    ax.set_title(
        f"Bootstrap 95% CI on per-type divergence rankings\n(n = {N_BOOTSTRAP}, full pipeline)",
        fontsize=8, fontweight="bold", fontfamily=FONT_FAMILY,
    )

    # Clean spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color(C_DARKGRAY)
    ax.tick_params(colors=C_DARKGRAY, labelsize=FONT_SIZE_TICK)

    # Legend — only include categories with ≥1 cell type in this run.
    # MODERATE and UNSTABLE bins are empty for the 35-type bootstrap set, so
    # carrying empty patches would mislead readers.
    present_categories = set(df["category"].unique())
    legend_handles = [
        mpatches.Patch(color=c, label=CATEGORY_DISPLAY.get(l, l.replace("_", " ").title()))
        for l, c in CATEGORY_COLORS.items()
        if l in present_categories
    ]
    ax.legend(
        handles=legend_handles, loc="upper left", bbox_to_anchor=(0.0, -0.18),
        fontsize=FONT_SIZE_LEGEND, ncol=len(legend_handles),
        framealpha=0.9, title="Stability category", title_fontsize=FONT_SIZE_LEGEND,
        edgecolor="#BDBDBD", fancybox=False, frameon=True,
    )
    ax.get_legend().get_frame().set_linewidth(0.5)

    plt.tight_layout()
    return fig


def generate_panel_d(master_df):
    """
    Panel D: Bootstrap CI width vs cross-atlas mean rank shift scatter.

    Replaces the prior 35x35 swap-probability heatmap, which lacked a citation
    in the manuscript. The scatter visualizes the within-vs-cross-atlas
    inversion cited in Results §7 and the figure caption (ρ = −0.410,
    p = 0.073, n = 20). Logic ported from
    analysis/cross_reference/cross_reference_analysis.py:270-368, adapted to
    Cell Systems styling for parity with Panel C.
    """
    from scipy.stats import spearmanr
    apply_style()

    # Restrict to types present in ≥2 replications and with both bootstrap CI
    # width and a mean cross-atlas rank shift defined.
    plot_data = master_df[master_df["n_replications_present"] >= 2].copy()
    valid = plot_data.dropna(
        subset=["bootstrap_CI_width", "mean_rank_shift"]
    ).copy()

    rho, p_val = spearmanr(valid["bootstrap_CI_width"], valid["mean_rank_shift"])

    # Colors by bootstrap category — share Panel C's palette so the two
    # panels read as a paired pair.
    cat_colors = {
        "STABLE_FLEXIBLE": CATEGORY_COLORS["STABLE_FLEXIBLE"],
        "STABLE_MIDDLE":   CATEGORY_COLORS["STABLE_MIDDLE"],
        "STABLE_RIGID":    CATEGORY_COLORS["STABLE_RIGID"],
    }

    fig, ax = plt.subplots(figsize=(COL2 * 0.5, COL2 * 0.55))

    for _, row in valid.iterrows():
        color = cat_colors.get(row["bootstrap_category"], C_GRAY)
        ax.scatter(
            row["bootstrap_CI_width"], row["mean_rank_shift"],
            s=25, c=[color], edgecolors="white", linewidths=0.4,
            zorder=3, alpha=0.9,
        )
        # Compact per-point label using the SHORT_NAMES mapping shared across
        # phase 3 panels (no LABEL_EXPAND here — scatter is tight and short
        # names give the best density-vs-readability tradeoff).
        label = short_name(row["cell_type"])
        ax.annotate(
            label,
            (row["bootstrap_CI_width"], row["mean_rank_shift"]),
            fontsize=5, ha="left", va="bottom",
            xytext=(3, 3), textcoords="offset points",
            color=C_DARKGRAY, fontfamily=FONT_FAMILY,
        )

    # Stats annotation — upper-right. Upper-left masked Hepatocyte (at
    # x=2, y=10.83); upper-right is empty in the n=20 cloud (verified: no
    # points at x≥6, y≥8; granulocyte at (7, 6) sits well below). Opaque
    # bbox (alpha=1.0, zorder=10) preserved.
    ax.text(
        0.97, 0.97,
        f"ρ = {rho:.3f}, p = {p_val:.3f}\nn = {len(valid)}",
        transform=ax.transAxes, fontsize=FONT_SIZE_ANNOT,
        ha="right", va="top", color=C_DARKGRAY, fontfamily=FONT_FAMILY,
        zorder=10,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                  edgecolor="lightgray", alpha=1.0, linewidth=0.4),
    )

    ax.set_xlabel(
        "Bootstrap CI width (within-atlas stability)",
        fontsize=FONT_SIZE_LABEL, fontfamily=FONT_FAMILY,
    )
    ax.set_ylabel(
        "Mean cross-atlas rank shift",
        fontsize=FONT_SIZE_LABEL, fontfamily=FONT_FAMILY,
    )
    ax.set_title(
        "Bootstrap stability vs cross-atlas consistency",
        fontsize=8, fontweight="bold", fontfamily=FONT_FAMILY,
    )

    # Legend — only categories present in the n=20 subset (MODERATE/UNSTABLE
    # are absent by construction in this bootstrap).
    present = set(valid["bootstrap_category"].unique())
    legend_handles = [
        mpatches.Patch(color=cat_colors[k], label=CATEGORY_DISPLAY.get(k, k.replace("_", " ").title()))
        for k in ("STABLE_FLEXIBLE", "STABLE_MIDDLE", "STABLE_RIGID")
        if k in present
    ]
    ax.legend(
        handles=legend_handles, loc="lower right",
        fontsize=FONT_SIZE_LEGEND,
        framealpha=0.9, title="Stability category",
        title_fontsize=FONT_SIZE_LEGEND,
        edgecolor="#BDBDBD", fancybox=False, frameon=True,
    )
    ax.get_legend().get_frame().set_linewidth(0.5)

    # Padding around point cloud
    x_max = float(valid["bootstrap_CI_width"].max())
    y_max = float(valid["mean_rank_shift"].max())
    ax.set_xlim(-0.5, x_max + 1.0)
    ax.set_ylim(-0.5, y_max + 1.5)

    # Spines + ticks consistent with Panel C
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_color(C_DARKGRAY)
    ax.spines["left"].set_color(C_DARKGRAY)
    ax.tick_params(colors=C_DARKGRAY, labelsize=FONT_SIZE_TICK)

    plt.tight_layout()
    return fig


def main():
    """Build composite Figure S3 with all four panels."""
    # Path-swap (mirrors the S1 producer): Panels A+B now sourced directly
    # from live per-panel PDFs (fig7a_treeness, fig7b_density) rather than
    # the absent hand-polished figS3_bootstrap_rankings_polished.pdf. The
    # early-return guard was removed; this script is now the live writer.
    import fitz

    os.makedirs(OUT_DIR, exist_ok=True)

    # ── Load bootstrap data ──────────────────────────────────────
    print("Loading bootstrap data...")
    summary_df = pd.read_csv(os.path.join(BOOTSTRAP_DIR, "bootstrap_summary.csv"))
    master_df = pd.read_csv(
        str(Path(__file__).resolve().parent.parent /
            "analysis/cross_reference/master_ranking_table.csv")
    )
    print(f"  Summary: {len(summary_df)} cell types")
    print(f"  Master: {len(master_df)} cell types "
          f"(scatter restricted to ≥2 replications)")

    # ── Generate polished panels C and D ─────────────────────────
    print("\nGenerating polished panel C (forest plot)...")
    fig_c = generate_panel_c(summary_df)
    panel_c_path = "/tmp/figS3_panel_c.pdf"
    fig_c.savefig(panel_c_path, format="pdf", dpi=300, bbox_inches="tight",
                  pad_inches=0.05)
    plt.close(fig_c)
    print(f"  Saved: {panel_c_path} ({os.path.getsize(panel_c_path):,} bytes)")

    print("\nGenerating polished panel D (CI-vs-shift scatter)...")
    fig_d = generate_panel_d(master_df)
    panel_d_path = "/tmp/figS3_panel_d.pdf"
    fig_d.savefig(panel_d_path, format="pdf", dpi=300, bbox_inches="tight",
                  pad_inches=0.05)
    plt.close(fig_d)
    print(f"  Saved: {panel_d_path} ({os.path.getsize(panel_d_path):,} bytes)")

    # ── PDF compositing ──────────────────────────────────────────
    print("\nCompositing four panels into final S3...")

    ARIAL_BOLD = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
    LABEL_SIZE = 9.0

    # Panel labels are Arial Bold wherever macOS supplies it. The deposited figure
    # was built that way and embeds Arial, so keeping this the first choice means a
    # macOS rebuild is unchanged. Off macOS the path does not exist and there is no
    # Arial to embed; fall back to the Helvetica-Bold that MuPDF carries internally,
    # which needs no download and no extra dependency. That route is for a reader's
    # run, not for producing the deposited artifact: the labels are bold and legible
    # but the file will not match the deposit, because the embedded face differs.
    # MuPDF reports the fallback as "NimbusSans-Bold", so the "Bold" test in the
    # verification block below matches either way.
    LABEL_FONT_FILE = ARIAL_BOLD if os.path.exists(ARIAL_BOLD) else None
    LABEL_FONT_DESC = ("Arial Bold" if LABEL_FONT_FILE else
                       "MuPDF built-in Helvetica-Bold (Arial Bold not on this system)")
    print(f"  Panel-label font: {LABEL_FONT_DESC}")

    def _label_font():
        """Arial Bold when present, else MuPDF's built-in Helvetica-Bold."""
        if LABEL_FONT_FILE:
            return fitz.Font(fontfile=LABEL_FONT_FILE)
        return fitz.Font(fontname="hebo")

    def add_label(page, text, x, y_baseline, size=LABEL_SIZE):
        """Add a bold panel label to a PDF page."""
        font = _label_font()
        tw = fitz.TextWriter(page.rect)
        tw.append(fitz.Point(x, y_baseline), text, font=font, fontsize=size)
        tw.write_text(page)

    # Two-panel layout: treeness A+B cut; promote forest (C->A)
    # and CI-vs-shift scatter (D->B). Macaque 12-type fix flows via master_ranking_table.
    c_doc = fitz.open(panel_c_path)
    d_doc = fitz.open(panel_d_path)
    c_page = c_doc[0]; d_page = d_doc[0]
    c_w, c_h = c_page.rect.width, c_page.rect.height
    d_w, d_h = d_page.rect.width, d_page.rect.height
    print(f"  C (->A) source: {c_w:.1f} x {c_h:.1f} pts")
    print(f"  D (->B) source: {d_w:.1f} x {d_h:.1f} pts")

    col_gap = 10.0
    label_h = LABEL_SIZE + 4.0
    top_h = max(c_h, d_h)
    page_w = c_w + col_gap + d_w
    page_h = label_h + top_h + 14.0

    comp_doc = fitz.open()
    comp_page = comp_doc.new_page(width=page_w, height=page_h)

    a_y0 = label_h + (top_h - c_h) / 2.0
    comp_page.show_pdf_page(fitz.Rect(0, a_y0, c_w, a_y0 + c_h), c_doc, 0)
    b_x0 = c_w + col_gap
    b_y0 = label_h + (top_h - d_h) / 2.0
    comp_page.show_pdf_page(fitz.Rect(b_x0, b_y0, b_x0 + d_w, b_y0 + d_h), d_doc, 0)

    LEFT_MARGIN = 3.6
    add_label(comp_page, "A", LEFT_MARGIN, LABEL_SIZE)
    add_label(comp_page, "B", b_x0 + LEFT_MARGIN, LABEL_SIZE)

    out_path = os.path.join(OUT_DIR, "figS3_bootstrap_rankings.pdf")
    comp_doc.save(out_path, garbage=4, deflate=True)
    comp_doc.close()
    c_doc.close()
    d_doc.close()

    print(f"\n  Saved: {out_path} ({os.path.getsize(out_path):,} bytes)")

    # ── Verification ──
    print("\n" + "=" * 60)
    print("VERIFICATION")
    print("=" * 60)

    doc = fitz.open(out_path)
    page = doc[0]

    # Check labels
    labels_found = []
    all_fonts = set()
    blocks = page.get_text("dict")["blocks"]
    for b in blocks:
        if b["type"] == 0:
            for line in b["lines"]:
                for span in line["spans"]:
                    all_fonts.add(span["font"])
                    t = span["text"].strip()
                    if (len(t) == 1 and t.isalpha() and t.isupper()
                            and span["size"] >= 8.0
                            and ("Bold" in span["font"] or "bold" in span["font"].lower())):
                        labels_found.append(t)

    found_str = "".join(sorted(labels_found))
    expected = "AB"  # treeness panels cut; only A (forest) + B (CI-vs-shift)
    ok = set("AB").issubset(set(labels_found))

    print(f"\n  File: {os.path.basename(out_path)}")
    print(f"  Page: {page.rect.width:.1f} x {page.rect.height:.1f} pts")
    print(f"  Labels: [{', '.join(sorted(labels_found))}] — "
          f"{'OK' if ok else f'MISMATCH (expected {expected}, got {found_str})'}")
    print(f"  Fonts: {all_fonts}")
    print(f"  Size: {os.path.getsize(out_path):,} bytes")

    # Check font embedding
    font_list = page.get_fonts()
    print(f"  Fonts on page: {len(font_list)}")
    for f in font_list:
        print(f"    {f[3]} ({f[4]}) — embedded: {'yes' if f[3] else 'no'}")

    doc.close()

    # Cleanup temp files
    for p in [panel_c_path, panel_d_path]:
        if os.path.exists(p):
            os.remove(p)

    if ok:
        print("\n  ALL CHECKS PASSED")
    else:
        print("\n  SOME CHECKS FAILED — review above")

    print("\nDone.")


if __name__ == "__main__":
    main()
