#!/usr/bin/env python3
"""
FW-024 correction: regenerate Figure S4 panel A from matched 15-type
Procrustes baseline.

Panel A: rank scatter (primary matched-PCA rank vs. CellHint rank), n = 15,
annotated with rho = -0.139, p = 0.621, one label per point.

Panel B: existing factor dot plot (systematic_factors.csv), unchanged logic.

Reads per-type residuals from the newly-persisted Level 0 CSVs:
    analysis/harmonized_replication/harmonized_residuals_primary_level0.csv
    analysis/harmonized_replication/harmonized_residuals_cellhint_level0.csv

Outputs (overwritten):
    figures/panels/figs2a_cellhint_rank.{pdf,png}
    figures/panels/figs2b_cellhint_residual.{pdf,png}
    figures/supplementary/figS2_cellhint_investigation.{pdf,png}
    figures/supplementary/figS5_cellhint_investigation_polished.pdf
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

BASE = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(BASE / "src"))

from cellwarp.figure_style import (  # noqa: E402
    COL1,
    COL2,
    C_BLUE,
    C_DARKGRAY,
    C_GRAY,
    C_LIGHTGRAY,
    C_TEAL,
    FONT_FAMILY,
    FONT_SIZE_ANNOT,
    FONT_SIZE_LABEL,
    FONT_SIZE_LEGEND,
    FONT_SIZE_PANEL,
    FONT_SIZE_TICK,
    add_panel_label,
    apply_style,
    clean_spine,
    lineage_color,
    save_figure,
    short_name,
)

apply_style()

PANELS = BASE / "figures" / "panels"
SUPP = BASE / "figures" / "supplementary"
SUBMISSION_SUPP = BASE / "figures" / "submission" / "supplementary"
REVIEW_DIR = BASE / "docs" / "submission" / "figures_for_review"
HARMON_DIR = BASE / "analysis" / "harmonized_replication"
FACTORS_CSV = BASE / "analysis" / "cellhint_investigation" / "systematic_factors.csv"


def _render_panel_a() -> None:
    primary = pd.read_csv(HARMON_DIR / "harmonized_residuals_primary_level0.csv")
    cellhint = pd.read_csv(HARMON_DIR / "harmonized_residuals_cellhint_level0.csv")
    merged = primary.merge(cellhint, on="cell_type", suffixes=("_primary", "_cellhint"))
    if len(merged) != 15:
        raise SystemExit(f"Expected n=15 after merge, got {len(merged)}")

    rho, p = spearmanr(merged["rank_primary"], merged["rank_cellhint"])
    if abs(rho - (-0.139)) > 1e-3 or abs(p - 0.621) > 1e-2:
        raise SystemExit(
            f"Level 0 correlation drift: rho={rho:.6f}, p={p:.6f} "
            "(expected rho=-0.139, p=0.621)"
        )

    # Per-point lineage colors (matches Fig 3A / S2C/D / S3A-B convention).
    # All 15 matched cell_type strings resolve through LINEAGE_MAP — no
    # fallback to default; verified in audit.
    point_colors = [lineage_color(ct) for ct in merged["cell_type"]]

    fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.85))
    ax.scatter(
        merged["rank_primary"],
        merged["rank_cellhint"],
        c=point_colors,
        s=30,
        edgecolors="white",
        linewidths=0.3,
        zorder=3,
    )

    max_r = int(max(merged["rank_primary"].max(), merged["rank_cellhint"].max())) + 1
    ax.plot([0, max_r], [0, max_r], color=C_LIGHTGRAY, linewidth=0.8, linestyle="--", zorder=1)

    for _, row in merged.iterrows():
        ax.annotate(
            short_name(row["cell_type"]),
            xy=(row["rank_primary"], row["rank_cellhint"]),
            xytext=(3, 3),
            textcoords="offset points",
            fontsize=FONT_SIZE_ANNOT - 1,
            color=C_DARKGRAY,
            zorder=4,
        )

    # Axis direction annotation matches S3C ("1 = most diverged, 35 = most
    # conserved"); scaled here to the 15-type matched set.
    ax.set_xlabel(
        "Primary rank (matched 15-type Procrustes)\n"
        "(1 = most diverged, 15 = most conserved)",
        fontsize=FONT_SIZE_LABEL,
    )
    ax.set_ylabel(
        "CellHint rank\n(1 = most diverged, 15 = most conserved)",
        fontsize=FONT_SIZE_LABEL,
    )
    # Annotation in top-right (data points cluster toward the diagonal; top-right
    # at (rank 14-16, 14-16) is free — "NK cell" at (6, 15) overlapped a
    # centered top annotation).
    ax.text(
        0.97,
        0.97,
        f"ρ = {rho:.3f}\np = {p:.3f}\nn = 15",
        transform=ax.transAxes,
        fontsize=FONT_SIZE_ANNOT,
        ha="right",
        va="top",
        color=C_DARKGRAY,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=C_LIGHTGRAY, alpha=0.9),
    )
    clean_spine(ax)
    save_figure(fig, PANELS / "figs2a_cellhint_rank")
    plt.close(fig)


def _render_panel_b() -> None:
    factors_df = pd.read_csv(FACTORS_CSV)

    # Drop "Absolute residual difference" — tautological per the source CSV's
    # own "Sanity check: do large rank diffs correspond to large residual
    # diffs?" interpretation. |residual diff| and |rank diff| are monotonically
    # coupled by rank-transform construction; testing one against the other is
    # circular and not a substantive driver.
    before_n = len(factors_df)
    factors_df = factors_df[
        factors_df["factor"] != "Absolute residual difference"
    ].reset_index(drop=True)
    after_n = len(factors_df)
    if before_n - after_n != 1:
        raise SystemExit(
            f"Expected to drop exactly 1 'Absolute residual difference' row, "
            f"dropped {before_n - after_n} (n {before_n} → {after_n})"
        )

    fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.65))
    y_pos = np.arange(len(factors_df))[::-1]
    for i, (_, row) in enumerate(factors_df.iterrows()):
        y = y_pos[i]
        rho_val = float(row["rho"])
        p_val = float(row["p_value"])
        color = C_TEAL if p_val < 0.05 else C_GRAY
        marker = "o" if p_val < 0.05 else "s"
        ax.plot(
            rho_val,
            y,
            marker,
            color=color,
            markersize=6,
            zorder=4,
            markeredgecolor="white",
            markeredgewidth=0.5,
        )

    ax.axvline(0, color=C_DARKGRAY, linewidth=0.5, linestyle="-", zorder=1)
    ax.set_yticks(y_pos)
    factor_labels = []
    for f in factors_df["factor"]:
        f = str(f).replace("abs(log2 cell count ratio)", "Cell count asymmetry")
        f = f.replace("log2(CellHint/primary cell count)", "Count direction")
        f = f.replace("CellHint tissue count", "Tissue count")
        if len(f) > 30:
            f = f[:27] + "..."
        factor_labels.append(f)
    ax.set_yticklabels(factor_labels, fontsize=FONT_SIZE_TICK)
    ax.set_xlabel("Spearman ρ", fontsize=FONT_SIZE_LABEL)
    ax.set_title("Factors associated with rank reversal", fontsize=FONT_SIZE_LABEL)

    legend_el = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=C_TEAL, markersize=5, label="p < 0.05"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor=C_GRAY, markersize=5, label="NS"),
    ]
    ax.legend(handles=legend_el, fontsize=FONT_SIZE_LEGEND, frameon=False, loc="lower right")

    clean_spine(ax)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    save_figure(fig, PANELS / "figs2b_cellhint_residual")
    plt.close(fig)


def _assemble_composite() -> Path:
    def load_panel(name: str):
        path = PANELS / f"{name}.png"
        return mpimg.imread(str(path))

    fig = plt.figure(figsize=(COL2, COL2 * 0.45))
    gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.3)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.imshow(load_panel("figs2a_cellhint_rank"))
    ax_a.axis("off")
    add_panel_label(ax_a, "A", x=-0.02, y=1.02)
    ax_b = fig.add_subplot(gs[0, 1])
    ax_b.imshow(load_panel("figs2b_cellhint_residual"))
    ax_b.axis("off")
    add_panel_label(ax_b, "B", x=-0.02, y=1.02)

    composite_stem = SUPP / "figS2_cellhint_investigation"
    save_figure(fig, composite_stem, tight=False)
    plt.close(fig)
    return composite_stem.with_suffix(".pdf")


def main() -> int:
    PANELS.mkdir(parents=True, exist_ok=True)
    SUPP.mkdir(parents=True, exist_ok=True)
    SUBMISSION_SUPP.mkdir(parents=True, exist_ok=True)
    REVIEW_DIR.mkdir(parents=True, exist_ok=True)

    _render_panel_a()
    _render_panel_b()
    composite_pdf = _assemble_composite()

    # Out-of-bundle producer-tree copy; the cellhint figure was cut from the
    # submission, so no packet mirror is materialized from it.
    polished = SUPP / "figS5_cellhint_investigation_polished.pdf"
    shutil.copy2(composite_pdf, polished)

    print(f"Wrote {composite_pdf.relative_to(BASE)}")
    print(f"Wrote {polished.relative_to(BASE)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
