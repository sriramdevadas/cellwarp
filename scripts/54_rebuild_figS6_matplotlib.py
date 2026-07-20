#!/usr/bin/env python3
"""Rebuild Figure S6 (CellMarker enrichment) from matplotlib source.

Two panels:
  Panel A — pooled enrichment (primary 4.49x, expression-matched 3.32x) with
            Clopper-Pearson 95% CIs as whiskers, K/N/n inline below each bar.
  Panel B — per-cell-type forest plot (6 cell types, top-50 method) with
            Clopper-Pearson 95% CIs. CD4+ T (k=0, FAIL) rendered as a
            distinct '×' marker without a CI whisker.

Data sources:
  output/validation/cellmarker_35type_rerun/cellmarker_35type_results.json
    .human_primary_35type            → Panel A bar 1 (k=34, n=500, K=257, N=16,959)
    .human_expression_matched_35type → Panel A bar 2 (k=34, n=500, K=56,  N=2,738)
    .per_cell_type_35type            → Panel B (6 cell types)

Replaces the PyMuPDF-overlay version from Part 4A which left the original
Unicode superscript glyphs underneath the ASCII e-notation. Uses matplotlib
mathtext for proper superscript rendering.
"""
from __future__ import annotations
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
from scipy.stats import beta

PROJECT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT / "src"))
from cellwarp.figure_style import (  # noqa
    C_BLUE, C_TEAL, C_DARKGRAY, C_GRAY, C_LIGHTGRAY, C_ORANGE,
    short_name,
)

OUT_PDFS = [
    PROJECT / "figures/supplementary/figS6_cellmarker_enrichment.pdf",
]

DATA_JSON = PROJECT / "output/validation/cellmarker_35type_rerun/cellmarker_35type_results.json"

for font in ["Arial", "Helvetica", "DejaVu Sans"]:
    try:
        fm.findfont(font, fallback_to_default=False); FONT = font; break
    except Exception:
        continue
else:
    FONT = "sans-serif"

plt.rcParams.update({
    "font.family": FONT, "font.size": 8,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})


def cp_fold_ci(k: int, n: int, K: int, N: int, alpha: float = 0.05):
    """Clopper-Pearson 95% CI on hypergeometric fold-enrichment.

    Returns (fold, ci_low, ci_high). Computes exact CP CI on the observed
    proportion k/n (Beta(k, n-k+1) / Beta(k+1, n-k) tails), then divides by
    the expected proportion K/N. For k=0 the lower bound is 0 by convention.
    """
    p_low = 0.0 if k == 0 else float(beta.ppf(alpha / 2, k, n - k + 1))
    p_high = 1.0 if k == n else float(beta.ppf(1 - alpha / 2, k + 1, n - k))
    exp = K / N
    fold = (k / n) / exp if exp > 0 else 0.0
    return fold, p_low / exp, p_high / exp


def main():
    data = json.loads(DATA_JSON.read_text())
    primary = data["human_primary_35type"]
    exp_matched = data["human_expression_matched_35type"]
    per_type = data["per_cell_type_35type"]

    # ── Panel A: pooled enrichment with Clopper-Pearson CIs ────────
    panelA = [
        {
            "label": "Primary\nenrichment",
            # Two-line K/n/N annotation: narrower footprint so the two bars'
            # triplets don't overlap at the bottom of the axes.
            "sub": f"K = {primary['n_cellmarker_genes']:,}; n = {primary['n_identity_genes']}\n"
                   f"N = {primary['n_background_genes']:,}",
            "k": primary["observed_overlap"], "n": primary["n_identity_genes"],
            "K": primary["n_cellmarker_genes"], "N": primary["n_background_genes"],
            "color": C_BLUE,
            "p_mathtext": r"$p = 2.10 \times 10^{-13}$",
        },
        {
            "label": "Expression-\nmatched",
            "sub": f"K = {exp_matched['n_cellmarker_in_universe']}; n = {exp_matched['n_identity_genes']}\n"
                   f"N = {exp_matched['universe_size']:,}",
            "k": exp_matched["observed_overlap"], "n": exp_matched["n_identity_genes"],
            "K": exp_matched["n_cellmarker_in_universe"], "N": exp_matched["universe_size"],
            "color": C_TEAL,
            "p_mathtext": r"$p = 1.15 \times 10^{-12}$",
        },
    ]
    for b in panelA:
        f, lo, hi = cp_fold_ci(b["k"], b["n"], b["K"], b["N"])
        b["fold"] = f; b["ci_lo"] = lo; b["ci_hi"] = hi

    # ── Panel B: per-cell-type forest with Clopper-Pearson CIs ─────
    # Per-type analysis uses primary background (N=16,959) per
    # cellmarker_35type_rerun.py:180 (hypergeom_enrichment uses n_genes).
    N_PER_TYPE = primary["n_background_genes"]
    for r in per_type:
        f, lo, hi = cp_fold_ci(
            r["overlap"], r["n_loading_genes"],
            r["n_cellmarker_markers"], N_PER_TYPE,
        )
        r["fold"] = f; r["ci_lo"] = lo; r["ci_hi"] = hi

    # Sort PASS types by fold descending; FAIL pinned to bottom.
    pass_types = sorted([r for r in per_type if r["pass"] == "PASS"],
                        key=lambda r: r["fold"], reverse=True)
    fail_types = [r for r in per_type if r["pass"] != "PASS"]
    ordered = pass_types + fail_types  # 5 PASS at top, 1 FAIL at bottom

    # ── Layout: A left, B right ───────────────────────────────────
    fig = plt.figure(figsize=(7.0, 3.2))
    gs = fig.add_gridspec(
        1, 2, left=0.07, right=0.98, top=0.88, bottom=0.20,
        wspace=0.45, width_ratios=[1.0, 1.4],
    )
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])

    # ── Render Panel A ─────────────────────────────────────────────
    x = list(range(len(panelA)))
    folds = [b["fold"] for b in panelA]
    err_low = [b["fold"] - b["ci_lo"] for b in panelA]
    err_high = [b["ci_hi"] - b["fold"] for b in panelA]
    colors = [b["color"] for b in panelA]
    axA.bar(x, folds, yerr=[err_low, err_high], color=colors,
            edgecolor=C_DARKGRAY, linewidth=0.5, capsize=4, width=0.55,
            error_kw={"elinewidth": 0.6, "ecolor": C_DARKGRAY})

    # p-value above each bar (above CI upper)
    for i, b in enumerate(panelA):
        y = b["ci_hi"] + 0.20
        axA.text(i, y, b["p_mathtext"], ha="center", va="bottom",
                 fontsize=7, color=b["color"], fontweight="bold")

    # Reference line + label
    axA.axhline(1.0, color=C_GRAY, linestyle="--", linewidth=0.8, zorder=0)
    axA.text(1.55, 1.0, "no enrichment\n(fold = 1.0)",
             ha="left", va="center", fontsize=6.5,
             color=C_GRAY, style="italic")

    axA.set_xticks(x)
    axA.set_xticklabels([b["label"] for b in panelA], fontsize=7.5)
    # Sample-size annotations directly below x-tick labels (figure coord
    # via blended transform: x in data, y in axes — slightly below 0).
    for i, b in enumerate(panelA):
        axA.annotate(
            b["sub"], xy=(i, 0), xytext=(0, -28),
            textcoords="offset points", ha="center", va="top",
            fontsize=5.5, color=C_DARKGRAY,
        )

    axA.set_xlim(-0.6, 2.3)
    axA.set_ylabel("Fold enrichment vs background", fontsize=8)
    # Drop "Test" xlabel (bars are self-labeled; sample-size sub-line below).
    axA.set_ylim(0, max(b["ci_hi"] for b in panelA) + 1.0)
    axA.spines["top"].set_visible(False)
    axA.spines["right"].set_visible(False)
    axA.set_title("Pooled enrichment", fontsize=8, loc="left", pad=2.0)

    # ── Render Panel B: horizontal forest plot ─────────────────────
    n_rows = len(ordered)
    y_positions = list(range(n_rows))[::-1]  # top row = idx 0

    # Fixed x-axis 0–500 fully accommodates all CIs (Hepatocyte ci_hi=465,
    # CD8+ T ci_hi=451) — previously dynamic max*1.05 clipped or compressed
    # the upper whiskers. 500 is a clean tick boundary at this scale.
    x_max = 500.0
    for yp, r in zip(y_positions, ordered):
        if r["pass"] == "PASS":
            axB.plot([r["ci_lo"], r["ci_hi"]], [yp, yp],
                     color=C_TEAL, linewidth=1.4, solid_capstyle="round",
                     zorder=2)
            axB.plot(r["fold"], yp, "o", color=C_TEAL, markersize=4.5,
                     markeredgecolor="white", markeredgewidth=0.4, zorder=3)
        else:
            # CD4+ T (k=0): distinct '×' marker at fold=0, no CI whisker.
            axB.plot(0, yp, marker="x", color=C_ORANGE, markersize=6,
                     markeredgewidth=1.2, zorder=3)
            axB.annotate(
                "overlap = 0", xy=(0, yp), xytext=(6, 0),
                textcoords="offset points", ha="left", va="center",
                fontsize=6.5, color=C_ORANGE, style="italic",
            )

    # Reference line at fold=1.0
    axB.axvline(1.0, color=C_GRAY, linestyle="--", linewidth=0.6, zorder=0)

    # y-tick labels: short cell-type names
    axB.set_yticks(y_positions)
    axB.set_yticklabels([short_name(r["cell_type"]) for r in ordered],
                        fontsize=7)
    axB.set_xlabel("Fold enrichment vs background", fontsize=8)
    axB.set_xlim(-15, x_max)
    axB.set_ylim(-0.7, n_rows - 0.3)
    axB.spines["top"].set_visible(False)
    axB.spines["right"].set_visible(False)
    axB.set_title("Per-cell-type enrichment (n=50 identity genes each)",
                  fontsize=8, loc="left", pad=2.0)

    # ── Panel letters ──────────────────────────────────────────────
    for ax, letter in [(axA, "A"), (axB, "B")]:
        ax.text(-0.18, 1.06, letter, transform=ax.transAxes,
                fontsize=10, fontweight="bold", va="bottom", ha="left")

    for out in OUT_PDFS:
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(f"  wrote: {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
