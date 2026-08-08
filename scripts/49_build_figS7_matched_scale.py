#!/usr/bin/env python3
"""Figure S7 rebuild — 6-type matched-scale human-vs-human vs human-vs-mouse
negative control, under confounded (cross-atlas) protocol.

After Part 4A:
  y-axis 0 to 1.1
  dashed reference line at obs/null = 1.0 ('permutation null')
  violin of each null distribution rescaled to obs/null (centered at 1.0)
  observed points marked for both analyses
  p-values: H×M 0.0035, H×H 0.0088; 1.92-fold ratio
  protocol-confounded annotation pointing to S1 Text for the protocol-
    controlled donor-split result.

Source nulls:
  H×M 6-type: output/phase2/null_distribution.npy (d=12.427, null_median=39.258)
  H×H 6-type v2: output/phase2/negative_control_v2/null_distribution_negctrl_v2.npy
                 (d=24.009, null_median=39.576)
"""
from __future__ import annotations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np

PROJECT = Path(__file__).resolve().parent.parent
NULL_HM = PROJECT / "output/phase2/null_distribution.npy"
NULL_HH = PROJECT / "output/phase2/negative_control_v2/null_distribution_negctrl_v2.npy"
OUT_SUPP = PROJECT / "figures/submission/supplementary"
OUT_LEGACY = PROJECT / "figures/supplementary"

# Observed values from source JSONs
HM_DIST = 12.427310325548715
HM_NULL_MED = 39.25758219423939
HH_DIST = 24.008688408425147
HH_NULL_MED = 39.576332731611004
HM_P = 0.0035
HH_P = 0.0088

for font in ["Arial", "Helvetica", "DejaVu Sans"]:
    try:
        fm.findfont(font, fallback_to_default=False)
        FONT = font
        break
    except Exception:
        continue
else:
    FONT = "sans-serif"
plt.rcParams.update({
    "font.family": FONT,
    "font.size": 8,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

C_BLUE = "#1f77b4"
C_ORANGE = "#ff7f0e"
C_GRAY = "#7f7f7f"
C_LIGHTGRAY = "#cccccc"
C_RED = "#d62728"


def main():
    OUT_SUPP.mkdir(parents=True, exist_ok=True)
    OUT_LEGACY.mkdir(parents=True, exist_ok=True)

    hm_null = np.load(NULL_HM)
    hh_null = np.load(NULL_HH)
    # Express each null distribution as obs/null ratio (null / own null_median)
    hm_norm = hm_null / HM_NULL_MED
    hh_norm = hh_null / HH_NULL_MED
    hm_obs = HM_DIST / HM_NULL_MED
    hh_obs = HH_DIST / HH_NULL_MED
    fold = hh_obs / hm_obs

    fig, ax = plt.subplots(figsize=(4.2, 3.7))

    # Violins for null distributions
    parts = ax.violinplot(
        [hm_norm, hh_norm], positions=[0, 1], widths=0.6, showmeans=False,
        showextrema=False, showmedians=False,
    )
    colors = [C_BLUE, C_ORANGE]
    for pc, c in zip(parts["bodies"], colors):
        pc.set_facecolor(c)
        pc.set_edgecolor(c)
        pc.set_alpha(0.25)

    # Observed values
    ax.scatter([0], [hm_obs], marker="D", s=50, color=C_BLUE,
               edgecolor="white", linewidth=0.8, zorder=10,
               label=f"H×M observed: obs/null = {hm_obs:.3f}")
    ax.scatter([1], [hh_obs], marker="D", s=50, color=C_ORANGE,
               edgecolor="white", linewidth=0.8, zorder=10,
               label=f"H×H observed: obs/null = {hh_obs:.3f}")

    # Permutation null reference
    ax.axhline(1.0, color=C_GRAY, linestyle="--", linewidth=0.8, zorder=1)
    ax.text(1.35, 1.005, "permutation null (obs/null = 1.0)",
            color=C_GRAY, fontsize=7, ha="right", va="bottom")

    # Observed annotations — both placed to the RIGHT of their markers,
    # mirroring layouts so the H×M label clears the y-axis tick column.
    ax.annotate(f"obs/null = {hm_obs:.3f}\n$p = {HM_P}$",
                xy=(0, hm_obs), xytext=(0 + 0.35, hm_obs - 0.04),
                ha="left", va="top", fontsize=7.5, color=C_BLUE, fontweight="bold",
                arrowprops=dict(arrowstyle="-", color=C_BLUE, lw=0.5))
    ax.annotate(f"obs/null = {hh_obs:.3f}\n$p = {HH_P}$",
                xy=(1, hh_obs), xytext=(1 + 0.35, hh_obs + 0.04),
                ha="left", va="bottom", fontsize=7.5, color=C_ORANGE, fontweight="bold",
                arrowprops=dict(arrowstyle="-", color=C_ORANGE, lw=0.5))

    # Fold-change annotation
    ax.annotate(f"{fold:.2f}-fold",
                xy=(0.5, (hm_obs + hh_obs) / 2), ha="center", va="center",
                fontsize=8, color="black", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="black", lw=0.5))
    # Connecting line
    ax.plot([0.05, 0.95], [hm_obs, hh_obs], color="black", linewidth=0.6, zorder=3)

    # Axes
    ax.set_xticks([0, 1])
    ax.set_xticklabels(
        ["Human × Mouse\n(6-type matched scale)",
         "Human × Human\n(6-type within-species)"],
        fontsize=8,
    )
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(0, 1.1)
    ax.set_ylabel("obs/null ratio", fontsize=8.5)
    # In-panel title removed (Cell Systems convention: panels carry no
    # in-figure titles; identification is via figure legend).
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    # Banner annotation
    ax.text(0.5, 0.03,
            "Protocol-confounded (different atlases);\n"
            "see S1 Text for the protocol-controlled donor-split result",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=7,
            color=C_GRAY, style="italic")

    fig.tight_layout()
    # Canonical write to OUT_SUPP only; OUT_LEGACY mirrors (.pdf + .png) are
    # materialized by scripts/build_submission_packet.py (R21 build script).
    base = OUT_SUPP / "figS4_matched_scale_control"
    base.parent.mkdir(parents=True, exist_ok=True)
    # PDF only: suppress the embedded /CreationDate so that re-running the
    # unchanged script in one environment reproduces the PDF byte-for-byte. The
    # PNG writer emits no timestamp (Software + dpi only), so it needs no
    # equivalent. Across environments the bytes still move with whichever
    # FreeType matplotlib links -- see reproduce/figure_script_map.md.
    for ext, meta in (("pdf", {"CreationDate": None}), ("png", None)):
        fig.savefig(f"{base}.{ext}", dpi=300, bbox_inches="tight", metadata=meta)
    print(f"  wrote: {base.name}.{{pdf,png}}")
    plt.close(fig)


if __name__ == "__main__":
    main()
