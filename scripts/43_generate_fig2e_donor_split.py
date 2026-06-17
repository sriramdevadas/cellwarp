#!/usr/bin/env python3
"""
Generate publication-quality Figure 3E panel: donor-split coherence hierarchy.

Reads shared-PCA donor-split results and creates a box plot showing:
  - Self-comparison floor (0.033, green dot)
  - Donor-split within-species distribution (blue box)
  - Cross-species matched distribution (orange box)
  - Permutation null reference line (gray dotted at 1.0)
  - Delta annotation

Style matches existing Figure 3 panels (Arial, 8pt labels, clean axes).

Outputs:
  figures/panels/fig2e_donor_split.png (300 dpi)
  figures/panels/fig2e_donor_split.pdf
"""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import pandas as pd

# Paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent
RESULTS_PATH = PROJECT_ROOT / "analysis" / "donor_split" / "donor_split_shared_pca_distributions.csv"
OUT_DIR = PROJECT_ROOT / "figures" / "panels"

# Reference values
SELF_COMPARISON = 0.033

# Try to use Arial/Helvetica
for font_name in ["Arial", "Helvetica", "DejaVu Sans"]:
    try:
        fm.findfont(font_name, fallback_to_default=False)
        FONT = font_name
        break
    except Exception:
        continue
else:
    FONT = "sans-serif"


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load data
    df = pd.read_csv(RESULTS_PATH)
    ws_ratios = df["ws_obs_null_ratio"].values
    cs_ratios = df["cs_obs_null_ratio"].values

    ws_median = np.median(ws_ratios)
    cs_median = np.median(cs_ratios)
    delta_median = np.median(df["delta_obs_null"].values)
    delta_lo = np.percentile(df["delta_obs_null"].values, 2.5)
    delta_hi = np.percentile(df["delta_obs_null"].values, 97.5)
    n_positive = (df["delta_obs_null"].values > 0).sum()

    # Create figure
    plt.rcParams.update({
        "font.family": FONT,
        "font.size": 8,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
    })

    # Full-width source panel for the new 3-row Fig 3 composite layout
    # (X-fig3-restructure). Panel E now occupies a full row, so the
    # three multi-line xtick labels have ample horizontal room.
    fig, ax = plt.subplots(figsize=(7.0, 3.8))

    # Colors matching existing Figure 3 palette
    green = "#2ca02c"
    blue = "#1f77b4"
    orange = "#ff7f0e"

    # Self-comparison point
    ax.scatter([1], [SELF_COMPARISON], s=80, color=green, zorder=5,
               edgecolors="black", linewidth=0.6)

    # Donor-split within-species box
    bp1 = ax.boxplot([ws_ratios], positions=[2], widths=0.45,
                     patch_artist=True, showfliers=False,
                     medianprops=dict(color="black", linewidth=1.2),
                     whiskerprops=dict(linewidth=0.8),
                     capprops=dict(linewidth=0.8))
    bp1["boxes"][0].set_facecolor(blue)
    bp1["boxes"][0].set_alpha(0.7)
    bp1["boxes"][0].set_edgecolor("black")
    bp1["boxes"][0].set_linewidth(0.8)

    # Cross-species matched box
    bp2 = ax.boxplot([cs_ratios], positions=[3], widths=0.45,
                     patch_artist=True, showfliers=False,
                     medianprops=dict(color="black", linewidth=1.2),
                     whiskerprops=dict(linewidth=0.8),
                     capprops=dict(linewidth=0.8))
    bp2["boxes"][0].set_facecolor(orange)
    bp2["boxes"][0].set_alpha(0.7)
    bp2["boxes"][0].set_edgecolor("black")
    bp2["boxes"][0].set_linewidth(0.8)

    # Permutation null reference line
    ax.axhline(1.0, color="gray", linestyle=":", linewidth=0.8, zorder=1)
    ax.text(3.55, 1.0, "Permutation\nnull", fontsize=6.5, color="gray",
            va="center", ha="left")

    # Delta annotation with bracket
    bracket_x = 3.55
    ax.annotate("", xy=(bracket_x, ws_median), xytext=(bracket_x, cs_median),
                arrowprops=dict(arrowstyle="<->", color="black", linewidth=0.8))
    ax.text(bracket_x + 0.12, (ws_median + cs_median) / 2,
            f"\u0394 = +{delta_median:.3f}\n100/100 positive",
            fontsize=6.5, va="center", ha="left",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                      edgecolor="gray", alpha=0.8, linewidth=0.5))

    # Labels
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels([
        f"Single TS\nself-comparison\n({SELF_COMPARISON:.3f})",
        f"Donor-split\nwithin-species\n(median {ws_median:.3f})",
        f"Cross-species\nmatched\n(median {cs_median:.3f})",
    ], fontsize=7, ha='center')
    ax.set_ylabel("obs/null ratio", fontsize=8)
    ax.set_ylim(-0.03, 1.12)
    ax.set_xlim(0.4, 4.3)

    # Clean up
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()

    # Save
    for ext in ["png", "pdf"]:
        path = OUT_DIR / f"fig2e_donor_split.{ext}"
        fig.savefig(path, dpi=300, bbox_inches="tight",
                    facecolor="white", edgecolor="none")
        print(f"Saved: {path}")

    plt.close(fig)

    # Text description
    print(f"\nFigure 3E: Donor-split coherence hierarchy.")
    print(f"  Self-comparison: {SELF_COMPARISON}")
    print(f"  Within-species (donor-split): median {ws_median:.3f}")
    print(f"  Cross-species (matched): median {cs_median:.3f}")
    print(f"  Delta: +{delta_median:.3f} [{delta_lo:+.3f}, {delta_hi:+.3f}]")
    print(f"  Positive: {n_positive}/100")


if __name__ == "__main__":
    main()
