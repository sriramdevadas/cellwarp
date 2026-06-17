#!/usr/bin/env python3
"""
CellWarp — Standalone figure script for expanded negative controls.

Reads pre-computed results from CSVs and generates publication figure.
Run after expanded_negative_controls.py has completed.

Usage:
    python negative_control_figure.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Paths
OUTPUT_DIR = PROJECT_ROOT / "analysis" / "expanded_negative_controls"
FIGURE_DIR = PROJECT_ROOT / "figures" / "supplementary"
CROSS_SPECIES_RESULTS = PROJECT_ROOT / "output" / "phase2" / "scaled_35types" / "procrustes_results_35.json"
CROSS_SPECIES_NULL = PROJECT_ROOT / "output" / "phase2" / "scaled_35types" / "null_distribution_35.npy"


def main():
    # Load results
    within_df = pd.read_csv(OUTPUT_DIR / "within_species_pairs.csv")
    self_df = pd.read_csv(OUTPUT_DIR / "self_comparison_results.csv")

    with open(CROSS_SPECIES_RESULTS) as f:
        xsp = json.load(f)
    xsp_distance = xsp["procrustes"]["distance"]
    xsp_null_median = xsp["permutation_test"]["null_distribution_summary"]["median"]
    xsp_ratio = xsp_distance / xsp_null_median
    xsp_null = np.load(CROSS_SPECIES_NULL)

    # Figure
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    rng = np.random.RandomState(42)

    colors = {
        "Permutation\nnull": "#bdbdbd",
        "Within-species\ntissue pairs": "#fc8d62",
        "Cross-species\n(primary)": "#e41a1c",
        "Self-comparison\n(random split)": "#66c2a5",
    }

    cat_order = [
        "Permutation\nnull",
        "Within-species\ntissue pairs",
        "Cross-species\n(primary)",
        "Self-comparison\n(random split)",
    ]

    # Data
    null_median = float(np.median(xsp_null))
    null_ratios = xsp_null / null_median
    null_sample = rng.choice(null_ratios, size=500, replace=False)

    datasets = {
        "Permutation\nnull": null_sample,
        "Within-species\ntissue pairs": within_df["obs_to_null_ratio"].values,
        "Cross-species\n(primary)": np.array([xsp_ratio]),
        "Self-comparison\n(random split)": self_df["obs_to_null_ratio"].values,
    }

    for i, cat in enumerate(cat_order):
        data = datasets[cat]

        if len(data) > 5:
            parts = ax.violinplot(
                data, positions=[i], showmedians=True,
                showextrema=False, widths=0.6,
            )
            for pc in parts["bodies"]:
                pc.set_facecolor(colors[cat])
                pc.set_alpha(0.6)
            parts["cmedians"].set_color("black")
            parts["cmedians"].set_linewidth(1.5)

        # Wider jitter (±0.12 → ±0.25) spreads bunched dots across the violin
        # width — particularly important for Self-comparison where values
        # cluster tightly near y ≈ 0.03 and were visually stacking at one x.
        jitter = rng.uniform(-0.25, 0.25, size=len(data))

        if cat == "Permutation\nnull":
            ax.scatter([i], [np.median(data)], color="black", s=30, zorder=5, marker="D")
        elif cat == "Cross-species\n(primary)":
            ax.scatter([i], data, color=colors[cat], s=200, zorder=5, marker="*",
                       edgecolors="black", linewidths=0.8)
        else:
            ax.scatter(i + jitter, data, color=colors[cat], s=50, alpha=0.75,
                       zorder=4, edgecolors="white", linewidths=0.4)

    ax.axhline(y=1.0, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.text(3.55, 1.0, "random", fontsize=8, color="gray", va="center")
    ax.axhline(y=xsp_ratio, color="#e41a1c", linestyle=":", linewidth=0.8, alpha=0.5)

    ax.set_xticks(range(len(cat_order)))
    ax.set_xticklabels(cat_order, fontsize=8, rotation=30, ha='right')
    ax.set_ylabel("Procrustes distance / null median\n(lower = more coherent)", fontsize=11)
    ax.set_title("Within-species negative control: coherence hierarchy", fontsize=12, fontweight="bold")

    ws_ratios = within_df["obs_to_null_ratio"].values
    n_as_strong = np.sum(ws_ratios <= xsp_ratio)
    frac = n_as_strong / len(ws_ratios) * 100
    ax.text(
        0.98, 0.97,
        f"Within-species: {n_as_strong}/{len(ws_ratios)} pairs ({frac:.0f}%)\n"
        f"more coherent than cross-species\n"
        f"(see Methods,\n"
        f" Expanded negative controls)",
        transform=ax.transAxes, fontsize=8, va="top", ha="right",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8),
    )

    # Place n-counts as superscripts on tick labels to avoid overlap with rotated labels
    current_labels = [t.get_text() for t in ax.get_xticklabels()]
    new_labels = []
    for i, cat in enumerate(cat_order):
        n = len(datasets[cat])
        new_labels.append(f"{cat}\nn={n}")
    ax.set_xticklabels(new_labels, fontsize=8, rotation=30, ha='right')

    ax.set_xlim(-0.5, len(cat_order) - 0.5)
    sns.despine(ax=ax)
    plt.tight_layout()

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_DIR / "negative_control_distributions.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURE_DIR / "negative_control_distributions.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved to {FIGURE_DIR / 'negative_control_distributions.pdf'}")


if __name__ == "__main__":
    main()
