#!/usr/bin/env python3
"""
Standalone figure script for the CellWarp simulation study.

Reads simulation_results.json and generates the 4-panel publication figure
(Fig S7). Can be re-run independently to adjust styling without re-running
the full simulation.

Usage:
    python simulation_figures.py                     # default paths
    python simulation_figures.py results.json        # custom input
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ── Paths ──
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[1]
FIG_DIR = PROJECT_ROOT / "figures" / "supplementary"

# ── Constants (must match simulation_study.py) ──
REAL_OBS_NULL = 0.522
POWER_N_TYPES = [15, 25, 35]
RECOVERY_N_CELLS = [50, 200, 500, 2000]


def load_results(path=None):
    if path is None:
        path = SCRIPT_DIR / "simulation_results.json"
    with open(path) as f:
        return json.load(f)


def make_figure(results, output_dir=None):
    """Generate the 4-panel Fig S7."""
    if output_dir is None:
        output_dir = FIG_DIR

    power = results["power_curve"]
    recovery = results["ranking_recovery"]
    stability = results["stability"]
    null_cal = results["null_calibration"]
    real_signal = results["calibration"]["estimated_real_signal"]

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 9,
        "axes.labelsize": 10,
        "axes.titlesize": 10,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7.5,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.5))

    type_colors = {15: "#1f77b4", 25: "#ff7f0e", 35: "#2ca02c"}
    cell_colors = {
        25: "#bcbd22", 50: "#d62728", 100: "#8c564b", 200: "#1f77b4",
        500: "#2ca02c", 1000: "#e377c2", 2000: "#9467bd", 5000: "#17becf",
    }

    # Panel A: Power
    ax = axes[0, 0]
    for nt in POWER_N_TYPES:
        sub = [r for r in power if r["n_types"] == nt]
        sigs = [r["signal_strength"] for r in sub]
        det = [r["detection_rate"] for r in sub]
        ax.plot(sigs, det, "o-", color=type_colors[nt], markersize=4,
                linewidth=1.5, label=f"n = {nt} types")

    ax.axvline(real_signal, color="0.45", ls="--", lw=1, alpha=0.7)
    ann_x_offset = 0.3 if real_signal < 5 else -0.3
    ann_ha = "left" if real_signal < 5 else "right"
    ax.annotate(f"real data\nobs/null = {REAL_OBS_NULL}",
                xy=(real_signal, 0.50), fontsize=7, color="0.35",
                ha=ann_ha, va="center")
    ax.axhline(0.05, color="red", ls=":", lw=0.8, alpha=0.5, label="α = 0.05")
    ax.set_xlabel("Signal strength")
    ax.set_ylabel("Detection rate (p < 0.05)")
    ax.set_ylim(-0.05, 1.08)
    ax.legend(loc="center right", frameon=False)
    ax.set_title("A.  Detection power", fontweight="bold", loc="left")

    # Panel B: Ranking recovery
    ax = axes[0, 1]
    for nc in RECOVERY_N_CELLS:
        sub = [r for r in recovery if r["n_cells_per_type"] == nc]
        if not sub:
            continue
        sigs = [r["signal_strength"] for r in sub]
        rhos = [r["mean_rho"] for r in sub]
        stds = [r["std_rho"] for r in sub]
        ax.errorbar(sigs, rhos, yerr=stds, fmt="o-", color=cell_colors[nc],
                    markersize=4, linewidth=1.5, capsize=2,
                    label=f"{nc} cells/type")

    ax.axvline(real_signal, color="0.45", ls="--", lw=1, alpha=0.7)
    ax.set_xlabel("Signal strength")
    ax.set_ylabel("Spearman ρ (planted vs recovered)")
    ax.set_ylim(-0.15, 1.08)
    ax.legend(loc="lower right", frameon=False)
    ax.set_title("B.  Ranking recovery", fontweight="bold", loc="left")

    # Panel C: Ranking stability
    ax = axes[1, 0]
    nc_vals = [r["n_cells_per_type"] for r in stability]
    m_rhos = [r["mean_rho"] for r in stability]
    ci_lo = [r["ci_lower"] for r in stability]
    ci_hi = [r["ci_upper"] for r in stability]

    ax.fill_between(nc_vals, ci_lo, ci_hi, alpha=0.20, color="#1f77b4")
    ax.plot(nc_vals, m_rhos, "o-", color="#1f77b4", markersize=5, linewidth=1.5)

    ax.axvline(200, color="0.45", ls="--", lw=1, alpha=0.7)
    ax.annotate("real sample size\n(~200 cells/type)", xy=(210, 0.15),
                fontsize=7, color="0.35", ha="left", va="bottom")

    ax.axhline(0.15, color="#d62728", ls=":", lw=1, alpha=0.7)
    ax.annotate("real replication ρ ≈ 0.15",
                xy=(max(nc_vals), 0.19), fontsize=7,
                color="#d62728", ha="right", va="bottom")

    ax.set_xlabel("Cells per type")
    ax.set_ylabel("Test-retest Spearman ρ")
    ax.set_xscale("log")
    ax.set_ylim(-0.15, 1.08)
    ax.set_title("C.  Ranking stability", fontweight="bold", loc="left")

    # Panel D: Null calibration
    ax = axes[1, 1]
    p_sorted = np.sort(null_cal["p_values"])
    n = len(p_sorted)
    expected = (np.arange(1, n + 1) - 0.5) / n

    ax.scatter(expected, p_sorted, s=6, alpha=0.35, color="#1f77b4",
               edgecolors="none")
    ax.plot([0, 1], [0, 1], "k--", lw=1, label="Ideal (uniform)")
    ax.set_xlabel("Expected quantile (uniform)")
    ax.set_ylabel("Observed p-value")
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.03, 1.03)
    ax.set_aspect("equal")

    rej = null_cal["rejection_rate_05"]
    ax.annotate(
        f"α = 0.05 rejection: {rej:.1%}\n(expected: 5.0%)",
        xy=(0.05, 0.92), fontsize=7, va="top",
    )
    ax.legend(loc="lower right", fontsize=7, frameon=False)
    ax.set_title("D.  Null calibration", fontweight="bold", loc="left")

    plt.tight_layout(h_pad=2.0, w_pad=1.5)

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(Path(output_dir) / f"figS7_simulation_study.{ext}",
                    bbox_inches="tight")
    fig.savefig(SCRIPT_DIR / "simulation_study_figure.png", bbox_inches="tight")
    plt.close(fig)
    print(f"Figures saved to {output_dir}/")


if __name__ == "__main__":
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    results = load_results(path)
    make_figure(results)
