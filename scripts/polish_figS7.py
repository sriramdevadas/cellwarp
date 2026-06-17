#!/usr/bin/env python3
"""
Generate Cell Systems-polished version of Figure S7 (simulation study).

Fixes applied vs. original:
  1. Font: DejaVuSans → Arial (all text)
  2. Axis label size: 10pt → 8pt (spec: 6-8pt at print)
  3. Panel title size: 10pt → 9pt bold (spec: 8-9pt)
  4. Figure width: 7.2" → 6.85" (174mm, 2-column max)
  5. Red-green colors replaced with colorblind-safe palette:
       Panel A: green (#2ca02c) → teal (#00897B)
       Panel A: red α-line → gray (#757575)
       Panel B: red 50-cells (#d62728) → dark orange (#E8820E)
       Panel B: green 500-cells (#2ca02c) → teal (#00897B)
       Panel C: red annotation (#d62728) → orange (#E8820E)
  6. pdf.fonttype=42 (TrueType, editable in Illustrator)

NO scientific content changes: same data, same axis ranges, same annotations.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ── Paths ──
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
SIM_DIR = PROJECT_ROOT / "analysis" / "simulation_study"
FIG_DIR = PROJECT_ROOT / "figures" / "supplementary"

# ── Constants (unchanged from original) ──
REAL_OBS_NULL = 0.522
POWER_N_TYPES = [15, 25, 35]
RECOVERY_N_CELLS = [50, 200, 500, 2000]

# ── Cell Systems color palette (colorblind-safe, no red-green) ──
C_BLUE = "#3574B0"
C_ORANGE = "#E8820E"
C_TEAL = "#00897B"
C_PURPLE = "#9467bd"
C_GRAY = "#757575"


def main():
    results_path = SIM_DIR / "simulation_results.json"
    with open(results_path) as f:
        results = json.load(f)

    power = results["power_curve"]
    recovery = results["ranking_recovery"]
    stability = results["stability"]
    null_cal = results["null_calibration"]
    real_signal = results["calibration"]["estimated_real_signal"]

    # ── Cell Systems rcParams ──
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.labelsize": 8,         # was 10 → 8 (spec: 6-8pt)
        "axes.titlesize": 9,         # was 10 → 9 (spec: 8-9pt)
        "xtick.labelsize": 7,        # was 8 → 7
        "ytick.labelsize": 7,        # was 8 → 7
        "legend.fontsize": 7,        # was 7.5 → 7
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "lines.linewidth": 1.0,
        "pdf.fonttype": 42,          # TrueType (editable text in PDF)
        "ps.fonttype": 42,
    })

    # 174mm / 25.4 = 6.85" (2-column max)
    fig, axes = plt.subplots(2, 2, figsize=(6.85, 6.2))

    # Colorblind-safe palettes (no red-green pairs)
    type_colors = {15: C_BLUE, 25: C_ORANGE, 35: C_TEAL}  # was green→teal
    cell_colors = {50: C_ORANGE, 200: C_BLUE, 500: C_TEAL, 2000: C_PURPLE}

    # ── Panel A: Power ──
    ax = axes[0, 0]
    for nt in POWER_N_TYPES:
        sub = [r for r in power if r["n_types"] == nt]
        sigs = [r["signal_strength"] for r in sub]
        det = [r["detection_rate"] for r in sub]
        ax.plot(sigs, det, "o-", color=type_colors[nt], markersize=4,
                linewidth=1.5, label=f"n = {nt} types")

    ax.axvline(real_signal, color="0.45", ls="--", lw=1, alpha=0.7)
    ann_ha = "left" if real_signal < 5 else "right"
    # Small x-offset so the label doesn't sit directly on the dashed line.
    ann_x = real_signal + (0.3 if ann_ha == "left" else -0.3)
    ax.annotate(f"real data signal\nobs/null = {REAL_OBS_NULL}",
                xy=(ann_x, 0.50), fontsize=6.5, color="0.35",
                ha=ann_ha, va="center")
    ax.axhline(0.05, color=C_GRAY, ls=":", lw=0.8, alpha=0.6,
               label="α = 0.05")                       # was red → gray
    ax.set_xlabel("Signal strength")
    ax.set_ylabel("Detection rate (p < 0.05)")
    ax.set_xlim(-0.5, 10.5)
    ax.set_ylim(-0.05, 1.08)
    ax.legend(loc="center right", frameon=False)
    ax.set_title("A.  Detection power", fontweight="bold", loc="left")

    # ── Panel B: Ranking recovery ──
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
    # Match Panel A's x-range so the two panels read on a common scale.
    ax.set_xlim(-0.5, 10.5)
    ax.set_ylim(-0.15, 1.08)
    # Legend nudged left from the axes-right edge.
    ax.legend(loc="lower right", bbox_to_anchor=(0.92, 0.02), frameon=False)
    ax.set_title("B.  Ranking recovery", fontweight="bold", loc="left")

    # ── Panel C: Ranking stability ──
    ax = axes[1, 0]
    nc_vals = [r["n_cells_per_type"] for r in stability]
    m_rhos = [r["mean_rho"] for r in stability]
    ci_lo = [r["ci_lower"] for r in stability]
    ci_hi = [r["ci_upper"] for r in stability]

    ax.fill_between(nc_vals, ci_lo, ci_hi, alpha=0.20, color=C_BLUE)
    ax.plot(nc_vals, m_rhos, "o-", color=C_BLUE, markersize=5, linewidth=1.5)

    ax.axvline(200, color="0.45", ls="--", lw=1, alpha=0.7)
    # Multiplicative x-offset (200 * 1.1) parallels Panel A's linear-axis
    # +0.3-unit shift on the log axis: same visual spacing in axes-fraction
    # terms (~2.7% of plot width). va/ha aligned to Panel A's convention.
    ax.annotate("real sample size\n(~200 cells/type)", xy=(200 * 1.1, 0.5),
                fontsize=6.5, color="0.35", ha="left", va="center")

    ax.axhline(0.15, color=C_ORANGE, ls=":", lw=1, alpha=0.7)  # was red→orange
    # Two-line annotation matching Panel A's "real data signal\nobs/null=..."
    # convention; sits just below the dashed line at the right edge.
    ax.annotate("real replication\nρ ≈ 0.15",
                xy=(max(nc_vals), 0.15), fontsize=6.5,
                color=C_ORANGE, ha="right", va="top")          # was red→orange

    ax.set_xlabel("Cells per type")
    ax.set_ylabel("Test-retest Spearman ρ")
    ax.set_xscale("log")
    ax.set_ylim(-0.15, 1.08)
    ax.set_title("C.  Ranking stability", fontweight="bold", loc="left")

    # ── Panel D: Null calibration ──
    ax = axes[1, 1]
    p_sorted = np.sort(null_cal["p_values"])
    n = len(p_sorted)
    expected = (np.arange(1, n + 1) - 0.5) / n

    ax.scatter(expected, p_sorted, s=6, alpha=0.35, color=C_BLUE,
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
        xy=(0.05, 0.92), fontsize=6.5, va="top",
    )
    ax.legend(loc="lower right", fontsize=7, frameon=False)
    ax.set_title("D.  Null calibration", fontweight="bold", loc="left")

    plt.tight_layout(h_pad=2.0, w_pad=1.5)

    # Save polished versions (alongside originals, not overwriting)
    out_stem = FIG_DIR / "figS7_simulation_study_polished"
    fig.savefig(str(out_stem) + ".pdf", format="pdf",
                bbox_inches="tight", pad_inches=0.05)
    fig.savefig(str(out_stem) + ".png", format="png", dpi=300,
                bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"Saved: {out_stem}.pdf and {out_stem}.png")


if __name__ == "__main__":
    main()
