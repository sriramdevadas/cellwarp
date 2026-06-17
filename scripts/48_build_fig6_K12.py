#!/usr/bin/env python3
"""Build Figure 6 composite (human-macaque extension, K12 design).

Three panels:
  (A) Null distribution histogram for human-macaque primary (Qu-only, 12 types).
      obs/null = 0.810, p = 0.0043, n = 12.
  (B) Null distribution histogram for no-immune sensitivity (7 types).
      obs/null = 0.733, p = 0.013, n = 7.
  (C) Hepatocyte rank-reversal slope chart across species pairs:
      human-mouse 12-type matched control (rank 12/12, most rigid) →
      human-macaque 12-type (rank 1/12, most flexible, 47.3% of SSR).
      Other 11 types shown as gray background lines.

Writes:
  figures/main/fig4_human_macaque.pdf            (canonical; Figure_4.pdf source)
  figures/main/fig4_human_macaque_polished.pdf  (identical; polishing step was
                                                  historically a cosmetic overlay)
  figures/panels/fig6a_macaque_primary.{pdf,png}
  figures/panels/fig6b_macaque_sensitivity.{pdf,png}
  figures/panels/fig6c_hepatocyte_rigidity.{pdf,png}
"""
from __future__ import annotations
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np

PROJECT = Path(__file__).resolve().parent.parent
OUT_MAIN = PROJECT / "figures/main"
OUT_PANELS = PROJECT / "figures/panels"
MACAQUE = PROJECT / "output/macaque_pipeline"

NULL_QU12 = MACAQUE / "null_distribution_qu12.npy"
NULL_QU7 = MACAQUE / "null_distribution_qu7_D1.npy"
M1 = MACAQUE / "m1_close_table1_summary.json"

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
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "pdf.fonttype": 42,  # embed TrueType
    "ps.fonttype": 42,
})

C_BLUE = "#1f77b4"
C_ORANGE = "#ff7f0e"
C_GRAY = "#7f7f7f"
C_LIGHTGRAY = "#cccccc"
C_RED = "#d62728"
C_DARKGRAY = "#404040"


def _null_panel(ax, null, obs, obs_null, p, n_types, n_perms, title, extra_inset=None):
    """Draw a null-distribution histogram with observed line + inset stats.

    Inset placed in whichever top corner is farther from the histogram bulk
    (observed is always left of bulk for our panels, so the inset goes top-left
    near the red line to avoid the main mass of bars). y-limit is extended 35%
    above the tallest bar so the inset box has clear vertical headroom.
    """
    counts, edges, _ = ax.hist(null, bins=50, color=C_LIGHTGRAY,
                               edgecolor=C_DARKGRAY, linewidth=0.4,
                               density=True)
    ax.axvline(obs, color=C_RED, linewidth=1.6, zorder=10)

    # p-value formatting per Part 4A convention (superscript / floor indication)
    if n_perms == 10_000 and p <= 1.5 / n_perms:
        p_str = rf"$p < 10^{{-4}}$"
    elif p < 1e-3:
        mantissa, exp = f"{p:.1e}".split("e")
        p_str = rf"$p = {mantissa} \times 10^{{{int(exp)}}}$"
    else:
        p_str = f"p = {p:.4f}"

    # Extra vertical headroom so inset box doesn't crash into the tallest bar
    ax.set_ylim(0, counts.max() * 1.45)
    # Extend x-axis leftward by ~10% so the "observed" label has room to sit
    # to the left of the red line without clipping the panel boundary.
    xmin_cur, xmax_cur = ax.get_xlim()
    x_span = xmax_cur - xmin_cur
    ax.set_xlim(xmin_cur - x_span * 0.10, xmax_cur)

    inset_lines = [
        f"obs/null = {obs_null:.3f}",
        p_str,
        f"n = {n_types} cell types",
        f"{n_perms:,} permutations",
    ]
    if extra_inset:
        inset_lines.append(extra_inset)
    # Inset: top-right of panel. Observed is on the far left (obs ≪ null bulk),
    # so right side is farthest from the red line; headroom above histogram max
    # keeps the box from overlapping any bar.
    ax.text(0.98, 0.97, "\n".join(inset_lines),
            transform=ax.transAxes, ha="right", va="top", fontsize=6.5,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=C_DARKGRAY, lw=0.5))

    # "observed" label placed to the LEFT of the red line, in the clear space
    # below the histogram bulk. Arrow points from label to red line.
    y_top = ax.get_ylim()[1]
    ax.annotate(f"observed\n({obs:.2f})",
                xy=(obs, y_top * 0.45),
                xytext=(obs - x_span * 0.04, y_top * 0.75),
                ha="right", fontsize=6.5, color=C_RED,
                arrowprops=dict(arrowstyle="->", color=C_RED, lw=0.6,
                                shrinkA=0, shrinkB=2))

    ax.set_xlabel("Procrustes distance", fontsize=8)
    ax.set_ylabel("Density", fontsize=8)
    ax.set_title(title, fontsize=8.5, pad=6)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def panel_A(ax):
    null = np.load(NULL_QU12)
    _null_panel(ax, null, obs=39.88, obs_null=0.810, p=0.0043,
                n_types=12, n_perms=10_000,
                title="A  Human-macaque primary (Qu only, 12 types)",
                extra_inset="13,927 three-way orthologs")


def panel_B(ax):
    null = np.load(NULL_QU7)
    # No extra_inset: "No-immune sensitivity (7 types)" title already says
    # what this panel is. Previous "excludes 5 immune types added in K12"
    # line widened the inset box enough to overlap the "observed" label.
    _null_panel(ax, null, obs=20.56, obs_null=0.733, p=0.0128,
                n_types=7, n_perms=10_000,
                title="B  No-immune sensitivity (7 types)")


def panel_C(ax):
    """Slope chart: hepatocyte rank reversal + 11 gray background lines."""
    summ = json.loads(M1.read_text())
    mac_res = {r["type"]: r["magnitude"]
               for r in summ["human_macaque_12type"].get("per_type_residuals_ranked", [])}
    # If the key name differs, reach into the full JSONs
    if not mac_res:
        mac_raw = json.loads((MACAQUE / "reconstruction_qu12_results.json").read_text())
        mac_res = {r["type"]: r["magnitude"] for r in mac_raw["per_type_residuals_ranked"]}
    hm_raw = json.loads((MACAQUE / "human_mouse_12type_control.json").read_text())
    mm_res = {r["type"]: r["magnitude"]
              for r in hm_raw["control_16959"]["per_type_residuals_ranked"]}

    types = sorted(set(mac_res) & set(mm_res))
    assert len(types) == 12, f"expected 12 types, got {len(types)}"
    # Compute ranks (1 = largest residual = most flexible)
    def ranks_from(res):
        order = sorted(res, key=lambda t: -res[t])
        return {t: i + 1 for i, t in enumerate(order)}
    mac_rank = ranks_from(mac_res)
    mm_rank = ranks_from(mm_res)

    # Hepatocyte first + last so it renders on top
    non_hep = [t for t in types if t != "hepatocyte"]
    for t in non_hep:
        ax.plot([0, 1], [mm_rank[t], mac_rank[t]],
                color=C_LIGHTGRAY, linewidth=0.8, marker="o",
                markerfacecolor=C_LIGHTGRAY, markeredgecolor=C_DARKGRAY,
                markeredgewidth=0.4, markersize=3, zorder=2)
    # Hepatocyte highlighted
    ax.plot([0, 1], [mm_rank["hepatocyte"], mac_rank["hepatocyte"]],
            color=C_RED, linewidth=1.6, marker="o",
            markerfacecolor=C_RED, markeredgecolor="white",
            markeredgewidth=0.6, markersize=6, zorder=10)
    # Endpoint labels — pushed clear of the red markers. Marker size 6 at
    # 7.5×3.6 fig = ~0.17 x-units radius, so xytext at ±0.04 from marker
    # center put the label's inside edge INSIDE the marker footprint and
    # obscured the "h" of "hepatocyte" / the "4" of "47.3% of SSR". New
    # offsets (±0.15) place the label in clear air beyond the marker.
    ax.annotate("hepatocyte",
                xy=(0, mm_rank["hepatocyte"]),
                xytext=(-0.15, mm_rank["hepatocyte"]),
                ha="right", va="center", fontsize=6.5,
                color=C_RED, fontweight="bold")
    ax.annotate("47.3% of SSR",
                xy=(1, mac_rank["hepatocyte"]),
                xytext=(1.15, mac_rank["hepatocyte"]),
                ha="left", va="center", fontsize=6.5,
                color=C_RED, fontweight="bold")

    ax.set_xticks([0, 1])
    # Compacted tick labels (2 lines, not 3) so they don't overrun each
    # other when the x-axis is stretched to accommodate the endpoint
    # labels. Panel C is also given a larger width ratio in save_composite.
    ax.set_xticklabels(["human-mouse\n(matched 12-type)",
                        "human-macaque\n(primary 12-type)"],
                       fontsize=6.5)
    # Widen x-axis so both endpoint labels fit in clear air past the
    # markers. The endpoint labels are offset ±0.15 from each marker and
    # extend outward another ~0.6 x-units.
    ax.set_xlim(-0.85, 1.85)
    # Rigid-at-top convention, matching Fig 5A and Fig 1E: most rigid
    # (rank 12) at top, most flexible (rank 1) at bottom. No axis inversion.
    # Lock y-ticks to valid rank values so the caption headroom doesn't
    # introduce a phantom tick.
    ax.set_yticks(list(range(1, 13)))
    # Leave headroom below the data (rank-1 end) for the gray-lines caption.
    ax.set_ylim(-0.6, 12.6)
    ax.set_ylabel("divergence rank (1 = most divergent)", fontsize=8)
    ax.set_title("C  Hepatocyte rank reversal across species pairs",
                 fontsize=8.5, pad=6)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    # Gray-lines caption placed BELOW the data (inside extra headroom added
    # below the rank-1 end) so it doesn't crash into the title or markers.
    ax.text(0.5, -0.3,
            "Gray: 11 other matched types "
            "(overall Spearman ρ = 0.147, p = 0.649)",
            ha="center", va="center", fontsize=6.0, color=C_DARKGRAY,
            style="italic", transform=ax.transData)


def save_individual_panels():
    for name, fn in [("fig6a_macaque_primary", panel_A),
                     ("fig6b_macaque_sensitivity", panel_B),
                     ("fig6c_hepatocyte_rigidity", panel_C)]:
        fig, ax = plt.subplots(figsize=(3.4, 3.0))
        fn(ax)
        fig.tight_layout()
        for ext in ("pdf", "png"):
            fig.savefig(OUT_PANELS / f"{name}.{ext}", dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  wrote panel: {name}.{{pdf,png}}")


def save_composite():
    # Wider + taller composite: previous 7.0×3.2 was too cramped — titles
    # sat flush against the inset, and Panel C's endpoint labels clipped.
    # Panel C gets a 1.3× width ratio so the widened xlim (needed for the
    # hepatocyte and 47.3% endpoint labels) does not compress the two
    # x-tick labels into each other.
    fig = plt.figure(figsize=(8.2, 3.6))
    gs = fig.add_gridspec(1, 3, wspace=0.40, width_ratios=[1.0, 1.0, 1.3],
                          left=0.05, right=0.98, top=0.88, bottom=0.22)
    panel_A(fig.add_subplot(gs[0, 0]))
    panel_B(fig.add_subplot(gs[0, 1]))
    panel_C(fig.add_subplot(gs[0, 2]))
    out_base = OUT_MAIN / "fig4_human_macaque"
    out_polished = OUT_MAIN / "fig4_human_macaque_polished"
    for out in (out_base, out_polished):
        for ext in ("pdf", "png"):
            fig.savefig(f"{out}.{ext}", dpi=300, bbox_inches="tight")
        print(f"  wrote: {out.name}.{{pdf,png}}")
    plt.close(fig)


def main():
    OUT_MAIN.mkdir(parents=True, exist_ok=True)
    OUT_PANELS.mkdir(parents=True, exist_ok=True)
    save_individual_panels()
    save_composite()


if __name__ == "__main__":
    main()
