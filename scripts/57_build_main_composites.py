#!/usr/bin/env python3
"""Build the panel-embedded main-figure composites under their CANONICAL names.

This is the tracked, self-contained producer for the five main composites that
embed pre-rendered panel PNGs. Each function is lifted verbatim (proven
vector-identical, CellWarp Stage H2b) from the previously-gitignored
scripts/assemble_figures.py -- the true producer of the shipped figures -- so
the canonical composites regenerate from tracked code. Figure 4 (human-macaque)
is a native vector composite built separately by scripts/48_build_fig6_K12.py
and is NOT produced here.

The source file used a scrambled internal numbering; this producer un-scrambles
it. Mapping (source function -> canonical output):
    assemble_fig1 -> fig1_global_coherence   (build_fig1)
    assemble_fig3 -> fig2_two_layer          (build_fig2)
    assemble_fig4 -> fig3_replication        (build_fig3)
    assemble_fig2 -> fig5_rigidity_ranking   (build_fig5)
    assemble_fig6 -> fig6_l1000_nulls        (build_fig6)

fig6's panels carry the post biology-review divergence / Procrustes-residual
labels. Styling comes from the shared src/cellwarp/figure_style.py (NOT
assemble_figures.py). Portable: all paths are Path(__file__)-relative.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg

sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))
from cellwarp.figure_style import apply_style, add_panel_label, save_figure, COL2

apply_style()

BASE = Path(__file__).parent.parent
PANELS = BASE / 'figures' / 'panels'
MAIN_OUT = BASE / 'figures' / 'main'
MAIN_OUT.mkdir(parents=True, exist_ok=True)


def load_panel(name):
    """Load a rendered panel PNG image (figures/panels/<name>.png)."""
    path = PANELS / f'{name}.png'
    if path.exists():
        return mpimg.imread(str(path))
    print(f"  WARNING: Panel {name}.png not found")
    return None


def embed_panel(ax, img, label=None):
    """Embed a panel image into an axes and optionally add a panel label."""
    if img is not None:
        ax.imshow(img)
    ax.axis('off')
    if label:
        add_panel_label(ax, label, x=-0.02, y=1.02)


def build_fig1():
    """Figure 1 (global coherence) <- assemble_fig1. 5 panels A-E:
    pipeline schematic + 1M null + lineage-stratified null + bootstrap + LOOCV."""
    print("Assembling Figure 1: Global Procrustes Coherence...")
    fig = plt.figure(figsize=(COL2, COL2 * 1.1))
    # Row 0 (Panel A) height ratio 0.65 so the slot aspect matches the
    # schematic's source aspect (~4.29); otherwise the slot is too wide.
    gs = gridspec.GridSpec(3, 2, figure=fig, hspace=0.35, wspace=0.3,
                           height_ratios=[0.65, 1, 1])
    embed_panel(fig.add_subplot(gs[0, :]), load_panel('fig1a_pipeline_schematic'), 'A')
    embed_panel(fig.add_subplot(gs[1, 0]), load_panel('fig1b_null_1M'), 'B')
    embed_panel(fig.add_subplot(gs[1, 1]), load_panel('fig1c_lineage_stratified'), 'C')
    embed_panel(fig.add_subplot(gs[2, 0]), load_panel('fig1d_bootstrap'), 'D')
    embed_panel(fig.add_subplot(gs[2, 1]), load_panel('fig1e_loocv'), 'E')
    save_figure(fig, MAIN_OUT / 'fig1_global_coherence', tight=False)


def build_fig2():
    """Figure 2 (two-layer geometric conservation) <- assemble_fig3. 4 panels A-D:
    ellipsoid heatmap + pre/post rotation + layer nulls + per-type scatter."""
    print("Assembling Figure 2: Two-Layer Geometric Conservation...")
    fig = plt.figure(figsize=(COL2, COL2 * 0.85))
    # width_ratios match source aspects: Panel A taller/wider, Panel C wide-short.
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.22, wspace=0.22,
                           width_ratios=[1.15, 1.0])
    embed_panel(fig.add_subplot(gs[0, 0]), load_panel('fig3a_ellipsoid_heatmap'), 'A')
    embed_panel(fig.add_subplot(gs[0, 1]), load_panel('fig3b_pre_post'), 'B')
    embed_panel(fig.add_subplot(gs[1, 0]), load_panel('fig3c_layer_nulls'), 'C')
    embed_panel(fig.add_subplot(gs[1, 1]), load_panel('fig3d_layer_scatter'), 'D')
    save_figure(fig, MAIN_OUT / 'fig2_two_layer', tight=False)


def build_fig3():
    """Figure 3 (replication across datasets) <- assemble_fig4. 5 panels A-E:
    Sun2023 + PanSci + CellHint nulls (row A/B/C), replication summary (D),
    donor-split within-species control (E)."""
    print("Assembling Figure 3: Replication Across Datasets...")
    fig = plt.figure(figsize=(7.2, 8.5))
    # 3-row layout: row 0 = A/B/C (3 equal panels), row 1 = D (full width),
    # row 2 = E (full width).
    gs_outer = gridspec.GridSpec(
        3, 1, figure=fig,
        height_ratios=[1.0, 1.0, 1.2],
        hspace=0.28,
        left=0.02, right=0.98, top=0.98, bottom=0.04,
    )
    gs_top = gs_outer[0].subgridspec(1, 3, wspace=0.30)
    embed_panel(fig.add_subplot(gs_top[0, 0]), load_panel('fig4a_sun2023_null'), 'A')
    embed_panel(fig.add_subplot(gs_top[0, 1]), load_panel('fig4b_pansci_null'), 'B')
    embed_panel(fig.add_subplot(gs_top[0, 2]), load_panel('fig4c_cellhint_null'), 'C')
    embed_panel(fig.add_subplot(gs_outer[1]), load_panel('fig4d_replication_summary'), 'D')
    embed_panel(fig.add_subplot(gs_outer[2]), load_panel('fig2e_donor_split'), 'E')
    save_figure(fig, MAIN_OUT / 'fig3_replication', tight=False)


def build_fig5():
    """Figure 5 (per-type divergence / rigidity ranking) <- assemble_fig2.
    2 panels A-B: ranking bar chart + cell-count confound (CellMarker demoted
    to Fig S6)."""
    print("Assembling Figure 5: Divergence / Rigidity Ranking...")
    fig = plt.figure(figsize=(7.2, 5.5))
    gs = gridspec.GridSpec(2, 2, figure=fig,
                           height_ratios=[2.2, 1.0],
                           hspace=0.22, wspace=0.22,
                           left=0.02, right=0.98, top=0.97, bottom=0.03)
    embed_panel(fig.add_subplot(gs[0, :]), load_panel('fig2a_rigidity_ranking'), 'A')
    embed_panel(fig.add_subplot(gs[1, :]), load_panel('fig2c_cellcount'), 'B')
    save_figure(fig, MAIN_OUT / 'fig5_rigidity_ranking', tight=False)


def build_fig6():
    """Figure 6 (L1000 robustness + mechanistic nulls) <- assemble_fig6.
    2 panels A-B; x-axes read 'divergence ranking correlation' (A) and
    'Spearman rho with Procrustes residual' (B) -- post biology-review labels."""
    print("Assembling Figure 6: L1000 & Mechanistic Nulls...")
    fig = plt.figure(figsize=(COL2, COL2 * 0.5))
    gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.3)
    embed_panel(fig.add_subplot(gs[0, 0]), load_panel('fig6a_l1000'), 'A')
    embed_panel(fig.add_subplot(gs[0, 1]), load_panel('fig6b_mechanistic_nulls'), 'B')
    save_figure(fig, MAIN_OUT / 'fig6_l1000_nulls', tight=False)


def main():
    build_fig1()
    build_fig2()
    build_fig3()
    build_fig5()
    build_fig6()


if __name__ == '__main__':
    main()
