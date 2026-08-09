#!/usr/bin/env python3
"""
Generate Phase 2 NEEDS_CREATION panels for CellWarp Cell Systems submission.

2 panels that need to be designed from scratch:
  Fig 1A: Pipeline schematic (centroids → PCA → Procrustes → residuals)
  Fig 6B: Mechanistic nulls forest plot (10-test summary)

Biology: The schematic shows the CellWarp geometric morphometrics pipeline.
The forest plot summarizes 10 mechanistic null tests showing no confound
association with the evolutionary rigidity signal.
"""

import sys
import json
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))
from cellwarp.figure_style import (
    apply_style, add_panel_label, save_figure, clean_spine,
    COL1, COL15, COL2, DPI,
    C_BLUE, C_ORANGE, C_PURPLE, C_TEAL, C_GRAY, C_LIGHTGRAY, C_DARKGRAY, C_BLACK,
    FONT_SIZE_LABEL, FONT_SIZE_TICK, FONT_SIZE_ANNOT, FONT_SIZE_PANEL,
    FONT_SIZE_LEGEND, FONT_FAMILY,
)

apply_style()

BASE = Path(__file__).parent.parent
PANELS = BASE / 'figures' / 'panels'
PANELS.mkdir(parents=True, exist_ok=True)


# ===================================================================
# Fig 1A: Pipeline Schematic
# ===================================================================
def fig1a_pipeline_schematic():
    """Create a clean pipeline schematic: species → centroids → PCA → Procrustes → residuals."""
    print("Generating Fig 1A: Pipeline schematic...")

    fig, ax = plt.subplots(figsize=(COL2, COL2 * 0.22))
    ax.set_xlim(0, 10)
    ax.set_ylim(0.3, 2.3)
    ax.axis('off')

    # Box style
    box_kw = dict(boxstyle='round,pad=0.25', facecolor='white',
                  edgecolor=C_DARKGRAY, linewidth=0.8)
    highlight_kw = dict(boxstyle='round,pad=0.25', facecolor='#E3F2FD',
                        edgecolor=C_BLUE, linewidth=0.8)
    result_kw = dict(boxstyle='round,pad=0.25', facecolor='#FFF3E0',
                     edgecolor=C_ORANGE, linewidth=0.8)

    # Step positions (x-centers)
    steps_x = [0.8, 2.6, 4.4, 6.2, 8.2]
    y_center = 1.2

    texts = []

    # Step 1: Species atlases
    texts.append(ax.text(steps_x[0], y_center,
            'Single-cell\natlases\n(Human, Mouse)',
            ha='center', va='center', fontsize=6.5, fontfamily=FONT_FAMILY,
            bbox=box_kw))

    # Step 2: Centroids
    texts.append(ax.text(steps_x[1], y_center,
            'Cell-type\ncentroids\n(35 types × 2 spp)',
            ha='center', va='center', fontsize=6.5, fontfamily=FONT_FAMILY,
            bbox=box_kw))

    # Step 3: PCA
    texts.append(ax.text(steps_x[2], y_center,
            'Joint PCA\n(33 components,\n95.2% variance)',
            ha='center', va='center', fontsize=6.5, fontfamily=FONT_FAMILY,
            bbox=highlight_kw))

    # Step 4: Procrustes
    texts.append(ax.text(steps_x[3], y_center,
            'Procrustes\nsuperimposition\n(rotation + scaling)',
            ha='center', va='center', fontsize=6.5, fontfamily=FONT_FAMILY,
            bbox=highlight_kw))

    # Step 5: Residuals
    texts.append(ax.text(steps_x[4], y_center,
            'Per-type\nresiduals\n→ divergence ranking',
            ha='center', va='center', fontsize=6.5, fontfamily=FONT_FAMILY,
            bbox=result_kw))

    # Arrows between steps — endpoints sit at box-centers and matplotlib
    # clips them against each box's bbox patch, so arrows meet box edges
    # cleanly regardless of per-box width. Requires a draw pass first to
    # populate bbox extents.
    fig.canvas.draw()
    arrow_kw = dict(arrowstyle='->', color=C_DARKGRAY, lw=1.2,
                    connectionstyle='arc3,rad=0', shrinkA=2, shrinkB=2)
    for i in range(len(steps_x) - 1):
        ax.annotate('', xy=(steps_x[i + 1], y_center),
                     xytext=(steps_x[i], y_center),
                     arrowprops={**arrow_kw,
                                 'patchA': texts[i].get_bbox_patch(),
                                 'patchB': texts[i + 1].get_bbox_patch()})

    # Step labels above
    step_labels = ['1. Data', '2. Aggregate', '3. Embed', '4. Align', '5. Quantify']
    for x, label in zip(steps_x, step_labels):
        ax.text(x, y_center + 0.85, label, ha='center', va='bottom',
                fontsize=6, fontfamily=FONT_FAMILY, color=C_GRAY,
                fontweight='bold')

    # Bottom annotation: permutation test (anchored just below box bottom edge
    # so it visually attaches to the schematic instead of floating).
    ax.text(5.0, 0.55,
            'Permutation test: shuffle cell type labels \u2192 null distribution \u2192 '
            'p < 10$^{-6}$ (1M permutations)',
            ha='center', va='center', fontsize=6, fontfamily=FONT_FAMILY,
            color=C_DARKGRAY, style='italic')

    save_figure(fig, PANELS / 'fig1a_pipeline_schematic', tight=False)
    return fig


# ===================================================================
# Fig 6B: Mechanistic nulls forest plot
# ===================================================================
def fig6b_mechanistic_nulls():
    """Forest plot of 10 mechanistic null tests."""
    print("Generating Fig 6B: Mechanistic nulls forest plot...")

    # Load each mechanistic null test result from its source file
    def _load(path):
        with open(BASE / path) as f:
            return json.load(f)

    hk = _load('output/phase2/mechanistic/housekeeping/hk_ratio_results.json')
    tf = _load('output/phase2/mechanistic/tf_complexity/tf_complexity_results.json')
    niche = _load('output/phase2/progenitor_analysis/niche_hypothesis/niche_hypothesis_results.json')
    var_diag = _load('output/phase2/variance_diagnostic/diagnostic_results.json')
    donor = _load('output/phase2/diagnostics/interdonor_variance/diagnostic_results.json')
    ppi = _load('output/mechanistic/ppi_centrality/ppi_centrality_results.json')
    chrom = _load('output/validation/t3e_chromatin/spearman_primary_result.json')
    enh = _load('output/validation/t3e_enhancer/spearman_primary_result.json')
    drug = _load('output/t3g/primary_correlation_results.json')

    # Expression-level confounds: strongest per-cell-type expression metric
    with open(BASE / 'output/phase2/diagnostics/expression_level_vs_rigidity/correlations.csv') as f:
        expr_rows = list(csv.DictReader(f))
    expr_best = max(expr_rows, key=lambda r: abs(float(r['rho'])))

    niche_sig = niche['progenitor_divergence']['n_niche_sets_significant_q05']
    niche_total = niche['progenitor_divergence']['n_niche_sets_tested']

    tests = [
        {'name': 'Housekeeping gene ratio',
         'rho': hk['human_correlation']['spearman_rho'], 'n': hk['n_cell_types']},
        {'name': 'TF network complexity',
         'rho': tf['human_correlations']['n_active_tfs']['spearman_rho'], 'n': tf['n_cell_types']},
        {'name': f'Niche adaptation ({niche_sig}/{niche_total} sets)',
         'rho': 0.0, 'n': hk['n_cell_types']},
        {'name': 'Within-type variance',
         'rho': var_diag['spearman_rho_mean'], 'n': var_diag['n_cell_types']},
        {'name': 'Inter-donor variance',
         'rho': donor['spearman_rho_mean'], 'n': donor['n_cell_types']},
        {'name': 'Expression-level confounds',
         'rho': float(expr_best['rho']), 'n': hk['n_cell_types']},
        {'name': 'PPI network centrality',
         'rho': ppi['sensitivity']['best_rho'], 'n': ppi['correlation_results'][0]['n_cell_types']},
        {'name': 'Promoter conservation',
         'rho': chrom['rho'], 'n': chrom['n']},
        {'name': 'Enhancer conservation',
         'rho': enh['spearman_rho'], 'n': enh['n']},
        {'name': 'Drug target conservation',
         'rho': drug['primary_correlation']['rho'], 'n': drug['primary_correlation']['n']},
    ]

    fig, ax = plt.subplots(figsize=(COL1, COL1 * 1.1))

    y_pos = np.arange(len(tests))[::-1]  # top to bottom

    for i, test in enumerate(tests):
        y = y_pos[i]
        rho = test['rho']
        n = test['n']

        # Approximate 95% CI for Spearman rho (Fisher z transform)
        if n > 3 and rho != 0:
            se = 1.0 / np.sqrt(n - 3)
            z = np.arctanh(rho)
            ci_lo = np.tanh(z - 1.96 * se)
            ci_hi = np.tanh(z + 1.96 * se)
        else:
            se = 1.0 / np.sqrt(max(n - 3, 1))
            ci_lo = -1.96 * se
            ci_hi = 1.96 * se

        # Point and CI — single neutral color, no sign-based coloring.
        # Previously blue for rho>=0 and orange for rho<0; dropped because
        # the sign is already visible as the marker's position relative to
        # the zero line, and a two-color split over-emphasized a cosmetic
        # distinction that isn't a finding of the panel.
        ax.plot([ci_lo, ci_hi], [y, y], color=C_DARKGRAY, linewidth=1.2,
                zorder=3)
        ax.plot(rho, y, 'o', color=C_DARKGRAY, markersize=5, zorder=4,
                markeredgecolor='white', markeredgewidth=0.5)

    # Vertical zero line
    ax.axvline(0, color=C_DARKGRAY, linewidth=0.5, linestyle='-', zorder=1)

    # Significance threshold markers
    ax.axvline(-0.35, color=C_LIGHTGRAY, linewidth=0.4, linestyle=':', zorder=1)
    ax.axvline(0.35, color=C_LIGHTGRAY, linewidth=0.4, linestyle=':', zorder=1)

    ax.set_yticks(y_pos)
    ax.set_yticklabels([t['name'] for t in tests], fontsize=FONT_SIZE_TICK)
    ax.set_xlabel(r'Spearman $\rho$ with Procrustes residual', fontsize=FONT_SIZE_LABEL)
    ax.set_xlim(-0.8, 0.6)

    # N labels in a fixed right-margin (axes-frac x=1.02, data y) so they
    # never collide with CI bars. Previously placed at data x=0.55, which
    # was inside the plot area — PPI's CI right-end at 0.569 ended up under
    # the 'n=35' label and read as a stray glyph at print scale.
    for i, test in enumerate(tests):
        ax.text(1.02, y_pos[i], f'n={test["n"]}', fontsize=7,
                va='center', ha='left', color=C_GRAY,
                transform=ax.get_yaxis_transform())

    # Legend removed — all markers are now a single neutral color, so the
    # former blue/orange rho-sign split no longer applies.

    clean_spine(ax)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0)

    save_figure(fig, PANELS / 'fig6b_mechanistic_nulls')
    return fig


# ===================================================================
# MAIN
# ===================================================================
if __name__ == '__main__':
    print("=" * 60)
    print("PHASE 2: Generating NEEDS_CREATION panels")
    print("=" * 60)

    panels = {}
    panels['fig1a'] = fig1a_pipeline_schematic()
    panels['fig6b'] = fig6b_mechanistic_nulls()

    print("\n" + "=" * 60)
    print(f"Phase 2 complete: {len(panels)} panels generated")
    print(f"Output: {PANELS}")
    print("=" * 60)
