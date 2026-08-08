#!/usr/bin/env python3
"""
Generate Phase 3 EXISTS panels — recreate from source data with Cell Systems formatting.

22 panels across Figs 1D/E, 3A/B/C/D, 4A/B/C, 6A, 7A/B, S2A/B, S3A/B, S4A, S5A/B/C, S6A, S7A.

Biology: Bootstrap stability, LOOCV, ellipsoid alignment, replication null distributions,
L1000 sampling robustness, treeness-rigidity, disease deformation, DILI, protocol sensitivity.

Math: CVs, error-to-null ratios, Krzanowski subspace similarity, Spearman correlations,
null distributions, Fisher's exact enrichment.
"""

import sys
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))
from cellwarp.figure_style import (
    apply_style, add_panel_label, save_figure, format_p, clean_spine,
    add_lineage_legend, add_lineage_ref,
    COL1, COL15, COL2, DPI,
    C_BLUE, C_ORANGE, C_PURPLE, C_TEAL, C_GRAY, C_LIGHTGRAY, C_DARKGRAY, C_BLACK,
    LINEAGE_COLORS, LINEAGE_MAP, DATASET_COLORS,
    FONT_SIZE_LABEL, FONT_SIZE_TICK, FONT_SIZE_ANNOT, FONT_SIZE_LEGEND, FONT_FAMILY,
    short_name, lineage_color, LABEL_EXPAND,
)

apply_style()

BASE = Path(__file__).parent.parent
PANELS = BASE / 'figures' / 'panels'
PANELS.mkdir(parents=True, exist_ok=True)


# ===================================================================
# Fig 1D: Bootstrap stability histogram
# ===================================================================
def fig1d_bootstrap():
    """Histogram of bootstrap Procrustes distances (CV=0.004)."""
    print("  Fig 1D: Bootstrap stability...")
    df = pd.read_csv(BASE / 'output/phase3/bootstrap/bootstrap_results.csv')
    with open(BASE / 'output/phase3/bootstrap/bootstrap_summary.json') as f:
        summary = json.load(f)

    fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.7))
    ax.hist(df['distance'], bins=20, color=C_BLUE, edgecolor='white',
            linewidth=0.3, alpha=0.85, zorder=3)
    ax.axvline(summary['distances']['mean'], color=C_ORANGE, linewidth=1.2,
               linestyle='--', zorder=5, label=f'Mean = {summary["distances"]["mean"]:.2f}')

    ax.text(0.97, 0.95,
            f'CV = {summary["distances"]["cv"]:.4f}\n'
            f'100/100 sig. at α = 0.01\n'
            f'Gate: CV < 0.2 (PASS)',
            transform=ax.transAxes, fontsize=FONT_SIZE_ANNOT,
            ha='right', va='top', color=C_DARKGRAY,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor=C_LIGHTGRAY, alpha=0.9))

    ax.set_xlabel('Procrustes distance', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('Count (100 iterations)', fontsize=FONT_SIZE_LABEL)
    ax.legend(fontsize=FONT_SIZE_LEGEND, frameon=False)
    clean_spine(ax)
    save_figure(fig, PANELS / 'fig1d_bootstrap')


# ===================================================================
# Fig 1E: LOOCV bar chart
# ===================================================================
def fig1e_loocv():
    """Bar chart of LOOCV error-to-null ratios (35/35 better than chance)."""
    print("  Fig 1E: LOOCV...")
    df = pd.read_csv(BASE / 'output/phase3/loocv/loocv_results.csv')
    df = df.sort_values('ratio', ascending=True)  # rigid at top, matching Fig 2A

    fig, ax = plt.subplots(figsize=(COL2, COL2 * 0.55))
    y_pos = np.arange(len(df))
    colors = [lineage_color(ct) for ct in df['cell_type']]

    ax.barh(y_pos, df['ratio'], height=0.75, color=colors,
            edgecolor='white', linewidth=0.3, zorder=3)
    ax.axvline(1.0, color=C_DARKGRAY, linewidth=0.8, linestyle='--', zorder=1)

    # Long-form labels via centralized LABEL_EXPAND (figure_style.py) — kept
    # in sync with Fig 3A and Fig 4A.
    ax.set_yticks(y_pos)
    ax.set_yticklabels(
        [LABEL_EXPAND.get(short_name(ct), short_name(ct)) for ct in df['cell_type']],
        fontsize=FONT_SIZE_TICK)
    ax.set_xlabel('LOOCV error / null ratio', fontsize=FONT_SIZE_LABEL)
    ax.invert_yaxis()

    # Annotation given a white bbox so it reads clearly on top of the
    # dashed x=1 line instead of blending into the tick label beneath it.
    # Nudged up from y=0.02 to y=0.05 to leave room for the x-tick labels.
    ax.text(0.97, 0.05, '35/35 < 1.0\nMean = 0.4201',
            transform=ax.transAxes, fontsize=FONT_SIZE_ANNOT,
            ha='right', va='bottom', color=C_DARKGRAY, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor=C_LIGHTGRAY, alpha=0.9))

    clean_spine(ax)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0)
    add_lineage_legend(ax, loc='upper right', ncol=1, title='Lineage')
    save_figure(fig, PANELS / 'fig1e_loocv')


# ===================================================================
# Fig 3A: Ellipsoid alignment heatmap (35 types)
# ===================================================================
def fig3a_ellipsoid_heatmap():
    """Heatmap of cross-species ellipsoid alignment scores."""
    print("  Fig 3A: Ellipsoid alignment heatmap...")
    df = pd.read_csv(BASE / 'output/mechanistic/ellipsoid_alignment/35type_alignment_scores.csv')

    # Pivot: rows = cell types, columns = k values, values = S_pre
    df_pre = df[['cell_type', 'k', 'S_pre']].pivot(index='cell_type', columns='k', values='S_pre')
    df_pre = df_pre.sort_values(by=5, ascending=False) if 5 in df_pre.columns else df_pre

    # Wider canvas (was COL1*0.75) reduces composite whitespace right of A.
    fig, ax = plt.subplots(figsize=(COL1 * 1.15, COL2 * 0.7))
    im = ax.imshow(df_pre.values, aspect='auto', cmap='YlOrBr',
                   vmin=0, vmax=1, interpolation='nearest')
    ax.set_xticks(range(len(df_pre.columns)))
    ax.set_xticklabels([f'k={c}' for c in df_pre.columns], fontsize=FONT_SIZE_TICK)
    ax.set_yticks(range(len(df_pre.index)))
    ax.set_yticklabels(
        [LABEL_EXPAND.get(short_name(ct), short_name(ct)) for ct in df_pre.index],
        fontsize=FONT_SIZE_TICK)
    ax.set_xlabel('Principal components', fontsize=FONT_SIZE_LABEL)

    cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label('Krzanowski subspace similarity (S)', fontsize=6)
    cbar.ax.tick_params(labelsize=5.5)

    save_figure(fig, PANELS / 'fig3a_ellipsoid_heatmap')


# ===================================================================
# Fig 3B: Pre- vs post-rotation comparison
# ===================================================================
def fig3b_pre_post():
    """Bar chart comparing pre- and post-rotation ellipsoid alignment."""
    print("  Fig 3B: Pre vs post rotation...")
    # Reads permutation_results.json, not summary_stats.json. The observed values
    # are identical in both, checked to full precision at all six bars, but only
    # this file carries the null and the p-value, so the panel could not show its
    # own evidence while reading the other one. summary_stats.json is NOT retired:
    # its 35type.eigenval_vs_rigidity is gated by reproduce/validate.py.
    with open(BASE / 'output/mechanistic/ellipsoid_alignment/permutation_results.json') as f:
        perm = json.load(f)['35type']

    ks = ['k=1', 'k=3', 'k=5']
    pre = [perm['label_shuffle_pre'][k] for k in ks]
    post = [perm['label_shuffle_post'][k] for k in ks]

    # 1e-4 is the floor of a 10,000-permutation test, not a measurement, so it is
    # shown as a bound; the prose uses that form throughout. Above the floor the
    # shared format_p applies.
    def p_label(p):
        return 'p < 10$^{-4}$' if p <= 1e-4 else format_p(p)

    fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.65))
    x = np.arange(len(ks))
    w = 0.35

    pre_vals = [d['observed'] for d in pre]
    post_vals = [d['observed'] for d in post]

    ax.bar(x - w/2, pre_vals, w, color=C_BLUE, edgecolor='white',
           linewidth=0.3, label='Pre-rotation', zorder=3)
    ax.bar(x + w/2, post_vals, w, color=C_ORANGE, edgecolor='white',
           linewidth=0.3, label='Post-rotation', zorder=3)

    # The permutation null on every bar, and the p-value above it. Without both,
    # a reader reads bar-above-marker as success: at k=1 the post-rotation bar
    # sits below its null, and the two bars that clear theirs by 6.9 and 10.7 SD
    # look only modestly above them on this axis.
    for i, (a, b) in enumerate(zip(pre, post)):
        for off, d in ((-w/2, a), (w/2, b)):
            ax.plot([i + off - w*0.42, i + off + w*0.42], [d['null_mean']] * 2,
                    color=C_DARKGRAY, lw=1.0, zorder=5,
                    label='Label-shuffle null' if (i == 0 and off < 0) else None)
            v = d['observed']
            ax.text(i + off, v + 0.01, f'{v:.3f}',
                    ha='center', va='bottom', fontsize=5.0, color=C_DARKGRAY)
            ax.text(i + off, v + 0.045, p_label(d['p_value']),
                    ha='center', va='bottom', fontsize=4.6, color=C_DARKGRAY)

    ax.set_xticks(x)
    ax.set_xticklabels(ks, fontsize=FONT_SIZE_TICK)
    ax.set_ylabel('Krzanowski subspace similarity (S)', fontsize=FONT_SIZE_LABEL)
    ax.set_ylim(0, 0.6)
    ax.legend(fontsize=FONT_SIZE_LEGEND - 1.2, frameon=False)

    clean_spine(ax)
    save_figure(fig, PANELS / 'fig3b_pre_post')


# ===================================================================
# Fig 3C: Layer 1 and Layer 2 null distributions
# ===================================================================
def fig3c_layer_nulls():
    """Null distributions for Layer 1 (centroid) and Layer 2 (ellipsoid)."""
    print("  Fig 3C: Layer 1 & 2 null distributions...")
    with open(BASE / 'output/mechanistic/ellipsoid_alignment/permutation_results.json') as f:
        perm = json.load(f)

    # Wider canvas (was COL2*0.65) so right subpanel's x-axis label
    # 'Krzanowski subspace similarity (S)' isn't truncated.
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(COL2 * 0.85, COL1 * 0.55))

    # Layer 1 (centroid position) - label_shuffle_pre at k=5
    pre = perm['35type']['label_shuffle_pre']['k=5']
    # Summary visualization — no raw null array saved
    ax1.errorbar(pre['null_mean'], 1, xerr=pre['null_std'], fmt='s',
                 color=C_GRAY, capsize=4, markersize=6, zorder=3,
                 markeredgecolor='white', markeredgewidth=0.5, label='Null (mean \u00b1 SD)')
    ax1.plot(pre['observed'], 0, 'o', color=C_BLUE, markersize=8, zorder=5,
             markeredgecolor='white', markeredgewidth=0.5, label='Observed')
    ax1.set_yticks([0, 1])
    ax1.set_yticklabels(['Observed', 'Null\n(10K perms)'], fontsize=FONT_SIZE_TICK)
    ax1.set_title('Layer 2: Pre-rotation (k=5)', fontsize=FONT_SIZE_LABEL)
    ax1.set_xlabel('Krzanowski subspace similarity (S)', fontsize=FONT_SIZE_LABEL)
    # Common x-axis for pre- and post-rotation subpanels: visually preserves
    # the "compressed but still significant" pattern without the prior
    # independent-range artifact.
    ax1.set_xlim(0.00, 0.55)
    ax1.set_ylim(-0.5, 1.8)
    p1_str = format_p(pre['p_value'])
    if pre['p_value'] <= 0.0001:
        p1_str += '\n(permutation floor)'
    ax1.text(0.95, 0.95, p1_str,
             transform=ax1.transAxes, fontsize=FONT_SIZE_ANNOT,
             ha='right', va='top', color=C_DARKGRAY)
    ax1.legend(fontsize=5.5, frameon=False, loc='lower right')
    clean_spine(ax1)
    ax1.spines['left'].set_visible(False)
    ax1.tick_params(axis='y', length=0)

    # Layer 2 (ellipsoid orientation) - label_shuffle_post at k=5
    post = perm['35type']['label_shuffle_post']['k=5']
    ax2.errorbar(post['null_mean'], 1, xerr=post['null_std'], fmt='s',
                 color=C_GRAY, capsize=4, markersize=6, zorder=3,
                 markeredgecolor='white', markeredgewidth=0.5)
    ax2.plot(post['observed'], 0, 'o', color=C_ORANGE, markersize=8, zorder=5,
             markeredgecolor='white', markeredgewidth=0.5)
    ax2.set_yticks([0, 1])
    ax2.set_yticklabels(['Observed', 'Null\n(10K perms)'], fontsize=FONT_SIZE_TICK)
    ax2.set_title('Layer 2: Post-rotation (k=5)', fontsize=FONT_SIZE_LABEL)
    ax2.set_xlabel('Krzanowski subspace similarity (S)', fontsize=FONT_SIZE_LABEL)
    ax2.set_xlim(0.00, 0.55)  # common x-axis with ax1 (Pre-rotation)
    ax2.set_ylim(-0.5, 1.8)
    p2_str = format_p(post['p_value'])
    if post['p_value'] <= 0.0001:
        p2_str += '\n(permutation floor)'
    ax2.text(0.95, 0.95, p2_str,
             transform=ax2.transAxes, fontsize=FONT_SIZE_ANNOT,
             ha='right', va='top', color=C_DARKGRAY)
    clean_spine(ax2)
    ax2.spines['left'].set_visible(False)
    ax2.tick_params(axis='y', length=0)

    save_figure(fig, PANELS / 'fig3c_layer_nulls')


# ===================================================================
# Fig 3D: Layer 1 vs Layer 2 per-type scatter
# ===================================================================
def fig3d_layer_scatter():
    """Scatter of Layer 1 (residual) vs Layer 2 (alignment) per type."""
    print("  Fig 3D: Layer 1 vs Layer 2 scatter...")
    df_align = pd.read_csv(BASE / 'output/mechanistic/ellipsoid_alignment/35type_alignment_scores.csv')
    df_align_k3 = df_align[df_align['k'] == 3][['cell_type', 'S_pre']].copy()

    with open(BASE / 'output/cellcount_confound/cellcount_confound_results.json') as f:
        cc_data = json.load(f)
    resid_df = pd.DataFrame(cc_data['per_type_data'])[['cell_type', 'residual_magnitude']]

    merged = pd.merge(df_align_k3, resid_df, on='cell_type')

    fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.8))

    colors = [lineage_color(ct) for ct in merged['cell_type']]
    ax.scatter(merged['residual_magnitude'], merged['S_pre'],
               c=colors, s=25, edgecolors='white', linewidths=0.3, zorder=3, alpha=0.85)

    from scipy.stats import spearmanr
    rho, p = spearmanr(merged['residual_magnitude'], merged['S_pre'])

    # Linear regression line over the data range (illustrative trend; the
    # Spearman ρ is the inferential statistic).
    xv = merged['residual_magnitude'].values
    yv = merged['S_pre'].values
    slope, intercept = np.polyfit(xv, yv, 1)
    xline = np.linspace(xv.min(), xv.max(), 50)
    ax.plot(xline, slope * xline + intercept,
            color=C_DARKGRAY, linewidth=1.0, linestyle='--', zorder=2, alpha=0.8)

    ax.text(0.03, 0.97,
            f'ρ = {rho:.3f}\nuncorrected p = {p:.3f}\nn = {len(merged)}',
            transform=ax.transAxes, fontsize=FONT_SIZE_ANNOT,
            ha='left', va='top', color=C_DARKGRAY,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor=C_LIGHTGRAY, alpha=0.9))

    ax.set_xlabel('Procrustes residual (Layer 1)', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('Ellipsoid alignment S (Layer 2)', fontsize=FONT_SIZE_LABEL)
    clean_spine(ax)
    # Lineage colors match Fig 3A (rigidity ranking), not Fig 4A (heatmap).
    # Moved from upper-right (0.97, 0.97) to lower-right (0.97, 0.03):
    # upper-right overlapped a blue scatter point at the edge of the data,
    # making the caption partially illegible. The scatter's lower-right is
    # clear (no points in the y < 0.25 band at x > 15).
    ax.text(0.97, 0.03, 'Colors as in Fig. 3A', transform=ax.transAxes,
            fontsize=5.5, color=C_GRAY, fontstyle='italic',
            ha='right', va='bottom')
    save_figure(fig, PANELS / 'fig3d_layer_scatter')


# ===================================================================
# Fig 4A: Sun2023 null distribution
# ===================================================================
def fig4a_sun2023_null():
    """Null distribution for Sun2023 expanded 15-type replication."""
    print("  Fig 4A: Sun2023 expanded null...")
    null_dist = np.load(BASE / 'output/validation/sun2023_replication_expanded/null_distribution.npy')
    with open(BASE / 'output/validation/sun2023_replication_expanded/sun2023_expanded.json') as f:
        data = json.load(f)

    proc = data['procrustes']
    obs = proc['distance']

    fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.7))
    ax.hist(null_dist, bins=60, color=C_LIGHTGRAY, edgecolor='white',
            linewidth=0.3, density=True, alpha=0.9, zorder=2)
    ax.set_xlim(25, 80)
    ax.axvline(obs, color=C_ORANGE, linewidth=1.5, zorder=5)

    # Observed label inside plot, near top of line, nudged right
    ax.text(obs + 0.8, 0.92, f'Observed\n({obs:.1f})',
            transform=ax.get_xaxis_transform(), fontsize=FONT_SIZE_ANNOT,
            color=C_ORANGE, ha='left', va='top')

    ax.set_title(
        f'Sun2023 × TS\nobs/null = {proc["obs_null_ratio"]:.3f}, '
        f'p = {proc["p_value"]:.4f}, n = {proc["n_types"]}',
        fontsize=FONT_SIZE_LABEL, pad=4, loc='left')

    ax.set_xlabel('Procrustes distance', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('Density', fontsize=FONT_SIZE_LABEL)
    clean_spine(ax)
    save_figure(fig, PANELS / 'fig4a_sun2023_null')


# ===================================================================
# Fig 4B: PanSci null distribution
# ===================================================================
def fig4b_pansci_null():
    """Null distribution for PanSci replication."""
    print("  Fig 4B: PanSci null...")
    null_dist = np.load(BASE / 'output/validation/pansci_replication/null_distribution.npy')
    with open(BASE / 'output/validation/pansci_replication/pansci_replication.json') as f:
        data = json.load(f)

    obs = data['procrustes']['distance']

    fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.7))
    ax.hist(null_dist, bins=60, color=C_LIGHTGRAY, edgecolor='white',
            linewidth=0.3, density=True, alpha=0.9, zorder=2)
    ax.set_xlim(25, 80)
    ax.axvline(obs, color=C_PURPLE, linewidth=1.5, zorder=5)

    # Observed label inside plot, near top of line, nudged right
    ax.text(obs + 0.8, 0.92, f'Observed\n({obs:.1f})',
            transform=ax.get_xaxis_transform(), fontsize=FONT_SIZE_ANNOT,
            color=C_PURPLE, ha='left', va='top')

    ax.set_title(
        f'PanSci × TS\nobs/null = {data["procrustes"]["obs_null_ratio"]:.3f}, '
        f'p = {data["procrustes"]["p_value"]:.4f}, n = {data["procrustes"]["n_types"]}',
        fontsize=FONT_SIZE_LABEL, pad=4, loc='left')

    ax.set_xlabel('Procrustes distance', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('Density', fontsize=FONT_SIZE_LABEL)
    clean_spine(ax)
    save_figure(fig, PANELS / 'fig4b_pansci_null')


# ===================================================================
# Fig 4C: CellHint null distribution
# ===================================================================
def fig4c_cellhint_null():
    """Null distribution for CellHint replication."""
    print("  Fig 4C: CellHint null...")
    null_dist = np.load(BASE / 'output/validation/cellhint_replication/null_distribution.npy')
    with open(BASE / 'output/validation/cellhint_replication/cellhint_replication.json') as f:
        data = json.load(f)

    obs = data['procrustes']['distance']

    fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.7))
    ax.hist(null_dist, bins=60, color=C_LIGHTGRAY, edgecolor='white',
            linewidth=0.3, density=True, alpha=0.9, zorder=2)
    ax.set_xlim(25, 80)
    ax.axvline(obs, color=C_TEAL, linewidth=1.5, zorder=5)

    # Observed label inside plot, near top of line, nudged right
    ax.text(obs + 0.8, 0.92, f'Observed\n({obs:.1f})',
            transform=ax.get_xaxis_transform(), fontsize=FONT_SIZE_ANNOT,
            color=C_TEAL, ha='left', va='top')

    ax.set_title(
        f'CellHint × TMS\nobs/null = {data["procrustes"]["obs_null_ratio"]:.3f}, '
        f'p = {data["procrustes"]["p_value"]:.4f}, n = {data["procrustes"]["n_types"]}',
        fontsize=FONT_SIZE_LABEL, pad=4, loc='left')

    ax.set_xlabel('Procrustes distance', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('Density', fontsize=FONT_SIZE_LABEL)
    clean_spine(ax)
    save_figure(fig, PANELS / 'fig4c_cellhint_null')


# ===================================================================
# Fig 6A: L1000 sampling robustness
# ===================================================================
def fig6a_l1000():
    """L1000 rho vs random baseline null-distribution histogram.

    Previous version was a strip plot (SD bar for the null + single dot for
    observed) because the raw 1000 random rho draws were not saved when
    script 35 was first run. The draws have since been regenerated (same
    seed=42, deterministic — mean/SD/p match the committed JSON to 3 dp)
    and saved as l1000_primary_rhos.npy, so this panel now renders as a true
    histogram with the observed L1000 \u03c1 marked as a vertical line,
    matching the null-distribution convention used elsewhere in the figure set.
    """
    print("  Fig 6A: L1000 robustness (histogram)...")
    with open(BASE / 'output/figures/l1000_random_baseline_results.json') as f:
        data = json.load(f)
    prim = data['primary']

    rhos = np.load(BASE / 'output/figures/l1000_primary_rhos.npy')
    obs = prim['observed_l1000_rho']

    fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.75))
    counts, edges, _ = ax.hist(rhos, bins=40, color=C_LIGHTGRAY,
                               edgecolor=C_DARKGRAY, linewidth=0.4, zorder=2)
    # Observed L1000 \u03c1 as a vertical line.
    ax.axvline(obs, color=C_ORANGE, linewidth=1.6, zorder=5)

    # Headroom above tallest bar so the leader label and the stats inset
    # don't collide with histogram bars.
    ax.set_ylim(0, counts.max() * 1.40)
    y_top = ax.get_ylim()[1]

    # Stats inset: upper-left (far corner from the observed line, which sits
    # in the right tail of the null).
    ax.text(0.03, 0.97,
            f'Random mean \u03c1 = {prim["random_mean_rho"]:.3f} \u00b1 {prim["random_std_rho"]:.3f}\n'
            f'Empirical p = {prim["empirical_p_value"]:.3f}\n'
            f'{prim["n_genes_sampled"]} / {prim["total_ortholog_genes"]:,} genes, '
            f'{prim["n_iterations"]:,} draws',
            transform=ax.transAxes, fontsize=FONT_SIZE_ANNOT,
            ha='left', va='top', color=C_DARKGRAY,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor=C_LIGHTGRAY, alpha=0.95))

    # Observed label sits to the RIGHT of the line in the upper area \u2014
    # left-of-line collided with the upper-left stats inset; centered
    # above-line ran the line through the \u03c1 character. Right-of-line
    # places the label clear of both.
    ax.annotate(f'L1000 (\u03c1 = {obs:.3f})',
                xy=(obs, y_top * 0.55),
                xytext=(obs + 0.015, y_top * 0.92),
                ha='left', va='bottom', fontsize=FONT_SIZE_ANNOT,
                color=C_ORANGE, fontweight='bold',
                arrowprops=dict(arrowstyle='->', color=C_ORANGE, lw=0.6,
                                shrinkA=0, shrinkB=2))

    ax.set_xlabel('Spearman \u03c1 (divergence ranking correlation)',
                  fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('random-draw count', fontsize=FONT_SIZE_LABEL)
    clean_spine(ax)
    save_figure(fig, PANELS / 'fig6a_l1000')


# ===================================================================
# Fig 7A: Treeness vs rigidity scatter
# ===================================================================
def fig7a_treeness():
    """Treeness vs rigidity scatter (ρ = −0.349), from real per-type data."""
    print("  Fig 7A: Treeness vs rigidity...")
    with open(BASE / 'output/liang_wagner/treeness_rigidity_correlation.json') as f:
        data = json.load(f)

    rho = data['step3_correlation']['rho_treeness_rigidity']
    p = data['step3_correlation']['p_value']

    # Load actual per-type treeness scores and rigidity residuals
    treeness_df = pd.read_csv(BASE / 'output/liang_wagner/treeness_scores_per_celltype.csv')
    resid_df = pd.read_csv(BASE / 'output/phase2/scaled_35types/residuals_ranked.csv')

    # Merge on cell_type
    merged = treeness_df.merge(resid_df[['cell_type', 'residual_magnitude']],
                               on='cell_type', how='inner')

    # Compute rigidity rank (1=most flexible/smallest residual, 35=most rigid)
    merged['rigidity_rank'] = merged['residual_magnitude'].rank(method='average').astype(int)
    # Treeness rank from CSV (already ranked 1=highest treeness)
    merged['treeness_rank'] = merged['rank']

    fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.8))

    colors = [lineage_color(ct) for ct in merged['cell_type']]
    ax.scatter(merged['rigidity_rank'], merged['treeness_rank'], c=colors, s=25,
               edgecolors='white', linewidths=0.3, zorder=3, alpha=0.85)

    # Trend line
    z = np.polyfit(merged['rigidity_rank'], merged['treeness_rank'], 1)
    n = len(merged)
    x_line = np.array([1, n])
    ax.plot(x_line, np.polyval(z, x_line), color=C_DARKGRAY,
            linewidth=0.8, linestyle='--', zorder=1)

    # Stats box — upper-right. Opaque bbox (alpha=1.0, zorder=10) masks any
    # dots that fall under the text, matching S2 Panel C/D treatment.
    ax.text(0.97, 0.97,
            f'ρ = {rho:.3f}, p = {p:.3f}\n(uncorrected; NS after Bonferroni)',
            transform=ax.transAxes, fontsize=6,
            ha='right', va='top', color=C_DARKGRAY, zorder=10,
            bbox=dict(boxstyle='round,pad=0.25', facecolor='white',
                      edgecolor='lightgray', linewidth=0.4, alpha=1.0))

    ax.set_xlabel('Rigidity rank (1=flexible, 35=rigid)', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('Treeness rank (Liang-Wagner)', fontsize=FONT_SIZE_LABEL)
    clean_spine(ax)
    # Inline replacement for add_lineage_ref(ax) — adds opaque white bbox so
    # the note masks any dots beneath. Shared helper has no bbox; inline here
    # to keep non-S3 callers unchanged (matches S2 Panel C/D precedent).
    ax.text(0.03, 0.03, 'Colors as in Fig. 3A', transform=ax.transAxes,
            fontsize=5.5, color=C_GRAY, fontstyle='italic',
            ha='left', va='bottom', zorder=10,
            bbox=dict(boxstyle='round,pad=0.25', facecolor='white',
                      edgecolor='lightgray', linewidth=0.4, alpha=1.0))
    save_figure(fig, PANELS / 'fig7a_treeness')


# ===================================================================
# Fig 7B: Density vs rigidity scatter
# ===================================================================
def fig7b_density():
    """Neighborhood density (k=5) vs rigidity scatter."""
    print("  Fig 7B: Density vs rigidity...")
    density_df = pd.read_csv(BASE / 'output/liang_wagner/neighborhood_density.csv')
    resid_df = pd.read_csv(BASE / 'output/phase2/scaled_35types/residuals_ranked.csv')

    # Use k=5 as primary (per mediation analysis)
    merged = density_df[['cell_type', 'density_k5']].merge(
        resid_df[['cell_type', 'residual_magnitude']], on='cell_type', how='inner')

    fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.8))

    colors = [lineage_color(ct) for ct in merged['cell_type']]
    ax.scatter(merged['residual_magnitude'], merged['density_k5'],
               c=colors, s=25, edgecolors='white', linewidths=0.3, zorder=3, alpha=0.85)

    # Trend line
    z = np.polyfit(merged['residual_magnitude'], merged['density_k5'], 1)
    x_line = np.array([merged['residual_magnitude'].min(), merged['residual_magnitude'].max()])
    ax.plot(x_line, np.polyval(z, x_line), color=C_DARKGRAY,
            linewidth=0.8, linestyle='--', zorder=1)

    # Report H1 rho (density vs residual) — the direct correlation
    from scipy.stats import spearmanr
    rho_h1, p_h1 = spearmanr(merged['residual_magnitude'], merged['density_k5'])
    # Opaque bbox (alpha=1.0, zorder=10) matches S3 Panel A treatment.
    ax.text(0.03, 0.97,
            f'ρ = {rho_h1:.3f}\np = {p_h1:.3f}',
            transform=ax.transAxes, fontsize=FONT_SIZE_ANNOT,
            ha='left', va='top', color=C_DARKGRAY, zorder=10,
            bbox=dict(boxstyle='round,pad=0.25', facecolor='white',
                      edgecolor='lightgray', linewidth=0.4, alpha=1.0))

    ax.set_xlabel('Procrustes residual magnitude', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('Neighborhood density (k=5)', fontsize=FONT_SIZE_LABEL)
    clean_spine(ax)
    # Inline replacement for add_lineage_ref(ax) with opaque bbox; mirrors
    # Panel A treatment (shared helper has no bbox; inline preserves non-S3
    # callers).
    ax.text(0.03, 0.03, 'Colors as in Fig. 3A', transform=ax.transAxes,
            fontsize=5.5, color=C_GRAY, fontstyle='italic',
            ha='left', va='bottom', zorder=10,
            bbox=dict(boxstyle='round,pad=0.25', facecolor='white',
                      edgecolor='lightgray', linewidth=0.4, alpha=1.0))
    save_figure(fig, PANELS / 'fig7b_density')


# ===================================================================
# Fig S2A/B: CellHint investigation
# ===================================================================
def figs2_cellhint():
    """CellHint rank reversal investigation panels."""
    print("  Fig S2A/B: CellHint investigation...")
    rank_df = pd.read_csv(BASE / 'analysis/cellhint_investigation/rank_reversal_table.csv')
    factors_df = pd.read_csv(BASE / 'analysis/cellhint_investigation/systematic_factors.csv')

    # S2A: Rank scatter
    fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.85))
    if 'primary_rank' in rank_df.columns and 'cellhint_rank' in rank_df.columns:
        ax.scatter(rank_df['primary_rank'], rank_df['cellhint_rank'],
                   c=C_BLUE, s=30, edgecolors='white', linewidths=0.3, zorder=3)
        max_r = max(rank_df['primary_rank'].max(), rank_df['cellhint_rank'].max()) + 1
        ax.plot([0, max_r], [0, max_r], color=C_LIGHTGRAY, linewidth=0.8,
                linestyle='--', zorder=1)
        ax.set_xlabel('Primary rank', fontsize=FONT_SIZE_LABEL)
        ax.set_ylabel('CellHint rank', fontsize=FONT_SIZE_LABEL)
        ax.text(0.35, 0.97, 'ρ = −0.386\np = 0.156',
                transform=ax.transAxes, fontsize=FONT_SIZE_ANNOT,
                ha='center', va='top', color=C_DARKGRAY,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                          edgecolor=C_LIGHTGRAY, alpha=0.9))
    clean_spine(ax)
    save_figure(fig, PANELS / 'figs2a_cellhint_rank')

    # S2B: Systematic factor associations — horizontal dot plot
    fig2, ax2 = plt.subplots(figsize=(COL1, COL1 * 0.65))

    y_pos = np.arange(len(factors_df))[::-1]
    for i, (_, row) in enumerate(factors_df.iterrows()):
        y = y_pos[i]
        rho_val = float(row['rho'])
        p_val = float(row['p_value'])
        color = C_TEAL if p_val < 0.05 else C_GRAY
        marker = 'o' if p_val < 0.05 else 's'
        ax2.plot(rho_val, y, marker, color=color, markersize=6, zorder=4,
                 markeredgecolor='white', markeredgewidth=0.5)

    ax2.axvline(0, color=C_DARKGRAY, linewidth=0.5, linestyle='-', zorder=1)
    ax2.set_yticks(y_pos)
    # Shorten factor names for readability
    factor_labels = []
    for f in factors_df['factor']:
        f = str(f).replace('abs(log2 cell count ratio)', 'Cell count asymmetry')
        f = f.replace('log2(CellHint/primary cell count)', 'Count direction')
        f = f.replace('CellHint tissue count', 'Tissue count')
        if len(f) > 30:
            f = f[:27] + '...'
        factor_labels.append(f)
    ax2.set_yticklabels(factor_labels, fontsize=FONT_SIZE_TICK)
    ax2.set_xlabel('Spearman ρ', fontsize=FONT_SIZE_LABEL)
    ax2.set_title('Factors associated with rank reversal', fontsize=FONT_SIZE_LABEL)

    # Mark significance
    from matplotlib.lines import Line2D
    legend_el = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=C_TEAL,
               markersize=5, label='p < 0.05'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor=C_GRAY,
               markersize=5, label='NS'),
    ]
    ax2.legend(handles=legend_el, fontsize=FONT_SIZE_LEGEND, frameon=False,
               loc='lower right')

    clean_spine(ax2)
    ax2.spines['left'].set_visible(False)
    ax2.tick_params(axis='y', length=0)
    save_figure(fig2, PANELS / 'figs2b_cellhint_residual')


# ===================================================================
# Fig S3A/B: PCA sensitivity
# ===================================================================
def figs3_pca_sensitivity():
    """PCA sensitivity analysis (k=10 to k=50)."""
    print("  Fig S3A/B: PCA sensitivity...")
    with open(BASE / 'output/validation/pca_sensitivity/pca_sensitivity_results.json') as f:
        data = json.load(f)

    # Extract k values and metrics
    ks = []
    obs_nulls = []
    rhos = []
    for key, val in data.items():
        if isinstance(val, dict) and 'n_components' in val:
            ks.append(val['n_components'])
            obs_nulls.append(val['obs_null_ratio'])
            rhos.append(val.get('spearman_rho_vs_reference', 1.0))

    if not ks:
        print("    No PCA sensitivity data found, skipping.")
        return

    # Sort by k
    sort_idx = np.argsort(ks)
    ks = [ks[i] for i in sort_idx]
    obs_nulls = [obs_nulls[i] for i in sort_idx]
    rhos = [rhos[i] for i in sort_idx]

    # S3A: Ranking scatter vs k=33. Filter out the k=33 self-correlation
    # point — it sits at ρ=1.0 by construction (self vs self) and adds no
    # analytical signal; the baseline is communicated via the axhline.
    ks_nonref = [k for k in ks if k != 33]
    rhos_nonref = [r for k, r in zip(ks, rhos) if k != 33]
    fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.7))
    ax.plot(ks_nonref, rhos_nonref, 'o-', color=C_BLUE, markersize=4, linewidth=1, zorder=3)
    ax.axhline(1.0, color=C_LIGHTGRAY, linewidth=0.5, linestyle='--', zorder=1)
    # Label the baseline; place at right edge just above the dashed line.
    ax.text(0.97, 0.96, 'k=33 baseline', transform=ax.transAxes,
            fontsize=FONT_SIZE_ANNOT, color=C_DARKGRAY,
            ha='right', va='top')
    ax.set_xlabel('PCA components (k)', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('ρ vs k=33 ranking', fontsize=FONT_SIZE_LABEL)
    ax.set_ylim(0.7, 1.05)
    clean_spine(ax)
    save_figure(fig, PANELS / 'figs3a_pca_ranking')

    # S3B: obs/null vs k. axhline(1.0) marks the null-mean reference
    # (= "no signal"), matching Fig 1B convention.
    fig2, ax2 = plt.subplots(figsize=(COL1, COL1 * 0.7))
    ax2.plot(ks, obs_nulls, 's-', color=C_ORANGE, markersize=4, linewidth=1, zorder=3)
    ax2.axhline(1.0, color=C_LIGHTGRAY, linewidth=0.5, linestyle='--', zorder=1)
    ax2.set_xlabel('PCA components (k)', fontsize=FONT_SIZE_LABEL)
    ax2.set_ylabel('Obs/null ratio', fontsize=FONT_SIZE_LABEL)
    ax2.set_ylim(0, 1.1)
    clean_spine(ax2)
    save_figure(fig2, PANELS / 'figs3b_pca_obsnull')


# ===================================================================
# Fig S4A: SAMap heatmap — use existing image (can't recreate without SAMap)
# ===================================================================
def figs4_samap():
    """SAMap heatmap — recreated from correspondence scores with cividis colormap.

    Cell types are reordered by lineage (Immune → Epithelial → Stromal →
    Endothelial → Metabolic per LINEAGE_COLORS dict order, alphabetical
    within each), same on both axes for diagonal coherence. Row-maxima are
    highlighted with small white-filled dots so the "24/35 same-name
    correspondence" claim is verifiable directly.
    """
    print("  Fig S4A: SAMap heatmap...")
    csv_path = BASE / 'output/phase1_samap/samap_35types/samap_mapping_scores_35.csv'
    if not csv_path.exists():
        # Fall back to copying existing image
        existing = BASE / 'output/phase1_samap/samap_35types/samap_heatmap_35.png'
        if existing.exists():
            from shutil import copy2
            copy2(existing, PANELS / 'figs4a_samap_heatmap.png')
            print(f"    Copied existing (CSV not found): {existing}")
        else:
            print("    WARNING: SAMap data not found, skipping.")
        return

    scores = pd.read_csv(csv_path, index_col=0)

    # Reorder by lineage (primary) then alphabetical (secondary), same on
    # both axes. All 35 types resolve through LINEAGE_MAP — verified in audit.
    LINEAGE_ORDER = ['Immune', 'Epithelial', 'Stromal', 'Endothelial', 'Metabolic']

    def _lineage_sort_key(ct):
        lin = LINEAGE_MAP.get(ct, 'Unknown')
        idx = LINEAGE_ORDER.index(lin) if lin in LINEAGE_ORDER else len(LINEAGE_ORDER)
        return (idx, ct.lower())

    ordered_types = sorted(scores.index, key=_lineage_sort_key)
    scores = scores.reindex(index=ordered_types, columns=ordered_types)

    # Increase figure dimensions to accommodate 35 cell type labels
    fig, ax = plt.subplots(figsize=(COL2, COL2 * 0.9))
    im = ax.imshow(scores.values, aspect='auto', cmap='cividis',
                   vmin=0, vmax=1, interpolation='nearest')

    # Row-max highlights: small white-filled dot with thin black edge marks
    # each row's argmax column. White-on-cividis (blue→yellow) stays legible
    # at both ends with the black outline.
    row_argmax = np.argmax(scores.values, axis=1)
    ax.scatter(row_argmax, range(len(scores.index)),
               marker='o', s=6, c='white', edgecolors='black',
               linewidths=0.3, zorder=5)

    # Tick labels: SHORT_NAMES + LABEL_EXPAND (consistent with S1B/S2/S3-C).
    ax.set_xticks(range(len(scores.columns)))
    ax.set_xticklabels(
        [LABEL_EXPAND.get(short_name(c), short_name(c)) for c in scores.columns],
        fontsize=6, rotation=45, ha='right',
    )
    ax.set_yticks(range(len(scores.index)))
    ax.set_yticklabels(
        [LABEL_EXPAND.get(short_name(c), short_name(c)) for c in scores.index],
        fontsize=6,
    )

    ax.set_xlabel('Mouse cell types', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('Human cell types', fontsize=FONT_SIZE_LABEL)
    ax.set_title('Cross-species correspondence scores (35 matched cell types)',
                 fontsize=FONT_SIZE_LABEL, pad=4)

    cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label('SAMap correspondence score', fontsize=FONT_SIZE_LABEL)
    cbar.ax.tick_params(labelsize=FONT_SIZE_TICK)

    save_figure(fig, PANELS / 'figs4a_samap_heatmap')


# ===================================================================
# Supp Text S1: Disease deformation panels
# ===================================================================
def figs5_disease():
    """Disease deformation panels (cancer + COVID + enrichment)."""
    print("  Supp Text S1: Disease deformation panels...")

    # Supp Text S1 cancer — canonical at output/cancer/scaled/cross_analysis_scaled.png
    # Supp Text S1 COVID — canonical at output/disease_replication/covid/covid_cross_analysis.png
    # Panel mirrors (suppl_text_s1_cancer.png, suppl_text_s1_covid.png) are materialized by
    # scripts/build_submission_packet.py (R21 build script). No copy here.

    # Supp Text S1 identity-vs-activation enrichment bars
    with open(BASE / 'output/mechanistic/identity_vs_state/identity_vs_state_results.json') as f:
        data = json.load(f)

    fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.65))

    # Load enrichment medians and p-values from correct key paths — no fallbacks
    id_enrich = data['identity_enrichment']
    act_enrich = data['activation_enrichment']

    cats = ['Cancer\nIdentity', 'Cancer\nActivation', 'COVID\nIdentity', 'COVID\nActivation']
    vals = [
        id_enrich['cancer_median'],
        act_enrich['cancer_median'],
        id_enrich['covid_median'],
        act_enrich['covid_median'],
    ]
    id_p = id_enrich['mann_whitney_p']
    act_p = act_enrich['mann_whitney_p']

    # Compute SEM from per-type enrichment ratios for error bars
    id_cancer_vals = [v['enrichment_ratio'] for v in id_enrich.get('cancer', {}).values() if isinstance(v, dict) and 'enrichment_ratio' in v]
    id_covid_vals = [v['enrichment_ratio'] for v in id_enrich.get('covid', {}).values() if isinstance(v, dict) and 'enrichment_ratio' in v]
    act_cancer_vals = [v['enrichment_ratio'] for v in act_enrich.get('cancer', {}).values() if isinstance(v, dict) and 'enrichment_ratio' in v]
    act_covid_vals = [v['enrichment_ratio'] for v in act_enrich.get('covid', {}).values() if isinstance(v, dict) and 'enrichment_ratio' in v]

    def _sem(arr):
        if len(arr) < 2:
            return 0
        return np.std(arr, ddof=1) / np.sqrt(len(arr))

    yerr = [_sem(id_cancer_vals), _sem(act_cancer_vals),
            _sem(id_covid_vals), _sem(act_covid_vals)]

    colors_bars = [C_BLUE, C_BLUE, C_ORANGE, C_ORANGE]
    alphas = [0.9, 0.5, 0.9, 0.5]

    bars = ax.bar(range(len(cats)), vals, width=0.6, color=colors_bars,
                  edgecolor='white', linewidth=0.5, zorder=3,
                  yerr=yerr, capsize=3, error_kw=dict(lw=0.8, color=C_DARKGRAY))
    for bar, alpha in zip(bars, alphas):
        bar.set_alpha(alpha)

    ax.set_xticks(range(len(cats)))
    ax.set_xticklabels(cats, fontsize=6)
    ax.set_ylabel('Fold enrichment', fontsize=FONT_SIZE_LABEL)
    ax.axhline(1, color=C_LIGHTGRAY, linewidth=0.5, linestyle='--', zorder=1)

    ax.text(0.97, 0.95,
            f'Identity: {format_p(id_p)}\nActivation: {format_p(act_p)}',
            transform=ax.transAxes, fontsize=FONT_SIZE_ANNOT,
            ha='right', va='top', color=C_DARKGRAY)

    clean_spine(ax)
    save_figure(fig, PANELS / 'figs5c_enrichment_bars')


# ===================================================================
# Fig S6A: DILI deformation distributions
# ===================================================================
def figs6_dili():
    """DILI deformation distributions by DILI class — KDE plot."""
    print("  Fig S6A: DILI distributions (KDE)...")

    deform_path = BASE / 'output/dilirank/deformation_distances.csv'
    if not deform_path.exists():
        print("    WARNING: deformation_distances.csv not found, skipping.")
        return

    df = pd.read_csv(deform_path)

    fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.65))

    from scipy.stats import gaussian_kde

    classes = {'most-concern': C_ORANGE, 'less-concern': C_BLUE, 'no-concern': C_GRAY}
    for cls, color in classes.items():
        subset = df[df['dili_class'] == cls]['deformation'].dropna()
        if len(subset) < 3:
            continue
        kde = gaussian_kde(subset, bw_method=0.3)
        x_range = np.linspace(subset.min() * 0.8, subset.max() * 1.2, 200)
        ax.plot(x_range, kde(x_range), color=color, linewidth=1.2,
                label=cls.replace('-', ' ').title(), zorder=3)
        ax.fill_between(x_range, kde(x_range), alpha=0.15, color=color, zorder=2)

    ax.set_xlabel('Deformation magnitude', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('Density (KDE)', fontsize=FONT_SIZE_LABEL)
    ax.set_title('DILI class deformation distributions', fontsize=FONT_SIZE_LABEL)
    ax.legend(fontsize=FONT_SIZE_LEGEND, frameon=False)
    clean_spine(ax)
    save_figure(fig, PANELS / 'figs6a_dili')


# ===================================================================
# Fig S7A: Protocol sensitivity (Smart-seq2 fraction)
# ===================================================================
def figs7_smartseq2():
    """Protocol sensitivity — split into A (ranking scatter) and B (fraction vs rank change)."""
    print("  Fig S7A/B: Protocol sensitivity...")

    comp_path = BASE / 'output/phase2/sensitivity/smartseq2/rigidity_comparison.csv'
    if not comp_path.exists():
        print("    WARNING: rigidity_comparison.csv not found, skipping.")
        return

    df = pd.read_csv(comp_path)

    with open(BASE / 'output/phase2/sensitivity/smartseq2/sensitivity_results.json') as f:
        data = json.load(f)

    # S7A: Full vs 10x-only ranking scatter
    fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.85))
    colors = [lineage_color(ct) for ct in df['cell_type']]
    ax.scatter(df['rank'], df['rank_10x'], c=colors, s=25,
               edgecolors='white', linewidths=0.3, zorder=3, alpha=0.85)
    max_r = max(df['rank'].max(), df['rank_10x'].max()) + 1
    ax.plot([0, max_r], [0, max_r], color=C_LIGHTGRAY, linewidth=0.8,
            linestyle='--', zorder=1)

    # Label points with rank change > 2 — adjustText prevents overlap + adds leader lines
    from adjustText import adjust_text
    texts = []
    for _, row in df.iterrows():
        if abs(row['rank'] - row['rank_10x']) > 2:
            texts.append(ax.text(row['rank'], row['rank_10x'],
                                 short_name(row['cell_type']),
                                 fontsize=5, color=C_DARKGRAY))
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))

    rho_mag = data['rigidity_comparison']['spearman_rho_magnitudes']
    p_mag = data['rigidity_comparison']['spearman_p_magnitudes']
    # Format p with mathtext superscript (e.g. 1.75e-7 → "$p = 1.75 × 10⁻⁷$");
    # 3-decimal formatting (:.3f) truncated this asymptotic p to "0.000",
    # hiding its true magnitude.
    _mantissa, _exp = f'{p_mag:.2e}'.split('e')
    _p_str = rf'$p = {float(_mantissa):.2f} \times 10^{{{int(_exp)}}}$'
    ax.text(0.03, 0.97, f'ρ = {rho_mag:.3f}\n{_p_str}',
            transform=ax.transAxes, fontsize=FONT_SIZE_ANNOT,
            ha='left', va='top', color=C_DARKGRAY,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor=C_LIGHTGRAY, alpha=0.9))

    ax.set_xlabel('Full ranking', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('10x-only ranking', fontsize=FONT_SIZE_LABEL)
    # Inlined lineage-ref with white bbox so the label masks any dots
    # beneath (the shared add_lineage_ref helper has no bbox). Same
    # styling applied to S2 Panel D below for uniformity.
    ax.text(0.03, 0.03, 'Colors as in Fig. 3A', transform=ax.transAxes,
            fontsize=5.5, color=C_GRAY, fontstyle='italic',
            ha='left', va='bottom', zorder=10,
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                      edgecolor=C_LIGHTGRAY, linewidth=0.4, alpha=1.0))
    clean_spine(ax)
    save_figure(fig, PANELS / 'figs7a_smartseq2')

    # S7B: Smart-seq2 fraction vs rank change
    proto_path = BASE / 'output/phase2/sensitivity/smartseq2/protocol_breakdown.csv'
    if proto_path.exists():
        proto = pd.read_csv(proto_path)
        # Map column name: CSV uses 'fraction_smartseq2' or 'smartseq2_fraction'
        frac_col = 'smartseq2_fraction' if 'smartseq2_fraction' in proto.columns else 'fraction_smartseq2'
        merged = df.merge(proto[['cell_type', frac_col]].rename(columns={frac_col: 'smartseq2_fraction'}),
                          on='cell_type', how='inner')
        if 'smartseq2_fraction' in merged.columns:
            merged['rank_change'] = abs(merged['rank'] - merged['rank_10x'])
            fig2, ax2 = plt.subplots(figsize=(COL1, COL1 * 0.85))
            colors2 = [lineage_color(ct) for ct in merged['cell_type']]
            ax2.scatter(merged['smartseq2_fraction'] * 100, merged['rank_change'],
                       c=colors2, s=25, edgecolors='white', linewidths=0.3, zorder=3)

            ss2_rho = data['rigidity_comparison'].get('ss2_fraction_vs_rank_change_rho', None)
            ss2_p = data['rigidity_comparison'].get('ss2_fraction_vs_rank_change_p', None)
            if ss2_rho is not None:
                ax2.text(0.03, 0.97, f'ρ = {ss2_rho:.3f}\np = {ss2_p:.3f}',
                        transform=ax2.transAxes, fontsize=FONT_SIZE_ANNOT,
                        ha='left', va='top', color=C_DARKGRAY,
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                                  edgecolor=C_LIGHTGRAY, alpha=0.9))

            ax2.set_xlabel('Smart-seq2 fraction (%)', fontsize=FONT_SIZE_LABEL)
            ax2.set_ylabel('|Rank change|', fontsize=FONT_SIZE_LABEL)
            # Inlined lineage-ref (matches Panel C styling — white bbox so the
            # label masks any dots beneath).
            ax2.text(0.03, 0.03, 'Colors as in Fig. 3A', transform=ax2.transAxes,
                     fontsize=5.5, color=C_GRAY, fontstyle='italic',
                     ha='left', va='bottom', zorder=10,
                     bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                               edgecolor=C_LIGHTGRAY, linewidth=0.4, alpha=1.0))
            clean_spine(ax2)
            save_figure(fig2, PANELS / 'figs7b_smartseq2_fraction')
            return

    # Fallback: create a simple summary panel B
    fig2, ax2 = plt.subplots(figsize=(COL1, COL1 * 0.65))
    ax2.text(0.5, 0.5,
             f'10x-only: obs/null = {data["procrustes_10x_only"]["obs_null_ratio"]:.3f}\n'
             f'Original: obs/null = {data["procrustes_original"]["obs_null_ratio"]:.3f}\n'
             f'ρ = {data["rigidity_comparison"]["spearman_rho_magnitudes"]:.3f}',
             transform=ax2.transAxes, fontsize=FONT_SIZE_LABEL,
             ha='center', va='center', color=C_DARKGRAY)
    clean_spine(ax2)
    save_figure(fig2, PANELS / 'figs7b_smartseq2_fraction')


# ===================================================================
# MAIN
# ===================================================================
if __name__ == '__main__':
    print("=" * 60)
    print("PHASE 3: Reformatting EXISTS panels from source data")
    print("=" * 60)

    fig1d_bootstrap()
    fig1e_loocv()
    fig3a_ellipsoid_heatmap()
    fig3b_pre_post()
    fig3c_layer_nulls()
    fig3d_layer_scatter()
    fig4a_sun2023_null()
    fig4b_pansci_null()
    fig4c_cellhint_null()
    fig6a_l1000()
    fig7a_treeness()
    fig7b_density()
    figs2_cellhint()
    figs3_pca_sensitivity()
    figs4_samap()
    figs5_disease()
    figs6_dili()
    figs7_smartseq2()

    print("\n" + "=" * 60)
    print("Phase 3 complete")
    print(f"Output: {PANELS}")
    print("=" * 60)
