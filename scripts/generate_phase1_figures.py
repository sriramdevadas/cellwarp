#!/usr/bin/env python3
"""
Generate Phase 1 DATA_ONLY panels for CellWarp Cell Systems submission.

13 panels from source data (JSON/NPY/CSV):
  Fig 1B: 1M null distribution histogram (C10 data)
  Fig 1C: Lineage-stratified null distribution
  Fig 2A: Rigidity ranking bar chart (ICONIC)
  Fig 2B: CellMarker enrichment
  Fig 2C: Cell count confound
  Fig 4D: Summary obs/null across 4 datasets
  Fig 4E: Human control comparison
  Fig 5A: Three-species Procrustes coherence
  Fig 5B: Sensitivity (13 RIRA-only types)
  Fig 5C: Hepatocyte rigidity across species
  Fig S1A: Independent-PCA null distribution (C1 data)
  Fig S1B: Joint vs independent PCA ranking scatter

Biology: Cross-species Procrustes analysis of cell type centroids
in transcriptomic space across human, mouse, and macaque.

Math: Permutation tests, Spearman rank correlations, fold-enrichment,
observed-to-null ratios, Fisher's exact tests.
"""

import sys
import json
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import FancyBboxPatch
from pathlib import Path
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))
from cellwarp.figure_style import (
    apply_style, add_panel_label, save_figure, format_p, clean_spine,
    add_lineage_legend, add_lineage_ref,
    COL1, COL15, COL2, DPI, MM_PER_INCH,
    C_BLUE, C_ORANGE, C_PURPLE, C_TEAL, C_GRAY, C_LIGHTGRAY, C_DARKGRAY, C_BLACK,
    LINEAGE_COLORS, LINEAGE_MAP, DATASET_COLORS,
    FONT_SIZE_LABEL, FONT_SIZE_TICK, FONT_SIZE_ANNOT, FONT_SIZE_PANEL,
    FONT_SIZE_LEGEND, FONT_FAMILY,
    short_name, lineage_color, LABEL_EXPAND,
)

apply_style()

BASE = Path(__file__).parent.parent
PANELS = BASE / 'figures' / 'panels'
PANELS.mkdir(parents=True, exist_ok=True)


# ===================================================================
# Fig 1B: 1M permutation null distribution
# ===================================================================
def _p_bound(p):
    """Render a permutation p as the tightest power-of-ten bound above it.

    A permutation p at the (k+1)/(n+1) floor is a bound, not a measurement,
    and a fixed-decimal format rounds it up into an assertion of equality:
    9.999e-05 through '%.4f' reads 'p = 0.0001' for a value that is strictly
    below 1e-4. This returns the bound the value actually satisfies, in the
    form the prose and the Fig 2 panels use (generate_phase3_figures.p_label,
    build_fig2c_bg.py). figure_style.format_p is left alone: it is shared with
    panels outside this figure.
    """
    exp = math.ceil(math.log10(p))
    if 10.0 ** exp <= p:            # p sits exactly on a power of ten
        exp += 1
    return f'p < 10$^{{{exp}}}$'


def fig1b_null_distribution():
    """Plot 1M null distribution histogram with observed distance marked."""
    print("Generating Fig 1B: 1M null distribution...")

    with open(BASE / 'analysis/permutation_1M/results_1M.json') as f:
        results = json.load(f)
    null_dist = np.load(BASE / 'analysis/permutation_1M/null_distribution_1M.npy')

    obs_dist = results['observed_procrustes_distance']
    obs_null = results['obs_null_ratio']
    p_str = _p_bound(results['p_value'])

    fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.75))

    # Histogram of null
    ax.hist(null_dist, bins=120, color=C_LIGHTGRAY, edgecolor='white',
            linewidth=0.3, density=True, alpha=0.9, zorder=2)

    # Observed line
    ax.axvline(obs_dist, color=C_BLUE, linewidth=1.5, linestyle='-', zorder=5)

    # Observed label inside plot, near top of line, nudged right (matches Fig 2's pattern).
    ax.text(obs_dist + 1.0, 0.92, f'Observed\n({obs_dist:.2f})',
            transform=ax.get_xaxis_transform(), fontsize=FONT_SIZE_ANNOT,
            color=C_BLUE, ha='left', va='top')

    # Stats text — truncate to 3dp to match manuscript (0.52296 → 0.522).
    # Anchored lower-left (matches fig1c's iter-3 layout) so the upper area is
    # free for the Observed label. Histogram bars sit at high x (x_frac > 0.6)
    # so this corner is clear.
    obs_null_3dp = int(obs_null * 1000) / 1000
    ax.text(0.06, 0.03,
            f'obs/null = {obs_null_3dp:.3f}\n{p_str}\nn = 1,000,000 perms',
            transform=ax.transAxes, fontsize=FONT_SIZE_ANNOT,
            ha='left', va='bottom', color=C_DARKGRAY,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor=C_LIGHTGRAY, alpha=0.9))

    ax.set_xlabel('Procrustes distance', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('Density', fontsize=FONT_SIZE_LABEL)
    clean_spine(ax)

    save_figure(fig, PANELS / 'fig1b_null_1M')
    return fig


# ===================================================================
# Fig 1C: Lineage-stratified null distribution
# ===================================================================
def fig1c_lineage_stratified():
    """Plot lineage-stratified permutation null vs global null."""
    print("Generating Fig 1C: Lineage-stratified null...")

    with open(BASE / 'output/validation/lineage_stratified/lineage_stratified_results.json') as f:
        data = json.load(f)

    # Load global null for overlay
    null_global = np.load(BASE / 'output/phase2/scaled_35types/null_distribution_35.npy')

    obs = data['observed_procrustes_distance']
    strat = data['stratified_null']

    # Load real stratified null distribution array
    strat_null = np.load(BASE / 'output/validation/lineage_stratified/null_distribution_stratified.npy')

    fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.75))

    # Global null
    ax.hist(null_global, bins=80, color=C_LIGHTGRAY, edgecolor='white',
            linewidth=0.3, density=True, alpha=0.7, label='Global null', zorder=2)

    # Stratified null (real data)
    ax.hist(strat_null, bins=80, color=C_ORANGE, edgecolor='white',
            linewidth=0.3, density=True, alpha=0.5, label='Within-lineage null', zorder=3)

    # Observed line
    ax.axvline(obs, color=C_BLUE, linewidth=1.5, linestyle='-', zorder=5)

    # Observed label inside plot, near top of line, nudged right (matches Fig 2's pattern).
    ax.text(obs + 1.0, 0.92, f'Observed\n({obs:.2f})',
            transform=ax.get_xaxis_transform(), fontsize=FONT_SIZE_ANNOT,
            color=C_BLUE, ha='left', va='top')

    # Truncate global obs/null to 3dp (matching fig1b: 0.52296→0.522, not rounded 0.523)
    _global_obs_null = int(data['global_null']['obs_null_ratio'] * 1000) / 1000
    # Stats inset placed in the lower-left empty zone — left of the orange
    # within-lineage histogram (which starts at x_frac ≈ 0.31) and below the
    # 'upper left' legend (which ends at y_frac ≈ 0.85). The right tail of
    # the gray Global histogram now occupies the upper-right of the panel,
    # so the previous (0.62, 0.97) center-top position clipped its left edge.
    ax.text(0.06, 0.78,
            f'Within-lineage:\nobs/null = {strat["obs_null_ratio"]:.3f}\n'
            f'{_p_bound(strat["p_value"])}\n\n'
            f'Global:\nobs/null = {_global_obs_null:.3f}',
            transform=ax.transAxes, fontsize=FONT_SIZE_ANNOT,
            ha='left', va='top', color=C_DARKGRAY, zorder=10,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor=C_LIGHTGRAY, alpha=1.0))

    ax.set_xlabel('Procrustes distance', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('Density', fontsize=FONT_SIZE_LABEL)
    # Legend outside plot above panel — upper-right placement covered the gray
    # histogram peak; outside-above keeps the data area uncluttered. ncol=2
    # so both entries sit on one row. Reorder so the legend reads
    # left-to-right matching the histogram x-ordering (orange Within-lineage
    # at x≈90 left of gray Global at x≈115).
    handles, labels = ax.get_legend_handles_labels()
    label_to_handle = dict(zip(labels, handles))
    desired_order = ['Within-lineage null', 'Global null']
    ax.legend([label_to_handle[l] for l in desired_order if l in label_to_handle],
              [l for l in desired_order if l in label_to_handle],
              fontsize=FONT_SIZE_LEGEND, frameon=False,
              loc='lower center', bbox_to_anchor=(0.5, 1.02), ncol=2)
    clean_spine(ax)

    save_figure(fig, PANELS / 'fig1c_lineage_stratified')
    return fig


# ===================================================================
# Fig 2A: Rigidity ranking bar chart — THE ICONIC FIGURE
# ===================================================================
def fig2a_rigidity_ranking():
    """Horizontal bar chart of 35 cell types by Procrustes residual."""
    print("Generating Fig 2A: Rigidity ranking (ICONIC)...")

    with open(BASE / 'output/cellcount_confound/cellcount_confound_results.json') as f:
        data = json.load(f)

    # Extract per-type data, sorted by residual (ascending = most rigid at bottom)
    types = data['per_type_data']
    # Already sorted descending by residual (flexible at top, rigid at bottom)
    # Reverse so most rigid (smallest residual) at bottom → reading top to bottom
    # gives flexible to rigid

    df = pd.DataFrame(types)
    df = df.sort_values('residual_magnitude', ascending=True)  # rigid at top

    # Widened canvas (was COL1) absorbs the whitespace that appeared between
    # bar area and lineage legend after Fig 3 went to the 2-panel layout.
    fig, ax = plt.subplots(figsize=(COL2 * 0.7, COL2 * 0.85))

    y_pos = np.arange(len(df))
    colors = [lineage_color(ct) for ct in df['cell_type']]

    bars = ax.barh(y_pos, df['residual_magnitude'], height=0.75,
                   color=colors, edgecolor='white', linewidth=0.3, zorder=3)

    # Long-form labels via centralized LABEL_EXPAND (figure_style.py) — kept
    # in sync with Fig 1E and Fig 4A.
    ax.set_yticks(y_pos)
    ax.set_yticklabels(
        [LABEL_EXPAND.get(short_name(ct), short_name(ct)) for ct in df['cell_type']],
        fontsize=6.5)
    ax.set_xlabel('Procrustes residual magnitude', fontsize=FONT_SIZE_LABEL)
    ax.invert_yaxis()  # most rigid at top

    # Add rank numbers on right side
    for i, (_, row) in enumerate(df.iterrows()):
        rank = len(df) - i  # rank 1 = most flexible, rank 35 = most rigid
        # Actually: rigid = smallest residual. Let me re-derive.
        # Sorted ascending: idx 0 = smallest residual = most rigid = rank 35
        rank = 35 - i
        ax.text(row['residual_magnitude'] + 0.15, i, str(rank),
                fontsize=5.5, va='center', ha='left', color=C_DARKGRAY)

    # Standalone lineage legend — place in upper-right empty space
    # (rigid types have short bars, leaving room at top-right)
    add_lineage_legend(ax, loc='upper right', ncol=1, title='Lineage')

    ax.set_xlim(0, max(df['residual_magnitude']) * 1.15)
    clean_spine(ax)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0)

    save_figure(fig, PANELS / 'fig2a_rigidity_ranking')
    return fig


# ===================================================================
# Fig 2B: CellMarker enrichment
# ===================================================================
def fig2b_cellmarker():
    """Bar chart showing CellMarker identity gene enrichment."""
    print("Generating Fig 2B: CellMarker enrichment...")

    with open(BASE / 'output/validation/cellmarker_35type_rerun/cellmarker_35type_results.json') as f:
        data = json.load(f)

    fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.6))

    # Two bars: primary enrichment and expression-matched
    categories = ['Primary\nenrichment', 'Expression-\nmatched']
    values = [data['global_enrichment']['enrichment'],
              data['expression_matched']['enrichment']]
    pvals = [data['global_enrichment']['p_value'],
             data['expression_matched']['p_value']]
    colors_bars = [C_BLUE, C_TEAL]

    # Approximate 95% CI: enrichment ≈ obs/exp, SE ≈ enrichment/sqrt(obs)
    obs_counts = [data['global_enrichment']['observed'],
                  data['expression_matched']['observed']]
    ci_half = [1.96 * v / np.sqrt(max(o, 1)) for v, o in zip(values, obs_counts)]

    bars = ax.bar(categories, values, width=0.55, color=colors_bars,
                  edgecolor='white', linewidth=0.5, zorder=3,
                  yerr=ci_half, capsize=3, error_kw=dict(lw=0.8, color=C_DARKGRAY))

    # Add exact p-value annotations with LaTeX superscripts
    p_labels = [r'$p = 2.1 \times 10^{-13}$', r'$p = 1.2 \times 10^{-12}$']
    for bar, p_label, ci in zip(bars, p_labels, ci_half):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + ci + 0.15,
                p_label, ha='center', va='bottom',
                fontsize=FONT_SIZE_ANNOT, color=C_DARKGRAY)

    # Reference line at enrichment = 1
    ax.axhline(1, color=C_DARKGRAY, linewidth=0.5, linestyle='--', zorder=1)
    ax.text(1.5, 1.15, 'No enrichment', fontsize=5.5, color=C_DARKGRAY,
            ha='right', va='bottom', alpha=0.7)

    ax.set_ylabel('Fold enrichment', fontsize=FONT_SIZE_LABEL)
    ax.set_ylim(0, max(values) * 1.3)
    clean_spine(ax)

    # Additional info — gap between bars, above "No enrichment" line
    ax.text(0.50, 0.28,
            f'{data["per_type_pass"]} types pass\n'
            f'{data["global_enrichment"]["observed"]} genes observed\n'
            f'{data["global_enrichment"]["expected"]:.1f} expected',
            transform=ax.transAxes, fontsize=FONT_SIZE_ANNOT,
            ha='center', va='bottom', color=C_DARKGRAY)

    save_figure(fig, PANELS / 'fig2b_cellmarker')
    return fig


# ===================================================================
# Fig 2C: Cell count confound null
# ===================================================================
def fig2c_cellcount():
    """Scatter plot of cell count vs rigidity showing no confound."""
    print("Generating Fig 2C: Cell count confound...")

    with open(BASE / 'output/cellcount_confound/cellcount_confound_results.json') as f:
        data = json.load(f)

    df = pd.DataFrame(data['per_type_data'])

    fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.75))

    colors = [lineage_color(ct) for ct in df['cell_type']]
    ax.scatter(df['min_cells'], df['residual_magnitude'],
               c=colors, s=25, edgecolors='white', linewidths=0.3,
               zorder=3, alpha=0.85)

    ax.set_xlabel('Min cell count (human, mouse)', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('Procrustes residual', fontsize=FONT_SIZE_LABEL)

    # Stats inset created BEFORE the extreme-label block so adjust_text can
    # repel labels from it. Lower-LEFT placement: only three unlabeled
    # points there (adventitial, granulocyte, pancreatic ductal); the
    # extreme labels Epithelial (UL), Stromal/HPC (UR), and the rigid
    # cluster CD8+T/Endothelial/Hepatocyte (LR) are all clear of LL.
    inset_artist = ax.text(0.03, 0.03,
            f'Spearman $\\rho$ = {data["spearman_rho"]:.3f}\n'
            f'         p = {data["spearman_p"]:.3f}\n'
            f' Partial $\\rho$ = {data["partial_rho"]:.3f}\n'
            f'         p = {data["partial_p"]:.3f}',
            transform=ax.transAxes, fontsize=FONT_SIZE_ANNOT,
            ha='left', va='bottom', color=C_DARKGRAY,
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                      edgecolor=C_LIGHTGRAY, alpha=0.9))

    clean_spine(ax)

    # Manuscript-canonical extreme labels for identifiability (top-3 most
    # rigid + bottom-3 most flexible per Fig 3 legend). Expanded forms match
    # Fig 1E / Fig 2A. adjustText handles non-overlap among labels AND
    # repels them from the stats inset (passed via `objects=`).
    from adjustText import adjust_text
    extreme_labels = [
        ('CD8-positive, alpha-beta T cell', 'CD8+ T cell'),       # rank 35
        ('endothelial cell',                'Endothelial cell'),  # rank 33
        ('hepatocyte',                      'Hepatocyte'),        # rank 32
        ('hematopoietic precursor cell',    'Hematopoietic precursor cell'),  # rank 3
        ('epithelial cell',                 'Epithelial cell'),   # rank 2
        ('stromal cell',                    'Stromal cell'),      # rank 1
    ]
    df_idx = df.set_index('cell_type')
    label_artists = []
    for ct, label in extreme_labels:
        if ct in df_idx.index:
            x = df_idx.loc[ct, 'min_cells']
            y = df_idx.loc[ct, 'residual_magnitude']
            label_artists.append(
                ax.text(x, y, label, fontsize=5.5, color=C_DARKGRAY,
                        zorder=5))
    adjust_text(label_artists, ax=ax, objects=[inset_artist],
                arrowprops=dict(arrowstyle='-', color=C_GRAY, lw=0.4),
                expand=(1.2, 1.5))
    # Lineage-color cross-reference text dropped — redundant with the Fig 3
    # caption (which describes the lineage palette for Panel A and notes
    # Panel B uses the same scheme). Removing it freed UR for the Stromal
    # label without crowding.
    save_figure(fig, PANELS / 'fig2c_cellcount')
    return fig


# ===================================================================
# Fig 4D: Summary obs/null across 4 datasets
# ===================================================================
def fig4d_replication_summary():
    """Grouped bar plot showing obs/null ratios across all replication
    attempts: primary + 4 successes + 2 failures (Andrews et al. and
    MCA × HCA). Successes use filled bars; failures use hatched bars.
    """
    print("Generating Fig 4D: Replication summary (all 7 datasets)...")

    def _fmt_p(p):
        """Bound form at a permutation floor, value form otherwise.

        Delegates the bound to the module-level _p_bound so this panel draws
        the same notation as Fig 1 and Fig 2 rather than the '1e-6' dialect it
        used before. The two hand-written branches it replaces chose between a
        fixed 1e-6 and a fixed 1e-4; _p_bound derives whichever power of ten
        the value actually sits below, which is the same answer for every p
        this panel receives and stays correct if one moves.
        """
        if p <= 1e-4:
            return _p_bound(p)
        return f'p = {p:.2g}'

    # Load all replication results from source JSONs (no hardcoding for
    # the primary set; failures pulled from their own JSONs).
    #
    # Each path is spelled ONCE. It is both the file this panel reads and the key
    # it joins the matched-n baselines on below, so a second spelling could drift
    # from the first and send a bar's baseline to the wrong arm without failing.
    P_PRIMARY = 'analysis/permutation_1M/results_1M.json'
    P_SUN2023 = 'output/validation/sun2023_replication_expanded/sun2023_expanded.json'
    P_PANSCI = 'output/validation/pansci_replication/pansci_replication.json'
    P_CELLHINT = 'output/validation/cellhint_replication/cellhint_replication.json'
    P_PANCENSUS = 'analysis/census_replication/replication_results.json'
    P_ANDREWS = 'output/validation/andrews_replication/andrews_replication_results.json'
    P_MCA_HCA = 'output/validation/t1a_replication/t1a_results.json'

    def _load(rel):
        with open(BASE / rel) as f:
            return json.load(f)

    primary = _load(P_PRIMARY)
    sun2023 = _load(P_SUN2023)
    pansci = _load(P_PANSCI)
    cellhint = _load(P_CELLHINT)
    pancensus = _load(P_PANCENSUS)
    andrews = _load(P_ANDREWS)
    mca_hca = _load(P_MCA_HCA)

    # Extract obs/null, n, and p from each source
    # (success, dataset_label, source path, obs_null, n, p)
    sources = [
        (True,  'Primary\n(TS×TMS)', P_PRIMARY,
         primary['obs_null_ratio'],
         primary['n_cell_types'], primary['p_value']),
        (True,  'Sun2023\n(10x v3)', P_SUN2023,
         sun2023['procrustes']['obs_null_ratio'],
         sun2023['procrustes']['n_types'], sun2023['procrustes']['p_value']),
        (True,  'PanSci\n(EasySci)', P_PANSCI,
         pansci['procrustes']['obs_null_ratio'],
         pansci['procrustes']['n_types'], pansci['procrustes']['p_value']),
        (True,  'CellHint\n(Human)', P_CELLHINT,
         cellhint['procrustes']['obs_null_ratio'],
         cellhint['procrustes']['n_types'], cellhint['procrustes']['p_value']),
        (True,  'pan-Census\n(pooled)', P_PANCENSUS,
         pancensus['permutation_test']['obs_null_ratio'],
         pancensus['n_cell_types'], pancensus['permutation_test']['p_value']),
        (False, 'Andrews\n(liver)', P_ANDREWS,
         andrews['obs_null_ratio'],
         andrews['n_types'], andrews['p_value']),
        (False, 'MCA × HCA\n(microwell)', P_MCA_HCA,
         mca_hca['t1a_procrustes']['obs_null_ratio'],
         mca_hca['t1a_procrustes']['n_types'], mca_hca['t1a_procrustes']['p_value']),
    ]

    success_flags = [s[0] for s in sources]
    datasets = [s[1] for s in sources]
    source_paths = [s[2] for s in sources]
    obs_null = [s[3] for s in sources]
    n_types = [s[4] for s in sources]
    p_values = [s[5] for s in sources]

    # ---- matched-n baselines, joined on the source path --------------------
    #
    # Each bar is read against the primary restricted to that bar's own matched
    # types, because the null median moves with type count and a replication's
    # obs/null is not comparable to the primary's 0.522 at 35. The baselines come
    # from analysis/ranking_replication/block2_matched_n.py.
    #
    # The join is on the source JSON path, and it is the path because nothing else
    # works. That producer's labels and this panel's differ three separate ways --
    # separator (space vs newline), multiplication sign (x vs U+00D7) and, for
    # MCAxHCA, the text itself -- so a label join matches nothing at all. The two
    # lists happen to be in the same order, so a positional join would be correct
    # today with nothing asserting it stays correct; the failure it would produce
    # after a reorder is six bars silently carrying each other's baselines. The
    # path is the one key both sides already hold identically, and every lookup
    # below must resolve to exactly one entry or the panel refuses to draw.
    with open(BASE / 'analysis/ranking_replication/block2_matched_n_results.json') as f:
        matched = json.load(f)

    by_path = {}
    for bar in matched['bars']:
        by_path.setdefault(bar['path'], []).append(bar)

    baselines = []
    for label, path in zip(datasets, source_paths):
        hits = by_path.get(path, [])
        if len(hits) != 1:
            raise ValueError(
                f"Fig 3 baseline join failed for {label.replace(chr(10), ' ')!r}: "
                f"{path} resolved to {len(hits)} entries in "
                f"block2_matched_n_results.json, expected exactly 1. "
                f"Paths that artifact carries: {sorted(by_path)}")
        baselines.append(hits[0])

    # The seventh arm is the primary, and it has no baseline to draw: its matched
    # type set is the whole primary set, so block2 hands it back its own value.
    # Derived from the artifact's own faithfulness gate rather than by matching a
    # label, so it stays right if the panel's bar order or wording changes.
    full_n = matched['gate']['n']
    draw_baseline = [b['n_inter'] != full_n for b in baselines]
    # Color palette: successes get distinct colors; failures rendered gray.
    success_colors = [C_BLUE, C_ORANGE, C_PURPLE, C_TEAL, C_GRAY]
    success_idx = 0
    bar_colors = []
    for ok in success_flags:
        if ok:
            bar_colors.append(success_colors[success_idx])
            success_idx += 1
        else:
            bar_colors.append('#D8D8D8')  # light gray for failures

    # Full-width source panel (matches the full-row slot in the new 3-row
    # Fig 3 composite layout — X-fig3-restructure). Bar labels and per-bar
    # n=/p= annotations render at native size without scaling.
    fig, ax = plt.subplots(figsize=(COL2, COL1 * 0.75))

    x = np.arange(len(datasets))
    bars = []
    for i, (xi, ratio, color, ok) in enumerate(zip(x, obs_null, bar_colors, success_flags)):
        if ok:
            bar = ax.bar(xi, ratio, width=0.6, color=color,
                         edgecolor='white', linewidth=0.5, zorder=3)
        else:
            bar = ax.bar(xi, ratio, width=0.6, color=color,
                         edgecolor=C_DARKGRAY, linewidth=0.7, zorder=3,
                         hatch='///')
        bars.append(bar[0])

    # Matched-n baseline on each replication bar: the primary's own obs/null when
    # it is restricted to that bar's matched types, drawn as a segment across the
    # bar rather than as a separate bar, so it reads as a level this bar is above
    # and not as a seventh measurement. Same idiom fig3b_pre_post uses for its
    # label-shuffle nulls, and for the same reason -- a reader who cannot see the
    # level a bar is meant to clear reads the bar's height as the result.
    #
    # Six segments across seven bars. The primary is its own baseline, and drawing
    # a segment 0.0001 from its own bar top would assert a comparison that is not
    # one; that 0.0001 is the 1,000,000-permutation bar against this baseline's
    # 10,000-permutation null median, not a gap. Fig 3's caption names it.
    for bar, base, draw in zip(bars, baselines, draw_baseline):
        if not draw:
            continue
        level = base['matched_n_primary']['obs_null']
        half = bar.get_width() * 0.46
        centre = bar.get_x() + bar.get_width() / 2
        ax.plot([centre - half, centre + half], [level] * 2,
                color=C_DARKGRAY, lw=1.1, solid_capstyle='butt', zorder=6,
                label='Matched-n primary baseline' if bar is bars[1] else None)

    # Add n and p labels above each bar
    for bar, n, p in zip(bars, n_types, p_values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                f'n={n}\n{_fmt_p(p)}', ha='center', va='bottom',
                fontsize=5.5, color=C_DARKGRAY)

    ax.set_xticks(x)
    ax.set_xticklabels(datasets, fontsize=6)
    for lbl, c, ok in zip(ax.get_xticklabels(), bar_colors, success_flags):
        if ok:
            lbl.set_color(c)
    ax.set_ylabel('Obs/null ratio', fontsize=FONT_SIZE_LABEL)
    ax.set_ylim(0, 1.15)

    # Reference lines: chance at 1.0; failure-cutoff annotation (≈0.82).
    # Annotation pinned to the far-left of the axes to avoid colliding with
    # the n=17, p=0.54 annotation above the MCA × HCA bar at the right edge.
    ax.axhline(1.0, color=C_LIGHTGRAY, linewidth=0.7, linestyle='--', zorder=1)
    ax.text(-0.5, 1.005, 'chance (obs/null = 1)',
            fontsize=5.0, color=C_GRAY, ha='left', va='bottom')

    # Upper left, but well below the chance line: that line's annotation is drawn
    # at data y = 1.005 starting left of the first bar, and an upper-left legend
    # lands on top of it. The four left-hand bars and their n=/p= blocks top out
    # near 0.62 and the tall bars are all at the right, so the band from about 0.7
    # to 0.95 above the left half is the one clear region on this panel.
    ax.legend(loc='upper left', bbox_to_anchor=(0.0, 0.82), fontsize=5.2,
              frameon=False, handlelength=1.4, borderpad=0.0, handletextpad=0.5)

    clean_spine(ax)
    save_figure(fig, PANELS / 'fig4d_replication_summary')
    return fig


# ===================================================================
# Fig 4E: Human control comparison
# ===================================================================
def fig4e_human_control():
    """Bar chart comparing cross-species vs within-species obs/null ratios."""
    print("Generating Fig 4E: Human control comparison...")

    # 6-type comparison from sun2023 replication
    with open(BASE / 'output/validation/sun2023_replication/sun2023_replication.json') as f:
        sun = json.load(f)

    # Within-species control from negctrl_v2
    with open(BASE / 'output/phase2/negative_control_v2/negctrl_v2_results.json') as f:
        negctrl = json.load(f)

    # Use primary analysis p-value (0.0035) for cross-species, not sun2023
    # comparison_b re-analysis (0.0023) — see v32.3 fix
    with open(BASE / 'output/phase2/procrustes_results.json') as f:
        primary = json.load(f)
    cross_obs_null = sun['comparison_b']['obs_null_ratio']
    cross_p = primary['permutation_test']['p_value']
    within_obs_null = negctrl['procrustes']['distance'] / negctrl['permutation_test']['null_distribution_summary']['median']
    within_p = negctrl['permutation_test']['p_value']
    fold_diff = within_obs_null / cross_obs_null

    fig, ax = plt.subplots(figsize=(COL1 * 0.85, COL1 * 0.7))

    # Individual cell type residuals as dots
    cross_types = sun['comparison_b']['cell_types']
    cross_resids = [sun['comparison_b']['per_type_residuals'][ct]['magnitude']
                    for ct in cross_types]
    within_resids = [negctrl['per_type_residuals'][ct]['magnitude']
                     for ct in cross_types if ct in negctrl.get('per_type_residuals', {})]

    # If per-type residuals not available for negctrl, show aggregate comparison
    if not within_resids:
        # Fallback: show the two aggregate values as larger dots with CI
        categories = ['Cross-species\n(H×M)', 'Within-species\n(H×H)']
        x_pos = [0, 1]
        ax.scatter(x_pos, [cross_obs_null, within_obs_null],
                   c=[C_BLUE, C_ORANGE], s=80, edgecolors='white',
                   linewidths=0.8, zorder=5)
        # Cross-species label: right of dot to clear y-axis tick collision
        ax.text(x_pos[0] + 0.20, cross_obs_null + 0.01, format_p(cross_p),
                ha='left', va='center', fontsize=FONT_SIZE_ANNOT, color=C_DARKGRAY)
        # Within-species label: above dot
        ax.text(x_pos[1], within_obs_null + 0.02, format_p(within_p),
                ha='center', va='bottom', fontsize=FONT_SIZE_ANNOT, color=C_DARKGRAY)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(categories, fontsize=FONT_SIZE_TICK)
    else:
        # Show individual cell type dots
        jitter = 0.08
        for i, (ct, r) in enumerate(zip(cross_types, cross_resids)):
            ax.scatter(0 + np.random.uniform(-jitter, jitter), r,
                       c=C_BLUE, s=30, edgecolors='white', linewidths=0.3, zorder=3)
        for i, r in enumerate(within_resids):
            ax.scatter(1 + np.random.uniform(-jitter, jitter), r,
                       c=C_ORANGE, s=30, edgecolors='white', linewidths=0.3, zorder=3)
        # Median lines
        ax.plot([-0.15, 0.15], [np.median(cross_resids)] * 2, color=C_BLUE,
                linewidth=1.5, zorder=4)
        ax.plot([0.85, 1.15], [np.median(within_resids)] * 2, color=C_ORANGE,
                linewidth=1.5, zorder=4)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['Cross-species\n(H×M)', 'Within-species\n(H×H)'],
                           fontsize=FONT_SIZE_TICK)
        ax.set_ylabel('Per-type residual magnitude', fontsize=FONT_SIZE_LABEL)

    # Summary stats — truncate fold to 2dp to match manuscript (1.9164 → 1.91)
    fold_diff_2dp = int(fold_diff * 100) / 100
    ax.text(0.95, 0.70,
            f'Obs/null: {cross_obs_null:.3f} vs {within_obs_null:.3f}\n'
            f'{fold_diff_2dp:.2f}\u00d7 difference\n'
            f'{len(cross_types)} cell types',
            transform=ax.transAxes, fontsize=FONT_SIZE_ANNOT,
            ha='right', va='top', color=C_DARKGRAY,
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                      edgecolor=C_LIGHTGRAY, alpha=0.9))

    clean_spine(ax)
    save_figure(fig, PANELS / 'fig4e_human_control')
    return fig


# ===================================================================
# Fig 5A: Three-species Procrustes coherence
# ===================================================================
def fig5a_macaque_primary():
    """Null distribution for three-species (macaque) primary analysis."""
    print("Generating Fig 5A: Macaque primary...")

    with open(BASE / 'output/macaque_pipeline/primary_procrustes_results.json') as f:
        data = json.load(f)

    obs = data['procrustes']['distance']
    perm = data['permutation_test']
    null_median = perm['null_median']

    # Summary visualization — no raw null array exists; show what we have
    fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.55))

    ax.plot(null_median, 1, 's', color=C_GRAY, markersize=8, zorder=4,
            markeredgecolor='white', markeredgewidth=0.5, label='Null median')
    ax.plot(obs, 0, 'o', color=C_PURPLE, markersize=9, zorder=5,
            markeredgecolor='white', markeredgewidth=0.5, label='Observed')

    # Connecting line to show distance
    ax.plot([obs, null_median], [0, 1], color=C_LIGHTGRAY, linewidth=0.8,
            linestyle=':', zorder=1)

    ax.set_yticks([0, 1])
    ax.set_yticklabels(['Observed', 'Null median\n(10K perms)'], fontsize=FONT_SIZE_TICK)
    ax.set_xlabel('Procrustes distance', fontsize=FONT_SIZE_LABEL)
    ax.set_ylim(-0.5, 1.8)

    ax.text(0.55, 0.95,
            f'3 species, {data["n_types"]} types\n'
            f'obs/null = {perm["obs_null_ratio"]:.3f}\n'
            f'{format_p(perm["p_value"])}\n'
            f'10K perms',
            transform=ax.transAxes, fontsize=FONT_SIZE_ANNOT,
            ha='left', va='top', color=C_DARKGRAY,
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                      edgecolor=C_LIGHTGRAY, alpha=0.9))

    # No separate legend — stats box is sufficient
    clean_spine(ax)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0)

    save_figure(fig, PANELS / 'fig5a_macaque_primary')
    return fig


# ===================================================================
# Fig 5B: Sensitivity (13 RIRA-only types)
# ===================================================================
def fig5b_macaque_sensitivity():
    """Null distribution for RIRA-only sensitivity analysis."""
    print("Generating Fig 5B: Macaque sensitivity...")

    with open(BASE / 'output/macaque_pipeline/sensitivity_procrustes_results.json') as f:
        data = json.load(f)

    obs = data['procrustes']['distance']
    perm = data['permutation_test']
    # Compute null median from obs/null ratio (raw array not saved)
    null_median = obs / perm['obs_null_ratio']

    # Summary visualization — no raw null array exists; show what we have
    fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.55))

    ax.plot(null_median, 1, 's', color=C_GRAY, markersize=8, zorder=4,
            markeredgecolor='white', markeredgewidth=0.5, label='Null median')
    ax.plot(obs, 0, 'o', color=C_TEAL, markersize=9, zorder=5,
            markeredgecolor='white', markeredgewidth=0.5, label='Observed')

    ax.plot([obs, null_median], [0, 1], color=C_LIGHTGRAY, linewidth=0.8,
            linestyle=':', zorder=1)

    ax.set_yticks([0, 1])
    ax.set_yticklabels(['Observed', 'Null median\n(10K perms)'], fontsize=FONT_SIZE_TICK)
    ax.set_xlabel('Procrustes distance', fontsize=FONT_SIZE_LABEL)
    ax.set_ylim(-0.5, 1.8)

    ax.text(0.97, 0.95,
            f'{len(data["cell_types"])} RIRA-only types\n'
            f'obs/null = {perm["obs_null_ratio"]:.3f}\n'
            f'{format_p(perm["p_value"])}',
            transform=ax.transAxes, fontsize=FONT_SIZE_ANNOT,
            ha='right', va='top', color=C_DARKGRAY,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor=C_LIGHTGRAY, alpha=0.9))

    clean_spine(ax)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0)

    save_figure(fig, PANELS / 'fig5b_macaque_sensitivity')
    return fig


# ===================================================================
# Fig 5C: Hepatocyte rigidity rank reversal across analyses
# ===================================================================
def fig5c_hepatocyte_rigidity():
    """Slope chart showing hepatocyte rank reversal: rigid in two-species,
    most flexible in three-species analysis.

    Biology: Hepatocyte ranks 32/35 (rigid) in the primary human-mouse
    analysis but 20/20 (most flexible) in the three-species macaque
    analysis. This reversal demonstrates atlas-conditionality of per-type
    rigidity rankings — a key limitation. Macaque hepatocyte pairing is
    LOW_CONFIDENCE (Qu et al., 2 donors).

    Math: Ranks converted to rigidity percentiles (100%=most rigid,
    0%=most flexible) for cross-analysis comparison on a common scale.
    All 20 matched types shown as context slopes; hepatocyte highlighted.
    """
    print("Generating Fig 5C: Hepatocyte rigidity rank reversal...")

    mac_df = pd.read_csv(BASE / 'output/macaque_pipeline/rigidity_ranking_comparison.csv')

    primary_total = 35
    mac_total = len(mac_df)  # 20

    # Convert ranks to rigidity percentiles (100% = most rigid, 0% = most flexible)
    # hm_rank: 1=most flexible, 35=most rigid → pct = hm_rank / primary_total * 100
    # mac_rank: 1=most rigid (smallest residual), N=most flexible (largest residual)
    #   → rigid equivalent = mac_total + 1 - mac_rank
    #   → pct = rigid_equiv / mac_total * 100
    mac_df['primary_pct'] = mac_df['hm_rank'] / primary_total * 100
    mac_df['mac_pct'] = (mac_total + 1 - mac_df['mac_rank']) / mac_total * 100

    fig, ax = plt.subplots(figsize=(COL1 * 0.85, COL1 * 0.7))

    x = [0, 1]

    # Background: all other matched types as faint context slopes
    for _, row in mac_df.iterrows():
        if row['cell_type'] == 'hepatocyte':
            continue
        lc = row['LOW_CONFIDENCE']
        ax.plot(x, [row['primary_pct'], row['mac_pct']],
                color='#E0E0E0' if lc else C_LIGHTGRAY,
                linewidth=0.5, linestyle=':' if lc else '-',
                alpha=0.35, zorder=1)

    # Hepatocyte: highlighted slope showing the reversal
    hep = mac_df[mac_df['cell_type'] == 'hepatocyte'].iloc[0]
    y_left = hep['primary_pct']   # ~91.4 (rigid)
    y_right = hep['mac_pct']      # ~5.0  (flexible)

    ax.plot(x, [y_left, y_right], color=C_ORANGE, linewidth=2.0, zorder=4)
    ax.scatter(x, [y_left, y_right], color=C_ORANGE, s=45, zorder=5,
               edgecolor='white', linewidth=0.8)

    # Combined label to the RIGHT of the left data point — single two-line annotation
    ax.text(0.06, y_left, f'Hepatocyte\nRank {int(hep["hm_rank"])}/{primary_total} (rigid)',
            fontsize=FONT_SIZE_ANNOT, color=C_BLACK,
            ha='left', va='center')

    # Right data point label
    ax.text(1.06, y_right, f'Rank {int(hep["mac_rank"])}/{mac_total}\n(most flexible)',
            fontsize=FONT_SIZE_ANNOT, color=C_DARKGRAY,
            ha='left', va='center')

    # LOW CONFIDENCE badge — centered between data points at y=55
    ax.text(0.5, 55, 'LOW CONFIDENCE',
            fontsize=8, color=C_ORANGE, ha='center', va='center',
            fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                      edgecolor=C_ORANGE, alpha=0.9, linewidth=0.5))

    # Median reference line
    ax.axhline(50, color=C_LIGHTGRAY, linewidth=0.5, linestyle='--', zorder=0)
    ax.text(1.02, 51, '50th', fontsize=5, color=C_GRAY, ha='left', va='bottom')

    # Axis formatting
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylim(-5, 115)
    ax.set_xticks(x)
    ax.set_xticklabels([f'Primary\n(H\u00d7M, n={primary_total})',
                         f'Three-species\n(n={mac_total})'],
                        fontsize=FONT_SIZE_TICK)
    ax.set_ylabel('Divergence percentile', fontsize=FONT_SIZE_LABEL)
    ax.set_title('Rank reversal across analyses', fontsize=FONT_SIZE_LABEL, pad=8)

    clean_spine(ax)
    save_figure(fig, PANELS / 'fig5c_hepatocyte_rigidity')
    return fig


# ===================================================================
# Fig S1A: Independent-PCA null distribution
# ===================================================================
def figs1a_indep_pca_null():
    """Plot independent-PCA null distribution (1M permutations)."""
    print("Generating Fig S1A: Independent-PCA null...")

    with open(BASE / 'analysis/independent_pca_sensitivity/independent_pca_results.json') as f:
        data = json.load(f)
    null_dist = np.load(BASE / 'analysis/independent_pca_sensitivity/null_distribution_independent_pca.npy')

    obs = data['procrustes']['distance']
    perm = data['permutation_test']

    fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.75))

    ax.hist(null_dist, bins=120, color=C_LIGHTGRAY, edgecolor='white',
            linewidth=0.3, density=True, alpha=0.9, zorder=2)
    ax.axvline(obs, color=C_ORANGE, linewidth=1.5, zorder=5)

    # Annotation — place in empty space to right of observed line
    ymax = ax.get_ylim()[1]
    ax.annotate(f'Observed ({obs:.1f})',
                xy=(obs, ymax * 0.85), xytext=(obs + 14, ymax * 0.55),
                fontsize=FONT_SIZE_ANNOT, color=C_ORANGE,
                arrowprops=dict(arrowstyle='->', color=C_ORANGE, lw=0.8),
                ha='center', va='top')

    # Place stat box in the gap between observed line and null distribution
    ax.text(0.42, 0.95,
            f'Independent PCA\nobs/null = {perm["obs_null_ratio"]:.3f}\n'
            r'$p < 10^{-6}$' + '\nn = 1,000,000 iterations',
            transform=ax.transAxes, fontsize=FONT_SIZE_ANNOT,
            ha='center', va='top', color=C_DARKGRAY,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor=C_LIGHTGRAY, alpha=0.9))

    ax.set_xlabel('Procrustes distance', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('Density', fontsize=FONT_SIZE_LABEL)
    clean_spine(ax)

    save_figure(fig, PANELS / 'figs1a_indep_pca_null')
    return fig


# ===================================================================
# Fig S1B: Joint vs independent PCA ranking scatter
# ===================================================================
def figs1b_ranking_scatter():
    """Scatter plot of joint-PCA vs independent-PCA rigidity rankings."""
    print("Generating Fig S1B: Joint vs independent PCA ranking scatter...")

    with open(BASE / 'analysis/independent_pca_sensitivity/independent_pca_results.json') as f:
        data = json.load(f)

    comp = data['comparison_to_joint_pca']
    residuals = data['residuals']

    # Build rank pairs
    cell_types = list(residuals.keys())
    joint_ranks = [residuals[ct]['joint_rank'] for ct in cell_types]
    indep_ranks = [residuals[ct]['rank'] for ct in cell_types]

    fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.85))

    colors = [lineage_color(ct) for ct in cell_types]
    ax.scatter(joint_ranks, indep_ranks, c=colors, s=25,
               edgecolors='white', linewidths=0.3, zorder=3, alpha=0.85)

    # Identity line
    ax.plot([0, 36], [0, 36], color=C_LIGHTGRAY, linewidth=0.8,
            linestyle='--', zorder=1)

    # Label extremes — manual offsets for crowded labels, with leader lines.
    # Offsets re-tuned after the LABEL_EXPAND swap (HPC/HSC/etc → multi-word
    # forms): the four lower-left cluster labels (HSC/HPC/Stromal/Epithelial)
    # all extend right of their data points and stack vertically at staggered
    # dy values to avoid mutual collision. dx values increase for labels with
    # smaller-jr (leftmost) data points so the label left-edge clears the
    # y-axis tick label column (PDF x ≈ 43–51).
    label_offsets = {
        'endothelial cell': (3, -8),                # top-right cluster
        # Lower-left cluster (jr=1-4, ir=1-5): all labels stack to the upper
        # right of the cluster. dx tuned per data jr so label-left clears the
        # y-axis tick column (PDF x ≈ 43-51). dy tuned per data ir so labels
        # stack vertically without mutual overlap (each ~10pt apart in PDF y).
        'hematopoietic stem cell': (15, 13),        # jr=4, ir=1 — stack bottom
        'hematopoietic precursor cell': (20, 18),   # jr=3, ir=2
        'stromal cell': (60, 24),                   # jr=1, ir=3 — dx bumped 28→60 to clear MSC-adipose dot at PDF (74.9, 131.3)
        'epithelial cell': (24, 26),                # jr=2, ir=5 — stack top
    }
    for ct in cell_types:
        jr = residuals[ct]['joint_rank']
        ir = residuals[ct]['rank']
        if abs(jr - ir) > 8 or jr >= 34 or ir >= 34 or jr <= 2 or ir <= 2:
            dx, dy = label_offsets.get(ct, (3, 3))
            ha = 'right' if dx < 0 else 'left'
            label = LABEL_EXPAND.get(short_name(ct), short_name(ct))
            ax.annotate(label, (jr, ir),
                        fontsize=5, ha=ha, va='bottom',
                        xytext=(dx, dy), textcoords='offset points',
                        color=C_DARKGRAY,
                        arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))

    ax.set_xlabel('Joint-PCA rank', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('Independent-PCA rank', fontsize=FONT_SIZE_LABEL)
    ax.set_xlim(0, 36)
    ax.set_ylim(0, 36)
    ax.set_aspect('equal')

    # Format p-value with mathtext (matches local convention in
    # scripts/48_build_fig6_K12.py and Arial-safe — Unicode superscripts
    # like U+2074 are missing from Arial).
    _p = comp["spearman_p"]
    _mantissa, _exp = f'{_p:.1e}'.split('e')
    _p_str = rf'$p = {float(_mantissa):.1f} \times 10^{{{int(_exp)}}}$'

    ax.text(0.03, 0.97,
            f'ρ = {comp["spearman_rho"]:.3f}\n'
            f'{_p_str}',
            transform=ax.transAxes, fontsize=FONT_SIZE_ANNOT,
            ha='left', va='top', color=C_DARKGRAY,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor=C_LIGHTGRAY, alpha=0.9))

    clean_spine(ax)
    save_figure(fig, PANELS / 'figs1b_ranking_scatter')
    return fig


# ===================================================================
# MAIN
# ===================================================================
if __name__ == '__main__':
    print("=" * 60)
    print("PHASE 1: Generating DATA_ONLY panels")
    print("=" * 60)

    panels = {}
    panels['fig1b'] = fig1b_null_distribution()
    panels['fig1c'] = fig1c_lineage_stratified()
    panels['fig2a'] = fig2a_rigidity_ranking()
    panels['fig2b'] = fig2b_cellmarker()
    panels['fig2c'] = fig2c_cellcount()
    panels['fig4d'] = fig4d_replication_summary()
    panels['fig4e'] = fig4e_human_control()
    panels['fig5a'] = fig5a_macaque_primary()
    panels['fig5b'] = fig5b_macaque_sensitivity()
    panels['fig5c'] = fig5c_hepatocyte_rigidity()
    panels['figs1a'] = figs1a_indep_pca_null()
    panels['figs1b'] = figs1b_ranking_scatter()

    print("\n" + "=" * 60)
    print(f"Phase 1 complete: {len(panels)} panels generated")
    print(f"Output: {PANELS}")
    print("=" * 60)
