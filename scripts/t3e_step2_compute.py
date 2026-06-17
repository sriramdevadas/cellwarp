#!/usr/bin/env python3
"""
T3-E Step 2d-2h: Conservation score computation, Spearman correlation,
sensitivity analysis, and visualization.

Prerequisite: gene_sets.json, tss_lookup.json, rigidity_scores.csv,
and both bigWig files must already exist.

Biology: phastCons scores at promoter windows of cell-type-specific genes
measure how much purifying selection has maintained the regulatory architecture
of genes that define each cell type's geometric fingerprint. If rigidity
reflects upstream regulatory constraint, rigid cell types should have higher
promoter conservation at their identity-defining loci.

Math: For Option A, we take the top-50 genes by Procrustes loading magnitude
per cell type — the genes most responsible for each type's position in
Procrustes space — and compute the mean phastCons score over their promoter
windows. For Option B, we use the shared top-200 variance genes weighted by
each cell type's expression level.
"""

import json
import os
import time
import traceback

import numpy as np
import pandas as pd
import pyBigWig
from scipy import stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# ============================================================
# Configuration
# ============================================================
BASE = str(Path(__file__).resolve().parent.parent)
OUT = os.path.join(BASE, 'output/validation/t3e_chromatin')

BW_PATHS = {
    'placental_20way': os.path.join(BASE, 'data/ucsc/phastCons_placental.bw'),
    '100way_vertebrate': os.path.join(BASE, 'data/ucsc/phastCons100way.bw'),
}
WINDOW_KBS = [1, 2, 5]

# Pre-registered thresholds
RHO_POS = 0.50
RHO_NULL = 0.35

CLEAN_T3C = {
    'natural killer cell', 'B cell', 'plasma cell',
    'endothelial cell', 'CD8-positive, alpha-beta T cell'
}


def load_precomputed():
    """Load gene sets, TSS lookup, rigidity scores, and centroids."""
    with open(os.path.join(OUT, 'gene_sets.json')) as f:
        gs = json.load(f)
    with open(os.path.join(OUT, 'tss_lookup.json')) as f:
        tss = json.load(f)
    rigidity = pd.read_csv(os.path.join(OUT, 'rigidity_scores.csv'))
    centroids = pd.read_csv(
        os.path.join(BASE, 'output/phase2/scaled_35types/centroids_human_35.csv'),
        index_col=0)
    return gs, tss, rigidity, centroids


def compute_scores(gs, tss, centroids):
    """Step 2d: Compute phastCons conservation scores.

    For each cell type × track × option × window:
      Option A: unweighted mean over top-50 loading genes
      Option B: expression-weighted mean over top-200 identity genes
    """
    print("\n=== Step 2d: Computing conservation scores ===")

    scores_path = os.path.join(OUT, 'conservation_scores.csv')

    # Check for completed results
    if os.path.exists(scores_path):
        existing = pd.read_csv(scores_path)
        if len(existing) >= 420:
            print(f"  Found complete results: {len(existing)} rows")
            return existing

    option_a = gs['option_a']  # dict: cell_type -> list of ensembl_ids
    option_b = gs['option_b']  # list of ensembl_ids
    cell_types = sorted(option_a.keys())

    all_rows = []

    for track_name, bw_path in BW_PATHS.items():
        if not os.path.exists(bw_path):
            print(f"  SKIP: {bw_path} not found")
            continue

        print(f"\n  Track: {track_name}")
        bw = pyBigWig.open(bw_path)
        chrom_sizes = bw.chroms()

        for ci, ct in enumerate(cell_types):
            for opt_label in ['A', 'B']:
                gene_list = option_a[ct] if opt_label == 'A' else option_b

                for wkb in WINDOW_KBS:
                    half = wkb * 1000
                    scores = []
                    weights = []
                    n_found = 0
                    n_extracted = 0
                    n_failed = 0

                    for eid in gene_list:
                        if eid not in tss:
                            n_failed += 1
                            continue
                        n_found += 1
                        info = tss[eid]
                        chrom = info['chrom']
                        tss_pos = info['tss']

                        start = max(0, tss_pos - half)
                        end = tss_pos + half

                        if chrom not in chrom_sizes:
                            n_failed += 1
                            continue
                        end = min(end, chrom_sizes[chrom])

                        try:
                            vals = bw.stats(chrom, start, end, type="mean")
                            if vals and vals[0] is not None:
                                scores.append(vals[0])
                                n_extracted += 1

                                # For Option B: get expression weight
                                if opt_label == 'B' and eid in centroids.columns:
                                    if ct in centroids.index:
                                        w = centroids.loc[ct, eid]
                                        weights.append(max(w, 0))
                                    else:
                                        weights.append(0)
                            else:
                                n_failed += 1
                        except Exception:
                            n_failed += 1

                    # Compute mean score
                    if opt_label == 'A':
                        mean_score = np.mean(scores) if scores else np.nan
                    else:
                        # Expression-weighted mean
                        if scores and weights:
                            s = np.array(scores)
                            w = np.array(weights)
                            if w.sum() > 0:
                                mean_score = np.sum(s * w) / np.sum(w)
                            else:
                                mean_score = np.mean(scores)
                        else:
                            mean_score = np.nan

                    all_rows.append({
                        'cell_type': ct,
                        'track': track_name,
                        'option': opt_label,
                        'window_kb': wkb,
                        'mean_phastCons': mean_score,
                        'n_genes_found': n_found,
                        'n_windows_extracted': n_extracted,
                        'n_windows_failed': n_failed
                    })

            # Progress
            if (ci + 1) % 5 == 0 or ci == len(cell_types) - 1:
                print(f"    {ci+1}/{len(cell_types)} cell types processed")

        bw.close()

        # Incremental save after each track
        pd.DataFrame(all_rows).to_csv(scores_path, index=False)

    df = pd.DataFrame(all_rows)
    df.to_csv(scores_path, index=False)
    print(f"\n  Saved: {scores_path} ({len(df)} rows)")

    # Low coverage check
    opt_a = df[df['option'] == 'A']
    low = opt_a[opt_a['n_genes_found'] < 40]
    if len(low) > 0:
        print(f"\n  WARNING: {len(low)} Option A entries with <80% coverage:")
        for _, r in low.iterrows():
            print(f"    {r['cell_type']} | {r['track']} | {r['window_kb']}kb: "
                  f"n_found={r['n_genes_found']}")
    else:
        print("  All Option A entries have >=80% gene coverage")

    return df


def primary_spearman(scores_df, rigidity_df):
    """Step 2f: Primary Spearman correlation.

    Configuration: placental_20way / Option A / ±2kb.
    """
    print("\n=== Step 2f: Primary Spearman correlation ===")

    primary = scores_df[
        (scores_df['track'] == 'placental_20way') &
        (scores_df['option'] == 'A') &
        (scores_df['window_kb'] == 2)
    ].copy()

    merged = primary.merge(rigidity_df, on='cell_type', how='inner')
    merged = merged.dropna(subset=['mean_phastCons'])
    n = len(merged)

    print(f"  n = {n} cell types in primary analysis")
    if n < 35:
        missing = set(rigidity_df['cell_type']) - set(merged['cell_type'])
        if missing:
            print(f"  Dropped: {missing}")

    # Spearman: phastCons vs rigidity_score (higher = more rigid)
    rho, p = stats.spearmanr(merged['mean_phastCons'], merged['rigidity_score'])

    # 95% CI via Fisher z
    z = np.arctanh(rho)
    se = 1.0 / np.sqrt(n - 3) if n > 3 else float('inf')
    ci_lo, ci_hi = np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)

    # Determine threshold
    if rho >= RHO_POS and p < 0.05:
        conclusion = "POSITIVE"
    elif rho >= RHO_NULL:
        conclusion = "TREND"
    else:
        conclusion = "NULL_TRIGGERED"

    print(f"\n  *** PRIMARY RESULT ***")
    print(f"  Spearman ρ = {rho:.4f}")
    print(f"  p-value    = {p:.6f}")
    print(f"  n          = {n}")
    print(f"  95% CI     = [{ci_lo:.4f}, {ci_hi:.4f}]")
    print(f"  Conclusion = {conclusion}")

    result = {
        'rho': round(float(rho), 6),
        'p_value': round(float(p), 6),
        'n': int(n),
        'ci_95_lower': round(float(ci_lo), 4),
        'ci_95_upper': round(float(ci_hi), 4),
        'track': 'placental_20way',
        'option': 'A',
        'window_kb': 2,
        'conclusion': conclusion,
    }

    # Save
    merged.to_csv(os.path.join(OUT, 'rigidity_conservation_merged.csv'), index=False)
    with open(os.path.join(OUT, 'spearman_primary_result.json'), 'w') as f:
        json.dump(result, f, indent=2)

    return result, merged


def sensitivity_analysis(scores_df, rigidity_df, merged_primary):
    """Step 2g: Sensitivity grid + partial correlation + rank check."""
    print("\n=== Step 2g: Sensitivity analysis ===")

    rows = []
    for track in scores_df['track'].unique():
        for opt in ['A', 'B']:
            for wkb in WINDOW_KBS:
                sub = scores_df[
                    (scores_df['track'] == track) &
                    (scores_df['option'] == opt) &
                    (scores_df['window_kb'] == wkb)
                ]
                m = sub.merge(rigidity_df, on='cell_type', how='inner')
                m = m.dropna(subset=['mean_phastCons'])
                n = len(m)

                if n < 5:
                    rows.append({'track': track, 'option': opt, 'window_kb': wkb,
                                 'rho': np.nan, 'p': np.nan, 'n': n,
                                 'conclusion': 'INSUFFICIENT_N'})
                    continue

                rho, p = stats.spearmanr(m['mean_phastCons'], m['rigidity_score'])

                if rho >= RHO_POS and p < 0.05:
                    conc = "POSITIVE"
                elif rho >= RHO_NULL:
                    conc = "TREND"
                else:
                    conc = "NULL_TRIGGERED"

                rows.append({'track': track, 'option': opt, 'window_kb': wkb,
                             'rho': round(float(rho), 4), 'p': round(float(p), 6),
                             'n': n, 'conclusion': conc})

    sens = pd.DataFrame(rows)
    sens.to_csv(os.path.join(OUT, 'sensitivity_table.csv'), index=False)

    print("\n  Sensitivity table:")
    print(f"  {'Track':<22} {'Opt':>3} {'Win':>4} {'ρ':>8} {'p':>10} {'Conclusion'}")
    print("  " + "-" * 65)
    for _, r in sens.iterrows():
        rho_s = f"{r['rho']:.4f}" if pd.notna(r['rho']) else "  N/A "
        p_s = f"{r['p']:.6f}" if pd.notna(r['p']) else "    N/A   "
        print(f"  {r['track']:<22} {r['option']:>3} {r['window_kb']:>4} "
              f"{rho_s:>8} {p_s:>10} {r['conclusion']}")

    # Partial correlation controlling for expression level
    print("\n  --- Partial correlation (controlling for mean expression) ---")
    centroids = pd.read_csv(
        os.path.join(BASE, 'output/phase2/scaled_35types/centroids_human_35.csv'),
        index_col=0)
    mean_expr = centroids.mean(axis=1)
    mean_expr.name = 'mean_expression'

    mp = merged_primary.set_index('cell_type')
    mp = mp.join(mean_expr, how='left').dropna(
        subset=['mean_expression', 'mean_phastCons'])
    n_p = len(mp)

    if n_p >= 5:
        from scipy.stats import rankdata

        x = rankdata(mp['mean_phastCons'].values)
        y = rankdata(mp['rigidity_score'].values)
        z = rankdata(mp['mean_expression'].values)

        # Residualize ranks on expression
        cx = np.polyfit(z, x, 1)
        cy = np.polyfit(z, y, 1)
        xr = x - np.polyval(cx, z)
        yr = y - np.polyval(cy, z)

        rho_partial = np.corrcoef(xr, yr)[0, 1]
        t_stat = rho_partial * np.sqrt((n_p - 3) / (1 - rho_partial ** 2 + 1e-12))
        p_partial = 2 * stats.t.sf(np.abs(t_stat), df=n_p - 3)

        print(f"  Partial ρ = {rho_partial:.4f}, p = {p_partial:.6f}, n = {n_p}")

        # Also check mean expression vs phastCons correlation directly
        rho_expr_cons, p_expr_cons = stats.spearmanr(
            mp['mean_expression'], mp['mean_phastCons'])
        print(f"  Expression-conservation ρ = {rho_expr_cons:.4f}, p = {p_expr_cons:.6f}")

        rho_expr_rig, p_expr_rig = stats.spearmanr(
            mp['mean_expression'], mp['rigidity_score'])
        print(f"  Expression-rigidity ρ = {rho_expr_rig:.4f}, p = {p_expr_rig:.6f}")
    else:
        rho_partial = np.nan
        p_partial = np.nan
        print("  Insufficient data for partial correlation")

    # Rank-based check
    print("\n  --- Rank-based vs continuous check ---")
    rho_rank, p_rank = stats.spearmanr(
        merged_primary['mean_phastCons'], merged_primary['rigidity_rank'])
    print(f"  ρ(phastCons, rigidity_rank) = {rho_rank:.4f}, p = {p_rank:.6f}")
    print(f"  (Negative ρ = higher conservation → lower rank → more rigid)")

    return sens, rho_partial, p_partial


def make_plots(merged_primary, primary_result, sens_df):
    """Step 2h: Scatter plot and sensitivity heatmaps."""
    print("\n=== Step 2h: Generating plots ===")

    rho = primary_result['rho']
    p = primary_result['p_value']
    n = primary_result['n']
    conc = primary_result['conclusion']

    # --- Scatter: primary result ---
    fig, ax = plt.subplots(figsize=(10, 8))

    is_clean = merged_primary['cell_type'].isin(CLEAN_T3C)
    colors = ['#e74c3c' if c else '#3498db' for c in is_clean]

    ax.scatter(merged_primary['mean_phastCons'],
               merged_primary['rigidity_score'],
               c=colors, s=60, alpha=0.75, edgecolors='black', linewidth=0.5,
               zorder=3)

    for _, row in merged_primary.iterrows():
        label = row['cell_type']
        if len(label) > 30:
            label = label[:27] + '...'
        ax.annotate(label, (row['mean_phastCons'], row['rigidity_score']),
                    fontsize=5.5, alpha=0.85, xytext=(4, 3),
                    textcoords='offset points')

    ax.set_xlabel('Mean phastCons at promoters\n'
                  '(placental 20way, Option A top-50 loading genes, ±2kb)',
                  fontsize=10)
    ax.set_ylabel('Procrustes rigidity score (higher = more rigid)', fontsize=10)
    ax.set_title(f'T3-E: Regulatory Sequence Conservation vs Procrustes Rigidity\n'
                 f'Spearman ρ = {rho:.3f}, p = {p:.4f}, n = {n} | {conc}',
                 fontsize=11, fontweight='bold')

    from matplotlib.lines import Line2D
    legend = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#e74c3c',
               markersize=8, label='T3-C clean types'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#3498db',
               markersize=8, label='Other cell types'),
    ]
    ax.legend(handles=legend, loc='best', fontsize=8)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    scatter_path = os.path.join(OUT, 'scatter_primary.png')
    plt.savefig(scatter_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: scatter_primary.png")

    # --- Sensitivity heatmaps ---
    for track, fname in [('placental_20way', 'sensitivity_heatmap.png'),
                          ('100way_vertebrate', 'sensitivity_heatmap_100way.png')]:
        tdata = sens_df[sens_df['track'] == track]
        if len(tdata) == 0:
            continue

        pivot = tdata.pivot_table(index='option', columns='window_kb', values='rho')

        fig, ax = plt.subplots(figsize=(7, 3.5))
        vmax = max(abs(pivot.values[~np.isnan(pivot.values)]).max(), 0.3) \
            if not np.all(np.isnan(pivot.values)) else 0.5

        sns.heatmap(pivot, annot=True, fmt='.3f', cmap='RdBu_r',
                    center=0, vmin=-vmax, vmax=vmax,
                    linewidths=0.5, ax=ax, cbar_kws={'label': 'Spearman ρ'})

        track_label = 'Placental 20way' if 'placental' in track else '100-way Vertebrate'
        ax.set_title(f'T3-E Sensitivity: Spearman ρ\nTrack: {track_label}',
                     fontsize=11)
        ax.set_xlabel('Promoter window (±kb)', fontsize=10)
        ax.set_ylabel('Scoring option', fontsize=10)
        yticklabels = []
        for label in pivot.index:
            if label == 'A':
                yticklabels.append('A (loading genes)')
            else:
                yticklabels.append('B (expr-weighted)')
        ax.set_yticklabels(yticklabels, rotation=0, fontsize=9)

        plt.tight_layout()
        plt.savefig(os.path.join(OUT, fname), dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {fname}")


def write_summary(primary_result, sens_df, scores_df, rho_partial, p_partial):
    """Write t3e_step2_summary.md."""
    print("\n=== Writing summary ===")

    rho = primary_result['rho']
    p = primary_result['p_value']
    n = primary_result['n']
    ci_lo = primary_result['ci_95_lower']
    ci_hi = primary_result['ci_95_upper']
    conc = primary_result['conclusion']

    lines = [
        "# T3-E Step 2: Regulatory Sequence Conservation vs Procrustes Rigidity",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M')}",
        "",
        "## 1. UCSC Tracks",
        "",
        "**Primary:** phastCons20way — 20 species (17 primates + treeshrew + mouse + dog).",
        "Includes Mus musculus (mm10). Primate-dominated alignment; conservation scores",
        "primarily reflect primate constraint with mouse as an outgroup anchor.",
        "",
        "**Sensitivity:** phastCons100way — 100 vertebrates (mammals + birds + reptiles +",
        "amphibians + fish). Broadest taxonomic scope.",
        "",
        "Both tracks verified via md5 checksum.",
        "",
        "## 2. Primary Spearman Result",
        "",
        "| Metric | Value |",
        "|--------|-------|",
        f"| Spearman ρ | {rho:.4f} |",
        f"| p-value | {p:.6f} |",
        f"| n (cell types) | {n} |",
        f"| 95% CI | [{ci_lo:.4f}, {ci_hi:.4f}] |",
        "| Track | placental_20way (phastCons20way) |",
        "| Option | A (cell-type-specific top-50 loading genes) |",
        "| Window | ±2kb |",
        "",
    ]

    if conc == "POSITIVE":
        lines.append(f"**Pre-registered threshold: POSITIVE** — ρ = {rho:.3f} ≥ 0.50 AND "
                     f"p = {p:.4f} < 0.05. Regulatory sequence conservation at identity-gene "
                     "promoters predicts Procrustes rigidity.")
    elif conc == "TREND":
        lines.append(f"**Pre-registered threshold: TREND** — ρ = {rho:.3f} in [0.35, 0.50). "
                     "Directionally consistent but below positive threshold.")
    else:
        lines.append(f"**Pre-registered threshold: 8TH NULL TRIGGERED** — ρ = {rho:.3f} < 0.35. "
                     "Close chromatin/regulatory sequence as proximate mechanism.")
    lines.append("")

    # Sensitivity
    lines.extend([
        "## 3. Sensitivity Analysis",
        "",
        "| Track | Option | Window (kb) | ρ | p | Conclusion |",
        "|-------|--------|-------------|------|-------|-----------|",
    ])
    for _, r in sens_df.iterrows():
        rho_s = f"{r['rho']:.4f}" if pd.notna(r['rho']) else "N/A"
        p_s = f"{r['p']:.6f}" if pd.notna(r['p']) else "N/A"
        lines.append(f"| {r['track']} | {r['option']} | {r['window_kb']} | "
                     f"{rho_s} | {p_s} | {r['conclusion']} |")
    lines.append("")

    # Check consistency
    contradictions = sens_df[sens_df['conclusion'] != conc]
    contradictions = contradictions[contradictions['conclusion'] != 'INSUFFICIENT_N']
    if len(contradictions) > 0:
        lines.append(f"**{len(contradictions)} sensitivity results contradict primary "
                     f"conclusion ({conc}).** Review required.")
    else:
        lines.append("**All sensitivity results are consistent with primary conclusion.**")
    lines.append("")

    # Data quality
    opt_a = scores_df[scores_df['option'] == 'A']
    low = opt_a[opt_a['n_genes_found'] < 40]
    lines.extend([
        "## 4. Data Quality Flags",
        "",
        f"- Total conservation scores computed: {len(scores_df)}",
        f"- Genes resolved via Ensembl REST: 658/670 (12 MT-excluded, 0 failed)",
        f"- Option A entries with <80% gene coverage: {len(low)}",
    ])
    if len(low) > 0:
        for _, r in low.iterrows():
            lines.append(f"  - {r['cell_type']} ({r['track']}, {r['window_kb']}kb): "
                         f"n_found={r['n_genes_found']}")
    lines.append("")

    # Partial correlation
    lines.extend([
        "## 5. Partial Correlation",
        "",
    ])
    if not np.isnan(rho_partial):
        lines.append(f"Partial ρ (controlling for mean expression level) = {rho_partial:.4f}, "
                     f"p = {p_partial:.6f}")
        if abs(rho_partial - rho) > 0.1:
            lines.append("**Partial correlation substantially changes the result.** "
                         "The phastCons signal is partially confounded with expression level.")
        else:
            lines.append("Partial correlation does not substantially change the conclusion.")
    else:
        lines.append("Partial correlation not computed (insufficient data).")
    lines.append("")

    # Files
    lines.extend([
        "## Files Generated",
        "",
        "- `conservation_scores.csv` — All conservation scores",
        "- `rigidity_conservation_merged.csv` — Primary analysis merged data (35 rows)",
        "- `spearman_primary_result.json` — Primary Spearman result",
        "- `sensitivity_table.csv` — Sensitivity grid (12 rows)",
        "- `identity_gene_tss_hg38.bed` — TSS coordinates (BED6)",
        "- `scatter_primary.png` — Primary result scatter plot",
        "- `sensitivity_heatmap.png` — Sensitivity heatmap (placental)",
        "- `sensitivity_heatmap_100way.png` — Sensitivity heatmap (100way)",
        "- `download_log.txt` — BigWig download sizes and checksums",
        "- `ucsc_track_availability.txt` — UCSC track availability report",
    ])

    summary_path = os.path.join(OUT, 't3e_step2_summary.md')
    with open(summary_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print(f"  Saved: {summary_path}")


# ============================================================
# Main
# ============================================================

def main():
    print("=" * 70)
    print("T3-E Step 2: Regulatory Sequence Conservation vs Procrustes Rigidity")
    print("=" * 70)

    gs, tss, rigidity, centroids = load_precomputed()

    # Step 2d
    scores_df = compute_scores(gs, tss, centroids)

    # Step 2f: Primary Spearman
    primary_result, merged = primary_spearman(scores_df, rigidity)

    # Step 2g: Sensitivity
    sens_df, rho_partial, p_partial = sensitivity_analysis(
        scores_df, rigidity, merged)

    # Step 2h: Plots
    make_plots(merged, primary_result, sens_df)

    # Summary
    write_summary(primary_result, sens_df, scores_df, rho_partial, p_partial)

    print("\n" + "=" * 70)
    print("T3-E Step 2 COMPLETE")
    print("=" * 70)


if __name__ == '__main__':
    main()
