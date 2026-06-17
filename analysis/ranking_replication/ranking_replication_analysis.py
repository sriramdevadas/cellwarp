"""
Ranking Replication Analysis
============================
Computes Spearman rank correlations between per-cell-type Procrustes
residual magnitudes (rigidity scores) from the primary 35-type analysis
and three independent replications: Sun2023, PanSci, and CellHint.

Biology: Higher Procrustes residuals = more diverged between species
(less conserved geometry). Lower residuals = more rigid/conserved.

Math: Spearman ρ on shared cell types tests whether the ordinal ranking
of cell-type rigidity is reproducible across independent datasets.
"""

import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

OUT_DIR = Path(__file__).parent
DATA_DIR = Path(__file__).resolve().parent.parent.parent

# ── Load data ──────────────────────────────────────────────────────────
primary = pd.read_csv(DATA_DIR / 'output/phase2/scaled_35types/residuals_ranked.csv')
sun2023 = pd.read_csv(DATA_DIR / 'output/validation/sun2023_replication_expanded/ranking_comparison.csv')
pansci = pd.read_csv(DATA_DIR / 'output/validation/pansci_replication/ranking_comparison.csv')
cellhint = pd.read_csv(DATA_DIR / 'output/validation/cellhint_replication/ranking_comparison.csv')

# Build primary lookup: cell_type -> (rank, residual)
primary_lookup = {}
for _, row in primary.iterrows():
    primary_lookup[row['cell_type']] = {
        'rank': int(row['rank']),
        'residual': row['residual_magnitude']
    }

def analyze_replication(name, df, resid_col):
    """Compute Spearman correlation and per-type rank comparison."""
    # Get shared cell types
    shared = [ct for ct in df['cell_type'] if ct in primary_lookup]
    n = len(shared)

    # Build comparison dataframe
    rows = []
    for ct in shared:
        primary_resid = primary_lookup[ct]['residual']
        primary_rank = primary_lookup[ct]['rank']
        repl_resid = df.loc[df['cell_type'] == ct, resid_col].values[0]
        rows.append({
            'cell_type': ct,
            'primary_residual': primary_resid,
            'primary_rank': primary_rank,
            'replication_residual': repl_resid,
        })

    comp = pd.DataFrame(rows)

    # Compute replication ranks (1 = highest residual = most diverged)
    comp['replication_rank'] = comp['replication_residual'].rank(ascending=False).astype(int)

    # Re-rank primary within shared set only
    comp['primary_rank_shared'] = comp['primary_residual'].rank(ascending=False).astype(int)

    # Spearman on residual magnitudes directly
    rho, pval = spearmanr(comp['primary_residual'], comp['replication_residual'])

    # Rank shifts
    comp['rank_shift'] = abs(comp['primary_rank_shared'] - comp['replication_rank'])
    comp = comp.sort_values('primary_rank_shared')

    return {
        'name': name,
        'n_shared': n,
        'rho': rho,
        'pval': pval,
        'comparison': comp,
        'shared_types': shared
    }

# ── Run analyses ───────────────────────────────────────────────────────
results = {
    'Sun2023': analyze_replication('Sun2023 (15 types)', sun2023, 'sun2023_residual'),
    'PanSci': analyze_replication('PanSci (16 types)', pansci, 'pansci_residual'),
    'CellHint': analyze_replication('CellHint (15 types)', cellhint, 'cellhint_residual'),
}

# ── Scatter plots ──────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

for idx, (key, res) in enumerate(results.items()):
    ax = axes[idx]
    comp = res['comparison']

    ax.scatter(comp['primary_rank_shared'], comp['replication_rank'],
               s=80, c='#2563eb', alpha=0.8, edgecolors='white', linewidth=0.5)

    # Label points
    for _, row in comp.iterrows():
        # Shorten long names
        label = row['cell_type']
        if len(label) > 20:
            label = label[:18] + '…'
        ax.annotate(label,
                    (row['primary_rank_shared'], row['replication_rank']),
                    fontsize=6.5, ha='left', va='bottom',
                    xytext=(3, 3), textcoords='offset points')

    # Perfect agreement line
    max_rank = res['n_shared']
    ax.plot([1, max_rank], [1, max_rank], 'k--', alpha=0.3, linewidth=1)

    ax.set_xlabel('Primary rank (within shared types)', fontsize=10)
    ax.set_ylabel(f'{key} rank', fontsize=10)
    ax.set_title(f'{res["name"]}\nρ = {res["rho"]:.3f}, p = {res["pval"]:.3f}, n = {res["n_shared"]}',
                 fontsize=11)
    ax.set_xlim(0, max_rank + 1)
    ax.set_ylim(0, max_rank + 1)
    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.invert_xaxis()
    # Rank 1 = most diverged, so top-left corner = most diverged in both

plt.tight_layout()
plt.savefig(OUT_DIR / 'ranking_scatter_all.png', dpi=200, bbox_inches='tight')
plt.close()
print(f"Saved scatter plot to {OUT_DIR / 'ranking_scatter_all.png'}")

# Individual scatter plots
for key, res in results.items():
    fig, ax = plt.subplots(figsize=(7, 7))
    comp = res['comparison']

    ax.scatter(comp['primary_rank_shared'], comp['replication_rank'],
               s=100, c='#2563eb', alpha=0.8, edgecolors='white', linewidth=0.5)

    for _, row in comp.iterrows():
        label = row['cell_type']
        ax.annotate(label,
                    (row['primary_rank_shared'], row['replication_rank']),
                    fontsize=7.5, ha='left', va='bottom',
                    xytext=(4, 4), textcoords='offset points')

    max_rank = res['n_shared']
    ax.plot([1, max_rank], [1, max_rank], 'k--', alpha=0.3, linewidth=1)
    ax.set_xlabel('Primary analysis rank (1 = most diverged)', fontsize=11)
    ax.set_ylabel(f'{key} replication rank', fontsize=11)
    ax.set_title(f'Rigidity ranking: Primary vs {res["name"]}\n'
                 f'Spearman ρ = {res["rho"]:.3f}, p = {res["pval"]:.3f}, n = {res["n_shared"]}',
                 fontsize=12)
    ax.set_xlim(0, max_rank + 1)
    ax.set_ylim(0, max_rank + 1)
    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.invert_xaxis()

    plt.tight_layout()
    fname = f'ranking_scatter_{key.lower()}.png'
    plt.savefig(OUT_DIR / fname, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {fname}")

# ── Combined stability analysis ────────────────────────────────────────
# For each cell type, collect rank shifts across all replications it appears in
all_types = set()
for res in results.values():
    all_types.update(res['comparison']['cell_type'].tolist())

stability = []
for ct in sorted(all_types):
    rank_shifts = []
    appearances = []
    for key, res in results.items():
        comp = res['comparison']
        row = comp[comp['cell_type'] == ct]
        if not row.empty:
            rank_shifts.append(row['rank_shift'].values[0])
            appearances.append(key)

    if rank_shifts:
        stability.append({
            'cell_type': ct,
            'n_datasets': len(rank_shifts),
            'mean_rank_shift': np.mean(rank_shifts),
            'max_rank_shift': max(rank_shifts),
            'datasets': ', '.join(appearances),
            'primary_residual': primary_lookup.get(ct, {}).get('residual', np.nan),
            'primary_rank_35': primary_lookup.get(ct, {}).get('rank', np.nan),
        })

stab_df = pd.DataFrame(stability).sort_values('mean_rank_shift')

# Classify stability
def classify(row):
    if row['n_datasets'] < 2:
        return 'insufficient data'
    if row['mean_rank_shift'] <= 2.0:
        return 'stable'
    elif row['mean_rank_shift'] <= 4.0:
        return 'moderate'
    else:
        return 'volatile'

stab_df['classification'] = stab_df.apply(classify, axis=1)
stab_df.to_csv(OUT_DIR / 'stability_classification.csv', index=False)

# ── Print results ──────────────────────────────────────────────────────
print("\n" + "="*70)
print("RANKING REPLICATION ANALYSIS — SUMMARY")
print("="*70)

for key, res in results.items():
    print(f"\n--- {res['name']} ---")
    print(f"  Shared cell types: {res['n_shared']}")
    print(f"  Spearman ρ = {res['rho']:.4f}")
    print(f"  p-value   = {res['pval']:.4f}")
    print(f"\n  Per-type comparison (ranked by primary):")
    comp = res['comparison'].sort_values('primary_rank_shared')
    for _, row in comp.iterrows():
        shift = int(row['rank_shift'])
        arrow = '→' if shift <= 2 else '↕' if shift <= 4 else '⇕'
        print(f"    {row['cell_type']:45s}  "
              f"P:{int(row['primary_rank_shared']):2d}  "
              f"R:{int(row['replication_rank']):2d}  "
              f"shift:{shift:2d} {arrow}")

print(f"\n--- Combined Stability Analysis ---")
print(f"{'Cell Type':45s} {'N':>3s} {'Mean Shift':>10s} {'Class':>12s}")
print("-" * 75)
for _, row in stab_df.iterrows():
    print(f"{row['cell_type']:45s} {row['n_datasets']:3.0f} "
          f"{row['mean_rank_shift']:10.1f} {row['classification']:>12s}")

# ── Generate markdown report ───────────────────────────────────────────
report_lines = []
report_lines.append("# Ranking Replication Analysis")
report_lines.append("")
report_lines.append("**Question:** Do per-cell-type rigidity rankings (Procrustes residual")
report_lines.append("magnitudes) replicate across independent datasets?")
report_lines.append("")
report_lines.append("**Method:** For each replication dataset, we identified cell types shared")
report_lines.append("with the primary 35-type analysis, re-ranked within the shared subset, and")
report_lines.append("computed Spearman rank correlation (ρ) between primary and replication rankings.")
report_lines.append("")

# Summary table
report_lines.append("## Per-Dataset Correlation Summary")
report_lines.append("")
report_lines.append("| Dataset | N shared types | Spearman ρ | p-value | Interpretation |")
report_lines.append("|---------|---------------|-----------|---------|----------------|")
for key, res in results.items():
    if res['pval'] < 0.05:
        interp = f"{'Positive' if res['rho'] > 0 else 'Negative'}, significant"
    else:
        interp = "Not significant"
    report_lines.append(f"| {res['name']} | {res['n_shared']} | {res['rho']:.3f} | {res['pval']:.4f} | {interp} |")
report_lines.append("")

# Per-dataset detail
for key, res in results.items():
    report_lines.append(f"### {res['name']}")
    report_lines.append("")
    report_lines.append(f"Shared cell types: **{res['n_shared']}** | Spearman ρ = **{res['rho']:.3f}** | p = **{res['pval']:.4f}**")
    report_lines.append("")
    report_lines.append(f"![Scatter plot](ranking_scatter_{key.lower()}.png)")
    report_lines.append("")
    report_lines.append("| Cell Type | Primary Rank | Replication Rank | Rank Shift |")
    report_lines.append("|-----------|:------------|:----------------|:----------|")
    comp = res['comparison'].sort_values('primary_rank_shared')
    for _, row in comp.iterrows():
        shift = int(row['rank_shift'])
        flag = ' **' if shift >= 5 else ''
        endflag = '**' if shift >= 5 else ''
        report_lines.append(f"| {row['cell_type']} | {int(row['primary_rank_shared'])} | "
                          f"{int(row['replication_rank'])} | {flag}{shift}{endflag} |")

    # Largest rank shifts
    top_shifts = comp.nlargest(3, 'rank_shift')
    report_lines.append("")
    report_lines.append("**Largest rank shifts:**")
    for _, row in top_shifts.iterrows():
        report_lines.append(f"- {row['cell_type']}: shifted {int(row['rank_shift'])} positions "
                          f"(primary #{int(row['primary_rank_shared'])} → replication #{int(row['replication_rank'])})")
    report_lines.append("")

# Combined stability
report_lines.append("## Combined Stability Analysis")
report_lines.append("")
report_lines.append("Cell types appearing in ≥2 replications are classified by mean absolute rank shift:")
report_lines.append("- **Stable** (≤2.0): ranking is reproducible across datasets")
report_lines.append("- **Moderate** (2.1–4.0): some ranking variability")
report_lines.append("- **Volatile** (>4.0): ranking is dataset-dependent")
report_lines.append("")
report_lines.append("| Cell Type | N datasets | Mean Rank Shift | Max Shift | Classification | Primary Rank (of 35) |")
report_lines.append("|-----------|:----------|:---------------|:---------|:--------------|:---------------------|")

for _, row in stab_df.iterrows():
    report_lines.append(f"| {row['cell_type']} | {int(row['n_datasets'])} | "
                      f"{row['mean_rank_shift']:.1f} | {int(row['max_rank_shift'])} | "
                      f"{row['classification']} | {int(row['primary_rank_35'])} |")

report_lines.append("")

# Key findings
stable_types = stab_df[(stab_df['classification'] == 'stable') & (stab_df['n_datasets'] >= 2)]
moderate_types = stab_df[(stab_df['classification'] == 'moderate') & (stab_df['n_datasets'] >= 2)]
volatile_types = stab_df[(stab_df['classification'] == 'volatile') & (stab_df['n_datasets'] >= 2)]

report_lines.append("## Key Findings")
report_lines.append("")

report_lines.append("### Correlation Results")
rhos = [res['rho'] for res in results.values()]
report_lines.append(f"- **Sun2023** (ρ = {results['Sun2023']['rho']:.3f}): "
                   f"{'Significant' if results['Sun2023']['pval'] < 0.05 else 'Not significant'} "
                   f"{'positive' if results['Sun2023']['rho'] > 0 else 'negative'} correlation")
report_lines.append(f"- **PanSci** (ρ = {results['PanSci']['rho']:.3f}): "
                   f"{'Significant' if results['PanSci']['pval'] < 0.05 else 'Not significant'} "
                   f"{'positive' if results['PanSci']['rho'] > 0 else 'negative'} correlation")
report_lines.append(f"- **CellHint** (ρ = {results['CellHint']['rho']:.3f}): "
                   f"{'Significant' if results['CellHint']['pval'] < 0.05 else 'Not significant'} "
                   f"{'positive' if results['CellHint']['rho'] > 0 else 'negative'} correlation")
report_lines.append("")

report_lines.append("### Stability Summary")
if len(stable_types) > 0:
    report_lines.append(f"- **Stable types** ({len(stable_types)}): " +
                       ', '.join(stable_types['cell_type'].tolist()))
if len(moderate_types) > 0:
    report_lines.append(f"- **Moderate types** ({len(moderate_types)}): " +
                       ', '.join(moderate_types['cell_type'].tolist()))
if len(volatile_types) > 0:
    report_lines.append(f"- **Volatile types** ({len(volatile_types)}): " +
                       ', '.join(volatile_types['cell_type'].tolist()))
report_lines.append("")

report_lines.append("### Interpretation")
report_lines.append("")
pos_count = sum(1 for r in rhos if r > 0)
neg_count = sum(1 for r in rhos if r < 0)
sig_count = sum(1 for res in results.values() if res['pval'] < 0.05)

if sig_count >= 2 and pos_count >= 2:
    report_lines.append("Per-type rigidity rankings show significant positive replication in the "
                       "majority of independent datasets, indicating that which cell types are "
                       "most/least conserved is a reproducible biological signal, not an artifact "
                       "of dataset-specific batch effects.")
elif pos_count >= 2:
    report_lines.append("Per-type rigidity rankings trend positive across replications, though "
                       "statistical significance varies with sample size. The direction of the "
                       "correlation is consistent, suggesting an underlying biological signal "
                       "modulated by dataset-specific factors (tissue representation, cell counts, "
                       "protocol differences).")
else:
    report_lines.append("Per-type rigidity rankings show mixed replication, with some datasets "
                       "agreeing and others diverging. This suggests that while global geometric "
                       "coherence is robust, the fine-grained ranking of individual cell types "
                       "is sensitive to dataset composition and processing pipeline differences.")

report_lines.append("")
report_lines.append("---")
report_lines.append(f"*Generated: 2026-04-05*")

with open(OUT_DIR / 'RANKING_REPLICATION_ANALYSIS.md', 'w') as f:
    f.write('\n'.join(report_lines))

print(f"\nSaved report to {OUT_DIR / 'RANKING_REPLICATION_ANALYSIS.md'}")
print("Done.")
