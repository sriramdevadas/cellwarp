"""
Cross-Reference: Bootstrap Stability vs Cross-Atlas Ranking Consistency
=======================================================================
Checks whether cell types that are bootstrap-stable (narrow CI within
the primary analysis) are also consistent across independent replications
(low rank shift across atlases).

Biology: If a cell type's rigidity rank is biologically meaningful,
it should be stable both under resampling (bootstrap) and when measured
in completely different atlases (cross-atlas replication).

Math: For each cell type present in multiple replications, we compute
the mean absolute rank shift from the primary analysis (within each
replication's shared type subset). We then correlate this cross-atlas
instability with the bootstrap CI width.
"""

import pandas as pd
import numpy as np
from scipy.stats import spearmanr, pearsonr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent.parent
OUT = Path(__file__).resolve().parent

# ── STEP 1: Load bootstrap results ───────────────────────────────────
bootstrap = pd.read_csv(BASE / 'analysis/bootstrap_rankings/bootstrap_summary.csv')
print(f"Bootstrap: {len(bootstrap)} types loaded")

# Build bootstrap lookup
boot_lookup = {}
for _, row in bootstrap.iterrows():
    boot_lookup[row['cell_type']] = {
        'median_rank': row['median_rank'],
        'ci_width': row['ci_width'],
        'category': row['category'],
        'ci_lower': row['ci_lower'],
        'ci_upper': row['ci_upper'],
    }

# ── STEP 2: Load primary analysis ────────────────────────────────────
primary = pd.read_csv(BASE / 'output/phase2/scaled_35types/residuals_ranked.csv')
primary_lookup = {}
for _, row in primary.iterrows():
    primary_lookup[row['cell_type']] = {
        'rank': int(row['rank']),
        'residual': row['residual_magnitude'],
    }

# ── STEP 3: Load all replication data ────────────────────────────────
# Each replication: dict of cell_type -> residual_magnitude
# We'll compute ranks within each replication's shared type subset

def load_replication_residuals(filepath, resid_col, ct_col='cell_type'):
    """Load replication data and return {cell_type: residual}."""
    df = pd.read_csv(filepath)
    return dict(zip(df[ct_col], df[resid_col]))

# Sun2023
sun2023_resids = load_replication_residuals(
    BASE / 'output/validation/sun2023_replication_expanded/ranking_comparison.csv',
    'sun2023_residual')

# PanSci
pansci_resids = load_replication_residuals(
    BASE / 'output/validation/pansci_replication/ranking_comparison.csv',
    'pansci_residual')

# CellHint
cellhint_resids = load_replication_residuals(
    BASE / 'output/validation/cellhint_replication/ranking_comparison.csv',
    'cellhint_residual')

# CellHint harmonized (12-type restricted set)
cellhint_harm = pd.read_csv(BASE / 'analysis/harmonized_replication/harmonized_residuals_cellhint.csv')
cellhint_harm_resids = dict(zip(cellhint_harm['cell_type'], cellhint_harm['residual_magnitude']))

# Pan-Census (22 types)
census = pd.read_csv(BASE / 'analysis/census_replication/ranking_comparison.csv')
census_resids = dict(zip(census['cell_type'], census['replication_residual']))

# Macaque (12 types, Qu-only canonical; was deprecated 20-type RIRA-mixed)
import json as _json
_qu12 = _json.load(open(BASE / 'output/macaque_pipeline/reconstruction_qu12_results.json'))
macaque_resids = {r['type']: r['magnitude'] for r in _qu12['per_type_residuals_ranked']}

# Mouse lemur (15 types)
lemur = pd.read_csv(BASE / 'analysis/mouse_lemur/per_type_residuals.csv')
lemur_resids = dict(zip(lemur['cell_type'], lemur['residual_magnitude']))


def compute_ranks_and_shifts(repl_resids, primary_lookup):
    """
    For a replication's residuals, compute:
    - primary_rank_within_subset: re-rank primary residuals among shared types
    - replication_rank: rank replication residuals (1 = highest = most flexible)
    - rank_shift: |primary_rank_within_subset - replication_rank|

    Returns dict: cell_type -> {repl_rank, primary_subset_rank, rank_shift}
    """
    shared = [ct for ct in repl_resids if ct in primary_lookup]
    if not shared:
        return {}

    rows = []
    for ct in shared:
        rows.append({
            'cell_type': ct,
            'primary_residual': primary_lookup[ct]['residual'],
            'repl_residual': repl_resids[ct],
        })
    df = pd.DataFrame(rows)

    # Rank 1 = highest residual = most flexible (descending)
    df['primary_subset_rank'] = df['primary_residual'].rank(ascending=False).astype(int)
    df['repl_rank'] = df['repl_residual'].rank(ascending=False).astype(int)
    df['rank_shift'] = abs(df['primary_subset_rank'] - df['repl_rank'])

    result = {}
    for _, row in df.iterrows():
        result[row['cell_type']] = {
            'repl_rank': int(row['repl_rank']),
            'n_types': len(shared),
            'primary_subset_rank': int(row['primary_subset_rank']),
            'rank_shift': int(row['rank_shift']),
        }
    return result


# Compute for each replication
replications = {
    'Sun2023': compute_ranks_and_shifts(sun2023_resids, primary_lookup),
    'PanSci': compute_ranks_and_shifts(pansci_resids, primary_lookup),
    'CellHint': compute_ranks_and_shifts(cellhint_resids, primary_lookup),
    'CellHint_harmonized': compute_ranks_and_shifts(cellhint_harm_resids, primary_lookup),
    'Pan_Census': compute_ranks_and_shifts(census_resids, primary_lookup),
    'Macaque': compute_ranks_and_shifts(macaque_resids, primary_lookup),
    'Mouse_lemur': compute_ranks_and_shifts(lemur_resids, primary_lookup),
}

repl_names = list(replications.keys())
print(f"\nReplications loaded:")
for name, data in replications.items():
    print(f"  {name}: {len(data)} types")

# ── STEP 4: Build master table ────────────────────────────────────────
all_types = sorted(primary_lookup.keys())

master_rows = []
for ct in all_types:
    row = {
        'cell_type': ct,
        'primary_rank': primary_lookup[ct]['rank'],
        'bootstrap_median_rank': boot_lookup.get(ct, {}).get('median_rank', np.nan),
        'bootstrap_CI_width': boot_lookup.get(ct, {}).get('ci_width', np.nan),
        'bootstrap_CI_lower': boot_lookup.get(ct, {}).get('ci_lower', np.nan),
        'bootstrap_CI_upper': boot_lookup.get(ct, {}).get('ci_upper', np.nan),
        'bootstrap_category': boot_lookup.get(ct, {}).get('category', ''),
    }

    # Add replication ranks and rank shifts
    shifts = []
    n_present = 0
    for rname in repl_names:
        rdata = replications[rname].get(ct)
        if rdata:
            row[f'{rname}_rank'] = rdata['repl_rank']
            row[f'{rname}_n_types'] = rdata['n_types']
            row[f'{rname}_rank_shift'] = rdata['rank_shift']
            shifts.append(rdata['rank_shift'])
            n_present += 1
        else:
            row[f'{rname}_rank'] = np.nan
            row[f'{rname}_n_types'] = np.nan
            row[f'{rname}_rank_shift'] = np.nan

    row['n_replications_present'] = n_present
    row['mean_rank_shift'] = np.mean(shifts) if shifts else np.nan
    row['sd_rank_shift'] = np.std(shifts, ddof=1) if len(shifts) > 1 else np.nan

    master_rows.append(row)

master = pd.DataFrame(master_rows).sort_values('primary_rank')

# Save master table (select key columns for the CSV)
csv_cols = ['cell_type', 'primary_rank', 'bootstrap_median_rank',
            'bootstrap_CI_width', 'bootstrap_category']
for rname in repl_names:
    csv_cols.append(f'{rname}_rank')
csv_cols += ['n_replications_present', 'mean_rank_shift', 'sd_rank_shift']

master[csv_cols].to_csv(OUT / 'master_ranking_table.csv', index=False)
print(f"\nSaved master_ranking_table.csv ({len(master)} types)")

# ── STEP 5: Identify convergent types ─────────────────────────────────

# For types present in ≥2 replications, compute quartiles of mean_rank_shift
multi_repl = master[master['n_replications_present'] >= 2].copy()

if len(multi_repl) > 0:
    sd_q25 = multi_repl['mean_rank_shift'].quantile(0.25)
    sd_q75 = multi_repl['mean_rank_shift'].quantile(0.75)
    print(f"\nMean rank shift quartiles (types in ≥2 replications):")
    print(f"  Q25 = {sd_q25:.2f}")
    print(f"  Q75 = {sd_q75:.2f}")

# (a) STABLE in bootstrap AND consistent across atlases
stable_boot = master[master['bootstrap_CI_width'] <= 3].copy()
stable_boot_consistent = stable_boot[
    (stable_boot['n_replications_present'] >= 2) &
    (stable_boot['mean_rank_shift'] <= sd_q25)
].copy()

# (b) STABLE in bootstrap BUT inconsistent across atlases
stable_boot_inconsistent = stable_boot[
    (stable_boot['n_replications_present'] >= 2) &
    (stable_boot['mean_rank_shift'] >= sd_q75)
].copy()

# (c) Specific types with zero rank shift in Pan-Census
zero_shift_census = ['myeloid leukocyte', 'plasma cell', 'B cell']

# (d) CD8+ T cell special check
cd8_name = 'CD8-positive, alpha-beta T cell'

print(f"\n{'='*70}")
print("CONVERGENT TYPES ANALYSIS")
print(f"{'='*70}")

print(f"\n(a) Bootstrap-stable (CI ≤ 3) AND cross-atlas consistent (mean shift ≤ Q25={sd_q25:.1f}):")
if len(stable_boot_consistent) > 0:
    for _, r in stable_boot_consistent.iterrows():
        print(f"  {r['cell_type']:45s}  CI={r['bootstrap_CI_width']:.0f}  "
              f"mean_shift={r['mean_rank_shift']:.1f}  cat={r['bootstrap_category']}")
else:
    print("  (none)")

print(f"\n(b) Bootstrap-stable (CI ≤ 3) BUT cross-atlas inconsistent (mean shift ≥ Q75={sd_q75:.1f}):")
if len(stable_boot_inconsistent) > 0:
    for _, r in stable_boot_inconsistent.iterrows():
        print(f"  {r['cell_type']:45s}  CI={r['bootstrap_CI_width']:.0f}  "
              f"mean_shift={r['mean_rank_shift']:.1f}  cat={r['bootstrap_category']}")
else:
    print("  (none)")

print(f"\n(c) Pan-Census zero-shift types — bootstrap CIs:")
for ct in zero_shift_census:
    bci = boot_lookup.get(ct, {})
    census_data = replications['Pan_Census'].get(ct)
    census_shift = census_data['rank_shift'] if census_data else 'N/A'
    print(f"  {ct:45s}  CI=[{bci.get('ci_lower','?')}, {bci.get('ci_upper','?')}]  "
          f"width={bci.get('ci_width','?')}  cat={bci.get('category','?')}  "
          f"census_shift={census_shift}")

print(f"\n(d) CD8+ T cell across replications:")
cd8_boot = boot_lookup.get(cd8_name, {})
print(f"  Bootstrap: median_rank={cd8_boot.get('median_rank')}, "
      f"CI=[{cd8_boot.get('ci_lower')}, {cd8_boot.get('ci_upper')}], "
      f"width={cd8_boot.get('ci_width')}, cat={cd8_boot.get('category')}")
for rname in repl_names:
    rdata = replications[rname].get(cd8_name)
    if rdata:
        print(f"  {rname:25s}: rank {rdata['repl_rank']}/{rdata['n_types']}  "
              f"(primary subset rank: {rdata['primary_subset_rank']}, shift: {rdata['rank_shift']})")
    else:
        print(f"  {rname:25s}: not present")

# ── STEP 6: Scatter plot — bootstrap CI width vs cross-atlas SD ───────
plot_data = master[master['n_replications_present'] >= 2].copy()

fig, ax = plt.subplots(figsize=(10, 8))

# Color by bootstrap category
colors = {
    'STABLE_RIGID': '#1e40af',     # dark blue
    'STABLE_MIDDLE': '#7c3aed',    # purple
    'STABLE_FLEXIBLE': '#dc2626',  # red
}

for _, row in plot_data.iterrows():
    color = colors.get(row['bootstrap_category'], '#6b7280')
    ax.scatter(row['bootstrap_CI_width'], row['mean_rank_shift'],
               s=80, c=color, alpha=0.8, edgecolors='white', linewidth=0.5,
               zorder=3)

    # Label each point
    label = row['cell_type']
    if 'CD8-positive' in label:
        label = 'CD8+ T'
    elif 'CD4-positive' in label:
        label = 'CD4+ T'
    elif 'mesenchymal stem cell of' in label:
        label = 'MSC adipose'
    elif label == 'mesenchymal stem cell':
        label = 'MSC'
    elif label == 'luminal epithelial cell of mammary gland':
        label = 'luminal epithelial'
    elif label == 'large intestine goblet cell':
        label = 'goblet cell'
    elif label == 'enterocyte of epithelium of large intestine':
        label = 'enterocyte'
    elif label == 'fibroblast of cardiac tissue':
        label = 'cardiac fibroblast'
    elif label == 'bladder urothelial cell':
        label = 'urothelial'
    elif label == 'pancreatic ductal cell':
        label = 'pancreatic ductal'
    elif label == 'pancreatic acinar cell':
        label = 'pancreatic acinar'
    elif label == 'non-classical monocyte':
        label = 'NC monocyte'
    elif label == 'intermediate monocyte':
        label = 'int. monocyte'
    elif label == 'classical monocyte':
        label = 'class. monocyte'
    elif label == 'mature NK T cell':
        label = 'NKT'
    elif label == 'myeloid dendritic cell':
        label = 'mDC'
    elif label == 'natural killer cell':
        label = 'NK'
    elif label == 'myeloid leukocyte':
        label = 'myeloid leuk.'
    elif label == 'hematopoietic precursor cell':
        label = 'HPC'
    elif label == 'hematopoietic stem cell':
        label = 'HSC'
    elif label == 'adventitial cell':
        label = 'adventitial'

    ax.annotate(label, (row['bootstrap_CI_width'], row['mean_rank_shift']),
                fontsize=7, ha='left', va='bottom',
                xytext=(4, 4), textcoords='offset points')

# Correlation
valid = plot_data.dropna(subset=['bootstrap_CI_width', 'mean_rank_shift'])
if len(valid) >= 3:
    rho, pval = spearmanr(valid['bootstrap_CI_width'], valid['mean_rank_shift'])
    r_pearson, p_pearson = pearsonr(valid['bootstrap_CI_width'], valid['mean_rank_shift'])
    ax.set_title(f'Bootstrap Stability vs Cross-Atlas Ranking Consistency\n'
                 f'Spearman ρ = {rho:.3f} (p = {pval:.3f}), '
                 f'Pearson r = {r_pearson:.3f} (p = {p_pearson:.3f})\n'
                 f'n = {len(valid)} types in ≥2 replications',
                 fontsize=11)

ax.set_xlabel('Bootstrap CI width (within-atlas stability)\n'
              '0 = perfectly stable rank, larger = more variable', fontsize=10)
ax.set_ylabel('Mean absolute rank shift across replications\n'
              '(cross-atlas instability)', fontsize=10)

# Legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor=colors['STABLE_RIGID'], label='STABLE_RIGID'),
    Patch(facecolor=colors['STABLE_MIDDLE'], label='STABLE_MIDDLE'),
    Patch(facecolor=colors['STABLE_FLEXIBLE'], label='STABLE_FLEXIBLE'),
]
ax.legend(handles=legend_elements, loc='upper left', fontsize=9)

ax.set_xlim(-0.5, max(valid['bootstrap_CI_width']) + 1)
ax.set_ylim(-0.5, max(valid['mean_rank_shift']) + 1.5)

plt.tight_layout()
plt.savefig(OUT / 'bootstrap_vs_crossatlas_stability.png', dpi=200, bbox_inches='tight')
plt.close()
print(f"\nSaved bootstrap_vs_crossatlas_stability.png")

# ── STEP 7: Generate convergent types summary ─────────────────────────
lines = []
lines.append("# Cross-Reference: Bootstrap Stability vs Cross-Atlas Consistency")
lines.append("")
lines.append(f"**Date:** 2026-04-06")
lines.append(f"**Purpose:** Check whether types that are bootstrap-stable (narrow CI)")
lines.append(f"are also consistent when measured in independent atlases.")
lines.append("")
lines.append("---")
lines.append("")

# Correlation result
lines.append("## Overall Correlation")
lines.append("")
lines.append(f"| Metric | Value |")
lines.append(f"|--------|-------|")
lines.append(f"| Types in ≥2 replications | {len(valid)} |")
lines.append(f"| Spearman ρ (CI width vs mean rank shift) | **{rho:.3f}** |")
lines.append(f"| Spearman p-value | {pval:.4f} |")
lines.append(f"| Pearson r | {r_pearson:.3f} |")
lines.append(f"| Pearson p-value | {p_pearson:.4f} |")
lines.append("")

if pval < 0.05 and rho > 0:
    lines.append("**Interpretation:** Significant positive correlation — bootstrap-unstable types")
    lines.append("are also cross-atlas-unstable. Within-atlas variability predicts cross-atlas variability.")
elif pval < 0.05 and rho < 0:
    lines.append("**Interpretation:** Significant negative correlation — bootstrap-stable types")
    lines.append("are paradoxically MORE variable across atlases. This would suggest that")
    lines.append("within-resampling stability and cross-atlas stability measure different things.")
elif abs(rho) < 0.2:
    lines.append("**Interpretation:** No meaningful correlation. Bootstrap stability (within-atlas")
    lines.append("resampling) and cross-atlas consistency measure largely independent properties.")
    lines.append("A type can be locked to a rank under resampling but shift when measured in")
    lines.append("a different atlas, or vice versa.")
else:
    lines.append(f"**Interpretation:** Weak {'positive' if rho > 0 else 'negative'} trend but not")
    lines.append(f"statistically significant (p = {pval:.3f}). The two stability measures are")
    lines.append(f"at best loosely coupled.")

lines.append("")
lines.append("---")
lines.append("")

# Category (a): stable both ways
lines.append("## (a) Stable in Bootstrap AND Consistent Across Atlases")
lines.append("")
lines.append(f"Criteria: bootstrap CI width ≤ 3, mean rank shift ≤ Q25 ({sd_q25:.1f}), "
             f"present in ≥2 replications.")
lines.append("")
if len(stable_boot_consistent) > 0:
    lines.append("| Cell Type | Primary Rank | Bootstrap CI | Category | Mean Shift | N Replications |")
    lines.append("|-----------|:------------|:------------|:---------|:----------|:--------------|")
    for _, r in stable_boot_consistent.sort_values('primary_rank').iterrows():
        lines.append(f"| {r['cell_type']} | {int(r['primary_rank'])} | "
                     f"[{r['bootstrap_CI_lower']:.0f}, {r['bootstrap_CI_upper']:.0f}] | "
                     f"{r['bootstrap_category']} | {r['mean_rank_shift']:.1f} | "
                     f"{int(r['n_replications_present'])} |")
else:
    lines.append("**No types meet both criteria.**")
lines.append("")

# Category (b): stable bootstrap, inconsistent across atlases
lines.append("## (b) Stable in Bootstrap BUT Inconsistent Across Atlases")
lines.append("")
lines.append(f"Criteria: bootstrap CI width ≤ 3, mean rank shift ≥ Q75 ({sd_q75:.1f}), "
             f"present in ≥2 replications.")
lines.append("")
if len(stable_boot_inconsistent) > 0:
    lines.append("| Cell Type | Primary Rank | Bootstrap CI | Category | Mean Shift | N Replications |")
    lines.append("|-----------|:------------|:------------|:---------|:----------|:--------------|")
    for _, r in stable_boot_inconsistent.sort_values('primary_rank').iterrows():
        lines.append(f"| {r['cell_type']} | {int(r['primary_rank'])} | "
                     f"[{r['bootstrap_CI_lower']:.0f}, {r['bootstrap_CI_upper']:.0f}] | "
                     f"{r['bootstrap_category']} | {r['mean_rank_shift']:.1f} | "
                     f"{int(r['n_replications_present'])} |")
else:
    lines.append("**No types meet both criteria.**")
lines.append("")

# Category (c): Pan-Census zero-shift types
lines.append("## (c) Pan-Census Zero-Shift Types")
lines.append("")
lines.append("These types had exactly zero rank shift in the Pan-Census (22-type) replication.")
lines.append("")
lines.append("| Cell Type | Primary Rank | Bootstrap CI | Width | Category | Census Rank Shift |")
lines.append("|-----------|:------------|:------------|:------|:---------|:-----------------|")
for ct in zero_shift_census:
    bci = boot_lookup.get(ct, {})
    mrow = master[master['cell_type'] == ct]
    if not mrow.empty:
        r = mrow.iloc[0]
        census_shift = r.get('Pan_Census_rank_shift', 'N/A')
        lines.append(f"| {ct} | {int(r['primary_rank'])} | "
                     f"[{bci.get('ci_lower','?')}, {bci.get('ci_upper','?')}] | "
                     f"{bci.get('ci_width','?')} | {bci.get('category','?')} | "
                     f"{int(census_shift) if not pd.isna(census_shift) else 'N/A'} |")
lines.append("")

# Category (d): CD8+ T cell
lines.append("## (d) CD8+ T Cell — Detailed Cross-Atlas Profile")
lines.append("")
cd8_row = master[master['cell_type'] == cd8_name].iloc[0]
lines.append(f"- **Primary rank:** {int(cd8_row['primary_rank'])} of 35 (most rigid)")
lines.append(f"- **Bootstrap:** median rank = {cd8_boot.get('median_rank')}, "
             f"CI = [{cd8_boot.get('ci_lower')}, {cd8_boot.get('ci_upper')}], "
             f"width = {cd8_boot.get('ci_width')}")
lines.append(f"- **Category:** {cd8_boot.get('category')}")
lines.append("")
lines.append("| Replication | Rank | N Types | Primary Subset Rank | Rank Shift |")
lines.append("|-------------|:-----|:--------|:-------------------|:----------|")
for rname in repl_names:
    rdata = replications[rname].get(cd8_name)
    if rdata:
        lines.append(f"| {rname} | {rdata['repl_rank']} | {rdata['n_types']} | "
                     f"{rdata['primary_subset_rank']} | {rdata['rank_shift']} |")
    else:
        lines.append(f"| {rname} | — | — | — | — |")

lines.append("")
cd8_shifts = [replications[rn][cd8_name]['rank_shift']
              for rn in repl_names if cd8_name in replications[rn]]
lines.append(f"**Mean rank shift:** {np.mean(cd8_shifts):.1f} "
             f"(SD = {np.std(cd8_shifts, ddof=1):.1f}, n = {len(cd8_shifts)})")
lines.append("")

if cd8_boot.get('ci_width') == 0:
    lines.append("**Key finding:** CD8+ T cell has bootstrap CI width = 0 (perfectly locked "
                 "at rank 35 under resampling) but shows substantial cross-atlas variability. "
                 "It is always among the most rigid types in the primary analysis bootstrap, "
                 "but its relative rank shifts considerably when measured in different atlases.")
lines.append("")

lines.append("---")
lines.append("")

# Full ranking table for types in ≥2 replications
lines.append("## Full Per-Type Summary (≥2 replications)")
lines.append("")
lines.append("| Cell Type | Primary | Boot CI | Boot Cat | Mean Shift | SD Shift | N Repl |")
lines.append("|-----------|:--------|:--------|:---------|:----------|:---------|:-------|")
for _, r in multi_repl.sort_values('mean_rank_shift').iterrows():
    sd_str = f"{r['sd_rank_shift']:.1f}" if not pd.isna(r['sd_rank_shift']) else "—"
    lines.append(f"| {r['cell_type']} | {int(r['primary_rank'])} | "
                 f"{r['bootstrap_CI_width']:.0f} | {r['bootstrap_category']} | "
                 f"{r['mean_rank_shift']:.1f} | {sd_str} | {int(r['n_replications_present'])} |")
lines.append("")

lines.append("---")
lines.append(f"*Generated: 2026-04-06*")

with open(OUT / 'convergent_types_summary.md', 'w') as f:
    f.write('\n'.join(lines))
print(f"Saved convergent_types_summary.md")

# ── STEP 8: Print final report ────────────────────────────────────────
print(f"\n{'='*70}")
print("CROSS-REFERENCE REPORT")
print(f"{'='*70}")
print(f"\nCorrelation (bootstrap CI width vs mean cross-atlas rank shift):")
print(f"  Spearman ρ = {rho:.3f}, p = {pval:.4f}")
print(f"  Pearson  r = {r_pearson:.3f}, p = {p_pearson:.4f}")
print(f"  n = {len(valid)} types")

print(f"\nTypes reliably RIGID everywhere:")
rigid_everywhere = multi_repl[
    (multi_repl['bootstrap_category'] == 'STABLE_RIGID') &
    (multi_repl['mean_rank_shift'] <= sd_q25)
]
if len(rigid_everywhere) > 0:
    for _, r in rigid_everywhere.iterrows():
        print(f"  {r['cell_type']:40s}  boot_CI={r['bootstrap_CI_width']:.0f}  "
              f"mean_shift={r['mean_rank_shift']:.1f}")
else:
    print("  (none)")

print(f"\nTypes reliably FLEXIBLE everywhere:")
flex_everywhere = multi_repl[
    (multi_repl['bootstrap_category'] == 'STABLE_FLEXIBLE') &
    (multi_repl['mean_rank_shift'] <= sd_q25)
]
if len(flex_everywhere) > 0:
    for _, r in flex_everywhere.iterrows():
        print(f"  {r['cell_type']:40s}  boot_CI={r['bootstrap_CI_width']:.0f}  "
              f"mean_shift={r['mean_rank_shift']:.1f}")
else:
    print("  (none)")

print(f"\nTypes stable in bootstrap but volatile across atlases:")
boot_stable_atlas_volatile = multi_repl[
    (multi_repl['bootstrap_CI_width'] <= 2) &
    (multi_repl['mean_rank_shift'] >= sd_q75)
]
if len(boot_stable_atlas_volatile) > 0:
    for _, r in boot_stable_atlas_volatile.iterrows():
        print(f"  {r['cell_type']:40s}  boot_CI={r['bootstrap_CI_width']:.0f}  "
              f"mean_shift={r['mean_rank_shift']:.1f}")
else:
    print("  (none)")

print("\nDone.")
