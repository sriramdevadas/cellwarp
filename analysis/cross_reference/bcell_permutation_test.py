"""
B cell cross-atlas ranking consistency permutation test.

Biology: Tests whether B cell's remarkably stable ranking across 7 independent
atlases/replications is significantly better than expected by chance.

Math: In each replication, types are ranked by residual magnitude. The primary
ranking is re-ranked within the subset of types present in that replication.
rank_shift = |primary_subset_rank - replication_rank|. We permute the
replication rank assignments (keeping which types are present fixed) and ask:
how often does ANY type with 7+ replications achieve a mean rank shift as low
as B cell's observed 3.14?
"""

import numpy as np
import pandas as pd
from scipy.stats import rankdata

np.random.seed(42)

# Load data
df = pd.read_csv("analysis/cross_reference/master_ranking_table.csv")

# Replication columns
rep_cols = [
    "Sun2023_rank", "PanSci_rank", "CellHint_rank",
    "CellHint_harmonized_rank", "Pan_Census_rank",
    "Macaque_rank", "Mouse_lemur_rank",
]

primary_rank = df["primary_rank"].values  # global rank (1=most flexible)
n_types = len(df)

# For each replication, build:
#   present_idx: which types are present
#   repl_ranks: their replication ranks
#   primary_subset_ranks: re-ranked primary among the subset
rep_data = []
for col in rep_cols:
    present_mask = df[col].notna().values
    present_idx = np.where(present_mask)[0]
    if len(present_idx) == 0:
        rep_data.append(None)
        continue

    repl_ranks = df[col].values[present_idx].astype(float)

    # Re-rank primary within subset: rank 1 = smallest primary_rank = most flexible
    primary_in_subset = primary_rank[present_idx]
    primary_subset_ranks = rankdata(primary_in_subset, method='ordinal').astype(float)

    rep_data.append({
        'present_idx': present_idx,
        'repl_ranks': repl_ranks,
        'primary_subset_ranks': primary_subset_ranks,
        'n_types_in_repl': len(present_idx),
    })

# Verify B cell observed stats
bcell_idx = df.index[df["cell_type"] == "B cell"][0]
bcell_primary = primary_rank[bcell_idx]
bcell_shifts = []
for rd in rep_data:
    if rd is None:
        continue
    loc = np.where(rd['present_idx'] == bcell_idx)[0]
    if len(loc) > 0:
        i = loc[0]
        shift = abs(rd['primary_subset_ranks'][i] - rd['repl_ranks'][i])
        bcell_shifts.append(shift)

bcell_observed_mean_shift = np.mean(bcell_shifts)
bcell_n_reps = len(bcell_shifts)

print(f"B cell: primary_rank={bcell_primary}, n_replications={bcell_n_reps}")
print(f"B cell observed per-replication shifts: {bcell_shifts}")
print(f"B cell observed mean shift: {bcell_observed_mean_shift:.4f}")
print(f"  (matches CSV mean_rank_shift={df.iloc[bcell_idx]['mean_rank_shift']:.4f})")
print()

# Count replications per type
n_reps_per_type = np.zeros(n_types, dtype=int)
for rd in rep_data:
    if rd is None:
        continue
    n_reps_per_type[rd['present_idx']] += 1

types_with_7plus = np.where(n_reps_per_type >= 7)[0]
print(f"Types with 7+ replications: {len(types_with_7plus)}")
for idx in types_with_7plus:
    print(f"  {df.iloc[idx]['cell_type']}: {n_reps_per_type[idx]} reps, "
          f"observed mean shift = {df.iloc[idx]['mean_rank_shift']:.2f}")
print()

# Permutation test
N_PERM = 10000
observed_threshold = bcell_observed_mean_shift

null_min_shifts = np.zeros(N_PERM)

# Also collect per-type null distributions for types with 7+ reps
type_null_shifts = {idx: [] for idx in types_with_7plus}

for p in range(N_PERM):
    shift_sums = np.zeros(n_types)
    shift_counts = np.zeros(n_types, dtype=int)

    for rd in rep_data:
        if rd is None:
            continue
        present_idx = rd['present_idx']
        psr = rd['primary_subset_ranks']

        # Shuffle replication ranks among the present types
        shuffled_repl_ranks = np.random.permutation(rd['repl_ranks'])

        # Compute |primary_subset_rank - shuffled_repl_rank|
        shifts = np.abs(psr - shuffled_repl_ranks)
        shift_sums[present_idx] += shifts
        shift_counts[present_idx] += 1

    # For types with 7+ reps, compute mean shift
    min_shift = np.inf
    for idx in types_with_7plus:
        if shift_counts[idx] >= 7:
            ms = shift_sums[idx] / shift_counts[idx]
            type_null_shifts[idx].append(ms)
            if ms < min_shift:
                min_shift = ms
    null_min_shifts[p] = min_shift

# P-value: fraction where ANY type with 7+ reps achieves mean shift ≤ observed
p_value = np.mean(null_min_shifts <= observed_threshold)

# Expected mean shift under null per type
print("Expected mean shift under null (per type with 7+ reps):")
for idx in types_with_7plus:
    arr = np.array(type_null_shifts[idx])
    obs = df.iloc[idx]['mean_rank_shift']
    print(f"  {df.iloc[idx]['cell_type']:45s}: null={arr.mean():.2f}±{arr.std():.2f}, "
          f"observed={obs:.2f}")
print()

overall_null_mean = np.mean([np.mean(type_null_shifts[idx]) for idx in types_with_7plus])
print(f"Overall expected mean shift under null: {overall_null_mean:.2f}")

# For a uniform permutation of n items, E[|i - sigma(i)|] ≈ n/3
# Compute per-replication expected
for rd in rep_data:
    if rd is not None:
        n = rd['n_types_in_repl']
        print(f"  Replication with {n} types: n/3 = {n/3:.1f}")
print()

# Results
print("=" * 60)
print("PERMUTATION TEST RESULTS")
print("=" * 60)
print(f"Observed B cell mean rank shift: {bcell_observed_mean_shift:.2f}")
print(f"  across {bcell_n_reps} replications")
print(f"Number of permutations: {N_PERM:,}")
print(f"P-value (any type with 7+ reps achieves ≤ {observed_threshold:.2f}): {p_value:.4f}")
print()
print(f"Null distribution of minimum mean shift (7+ reps):")
print(f"  Mean: {null_min_shifts.mean():.2f}")
print(f"  SD:   {null_min_shifts.std():.2f}")
print(f"  5th percentile:  {np.percentile(null_min_shifts, 5):.2f}")
print(f"  Median:          {np.median(null_min_shifts):.2f}")
print(f"  95th percentile: {np.percentile(null_min_shifts, 95):.2f}")
print()

# B cell type-specific p-value
bcell_null = np.array(type_null_shifts[bcell_idx])
bcell_type_p = np.mean(bcell_null <= bcell_observed_mean_shift)
print(f"B cell type-specific p-value: {bcell_type_p:.4f}")
print(f"  (fraction of null where B cell alone achieves ≤ {bcell_observed_mean_shift:.2f})")
print(f"B cell null: mean={bcell_null.mean():.2f}, sd={bcell_null.std():.2f}")
print(f"B cell z-score: {(bcell_observed_mean_shift - bcell_null.mean()) / bcell_null.std():.2f}")
print()

# All types with 7+ reps: type-specific p-values
print("Type-specific p-values (all types with 7+ reps):")
for idx in types_with_7plus:
    arr = np.array(type_null_shifts[idx])
    obs = df.iloc[idx]['mean_rank_shift']
    tp = np.mean(arr <= obs)
    z = (obs - arr.mean()) / arr.std()
    print(f"  {df.iloc[idx]['cell_type']:45s}: p={tp:.4f}, z={z:.2f}, "
          f"obs={obs:.2f}, null={arr.mean():.2f}")
print()

# Percentile of B cell in null min distribution
bcell_percentile = np.mean(null_min_shifts <= bcell_observed_mean_shift) * 100
print(f"B cell's mean shift ({bcell_observed_mean_shift:.2f}) falls at the "
      f"{bcell_percentile:.1f}th percentile of the null min distribution")
