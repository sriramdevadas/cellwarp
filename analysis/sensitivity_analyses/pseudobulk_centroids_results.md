# Pseudo-bulk centroid-definition sensitivity (Layer-1, human-mouse, 35 types)

Seed 42, 10,000 permutations, target_sum 1e+04, joint PCA at 95% variance.
Everything except the centroid definition is held fixed: same 16,959 orthologs,
same 35 cell types, same cells, same Procrustes/permutation engine.

| centroid definition | PCs | Procrustes d | obs/null (global) | p | obs/null (lineage-stratified) | p (lineage) | rank rho vs deposited |
|---|---|---|---|---|---|---|---|
| deposited mean-of-log1p (anchor) | 33 | 61.153 | 0.5222 | 1.00e-04 | 0.6683 | 1.00e-04 | 1.000 |
| mean-of-log1p rebuilt from raw (anchor) | 33 | 61.153 | 0.5222 | 1.00e-04 | 0.6683 | 1.00e-04 | 1.000 |
| A. equal-cell-weight  log1p(mean CP10K) | 36 | 71.771 | 0.4866 | 1.00e-04 | 0.6498 | 1.00e-04 | 0.827 |
| B. aggregate pseudo-bulk  log1p(CP10K(sum raw)) | 45 | 98.538 | 0.6200 | 1.00e-04 | 0.7991 | 1.00e-04 | 0.496 |

Deposited baseline: global obs/null = 0.5222, lineage-stratified obs/null = 0.6683, 33 PCs, Procrustes d = 61.153.

## Verdict

**equal_cell_weight: PASS** (obs/null = 0.4866 < 0.8, p = 1.00e-04 < 0.01; delta vs 0.522 baseline = -0.0356; lineage-stratified obs/null = 0.6498, delta = -0.0184; retained PCs 36 vs 33 (CHANGED))

**aggregate_pseudobulk: PASS** (obs/null = 0.6200 < 0.8, p = 1.00e-04 < 0.01; delta vs 0.522 baseline = +0.0978; lineage-stratified obs/null = 0.7991, delta = +0.1308; retained PCs 45 vs 33 (CHANGED))

## Fidelity gates

Both anchors were pushed through the identical code path before the pseudo-bulk flavors:

- `deposited` (the shipped centroid CSVs) returns obs/null = 0.5222043226858 and lineage-stratified 0.6682720235159, bit-identical to the published 0.5222043226858 / 0.6682720235159.
- `recon_mean_log1p` rebuilds the primary definition from the same raw counts this script reads, and lands at 0.522204353 (delta +3.01e-08); centroid values agree with the deposited CSVs to max|diff| = 1.00e-05 (float32 round-off). The raw-count provenance and the streaming accumulator are therefore the same object the deposit was built from.
- Cells: 62,708 human / 59,745 mouse, 0 zero-library cells.

## Caveats

- **p-values are at the permutation floor.** 0 of 10,000 permutations reached the observed distance for every flavor and for both nulls, so p = 1/(10,000+1) = 1.0e-04 is a bound, not an estimate. The obs/null ratio, not p, is the quantity being compared here.
- **The retained-PC count is not stable across centroid definitions.** 95% variance is a relative criterion, and the pseudo-bulk definitions redistribute variance across the joint spectrum: 33 PCs (primary) -> 36 (equal_cell_weight), 45 (aggregate_pseudobulk). Procrustes distances are consequently not comparable across rows in absolute terms; each row's obs/null is computed against its own permutation null in its own space, which is what makes the ratios comparable.
- **Per-type rigidity ordering is more sensitive than configuration-level coherence.** Spearman rho of per-type residual magnitude against the deposited ranking is 0.827 (p = 9.3e-10) for A and 0.496 (p = 2.4e-03) for B. Library-size weighting (B) reshuffles which types look rigid more than it weakens the global result.
