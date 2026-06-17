# Data Provenance — Donor-Stability Resampling

Human-readable mirror of the machine-readable records
`donor_stability/phase0_summary.json` (cell/donor structure) and the validity
block of `donor_stability/donor_stability_results.json` (deposit correspondence).
Provenance only.

## Source

The donor-stability analysis re-acquires cell-level data from
**CZ CELLxGENE Census, version 2025-11-08**, for the same **Tabula Sapiens**
(human) and **Tabula Muris Senis** (mouse) collections used to build the
deposit's 35-type centroids. Cells are selected with
`is_primary_data == True and disease == 'normal'`, restricted to the 35 matched
cell types, normalized per cell (`normalize_total(target_sum=1e4)` + `log1p`),
and aggregated to mean centroids — matching the deposit pipeline exactly. The
per-(type, donor, protocol) aggregates and within-donor cell-split half-centroids
are regenerable from Census via `donor_stability/pull_aggregate.py` (the metadata
map is produced by `donor_stability/phase0_metadata.py`); these aggregate files
are large and are not redistributed (see `.gitignore`).

## Cell and donor structure (per `phase0_summary.json`)

| Atlas | Cells | Donors | Datasets | Assays |
|---|---|---|---|---|
| Human (Tabula Sapiens) | 804,639 | 24 | 1 | 10x 3' v3, 10x 5' v2, Smart-seq2, Smart-seq3 |
| Mouse (Tabula Muris Senis) | 187,553 | 44 | 2 | 10x 3' v2, Smart-seq2 |

All 35 cell types are donor-powered in both species (human: 32 types with ≥ 4
donors, 29 with ≥ 6; mouse: all 35 with ≥ 6), so a per-species donor-split is
feasible in both atlases.

## Deposit-centroid correspondence (validity, per `donor_stability_results.json`)

Re-pulling the cells from Census and re-deriving centroids reproduces the
deposited 35-type centroids and the per-gene conservation score C:

- median per-type human centroid correlation, fresh-pull vs deposit = **1.000**
- per-gene C correlation, fresh-pull vs deposit = **0.989** (Spearman and Pearson)
- vectorized per-gene Pearson vs the frozen gate library, max\|Δ\| = **1.0e-15**
- fresh-pull 35-type Procrustes obs/null = **0.521** (deposit anchor 0.522)

This faithful re-acquisition is what underwrites the donor-split, cell-sampling
ceiling, cap, and cross-protocol comparisons that follow; the donor resampling is
the only operation introduced on top of the deposit's frozen quantities.
