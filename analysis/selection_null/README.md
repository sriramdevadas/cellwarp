# Selection-aware null — C-selection circularity control (S1 Text §10)

The pre-registered control behind **S1 Text §10**, which tests the sharpest constructible
circularity objection to **Results §5**: that selecting genes on the per-gene conservation
score *C* manufactures the conserved-gene geometry (conserved-quartile obs/null = 0.384)
rather than measuring it.

The pipeline is not modified. `selection_null.py` imports the published code
(`gate_lib.per_gene_corr` for *C*; `cellwarp.procrustes.pca_reduce_centroids ->
procrustes_align -> permutation_test` for obs/null) and only re-selects its inputs.

## Pre-registered criteria

Fixed in `selection_null.py` before any draw was run, and written verbatim into each summary
JSON under `preregistered_criteria_SURFACED_not_declared`:

```
PASS_if : "real conserved obs/null < sigma-null 1st percentile (equivalently z <= -3)"
FAIL_if : "real falls within / near the sigma-null distribution"
```

The script surfaces the arithmetic; it declares no verdict. Design date recorded in the
script docstring and in each summary's `design` field: 2026-06-28. There is no separate
frozen analysis plan for this control, unlike the five `docs/preregistration_*.md` documents;
the conditions live in the code that computes them.

## Design

Per draw, with the selection procedure held fixed and only the cross-species correspondence
destroyed:

1. Draw sigma over the 35 cell types — a **derangement** (primary; no type maps to itself) or
   a **full permutation** (`--mode labelshuffle`, cross-check).
2. Permute the mouse centroid rows by sigma.
3. Recompute *C_sigma* and re-select the top quartile (3,985 genes; validity is
   permutation-invariant, so every draw selects exactly 3,985).
4. Re-run the unmodified obs/null pipeline under the sigma pairing (joint PCA, Procrustes,
   2,000-permutation inner label-shuffle null, median denominator).

N = 1,000 draws per mode. Sigma seed 12345; the inner pipeline keeps its published seed 42.

## Inputs

- `output/phase2/scaled_35types/centroids_human_35.csv` and `centroids_mouse_35.csv` — the
  35-type × 16,959-ortholog matched centroids, loaded via `gate_lib.load_centroids()`.
- `analysis/conserved_contribution/gate_lib.py` — `per_gene_corr`, `load_centroids`.
- `src/cellwarp/procrustes.py` — the unmodified obs/null pipeline.

No `.h5ad` and no network access: the whole control runs from the deposited centroids.

## Outputs (`outputs/`)

| file | contents |
|---|---|
| `selection_null_summary_{derangement,labelshuffle}.json` | substrate counts, real values, sigma-null moments and percentiles, z, and the pre-registered arithmetic |
| `sigma_null_draws_{...}.csv` | one row per draw: `draw_id,ratio,n_conserved,n_valid,q75,distance,null_median,n_fixed_points` |
| `sigma_perms_{...}.npy` | the (1000, 35) int64 sigma array, for exact replay |

## Results

Real conserved obs/null **0.3843**; full-space **0.5219**.

| mode | sigma-null mean ± sd | 1st pct | z | draws ≤ real |
|---|---|---|---|---|
| derangement (primary) | 0.9907 ± 0.0206 | 0.9273 | −29.49 | 0 / 1000 |
| label-shuffle (cross-check) | 0.9833 ± 0.0236 | 0.9165 | −25.34 | 0 / 1000 |

Both pre-registered PASS conditions hold in both modes. Mechanism: a derangement collapses
the selection threshold (mean Q75 of *C* falls from 0.5919 to 0.078 under derangement, 0.090
under label-shuffle), so with the correspondence destroyed the re-selection has no genuine
signal to grab and obs/null returns to ~1.0.

## Scope

Addresses whether selecting on *C* manufactures the configuration result. It does **not**
address the ontology type-matching upstream of the full-space 0.522 — that is what the
marker-similarity-stratified null (`analysis/sensitivity_analyses/markernull.py`, S5 Fig and
S10 Table) speaks to. It does not bear on the Layer-2 covariance analysis.

## Reproduce

```bash
.venv/bin/python analysis/selection_null/repro_baseline.py                       # ~12 s
.venv/bin/python analysis/selection_null/selection_null.py --mode derangement    # ~52 s
.venv/bin/python analysis/selection_null/selection_null.py --mode labelshuffle   # ~53 s
```

Defaults are the published settings (`--n 1000 --nperm 2000 --workers 8 --seed 12345`), and
`--out` defaults to `outputs/` beside this file, so a re-run overwrites the deposited results.
Pass `--out` elsewhere to keep them.

Re-run against this tree, the control is deterministic and reproduces the deposited artifacts
exactly: both `sigma_null_draws_*.csv` and both `sigma_perms_*.npy` come back byte-identical,
and both summary JSONs match on every field except `runtime_sec`.
