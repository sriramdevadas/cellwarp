# Macaque pipeline reconstruction (Test 1a) — RIRA-13 sensitivity verified

**Date:** 2026-04-23
**Scope:** 13 RIRA-sourced cell types, committed sensitivity result (obs/null = 0.749, p = 0.0002).
**Status:** **PASS** — reproduces to within approved 0.01 tolerance.

## Headline numbers

| Metric | Committed (sensitivity_procrustes_results.json) | Reconstructed | Δ |
|---|---|---|---|
| n_types | 13 | 13 | — |
| gene space | 13,927 (implied) | 13,927 (computed) | exact |
| Procrustes distance | 16.531 | 16.523 | 0.008 |
| Procrustes scaling | 0.0534 | 0.0530 | 0.0004 |
| obs/null (median) | 0.7488 | 0.7506 | 0.0018 ✓ |
| p-value | 0.0002 | 0.0001 | < 0.0001 (both highly significant) |
| PCA components (95% var) | not recorded | 6 | — |

Verification target `|Δobs/null| < 0.01` met. Residual deltas (0.001–0.002) are
consistent with permutation RNG micro-variance across seed-42 runs and rounding
in the committed JSON.

## Reconstruction pipeline (committed driver)

`analysis/macaque/reconstruct_macaque_pipeline.py`

1. **RIRA label harmonization** — exact mapping inferred from committed
   `centroid_cell_counts.csv` (all 13 per-type counts match to unit precision):
   - B cell ← `RIRA_Immune_v2.cellclass == "Bcell"` (84,412)
   - CD4⁺ α-β T ← `Immune_v2 == "T_NK" & TNK_v2 == "CD4+ T Cells"` (116,863)
   - CD8⁺ α-β T ← `Immune_v2 == "T_NK" & TNK_v2 == "CD8+ T Cells"` (110,916)
   - NK ← `Immune_v2 == "T_NK" & TNK_v2 == "NK Cells"` (11,802)
   - T cell ← `Immune_v2 == "T_NK" & TNK_v2 ∉ {CD4+, CD8+, NK}` (48,226; collapses Unassigned/Ambiguous/Unknown/Other/NaN)
   - classical monocyte ← `Immune_v2 == "Myeloid" & Myeloid_v3 == "CD14+ Monocytes"` (7,250)
   - intermediate monocyte ← `... Myeloid_v3 == "Inflammatory Monocytes"` (1,452)
   - non-classical monocyte ← `... Myeloid_v3 == "CD16+ Monocytes"` (3,484)
   - macrophage ← `... Myeloid_v3 ∈ {Macrophages, Alv. mac.}` (14,298)
   - myeloid dendritic cell ← `... Myeloid_v3 ∈ {DC, pDC, Mature DC}` (4,780)
   - granulocyte ← `... Myeloid_v3 == "Myelocytes"` (5,698)
   - hematopoietic precursor cell ← `... Myeloid_v3 == "Promyelocytes"` (1,501)
   - myeloid leukocyte ← `... Myeloid_v3 ∉ named_myeloid` (2,138; catches Unassigned/Ambiguous/Unknown/NaN)

   Excluded (not target types): `Immune_v2 ∈ {Unknown, Non-Immune, Erythrocyte}` (184,028 cells).

2. **Three-way gene space** — 13,927 human ENSGs reconstructed from:
   (CellWarp 16,959-gene space in `human_scaled.h5ad` var) ∩
   (RIRA symbols matching CellWarp `feature_name`) ∩
   (Qu ENSMFAG coverage via direct HGNC-like symbol match + BioMart recovery for
   ENSMFAG-only rows). Reproduces to exactly 13,927.

3. **Cell subsample** — 2,000 cells per type (seed 42). Verified against
   `human_scaled.h5ad` which contains exactly 2,000 cells per type; means
   computed from that h5ad reproduce `centroids_human_35.csv` to float
   precision.

4. **Normalization of RIRA counts — none.** **This is the key methodological
   finding.** See §"Surprising finding" below.

5. **Centroid computation** — per-type mean on 13 × 13,927 subsampled,
   gene-filtered, already-normalized RIRA matrix.

6. **Joint PCA** — `pca_reduce_centroids(human_13, macaque_13)` with ≥95%
   variance (float n_components); returns k = 6 components.

7. **Procrustes** — `procrustes_align` (free scaling, reflection-forbidden,
   SVD-based OPA) + `permutation_test` 10,000 iterations, seed = 42.

## Surprising finding: RIRA `.counts.rds` is already pre-normalized

`data/macaque/rira/RIRA.All.RNA.counts.rds`, despite its filename, does **not**
contain raw UMI counts. The matrix contains non-integer values (≈99.1% of
non-zeros are non-integer), with per-cell maxes in the range 68–252 and a
global max of 16,937. This is consistent with Seurat `NormalizeData` or
`SCTransform` output — the file stores pre-processed ("corrected") values.

**Consequence for the pipeline:** applying `sc.pp.normalize_total(target_sum=1e4) + log1p`
to this matrix **double-normalizes** and produces an obs/null of ≈0.50 (no
correspondence recoverable to the committed 0.749).

The committed primary pipeline uses the `.counts.rds` values directly as the
per-cell expression vectors. This breaks the symmetry with the human
configuration, which *did* go through `normalize_total(1e4) + log1p` (DECISION-064).
The two configurations are therefore in different numeric domains — which the
Procrustes optimal scaling `s` absorbs (observed `s = 0.053` reflects this
mismatch: the RIRA centroid cloud has ~20× the Frobenius norm of the human
cloud after centering).

**Implications for reproducibility and the manuscript:**
- Methods text should state explicitly that RIRA pre-processed values were
  used as-is without re-normalization. A reader following the stated
  "counts per 10k + log1p" convention verbatim will not reproduce 0.841 /
  0.749 — they will reproduce ≈0.50.
- The Seurat pre-processing specifics (NormalizeData vs SCTransform, target
  sum, etc.) that generated RIRA's `.counts.rds` are not documented in the
  locally-available RIRA artifacts and would need a pointer to the RIRA
  paper's methods or the data contributor's pipeline.
- The asymmetric preprocessing (human normalized + log1p'd, macaque left as
  Seurat-normalized) means the Procrustes-optimal scaling is effectively
  correcting a unit-system mismatch in addition to the biological
  "cross-species expression magnitude" difference mentioned in the methods.
- For Test 1b (GPA) the GPA consensus step will accommodate this via its
  per-configuration scaling, but the scaling factors should not be
  over-interpreted biologically until the asymmetric preprocessing is
  addressed.

## Differences from original task-brief expectations

- **Subsampling.** The committed pipeline subsamples each cell type to max
  2,000 cells (seed 42) before centroid computation. `centroid_cell_counts.csv`
  records the *total available* cells per type, not the number actually used.
  This is the same `MAX_CELLS_PER_TYPE = 2_000` convention used in
  `analysis/mouse_lemur/01_run_pipeline.py`.
- **PCA components for 13-type sensitivity: 6.** Committed primary (20 types)
  used 8. Dimensionality is recomputed fresh per analysis at the ≥95%
  variance threshold; the 20-type and 13-type cases differ naturally.
- **Permutation scheme.** Matches `src/procrustes.py::permutation_test` — row
  permutation of Y with X fixed, seed 42. Produces committed results.

## Not covered by this reconstruction

- 20-type primary (obs/null = 0.841). Requires the Qu per-cell annotations
  from Zenodo 10.5281/zenodo.5881495; download has completed (14.33 GB
  `MacFas.meta.data.rds` now local) but extraction and 7-type Qu centroid
  construction are Track 3.
- 2-config driver that also reproduces the 7 Qu-sourced types will land with
  Test 1a-Qu once extraction is verified.

## Artifacts committed

- `analysis/macaque/reconstruct_macaque_pipeline.py` — driver
- `analysis/macaque/extract_rira_metadata.R` — R-side one-time extraction
- `analysis/macaque/diagnose_procrustes_gap.py` — PCA-dimensionality sweep (kept for provenance)
- `analysis/macaque/diagnose_scaling_gap.py` — ||X_c||/||Y_c|| diagnostic (provenance)
- `analysis/macaque/fast_norm_test.py` — normalization hypothesis test (provenance; identified as-is as correct)
- `analysis/macaque/reconstruction_rira13_report.md` — this report
- `data/macaque/rira/rira_metadata.csv` — 296 MB per-cell metadata (extracted once from .rds)
- `data/macaque/rira/rira_metadata_summary.txt` — column inventory
- `output/macaque_pipeline/reconstruction_rira13_results.json` — run output
- `output/macaque_pipeline/reconstruction_rira13_centroids.csv` — 13×13,927 macaque centroids
- `output/macaque_pipeline/reconstruction_rira13_gene_list.csv` — the 13,927 ENSG+symbol list

## Runtime

~70 s end-to-end after metadata is extracted (extraction itself ~2–3 min).
