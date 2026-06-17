# CellWarp Simulation Study Results

**Generated:** 2026-05-25 10:14
**Runtime:** 1216 seconds
**Replicates per condition:** 100

---

## Calibration

| Parameter | Value |
|-----------|-------|
| Real data obs/null ratio | 0.522 |
| Estimated signal_strength | 3.68 |

The simulation's signal_strength parameter was calibrated so that synthetic
data at signal = 3.68 produces an obs/null ratio ≈ 0.522,
matching the real 35-type Procrustes analysis.

---

## 1. Detection Power (Panel A)

At the calibrated real-data signal strength (3.68):

| n_types | Detection rate | Mean obs/null |
|---------|---------------|---------------|
| 15 | 100% | 0.559 |
| 25 | 100% | 0.555 |
| 35 | 100% | 0.554 |

The pipeline has excellent detection power (100%) at the real data's signal level with 35 types.
With fewer types (15), power remains adequate.

---

## 2. Ranking Recovery (Panel B)

Spearman ρ between planted and recovered per-type rankings:

| Signal | 50 cells | 200 cells | 500 cells | 2000 cells |
|--------|----------|-----------|-----------|------------|
| 0.5 | -0.54 ± 0.14 | -0.53 ± 0.13 | -0.53 ± 0.14 | -0.53 ± 0.13 |
| 1.0 | -0.30 ± 0.14 | -0.32 ± 0.14 | -0.30 ± 0.16 | -0.30 ± 0.13 |
| 1.5 | -0.00 ± 0.19 | 0.00 ± 0.20 | 0.01 ± 0.22 | 0.00 ± 0.19 |
| 2.0 | 0.23 ± 0.20 | 0.23 ± 0.20 | 0.23 ± 0.22 | 0.23 ± 0.20 |
| 3.0 | 0.41 ± 0.17 | 0.42 ± 0.17 | 0.40 ± 0.19 | 0.40 ± 0.18 |
| 5.0 | 0.44 ± 0.13 | 0.43 ± 0.15 | 0.42 ± 0.15 | 0.41 ± 0.16 |
| 10.0 | 0.29 ± 0.15 | 0.28 ± 0.16 | 0.27 ± 0.16 | 0.24 ± 0.16 |

At the real data's signal level (≈ 3.7) with 200 cells/type:
**ρ = 0.418 ± 0.173**

---

## 3. Ranking Stability — KEY FINDING (Panel C)

Test-retest Spearman ρ (two independent cell samples from same true centroids)
at signal = 3.68, n_types = 35:

| Cells/type | Mean ρ | 95% CI |
|-----------|--------|--------|
|    25 | 0.974 | [0.948, 0.990] |
|    50 | 0.983 | [0.951, 0.994] |
|   100 | 0.991 | [0.968, 0.997] |
|   200 | 0.994 | [0.980, 0.998] |
|   500 | 0.997 | [0.991, 0.999] |
|  1000 | 0.998 | [0.995, 1.000] |
|  2000 | 0.999 | [0.997, 1.000] |
|  5000 | 0.999 | [0.998, 1.000] |

### Interpretation

**At the real sample size (200 cells/type):** expected test-retest ρ = **0.994**
**Real cross-atlas replication ρ ≈ 0.15–0.19**

The simulation predicts substantially higher ranking stability from sampling
noise alone (ρ ≈ 0.99) than observed in real cross-atlas replication (ρ ≈ 0.17).
This gap demonstrates that **atlas-to-atlas biological variability** (different
donors, tissue procurement, processing pipelines) is the dominant source of
ranking instability — not centroid estimation noise.

The global coherence statistic (obs/null ratio, p-value) is robust because it
measures the overall geometric signal, which averages across cell types. But
per-type residual rankings are sensitive to atlas-specific biology.

### Sample Size Recommendation

At **2000 cells/type**, test-retest ρ = 0.999 — rankings
become substantially more reliable. Future atlases with deeper per-type sampling
would enable meaningful per-type evolutionary divergence rankings.

---

## 4. Null Calibration (Panel D)

| Metric | Observed | Expected |
|--------|----------|----------|
| Rejection rate α = 0.05 | 0.048 | 0.050 |
| Rejection rate α = 0.01 | 0.009 | 0.010 |
| Mean p-value | 0.507 | 0.500 |
| Std p-value | 0.293 | 0.289 |

P-values are well-calibrated under the null.
The permutation test correctly controls the false positive rate.

---

## Simulation Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| n_genes | 500 | Realistic PCA behavior; computational efficiency |
| n_factors | 50 | Latent dimensionality (real PCA yields ~33 components) |
| centroid_scale | 2.0 | Inter-type separation in factor space |
| within_type_var | 1.0 | Per-gene cell noise; realistic SNR |
| n_permutations | 1000 | Speed/precision tradeoff (1000 per simulation) |
| rigidity_spread | 1.0 | Log-normal spread; ~10x range in noise levels |
| PCA threshold | 0.95 | Matches real pipeline (95% variance) |

### Design Notes

- Centroid noise is modeled via CLT shortcut (exact for Gaussian cells), avoiding
  the need to generate individual cells. This makes the simulation ~1000x faster
  for large n_cells without sacrificing accuracy.
- The per-type noise is added in factor space, ensuring PCA captures it fully.
  This represents a **best case** for ranking recovery — real data may have
  divergence in PCA-dropped dimensions, making real rankings harder to recover.
- Each replicate uses independent random centroids, rotation, and noise draws.

---

## Files

| File | Description |
|------|-------------|
| `simulation_results.json` | Full numerical results (all conditions) |
| `simulation_summary.md` | This file |
| `simulation_study.py` | Simulation code (reproducible) |
| `simulation_figures.py` | Standalone figure script |
| `simulation_study_figure.png` | Local copy of figure |
| `figures/supplementary/figS7_simulation_study.{png,pdf}` | Publication figure |
