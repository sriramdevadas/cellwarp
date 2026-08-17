> **Historical snapshot (as of 2026-04-06).** The reproducibility gaps recorded below reflect the repository at that time. The statsmodels item (M1) has since been resolved: `statsmodels==0.14.6` is now pinned in `requirements.txt`, and the Key Resources Table referenced as its source has been retired in the PLOS reformat. Other snapshot gaps have likewise been closed since: the Census version pin is now present in `08_scaled_procrustes.py` and the other two flagged scripts (C1); the code and data repository URLs are now live Zenodo DOIs (C2); and the README permutation count (m6) and the 0.05^4 multiple-comparison note (m10) have been corrected.

# Reproducibility Audit v2: Full Manuscript

**Audit mode:** Self-audit against an external-reproducer standard: every quantitative claim re-derived from the manuscript and code alone, without recourse to working notes or undeposited intermediates. Not an independent review – the paper is single-author, and the audit was carried out by the author with AI assistance (see Use of generative AI in Methods).
**Date:** 2026-04-06
**Scope:** All quantitative claims in manuscript_combined.txt, full codebase review

---

## EXECUTIVE SUMMARY

**Overall grade: B (Mostly Reproducible)**

The manuscript is unusually thorough for a computational biology paper. Software versions are pinned, random seeds are documented and enforced, the core algorithm is described at implementation level, and the README provides a clear execution order. However, several gaps prevent an "A" rating: (1) three scripts that query CZ CELLxGENE Census omit the version pin, including the primary analysis script; (2) BioMart ortholog queries hit a live (unversioned) API endpoint; (3) no automated pipeline runner exists; (4) two software dependencies listed in the manuscript are missing from or unpinned in requirements.txt; (5) the code repository URL is a placeholder. None of these are fatal – cached intermediate files and the archived ortholog CSV provide de facto reproducibility for reruns – but a from-scratch reproduction on a new machine would require manual intervention.

---

## 1. DATA SOURCES AND ACCESSIBILITY

### 1.1 Primary atlases

| Dataset | Accession in Manuscript | Accessible? | Version Pinned? |
|---|---|---|---|
| Tabula Sapiens | Census v2025-11-08; figshare 10.6084/m9.figshare.14267219 | Yes | Yes (Census version in Methods + README) |
| Tabula Muris Senis | Census v2025-11-08; figshare 10.6084/m9.figshare.12654728 | Yes | Yes |
| Sun et al. 2023 | OMIX002605 / GSA CRA007207 | Yes | N/A (static deposit) |
| PanSci | GEO GSE247719 | Yes | N/A (static deposit) |
| CellHint Human | Census v2025-11-08 | Yes | Yes |
| Tabula Microcebus | "CZ CELLxGENE Discover" | **Partial** – no DOI/GEO; only reference [28] | **Gap: no version pin or permanent identifier** |
| RIRA macaque | GEO GSE277821 | Yes | N/A |
| Qu et al. macaque | GEO GSE196791 | Yes | N/A |
| Cancer/COVID | Census v2025-11-08 | Yes | Yes |

### 1.2 Reference databases

| Resource | Version in Manuscript | Pinned in Code? |
|---|---|---|
| Ensembl BioMart | Release 115, GRCh38.p14; accessed 12 Mar 2026 | **No** – `src/data_loader.py` queries live `http://www.ensembl.org`, not an archive URL |
| CellMarker 2.0 | Downloaded 2026-03-16 | **No version number** – only download date |
| STRING | v12.0 | Not verified in code |
| MSigDB (niche gene sets) | **Not specified** | **Gap** |
| DoRothEA | **Not specified** | **Gap** |
| ChEMBL/DrugBank | **Not specified** | **Gap** |
| ENCODE (H3K27ac) | **Not specified** – no experiment IDs | **Gap** |
| DILIrank | FDA LTKB (URL given) | Static database |
| L1000 | LINCS; 978 genes | Static gene set |

**Mitigating factor:** The ortholog table is cached as `data/phase1/orthologs_human_mouse.csv` after the first BioMart query. The manuscript states "The full ortholog table is archived in the repository as a CSV file and constitutes the reproducibility anchor." This provides de facto reproducibility if the cached file is distributed with the code. The macaque ortholog table is similarly cached at `data/macaque/biomart_macaque_human_orthologs.csv`.

### 1.3 Data Availability Statement

The manuscript states:
- Code: `[REPOSITORY URL – to be inserted upon submission]` ← **placeholder**
- Data: `[DATA REPOSITORY – to be inserted upon submission]` ← **placeholder**
- DILI pre-registration archived at git commit 415c806 ← verifiable

**Finding:** Code and data repository URLs are not yet filled in. This is expected pre-submission but must be resolved.

---

## 2. METHODS SPECIFICATION

### 2.1 Core pipeline parameters

| Parameter | Manuscript | Code | Match? |
|---|---|---|---|
| Random seed | RANDOM_SEED=42 | `RANDOM_SEED = 42` in `src/procrustes.py` and all scripts | **Yes** |
| PCA threshold | 95% cumulative variance → 33 components | `PCA_VARIANCE_THRESHOLD = 0.95` in `src/procrustes.py` | **Yes** |
| Normalization | CP10K + log1p | `TARGET_SUM = 10000` in `src/qc.py` | **Yes** |
| Cell cap | 2,000 per type per species | `MAX_CELLS_PER_TYPE = 2000` in `src/data_loader.py` | **Yes** |
| Min cells | 500 per type per species | `MIN_CELLS_PER_TYPE = 500` in `src/data_loader.py` | **Yes** |
| Primary permutations | 1,000,000 | `N_PERMUTATIONS = 1_000_000` in `scripts/permutation_1M.py` | **Yes** |
| Replication permutations | 10,000 | `N_PERMUTATIONS = 10_000` in `src/procrustes.py` | **Yes** |
| Bootstrap iterations | 100 (robustness), 1,000 (ranking) | `N_BOOTSTRAP = 100` in `scripts/07_bootstrap.py` | **Yes** (robustness); ranking script not verified |
| Bootstrap subsampling | 50% | `SUBSAMPLE_FRACTION = 0.5` | **Yes** |
| Alpha | 0.01 (primary), 0.05 (secondary) | Hardcoded in analysis scripts | **Yes** |
| Census version | 2025-11-08 | See Section 3.1 | **Partial** |
| SVD rotation | R = VDU^T, det(R) = +1 | Implemented in `src/procrustes.py` | **Yes** |
| Ortholog filter | Strict 1:1, BioMart release 115 | `src/data_loader.py` (1:1 filter present) | **Yes** |

### 2.2 Parameters NOT documented in manuscript

| Parameter | Value in Code | Impact |
|---|---|---|
| QC min genes/cell | 200 (`src/qc.py`) | Low – standard default |
| HVG count for QC UMAP | 2,000 (`src/qc.py`) | Low – QC visualization only, not used in Procrustes |
| UMAP neighbors | 15 (`src/qc.py`) | Low – visualization only |
| BioMart retry logic | 3 retries, 30s delay | None – operational parameter |
| Bootstrap permutations per iter | 1,000 (reduced from 10K) | **Moderate** – manuscript says "1,000 permutations per iteration" but this is a speed optimization not stated in results |

### 2.3 Reimplementation feasibility

The Methods section describes the Procrustes algorithm at mathematical detail sufficient for reimplementation: centering, SVD decomposition, reflection correction (det check), scaling, distance formula. The permutation null (label shuffle) is clearly specified. A competent computational biologist could reimplement the core pipeline from the manuscript alone without seeing the code.

**However:** Some secondary analyses are described only at a high level:
- Krzanowski S statistic for ellipsoid alignment: formula not given (reference to Krzanowski 1979 suffices for an expert)
- Liang-Wagner treeness: delta statistic described by reference only
- CellMarker enrichment: the shift from top-20 to top-50 genes is documented as data-informed, which is transparent but introduces a researcher degree of freedom

---

## 3. CODEBASE REPRODUCIBILITY

### 3.1 Census version pinning

**CRITICAL FINDING:** Three scripts call `cellxgene_census.open_soma()` without specifying `census_version="2025-11-08"`:

| Script | Role | Severity |
|---|---|---|
| `scripts/08_scaled_procrustes.py` (line 215) | **Primary 35-type analysis** | **HIGH** – this is the script that produces the headline result |
| `scripts/33_cellhint_replication.py` (line 324) | CellHint replication | Medium |
| `scripts/15_hsc_10x_validation.py` (line 116) | HSC validation | Low |
| `scripts/08_scaled_procrustes.py` (line 215) | (also used by Census for 35-type download) | **HIGH** |

Without a version pin, `open_soma()` uses the latest Census version, which may differ from v2025-11-08 and could change cell counts, annotations, or available datasets. The manuscript claims "Census version 2025-11-08" for all analyses.

**Mitigating factor:** The script saves normalized `.h5ad` checkpoint files. If these are distributed or previously generated, the Census call is skipped on rerun. However, a from-scratch reproduction would not be deterministic.

All other Census-calling scripts (17 checked) correctly pin to `"2025-11-08"`.

### 3.2 Random seed coverage

| Component | Seed Set? | Method |
|---|---|---|
| Permutation tests | **Yes** | `np.random.RandomState(42)` or `np.random.default_rng(42)` |
| PCA in `src/procrustes.py` | **Yes** | `PCA(random_state=RANDOM_SEED)` |
| PCA in `src/qc.py` | **No** (implicit) | `sc.tl.pca()` uses ARPACK solver which is deterministic, but no explicit `random_state` |
| Cell subsampling | **Yes** | `np.random.default_rng(42)` in `data_loader.py` |
| Bootstrap | **Yes** | Iteration-indexed seeds (`i` for subsampling, `i+1000` for permutations) |
| LOOCV | **Partial** | Hardcoded `random_state=42` literal in `08_loocv.py` instead of imported constant |
| UMAP | **Yes** | `random_state=42` in `src/qc.py` |

The implicit PCA seed in `src/qc.py` is a minor concern – Scanpy's ARPACK solver is deterministic for `n_comps < n_features`, but this is not explicitly documented or guaranteed across platforms.

### 3.3 Software dependency pinning

| Package | Manuscript Version | requirements.txt | Match? |
|---|---|---|---|
| Python | 3.12.12 | Not in requirements.txt (expected) | N/A |
| NumPy | 2.4.3 | `numpy==2.4.3` | **Yes** |
| SciPy | 1.17.1 | `scipy==1.17.1` | **Yes** |
| pandas | 2.3.3 | `pandas==2.3.3` | **Yes** |
| Scanpy | 1.12 | `scanpy==1.12` | **Yes** |
| AnnData | 0.12.10 | `anndata==0.12.10` | **Yes** |
| scikit-learn | 1.8.0 | `scikit-learn==1.8.0` | **Yes** |
| statsmodels | 0.14.6 | **MISSING** | **Gap** |
| cellxgene-census | 1.17.0 | `cellxgene-census==1.17.0` | **Yes** |
| SAMap | 1.0.14 | `samap==1.0.14` | **Yes** |
| Matplotlib | 3.10.8 | `matplotlib==3.10.8` | **Yes** |
| Seaborn | 0.13.2 | `seaborn==0.13.2` | **Yes** |
| pybiomart | Not in manuscript | `pybiomart>=0.2` (unpinned) | **Gap** – used for ortholog queries |
| gseapy | Not in manuscript | `gseapy>=1.0` (unpinned) | **Gap** – used for GO enrichment |

**Finding:** `statsmodels==0.14.6` is listed in the manuscript Key Resources Table but missing from `requirements.txt`. It is imported in at least 2 scripts. A fresh `pip install -r requirements.txt` would not install it (it may come as a transitive dependency, but this is not guaranteed).

### 3.4 Execution order and automation

- **README.md** provides a 7-step execution table with numbered scripts, descriptions, and approximate runtimes. This covers the primary result pipeline.
- **No Makefile, Snakemake, or workflow runner** exists. Reproduction requires manual execution of ~7+ scripts in sequence.
- The README pipeline covers only the primary result (steps 1-7). The full set of results in the manuscript requires running 80+ additional scripts (replications, negative controls, sensitivity analyses, disease analyses, mechanistic nulls, etc.) with **no documented execution order** for these additional analyses.
- Script numbering is inconsistent (e.g., 01, 02, 04, 04b, 04c, 06, 07, 08, 08b, 09, 10, 11, 12, 13, 14, ..., 35). Some numbers are reused with different prefixes. Several scripts have non-numeric names (e.g., `test_*`, `step_*`, `v1_*`).

### 3.5 Code duplication risk

The Procrustes distance function is implemented in three places:
1. `src/procrustes.py` – canonical module
2. `scripts/permutation_1M.py` – standalone copy with its own determinant-sign handling (written as overflow-safe; the matrix is orthogonal, so no overflow arises)
3. `scripts/permutation_1M_independent_pca.py` – another standalone copy

The standalone copies include defensive numerical handling (`np.clip` on determinants) that the canonical module does not. If these diverge, results could differ silently. The 1M scripts also contain hardcoded assertions against expected observed distances (`~61.153` and `~52.716`), which provide a correctness check but make the scripts brittle to upstream changes.

### 3.6 External API dependencies

| Service | Used For | Reproducibility Risk |
|---|---|---|
| Ensembl BioMart (live) | Ortholog mapping | Medium – cached CSV mitigates |
| Enrichr (via gseapy) | GO enrichment | Medium – library version pinned to "GO_Biological_Process_2023" but server content not version-controlled |
| CZ CELLxGENE Census API | Primary data access | Low if version pinned; **High for 3 unpinned scripts** |

---

## 4. PER-CLAIM REPRODUCIBILITY ASSESSMENT

### 4.1 Primary result (obs/null = 0.522, p < 10^-6)

| Criterion | Status |
|---|---|
| Data source specified | Yes (Census v2025-11-08) |
| Method reimplementable from text | Yes |
| Parameters reported | Yes (33 PCs, 1M permutations, seed=42) |
| Software versions pinned | Yes (in manuscript + requirements.txt) |
| Exact number reproducible | **Likely yes from cached data; uncertain from scratch** due to unpinned Census in `08_scaled_procrustes.py` |

### 4.2 Bootstrap (CV = 0.004, 100/100 significant)

| Criterion | Status |
|---|---|
| Parameters reported | Yes (100 iterations, 50% subsampling, 1,000 permutations/iter) |
| Seeds reported | Yes (seed=42 global; iteration-indexed subseed) |
| Exact number reproducible | Yes if using saved .h5ad files |

### 4.3 LOOCV (35/35 better than chance, mean ratio 0.4201)

| Criterion | Status |
|---|---|
| Method described | Yes (PCA refit per fold, 95% threshold) |
| Parameters reported | Yes |
| Minor code issue | Hardcoded literal `random_state=42` instead of shared constant |
| Exact number reproducible | Yes |

### 4.4 Replication datasets (Sun2023, PanSci, CellHint, Pan-Census)

| Criterion | Status |
|---|---|
| Data sources | Yes for Sun2023 (OMIX/GSA), PanSci (GEO), CellHint (Census) |
| Pan-Census | "22 types from 15 independent datasets" – individual dataset IDs not tracked |
| Parameters | Yes (10,000 permutations, seed=42) |
| Exact numbers reproducible | Likely yes from cached data; Pan-Census dataset provenance incomplete |

### 4.5 Three-species extension (macaque, mouse lemur)

| Criterion | Status |
|---|---|
| Macaque data | RIRA (GEO GSE277821) + Qu et al. (GEO GSE196791) – both accessible |
| Mouse lemur data | "CZ CELLxGENE Discover" – **no specific accession ID** |
| Ortholog mapping | BioMart release 115, macaque + mouse lemur cross-references |
| Exact numbers reproducible | **Uncertain for mouse lemur** due to missing data accession |

### 4.6 Negative controls

| Criterion | Status |
|---|---|
| Human-vs-human (6-type) | Reproducible – Census v2025-11-08 pinned in `04c_negative_control.py` and `09_negative_control_v2.py` |
| Within-species tissue pairs (24 pairs) | Method described but 24 specific tissue pairs not enumerated in manuscript |
| Self-comparison (50 splits) | Method described; seed not stated for split selection |

### 4.7 Two-layer decomposition (Layer 1 p=0.0001, Layer 2 p=0.0001)

| Criterion | Status |
|---|---|
| Krzanowski S statistic | Referenced (Krzanowski 1979) but formula not reproduced in Methods |
| Centroid-optimal rotation | Described conceptually |
| Exact reimplementation | Would require reading the cited reference |

### 4.8 Ten mechanistic null tests

| Criterion | Status |
|---|---|
| Housekeeping, TF complexity, variance, etc. | Gene sources described but database versions not pinned (MSigDB, DoRothEA, ChEMBL, DrugBank) |
| STRING PPI | v12.0 specified, 15,977/16,959 genes mapped |
| phastCons | "UCSC 100-vertebrate alignment" – no download date |
| ENCODE H3K27ac | No experiment IDs given |
| Reproducible | **Partially** – Spearman correlations are simple to verify if gene sets are available, but gene set provenance is incomplete for 5/10 tests |

### 4.9 CellMarker validation (4.49-fold, p = 2.10 × 10^-13)

| Criterion | Status |
|---|---|
| CellMarker 2.0 | Downloaded 2026-03-16; no database version number |
| Method | Hypergeometric test, fully described |
| Researcher degree of freedom | Top-20 → top-50 gene threshold change documented as data-informed |
| Reproducible | **Partially** – CellMarker database content at download date is not archived |

### 4.10 Cancer/COVID deformation

| Criterion | Status |
|---|---|
| Data | Census v2025-11-08 (pinned in relevant scripts) |
| Parameters | 14 types, enrichment method described |
| Individual study accessions | **Not tracked** ("can be recovered from Census metadata") |
| Reproducible | Yes for main results; study-level provenance incomplete |

### 4.11 DILI analysis

| Criterion | Status |
|---|---|
| Pre-registration | Git commit 415c806 – verifiable |
| Data | DILIrank (FDA URL given), L1000 (LINCS) |
| Hard Abort framework | Described in detail with numbered tests |
| Reproducible | **Yes** – pre-registered, well-documented |

### 4.12 Simulation study

| Criterion | Status |
|---|---|
| Parameters | 500 genes, 50 latent factors, 100 replicates, type counts 15/25/35, cell counts 50-2000, signal 0.5-10.0 |
| Noise model | "Log-normal noise with 10-fold range, CLT centroid shortcut" |
| Reproducible | **Yes** – fully parameterized, deterministic simulation |

---

## 5. NUMERICAL DISCREPANCIES NOTED

| Claim | Manuscript | Issue |
|---|---|---|
| "0.05^4 = 0.000125" | Line ~110 (multiple comparison note) | 0.05^4 = 6.25 × 10^-6, not 0.000125. If "four" means primary + three replications, this is correct for the count but wrong for the arithmetic. If three independent tests: 0.05^3 = 0.000125. Likely the exponent should be 3, not 4. |
| obs/null in README vs manuscript | README: "obs/null = 0.522, p < 10^-6 (35 cell types, 10,000 permutations)" | The primary result used 1,000,000 permutations, not 10,000. The README's "10,000" matches the default module constant, not the actual primary analysis. |

---

## 6. CATEGORIZED FINDINGS

### 6.1 CRITICAL (would prevent reproduction)

| # | Finding | Impact |
|---|---|---|
| C1 | `08_scaled_procrustes.py` calls `open_soma()` without `census_version="2025-11-08"` | Primary result may not reproduce from scratch on a future Census release |
| C2 | Code/data repository URLs are placeholders | Reviewer cannot access code |

### 6.2 MAJOR (would complicate reproduction)

| # | Finding | Impact |
|---|---|---|
| M1 | `statsmodels==0.14.6` listed in manuscript but missing from `requirements.txt` | `pip install -r requirements.txt` would not install it; scripts importing statsmodels would fail |
| M2 | No execution order documented for ~80 scripts beyond the primary 7-step pipeline | Reproducing replication/sensitivity/disease analyses requires guessing script order |
| M3 | Ensembl BioMart queried via live API (not versioned archive) in `src/data_loader.py` | From-scratch ortholog mapping could differ if Ensembl updates release 115 annotations |
| M4 | Database versions unspecified for 5/10 mechanistic null tests (MSigDB, DoRothEA, ChEMBL, DrugBank, ENCODE) | Cannot guarantee identical gene sets |
| M5 | Tabula Microcebus has no specific data accession beyond "CZ CELLxGENE Discover" | Mouse lemur results cannot be independently sourced |
| M6 | Two Census-calling scripts unpinned (`33_cellhint_replication.py`, `15_hsc_10x_validation.py`) | Replication results may not reproduce from scratch |

### 6.3 MINOR (cosmetic or low-impact)

| # | Finding | Impact |
|---|---|---|
| m1 | `pybiomart>=0.2` and `gseapy>=1.0` unpinned in requirements.txt | Could introduce behavioral differences; both are noted as "not listed in manuscript methods" |
| m2 | `src/qc.py` PCA lacks explicit `random_state` (relies on implicit ARPACK determinism) | Unlikely to affect results; platform-dependent in theory |
| m3 | `08_loocv.py` uses literal `42` instead of imported `RANDOM_SEED` constant | Functionally identical; maintenance concern only |
| m4 | Code duplication of Procrustes function across 3 files with slight divergence | Could silently produce different results if one copy is patched and others are not |
| m5 | No Makefile/Snakemake workflow for automated end-to-end reproduction | Manual execution required; error-prone but documented in README |
| m6 | README states "10,000 permutations" for primary result; manuscript and actual analysis used 1,000,000 | Inconsistency in documentation |
| m7 | CellMarker 2.0 archived by download date (2026-03-16) but database version not recorded | Minor – download date is sufficient for most practical purposes |
| m8 | No OS/hardware specification in manuscript (runs on Apple Silicon / arm64 macOS; recorded in reproduce/environment_ground_truth.txt) | Computational results should be platform-independent; no GPU-dependent steps |
| m9 | Pan-Census "22 types from 15 independent datasets" – individual study accessions not enumerated | Provenance incomplete but recoverable from Census metadata |
| m10 | 0.05^4 arithmetic error in multiple comparison justification | Does not affect any reported result |

---

## 7. REPRODUCIBILITY SCORECARD

| Dimension | Score | Notes |
|---|---|---|
| **Data accessibility** | 8/10 | All major datasets have accessions; mouse lemur and Pan-Census provenance gaps |
| **Method specification** | 9/10 | Core algorithm described at reimplementation level; some secondary analyses by reference only |
| **Parameter documentation** | 9/10 | All critical parameters reported; 5 database versions missing for mechanistic nulls |
| **Software pinning** | 8/10 | Core packages pinned with ==; statsmodels missing; 2 utility packages unpinned |
| **Random seed discipline** | 9/10 | Seed=42 consistently used; one implicit PCA seed; one hardcoded literal |
| **Code-manuscript consistency** | 7/10 | Census version pin missing in 3 scripts including primary; README permutation count discrepancy |
| **Execution documentation** | 6/10 | Primary pipeline documented; 80+ secondary scripts undocumented |
| **Automation** | 4/10 | No Makefile, Snakemake, or CI; manual execution only |

**Weighted overall: B (Mostly Reproducible)**

---

## 8. RECOMMENDATIONS (in priority order)

1. **Pin Census version in `08_scaled_procrustes.py`** – change `open_soma()` to `open_soma(census_version="2025-11-08")`. Same for `33_cellhint_replication.py` and `15_hsc_10x_validation.py`.
2. **Add `statsmodels==0.14.6`** to `requirements.txt`.
3. **Fill in repository URL** in the Data Availability section before submission.
4. **Fix README** to say "1,000,000 permutations" for the primary result.
5. **Use Ensembl archive URL** (`https://mar2026.archive.ensembl.org`) in `src/data_loader.py` instead of the live endpoint, or document that the cached CSV is the reproducibility anchor.
6. **Add a `SCRIPTS.md`** or second README section documenting execution order for all analyses beyond the primary 7-step pipeline.
7. **Specify database versions** for MSigDB, DoRothEA, ChEMBL, DrugBank, and ENCODE experiment IDs.
8. **Add a Tabula Microcebus data accession** (DOI or CZ CELLxGENE dataset ID).
9. **Fix the 0.05^4 arithmetic** in the multiple comparison note (should be 0.05^3 = 0.000125 for three replications, or 0.05^4 = 6.25 × 10^-6 for four).
10. **Consider adding a Makefile or Snakemake workflow** for end-to-end automated reproduction.

---

## APPENDIX A: Census Version Pin Audit (All Scripts)

Scripts that correctly pin `census_version="2025-11-08"`:
- `01_download_data.py` ✓
- `04c_negative_control.py` ✓
- `09_negative_control_v2.py` ✓
- `10_cancer_download.py` ✓
- `12_cancer_scaled.py` ✓ (3 call sites)
- `13_covid_procrustes.py` ✓ (4 call sites)
- `13_replication_inventory.py` ✓
- `13_t1a_hca_download.py` ✓
- `30_fourth_replication_candidates.py` ✓
- `30b_candidate_deep_dive.py` ✓
- `30c_hepatocyte_source_scout.py` ✓
- `31_andrews_replication.py` ✓
- `32_ts2_verification.py` ✓
- `disease_inventory_covid.py` ✓
- `disease_inventory_lupus.py` ✓
- `t1a_hca_feasibility.py` ✓
- `t1a_hca_curated_set.py` ✓
- `t1a_hca_curated_nohcl.py` ✓
- `08_cell_type_inventory.py` ✓

Scripts that call `open_soma()` **without** version pin:
- `08_scaled_procrustes.py` ✗ ← **PRIMARY ANALYSIS**
- `33_cellhint_replication.py` ✗
- `15_hsc_10x_validation.py` ✗

## APPENDIX B: Procrustes Distance Function Locations

| Location | Determinant Handling | Notes |
|---|---|---|
| `src/procrustes.py` | Standard `np.linalg.det()` | Canonical module |
| `scripts/permutation_1M.py` | `np.clip()` overflow-safe | Standalone; hardcoded assertion on expected distance |
| `scripts/permutation_1M_independent_pca.py` | `np.clip()` overflow-safe | Standalone; hardcoded assertion on expected distance |

The overflow-safe handling in the 1M scripts is a defensive improvement not present in the canonical module. It cannot cause divergence: the matrix whose determinant is taken is V @ U.T from an SVD, which is orthogonal at any k, so its determinant is ±1 and no overflow arises in either implementation. The inconsistency is stylistic rather than numerical.
