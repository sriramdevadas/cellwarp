# Reproducing the CellWarp Paper

## Requirements
- Python 3.12 (the project pins `>=3.12,<3.13`)
- ~6 GB disk space for core data (Tier 1, downloaded automatically)
- ~15 GB additional disk space for supplementary data (Tier 2)
- Estimated runtime: ~4 hours total (~1 hour core + ~3 hours supplementary)
- Internet connection required for initial data download

## Fast path (no download, a few minutes)

For a turnkey check that the headline result reproduces with **no atlas
download**, run the deposited-centroids demo:

```bash
# from a fresh clone -- complete, self-contained sequence:
python3.12 -m venv .venv && source .venv/bin/activate # 1. create + activate env (Python 3.12 required)
pip install -e .                                     # 2. install cellwarp (base, no compiler)
python reproduce/fast_path.py                        # 3. run (a few minutes, no network)
```

(Already have the environment from the Steps below? Just run
`python reproduce/fast_path.py`. Prefer not to activate? Use the venv
interpreter directly: `.venv/bin/python reproduce/fast_path.py`.)

It loads the deposited PCA centroids
(`output/phase2/scaled_35types/pca_centroids_35.npz`), runs the 1,000,000-
permutation Procrustes test (seed 42), and prints PASS/FAIL against the
published primary result: **obs/null = 0.522 (median denominator), p < 1e-6**.
Runtime a few minutes, no network required, read-only (it does not overwrite any
committed output). The fast path needs only the base install:
`pip install -e .` pulls just `numpy/scipy/scikit-learn/pandas/anndata/h5py`
as prebuilt wheels (no C compiler, no SAMap). The full download-based
pipeline below (`pip install -e ".[lock]"`) reproduces the complete set of
analyses.

## Steps

> Requires **Python 3.12** (`>=3.12,<3.13`); the stock `python3` on Ubuntu 22.04 / macOS is older and the install will fail against it. Get 3.12 first (Ubuntu: deadsnakes PPA; macOS: `brew install python@3.12`), or skip host Python entirely with the Docker image (main [README → Reproduce in Docker](../README.md#reproduce-in-docker)).

```bash
git clone https://github.com/sriramdevadas/cellwarp.git
cd cellwarp
python3.12 -m venv .venv
source .venv/bin/activate
pip install -e ".[lock]"

# Run full reproduction (Tier 1 core + Tier 2 supplementary):
bash reproduce/run_all.sh

# Or run just the core pipeline (Tier 1 only):
# Stop run_all.sh after "TIER 1 COMPLETE" with Ctrl+C,
# or run the 8 core scripts individually (see below).
```

> **Slim-image note:** the full `.[lock]` install (and `cellwarp[samap]`)
> builds `hnswlib` (a SAMap dependency) from source; on a minimal image
> without build tools, run `apt-get install -y build-essential` first, or
> use the full `python:3.12` image. The fast-path base install above needs
> none of this.

## What each tier does

**Tier 1 -- Core pipeline (main result):**
Downloads primary human/mouse atlas data from CELLxGENE Census, runs QC
and normalization, identifies qualifying cell types, executes 35-type
Procrustes analysis with 10,000 permutations, and generates main results.

**Tier 2 -- Supplementary analyses:**
Runs all supplementary analyses: independent PCA sensitivity check,
simulation study, parameter and protocol sensitivity, expanded negative
controls, bootstrap ranking stability, CellHint investigation, SAMap
validation, CellMarker validation, biological predictors, cross-atlas
replication, disease deformation, and CPC1 driver gene extraction.

## Validation

After the pipeline completes, `reproduce/validate.py` automatically
compares all key output statistics against values reported in the
manuscript and prints a pass/fail summary.

## Figure-to-script mapping

See `reproduce/figure_script_map.md` for a complete table showing which
script generates each figure and table in the paper.

## SAMap validation (optional, requires PyTorch)

SAMap has heavy dependencies (including PyTorch) and is excluded from
the default install. To run the SAMap validation (Figure S5):

```bash
pip install cellwarp[samap]
# SAMap step will then run automatically in reproduce/run_all.sh
```

## Seed and determinism

All stochastic operations use seed 42. On the same platform with the same
library versions, results reproduce to ~6 significant figures. Across
different platforms or library versions, two kinds of small variation can
occur:

- **Floating-point / BLAS differences** shift numeric results around the
  6th-7th significant figure (e.g. obs/null ratios) - well within reported
  precision.
- **Stochastic / environment-sensitive steps** - Leiden cell-type assignment
  (e.g. the Sun2023 replication) and SAMap graph mapping - can vary slightly
  across platforms and library versions, occasionally nudging a per-type cell
  count near a hard threshold or a reported correlation in its third decimal.

In all cases, analyzed cell-type membership, every headline result, and all
scientific conclusions are unaffected.
