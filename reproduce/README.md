# Reproducing the CellWarp Paper

## Requirements
- Python 3.12 (the project pins `>=3.12,<3.13`)
- ~6 GB disk space for core data (Tier 1, downloaded automatically)
- additional disk space for the optional Tier-2 datasets: see
  [DATA_SOURCES.md](../DATA_SOURCES.md) for per-dataset sizes
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
# Stop run_all.sh after "TIER 1 COMPLETE" with Ctrl+C, or run the eight
# Tier-1 steps individually -- they are the [1/8]..[8/8] steps at the top
# of reproduce/run_all.sh.
```

> **Slim-image note:** the full `.[lock]` install (and `cellwarp[samap]`)
> builds `hnswlib` (a SAMap dependency) from source; on a minimal image
> without build tools, run `apt-get install -y build-essential` first, or
> use the full `python:3.12` image. The fast-path base install above needs
> none of this.

## What each tier does

**Tier 1 -- Core pipeline (main result):** eight steps, banners `[1/8]`
through `[8/8]`. Downloads primary human/mouse atlas data from CELLxGENE
Census, runs QC and normalization, identifies qualifying cell types, then runs
the 35-type Procrustes analysis at 10,000 permutations (step 4) followed by a
separate 1,000,000-permutation test on the same comparison (step 5), and
finishes with GO enrichment, bootstrap robustness, and LOOCV.

**Tier 2 -- Supplementary analyses:** everything after `TIER 1 COMPLETE`.
Independent PCA sensitivity, simulation study, parameter and protocol
sensitivity, expanded negative controls, bootstrap ranking stability, CellHint
investigation, SAMap validation, CellMarker validation, biological predictors,
cross-atlas replication, disease deformation, and CPC1 driver gene extraction.
`reproduce/run_all.sh` is the authority for what runs and in what order. Some
of these belong to analyses this repository holds but the paper does not
report; `SCOPE.md` classifies every script, analysis directory, and output
directory on that boundary.

Steps whose optional inputs are absent are skipped by the `require_data` guard
where one is present; not every step is guarded, so a standalone Tier-2 run on
a fresh clone can error rather than skip. Run Tier 1 first -- several Tier-2
steps read intermediates it produces.

## Validation

After the pipeline completes, `reproduce/validate.py` reads the tracked
output artifacts and checks each recorded statistic against its expected
value, then prints a pass/fail summary. Each check carries a `paper_ref`
naming where the value appears; some are retained tooling for analyses the
current paper does not report.

## Figure-to-script mapping

See `reproduce/figure_script_map.md` for a complete table showing which
script generates each figure and table in the paper.

## SAMap validation (optional, requires PyTorch)

SAMap has heavy dependencies (including PyTorch) and is excluded from
the default install. To enable the SAMap step:

```bash
pip install cellwarp[samap]
# SAMap step will then run automatically in reproduce/run_all.sh
```

`run_all.sh` gates that step on an `import samap` probe and skips with a
message when the package is absent, so install into the same environment you
run `run_all.sh` from.

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

## Counting words in the submitted texts

There are three different word counts of the manuscript in circulation, and
only one of them is asserted by anything. Post-Tier-1 values:

| quantity | value | who computes it | gated? |
|---|---|---|---|
| `EXPECTED_JOINED_WORDS` | 13612 | `docs/submission/plosone/build_manuscript_docx.py` | **yes** — the builder aborts on mismatch |
| `wc -w` on `manuscript_combined.txt` | 13775 | shell | no |
| `wc -w` on `S1_Text.txt` / `S2_Text.txt` | 4947 / 975 | shell | no |

`EXPECTED_JOINED_WORDS` is smaller because it counts the **158 content lines
joined**, i.e. the manuscript body as the DOCX renders it: the section banners,
the `====` rules and the blank structural lines that `wc` sees are not content
and are not counted. It is not a "wrong" `wc`; it is a different object. Do not
reconcile the two, and do not copy a `wc` number into the builder.

`wc -w` is in principle both locale- and implementation-dependent: GNU
coreutils `wc` splits on `iswspace()` in a UTF-8 locale and on `isspace()` over
bytes under `LC_ALL=C`, which can disagree on a file carrying non-ASCII
characters, and these texts carry plenty (82 en dashes, 42 rho, 59 minus
signs). **On the reference machine that divergence does not appear.** Measured
against these files, macOS BSD `/usr/bin/wc` returns 13775 identically under
`LC_ALL=C`, `POSIX`, `en_US.UTF-8`, `UTF-8` and `en_US.ISO8859-1`, and Python's
`str.split()` agrees. The only whitespace code points in the file are U+0020
and U+000A, and there is no U+00A0, which is why the locales cannot disagree
here. A C-versus-UTF-8 split was reported during the Tier-1 pass (13606 against
13524 at the pre-edit HEAD); the 13606 reproduces exactly under every locale
above and 13524 reproduces under none, so it is recorded here as unconfirmed
rather than as a property of these files. GNU `wc` is not installed on the
reference machine and was not tested.

Practical consequence: if a submission form wants a word count, say which tool
produced it. The values above are BSD `wc -w`. Only `EXPECTED_JOINED_WORDS` is
a gate, and it is the only one that must not drift silently.
