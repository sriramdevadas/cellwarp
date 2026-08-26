# CellWarp

Cross-species geometric morphometrics in transcriptomic space, inspired by D'Arcy Thompson's *On Growth and Form*.

> **Reading this from a Zenodo archive?** That archive is a **frozen snapshot** of one commit, not
> the living project. The current version is at
> **<https://github.com/sriramdevadas/cellwarp>**, and the Zenodo record's own concept DOI
> [10.5281/zenodo.20735611](https://doi.org/10.5281/zenodo.20735611) always resolves to the newest
> deposited version. If anything below does not work, check the repository before working around it:
> reader-path fixes land there first and reach a snapshot only at the next deposit.
>
> Two things that bite only on the archive route. **Versions 1 and 2 of the zip have no wrapper
> directory**, so all 1,094 files land wherever you unpack them — unpack those into an empty
> directory rather than into `~/Downloads`; version 3 onward carries a wrapper and does not have
> this problem. And extracting needs `unzip`, which a minimal Linux image (including stock
> `ubuntu:24.04`) does not ship: `sudo apt-get install -y unzip` first, or use
> `python3 -m zipfile -e <file>.zip .`, which needs nothing beyond the Python you already have.

## Quick Start

```bash
git clone https://github.com/sriramdevadas/cellwarp.git && cd cellwarp
```

**Verify the deposit in Docker** (no host Python needed): it builds Ubuntu 22.04 + Python 3.12 and runs the four reproduction gates plus the headline fast-path. See [Reproduce in Docker](#reproduce-in-docker).

### Install prerequisites first

Every local path in this file — the gates, the fast path and the full reproduction — opens with
`python3.12 -m venv .venv`. **On a stock Linux image that command fails**, because `python3.12-venv`
is a package separate from the interpreter. Install prerequisites before either block below, not
after one has failed.

Two package sets, and which one you need depends on where you are heading:

- **Gates, fast path, library** (the `.[dev]` block just below): `python3.12` and `python3.12-venv`.
  Nothing is compiled, so no compiler and no headers.
- **Full reproduction** ([Setup](#setup-requires-python-312), the `.[lock,dev]` block): the same two
  **plus** `python3.12-dev` and a compiler, because `[lock]` installs SAMap and its dependency
  `hnswlib` is built from source.

**Ubuntu 24.04** ships Python 3.12.3, which satisfies the pin, so no PPA is needed:

```bash
sudo apt-get update
sudo apt-get install -y python3.12 python3.12-venv        # gates, fast path, library
sudo apt-get install -y python3.12-dev build-essential    # additionally, for the full reproduction
```

**Ubuntu 22.04** stocks Python 3.10, which fails the pin. Add the deadsnakes PPA first:

```bash
sudo add-apt-repository -y ppa:deadsnakes/ppa && sudo apt-get update
sudo apt-get install -y python3.12 python3.12-venv        # gates, fast path, library
sudo apt-get install -y python3.12-dev build-essential    # additionally, for the full reproduction
```

**macOS** needs neither a `-dev` package nor `build-essential`: the framework Python ships `venv` and
its own headers, and the Xcode command line tools supply the compiler. One line covers both paths:

```bash
brew install python@3.12
```

**RHEL-family (Fedora, Rocky, AlmaLinux) — untested here.** No run of this project has been made on
these distributions; the equivalent packages are named below, but unlike the Ubuntu and macOS lines
they have not been executed:

```bash
sudo dnf install -y python3.12 python3.12-devel gcc gcc-c++
```

**Verify before continuing.** Both lines are platform-independent and must print `OK`; a line that
prints nothing is the prerequisite you are missing, and it is far cheaper to see it here than at the
`venv` step:

```bash
python3.12 -m venv /tmp/cw-check && rm -rf /tmp/cw-check && echo "OK: venv works"
gcc --version >/dev/null 2>&1 && python3.12-config --includes >/dev/null 2>&1 && echo "OK: compiler and headers"
```

The first line is required for every path; the second only for the full reproduction.

> **If you already ran `python3.12 -m venv .venv` and it failed**, delete the directory before
> retrying: `rm -rf .venv`. A failed `venv` leaves `.venv/bin/python` behind while `pip` and
> `activate` are absent, so retrying over it fails a second and less obvious way. The error message
> says the same thing — "recreate your virtual environment" — and this is what it means.

### Install locally

With the prerequisites above in place (and their verification printing `OK`):

```bash
python3.12 -m venv .venv
source .venv/bin/activate          # or call .venv/bin/python directly (the gates assume it)
pip install -e ".[dev]"
```

**Which install?** `.[dev]` is enough for the four reproduction gates, the test suite and the
library snippet below; it is what the Docker image installs. **Measured, not assumed** — in a clean
Python 3.12 environment carrying nothing but `pip install -e ".[dev]"`: Gate 1 232/232, Gate 2 195
passed, Gate 3 30/30 pairs, Gate 4 3 of 3. That sentence was false until 2026-08-25, when the first
`docker build` ever executed failed on a missing `openpyxl`; the fix was to add it to `[dev]`, where
a test dependency belongs, rather than to the image. Reproducing the paper end to end
additionally needs `[lock]`, which pins the exact versions that produced the published numbers:
`pip install -e ".[lock,dev]"`, as in [Setup](#setup-requires-python-312). `[lock]` on its own is
not sufficient for the gates, because Gate 2 is `pytest` and `[lock]` does not install it.

If `pip install` fails, you're on the wrong Python. The pin is two-sided, so newer interpreters are rejected as well as older ones: 3.11 and below give `No matching distribution found for numpy…`, and 3.13 and above give `Package 'cellwarp' requires a different Python`. Use 3.12; [Setup](#setup-requires-python-312) has the detail.

Once installed, use the library on the deposited centroids:

```python
import numpy as np
from cellwarp.procrustes import procrustes_align, permutation_test

d = np.load("output/phase2/scaled_35types/pca_centroids_35.npz")
human, mouse = d["human"], d["mouse"]   # 35 matched cell-type centroids in 33-D PCA space

result = procrustes_align(human, mouse)              # align mouse → human (rotation, scaling, translation)
p_value, null_dist = permutation_test(human, mouse)  # geometric-coherence significance vs a label-permutation null
```

To reproduce the paper, see [Reproducing the Paper](#reproducing-the-paper), the full native-host pipeline (atlas download + all analyses, including the conda / `environment.yml` route), or [Reproduce in Docker](#reproduce-in-docker) for gate + fast-path verification without the large data download.

## Reproducing the Paper

This is the full reproduction (every analysis plus the CELLxGENE Census atlas download driven by `reproduce/run_all.sh`) on a native host with Python 3.12. For a quick, host-Python-independent **verification** of the gates and the headline fast-path (no large data download), use [Reproduce in Docker](#reproduce-in-docker) instead.

### Setup (requires Python 3.12)

This project pins **Python `>=3.12,<3.13`**. The stock `python3` on Ubuntu 22.04 (3.10) and macOS (system Python) is too old, and a current Homebrew default (3.14) is too new; the install **will fail** against either. **Do [Install prerequisites first](#install-prerequisites-first) before the block below.** Those commands are not repeated here, and on a stock image the block's first command fails without them. From `python3.12 -m venv` onward the block is kept contiguous, so the `pip install` line never runs against the wrong interpreter:

> **Read this before you run the block below.** It ends with `bash reproduce/run_all.sh`, which
> is a multi-hour job with a large memory footprint, and both of the things that stop it are
> described *after* the block rather than before it:
>
> - **Memory: the peak is 63.95 GiB**, in `[4/8]`, and it is inside Tier 1, so stopping after
>   `TIER 1 COMPLETE` does not avoid it. Provision against the **system-wide** figure, not the
>   per-process one: when `[4/8]` peaked the machine was using **65.14 GiB**. **128 GiB is the
>   recommendation and 64 GB is not** — that is measured, not cautioned: a 64 GiB machine would
>   have been OOM-killed on the third full run. A 32 GB instance is OOM-killed after about five
>   minutes. The full measurements, with the instance they were taken on, are in the
>   **Requirements** block further down this section.
> - **Build tools on Linux:** the `pip install -e ".[lock,dev]"` line builds `hnswlib` from
>   source and fails without a compiler and the CPython headers.
>   [Install prerequisites first](#install-prerequisites-first) installs them; do not skip it.
>
> Wall-clock for the whole pipeline has been measured at **3 h 09 m** and **9 h 40 m** on the same
> instance type. **Plan by where you are running, because that spread is one stage.**
> `[S23]` (`13_covid_procrustes.py`) issues a CELLxGENE Census query per cell-type × tissue with a
> 600 s timeout: **33 minutes from a cloud instance with no timeout firing, 6 h 29 m from a home
> connection with 14 of them firing.** Everything else is local compute and barely moves. An
> afternoon on a well-connected host; most of a day otherwise.
>
> **If you only need the gates, you do not need any of this.**
> [Reproduce in Docker](#reproduce-in-docker) runs all four reproduction gates and the headline
> fast-path in minutes, on any machine, with no large download and no memory requirement worth
> naming. Read that section's **Scope** paragraph before deciding: Docker verifies the gates and the
> fast path, **not** the pipeline that regenerates the numbers. If that is what you came to check,
> the block below is the one you want.

```bash
git clone https://github.com/sriramdevadas/cellwarp.git
cd cellwarp

# Prerequisites are NOT in this block: do "Install prerequisites first" above, and run its
# two verification lines, before pasting this. The next command fails on a stock image without them.
python3.12 -m venv .venv
source .venv/bin/activate          # or invoke .venv/bin/python directly (the gates assume .venv/bin/python)

# Recommended: [lock] pins the exact versions that produced the published results,
# and [dev] is required alongside it because Gate 2 is pytest, which [lock] does not install.
pip install -e ".[lock,dev]"
# Alternatives (same venv):
#   pip install -e ".[reproduce,dev]"  # bounded version ranges; numbers may differ slightly
#   pip install -r requirements.txt    # manuscript-anchored snapshot
#   conda env create -f environment.yml && conda activate cellwarp   # mirrors [lock] and includes pytest; Python 3.12.12

bash reproduce/run_all.sh

# A successful run leaves Gate 3 red until you refresh the packet mirrors.
# This is required, not optional -- see the note directly below the block.
python scripts/build_submission_packet.py --rebuild
```

**After a successful run, Gate 3 and part of Gate 2 fail until you rebuild the packet.** This is
not a failed reproduction: `run_all.sh` rewrites tracked canonical artifacts, and their byte-copies
under `docs/submission/figures_for_review/` are written only by
`scripts/build_submission_packet.py --rebuild`. Until you run it, the canonical and its mirror are
out of sync and `--verify` reports the drift. It is **one pair of thirty**, the same one on each
of the three full runs: `figS3_bootstrap_rankings.pdf` against `Figure_S3.pdf`. On the pristine archive all 30 pairs match,
which is why the gates are green *before* you reproduce anything.
[reproduce/README.md](reproduce/README.md#after-a-full-run-rebuild-the-packet-before-the-gates-mean-anything)
gives the full stage-by-stage account of which pairs drift and why.

`build-essential` and `python3.12-dev` are needed because `[lock]` installs SAMap, whose dependency `hnswlib` is built from source. This bites on Linux only. Neither the fast path nor the `.[dev]` gate install needs them. It fails two different ways and the message tells you which package is missing:

- **No compiler at all** (`build-essential` missing):
  `RuntimeError: Unsupported compiler -- at least C++11 support is needed!`
  This reads as though a compiler was found and judged inadequate. There is none. Install
  `build-essential`; do not go looking for a newer `g++`.
- **Compiler present, headers missing** (`python3.12-dev` missing):
  `fatal error: Python.h: No such file or directory`, then
  `ERROR: Could not build wheels for hnswlib`, exit 1.

Either way it fails at `pip install`, before anything has run.

If you see `No matching distribution found for numpy…`, you're on the wrong Python: that error means an interpreter older than 3.12; use 3.12.

The pin is two-sided, so a **newer** Python fails too, with a different message: `Package 'cellwarp' requires a different Python`. Homebrew now installs 3.14 by default, so on a fresh machine you are likelier to hit the ceiling than the floor. Measured with `pip install --dry-run`: exit 0 on 3.12.13, exit 1 on 3.14.3. There is no 3.13+ fallback, so install 3.12 explicitly (`brew install python@3.12` on macOS, the deadsnakes PPA on Ubuntu) and create the venv with `python3.12` as shown above.

**Requirements:** Python 3.12 (`>=3.12,<3.13`), ~6 GB disk (core), internet for initial data download, and enough memory for the tier you are running:

| what you are running | peak resident | leave yourself |
|---|---|---|
| fast path, or the Docker gate run | negligible | any machine |
| Tier 1 (`[1/8]`–`[8/8]`) | **63.95 GiB** (**65.14 GiB** system-wide) | 128 GiB recommended |
| full pipeline | **63.95 GiB**, the maximum is in Tier 1 | 128 GiB recommended |

**Three measurements, not an estimate: 51.4, 58.9 and 63.95 GiB** — same instance type, same commit, a 24% spread, the largest being the most recent (2026-08-26). Provision against the **system-wide** figure rather than the per-process one: at the moment `[4/8]` peaked at 63.95 GiB the machine was using **65.14 GiB**. **64 GB is not enough, and that is now measured rather than cautioned** — a 64 GiB machine would have been OOM-killed on the third run, not merely left tight; a 32 GB instance was OOM-killed at 30.2 GiB anon-rss. The peak is `[4/8]`, `scripts/08_scaled_procrustes.py`, which downloads 992,192 cells in order to keep 140,000. Because that step is **inside Tier 1**, stopping after `TIER 1 COMPLETE` does not avoid it. One other stage exceeds 40 GiB — `[S11]`, `scripts/33_cellhint_replication.py`, at **41.18 GiB** — but that one is in **Tier 2**, so stopping at `TIER 1 COMPLETE` does avoid it. Measured on AWS `r6i.4xlarge` (128 GiB) under Ubuntu 24.04.4: these are peaks to have headroom above, not thresholds to sit on. See [reproduce/README.md](reproduce/README.md) for the full measurement note.

**Tested on.** Fast path: macOS 15 on Apple silicon, and Ubuntu 24.04 on x86-64. Full pipeline: Ubuntu 24.04 on x86-64 only, three times, end to end with four green gates — most recently 2026-08-26 on AWS `r6i.4xlarge`, AMI `ami-052355af2a014bd2c` (`ubuntu-noble-24.04-amd64-server-20260714`, booting as Ubuntu 24.04.4), 3 h 09 m 10 s, exit 0. The full pipeline on macOS is untested. Windows is untested throughout, natively and in Docker alike. Nothing is known to be wrong with either; nothing has been checked.

The pipeline downloads human/mouse atlas data from CELLxGENE Census, runs QC, executes the 35-type Procrustes analysis with permutation testing, and validates all supplementary analyses. After completion, `reproduce/validate.py` checks key statistics against the manuscript values, and the frozen submission text is pinned by `reproduce/MANUSCRIPT_MD5`: the manuscript (`docs/submission/plosone/manuscript_combined.txt`) and both supporting-information texts (`S1_Text.txt`, `S2_Text.txt`), verified with `md5sum -c reproduce/MANUSCRIPT_MD5`, Gate 4.

## Reproduce in Docker

**Docker or native is a choice of scope, not of platform.** Both run on any host: the image needs only Docker, the native path only Python 3.12. Docker covers the four gates and the fast path with no large download; native covers the full pipeline. All four combinations are valid, so choose by what you want to check rather than by where you are running it.

A host-Python-independent **verification** path: it runs the four reproduction gates and the no-download fast-path; it does **not** run the full analysis pipeline (see Scope). The image is Ubuntu 22.04 + Python 3.12, and the gates run at build time, so a green build itself certifies that step:

```bash
git clone https://github.com/sriramdevadas/cellwarp.git
cd cellwarp
docker build -t cellwarp .        # installs Python 3.12 + deps, runs the 4 gates in-build
docker run  --rm cellwarp         # re-runs the 4 gates + the no-download fast-path
```

**What it costs, measured 2026-08-25.** A cold `docker build --no-cache` took **2 m 07 s** and the
image is **3.69 GB**; `docker run` took **5 m 28 s**, of which the 1,000,000-permutation fast path is
**5 m 22 s**. Both were measured on **arm64 macOS under emulation**, which is the slow case: the fast
path takes **1 m 49 s** natively on the same machine, so an amd64 reviewer should expect materially
less than these figures rather than more. The 3.69 GB is on top of a ~101 MB clone.

The image targets `linux/amd64` (deadsnakes ships Python 3.12 for Ubuntu 22.04 on amd64 only, so the platform is pinned for a reproducible build). On x86-64 hosts it builds natively; on Apple Silicon / ARM it builds under emulation automatically (no extra flags); `docker build`/`docker run` work exactly as written. The gates are arch-robust (md5 and packet checks are arch-independent; numeric checks carry tolerance).

**Scope:** the image reproduces the four gates (`reproduce/validate.py`, `pytest`, `build_submission_packet.py --verify`, `md5sum -c reproduce/MANUSCRIPT_MD5` over the manuscript and both SI texts) and the no-download fast-path (`obs/null ≈ 0.522` from deposited centroids). The **full paper reproduction**, which runs all analyses plus the external CELLxGENE Census atlas download via `reproduce/run_all.sh` (see [DATA_SOURCES.md](DATA_SOURCES.md)), is the native-host path only (see [Reproducing the Paper](#reproducing-the-paper)); there is no full-pipeline Docker route.

See [reproduce/README.md](reproduce/README.md) for detailed instructions.

## Repository Structure

```
src/cellwarp/       Library modules (Procrustes, data loading, QC, enrichment)
scripts/            Pipeline scripts (01-08 core, numbered analysis scripts)
analysis/           Supplementary analysis directories
reproduce/          Reproduction pipeline (run_all.sh, validate.py, config.py)
tests/              Unit and integration tests
docs/               Manuscript materials (current PLOS ONE submission: docs/submission/plosone/)
figures/            Paper figures
data/               Downloaded data (gitignored)
output/             Analysis outputs (gitignored)
```

## Data Sources

| Source | Version/Accession |
|---|---|
| CZ CELLxGENE Census | `census_version="2025-11-08"` |
| Tabula Sapiens (human) | Census v2025-11-08; figshare 10.6084/m9.figshare.14267219 |
| Tabula Muris Senis (mouse) | Census v2025-11-08; figshare 10.6084/m9.figshare.12654728 |
| Ensembl BioMart | Release 115 (GRCh38.p14) |

Primary atlas data (Tabula Sapiens, Tabula Muris Senis, CellHint Human) is
downloaded programmatically from CELLxGENE Census on first run. Optional
replication and validation datasets (Sun2023, PanSci, T3-E phastCons +
H3K27ac, DILIrank) require manual fetch from external sources; see
[DATA_SOURCES.md](DATA_SOURCES.md) for instructions. The pipeline
gracefully skips analyses whose optional data is absent.

## Deposit Artifacts

The repository tracks a small number of small CSV files that document specific
data subsets used in the analysis. These are reproducibility anchors that
travel with the deposit; they are not themselves analysis inputs that must be
regenerated.

| Artifact | Path | Rows | Purpose |
|---|---|---|---|
| Pan-Census replication manifest | `data/replication/pan_census_manifest.csv` | 15 | Lists the 15 CELLxGENE Census datasets (9 mouse + 6 human) used in the pan-Census replication (S2 Fig, panel F). Schema: `dataset_id`, `species`, `collection_id`, `collection_name`, `dataset_title`, `dataset_version_id`, `dataset_total_cell_count`, `collection_doi`, `citation`, `census_version`. |
| Tabula Microcebus metadata | `data/replication/tabula_microcebus_metadata.csv` | 1 | Records the CELLxGENE Discover deposit anchors for the Tabula Microcebus mouse lemur atlas (Ezran et al. 2025) used in the human–mouse lemur Procrustes comparison (Fig 1D). Schema: `dataset_id`, `species`, `collection_id`, `collection_name`, `dataset_title`, `dataset_total_cell_count`, `collection_doi`, `source`, `download_date`, `analysis_assay_filter`. Discover does not version-pin in the manner of Census, so download date is the reproducibility anchor. |

This table will grow as additional manifest CSVs are deposited.

## Figure-to-Script Mapping

See [reproduce/figure_script_map.md](reproduce/figure_script_map.md) for a complete
table mapping every figure and table in the paper to its generating script.

## Running Tests

```bash
pip install -e ".[dev]"     # already present if you installed ".[lock,dev]" above
pytest tests/ -v
```

195 tests covering Procrustes alignment, permutation calibration, end-to-end
integration, manifest sync, and the submission packet. The count is larger than
the number of test functions because `tests/test_table1_callouts.py`
parametrizes two checks over all 64 numbered S13 Table rows.

## Citation

If you use this work, please cite it according to the metadata in `CITATION.cff` (machine-readable Citation File Format) at the repository root. The analysis code is archived at Zenodo: https://doi.org/10.5281/zenodo.20735611 (MIT License). The analysis outputs are archived at Zenodo: https://doi.org/10.5281/zenodo.20735639 (CC0 1.0). A formal manuscript citation will be added to this section after journal acceptance.

## Key References

- Thompson, D'Arcy W. *On Growth and Form* (1917)
- Tarashansky et al. (2021) "Mapping single-cell atlases throughout Metazoa" *eLife*

## Notes for Reviewers and Reproducers

**Reviewer entry-point:** [`CROSSWALK.md`](CROSSWALK.md) maps every Methods subsection,
figure, and numerical claim to its generating code and output file.

**AI use disclosure:** Generative AI tools (Claude, Anthropic) were used to assist with code development, literature search, and manuscript drafting. All AI-generated content was reviewed, verified, and revised by the author. The author takes full responsibility for the accuracy and integrity of the work. See `docs/declarations.txt` for the full disclosure statement and `docs/REPRODUCIBILITY_AUDIT_v2.md` for a reproducibility audit of an earlier state of this repository (historical snapshot, 2026-04-06; the current reproduction path is the four gates described above).

**Tracked files under `data/` and `output/`:** This repository tracks several files inside `data/` and `output/` directories despite those paths appearing in `.gitignore`. These tracked files (centroid CSVs, null distribution arrays, intermediate analysis results, etc.) are reproducibility anchors that travel with the deposit. Verify the tracked set with `git ls-files data/ output/`. If you need to delete and re-add one of these files, use `git add -f` to override the gitignore pattern. **Both commands are repository-only**: a Zenodo archive is not a git repository, so neither runs there. From an archive the equivalent check is that the files are simply present — `find data output -type f | wc -l` counts what the snapshot carries, and every file it lists is by construction one the deposit intended to include, because the archive is built from `git ls-files` at the release commit. There is nothing to reconcile against, and nothing to re-add.

## License

MIT. See [LICENSE](LICENSE).
