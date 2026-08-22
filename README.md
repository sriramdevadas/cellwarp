# CellWarp

Cross-species geometric morphometrics in transcriptomic space, inspired by D'Arcy Thompson's *On Growth and Form*.

## Quick Start

```bash
git clone https://github.com/sriramdevadas/cellwarp.git && cd cellwarp
```

**Verify the deposit in Docker** (no host Python needed): it builds Ubuntu 22.04 + Python 3.12 and runs the four reproduction gates plus the headline fast-path. See [Reproduce in Docker](#reproduce-in-docker).

**Or install locally**, which requires Python 3.12 (the project pins `>=3.12,<3.13`; Ubuntu 24.04 already ships 3.12, Ubuntu 22.04 needs the deadsnakes PPA, and macOS needs `brew install python@3.12` — [Setup](#setup-requires-python-312) has the per-distribution commands):

```bash
python3.12 -m venv .venv
source .venv/bin/activate          # or call .venv/bin/python directly (the gates assume it)
pip install -e ".[dev]"
```

**Which install?** `.[dev]` is enough for the four reproduction gates, the test suite and the
library snippet below; it is what the Docker image installs. Reproducing the paper end to end
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

This project pins **Python `>=3.12,<3.13`**. The stock `python3` on Ubuntu 22.04 (3.10) and macOS (system Python) is too old, and a current Homebrew default (3.14) is too new; the install **will fail** against either. Get 3.12 first if needed, then create the venv with 3.12 and install, kept as one contiguous block so the `pip install` line never runs against the wrong interpreter:

```bash
git clone https://github.com/sriramdevadas/cellwarp.git
cd cellwarp

# Get Python 3.12 and the build prerequisites first:
#   Ubuntu 24.04 (ships Python 3.12.3; no PPA needed):
#     sudo apt-get update
#     sudo apt-get install -y python3.12 python3.12-venv python3.12-dev build-essential
#   Ubuntu 22.04 (stock python3 is 3.10, which fails the pin; add deadsnakes first):
#     sudo add-apt-repository -y ppa:deadsnakes/ppa && sudo apt-get update
#     sudo apt-get install -y python3.12 python3.12-venv python3.12-dev build-essential
#   macOS (the framework Python ships its own headers; no -dev package to install):
#     brew install python@3.12
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
```

`build-essential` and `python3.12-dev` are needed because `[lock]` installs SAMap, whose dependency `hnswlib` is built from source: without the CPython headers the install ends `fatal error: Python.h: No such file or directory` / `ERROR: Could not build wheels for hnswlib`, exit 1, before anything has run. This bites on Linux only. Neither the fast path nor the `.[dev]` gate install needs them.

If you see `No matching distribution found for numpy…`, you're on the wrong Python: that error means an interpreter older than 3.12; use 3.12.

The pin is two-sided, so a **newer** Python fails too, with a different message: `Package 'cellwarp' requires a different Python`. Homebrew now installs 3.14 by default, so on a fresh machine you are likelier to hit the ceiling than the floor. Measured with `pip install --dry-run`: exit 0 on 3.12.13, exit 1 on 3.14.3. There is no 3.13+ fallback, so install 3.12 explicitly (`brew install python@3.12` on macOS, the deadsnakes PPA on Ubuntu) and create the venv with `python3.12` as shown above.

**Requirements:** Python 3.12 (`>=3.12,<3.13`), ~6 GB disk (core), internet for initial data download, and enough memory for the tier you are running:

| what you are running | peak resident | leave yourself |
|---|---|---|
| fast path, or the Docker gate run | negligible | any machine |
| Tier 1 (`[1/8]`–`[8/8]`) | **58.9 GiB** | 64 GB minimum, 128 GB comfortable |
| full pipeline | **58.9 GiB** — the maximum is in Tier 1 | 64 GB minimum, 128 GB comfortable |

The peak is `[4/8]`, `scripts/08_scaled_procrustes.py`, which downloads 992,192 cells in order to keep 140,000. Because that step is **inside Tier 1**, stopping after `TIER 1 COMPLETE` does not avoid it; a 32 GB instance was OOM-killed there at 30.2 GiB. Measured on one platform (AWS `r6i.4xlarge`, 128 GiB, Ubuntu 24.04.4): these are peaks to have headroom above, not thresholds to sit on. See [reproduce/README.md](reproduce/README.md) for the full measurement note.

The pipeline downloads human/mouse atlas data from CELLxGENE Census, runs QC, executes the 35-type Procrustes analysis with permutation testing, and validates all supplementary analyses. After completion, `reproduce/validate.py` checks key statistics against the manuscript values, and the frozen submission text is pinned by `reproduce/MANUSCRIPT_MD5`: the manuscript (`docs/submission/plosone/manuscript_combined.txt`) and both supporting-information texts (`S1_Text.txt`, `S2_Text.txt`), verified with `md5sum -c reproduce/MANUSCRIPT_MD5`, Gate 4.

## Reproduce in Docker

A host-Python-independent **verification** path: it runs the four reproduction gates and the no-download fast-path; it does **not** run the full analysis pipeline (see Scope). The image is Ubuntu 22.04 + Python 3.12, and the gates run at build time, so a green build itself certifies that step:

```bash
git clone https://github.com/sriramdevadas/cellwarp.git
cd cellwarp
docker build -t cellwarp .        # installs Python 3.12 + deps, runs the 4 gates in-build
docker run  --rm cellwarp         # re-runs the 4 gates + the no-download fast-path
```

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

**Tracked files under `data/` and `output/`:** This repository tracks several files inside `data/` and `output/` directories despite those paths appearing in `.gitignore`. These tracked files (centroid CSVs, null distribution arrays, intermediate analysis results, etc.) are reproducibility anchors that travel with the deposit. Verify the tracked set with `git ls-files data/ output/`. If you need to delete and re-add one of these files, use `git add -f` to override the gitignore pattern.

## License

MIT. See [LICENSE](LICENSE).
