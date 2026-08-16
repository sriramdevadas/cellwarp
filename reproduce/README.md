# Reproducing the CellWarp Paper

## Requirements
- Python 3.12 (the project pins `>=3.12,<3.13`)
- ~6 GB disk space for core data (Tier 1, downloaded automatically)
- additional disk space for the optional Tier-2 datasets: see
  [DATA_SOURCES.md](../DATA_SOURCES.md) for per-dataset sizes
- Estimated runtime: **budget a day, not an afternoon.** A measured cold run took
  **~9 h 40 m**: Tier 1 31 min (6.1 GB download included), the rest of Tier 2 87 min,
  out-of-pipeline producers ~70 min — and **`scripts/13_covid_procrustes.py` alone
  6 h 29 m**. That one step issues a CZ CELLxGENE Census query per (cell type × tissue)
  combination across 20 cell types with `CENSUS_QUERY_TIMEOUT = 600`; in the measured run
  14 of those timeouts fired, which is over two hours of waiting on its own. It prints
  nothing while a query is outstanding. **It has not hung.** Timings scale with Census
  responsiveness, so treat them as an order of magnitude, not a promise.
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
pip install -e ".[lock,dev]"   # [dev] is required: Gate 2 is pytest, which [lock] does not install

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

**What Gate 1 does and does not establish.** Every artifact it reads is a tracked file in
this repository, so it returns a full pass on a freshly unpacked archive with no pipeline
step run at all. It checks that the values reported in the manuscript match the artifacts
deposited beside them — a consistency check, not evidence that the code regenerates those
artifacts. To use it as a reproduction test, run the pipeline first and then re-run it, so
the artifacts it reads are ones this machine produced.

## After a full run: rebuild the packet before the gates mean anything

`run_all.sh` rewrites tracked canonical artifacts. Six of them have byte-copies in the
submission packet (`docs/submission/figures_for_review/` and one panel mirror), and nothing
but `scripts/build_submission_packet.py --rebuild` writes those copies. So a **successful**
reproduction leaves the canonical and its mirror out of sync, and both Gate 3 and the
`tests/test_submission_packet_consistency.py` half of Gate 2 fail — because the run worked,
not because it did not.

```bash
python scripts/build_submission_packet.py --rebuild   # after run_all.sh, before the gates
```

Expect exactly these five pairs to drift on a fresh run: `figS1_pipeline_validation.pdf`,
`figS2_parameter_protocol_sensitivity.pdf`, `figS3_bootstrap_rankings.pdf`,
`table_S1.xlsx`, and `covid_cross_analysis.png`. On the pristine archive all 30 pairs
match, which is why the gates are green before you reproduce anything.

`table_S2.xlsx` was the sixth until `scripts/create_table_S2.py` was given the fixed-epoch
stamp described below; it now regenerates byte-for-byte and its pair no longer drifts.
`[S29c]` also refreshes the mirrors mid-run, so the pairs listed above drift only if a
later stage rewrites them.

### Why a spreadsheet pair drifts: the writer's clock, not the content

`table_S1.xlsx` drifts for a different reason from the four figure pairs.
`openpyxl` stamps wall-clock time into every workbook it writes — `dcterms:created` and
`dcterms:modified` in `docProps/core.xml`, and an mtime on every zip entry — so
regenerating one yields a new md5 **even when every cell is identical**. An md5 over such a
file pins the moment it was written rather than its content, and no content fix removes the
drift.

That is a property of the writer, not of the format. `scripts/table1_formatting.py`
normalises it: it writes `dcterms:modified = 2026-01-01T00:00:00Z` and a fixed date on
every zip entry, and is byte-idempotent as a result — consecutive runs produce identical
bytes. `scripts/create_table_S2.py` now does the same, via a
`normalize_xlsx_timestamps()` that `edit_table_s2()` imports rather than reimplements, so
whichever of the two writes last applies the stamp. Compare a workbook from a writer that
does this by md5; compare one that does not, such as `table_S1.xlsx`, cell by cell.

## The supplementary tables `run_all.sh` does not finish

`run_all.sh` generates three supplementary tables — `[S28]` `table_S1.xlsx`, `[S29]`
`table_S2.xlsx`, `[S30]` `Table_S6_CPC1_driver_genes.xlsx` — and finishes none of them.

**A required hand step is missing from the pipeline.**
`scripts/46_synthesis_pass_supplementary_table_edits.py` post-processes six deposited
tables — `table_S1.xlsx`, `table_S2.xlsx`, `table_S3.csv`, `table_S4.csv`, `table_S5.csv`
and `Table_S6_CPC1_driver_genes.xlsx` — plus `docs/submission/key_resources_table.md`. It
is referenced nowhere in `run_all.sh`, and its own closing line asks for a packet refresh
that `run_all.sh` also does not perform. So `[S28]` and `[S29]` write the create-stage
output only, and a reader's `table_S1.xlsx` and `table_S2.xlsx` will differ from the
deposited copies. `table_S3.csv`, `table_S4.csv` and `table_S5.csv` have no generator at
all — `run_all.sh` never writes them, so those three remain as deposited.

**`[S30]` cannot run from the deposit.** `scripts/generate_table_S6.py` opens
`data/phase2_scaled/human_scaled.h5ad`. `.gitignore` excludes `data/` wholesale, no file
under `data/phase2_scaled/` is tracked, and the archive does not carry it. Under
`set -euo pipefail` (line 2) the run stops there, before `reproduce/validate.py` — so fetch
the Tier-2 inputs first, or skip the stage. `Table_S6_CPC1_driver_genes.xlsx` is therefore
**untested against the deposit**, which is not the same as passing.

## Analyses `run_all.sh` does not run

The manuscript states that the conserved-contribution analyses in Results section 5 "are
run outside the automated reproduction script and must be invoked directly". They are, in
this order — each step reads the previous one's output:

```bash
python analysis/conserved_contribution/run_gate.py            # gate_results.json, gene_conservation_core.csv
python analysis/conserved_contribution/run_robustness.py      # robustness_results.json
python analysis/conserved_contribution/highN_tf_pvalues.py    # highN_tf_pvalues.json
python analysis/conserved_contribution/breadth_sensitivity.py # breadth_sensitivity_results.json
python analysis/conserved_contribution/make_table_s11.py      # S11 Table
python analysis/conserved_contribution/make_figure7.py        # Fig 5
```

The last two need `donor_stability/agg_{human,mouse}_cap10000.npz`, which are gitignored
and not redistributed, and they stop cleanly when the files are absent. `run_all.sh`
rebuilds no deposited figure at all; `reproduce/figure_script_map.md` (Known gaps) lists
the eight producers that write one and the order the S2 Fig chain requires.

## Which manifest is authoritative

**`pyproject.toml`'s `[lock]` extra is the instruction; `requirements.txt` is a record.**
`[lock]` pins the 27 direct dependencies and lets pip resolve the rest. `requirements.txt`
is a 196-package freeze taken from one machine on 2026-05-18, and installing `.[lock]`
today resolves 59 of those packages to different versions and omits 6 — including
`leidenalg`, which "Seed and determinism" below names as one of the two
environment-sensitive steps. Reproduce from `.[lock]`; read `requirements.txt` to see what
the authoring environment happened to contain, and do not `pip install -r` it expecting the
two to agree.

## Figure bytes record the matplotlib version

Figure metadata records the matplotlib version. PNGs carry it in a `Software` chunk and
PDFs in `/Producer`, so rebuilding a figure under a different matplotlib gives a file whose
bytes differ while every pixel is identical.

Fourteen deposited figures were built with matplotlib 3.10.9. The manifests pin 3.10.8, and
the other 190 deposited PNGs were built with the pinned version. Rebuilding any of the
fourteen under the pinned environment reproduces them pixel for pixel — zero differing
pixels, maximum delta 0.0000 — and changes only that string. **For these files, compare
pixels rather than md5.**

The fourteen are: `Fig1_configuration_conserved.png` through
`Fig5_conserved_identity_genes.png` and `Fig2C_bg_replication.png` under
`docs/submission/plosone/figures/`; `figures/main/fig7_conserved_contribution.png`;
`docs/supplementary_materials/figure_S8_markernull.png`;
`figures/supplementary/negative_control_distributions.png`; and
`fig1a_pipeline_schematic.png`, `fig1b_null_1M.png`, `fig1c_lineage_stratified.png`,
`fig3b_pre_post.png` and `fig4d_replication_summary.png` under `figures/panels/`. Fig 5 is
a byte copy of `fig7_conserved_contribution.png` and carries that file's chunk unchanged,
which is why it alone matches a 3.10.9 rebuild.

The five submission TIFFs are unaffected: the TIFF writer records no version string, and
all five regenerate byte-identically under the pinned environment.

## Two interpreters: the gates and the DOCX build

The four gates and the DOCX build need not run under the same interpreter, and
the final submission sequence touches both.

| what | interpreter | needs |
|---|---|---|
| Gate 1 `reproduce/validate.py` | `.venv` | core deps |
| Gate 2 `pytest -q` | `.venv` | `[dev]` |
| Gate 3 `scripts/build_submission_packet.py --verify` | `.venv` | core deps |
| Gate 4 `md5sum -c reproduce/MANUSCRIPT_MD5` | no interpreter | — |
| `docs/submission/plosone/build_manuscript_docx.py` | **whichever one can `import docx`** | `python-docx` |

`build_manuscript_docx.py` imports `python-docx`. The rule is a condition, not a
named interpreter: **check that `import docx` succeeds in the interpreter you
are about to use**, and run the builder there. Stated that way it stays true
however the environments on a particular machine drift. On the reference machine
the miniforge base python satisfies it; a given `.venv` may or may not, according
to when and from what it was built. The builder fails fast and legibly when it is
missing:

```
ERROR: python-docx is required (pip install 'python-docx>=1.1'): No module named 'docx'
```

Where it is declared: commit `3f1d326` added `python-docx` to `pyproject.toml`'s
`[dev]` and `[lock]` extras and to `requirements.txt` and `environment.yml`,
precisely because no documented install read the `[reproduce]` extra it had been
declared in alone. So every install path the documentation names now declares
it. Text here previously said the opposite; it was left behind by that commit
and is corrected rather than restated, since a claim about what some environment
contains goes stale and the import condition above does not.

Gate 3 stays pointed at `scripts/build_submission_packet.py --verify` rather
than at `docs/submission/plosone/figures/build_submission_tiffs.py`. That was
deferred on cost, not on whether `pymupdf` can be imported, so finding `pymupdf`
present in a given interpreter does not reopen it. Noted separately because it
is a real gap and not the reason: `[dev]`, the only extra the Dockerfile
installs, declares `python-docx` but not `pymupdf`.

Why this is easy to miss: the DOCX is gitignored and no gate reads it, so a
**broken DOCX build leaves all four gates green**. Green gates are not evidence
the submission document builds. Run the builder explicitly and check its exit
code before submitting.

## Figure-to-script mapping

See `reproduce/figure_script_map.md` for a complete table showing which
script generates each figure and table in the paper.

## SAMap validation

**`.[lock]` already installs SAMap, so the SAMap step runs by default.** `[lock]` pins
`samap==1.0.14`, `import samap` succeeds, and `run_all.sh`'s probe therefore takes the
SAMap branch: step S27 executes (~4 min in the measured run). Earlier text here said SAMap
was excluded from the default install; that was true of `pip install -e .` (the fast-path
base install) and never of `.[lock]`, which is what the Steps section prescribes.

`torch` is **not** installed by `.[lock]`, despite `environment_ground_truth.txt` recording
`torch==2.10.0`; `samap` imports and the step runs without it. To add SAMap to an
environment built some other way:

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
| `EXPECTED_JOINED_WORDS` | 14952 | `docs/submission/plosone/build_manuscript_docx.py` | **yes** — the builder aborts on mismatch |
| `wc -w` on `manuscript_combined.txt` | 15117 | shell | no |
| `wc -w` on `S1_Text.txt` / `S2_Text.txt` | 5979 / 975 | shell | no |

`EXPECTED_JOINED_WORDS` is smaller because it counts the **161 content lines
joined**, i.e. the manuscript body as the DOCX renders it: the section banners,
the `====` rules and the blank structural lines that `wc` sees are not content
and are not counted. It is not a "wrong" `wc`; it is a different object. Do not
reconcile the two, and do not copy a `wc` number into the builder.

### `wc -w` is GNU-vs-BSD dependent, and the rule is exact

**GNU `wc -w` under `LC_ALL=C` counts zero words for a token made entirely of
non-ASCII characters.** In the byte locale those bytes are not printable, and
GNU `wc` excludes non-printable characters when deciding what a word is, so a
standalone `ρ` or `≈` contributes nothing. Under a UTF-8 locale the same bytes
decode to one printable wide character and the token counts as one word. BSD
`wc` does not apply the printability filter at all, which is why it returns the
same number under every locale.

The rule is exact, not approximate: the two counts differ by precisely the
number of whitespace-delimited tokens composed **only** of non-ASCII
characters. For the post-edit manuscript that is 77 tokens — the spaced
operators, 31 `ρ`, 13 `×`, 9 `≈`, 8 `→`, 5 `∈`, 5 `≤`, 4 `≥`, 1 `α`, 1 `—`.

| count | value | how |
|---|---|---|
| BSD `wc -w`, any locale; GNU under UTF-8; Python `str.split()` | **15117** | every token counts |
| GNU `wc -w` under `LC_ALL=C` | **15040** | 15117 − 77 all-non-ASCII tokens |
| `EXPECTED_JOINED_WORDS` | **14952** | the gate; see above |

A token that *mixes* ASCII and non-ASCII — `human–mouse`, `50–2,000`, `ρ = 0.45`
once the `=` is its own token — still contains printable ASCII and still counts
as one under both. That is the whole reason en dashes are irrelevant here: all
89 of them sit inside mixed tokens and none is ever a bare token.

**The container runs a GNU userland.** A word count taken inside it with `LANG`
unset gets the *lower* number, 15040, on a correct tree. That is not drift and
not a corrupted file; it is this rule. Neither 15117 nor 15040 is
`EXPECTED_JOINED_WORDS`, and neither is asserted by any gate.

Which of the two a GNU `wc` returns depends on how it was launched, not only
on the locale variables: CPython coerces the legacy C locale and exports
`LC_CTYPE=C.UTF-8` to child processes, so `wc -w` invoked from a shell with
`LANG` unset returns 15040 while the same `wc -w` invoked through Python
returns 15117. Anything run through the pipeline is Python-launched.

GNU `wc` is not installed on the macOS reference machine, so the 15040 above is
the arithmetic prediction of the rule, and it matches the GNU measurement taken
during the Tier-1 pass exactly.

#### A coincidence, recorded so nobody chases it twice

At the pre-edit HEAD the two manuscript counts were 13606 and 13524, a gap of
**82** — and the manuscript contained exactly **82** en dashes. Two people
independently chased the en dash as the cause. It is not. The gap was 82
because the pre-edit manuscript happened to contain 82 all-non-ASCII tokens as
well; the edits moved that to 77 while leaving the en dashes at 82 as of that
HEAD (they are 89 at `453c804`), which separates the two numbers. The other two
texts refute the en-dash reading outright, at that same HEAD:

| file (at the pre-edit HEAD) | all-non-ASCII tokens = the gap | en dashes |
|---|---|---|
| `manuscript_combined.txt` | 82 | 82 ← the coincidence |
| `S1_Text.txt` | 39 | 18 |
| `S2_Text.txt` | 5 | 0 |

Practical consequence: if a submission form wants a word count, say which tool
and which locale produced it.
