#!/usr/bin/env bash
set -euo pipefail

# ──────────────────────────────────────────────────────────────────────
# require_data: gate a pipeline step on the presence of a gitignored
# raw-data file. See DATA_SOURCES.md for which datasets are excluded
# from git and how to fetch them.
#
# Usage: require_data <sentinel_path> <skip_message> && python3 scripts/...
# Returns 0 (proceed) if sentinel exists; returns 1 (skip) otherwise.
# Prints a SKIPPED message on the skip path so the log is self-documenting.
# ──────────────────────────────────────────────────────────────────────
require_data() {
    local sentinel="$1"
    local skip_msg="$2"
    if [ -f "$sentinel" ]; then
        return 0
    else
        echo "  SKIPPED: $skip_msg"
        return 1
    fi
}

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$REPO_ROOT"

echo "=== CellWarp Full Reproduction Pipeline ==="
echo "Start: $(date)"
echo "Python: $(python3 --version)"
echo ""

# ── Preflight checks ──────────────────────────────────────────────

# Verify cellwarp is installed
python3 -c "import cellwarp; print('cellwarp package: OK')"
python3 -c "import scanpy, anndata, scipy, numpy; print('Dependencies: OK')"

# Fail loudly if any hardcoded user paths remain in code
# Matches a literal /Users/ and Path.home(); the latter resolves the repository
# from $HOME, which silently writes outside the tree the deposit was unpacked in.
# docs/ is scanned because a figure producer lives there.
if grep -rqE "/Users/|Path\.home\(\)" --include="*.py" --include="*.R" "$REPO_ROOT/scripts" "$REPO_ROOT/src" \
   "$REPO_ROOT/analysis" "$REPO_ROOT/reproduce" "$REPO_ROOT/docs" 2>/dev/null; then
    echo "ERROR: Hardcoded user paths detected (/Users/ or Path.home()). These must be fixed."
    grep -rnE "/Users/|Path\.home\(\)" --include="*.py" --include="*.R" "$REPO_ROOT/scripts" "$REPO_ROOT/src" \
       "$REPO_ROOT/analysis" "$REPO_ROOT/reproduce" "$REPO_ROOT/docs"
    exit 1
fi

echo "Preflight checks passed."
echo ""

# ── TIER 1: Core pipeline (main result) ──────────────────────────
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  TIER 1: Core Pipeline"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

echo "[1/8] Downloading data from CELLxGENE Census..."
python3 scripts/01_download_data.py

echo "[2/8] QC and normalization..."
python3 scripts/02_qc_and_normalize.py

echo "[3/8] Cell type inventory..."
python3 scripts/08_cell_type_inventory.py

echo "[4/8] Primary 35-type Procrustes analysis (10K permutations)..."
python3 scripts/08_scaled_procrustes.py

echo "[5/8] 1M permutation test..."
python3 scripts/permutation_1M.py

echo "[6/8] GO enrichment..."
python3 scripts/06_go_enrichment.py

echo "[7/8] Bootstrap robustness (100 iterations)..."
python3 scripts/07_bootstrap.py

echo "[8/8] LOOCV..."
python3 scripts/08_loocv.py

echo ""
echo "TIER 1 COMPLETE"
echo ""

# ── TIER 2: Supplementary analyses ───────────────────────────────
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  TIER 2: Supplementary Analyses"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# ── Validation & Controls ─────────────────────────────────────────

echo "[S1] Independent PCA sensitivity (Figure S1A-B)..."
python3 analysis/independent_pca_sensitivity/run_independent_pca.py

echo "[S2] Independent PCA 1M permutation..."
python3 scripts/permutation_1M_independent_pca.py

echo "[S3] Lineage-stratified permutation test (Figure 1C)..."
python3 scripts/test_lineage_stratified_permutation.py

echo "[S4] Human-vs-human negative control (S4 Fig)..."
python3 scripts/test_35type_human_control.py

echo "[S5] Expanded negative controls (Figure S2E)..."
python3 analysis/expanded_negative_controls/expanded_negative_controls.py
python3 analysis/expanded_negative_controls/negative_control_figure.py

# ── Simulation ────────────────────────────────────────────────────

echo "[S6] Simulation study (Figure S1C-F, Table S2 sheet 1)..."
python3 analysis/simulation_study/simulation_study.py
python3 analysis/simulation_study/simulation_figures.py

# ── Parameter & Protocol Sensitivity ──────────────────────────────

echo "[S7] PCA k-sensitivity (Figure S2A-B)..."
python3 scripts/17_pca_sensitivity.py
python3 scripts/18_pca_sensitivity_v2.py

echo "[S8] Smart-seq2 protocol sensitivity (Figure S2C-D)..."
python3 scripts/14_smartseq2_sensitivity.py

# ── Replication ───────────────────────────────────────────────────

echo "[S9] Sun2023 replication (Figure 3A)..."
require_data "data/replication/sun2023/extracted/YC-Liver/matrix.mtx.gz" \
    "Sun2023 raw matrices not present. See DATA_SOURCES.md to fetch from OMIX002605." \
    && python3 scripts/16_sun2023_replication.py

echo "[S10] PanSci replication (Figure 3B)..."
require_data "data/replication/pansci/lung_genecount.mtx.gz" \
    "PanSci raw matrices not present. See DATA_SOURCES.md to fetch from GSE247719." \
    && python3 scripts/pansci_replication.py

echo "[S11] CellHint replication (Figure 3C)..."
python3 scripts/33_cellhint_replication.py

echo "[S12] Cross-atlas ranking replication (Table S1 sheet 2)..."
python3 analysis/ranking_replication/ranking_replication_analysis.py

echo "[S13] Harmonized CellHint replication (Table S4)..."
python3 analysis/harmonized_replication/harmonized_replication.py

# ── Bootstrap Rankings ────────────────────────────────────────────

echo "[S14] Bootstrap ranking analysis (Figure S3C-D, Table S2 sheet 2)..."
python3 analysis/bootstrap_rankings/bootstrap_ranking_analysis.py

# ── Biological & Mechanistic ─────────────────────────────────────

echo "[S15] Treeness-rigidity analysis -- EXPLORATORY; cut from the paper (Table 1 marks T55-T57 intentionally absent), retained for transparency and reproducibility; not a manuscript figure..."
python3 scripts/20_liang_wagner_treeness.py

echo "[S16] CellMarker validation (Figure S6)..."
python3 scripts/cellmarker_35type_rerun.py

echo "[S17] Cell count confound (Figure 6B)..."
python3 scripts/confound_cellcount_rigidity.py

echo "[S18] Biological predictors (Table S1 sheet 1)..."
python3 analysis/biological_predictors/biological_predictors.py

echo "[S19] Two-layer ellipsoid alignment (Figure 4A-D)..."
python3 scripts/t3b_ellipsoid_alignment.py
python3 scripts/layer3_permutation_test.py

echo "[S20] L1000 random baseline (Figure 7A)..."
python3 scripts/35_l1000_random_baseline.py

echo "[S21] Mechanistic null tests (Figure 7B)..."
python3 scripts/12_housekeeping_ratio.py
python3 scripts/13_tf_complexity.py
python3 scripts/12_niche_hypothesis.py
python3 scripts/12_variance_diagnostic.py
python3 scripts/16_interdonor_variance.py
python3 scripts/19_ppi_centrality.py
require_data "data/ucsc/phastCons_placental.bw" \
    "phastCons bigWig not present. See DATA_SOURCES.md to fetch UCSC phastCons tracks for hg38." \
    && python3 scripts/t3e_step2_compute.py
require_data "data/h3k27ac/SENTINEL_FETCHED" \
    "H3K27ac bigWigs not present. See DATA_SOURCES.md to fetch from ENCODE/GEO (script 03b downloads to data/h3k27ac/)." \
    && python3 scripts/t3e_step3b_enhancer.py
python3 scripts/diagnostic_expression_vs_rigidity.py
python3 scripts/16_ribosomal_confound_test.py

# ── Disease Deformation ───────────────────────────────────────────

echo "[S22] Cancer deformation (Figure S4A)..."
python3 scripts/12_cancer_scaled.py

echo "[S23] COVID-19 deformation (Figure S4B)..."
python3 scripts/13_covid_procrustes.py

echo "[S24] Identity vs state analysis (Figure S4C)..."
python3 scripts/identity_vs_state_analysis.py

echo "[S25] DILI deformation (Figure S4D)..."
require_data "data/dilirank/dilirank_v2.xlsx" \
    "DILIrank input dilirank_v2.xlsx not present. See DATA_SOURCES.md to fetch DILIrank v2." \
    && require_data "data/dilirank/lincs_l2_epsilon.gctx" \
        "LINCS L1000 epsilon GCTX not present. See DATA_SOURCES.md." \
    && require_data "data/dilirank/lincs_l2_delta.gctx" \
        "LINCS L1000 delta GCTX not present. See DATA_SOURCES.md." \
    && require_data "/tmp/lincs_sig_info_phase1.txt" \
        "LINCS sig_info Phase I metadata not present at /tmp/lincs_sig_info_phase1.txt. This file is generated outside the pipeline (extracted from LINCS GSE92742 metadata); S25 (DILI) is an optional extended analysis and is skipped without it." \
    && python3 scripts/24_dilirank_analysis.py

# ── CellHint Investigation ────────────────────────────────────────

echo "[S26] CellHint rank reversal investigation (Figure S5, Tables S3-S4)..."
python3 analysis/cellhint_investigation/investigate_rank_reversal.py

# ── SAMap Validation ──────────────────────────────────────────────

echo "[S27] SAMap validation (Figure S6) — requires cellwarp[samap]..."
if python3 -c "import samap" 2>/dev/null; then
    python3 scripts/34_samap_35types.py
else
    echo "  SKIPPED: samap not installed. Install with: pip install cellwarp[samap]"
fi

# ── Tables ────────────────────────────────────────────────────────

echo "[S28] Generate Table S1..."
python3 scripts/create_table_S1.py

echo "[S29] Generate Table S2..."
python3 scripts/create_table_S2.py

# [S29b] table_S2.xlsx is finished by a post-processor, so [S29] alone does not
# reproduce the deposited file. Only S2's post-processor runs here, deliberately.
#
# scripts/46_synthesis_pass_supplementary_table_edits.py also carries editors for
# S1, S3, S4, S5, S6 and key_resources_table.md, and running its main() from this
# pipeline is RULED OUT, not merely avoided:
#   - edit_table_s1 would move table_S1.xlsx, which is under a standing
#     do-not-regenerate ruling: it appends explanatory prose and changes the macaque
#     three-species value from 0.81 to 0.811, neither of which is corroborated by any
#     submitted text. Its cell-count note also cited "Figure 6B" until 184f5ff --
#     this manuscript has five figure captions (Fig 1-5) and zero occurrences of
#     "Fig 6"/"Figure 6" in manuscript_combined.txt, S1_Text.txt or S2_Text.txt, so
#     the reference pointed into the superseded seven-figure version. If that removal
#     is ever reverted, running main() here would inject it into every reproduction.
#   - edit_table_s3 no longer edits anything (its legend was false; see that
#     function's docstring), so calling it would be pointless here.
# Do NOT "complete" this by calling the script whole. The restriction is the point.
echo "[S29b] Finish Table S2 (post-processor)..."
python3 -c "
import importlib.util, pathlib
p = pathlib.Path('scripts/46_synthesis_pass_supplementary_table_edits.py')
spec = importlib.util.spec_from_file_location('s46', p)
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
m.edit_table_s2(m.CANON / 'table_S2.xlsx')
"

# The canonical just moved, and its submission-packet mirror is a frozen byte copy
# that nothing else writes. Refresh it here so Gate 3 and the packet half of Gate 2
# reflect this run. Later stages dirty further pairs; reproduce/README.md covers
# re-running this after the full pipeline.
echo "[S29c] Refresh submission-packet mirrors..."
python3 scripts/build_submission_packet.py --rebuild

echo "[S30] Generate Table S6 (CPC1 driver genes)..."
python3 scripts/generate_table_S6.py

# ── Figure Assembly ───────────────────────────────────────────────

echo "[S31] Generating figure panels..."
python3 scripts/generate_phase1_figures.py
python3 scripts/generate_phase2_figures.py
python3 scripts/generate_phase3_figures.py
python3 scripts/composite_figS3.py

# The deposited figures are assembled outside this pipeline. None of the eight
# producers that write them is invoked here, so this step deliberately runs
# nothing and says so; the previous line claimed the assembly it did not do.
echo "[S32] Deposited figures are NOT rebuilt by this pipeline (see reproduce/figure_script_map.md, Known gaps)."

echo ""
echo "TIER 2 COMPLETE"
echo ""

# ── Validation ───────────────────────────────────────────────────
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Validating Outputs Against Paper"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
python3 reproduce/validate.py

echo ""
echo "=== Pipeline complete. End: $(date) ==="
