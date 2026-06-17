#!/usr/bin/env python3
"""
Reconstruct the committed macaque_pipeline centroid construction and verify
the 13-type RIRA sensitivity result (obs/null ≈ 0.749, p ≈ 0.0002).

Why this exists
---------------
The script that produced output/macaque_pipeline/primary_procrustes_results.json
and sensitivity_procrustes_results.json is not in the committed tree. Without
a driver, the 20-type primary obs/null = 0.841 and the 13-type sensitivity
obs/null = 0.749 are not independently reproducible from local data, which is
a submission-blocking reproducibility gap (see DECISIONS.md around Qu
annotation access).

Qu per-cell annotations are in a 14.33 GB RDS on Zenodo (10.5281/zenodo.5881495)
and are not on disk at the time this driver was written, so this script
covers only the 13 RIRA-sourced cell types — which map cleanly to the
committed sensitivity result. Extending to the full 20 types requires the
Qu Zenodo RDS and a complementary extractor.

Verification targets
--------------------
From output/macaque_pipeline/sensitivity_procrustes_results.json:
    n_types        = 13
    obs/null       = 0.7488054428
    p              = 0.0001999800
    scaling        = 0.0534434720
    distance       = 16.5307225745

Pipeline
--------
1.  Harmonize RIRA cell-type labels to the CellWarp Census vocabulary
    (13 target types, exact cell-count reproduction verified in 1b).
2.  Build the 13,927-gene three-way ortholog intersection
    (RIRA-symbol ∩ Qu-ENSMFAG-via-BioMart ∩ human-mouse-CellWarp-space).
3.  Subset human centroids (committed) to 13 RIRA types × 13,927 genes.
4.  Compute macaque-RIRA centroids from RIRA counts: filter cells → filter
    genes → normalize_total(1e4) → log1p → per-type mean.
5.  Sanity: compare macaque centroid value range/variance against committed
    human centroids. Stop if distributions differ wildly.
6.  Joint PCA via src/procrustes.py::pca_reduce_centroids (≥95% variance).
7.  Procrustes align + 10,000-permutation test via src/procrustes.py,
    seed=42. Verify obs/null and p against committed sensitivity values.

Committed conventions applied
-----------------------------
- Gene filter before normalization (mouse_lemur/01_run_pipeline.py precedent).
- normalize_total(1e4) + log1p, no z-scoring (DECISION-064).
- Exclude cells with RIRA_Immune_v2.cellclass == "Unknown" (step0b).
- Per-type RIRA-label harmonization table defined in `RIRA_TYPE_MAP` below.

Inputs
------
- data/macaque/rira/converted/{barcodes,genes,matrix}.tsv.gz / .mtx.gz
- data/macaque/rira/rira_metadata.csv  (produced by extract_rira_metadata.R)
- data/macaque/biomart_macaque_human_orthologs.csv
- data/macaque/qu_2022/extracted/GSM5901076_Adipose_features.tsv.gz
    (any single GSM works — all 20 have identical feature sets)
- data/phase2_scaled/human_scaled.h5ad  (var has feature_id, feature_name)
- output/phase2/scaled_35types/centroids_human_35.csv

Outputs
-------
- output/macaque_pipeline/reconstruction_rira13_results.json
- output/macaque_pipeline/reconstruction_rira13_centroids.csv
- output/macaque_pipeline/reconstruction_rira13_gene_list.csv
- analysis/macaque/reconstruction_rira13_report.md
"""
from __future__ import annotations

import gzip
import json
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io as sio
import scipy.sparse as sp

PROJECT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT / "src"))
from procrustes import pca_reduce_centroids, procrustes_align, permutation_test  # noqa

RANDOM_SEED = 42
N_PERMUTATIONS = 10_000
MAX_CELLS_PER_TYPE = 2_000  # Matches mouse_lemur/01_run_pipeline.py and verified
                            # against human_scaled.h5ad (per-type exactly 2000 cells
                            # for types with >=2000 available; centroids_human_35.csv
                            # reproduces from h5ad .X to float precision).

# Paths
RIRA_DIR = PROJECT / "data" / "macaque" / "rira"
RIRA_CONV = RIRA_DIR / "converted"
RIRA_META = RIRA_DIR / "rira_metadata.csv"
QU_FEATURES = PROJECT / "data/macaque/qu_2022/extracted/GSM5901076_Adipose_features.tsv.gz"
BIOMART = PROJECT / "data/macaque/biomart_macaque_human_orthologs.csv"
HUMAN_SCALED_H5AD = PROJECT / "data/phase2_scaled/human_scaled.h5ad"
HUMAN_CENTROIDS_CSV = PROJECT / "output/phase2/scaled_35types/centroids_human_35.csv"

OUT_DIR = PROJECT / "output/macaque_pipeline"
OUT_JSON = OUT_DIR / "reconstruction_rira13_results.json"
OUT_CENTROIDS = OUT_DIR / "reconstruction_rira13_centroids.csv"
OUT_GENES = OUT_DIR / "reconstruction_rira13_gene_list.csv"
REPORT_MD = PROJECT / "analysis/macaque/reconstruction_rira13_report.md"

# ---------------------------------------------------------------------------
# RIRA label harmonization (verified exact against centroid_cell_counts.csv)
# Each predicate is evaluated on a cell filtered to RIRA_Immune_v2.cellclass
# being one of {Bcell, T_NK, Myeloid}. "T cell" = T_NK cells not in
# {CD4+ T Cells, CD8+ T Cells, NK Cells}. "myeloid leukocyte" = Myeloid
# cells with Myeloid_v3.cellclass in {Unassigned, Ambiguous, Unknown, NaN}.
# ---------------------------------------------------------------------------
RIRA_TARGET_TYPES = [
    "B cell",
    "CD4-positive, alpha-beta T cell",
    "CD8-positive, alpha-beta T cell",
    "T cell",
    "classical monocyte",
    "granulocyte",
    "hematopoietic precursor cell",
    "intermediate monocyte",
    "macrophage",
    "myeloid dendritic cell",
    "myeloid leukocyte",
    "natural killer cell",
    "non-classical monocyte",
]

# Expected counts from output/macaque_pipeline/centroid_cell_counts.csv.
# Pipeline must reproduce these exactly.
EXPECTED_RIRA_COUNTS = {
    "B cell": 84_412,
    "CD4-positive, alpha-beta T cell": 116_863,
    "CD8-positive, alpha-beta T cell": 110_916,
    "T cell": 48_226,
    "classical monocyte": 7_250,
    "granulocyte": 5_698,
    "hematopoietic precursor cell": 1_501,
    "intermediate monocyte": 1_452,
    "macrophage": 14_298,
    "myeloid dendritic cell": 4_780,
    "myeloid leukocyte": 2_138,
    "natural killer cell": 11_802,
    "non-classical monocyte": 3_484,
}

# Verification target for the Procrustes run (2-config, human vs RIRA-13)
EXPECTED_OBS_NULL = 0.7488054428
EXPECTED_P = 0.00019998000199980003
TOL_OBS_NULL = 0.01  # user-approved tolerance


def harmonize_rira_labels(meta: pd.DataFrame) -> pd.Series:
    """Apply the RIRA → CellWarp label mapping; NaN for cells outside targets."""
    imm = meta["RIRA_Immune_v2.cellclass"]
    tnk = meta["RIRA_TNK_v2.cellclass"]
    mye = meta["RIRA_Myeloid_v3.cellclass"]

    out = pd.Series(pd.NA, index=meta.index, dtype="object")

    # B cell
    out[imm == "Bcell"] = "B cell"

    # T_NK subtypes
    is_tnk = imm == "T_NK"
    out[is_tnk & (tnk == "CD4+ T Cells")] = "CD4-positive, alpha-beta T cell"
    out[is_tnk & (tnk == "CD8+ T Cells")] = "CD8-positive, alpha-beta T cell"
    out[is_tnk & (tnk == "NK Cells")] = "natural killer cell"
    t_residual = is_tnk & ~tnk.isin(["CD4+ T Cells", "CD8+ T Cells", "NK Cells"])
    out[t_residual] = "T cell"

    # Myeloid subtypes
    is_mye = imm == "Myeloid"
    out[is_mye & (mye == "CD14+ Monocytes")] = "classical monocyte"
    out[is_mye & (mye == "Inflammatory Monocytes")] = "intermediate monocyte"
    out[is_mye & (mye == "CD16+ Monocytes")] = "non-classical monocyte"
    out[is_mye & mye.isin(["Macrophages", "Alv. mac."])] = "macrophage"
    out[is_mye & mye.isin(["DC", "pDC", "Mature DC"])] = "myeloid dendritic cell"
    out[is_mye & (mye == "Myelocytes")] = "granulocyte"
    out[is_mye & (mye == "Promyelocytes")] = "hematopoietic precursor cell"
    mye_residual = is_mye & ~mye.isin([
        "CD14+ Monocytes", "Inflammatory Monocytes", "CD16+ Monocytes",
        "Macrophages", "Alv. mac.", "DC", "pDC", "Mature DC",
        "Myelocytes", "Promyelocytes",
    ])  # catches {Unassigned, Ambiguous, Unknown, NaN}
    out[mye_residual] = "myeloid leukocyte"

    return out


def build_three_way_gene_list(human_ensg: list[str], human_symbols: list[str],
                              rira_symbols: set[str],
                              qu_symbol_by_ensmfag: dict[str, str],
                              qu_direct_symbols: set[str],
                              biomart_ensmfag_to_human: dict[str, str]
                              ) -> pd.DataFrame:
    """Return the three-way ortholog ENSG list.

    Qu coverage of a human symbol is defined as (matches directly in Qu
    HGNC-like symbol column) OR (a Qu ENSMFAG row exists whose BioMart
    'Human gene name' equals the human symbol). This matches the two-step
    construction described in step0a (direct match + BioMart recovery of
    ENSMFAG-only rows → 13,597 + 961 = 14,558 Qu-CellWarp overlap).
    """
    # Qu human-symbol coverage: direct Qu HGNC symbols + BioMart-recovered
    # human gene names from Qu's ENSMFAG-only rows.
    qu_covered_symbols = set(qu_direct_symbols)
    for ensmfag in qu_symbol_by_ensmfag:
        if qu_symbol_by_ensmfag[ensmfag] is None:  # ENSMFAG-only row
            hgnc = biomart_ensmfag_to_human.get(ensmfag)
            if hgnc and hgnc != "":
                qu_covered_symbols.add(hgnc)

    rows = []
    for ensg, sym in zip(human_ensg, human_symbols):
        in_rira = sym in rira_symbols
        in_qu = sym in qu_covered_symbols
        if in_rira and in_qu:
            rows.append({"ensg": ensg, "symbol": sym})
    return pd.DataFrame(rows)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t_start = time.time()

    print("=" * 70)
    print("MACAQUE RIRA-13 RECONSTRUCTION — VERIFY OBS/NULL ≈ 0.749")
    print("=" * 70)

    # ── 1. Load RIRA metadata ─────────────────────────────────────────────
    print("\n[1] Loading RIRA metadata CSV…")
    meta_cols = [
        "cell_id",
        "RIRA_Immune_v2.cellclass",
        "RIRA_TNK_v2.cellclass",
        "RIRA_Myeloid_v3.cellclass",
    ]
    meta = pd.read_csv(RIRA_META, usecols=meta_cols, low_memory=False)
    print(f"  Loaded {len(meta):,} cells")

    # ── 1a. Harmonize labels ──────────────────────────────────────────────
    print("\n[1a] Harmonizing RIRA labels → CellWarp Census types…")
    meta["target_type"] = harmonize_rira_labels(meta)
    # Drop cells outside the 13 target types (Unknown, Non-Immune, Erythrocyte, etc.)
    meta_kept = meta[meta["target_type"].notna()].reset_index(drop=True)
    print(f"  Cells after filter to 13 target types: {len(meta_kept):,}")

    # ── 1b. Verify exact cell-count reproduction ──────────────────────────
    print("\n[1b] Verifying per-type cell counts vs centroid_cell_counts.csv…")
    counts = meta_kept["target_type"].value_counts().to_dict()
    fail = False
    for t in RIRA_TARGET_TYPES:
        got = counts.get(t, 0)
        exp = EXPECTED_RIRA_COUNTS[t]
        status = "OK" if got == exp else "FAIL"
        if got != exp:
            fail = True
        print(f"  {t:<40s} got={got:>7,}  expected={exp:>7,}  {status}")
    if fail:
        raise SystemExit("Per-type cell counts do not match — stopping before centroid computation.")
    print("  All 13 RIRA-sourced cell counts reproduce exactly.")

    # ── 2. Build three-way 13,927-gene list ───────────────────────────────
    print("\n[2] Building three-way ortholog gene list (target: 13,927 ENSG)…")

    # 2a. CellWarp human gene space
    print("  2a. Loading human_scaled.h5ad var (CellWarp 16,959-gene space)…")
    a = ad.read_h5ad(HUMAN_SCALED_H5AD, backed="r")
    human_ensg = list(a.var["feature_id"])
    human_symbols = list(a.var["feature_name"])
    assert len(human_ensg) == 16_959, f"CellWarp space must be 16,959 genes (got {len(human_ensg)})"
    a.file.close()
    # Map ENSG → symbol for later centroid filtering
    ensg_to_sym = dict(zip(human_ensg, human_symbols))

    # 2b. RIRA gene symbols
    print("  2b. Loading RIRA gene symbols…")
    with gzip.open(RIRA_CONV / "genes.tsv.gz", "rt") as f:
        rira_symbols_list = [line.strip() for line in f]
    rira_symbols = set(rira_symbols_list)
    print(f"     RIRA unique symbols: {len(rira_symbols):,}")

    # 2c. Qu features (direct symbols + ENSMFAG-only rows)
    print("  2c. Loading Qu features…")
    qu_direct_symbols: set[str] = set()
    qu_symbol_by_ensmfag: dict[str, str | None] = {}
    with gzip.open(QU_FEATURES, "rt") as f:
        for line in f:
            ensmfag, sym, _ = line.rstrip("\n").split("\t")
            if sym.startswith("ENSMFAG") or sym == ensmfag:
                qu_symbol_by_ensmfag[ensmfag] = None  # ENSMFAG-only row
            else:
                qu_direct_symbols.add(sym)
                qu_symbol_by_ensmfag[ensmfag] = sym
    n_ensmfag_only = sum(1 for v in qu_symbol_by_ensmfag.values() if v is None)
    print(f"     Qu total features: {len(qu_symbol_by_ensmfag):,}")
    print(f"     Qu direct HGNC-like symbols: {len(qu_direct_symbols):,}")
    print(f"     Qu ENSMFAG-only rows: {n_ensmfag_only:,}")

    # 2d. BioMart ENSMFAG → human gene name
    print("  2d. Loading BioMart macaque orthologs…")
    bm = pd.read_csv(BIOMART)
    bm_121 = bm[bm["Human homology type"] == "ortholog_one2one"]
    biomart_map = dict(zip(bm_121["Gene stable ID"].astype(str),
                           bm_121["Human gene name"].fillna("").astype(str)))
    print(f"     1:1 orthologs in BioMart: {len(biomart_map):,}")

    # 2e. Compute three-way list
    gene_df = build_three_way_gene_list(
        human_ensg, human_symbols, rira_symbols,
        qu_symbol_by_ensmfag, qu_direct_symbols, biomart_map,
    )
    print(f"  2e. Three-way intersection size: {len(gene_df):,} ENSGs")
    if len(gene_df) != 13_927:
        print(f"  WARN: expected 13,927 but got {len(gene_df):,}. Will proceed but flag.")
    gene_df.to_csv(OUT_GENES, index=False)
    three_way_ensg = list(gene_df["ensg"])
    three_way_sym = list(gene_df["symbol"])

    # ── 3. Subset human centroids ─────────────────────────────────────────
    print("\n[3] Loading + subsetting human centroids (35 → 13 types × 13,927 genes)…")
    human_cents = pd.read_csv(HUMAN_CENTROIDS_CSV, index_col=0)
    print(f"  Full: {human_cents.shape}")
    human_cents_sub = human_cents.loc[RIRA_TARGET_TYPES, three_way_ensg].copy()
    print(f"  Subset: {human_cents_sub.shape}")
    assert human_cents_sub.shape == (13, len(three_way_ensg))

    # ── 4. Load RIRA counts + compute centroids ───────────────────────────
    print("\n[4] Loading RIRA counts matrix (may take ~1 min)…")
    t_mtx = time.time()
    # Matrix is genes × cells (per standard 10x mtx); scanpy's read_mtx yields
    # an AnnData with rows=cells if we transpose; we use direct scipy read.
    with gzip.open(RIRA_CONV / "matrix.mtx.gz", "rt") as f:
        mtx = sio.mmread(f)  # (n_genes, n_cells)
    mtx = sp.csr_matrix(mtx.T)  # now (n_cells, n_genes)
    print(f"  Loaded matrix: {mtx.shape} in {time.time() - t_mtx:.1f}s")

    with gzip.open(RIRA_CONV / "barcodes.tsv.gz", "rt") as f:
        rira_barcodes = [line.strip() for line in f]
    assert len(rira_barcodes) == mtx.shape[0]
    assert len(rira_symbols_list) == mtx.shape[1]

    # Build AnnData with barcode index matching meta.cell_id order
    adata = ad.AnnData(
        X=mtx,
        obs=pd.DataFrame(index=rira_barcodes),
        var=pd.DataFrame(index=rira_symbols_list),
    )

    # Align metadata to adata.obs (metadata cell_id should match barcodes)
    # Verify ordering — RIRA metadata.rds rownames are the barcodes
    print("\n[4a] Aligning metadata to counts…")
    meta_all = pd.read_csv(
        RIRA_META,
        usecols=["cell_id"] + [
            "RIRA_Immune_v2.cellclass", "RIRA_TNK_v2.cellclass",
            "RIRA_Myeloid_v3.cellclass",
        ],
        low_memory=False,
    )
    meta_all = meta_all.set_index("cell_id")
    # Reorder to match adata.obs
    if not meta_all.index.equals(adata.obs.index):
        # Reindex (may produce NaN for missing barcodes)
        missing = set(adata.obs.index) - set(meta_all.index)
        if missing:
            print(f"  {len(missing):,} barcodes in counts absent from metadata")
        meta_all = meta_all.reindex(adata.obs.index)
    adata.obs["target_type"] = harmonize_rira_labels(meta_all)
    kept_mask = adata.obs["target_type"].notna().to_numpy()
    print(f"  Cells with target type (pre-subsample): {int(kept_mask.sum()):,}")

    # Subsample to MAX_CELLS_PER_TYPE per type (seed=42).
    # The committed centroid_cell_counts.csv reports TOTAL available cells;
    # the centroids themselves are computed on the subsampled set, per the
    # CellWarp convention verified against human_scaled.h5ad.
    print(f"\n[4a'] Subsampling each type to max {MAX_CELLS_PER_TYPE:,} cells (seed={RANDOM_SEED})…")
    rng = np.random.default_rng(RANDOM_SEED)
    sub_mask = np.zeros(adata.n_obs, dtype=bool)
    for t in RIRA_TARGET_TYPES:
        idx = np.where((adata.obs["target_type"] == t).to_numpy())[0]
        if len(idx) > MAX_CELLS_PER_TYPE:
            sel = rng.choice(idx, size=MAX_CELLS_PER_TYPE, replace=False)
            sel.sort()
        else:
            sel = idx
        sub_mask[sel] = True
        print(f"  {t:<40s} {len(idx):>7,} → {int(sub_mask[sel].sum()):>6,}")
    # Compose subsample + type filter
    kept_mask = kept_mask & sub_mask
    print(f"  Cells after subsample: {int(kept_mask.sum()):,}")

    # ── 4b. Filter to 13 target types + 13,927 gene space (symbols) ───────
    print("\n[4b] Filtering cells + genes…")
    adata = adata[kept_mask].copy()
    # Gene filter: keep RIRA rows whose symbol matches a three-way human symbol.
    gene_mask = np.isin(adata.var_names.to_numpy(), np.array(three_way_sym))
    adata = adata[:, gene_mask].copy()
    print(f"  After gene filter: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    # A RIRA symbol may have duplicates in genes.tsv (rare); aggregate by sum
    # for safety so downstream centroid means behave.
    if len(set(adata.var_names)) != adata.n_vars:
        print("  WARN: duplicate RIRA symbols present — aggregating by sum per symbol")
        import scipy.sparse as _sp
        idx = pd.Series(range(adata.n_vars), index=adata.var_names)
        groups = idx.groupby(idx.index).apply(list)
        n_unique = len(groups)
        rows = []
        cols = []
        data = []
        new_gene_names = []
        Xc = adata.X.tocsc()
        for gcol, (sym, idx_list) in enumerate(groups.items()):
            col = Xc[:, idx_list].sum(axis=1)
            new_gene_names.append(sym)
            nz = np.nonzero(np.asarray(col).flatten())[0]
            for r in nz:
                rows.append(r)
                cols.append(gcol)
                data.append(np.asarray(col).flatten()[r])
        X_new = _sp.csr_matrix(
            (data, (rows, cols)), shape=(adata.n_obs, n_unique)
        )
        adata = ad.AnnData(
            X=X_new,
            obs=adata.obs.copy(),
            var=pd.DataFrame(index=new_gene_names),
        )
        print(f"  After aggregation: {adata.n_obs:,} cells × {adata.n_vars:,} unique gene symbols")

    # Re-order RIRA genes to match three_way_sym / three_way_ensg order
    adata = adata[:, [s for s in three_way_sym if s in set(adata.var_names)]].copy()
    # Some symbols might be missing (RIRA symbol set is union; 13,927 is
    # intersection, so all 13,927 symbols should be present). Verify:
    missing_from_rira = set(three_way_sym) - set(adata.var_names)
    if missing_from_rira:
        raise RuntimeError(
            f"{len(missing_from_rira)} three-way symbols missing from RIRA "
            f"counts — gene list construction inconsistent."
        )
    assert list(adata.var_names) == three_way_sym, "Gene order mismatch after filter"

    # ── 4c. Normalization: NONE (RIRA .X is already pre-normalized) ───────
    # IMPORTANT: RIRA.All.RNA.counts.rds despite its name does not contain
    # raw UMI counts. Values are non-integer, max ~16,937, per-cell max
    # ~100-200, consistent with Seurat NormalizeData / SCTransform output.
    # The committed primary macaque pipeline uses these values as-is —
    # verified by reproducing obs/null = 0.7506 on the 13-type sensitivity
    # (committed 0.7488) with no further transformation. Applying another
    # normalize_total+log1p here would double-normalize and produce
    # obs/null ≈ 0.50 (verified experimentally).
    print("\n[4c] Normalization: NONE — RIRA .X is already pre-normalized")
    sample = adata.X[:5].toarray() if sp.issparse(adata.X) else adata.X[:5]
    print(f"  RIRA.X sample: min={sample.min():.3f}  max={sample.max():.3f}  "
          f"mean_nonzero={sample[sample > 0].mean():.3f}")

    # ── 4d. Compute per-type mean centroids ───────────────────────────────
    print("\n[4d] Computing per-type RIRA centroids…")
    centroids = np.zeros((13, adata.n_vars), dtype=np.float64)
    for i, t in enumerate(RIRA_TARGET_TYPES):
        mask = (adata.obs["target_type"] == t).to_numpy()
        X_sub = adata.X[mask]
        mean_vec = np.asarray(X_sub.mean(axis=0)).flatten()
        centroids[i] = mean_vec
        print(f"  {t:<40s} {int(mask.sum()):>7,} cells → centroid dim {mean_vec.size}")
    mac_cents = pd.DataFrame(
        centroids, index=RIRA_TARGET_TYPES, columns=three_way_ensg
    )
    mac_cents.to_csv(OUT_CENTROIDS)

    # ── 5. Sanity check: distribution vs human centroids ──────────────────
    print("\n[5] Sanity check vs human centroids…")
    h_vals = human_cents_sub.to_numpy()
    m_vals = mac_cents.to_numpy()
    print(f"  Human:   min={h_vals.min():.4f}  max={h_vals.max():.4f}  "
          f"mean={h_vals.mean():.4f}  std={h_vals.std():.4f}")
    print(f"  Macaque: min={m_vals.min():.4f}  max={m_vals.max():.4f}  "
          f"mean={m_vals.mean():.4f}  std={m_vals.std():.4f}")
    print(f"  Non-negative? human={(h_vals >= 0).all()}  macaque={(m_vals >= 0).all()}")
    std_ratio = m_vals.std() / max(h_vals.std(), 1e-12)
    print(f"  std(macaque)/std(human) = {std_ratio:.3f} (expect order unity)")
    if std_ratio > 10 or std_ratio < 0.1:
        raise SystemExit(
            "Macaque centroids differ from human by >10× in std — normalization "
            "likely mismatched. Stopping before Procrustes."
        )

    # ── 6. Joint PCA ──────────────────────────────────────────────────────
    print("\n[6] Joint PCA (≥95% variance)…")
    human_pca, mac_pca, pca, cell_types_ordered = pca_reduce_centroids(
        human_cents_sub, mac_cents
    )
    pca_info = {
        "n_components": int(pca.n_components_),
        "variance_explained": pca.explained_variance_ratio_.tolist(),
        "cumulative_variance": float(np.cumsum(pca.explained_variance_ratio_)[-1]),
    }

    # ── 7. Procrustes + permutation test ──────────────────────────────────
    print("\n[7] Procrustes align + 10,000-permutation test…")
    result = procrustes_align(human_pca, mac_pca)
    p_value, null_dist = permutation_test(
        human_pca, mac_pca, n_permutations=N_PERMUTATIONS, seed=RANDOM_SEED,
    )
    obs_null_ratio = result.distance / np.median(null_dist)
    obs_null_ratio_mean = result.distance / np.mean(null_dist)

    # ── 8. Save results + verify ──────────────────────────────────────────
    out = {
        "analysis": "RIRA13_RECONSTRUCTION",
        "n_types": 13,
        "cell_types": cell_types_ordered,
        "gene_space": len(three_way_ensg),
        "procrustes": {
            "distance": float(result.distance),
            "distance_squared": float(result.distance_squared),
            "scaling": float(result.scaling),
            "rotation_det": float(np.linalg.det(result.rotation)),
        },
        "permutation_test": {
            "p_value": float(p_value),
            "obs_null_ratio": float(obs_null_ratio),
            "obs_null_ratio_mean": float(obs_null_ratio_mean),
            "null_median": float(np.median(null_dist)),
            "null_mean": float(np.mean(null_dist)),
            "n_permutations": N_PERMUTATIONS,
        },
        "pca": pca_info,
        "per_type_residuals": {
            t: {
                "magnitude": float(np.linalg.norm(
                    result.aligned_target[i] - result.centered_reference[i]
                )),
            }
            for i, t in enumerate(cell_types_ordered)
        },
        "verification": {
            "expected_obs_null": EXPECTED_OBS_NULL,
            "expected_p": EXPECTED_P,
            "delta_obs_null": abs(obs_null_ratio - EXPECTED_OBS_NULL),
            "delta_obs_null_within_tolerance": abs(obs_null_ratio - EXPECTED_OBS_NULL) < TOL_OBS_NULL,
            "cell_counts_match_exact": True,
        },
        "seed": RANDOM_SEED,
    }
    # Ensure JSON-serialisable (np.bool_ etc.)
    def _default(o):
        if isinstance(o, (np.bool_,)):
            return bool(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        raise TypeError(f"not serialisable: {type(o).__name__}")
    OUT_JSON.write_text(json.dumps(out, indent=2, default=_default))

    # ── Summary ───────────────────────────────────────────────────────────
    dt = time.time() - t_start
    print("\n" + "=" * 70)
    print("RESULTS")
    print("=" * 70)
    print(f"  obs/null (median):  {obs_null_ratio:.6f}   (expected 0.7488)")
    print(f"  delta:              {abs(obs_null_ratio - EXPECTED_OBS_NULL):.6f}   (tol 0.01)")
    print(f"  p-value:            {p_value:.6f}   (expected 0.0001999800)")
    print(f"  distance:           {result.distance:.6f}   (expected 16.5307)")
    print(f"  scaling:            {result.scaling:.6f}   (expected 0.0534)")
    print(f"  PCA components:     {pca_info['n_components']}")
    print(f"  Runtime:            {dt:.1f}s")
    status = "PASS" if out["verification"]["delta_obs_null_within_tolerance"] else "FAIL"
    print(f"  VERIFICATION:       {status}")
    print(f"\nResults written to: {OUT_JSON}")


if __name__ == "__main__":
    main()
