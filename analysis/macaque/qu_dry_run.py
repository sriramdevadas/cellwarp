#!/usr/bin/env python3
"""Dry run for Qu-only pipeline — exercises everything except cell-type
filtering and centroid computation, so that when MacFas.meta.data.rds
finishes downloading, only the label-mapping step remains.

Steps:
  1. Inspect barcode format in all 20 per-tissue 10x bundles (length,
     suffixes, duplicates).
  2. Verify Qu feature lists are identical across tissues.
  3. Build Qu-ENSMFAG → human-ENSG mapping to the 13,927-gene space.
     Sanity: the mapping should cover Qu's overlap with three-way (13,927
     of which 14,558-subset-minus-RIRA-gap = 13,927 guaranteed if mapping
     is sound).
  4. Aggregate all 20 Qu 10x matrices → one sparse matrix in the 13,927-gene
     space, with per-cell 'tissue' annotation.
  5. Sanity-check value range on the aggregated matrix (should be pure
     integer UMIs, matching the per-tissue checks earlier).
  6. Apply normalize_total(1e4) + log1p. Report post-normalization stats.
  7. Stop — centroid computation requires cell-type annotations we don't
     have yet.

Output artifact:
  output/macaque_pipeline/qu_dry_run_aggregate.h5ad
    (aggregate of all 230,882 Qu cells × 13,927 ortholog genes, RAW counts
    pre-normalization — reusable by the real pipeline once annotations arrive)
"""
from __future__ import annotations
import gzip, glob, re, sys, time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io as sio
import scipy.sparse as sp

PROJECT = Path(__file__).resolve().parent.parent.parent
QU_DIR = PROJECT / "data/macaque/qu_2022/extracted"
BIOMART = PROJECT / "data/macaque/biomart_macaque_human_orthologs.csv"
THREEWAY_GENELIST = PROJECT / "output/macaque_pipeline/reconstruction_rira13_gene_list.csv"
OUT_AGG = PROJECT / "output/macaque_pipeline/qu_dry_run_aggregate.h5ad"


def discover_tissues():
    matrix_files = sorted(glob.glob(str(QU_DIR / "GSM*_matrix.mtx.gz")))
    tissues = []
    for m in matrix_files:
        stem = Path(m).name
        match = re.match(r"(GSM\d+)_(.+)_matrix\.mtx\.gz", stem)
        if match:
            tissues.append({
                "gsm": match.group(1),
                "tissue": match.group(2),
                "matrix": m,
                "barcodes": str(QU_DIR / f"{match.group(1)}_{match.group(2)}_barcodes.tsv.gz"),
                "features": str(QU_DIR / f"{match.group(1)}_{match.group(2)}_features.tsv.gz"),
            })
    return tissues


def inspect_barcodes(tissues):
    print("\n[1] Barcode format inspection")
    print(f"{'tissue':<12} {'n_cells':>8}  {'len':<4}  {'suffix':<6}  {'sample':<25}")
    for t in tissues:
        with gzip.open(t["barcodes"], "rt") as f:
            bcs = [line.strip() for line in f]
        # Format signature
        sample_bc = bcs[0]
        suffix = "-1" if sample_bc.endswith("-1") else "none"
        lengths = set(len(bc) for bc in bcs[:200])
        dup_count = len(bcs) - len(set(bcs))
        print(f"  {t['tissue']:<12} {len(bcs):>8,}  {list(lengths)}  {suffix:<6}  {sample_bc:<25}  dupes={dup_count}")


def verify_feature_consistency(tissues):
    print("\n[2] Feature list consistency across tissues")
    sets = []
    for t in tissues[:5]:
        with gzip.open(t["features"], "rt") as f:
            s = set(line.rstrip("\n").split("\t")[0] for line in f)
        sets.append(s)
    all_equal = all(s == sets[0] for s in sets)
    print(f"  First 5 tissues have identical feature sets: {all_equal}  (n_features = {len(sets[0]):,})")
    # Check all 20 against tissue 0
    with gzip.open(tissues[0]["features"], "rt") as f:
        rows0 = [line.rstrip("\n").split("\t") for line in f]
    for t in tissues:
        with gzip.open(t["features"], "rt") as f:
            rows = [line.rstrip("\n").split("\t") for line in f]
        if rows != rows0:
            print(f"  WARN: {t['tissue']} has different feature list")
    return rows0  # list of (ensmfag, symbol, type)


def build_qu_to_human_map(qu_features: list[list[str]]) -> dict[str, str]:
    """Map Qu ENSMFAG feature → human ENSG in the 13,927 three-way space.

    Reuses the three-way gene list built by the RIRA reconstruction (each row
    is (ensg, symbol)). For each Qu feature:
      - If its HGNC-like symbol (col 2, when not ENSMFAG) matches a three-way
        symbol: direct map to that ENSG.
      - Else use BioMart ENSMFAG → human gene name and match to three-way.
    Returns dict ENSMFAG → human_ENSG (only for features that are in the
    three-way set).
    """
    threeway = pd.read_csv(THREEWAY_GENELIST)
    sym_to_ensg = dict(zip(threeway["symbol"], threeway["ensg"]))
    three_symbols = set(sym_to_ensg.keys())
    print(f"  Three-way symbols loaded: {len(three_symbols):,}")

    biomart = pd.read_csv(BIOMART)
    bm_121 = biomart[biomart["Human homology type"] == "ortholog_one2one"]
    biomart_map = dict(zip(bm_121["Gene stable ID"].astype(str),
                           bm_121["Human gene name"].fillna("").astype(str)))

    qu_to_human: dict[str, str] = {}
    n_direct = n_biomart = 0
    for ensmfag, sym, _ in qu_features:
        human_sym = None
        if not sym.startswith("ENSMFAG") and sym != ensmfag and sym != "":
            if sym in three_symbols:
                human_sym = sym
                n_direct += 1
        if human_sym is None:
            bm_sym = biomart_map.get(ensmfag, "")
            if bm_sym in three_symbols:
                human_sym = bm_sym
                n_biomart += 1
        if human_sym is not None:
            qu_to_human[ensmfag] = sym_to_ensg[human_sym]
    print(f"  Qu features in three-way: {len(qu_to_human):,}  "
          f"(direct={n_direct:,}, biomart={n_biomart:,})")
    return qu_to_human


def aggregate(tissues, qu_to_human, qu_features):
    """Load each tissue's matrix, gene-filter, vstack."""
    print("\n[4] Aggregate Qu 10x matrices")
    # Build column index in three-way gene space
    # Preserve three-way human ENSG order (canonical) for compatibility with
    # the RIRA reconstruction centroids.
    threeway = pd.read_csv(THREEWAY_GENELIST)
    ensg_order = list(threeway["ensg"])
    ensg_to_col = {e: i for i, e in enumerate(ensg_order)}

    # Map each Qu feature row → column in three-way, or -1 if excluded
    n_qu_feat = len(qu_features)
    qu_row_to_col = np.full(n_qu_feat, -1, dtype=np.int64)
    for r, (ensmfag, _sym, _typ) in enumerate(qu_features):
        target_ensg = qu_to_human.get(ensmfag)
        if target_ensg is not None:
            qu_row_to_col[r] = ensg_to_col[target_ensg]

    # Build column-collapse matrix: shape (n_qu_feat, n_three_way).
    # For each Qu feature that maps to a three-way ENSG, this matrix sums
    # its counts into that ENSG column. Resulting matrix dimensions are
    # exactly n_three_way — no duplicates.
    kept_rows = np.where(qu_row_to_col >= 0)[0]
    kept_cols = qu_row_to_col[kept_rows]
    n_ensg = len(ensg_order)
    collapse = sp.coo_matrix(
        (np.ones(len(kept_rows), dtype=np.float64),
         (kept_rows, kept_cols)),
        shape=(n_qu_feat, n_ensg),
    ).tocsr()
    # Count how many Qu rows collapse per ENSG (sanity)
    fan_in = np.bincount(kept_cols, minlength=n_ensg)
    n_dup_ensg = int((fan_in > 1).sum())
    print(f"  Human ENSGs with >1 Qu-feature (summed): {n_dup_ensg}  "
          f"(max fan-in = {int(fan_in.max())})")

    # Per-tissue list of sparse matrices
    parts = []
    obs_rows = []
    n_cells_total = 0
    for t in tissues:
        t_start = time.time()
        with gzip.open(t["matrix"], "rt") as f:
            mat = sio.mmread(f)
        X = sp.csr_matrix(mat.T)  # (cells, n_qu_feat)
        with gzip.open(t["barcodes"], "rt") as f:
            bcs = [line.strip() for line in f]
        assert X.shape[0] == len(bcs)
        # X (cells × n_qu_feat) @ collapse (n_qu_feat × n_ensg) → cells × n_ensg
        X_filt = X @ collapse
        n_cells_tissue = X_filt.shape[0]
        parts.append(X_filt)
        obs_rows.extend([(bc, t["gsm"], t["tissue"]) for bc in bcs])
        n_cells_total += n_cells_tissue
        t_el = time.time() - t_start
        print(f"  {t['tissue']:<12} {n_cells_tissue:>8,} cells × {X_filt.shape[1]:,} genes  ({t_el:.1f}s)")

    X_all = sp.vstack(parts)
    print(f"  Total: {X_all.shape}  nnz={X_all.nnz:,}")
    obs_df = pd.DataFrame(obs_rows, columns=["barcode", "gsm", "tissue"])
    obs_df.index = obs_df["barcode"].astype(str) + "_" + obs_df["tissue"]
    # We include tissue tag in the AnnData row index to guarantee uniqueness
    # across tissues even when the same barcode repeats.

    # Verify sanity of values (should still be pure integers)
    sample = X_all.data[:100000]
    int_frac = (sample == sample.astype(int)).mean()
    print(f"  Sample integer fraction: {int_frac:.4%}  max={X_all.data.max():.2f}  "
          f"negatives={int((X_all.data < 0).sum())}")
    return X_all, obs_df, ensg_order


def main():
    t0 = time.time()
    print("=" * 70)
    print("Qu dry run — aggregate raw 10x bundles into three-way gene space")
    print("=" * 70)

    tissues = discover_tissues()
    print(f"Discovered {len(tissues)} tissues")

    inspect_barcodes(tissues)

    print("\n[2] Verify feature consistency")
    qu_features = verify_feature_consistency(tissues)

    print("\n[3] Build Qu → human ENSG map")
    qu_to_human = build_qu_to_human_map(qu_features)

    X_all, obs_df, ensg_order = aggregate(tissues, qu_to_human, qu_features)

    # Save pre-normalization aggregate
    adata = ad.AnnData(
        X=X_all,
        obs=obs_df,
        var=pd.DataFrame(index=ensg_order),
    )
    # Need unique obs index — already built above
    OUT_AGG.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(OUT_AGG, compression="gzip")
    print(f"\n  Wrote: {OUT_AGG}")

    # Post-norm stats (dry run)
    print("\n[6] Applying normalize_total(1e4) + log1p to full aggregate…")
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)
    print(f"  Post-norm sample: min={adata_norm.X.data.min():.4f}  "
          f"max={adata_norm.X.data.max():.4f}  "
          f"mean(nonzero)={adata_norm.X.data.mean():.4f}")

    dt = time.time() - t0
    print(f"\nTotal dry-run time: {dt:.1f}s")
    print("\nReady for next step: filter to 7 target types once Qu @meta.data arrives.")


if __name__ == "__main__":
    main()
