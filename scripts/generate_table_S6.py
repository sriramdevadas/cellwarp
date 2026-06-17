#!/usr/bin/env python3
"""
Generate Supplementary Table S6: CPC1 driver genes for all 35 matched cell types.

CPC1 = Common Principal Component 1 (Krzanowski 1984): the dominant shared axis of
within-type variation across species. Computed by weighting both species' within-type
covariance matrices (in 33-D PCA space) and eigendecomposing the sum, then projecting
the top eigenvector back to gene space.

This is the same computation as t3b_ellipsoid_alignment.py Step 10, extended from
top-5 to top-20 gene loadings for the supplementary table.

Biology: Within each cell type, cells vary along certain gene axes. CPC1 captures
the dominant shared axis. In 25 of 35 types, the highest-loading gene is a ribosomal
protein (RPL/RPS), reflecting conserved translational state variation. In 10 types,
CPC1 instead reflects cell-type-specific biology.

Output:
  - Table_S6_CPC1_driver_genes.xlsx (two sheets: summary + full loadings)
  - Table_S6_validation.txt (consistency checks)
"""

from __future__ import annotations

import re
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.decomposition import PCA

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT = Path(__file__).resolve().parent.parent
DATA_SCALED = PROJECT / "data" / "phase2_scaled"
OUTPUT_SCALED = PROJECT / "output" / "phase2" / "scaled_35types"
OUTPUT_DIR = PROJECT / "output" / "supplementary"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

SEED = 42
TOP_N = 20
N_PCA = 33  # Same as primary analysis


def is_ribosomal(symbol: str) -> bool:
    """Check if a gene symbol is a ribosomal protein gene.

    Handles both human (RPL/RPS) and mouse (Rpl/Rps) nomenclature.
    """
    return bool(re.match(r"^(RPL|RPS|Rpl|Rps)\d", symbol))


def main():
    t_start = time.time()

    # ------------------------------------------------------------------
    # Step 0: Load data
    # ------------------------------------------------------------------
    print("Loading data files...")
    h_adata = ad.read_h5ad(DATA_SCALED / "human_scaled.h5ad")
    m_adata = ad.read_h5ad(DATA_SCALED / "mouse_scaled.h5ad")

    # Gene name mappings (Ensembl ID → symbol)
    ensembl_to_human = dict(zip(h_adata.var_names, h_adata.var["feature_name"]))
    ensembl_to_mouse = dict(zip(m_adata.var_names, m_adata.var["feature_name"]))

    # Rigidity rankings
    ranks_df = pd.read_csv(OUTPUT_SCALED / "residuals_ranked.csv")
    rank_map = dict(zip(ranks_df["cell_type"], ranks_df["rank"]))

    # ------------------------------------------------------------------
    # Step 1: Reconstruct PCA model (identical to t3b Step 0)
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("Reconstructing PCA model from centroids")
    print("=" * 70)

    saved = np.load(OUTPUT_SCALED / "pca_centroids_35.npz", allow_pickle=True)
    cell_types = list(saved["cell_types"])
    assert len(cell_types) == 35

    hc = pd.read_csv(OUTPUT_SCALED / "centroids_human_35.csv", index_col=0)
    mc = pd.read_csv(OUTPUT_SCALED / "centroids_mouse_35.csv", index_col=0)
    hc = hc.loc[cell_types]
    mc = mc.loc[cell_types]
    gene_ids = list(hc.columns)  # Ensembl IDs

    combined = np.vstack([hc.values, mc.values])  # (70, 16959)
    pca = PCA(n_components=N_PCA, svd_solver="full", random_state=SEED)
    pca.fit_transform(combined)

    W = pca.components_  # (33, 16959)
    pca_mean = pca.mean_  # (16959,)
    n_genes = len(gene_ids)

    print(f"  PCA: {pca.n_components_} components, "
          f"{sum(pca.explained_variance_ratio_) * 100:.1f}% variance")
    print(f"  Gene space: {n_genes} genes")

    # ------------------------------------------------------------------
    # Step 2: Compute CPC1 loadings for all 35 types (Krzanowski approach)
    # Identical to t3b Step 10, but extracting top 20 instead of top 5
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("Computing CPC1 loadings (Krzanowski approach, top 20)")
    print("=" * 70)

    summary_rows = []
    full_rows = []

    for ct in cell_types:
        # --- Project cells to PCA space and compute covariance ---
        species_cov = {}
        species_n = {}
        for sp_label, adata in [("human", h_adata), ("mouse", m_adata)]:
            mask = adata.obs["cell_type"].values == ct
            X = adata[mask][:, gene_ids].X
            if sp.issparse(X):
                X = X.toarray()
            X = np.asarray(X, dtype=np.float64)
            X_pca = (X - pca_mean) @ W.T  # (n_cells, 33)

            centroid = X_pca.mean(axis=0)
            centered = X_pca - centroid
            n = centered.shape[0]
            cov = (centered.T @ centered) / (n - 1)

            species_cov[sp_label] = cov
            species_n[sp_label] = n

        # --- Weighted sum of covariance matrices → CPC1 ---
        n_h, n_m = species_n["human"], species_n["mouse"]
        weighted_cov = n_h * species_cov["human"] + n_m * species_cov["mouse"]

        eigvals, eigvecs = np.linalg.eigh(weighted_cov)
        idx = np.argsort(eigvals)[::-1]
        eigvecs = eigvecs[:, idx]

        cpc1 = eigvecs[:, 0]  # (33,)

        # --- Variance fractions in each species ---
        var_h = cpc1 @ species_cov["human"] @ cpc1
        var_m = cpc1 @ species_cov["mouse"] @ cpc1
        total_h = np.trace(species_cov["human"])
        total_m = np.trace(species_cov["mouse"])
        frac_h = float(var_h / total_h) if total_h > 0 else 0.0
        frac_m = float(var_m / total_m) if total_m > 0 else 0.0

        # --- Project CPC1 back to gene space ---
        cpc1_genes = cpc1 @ W  # (16959,)
        abs_loadings = np.abs(cpc1_genes)
        top_idx = np.argsort(abs_loadings)[::-1][:TOP_N]

        # --- Extract top 20 genes with both human and mouse names ---
        top20_human = []
        top20_mouse = []
        n_ribo_human = 0
        n_ribo_mouse = 0

        for rank_i, g_idx in enumerate(top_idx, 1):
            ens_id = gene_ids[g_idx]
            h_symbol = ensembl_to_human.get(ens_id, ens_id)
            m_symbol = ensembl_to_mouse.get(ens_id, ens_id)
            loading = float(cpc1_genes[g_idx])
            abs_load = float(abs_loadings[g_idx])
            ribo_h = is_ribosomal(h_symbol)
            ribo_m = is_ribosomal(m_symbol)
            if ribo_h:
                n_ribo_human += 1
            if ribo_m:
                n_ribo_mouse += 1

            top20_human.append({
                "gene_name": h_symbol, "ensembl_id": ens_id,
                "loading_rank": rank_i, "loading_value": loading,
                "abs_loading": abs_load, "is_ribosomal": ribo_h,
            })
            top20_mouse.append({
                "gene_name": m_symbol, "ensembl_id": ens_id,
                "loading_rank": rank_i, "loading_value": loading,
                "abs_loading": abs_load, "is_ribosomal": ribo_m,
            })

        ribo_frac_h = n_ribo_human / TOP_N
        ribo_frac_m = n_ribo_mouse / TOP_N

        # Classification: rank-1 gene is ribosomal → "ribosomal-dominated"
        rank1_symbol = ensembl_to_human.get(gene_ids[top_idx[0]], "")
        classification = "ribosomal-dominated" if is_ribosomal(rank1_symbol) else "cell-type-specific"

        h_top5 = ", ".join(g["gene_name"] for g in top20_human[:5])
        m_top5 = ", ".join(g["gene_name"] for g in top20_mouse[:5])

        # Gene overlap: CPC1 is a shared direction, so the same 20 orthologs
        # appear in both species' top-20 (overlap = 20 by construction).
        overlap = TOP_N

        print(f"  {ct}: var_h={frac_h:.1%}, var_m={frac_m:.1%}, "
              f"rank1={rank1_symbol}, class={classification}")

        summary_rows.append({
            "cell_type": ct,
            "rank": rank_map[ct],
            "human_CPC1_var_explained": round(frac_h * 100, 2),
            "mouse_CPC1_var_explained": round(frac_m * 100, 2),
            "human_top5_genes": h_top5,
            "mouse_top5_genes": m_top5,
            "human_ribosomal_frac_top20": round(ribo_frac_h, 2),
            "mouse_ribosomal_frac_top20": round(ribo_frac_m, 2),
            "classification": classification,
            "human_mouse_CPC1_genes_overlap": overlap,
        })

        # Full loadings rows (human and mouse naming)
        for g in top20_human:
            full_rows.append({
                "cell_type": ct,
                "species": "human",
                "gene_name": g["gene_name"],
                "loading_rank": g["loading_rank"],
                "loading_value": round(g["loading_value"], 6),
                "abs_loading": round(g["abs_loading"], 6),
                "is_ribosomal": g["is_ribosomal"],
            })
        for g in top20_mouse:
            full_rows.append({
                "cell_type": ct,
                "species": "mouse",
                "gene_name": g["gene_name"],
                "loading_rank": g["loading_rank"],
                "loading_value": round(g["loading_value"], 6),
                "abs_loading": round(g["abs_loading"], 6),
                "is_ribosomal": g["is_ribosomal"],
            })

    # ------------------------------------------------------------------
    # Build DataFrames
    # ------------------------------------------------------------------
    summary_df = pd.DataFrame(summary_rows)
    summary_df = summary_df.sort_values("rank").reset_index(drop=True)

    full_df = pd.DataFrame(full_rows)

    # ------------------------------------------------------------------
    # Save Excel
    # ------------------------------------------------------------------
    xlsx_path = OUTPUT_DIR / "Table_S6_CPC1_driver_genes.xlsx"
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        summary_df.to_excel(writer, sheet_name="CPC1_summary", index=False)
        full_df.to_excel(writer, sheet_name="CPC1_full_loadings", index=False)

    print(f"\n  Saved: {xlsx_path}")
    print(f"  Sheet 1: {len(summary_df)} rows (CPC1_summary)")
    print(f"  Sheet 2: {len(full_df)} rows (CPC1_full_loadings)")

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("VALIDATION")
    print("=" * 70)

    n_ribo_final = sum(1 for _, r in summary_df.iterrows()
                       if r["classification"] == "ribosomal-dominated")
    n_specific_final = sum(1 for _, r in summary_df.iterrows()
                          if r["classification"] == "cell-type-specific")
    print(f"  Ribosomal-dominated: {n_ribo_final}")
    print(f"  Cell-type-specific:  {n_specific_final}")

    specific_types = summary_df[summary_df["classification"] == "cell-type-specific"]
    print(f"\n  Cell-type-specific types ({len(specific_types)}):")
    for _, row in specific_types.iterrows():
        top3 = ", ".join(row["human_top5_genes"].split(", ")[:3])
        print(f"    rank {row['rank']:2d}: {row['cell_type']:<50} top3: {top3}")

    discrepancies = []
    if n_ribo_final != 25:
        discrepancies.append(
            f"DISCREPANCY: {n_ribo_final} ribosomal-dominated (expected 25)")
    if n_specific_final != 10:
        discrepancies.append(
            f"DISCREPANCY: {n_specific_final} cell-type-specific (expected 10)")

    # ------------------------------------------------------------------
    # Save validation log
    # ------------------------------------------------------------------
    log_lines = [
        "Table S6 Validation Log",
        "=" * 50,
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}",
        f"Random seed: {SEED}",
        f"Census version: 2025-11-08 (pinned)",
        f"Gene space: {n_genes} shared 1:1 orthologs",
        f"Top-N genes: {TOP_N}",
        f"PCA dimensions: {N_PCA} components ({sum(pca.explained_variance_ratio_) * 100:.1f}% variance)",
        f"Classification criterion: rank-1 CPC1 gene is ribosomal protein (RPL/RPS)",
        "",
        "METHOD",
        "  CPC1 = Common Principal Component 1 (Krzanowski 1984).",
        "  Within-type covariance matrices computed in 33-D PCA space",
        "  (fitted on 70 centroids = 35 types x 2 species).",
        "  Weighted sum of human and mouse covariance eigendecomposed;",
        "  top eigenvector projected to gene space via PCA loading matrix.",
        "  Gene loadings are shared across species (CPC is a common direction).",
        "  Variance fractions differ per species (how much of each species'",
        "  within-type variance lies along the shared CPC1 axis).",
        "",
        "COUNTS",
        f"  Ribosomal-dominated: {n_ribo_final}",
        f"  Cell-type-specific:  {n_specific_final}",
        f"  Expected:            25 ribosomal / 10 specific",
        "",
        "CELL-TYPE-SPECIFIC TYPES (with top 3 genes):",
    ]
    for _, row in specific_types.iterrows():
        top3 = ", ".join(row["human_top5_genes"].split(", ")[:3])
        log_lines.append(
            f"  rank {row['rank']:2d}: {row['cell_type']:<50} {top3}")

    if discrepancies:
        log_lines.append("")
        log_lines.append("DISCREPANCIES:")
        for d in discrepancies:
            log_lines.append(f"  {d}")
    else:
        log_lines.append("")
        log_lines.append("NO DISCREPANCIES — matches manuscript claims (25/10 split).")

    log_lines.append("")
    log_lines.append("NOTE: CPC1 is a shared direction, so human and mouse top-20")
    log_lines.append("genes are the same orthologs (overlap = 20 by construction).")
    log_lines.append("Human and mouse columns use species-appropriate gene nomenclature.")

    log_path = OUTPUT_DIR / "Table_S6_validation.txt"
    with open(log_path, "w") as f:
        f.write("\n".join(log_lines) + "\n")
    print(f"\n  Saved: {log_path}")

    elapsed = time.time() - t_start
    print(f"\n  Done in {elapsed:.0f}s")


if __name__ == "__main__":
    main()
