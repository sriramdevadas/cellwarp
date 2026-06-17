#!/usr/bin/env python3
"""
CellWarp — Hepatocyte Zonation Diagnostic

Targeted diagnostic to determine whether the hepatocyte rank reversal in T3-C
(aggregate rank 4 → tissue-stratified rank 28, Δ=24) reflects a genuine finding
or a Sun2023 Alb-rescue annotation artifact.

Biology
-------
Hepatic zonation is a fundamental organizational principle of the liver.
Hepatocytes near the central vein (pericentral zone) express Cyp2e1, Glul,
Cyp1a2, Axin2. Hepatocytes near the portal triad (periportal zone) express
Ass1, Pck1, Hal, Sds, Arg1. If the Sun2023 Alb-rescue annotation selectively
captures pericentral hepatocytes (high Alb expressors), the Sun2023 centroid
will be pericentral-shifted relative to Tabula, producing a centroid mismatch
that inflates the cross-species distance.

Pipeline
--------
  Task 0: Locate hepatocytes across datasets
  Task 1: Zonation marker expression comparison
  Task 2: Centroid shift in PCA space
  Task 3: Alb-rescue bias check
  Task 4: Verdict and documentation

Output
------
  output/validation/hepatocyte_zonation_diagnostic/
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

OUTPUT_DIR = Path("output/validation/hepatocyte_zonation_diagnostic")
TABULA_HUMAN_PATH = Path("data/phase2_scaled/human_scaled.h5ad")
TABULA_MOUSE_PATH = Path("data/phase2_scaled/mouse_scaled.h5ad")
SUN2023_PATH = Path("data/replication/sun2023/sun2023_yc.h5ad")
ORTHOLOG_PATH = Path("data/phase1/orthologs_human_mouse.csv")

# Primary PCA centroids (for reconstructing PCA model)
CENTROIDS_HUMAN_PATH = Path("output/phase2/scaled_35types/centroids_human_35.csv")
CENTROIDS_MOUSE_PATH = Path("output/phase2/scaled_35types/centroids_mouse_35.csv")
PCA_CENTROIDS_PATH = Path("output/phase2/scaled_35types/pca_centroids_35.npz")
RESIDUALS_RANKED_PATH = Path("output/phase2/scaled_35types/residuals_ranked.csv")

# T3-C tissue-stratified distances
T3C_RANKING_PATH = Path("output/validation/t3c_tissue_stratified/tissue_stratified_ranking.csv")

RANDOM_SEED = 42
PCA_VARIANCE_THRESHOLD = 0.95

# Zonation marker genes (human Ensembl IDs — all confirmed present in ortholog space)
PERICENTRAL_GENES = {
    "CYP2E1": "ENSG00000130649",
    "GLUL": "ENSG00000135821",
    "CYP1A2": "ENSG00000140505",
    "AXIN2": "ENSG00000168646",
}
PERIPORTAL_GENES = {
    "ASS1": "ENSG00000130707",
    "PCK1": "ENSG00000124253",
    "HAL": "ENSG00000084110",
    "SDS": "ENSG00000135094",
    "ARG1": "ENSG00000118520",
}
ALB_GENE = {"ALB": "ENSG00000163631"}


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    ortho = pd.read_csv(ORTHOLOG_PATH)
    human_to_name = dict(zip(ortho["human_ensembl_id"], ortho["human_gene_name"]))

    # ==================================================================
    # TASK 0: Locate hepatocyte data across all datasets
    # ==================================================================
    print("=" * 72)
    print("TASK 0: Locate hepatocyte data across all datasets")
    print("=" * 72)

    # --- Tabula Human ---
    h = ad.read_h5ad(TABULA_HUMAN_PATH)
    h_hep = h[h.obs["cell_type"] == "hepatocyte"]
    h_hep_tissues = h_hep.obs["tissue_general"].value_counts()
    print(f"\nTabula Human hepatocytes: {h_hep.n_obs} cells")
    print(f"  Tissues: {dict(h_hep_tissues)}")
    has_alb_rescue_flag_h = "alb_rescued" in h.obs.columns
    print(f"  Alb-rescue flag in .obs: {has_alb_rescue_flag_h}")

    # --- Tabula Mouse ---
    m = ad.read_h5ad(TABULA_MOUSE_PATH)
    m_hep = m[m.obs["cell_type"] == "hepatocyte"]
    m_hep_tissues = m_hep.obs["tissue_general"].value_counts()
    print(f"\nTabula Mouse hepatocytes: {m_hep.n_obs} cells")
    print(f"  Tissues: {dict(m_hep_tissues)}")
    has_alb_rescue_flag_m = "alb_rescued" in m.obs.columns
    print(f"  Alb-rescue flag in .obs: {has_alb_rescue_flag_m}")

    # --- Sun2023 ---
    s = ad.read_h5ad(SUN2023_PATH)
    s_hep = s[s.obs["cell_type"] == "hepatocyte"]
    s_hep_tissues = s_hep.obs["tissue"].value_counts()
    print(f"\nSun2023 hepatocytes: {s_hep.n_obs} cells")
    print(f"  Tissues: {dict(s_hep_tissues)}")
    has_alb_rescue_flag_s = "alb_rescued" in s.obs.columns
    print(f"  Alb-rescue flag in .obs: {has_alb_rescue_flag_s}")
    print(f"  obs columns: {list(s.obs.columns)}")

    # Confirm all liver-only
    all_liver = True
    for name, tissues in [("Tabula Human", h_hep_tissues),
                           ("Tabula Mouse", m_hep_tissues),
                           ("Sun2023", s_hep_tissues)]:
        non_liver = {k: v for k, v in tissues.items() if "liver" not in k.lower()}
        if non_liver:
            print(f"  WARNING: {name} has non-liver hepatocytes: {non_liver}")
            all_liver = False
    if all_liver:
        print("\n  CONFIRMED: All hepatocytes are liver-only across all datasets.")

    task0_results = {
        "tabula_human": {"n_cells": int(h_hep.n_obs), "tissues": dict(h_hep_tissues)},
        "tabula_mouse": {"n_cells": int(m_hep.n_obs), "tissues": dict(m_hep_tissues)},
        "sun2023": {"n_cells": int(s_hep.n_obs), "tissues": dict(s_hep_tissues)},
        "all_liver_only": all_liver,
    }

    # ==================================================================
    # TASK 1: Zonation marker gene expression comparison
    # ==================================================================
    print(f"\n{'='*72}")
    print("TASK 1: Zonation marker gene expression comparison")
    print(f"{'='*72}")

    # Build gene index for the shared gene space
    var_names = list(h.var_names)  # Same across all 3 datasets
    gene_idx = {g: i for i, g in enumerate(var_names)}

    # Check which zonation markers are present
    print("\nZonation marker availability in 16,959 ortholog gene space:")
    all_markers = {}
    for category, genes in [("Pericentral", PERICENTRAL_GENES),
                             ("Periportal", PERIPORTAL_GENES),
                             ("ALB", ALB_GENE)]:
        for name, ensembl_id in genes.items():
            present = ensembl_id in gene_idx
            all_markers[name] = {"ensembl_id": ensembl_id, "present": present, "category": category}
            status = "PRESENT" if present else "ABSENT"
            print(f"  {category:<12} {name:<10} ({ensembl_id}) — {status}")

    # Compute mean expression per marker per dataset (hepatocytes only)
    def get_mean_expression(adata_subset, ensembl_ids):
        """Get mean expression for a list of Ensembl IDs across cells."""
        result = {}
        X = adata_subset.X
        for name, eid in ensembl_ids.items():
            if eid in gene_idx:
                idx = gene_idx[eid]
                vals = X[:, idx]
                if hasattr(vals, 'toarray'):
                    vals = vals.toarray().flatten()
                else:
                    vals = np.asarray(vals).flatten()
                result[name] = {
                    "mean": float(np.mean(vals)),
                    "median": float(np.median(vals)),
                    "pct_nonzero": float(np.mean(vals > 0) * 100),
                }
        return result

    print("\nPericentral marker expression in hepatocytes:")
    h_peri_c = get_mean_expression(h_hep, PERICENTRAL_GENES)
    m_peri_c = get_mean_expression(m_hep, PERICENTRAL_GENES)
    s_peri_c = get_mean_expression(s_hep, PERICENTRAL_GENES)

    print(f"  {'Gene':<10} {'Tab Human mean':>15} {'Tab Mouse mean':>15} {'Sun2023 mean':>15}")
    print(f"  {'-'*58}")
    for gene in PERICENTRAL_GENES:
        h_val = h_peri_c.get(gene, {}).get("mean", float("nan"))
        m_val = m_peri_c.get(gene, {}).get("mean", float("nan"))
        s_val = s_peri_c.get(gene, {}).get("mean", float("nan"))
        print(f"  {gene:<10} {h_val:>15.4f} {m_val:>15.4f} {s_val:>15.4f}")

    print("\nPeriportal marker expression in hepatocytes:")
    h_peri_p = get_mean_expression(h_hep, PERIPORTAL_GENES)
    m_peri_p = get_mean_expression(m_hep, PERIPORTAL_GENES)
    s_peri_p = get_mean_expression(s_hep, PERIPORTAL_GENES)

    print(f"  {'Gene':<10} {'Tab Human mean':>15} {'Tab Mouse mean':>15} {'Sun2023 mean':>15}")
    print(f"  {'-'*58}")
    for gene in PERIPORTAL_GENES:
        h_val = h_peri_p.get(gene, {}).get("mean", float("nan"))
        m_val = m_peri_p.get(gene, {}).get("mean", float("nan"))
        s_val = s_peri_p.get(gene, {}).get("mean", float("nan"))
        print(f"  {gene:<10} {h_val:>15.4f} {m_val:>15.4f} {s_val:>15.4f}")

    # Compute aggregate pericentral vs periportal scores
    def aggregate_zone_score(zone_expr):
        """Mean of mean expression across all markers in a zone."""
        vals = [v["mean"] for v in zone_expr.values()]
        return np.mean(vals) if vals else 0.0

    tab_h_pc = aggregate_zone_score(h_peri_c)
    tab_h_pp = aggregate_zone_score(h_peri_p)
    tab_m_pc = aggregate_zone_score(m_peri_c)
    tab_m_pp = aggregate_zone_score(m_peri_p)
    sun_pc = aggregate_zone_score(s_peri_c)
    sun_pp = aggregate_zone_score(s_peri_p)

    tab_h_ratio = tab_h_pc / tab_h_pp if tab_h_pp > 0 else float("inf")
    tab_m_ratio = tab_m_pc / tab_m_pp if tab_m_pp > 0 else float("inf")
    sun_ratio = sun_pc / sun_pp if sun_pp > 0 else float("inf")

    print(f"\nAggregate zonation scores (mean of marker means):")
    print(f"  {'Dataset':<18} {'Pericentral':>12} {'Periportal':>12} {'PC/PP ratio':>12}")
    print(f"  {'-'*56}")
    print(f"  {'Tabula Human':<18} {tab_h_pc:>12.4f} {tab_h_pp:>12.4f} {tab_h_ratio:>12.4f}")
    print(f"  {'Tabula Mouse':<18} {tab_m_pc:>12.4f} {tab_m_pp:>12.4f} {tab_m_ratio:>12.4f}")
    print(f"  {'Sun2023':<18} {sun_pc:>12.4f} {sun_pp:>12.4f} {sun_ratio:>12.4f}")

    pericentral_shifted = sun_ratio > max(tab_h_ratio, tab_m_ratio)
    print(f"\n  Sun2023 pericentral-shifted vs both Tabula: {pericentral_shifted}")
    if pericentral_shifted:
        print(f"  Sun2023 PC/PP ratio ({sun_ratio:.3f}) > Tabula max ({max(tab_h_ratio, tab_m_ratio):.3f})")

    # Detection rates
    print(f"\nDetection rates (% cells with expression > 0):")
    print(f"  {'Gene':<10} {'Tab Human %':>12} {'Tab Mouse %':>12} {'Sun2023 %':>12}")
    print(f"  {'-'*50}")
    for gene in list(PERICENTRAL_GENES) + list(PERIPORTAL_GENES) + list(ALB_GENE):
        h_pct = {**h_peri_c, **h_peri_p, **get_mean_expression(h_hep, ALB_GENE)}.get(gene, {}).get("pct_nonzero", float("nan"))
        m_pct = {**m_peri_c, **m_peri_p, **get_mean_expression(m_hep, ALB_GENE)}.get(gene, {}).get("pct_nonzero", float("nan"))
        s_pct = {**s_peri_c, **s_peri_p, **get_mean_expression(s_hep, ALB_GENE)}.get(gene, {}).get("pct_nonzero", float("nan"))
        print(f"  {gene:<10} {h_pct:>12.1f} {m_pct:>12.1f} {s_pct:>12.1f}")

    task1_results = {
        "pericentral": {
            "tabula_human": {k: v for k, v in h_peri_c.items()},
            "tabula_mouse": {k: v for k, v in m_peri_c.items()},
            "sun2023": {k: v for k, v in s_peri_c.items()},
        },
        "periportal": {
            "tabula_human": {k: v for k, v in h_peri_p.items()},
            "tabula_mouse": {k: v for k, v in m_peri_p.items()},
            "sun2023": {k: v for k, v in s_peri_p.items()},
        },
        "aggregate_ratios": {
            "tabula_human_pc_pp": float(tab_h_ratio),
            "tabula_mouse_pc_pp": float(tab_m_ratio),
            "sun2023_pc_pp": float(sun_ratio),
        },
        "pericentral_shifted": pericentral_shifted,
    }

    # Save zonation marker table
    zon_rows = []
    for gene in list(PERICENTRAL_GENES) + list(PERIPORTAL_GENES) + list(ALB_GENE):
        zone = "pericentral" if gene in PERICENTRAL_GENES else ("periportal" if gene in PERIPORTAL_GENES else "alb")
        eid = {**PERICENTRAL_GENES, **PERIPORTAL_GENES, **ALB_GENE}[gene]
        all_expr = {
            "tabula_human": {**h_peri_c, **h_peri_p, **get_mean_expression(h_hep, ALB_GENE)},
            "tabula_mouse": {**m_peri_c, **m_peri_p, **get_mean_expression(m_hep, ALB_GENE)},
            "sun2023": {**s_peri_c, **s_peri_p, **get_mean_expression(s_hep, ALB_GENE)},
        }
        zon_rows.append({
            "gene": gene, "ensembl_id": eid, "zone": zone,
            "tabula_human_mean": all_expr["tabula_human"].get(gene, {}).get("mean"),
            "tabula_mouse_mean": all_expr["tabula_mouse"].get(gene, {}).get("mean"),
            "sun2023_mean": all_expr["sun2023"].get(gene, {}).get("mean"),
            "tabula_human_pct_detected": all_expr["tabula_human"].get(gene, {}).get("pct_nonzero"),
            "tabula_mouse_pct_detected": all_expr["tabula_mouse"].get(gene, {}).get("pct_nonzero"),
            "sun2023_pct_detected": all_expr["sun2023"].get(gene, {}).get("pct_nonzero"),
        })
    zon_df = pd.DataFrame(zon_rows)
    zon_df.to_csv(OUTPUT_DIR / "zonation_marker_expression.csv", index=False)

    # ==================================================================
    # TASK 2: Centroid shift in PCA space
    # ==================================================================
    print(f"\n{'='*72}")
    print("TASK 2: Centroid shift in PCA space")
    print(f"{'='*72}")

    # 2a: Reconstruct PCA model from primary centroids
    from sklearn.decomposition import PCA

    h_cent = pd.read_csv(CENTROIDS_HUMAN_PATH, index_col=0)
    m_cent = pd.read_csv(CENTROIDS_MOUSE_PATH, index_col=0)
    cell_types = sorted(h_cent.index.tolist())
    gene_list = list(h_cent.columns)

    combined = np.vstack([
        h_cent.loc[cell_types].values,
        m_cent.loc[cell_types].values,
    ])
    pca = PCA(n_components=PCA_VARIANCE_THRESHOLD, svd_solver="full", random_state=RANDOM_SEED)
    pca.fit(combined)
    n_components = pca.n_components_
    print(f"  PCA reconstructed: {n_components} components, "
          f"{sum(pca.explained_variance_ratio_)*100:.1f}% variance")

    # Verify against saved centroids
    saved = np.load(PCA_CENTROIDS_PATH)
    recon_pca = pca.transform(combined)
    max_diff = np.max(np.abs(recon_pca[:35] - saved["human"]))
    print(f"  PCA reconstruction max diff vs saved: {max_diff:.2e}")

    # Compute hepatocyte centroids in gene space
    gene_list_idx = {g: i for i, g in enumerate(gene_list)}

    def compute_centroid_gene_space(adata_subset):
        """Compute mean expression vector in the 16,959 gene space."""
        X = adata_subset.X
        if hasattr(X, 'toarray'):
            mean_vec = np.asarray(X.mean(axis=0)).flatten()
        else:
            mean_vec = np.mean(X, axis=0).flatten()
        return mean_vec

    # Hepatocyte centroids
    h_hep_centroid = compute_centroid_gene_space(h_hep)
    m_hep_centroid = compute_centroid_gene_space(m_hep)
    s_hep_centroid = compute_centroid_gene_space(s_hep)

    # 2a: Project into PCA space
    h_hep_pca = pca.transform(h_hep_centroid.reshape(1, -1))[0]
    m_hep_pca = pca.transform(m_hep_centroid.reshape(1, -1))[0]
    s_hep_pca = pca.transform(s_hep_centroid.reshape(1, -1))[0]

    print(f"\n  Hepatocyte centroids projected into {n_components}-D PCA space")

    # 2b: Pairwise Euclidean distances
    d_hm = np.linalg.norm(h_hep_pca - m_hep_pca)
    d_hs = np.linalg.norm(h_hep_pca - s_hep_pca)
    d_ms = np.linalg.norm(m_hep_pca - s_hep_pca)

    print(f"\n  Pairwise Euclidean distances in {n_components}-D PCA space:")
    print(f"    Tabula human vs Tabula mouse:   {d_hm:.3f}  (aggregate cross-species residual)")
    print(f"    Tabula human vs Sun2023 mouse:  {d_hs:.3f}  (T3-C tissue-matched distance → rank 28)")
    print(f"    Sun2023 mouse vs Tabula mouse:  {d_ms:.3f}  (how different are the two mouse centroids)")

    # 2c: Key diagnostic interpretation
    print(f"\n  DIAGNOSTIC:")
    print(f"    d(Sun2023 mouse, Tabula mouse) = {d_ms:.3f}")
    print(f"    d(Tabula human, Tabula mouse)  = {d_hm:.3f}")
    ratio_ms_hm = d_ms / d_hm
    print(f"    Ratio: {ratio_ms_hm:.3f}")
    if d_ms > d_hm:
        print(f"    *** Sun2023 mouse hepatocyte is FARTHER from Tabula mouse than Tabula human is! ***")
        print(f"    This strongly suggests an annotation artifact.")
    elif d_ms > d_hm * 0.5:
        print(f"    Sun2023 mouse is substantially different from Tabula mouse (>{d_hm*0.5:.1f}).")
        print(f"    Partial annotation artifact likely.")
    else:
        print(f"    Sun2023 mouse is close to Tabula mouse (ratio < 0.5).")
        print(f"    Rank 28 result appears genuine.")

    # Also compute distances for OTHER types to contextualize
    # How far are Sun2023 vs Tabula mouse centroids for well-replicating types?
    print(f"\n  Context: inter-mouse centroid distances for other T3-C types")

    # Load T3-C per-tissue distances for comparison
    t3c_df = pd.read_csv(T3C_RANKING_PATH)
    residuals_df = pd.read_csv(RESIDUALS_RANKED_PATH)

    print(f"    Hepatocyte d(Sun2023, Tabula mouse) = {d_ms:.3f}")
    print(f"    (For context: T3-C tissue-stratified hepatocyte score = "
          f"{t3c_df[t3c_df['cell_type']=='hepatocyte']['tissue_stratified_score'].values[0]:.3f})")
    print(f"    Aggregate hepatocyte residual = "
          f"{residuals_df[residuals_df['cell_type']=='hepatocyte']['residual_magnitude'].values[0]:.3f}")

    # 2d: Direction of Sun2023 vs Tabula mouse shift
    print(f"\n  Top 10 genes driving Sun2023 vs Tabula mouse hepatocyte centroid shift:")
    shift_vec = s_hep_centroid - m_hep_centroid
    abs_shift = np.abs(shift_vec)
    top_idx = np.argsort(abs_shift)[::-1][:10]

    shift_gene_rows = []
    for rank_i, idx in enumerate(top_idx):
        eid = gene_list[idx]
        gene_name = human_to_name.get(eid, eid)
        loading = shift_vec[idx]
        zone = ""
        if eid in PERICENTRAL_GENES.values():
            zone = "PERICENTRAL"
        elif eid in PERIPORTAL_GENES.values():
            zone = "PERIPORTAL"
        elif eid == ALB_GENE["ALB"]:
            zone = "ALB"
        shift_gene_rows.append({
            "rank": rank_i + 1,
            "gene_name": gene_name,
            "ensembl_id": eid,
            "shift_loading": float(loading),
            "abs_loading": float(abs_shift[idx]),
            "sun2023_mean": float(s_hep_centroid[idx]),
            "tabula_mouse_mean": float(m_hep_centroid[idx]),
            "zone_marker": zone,
        })
        flag = f" *** {zone} ***" if zone else ""
        direction = "↑ Sun2023" if loading > 0 else "↓ Sun2023"
        print(f"    {rank_i+1:>2}. {gene_name:<15} ({eid}) shift={loading:>+.4f}  "
              f"Sun={s_hep_centroid[idx]:.3f}  TabM={m_hep_centroid[idx]:.3f}  {direction}{flag}")

    shift_df = pd.DataFrame(shift_gene_rows)
    shift_df.to_csv(OUTPUT_DIR / "centroid_shift_top_genes.csv", index=False)

    # Check if any zonation markers are in the top 50
    print(f"\n  Zonation markers in top 50 shift genes:")
    top_50_idx = np.argsort(abs_shift)[::-1][:50]
    top_50_eids = set(gene_list[i] for i in top_50_idx)
    for name, eid in {**PERICENTRAL_GENES, **PERIPORTAL_GENES, **ALB_GENE}.items():
        if eid in top_50_eids:
            rank_in_shift = np.where(np.argsort(abs_shift)[::-1] == list(gene_list).index(eid))[0][0] + 1
            loading = shift_vec[list(gene_list).index(eid)]
            print(f"    {name}: rank {rank_in_shift}, loading={loading:+.4f}")
        else:
            print(f"    {name}: NOT in top 50")

    task2_results = {
        "distances": {
            "tabula_human_vs_tabula_mouse": float(d_hm),
            "tabula_human_vs_sun2023_mouse": float(d_hs),
            "sun2023_mouse_vs_tabula_mouse": float(d_ms),
        },
        "ratio_ms_to_hm": float(ratio_ms_hm),
        "n_pca_components": int(n_components),
    }

    # Save distance table
    dist_rows = [
        {"comparison": "Tabula human vs Tabula mouse", "distance": float(d_hm),
         "interpretation": "aggregate cross-species residual"},
        {"comparison": "Tabula human vs Sun2023 mouse", "distance": float(d_hs),
         "interpretation": "T3-C tissue-matched distance (rank 28)"},
        {"comparison": "Sun2023 mouse vs Tabula mouse", "distance": float(d_ms),
         "interpretation": "inter-mouse centroid divergence (artifact check)"},
    ]
    pd.DataFrame(dist_rows).to_csv(OUTPUT_DIR / "centroid_distances.csv", index=False)

    # ==================================================================
    # TASK 3: Alb-rescue bias check
    # ==================================================================
    print(f"\n{'='*72}")
    print("TASK 3: Alb-rescue bias check")
    print(f"{'='*72}")

    # Check how Sun2023 hepatocytes were annotated
    # No explicit flag exists — reconstruct from Leiden clusters
    s_hep_obs = s_hep.obs.copy()
    s_hep_leiden = s_hep_obs["leiden"].value_counts()
    print(f"\n  Sun2023 hepatocytes ({s_hep.n_obs} cells) by Leiden cluster:")
    for cluster, count in s_hep_leiden.items():
        print(f"    Cluster {cluster}: {count} cells")

    # The Alb-rescue logic: tissue=="liver" AND Alb>0 → hepatocyte
    # Check how many have Alb > 0
    alb_idx = gene_idx[ALB_GENE["ALB"]]
    s_alb_vals = s_hep.X[:, alb_idx]
    if hasattr(s_alb_vals, 'toarray'):
        s_alb_vals = s_alb_vals.toarray().flatten()
    else:
        s_alb_vals = np.asarray(s_alb_vals).flatten()
    n_alb_positive = (s_alb_vals > 0).sum()
    n_alb_zero = (s_alb_vals == 0).sum()

    print(f"\n  ALB expression in Sun2023 hepatocytes:")
    print(f"    ALB > 0: {n_alb_positive} cells ({n_alb_positive/s_hep.n_obs*100:.1f}%)")
    print(f"    ALB = 0: {n_alb_zero} cells ({n_alb_zero/s_hep.n_obs*100:.1f}%)")

    if n_alb_zero == 0:
        print(f"    *** ALL Sun2023 hepatocytes have ALB > 0 ***")
        print(f"    This means ALL cells could have come from Alb-rescue.")
        print(f"    Some may also have been annotated via initial marker scoring (Alb+Ttr+Apoa1),")
        print(f"    but there is no way to distinguish — both pathways require Alb detection.")
        all_alb_selected = True
    else:
        print(f"    {n_alb_zero} cells have ALB=0, indicating they came from marker scoring")
        print(f"    (Ttr/Apoa1 detection without Alb) or cluster-level annotation, not Alb-rescue.")
        all_alb_selected = False

    # Compare ALB expression levels
    print(f"\n  ALB expression comparison:")
    h_alb = h_hep.X[:, alb_idx]
    if hasattr(h_alb, 'toarray'):
        h_alb = h_alb.toarray().flatten()
    else:
        h_alb = np.asarray(h_alb).flatten()
    m_alb = m_hep.X[:, alb_idx]
    if hasattr(m_alb, 'toarray'):
        m_alb = m_alb.toarray().flatten()
    else:
        m_alb = np.asarray(m_alb).flatten()

    h_alb_mean = np.mean(h_alb)
    m_alb_mean = np.mean(m_alb)
    s_alb_mean = np.mean(s_alb_vals)

    h_alb_pct = np.mean(h_alb > 0) * 100
    m_alb_pct = np.mean(m_alb > 0) * 100
    s_alb_pct = np.mean(s_alb_vals > 0) * 100

    print(f"    {'Dataset':<18} {'ALB mean':>10} {'ALB % detected':>16}")
    print(f"    {'-'*46}")
    print(f"    {'Tabula Human':<18} {h_alb_mean:>10.4f} {h_alb_pct:>15.1f}%")
    print(f"    {'Tabula Mouse':<18} {m_alb_mean:>10.4f} {m_alb_pct:>15.1f}%")
    print(f"    {'Sun2023':<18} {s_alb_mean:>10.4f} {s_alb_pct:>15.1f}%")

    # Check >2x criterion
    alb_ratio_vs_h = s_alb_mean / h_alb_mean if h_alb_mean > 0 else float("inf")
    alb_ratio_vs_m = s_alb_mean / m_alb_mean if m_alb_mean > 0 else float("inf")
    print(f"\n    Sun2023 ALB / Tabula Human ALB: {alb_ratio_vs_h:.3f}x")
    print(f"    Sun2023 ALB / Tabula Mouse ALB: {alb_ratio_vs_m:.3f}x")

    alb_bias = alb_ratio_vs_m > 2.0
    if alb_bias:
        print(f"    *** Sun2023 ALB is >{alb_ratio_vs_m:.1f}x higher than Tabula Mouse — SELECTION BIAS ***")
    else:
        print(f"    Sun2023 ALB is within 2x of Tabula Mouse — no strong bias signal.")

    # Also check the cluster composition — how many clusters contribute to hepatocytes?
    # This helps distinguish "one big hepatocyte cluster" from "scattered Alb+ cells"
    n_clusters = s_hep_leiden.shape[0]
    largest_cluster = s_hep_leiden.iloc[0]
    largest_cluster_pct = largest_cluster / s_hep.n_obs * 100
    print(f"\n  Hepatocyte cluster structure:")
    print(f"    Number of Leiden clusters contributing: {n_clusters}")
    print(f"    Largest cluster: {largest_cluster} cells ({largest_cluster_pct:.1f}%)")
    if n_clusters > 5:
        print(f"    *** Hepatocytes come from {n_clusters} clusters — suggests scattered Alb+ rescue ***")
        scattered = True
    else:
        print(f"    Hepatocytes concentrated in {n_clusters} cluster(s) — suggests genuine cell type cluster")
        scattered = False

    task3_results = {
        "all_alb_selected": all_alb_selected,
        "n_alb_positive": int(n_alb_positive),
        "n_alb_zero": int(n_alb_zero),
        "alb_mean_tabula_human": float(h_alb_mean),
        "alb_mean_tabula_mouse": float(m_alb_mean),
        "alb_mean_sun2023": float(s_alb_mean),
        "alb_ratio_sun_vs_tab_mouse": float(alb_ratio_vs_m),
        "alb_bias_gt2x": alb_bias,
        "n_leiden_clusters": int(n_clusters),
        "scattered_across_clusters": scattered,
    }

    # Save Alb comparison
    alb_rows = [
        {"dataset": "Tabula Human", "alb_mean": float(h_alb_mean), "alb_pct_detected": float(h_alb_pct),
         "n_hepatocytes": int(h_hep.n_obs)},
        {"dataset": "Tabula Mouse", "alb_mean": float(m_alb_mean), "alb_pct_detected": float(m_alb_pct),
         "n_hepatocytes": int(m_hep.n_obs)},
        {"dataset": "Sun2023", "alb_mean": float(s_alb_mean), "alb_pct_detected": float(s_alb_pct),
         "n_hepatocytes": int(s_hep.n_obs)},
    ]
    pd.DataFrame(alb_rows).to_csv(OUTPUT_DIR / "alb_expression_comparison.csv", index=False)

    # ==================================================================
    # TASK 4: Verdict and documentation
    # ==================================================================
    print(f"\n{'='*72}")
    print("TASK 4: Verdict")
    print(f"{'='*72}")

    # Apply decision logic
    centroid_far = d_ms > d_hm * 0.5  # Sun2023 mouse is >50% of cross-species distance from Tabula mouse

    print(f"\n  Decision criteria:")
    print(f"    1. Sun2023 mouse centroid far from Tabula mouse? "
          f"d={d_ms:.3f} vs d_hm={d_hm:.3f} (ratio={ratio_ms_hm:.3f})")
    print(f"       → {'YES' if centroid_far else 'NO'} (threshold: ratio > 0.5)")
    print(f"    2. Pericentral marker shift confirmed? → {'YES' if pericentral_shifted else 'NO'}")
    print(f"    3. ALB mean > 2x in Sun2023 vs Tabula Mouse? → {'YES' if alb_bias else 'NO'} "
          f"(ratio={alb_ratio_vs_m:.3f})")
    print(f"    4. All cells Alb-selected? → {'YES' if all_alb_selected else 'NO'}")
    print(f"    5. Scattered across many clusters? → {'YES' if scattered else 'NO'}")

    # Determine verdict
    n_positive = sum([centroid_far, pericentral_shifted, alb_bias, all_alb_selected, scattered])

    if centroid_far and (pericentral_shifted or alb_bias):
        verdict = "ZONATION_ARTIFACT"
        explanation = (
            f"Sun2023 mouse hepatocyte centroid is far from Tabula mouse "
            f"(d={d_ms:.3f}, ratio={ratio_ms_hm:.3f} of cross-species distance). "
            f"{'Pericentral marker shift confirmed. ' if pericentral_shifted else ''}"
            f"{'ALB >2x higher in Sun2023. ' if alb_bias else ''}"
            f"The rank 28 result is NOT biological — it reflects Alb-rescue annotation bias "
            f"in Sun2023 creating a pericentral-shifted centroid."
        )
    elif not centroid_far:
        verdict = "GENUINE"
        explanation = (
            f"Sun2023 mouse hepatocyte centroid is close to Tabula mouse "
            f"(d={d_ms:.3f}, ratio={ratio_ms_hm:.3f} of cross-species distance). "
            f"The rank 28 result is real — hepatocyte is genuinely less rigid "
            f"than the aggregate ranking suggested."
        )
    else:
        verdict = "AMBIGUOUS"
        explanation = (
            f"Mixed evidence: centroid far={'YES' if centroid_far else 'NO'}, "
            f"pericentral shift={'YES' if pericentral_shifted else 'NO'}, "
            f"ALB bias={'YES' if alb_bias else 'NO'}, "
            f"all Alb-selected={'YES' if all_alb_selected else 'NO'}. "
            f"Positive criteria: {n_positive}/5. Flag to advisor with all numbers."
        )

    print(f"\n  *** VERDICT: {verdict} ***")
    print(f"  {explanation}")

    # Additional diagnostic: what drives the shift?
    print(f"\n  Additional analysis — shift decomposition:")
    # Check mitochondrial gene contribution
    mt_genes_in_shift = []
    non_mt_genes_in_shift = []
    for idx_i in range(len(gene_list)):
        gname = human_to_name.get(gene_list[idx_i], gene_list[idx_i])
        if isinstance(gname, str) and gname.startswith("MT-"):
            mt_genes_in_shift.append(abs(shift_vec[idx_i]))
        else:
            non_mt_genes_in_shift.append(abs(shift_vec[idx_i]))
    mt_l2 = np.sqrt(sum(x**2 for x in mt_genes_in_shift))
    total_l2 = np.linalg.norm(shift_vec)
    mt_pct_shift = (mt_l2 / total_l2) * 100 if total_l2 > 0 else 0
    print(f"    Mitochondrial genes contribute {mt_pct_shift:.1f}% of shift L2 norm")
    print(f"    (Sun2023 has MT expression; Tabula Mouse scaled data has MT=0)")

    # Check mean expression levels overall
    h_total_mean = np.mean(h_hep_centroid)
    m_total_mean = np.mean(m_hep_centroid)
    s_total_mean = np.mean(s_hep_centroid)
    print(f"\n    Global mean expression across 16,959 genes:")
    print(f"      Tabula Human: {h_total_mean:.4f}")
    print(f"      Tabula Mouse: {m_total_mean:.4f}")
    print(f"      Sun2023:      {s_total_mean:.4f}")
    depth_ratio = s_total_mean / m_total_mean if m_total_mean > 0 else 0
    print(f"      Sun2023 / Tabula Mouse ratio: {depth_ratio:.3f}")
    if depth_ratio < 0.5:
        print(f"      *** Sun2023 hepatocytes have {depth_ratio:.1%} of Tabula Mouse expression depth ***")
        print(f"      This is a GENE DETECTION DEPTH issue, not a zonation issue.")

    # Save final summary (convert numpy types)
    def to_native(obj):
        """Recursively convert numpy types to native Python for JSON."""
        if isinstance(obj, dict):
            return {k: to_native(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [to_native(v) for v in obj]
        elif isinstance(obj, (np.integer,)):
            return int(obj)
        elif isinstance(obj, (np.floating,)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.bool_):
            return bool(obj)
        return obj

    final_summary = to_native({
        "task": "Hepatocyte Zonation Diagnostic",
        "task0": task0_results,
        "task1": task1_results,
        "task2": task2_results,
        "task3": task3_results,
        "additional": {
            "mt_pct_of_shift": mt_pct_shift,
            "global_mean_tabula_human": h_total_mean,
            "global_mean_tabula_mouse": m_total_mean,
            "global_mean_sun2023": s_total_mean,
            "depth_ratio_sun_vs_tabm": depth_ratio,
        },
        "verdict": verdict,
        "explanation": explanation,
        "decision_criteria": {
            "centroid_far": centroid_far,
            "pericentral_shifted": pericentral_shifted,
            "alb_bias_gt2x": alb_bias,
            "all_alb_selected": all_alb_selected,
            "scattered_clusters": scattered,
            "n_positive_of_5": n_positive,
        },
    })

    with open(OUTPUT_DIR / "diagnostic_summary.json", "w") as f:
        json.dump(final_summary, f, indent=2)

    print(f"\n  All outputs saved to {OUTPUT_DIR}/")
    print(f"  - zonation_marker_expression.csv")
    print(f"  - centroid_distances.csv")
    print(f"  - centroid_shift_top_genes.csv")
    print(f"  - alb_expression_comparison.csv")
    print(f"  - diagnostic_summary.json")

    print(f"\n{'='*72}")
    print("HEPATOCYTE ZONATION DIAGNOSTIC COMPLETE")
    print(f"{'='*72}")

    return final_summary


if __name__ == "__main__":
    main()
