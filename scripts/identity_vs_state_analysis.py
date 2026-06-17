#!/usr/bin/env python3
"""
Identity vs Activation State Analysis
======================================
Tests whether cancer and COVID-19 Procrustes deformation differences
reflect cell-type identity transformation vs immune activation state.

Biology: Cancer deformation (normal→tumor) might reshape the geometric
identity of cell types (changing *what* a cell is), while COVID-19
deformation (normal→infected) might primarily shift cells along
activation axes (changing *how active* a cell is) without changing
identity. If true, projecting out activation-related genes should
collapse the COVID signal more than the cancer signal.

Math: We reconstruct gene-space deformation vectors from stored
Procrustes residuals in PCA space, partition genes into activation vs
identity sets, and compare enrichment ratios and residual cross-axis
correlations.
"""

import json
import os
import sys

import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# ─── Paths ──────────────────────────────────────────────────────────
OUT_DIR = "output/mechanistic/identity_vs_state"
CANCER_PCA = "output/cancer/scaled/pca_cancer_scaled.npz"
CANCER_JSON = "output/cancer/scaled/cancer_scaled_results.json"
COVID_PCA = "output/disease_replication/covid/pca_covid.npz"
COVID_JSON = "output/disease_replication/covid/covid_procrustes_results.json"
CANCER_CENTROIDS = "output/cancer/scaled/centroids_normal_scaled.csv"
COVID_CENTROIDS = "output/disease_replication/covid/centroids_normal.csv"
ORTHOLOG_MAP = "data/phase1/orthologs_human_mouse.csv"

os.makedirs(OUT_DIR, exist_ok=True)


# ═══════════════════════════════════════════════════════════════════
# STEP 1: Reconstruct gene-space deformation vectors from stored
#         Procrustes residuals
# ═══════════════════════════════════════════════════════════════════

def reconstruct_gene_vectors(pca_path, results_json_path):
    """Reconstruct full gene-space deformation vectors from stored Procrustes residuals.

    The results JSON stores per-cell-type residual vectors in PCA space
    (aligned_target - centered_reference after Procrustes rotation/scaling).
    We project these back to gene space via residual_pca @ pca_components.

    Returns: (gene_vectors, cell_types, pca_components)
        gene_vectors: ndarray shape (n_cell_types, n_genes)
        cell_types: list of cell type names
        pca_components: ndarray shape (n_components, n_genes)
    """
    with open(results_json_path) as f:
        results = json.load(f)

    pca_data = np.load(pca_path, allow_pickle=True)
    components = pca_data["components"]  # (n_comp, n_genes)

    cell_types = results["cell_types"]
    residuals = results["residuals"]

    # Build matrix of PCA-space residual vectors
    n_types = len(cell_types)
    n_comp = components.shape[0]
    vector_pca = np.zeros((n_types, n_comp))

    for i, ct in enumerate(cell_types):
        vector_pca[i] = np.array(residuals[ct]["vector_pca"])

    # Project to gene space
    vector_gene = vector_pca @ components  # (n_types, n_genes)
    return vector_gene, cell_types, components


def verify_top_genes(vector_gene, gene_symbols, cell_types, results_json_path, label):
    """Verify reconstructed vectors match stored top-3 genes per cell type."""
    with open(results_json_path) as f:
        results = json.load(f)
    stored = results["top_genes_per_cell_type"]

    mismatches = []
    for i, ct in enumerate(cell_types):
        if ct not in stored:
            continue
        abs_loadings = np.abs(vector_gene[i])
        top_idx = np.argsort(abs_loadings)[::-1][:3]
        reconstructed_top3 = [gene_symbols[j] for j in top_idx]
        stored_top3 = [g["gene"] for g in stored[ct][:3]]

        if reconstructed_top3 != stored_top3:
            mismatches.append({
                "cell_type": ct,
                "reconstructed": reconstructed_top3,
                "stored": stored_top3,
            })

    if mismatches:
        print(f"\n  WARNING {label}: {len(mismatches)} cell types have top-3 gene mismatches:")
        for m in mismatches:
            print(f"    {m['cell_type']}: reconstructed={m['reconstructed']}, stored={m['stored']}")
            overlap = len(set(m["reconstructed"]) & set(m["stored"]))
            if overlap < 2:
                print(f"    CRITICAL: <2 overlapping top-3 genes — STOPPING")
                sys.exit(1)
        print("  All mismatches are minor rank swaps — proceeding.")
    else:
        print(f"  {label}: All top-3 genes verified.")


# Load gene name mapping
centroids_cols = pd.read_csv(CANCER_CENTROIDS, nrows=0).columns.tolist()
gene_ensg = centroids_cols[1:]
orth = pd.read_csv(ORTHOLOG_MAP)
ensg2sym = dict(zip(orth["human_ensembl_id"], orth["human_gene_name"]))
gene_symbols = [ensg2sym[g] for g in gene_ensg]
n_genes = len(gene_symbols)

print(f"Gene space: {n_genes} genes")

# Load results JSONs
with open(CANCER_JSON) as f:
    cancer_results = json.load(f)
with open(COVID_JSON) as f:
    covid_results = json.load(f)

# Reconstruct gene-space vectors
print("\n=== Step 1: Reconstructing gene-space deformation vectors ===")
cancer_vec, cancer_types, cancer_comp = reconstruct_gene_vectors(CANCER_PCA, CANCER_JSON)
covid_vec, covid_types, covid_comp = reconstruct_gene_vectors(COVID_PCA, COVID_JSON)

print(f"  Cancer: {cancer_vec.shape} ({len(cancer_types)} types x {n_genes} genes)")
print(f"  COVID:  {covid_vec.shape} ({len(covid_types)} types x {n_genes} genes)")

# Verify reconstruction
verify_top_genes(cancer_vec, gene_symbols, cancer_types, CANCER_JSON, "Cancer")
verify_top_genes(covid_vec, gene_symbols, covid_types, COVID_JSON, "COVID")


# ═══════════════════════════════════════════════════════════════════
# STEP 2: Download MSigDB Hallmark gene sets & define activation set
# ═══════════════════════════════════════════════════════════════════

print("\n=== Step 2: Downloading MSigDB Hallmark gene sets ===")

ACTIVATION_SETS = [
    "Interferon Alpha Response",
    "Interferon Gamma Response",
    "TNF-alpha Signaling via NF-kB",
    "Inflammatory Response",
    "IL-6/JAK/STAT3 Signaling",
]

hallmark = gp.get_library("MSigDB_Hallmark_2020", organism="Human")

activation_genes_raw = set()
for gs_name in ACTIVATION_SETS:
    if gs_name in hallmark:
        genes = hallmark[gs_name]
        activation_genes_raw.update(genes)
        print(f"  {gs_name}: {len(genes)} genes")
    else:
        for k in hallmark:
            if gs_name in k:
                genes = hallmark[k]
                activation_genes_raw.update(genes)
                print(f"  {gs_name} (matched as {k}): {len(genes)} genes")
                break
        else:
            print(f"  WARNING: {gs_name} not found in Hallmark library")

print(f"\n  Total unique activation genes (raw): {len(activation_genes_raw)}")

# Map to our gene space
gene_set = set(gene_symbols)
activation_genes_mapped = activation_genes_raw & gene_set
activation_idx = np.array([i for i, g in enumerate(gene_symbols) if g in activation_genes_mapped])
non_activation_idx = np.array([i for i in range(n_genes) if i not in set(activation_idx)])

print(f"  Mapped to our {n_genes}-gene space: {len(activation_genes_mapped)} activation genes")
print(f"  Non-activation genes: {len(non_activation_idx)}")

# Binary indicator
activation_indicator = np.zeros(n_genes)
activation_indicator[activation_idx] = 1.0


# ═══════════════════════════════════════════════════════════════════
# STEP 3: Compute activation enrichment per cell type per disease
# ═══════════════════════════════════════════════════════════════════

print("\n=== Step 3: Activation enrichment scores ===")


def compute_enrichment(deform_vec, target_idx, cell_types):
    """Enrichment = mean |target genes| / mean |all genes|."""
    results = {}
    for i, ct in enumerate(cell_types):
        abs_load = np.abs(deform_vec[i])
        mean_target = abs_load[target_idx].mean()
        mean_all = abs_load.mean()
        ratio = mean_target / mean_all
        results[ct] = {
            "mean_target": float(mean_target),
            "mean_all": float(mean_all),
            "enrichment_ratio": float(ratio),
        }
    return results


cancer_enrichment = compute_enrichment(cancer_vec, activation_idx, cancer_types)
covid_enrichment = compute_enrichment(covid_vec, activation_idx, covid_types)

print("\n  Cancer activation enrichment per cell type:")
for ct, v in cancer_enrichment.items():
    print(f"    {ct:45s}: {v['enrichment_ratio']:.3f}")

print("\n  COVID activation enrichment per cell type:")
for ct, v in covid_enrichment.items():
    print(f"    {ct:45s}: {v['enrichment_ratio']:.3f}")

cancer_ratios = [v["enrichment_ratio"] for v in cancer_enrichment.values()]
covid_ratios = [v["enrichment_ratio"] for v in covid_enrichment.values()]

mw_stat, mw_p = stats.mannwhitneyu(cancer_ratios, covid_ratios, alternative="two-sided")
print(f"\n  Mann-Whitney U: statistic={mw_stat:.1f}, p={mw_p:.4f}")
print(f"  Cancer median enrichment: {np.median(cancer_ratios):.3f}")
print(f"  COVID median enrichment:  {np.median(covid_ratios):.3f}")


# ═══════════════════════════════════════════════════════════════════
# STEP 4: Project out activation component, recompute cross-axis ρ
# ═══════════════════════════════════════════════════════════════════

print("\n=== Step 4: Project out activation component ===")


def project_out_activation(deform_vec, indicator):
    """Remove activation-correlated magnitude from each deformation vector.

    For each cell type: regress |loading_i| ~ beta0 + beta1 * indicator_i.
    Residual = original loading with magnitude reduced by beta1 for activation genes.
    """
    residual_vec = np.zeros_like(deform_vec)
    for i in range(deform_vec.shape[0]):
        abs_load = np.abs(deform_vec[i])
        X = np.column_stack([np.ones(len(indicator)), indicator])
        beta = np.linalg.lstsq(X, abs_load, rcond=None)[0]
        predicted_activation = beta[1] * indicator
        sign = np.sign(deform_vec[i])
        residual_abs = np.maximum(abs_load - predicted_activation, 0)
        residual_vec[i] = sign * residual_abs
    return residual_vec


cancer_residual = project_out_activation(cancer_vec, activation_indicator)
covid_residual = project_out_activation(covid_vec, activation_indicator)

# Deformation scores (L2 norm)
cancer_scores_orig = np.linalg.norm(cancer_vec, axis=1)
cancer_scores_resid = np.linalg.norm(cancer_residual, axis=1)
covid_scores_orig = np.linalg.norm(covid_vec, axis=1)
covid_scores_resid = np.linalg.norm(covid_residual, axis=1)

print("\n  Cancer deformation scores (original -> residual):")
for i, ct in enumerate(cancer_types):
    pct = (1 - cancer_scores_resid[i] / cancer_scores_orig[i]) * 100
    print(f"    {ct:45s}: {cancer_scores_orig[i]:.3f} -> {cancer_scores_resid[i]:.3f}  ({pct:+.1f}%)")

print("\n  COVID deformation scores (original -> residual):")
for i, ct in enumerate(covid_types):
    pct = (1 - covid_scores_resid[i] / covid_scores_orig[i]) * 100
    print(f"    {ct:45s}: {covid_scores_orig[i]:.3f} -> {covid_scores_resid[i]:.3f}  ({pct:+.1f}%)")

# Cross-axis Spearman: disease deformation scores vs cross-species residuals
cancer_xs_data = cancer_results["cross_analysis_correlation"]["matched_types"]
cancer_xs_df = pd.DataFrame(cancer_xs_data)
covid_xs_data = covid_results["cross_analysis_correlation"]["matched_types"]
covid_xs_df = pd.DataFrame(covid_xs_data)

# Original ρ values
cancer_orig_rho = cancer_results["cross_analysis_correlation"]["spearman_rho"]
cancer_orig_p = cancer_results["cross_analysis_correlation"]["spearman_p"]
covid_orig_rho = covid_results["cross_analysis_correlation"]["spearman_rho"]
covid_orig_p = covid_results["cross_analysis_correlation"]["spearman_p"]

# Map residual scores to cell types
cancer_type_to_resid = {ct: float(cancer_scores_resid[i]) for i, ct in enumerate(cancer_types)}
covid_type_to_resid = {ct: float(covid_scores_resid[i]) for i, ct in enumerate(covid_types)}

# Cancer residual cross-axis
cancer_key = "cancer_type"
cancer_xs_df["resid_deformation"] = cancer_xs_df[cancer_key].map(cancer_type_to_resid)
cancer_valid = cancer_xs_df.dropna(subset=["resid_deformation"])
cancer_resid_rho, cancer_resid_p = stats.spearmanr(
    cancer_valid["resid_deformation"], cancer_valid["xs_residual"]
)

# COVID residual cross-axis
covid_key = "covid_type"
covid_xs_df["resid_deformation"] = covid_xs_df[covid_key].map(covid_type_to_resid)
covid_valid = covid_xs_df.dropna(subset=["resid_deformation"])
covid_resid_rho, covid_resid_p = stats.spearmanr(
    covid_valid["resid_deformation"], covid_valid["xs_residual"]
)

print(f"\n  Cross-axis Spearman rho (original -> residual):")
print(f"    Cancer: rho={cancer_orig_rho:.4f} (p={cancer_orig_p:.4f}) -> rho={cancer_resid_rho:.4f} (p={cancer_resid_p:.4f})")
print(f"    COVID:  rho={covid_orig_rho:.4f} (p={covid_orig_p:.4f}) -> rho={covid_resid_rho:.4f} (p={covid_resid_p:.4f})")

cancer_pct_reduction = np.mean(1 - cancer_scores_resid / cancer_scores_orig) * 100
covid_pct_reduction = np.mean(1 - covid_scores_resid / covid_scores_orig) * 100
print(f"\n  Mean magnitude reduction:")
print(f"    Cancer: {cancer_pct_reduction:.1f}%")
print(f"    COVID:  {covid_pct_reduction:.1f}%")


# ═══════════════════════════════════════════════════════════════════
# STEP 5: Identity proxy — cell-type-specific genes
# ═══════════════════════════════════════════════════════════════════

print("\n=== Step 5: Identity proxy — cell-type-specific genes ===")

# Load normal centroids and compute cross-cell-type variance
cancer_normal = pd.read_csv(CANCER_CENTROIDS, index_col=0)
covid_normal = pd.read_csv(COVID_CENTROIDS, index_col=0)
all_normal = pd.concat([cancer_normal, covid_normal], axis=0)

gene_variance = all_normal.var(axis=0).values
identity_k = 500
identity_idx = np.argsort(gene_variance)[::-1][:identity_k]
identity_genes = [gene_symbols[i] for i in identity_idx]

print(f"  Top {identity_k} identity genes (by cross-cell-type variance):")
print(f"  Examples: {identity_genes[:20]}")

overlap = set(identity_genes) & activation_genes_mapped
print(f"  Overlap with activation genes: {len(overlap)}")

identity_idx_arr = np.array(identity_idx)

cancer_identity = compute_enrichment(cancer_vec, identity_idx_arr, cancer_types)
covid_identity = compute_enrichment(covid_vec, identity_idx_arr, covid_types)

print("\n  Cancer identity enrichment per cell type:")
for ct, v in cancer_identity.items():
    print(f"    {ct:45s}: {v['enrichment_ratio']:.3f}")

print("\n  COVID identity enrichment per cell type:")
for ct, v in covid_identity.items():
    print(f"    {ct:45s}: {v['enrichment_ratio']:.3f}")

cancer_id_ratios = [v["enrichment_ratio"] for v in cancer_identity.values()]
covid_id_ratios = [v["enrichment_ratio"] for v in covid_identity.values()]

mw_id_stat, mw_id_p = stats.mannwhitneyu(cancer_id_ratios, covid_id_ratios, alternative="two-sided")
print(f"\n  Mann-Whitney U (identity): statistic={mw_id_stat:.1f}, p={mw_id_p:.4f}")
print(f"  Cancer median identity enrichment: {np.median(cancer_id_ratios):.3f}")
print(f"  COVID median identity enrichment:  {np.median(covid_id_ratios):.3f}")


# ═══════════════════════════════════════════════════════════════════
# FIGURE 1: Side-by-side bar chart of activation enrichment
# ═══════════════════════════════════════════════════════════════════

print("\n=== Generating figures ===")

fig, axes = plt.subplots(1, 2, figsize=(16, 7), sharey=True)

cancer_df = pd.DataFrame([
    {"cell_type": ct, "enrichment": v["enrichment_ratio"]}
    for ct, v in cancer_enrichment.items()
]).sort_values("enrichment", ascending=True)

axes[0].barh(range(len(cancer_df)), cancer_df["enrichment"].values, color="#E74C3C", alpha=0.8)
axes[0].set_yticks(range(len(cancer_df)))
axes[0].set_yticklabels(cancer_df["cell_type"].values, fontsize=9)
axes[0].axvline(1.0, color="gray", linestyle="--", alpha=0.5, label="No enrichment")
axes[0].set_xlabel("Activation Enrichment Ratio")
axes[0].set_title(f"Cancer (n={len(cancer_types)} types)\nMedian={np.median(cancer_ratios):.2f}")
axes[0].legend(fontsize=8)

covid_df = pd.DataFrame([
    {"cell_type": ct, "enrichment": v["enrichment_ratio"]}
    for ct, v in covid_enrichment.items()
]).sort_values("enrichment", ascending=True)

axes[1].barh(range(len(covid_df)), covid_df["enrichment"].values, color="#3498DB", alpha=0.8)
axes[1].set_yticks(range(len(covid_df)))
axes[1].set_yticklabels(covid_df["cell_type"].values, fontsize=9)
axes[1].axvline(1.0, color="gray", linestyle="--", alpha=0.5, label="No enrichment")
axes[1].set_xlabel("Activation Enrichment Ratio")
axes[1].set_title(f"COVID-19 (n={len(covid_types)} types)\nMedian={np.median(covid_ratios):.2f}")
axes[1].legend(fontsize=8)

fig.suptitle(
    f"Activation Gene Enrichment in Disease Deformation Vectors\n"
    f"Mann-Whitney U p={mw_p:.4f}  |  {len(activation_genes_mapped)} activation genes from 5 Hallmark sets",
    fontsize=12, fontweight="bold",
)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/activation_enrichment_bars.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Saved: activation_enrichment_bars.png")
print("    Horizontal bar chart comparing activation gene enrichment ratios")
print("    per cell type for cancer (red) vs COVID (blue). Dashed line at 1.0 = no enrichment.")


# ═══════════════════════════════════════════════════════════════════
# FIGURE 2: Scatter plot — original vs residual deformation scores
# ═══════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(9, 8))

ax.scatter(cancer_scores_orig, cancer_scores_resid,
           c="#E74C3C", s=80, alpha=0.8, edgecolors="white", linewidth=0.5,
           label=f"Cancer (n={len(cancer_types)})", zorder=3)
for i, ct in enumerate(cancer_types):
    short = ct[:20] + "..." if len(ct) > 20 else ct
    ax.annotate(short, (cancer_scores_orig[i], cancer_scores_resid[i]),
                fontsize=6, alpha=0.7, xytext=(3, 3), textcoords="offset points")

ax.scatter(covid_scores_orig, covid_scores_resid,
           c="#3498DB", s=80, alpha=0.8, edgecolors="white", linewidth=0.5,
           label=f"COVID (n={len(covid_types)})", zorder=3)
for i, ct in enumerate(covid_types):
    short = ct[:20] + "..." if len(ct) > 20 else ct
    ax.annotate(short, (covid_scores_orig[i], covid_scores_resid[i]),
                fontsize=6, alpha=0.7, xytext=(3, 3), textcoords="offset points")

max_val = max(cancer_scores_orig.max(), covid_scores_orig.max(),
              cancer_scores_resid.max(), covid_scores_resid.max()) * 1.1
ax.plot([0, max_val], [0, max_val], "k--", alpha=0.3, label="y=x (no reduction)")
ax.set_xlabel("Original Deformation Score (L2 norm)")
ax.set_ylabel("Residual Deformation Score (activation projected out)")
ax.set_title(
    f"Deformation Before vs After Activation Projection\n"
    f"Cancer mean reduction: {cancer_pct_reduction:.1f}%  |  "
    f"COVID mean reduction: {covid_pct_reduction:.1f}%"
)
ax.legend()
ax.set_xlim(0, max_val)
ax.set_ylim(0, max_val)
ax.set_aspect("equal")
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/original_vs_residual_scatter.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Saved: original_vs_residual_scatter.png")
print("    Scatter: each point is a cell type. X=original, Y=residual magnitude.")
print("    Points below y=x lost deformation magnitude from activation projection.")


# ═══════════════════════════════════════════════════════════════════
# FIGURE 3: Summary table
# ═══════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(10, 3.5))
ax.axis("off")

table_data = [
    ["", "Cancer", "Cancer", "COVID-19", "COVID-19"],
    ["", "Original", "Residual", "Original", "Residual"],
    ["Cross-axis Spearman rho",
     f"{cancer_orig_rho:.4f}", f"{cancer_resid_rho:.4f}",
     f"{covid_orig_rho:.4f}", f"{covid_resid_rho:.4f}"],
    ["p-value",
     f"{cancer_orig_p:.4f}", f"{cancer_resid_p:.4f}",
     f"{covid_orig_p:.4f}", f"{covid_resid_p:.4f}"],
    ["Mean deform. reduction", "—", f"{cancer_pct_reduction:.1f}%",
     "—", f"{covid_pct_reduction:.1f}%"],
    ["Activation enrichment (med.)",
     f"{np.median(cancer_ratios):.3f}", "—",
     f"{np.median(covid_ratios):.3f}", "—"],
    ["Identity enrichment (med.)",
     f"{np.median(cancer_id_ratios):.3f}", "—",
     f"{np.median(covid_id_ratios):.3f}", "—"],
    [f"Activation MW-U p", f"{mw_p:.4f}", "", "", ""],
    [f"Identity MW-U p", f"{mw_id_p:.4f}", "", "", ""],
]

table = ax.table(cellText=table_data, loc="center", cellLoc="center")
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1.2, 1.5)

for j in range(5):
    table[0, j].set_facecolor("#E8E8E8")
    table[0, j].set_text_props(fontweight="bold")
    table[1, j].set_facecolor("#F0F0F0")
    table[1, j].set_text_props(fontweight="bold")

ax.set_title("Identity vs Activation State: Summary Statistics",
             fontsize=12, fontweight="bold", pad=20)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/summary_table.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Saved: summary_table.png")


# ═══════════════════════════════════════════════════════════════════
# Save numerical results
# ═══════════════════════════════════════════════════════════════════

results = {
    "activation_genes": {
        "n_hallmark_sets": len(ACTIVATION_SETS),
        "sets_used": ACTIVATION_SETS,
        "n_raw": len(activation_genes_raw),
        "n_mapped": len(activation_genes_mapped),
        "n_total_genes": n_genes,
    },
    "identity_genes": {
        "n_identity": identity_k,
        "method": "top 500 genes by variance across cell type centroids",
        "overlap_with_activation": len(overlap),
        "examples": identity_genes[:20],
    },
    "activation_enrichment": {
        "cancer": cancer_enrichment,
        "covid": covid_enrichment,
        "mann_whitney_U": float(mw_stat),
        "mann_whitney_p": float(mw_p),
        "cancer_median": float(np.median(cancer_ratios)),
        "covid_median": float(np.median(covid_ratios)),
    },
    "identity_enrichment": {
        "cancer": {ct: v["enrichment_ratio"] for ct, v in cancer_identity.items()},
        "covid": {ct: v["enrichment_ratio"] for ct, v in covid_identity.items()},
        "mann_whitney_U": float(mw_id_stat),
        "mann_whitney_p": float(mw_id_p),
        "cancer_median": float(np.median(cancer_id_ratios)),
        "covid_median": float(np.median(covid_id_ratios)),
    },
    "cross_axis_spearman": {
        "cancer_original": {"rho": float(cancer_orig_rho), "p": float(cancer_orig_p)},
        "cancer_residual": {"rho": float(cancer_resid_rho), "p": float(cancer_resid_p)},
        "covid_original": {"rho": float(covid_orig_rho), "p": float(covid_orig_p)},
        "covid_residual": {"rho": float(covid_resid_rho), "p": float(covid_resid_p)},
    },
    "deformation_magnitude": {
        "cancer_mean_reduction_pct": float(cancer_pct_reduction),
        "covid_mean_reduction_pct": float(covid_pct_reduction),
        "cancer_per_type": {ct: {"original": float(cancer_scores_orig[i]),
                                  "residual": float(cancer_scores_resid[i])}
                            for i, ct in enumerate(cancer_types)},
        "covid_per_type": {ct: {"original": float(covid_scores_orig[i]),
                                 "residual": float(covid_scores_resid[i])}
                           for i, ct in enumerate(covid_types)},
    },
}

with open(f"{OUT_DIR}/identity_vs_state_results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"\n  Saved: {OUT_DIR}/identity_vs_state_results.json")

summary_df = pd.DataFrame({
    "Disease": ["Cancer", "Cancer", "COVID-19", "COVID-19"],
    "Condition": ["Original", "Activation projected out", "Original", "Activation projected out"],
    "Spearman_rho": [cancer_orig_rho, cancer_resid_rho, covid_orig_rho, covid_resid_rho],
    "p_value": [cancer_orig_p, cancer_resid_p, covid_orig_p, covid_resid_p],
})
summary_df.to_csv(f"{OUT_DIR}/spearman_summary.csv", index=False)
print(f"  Saved: {OUT_DIR}/spearman_summary.csv")


# ═══════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ═══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("IDENTITY vs ACTIVATION STATE ANALYSIS — FINAL SUMMARY")
print("=" * 70)

print(f"""
ACTIVATION ENRICHMENT:
  Cancer median enrichment: {np.median(cancer_ratios):.3f}
  COVID  median enrichment: {np.median(covid_ratios):.3f}
  Mann-Whitney U p = {mw_p:.4f}
  Interpretation: {'COVID significantly more activation-enriched' if mw_p < 0.05 and np.median(covid_ratios) > np.median(cancer_ratios) else 'Cancer significantly more activation-enriched' if mw_p < 0.05 and np.median(cancer_ratios) > np.median(covid_ratios) else 'No significant difference in activation enrichment'}

IDENTITY ENRICHMENT:
  Cancer median enrichment: {np.median(cancer_id_ratios):.3f}
  COVID  median enrichment: {np.median(covid_id_ratios):.3f}
  Mann-Whitney U p = {mw_id_p:.4f}
  Interpretation: {'Cancer significantly more identity-enriched' if mw_id_p < 0.05 and np.median(cancer_id_ratios) > np.median(covid_id_ratios) else 'COVID significantly more identity-enriched' if mw_id_p < 0.05 and np.median(covid_id_ratios) > np.median(cancer_id_ratios) else 'No significant difference in identity enrichment'}

CROSS-AXIS SPEARMAN (disease deformation <-> cross-species residual):
  Cancer: rho={cancer_orig_rho:.4f} (p={cancer_orig_p:.4f}) -> rho={cancer_resid_rho:.4f} (p={cancer_resid_p:.4f}) after projection
  COVID:  rho={covid_orig_rho:.4f} (p={covid_orig_p:.4f}) -> rho={covid_resid_rho:.4f} (p={covid_resid_p:.4f}) after projection

DEFORMATION MAGNITUDE REDUCTION:
  Cancer mean: {cancer_pct_reduction:.1f}%
  COVID  mean: {covid_pct_reduction:.1f}%

PREDICTIONS (identity vs activation hypothesis):
  (a) COVID activation enrichment > cancer:  {'YES' if np.median(covid_ratios) > np.median(cancer_ratios) else 'NO'}
  (b) COVID loses more magnitude on projection: {'YES' if covid_pct_reduction > cancer_pct_reduction else 'NO'}
  (c) Cancer has higher identity enrichment:  {'YES' if np.median(cancer_id_ratios) > np.median(covid_id_ratios) else 'NO'}
""")
