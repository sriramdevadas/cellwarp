#!/usr/bin/env python3
"""
Ribosomal confound test for sequence-expression anticorrelation.

Biology: The negative correlation between protein sequence conservation and
expression divergence (rho = -0.234) could be inflated by ribosomal proteins,
which are both extremely sequence-conserved AND highly expressed (thus having
large Procrustes loading magnitudes simply due to scale). This script tests
whether the anticorrelation survives after removing ribosomal and other
high-expression, sequence-constrained gene families.

Math: We recompute Spearman correlations between protein percent identity and
mean expression divergence (Procrustes loading magnitude) under four conditions:
  1. All 16,324 genes (original result)
  2. Excluding ribosomal protein genes (RPL/RPS/MRPL/MRPS)
  3. Excluding ribosomal + histones + heat shock proteins
  4. Expression-level normalization: divergence / mean_expression_level

Output: output/phase2/dnds/confound_analysis/
"""

import numpy as np
import pandas as pd
from scipy import stats
import json
import re

# ── Load data ────────────────────────────────────────────────────────────────

# Merged sequence identity + expression divergence (16,324 genes with both)
merged = pd.read_csv("output/phase2/dnds/merged_seq_expr_divergence.csv")

# Gene loadings (has gene names)
loadings = pd.read_csv("output/phase2/dnds/gene_loadings_all.csv")
gene_names = loadings[["ensembl_gene_id", "gene_name"]].drop_duplicates()

# Merge gene names into the analysis set
df = merged.merge(gene_names, left_on="human_ensembl_id", right_on="ensembl_gene_id",
                  how="left")
print(f"Analysis set: {len(df)} genes with both sequence identity and expression divergence")
print(f"  Genes with names: {df['gene_name'].notna().sum()}")

# ── Define gene family filters ───────────────────────────────────────────────

def is_ribosomal(name):
    """Ribosomal protein genes: RPL*, RPS*, MRPL*, MRPS*, RPP*."""
    if pd.isna(name):
        return False
    return bool(re.match(r'^(RPL|RPS|MRPL|MRPS|RPP)\d', name))

def is_histone(name):
    """Histone genes: HIST*, H1*, H2*, H3*, H4*."""
    if pd.isna(name):
        return False
    return bool(re.match(r'^(HIST|H1[A-Z]|H2[A-Z]|H3[A-Z]|H4[A-Z])', name))

def is_hsp(name):
    """Heat shock proteins: HSP*, HSPA*, HSPB*, HSPC*, HSPD*, HSPE*, HSPH*."""
    if pd.isna(name):
        return False
    return bool(re.match(r'^HSP[A-Z0-9]', name))

# Flag each gene
df["is_ribosomal"] = df["gene_name"].apply(is_ribosomal)
df["is_histone"] = df["gene_name"].apply(is_histone)
df["is_hsp"] = df["gene_name"].apply(is_hsp)
df["is_any_confound"] = df["is_ribosomal"] | df["is_histone"] | df["is_hsp"]

# ── STEP 1: Characterize ribosomal genes ─────────────────────────────────────

ribo = df[df["is_ribosomal"]]
non_ribo = df[~df["is_ribosomal"]]

print("\n" + "="*70)
print("STEP 1: Ribosomal gene characterization")
print("="*70)
print(f"\nRibosomal genes in analysis: {len(ribo)} / {len(df)} "
      f"({100*len(ribo)/len(df):.1f}%)")

# Breakdown by prefix
for prefix in ["RPL", "RPS", "MRPL", "MRPS", "RPP"]:
    n = sum(1 for name in ribo["gene_name"] if name.startswith(prefix))
    if n > 0:
        print(f"  {prefix}*: {n}")

print(f"\nMean protein identity:")
print(f"  Ribosomal:     {ribo['pct_identity'].mean():.2f}% "
      f"(SD {ribo['pct_identity'].std():.2f})")
print(f"  Non-ribosomal: {non_ribo['pct_identity'].mean():.2f}% "
      f"(SD {non_ribo['pct_identity'].std():.2f})")

# Mann-Whitney U test for identity difference
u_stat, u_p = stats.mannwhitneyu(ribo["pct_identity"], non_ribo["pct_identity"],
                                  alternative="greater")
print(f"  Mann-Whitney U (ribo > non-ribo): p = {u_p:.2e}")

print(f"\nMean expression divergence (Procrustes loading):")
print(f"  Ribosomal:     {ribo['mean_expression_divergence'].mean():.4f} "
      f"(SD {ribo['mean_expression_divergence'].std():.4f})")
print(f"  Non-ribosomal: {non_ribo['mean_expression_divergence'].mean():.4f} "
      f"(SD {non_ribo['mean_expression_divergence'].std():.4f})")

u_stat2, u_p2 = stats.mannwhitneyu(ribo["mean_expression_divergence"],
                                     non_ribo["mean_expression_divergence"],
                                     alternative="greater")
print(f"  Mann-Whitney U (ribo > non-ribo): p = {u_p2:.2e}")

# List top ribosomal genes by expression divergence
ribo_sorted = ribo.sort_values("mean_expression_divergence", ascending=False)
print(f"\nTop 10 ribosomal genes by expression divergence:")
for _, row in ribo_sorted.head(10).iterrows():
    print(f"  {row['gene_name']:10s}  identity={row['pct_identity']:.1f}%  "
          f"expr_div={row['mean_expression_divergence']:.4f}")

# ── Also characterize histones and HSPs ──────────────────────────────────────

hist = df[df["is_histone"]]
hsp = df[df["is_hsp"]]

print(f"\nOther constrained families in analysis set:")
print(f"  Histones: {len(hist)} genes, mean identity = "
      f"{hist['pct_identity'].mean():.2f}%" if len(hist) > 0 else "  Histones: 0 genes")
print(f"  HSPs:     {len(hsp)} genes, mean identity = "
      f"{hsp['pct_identity'].mean():.2f}%" if len(hsp) > 0 else "  HSPs: 0 genes")
print(f"  Total confound genes: {df['is_any_confound'].sum()}")

# ── STEP 2: Correlation excluding ribosomal genes ────────────────────────────

print("\n" + "="*70)
print("STEP 2: Spearman correlation EXCLUDING ribosomal genes")
print("="*70)

rho_noribo, p_noribo = stats.spearmanr(non_ribo["pct_identity"],
                                         non_ribo["mean_expression_divergence"])
print(f"\nN genes: {len(non_ribo)}")
print(f"Spearman rho: {rho_noribo:.4f}")
print(f"P-value: {p_noribo:.2e}")
print(f"Direction: {'NEGATIVE (conserved seq → more expr divergence)' if rho_noribo < 0 else 'POSITIVE'}")

# ── STEP 3: Correlation excluding all three families ─────────────────────────

print("\n" + "="*70)
print("STEP 3: Spearman correlation EXCLUDING ribosomal + histones + HSP")
print("="*70)

clean = df[~df["is_any_confound"]]
rho_clean, p_clean = stats.spearmanr(clean["pct_identity"],
                                       clean["mean_expression_divergence"])
print(f"\nN genes: {len(clean)} (removed {len(df) - len(clean)} confound genes)")
print(f"Spearman rho: {rho_clean:.4f}")
print(f"P-value: {p_clean:.2e}")
print(f"Direction: {'NEGATIVE' if rho_clean < 0 else 'POSITIVE'}")

# ── STEP 4: Expression-level normalization ───────────────────────────────────

print("\n" + "="*70)
print("STEP 4: Expression-level normalization")
print("="*70)

# Load centroid expression levels
centroids_h = pd.read_csv("output/phase2/centroids_human.csv", index_col=0)
centroids_m = pd.read_csv("output/phase2/centroids_mouse.csv", index_col=0)

# Mean expression across all cell types and both species
mean_expr_h = centroids_h.mean(axis=0)  # mean across 6 cell types per gene
mean_expr_m = centroids_m.mean(axis=0)

# Average human + mouse expression per gene
mean_expr = (mean_expr_h + mean_expr_m) / 2
mean_expr = mean_expr.to_frame("mean_expression_level")
mean_expr.index.name = "ensembl_gene_id"

# Merge expression levels
df_expr = df.merge(mean_expr, left_on="human_ensembl_id", right_index=True, how="left")

# Filter: only genes with non-zero expression
df_expr = df_expr[df_expr["mean_expression_level"] > 0].copy()
print(f"Genes with expression level > 0: {len(df_expr)}")

# Compute normalized divergence
df_expr["normalized_divergence"] = (df_expr["mean_expression_divergence"] /
                                     df_expr["mean_expression_level"])

# Handle inf/nan
df_expr = df_expr.replace([np.inf, -np.inf], np.nan).dropna(
    subset=["normalized_divergence"])
print(f"Genes after removing inf/nan: {len(df_expr)}")

rho_norm, p_norm = stats.spearmanr(df_expr["pct_identity"],
                                     df_expr["normalized_divergence"])
print(f"\nSpearman rho (identity vs normalized divergence): {rho_norm:.4f}")
print(f"P-value: {p_norm:.2e}")
print(f"Direction: {'NEGATIVE' if rho_norm < 0 else 'POSITIVE'}")

# Also: expression-normalized, EXCLUDING ribosomal
df_expr_noribo = df_expr[~df_expr["is_ribosomal"]]
rho_norm_noribo, p_norm_noribo = stats.spearmanr(
    df_expr_noribo["pct_identity"], df_expr_noribo["normalized_divergence"])
print(f"\nExpression-normalized, ALSO excluding ribosomal:")
print(f"  N = {len(df_expr_noribo)}, rho = {rho_norm_noribo:.4f}, p = {p_norm_noribo:.2e}")

# ── STEP 5: Summary report ──────────────────────────────────────────────────
# Use seq_divergence convention throughout (to match original rho = -0.234)

print("\n" + "="*70)
print("STEP 5: SUMMARY — Robustness of sequence-expression anticorrelation")
print("="*70)

# Recompute all-genes with seq_divergence (to match original convention)
rho_all_sd, p_all_sd = stats.spearmanr(df["seq_divergence"],
                                         df["mean_expression_divergence"])
rho_noribo_sd = -rho_noribo   # flip sign: pct_identity → seq_divergence
rho_clean_sd = -rho_clean
rho_norm_sd = -rho_norm
rho_norm_noribo_sd = -rho_norm_noribo

print(f"\nAll correlations use seq_divergence convention (matching original rho = -0.234)")
print(f"Negative rho = conserved sequence → MORE expression divergence")

results = {
    "all_genes": {
        "n": len(df), "rho": round(rho_all_sd, 4), "p": float(p_all_sd),
        "note": "Original (all 16,324 genes)"
    },
    "excl_ribosomal": {
        "n": len(non_ribo), "rho": round(rho_noribo_sd, 4), "p": float(p_noribo),
        "note": f"Excluded {len(ribo)} ribosomal genes"
    },
    "excl_ribo_hist_hsp": {
        "n": len(clean), "rho": round(rho_clean_sd, 4), "p": float(p_clean),
        "note": f"Excluded {len(df) - len(clean)} genes (ribo+hist+HSP)"
    },
    "expr_normalized": {
        "n": len(df_expr), "rho": round(rho_norm_sd, 4), "p": float(p_norm),
        "note": "divergence / mean_expression_level"
    },
    "expr_normalized_excl_ribo": {
        "n": len(df_expr_noribo), "rho": round(rho_norm_noribo_sd, 4),
        "p": float(p_norm_noribo),
        "note": "Normalized + excluded ribosomal"
    }
}

print(f"\n{'Condition':<35s} {'N':>6s} {'rho':>8s} {'p-value':>12s} {'Meaning':>15s}")
print("-" * 80)
for key, v in results.items():
    if v["rho"] < 0:
        meaning = "conserved→MORE"
    else:
        meaning = "conserved→LESS"
    label = key.replace("_", " ").title()
    print(f"{label:<35s} {v['n']:>6d} {v['rho']:>8.4f} {v['p']:>12.2e} {meaning:>15s}")

# Nuanced verdict based on actual pattern
ribo_exclusion_same_sign = (np.sign(results["all_genes"]["rho"]) ==
                             np.sign(results["excl_ribosomal"]["rho"]))
norm_flips_sign = (np.sign(results["all_genes"]["rho"]) !=
                    np.sign(results["expr_normalized"]["rho"]))

print()
if ribo_exclusion_same_sign and not norm_flips_sign:
    verdict = ("FULLY ROBUST: The anticorrelation persists across all conditions "
               "including gene family exclusion and expression normalization.")
elif ribo_exclusion_same_sign and norm_flips_sign:
    verdict = (
        "RIBOSOMAL CONFOUND NEGATIVE — EXPRESSION LEVEL CONFOUND POSITIVE. "
        "Excluding ribosomal genes barely changes the correlation "
        f"(rho {results['all_genes']['rho']:.3f} → {results['excl_ribosomal']['rho']:.3f}), "
        "so the anticorrelation is NOT driven by ribosomal proteins specifically. "
        "However, expression-level normalization REVERSES the sign "
        f"(rho → {results['expr_normalized']['rho']:+.3f}). "
        "This means the raw anticorrelation reflects a scale confound: "
        "highly conserved housekeeping genes are highly expressed, and highly "
        "expressed genes have larger absolute Procrustes loadings. "
        "After normalizing for expression level, conserved genes show LESS "
        "proportional expression divergence — the intuitive direction. "
        "The original finding (rho = -0.234) should be reported with this caveat."
    )
elif not ribo_exclusion_same_sign:
    verdict = (
        "RIBOSOMAL CONFOUND DETECTED: Removing ribosomal proteins changes "
        "the correlation direction. The original result is driven by this "
        "specific gene family."
    )
else:
    verdict = "Mixed results — see individual conditions."

print(f"VERDICT: {verdict}")

# ── Save outputs ─────────────────────────────────────────────────────────────

outdir = "output/phase2/dnds/confound_analysis"

# Save full annotated data
df.to_csv(f"{outdir}/genes_annotated_confound.csv", index=False)

# Save ribosomal gene list
ribo[["gene_name", "human_ensembl_id", "pct_identity",
      "mean_expression_divergence"]].to_csv(
    f"{outdir}/ribosomal_genes.csv", index=False)

# Save summary JSON
summary = {
    "verdict": verdict,
    "results": results,
    "ribosomal_stats": {
        "n_ribosomal": int(len(ribo)),
        "n_histone": int(len(hist)),
        "n_hsp": int(len(hsp)),
        "ribo_mean_identity": round(float(ribo["pct_identity"].mean()), 2),
        "nonribo_mean_identity": round(float(non_ribo["pct_identity"].mean()), 2),
        "ribo_mean_expr_div": round(float(ribo["mean_expression_divergence"].mean()), 4),
        "nonribo_mean_expr_div": round(float(non_ribo["mean_expression_divergence"].mean()), 4),
    }
}
with open(f"{outdir}/confound_summary.json", "w") as f:
    json.dump(summary, f, indent=2)

print(f"\nOutputs saved to {outdir}/")
print(f"  genes_annotated_confound.csv — full gene set with confound flags")
print(f"  ribosomal_genes.csv — ribosomal gene details")
print(f"  confound_summary.json — results summary")
