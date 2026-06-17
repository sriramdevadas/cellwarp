#!/usr/bin/env python3
# DIAGNOSTIC ONLY — not part of canonical validation pipeline
"""
Diagnostic: Annotation Granularity vs Procrustes Rigidity

Tests the adversarial claim that "flexible" cell types (large Procrustes
residuals) are simply types with broad, ambiguous annotations — i.e., that
annotation granularity drives the rigidity ranking rather than biology.

Biology: If "stromal cell" ranks as flexible merely because it is a catch-all
label encompassing heterogeneous subtypes, while "CD8-positive, alpha-beta
T cell" ranks as rigid because it denotes a narrow, homogeneous population,
then rigidity would be an artifact of annotation precision, not a signal of
evolutionary constraint.

Math: Spearman rank correlation between annotation specificity proxy and
Procrustes residual magnitude (higher residual = less rigid). If broad
annotations predict large residuals, ρ will be positive and significant.

Proxies used (multiple, for robustness):
  1. Manual expert classification: each of 35 cell types classified as
     "broad" (umbrella label covering known heterogeneous subtypes) or
     "specific" (narrowly defined, relatively homogeneous population).
     Tested with Mann-Whitney U (rank-sum) test.
  2. Label word count: number of whitespace-delimited tokens in the cell
     type label. More specific annotations tend to use more qualifying
     words.
  3. Label character length: total characters in the label string. Crude
     proxy, included for completeness.
  4. Cell Ontology depth: hierarchical depth in CL ontology (if available
     via OLS API). Deeper terms are more specific.

Outputs:
  - Console summary of all correlations and the rank-sum test
  - Scatter plot of word count vs residual magnitude
  - Box plot of broad vs specific residual distributions
  - CSV of per-cell-type annotation metrics
"""

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import os
import json

# -------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------
DATA_DIR = "output/phase2/scaled_35types"
OUT_DIR = "output/phase2/diagnostics/annotation_granularity"
os.makedirs(OUT_DIR, exist_ok=True)

# -------------------------------------------------------------------
# 1. Load residual data
# -------------------------------------------------------------------
print("=" * 70)
print("DIAGNOSTIC: Annotation Granularity vs Procrustes Rigidity")
print("=" * 70)

residuals = pd.read_csv(f"{DATA_DIR}/residuals_ranked.csv")
print(f"Loaded {residuals.shape[0]} cell types with residual magnitudes")
print()

# -------------------------------------------------------------------
# 2. Manual expert classification: broad vs specific
# -------------------------------------------------------------------
# Classification rationale:
#   "Broad" = umbrella label that subsumes multiple known subtypes or
#   encompasses heterogeneous populations across tissues/lineages.
#   "Specific" = narrowly defined cell type with a clear identity,
#   limited heterogeneity, and often tissue-restricted.
#
# This classification was made by examining each label against the
# Cell Ontology hierarchy and known cell biology.

BROAD_TYPES = {
    "stromal cell",           # Umbrella: fibroblasts, pericytes, MSCs, etc.
    "epithelial cell",        # Umbrella: many tissue-specific subtypes
    "T cell",                 # Umbrella: CD4+, CD8+, gamma-delta, etc.
    "endothelial cell",       # Present in all vascularized tissues
    "fibroblast",             # Broad mesenchymal — many tissue variants
    "monocyte",               # Umbrella for classical/non-classical/intermediate
    "granulocyte",            # Umbrella: neutrophils, eosinophils, basophils
    "macrophage",             # Tissue-resident vs recruited, many subtypes
    "B cell",                 # Naive, memory, germinal center, etc.
    "natural killer cell",    # CD56bright/dim, tissue variants
    "basal cell",             # Found in many epithelia (skin, lung, prostate)
    "myeloid leukocyte",      # Very broad myeloid category
    "neutrophil",             # Relatively defined but recent subtype discovery
    "mesenchymal stem cell",  # Broad progenitor category
    "plasma cell",            # Post-GC differentiation but subtypes exist
}

SPECIFIC_TYPES = {
    "CD8-positive, alpha-beta T cell",          # Narrow T cell subset
    "CD4-positive, alpha-beta T cell",          # Narrow T cell subset
    "non-classical monocyte",                   # Defined monocyte subset (CD14lo CD16+)
    "classical monocyte",                       # Defined monocyte subset (CD14+ CD16-)
    "intermediate monocyte",                    # Defined monocyte subset (CD14+ CD16+)
    "hepatocyte",                               # Single cell type, tissue-specific
    "smooth muscle cell",                       # Defined contractile cell
    "pancreatic ductal cell",                   # Tissue- and structure-specific
    "adventitial cell",                         # Specific perivascular cell
    "mature NK T cell",                         # Narrow lymphocyte subset
    "myeloid dendritic cell",                   # Defined DC subset
    "luminal epithelial cell of mammary gland", # Tissue- and layer-specific
    "large intestine goblet cell",              # Tissue- and function-specific
    "enterocyte of epithelium of large intestine",  # Tissue- and function-specific
    "fibroblast of cardiac tissue",             # Tissue-specific fibroblast
    "mesenchymal stem cell of adipose tissue",  # Tissue-specific MSC
    "hematopoietic stem cell",                  # Defined progenitor
    "hematopoietic precursor cell",             # Defined progenitor
    "pancreatic acinar cell",                   # Tissue- and function-specific
    "bladder urothelial cell",                  # Tissue-specific epithelial
}

# Verify all 35 types are classified
all_types = set(residuals["cell_type"].values)
classified = BROAD_TYPES | SPECIFIC_TYPES
assert classified == all_types, (
    f"Classification mismatch!\n"
    f"  Missing: {all_types - classified}\n"
    f"  Extra: {classified - all_types}"
)
print(f"Manual classification: {len(BROAD_TYPES)} broad, {len(SPECIFIC_TYPES)} specific")

# -------------------------------------------------------------------
# 3. Compute annotation proxies
# -------------------------------------------------------------------
print()
print("-" * 70)
print("Computing annotation specificity proxies")
print("-" * 70)

df = residuals.copy()
df["word_count"] = df["cell_type"].apply(lambda x: len(x.split()))
df["char_length"] = df["cell_type"].apply(len)
df["is_broad"] = df["cell_type"].apply(lambda x: x in BROAD_TYPES)
df["specificity_class"] = df["is_broad"].map({True: "broad", False: "specific"})

# Invert rank so higher = more rigid (rank 1 = least rigid in original)
# Actually, residual_magnitude directly captures flexibility — use it.
# Rigidity rank: rank 1 = most flexible (largest residual) in the CSV.

print("\nAnnotation metrics by cell type:")
print(df[["rank", "cell_type", "residual_magnitude", "word_count",
          "char_length", "specificity_class"]].to_string(index=False))

# -------------------------------------------------------------------
# 4. Spearman correlations
# -------------------------------------------------------------------
print()
print("-" * 70)
print("Spearman correlations: annotation proxy vs residual magnitude")
print("-" * 70)

results = {}

# Word count vs residual magnitude
rho_wc, p_wc = stats.spearmanr(df["word_count"], df["residual_magnitude"])
results["word_count"] = {"rho": rho_wc, "p": p_wc}
print(f"  Word count vs residual:  ρ = {rho_wc:+.4f}, p = {p_wc:.4f}")

# Char length vs residual magnitude
rho_cl, p_cl = stats.spearmanr(df["char_length"], df["residual_magnitude"])
results["char_length"] = {"rho": rho_cl, "p": p_cl}
print(f"  Char length vs residual: ρ = {rho_cl:+.4f}, p = {p_cl:.4f}")

# Note: negative ρ would mean more specific labels (longer names) have
# LARGER residuals — the opposite of the confound hypothesis.
# The confound predicts POSITIVE ρ between specificity and rigidity,
# i.e., NEGATIVE ρ between word_count and residual_magnitude
# (specific labels → small residuals → rigid).

# -------------------------------------------------------------------
# 5. Mann-Whitney U test: broad vs specific
# -------------------------------------------------------------------
print()
print("-" * 70)
print("Mann-Whitney U test: broad vs specific cell types")
print("-" * 70)

broad_resid = df.loc[df["is_broad"], "residual_magnitude"].values
specific_resid = df.loc[~df["is_broad"], "residual_magnitude"].values

print(f"  Broad types (n={len(broad_resid)}):")
print(f"    Median residual: {np.median(broad_resid):.3f}")
print(f"    Mean residual:   {np.mean(broad_resid):.3f}")
print(f"    Range: [{np.min(broad_resid):.3f}, {np.max(broad_resid):.3f}]")
print(f"  Specific types (n={len(specific_resid)}):")
print(f"    Median residual: {np.median(specific_resid):.3f}")
print(f"    Mean residual:   {np.mean(specific_resid):.3f}")
print(f"    Range: [{np.min(specific_resid):.3f}, {np.max(specific_resid):.3f}]")

U_stat, p_mw = stats.mannwhitneyu(
    broad_resid, specific_resid, alternative="two-sided"
)
# Compute rank-biserial correlation as effect size
n1, n2 = len(broad_resid), len(specific_resid)
r_rb = 1 - (2 * U_stat) / (n1 * n2)
results["broad_vs_specific"] = {"U": U_stat, "p": p_mw, "r_rankbiserial": r_rb}

print(f"\n  Mann-Whitney U = {U_stat:.1f}, p = {p_mw:.4f}")
print(f"  Rank-biserial r = {r_rb:+.4f}")
direction = "broad > specific" if np.median(broad_resid) > np.median(specific_resid) else "specific > broad"
print(f"  Direction: {direction} residuals")

# -------------------------------------------------------------------
# 6. Cell Ontology depth via OLS API (attempt)
# -------------------------------------------------------------------
print()
print("-" * 70)
print("Attempting Cell Ontology lookup via OLS API...")
print("-" * 70)

# Map cell type labels to Cell Ontology (CL) term IDs.
# These were manually curated by searching the OLS for each label.
CL_TERM_MAP = {
    "stromal cell": "CL:0000499",
    "epithelial cell": "CL:0000066",
    "hematopoietic precursor cell": "CL:0008001",
    "hematopoietic stem cell": "CL:0000037",
    "pancreatic acinar cell": "CL:0002064",
    "basal cell": "CL:0000646",
    "T cell": "CL:0000084",
    "neutrophil": "CL:0000775",
    "fibroblast of cardiac tissue": "CL:0002548",
    "myeloid leukocyte": "CL:0000766",
    "mesenchymal stem cell of adipose tissue": "CL:0002570",
    "plasma cell": "CL:0000786",
    "mesenchymal stem cell": "CL:0000134",
    "CD4-positive, alpha-beta T cell": "CL:0000624",
    "classical monocyte": "CL:0000860",
    "macrophage": "CL:0000235",
    "B cell": "CL:0000236",
    "luminal epithelial cell of mammary gland": "CL:0002326",
    "large intestine goblet cell": "CL:1000320",
    "enterocyte of epithelium of large intestine": "CL:0002071",
    "myeloid dendritic cell": "CL:0000782",
    "monocyte": "CL:0000576",
    "natural killer cell": "CL:0000623",
    "intermediate monocyte": "CL:0002393",
    "mature NK T cell": "CL:0000814",
    "adventitial cell": "CL:0002503",
    "granulocyte": "CL:0000094",
    "fibroblast": "CL:0000057",
    "bladder urothelial cell": "CL:1001428",
    "pancreatic ductal cell": "CL:0002079",
    "smooth muscle cell": "CL:0000192",
    "hepatocyte": "CL:0000182",
    "endothelial cell": "CL:0000115",
    "non-classical monocyte": "CL:0000875",
    "CD8-positive, alpha-beta T cell": "CL:0000625",
}

cl_depth_success = False
try:
    import urllib.request
    import urllib.error

    def get_cl_descendants(cl_id, timeout=10):
        """Query OLS4 API for number of descendants of a CL term."""
        # OLS4 API endpoint for term hierarchy
        encoded_iri = urllib.parse.quote_plus(
            urllib.parse.quote_plus(f"http://purl.obolibrary.org/obo/{cl_id.replace(':', '_')}")
        )
        url = (f"https://www.ebi.ac.uk/ols4/api/ontologies/cl/terms/"
               f"{encoded_iri}/hierarchicalDescendants?size=1")
        req = urllib.request.Request(url)
        req.add_header("Accept", "application/json")
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = json.loads(resp.read().decode())
                # Total number of descendants from pagination info
                total = data.get("page", {}).get("totalElements", 0)
                return total
        except (urllib.error.URLError, urllib.error.HTTPError, Exception) as e:
            return None

    import urllib.parse

    # Test with one term first
    test_id = CL_TERM_MAP["hepatocyte"]
    test_result = get_cl_descendants(test_id)

    if test_result is not None:
        print(f"  OLS API accessible. Test: hepatocyte has {test_result} descendants.")
        print("  Querying all 35 terms (this may take ~1 minute)...")

        descendant_counts = {}
        for ct, cl_id in CL_TERM_MAP.items():
            count = get_cl_descendants(cl_id)
            if count is not None:
                descendant_counts[ct] = count
            else:
                print(f"    WARNING: Failed to get descendants for {ct} ({cl_id})")

        if len(descendant_counts) >= 30:  # Accept if most succeeded
            df["cl_descendants"] = df["cell_type"].map(descendant_counts)
            valid = df.dropna(subset=["cl_descendants"])
            rho_cl_desc, p_cl_desc = stats.spearmanr(
                valid["cl_descendants"], valid["residual_magnitude"]
            )
            results["cl_descendants"] = {
                "rho": rho_cl_desc, "p": p_cl_desc,
                "n_terms": len(valid)
            }
            print(f"\n  CL descendants vs residual (n={len(valid)}):")
            print(f"    ρ = {rho_cl_desc:+.4f}, p = {p_cl_desc:.4f}")
            print(f"    (Positive ρ = more descendants → larger residual → less rigid)")
            cl_depth_success = True
        else:
            print(f"  Only {len(descendant_counts)}/35 terms succeeded. Skipping CL proxy.")
    else:
        print("  OLS API not accessible. Skipping Cell Ontology proxy.")

except ImportError:
    print("  urllib not available. Skipping Cell Ontology proxy.")
except Exception as e:
    print(f"  CL lookup failed: {e}. Skipping Cell Ontology proxy.")

# -------------------------------------------------------------------
# 7. Summary and interpretation
# -------------------------------------------------------------------
print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)

print("\nProxy correlations with residual magnitude:")
print(f"  {'Proxy':<30s} {'ρ':>8s} {'p-value':>10s} {'|ρ|≥0.3?':>10s} {'p≤0.05?':>10s}")
print(f"  {'-'*30} {'-'*8} {'-'*10} {'-'*10} {'-'*10}")

confound_detected = False
for name, res in results.items():
    if name == "broad_vs_specific":
        continue
    rho = res["rho"]
    p = res["p"]
    flag_rho = "YES" if abs(rho) >= 0.3 else "no"
    flag_p = "YES" if p <= 0.05 else "no"
    if abs(rho) >= 0.3 or p <= 0.05:
        confound_detected = True
    print(f"  {name:<30s} {rho:>+8.4f} {p:>10.4f} {flag_rho:>10s} {flag_p:>10s}")

print(f"\nBroad vs specific (Mann-Whitney U):")
mw = results["broad_vs_specific"]
print(f"  U = {mw['U']:.1f}, p = {mw['p']:.4f}, rank-biserial r = {mw['r_rankbiserial']:+.4f}")
flag_mw = "YES" if mw["p"] <= 0.05 else "no"
print(f"  Significant at α=0.05? {flag_mw}")
if mw["p"] <= 0.05:
    confound_detected = True

# -------------------------------------------------------------------
# Interpretation
# -------------------------------------------------------------------
print()
print("-" * 70)
if confound_detected:
    print("RESULT: Annotation granularity confound detected.")
    print("  At least one proxy shows |ρ| ≥ 0.3 or p ≤ 0.05.")
    print("  RECOMMENDATION: ADVISOR must decide on limitation language")
    print("  or analytical correction (e.g., restrict to 'specific' types")
    print("  and re-run Procrustes to verify ranking stability).")
else:
    print("RESULT: No evidence that annotation granularity drives rigidity ranking.")
    print("  All proxies show |ρ| < 0.3 and p > 0.05.")
    print("  RECOMMENDATION: Add one sentence to Discussion noting this was tested.")
print("-" * 70)

# -------------------------------------------------------------------
# 8. Figures
# -------------------------------------------------------------------
print()
print("Generating figures...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# (a) Word count vs residual magnitude
ax = axes[0]
colors = ["#d62728" if b else "#1f77b4" for b in df["is_broad"]]
ax.scatter(df["word_count"], df["residual_magnitude"], c=colors, s=60, alpha=0.7, edgecolors="k", linewidths=0.5)
# Add cell type labels for extremes
for _, row in df.iterrows():
    if row["residual_magnitude"] > 14 or row["word_count"] >= 6:
        ax.annotate(
            row["cell_type"], (row["word_count"], row["residual_magnitude"]),
            fontsize=6, ha="left", va="bottom", rotation=15,
            xytext=(3, 3), textcoords="offset points"
        )
rho_wc = results["word_count"]["rho"]
p_wc = results["word_count"]["p"]
ax.set_xlabel("Label word count", fontsize=12)
ax.set_ylabel("Procrustes residual magnitude", fontsize=12)
ax.set_title(f"Word Count vs Rigidity\nSpearman ρ = {rho_wc:+.3f}, p = {p_wc:.3f}", fontsize=11)
ax.legend(
    handles=[
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#d62728", markersize=8, label="Broad"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#1f77b4", markersize=8, label="Specific"),
    ],
    loc="upper right", fontsize=9
)

# (b) Box plot: broad vs specific
ax = axes[1]
data_box = [broad_resid, specific_resid]
bp = ax.boxplot(data_box, tick_labels=["Broad\n(n={})".format(len(broad_resid)),
                                       "Specific\n(n={})".format(len(specific_resid))],
                patch_artist=True, widths=0.5)
bp["boxes"][0].set_facecolor("#d62728")
bp["boxes"][0].set_alpha(0.4)
bp["boxes"][1].set_facecolor("#1f77b4")
bp["boxes"][1].set_alpha(0.4)
# Overlay individual points
for i, (data, color) in enumerate(zip(data_box, ["#d62728", "#1f77b4"])):
    x = np.random.normal(i + 1, 0.04, size=len(data))
    ax.scatter(x, data, c=color, s=30, alpha=0.6, edgecolors="k", linewidths=0.5, zorder=3)

p_mw = results["broad_vs_specific"]["p"]
ax.set_ylabel("Procrustes residual magnitude", fontsize=12)
ax.set_title(f"Broad vs Specific Annotations\nMann-Whitney p = {p_mw:.3f}", fontsize=11)

plt.tight_layout()
fig_path = f"{OUT_DIR}/annotation_granularity_vs_rigidity.png"
plt.savefig(fig_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved: {fig_path}")

# -------------------------------------------------------------------
# 9. Save data
# -------------------------------------------------------------------
cols_save = ["rank", "cell_type", "residual_magnitude", "pct_of_ssr",
             "word_count", "char_length", "specificity_class"]
if "cl_descendants" in df.columns:
    cols_save.append("cl_descendants")
csv_path = f"{OUT_DIR}/annotation_metrics.csv"
df[cols_save].to_csv(csv_path, index=False)
print(f"  Saved: {csv_path}")

# Save summary JSON
summary = {
    "proxy_used": "Multiple: manual broad/specific classification, label word count, label character length"
        + (", Cell Ontology descendant count" if cl_depth_success else ""),
    "n_cell_types": len(df),
    "n_broad": int(df["is_broad"].sum()),
    "n_specific": int((~df["is_broad"]).sum()),
    "correlations": {k: {kk: round(float(vv), 4) for kk, vv in v.items()} for k, v in results.items()},
    "confound_detected": confound_detected,
}
json_path = f"{OUT_DIR}/annotation_granularity_summary.json"
with open(json_path, "w") as f:
    json.dump(summary, f, indent=2)
print(f"  Saved: {json_path}")

print()
print("DIAGNOSTIC: Annotation Granularity vs Rigidity — COMPLETE")
print("=" * 70)
