#!/usr/bin/env python3
"""
Thread 4B: Transcription factor network complexity vs Procrustes rigidity.

Biology:
    Transcription factors (TFs) orchestrate gene expression programs. The
    DoRothEA database catalogs TF-target regulatory interactions curated from
    literature and ChIP-seq experiments (confidence levels A-C = highest
    quality). Using decoupler's univariate linear model (ULM), we estimate
    TF activity from expression data: for each TF, how well do its known
    target genes' expression match the expected activation pattern?

    Hypothesis: cell types with more complex TF regulatory networks —
    more active TFs, broader regulatory control — are more evolutionarily
    rigid (lower Procrustes residual). Complex networks resist evolutionary
    drift because altering one TF's activity disrupts many downstream targets
    simultaneously (pleiotropic constraint at the regulatory level).

Math:
    ULM fits y_i = β₀ + β₁x_i + ε for each TF, where y is gene expression
    across all genes and x is the TF's target weight vector. The t-statistic
    of β₁ is the TF activity score (positive = active, negative = repressed).

    Complexity metrics per cell type:
      a) n_active_tfs: |{TF : |activity| > 1.0}|  (z-score threshold)
      b) mean_activity: mean(|activity|) across all TFs
      c) activity_entropy: H(p) where p_i = |activity_i| / Σ|activity|
         High entropy = many TFs active at moderate levels = complex network
         Low entropy = few dominant TFs = simple/concentrated regulation

Data:
    DoRothEA regulon: levels A+B+C from OmniPath via decoupler (v2.1.4)
    Centroids: output/phase2/scaled_35types/centroids_{human,mouse}_35.csv
    Residuals: output/phase2/scaled_35types/procrustes_results_35.json
    Annotations: output/phase2/developmental_constraint/developmental_annotations.csv
    Ortholog map: data/phase1/orthologs_human_mouse.csv (Ensembl→symbol)

Output:
    output/phase2/mechanistic/tf_complexity/
"""

import json
from pathlib import Path

import decoupler as dc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# ── Paths ────────────────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "output" / "phase2" / "mechanistic" / "tf_complexity"
OUT.mkdir(parents=True, exist_ok=True)

CENTROIDS_H = ROOT / "output" / "phase2" / "scaled_35types" / "centroids_human_35.csv"
CENTROIDS_M = ROOT / "output" / "phase2" / "scaled_35types" / "centroids_mouse_35.csv"
RESULTS_JSON = ROOT / "output" / "phase2" / "scaled_35types" / "procrustes_results_35.json"
ANNOTATIONS = ROOT / "output" / "phase2" / "developmental_constraint" / "developmental_annotations.csv"
ORTHOLOGS = ROOT / "data" / "phase1" / "orthologs_human_mouse.csv"


def load_ortholog_maps():
    """Build Ensembl→symbol and symbol→Ensembl maps from ortholog table."""
    orth = pd.read_csv(ORTHOLOGS)
    ens_to_sym = dict(zip(orth["human_ensembl_id"], orth["human_gene_name"]))
    sym_to_ens = dict(zip(orth["human_gene_name"], orth["human_ensembl_id"]))
    return ens_to_sym, sym_to_ens


def load_centroids_as_symbols(path, ens_to_sym):
    """Load centroid CSV and rename Ensembl columns to gene symbols.

    Multiple Ensembl IDs can map to the same symbol (e.g. pseudogenes).
    We average expression across duplicates to produce unique column names.
    """
    df = pd.read_csv(path)
    gene_cols = [c for c in df.columns if c.startswith("ENSG")]
    # Map Ensembl → symbol, drop unmapped
    rename = {ens: ens_to_sym[ens] for ens in gene_cols if ens in ens_to_sym}
    df_renamed = df[["cell_type"] + list(rename.keys())].rename(columns=rename)
    df_renamed = df_renamed.set_index("cell_type")
    # Average duplicate symbols
    n_dups = df_renamed.columns.duplicated().sum()
    if n_dups > 0:
        df_renamed = df_renamed.T.groupby(level=0).mean().T
    return df_renamed


def load_residuals():
    """Load per-cell-type Procrustes residual magnitudes."""
    with open(RESULTS_JSON) as f:
        data = json.load(f)
    residuals = {ct: info["magnitude"] for ct, info in data["residuals"].items()}
    return pd.Series(residuals, name="residual_magnitude")


def load_annotations():
    """Load developmental annotations for cell category coloring."""
    return pd.read_csv(ANNOTATIONS)


def assign_category(row):
    """Assign broad cell category from lineage annotations."""
    lin = row.get("lineage", "")
    if lin == "hematopoietic":
        return "immune"
    elif lin == "epithelial":
        return "epithelial"
    elif lin in ("mesenchymal", "endothelial"):
        return "stromal"
    else:
        return "other"


def compute_complexity_metrics(activity_df):
    """
    Compute TF network complexity metrics from activity score matrix.

    Parameters
    ----------
    activity_df : DataFrame
        Cell types × TFs, values are ULM activity scores (t-statistics).

    Returns
    -------
    DataFrame with columns: n_active_tfs, mean_activity, activity_entropy
    """
    abs_act = activity_df.abs()

    # a) Number of TFs with |activity| > 1.0
    n_active = (abs_act > 1.0).sum(axis=1)

    # b) Mean absolute activity across all TFs
    mean_act = abs_act.mean(axis=1)

    # c) Shannon entropy of the absolute activity distribution
    #    Normalize to probability distribution per cell type
    row_sums = abs_act.sum(axis=1)
    probs = abs_act.div(row_sums, axis=0)
    # Entropy: -Σ p*log(p), handling zeros
    log_probs = np.log2(probs.replace(0, np.nan))
    entropy = -(probs * log_probs).sum(axis=1)

    return pd.DataFrame({
        "n_active_tfs": n_active,
        "mean_activity": mean_act,
        "activity_entropy": entropy,
    })


def correlate_and_report(metrics_df, residuals, species):
    """Compute Spearman correlations for each metric vs residual."""
    common = metrics_df.index.intersection(residuals.index)
    res = residuals.loc[common]
    results = {}
    metric_names = {
        "n_active_tfs": "N active TFs (|score|>1)",
        "mean_activity": "Mean |activity|",
        "activity_entropy": "Activity entropy",
    }
    for col in ["n_active_tfs", "mean_activity", "activity_entropy"]:
        vals = metrics_df.loc[common, col]
        rho, pval = stats.spearmanr(vals, res)
        results[col] = {"rho": rho, "pval": pval}
        sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "n.s."
        print(f"  {metric_names[col]:30s}: ρ = {rho:+.4f}, p = {pval:.4f} ({sig})")
    return results


def plot_metric_scatter(
    metric_values, residuals, categories, metric_name, metric_label,
    species, rho, pval, out_path,
):
    """Scatter plot: one TF complexity metric vs residual magnitude."""
    fig, ax = plt.subplots(figsize=(10, 8))
    colors = {
        "immune": "#2196F3",
        "epithelial": "#4CAF50",
        "stromal": "#FF9800",
        "other": "#9E9E9E",
    }
    common = metric_values.index.intersection(residuals.index).intersection(categories.index)
    mv = metric_values.loc[common]
    res = residuals.loc[common]
    cats = categories.loc[common]

    for cat, color in colors.items():
        mask = cats == cat
        if mask.sum() == 0:
            continue
        ax.scatter(
            mv.loc[mask], res.loc[mask],
            c=color, label=cat, s=60, alpha=0.8, edgecolors="white", linewidth=0.5,
        )

    for ct in common:
        label = ct if len(ct) <= 25 else ct[:22] + "..."
        ax.annotate(
            label, (mv.loc[ct], res.loc[ct]),
            fontsize=5.5, alpha=0.75, xytext=(4, 4), textcoords="offset points",
        )

    ax.set_xlabel(metric_label, fontsize=11)
    ax.set_ylabel("Procrustes Residual Magnitude\n(evolutionary divergence)", fontsize=11)
    ax.set_title(
        f"{metric_name} vs Evolutionary Rigidity — {species}\n"
        f"Spearman ρ = {rho:+.3f}, p = {pval:.4f}",
        fontsize=12,
    )
    ax.legend(title="Cell category", fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def main():
    print("=" * 70)
    print("Thread 4B: TF Network Complexity vs Procrustes Rigidity")
    print("=" * 70)

    # ── Step 1: Load DoRothEA regulon ────────────────────────────────────
    print("\n── Step 1: Loading DoRothEA regulon via decoupler ──")
    regulon = dc.op.dorothea(organism="human", levels=["A", "B", "C"])
    n_pairs = len(regulon)
    n_tfs = regulon["source"].nunique()
    n_targets = regulon["target"].nunique()
    print(f"  Total TF-target pairs: {n_pairs}")
    print(f"  Unique TFs: {n_tfs}")
    print(f"  Unique target genes: {n_targets}")

    # ── Step 2: Filter regulon to our gene space ─────────────────────────
    print("\n── Step 2: Filtering regulon to our 16,959-gene ortholog space ──")
    ens_to_sym, sym_to_ens = load_ortholog_maps()
    our_symbols = set(ens_to_sym.values())

    # Filter: keep only target genes present in our space
    regulon_filt = regulon[regulon["target"].isin(our_symbols)].copy()
    # Also filter: keep only TFs that are themselves genes in our space (for interpretability)
    # But TFs don't need to be in expression matrix for ULM — only targets do
    n_pairs_filt = len(regulon_filt)
    n_tfs_filt = regulon_filt["source"].nunique()
    n_targets_filt = regulon_filt["target"].nunique()
    print(f"  After filtering targets to ortholog space:")
    print(f"    TF-target pairs: {n_pairs_filt} ({n_pairs_filt}/{n_pairs} = {100*n_pairs_filt/n_pairs:.1f}%)")
    print(f"    Unique TFs: {n_tfs_filt} ({n_tfs_filt}/{n_tfs})")
    print(f"    Unique targets: {n_targets_filt} ({n_targets_filt}/{n_targets})")

    # ── Step 3: Load centroids and run ULM ───────────────────────────────
    print("\n── Step 3: Loading centroids and running TF activity inference ──")
    centroids_h = load_centroids_as_symbols(CENTROIDS_H, ens_to_sym)
    centroids_m = load_centroids_as_symbols(CENTROIDS_M, ens_to_sym)
    print(f"  Human centroids: {centroids_h.shape[0]} cell types × {centroids_h.shape[1]} genes (symbols)")
    print(f"  Mouse centroids: {centroids_m.shape[0]} cell types × {centroids_m.shape[1]} genes (symbols)")

    print("\n  Running ULM on human centroids...")
    activity_h, pvals_h = dc.mt.ulm(centroids_h, regulon_filt)
    print(f"  Human TF activity matrix: {activity_h.shape[0]} cell types × {activity_h.shape[1]} TFs")

    print("  Running ULM on mouse centroids...")
    activity_m, pvals_m = dc.mt.ulm(centroids_m, regulon_filt)
    print(f"  Mouse TF activity matrix: {activity_m.shape[0]} cell types × {activity_m.shape[1]} TFs")

    # ── Step 4: Compute complexity metrics ───────────────────────────────
    print("\n── Step 4: Computing TF network complexity metrics ──")
    metrics_h = compute_complexity_metrics(activity_h)
    metrics_m = compute_complexity_metrics(activity_m)

    print(f"\n  Human metrics summary:")
    print(f"    N active TFs:     mean={metrics_h['n_active_tfs'].mean():.1f}, "
          f"range=[{metrics_h['n_active_tfs'].min()}, {metrics_h['n_active_tfs'].max()}]")
    print(f"    Mean |activity|:  mean={metrics_h['mean_activity'].mean():.3f}, "
          f"range=[{metrics_h['mean_activity'].min():.3f}, {metrics_h['mean_activity'].max():.3f}]")
    print(f"    Entropy:          mean={metrics_h['activity_entropy'].mean():.3f}, "
          f"range=[{metrics_h['activity_entropy'].min():.3f}, {metrics_h['activity_entropy'].max():.3f}]")

    print(f"\n  Top 5 most complex (by n_active_tfs, human):")
    for ct, n in metrics_h["n_active_tfs"].sort_values(ascending=False).head(5).items():
        print(f"    {ct}: {n} active TFs")
    print(f"  Top 5 least complex (by n_active_tfs, human):")
    for ct, n in metrics_h["n_active_tfs"].sort_values().head(5).items():
        print(f"    {ct}: {n} active TFs")

    # ── Step 5: Correlate with Procrustes residuals ──────────────────────
    print("\n── Step 5: Correlating with Procrustes residual magnitude ──")
    residuals = load_residuals()
    print(f"  Loaded residuals for {len(residuals)} cell types")

    # Load annotations
    annot = load_annotations()
    annot["category"] = annot.apply(assign_category, axis=1)
    cat_map = annot.set_index("cell_type")["category"]

    print(f"\n  Human correlations:")
    corr_h = correlate_and_report(metrics_h, residuals, "Human")

    # Plot each metric
    metric_info = {
        "n_active_tfs": ("N Active TFs", "Number of active TFs (|activity score| > 1.0)"),
        "mean_activity": ("Mean TF Activity", "Mean absolute TF activity score"),
        "activity_entropy": ("TF Activity Entropy", "Shannon entropy of |activity| distribution (bits)"),
    }
    for col, (name, label) in metric_info.items():
        rho = corr_h[col]["rho"]
        pval = corr_h[col]["pval"]
        plot_metric_scatter(
            metrics_h[col], residuals, cat_map,
            name, label, "Human", rho, pval,
            OUT / f"{col}_vs_residual_human.png",
        )

    # ── Step 6: Mouse replication ────────────────────────────────────────
    print(f"\n── Step 6: Mouse replication ──")
    print(f"\n  Mouse correlations:")
    corr_m = correlate_and_report(metrics_m, residuals, "Mouse")

    for col, (name, label) in metric_info.items():
        rho = corr_m[col]["rho"]
        pval = corr_m[col]["pval"]
        plot_metric_scatter(
            metrics_m[col], residuals, cat_map,
            name, label, "Mouse", rho, pval,
            OUT / f"{col}_vs_residual_mouse.png",
        )

    # ── Cross-species TF activity agreement ──────────────────────────────
    print(f"\n── Cross-species TF activity agreement ──")
    common_tfs = activity_h.columns.intersection(activity_m.columns)
    common_cts = activity_h.index.intersection(activity_m.index)
    print(f"  Shared TFs: {len(common_tfs)}, shared cell types: {len(common_cts)}")

    # Per-metric cross-species agreement
    for col in ["n_active_tfs", "mean_activity", "activity_entropy"]:
        common_idx = metrics_h.index.intersection(metrics_m.index)
        rho, pval = stats.spearmanr(
            metrics_h.loc[common_idx, col], metrics_m.loc[common_idx, col]
        )
        print(f"  {col:25s}: human-mouse ρ = {rho:+.4f}, p = {pval:.4e}")

    # ── Save results ─────────────────────────────────────────────────────
    print("\n── Saving results ──")

    # Combined CSV
    common = metrics_h.index.intersection(residuals.index)
    results_df = pd.DataFrame({
        "cell_type": common,
        "n_active_tfs_human": metrics_h.loc[common, "n_active_tfs"].values,
        "n_active_tfs_mouse": metrics_m.loc[common, "n_active_tfs"].values,
        "mean_activity_human": metrics_h.loc[common, "mean_activity"].values,
        "mean_activity_mouse": metrics_m.loc[common, "mean_activity"].values,
        "entropy_human": metrics_h.loc[common, "activity_entropy"].values,
        "entropy_mouse": metrics_m.loc[common, "activity_entropy"].values,
        "residual_magnitude": residuals.loc[common].values,
    })
    results_df = results_df.merge(
        annot[["cell_type", "lineage", "progenitor"]], on="cell_type", how="left",
    )
    results_df = results_df.sort_values("residual_magnitude")
    csv_path = OUT / "tf_complexity_vs_residual.csv"
    results_df.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path}")

    # Save TF activity matrices
    activity_h.to_csv(OUT / "tf_activity_human.csv")
    activity_m.to_csv(OUT / "tf_activity_mouse.csv")
    print(f"  Saved: {OUT / 'tf_activity_human.csv'}")
    print(f"  Saved: {OUT / 'tf_activity_mouse.csv'}")

    # JSON summary
    summary = {
        "analysis": "Thread 4B: TF network complexity vs Procrustes rigidity",
        "dorothea_levels": ["A", "B", "C"],
        "total_tf_target_pairs": n_pairs,
        "filtered_tf_target_pairs": n_pairs_filt,
        "n_tfs_total": n_tfs,
        "n_tfs_after_filter": n_tfs_filt,
        "n_targets_in_ortholog_space": n_targets_filt,
        "n_cell_types": len(common),
        "n_tfs_in_activity_matrix": int(activity_h.shape[1]),
        "human_correlations": {
            col: {"spearman_rho": round(corr_h[col]["rho"], 4),
                  "p_value": round(corr_h[col]["pval"], 4),
                  "significant_at_005": bool(corr_h[col]["pval"] < 0.05)}
            for col in corr_h
        },
        "mouse_correlations": {
            col: {"spearman_rho": round(corr_m[col]["rho"], 4),
                  "p_value": round(corr_m[col]["pval"], 4),
                  "significant_at_005": bool(corr_m[col]["pval"] < 0.05)}
            for col in corr_m
        },
    }
    json_path = OUT / "tf_complexity_results.json"
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  Saved: {json_path}")

    # ── Summary ──────────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"\nDoRothEA regulon: {n_tfs_filt} TFs, {n_pairs_filt} pairs (levels A+B+C)")
    print(f"ULM activity matrix: {activity_h.shape[0]} cell types × {activity_h.shape[1]} TFs")

    print(f"\n  {'Metric':<30s} {'Human ρ':>10s} {'p':>10s} {'Mouse ρ':>10s} {'p':>10s}")
    print(f"  {'-'*30} {'-'*10} {'-'*10} {'-'*10} {'-'*10}")
    for col in ["n_active_tfs", "mean_activity", "activity_entropy"]:
        rh = corr_h[col]["rho"]
        ph = corr_h[col]["pval"]
        rm = corr_m[col]["rho"]
        pm = corr_m[col]["pval"]
        print(f"  {col:<30s} {rh:>+10.4f} {ph:>10.4f} {rm:>+10.4f} {pm:>10.4f}")

    print("\n── Interpretation ──")
    # Check if any metric is significant in human
    any_sig_h = any(corr_h[col]["pval"] < 0.05 for col in corr_h)
    any_sig_m = any(corr_m[col]["pval"] < 0.05 for col in corr_m)

    if any_sig_h:
        sig_metrics = [col for col in corr_h if corr_h[col]["pval"] < 0.05]
        for col in sig_metrics:
            rho = corr_h[col]["rho"]
            if rho < 0:
                print(f"HYPOTHESIS SUPPORTED ({col}): Cell types with greater TF network")
                print(f"complexity have LOWER Procrustes residuals (more rigid). This is")
                print(f"consistent with complex regulatory networks resisting evolutionary drift.")
            else:
                print(f"UNEXPECTED POSITIVE ({col}): More complex TF networks → MORE diverged.")
                print(f"This contradicts the constraint hypothesis. Possible explanation: complex")
                print(f"networks have more parameters that CAN change, increasing divergence capacity.")
        if any_sig_m:
            print(f"\nMouse replication: {'CONFIRMED' if all(np.sign(corr_m[col]['rho']) == np.sign(corr_h[col]['rho']) for col in sig_metrics) else 'MIXED'}.")
        else:
            print(f"\nMouse: Not significant — human result does not replicate.")
    else:
        print("HYPOTHESIS NOT SUPPORTED: No TF complexity metric significantly correlates")
        print("with Procrustes residual magnitude in human. TF regulatory network complexity,")
        print("as measured by DoRothEA activity inference, does not explain cross-species")
        print("rigidity differences.")
        if not any_sig_m:
            print("\nMouse: Also not significant (consistent null in both species).")
        else:
            print(f"\nMouse: Some metrics significant — inconsistent with human null.")

    print(f"\nAll outputs saved to: {OUT}")
    print("=" * 70)


if __name__ == "__main__":
    main()
