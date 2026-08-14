"""
CellWarp — L1000 Random Baseline Analysis

Establishes whether the observed ρ=0.852 Spearman correlation between
L1000-landmark rigidity ranking and full-space rigidity ranking is
significantly above what random gene subsets of equal size produce.

Biology
-------
The L1000 platform uses 978 "landmark" genes chosen by the Broad Institute
to predict whole-transcriptome variance. A reviewer argues that because these
genes were selected for transcriptomic predictability, high ρ is trivially
expected. This script tests that claim by drawing random gene subsets from the
ortholog space and asking: does an arbitrary ~907-gene subset produce a
similarly high ρ?

Math
----
For each of N iterations:
    1. Sample 907 genes uniformly at random (without replacement) from the
       16,959-gene ortholog space.
    2. Subset human and mouse centroid matrices to those genes.
    3. PCA-reduce the subsetted centroids (≥95% variance).
    4. Procrustes alignment (rotation + scaling).
    5. Compute per-cell-type residual magnitudes → rigidity ranking.
    6. Spearman ρ between random-subset ranking and full-space ranking.

The distribution of ρ values across iterations is the null distribution.
The observed L1000 ρ = 0.852 is compared against this null.

Output
------
- Console summary with mean, SD, percentiles, z-score, empirical p
- Histogram: output/figures/l1000_random_baseline.pdf
- JSON summary: output/figures/l1000_random_baseline_results.json
"""

import json
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# Add src/ to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from cellwarp.procrustes import pca_reduce_centroids, procrustes_align

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT = Path(__file__).resolve().parent.parent
CENTROIDS_HUMAN = PROJECT / "output/phase2/scaled_35types/centroids_human_35.csv"
CENTROIDS_MOUSE = PROJECT / "output/phase2/scaled_35types/centroids_mouse_35.csv"
FULL_RANKING = PROJECT / "output/phase2/scaled_35types/residuals_ranked.csv"
SENSITIVITY_RESULT = PROJECT / "output/landmark_sensitivity/sensitivity_gate_result.json"
OUTPUT_DIR = PROJECT / "output/figures"

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
N_ITERATIONS = 1000
OBSERVED_RHO = 0.8515  # From sensitivity_gate_result.json
RANDOM_SEED = 42
SECONDARY_GENE_COUNTS = [500, 2000]

# ---------------------------------------------------------------------------
# Suppress verbose printing from procrustes module during iterations
# ---------------------------------------------------------------------------
import io
import contextlib


def silent_pipeline(human_full, mouse_full, gene_indices, full_residuals):
    """
    Run CellWarp pipeline on a gene subset and return Spearman ρ.

    Silences all print output from the procrustes module.
    Returns ρ (Spearman correlation of residual magnitudes vs full-space).
    """
    gene_cols = [human_full.columns[i] for i in gene_indices]
    human_sub = human_full[gene_cols]
    mouse_sub = mouse_full[gene_cols]

    with contextlib.redirect_stdout(io.StringIO()):
        human_pca, mouse_pca, pca_model, cell_types = pca_reduce_centroids(
            human_sub, mouse_sub, variance_threshold=0.95
        )
        result = procrustes_align(human_pca, mouse_pca)

    # Compute residual magnitudes
    residual_mags = {}
    for i, ct in enumerate(cell_types):
        r = result.aligned_target[i] - result.centered_reference[i]
        residual_mags[ct] = np.linalg.norm(r)

    # Build DataFrame aligned to full ranking
    sub_df = pd.DataFrame(
        [{"cell_type": ct, "residual_magnitude_sub": mag}
         for ct, mag in residual_mags.items()]
    )

    merged = pd.merge(full_residuals, sub_df, on="cell_type")
    rho, _ = stats.spearmanr(
        merged["residual_magnitude"],
        merged["residual_magnitude_sub"],
    )
    return rho


def run_random_baseline(human_full, mouse_full, full_residuals, n_genes,
                        n_iterations, seed, label=""):
    """
    Run n_iterations of random gene sampling and return array of ρ values.
    """
    rng = np.random.RandomState(seed)
    n_total_genes = human_full.shape[1]
    rhos = np.zeros(n_iterations)

    t0 = time.time()
    for i in range(n_iterations):
        gene_idx = rng.choice(n_total_genes, size=n_genes, replace=False)
        rhos[i] = silent_pipeline(human_full, mouse_full, gene_idx, full_residuals)

        if (i + 1) % 100 == 0:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed
            eta = (n_iterations - i - 1) / rate
            print(f"  {label}Iteration {i+1}/{n_iterations} "
                  f"(ρ={rhos[i]:.3f}, elapsed={elapsed:.0f}s, ETA={eta:.0f}s)")

    elapsed = time.time() - t0
    print(f"  {label}Completed {n_iterations} iterations in {elapsed:.1f}s")
    return rhos


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("L1000 RANDOM BASELINE — Null Distribution for Rigidity Ranking ρ")
    print("=" * 70)

    # Load data
    print("\n--- Loading data ---")
    human_full = pd.read_csv(CENTROIDS_HUMAN, index_col=0)
    mouse_full = pd.read_csv(CENTROIDS_MOUSE, index_col=0)
    full_ranking = pd.read_csv(FULL_RANKING)
    full_residuals = full_ranking[["cell_type", "residual_magnitude"]]

    with open(SENSITIVITY_RESULT) as f:
        sens = json.load(f)

    n_l1000 = sens["landmark_genes_in_centroids"]  # 907
    observed_rho = sens["spearman_rho"]             # 0.8515

    print(f"  Human centroids: {human_full.shape}")
    print(f"  Mouse centroids: {mouse_full.shape}")
    print(f"  Full-space ranking: {len(full_residuals)} cell types")
    print(f"  L1000 gene count in centroids: {n_l1000}")
    print(f"  Observed L1000 ρ: {observed_rho}")
    print(f"  Total ortholog genes: {human_full.shape[1]}")

    # -----------------------------------------------------------------------
    # Primary analysis: 907 genes × 1000 iterations
    # -----------------------------------------------------------------------
    print(f"\n--- Primary analysis: {n_l1000} random genes × "
          f"{N_ITERATIONS} iterations ---")

    rhos_primary = run_random_baseline(
        human_full, mouse_full, full_residuals,
        n_genes=n_l1000,
        n_iterations=N_ITERATIONS,
        seed=RANDOM_SEED,
        label=f"[{n_l1000} genes] ",
    )

    # Compute statistics
    mean_rho = np.mean(rhos_primary)
    std_rho = np.std(rhos_primary)
    p5 = np.percentile(rhos_primary, 5)
    p95 = np.percentile(rhos_primary, 95)
    z_score = (observed_rho - mean_rho) / std_rho if std_rho > 0 else float("inf")
    empirical_p = np.sum(rhos_primary >= observed_rho) / N_ITERATIONS

    # -----------------------------------------------------------------------
    # Secondary analysis: 500 and 2000 genes
    # -----------------------------------------------------------------------
    secondary_results = {}
    for n_genes in SECONDARY_GENE_COUNTS:
        print(f"\n--- Secondary analysis: {n_genes} random genes × "
              f"{N_ITERATIONS} iterations ---")
        rhos_sec = run_random_baseline(
            human_full, mouse_full, full_residuals,
            n_genes=n_genes,
            n_iterations=N_ITERATIONS,
            seed=RANDOM_SEED + n_genes,  # Different seed for independence
            label=f"[{n_genes} genes] ",
        )
        secondary_results[n_genes] = {
            "rhos": rhos_sec,
            "mean": float(np.mean(rhos_sec)),
            "std": float(np.std(rhos_sec)),
            "p5": float(np.percentile(rhos_sec, 5)),
            "p95": float(np.percentile(rhos_sec, 95)),
        }

    # -----------------------------------------------------------------------
    # Report
    # -----------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("RESULTS")
    print("=" * 70)

    print(f"\nInput files:")
    print(f"  Human centroids:    {CENTROIDS_HUMAN}")
    print(f"  Mouse centroids:    {CENTROIDS_MOUSE}")
    print(f"  Full-space ranking: {FULL_RANKING}")
    print(f"  L1000 sensitivity:  {SENSITIVITY_RESULT}")

    print(f"\nPrimary analysis ({n_l1000} genes, {N_ITERATIONS} iterations):")
    print(f"  Random mean ρ:      {mean_rho:.3f}")
    print(f"  Random SD ρ:        {std_rho:.3f}")
    print(f"  5th percentile:     {p5:.3f}")
    print(f"  95th percentile:    {p95:.3f}")
    print(f"  Observed L1000 ρ:   {observed_rho:.3f}")
    print(f"  Z-score:            {z_score:.3f}")
    print(f"  Empirical p-value:  {empirical_p:.3f}")

    print(f"\nSecondary analysis (scaling with gene count):")
    for n_genes in SECONDARY_GENE_COUNTS:
        sr = secondary_results[n_genes]
        print(f"  {n_genes} genes: mean ρ = {sr['mean']:.3f} ± {sr['std']:.3f} "
              f"(5th={sr['p5']:.3f}, 95th={sr['p95']:.3f})")

    # Characterization
    print()
    if empirical_p < 0.05:
        verdict = "is"
        flag = ("L1000 FINDING HOLDS — add random baseline to methods")
    else:
        verdict = "is not"
        flag = ("L1000 FINDING WEAKENED — Language Ruling 3 must be revised, "
                "flagged for review")

    characterization = (
        f"The L1000 ρ={observed_rho:.3f} {verdict} significantly above the "
        f"random baseline (random mean ρ={mean_rho:.3f} ± {std_rho:.3f}, "
        f"empirical p={empirical_p:.3f})"
    )
    print(f"CHARACTERIZATION: {characterization}")
    print(f"\nFLAG: {flag}")

    # -----------------------------------------------------------------------
    # Plot
    # -----------------------------------------------------------------------
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Panel A: Primary histogram (907 genes)
    ax = axes[0]
    ax.hist(rhos_primary, bins=40, color="#4C72B0", edgecolor="white",
            alpha=0.8, density=True)
    ax.axvline(observed_rho, color="red", linewidth=2, linestyle="--",
               label=f"L1000 ρ = {observed_rho:.3f}")
    ax.axvline(mean_rho, color="black", linewidth=1, linestyle="-",
               label=f"Random mean = {mean_rho:.3f}")
    ax.set_xlabel("Spearman ρ (subset vs full-space ranking)", fontsize=11)
    ax.set_ylabel("Density", fontsize=11)
    ax.set_title(f"Random {n_l1000}-gene subsets (n={N_ITERATIONS})", fontsize=12)
    ax.legend(fontsize=9, loc="upper left")

    # Add text box with stats
    stats_text = (
        f"mean={mean_rho:.3f} ± {std_rho:.3f}\n"
        f"z = {z_score:.1f}\n"
        f"p = {empirical_p:.3f}"
    )
    ax.text(0.98, 0.95, stats_text, transform=ax.transAxes, fontsize=9,
            verticalalignment="top", horizontalalignment="right",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow",
                      edgecolor="gray", alpha=0.9))

    # Panel B: 500 genes
    ax = axes[1]
    sr500 = secondary_results[500]
    ax.hist(sr500["rhos"], bins=40, color="#55A868", edgecolor="white",
            alpha=0.8, density=True)
    ax.axvline(observed_rho, color="red", linewidth=2, linestyle="--",
               label=f"L1000 ρ = {observed_rho:.3f}")
    ax.axvline(sr500["mean"], color="black", linewidth=1, linestyle="-",
               label=f"Random mean = {sr500['mean']:.3f}")
    ax.set_xlabel("Spearman ρ", fontsize=11)
    ax.set_title(f"Random 500-gene subsets (n={N_ITERATIONS})", fontsize=12)
    ax.legend(fontsize=9, loc="upper left")

    # Panel C: 2000 genes
    ax = axes[2]
    sr2000 = secondary_results[2000]
    ax.hist(sr2000["rhos"], bins=40, color="#C44E52", edgecolor="white",
            alpha=0.8, density=True)
    ax.axvline(observed_rho, color="red", linewidth=2, linestyle="--",
               label=f"L1000 ρ = {observed_rho:.3f}")
    ax.axvline(sr2000["mean"], color="black", linewidth=1, linestyle="-",
               label=f"Random mean = {sr2000['mean']:.3f}")
    ax.set_xlabel("Spearman ρ", fontsize=11)
    ax.set_title(f"Random 2000-gene subsets (n={N_ITERATIONS})", fontsize=12)
    ax.legend(fontsize=9, loc="upper left")

    fig.suptitle("L1000 Landmark Gene Random Baseline Test", fontsize=14,
                 fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "l1000_random_baseline.pdf",
                bbox_inches="tight", dpi=300)
    fig.savefig(OUTPUT_DIR / "l1000_random_baseline.png",
                bbox_inches="tight", dpi=150)
    print(f"\nPlot saved: {OUTPUT_DIR / 'l1000_random_baseline.pdf'}")

    # Text description of plot
    print(f"\nPlot description: Three-panel figure showing histograms of "
          f"Spearman ρ values from random gene subsets. Panel A ({n_l1000} "
          f"genes, matching L1000 count): distribution centered at "
          f"ρ={mean_rho:.3f}, observed L1000 ρ={observed_rho:.3f} shown as "
          f"red dashed line. Panel B (500 genes): distribution centered at "
          f"ρ={sr500['mean']:.3f}. Panel C (2000 genes): distribution "
          f"centered at ρ={sr2000['mean']:.3f}.")

    # -----------------------------------------------------------------------
    # Save JSON results
    # -----------------------------------------------------------------------
    results_json = {
        "analysis": "L1000 Random Baseline Test",
        "date": "2026-03-18",
        "input_files": {
            "human_centroids": str(CENTROIDS_HUMAN.relative_to(PROJECT)),
            "mouse_centroids": str(CENTROIDS_MOUSE.relative_to(PROJECT)),
            "full_space_ranking": str(FULL_RANKING.relative_to(PROJECT)),
            "sensitivity_gate_result": str(SENSITIVITY_RESULT.relative_to(PROJECT)),
        },
        "primary": {
            "n_genes_sampled": int(n_l1000),
            "n_iterations": int(N_ITERATIONS),
            "n_cell_types": int(len(full_residuals)),
            "total_ortholog_genes": int(human_full.shape[1]),
            "random_seed": int(RANDOM_SEED),
            "random_mean_rho": round(float(mean_rho), 3),
            "random_std_rho": round(float(std_rho), 3),
            "percentile_5": round(float(p5), 3),
            "percentile_95": round(float(p95), 3),
            "observed_l1000_rho": round(float(observed_rho), 3),
            "z_score": round(float(z_score), 3),
            "empirical_p_value": round(float(empirical_p), 3),
            "characterization": characterization,
            "flag": flag,
        },
        "secondary": {
            str(n_genes): {
                "n_genes": int(n_genes),
                "mean_rho": round(float(sr["mean"]), 3),
                "std_rho": round(float(sr["std"]), 3),
                "percentile_5": round(float(sr["p5"]), 3),
                "percentile_95": round(float(sr["p95"]), 3),
            }
            for n_genes, sr in secondary_results.items()
        },
    }

    with open(OUTPUT_DIR / "l1000_random_baseline_results.json", "w") as f:
        json.dump(results_json, f, indent=2)
    print(f"Results saved: {OUTPUT_DIR / 'l1000_random_baseline_results.json'}")

    print("\nDone.")


if __name__ == "__main__":
    main()
