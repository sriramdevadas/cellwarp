#!/usr/bin/env python3
"""
CellWarp — Bootstrap Confidence Intervals on Per-Type Rigidity Rankings

Addresses reviewer criticism of ranking instability by bootstrapping the FULL
Procrustes pipeline 1000 times and computing confidence intervals on each cell
type's residual-magnitude rank.

Biology
-------
The 35-cell-type Procrustes analysis produces per-type residual magnitudes that
rank cell types from most rigid (small residual, well-explained by the global
cross-species transformation) to most flexible (large residual, cell-type-specific
divergence). Reviewers noted that point rankings without uncertainty estimates
could be misleading — small differences in residual magnitude might not reflect
real biological differences.

By resampling cells within each cell type and re-running the full pipeline, we
test how stable these rankings are to sampling variation. Stable rankings
indicate robust biological signals; unstable rankings indicate that apparent
rank differences are within noise.

Math
----
For b = 1..1000 bootstrap iterations:
  1. For each cell type t in each species s: resample n_t cells WITH REPLACEMENT
     from the original n_t cells (classic nonparametric bootstrap).
  2. Recompute centroids: mu_t^(b) = (1/n_b) sum x_i
  3. Rerun PCA on resampled centroids (fresh PCA, not projection onto original
     PC space — see design decision below).
  4. Rerun Procrustes alignment (mouse -> human).
  5. Compute per-type residual magnitudes ||r_i^(b)||.
  6. Rank the 35 types by residual magnitude (1 = most diverged/flexible).
  7. Store the full ranking vector.

DESIGN DECISION — Fresh PCA vs projection:
  We refit PCA each iteration (fresh PCA) rather than projecting onto the
  original PC space. Rationale: the resampled centroids shift positions, so
  the principal axes of variation also shift. Refitting PCA tests the stability
  of the FULL pipeline end-to-end, which is what reviewers care about. Projecting
  onto the original PC space would understate uncertainty by holding the coordinate
  system fixed. Fresh PCA is the more conservative, more honest choice.

Inputs:
    data/phase2_scaled/human_scaled.h5ad
    data/phase2_scaled/mouse_scaled.h5ad

Outputs (all in analysis/bootstrap_rankings/):
    bootstrap_rankings_raw.csv       — 1000 rows x 35 columns (rank per type per iteration)
    bootstrap_summary.csv            — per-type: median rank, 95% CI, SD, category
    bootstrap_results.md             — human-readable summary
    pairwise_swap_matrix.csv         — 35x35 pairwise rank-swap probabilities

Figures (in figures/supplementary/):
    bootstrap_forest_plot.pdf        — Panel A: forest plot of rank CIs
    bootstrap_swap_heatmap.pdf       — Panel B: pairwise swap probability heatmap

Usage:
    python analysis/bootstrap_rankings/bootstrap_ranking_analysis.py
"""

from __future__ import annotations

import contextlib
import io
import sys
import time
import warnings
from multiprocessing import Pool, cpu_count
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
import anndata as ad

from cellwarp.procrustes import (
    PCA_VARIANCE_THRESHOLD,
    _procrustes_distance,
    pca_reduce_centroids,
    procrustes_align,
)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

N_BOOTSTRAP = 1000
RANDOM_SEED = 42
CELL_TYPE_COL = "cell_type"
N_TYPES = 35

# Stability classification thresholds
RANK_TOP_THIRD = 12        # ranks 1-12 = top third (most diverged/flexible)
RANK_BOTTOM_THIRD = 24     # ranks 24-35 = bottom third (most rigid/conserved)
STABLE_CI_WIDTH = 10       # CI width <= 10 = stable
UNSTABLE_CI_WIDTH = 15     # CI width > 15 = unstable

# Paths
HUMAN_DATA = PROJECT_ROOT / "data" / "phase2_scaled" / "human_scaled.h5ad"
MOUSE_DATA = PROJECT_ROOT / "data" / "phase2_scaled" / "mouse_scaled.h5ad"
OUTPUT_DIR = PROJECT_ROOT / "analysis" / "bootstrap_rankings"
FIGURE_DIR = PROJECT_ROOT / "figures" / "supplementary"


# ---------------------------------------------------------------------------
# Bootstrap helpers
# ---------------------------------------------------------------------------


def _resample_and_centroid(
    X_data: np.ndarray,
    cell_type_labels: np.ndarray,
    gene_names: list[str],
    sorted_types: list[str],
    rng: np.random.RandomState,
) -> pd.DataFrame:
    """
    Resample cells with replacement within each type and compute centroids.

    Operates on raw numpy arrays to avoid AnnData pickling overhead in
    multiprocessing. Combines resampling + centroid computation in one pass
    to minimize memory allocation.

    Args:
        X_data: Dense expression matrix (n_cells, n_genes).
        cell_type_labels: Array of cell type strings per cell.
        gene_names: Gene identifiers (column names).
        sorted_types: Sorted list of cell type names.
        rng: Seeded RandomState.

    Returns:
        DataFrame (n_types, n_genes) of resampled centroids.
    """
    centroids = {}
    for ct in sorted_types:
        ct_idx = np.where(cell_type_labels == ct)[0]
        resampled_idx = rng.choice(ct_idx, size=len(ct_idx), replace=True)
        centroids[ct] = X_data[resampled_idx].mean(axis=0)
    df = pd.DataFrame(centroids, index=gene_names).T
    df.index.name = "cell_type"
    return df


# Global variables set by _init_worker for multiprocessing
_worker_human_X = None
_worker_human_labels = None
_worker_mouse_X = None
_worker_mouse_labels = None
_worker_gene_names = None
_worker_sorted_types = None


def _init_worker(h_X, h_labels, m_X, m_labels, gene_names, sorted_types):
    """Initialize worker process with shared data (avoids pickling per task)."""
    global _worker_human_X, _worker_human_labels
    global _worker_mouse_X, _worker_mouse_labels
    global _worker_gene_names, _worker_sorted_types
    _worker_human_X = h_X
    _worker_human_labels = h_labels
    _worker_mouse_X = m_X
    _worker_mouse_labels = m_labels
    _worker_gene_names = gene_names
    _worker_sorted_types = sorted_types


def run_single_bootstrap(i: int) -> dict:
    """
    Execute one bootstrap iteration of the full Procrustes pipeline.

    Uses global worker variables initialized by _init_worker to avoid
    pickling large arrays for each task. Only the iteration index is sent.

    Pipeline per iteration:
      1. Resample cells with replacement within each type
      2. Recompute centroids from resampled cells
      3. Fresh PCA on resampled centroids
      4. Procrustes alignment (mouse -> human)
      5. Compute per-type residual magnitudes
      6. Rank by residual magnitude (1 = most diverged)

    Args:
        i: Bootstrap iteration index.

    Returns:
        Dict with iteration number, per-type rankings, residual magnitudes,
        Procrustes distance, and number of PCA components.
    """
    rng = np.random.RandomState(RANDOM_SEED + i)

    # 1-2. Resample cells and compute centroids (combined for efficiency)
    h_cent = _resample_and_centroid(
        _worker_human_X, _worker_human_labels, _worker_gene_names,
        _worker_sorted_types, rng,
    )
    m_cent = _resample_and_centroid(
        _worker_mouse_X, _worker_mouse_labels, _worker_gene_names,
        _worker_sorted_types, rng,
    )

    # 3. Fresh PCA on resampled centroids (silent)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        with contextlib.redirect_stdout(io.StringIO()):
            h_pca, m_pca, pca_model, cell_types = pca_reduce_centroids(
                h_cent, m_cent, PCA_VARIANCE_THRESHOLD
            )

    # 4. Procrustes alignment (silent — use internal helper logic)
    n, k = h_pca.shape
    X_c = h_pca - h_pca.mean(axis=0)
    Y_c = m_pca - m_pca.mean(axis=0)

    M = X_c.T @ Y_c
    U, sigma, Vt = np.linalg.svd(M)
    V = Vt.T

    d = np.linalg.det(V @ U.T)
    D_diag = np.ones(k)
    D_diag[-1] = np.sign(d)
    D = np.diag(D_diag)

    R = V @ D @ U.T
    ss_Y = np.sum(Y_c ** 2)
    trace_sigma_D = np.sum(sigma * D_diag)
    s = trace_sigma_D / ss_Y

    Y_aligned = s * (Y_c @ R)

    proc_dist = float(np.sqrt(np.sum((X_c - Y_aligned) ** 2)))

    # 5. Per-type residual magnitudes
    residual_mags = {}
    for j, ct in enumerate(cell_types):
        r = Y_aligned[j] - X_c[j]
        residual_mags[ct] = float(np.linalg.norm(r))

    # 6. Rank by residual magnitude (1 = most diverged/flexible)
    sorted_types = sorted(residual_mags.keys(), key=lambda x: residual_mags[x], reverse=True)
    rankings = {ct: rank + 1 for rank, ct in enumerate(sorted_types)}

    return {
        "iteration": i,
        "rankings": rankings,
        "residual_mags": residual_mags,
        "procrustes_distance": proc_dist,
        "n_pca_components": pca_model.n_components_,
    }


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------


def classify_stability(row: pd.Series) -> str:
    """
    Classify a cell type's ranking stability based on its bootstrap CI.

    Categories:
      STABLE_FLEXIBLE: 95% CI entirely within top third (ranks 1-12)
                       These types are reliably the most diverged.
      STABLE_RIGID:    95% CI entirely within bottom third (ranks 24-35)
                       These types are reliably the most conserved.
      STABLE_MIDDLE:   CI width <= 10 ranks, centered in middle third
      UNSTABLE:        CI width > 15 ranks, or spans flexible-to-rigid

    Note on naming convention: in the residual ranking, rank 1 = LARGEST
    residual = most FLEXIBLE (diverged). Rank 35 = smallest residual =
    most RIGID (conserved). This matches the D'Arcy Thompson framing:
    "rigid" types are well-explained by the global transformation.

    Args:
        row: Series with ci_lower, ci_upper, ci_width fields.

    Returns:
        String classification label.
    """
    lower = row["ci_lower"]
    upper = row["ci_upper"]
    width = row["ci_width"]

    # CI entirely in top third (ranks 1-12) = stably flexible/diverged
    if upper <= RANK_TOP_THIRD:
        return "STABLE_FLEXIBLE"
    # CI entirely in bottom third (ranks 24-35) = stably rigid/conserved
    if lower >= RANK_BOTTOM_THIRD:
        return "STABLE_RIGID"
    # Narrow CI in middle
    if width <= STABLE_CI_WIDTH:
        return "STABLE_MIDDLE"
    # Wide CI = unstable
    if width > UNSTABLE_CI_WIDTH:
        return "UNSTABLE"
    # Moderate width, not clearly in one region
    if width <= UNSTABLE_CI_WIDTH:
        return "MODERATE"
    return "UNSTABLE"


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def plot_forest(summary_df: pd.DataFrame, output_path: Path) -> str:
    """
    Panel A: Forest plot of bootstrap ranking CIs per cell type.

    Each cell type is shown as a point (median rank) with horizontal bar
    (95% CI). Types are ordered by median rank on the y-axis. Color-coded
    by stability category.

    Args:
        summary_df: DataFrame with columns median_rank, ci_lower, ci_upper, category.
        output_path: Path to save the figure.

    Returns:
        Text description of the plot for terminal output.
    """
    # Sort by median rank (most flexible at top)
    df = summary_df.sort_values("median_rank").reset_index()

    category_colors = {
        "STABLE_FLEXIBLE": "#d62728",   # red — reliably diverged
        "STABLE_RIGID": "#2ca02c",      # green — reliably conserved
        "STABLE_MIDDLE": "#1f77b4",     # blue — stable middle
        "MODERATE": "#ff7f0e",          # orange — moderate uncertainty
        "UNSTABLE": "#7f7f7f",          # gray — unstable
    }

    fig, ax = plt.subplots(figsize=(10, 14))

    for idx, row in df.iterrows():
        color = category_colors.get(row["category"], "#7f7f7f")
        y = len(df) - idx - 1  # flip so rank 1 is at top

        # CI bar
        ax.plot(
            [row["ci_lower"], row["ci_upper"]], [y, y],
            color=color, linewidth=2.5, solid_capstyle="round", zorder=2,
        )
        # Median point
        ax.plot(
            row["median_rank"], y,
            "o", color=color, markersize=7, zorder=3,
            markeredgecolor="white", markeredgewidth=0.5,
        )
        # Label
        ax.text(
            -0.5, y, row["cell_type"],
            ha="right", va="center", fontsize=8,
        )

    # Reference lines for thirds
    ax.axvline(RANK_TOP_THIRD + 0.5, color="gray", linestyle=":", alpha=0.5, linewidth=1)
    ax.axvline(RANK_BOTTOM_THIRD - 0.5, color="gray", linestyle=":", alpha=0.5, linewidth=1)
    ax.text(
        RANK_TOP_THIRD / 2, len(df) + 0.5, "FLEXIBLE\n(diverged)",
        ha="center", va="bottom", fontsize=8, color="#d62728", alpha=0.7,
    )
    ax.text(
        (RANK_BOTTOM_THIRD + N_TYPES) / 2, len(df) + 0.5, "RIGID\n(conserved)",
        ha="center", va="bottom", fontsize=8, color="#2ca02c", alpha=0.7,
    )

    ax.set_xlim(0, N_TYPES + 1)
    ax.set_ylim(-1, len(df) + 1.5)
    ax.set_xlabel("Residual Magnitude Rank (1 = most diverged, 35 = most conserved)", fontsize=11)
    ax.set_yticks([])
    ax.set_title(
        "Bootstrap 95% CI on Per-Type Rigidity Rankings\n"
        f"(n = {N_BOOTSTRAP} bootstrap iterations, full pipeline re-run)",
        fontsize=12, fontweight="bold",
    )

    # Legend
    legend_handles = [
        mpatches.Patch(color=c, label=l.replace("_", " ").title())
        for l, c in category_colors.items()
    ]
    ax.legend(
        handles=legend_handles, loc="lower right", fontsize=9,
        framealpha=0.9, title="Stability Category",
    )

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Count categories for description
    cats = summary_df["category"].value_counts()
    desc = (
        f"Forest plot: 35 cell types ordered by median bootstrap rank. "
        f"Horizontal bars show 95% CIs from {N_BOOTSTRAP} bootstrap iterations. "
        f"Categories: {', '.join(f'{k}={v}' for k, v in sorted(cats.items()))}. "
        f"Dashed lines mark rank boundaries between thirds."
    )
    return desc


def plot_swap_heatmap(swap_matrix: pd.DataFrame, output_path: Path) -> str:
    """
    Panel B: Heatmap of pairwise rank-swap probability.

    For each pair of cell types (i, j), shows the fraction of bootstraps
    where type i was ranked higher (more diverged) than type j. Values near
    0.5 indicate the ordering is a coin flip; values near 0 or 1 indicate
    reliable pairwise ordering.

    Args:
        swap_matrix: 35x35 DataFrame of swap probabilities.
        output_path: Path to save the figure.

    Returns:
        Text description of the plot.
    """
    # Order by original rank (median bootstrap rank)
    n = len(swap_matrix)

    # Custom colormap: blue (0) -> white (0.5) -> red (1)
    cmap = LinearSegmentedColormap.from_list(
        "swap", ["#2166ac", "#f7f7f7", "#b2182b"], N=256
    )

    fig, ax = plt.subplots(figsize=(14, 12))

    im = ax.imshow(swap_matrix.values, cmap=cmap, vmin=0, vmax=1, aspect="equal")

    # Labels
    ax.set_xticks(range(n))
    ax.set_xticklabels(swap_matrix.columns, rotation=90, fontsize=6.5, ha="center")
    ax.set_yticks(range(n))
    ax.set_yticklabels(swap_matrix.index, fontsize=6.5)

    # Colorbar
    cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label("P(row type ranked more diverged than column type)", fontsize=10)

    ax.set_title(
        "Pairwise Rank-Swap Probability Matrix\n"
        f"(n = {N_BOOTSTRAP} bootstrap iterations)",
        fontsize=12, fontweight="bold",
    )

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Count coin-flip pairs (0.4 < p < 0.6)
    upper_tri = swap_matrix.values[np.triu_indices(n, k=1)]
    n_coin_flip = int(np.sum((upper_tri > 0.4) & (upper_tri < 0.6)))
    n_total_pairs = len(upper_tri)
    n_reliable = int(np.sum((upper_tri > 0.9) | (upper_tri < 0.1)))

    desc = (
        f"Heatmap: {n}x{n} pairwise rank-swap probabilities. "
        f"Of {n_total_pairs} pairs, {n_reliable} ({n_reliable/n_total_pairs*100:.0f}%) "
        f"have reliable ordering (P>0.9 or P<0.1), "
        f"{n_coin_flip} ({n_coin_flip/n_total_pairs*100:.0f}%) are coin flips (0.4<P<0.6)."
    )
    return desc


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    t_start = time.time()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("BOOTSTRAP RANKING CONFIDENCE INTERVALS (35 CELL TYPES)")
    print("=" * 70)

    # ------------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------------
    print("\n  Loading 35-type normalized data...")
    human = ad.read_h5ad(HUMAN_DATA)
    mouse = ad.read_h5ad(MOUSE_DATA)
    print(f"  Human: {human.n_obs:,} cells x {human.n_vars:,} genes")
    print(f"  Mouse: {mouse.n_obs:,} cells x {mouse.n_vars:,} genes")
    n_types = len(human.obs[CELL_TYPE_COL].unique())
    print(f"  Cell types: {n_types}")

    # ------------------------------------------------------------------
    # Bootstrap loop
    # ------------------------------------------------------------------
    print(f"\n  Running {N_BOOTSTRAP} bootstrap iterations (full pipeline)...")
    print(f"  Seed: {RANDOM_SEED}, resampling: WITH REPLACEMENT")
    print(f"  PCA: fresh per iteration (not projecting onto original PC space)")

    # Pre-extract dense matrices and labels to avoid pickling AnnData
    # (AnnData pickling is extremely slow — ~1.4GB per object per task)
    print("  Pre-extracting expression matrices (dense)...")
    import scipy.sparse as sp
    human_X = human.X.toarray() if sp.issparse(human.X) else np.asarray(human.X)
    mouse_X = mouse.X.toarray() if sp.issparse(mouse.X) else np.asarray(mouse.X)
    human_labels = human.obs[CELL_TYPE_COL].values
    mouse_labels = mouse.obs[CELL_TYPE_COL].values
    gene_names = human.var_names.tolist()
    sorted_types = sorted(human.obs[CELL_TYPE_COL].unique())
    print(f"  Human: {human_X.shape}, Mouse: {mouse_X.shape}")

    # Free AnnData objects to save memory
    del human, mouse

    n_workers = max(1, cpu_count() - 1)
    print(f"  Workers: {n_workers} (multiprocessing)")

    # Use multiprocessing with initializer to share data across workers
    # (each worker gets a copy via fork, avoiding per-task pickling)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)

        try:
            with Pool(
                processes=n_workers,
                initializer=_init_worker,
                initargs=(human_X, human_labels, mouse_X, mouse_labels,
                          gene_names, sorted_types),
            ) as pool:
                results = []
                for j, result in enumerate(pool.imap_unordered(run_single_bootstrap, range(N_BOOTSTRAP))):
                    results.append(result)
                    if (j + 1) % 50 == 0 or j == 0:
                        elapsed = time.time() - t_start
                        rate = (j + 1) / elapsed
                        eta = (N_BOOTSTRAP - j - 1) / rate if rate > 0 else 0
                        print(
                            f"    [{j + 1:>4}/{N_BOOTSTRAP}] "
                            f"d={result['procrustes_distance']:.4f}, "
                            f"PC={result['n_pca_components']}, "
                            f"{elapsed:.0f}s elapsed, ~{eta:.0f}s remaining"
                        )
        except Exception as e:
            # Fallback to sequential if multiprocessing fails
            print(f"  Multiprocessing failed ({e}), falling back to sequential...")
            # Set globals for sequential execution
            _init_worker(human_X, human_labels, mouse_X, mouse_labels,
                         gene_names, sorted_types)
            results = []
            for j in range(N_BOOTSTRAP):
                result = run_single_bootstrap(j)
                results.append(result)
                if (j + 1) % 50 == 0 or j == 0:
                    elapsed = time.time() - t_start
                    rate = (j + 1) / elapsed
                    eta = (N_BOOTSTRAP - j - 1) / rate if rate > 0 else 0
                    print(
                        f"    [{j + 1:>4}/{N_BOOTSTRAP}] "
                        f"d={result['procrustes_distance']:.4f}, "
                        f"PC={result['n_pca_components']}, "
                        f"{elapsed:.0f}s elapsed, ~{eta:.0f}s remaining"
                    )

    t_bootstrap = time.time() - t_start
    print(f"\n  Bootstrap complete: {t_bootstrap:.0f}s ({t_bootstrap/60:.1f} min)")

    # ------------------------------------------------------------------
    # Assemble ranking matrix
    # ------------------------------------------------------------------
    print("\n  Assembling ranking matrix...")

    # Get cell type order from first result
    cell_types = sorted(results[0]["rankings"].keys())

    # Build rankings matrix: (N_BOOTSTRAP x N_TYPES)
    ranking_matrix = np.zeros((N_BOOTSTRAP, len(cell_types)), dtype=int)
    residual_matrix = np.zeros((N_BOOTSTRAP, len(cell_types)))
    distances = np.zeros(N_BOOTSTRAP)
    n_pcs = np.zeros(N_BOOTSTRAP, dtype=int)

    for res in results:
        i = res["iteration"]
        for j, ct in enumerate(cell_types):
            ranking_matrix[i, j] = res["rankings"][ct]
            residual_matrix[i, j] = res["residual_mags"][ct]
        distances[i] = res["procrustes_distance"]
        n_pcs[i] = res["n_pca_components"]

    # Save raw rankings
    raw_df = pd.DataFrame(ranking_matrix, columns=cell_types)
    raw_df.index.name = "iteration"
    raw_df.to_csv(OUTPUT_DIR / "bootstrap_rankings_raw.csv")
    print(f"  Saved: bootstrap_rankings_raw.csv ({raw_df.shape})")

    # ------------------------------------------------------------------
    # Compute summary statistics per cell type
    # ------------------------------------------------------------------
    print("\n  Computing summary statistics...")

    # Load original rankings once (not inside the loop)
    original_rank_df = pd.read_csv(
        PROJECT_ROOT / "output" / "phase2" / "scaled_35types" / "residuals_ranked.csv",
        index_col=0,
    )

    summary_rows = []
    for j, ct in enumerate(cell_types):
        ranks = ranking_matrix[:, j]
        mags = residual_matrix[:, j]

        median_rank = float(np.median(ranks))
        mean_rank = float(np.mean(ranks))
        ci_lower = float(np.percentile(ranks, 2.5))
        ci_upper = float(np.percentile(ranks, 97.5))
        ci_width = ci_upper - ci_lower
        sd_rank = float(np.std(ranks))
        frac_top10 = float(np.mean(ranks <= 10))
        frac_bottom10 = float(np.mean(ranks >= 26))

        orig_row = original_rank_df[original_rank_df["cell_type"] == ct]
        original_rank = int(orig_row.index[0]) if len(orig_row) > 0 else -1

        summary_rows.append({
            "cell_type": ct,
            "original_rank": original_rank,
            "median_rank": median_rank,
            "mean_rank": mean_rank,
            "ci_lower": ci_lower,
            "ci_upper": ci_upper,
            "ci_width": ci_width,
            "sd_rank": sd_rank,
            "frac_top10_rigid": frac_top10,
            "frac_bottom10_flexible": frac_bottom10,
            "mean_residual_mag": float(np.mean(mags)),
            "sd_residual_mag": float(np.std(mags)),
        })

    summary_df = pd.DataFrame(summary_rows)

    # Classify stability
    summary_df["category"] = summary_df.apply(classify_stability, axis=1)

    # Sort by median rank
    summary_df = summary_df.sort_values("median_rank").reset_index(drop=True)

    summary_df.to_csv(OUTPUT_DIR / "bootstrap_summary.csv", index=False)
    print(f"  Saved: bootstrap_summary.csv")

    # ------------------------------------------------------------------
    # Pairwise swap matrix
    # ------------------------------------------------------------------
    print("\n  Computing pairwise rank-swap probabilities...")

    # Order types by median rank for display
    ordered_types = summary_df["cell_type"].tolist()
    type_to_col = {ct: j for j, ct in enumerate(cell_types)}

    swap_matrix = np.zeros((len(ordered_types), len(ordered_types)))
    for a_idx, ct_a in enumerate(ordered_types):
        col_a = type_to_col[ct_a]
        for b_idx, ct_b in enumerate(ordered_types):
            col_b = type_to_col[ct_b]
            # P(type a ranked more diverged = smaller rank number = larger residual)
            swap_matrix[a_idx, b_idx] = float(
                np.mean(ranking_matrix[:, col_a] < ranking_matrix[:, col_b])
            )

    swap_df = pd.DataFrame(swap_matrix, index=ordered_types, columns=ordered_types)
    swap_df.to_csv(OUTPUT_DIR / "pairwise_swap_matrix.csv")
    print(f"  Saved: pairwise_swap_matrix.csv")

    # ------------------------------------------------------------------
    # Generate figures
    # ------------------------------------------------------------------
    print("\n  Generating figures...")

    forest_desc = plot_forest(
        summary_df,
        FIGURE_DIR / "bootstrap_forest_plot.pdf",
    )
    print(f"  Saved: figures/supplementary/bootstrap_forest_plot.pdf")
    print(f"  {forest_desc}")

    swap_desc = plot_swap_heatmap(
        swap_df,
        FIGURE_DIR / "bootstrap_swap_heatmap.pdf",
    )
    print(f"  Saved: figures/supplementary/bootstrap_swap_heatmap.pdf")
    print(f"  {swap_desc}")

    # ------------------------------------------------------------------
    # Classification summary
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("STABILITY CLASSIFICATION")
    print("=" * 70)

    category_counts = summary_df["category"].value_counts()
    for cat in ["STABLE_FLEXIBLE", "STABLE_RIGID", "STABLE_MIDDLE", "MODERATE", "UNSTABLE"]:
        count = category_counts.get(cat, 0)
        types_in_cat = summary_df[summary_df["category"] == cat]["cell_type"].tolist()
        print(f"\n  {cat}: {count} types")
        for ct in types_in_cat:
            row = summary_df[summary_df["cell_type"] == ct].iloc[0]
            print(
                f"    {ct:<50} median={row['median_rank']:.0f}  "
                f"CI=[{row['ci_lower']:.0f}, {row['ci_upper']:.0f}]  "
                f"width={row['ci_width']:.0f}"
            )

    n_stable = int((summary_df["ci_width"] <= STABLE_CI_WIDTH).sum())
    n_unstable = int((summary_df["ci_width"] > UNSTABLE_CI_WIDTH).sum())

    # ------------------------------------------------------------------
    # Summary table
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("PER-TYPE RANKING SUMMARY (ordered by median rank)")
    print("=" * 70)

    print(f"\n  {'Cell Type':<50} {'Orig':>5} {'Med':>5} {'CI':>12} {'W':>4} {'SD':>5} {'Cat':<18}")
    print(f"  {'-' * 102}")
    for _, row in summary_df.iterrows():
        print(
            f"  {row['cell_type']:<50} {row['original_rank']:>5} "
            f"{row['median_rank']:>5.0f} "
            f"[{row['ci_lower']:>4.0f}, {row['ci_upper']:>4.0f}] "
            f"{row['ci_width']:>4.0f} {row['sd_rank']:>5.1f} "
            f"{row['category']:<18}"
        )

    # ------------------------------------------------------------------
    # Write human-readable results
    # ------------------------------------------------------------------
    md_path = OUTPUT_DIR / "bootstrap_results.md"
    with open(md_path, "w") as f:
        f.write("# Bootstrap Ranking Confidence Intervals\n\n")
        f.write(f"**Date:** 2026-04-05\n")
        f.write(f"**Bootstrap iterations:** {N_BOOTSTRAP}\n")
        f.write(f"**Random seed:** {RANDOM_SEED}\n")
        f.write(f"**Resampling:** With replacement, same n per type\n")
        f.write(f"**PCA:** Fresh per iteration (not projecting onto original PC space)\n")
        f.write(f"**Runtime:** {time.time() - t_start:.0f}s\n\n")

        f.write("## Key Finding\n\n")
        f.write(f"Of 35 cell types:\n")
        f.write(f"- **{n_stable}** have stable rankings (CI width <= {STABLE_CI_WIDTH})\n")
        f.write(f"- **{n_unstable}** have unstable rankings (CI width > {UNSTABLE_CI_WIDTH})\n\n")

        f.write("## Stability Categories\n\n")
        for cat in ["STABLE_FLEXIBLE", "STABLE_RIGID", "STABLE_MIDDLE", "MODERATE", "UNSTABLE"]:
            count = category_counts.get(cat, 0)
            types_in_cat = summary_df[summary_df["category"] == cat]["cell_type"].tolist()
            f.write(f"### {cat} ({count} types)\n\n")
            if types_in_cat:
                for ct in types_in_cat:
                    row = summary_df[summary_df["cell_type"] == ct].iloc[0]
                    f.write(
                        f"- **{ct}**: median rank {row['median_rank']:.0f}, "
                        f"95% CI [{row['ci_lower']:.0f}, {row['ci_upper']:.0f}], "
                        f"width {row['ci_width']:.0f}\n"
                    )
            else:
                f.write("(none)\n")
            f.write("\n")

        f.write("## Design Decisions\n\n")
        f.write("1. **Fresh PCA per iteration** (not projecting onto original PC space): "
                "The resampled centroids shift positions, so the principal axes of variation "
                "also shift. Refitting PCA tests the stability of the FULL pipeline end-to-end. "
                "Projecting onto the original PC space would understate uncertainty by holding "
                "the coordinate system fixed.\n\n")
        f.write("2. **Resampling with replacement** (classic nonparametric bootstrap): "
                "This tests sensitivity to which specific cells are sampled. The expected "
                "fraction of unique cells per type per iteration is ~63.2%.\n\n")
        f.write("3. **Same n per type**: Each bootstrap iteration preserves the original "
                "cell count per type, avoiding confounding sample-size effects with "
                "biological variation.\n\n")

        f.write("## Procrustes Distance Stability\n\n")
        f.write(f"- Mean distance: {np.mean(distances):.4f}\n")
        f.write(f"- SD: {np.std(distances):.4f}\n")
        f.write(f"- CV: {np.std(distances)/np.mean(distances):.4f}\n")
        f.write(f"- Range: [{np.min(distances):.4f}, {np.max(distances):.4f}]\n\n")

        f.write("## PCA Components\n\n")
        f.write(f"- Mean components: {np.mean(n_pcs):.1f}\n")
        f.write(f"- Range: [{np.min(n_pcs)}, {np.max(n_pcs)}]\n\n")

        f.write("## Full Ranking Table\n\n")
        f.write("| Cell Type | Orig Rank | Median | 95% CI | Width | SD | Category |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        for _, row in summary_df.iterrows():
            f.write(
                f"| {row['cell_type']} | {row['original_rank']} | "
                f"{row['median_rank']:.0f} | "
                f"[{row['ci_lower']:.0f}, {row['ci_upper']:.0f}] | "
                f"{row['ci_width']:.0f} | {row['sd_rank']:.1f} | "
                f"{row['category']} |\n"
            )

        f.write("\n## Pairwise Ordering Reliability\n\n")
        upper_tri = swap_df.values[np.triu_indices(len(swap_df), k=1)]
        n_reliable = int(np.sum((upper_tri > 0.9) | (upper_tri < 0.1)))
        n_coin_flip = int(np.sum((upper_tri > 0.4) & (upper_tri < 0.6)))
        n_total = len(upper_tri)
        f.write(f"- Total pairs: {n_total}\n")
        f.write(f"- Reliable ordering (P>0.9 or P<0.1): {n_reliable} ({n_reliable/n_total*100:.0f}%)\n")
        f.write(f"- Coin flip (0.4<P<0.6): {n_coin_flip} ({n_coin_flip/n_total*100:.0f}%)\n")

    print(f"\n  Saved: bootstrap_results.md")

    # ------------------------------------------------------------------
    # Final summary
    # ------------------------------------------------------------------
    t_total = time.time() - t_start

    print("\n" + "=" * 70)
    print("BOOTSTRAP RANKING ANALYSIS — COMPLETE")
    print("=" * 70)

    print(f"\n  KEY FINDING: {n_stable}/{N_TYPES} cell types have stable rankings (CI width <= {STABLE_CI_WIDTH})")
    print(f"               {n_unstable}/{N_TYPES} cell types are unstable (CI width > {UNSTABLE_CI_WIDTH})")

    print(f"\n  Procrustes distance: mean={np.mean(distances):.4f}, CV={np.std(distances)/np.mean(distances):.4f}")
    print(f"  PCA components: mean={np.mean(n_pcs):.1f}")

    print(f"\n  Output files:")
    for f in sorted(OUTPUT_DIR.iterdir()):
        if f.is_file():
            size = f.stat().st_size
            if size > 1024 * 1024:
                s = f"{size / 1024 / 1024:.1f} MB"
            elif size > 1024:
                s = f"{size / 1024:.1f} KB"
            else:
                s = f"{size} B"
            print(f"    {f.name:<40} {s:>10}")

    print(f"\n  Figures:")
    print(f"    figures/supplementary/bootstrap_forest_plot.pdf")
    print(f"    figures/supplementary/bootstrap_swap_heatmap.pdf")

    print(f"\n  Total runtime: {t_total:.0f}s ({t_total/60:.1f} min)")


if __name__ == "__main__":
    main()
