#!/usr/bin/env python3
"""
CellWarp — Cancer Procrustes Pipeline (Thread 1)

Applies the same Procrustes framework used for cross-species comparison to
normal colon vs colorectal cancer (CRC) tumor tissue. The condition variable
is disease (normal vs tumor) instead of organism (human vs mouse).

Biology
-------
Cancer reshapes the transcriptomic landscape of the colon. By treating normal
and tumor cell type centroids as two "constellations" in gene expression space,
Procrustes analysis asks: does the normal-to-tumor transformation follow a
coherent geometric deformation? Which cell types are most deformed? And
critically — do evolutionarily rigid cell types (from cross-species analysis)
also resist tumor reprogramming?

Math
----
Identical to cross-species pipeline (src/procrustes.py):
  min_{R, s, t}  ||X_normal - (s * X_tumor * R + 1*t^T)||_F^2

Per-donor centroids are computed first, then averaged across donors. This
controls for donor imbalance between conditions (232 normal vs 62 tumor donors).

Steps
-----
1. Load data (colon_normal.h5ad, colon_tumor.h5ad)
2. Compute per-donor centroids, then condition-level centroids
3. PCA on combined centroids (95% variance)
4. Procrustes alignment (tumor → normal)
5. Permutation test (10,000 iterations)
6. Per-cell-type deformation scores
7. Top genes per cell type via PCA loading projection
8. Cross-analysis correlation with cross-species residuals
9. Summary print
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from cellwarp.procrustes import (
    RANDOM_SEED,
    N_PERMUTATIONS,
    PCA_VARIANCE_THRESHOLD,
    N_TOP_GENES,
    procrustes_align,
    permutation_test,
    compute_residual_vectors,
    map_residuals_to_genes,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DATA_DIR = Path("data/cancer")
OUTPUT_DIR = Path("output/cancer")
CROSS_SPECIES_RESULTS = Path("output/phase2/scaled_35types/procrustes_results_35.json")
ORTHOLOG_GENES_PATH = Path("data/phase1/human_aligned.h5ad")

# Exclude "Other" — not a coherent cell type for geometric analysis
EXCLUDE_TYPES = {"Other"}

# Minimum cells per coarse type per condition
MIN_CELLS = 500


# ---------------------------------------------------------------------------
# Step 1 — Load data
# ---------------------------------------------------------------------------


def load_data() -> tuple:
    """Load normal and tumor h5ad files. Print cell counts and cell type breakdown."""
    import anndata as ad

    print("=" * 70)
    print("STEP 1 — Load Data")
    print("=" * 70)

    normal = ad.read_h5ad(DATA_DIR / "colon_normal.h5ad")
    tumor = ad.read_h5ad(DATA_DIR / "colon_tumor.h5ad")

    print(f"\n  Normal: {normal.n_obs:,} cells × {normal.n_vars:,} genes")
    print(f"  Tumor:  {tumor.n_obs:,} cells × {tumor.n_vars:,} genes")

    # Confirm gene spaces match
    assert normal.n_vars == tumor.n_vars, "Gene count mismatch"
    assert (normal.var.index == tumor.var.index).all(), "Gene order mismatch"
    print(f"  Gene space: identical ({normal.n_vars:,} genes)")

    # Check overlap with cross-species ortholog space
    xs = ad.read_h5ad(ORTHOLOG_GENES_PATH)
    xs_genes = set(xs.var.index)
    cancer_genes = set(normal.var.index)
    overlap = xs_genes & cancer_genes
    print(f"\n  Cross-species ortholog genes: {len(xs_genes):,}")
    print(f"  Cancer genes: {len(cancer_genes):,}")
    print(f"  Overlap: {len(overlap):,} ({len(overlap)/len(xs_genes)*100:.1f}%)")
    del xs

    # Filter cancer data to cross-species gene space for comparability
    shared_genes = sorted(overlap)
    normal = normal[:, normal.var.index.isin(shared_genes)].copy()
    tumor = tumor[:, tumor.var.index.isin(shared_genes)].copy()
    print(f"  After filtering to shared ortholog space: {normal.n_vars:,} genes")

    # Cell type breakdown
    for condition, adata in [("Normal", normal), ("Tumor", tumor)]:
        print(f"\n  {condition} — coarse cell type breakdown:")
        counts = adata.obs["coarse_cell_type"].value_counts()
        for ct, n in counts.items():
            marker = " [EXCLUDE]" if ct in EXCLUDE_TYPES else ""
            print(f"    {ct:<25} {n:>6,}{marker}")

    return normal, tumor, shared_genes


# ---------------------------------------------------------------------------
# Step 2 — Compute centroids (per-donor, then average)
# ---------------------------------------------------------------------------


def compute_donor_centroids(
    normal, tumor, shared_genes: list[str]
) -> tuple[pd.DataFrame, pd.DataFrame, list[str], list[str]]:
    """
    Compute per-cell-type centroids using donor-level averaging.

    Biology: Averaging per-donor first, then across donors, prevents a single
    donor with many cells from dominating the centroid. This is critical when
    donor counts differ between conditions (232 normal vs 62 tumor).

    Math: For cell type t in condition c:
      μ_{t,d} = mean(cells from donor d of type t)  [per-donor centroid]
      μ_{t,c} = mean(μ_{t,d} for all donors d)      [condition centroid]
    """
    print("\n" + "=" * 70)
    print("STEP 2 — Compute Per-Donor Centroids")
    print("=" * 70)

    # Identify valid cell types (≥500 in BOTH conditions, not excluded)
    normal_counts = normal.obs["coarse_cell_type"].value_counts()
    tumor_counts = tumor.obs["coarse_cell_type"].value_counts()

    all_types = sorted(
        set(normal_counts.index) & set(tumor_counts.index) - EXCLUDE_TYPES
    )

    passed = []
    dropped = []
    for ct in all_types:
        n_norm = normal_counts.get(ct, 0)
        n_tum = tumor_counts.get(ct, 0)
        if n_norm >= MIN_CELLS and n_tum >= MIN_CELLS:
            passed.append(ct)
        else:
            dropped.append(ct)

    print(f"\n  Cell types passing ≥{MIN_CELLS} gate in both conditions ({len(passed)}):")
    for ct in passed:
        print(f"    {ct:<25} normal={normal_counts[ct]:>5,}  tumor={tumor_counts[ct]:>5,}")

    if dropped:
        print(f"\n  Dropped ({len(dropped)}):")
        for ct in dropped:
            n = normal_counts.get(ct, 0)
            t = tumor_counts.get(ct, 0)
            print(f"    {ct:<25} normal={n:>5,}  tumor={t:>5,}")
    else:
        print(f"  No cell types dropped.")

    # Compute per-donor centroids for each condition
    normal_centroids = _donor_averaged_centroids(normal, passed, shared_genes, "Normal")
    tumor_centroids = _donor_averaged_centroids(tumor, passed, shared_genes, "Tumor")

    # Save
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    normal_centroids.to_csv(OUTPUT_DIR / "centroids_normal.csv")
    tumor_centroids.to_csv(OUTPUT_DIR / "centroids_tumor.csv")
    print(f"\n  Saved: {OUTPUT_DIR / 'centroids_normal.csv'}")
    print(f"  Saved: {OUTPUT_DIR / 'centroids_tumor.csv'}")

    return normal_centroids, tumor_centroids, passed, dropped


def _donor_averaged_centroids(
    adata, cell_types: list[str], gene_names: list[str], condition_label: str
) -> pd.DataFrame:
    """Compute mean-of-donor-means centroids for each cell type."""
    print(f"\n  {condition_label} — donor-averaged centroids:")
    centroids = {}

    for ct in cell_types:
        mask = adata.obs["coarse_cell_type"] == ct
        ct_data = adata[mask]
        donors = ct_data.obs["donor_id"].unique()
        n_donors = len(donors)

        # Compute per-donor centroids
        donor_centroids = []
        for donor in donors:
            donor_mask = ct_data.obs["donor_id"] == donor
            donor_cells = ct_data[donor_mask]
            if donor_cells.n_obs > 0:
                mean_vec = np.asarray(donor_cells.X.mean(axis=0)).flatten()
                donor_centroids.append(mean_vec)

        # Average across donors
        centroid = np.mean(donor_centroids, axis=0)
        centroids[ct] = centroid
        n_cells = mask.sum()
        print(
            f"    {ct:<25} {n_cells:>5,} cells, {n_donors:>4} donors → "
            f"centroid ({len(centroid):,} genes)"
        )

    df = pd.DataFrame(centroids, index=gene_names).T
    df.index.name = "cell_type"
    return df


# ---------------------------------------------------------------------------
# Step 3 — PCA on combined centroids
# ---------------------------------------------------------------------------


def run_pca(
    normal_centroids: pd.DataFrame,
    tumor_centroids: pd.DataFrame,
    cell_types: list[str],
) -> tuple[np.ndarray, np.ndarray, PCA]:
    """
    PCA on stacked normal + tumor centroids. Retain 95% variance.

    Math: Stack (2n × G) matrix, fit PCA. With 2n points in G-dim space,
    max rank = 2n - 1. PCA captures dominant axes of variation across
    cell types and conditions simultaneously.
    """
    print("\n" + "=" * 70)
    print("STEP 3 — PCA on Combined Centroids")
    print("=" * 70)

    normal_mat = normal_centroids.loc[cell_types].values
    tumor_mat = tumor_centroids.loc[cell_types].values
    combined = np.vstack([normal_mat, tumor_mat])

    n_max = min(combined.shape[0] - 1, combined.shape[1])

    pca = PCA(
        n_components=PCA_VARIANCE_THRESHOLD,
        svd_solver="full",
        random_state=RANDOM_SEED,
    )
    combined_pca = pca.fit_transform(combined)
    n_components = pca.n_components_

    n_types = len(cell_types)
    normal_pca = combined_pca[:n_types]
    tumor_pca = combined_pca[n_types:]

    cumvar = np.cumsum(pca.explained_variance_ratio_)
    print(f"\n  PCA on {2 * n_types} centroids ({combined.shape[1]:,} genes):")
    print(f"  Components retained: {n_components} / {n_max} possible")
    print(f"  Cumulative variance explained: {cumvar[-1] * 100:.1f}%")
    for i, (var, cum) in enumerate(zip(pca.explained_variance_ratio_, cumvar)):
        print(f"    PC{i + 1}: {var * 100:.2f}%  (cumulative: {cum * 100:.1f}%)")

    # Save PCA
    np.savez(
        OUTPUT_DIR / "pca_cancer.npz",
        normal_pca=normal_pca,
        tumor_pca=tumor_pca,
        components=pca.components_,
        explained_variance_ratio=pca.explained_variance_ratio_,
        mean=pca.mean_,
    )
    print(f"\n  Saved: {OUTPUT_DIR / 'pca_cancer.npz'}")

    return normal_pca, tumor_pca, pca


# ---------------------------------------------------------------------------
# Step 4 — Procrustes alignment
# ---------------------------------------------------------------------------


def run_procrustes(
    normal_pca: np.ndarray, tumor_pca: np.ndarray
):
    """Align tumor onto normal using Procrustes (identical to cross-species)."""
    print("\n" + "=" * 70)
    print("STEP 4 — Procrustes Alignment (tumor → normal)")
    print("=" * 70)

    result = procrustes_align(normal_pca, tumor_pca)

    # Verify proper rotation
    det = np.linalg.det(result.rotation)
    assert abs(det - 1.0) < 1e-6, f"Rotation determinant = {det}, expected +1.0"
    print(f"  Rotation determinant verified: +1.0 (proper rotation)")

    return result


# ---------------------------------------------------------------------------
# Step 5 — Permutation test
# ---------------------------------------------------------------------------


def run_permutation_test(
    normal_pca: np.ndarray, tumor_pca: np.ndarray
) -> tuple[float, np.ndarray]:
    """10,000-iteration permutation test for Procrustes significance."""
    print("\n" + "=" * 70)
    print("STEP 5 — Permutation Test (10,000 iterations)")
    print("=" * 70)

    p_value, null_dist = permutation_test(
        normal_pca, tumor_pca, n_permutations=N_PERMUTATIONS, seed=RANDOM_SEED
    )

    # Save null distribution
    np.save(OUTPUT_DIR / "null_distribution.npy", null_dist)
    print(f"\n  Saved: {OUTPUT_DIR / 'null_distribution.npy'}")

    # Compute obs/null ratio
    from cellwarp.procrustes import _procrustes_distance
    observed = _procrustes_distance(normal_pca, tumor_pca)
    null_median = np.median(null_dist)
    ratio = observed / null_median
    print(f"\n  Obs/null ratio: {ratio:.3f}")

    return p_value, null_dist


# ---------------------------------------------------------------------------
# Step 6 — Deformation scores
# ---------------------------------------------------------------------------


def compute_deformation_scores(
    result, cell_types: list[str]
) -> dict[str, float]:
    """
    Deformation score = ||aligned_tumor_centroid - normal_centroid|| per cell type.

    Biology: Larger deformation scores indicate cell types whose transcriptomic
    programs are more disrupted by cancer. The ranking reveals which cell types
    in the tumor microenvironment undergo the most dramatic reprogramming.
    """
    print("\n" + "=" * 70)
    print("STEP 6 — Per-Cell-Type Deformation Scores")
    print("=" * 70)

    residuals = compute_residual_vectors(result, cell_types)

    # Rank by deformation score
    scores = {ct: float(np.linalg.norm(residuals[ct])) for ct in cell_types}
    ranked = sorted(scores.items(), key=lambda x: x[1], reverse=True)

    print(f"\n  RANKED DEFORMATION TABLE:")
    print(f"  {'Rank':<6} {'Cell Type':<25} {'Deformation Score':>18}")
    print(f"  {'-' * 52}")
    for i, (ct, score) in enumerate(ranked, 1):
        print(f"  {i:<6} {ct:<25} {score:>18.4f}")

    return scores, residuals


# ---------------------------------------------------------------------------
# Step 7 — Top genes per cell type
# ---------------------------------------------------------------------------


def get_top_genes(
    residuals: dict, pca_model: PCA, gene_names: list[str]
) -> dict:
    """Project residuals back to gene space and report top 20 genes."""
    print("\n" + "=" * 70)
    print("STEP 7 — Top Genes per Cell Type (PCA Loading Projection)")
    print("=" * 70)

    top_genes = map_residuals_to_genes(
        residuals, pca_model, gene_names, n_top=N_TOP_GENES
    )
    return top_genes


# ---------------------------------------------------------------------------
# Step 8 — Cross-analysis correlation
# ---------------------------------------------------------------------------


def cross_analysis_correlation(
    cancer_scores: dict[str, float],
) -> dict:
    """
    Correlate cancer deformation scores with cross-species residual magnitudes.

    Biology: If evolutionarily rigid cell types (small cross-species residual)
    also resist tumor reprogramming (small cancer deformation), it suggests a
    shared constraint — perhaps chromatin architecture or regulatory network
    topology — that limits both evolutionary and oncogenic rewiring.

    Math: Spearman rank correlation between cancer deformation score and
    cross-species residual magnitude for matched cell types.
    """
    print("\n" + "=" * 70)
    print("STEP 8 — Cross-Analysis Correlation (CRITICAL)")
    print("=" * 70)

    # Load cross-species results
    with open(CROSS_SPECIES_RESULTS) as f:
        xs_data = json.load(f)

    xs_residuals = {
        ct: xs_data["residuals"][ct]["magnitude"]
        for ct in xs_data["cell_types"]
    }

    # Map cancer coarse types to cross-species type names
    # Cancer types are coarse categories; cross-species are fine-grained Census labels
    type_mapping = {
        "B cell": "B cell",
        "Endothelial cell": "endothelial cell",
        "Epithelial cell": "epithelial cell",
        "Fibroblast": "fibroblast",
        "Macrophage": "macrophage",
        "Mast cell": None,  # no mast cell in 35-type set
        "NK cell": "natural killer cell",
        "Smooth muscle cell": "smooth muscle cell",
        "T cell": "T cell",
    }

    # Build matched pairs
    matched = []
    print(f"\n  Matching cancer coarse types to cross-species types:")
    for cancer_type, xs_type in sorted(type_mapping.items()):
        if cancer_type not in cancer_scores:
            continue
        if xs_type is None:
            print(f"    {cancer_type:<25} → [no match in cross-species]")
            continue
        if xs_type not in xs_residuals:
            print(f"    {cancer_type:<25} → {xs_type} [NOT FOUND in 35-type set]")
            continue

        cancer_score = cancer_scores[cancer_type]
        xs_score = xs_residuals[xs_type]
        matched.append({
            "cancer_type": cancer_type,
            "xs_type": xs_type,
            "cancer_deformation": cancer_score,
            "xs_residual": xs_score,
        })
        print(
            f"    {cancer_type:<25} → {xs_type:<30} "
            f"cancer={cancer_score:.4f}  xs={xs_score:.4f}"
        )

    if len(matched) < 3:
        print(f"\n  WARNING: Only {len(matched)} matched types — too few for correlation.")
        return {"n_matched": len(matched), "insufficient": True}

    # Spearman correlation
    cancer_vals = [m["cancer_deformation"] for m in matched]
    xs_vals = [m["xs_residual"] for m in matched]

    rho, p_value = stats.spearmanr(cancer_vals, xs_vals)

    print(f"\n  Spearman correlation (n={len(matched)}):")
    print(f"    ρ = {rho:.4f}")
    print(f"    p = {p_value:.4f}")
    print(f"    Significant at α=0.05: {'YES' if p_value < 0.05 else 'NO'}")

    # Interpretation
    print(f"\n  INTERPRETATION:")
    if p_value < 0.05 and rho > 0:
        print(
            "    POSITIVE CORRELATION: Cell types with large cross-species divergence\n"
            "    also show large cancer deformation. Evolutionary flexibility and\n"
            "    oncogenic plasticity share a common axis."
        )
    elif p_value < 0.05 and rho < 0:
        print(
            "    NEGATIVE CORRELATION: Evolutionarily rigid cell types also resist\n"
            "    tumor reprogramming. A shared constraint (chromatin? regulatory\n"
            "    architecture?) limits both evolutionary and oncogenic rewiring."
        )
    else:
        print(
            "    NO SIGNIFICANT CORRELATION: Evolutionary rigidity and cancer\n"
            "    deformation appear to be independent axes. Cross-species divergence\n"
            "    does not predict susceptibility to tumor reprogramming."
        )

    return {
        "n_matched": len(matched),
        "matched_types": matched,
        "spearman_rho": float(rho),
        "spearman_p": float(p_value),
        "significant_005": bool(p_value < 0.05),
    }


# ---------------------------------------------------------------------------
# Save all results
# ---------------------------------------------------------------------------


def save_all_results(
    result,
    p_value: float,
    null_dist: np.ndarray,
    scores: dict,
    residuals: dict,
    top_genes: dict,
    cell_types: list[str],
    pca_model: PCA,
    cross_corr: dict,
    dropped_types: list[str],
) -> None:
    """Save comprehensive results JSON."""
    print("\n  Saving results...")

    from cellwarp.procrustes import _procrustes_distance
    observed = _procrustes_distance(
        result.centered_reference,
        # Reconstruct unaligned tumor in centered coords
        result.centered_reference + np.array([residuals[ct] for ct in cell_types]),
    )
    # Use the actual observed distance
    observed = result.distance

    null_median = float(np.median(null_dist))

    results_dict = {
        "procrustes": {
            "distance": float(result.distance),
            "distance_squared": float(result.distance_squared),
            "scaling": float(result.scaling),
            "rotation_det": float(np.linalg.det(result.rotation)),
        },
        "permutation_test": {
            "p_value": float(p_value),
            "n_permutations": len(null_dist),
            "null_median": null_median,
            "obs_null_ratio": float(result.distance / null_median),
            "significant_at_001": bool(p_value < 0.01),
        },
        "pca": {
            "n_components": int(pca_model.n_components_),
            "cumulative_variance": float(
                np.sum(pca_model.explained_variance_ratio_)
            ),
            "per_component_variance": pca_model.explained_variance_ratio_.tolist(),
        },
        "cell_types": cell_types,
        "dropped_types": dropped_types,
        "deformation_scores": {
            ct: float(scores[ct]) for ct in cell_types
        },
        "deformation_ranking": [
            {"rank": i + 1, "cell_type": ct, "score": float(scores[ct])}
            for i, (ct, _) in enumerate(
                sorted(scores.items(), key=lambda x: x[1], reverse=True)
            )
        ],
        "residuals": {
            ct: {
                "vector_pca": residuals[ct].tolist(),
                "magnitude": float(np.linalg.norm(residuals[ct])),
            }
            for ct in cell_types
        },
        "top_genes_per_cell_type": {
            ct: top_genes[ct][["gene", "loading", "abs_loading", "rank"]]
            .to_dict(orient="records")
            for ct in cell_types
        },
        "cross_analysis_correlation": cross_corr,
        "random_seed": RANDOM_SEED,
    }

    output_path = OUTPUT_DIR / "cancer_procrustes_results.json"
    tmp_path = output_path.with_suffix(".tmp")
    with open(tmp_path, "w") as f:
        json.dump(results_dict, f, indent=2)
    tmp_path.rename(output_path)
    print(f"  Saved: {output_path}")


# ---------------------------------------------------------------------------
# Step 9 — Summary
# ---------------------------------------------------------------------------


def print_summary(
    result,
    p_value: float,
    null_dist: np.ndarray,
    scores: dict,
    cell_types: list[str],
    cross_corr: dict,
    dropped_types: list[str],
    pca_model: PCA,
) -> None:
    """Print analyst-format summary of all results."""
    print("\n" + "=" * 70)
    print("CANCER PROCRUSTES PIPELINE — SUMMARY")
    print("=" * 70)

    null_median = float(np.median(null_dist))
    obs_null_ratio = result.distance / null_median

    print(f"""
1. WHAT WAS DONE
   Applied Procrustes analysis to normal colon vs CRC tumor tissue.
   {len(cell_types)} coarse cell types used as geometric landmarks.
   Per-donor centroids averaged to control for donor imbalance.
   Dropped types: {dropped_types if dropped_types else 'none'}.
   Gene space: {pca_model.mean_.shape[0]:,} genes (shared ortholog space).

2. KEY NUMBERS
   Procrustes distance:     {result.distance:.4f}
   Scaling factor:          {result.scaling:.6f}
   SSR:                     {result.distance_squared:.4f}
   Permutation p-value:     {p_value:.6f}
   Null median:             {null_median:.4f}
   Obs/null ratio:          {obs_null_ratio:.3f}
   PCA components retained: {pca_model.n_components_}
   Significant at α=0.01:   {'YES' if p_value < 0.01 else 'NO'}

3. RANKED DEFORMATION TABLE""")

    ranked = sorted(scores.items(), key=lambda x: x[1], reverse=True)
    print(f"   {'Rank':<6} {'Cell Type':<25} {'Score':>10}")
    print(f"   {'-' * 44}")
    for i, (ct, score) in enumerate(ranked, 1):
        print(f"   {i:<6} {ct:<25} {score:>10.4f}")

    print(f"""
4. CROSS-ANALYSIS CORRELATION
   Matched types:           {cross_corr.get('n_matched', 'N/A')}
   Spearman ρ:              {cross_corr.get('spearman_rho', 'N/A')}
   p-value:                 {cross_corr.get('spearman_p', 'N/A')}""")

    if cross_corr.get("spearman_p") is not None:
        rho = cross_corr["spearman_rho"]
        p = cross_corr["spearman_p"]
        print(f"""
5. BIOLOGICAL INTERPRETATION""")
        if p < 0.05 and rho > 0:
            print(
                "   Cell types that diverge most between species ALSO deform most in\n"
                "   cancer. This suggests a shared axis of transcriptomic plasticity:\n"
                "   cell types with flexible expression programs are vulnerable to both\n"
                "   evolutionary drift and oncogenic reprogramming."
            )
        elif p < 0.05 and rho < 0:
            print(
                "   Evolutionarily rigid cell types resist tumor reprogramming.\n"
                "   A shared constraint — likely chromatin architecture or regulatory\n"
                "   network topology — limits both evolutionary and oncogenic rewiring."
            )
        else:
            print(
                "   Evolutionary rigidity and cancer deformation are INDEPENDENT.\n"
                "   Cross-species divergence does not predict vulnerability to tumor\n"
                "   reprogramming. Cancer deformation is driven by TME-specific\n"
                "   pressures, not intrinsic plasticity constraints."
            )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    print("\n" + "#" * 70)
    print("# CellWarp — Cancer Procrustes Pipeline (Thread 1)")
    print("#" * 70 + "\n")

    # Step 1 — Load data
    normal, tumor, shared_genes = load_data()

    # Step 2 — Compute centroids
    normal_centroids, tumor_centroids, cell_types, dropped = (
        compute_donor_centroids(normal, tumor, shared_genes)
    )

    # Free memory
    del normal, tumor

    # Step 3 — PCA
    normal_pca, tumor_pca, pca_model = run_pca(
        normal_centroids, tumor_centroids, cell_types
    )

    # Step 4 — Procrustes alignment
    result = run_procrustes(normal_pca, tumor_pca)

    # Step 5 — Permutation test
    p_value, null_dist = run_permutation_test(normal_pca, tumor_pca)

    # Step 6 — Deformation scores
    scores, residuals = compute_deformation_scores(result, cell_types)

    # Step 7 — Top genes
    gene_names = list(normal_centroids.columns)

    # Need feature_name mapping for readable gene names
    import anndata as ad
    ref = ad.read_h5ad(DATA_DIR / "colon_normal.h5ad")
    ensembl_to_name = dict(zip(ref.var["feature_id"], ref.var["feature_name"]))
    del ref

    # Map gene names: use symbols where available, else keep Ensembl ID
    readable_genes = [ensembl_to_name.get(g, g) for g in gene_names]
    top_genes = get_top_genes(residuals, pca_model, readable_genes)

    # Step 8 — Cross-analysis correlation
    cross_corr = cross_analysis_correlation(scores)

    # Save all results
    save_all_results(
        result, p_value, null_dist, scores, residuals, top_genes,
        cell_types, pca_model, cross_corr, dropped,
    )

    # Step 9 — Summary
    print_summary(
        result, p_value, null_dist, scores, cell_types,
        cross_corr, dropped, pca_model,
    )

    print("\n  All outputs saved to output/cancer/")
    print("  Done.\n")


if __name__ == "__main__":
    main()
