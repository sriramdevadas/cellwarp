"""
SAMap validation module for CellWarp Phase 1.

Runs SAMap (Tarashansky et al. 2021, eLife) to independently identify cell type
correspondences between human and mouse single-cell atlases, then compares
SAMap's pairings to our manually curated 6-type pairing list.

Biology: SAMap uses a gene homology graph (normally from BLAST) plus mutual
nearest-neighbor manifold alignment to identify cross-species cell type mappings.
We provide a pre-computed gene homology graph built from BioMart 1:1 orthologs
instead of BLAST, because our datasets already share the same ortholog-aligned
gene space.

Math: SAMap's alignment score for a cell type pair (A, B) is the average
cross-species edge weight between cells of type A in species 1 and cells of
type B in species 2 after iterative manifold alignment. Higher scores indicate
stronger correspondence.
"""

import os
from pathlib import Path
from typing import Optional

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
import seaborn as sns
from samap.analysis import get_mapping_scores
from samap.mapping import SAMAP


def _apply_compat_patches():
    """
    Apply all compatibility patches needed for SAMap to work with modern
    numpy (>=2.0), scipy (>=1.14), and anndata (>=0.10).

    SAMap 1.0.x was written for numpy <2, scipy <1.14, anndata <0.10.
    Three known incompatibilities:

    1. np.in1d removed in numpy 2.0 (replaced by np.isin)
    2. Sparse matrix .A property removed in scipy 1.14
    3. anndata >=0.10 rejects LIL matrices in .obsp
    """
    # 1. numpy: restore np.in1d
    if not hasattr(np, "in1d"):
        np.in1d = np.isin

    # 2. scipy: restore .A property on all sparse matrix classes
    for cls in (sp.csr_matrix, sp.csc_matrix, sp.coo_matrix, sp.lil_matrix,
                sp.bsr_matrix, sp.dia_matrix, sp.dok_matrix):
        if not hasattr(cls, "A"):
            cls.A = property(lambda self: self.toarray())

    # 3. scipy sparse: redirect float16 to float32 (unsupported in modern scipy)
    for cls in (sp.csr_matrix, sp.csc_matrix, sp.coo_matrix):
        if not hasattr(cls, "_orig_astype"):
            cls._orig_astype = cls.astype

            def _patched_astype(self, dtype, casting="unsafe", copy=True,
                                _orig=cls._orig_astype):
                dt = np.dtype(dtype)
                if dt == np.float16:
                    dt = np.float32
                return _orig(self, dt, casting=casting, copy=copy)

            cls.astype = _patched_astype

    # 4. pandas: restore Series.nonzero() for scipy sparse indexing
    if not hasattr(pd.Series, "nonzero"):
        pd.Series.nonzero = lambda self: (self.values.nonzero()[0],)

    # 4. anndata: auto-convert LIL/COO to CSR in obsp
    try:
        from anndata._core.aligned_mapping import PairwiseArrays

        _orig_validate = PairwiseArrays._validate_value

        def _patched_validate(self, val, key):
            if sp.issparse(val) and not isinstance(val, (sp.csr_matrix, sp.csc_matrix)):
                val = val.tocsr()
            return _orig_validate(self, val, key)

        PairwiseArrays._validate_value = _patched_validate
    except (ImportError, AttributeError):
        pass  # anndata version doesn't need this patch

# Our 6 manually curated cell type pairings (human <-> mouse).
# Both species use the same cell_type labels after ortholog alignment.
MANUAL_PAIRINGS = [
    ("B cell", "B cell"),
    ("CD4-positive, alpha-beta T cell", "CD4-positive, alpha-beta T cell"),
    ("CD8-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell"),
    ("endothelial cell", "endothelial cell"),
    ("hepatocyte", "hepatocyte"),
    ("macrophage", "macrophage"),
]


def build_ortholog_gnnm(
    human_adata: ad.AnnData,
    mouse_adata: ad.AnnData,
    hu_id: str = "hu",
    mo_id: str = "mo",
) -> tuple:
    """
    Build a gene-gene homology graph (gnnm) from 1:1 orthologs.

    Since both datasets are already aligned to the same gene space (human
    Ensembl IDs), each gene in the human dataset maps 1:1 to the same-indexed
    gene in the mouse dataset. We represent this as a bipartite graph where
    each ortholog pair has edge weight 1.0.

    Parameters
    ----------
    human_adata : AnnData
        Human dataset with Ensembl IDs as var_names.
    mouse_adata : AnnData
        Mouse dataset with the same Ensembl IDs as var_names.
    hu_id : str
        Species prefix for human genes.
    mo_id : str
        Species prefix for mouse genes.

    Returns
    -------
    gnnm : scipy.sparse.csr_matrix
        Sparse adjacency matrix of gene-gene homology.
    gns : np.ndarray
        Array of all gene names (hu_* then mo_*).
    gns_dict : dict
        Dict mapping species IDs to their gene name arrays.
    """
    # Get shared genes (should be identical since both are ortholog-aligned)
    shared_genes = np.intersect1d(human_adata.var_names, mouse_adata.var_names)
    n_genes = len(shared_genes)
    print(f"  Building gnnm from {n_genes} shared 1:1 orthologs")

    # Create prefixed gene names
    hu_genes = np.array([f"{hu_id}_{g}" for g in shared_genes])
    mo_genes = np.array([f"{mo_id}_{g}" for g in shared_genes])
    gns = np.concatenate([hu_genes, mo_genes])

    # Build bipartite adjacency: hu_gene_i <-> mo_gene_i with weight 1.0
    total = 2 * n_genes
    rows = np.concatenate([np.arange(n_genes), np.arange(n_genes, total)])
    cols = np.concatenate([np.arange(n_genes, total), np.arange(n_genes)])
    vals = np.ones(2 * n_genes, dtype=np.float64)
    gnnm = sp.csr_matrix((vals, (rows, cols)), shape=(total, total))

    gns_dict = {hu_id: hu_genes, mo_id: mo_genes}
    return gnnm, gns, gns_dict


def run_samap(
    human_path: str,
    mouse_path: str,
    hu_id: str = "hu",
    mo_id: str = "mo",
    n_iters: int = 3,
    ncpus: Optional[int] = None,
) -> SAMAP:
    """
    Run SAMap alignment between human and mouse datasets.

    Loads raw count data, builds 1:1 ortholog gnnm, runs SAMap's iterative
    manifold alignment, and returns the fitted SAMAP object.

    Parameters
    ----------
    human_path : str
        Path to human .h5ad file (raw counts, ortholog-aligned).
    mouse_path : str
        Path to mouse .h5ad file (raw counts, ortholog-aligned).
    hu_id, mo_id : str
        Species identifiers for SAMap.
    n_iters : int
        Number of SAMap iterations (default 3).
    ncpus : int, optional
        Number of CPUs. Defaults to os.cpu_count().

    Returns
    -------
    sm : SAMAP
        Fitted SAMAP object with alignment results.
    """
    if ncpus is None:
        ncpus = os.cpu_count() or 4

    # Patch SAMap compatibility with modern numpy/scipy/anndata
    _apply_compat_patches()

    # Load data to build gnnm
    print("Loading data for gnnm construction...")
    hu_ad = ad.read_h5ad(human_path)
    mo_ad = ad.read_h5ad(mouse_path)
    gnnm_tuple = build_ortholog_gnnm(hu_ad, mo_ad, hu_id, mo_id)
    del hu_ad, mo_ad  # free memory

    # Initialize SAMAP — it will run SAM preprocessing on each dataset
    print("Initializing SAMap (SAM preprocessing per species)...")
    sm = SAMAP(
        sams={hu_id: human_path, mo_id: mouse_path},
        gnnm=gnnm_tuple,
        keys={hu_id: "cell_type", mo_id: "cell_type"},
        save_processed=False,
    )

    # Run iterative alignment
    print(f"Running SAMap alignment ({n_iters} iterations)...")
    sm.run(NUMITERS=n_iters, ncpus=ncpus, umap=False)

    return sm


def get_cell_type_scores(sm: SAMAP, hu_id: str = "hu", mo_id: str = "mo") -> pd.DataFrame:
    """
    Extract cell type mapping scores from a fitted SAMAP object.

    Returns a DataFrame where rows are human cell types and columns are mouse
    cell types, with values being SAMap alignment scores (higher = stronger
    correspondence).

    Parameters
    ----------
    sm : SAMAP
        Fitted SAMAP object.
    hu_id, mo_id : str
        Species identifiers.

    Returns
    -------
    scores_df : pd.DataFrame
        Pairwise mapping scores (human rows x mouse columns).
    """
    keys = {hu_id: "cell_type", mo_id: "cell_type"}
    D, A = get_mapping_scores(sm, keys=keys, n_top=0)
    # A is a full pairwise table. Filter to human-row x mouse-col
    hu_types = [c for c in A.index if c.startswith(f"{hu_id}_")]
    mo_types = [c for c in A.columns if c.startswith(f"{mo_id}_")]
    scores = A.loc[hu_types, mo_types].copy()

    # Strip species prefixes for readability
    scores.index = [x.replace(f"{hu_id}_", "") for x in scores.index]
    scores.columns = [x.replace(f"{mo_id}_", "") for x in scores.columns]

    return scores


def compare_pairings(
    scores_df: pd.DataFrame,
    manual_pairings: Optional[list] = None,
) -> dict:
    """
    Compare SAMap's mapping scores against our manual cell type pairings.

    For each human cell type, SAMap's top-scoring mouse cell type is compared
    to our manual pairing. A pairing is 'confirmed' if SAMap's top hit matches
    our assignment.

    Parameters
    ----------
    scores_df : pd.DataFrame
        Pairwise mapping scores (human x mouse).
    manual_pairings : list of (str, str)
        List of (human_type, mouse_type) manual pairings.

    Returns
    -------
    results : dict with keys:
        'n_confirmed' : int — number of confirmed pairings
        'n_total' : int — total manual pairings
        'details' : list of dict — per-pairing details
        'all_confirmed' : bool — True if all pairings confirmed
    """
    if manual_pairings is None:
        manual_pairings = MANUAL_PAIRINGS

    details = []
    n_confirmed = 0

    for hu_type, mo_type in manual_pairings:
        if hu_type not in scores_df.index:
            details.append({
                "human_type": hu_type,
                "mouse_type_manual": mo_type,
                "mouse_type_samap": "N/A",
                "score_manual": np.nan,
                "score_top": np.nan,
                "confirmed": False,
                "note": "Human type not found in SAMap scores",
            })
            continue

        row = scores_df.loc[hu_type]
        top_mouse = row.idxmax()
        top_score = row.max()
        manual_score = row.get(mo_type, np.nan)
        confirmed = top_mouse == mo_type

        if confirmed:
            n_confirmed += 1

        # Rank of manual pairing
        rank = int((row >= manual_score).sum()) if not np.isnan(manual_score) else -1

        details.append({
            "human_type": hu_type,
            "mouse_type_manual": mo_type,
            "mouse_type_samap": top_mouse,
            "score_manual": float(manual_score),
            "score_top": float(top_score),
            "rank": rank,
            "confirmed": confirmed,
        })

    return {
        "n_confirmed": n_confirmed,
        "n_total": len(manual_pairings),
        "fraction_confirmed": n_confirmed / len(manual_pairings),
        "details": details,
        "all_confirmed": n_confirmed == len(manual_pairings),
    }


def plot_mapping_heatmap(
    scores_df: pd.DataFrame,
    manual_pairings: Optional[list] = None,
    output_path: Optional[str] = None,
) -> plt.Figure:
    """
    Plot a heatmap of SAMap cell type mapping scores.

    Manual pairings are highlighted with a black border on the diagonal.

    Parameters
    ----------
    scores_df : pd.DataFrame
        Pairwise mapping scores.
    manual_pairings : list of (str, str), optional
        Manual pairings to highlight.
    output_path : str, optional
        Path to save figure.

    Returns
    -------
    fig : matplotlib Figure
    """
    if manual_pairings is None:
        manual_pairings = MANUAL_PAIRINGS

    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(
        scores_df,
        annot=True,
        fmt=".3f",
        cmap="YlOrRd",
        ax=ax,
        linewidths=0.5,
        square=True,
    )
    ax.set_title("SAMap Cell Type Mapping Scores\n(Human rows x Mouse columns)")
    ax.set_xlabel("Mouse cell types")
    ax.set_ylabel("Human cell types")

    # Highlight manual pairings
    for hu_type, mo_type in manual_pairings:
        if hu_type in scores_df.index and mo_type in scores_df.columns:
            i = list(scores_df.index).index(hu_type)
            j = list(scores_df.columns).index(mo_type)
            ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=False, edgecolor="black", lw=3))

    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"  Saved heatmap to {output_path}")
    return fig


def save_results(
    scores_df: pd.DataFrame,
    comparison: dict,
    output_dir: str,
    manual_pairings: Optional[list] = None,
) -> None:
    """
    Save all SAMap validation outputs: scores CSV, comparison report, heatmap.

    Parameters
    ----------
    scores_df : pd.DataFrame
        Pairwise mapping scores.
    comparison : dict
        Output of compare_pairings().
    output_dir : str
        Directory to save outputs.
    manual_pairings : list, optional
        Manual pairings for heatmap annotation.
    """
    os.makedirs(output_dir, exist_ok=True)

    # 1. Save raw scores
    scores_path = os.path.join(output_dir, "samap_mapping_scores.csv")
    scores_df.to_csv(scores_path)
    print(f"  Saved mapping scores to {scores_path}")

    # 2. Save comparison report
    report_path = os.path.join(output_dir, "samap_comparison_report.txt")
    with open(report_path, "w") as f:
        f.write("SAMap Validation Report — CellWarp Phase 1\n")
        f.write("=" * 50 + "\n\n")
        f.write(
            f"Confirmed: {comparison['n_confirmed']}/{comparison['n_total']} "
            f"({comparison['fraction_confirmed']:.0%})\n\n"
        )
        f.write("Per-pairing details:\n")
        f.write("-" * 50 + "\n")
        for d in comparison["details"]:
            status = "CONFIRMED" if d["confirmed"] else "MISMATCH"
            f.write(f"\n  Human: {d['human_type']}\n")
            f.write(f"  Manual mouse pairing: {d['mouse_type_manual']}\n")
            f.write(f"  SAMap top match:      {d['mouse_type_samap']}\n")
            f.write(f"  Score (manual pair):  {d['score_manual']:.4f}\n")
            f.write(f"  Score (SAMap top):    {d['score_top']:.4f}\n")
            rank = d.get("rank", "N/A")
            f.write(f"  Rank of manual pair:  {rank}\n")
            f.write(f"  Status: {status}\n")

        f.write("\n" + "=" * 50 + "\n")
        if comparison["all_confirmed"]:
            f.write("RESULT: ALL PAIRINGS CONFIRMED — Phase 1 SAMap gate PASSED\n")
        else:
            f.write(
                f"RESULT: {comparison['n_confirmed']}/{comparison['n_total']} confirmed "
                "— review mismatches before proceeding\n"
            )
    print(f"  Saved comparison report to {report_path}")

    # 3. Save heatmap
    heatmap_path = os.path.join(output_dir, "samap_heatmap.png")
    plot_mapping_heatmap(scores_df, manual_pairings, output_path=heatmap_path)
    plt.close("all")

    # 4. Save comparison as CSV
    details_path = os.path.join(output_dir, "samap_pairing_details.csv")
    pd.DataFrame(comparison["details"]).to_csv(details_path, index=False)
    print(f"  Saved pairing details to {details_path}")
