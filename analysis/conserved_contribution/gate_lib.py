"""
Conserved-contribution gate — shared library.

Implements the pre-registered analysis plan
(docs/preregistration_conserved_contribution_2026-06-05.md):
  - core quantity C_g = per-gene human-vs-mouse Pearson r across 35 centroids
    (non-circular; computed directly on centroids, independent of PCA/Procrustes)
  - expression-matched background sampler
  - conserved/divergent set definitions
  - curated lineage-defining TF positive-control list

Reuses cellwarp.procrustes routines for the secondary geometry attribution.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

SCALED35 = ROOT / "output/phase2/scaled_35types"
ORTHOLOGS = ROOT / "data/phase1/orthologs_human_mouse.csv"
CELLMARKER_H = ROOT / "data/validation/cellmarker/cellmarker_human_filtered.csv"
TF_ACTIVITY = ROOT / "output/phase2/mechanistic/tf_complexity/tf_activity_human.csv"
GENE_CONS_TABLE = ROOT / "analysis/gene_conservation/gene_conservation_table.csv"

SEED = 42

# Pre-registered canonical lineage-defining / master TFs (see the analysis plan).
POSITIVE_CONTROL_TFS = sorted(set([
    # T / NKT / thymocyte
    "TCF7", "LEF1", "GATA3", "TBX21", "RUNX3", "FOXP3", "ETS1", "BCL11B", "ZEB1", "TOX",
    # NK
    "EOMES", "ID2", "NFIL3",
    # B / plasma
    "PAX5", "EBF1", "POU2AF1", "SPIB", "IRF4", "PRDM1", "XBP1", "BACH2",
    # macrophage / monocyte / DC / granulocyte / neutrophil / microglia / myeloid
    "SPI1", "CEBPA", "CEBPB", "CEBPE", "CEBPD", "MAFB", "IRF8", "KLF4", "BATF3",
    # HSC / hematopoietic precursor
    "GATA1", "GATA2", "TAL1", "RUNX1", "GFI1B", "MEIS1", "HLF", "KLF1",
    # hepatocyte
    "HNF4A", "FOXA1", "FOXA2", "ONECUT1", "NR1H4", "HNF1A",
    # endothelial / vein endothelial
    "ERG", "FLI1", "SOX17", "SOX18", "KLF2", "FOXF1",
    # smooth muscle / cardiac myocyte / fibroblast / stromal / mesenchymal / adventitial
    "MYOCD", "SRF", "GATA4", "NKX2-5", "TBX5", "MEF2C", "TCF21", "PRRX1", "TWIST2",
    # epithelial / basal / urothelial / luminal mammary / enterocyte / goblet / acinar / ductal
    "TP63", "GRHL2", "ELF3", "EHF", "KLF5", "CDX2", "HNF1B", "SPDEF", "PTF1A",
    "RBPJL", "FOXA3", "GATA6", "SOX9", "ASCL2",
]))


def load_centroids():
    """35-type human & mouse centroids, aligned to shared 16,959 genes."""
    h = pd.read_csv(SCALED35 / "centroids_human_35.csv", index_col=0)
    m = pd.read_csv(SCALED35 / "centroids_mouse_35.csv", index_col=0)
    types = sorted(set(h.index) & set(m.index))
    genes = [g for g in h.columns if g in set(m.columns)]
    return h.loc[types, genes], m.loc[types, genes]


def ortholog_maps():
    """Return (sym->ensembl, ensembl->sym) for human, one2one only."""
    o = pd.read_csv(ORTHOLOGS)
    o = o[o.orthology_type == "ortholog_one2one"]
    ens2sym = dict(zip(o.human_ensembl_id, o.human_gene_name))
    sym2ens = {}
    for ens, sym in zip(o.human_ensembl_id, o.human_gene_name):
        sym2ens.setdefault(sym, ens)  # first wins (1:1 anyway)
    return sym2ens, ens2sym


def per_gene_corr(H, M, method="pearson"):
    """Per-gene cross-species correlation across cell types (rows=types)."""
    n = H.shape[1]
    r = np.full(n, np.nan)
    if method == "pearson":
        for j in range(n):
            a, b = H[:, j], M[:, j]
            if np.std(a) > 0 and np.std(b) > 0:
                r[j] = np.corrcoef(a, b)[0, 1]
    else:
        for j in range(n):
            a, b = H[:, j], M[:, j]
            if np.std(a) > 0 and np.std(b) > 0:
                r[j] = stats.spearmanr(a, b)[0]
    return r


def build_gene_table():
    """Core table: gene_id, symbol, C_pearson, C_spearman, mean_expression."""
    h, m = load_centroids()
    genes = list(h.columns)
    H, M = h.values, m.values
    cp = per_gene_corr(H, M, "pearson")
    csp = per_gene_corr(H, M, "spearman")
    _, ens2sym = ortholog_maps()
    df = pd.DataFrame({
        "gene_id": genes,
        "symbol": [ens2sym.get(g, g) for g in genes],
        "C_pearson": cp,
        "C_spearman": csp,
        "mean_expression": H.mean(axis=0),
        "max_expression": H.max(axis=0),
    })
    return df, h, m


# ---- expression-matched background ----------------------------------------

def expr_bins(mean_expr: np.ndarray, n_bins=20):
    """Equal-frequency bin index per gene by mean expression (rank-based)."""
    ranks = stats.rankdata(mean_expr, method="ordinal")
    return np.minimum((ranks - 1) * n_bins // len(ranks), n_bins - 1)


def matched_draws(target_idx, bins, n_draws, rng, pool_idx=None):
    """Yield n_draws index arrays matching the per-bin composition of target_idx.

    target_idx, pool_idx: positional indices into the bins array.
    Sampling is without replacement within a draw, from the full pool in each bin.
    """
    if pool_idx is None:
        pool_idx = np.arange(len(bins))
    pool_idx = np.asarray(pool_idx)
    # per-bin pools
    bin_pool = {b: pool_idx[bins[pool_idx] == b] for b in np.unique(bins[pool_idx])}
    tb = bins[target_idx]
    need = {b: int((tb == b).sum()) for b in np.unique(tb)}
    out = []
    for _ in range(n_draws):
        pick = []
        for b, k in need.items():
            pool = bin_pool[b]
            if k >= len(pool):
                pick.append(pool)
            else:
                pick.append(rng.choice(pool, size=k, replace=False))
        out.append(np.concatenate(pick))
    return out


def emp_p_greater(obs, null):
    null = np.asarray(null)
    return (np.sum(null >= obs) + 1) / (len(null) + 1)


def emp_p_two_sided(obs, null):
    null = np.asarray(null)
    c = min((np.sum(null >= obs) + 1), (np.sum(null <= obs) + 1))
    return min(1.0, 2 * c / (len(null) + 1))
