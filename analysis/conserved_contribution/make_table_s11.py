"""
Build Supplementary Table S11 (per-gene cross-species conservation score C).

Derives the deposited supplementary table from the verified per-gene table
(gene_conservation_core.csv), adding Tau specificity (computed from the human
35-type centroids, as in the gate) and the conserved/divergent-quartile
assignment (top/bottom quartile of valid C, frozen thresholds from
gate_results.json). Output: docs/supplementary_materials/table_S11_gene_conservation.csv.

Additive columns (existing columns/values unchanged):
  master_tf_flag           membership in the 73-gene pre-registered master-TF set
  procrustes_contribution  the within-data axis-loading metric (the rho=0.27
                           non-circularity contrast vs C)
  donor_split_C_A / _B     per-gene C in the two halves of a representative donor
                           split (split 0; BOTH atlases' donors partitioned), the
                           per-gene basis of the Figure 2D cross-half agreement.
                           Source aggregates are the large, gitignored donor-
                           resampling inputs in donor_stability/ (same inputs as
                           Figure 2D); left NaN if those aggregates are absent.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import gate_lib as G

core = pd.read_csv(HERE / "gene_conservation_core.csv")
gate = json.load(open(HERE / "gate_results.json"))
q75 = gate["thresholds"]["q75"]
q25 = gate["thresholds"]["q25"]

# Tau specificity from the human 35-type centroids (frozen, as in the gate).
h, _ = G.load_centroids()
Hc = h[core["gene_id"]].values  # 35 x n, in the deposit gene order
tau = np.array([np.sum(1 - Hc[:, j] / Hc[:, j].max()) / (Hc.shape[0] - 1)
                if Hc[:, j].max() > 0 else np.nan for j in range(len(core))])

C = core["C_pearson"].values
quartile = np.where(np.isnan(C), "",
            np.where(C >= q75, "conserved",
             np.where(C <= q25, "divergent", "intermediate")))

# ── additive columns (Stage 1 / cb m4) ──────────────────────────────────────
# master-TF membership: symbol in the 73-gene pre-registered positive-control set
master_tf_flag = core["symbol"].isin(set(G.POSITIVE_CONTROL_TFS)).astype(int)

# within-data axis-loading metric (the non-circularity contrast, rho=0.27 vs C)
_gc = pd.read_csv(G.GENE_CONS_TABLE)[["gene_id", "procrustes_contribution"]]
procrustes_contribution = core.merge(_gc, on="gene_id", how="left")["procrustes_contribution"]

# per-gene C in the two halves of a representative donor split (split 0), the
# per-gene basis of the Figure 2D cross-half reproducibility. BOTH atlases'
# donor sets are partitioned (human + mouse), matching donor_stability.
DONOR = HERE / "donor_stability"
def _donor_split_CA_CB(split=0, cap=10000):
    dh = np.load(DONOR / f"agg_human_cap{cap}.npz", allow_pickle=True)
    dm = np.load(DONOR / f"agg_mouse_cap{cap}.npz", allow_pickle=True)
    genes = list(dh["genes"]); NG = len(genes)
    def load(d, tag):
        g = pd.read_csv(DONOR / f"agg_{tag}_cap{cap}_groups.csv")
        return d["gsums"], g["type_idx"].values, g["donor"].astype(str).values, g["count"].values.astype(float)
    hg, ht, hdn, hc = load(dh, "human")
    mg, mt, mdn, mc = load(dm, "mouse")
    def centroid(gsums, tidx, donor, count, donor_set):
        sel = np.isin(donor, np.array(list(donor_set)))
        H = np.full((35, NG), np.nan)
        for t in range(35):
            m = sel & (tidx == t); c = count[m].sum()
            if c > 0:
                H[t] = gsums[m].sum(0) / c
        return H
    def cvec(H, M):
        ok = ~(np.isnan(H).any(1) | np.isnan(M).any(1))
        H, M = H[ok], M[ok]
        Hc = H - H.mean(0); Mc = M - M.mean(0)
        den = np.sqrt((Hc ** 2).sum(0) * (Mc ** 2).sum(0))
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(den == 0, np.nan, (Hc * Mc).sum(0) / den)
    rs = np.random.default_rng(100 + split)
    hp = rs.permutation(np.unique(hdn)); mp = rs.permutation(np.unique(mdn))
    hA, hB = set(hp[:len(hp) // 2]), set(hp[len(hp) // 2:])
    mA, mB = set(mp[:len(mp) // 2]), set(mp[len(mp) // 2:])
    CA = cvec(centroid(hg, ht, hdn, hc, hA), centroid(mg, mt, mdn, mc, mA))
    CB = cvec(centroid(hg, ht, hdn, hc, hB), centroid(mg, mt, mdn, mc, mB))
    return dict(zip(genes, CA)), dict(zip(genes, CB))

try:
    _cA, _cB = _donor_split_CA_CB(split=0)
except FileNotFoundError as _exc:
    # D67: this used to fall back to an all-NaN pair and carry on to exit 0. The
    # result was an S11 Table that looked built -- 16,959 rows, every other column
    # populated -- with both donor-split columns empty and nothing but a stdout
    # warning to say so. On a fresh clone the aggregates are always absent
    # (donor_stability/agg_*.npz and the agg_*_groups.csv beside them are
    # untracked), so that path was the one a reader would actually take.
    #
    # No fallback: the deposited table has these columns populated (15,959 and
    # 15,928 non-null of 16,959), so an absent aggregate means this run cannot
    # reproduce the deposited artifact, and it should say so rather than emit a
    # gutted one.
    raise SystemExit(
        "ERROR: cannot compute donor_split_C_A / donor_split_C_B.\n"
        "  missing: %s\n"
        "  wanted under: %s\n"
        "    agg_human_cap10000.npz, agg_mouse_cap10000.npz,\n"
        "    agg_human_cap10000_groups.csv, agg_mouse_cap10000_groups.csv\n"
        "  These aggregates are untracked, so a fresh clone does not carry them.\n"
        "  Regenerate them with donor_stability/run_donor_stability.py before\n"
        "  building S11 Table. Refusing to write the table with both donor-split\n"
        "  columns empty."
        % (getattr(_exc, "filename", None) or _exc, DONOR)
    ) from _exc
donor_split_C_A = core["gene_id"].map(_cA)
donor_split_C_B = core["gene_id"].map(_cB)

out = pd.DataFrame({
    "gene_id": core["gene_id"],
    "symbol": core["symbol"],
    "C_pearson": core["C_pearson"],
    "C_spearman": core["C_spearman"],
    "mean_expression": core["mean_expression"],
    "tau": tau,
    "quartile": quartile,
    "master_tf_flag": master_tf_flag,
    "procrustes_contribution": procrustes_contribution,
    "donor_split_C_A": donor_split_C_A,
    "donor_split_C_B": donor_split_C_B,
})
dst = ROOT / "docs/supplementary_materials/table_S11_gene_conservation.csv"
out.to_csv(dst, index=False)
n_valid = int(core["C_pearson"].notna().sum())
print(f"wrote {dst.relative_to(ROOT)}: {len(out)} rows ({n_valid} with defined C); "
      f"conserved={int((quartile=='conserved').sum())} divergent={int((quartile=='divergent').sum())}")
