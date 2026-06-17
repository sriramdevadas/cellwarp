"""
Pull cell-level TS/TMS from Census (35 types, 16,959 orthologs, with donor+assay),
normalize (CP10k+log1p) per cell, and aggregate to:
  - per-(type,donor,protocol) summed log-expression + cell count   [donor-split, caps, x-protocol]
  - within-donor cell-split half-centroids, K reps                 [cell-bootstrap ceiling]
Memory-efficient: one type at a time; checkpoints per (species,type).
"""
from __future__ import annotations
import sys, json, pickle, time
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import cellxgene_census

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from cellwarp.data_loader import (HUMAN_ORGANISM, MOUSE_ORGANISM, HUMAN_COLLECTION,
                                  MOUSE_COLLECTION, get_dataset_ids_for_collection)

HERE = Path(__file__).resolve().parent
CKPT = HERE / "checkpoints"; CKPT.mkdir(exist_ok=True)
CENSUS_VERSION = "2025-11-08"
CAP = 10000
CAPS = [500, 2000, 10000]
SEED = 42
K_CS = 20

gene_order = list(pd.read_csv(ROOT / "output/phase2/scaled_35types/centroids_human_35.csv",
                              index_col=0).columns)  # 16,959 human ENSG, deposit order
TYPES = list(pd.read_csv(ROOT / "output/phase2/scaled_35types/centroids_human_35.csv",
                         index_col=0).index)
NG = len(gene_order)
orth = pd.read_csv(ROOT / "data/phase1/orthologs_human_mouse.csv")
orth = orth[orth.orthology_type == "ortholog_one2one"]
h2m = dict(zip(orth.human_ensembl_id, orth.mouse_ensembl_id))
mouse_fids = [h2m[g] for g in gene_order]


def proto_of(assay: str) -> str:
    return "10x" if "10x" in assay else "SS"


def process_species(census, tag, organism, collection, species_target_fids):
    ds = get_dataset_ids_for_collection(census, collection)
    var = cellxgene_census.get_var(census, organism, column_names=["feature_id", "soma_joinid"])
    fid2coord = dict(zip(var.feature_id, var.soma_joinid))
    present = [f for f in species_target_fids if f in fid2coord]
    assert len(present) == NG, f"{tag}: only {len(present)}/{NG} ortholog genes in var"
    var_coords = [int(fid2coord[f]) for f in species_target_fids]
    idsf = ", ".join(f"'{d}'" for d in ds)
    base = f"is_primary_data == True and disease == 'normal' and dataset_id in [{idsf}]"

    for ti, t in enumerate(TYPES):
        cpath = CKPT / f"{tag}_{ti:02d}.pkl"
        if cpath.exists():
            print(f"  [skip] {tag} {t}", flush=True)
            continue
        t_start = time.time()
        vf = f"cell_type == '{t.replace(chr(39), chr(39)*2)}' and " + base
        # sequential server-side scan (fast) over the whole type, genes pre-restricted
        ad = cellxgene_census.get_anndata(census, organism, obs_value_filter=vf,
                                          var_coords=var_coords,
                                          obs_column_names=["donor_id", "assay"],
                                          var_column_names=["feature_id"])
        n = ad.n_obs
        rng = np.random.default_rng(SEED)
        if n > CAP:
            keep = np.sort(rng.choice(n, CAP, replace=False))
            ad = ad[keep].copy()
        # normalize per cell
        sc.pp.normalize_total(ad, target_sum=1e4)
        sc.pp.log1p(ad)
        # align genes to gene_order
        fid2idx = {f: i for i, f in enumerate(ad.var.feature_id.values)}
        perm = np.array([fid2idx[f] for f in species_target_fids], dtype=int)
        X = ad.X.tocsc()[:, perm].tocsr().astype(np.float32)  # cells x NG (gene_order)
        donors = ad.obs.donor_id.astype(str).values
        protos = np.array([proto_of(a) for a in ad.obs.assay.astype(str).values])
        npull = X.shape[0]

        # per-(donor,proto) sums + counts at each cell-count cap (power analysis)
        groups_by_cap = {}
        for cap in CAPS:
            if npull > cap:
                rcap = np.random.default_rng(SEED + cap)
                rows_cap = np.sort(rcap.choice(npull, cap, replace=False))
            else:
                rows_cap = np.arange(npull)
            dcap = donors[rows_cap]; pcap = protos[rows_cap]; Xc = X[rows_cap]
            g = {}
            for d in np.unique(dcap):
                for p in ("10x", "SS"):
                    m = (dcap == d) & (pcap == p)
                    c = int(m.sum())
                    if c:
                        g[(d, p)] = (np.asarray(Xc[m].sum(0)).ravel().astype(np.float32), c)
            groups_by_cap[cap] = g
        groups = groups_by_cap[CAP]  # back-compat alias (full cap)
        # within-donor cell-split half-centroids
        csA = np.zeros((K_CS, NG), np.float32); csB = np.zeros((K_CS, NG), np.float32)
        cntA = np.zeros(K_CS, int); cntB = np.zeros(K_CS, int)
        for r in range(K_CS):
            rr = np.random.default_rng(1000 + r)
            ai, bi = [], []
            for d in np.unique(donors):
                rows = np.where(donors == d)[0]
                rr.shuffle(rows)
                h = len(rows) // 2
                ai.append(rows[:h]); bi.append(rows[h:])
            ai = np.concatenate(ai) if ai else np.array([], int)
            bi = np.concatenate(bi) if bi else np.array([], int)
            if len(ai): csA[r] = np.asarray(X[ai].sum(0)).ravel(); cntA[r] = len(ai)
            if len(bi): csB[r] = np.asarray(X[bi].sum(0)).ravel(); cntB[r] = len(bi)
        pickle.dump(dict(type=t, n_total=n, n_pulled=X.shape[0],
                         groups_by_cap=groups_by_cap, csA=csA, csB=csB, cntA=cntA, cntB=cntB),
                    open(cpath, "wb"))
        print(f"  [done] {tag} {t}: pulled {X.shape[0]}/{n} cells, "
              f"{len(np.unique(donors))} donors, protos={sorted(set(protos))}, "
              f"{time.time()-t_start:.0f}s", flush=True)


def combine(tag):
    parts = [pickle.load(open(CKPT / f"{tag}_{ti:02d}.pkl", "rb")) for ti in range(len(TYPES))]
    csA = np.zeros((K_CS, len(TYPES), NG), np.float32); csB = np.zeros_like(csA)
    cntA = np.zeros((K_CS, len(TYPES)), int); cntB = np.zeros((K_CS, len(TYPES)), int)
    for ti, part in enumerate(parts):
        csA[:, ti] = part["csA"]; csB[:, ti] = part["csB"]
        cntA[:, ti] = part["cntA"]; cntB[:, ti] = part["cntB"]
    np.savez_compressed(HERE / f"agg_{tag}_cs.npz", csA=csA, csB=csB, cntA=cntA, cntB=cntB,
                        genes=np.array(gene_order), types=np.array(TYPES))
    # per-cap group tables
    for cap in CAPS:
        gmeta, gsums = [], []
        for ti, part in enumerate(parts):
            for (d, p), (s, c) in part["groups_by_cap"][cap].items():
                gmeta.append((ti, part["type"], d, p, c)); gsums.append(s)
        gmeta = pd.DataFrame(gmeta, columns=["type_idx", "type", "donor", "proto", "count"])
        gsums = np.array(gsums, np.float32)
        np.savez_compressed(HERE / f"agg_{tag}_cap{cap}.npz", gsums=gsums,
                            genes=np.array(gene_order), types=np.array(TYPES))
        gmeta.to_csv(HERE / f"agg_{tag}_cap{cap}_groups.csv", index=False)
        print(f"  combined {tag} cap{cap}: {len(gmeta)} groups, sums {gsums.shape}")


if __name__ == "__main__":
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        process_species(census, "human", HUMAN_ORGANISM, HUMAN_COLLECTION, gene_order)
        process_species(census, "mouse", MOUSE_ORGANISM, MOUSE_COLLECTION, mouse_fids)
    combine("human"); combine("mouse")
    print("PULL+AGGREGATE COMPLETE")
