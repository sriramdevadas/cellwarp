"""Pull METADATA ONLY from Census to map donor structure per type/species."""
from __future__ import annotations
import sys, json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
import cellxgene_census
from cellwarp.data_loader import (HUMAN_COLLECTION, MOUSE_COLLECTION, HUMAN_ORGANISM,
                                  MOUSE_ORGANISM, get_dataset_ids_for_collection)

HERE = Path(__file__).resolve().parent
CENSUS_VERSION = "2025-11-08"
TYPES = list(pd.read_csv(ROOT / "output/phase2/scaled_35types/centroids_human_35.csv",
                         index_col=0).index)
COLS = ["cell_type", "donor_id", "dataset_id", "assay", "tissue", "tissue_general", "sex"]

def pull(census, organism, dataset_ids):
    names = ", ".join(f"'{t}'" for t in TYPES)
    ids = ", ".join(f"'{d}'" for d in dataset_ids)
    vf = (f"cell_type in [{names}] and is_primary_data == True and "
          f"disease == 'normal' and dataset_id in [{ids}]")
    obs = cellxgene_census.get_obs(census, organism, value_filter=vf, column_names=COLS)
    return obs

summary = {}
with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
    for tag, organism, coll in [("human", HUMAN_ORGANISM, HUMAN_COLLECTION),
                                ("mouse", MOUSE_ORGANISM, MOUSE_COLLECTION)]:
        print("=" * 70, f"\n{tag.upper()} — {coll}")
        ds = get_dataset_ids_for_collection(census, coll)
        obs = pull(census, organism, ds)
        print(f"  total cells: {len(obs):,}; assays: {sorted(obs.assay.unique())}; "
              f"datasets: {obs.dataset_id.nunique()}; donors: {obs.donor_id.nunique()}; "
              f"tissues: {obs.tissue.nunique()}")
        obs.to_parquet(HERE / f"meta_{tag}.parquet")
        # per type donor counts
        rows = []
        for t in TYPES:
            sub = obs[obs.cell_type == t]
            dc = sub.groupby("donor_id").size()
            dc = dc[dc > 0]
            rows.append(dict(cell_type=t, n_cells=len(sub), n_donors=int((dc > 0).sum()),
                             n_donors_ge50=int((dc >= 50).sum()),
                             n_donors_ge100=int((dc >= 100).sum()),
                             median_cells_per_donor=float(dc.median()) if len(dc) else 0,
                             n_datasets=int(sub.dataset_id.nunique()),
                             n_assays=int(sub.assay.nunique())))
        tab = pd.DataFrame(rows).sort_values("n_donors", ascending=False)
        tab.to_csv(HERE / f"donor_structure_{tag}.csv", index=False)
        print(tab.to_string(index=False))
        summary[tag] = dict(total_cells=int(len(obs)), n_donors=int(obs.donor_id.nunique()),
                            assays=sorted(obs.assay.unique().tolist()),
                            n_datasets=int(obs.dataset_id.nunique()),
                            donor_powered_types_ge4donors=int((tab.n_donors >= 4).sum()),
                            donor_powered_types_ge6donors=int((tab.n_donors >= 6).sum()))
json.dump(summary, open(HERE / "phase0_summary.json", "w"), indent=2)
print("\nSUMMARY:", json.dumps(summary, indent=2))
