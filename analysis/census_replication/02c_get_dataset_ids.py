"""Print full dataset IDs for the greedy set cover selection."""
import cellxgene_census
import pandas as pd
from pathlib import Path

OUT = Path(__file__).parent
TARGET = sorted(pd.read_csv(OUT / "passing_cell_types.csv")["cell_type"].tolist())
EXCLUDED_NAMES = ["Tabula Sapiens", "Tabula Muris Senis", "CellHint", "PanSci", "EasySci"]

with cellxgene_census.open_soma(census_version="2025-11-08") as census:
    ds_df = census["census_info"]["datasets"].read().concat().to_pandas()
    excluded = set()
    for n in EXCLUDED_NAMES:
        excluded.update(ds_df.loc[ds_df["collection_name"].str.contains(n, case=False, na=False), "dataset_id"])

    for org, label in [("Mus musculus", "MOUSE"), ("Homo sapiens", "HUMAN")]:
        names = ", ".join(f"'{ct}'" for ct in TARGET)
        obs = cellxgene_census.get_obs(
            census, org,
            value_filter=f"cell_type in [{names}] and is_primary_data == True and disease == 'normal'",
            column_names=["cell_type", "dataset_id"],
        )
        obs = obs[~obs["dataset_id"].isin(excluded)]

        ds_stats = []
        for ds_id, grp in obs.groupby("dataset_id", observed=True):
            ct_counts = grp["cell_type"].value_counts()
            passing = [ct for ct in ct_counts.index if ct in TARGET and ct_counts[ct] >= 200]
            ds_stats.append({"dataset_id": ds_id, "total": len(grp), "n_types": len(passing),
                             "passing": passing})
        ds_stats = pd.DataFrame(ds_stats).sort_values("n_types", ascending=False)

        covered = set()
        selected = []
        remaining = ds_stats.copy()
        while len(covered) < len(TARGET) and len(remaining) > 0:
            remaining["new"] = remaining["passing"].apply(lambda x: len(set(x) - covered))
            remaining = remaining.sort_values("new", ascending=False)
            best = remaining.iloc[0]
            if best["new"] == 0:
                break
            selected.append(best)
            covered.update(best["passing"])
            remaining = remaining.iloc[1:]

        print(f"\n{label}_DATASETS = [")
        for s in selected:
            title = ds_df.loc[ds_df["dataset_id"] == s["dataset_id"], "dataset_title"].values
            t = title[0][:50] if len(title) > 0 else "?"
            print(f'    "{s["dataset_id"]}",  # {t}, {s["total"]:,} cells, {s["n_types"]} types')
        print("]")
        missing = set(TARGET) - covered
        if missing:
            print(f"Missing: {sorted(missing)}")
