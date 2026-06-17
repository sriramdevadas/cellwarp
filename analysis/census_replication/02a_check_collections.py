"""Quick check: which primary cell types are in each cross-species collection?"""
import cellxgene_census
import pandas as pd
from pathlib import Path

OUT_DIR = Path(__file__).parent

PRIMARY_TYPES = [
    "stromal cell", "epithelial cell", "hematopoietic precursor cell",
    "hematopoietic stem cell", "pancreatic acinar cell", "basal cell",
    "T cell", "neutrophil", "fibroblast of cardiac tissue",
    "myeloid leukocyte", "mesenchymal stem cell of adipose tissue",
    "plasma cell", "mesenchymal stem cell",
    "CD4-positive, alpha-beta T cell", "classical monocyte",
    "macrophage", "B cell",
    "luminal epithelial cell of mammary gland",
    "large intestine goblet cell",
    "enterocyte of epithelium of large intestine",
    "myeloid dendritic cell", "monocyte", "natural killer cell",
    "intermediate monocyte", "mature NK T cell", "adventitial cell",
    "granulocyte", "fibroblast", "bladder urothelial cell",
    "pancreatic ductal cell", "smooth muscle cell", "hepatocyte",
    "endothelial cell", "non-classical monocyte",
    "CD8-positive, alpha-beta T cell",
]

TARGET_COLLECTIONS = [
    "A single-cell atlas of human and mouse white adipose tissue",
    "Single-cell roadmap of human gonadal development",
    "Human skeletal muscle ageing atlas",
    "Cellular development and evolution of the mammalian cerebellum",
    "Human and mouse dermal fibroblast atlas",
    "Single cell analysis of mouse and human prostate reveals novel cell types and conserved programs of prostate cancer and inflammation",
]

with cellxgene_census.open_soma(census_version="2025-11-08") as census:
    datasets_df = census["census_info"]["datasets"].read().concat().to_pandas()

    for coll_name in TARGET_COLLECTIONS:
        mask = datasets_df["collection_name"].str.contains(coll_name[:40], case=False, na=False)
        ds_ids = datasets_df.loc[mask, "dataset_id"].tolist()
        actual_name = datasets_df.loc[mask, "collection_name"].iloc[0] if mask.any() else coll_name

        print(f"\n{'='*70}")
        print(f"COLLECTION: {actual_name[:70]}")
        print(f"  Datasets: {len(ds_ids)}")

        if not ds_ids:
            print("  NO DATASETS FOUND")
            continue

        ds_filter = ", ".join(f"'{d}'" for d in ds_ids)

        for organism, label in [("Homo sapiens", "Human"), ("Mus musculus", "Mouse")]:
            obs = cellxgene_census.get_obs(
                census, organism,
                value_filter=f"dataset_id in [{ds_filter}] and is_primary_data == True and disease == 'normal'",
                column_names=["cell_type"],
            )
            if len(obs) == 0:
                print(f"\n  {label}: 0 cells")
                continue

            ct_counts = obs["cell_type"].value_counts()
            primary_matches = [ct for ct in ct_counts.index if ct in PRIMARY_TYPES]

            print(f"\n  {label}: {len(obs):,} cells, {len(ct_counts)} types, "
                  f"{len(primary_matches)} overlap with primary 35")

            # Show overlapping types with counts
            for ct in sorted(primary_matches):
                n = ct_counts[ct]
                gate = "ok" if n >= 500 else "LOW"
                print(f"    {ct:<50} {n:>6,} {gate}")

print("\nDone.")
