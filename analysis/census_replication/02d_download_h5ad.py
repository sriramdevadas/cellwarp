"""Download h5ad files directly from CELLxGENE, skipping files > 5GB."""
import re
import subprocess
import cellxgene_census
from pathlib import Path

DATA = Path(__file__).parent / "h5ad_cache"
DATA.mkdir(exist_ok=True)

# Mouse datasets (skip 43GB embryonic)
MOUSE_DS = [
    "7b6bab5a-f9c4-4a56-9ed4-3b9079b14867",  # Mouse all cells 194K — hepatocyte, granulocyte, T, fibroblast
    "58b01044-c5e5-4b0f-8a2d-6ebf951e01ff",  # CNS borders 123K
    "4c4cfb38-c2af-4524-8ef8-bbcf1b6e2670",  # mouse limb 23K
    "a2da8d7b-54a8-47d1-a0d3-aafcd0535f00",  # white adipose 78K
    "731e0ae7-e600-470f-a6dc-8c35c28d6c3d",  # perinatal heart 16K
    "047d57f2-4d14-45de-aa98-336c6f583750",  # lung 17K
    "6c6b4c47-096d-4084-97e7-714ee10c556c",  # Brca1 mammary 27K
    "49e4ffcc-5444-406d-bdee-577127404ba8",  # pancreas 6K
    "25818bf7-e2a7-41ec-8ff2-bc369c0ff4f5",  # mouse kidney 203K
]
HUMAN_DS = [
    "2adb1f8a-a6b1-4909-8ee8-484814e2d4bf",  # cell landscape 393K
    "b61a921b-7fa3-4b42-b455-aaaf32447920",  # airways 38K
    "fd072bc3-2dfb-46f8-b4e3-467cb3223182",  # developmental 350K
    "37b21763-7f0f-41ae-9001-60bad6e2841d",  # pancreas 21K
    "ee195b7d-184d-4dfa-9b1c-51a7e601ac11",  # duodenum 2K
    "65badd7a-9262-4fd1-9ce2-eb5dc0ca8039",  # cardiac 19K
]

ALL = MOUSE_DS + HUMAN_DS

with cellxgene_census.open_soma(census_version="2025-11-08") as c:
    ds_df = c["census_info"]["datasets"].read().concat().to_pandas()

MAX_SIZE_GB = 5

for ds_id in ALL:
    matches = ds_df[ds_df["dataset_id"] == ds_id]
    if len(matches) == 0:
        raise ValueError(
            f"dataset_id {ds_id!r} not found in CELLxGENE Census "
            f"(version pinned in this script). This usually means the "
            f"UUID literal is wrong (typo, transcription error, or stale). "
            f"Check the source list and verify against "
            f"https://cellxgene.cziscience.com/."
        )
    if len(matches) > 1:
        raise ValueError(
            f"dataset_id {ds_id!r} matched {len(matches)} Census records; "
            f"expected exactly 1. Census may have changed schema."
        )
    row = matches

    citation = row.iloc[0].get("citation", "")
    urls = re.findall(r'https://datasets\.cellxgene\.cziscience\.com/[^\s]+\.h5ad', citation)
    vid = row.iloc[0].get("dataset_version_id", "")
    url = urls[0] if urls else f"https://datasets.cellxgene.cziscience.com/{vid}.h5ad"
    title = row.iloc[0].get("dataset_title", "?")[:50]

    outpath = DATA / f"{ds_id}.h5ad"
    if outpath.exists():
        mb = outpath.stat().st_size / 1e6
        print(f"EXISTS: {ds_id[:12]} ({mb:.0f} MB) {title}")
        continue

    # Check file size first
    print(f"CHECKING: {ds_id[:12]} {title}")
    result = subprocess.run(
        ["curl", "-sI", "-L", url],
        capture_output=True, text=True
    )
    size_header = [l for l in result.stdout.split('\n') if 'Content-Length' in l]
    if size_header:
        size_bytes = int(size_header[-1].split(':')[1].strip())
        size_gb = size_bytes / 1e9
        print(f"  Size: {size_gb:.1f} GB")
        if size_gb > MAX_SIZE_GB:
            print(f"  SKIP: too large (> {MAX_SIZE_GB} GB)")
            continue
    else:
        print("  Size: unknown (proceeding)")

    print(f"  Downloading...", flush=True)
    subprocess.run(["curl", "-L", "-o", str(outpath), "-#", url], check=True)
    mb = outpath.stat().st_size / 1e6
    print(f"  Saved: {mb:.0f} MB")

print("\nDone. Files in:", DATA)
