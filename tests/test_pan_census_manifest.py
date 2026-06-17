"""
Test that data/replication/pan_census_manifest.csv stays in sync with
the canonical UUID lists in
analysis/census_replication/02_run_replication.py.

Drift between the manifest and the script would mean the deposit's
documented data sources don't match the analysis's actual data sources.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
MANIFEST = REPO_ROOT / "data" / "replication" / "pan_census_manifest.csv"
SCRIPT = REPO_ROOT / "analysis" / "census_replication" / "02_run_replication.py"

UUID_RE = re.compile(
    r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}"
)

EXPECTED_COLUMNS = [
    "dataset_id",
    "species",
    "collection_id",
    "collection_name",
    "dataset_title",
    "dataset_version_id",
    "dataset_total_cell_count",
    "collection_doi",
    "citation",
    "census_version",
]


def _extract_script_uuids() -> tuple[list[str], list[str]]:
    """Return (mouse_uuids, human_uuids) from the script's MOUSE_FILES/HUMAN_FILES."""
    src = SCRIPT.read_text()
    mouse_block = re.search(
        r"MOUSE_FILES\s*=\s*\[(.*?)\]", src, re.DOTALL
    ).group(1)
    human_block = re.search(
        r"HUMAN_FILES\s*=\s*\[(.*?)\]", src, re.DOTALL
    ).group(1)
    return UUID_RE.findall(mouse_block), UUID_RE.findall(human_block)


def test_manifest_exists() -> None:
    assert MANIFEST.exists(), f"missing: {MANIFEST}"


def test_manifest_shape_and_columns() -> None:
    df = pd.read_csv(MANIFEST)
    assert df.shape == (15, 10), (
        f"expected (15, 10), got {df.shape}"
    )
    assert list(df.columns) == EXPECTED_COLUMNS, (
        f"column mismatch: expected {EXPECTED_COLUMNS}, got {list(df.columns)}"
    )


def test_manifest_no_missing_values() -> None:
    df = pd.read_csv(MANIFEST)
    na_per_col = df.isna().sum()
    cols_with_na = na_per_col[na_per_col > 0]
    assert len(cols_with_na) == 0, (
        f"manifest has missing values:\n{cols_with_na}"
    )


def test_manifest_uuids_match_script() -> None:
    """Every UUID in the script must appear in the manifest, and vice versa."""
    mouse_uuids, human_uuids = _extract_script_uuids()
    assert len(mouse_uuids) == 9, (
        f"expected 9 mouse UUIDs in script, got {len(mouse_uuids)}"
    )
    assert len(human_uuids) == 6, (
        f"expected 6 human UUIDs in script, got {len(human_uuids)}"
    )
    script_uuids = set(mouse_uuids) | set(human_uuids)

    df = pd.read_csv(MANIFEST)
    manifest_uuids = set(df["dataset_id"])

    only_in_script = sorted(script_uuids - manifest_uuids)
    only_in_manifest = sorted(manifest_uuids - script_uuids)
    assert not only_in_script, (
        f"UUIDs in script but missing from manifest: {only_in_script}"
    )
    assert not only_in_manifest, (
        f"UUIDs in manifest but absent from script: {only_in_manifest}"
    )


def test_species_assignments_match_script_blocks() -> None:
    """A UUID's species column must match which list block it came from."""
    mouse_uuids, human_uuids = _extract_script_uuids()
    df = pd.read_csv(MANIFEST).set_index("dataset_id")

    for u in mouse_uuids:
        assert df.loc[u, "species"] == "Mus musculus", (
            f"{u} is in MOUSE_FILES but manifest says "
            f"species={df.loc[u, 'species']!r}"
        )
    for u in human_uuids:
        assert df.loc[u, "species"] == "Homo sapiens", (
            f"{u} is in HUMAN_FILES but manifest says "
            f"species={df.loc[u, 'species']!r}"
        )


def test_census_version_uniform() -> None:
    """All rows must declare the same Census version (the analysis pin)."""
    df = pd.read_csv(MANIFEST)
    versions = df["census_version"].unique()
    assert len(versions) == 1, (
        f"expected 1 unique census_version, got {len(versions)}: {versions}"
    )
    assert versions[0] == "2025-11-08", (
        f"expected census_version='2025-11-08', got {versions[0]!r}"
    )


def test_collection_count() -> None:
    """
    The 15 pan-Census UUIDs occupy 15 distinct collection_ids in
    Census v2025-11-08 (the deposit-pinned version). Each dataset is
    in its own collection at this Census version.

    Note: an earlier Census version (v2025-01-30) returned 13
    collections because two of the 15 datasets were not yet indexed
    by Census at that snapshot. The manuscript's original "13" claim
    reflected that earlier state and is corrected to "15" in the
    Phase 4 manuscript-edit batch. See Phase 3c-finding-2 for the
    full diagnosis.
    """
    df = pd.read_csv(MANIFEST)
    n_collections = df["collection_id"].nunique()
    assert n_collections == 15, (
        f"expected 15 distinct collections (one per dataset at "
        f"Census v2025-11-08), got {n_collections}"
    )
