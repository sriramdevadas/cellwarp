"""
Test that data/replication/tabula_microcebus_metadata.csv stays in sync
with the Tabula Microcebus deposit anchors in the manuscript's Data
and code availability statement, and with the
analysis/mouse_lemur/feasibility_check.md artifact.

The manuscript is the canonical source-of-truth for deposit anchors
(collection_id, dataset_id, download date). The CSV must conform to it.
This test fails loudly if either the CSV or the manuscript drifts and
the other doesn't.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
MANIFEST = REPO_ROOT / "data" / "replication" / "tabula_microcebus_metadata.csv"
MANUSCRIPT = REPO_ROOT / "docs" / "submission" / "manuscript_combined.txt"
FEASIBILITY = REPO_ROOT / "analysis" / "mouse_lemur" / "feasibility_check.md"

EXPECTED_COLUMNS = [
    "dataset_id",
    "species",
    "collection_id",
    "collection_name",
    "dataset_title",
    "dataset_total_cell_count",
    "collection_doi",
    "source",
    "download_date",
    "analysis_assay_filter",
]


def test_manifest_exists() -> None:
    assert MANIFEST.exists(), f"missing: {MANIFEST}"


def test_manifest_shape_and_columns() -> None:
    df = pd.read_csv(MANIFEST)
    assert df.shape == (1, 10), f"expected (1, 10), got {df.shape}"
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


def test_manifest_anchors_present_in_manuscript() -> None:
    """
    The four deposit-anchor identifiers (collection_id, dataset_id,
    collection_doi, download_date) in the CSV must each appear in
    the manuscript text. The manuscript is the canonical source.
    """
    df = pd.read_csv(MANIFEST)
    row = df.iloc[0]
    manuscript = MANUSCRIPT.read_text(encoding='utf-8')
    anchor_fields = ["collection_id", "dataset_id", "collection_doi", "download_date"]
    missing = [
        field for field in anchor_fields
        if str(row[field]) not in manuscript
    ]
    assert not missing, (
        f"manifest values absent from manuscript "
        f"(canonical source-of-truth): {missing}\n"
        f"This means either the manifest is wrong or the manuscript "
        f"changed without updating the manifest."
    )


def test_manifest_anchors_present_in_feasibility_check() -> None:
    """
    Same anchors must also appear in feasibility_check.md, which is
    the analysis-side record of the dataset's provenance.
    """
    df = pd.read_csv(MANIFEST)
    row = df.iloc[0]
    feas = FEASIBILITY.read_text()
    anchor_fields = ["collection_id", "dataset_id", "collection_doi", "download_date"]
    missing = [
        field for field in anchor_fields
        if str(row[field]) not in feas
    ]
    assert not missing, (
        f"manifest values absent from analysis/mouse_lemur/feasibility_check.md: "
        f"{missing}\n"
        f"This usually means the analysis-side record was rebuilt "
        f"without updating the deposit manifest."
    )


def test_species_is_microcebus_murinus() -> None:
    """The dataset is the mouse lemur atlas; species must be Microcebus murinus."""
    df = pd.read_csv(MANIFEST)
    assert df.iloc[0]["species"] == "Microcebus murinus", (
        f"expected species='Microcebus murinus', got {df.iloc[0]['species']!r}"
    )


def test_source_is_cellxgene_discover() -> None:
    """
    Per the manuscript's primary-datasets paragraph, Tabula Microcebus is
    accessed via CELLxGENE Discover (not Census). source column must reflect
    this.
    """
    df = pd.read_csv(MANIFEST)
    assert df.iloc[0]["source"] == "CELLxGENE Discover", (
        f"expected source='CELLxGENE Discover', got {df.iloc[0]['source']!r}"
    )


def test_assay_filter_matches_manuscript() -> None:
    """
    Per the manuscript's mouse-lemur Methods paragraph, the analysis
    filtered to 10x 3' v2.
    The manifest's analysis_assay_filter must match what the
    manuscript says was used.
    """
    df = pd.read_csv(MANIFEST)
    assay = df.iloc[0]["analysis_assay_filter"]
    manuscript = MANUSCRIPT.read_text(encoding='utf-8')
    assert assay in manuscript, (
        f"analysis_assay_filter={assay!r} not found in manuscript text"
    )
