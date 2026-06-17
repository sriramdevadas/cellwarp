#!/usr/bin/env python3
"""
CellWarp — Scaled Cancer Procrustes Pipeline (Thread 1 expansion)

Expands the cancer analysis from 8 coarse cell types to a fine-grained set that
maximizes overlap with the 35-type cross-species analysis. Goal: get n≥15 matched
cell types for a definitive Spearman correlation between cancer deformation and
cross-species rigidity (currently ρ=0.679, p=0.094, n=7 — underpowered).

Biology
-------
The original cancer analysis pooled all T cells, all monocyte-lineage cells, etc.
into coarse categories. By resolving CD4+ T vs CD8+ T, macrophage vs monocyte,
goblet vs enterocyte, etc., we can match more cell types to the 35-type cross-species
set and test the evolutionary-rigidity-predicts-cancer-resistance hypothesis with
adequate statistical power.

Math
----
Identical to scripts/11_cancer_procrustes.py — Procrustes alignment of condition-level
centroids (per-donor averaged) in PCA space, permutation test, residual deformation
scores. The only difference is the cell type mapping: fine-grained categories that
map 1:1 to the 35-type cross-species set wherever possible.

Steps
-----
1. Load 35-type cross-species residuals — get the target cell type list
2. Query CELLxGENE Census for colon tissue inventory (normal + adenocarcinoma)
3. Build fine-grained coarse mapping maximizing overlap with 35-type set — PAUSE
4. Download ≤3,000 cells per type per condition, ≥500 gate, ortholog gene space
5. Run cancer Procrustes pipeline (centroids → PCA → Procrustes → permutation)
6. Cross-analysis Spearman with matched 35-type residuals + scatter plot
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import signal
from contextlib import contextmanager

from cellwarp.procrustes import (
    RANDOM_SEED,
    N_PERMUTATIONS,
    PCA_VARIANCE_THRESHOLD,
    N_TOP_GENES,
    procrustes_align,
    permutation_test,
    compute_residual_vectors,
    map_residuals_to_genes,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DATA_DIR = Path("data/cancer")
OUTPUT_DIR = Path("output/cancer/scaled")
CROSS_SPECIES_RESULTS = Path("output/phase2/scaled_35types/procrustes_results_35.json")
ORTHOLOG_CSV = Path("data/phase1/orthologs_human_mouse.csv")
ORTHOLOG_H5AD = Path("data/phase1/human_aligned.h5ad")

# CELLxGENE disease labels for colon tissue
DISEASE_NORMAL = "normal"
DISEASE_TUMOR = "colon adenocarcinoma"
TISSUE_GENERAL = ["large intestine"]

# Minimum cells per fine type per condition for inclusion
MIN_CELLS_FINE = 200   # for inventory (what's available)
MIN_CELLS_GATE = 500   # for Procrustes inclusion
MAX_CELLS_PER_TYPE = 3_000

EXCLUDE_TYPES = {"Other"}

# Census query timeout (seconds) — fail fast if Census hangs
# 600s for expression data downloads (3k cells × 60k genes is large);
# 300s was too aggressive, causing cascading failures on slow nights
CENSUS_QUERY_TIMEOUT = 1200


@contextmanager
def census_timeout(seconds: int = CENSUS_QUERY_TIMEOUT):
    """Raise TimeoutError if a Census API call exceeds *seconds*."""
    def _handler(signum, frame):
        raise TimeoutError(f"Census query timed out after {seconds}s")
    prev = signal.signal(signal.SIGALRM, _handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, prev)


# ---------------------------------------------------------------------------
# Fine-grained coarse mapping rules
# ---------------------------------------------------------------------------
# Unlike the original 8-category COARSE_RULES in src/cancer_loader.py, these
# resolve subtypes to match the 35-type cross-species set. Order matters:
# first match wins. More specific patterns come before general ones.

FINE_RULES: list[tuple[str, list[str]]] = [
    # Mast cell before T cell (substring safety)
    ("mast cell", ["mast cell"]),

    # Specific epithelial subtypes BEFORE general epithelial
    ("large intestine goblet cell", [
        "goblet",
    ]),
    ("enterocyte of epithelium of large intestine", [
        "enterocyte",
    ]),
    ("epithelial cell", [
        "epithelial", "colonocyte", "paneth", "tuft", "enteroendocrine",
        "transit amplifying", "stem cell of colon", "intestinal crypt stem cell",
        "colon epithelial",
    ]),

    # Myofibroblast before T cell — "myofibroblast cell" contains substring
    # "t cell", so must match fibroblast first (ISSUE-054 fix)
    ("fibroblast", ["myofibroblast"]),

    # T cell subtypes BEFORE general T cell
    ("CD4-positive, alpha-beta T cell", [
        "cd4",
    ]),
    ("CD8-positive, alpha-beta T cell", [
        "cd8",
    ]),
    # regulatory T cells are CD4+ subset
    ("CD4-positive, alpha-beta T cell", [
        "regulatory t", "treg",
    ]),
    # NKT cells are T cells (express TCR) — ISSUE-055 / DECISION-070
    ("T cell", [
        "nk t cell", "nkt cell",
    ]),
    ("T cell", [
        "t cell", "gamma-delta", "thymocyte",
    ]),

    # B cell subtypes
    ("plasma cell", [
        "plasma cell", "plasmablast", "plasma blast",
    ]),
    ("B cell", [
        "b cell",
    ]),

    # Myeloid subtypes — macrophage vs monocyte vs dendritic separately
    ("macrophage", [
        "macrophage",
    ]),
    ("monocyte", [
        "monocyte",
    ]),
    ("myeloid dendritic cell", [
        "dendritic",
    ]),

    # NK cell
    ("natural killer cell", [
        "natural killer", "nk cell",
    ]),

    # Stromal / structural
    ("fibroblast", [
        "fibroblast", "caf", "cancer associated fibroblast",
    ]),
    ("stromal cell", [
        "stromal",
    ]),
    ("endothelial cell", [
        "endothelial",
    ]),
    ("smooth muscle cell", [
        "smooth muscle",
    ]),

    # Granulocyte / neutrophil
    ("neutrophil", [
        "neutrophil",
    ]),
    ("granulocyte", [
        "granulocyte",
    ]),
]


def _classify_fine(raw_label: str) -> str:
    """Classify a raw cell type label using fine-grained rules."""
    label_lower = raw_label.lower()
    for fine_name, keywords in FINE_RULES:
        for kw in keywords:
            if kw in label_lower:
                return fine_name
    return "Other"


def audit_substring_collisions(raw_labels: list[str]) -> bool:
    """
    Check all raw Census labels against ALL FINE_RULES keywords simultaneously.
    Flag any label that matches keywords from multiple *different* fine categories.

    Biology: CELLxGENE cell type labels often contain overlapping substrings
    (e.g., "myofibroblast cell" contains "t cell"). If FINE_RULES ordering is
    wrong, a label can be silently misclassified. This audit catches such cases.

    Math: For each raw label, iterate ALL (fine_name, keyword) pairs and collect
    every match. A collision exists when a label matches keywords belonging to
    two or more distinct fine categories (regardless of rule order).

    Args:
        raw_labels: All unique raw cell_type strings from Census.

    Returns:
        True if any collisions found, False otherwise.
    """
    print("\n" + "=" * 70)
    print("COLLISION AUDIT — Checking FINE_RULES for substring overlaps")
    print("=" * 70)

    collisions = []
    for raw in sorted(set(raw_labels)):
        raw_lower = raw.lower()
        matched_rules: list[tuple[str, str]] = []
        seen_rules: set[int] = set()  # track rule index to avoid duplicates

        for rule_idx, (fine_name, keywords) in enumerate(FINE_RULES):
            if rule_idx in seen_rules:
                continue
            for kw in keywords:
                if kw in raw_lower:
                    matched_rules.append((fine_name, kw))
                    seen_rules.add(rule_idx)
                    break  # one keyword per rule is enough

        if len(matched_rules) > 1:
            # Check if matches map to DIFFERENT fine categories
            fine_labels = set(r[0] for r in matched_rules)
            if len(fine_labels) > 1:
                collisions.append((raw, matched_rules))

    if collisions:
        print(f"\n  *** {len(collisions)} COLLISION(S) FOUND ***")
        for raw, rules in collisions:
            first_match = _classify_fine(raw)
            print(f"\n  Label: '{raw}'")
            print(f"    First-match result: {first_match}")
            print(f"    All matching rules:")
            for fine_name, kw in rules:
                is_winner = " <-- WINNER (first match)" if fine_name == first_match else ""
                print(f"      '{kw}' -> {fine_name}{is_winner}")
        print(f"\n  Review above — rule ordering prevents bugs, but verify intent.")
        return True
    else:
        print(f"\n  No substring collisions found across {len(set(raw_labels))} raw labels.")
        print(f"  All labels match at most one FINE_RULES category.")
        return False


# ---------------------------------------------------------------------------
# Cross-species type mapping for Spearman correlation
# ---------------------------------------------------------------------------
# Maps each fine cancer category to its cross-species counterpart in the 35-type
# set. For types that split differently, we can average cross-species subtypes.

CANCER_TO_XS_MAP = {
    "B cell": ["B cell"],
    "CD4-positive, alpha-beta T cell": ["CD4-positive, alpha-beta T cell"],
    "CD8-positive, alpha-beta T cell": ["CD8-positive, alpha-beta T cell"],
    "T cell": ["T cell"],
    "endothelial cell": ["endothelial cell"],
    "enterocyte of epithelium of large intestine": ["enterocyte of epithelium of large intestine"],
    "epithelial cell": ["epithelial cell"],
    "fibroblast": ["fibroblast"],
    "granulocyte": ["granulocyte"],
    "large intestine goblet cell": ["large intestine goblet cell"],
    "macrophage": ["macrophage"],
    "mast cell": None,  # no match in 35-type set
    "monocyte": ["monocyte"],
    "myeloid dendritic cell": ["myeloid dendritic cell"],
    "natural killer cell": ["natural killer cell"],
    "neutrophil": ["neutrophil"],
    "plasma cell": ["plasma cell"],
    "smooth muscle cell": ["smooth muscle cell"],
    "stromal cell": ["stromal cell"],
}


# ===================================================================
# STEP 1 — Load 35-type cross-species residuals
# ===================================================================


def load_cross_species_residuals() -> dict[str, float]:
    """Load the 35-type cross-species residual magnitudes."""
    print("=" * 70)
    print("STEP 1 — Load 35-Type Cross-Species Residuals")
    print("=" * 70)

    with open(CROSS_SPECIES_RESULTS) as f:
        xs_data = json.load(f)

    xs_residuals = {
        ct: xs_data["residuals"][ct]["magnitude"]
        for ct in xs_data["cell_types"]
    }

    print(f"\n  35 cross-species cell types and residual magnitudes:")
    for i, (ct, mag) in enumerate(
        sorted(xs_residuals.items(), key=lambda x: x[1], reverse=True), 1
    ):
        print(f"    {i:>2}. {ct:<50} {mag:.4f}")

    return xs_residuals


# ===================================================================
# STEP 2 — Query CELLxGENE Census inventory
# ===================================================================


def query_inventory() -> pd.DataFrame:
    """
    Query CELLxGENE Census for all cell_type labels in colon tissue
    (normal + colon adenocarcinoma) with ≥200 cells in BOTH conditions.

    Returns DataFrame with columns: cell_type, normal_count, tumor_count.
    """
    import cellxgene_census

    print("\n" + "=" * 70)
    print("STEP 2 — CELLxGENE Census Inventory (colon tissue)")
    print("=" * 70)

    print(f"\n  Opening Census (2025-11-08)...")
    census = cellxgene_census.open_soma(census_version="2025-11-08")

    tissue_str = ", ".join(f"'{v}'" for v in TISSUE_GENERAL)
    results = {}

    for condition, disease_label in [
        ("normal", DISEASE_NORMAL),
        ("tumor", DISEASE_TUMOR),
    ]:
        value_filter = (
            f"is_primary_data == True and "
            f"tissue_general in [{tissue_str}] and "
            f"disease == '{disease_label}'"
        )
        print(f"  Querying {condition} (disease='{disease_label}')...")
        with census_timeout():
            obs_df = cellxgene_census.get_obs(
                census,
                "Homo sapiens",
                value_filter=value_filter,
                column_names=["cell_type"],
            )
        counts = obs_df["cell_type"].value_counts().to_dict()
        results[condition] = counts
        print(f"    Found {len(obs_df):,} cells, {len(counts)} distinct cell types")

    census.close()

    # Merge into DataFrame
    all_types = sorted(set(list(results["normal"].keys()) + list(results["tumor"].keys())))
    rows = []
    for ct in all_types:
        rows.append({
            "cell_type": ct,
            "normal_count": results["normal"].get(ct, 0),
            "tumor_count": results["tumor"].get(ct, 0),
        })

    df = pd.DataFrame(rows)

    # Print full inventory
    print(f"\n  FULL INVENTORY — all cell types in colon tissue:")
    print(f"  {'cell_type':<55} {'normal':>8} {'tumor':>8}")
    print(f"  {'-' * 75}")
    for _, row in df.sort_values("normal_count", ascending=False).iterrows():
        flag = ""
        if row["normal_count"] >= MIN_CELLS_FINE and row["tumor_count"] >= MIN_CELLS_FINE:
            flag = " *"
        print(
            f"  {str(row['cell_type']):<55} "
            f"{row['normal_count']:>8,} "
            f"{row['tumor_count']:>8,}{flag}"
        )
    print(f"  {'-' * 75}")
    print(f"  * = ≥{MIN_CELLS_FINE} cells in BOTH conditions")

    # Filter to types with ≥200 in both
    mask = (df["normal_count"] >= MIN_CELLS_FINE) & (df["tumor_count"] >= MIN_CELLS_FINE)
    df_eligible = df[mask].copy().reset_index(drop=True)
    print(f"\n  {len(df_eligible)} cell types with ≥{MIN_CELLS_FINE} cells in BOTH conditions")

    return df, df_eligible


# ===================================================================
# STEP 3 — Build fine-grained mapping and show proposed mapping
# ===================================================================


def build_fine_mapping(
    inventory_df: pd.DataFrame,
    xs_residuals: dict[str, float],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Apply FINE_RULES to map raw Census labels to fine categories, then aggregate
    and check which types pass the ≥500 gate. Print proposed mapping for user
    confirmation.

    Returns:
        mapping_df: per-raw-label mapping
        agg_df: aggregated fine-type counts with gate status and xs match info
    """
    print("\n" + "=" * 70)
    print("STEP 3 — Fine-Grained Mapping (maximizing 35-type overlap)")
    print("=" * 70)

    xs_types = set(xs_residuals.keys())

    # Apply fine classification
    rows = []
    for _, row in inventory_df.iterrows():
        raw = row["cell_type"]
        fine = _classify_fine(raw)
        rows.append({
            "raw_label": raw,
            "fine_label": fine,
            "normal_count": row["normal_count"],
            "tumor_count": row["tumor_count"],
        })

    mapping_df = pd.DataFrame(rows).sort_values(["fine_label", "raw_label"]).reset_index(drop=True)

    # AUDIT: flag every raw label containing "t cell" substring and show what it mapped to
    # This catches future substring collision bugs (ISSUE-054 safeguard)
    tcell_substr_labels = mapping_df[
        mapping_df["raw_label"].str.lower().str.contains("t cell", na=False)
    ][["raw_label", "fine_label"]].drop_duplicates()
    if len(tcell_substr_labels) > 0:
        print(f"\n  AUDIT — raw labels containing 't cell' substring:")
        for _, r in tcell_substr_labels.iterrows():
            flag = "" if r["fine_label"] in (
                "T cell", "CD4-positive, alpha-beta T cell",
                "CD8-positive, alpha-beta T cell", "mast cell",
            ) else " *** CHECK ***"
            print(f"    {r['raw_label']:<55} → {r['fine_label']}{flag}")
        print()

    # Print per-raw-label mapping
    print(f"\n  RAW LABEL → FINE CATEGORY MAPPING:")
    print(f"  {'raw_label':<55} {'fine_label':<40} {'normal':>7} {'tumor':>7}")
    print(f"  {'-' * 112}")

    current_fine = None
    for _, row in mapping_df.iterrows():
        if row["fine_label"] != current_fine:
            if current_fine is not None:
                print()
            current_fine = row["fine_label"]
        xs_match = " [XS]" if row["fine_label"] in xs_types else ""
        print(
            f"  {str(row['raw_label']):<55} "
            f"{row['fine_label']:<40} "
            f"{row['normal_count']:>7,} "
            f"{row['tumor_count']:>7,}{xs_match}"
        )

    # Aggregate to fine level
    agg = (
        mapping_df.groupby("fine_label")[["normal_count", "tumor_count"]]
        .sum()
        .reset_index()
    )
    agg["passes_500"] = (
        (agg["normal_count"] >= MIN_CELLS_GATE)
        & (agg["tumor_count"] >= MIN_CELLS_GATE)
    )
    agg["in_xs_35"] = agg["fine_label"].isin(xs_types)
    agg["xs_magnitude"] = agg["fine_label"].map(xs_residuals)

    print(f"\n\n  AGGREGATED FINE-GRAINED COUNTS:")
    print(
        f"  {'fine_label':<45} {'normal':>8} {'tumor':>8} "
        f"{'gate':>6} {'in_xs':>6} {'xs_resid':>9}"
    )
    print(f"  {'-' * 90}")
    for _, row in agg.sort_values("fine_label").iterrows():
        gate = "PASS" if row["passes_500"] else "FAIL"
        xs = "YES" if row["in_xs_35"] else "no"
        xs_mag = f"{row['xs_magnitude']:.2f}" if pd.notna(row["xs_magnitude"]) else "-"
        print(
            f"  {row['fine_label']:<45} "
            f"{row['normal_count']:>8,} "
            f"{row['tumor_count']:>8,} "
            f"{gate:>6} "
            f"{xs:>6} "
            f"{xs_mag:>9}"
        )

    n_pass = agg["passes_500"].sum()
    n_pass_xs = (agg["passes_500"] & agg["in_xs_35"]).sum()
    n_total = len(agg)
    print(f"\n  {n_pass}/{n_total} fine types pass ≥500 gate")
    print(f"  {n_pass_xs} of those match a 35-type cross-species type")
    print(f"  (Previous analysis: n=7 matched types. Target: n≥15)")

    # Show the proposed cross-analysis mapping
    print(f"\n\n  PROPOSED CROSS-ANALYSIS MAPPING (cancer → cross-species):")
    print(
        f"  {'cancer_fine_type':<45} {'xs_match':<45} {'status':>8}"
    )
    print(f"  {'-' * 100}")

    n_matched = 0
    for _, row in agg[agg["passes_500"]].sort_values("fine_label").iterrows():
        ft = row["fine_label"]
        if ft in EXCLUDE_TYPES:
            continue
        xs_list = CANCER_TO_XS_MAP.get(ft)
        if xs_list is None:
            status = "NO MATCH"
            xs_name = "-"
        elif all(x in xs_types for x in xs_list):
            status = "MATCHED"
            xs_name = " + ".join(xs_list)
            n_matched += 1
        else:
            missing = [x for x in xs_list if x not in xs_types]
            status = "MISSING"
            xs_name = f"{xs_list} (missing: {missing})"
        print(f"  {ft:<45} {xs_name:<45} {status:>8}")

    print(f"\n  → {n_matched} types will be used for cross-analysis Spearman")

    return mapping_df, agg


# ===================================================================
# STEP 4 — Download cells
# ===================================================================


def download_scaled_data(
    mapping_df: pd.DataFrame,
    valid_fine_types: list[str],
) -> tuple:
    """
    Download cells from CELLxGENE Census for the fine-grained cell types.
    Uses the same memory-efficient streaming approach as src/cancer_loader.py.

    Returns: normal_adata, tumor_adata
    """
    import anndata as ad
    import cellxgene_census
    import scanpy as sc

    print("\n" + "=" * 70)
    print("STEP 4 — Download Data (fine-grained cell types)")
    print("=" * 70)

    # Build raw→fine mapping dict
    raw_to_fine = dict(zip(mapping_df["raw_label"], mapping_df["fine_label"]))

    # Load ortholog gene IDs
    from cellwarp.cancer_loader import load_ortholog_gene_ids
    ortholog_ids = load_ortholog_gene_ids(ORTHOLOG_CSV)

    OBS_COLUMNS = [
        "cell_type", "tissue", "tissue_general", "disease",
        "donor_id", "dataset_id", "assay", "sex", "development_stage",
        "is_primary_data",
    ]
    VAR_COLUMNS = ["feature_id", "feature_name"]

    print(f"\n  Opening Census (2025-11-08)...")
    census = cellxgene_census.open_soma(census_version="2025-11-08")

    tissue_str = ", ".join(f"'{v}'" for v in TISSUE_GENERAL)
    rng = np.random.default_rng(RANDOM_SEED)

    results = {}
    for condition, disease_label in [
        ("normal", DISEASE_NORMAL),
        ("tumor", DISEASE_TUMOR),
    ]:
        value_filter = (
            f"is_primary_data == True and "
            f"tissue_general in [{tissue_str}] and "
            f"disease == '{disease_label}'"
        )

        # Phase 1: Fetch obs metadata only
        print(f"\n  [{condition}] Fetching obs metadata...")
        with census_timeout():
            obs_df = cellxgene_census.get_obs(
                census, "Homo sapiens",
                value_filter=value_filter,
                column_names=list(dict.fromkeys(["soma_joinid"] + OBS_COLUMNS)),
            )
        print(f"    Metadata: {len(obs_df):,} cells")

        # Apply fine mapping
        obs_df["fine_cell_type"] = obs_df["cell_type"].map(raw_to_fine)
        obs_df["fine_cell_type"] = obs_df["fine_cell_type"].fillna("Other")
        obs_df = obs_df[obs_df["fine_cell_type"].isin(valid_fine_types)].copy()
        print(f"    After fine type filter: {len(obs_df):,} cells")

        # Build fine→raw mapping for value_filter downloads
        fine_to_raw: dict[str, list[str]] = {}
        for raw, fine in raw_to_fine.items():
            if fine in valid_fine_types:
                fine_to_raw.setdefault(fine, []).append(raw)

        # Phase 2+3: Download per raw cell type using value_filter (fast indexed
        # lookup), then subsample per fine type. This replaces the obs_coords
        # approach which was too slow (scattered joinid reads >600s per query).
        print(f"    Downloading expression data per raw cell type (value_filter)...")
        chunks = []
        skipped_timeout: list[str] = []
        for ct in sorted(fine_to_raw.keys()):
            raw_types = fine_to_raw[ct]
            # Count available cells for this fine type
            ct_total = obs_df.loc[obs_df["fine_cell_type"] == ct].shape[0]
            ct_target = min(ct_total, MAX_CELLS_PER_TYPE)
            print(f"      {ct} ({ct_total:,} cells, target {ct_target:,}):", flush=True)

            ct_chunks = []
            ct_cells_so_far = 0
            skipped_raw: list[str] = []
            for raw_type in sorted(raw_types):
                # Count cells of this raw type in this condition
                raw_count = obs_df.loc[obs_df["cell_type"] == raw_type].shape[0]
                if raw_count == 0:
                    continue
                # Early stop: already have enough cells for this fine type
                if ct_cells_so_far >= MAX_CELLS_PER_TYPE:
                    print(f"        {raw_type} ({raw_count:,})... SKIP (already have {ct_cells_so_far:,})", flush=True)
                    continue
                # Escape single quotes in cell type name for value_filter
                raw_escaped = raw_type.replace("'", "\\'")
                raw_filter = (
                    f"is_primary_data == True and "
                    f"tissue_general in [{tissue_str}] and "
                    f"disease == '{disease_label}' and "
                    f"cell_type == '{raw_escaped}'"
                )
                print(f"        {raw_type} ({raw_count:,})...", end="", flush=True)
                try:
                    with census_timeout():
                        raw_chunk = cellxgene_census.get_anndata(
                            census=census,
                            organism="Homo sapiens",
                            obs_value_filter=raw_filter,
                            obs_column_names=OBS_COLUMNS,
                            var_column_names=VAR_COLUMNS,
                        )
                    ct_chunks.append(raw_chunk)
                    ct_cells_so_far += raw_chunk.n_obs
                    print(f" → {raw_chunk.n_obs:,} cells", flush=True)
                except TimeoutError:
                    print(f" *** TIMEOUT ({CENSUS_QUERY_TIMEOUT}s) ***", flush=True)
                    skipped_raw.append(raw_type)
                    # Reopen Census — SIGALRM may corrupt SOMA handle
                    try:
                        census.close()
                    except Exception:
                        pass
                    census = cellxgene_census.open_soma(census_version="2025-11-08")
                    # Continue to next raw type (don't skip entire fine type)

            if len(ct_chunks) == 0:
                skipped_timeout.append(ct)
                print(f"        → SKIPPED (no data downloaded)", flush=True)
                continue
            if skipped_raw:
                print(f"        NOTE: {len(skipped_raw)} raw subtype(s) timed out: {skipped_raw}", flush=True)
                print(f"        Using {ct_cells_so_far:,} cells from successful downloads", flush=True)

            # Concatenate raw type chunks for this fine type
            if len(ct_chunks) == 1:
                ct_adata = ct_chunks[0]
            else:
                ct_adata = ad.concat(ct_chunks, merge="same")
            del ct_chunks

            # Subsample to MAX_CELLS_PER_TYPE
            if ct_adata.n_obs > MAX_CELLS_PER_TYPE:
                idx = rng.choice(ct_adata.n_obs, size=MAX_CELLS_PER_TYPE, replace=False)
                idx.sort()
                ct_adata = ct_adata[idx].copy()
                print(f"        → subsampled to {ct_adata.n_obs:,} cells", flush=True)

            ct_adata.obs["coarse_cell_type"] = ct
            chunks.append(ct_adata)
            print(f"        → {ct_adata.n_obs:,} cells × {ct_adata.n_vars:,} genes", flush=True)

        if skipped_timeout:
            print(f"\n    WARNING: {len(skipped_timeout)} type(s) skipped due to timeout:")
            for st in skipped_timeout:
                print(f"      - {st}")
            print(f"    These types will be excluded from Procrustes analysis.")

        # Phase 4: Concatenate and filter to ortholog space
        print(f"    Concatenating {len(chunks)} chunks...")
        adata = ad.concat(chunks, merge="same")
        del chunks
        print(f"    Combined: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

        gene_mask = adata.var["feature_id"].isin(ortholog_ids)
        adata = adata[:, gene_mask].copy()
        adata.var.index = adata.var["feature_id"].values
        shared_ids = sorted(set(adata.var.index) & set(ortholog_ids))
        adata = adata[:, shared_ids].copy()
        print(f"    Ortholog-filtered: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

        # QC + normalize
        n_before = adata.n_obs
        sc.pp.filter_cells(adata, min_genes=200)
        n_after = adata.n_obs
        print(f"    QC (<200 genes): {n_before:,} → {n_after:,} ({n_before - n_after:,} removed)")
        sc.pp.normalize_total(adata, target_sum=10_000)
        sc.pp.log1p(adata)
        print(f"    Normalized: 10k counts + log1p")

        results[condition] = adata

    census.close()

    normal_adata = results["normal"]
    tumor_adata = results["tumor"]

    # Save with atomic writes
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    from cellwarp.cancer_loader import save_h5ad_atomic
    save_h5ad_atomic(normal_adata, DATA_DIR / "colon_normal_scaled.h5ad")
    save_h5ad_atomic(tumor_adata, DATA_DIR / "colon_tumor_scaled.h5ad")

    # Print summary
    print(f"\n  DOWNLOAD SUMMARY:")
    print(f"  {'Fine Type':<45} {'Normal':>8} {'Tumor':>8}")
    print(f"  {'-' * 65}")
    all_types = sorted(set(
        list(normal_adata.obs["coarse_cell_type"].unique()) +
        list(tumor_adata.obs["coarse_cell_type"].unique())
    ))
    for ct in all_types:
        n = (normal_adata.obs["coarse_cell_type"] == ct).sum()
        t = (tumor_adata.obs["coarse_cell_type"] == ct).sum()
        print(f"  {ct:<45} {n:>8,} {t:>8,}")
    print(f"  {'-' * 65}")
    print(f"  {'TOTAL':<45} {normal_adata.n_obs:>8,} {tumor_adata.n_obs:>8,}")
    print(f"  Genes: {normal_adata.n_vars:,}")

    # Donor counts
    print(f"\n  Donors: normal={normal_adata.obs['donor_id'].nunique()}, "
          f"tumor={tumor_adata.obs['donor_id'].nunique()}")

    print(f"\n  CHECKPOINT: Download complete.")
    return normal_adata, tumor_adata


# ===================================================================
# STEP 5 — Cancer Procrustes pipeline (identical to script 11)
# ===================================================================


def run_cancer_procrustes(
    normal_adata, tumor_adata
) -> tuple:
    """
    Run the full Procrustes pipeline: donor centroids → PCA → Procrustes →
    permutation test → deformation scores → top genes.

    Returns: cell_types, scores, residuals, result, p_value, null_dist, pca_model, gene_names, top_genes
    """
    import anndata as ad

    print("\n" + "=" * 70)
    print("STEP 5 — Cancer Procrustes Pipeline")
    print("=" * 70)

    shared_genes = list(normal_adata.var.index)

    # Identify valid cell types (≥500 in BOTH conditions)
    normal_counts = normal_adata.obs["coarse_cell_type"].value_counts()
    tumor_counts = tumor_adata.obs["coarse_cell_type"].value_counts()

    all_types = sorted(
        set(normal_counts.index) & set(tumor_counts.index) - EXCLUDE_TYPES
    )
    passed = []
    dropped = []
    for ct in all_types:
        n_n = normal_counts.get(ct, 0)
        n_t = tumor_counts.get(ct, 0)
        if n_n >= MIN_CELLS_GATE and n_t >= MIN_CELLS_GATE:
            passed.append(ct)
        else:
            dropped.append(ct)

    print(f"\n  Cell types passing ≥{MIN_CELLS_GATE} gate ({len(passed)}):")
    for ct in passed:
        print(f"    {ct:<45} normal={normal_counts[ct]:>5,}  tumor={tumor_counts[ct]:>5,}")
    if dropped:
        print(f"  Dropped ({len(dropped)}):")
        for ct in dropped:
            n = normal_counts.get(ct, 0)
            t = tumor_counts.get(ct, 0)
            print(f"    {ct:<45} normal={n:>5,}  tumor={t:>5,}")

    cell_types = passed

    # --- Per-donor centroids ---
    print(f"\n  Computing per-donor centroids...")
    normal_centroids = _donor_centroids(normal_adata, cell_types, shared_genes, "Normal")
    tumor_centroids = _donor_centroids(tumor_adata, cell_types, shared_genes, "Tumor")

    # Save centroids
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    normal_centroids.to_csv(OUTPUT_DIR / "centroids_normal_scaled.csv")
    tumor_centroids.to_csv(OUTPUT_DIR / "centroids_tumor_scaled.csv")

    # --- PCA ---
    print(f"\n  Running PCA on combined centroids...")
    normal_mat = normal_centroids.loc[cell_types].values
    tumor_mat = tumor_centroids.loc[cell_types].values
    combined = np.vstack([normal_mat, tumor_mat])

    pca = PCA(
        n_components=PCA_VARIANCE_THRESHOLD,
        svd_solver="full",
        random_state=RANDOM_SEED,
    )
    combined_pca = pca.fit_transform(combined)
    n_types = len(cell_types)
    normal_pca = combined_pca[:n_types]
    tumor_pca = combined_pca[n_types:]

    cumvar = np.cumsum(pca.explained_variance_ratio_)
    print(f"    PCA: {pca.n_components_} components, {cumvar[-1]*100:.1f}% variance")

    np.savez(
        OUTPUT_DIR / "pca_cancer_scaled.npz",
        normal_pca=normal_pca, tumor_pca=tumor_pca,
        components=pca.components_,
        explained_variance_ratio=pca.explained_variance_ratio_,
        mean=pca.mean_,
    )

    # --- Procrustes ---
    print(f"\n  Running Procrustes alignment (tumor → normal)...")
    result = procrustes_align(normal_pca, tumor_pca)
    det = np.linalg.det(result.rotation)
    assert abs(det - 1.0) < 1e-6, f"Rotation det={det}, expected +1.0"
    print(f"    Distance: {result.distance:.4f}, Scaling: {result.scaling:.6f}")

    # --- Permutation test ---
    print(f"\n  Running permutation test ({N_PERMUTATIONS:,} iterations)...")
    p_value, null_dist = permutation_test(
        normal_pca, tumor_pca, n_permutations=N_PERMUTATIONS, seed=RANDOM_SEED
    )
    np.save(OUTPUT_DIR / "null_distribution_scaled.npy", null_dist)
    null_median = float(np.median(null_dist))
    obs_null = result.distance / null_median
    print(f"    p={p_value:.6f}, obs/null={obs_null:.3f}")

    # --- Deformation scores ---
    print(f"\n  Computing deformation scores...")
    residuals = compute_residual_vectors(result, cell_types)
    scores = {ct: float(np.linalg.norm(residuals[ct])) for ct in cell_types}

    ranked = sorted(scores.items(), key=lambda x: x[1], reverse=True)
    print(f"\n  DEFORMATION RANKING:")
    print(f"  {'Rank':<6} {'Cell Type':<45} {'Score':>10}")
    print(f"  {'-' * 64}")
    for i, (ct, score) in enumerate(ranked, 1):
        print(f"  {i:<6} {ct:<45} {score:>10.4f}")

    # --- Top genes ---
    # Map Ensembl to gene names
    ref = ad.read_h5ad(DATA_DIR / "colon_normal_scaled.h5ad")
    ensembl_to_name = dict(zip(ref.var["feature_id"], ref.var["feature_name"]))
    del ref
    readable_genes = [ensembl_to_name.get(g, g) for g in shared_genes]
    top_genes = map_residuals_to_genes(residuals, pca, readable_genes, n_top=N_TOP_GENES)

    return (cell_types, scores, residuals, result, p_value, null_dist,
            pca, shared_genes, top_genes, dropped)


def _donor_centroids(
    adata, cell_types: list[str], gene_names: list[str], label: str
) -> pd.DataFrame:
    """Compute mean-of-donor-means centroids."""
    centroids = {}
    for ct in cell_types:
        mask = adata.obs["coarse_cell_type"] == ct
        ct_data = adata[mask]
        donors = ct_data.obs["donor_id"].unique()
        donor_means = []
        for d in donors:
            d_mask = ct_data.obs["donor_id"] == d
            d_cells = ct_data[d_mask]
            if d_cells.n_obs > 0:
                donor_means.append(np.asarray(d_cells.X.mean(axis=0)).flatten())
        centroids[ct] = np.mean(donor_means, axis=0)
        print(f"    {label} {ct:<45} {mask.sum():>5,} cells, {len(donors):>3} donors")
    df = pd.DataFrame(centroids, index=gene_names).T
    df.index.name = "cell_type"
    return df


# ===================================================================
# STEP 6 — Cross-analysis Spearman correlation + scatter plot
# ===================================================================


def cross_analysis_spearman(
    cancer_scores: dict[str, float],
    xs_residuals: dict[str, float],
) -> dict:
    """
    Match cancer fine types to cross-species 35-type residuals. Compute
    Spearman correlation. Save scatter plot.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    print("\n" + "=" * 70)
    print("STEP 6 — Cross-Analysis Spearman Correlation")
    print("=" * 70)

    matched = []
    print(f"\n  MATCH TABLE:")
    print(
        f"  {'cancer_type':<45} {'xs_type':<45} "
        f"{'cancer_def':>10} {'xs_resid':>10} {'quality':>8}"
    )
    print(f"  {'-' * 122}")

    for cancer_type in sorted(cancer_scores.keys()):
        xs_list = CANCER_TO_XS_MAP.get(cancer_type)
        if xs_list is None:
            print(
                f"  {cancer_type:<45} {'[no match]':<45} "
                f"{cancer_scores[cancer_type]:>10.4f} {'-':>10} {'SKIP':>8}"
            )
            continue

        # Check all xs types exist
        missing = [x for x in xs_list if x not in xs_residuals]
        if missing:
            print(
                f"  {cancer_type:<45} {str(xs_list):<45} "
                f"{cancer_scores[cancer_type]:>10.4f} {'-':>10} {'MISSING':>8}"
            )
            continue

        # Use mean if multiple xs matches
        xs_mag = np.mean([xs_residuals[x] for x in xs_list])
        xs_label = " + ".join(xs_list)
        quality = "exact" if len(xs_list) == 1 else "mean"

        matched.append({
            "cancer_type": cancer_type,
            "xs_type": xs_label,
            "cancer_deformation": cancer_scores[cancer_type],
            "xs_residual": float(xs_mag),
            "match_quality": quality,
        })
        print(
            f"  {cancer_type:<45} {xs_label:<45} "
            f"{cancer_scores[cancer_type]:>10.4f} {xs_mag:>10.4f} {quality:>8}"
        )

    n_matched = len(matched)
    print(f"\n  n = {n_matched} matched types (target was ≥15, previous was 7)")

    if n_matched < 3:
        print(f"  WARNING: Too few matches for correlation.")
        return {"n_matched": n_matched, "insufficient": True}

    # Spearman
    cancer_vals = [m["cancer_deformation"] for m in matched]
    xs_vals = [m["xs_residual"] for m in matched]

    rho, p_value = stats.spearmanr(cancer_vals, xs_vals)

    print(f"\n  Spearman ρ = {rho:.4f}")
    print(f"  p-value   = {p_value:.6f}")
    print(f"  n         = {n_matched}")
    print(f"  Significant at α=0.05: {'YES' if p_value < 0.05 else 'NO'}")
    print(f"  Significant at α=0.01: {'YES' if p_value < 0.01 else 'NO'}")

    # Interpretation
    print(f"\n  INTERPRETATION:")
    if p_value < 0.05 and rho > 0:
        print(
            "  POSITIVE CORRELATION: Evolutionarily flexible cell types are also\n"
            "  more deformed by cancer. Shared plasticity axis."
        )
    elif p_value < 0.05 and rho < 0:
        print(
            "  NEGATIVE CORRELATION: Evolutionarily rigid cell types resist tumor\n"
            "  reprogramming. Shared constraint on rewiring."
        )
    else:
        print(
            "  NO SIGNIFICANT CORRELATION: Evolutionary rigidity and cancer\n"
            "  deformation appear independent."
        )

    # ---------------------------------------------------------------
    # SENSITIVITY 1: Spearman without "T cell" (ambiguous match)
    # ---------------------------------------------------------------
    # "T cell" in cancer is a catch-all for gamma-delta T, NKT, and
    # unresolved T cells. Its cross-species counterpart "T cell" is
    # similarly heterogeneous. Flag whether removing it changes the result.
    matched_no_tcell = [m for m in matched if m["cancer_type"] != "T cell"]
    if len(matched_no_tcell) >= 3:
        c_vals_nt = [m["cancer_deformation"] for m in matched_no_tcell]
        x_vals_nt = [m["xs_residual"] for m in matched_no_tcell]
        rho_nt, p_nt = stats.spearmanr(c_vals_nt, x_vals_nt)
        print(f"\n  SENSITIVITY — Without 'T cell' (ambiguous match):")
        print(f"    Spearman ρ = {rho_nt:.4f}, p = {p_nt:.6f}, n = {len(matched_no_tcell)}")
        print(f"    Direction stable: {'YES' if (rho_nt > 0) == (rho > 0) else 'NO'}")
        delta_rho = abs(rho_nt - rho)
        print(f"    |Δρ| = {delta_rho:.4f} ({'sensitive' if delta_rho > 0.15 else 'robust'})")
    else:
        rho_nt, p_nt = None, None
        print(f"\n  SENSITIVITY — Cannot compute without T cell (too few remaining).")

    # ---------------------------------------------------------------
    # SENSITIVITY 2: Enterocyte low tumor count flag
    # ---------------------------------------------------------------
    # Enterocyte has only 1,347 tumor cells (vs 19,936 normal). The
    # centroid may be less reliable. Check if removing it changes result.
    entero_type = "enterocyte of epithelium of large intestine"
    matched_no_entero = [m for m in matched if m["cancer_type"] != entero_type]
    if len(matched_no_entero) >= 3:
        c_vals_ne = [m["cancer_deformation"] for m in matched_no_entero]
        x_vals_ne = [m["xs_residual"] for m in matched_no_entero]
        rho_ne, p_ne = stats.spearmanr(c_vals_ne, x_vals_ne)
        print(f"\n  SENSITIVITY — Without enterocyte (low tumor n=1,347):")
        print(f"    Spearman ρ = {rho_ne:.4f}, p = {p_ne:.6f}, n = {len(matched_no_entero)}")
        print(f"    Direction stable: {'YES' if (rho_ne > 0) == (rho > 0) else 'NO'}")
        delta_rho_e = abs(rho_ne - rho)
        print(f"    |Δρ| = {delta_rho_e:.4f} ({'sensitive' if delta_rho_e > 0.15 else 'robust'})")
    else:
        rho_ne, p_ne = None, None

    # ---------------------------------------------------------------
    # SENSITIVITY 3: Without BOTH T cell and enterocyte
    # ---------------------------------------------------------------
    matched_strict = [
        m for m in matched
        if m["cancer_type"] not in {"T cell", entero_type}
    ]
    if len(matched_strict) >= 3:
        c_vals_s = [m["cancer_deformation"] for m in matched_strict]
        x_vals_s = [m["xs_residual"] for m in matched_strict]
        rho_s, p_s = stats.spearmanr(c_vals_s, x_vals_s)
        print(f"\n  SENSITIVITY — Without T cell AND enterocyte (strictest):")
        print(f"    Spearman ρ = {rho_s:.4f}, p = {p_s:.6f}, n = {len(matched_strict)}")
        print(f"    Direction stable: {'YES' if (rho_s > 0) == (rho > 0) else 'NO'}")
    else:
        rho_s, p_s = None, None

    # ---------------------------------------------------------------
    # Enterocyte cell count warning
    # ---------------------------------------------------------------
    print(f"\n  NOTE on enterocyte:")
    print(f"    Tumor cell count (1,347) is low relative to normal (19,936).")
    print(f"    This 15:1 imbalance means the tumor centroid relies on fewer")
    print(f"    donors and may be noisier. Per-donor averaging mitigates this,")
    print(f"    but enterocyte deformation score should be interpreted with")
    print(f"    caution. Sensitivity analysis above quantifies the impact.")

    # Scatter plot
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.scatter(xs_vals, cancer_vals, s=60, alpha=0.8, edgecolors="black", linewidths=0.5)

    for m in matched:
        # Shorten label for plot
        label = m["cancer_type"]
        if len(label) > 30:
            label = label[:27] + "..."
        ax.annotate(
            label,
            (m["xs_residual"], m["cancer_deformation"]),
            fontsize=7,
            xytext=(5, 5),
            textcoords="offset points",
        )

    # Regression line
    z = np.polyfit(xs_vals, cancer_vals, 1)
    x_line = np.linspace(min(xs_vals) * 0.95, max(xs_vals) * 1.05, 100)
    ax.plot(x_line, np.polyval(z, x_line), "--", color="red", alpha=0.5, linewidth=1)

    ax.set_xlabel("Cross-Species Residual Magnitude (35-type)", fontsize=11)
    ax.set_ylabel("Cancer Deformation Score (tumor vs normal)", fontsize=11)
    ax.set_title(
        f"Evolutionary Rigidity vs Cancer Deformation\n"
        f"Spearman ρ={rho:.3f}, p={p_value:.4f}, n={n_matched}",
        fontsize=12,
    )
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    plot_path = OUTPUT_DIR / "cross_analysis_scaled.png"
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)
    print(f"\n  Saved scatter plot: {plot_path}")
    print(f"  Plot description: Scatter of {n_matched} cell types. X-axis is cross-species")
    print(f"  residual magnitude, Y-axis is cancer deformation score. Red dashed regression line.")
    if rho > 0:
        print(f"  Points trend upward — flexible types deform more in both contexts.")
    elif rho < 0:
        print(f"  Points trend downward — rigid types resist deformation in both contexts.")

    # Save match table as CSV
    match_df = pd.DataFrame(matched)
    match_df.to_csv(OUTPUT_DIR / "cross_analysis_match_table.csv", index=False)
    print(f"  Saved match table: {OUTPUT_DIR / 'cross_analysis_match_table.csv'}")

    sensitivity = {}
    if rho_nt is not None:
        sensitivity["without_T_cell"] = {
            "rho": float(rho_nt), "p": float(p_nt),
            "n": len(matched_no_tcell),
        }
    if rho_ne is not None:
        sensitivity["without_enterocyte"] = {
            "rho": float(rho_ne), "p": float(p_ne),
            "n": len(matched_no_entero),
        }
    if rho_s is not None:
        sensitivity["without_both"] = {
            "rho": float(rho_s), "p": float(p_s),
            "n": len(matched_strict),
        }

    return {
        "n_matched": n_matched,
        "matched_types": matched,
        "spearman_rho": float(rho),
        "spearman_p": float(p_value),
        "significant_005": bool(p_value < 0.05),
        "significant_001": bool(p_value < 0.01),
        "sensitivity": sensitivity,
    }


# ===================================================================
# Save all results
# ===================================================================


def save_results(
    result, p_value, null_dist, scores, residuals, top_genes,
    cell_types, pca_model, cross_corr, dropped,
):
    """Save comprehensive results JSON."""
    from cellwarp.procrustes import _procrustes_distance

    null_median = float(np.median(null_dist))

    results_dict = {
        "procrustes": {
            "distance": float(result.distance),
            "distance_squared": float(result.distance_squared),
            "scaling": float(result.scaling),
            "rotation_det": float(np.linalg.det(result.rotation)),
        },
        "permutation_test": {
            "p_value": float(p_value),
            "n_permutations": len(null_dist),
            "null_median": null_median,
            "obs_null_ratio": float(result.distance / null_median),
            "significant_at_001": bool(p_value < 0.01),
        },
        "pca": {
            "n_components": int(pca_model.n_components_),
            "cumulative_variance": float(np.sum(pca_model.explained_variance_ratio_)),
            "per_component_variance": pca_model.explained_variance_ratio_.tolist(),
        },
        "cell_types": cell_types,
        "dropped_types": dropped,
        "deformation_scores": {ct: float(scores[ct]) for ct in cell_types},
        "deformation_ranking": [
            {"rank": i + 1, "cell_type": ct, "score": float(scores[ct])}
            for i, (ct, _) in enumerate(
                sorted(scores.items(), key=lambda x: x[1], reverse=True)
            )
        ],
        "residuals": {
            ct: {
                "vector_pca": residuals[ct].tolist(),
                "magnitude": float(np.linalg.norm(residuals[ct])),
            }
            for ct in cell_types
        },
        "top_genes_per_cell_type": {
            ct: top_genes[ct][["gene", "loading", "abs_loading", "rank"]]
            .to_dict(orient="records")
            for ct in cell_types
        },
        "cross_analysis_correlation": cross_corr,
        "random_seed": RANDOM_SEED,
    }

    output_path = OUTPUT_DIR / "cancer_scaled_results.json"
    tmp_path = output_path.with_suffix(".tmp")
    with open(tmp_path, "w") as f:
        json.dump(results_dict, f, indent=2)
    tmp_path.rename(output_path)
    print(f"\n  Saved: {output_path}")


# ===================================================================
# Summary
# ===================================================================


def print_summary(
    result, p_value, null_dist, scores, cell_types,
    cross_corr, dropped, pca_model,
):
    """Print analyst-format summary."""
    null_median = float(np.median(null_dist))
    obs_null = result.distance / null_median

    print("\n" + "=" * 70)
    print("SCALED CANCER PROCRUSTES — SUMMARY")
    print("=" * 70)

    print(f"""
1. WHAT WAS DONE
   Scaled cancer Procrustes with fine-grained cell types (not pooled).
   {len(cell_types)} fine cell types as geometric landmarks.
   Per-donor centroids averaged to control for donor imbalance.
   Dropped: {dropped if dropped else 'none'}.
   Gene space: {pca_model.mean_.shape[0]:,} ortholog genes.

2. KEY NUMBERS
   Procrustes distance:     {result.distance:.4f}
   Scaling factor:          {result.scaling:.6f}
   Permutation p-value:     {p_value:.6f}
   Null median:             {null_median:.4f}
   Obs/null ratio:          {obs_null:.3f}
   PCA components:          {pca_model.n_components_}
   Significant at α=0.01:   {'YES' if p_value < 0.01 else 'NO'}

3. DEFORMATION RANKING""")

    ranked = sorted(scores.items(), key=lambda x: x[1], reverse=True)
    print(f"   {'Rank':<6} {'Cell Type':<45} {'Score':>10}")
    print(f"   {'-' * 64}")
    for i, (ct, score) in enumerate(ranked, 1):
        print(f"   {i:<6} {ct:<45} {score:>10.4f}")

    print(f"""
4. CROSS-ANALYSIS (CRITICAL)
   Matched types:   {cross_corr.get('n_matched', 'N/A')}
   Spearman ρ:      {cross_corr.get('spearman_rho', 'N/A')}
   p-value:         {cross_corr.get('spearman_p', 'N/A')}
   Significant:     {cross_corr.get('significant_005', 'N/A')}

5. COMPARISON WITH ORIGINAL 8-TYPE ANALYSIS
   Original: ρ=0.679, p=0.094, n=7
   Scaled:   ρ={cross_corr.get('spearman_rho', '?')}, p={cross_corr.get('spearman_p', '?')}, n={cross_corr.get('n_matched', '?')}""")

    sens = cross_corr.get("sensitivity", {})
    if sens:
        print(f"\n6. SENSITIVITY ANALYSIS")
        if "without_T_cell" in sens:
            s = sens["without_T_cell"]
            print(f"   Without T cell (ambiguous):    ρ={s['rho']:.4f}, p={s['p']:.6f}, n={s['n']}")
        if "without_enterocyte" in sens:
            s = sens["without_enterocyte"]
            print(f"   Without enterocyte (low tumor): ρ={s['rho']:.4f}, p={s['p']:.6f}, n={s['n']}")
        if "without_both" in sens:
            s = sens["without_both"]
            print(f"   Without both (strictest):       ρ={s['rho']:.4f}, p={s['p']:.6f}, n={s['n']}")

    if cross_corr.get("spearman_p") is not None:
        p = cross_corr["spearman_p"]
        rho = cross_corr["spearman_rho"]
        if p < 0.05 and rho > 0:
            print(
                "\n   POSITIVE CORRELATION CONFIRMED with adequate power.\n"
                "   Evolutionary flexibility predicts cancer vulnerability."
            )
        elif p < 0.05 and rho < 0:
            print(
                "\n   NEGATIVE CORRELATION CONFIRMED with adequate power.\n"
                "   Evolutionary rigidity predicts cancer resistance."
            )
        elif p >= 0.05:
            print(
                "\n   NOT SIGNIFICANT even with expanded n.\n"
                "   Evolutionary rigidity and cancer deformation are independent."
            )


# ===================================================================
# Main — two-phase execution
# ===================================================================


def main():
    print("\n" + "#" * 70)
    print("# CellWarp — Scaled Cancer Procrustes Pipeline")
    print("# (Fine-grained types, maximizing 35-type overlap)")
    print("#" * 70 + "\n")

    # Phase A: Inventory and mapping (steps 1-3)
    xs_residuals = load_cross_species_residuals()
    full_inventory_df, inventory_df = query_inventory()

    # Step 2.5 — Collision audit on ALL raw labels (before mapping)
    all_raw_labels = full_inventory_df["cell_type"].tolist()
    audit_substring_collisions(all_raw_labels)

    mapping_df, agg_df = build_fine_mapping(inventory_df, xs_residuals)

    # Get valid types (pass ≥500 gate, not Other)
    valid_types = sorted(
        agg_df.loc[agg_df["passes_500"] & ~agg_df["fine_label"].isin(EXCLUDE_TYPES), "fine_label"]
        .tolist()
    )

    print("\n" + "#" * 70)
    print("# PAUSED — Waiting for user confirmation of mapping")
    print(f"# {len(valid_types)} fine types pass gate. Review above and confirm.")
    print("#" * 70)
    print("\nTo continue, call main_phase_b() with the mapping_df and valid_types.")
    print("Or re-run with --continue flag after confirming.\n")

    # Check for --continue flag
    if "--continue" in sys.argv:
        print("  --continue flag detected. Proceeding with download + analysis...\n")
        main_phase_b(mapping_df, valid_types, xs_residuals)
    else:
        # Save state for phase B
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        mapping_df.to_csv(OUTPUT_DIR / "fine_mapping.csv", index=False)
        agg_df.to_csv(OUTPUT_DIR / "fine_aggregated.csv", index=False)
        pd.Series(valid_types).to_csv(OUTPUT_DIR / "valid_types.csv", index=False, header=["fine_type"])
        print(f"  Saved mapping state to {OUTPUT_DIR}/")
        print(f"  Re-run with: python scripts/12_cancer_scaled.py --continue")


def main_phase_b(mapping_df, valid_types, xs_residuals):
    """Phase B: download, Procrustes, cross-analysis (steps 4-6)."""

    # If called from --continue, reload saved state
    if mapping_df is None:
        mapping_df = pd.read_csv(OUTPUT_DIR / "fine_mapping.csv")
        valid_types = pd.read_csv(OUTPUT_DIR / "valid_types.csv")["fine_type"].tolist()
        xs_residuals = load_cross_species_residuals()

    # Step 4 — Download
    normal_adata, tumor_adata = download_scaled_data(mapping_df, valid_types)

    # Step 5 — Procrustes
    (cell_types, scores, residuals, result, p_value, null_dist,
     pca_model, gene_names, top_genes, dropped) = run_cancer_procrustes(
        normal_adata, tumor_adata
    )
    del normal_adata, tumor_adata  # free memory

    # Step 6 — Cross-analysis
    cross_corr = cross_analysis_spearman(scores, xs_residuals)

    # Save results
    save_results(
        result, p_value, null_dist, scores, residuals, top_genes,
        cell_types, pca_model, cross_corr, dropped,
    )

    # Summary
    print_summary(
        result, p_value, null_dist, scores, cell_types,
        cross_corr, dropped, pca_model,
    )

    print(f"\n  All outputs saved to {OUTPUT_DIR}/")
    print("  Done.\n")


if __name__ == "__main__":
    main()
