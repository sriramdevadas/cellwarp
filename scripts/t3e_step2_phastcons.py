#!/usr/bin/env python3
"""
T3-E Step 2: Regulatory Sequence Conservation vs Procrustes Rigidity

Tests whether cell types with more evolutionarily constrained regulatory
sequences at their identity-gene loci show higher cross-species geometric
rigidity (lower Procrustes residual magnitude).

Uses phastCons conservation scores (UCSC, hg38 coordinates) at promoter
windows around identity/loading genes, correlated with Procrustes rigidity
across 35 Tabula cell types.

Biology: phastCons measures the probability that a nucleotide is in a
conserved element, from a multiple alignment of vertebrate/mammal genomes.
High phastCons at a promoter = regulatory architecture under purifying
selection across species.

Math: Spearman rank correlation between per-cell-type mean phastCons
(at promoter windows) and Procrustes residual magnitude. Two scoring
approaches: Option A (cell-type-specific top-50 loading genes) and
Option B (shared top-200 identity genes, expression-weighted).
"""

import json
import os
import sys
import time
import traceback

import numpy as np
import pandas as pd
import requests
from scipy import stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# ============================================================
# Configuration
# ============================================================
BASE_DIR = str(Path(__file__).resolve().parent.parent)
OUTPUT_DIR = os.path.join(BASE_DIR, 'output/validation/t3e_chromatin')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Input data
CENTROIDS_HUMAN = os.path.join(BASE_DIR, 'output/phase2/scaled_35types/centroids_human_35.csv')
RESIDUALS_RANKED = os.path.join(BASE_DIR, 'output/phase2/scaled_35types/residuals_ranked.csv')
GENE_LOADINGS = os.path.join(BASE_DIR, 'output/phase2/dnds/gene_loadings_all.csv')
ORTHOLOGS = os.path.join(BASE_DIR, 'data/phase1/orthologs_human_mouse.csv')

# BigWig files
BW_PLACENTAL = os.path.join(BASE_DIR, 'data/ucsc/phastCons_placental.bw')
BW_100WAY = os.path.join(BASE_DIR, 'data/ucsc/phastCons100way.bw')

# Parameters
TOP_N_A = 50   # Top genes per cell type by loading magnitude (Option A)
TOP_N_B = 200  # Top genes by variance across centroids (Option B)
WINDOW_SIZES_KB = [1, 2, 5]

# T3-C clean types (highlight in plots)
CLEAN_T3C = {
    'natural killer cell', 'B cell', 'plasma cell',
    'endothelial cell', 'CD8-positive, alpha-beta T cell'
}

# Pre-registered thresholds
RHO_POSITIVE = 0.50
RHO_NULL = 0.35

# Output files
TSS_BED = os.path.join(OUTPUT_DIR, 'identity_gene_tss_hg38.bed')
SCORES_CSV = os.path.join(OUTPUT_DIR, 'conservation_scores.csv')
MERGED_CSV = os.path.join(OUTPUT_DIR, 'rigidity_conservation_merged.csv')
PRIMARY_JSON = os.path.join(OUTPUT_DIR, 'spearman_primary_result.json')
SENSITIVITY_CSV = os.path.join(OUTPUT_DIR, 'sensitivity_table.csv')
DOWNLOAD_LOG = os.path.join(OUTPUT_DIR, 'download_log.txt')


# ============================================================
# Step 2b: Identify gene sets + retrieve TSS coordinates
# ============================================================

def load_gene_sets():
    """Identify Option A (per-type top-50 loading genes) and Option B
    (shared top-200 identity genes) gene sets.

    Returns:
        option_a_genes: dict {cell_type: list of ensembl_gene_ids (top 50)}
        option_b_genes: list of ensembl_gene_ids (top 200 by variance)
        gene_id_to_name: dict mapping ensembl_id to gene_name
        centroids_df: human centroid DataFrame for Option B weighting
    """
    print("\n=== Step 2b: Loading gene sets ===")

    # Load gene loadings (Option A)
    loadings = pd.read_csv(GENE_LOADINGS)
    gene_id_to_name = dict(zip(loadings['ensembl_gene_id'], loadings['gene_name']))

    # Cell type columns (everything except first two and last)
    cell_type_cols = [c for c in loadings.columns
                      if c not in ('ensembl_gene_id', 'gene_name',
                                   'mean_expression_divergence')]

    # Option A: per-cell-type top 50 by loading magnitude
    option_a = {}
    for ct in cell_type_cols:
        top_genes = (loadings[['ensembl_gene_id', ct]]
                     .sort_values(ct, ascending=False)
                     .head(TOP_N_A)['ensembl_gene_id'].tolist())
        option_a[ct] = top_genes
    print(f"  Option A: {len(cell_type_cols)} cell types × {TOP_N_A} genes each")

    # Option B: top 200 by variance across human centroids
    print("  Loading human centroids for Option B...")
    centroids = pd.read_csv(CENTROIDS_HUMAN, index_col=0)
    # Columns are ensembl gene IDs, rows are cell types
    gene_variance = centroids.var(axis=0)  # variance across 35 cell types
    top_200 = gene_variance.sort_values(ascending=False).head(TOP_N_B).index.tolist()
    print(f"  Option B: top {TOP_N_B} genes by centroid variance")

    # Collect all unique gene IDs needed
    all_genes = set(top_200)
    for genes in option_a.values():
        all_genes.update(genes)
    print(f"  Total unique genes needing TSS: {len(all_genes)}")

    return option_a, top_200, gene_id_to_name, centroids, list(all_genes)


def fetch_tss_coordinates(ensembl_ids, gene_id_to_name):
    """Retrieve TSS coordinates from Ensembl REST API (batch POST).

    For + strand genes: TSS = gene start (most 5' position).
    For - strand genes: TSS = gene end (most 3' position = most 5' in
    transcript sense).

    Returns DataFrame with columns: ensembl_id, gene_name, chrom, tss, strand
    """
    print(f"\n=== Fetching TSS coordinates for {len(ensembl_ids)} genes ===")

    # Check for cached results
    if os.path.exists(TSS_BED):
        existing = pd.read_csv(TSS_BED, sep='\t', header=None,
                               names=['chrom', 'start', 'end', 'gene_name',
                                      'score', 'strand'],
                               comment='#')
        if len(existing) > 0:
            print(f"  Found cached TSS BED with {len(existing)} entries")
            # Build lookup from gene name to coordinates
            return existing

    url = "https://rest.ensembl.org/lookup/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    results = {}
    batch_size = 200
    ids_list = list(ensembl_ids)

    for i in range(0, len(ids_list), batch_size):
        batch = ids_list[i:i + batch_size]
        payload = {"ids": batch}

        success = False
        for attempt in range(3):
            try:
                resp = requests.post(url, json=payload, headers=headers, timeout=30)
                if resp.status_code == 200:
                    data = resp.json()
                    for eid, info in data.items():
                        if info and isinstance(info, dict) and 'seq_region_name' in info:
                            chrom = str(info['seq_region_name'])
                            strand = info.get('strand', 0)
                            start = info.get('start', 0)
                            end = info.get('end', 0)

                            # Valid chromosomes only
                            if chrom in [str(x) for x in range(1, 23)] + ['X', 'Y']:
                                tss = start if strand == 1 else end
                                results[eid] = {
                                    'ensembl_id': eid,
                                    'chrom': f'chr{chrom}',
                                    'tss': tss,
                                    'strand': '+' if strand == 1 else '-',
                                    'gene_name': gene_id_to_name.get(eid, eid)
                                }
                    success = True
                    break
                elif resp.status_code == 429:
                    wait = float(resp.headers.get('Retry-After', 2 ** (attempt + 1)))
                    print(f"  Rate limited, waiting {wait}s...")
                    time.sleep(wait)
                else:
                    print(f"  HTTP {resp.status_code} on attempt {attempt+1}")
                    time.sleep(2 ** (attempt + 1))
            except Exception as e:
                print(f"  Error on attempt {attempt+1}: {e}")
                time.sleep(2 ** (attempt + 1))

        if not success:
            print(f"  WARNING: Batch {i//batch_size + 1} failed after 3 retries")

        batch_num = i // batch_size + 1
        total_batches = (len(ids_list) - 1) // batch_size + 1
        print(f"  Batch {batch_num}/{total_batches}: "
              f"{len(results)}/{i + len(batch)} resolved")

        if i + batch_size < len(ids_list):
            time.sleep(0.5)

    # Filter out MT genes (task says exclude)
    results = {k: v for k, v in results.items() if v['chrom'] != 'chrMT'}

    print(f"  Successfully resolved: {len(results)}/{len(ids_list)} genes")
    failed = set(ids_list) - set(results.keys())
    if failed:
        print(f"  Failed/MT-excluded: {len(failed)} genes")

    # Save as BED6
    bed_rows = []
    for eid, info in results.items():
        # BED uses 0-based half-open coordinates
        # For the BED file, store TSS ± 2kb as the default window
        tss_0based = info['tss'] - 1  # Convert to 0-based
        bed_rows.append([
            info['chrom'],
            max(0, tss_0based - 2000),
            tss_0based + 2000,
            info['gene_name'],
            0,
            info['strand']
        ])

    bed_df = pd.DataFrame(bed_rows,
                           columns=['chrom', 'start', 'end', 'gene_name',
                                    'score', 'strand'])
    bed_df.to_csv(TSS_BED, sep='\t', header=False, index=False)
    print(f"  Saved TSS BED: {TSS_BED} ({len(bed_df)} rows)")

    return results


# ============================================================
# Step 2d: Conservation score computation
# ============================================================

def compute_conservation_scores(tss_data, option_a_genes, option_b_genes,
                                 centroids_df, gene_id_to_name):
    """Compute mean phastCons at promoter windows for each cell type.

    Option A: Mean phastCons over promoter windows of cell-type-specific
    top-50 loading genes (unweighted).

    Option B: Expression-weighted mean phastCons over shared top-200
    identity genes. Score = sum(phastCons_i * expr_i) / sum(expr_i).

    Computes for 2 tracks × 2 options × 3 window sizes = 12 scores per cell type.
    Total: 35 × 12 = 420 scores.

    Returns DataFrame with all scores.
    """
    print("\n=== Step 2d: Computing conservation scores ===")

    # Check for partial results
    if os.path.exists(SCORES_CSV):
        existing = pd.read_csv(SCORES_CSV)
        if len(existing) > 0:
            print(f"  Found partial results: {len(existing)} rows")
            completed_combos = set(zip(existing['cell_type'], existing['track'],
                                        existing['option'],
                                        existing['window_kb'].astype(str)))
        else:
            completed_combos = set()
    else:
        existing = pd.DataFrame()
        completed_combos = set()

    tracks = {
        'placental_20way': BW_PLACENTAL,
        '100way_vertebrate': BW_100WAY
    }

    # Build TSS lookup: ensembl_id -> {chrom, tss, strand}
    tss_lookup = {}
    if isinstance(tss_data, dict):
        tss_lookup = tss_data
    elif isinstance(tss_data, pd.DataFrame):
        # Rebuild from BED
        orthologs = pd.read_csv(ORTHOLOGS)
        name_to_id = dict(zip(orthologs['human_gene_name'],
                               orthologs['human_ensembl_id']))
        for _, row in tss_data.iterrows():
            gname = row['gene_name']
            eid = name_to_id.get(gname, gname)
            tss_lookup[eid] = {
                'chrom': row['chrom'],
                'tss': (row['start'] + row['end']) // 2,  # Center of BED interval
                'strand': row['strand'],
                'gene_name': gname
            }

    # Get cell type list from gene loadings
    loadings = pd.read_csv(GENE_LOADINGS)
    cell_type_cols = [c for c in loadings.columns
                      if c not in ('ensembl_gene_id', 'gene_name',
                                   'mean_expression_divergence')]

    # Build gene name to ensembl_id reverse lookup
    name_to_ensembl = {}
    orthologs = pd.read_csv(ORTHOLOGS)
    for _, row in orthologs.iterrows():
        name_to_ensembl[row['human_gene_name']] = row['human_ensembl_id']

    all_rows = []
    if len(existing) > 0:
        all_rows = existing.to_dict('records')

    for track_name, bw_path in tracks.items():
        if not os.path.exists(bw_path):
            print(f"  SKIP: {bw_path} not found (download pending)")
            continue

        print(f"\n  Track: {track_name}")
        try:
            bw = __import__('pyBigWig').open(bw_path)
        except Exception as e:
            print(f"  ERROR opening {bw_path}: {e}")
            continue

        chrom_sizes = bw.chroms()

        for ct in cell_type_cols:
            for option_label, gene_set in [('A', option_a_genes.get(ct, [])),
                                            ('B', option_b_genes)]:
                for wkb in WINDOW_SIZES_KB:
                    combo = (ct, track_name, option_label, str(wkb))
                    if combo in completed_combos:
                        continue

                    scores = []
                    n_found = 0
                    n_extracted = 0
                    n_failed = 0

                    if option_label == 'A':
                        # Option A: unweighted mean over top-50 loading genes
                        for eid in gene_set:
                            if eid not in tss_lookup:
                                n_failed += 1
                                continue
                            n_found += 1
                            info = tss_lookup[eid]
                            chrom = info['chrom']
                            tss = info['tss']

                            half = wkb * 1000
                            start = max(0, tss - half)
                            end = tss + half
                            if chrom in chrom_sizes:
                                end = min(end, chrom_sizes[chrom])
                            else:
                                n_failed += 1
                                continue

                            try:
                                vals = bw.stats(chrom, start, end, type="mean")
                                if vals and vals[0] is not None:
                                    scores.append(vals[0])
                                    n_extracted += 1
                                else:
                                    n_failed += 1
                            except Exception:
                                n_failed += 1

                        mean_score = np.mean(scores) if scores else np.nan

                    else:
                        # Option B: expression-weighted mean
                        weights = []
                        phastcons_vals = []

                        for eid in gene_set:
                            if eid not in tss_lookup:
                                n_failed += 1
                                continue
                            n_found += 1

                            # Get expression weight from centroids
                            if eid in centroids_df.columns:
                                expr = centroids_df.loc[ct, eid] if ct in centroids_df.index else 0
                            else:
                                expr = 0
                            if expr <= 0:
                                n_found += 1  # gene found but not expressed
                                continue

                            info = tss_lookup[eid]
                            chrom = info['chrom']
                            tss = info['tss']

                            half = wkb * 1000
                            start = max(0, tss - half)
                            end = tss + half
                            if chrom in chrom_sizes:
                                end = min(end, chrom_sizes[chrom])
                            else:
                                n_failed += 1
                                continue

                            try:
                                vals = bw.stats(chrom, start, end, type="mean")
                                if vals and vals[0] is not None:
                                    phastcons_vals.append(vals[0])
                                    weights.append(expr)
                                    n_extracted += 1
                                else:
                                    n_failed += 1
                            except Exception:
                                n_failed += 1

                        if phastcons_vals and weights:
                            w = np.array(weights)
                            v = np.array(phastcons_vals)
                            mean_score = np.sum(v * w) / np.sum(w)
                        else:
                            mean_score = np.nan

                    row = {
                        'cell_type': ct,
                        'track': track_name,
                        'option': option_label,
                        'window_kb': wkb,
                        'mean_phastCons': mean_score,
                        'n_genes_found': n_found,
                        'n_windows_extracted': n_extracted,
                        'n_windows_failed': n_failed
                    }
                    all_rows.append(row)

            # Save incrementally after each cell type
            pd.DataFrame(all_rows).to_csv(SCORES_CSV, index=False)

        bw.close()

    scores_df = pd.DataFrame(all_rows)
    scores_df.to_csv(SCORES_CSV, index=False)
    print(f"\n  Saved: {SCORES_CSV} ({len(scores_df)} rows)")

    # Flag low coverage
    option_a_scores = scores_df[scores_df['option'] == 'A']
    low_cov = option_a_scores[option_a_scores['n_genes_found'] < 40]
    if len(low_cov) > 0:
        print(f"\n  WARNING: {len(low_cov)} entries with <80% gene coverage (n<40):")
        for _, r in low_cov.iterrows():
            print(f"    {r['cell_type']} | {r['track']} | {r['window_kb']}kb: "
                  f"n={r['n_genes_found']}")

    return scores_df


# ============================================================
# Step 2e: Rigidity score retrieval
# ============================================================

def load_rigidity_scores():
    """Load existing Procrustes rigidity scores (residual magnitudes).

    Lower residual = more rigid. We create a rigidity_score = -residual
    so that higher = more rigid (matching the hypothesis direction).

    Returns DataFrame with cell_type, residual_magnitude, rigidity_score,
    rigidity_rank (1=most rigid).
    """
    print("\n=== Step 2e: Loading rigidity scores ===")
    df = pd.read_csv(RESIDUALS_RANKED)

    # rigidity_rank: 1 = most rigid (lowest residual = CD8+ T = rank 35 in CSV)
    df['rigidity_rank'] = df['rank'].max() + 1 - df['rank']
    df['rigidity_score'] = -df['residual_magnitude']  # higher = more rigid

    print(f"  Loaded {len(df)} cell types")
    print(f"  Most rigid (rank 1): {df.loc[df['rigidity_rank']==1, 'cell_type'].values[0]} "
          f"(residual={df.loc[df['rigidity_rank']==1, 'residual_magnitude'].values[0]:.3f})")
    print(f"  Least rigid (rank 35): {df.loc[df['rigidity_rank']==35, 'cell_type'].values[0]} "
          f"(residual={df.loc[df['rigidity_rank']==35, 'residual_magnitude'].values[0]:.3f})")

    rigidity = df[['cell_type', 'residual_magnitude', 'rigidity_score', 'rigidity_rank']]
    return rigidity


# ============================================================
# Step 2f: Spearman correlation — primary analysis
# ============================================================

def run_primary_spearman(scores_df, rigidity_df):
    """Run primary Spearman: placental track, Option A, 2kb window.

    Reports ρ, p-value, 95% CI (Fisher z), and pre-registered threshold.
    """
    print("\n=== Step 2f: Primary Spearman correlation ===")

    # Filter primary configuration
    primary = scores_df[
        (scores_df['track'] == 'placental_20way') &
        (scores_df['option'] == 'A') &
        (scores_df['window_kb'] == 2)
    ].copy()

    if len(primary) == 0:
        print("  ERROR: No primary scores found")
        return None

    # Merge with rigidity
    merged = primary.merge(rigidity_df, on='cell_type', how='inner')
    n = len(merged)
    print(f"  Merged: n={n} cell types")

    if n < 35:
        missing = set(rigidity_df['cell_type']) - set(merged['cell_type'])
        print(f"  WARNING: {35-n} cell types dropped: {missing}")

    # Drop NaN conservation scores
    merged_clean = merged.dropna(subset=['mean_phastCons'])
    if len(merged_clean) < n:
        print(f"  WARNING: {n - len(merged_clean)} cell types with NaN phastCons dropped")
    n = len(merged_clean)

    # Save merged for inspection
    merged_clean.to_csv(MERGED_CSV, index=False)
    print(f"  Saved merged data: {MERGED_CSV}")

    # Spearman: conservation vs rigidity_score (higher = more rigid)
    rho, p = stats.spearmanr(merged_clean['mean_phastCons'],
                              merged_clean['rigidity_score'])

    # 95% CI via Fisher z-transform
    z = np.arctanh(rho)
    se = 1.0 / np.sqrt(n - 3) if n > 3 else np.inf
    z_lo = z - 1.96 * se
    z_hi = z + 1.96 * se
    ci_lo = np.tanh(z_lo)
    ci_hi = np.tanh(z_hi)

    # Determine threshold
    if rho >= RHO_POSITIVE and p < 0.05:
        conclusion = "POSITIVE"
        threshold_msg = f"ρ={rho:.3f} ≥ {RHO_POSITIVE} AND p={p:.4f} < 0.05 → POSITIVE"
    elif rho >= RHO_NULL:
        conclusion = "TREND"
        threshold_msg = (f"ρ={rho:.3f} in [{RHO_NULL}, {RHO_POSITIVE}) → "
                         "TREND (underpowered or partial)")
    else:
        conclusion = "NULL_TRIGGERED"
        threshold_msg = (f"ρ={rho:.3f} < {RHO_NULL} → "
                         "TRIGGERED — 8th null, close chromatin hypothesis")

    print(f"\n  PRIMARY RESULT:")
    print(f"  Spearman ρ = {rho:.4f}")
    print(f"  p-value    = {p:.6f}")
    print(f"  n          = {n}")
    print(f"  95% CI     = [{ci_lo:.4f}, {ci_hi:.4f}]")
    print(f"  {threshold_msg}")

    result = {
        'rho': round(rho, 6),
        'p_value': round(p, 6),
        'n': n,
        'ci_95_lower': round(ci_lo, 4),
        'ci_95_upper': round(ci_hi, 4),
        'track': 'placental_20way',
        'option': 'A',
        'window_kb': 2,
        'conclusion': conclusion,
        'threshold_message': threshold_msg
    }

    with open(PRIMARY_JSON, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"  Saved: {PRIMARY_JSON}")

    return result, merged_clean


# ============================================================
# Step 2g: Sensitivity analysis
# ============================================================

def run_sensitivity_analysis(scores_df, rigidity_df, merged_primary):
    """Run full 12-combination sensitivity grid + partial correlation
    + rank-based check.
    """
    print("\n=== Step 2g: Sensitivity analysis ===")

    tracks = scores_df['track'].unique()
    options = scores_df['option'].unique()
    windows = sorted(scores_df['window_kb'].unique())

    results = []
    for track in tracks:
        for opt in options:
            for wkb in windows:
                subset = scores_df[
                    (scores_df['track'] == track) &
                    (scores_df['option'] == opt) &
                    (scores_df['window_kb'] == wkb)
                ].copy()

                merged = subset.merge(rigidity_df, on='cell_type', how='inner')
                merged_clean = merged.dropna(subset=['mean_phastCons'])
                n = len(merged_clean)

                if n < 5:
                    results.append({
                        'track': track, 'option': opt, 'window_kb': wkb,
                        'rho': np.nan, 'p': np.nan, 'n': n,
                        'conclusion': 'INSUFFICIENT_N'
                    })
                    continue

                rho, p = stats.spearmanr(merged_clean['mean_phastCons'],
                                          merged_clean['rigidity_score'])

                if rho >= RHO_POSITIVE and p < 0.05:
                    conc = "POSITIVE"
                elif rho >= RHO_NULL:
                    conc = "TREND"
                else:
                    conc = "NULL_TRIGGERED"

                results.append({
                    'track': track, 'option': opt, 'window_kb': wkb,
                    'rho': round(rho, 4), 'p': round(p, 6), 'n': n,
                    'conclusion': conc
                })

    sens_df = pd.DataFrame(results)
    sens_df.to_csv(SENSITIVITY_CSV, index=False)
    print(f"  Sensitivity table ({len(sens_df)} rows): {SENSITIVITY_CSV}")
    print(sens_df.to_string(index=False))

    # Check for contradictions with primary
    primary_row = sens_df[
        (sens_df['track'] == 'placental_20way') &
        (sens_df['option'] == 'A') &
        (sens_df['window_kb'] == 2)
    ]
    if len(primary_row) > 0:
        primary_conc = primary_row.iloc[0]['conclusion']
        contradictions = sens_df[sens_df['conclusion'] != primary_conc]
        if len(contradictions) > 0:
            print(f"\n  WARNING: {len(contradictions)} sensitivity results "
                  f"contradict primary conclusion ({primary_conc}):")
            for _, r in contradictions.iterrows():
                print(f"    {r['track']} / Option {r['option']} / {r['window_kb']}kb: "
                      f"ρ={r['rho']}, conclusion={r['conclusion']}")

    # Partial correlation controlling for expression level
    print("\n  Partial correlation (controlling for mean expression):")
    if merged_primary is not None and len(merged_primary) > 5:
        # Load centroids to compute mean expression per cell type
        centroids = pd.read_csv(CENTROIDS_HUMAN, index_col=0)
        mean_expr = centroids.mean(axis=1)  # mean across all genes per cell type
        mean_expr.name = 'mean_expression'

        partial_df = merged_primary.set_index('cell_type').join(
            mean_expr, how='left'
        ).dropna(subset=['mean_expression', 'mean_phastCons'])

        if len(partial_df) >= 5:
            # Partial Spearman via residualization
            from scipy.stats import rankdata
            x = rankdata(partial_df['mean_phastCons'])
            y = rankdata(partial_df['rigidity_score'])
            z = rankdata(partial_df['mean_expression'])

            # Residualize x and y on z
            from numpy.polynomial import polynomial as P
            # Simple linear residualization
            slope_xz = np.polyfit(z, x, 1)
            slope_yz = np.polyfit(z, y, 1)
            x_resid = x - np.polyval(slope_xz, z)
            y_resid = y - np.polyval(slope_yz, z)

            rho_partial = np.corrcoef(x_resid, y_resid)[0, 1]
            # p-value approximation for partial correlation
            n_partial = len(partial_df)
            t_stat = rho_partial * np.sqrt((n_partial - 3) / (1 - rho_partial**2))
            p_partial = 2 * stats.t.sf(np.abs(t_stat), df=n_partial - 3)

            print(f"    Partial ρ (controlling for expression) = {rho_partial:.4f}")
            print(f"    p = {p_partial:.6f}, n = {n_partial}")

    # Rank-based vs continuous
    print("\n  Rank-based correlation check:")
    if merged_primary is not None:
        rho_rank, p_rank = stats.spearmanr(
            merged_primary['mean_phastCons'],
            merged_primary['rigidity_rank']
        )
        # Note: rigidity_rank is inverted (1=most rigid), so sign flips
        print(f"    ρ(phastCons, rigidity_rank) = {rho_rank:.4f}, p = {p_rank:.6f}")
        print(f"    (Negative = higher conservation → lower rank number → more rigid)")

    return sens_df


# ============================================================
# Step 2h: Visualization
# ============================================================

def make_plots(merged_primary, sens_df, primary_result):
    """Generate scatter plot (primary) + sensitivity heatmaps."""
    print("\n=== Step 2h: Generating plots ===")

    if merged_primary is None or primary_result is None:
        print("  Skipping plots: no primary result")
        return

    # --- Plot 1: Primary scatter ---
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))

    is_clean = merged_primary['cell_type'].isin(CLEAN_T3C)
    colors = ['#e74c3c' if c else '#3498db' for c in is_clean]
    sizes = [80 if c else 50 for c in is_clean]

    ax.scatter(merged_primary['mean_phastCons'],
               merged_primary['rigidity_score'],
               c=colors, s=sizes, alpha=0.7, edgecolors='black', linewidth=0.5)

    # Label all points
    for _, row in merged_primary.iterrows():
        label = row['cell_type']
        if len(label) > 25:
            label = label[:22] + '...'
        fontsize = 6.5
        ax.annotate(label,
                    (row['mean_phastCons'], row['rigidity_score']),
                    fontsize=fontsize, alpha=0.8,
                    xytext=(5, 3), textcoords='offset points')

    rho = primary_result['rho']
    p = primary_result['p_value']
    n = primary_result['n']
    conc = primary_result['conclusion']

    ax.set_xlabel('Mean phastCons at promoters\n(placental 20way, Option A, ±2kb)',
                  fontsize=11)
    ax.set_ylabel('Procrustes rigidity score\n(higher = more rigid)', fontsize=11)
    ax.set_title(f'T3-E: Regulatory Sequence Conservation vs Procrustes Rigidity\n'
                 f'Spearman ρ = {rho:.3f}, p = {p:.4f}, n = {n} | {conc}',
                 fontsize=12)

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#e74c3c',
               markersize=10, label='T3-C clean types (5)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#3498db',
               markersize=8, label='Other cell types'),
    ]
    ax.legend(handles=legend_elements, loc='best', fontsize=9)

    plt.tight_layout()
    scatter_path = os.path.join(OUTPUT_DIR, 'scatter_primary.png')
    plt.savefig(scatter_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {scatter_path}")

    # --- Plot 2 & 3: Sensitivity heatmaps ---
    for track_filter, fname in [
        ('placental_20way', 'sensitivity_heatmap.png'),
        ('100way_vertebrate', 'sensitivity_heatmap_100way.png')
    ]:
        track_data = sens_df[sens_df['track'] == track_filter]
        if len(track_data) == 0:
            print(f"  Skipping {fname}: no data for {track_filter}")
            continue

        # Pivot: rows=option, cols=window_kb
        pivot = track_data.pivot_table(index='option', columns='window_kb',
                                        values='rho')

        fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        vmax = max(abs(pivot.min().min()), abs(pivot.max().max()), 0.5)
        sns.heatmap(pivot, annot=True, fmt='.3f', cmap='RdBu_r',
                    center=0, vmin=-vmax, vmax=vmax,
                    linewidths=0.5, ax=ax)
        ax.set_title(f'Sensitivity: Spearman ρ by Option × Window\n'
                     f'Track: {track_filter}', fontsize=11)
        ax.set_xlabel('Window size (kb)', fontsize=10)
        ax.set_ylabel('Scoring option', fontsize=10)
        ax.set_yticklabels(['A (loading genes)', 'B (expr-weighted)'],
                           rotation=0, fontsize=9)

        plt.tight_layout()
        heatmap_path = os.path.join(OUTPUT_DIR, fname)
        plt.savefig(heatmap_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {heatmap_path}")


# ============================================================
# Summary report
# ============================================================

def write_summary(primary_result, sens_df, scores_df):
    """Write the T3-E Step 2 summary markdown report."""
    print("\n=== Writing summary report ===")

    summary_path = os.path.join(OUTPUT_DIR, 't3e_step2_summary.md')

    # Gather info
    rho = primary_result['rho'] if primary_result else None
    p = primary_result['p_value'] if primary_result else None
    n = primary_result['n'] if primary_result else None
    conc = primary_result['conclusion'] if primary_result else 'NOT_RUN'

    lines = [
        "# T3-E Step 2: Regulatory Sequence Conservation vs Procrustes Rigidity",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M')}",
        "",
        "## 1. UCSC Track Used",
        "",
        "**Primary track:** phastCons20way (20 species: 17 primates + treeshrew + mouse + dog)",
        "- Includes Mus musculus (mm10) — confirmed",
        "- Source: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons20way/",
        "- Note: primate-dominated alignment. Conservation scores primarily reflect",
        "  primate constraint with mouse as an outgroup anchor.",
        "",
        "**Sensitivity track:** phastCons100way (100 vertebrates)",
        "- Includes Mus musculus (mm10), rat, zebrafish, and 97 other species",
        "- Broadest taxonomic scope for calibrating deep conservation",
        "",
    ]

    if primary_result:
        lines.extend([
            "## 2. Primary Spearman Result",
            "",
            f"| Metric | Value |",
            f"|--------|-------|",
            f"| Spearman ρ | {rho:.4f} |",
            f"| p-value | {p:.6f} |",
            f"| n (cell types) | {n} |",
            f"| 95% CI | [{primary_result['ci_95_lower']:.4f}, {primary_result['ci_95_upper']:.4f}] |",
            f"| Track | placental_20way |",
            f"| Option | A (cell-type-specific top-50 loading genes) |",
            f"| Window | ±2kb |",
            "",
            f"**Pre-registered threshold triggered:** {conc}",
            "",
        ])

        if conc == 'POSITIVE':
            lines.append("ρ ≥ 0.50 AND p < 0.05 → **POSITIVE** — regulatory sequence "
                         "conservation at identity-gene promoters predicts Procrustes rigidity.")
        elif conc == 'TREND':
            lines.append(f"ρ = {rho:.3f} in [0.35, 0.50) → **TREND** — directionally "
                         "consistent but below the positive threshold. Underpowered or partial signal.")
        else:
            lines.append(f"ρ = {rho:.3f} < 0.35 → **8TH NULL TRIGGERED** — close chromatin "
                         "accessibility/regulatory sequence as proximate mechanism for rigidity.")
        lines.append("")

    if sens_df is not None and len(sens_df) > 0:
        lines.extend([
            "## 3. Sensitivity Analysis",
            "",
            "| Track | Option | Window (kb) | ρ | p | Conclusion |",
            "|-------|--------|-------------|-----|-----|-----------|",
        ])
        for _, r in sens_df.iterrows():
            rho_s = f"{r['rho']:.4f}" if pd.notna(r['rho']) else "N/A"
            p_s = f"{r['p']:.6f}" if pd.notna(r['p']) else "N/A"
            lines.append(f"| {r['track']} | {r['option']} | {r['window_kb']} | "
                         f"{rho_s} | {p_s} | {r['conclusion']} |")
        lines.append("")

        if primary_result:
            primary_conc = conc
            contradictions = sens_df[sens_df['conclusion'] != primary_conc]
            if len(contradictions) > 0:
                lines.extend([
                    f"**{len(contradictions)} sensitivity results contradict "
                    f"primary conclusion ({primary_conc}).** These require explicit discussion.",
                    ""
                ])
            else:
                lines.extend([
                    "**All sensitivity results are consistent with the primary conclusion.**",
                    ""
                ])

    if scores_df is not None:
        option_a = scores_df[scores_df['option'] == 'A']
        low_cov = option_a[option_a['n_genes_found'] < 40]
        lines.extend([
            "## 4. Data Quality Flags",
            "",
            f"- Total conservation scores computed: {len(scores_df)}",
            f"- Option A entries with <80% gene coverage (n<40): {len(low_cov)}",
        ])
        if len(low_cov) > 0:
            for _, r in low_cov.iterrows():
                lines.append(f"  - {r['cell_type']} ({r['track']}, {r['window_kb']}kb): "
                             f"n={r['n_genes_found']}")
        lines.append("")

    lines.extend([
        "## 5. Partial Correlation",
        "",
        "See console output for partial correlation controlling for mean expression level.",
        "If the partial correlation changes the conclusion, it means the phastCons signal",
        "is confounded with expression level — the same genes that are highly expressed",
        "tend to be more conserved (known phenomenon).",
        "",
        "## Files Generated",
        "",
        f"- `{os.path.basename(SCORES_CSV)}` — All conservation scores (420 rows target)",
        f"- `{os.path.basename(MERGED_CSV)}` — Primary analysis merged data (35 rows)",
        f"- `{os.path.basename(PRIMARY_JSON)}` — Primary Spearman result",
        f"- `{os.path.basename(SENSITIVITY_CSV)}` — Sensitivity grid (12 rows)",
        f"- `{os.path.basename(TSS_BED)}` — TSS coordinates (BED6)",
        "- `scatter_primary.png` — Primary result scatter plot",
        "- `sensitivity_heatmap.png` — Sensitivity heatmap (placental track)",
        "- `sensitivity_heatmap_100way.png` — Sensitivity heatmap (100way track)",
    ])

    with open(summary_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print(f"  Saved: {summary_path}")


# ============================================================
# Main
# ============================================================

def main():
    print("=" * 70)
    print("T3-E Step 2: Regulatory Sequence Conservation vs Procrustes Rigidity")
    print("=" * 70)

    # Step 2b: Load gene sets
    option_a_genes, option_b_genes, gene_id_to_name, centroids_df, all_genes = \
        load_gene_sets()

    # Step 2b: Get TSS coordinates
    tss_data = fetch_tss_coordinates(all_genes, gene_id_to_name)

    # Step 2d: Compute conservation scores
    scores_df = compute_conservation_scores(
        tss_data, option_a_genes, option_b_genes,
        centroids_df, gene_id_to_name
    )

    # Step 2e: Load rigidity
    rigidity_df = load_rigidity_scores()

    # Step 2f: Primary Spearman
    primary_result, merged_primary = run_primary_spearman(scores_df, rigidity_df)

    # Step 2g: Sensitivity
    sens_df = run_sensitivity_analysis(scores_df, rigidity_df, merged_primary)

    # Step 2h: Visualization
    make_plots(merged_primary, sens_df, primary_result)

    # Summary report
    write_summary(primary_result, sens_df, scores_df)

    print("\n" + "=" * 70)
    print("T3-E Step 2 COMPLETE")
    print("=" * 70)


if __name__ == '__main__':
    main()
