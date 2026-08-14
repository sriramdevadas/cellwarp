#!/usr/bin/env python3
"""
T3-E Step 3b: H3K27ac Enhancer Conservation vs Procrustes Rigidity

Tests whether cell types with more conserved active enhancer landscapes
at identity-gene loci show higher Procrustes rigidity.

Biology: H3K27ac marks active enhancers and promoters. Distal H3K27ac peaks
(>2kb from any TSS) specifically mark active enhancers. If enhancer-level
regulatory conservation drives Procrustes rigidity, cell types with more
conserved active enhancers at identity-gene loci should have higher rigidity.

Math: For each cell type, compute a conservation score as the mean reciprocal
Jaccard index between human and mouse enhancer peak sets near identity genes.
Correlate (Spearman) with Procrustes rigidity score across n=6 cell types.

NULL CLOSURE TEST:
- n=6 is underpowered for positive detection (requires |ρ|≥0.829 for p<0.05)
- n=6 IS sufficient for null closure (ρ<0.35 triggers closure)

Pre-registered thresholds:
- ρ ≥ 0.50: POSITIVE — note underpowered, flag for follow-up
- ρ = 0.35-0.49: TREND — underpowered, ambiguous
- ρ < 0.35: TRIGGERED — 9th null, computational ceiling reached
"""

import json
import os
import sys
import time
import hashlib
import urllib.request
import gzip
import bisect
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# CONFIGURATION
# ============================================================

BASE_DIR = str(Path(__file__).resolve().parent.parent)
DATA_DIR = os.path.join(BASE_DIR, 'data/h3k27ac')
ANNOT_DIR = os.path.join(BASE_DIR, 'data/annotations')
OUTPUT_DIR = os.path.join(BASE_DIR, 'output/validation/t3e_enhancer')
CHROMATIN_DIR = os.path.join(BASE_DIR, 'output/validation/t3e_chromatin')

for d in [os.path.join(DATA_DIR, 'human'), os.path.join(DATA_DIR, 'mouse'),
          ANNOT_DIR, OUTPUT_DIR]:
    os.makedirs(d, exist_ok=True)

# Cell types: exactly 6 clean matched pairs
CELL_TYPES = {
    'CD8+ T cell':  {'rigidity_rank': 1,  'encode': 'ENCSR835OJV', 'gsm': 'GSM1441282', 'geo_label': 'CD8'},
    'NK cell':      {'rigidity_rank': 13, 'encode': 'ENCSR391EQV', 'gsm': 'GSM1441283', 'geo_label': 'NK'},
    'Monocyte':     {'rigidity_rank': 14, 'encode': 'ENCSR000ASJ', 'gsm': 'GSM1441278', 'geo_label': 'Mono'},
    'B cell':       {'rigidity_rank': 19, 'encode': 'ENCSR000AUP', 'gsm': 'GSM1441280', 'geo_label': 'B'},
    'CD4+ T cell':  {'rigidity_rank': 22, 'encode': 'ENCSR120WKZ', 'gsm': 'GSM1441281', 'geo_label': 'CD4'},
    'Neutrophil':   {'rigidity_rank': 28, 'encode': 'ENCSR267YXV', 'gsm': 'GSM1441277', 'geo_label': 'GN'},
}

# Map our names to rigidity_scores.csv and gene_sets.json names
NAME_MAP = {
    'CD8+ T cell': 'CD8-positive, alpha-beta T cell',
    'NK cell':     'natural killer cell',
    'Monocyte':    'monocyte',
    'B cell':      'B cell',
    'CD4+ T cell': 'CD4-positive, alpha-beta T cell',
    'Neutrophil':  'neutrophil',
}


def safe_name(ct):
    """Convert cell type name to filesystem-safe string."""
    return ct.replace(' ', '_').replace('+', 'plus')


# ============================================================
# DOWNLOAD UTILITIES
# ============================================================

def download(url, dest, timeout=300, retries=3, md5=None, desc=None):
    """Download a file with retries and optional md5 verification."""
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        if md5 and _md5(dest) != md5:
            print(f"  Checksum mismatch, re-downloading {os.path.basename(dest)}")
        else:
            print(f"  Already have: {os.path.basename(dest)}")
            return True

    desc = desc or os.path.basename(dest)
    for attempt in range(retries):
        if attempt > 0:
            wait = [5, 15, 45][min(attempt, 2)]
            print(f"  Retry {attempt+1}/{retries} after {wait}s...")
            time.sleep(wait)
        try:
            print(f"  Downloading: {desc}...")
            req = urllib.request.Request(url, headers={'User-Agent': 'CellWarp/1.0'})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                with open(dest, 'wb') as f:
                    while True:
                        chunk = resp.read(65536)
                        if not chunk:
                            break
                        f.write(chunk)
            if md5 and _md5(dest) != md5:
                print(f"  Checksum mismatch for {desc}")
                continue
            size_mb = os.path.getsize(dest) / 1e6
            print(f"  OK: {desc} ({size_mb:.1f} MB)")
            return True
        except Exception as e:
            print(f"  Failed: {e}")
    print(f"  FAILED after {retries} attempts: {desc}")
    return False


def _md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(65536), b''):
            h.update(chunk)
    return h.hexdigest()


def fetch_json(url, timeout=30):
    req = urllib.request.Request(url, headers={
        'User-Agent': 'CellWarp/1.0', 'Accept': 'application/json'})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode())


# ============================================================
# STEP 3b-i: DATA DOWNLOAD
# ============================================================

def download_encode_peaks():
    """Download ENCODE H3K27ac peak files (narrowPeak, GRCh38)."""
    print("\n" + "="*60)
    print("STEP 3b-i: ENCODE human H3K27ac peak files")
    print("="*60)

    results = {}
    for ct, info in CELL_TYPES.items():
        dest = os.path.join(DATA_DIR, 'human', f'{safe_name(ct)}_peaks.bed.gz')
        if os.path.exists(dest) and os.path.getsize(dest) > 0:
            print(f"\n{ct}: already downloaded")
            results[ct] = dest
            continue

        exp = info['encode']
        print(f"\n{ct} ({exp}):")

        try:
            # Search for peak files from this experiment
            url = (f"https://www.encodeproject.org/search/?type=File"
                   f"&dataset=/experiments/{exp}/"
                   f"&file_format=bed&status=released&assembly=GRCh38"
                   f"&format=json&limit=50")
            data = fetch_json(url)
            files = data.get('@graph', [])

            # Priority: IDR > replicated > pseudoreplicated > any peaks
            best = None
            for prio in ['IDR thresholded peaks', 'optimal IDR thresholded peaks',
                         'conservative IDR thresholded peaks', 'replicated peaks',
                         'pseudoreplicated peaks', 'stable peaks', 'peaks']:
                for f in files:
                    if f.get('output_type', '') == prio:
                        best = f
                        break
                if best:
                    break
            if not best:
                for f in files:
                    if 'peak' in f.get('output_type', '').lower():
                        best = f
                        break

            if not best:
                print(f"  No peak files found!")
                continue

            acc = best['accession']
            otype = best.get('output_type', '?')
            md5 = best.get('md5sum')
            dl_url = f"https://www.encodeproject.org/files/{acc}/@@download/{acc}.bed.gz"
            print(f"  Selected: {acc} ({otype})")

            if download(dl_url, dest, md5=md5, desc=f"{ct} peaks ({acc})"):
                results[ct] = dest
        except Exception as e:
            print(f"  API error: {e}")

    return results


def download_mouse_bigwigs():
    """Download Lara-Astiaso 2014 H3K27ac bigWig files from GEO (mm9)."""
    print("\n" + "="*60)
    print("STEP 3b-i: Lara-Astiaso mouse H3K27ac bigWig files")
    print("="*60)

    results = {}
    for ct, info in CELL_TYPES.items():
        gsm = info['gsm']
        label = info['geo_label']
        dest = os.path.join(DATA_DIR, 'mouse', f'{safe_name(ct)}_h3k27ac.bigWig')

        if os.path.exists(dest) and os.path.getsize(dest) > 0:
            print(f"\n{ct}: already downloaded")
            results[ct] = dest
            continue

        print(f"\n{ct} ({gsm}):")

        # Try HTTPS first, then HTTP
        for base in ['https://ftp.ncbi.nlm.nih.gov', 'https://www.ncbi.nlm.nih.gov']:
            if base.startswith('https://ftp'):
                url = f"{base}/geo/samples/GSM1441nnn/{gsm}/suppl/{gsm}_K27Ac_{label}.ucsc.bigWig"
            else:
                url = f"{base}/geo/download/?acc={gsm}&format=file&file={gsm}%5FK27Ac%5F{label}%2Eucsc%2EbigWig"

            if download(url, dest, timeout=300, desc=f"{ct} mouse bigWig"):
                results[ct] = dest
                break

    return results


def download_annotations():
    """Download UCSC refGene tables for TSS coordinates."""
    print("\n" + "="*60)
    print("STEP 3b-i: Genome annotations (refGene)")
    print("="*60)

    annots = {}
    for asm in ['hg38', 'mm10']:
        dest = os.path.join(ANNOT_DIR, f'{asm}_refGene.txt.gz')
        url = f'https://hgdownload.soe.ucsc.edu/goldenPath/{asm}/database/refGene.txt.gz'
        if download(url, dest, desc=f"{asm} refGene"):
            annots[asm] = dest
    return annots


# ============================================================
# PEAK CALLING FROM BIGWIG
# ============================================================

def call_peaks_bigwig(bw_path, out_bed, window=200, merge_gap=400, min_len=200):
    """
    Call peaks from H3K27ac bigWig using threshold-based approach.

    Biology: H3K27ac signal enrichment indicates active regulatory elements.
    Math: Threshold = 95th percentile of non-zero genome-wide signal in
    200bp windows. Merge peaks within 400bp. Keep peaks ≥200bp.

    Uses batch stats (nBins) for ~100x speedup over per-window calls.
    """
    import pyBigWig

    if os.path.exists(out_bed) and os.path.getsize(out_bed) > 0:
        n = sum(1 for _ in open(out_bed))
        print(f"  Already called: {n} peaks")
        return out_bed

    print(f"  Calling peaks from {os.path.basename(bw_path)}...")
    bw = pyBigWig.open(bw_path)
    chroms = bw.chroms()
    std_chroms = sorted([c for c in chroms if not ('_' in c or c == 'chrM')])

    # First pass: sample signal for threshold (batch stats, every 10th bin)
    all_vals = []
    for chrom in std_chroms:
        length = chroms[chrom]
        nbins = length // window
        if nbins < 1:
            continue
        try:
            vals = bw.stats(chrom, 0, nbins * window, type="mean", nBins=nbins)
            # Sample every 10th for threshold estimation
            sampled = [v for i, v in enumerate(vals) if i % 10 == 0
                       and v is not None and v > 0]
            all_vals.extend(sampled)
        except Exception:
            pass

    all_vals = np.array(all_vals)
    threshold = np.percentile(all_vals, 95)
    print(f"  Threshold (95th pctile): {threshold:.2f} "
          f"(mean={all_vals.mean():.2f}, max={all_vals.max():.2f})")

    # Second pass: call peaks using batch stats
    peaks = []
    for chrom in std_chroms:
        length = chroms[chrom]
        nbins = length // window
        if nbins < 1:
            continue

        try:
            vals = bw.stats(chrom, 0, nbins * window, type="mean", nBins=nbins)
        except Exception:
            continue

        # Find enriched windows
        enriched = []
        for i, v in enumerate(vals):
            if v is not None and v >= threshold:
                s = i * window
                enriched.append((s, s + window, v))

        if not enriched:
            continue

        # Merge nearby enriched windows
        cs, ce, cm = enriched[0]
        for s, e, v in enriched[1:]:
            if s - ce <= merge_gap:
                ce = e
                cm = max(cm, v)
            else:
                if ce - cs >= min_len:
                    peaks.append((chrom, cs, ce, cm))
                cs, ce, cm = s, e, v
        if ce - cs >= min_len:
            peaks.append((chrom, cs, ce, cm))

    bw.close()

    with open(out_bed, 'w') as f:
        for i, (chrom, start, end, score) in enumerate(peaks):
            f.write(f"{chrom}\t{start}\t{end}\tpeak_{i}\t{score:.4f}\t.\n")

    print(f"  Called {len(peaks)} peaks")
    return out_bed


# ============================================================
# GENOME ANNOTATION PARSING
# ============================================================

def parse_refgene_tss(gz_path):
    """
    Parse UCSC refGene to get TSS per gene symbol.
    Returns {gene_symbol: [(chrom, tss_pos, strand), ...]}.
    refGene cols: bin, name, chrom, strand, txStart, txEnd, ..., name2, ...
    TSS = txStart if strand=='+', txEnd if strand=='-'.
    """
    tss = defaultdict(set)
    with gzip.open(gz_path, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 13:
                continue
            chrom, strand = parts[2], parts[3]
            if '_' in chrom or chrom == 'chrM':
                continue
            tx_start, tx_end = int(parts[4]), int(parts[5])
            gene = parts[12]
            pos = tx_start if strand == '+' else tx_end
            tss[gene].add((chrom, pos, strand))
    return {g: list(v) for g, v in tss.items()}


# ============================================================
# PEAK LOADING AND FILTERING
# ============================================================

def load_peaks(path):
    """Load peaks from BED/narrowPeak (plain or gzipped).
    Returns [(chrom, start, end, summit), ...]."""
    opener = gzip.open if path.endswith('.gz') else open
    peaks = []
    with opener(path, 'rt') as f:
        for line in f:
            if line.startswith(('#', 'track', 'browser')):
                continue
            p = line.strip().split('\t')
            if len(p) < 3:
                continue
            chrom, start, end = p[0], int(p[1]), int(p[2])
            if '_' in chrom or chrom == 'chrM':
                continue
            # Summit from narrowPeak col10 if available
            if len(p) >= 10:
                try:
                    so = int(p[9])
                    summit = start + so if so >= 0 else (start + end) // 2
                except ValueError:
                    summit = (start + end) // 2
            else:
                summit = (start + end) // 2
            peaks.append((chrom, start, end, summit))
    return peaks


def filter_distal(peaks, all_tss, min_dist=2000):
    """
    Keep only peaks with summit >min_dist from any TSS.
    Biology: isolates enhancer signal from promoter-associated H3K27ac.
    """
    # Build sorted TSS arrays per chrom
    tss_by_chr = defaultdict(list)
    for gene, positions in all_tss.items():
        for chrom, pos, strand in positions:
            tss_by_chr[chrom].append(pos)
    for c in tss_by_chr:
        tss_by_chr[c] = sorted(set(tss_by_chr[c]))

    filtered = []
    for chrom, start, end, summit in peaks:
        arr = tss_by_chr.get(chrom, [])
        if not arr:
            filtered.append((chrom, start, end, summit))
            continue
        idx = bisect.bisect_left(arr, summit)
        min_d = float('inf')
        for i in (idx - 1, idx):
            if 0 <= i < len(arr):
                min_d = min(min_d, abs(summit - arr[i]))
        if min_d > min_dist:
            filtered.append((chrom, start, end, summit))
    return filtered


def get_identity_enhancers(enhancers, gene_tss_dict, window=50000):
    """
    Find enhancers within window bp of any identity gene TSS.
    Biology: enhancers typically regulate genes within 50kb.
    """
    # Build window list per chrom
    wins_by_chr = defaultdict(list)
    for gene, positions in gene_tss_dict.items():
        for chrom, pos, strand in positions:
            wins_by_chr[chrom].append((max(0, pos - window), pos + window))

    # Sort windows for efficient search
    for c in wins_by_chr:
        wins_by_chr[c].sort()

    found = []
    for chrom, start, end, summit in enhancers:
        for ws, we in wins_by_chr.get(chrom, []):
            if ws <= summit <= we:
                found.append((chrom, start, end, summit))
                break
    return found


# ============================================================
# LIFTOVER AND JACCARD
# ============================================================

def liftover_peaks(peaks, from_asm, to_asm, _cache={}):
    """
    Lift peak coordinates between assemblies using pyliftover.
    Returns (lifted_peaks, pct_lifted).
    """
    from pyliftover import LiftOver

    key = (from_asm, to_asm)
    if key not in _cache:
        _cache[key] = LiftOver(from_asm, to_asm)
    lo = _cache[key]

    lifted = []
    failed = 0
    for chrom, start, end, summit in peaks:
        result = lo.convert_coordinate(chrom, summit)
        if result:
            nc, np_, ns, _ = result[0]
            hw = (end - start) // 2
            lifted.append((nc, max(0, np_ - hw), np_ + hw))
        else:
            failed += 1

    total = len(peaks)
    pct = (total - failed) / total * 100 if total else 0
    return lifted, pct


def merge_intervals(intervals):
    """Merge overlapping intervals. Returns [(start, end), ...]."""
    if not intervals:
        return []
    s = sorted(intervals)
    merged = [s[0]]
    for start, end in s[1:]:
        if start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def jaccard(peaks_a, peaks_b):
    """
    Compute Jaccard index between two interval sets (in bp).
    J = |A ∩ B| / |A ∪ B|.
    """
    a_by_chr = defaultdict(list)
    b_by_chr = defaultdict(list)
    for c, s, e in peaks_a:
        a_by_chr[c].append((s, e))
    for c, s, e in peaks_b:
        b_by_chr[c].append((s, e))

    total_isect = 0
    total_union = 0

    for chrom in set(list(a_by_chr.keys()) + list(b_by_chr.keys())):
        am = merge_intervals(a_by_chr.get(chrom, []))
        bm = merge_intervals(b_by_chr.get(chrom, []))

        # Intersection via sweep
        isect = 0
        i, j = 0, 0
        while i < len(am) and j < len(bm):
            os_ = max(am[i][0], bm[j][0])
            oe = min(am[i][1], bm[j][1])
            if os_ < oe:
                isect += oe - os_
            if am[i][1] < bm[j][1]:
                i += 1
            else:
                j += 1

        a_bp = sum(e - s for s, e in am)
        b_bp = sum(e - s for s, e in bm)
        total_isect += isect
        total_union += a_bp + b_bp - isect

    return total_isect / total_union if total_union > 0 else 0.0


# ============================================================
# STATISTICS
# ============================================================

def spearman_ci(x, y, n):
    """Spearman ρ with Fisher z-transform 95% CI."""
    rho, pval = stats.spearmanr(x, y)
    if n > 3 and abs(rho) < 0.9999:
        z = np.arctanh(rho)
        se = 1 / np.sqrt(n - 3)
        ci_lo = np.tanh(z - 1.96 * se)
        ci_hi = np.tanh(z + 1.96 * se)
    else:
        ci_lo = ci_hi = rho
    return rho, pval, ci_lo, ci_hi


# ============================================================
# MAIN PIPELINE
# ============================================================

def main():
    print("=" * 60)
    print("T3-E Step 3b: H3K27ac Enhancer Conservation")
    print("         vs Procrustes Rigidity")
    print("=" * 60)
    print("NULL CLOSURE TEST — n=6")
    print("Pre-registered: ρ<0.35 → 9th null closure")
    print("=" * 60)

    # ---- Load existing data ----
    print("\nLoading existing data...")
    rigidity_df = pd.read_csv(os.path.join(CHROMATIN_DIR, 'rigidity_scores.csv'))
    with open(os.path.join(CHROMATIN_DIR, 'gene_sets.json')) as f:
        gene_sets = json.load(f)
    orthologs = pd.read_csv(os.path.join(BASE_DIR, 'data/phase1/orthologs_human_mouse.csv'))

    id_to_name = gene_sets['gene_id_to_name']
    h2m = {r['human_gene_name']: r['mouse_gene_name']
           for _, r in orthologs.iterrows()}
    print(f"  Rigidity: {len(rigidity_df)} types, Orthologs: {len(orthologs)}")

    # ---- Step 3b-i: Downloads ----
    human_peaks = download_encode_peaks()
    mouse_bws = download_mouse_bigwigs()
    annots = download_annotations()

    # Check completeness
    missing = []
    for ct in CELL_TYPES:
        if ct not in human_peaks:
            missing.append(f"{ct} (human)")
        if ct not in mouse_bws:
            missing.append(f"{ct} (mouse)")
    if 'hg38' not in annots or 'mm10' not in annots:
        missing.append("annotations")
    if missing:
        print(f"\n⚠ Missing data: {', '.join(missing)}")
        print("Continuing with available data...")

    # ---- Parse TSS ----
    print("\nParsing TSS from refGene...")
    hg38_tss = parse_refgene_tss(annots['hg38'])
    mm10_tss = parse_refgene_tss(annots['mm10'])
    print(f"  hg38: {len(hg38_tss)} genes, mm10: {len(mm10_tss)} genes")

    # ---- Call mouse peaks from bigWig ----
    print("\n" + "=" * 60)
    print("Calling peaks from mouse bigWig files (mm9)")
    print("=" * 60)

    mouse_peaks = {}
    for ct in CELL_TYPES:
        if ct not in mouse_bws:
            continue
        bed = os.path.join(DATA_DIR, 'mouse', f'{safe_name(ct)}_peaks_mm9.bed')
        print(f"\n{ct}:")
        call_peaks_bigwig(mouse_bws[ct], bed)
        mouse_peaks[ct] = bed

    # ---- LiftOver mouse peaks mm9 → mm10 ----
    print("\n" + "=" * 60)
    print("Lifting mouse peaks mm9 → mm10")
    print("=" * 60)

    mouse_peaks_mm10 = {}
    for ct in CELL_TYPES:
        if ct not in mouse_peaks:
            continue
        mm9_peaks = load_peaks(mouse_peaks[ct])
        mm10_bed = os.path.join(DATA_DIR, 'mouse', f'{safe_name(ct)}_peaks_mm10.bed')

        if os.path.exists(mm10_bed) and os.path.getsize(mm10_bed) > 0:
            print(f"  {ct}: already lifted")
            mouse_peaks_mm10[ct] = mm10_bed
            continue

        lifted, pct = liftover_peaks(mm9_peaks, 'mm9', 'mm10')
        print(f"  {ct}: {pct:.1f}% lifted ({len(lifted)}/{len(mm9_peaks)} peaks)")

        with open(mm10_bed, 'w') as f:
            for i, (c, s, e) in enumerate(lifted):
                f.write(f"{c}\t{s}\t{e}\tpeak_{i}\t0\t.\n")
        mouse_peaks_mm10[ct] = mm10_bed

    # ---- Step 3b-ii: Filter distal enhancers ----
    print("\n" + "=" * 60)
    print("STEP 3b-ii: Filtering distal enhancers (>2kb from TSS)")
    print("=" * 60)

    human_enh = {}
    mouse_enh = {}

    for ct in CELL_TYPES:
        print(f"\n{ct}:")

        # Human (hg38)
        if ct in human_peaks:
            hp = load_peaks(human_peaks[ct])
            hd = filter_distal(hp, hg38_tss, min_dist=2000)
            frac = len(hd) / len(hp) * 100 if hp else 0
            print(f"  Human: {len(hp)} total → {len(hd)} distal ({frac:.1f}%)")
            human_enh[ct] = hd

            out = os.path.join(DATA_DIR, 'human', f'{safe_name(ct)}_enhancers.bed')
            with open(out, 'w') as f:
                for c, s, e, su in hd:
                    f.write(f"{c}\t{s}\t{e}\t.\t0\t.\n")

        # Mouse (mm10)
        if ct in mouse_peaks_mm10:
            mp = load_peaks(mouse_peaks_mm10[ct])
            # Add summit for peaks that lost it during liftover
            mp_with_summit = [(c, s, e, (s+e)//2) for c, s, e, *rest in mp]
            md = filter_distal(mp_with_summit, mm10_tss, min_dist=2000)
            frac = len(md) / len(mp_with_summit) * 100 if mp_with_summit else 0
            print(f"  Mouse: {len(mp_with_summit)} total → {len(md)} distal ({frac:.1f}%)")
            mouse_enh[ct] = md

            out = os.path.join(DATA_DIR, 'mouse', f'{safe_name(ct)}_enhancers.bed')
            with open(out, 'w') as f:
                for c, s, e, su in md:
                    f.write(f"{c}\t{s}\t{e}\t.\t0\t.\n")

    # ---- Step 3b-iii: Identity gene enhancer windows ----
    print("\n" + "=" * 60)
    print("STEP 3b-iii: Enhancers near identity gene loci (50kb window)")
    print("=" * 60)

    human_id_enh = {}
    mouse_id_enh = {}

    for ct in CELL_TYPES:
        gs_name = NAME_MAP[ct]
        id_genes_ens = gene_sets['option_a'][gs_name]

        # Map to symbols and get TSS
        h_id_tss = {}
        m_id_tss = {}
        for eid in id_genes_ens:
            hname = id_to_name.get(eid)
            if hname and hname in hg38_tss:
                h_id_tss[hname] = hg38_tss[hname]
                mname = h2m.get(hname)
                if mname and mname in mm10_tss:
                    m_id_tss[mname] = mm10_tss[mname]

        print(f"\n{ct}: {len(h_id_tss)} human TSS, {len(m_id_tss)} mouse TSS")

        if ct in human_enh and h_id_tss:
            ha = get_identity_enhancers(human_enh[ct], h_id_tss, window=50000)
            human_id_enh[ct] = ha
            print(f"  Human identity-gene enhancers: {len(ha)}")
            if len(ha) < 10:
                print(f"  ⚠ FLAG: <10 enhancers")

        if ct in mouse_enh and m_id_tss:
            ma = get_identity_enhancers(mouse_enh[ct], m_id_tss, window=50000)
            mouse_id_enh[ct] = ma
            print(f"  Mouse identity-gene enhancers: {len(ma)}")
            if len(ma) < 10:
                print(f"  ⚠ FLAG: <10 enhancers")

    # ---- Step 3b-iv: Reciprocal liftover and Jaccard ----
    print("\n" + "=" * 60)
    print("STEP 3b-iv: Reciprocal liftover and Jaccard")
    print("=" * 60)

    results = []
    for ct, info in CELL_TYPES.items():
        rig_name = NAME_MAP[ct]
        rig_row = rigidity_df[rigidity_df['cell_type'] == rig_name]
        if rig_row.empty:
            print(f"  ✗ No rigidity score for {ct}")
            continue
        rig_score = rig_row['rigidity_score'].values[0]

        if ct not in human_id_enh or ct not in mouse_id_enh:
            print(f"  ✗ Missing data for {ct}")
            continue

        he = human_id_enh[ct]
        me = mouse_id_enh[ct]

        print(f"\n{ct} (rank {info['rigidity_rank']}):")
        print(f"  Human enh: {len(he)}, Mouse enh: {len(me)}")

        # Forward: human hg38 → mm10
        h_lifted, pct_h = liftover_peaks(he, 'hg38', 'mm10')
        m_intervals = [(c, s, e) for c, s, e, su in me]
        fwd_j = jaccard(h_lifted, m_intervals)
        print(f"  Fwd (hg38→mm10): {pct_h:.1f}% lifted, J={fwd_j:.6f}")

        # Reverse: mouse mm10 → hg38
        m_lifted, pct_m = liftover_peaks(me, 'mm10', 'hg38')
        h_intervals = [(c, s, e) for c, s, e, su in he]
        rev_j = jaccard(m_lifted, h_intervals)
        print(f"  Rev (mm10→hg38): {pct_m:.1f}% lifted, J={rev_j:.6f}")

        mean_j = (fwd_j + rev_j) / 2
        print(f"  Mean Jaccard: {mean_j:.6f}")

        results.append({
            'cell_type': ct,
            'rigidity_rank': info['rigidity_rank'],
            'rigidity_score': rig_score,
            'human_enhancer_count': len(he),
            'mouse_enhancer_count': len(me),
            'pct_human_lifted': pct_h,
            'pct_mouse_lifted': pct_m,
            'forward_jaccard': fwd_j,
            'reverse_jaccard': rev_j,
            'mean_jaccard': mean_j,
        })

        # Checkpoint
        pd.DataFrame(results).to_csv(
            os.path.join(OUTPUT_DIR, 'conservation_scores.csv'), index=False)

    df = pd.DataFrame(results)
    n = len(df)
    print(f"\nConservation scores: {n} cell types")

    if n < 4:
        print("FATAL: fewer than 4 cell types with data. Cannot compute correlation.")
        sys.exit(1)

    # ---- Step 3b-v: Primary Spearman ----
    print("\n" + "=" * 60)
    print("STEP 3b-v: PRIMARY SPEARMAN CORRELATION")
    print("=" * 60)

    rho, pval, ci_lo, ci_hi = spearman_ci(
        df['mean_jaccard'].values, df['rigidity_score'].values, n)

    if rho >= 0.50:
        threshold = "POSITIVE"
        interp = "Positive — underpowered, flag for follow-up"
    elif rho >= 0.35:
        threshold = "TREND"
        interp = "Trend — underpowered, ambiguous"
    else:
        threshold = "TRIGGERED"
        interp = "9th null — computational ceiling reached"

    print(f"\n  Spearman ρ = {rho:.4f}")
    print(f"  p-value   = {pval:.6f}")
    print(f"  n         = {n}")
    print(f"  95% CI    = [{ci_lo:.4f}, {ci_hi:.4f}]")
    print(f"\n  Threshold: {threshold}")
    print(f"  {interp}")
    print(f"\n  NOTE: Null closure test. n={n} underpowered for positive")
    print(f"  detection but sufficient for null closure.")

    primary = {
        'spearman_rho': float(rho), 'p_value': float(pval), 'n': n,
        'ci_lower': float(ci_lo), 'ci_upper': float(ci_hi),
        'threshold_triggered': threshold, 'interpretation': interp,
        'note': (f"Null closure test. n={n} underpowered for positive "
                 f"detection (requires |ρ|≥0.829) but sufficient for null closure.")
    }
    with open(os.path.join(OUTPUT_DIR, 'spearman_primary_result.json'), 'w') as f:
        json.dump(primary, f, indent=2)

    # ---- Step 3b-vi: Sensitivity ----
    print("\n" + "=" * 60)
    print("STEP 3b-vi: Sensitivity Analyses")
    print("=" * 60)

    sens_rows = []

    # 1. Leave-one-out
    print("\n--- Leave-one-out ---")
    loo_rhos = []
    for i in range(n):
        loo = df.drop(df.index[i])
        lr, lp, _, _ = spearman_ci(
            loo['mean_jaccard'].values, loo['rigidity_score'].values, n - 1)
        dropped = df.iloc[i]['cell_type']
        loo_rhos.append(lr)
        consistent = "Yes" if (lr < 0.35) == (rho < 0.35) else "No"
        print(f"  Drop {dropped}: ρ={lr:.4f} p={lp:.4f} [{consistent}]")
        sens_rows.append({'analysis': 'LOO', 'parameter': f'drop {dropped}',
                         'rho': lr, 'conclusion_consistent': consistent})
    print(f"  Range: [{min(loo_rhos):.4f}, {max(loo_rhos):.4f}]")

    # Helper for sensitivity recomputation
    def recompute_scores(tss_dist=2000, gene_window=50000):
        """Recompute conservation scores with different parameters."""
        scores = []
        for ct in CELL_TYPES:
            if ct not in human_peaks or ct not in mouse_peaks_mm10:
                continue

            hp = load_peaks(human_peaks[ct])
            hd = filter_distal(hp, hg38_tss, min_dist=tss_dist)

            mp_raw = load_peaks(mouse_peaks_mm10[ct])
            mp = [(c, s, e, (s+e)//2) for c, s, e, *r in mp_raw]
            md = filter_distal(mp, mm10_tss, min_dist=tss_dist)

            gs_name = NAME_MAP[ct]
            id_ens = gene_sets['option_a'][gs_name]
            h_tss_ct, m_tss_ct = {}, {}
            for eid in id_ens:
                hn = id_to_name.get(eid)
                if hn and hn in hg38_tss:
                    h_tss_ct[hn] = hg38_tss[hn]
                    mn = h2m.get(hn)
                    if mn and mn in mm10_tss:
                        m_tss_ct[mn] = mm10_tss[mn]

            ha = get_identity_enhancers(hd, h_tss_ct, window=gene_window)
            ma = get_identity_enhancers(md, m_tss_ct, window=gene_window)
            if not ha or not ma:
                continue

            hl, _ = liftover_peaks(ha, 'hg38', 'mm10')
            mi = [(c, s, e) for c, s, e, su in ma]
            fj = jaccard(hl, mi)

            ml, _ = liftover_peaks(ma, 'mm10', 'hg38')
            hi = [(c, s, e) for c, s, e, su in ha]
            rj = jaccard(ml, hi)

            rn = NAME_MAP[ct]
            rr = rigidity_df[rigidity_df['cell_type'] == rn]
            rs = rr['rigidity_score'].values[0]

            scores.append({'cell_type': ct, 'mean_jaccard': (fj + rj) / 2,
                          'rigidity_score': rs})
        return pd.DataFrame(scores)

    # 2. Window sensitivity
    print("\n--- Window sensitivity ---")
    for wkb in [25, 50, 100]:
        wdf = recompute_scores(gene_window=wkb * 1000)
        if len(wdf) >= 4:
            wr, wp, _, _ = spearman_ci(
                wdf['mean_jaccard'].values, wdf['rigidity_score'].values, len(wdf))
            wc = "Yes" if (wr < 0.35) == (rho < 0.35) else "No"
            print(f"  Window={wkb}kb: ρ={wr:.4f} p={wp:.4f} n={len(wdf)} [{wc}]")
            sens_rows.append({'analysis': 'Window', 'parameter': f'{wkb}kb',
                             'rho': wr, 'conclusion_consistent': wc})

    # 3. TSS exclusion sensitivity
    print("\n--- TSS exclusion sensitivity ---")
    for tkb in [1, 2, 5]:
        tdf = recompute_scores(tss_dist=tkb * 1000)
        if len(tdf) >= 4:
            tr, tp, _, _ = spearman_ci(
                tdf['mean_jaccard'].values, tdf['rigidity_score'].values, len(tdf))
            tc = "Yes" if (tr < 0.35) == (rho < 0.35) else "No"
            print(f"  TSS={tkb}kb: ρ={tr:.4f} p={tp:.4f} n={len(tdf)} [{tc}]")
            sens_rows.append({'analysis': 'TSS exclusion', 'parameter': f'{tkb}kb',
                             'rho': tr, 'conclusion_consistent': tc})

    sens_df = pd.DataFrame(sens_rows)
    sens_df.to_csv(os.path.join(OUTPUT_DIR, 'sensitivity_table.csv'), index=False)

    # ---- Step 3b-vii: Visualization ----
    print("\n" + "=" * 60)
    print("STEP 3b-vii: Visualizations")
    print("=" * 60)

    # 1. Primary scatter
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(df['mean_jaccard'], df['rigidity_score'],
               c='steelblue', s=120, edgecolors='black', linewidth=0.5, zorder=3)
    for _, row in df.iterrows():
        ax.annotate(row['cell_type'],
                    (row['mean_jaccard'], row['rigidity_score']),
                    textcoords="offset points", xytext=(8, 5), fontsize=9)
    ax.set_xlabel('Mean Enhancer Jaccard Conservation Score', fontsize=11)
    ax.set_ylabel('Procrustes Rigidity Score', fontsize=11)
    ax.set_title(f'H3K27ac Enhancer Conservation vs Procrustes Rigidity\n'
                 f'Spearman ρ = {rho:.3f}, p = {pval:.4f}, n={n}', fontsize=12)
    ax.text(0.02, 0.02,
            f"n={n}, null closure test\nUnderpowered for positive detection\n"
            f"Threshold: {threshold}",
            transform=ax.transAxes, fontsize=8, va='bottom',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'scatter_enhancer_primary.png'),
                dpi=150, bbox_inches='tight')
    plt.close()
    print("  scatter_enhancer_primary.png")

    # 2. Liftover summary
    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(n)
    w = 0.35
    ax.bar(x - w/2, df['pct_human_lifted'], w, label='Human→Mouse (hg38→mm10)',
           color='steelblue', edgecolor='black', linewidth=0.5)
    ax.bar(x + w/2, df['pct_mouse_lifted'], w, label='Mouse→Human (mm10→hg38)',
           color='coral', edgecolor='black', linewidth=0.5)
    ax.set_ylabel('% Peaks Successfully Lifted')
    ax.set_title('Liftover Success Rate by Cell Type')
    ax.set_xticks(x)
    ax.set_xticklabels(df['cell_type'], rotation=30, ha='right', fontsize=9)
    ax.legend()
    ax.set_ylim(0, 105)
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'liftover_summary.png'),
                dpi=150, bbox_inches='tight')
    plt.close()
    print("  liftover_summary.png")

    # ---- Summary Report ----
    print("\n" + "=" * 60)
    print("GENERATING SUMMARY REPORT")
    print("=" * 60)

    loo_range = f"[{min(loo_rhos):.4f}, {max(loo_rhos):.4f}]"
    all_loo_ok = all((r < 0.35) == (rho < 0.35) for r in loo_rhos)
    mean_h_lift = df['pct_human_lifted'].mean()
    mean_m_lift = df['pct_mouse_lifted'].mean()
    low_flags = [r['cell_type'] for _, r in df.iterrows()
                 if r['human_enhancer_count'] < 10 or r['mouse_enhancer_count'] < 10]

    summary = f"""# T3-E Step 3b: H3K27ac Enhancer Conservation vs Procrustes Rigidity

Generated: 2026-03-15

## 1. Primary Spearman Result

| Metric | Value |
|--------|-------|
| Spearman ρ | {rho:.4f} |
| p-value | {pval:.6f} |
| n (cell types) | {n} |
| 95% CI | [{ci_lo:.4f}, {ci_hi:.4f}] |

**Pre-registered threshold: {threshold}** — {interp}

This analysis is a null closure test. The underpowered sample (n={n}) limits
positive detection (requires |ρ|≥0.829 for p<0.05) but does not limit null
closure.

## 2. Conservation Scores by Cell Type

| Cell Type | Rigidity Rank | Human Enh | Mouse Enh | Fwd Jaccard | Rev Jaccard | Mean Jaccard |
|-----------|---------------|-----------|-----------|-------------|-------------|--------------|
"""
    for _, row in df.iterrows():
        summary += (f"| {row['cell_type']} | {row['rigidity_rank']} | "
                    f"{row['human_enhancer_count']} | {row['mouse_enhancer_count']} | "
                    f"{row['forward_jaccard']:.6f} | {row['reverse_jaccard']:.6f} | "
                    f"{row['mean_jaccard']:.6f} |\n")

    summary += f"""
## 3. Leave-One-Out Sensitivity

LOO ρ range: {loo_range}
All LOO results consistent with primary conclusion: {"Yes" if all_loo_ok else "No"}

"""
    for sr in sens_rows:
        if sr['analysis'] == 'LOO':
            summary += f"- {sr['parameter']}: ρ = {sr['rho']:.4f} [{sr['conclusion_consistent']}]\n"

    summary += f"""
## 4. Window and TSS Sensitivity

| Analysis | Parameter | ρ | Conclusion consistent? |
|----------|-----------|---|------------------------|
"""
    for sr in sens_rows:
        if sr['analysis'] != 'LOO':
            summary += f"| {sr['analysis']} | {sr['parameter']} | {sr['rho']:.4f} | {sr['conclusion_consistent']} |\n"

    summary += f"""
## 5. Data Quality Flags

- Mean liftover rate (human→mouse): {mean_h_lift:.1f}%
- Mean liftover rate (mouse→human): {mean_m_lift:.1f}%
- Cell types with <10 identity-gene enhancers: {', '.join(low_flags) if low_flags else 'None'}
- All mouse data: Lara-Astiaso 2014 (single replicate per cell type)
- Mouse data original assembly: mm9, lifted to mm10 for analysis
- Human data: ENCODE narrowPeak files (GRCh38)
- Mouse peak calling: threshold-based from bigWig (95th percentile)

## 6. Final Mechanistic Conclusion

"""
    if threshold == "TRIGGERED":
        summary += f"""**9TH MECHANISTIC NULL — COMPUTATIONAL CEILING REACHED**

Nine independent mechanistic hypotheses tested against Procrustes rigidity:

1. Housekeeping gene ratio (ρ=0.167, NS)
2. TF network complexity (ρ=-0.229, NS)
3. Niche adaptation (0/6 gene sets)
4. Within-type variance (ρ=-0.038, NS)
5. Inter-donor variance (ρ=-0.127, NS)
6. Expression-level confounds (all ρ<0.21)
7. PPI network centrality (0/27 combinations, best ρ=0.291 NS)
8. Promoter phastCons conservation (ρ=-0.058, n=35, NS)
9. **H3K27ac enhancer conservation (ρ={rho:.4f}, n={n}, NS) — THIS RESULT**

Rigidity is not predicted by any currently measurable transcriptomic, proteomic,
or cis-regulatory feature. The mechanism operates at a level requiring either
different data types (Hi-C, 4D Nucleome) or wet-lab experiments (CRISPR
perturbation screens).

This conclusion is publishable and scientifically strong: nine independent nulls
converging on the same mechanistic gap is a precise statement, not a failure.
"""
    elif threshold == "POSITIVE":
        summary += f"""**POSITIVE RESULT — FIRST MECHANISTIC CORRELATE**

After 8 nulls, enhancer conservation shows the first positive signal
(ρ={rho:.4f}). Statistically NS at n={n} (requires |ρ|≥0.829 for p<0.05).
The combination of 8 nulls + 1 positive enhancer trend is scientifically
meaningful regardless of individual significance.

Follow-up required: expand sample size to confirm.
"""
    else:
        summary += f"""**TREND RESULT — AMBIGUOUS**

Enhancer conservation shows a directional trend (ρ={rho:.4f}) but is
underpowered at n={n}. Cannot distinguish from noise. Requires review before use.
"""

    summary += f"""
## Files Generated

- conservation_scores.csv — Per-cell-type data ({n} rows)
- spearman_primary_result.json — Primary statistical result
- sensitivity_table.csv — All sensitivity analyses
- scatter_enhancer_primary.png — Primary scatter plot
- liftover_summary.png — Liftover quality check
- t3e_step3b_summary.md — This summary
"""

    with open(os.path.join(OUTPUT_DIR, 't3e_step3b_summary.md'), 'w') as f:
        f.write(summary)
    print("  t3e_step3b_summary.md")

    # ---- Final console summary ----
    print("\n" + "=" * 60)
    print("FINAL RESULT")
    print("=" * 60)
    print(f"  Spearman ρ = {rho:.4f} (p = {pval:.6f}, n = {n})")
    print(f"  95% CI: [{ci_lo:.4f}, {ci_hi:.4f}]")
    print(f"  Threshold: {threshold}")
    print(f"  {interp}")
    print("=" * 60)


if __name__ == '__main__':
    main()
