# T3-E Step 3a: H3K27ac Enhancer ChIP-seq Feasibility Report

Generated: 2026-03-15

## Purpose

Assess whether sufficient publicly available H3K27ac ChIP-seq data exists in
matched human/mouse primary cell types to test whether enhancer-level regulatory
conservation at identity-gene loci predicts Procrustes rigidity across cell types.

This follows the 8th mechanistic null (promoter-level phastCons conservation,
rho=-0.058, p=0.740, n=35). Enhancer conservation is the last computationally
testable mechanistic hypothesis before the analysis reaches the experimental ceiling.

## Methodology

Queried the ENCODE REST API (encodeproject.org) for ALL released H3K27ac ChIP-seq
experiments classified as "primary cell" in both Homo sapiens (49 experiments) and
Mus musculus (6 experiments). Also searched:
- ENCODE tissue-level H3K27ac (human: 296 experiments, mouse: 91 experiments)
- ENCODE in vitro differentiated H3K27ac (27 experiments, human only)
- ENCODE Mint-ChIP-seq H3K27ac (173 experiments)
- Roadmap Epigenomics consolidated H3K27ac narrowPeak files (63 human epigenomes)
- Lara-Astiaso et al. 2014, Science (GSE60103) — mouse hematopoietic H3K27ac (16 cell types)

---

## Section 1: Per-Cell-Type Coverage Table

### ENCODE Primary Cell H3K27ac Inventory

**Human primary cell H3K27ac (49 experiments total):**

| Cell Type (ENCODE term) | Experiments | Reps | Relevant to our 35? |
|---|---|---|---|
| CD8-positive, alpha-beta T cell | ENCSR835OJV | 2 | Yes (rank 1) |
| CD8-positive, alpha-beta T cell | ENCSR007HLH | 1 | Yes — FLAG single rep |
| Naive CD8+ T cell | ENCSR810EPZ (memory) | 1 | Subtype |
| Naive CD4+ T cell | ENCSR120WKZ | 2 | Yes (rank 22, via naive) |
| Naive CD4+ T cell | ENCSR674VPA | 1 | FLAG single rep |
| CD4+ T cell | ENCSR546SDM | 1 | FLAG single rep |
| CD4+ T cell | ENCSR561KOM | 1 | FLAG single rep |
| CD4+ T cell (treated) | ENCSR138DOM | 2 | EXCLUDED (PMA/ionomycin) |
| CD4+ memory T cell | ENCSR724GUS | 2 | Subtype |
| CD4+ memory T cell | ENCSR314BEX | 1 | FLAG single rep |
| Effector memory CD4+ T | ENCSR892HPQ | 2 | Subtype |
| Treg (CD4+CD25+) | ENCSR577GVS | 1 | Not in our 35 |
| Th17 (treated) | ENCSR041UZZ | 1 | EXCLUDED (PMA/ionomycin) |
| T-cell (generic) | ENCSR222QLW | 1 | Rank 29, FLAG single rep |
| B cell | ENCSR000AUP | 2 | Yes (rank 19) |
| B cell | ENCSR191ZQT | 1 | FLAG single rep |
| Natural killer cell | ENCSR391EQV | 1 | Yes (rank 13), FLAG single rep |
| CD14-positive monocyte | ENCSR000ASJ | 2 | Yes (rank 14) |
| CD14-positive monocyte | ENCSR102GGG | 1 | FLAG single rep |
| CD14-positive monocyte | ENCSR065SKV | 1 | FLAG single rep |
| CD14-positive monocyte | ENCSR012PII | 1 | FLAG single rep |
| Neutrophil | ENCSR267YXV | 3 | Yes (rank 28) |
| Endothelial (HUVEC) | ENCSR000ALB | 3 | Yes (rank 3), umbilical vein |
| Brain microvasc. endothelial | ENCSR111QCU | 2 | Subtype of rank 3 |
| Fibroblast of dermis | ENCSR000APN | 2 | Yes (rank 8) |
| Fibroblast of lung | ENCSR000AMR | 2 | Variant of rank 8 |
| Foreskin fibroblast | ENCSR917QEH, ENCSR822ZIG, ENCSR108NVQ | 2,1,1 | Variant |
| Fibroblast of breast | ENCSR953ULK, ENCSR579YLO | 1,1 | Variant |
| Mammary epithelial cell | ENCSR000ALW | 2 | Rank 18 proxy |
| Epithelial cell of prostate | ENCSR910PDW | 2 | Rank 34 proxy |
| Foreskin keratinocyte | ENCSR666TFS, ENCSR736ZEG | 2,2 | Not in our 35 |
| Foreskin keratinocyte (treated) | ENCSR637ISS, ENCSR709ABP | 2,2 | EXCLUDED (calcium) |
| Keratinocyte | ENCSR000ALK | 3 | Not in our 35 |
| Foreskin melanocyte | ENCSR227FYJ | 2 | Not in our 35 |
| Skeletal muscle myoblast | ENCSR000ANF | 2 | Not in our 35 |
| Osteoblast | ENCSR000APH | 2 | Not in our 35 |
| Astrocyte | ENCSR000AOQ | 2 | Not in our 35 |
| CMP CD34+ | ENCSR891KSP | 2 | HSC proxy (rank 32) |
| CMP CD34+ | ENCSR620AZM | 1 | FLAG single rep |
| PBMC | ENCSR615HXA, ENCSR105EMQ, ENCSR156XNC | 1,1,1 | Mixed — not usable |
| Neurosphere | ENCSR379WXM | 1 | Not in our 35 |

**Mouse primary cell H3K27ac (6 experiments total):**

| Cell Type (ENCODE term) | Experiment | Reps | Notes |
|---|---|---|---|
| Bone marrow macrophage | ENCSR000CFD | 2 | Rank 20, Bruce4 strain, adult 2mo |
| Macrophage | ENCSR933GNW | 2 | Rank 20, C57BL/6J |
| Embryonic fibroblast | ENCSR000CDI | 2 | Rank 8, FLAG: embryonic (E13.5) |
| Activated Treg | ENCSR585IGO | 2 | EXCLUDED (activated) |
| Treg | ENCSR896QWJ | 1 | Not in our 35, single rep |
| Inflammation-experienced Treg | ENCSR528PQK | 2 | EXCLUDED (inflammation) |

### ENCODE Tissue H3K27ac — Relevant Adult Tissues

**Human tissue (296 total, relevant subset):**

| Tissue | Exps with ≥2 reps | Cell type proxy | Rank |
|---|---|---|---|
| Liver (right lobe/liver) | ENCSR864OOO (2r), ENCSR779ISJ (2r), ENCSR119XNK (2r), ENCSR981UJA (2r) | Hepatocyte (60-80% of parenchyma) | 4 |
| Spleen | ENCSR593INW (4r), ENCSR726HTS (4r), ENCSR668GBL (4r), etc. | Mixed B/T/myeloid — NOT usable | — |
| Pancreas | ENCSR245GEV (2r), ENCSR868ZOR (2r), ENCSR029SIG (4r) | Mixed acinar/ductal/islet — NOT usable | — |
| Heart | Multiple 2-rep exps | Not in our 35 (cardiomyocyte dropped) | — |
| Arteries/aorta | ENCSR015GFK (2r), ENCSR355GNZ (2r), etc. | Endothelial mix with SMC/fibroblast — NOT usable | — |
| Small/large intestine | Multiple | Enterocyte/goblet mix — NOT usable | — |
| Thymus | ENCSR303IKJ (1r, child 3y) | T cell enriched but mixed — NOT usable | — |

**Mouse tissue (91 total, relevant adult subset):**

| Tissue | Experiment | Reps | Cell type proxy |
|---|---|---|---|
| Liver | ENCSR000CDH | 2 | Hepatocyte, adult 2mo |
| Spleen | ENCSR000CDJ | 2 | Mixed — NOT usable |
| Bone marrow | ENCSR000CCL | 2 | HSC/myeloid mix — NOT usable |
| Thymus | ENCSR000CCH | 2 | T cell mix — NOT usable |
| Heart | ENCSR000CDF | 2 | Not in our 35 |
| Small intestine | ENCSR000CCQ | 2 | Mixed — NOT usable |
| Kidney | ENCSR000CDG | 2 | Not in our 35 |

All other mouse tissue experiments are embryonic (E10.5-E16.5) or postnatal (P0).

### Mint-ChIP-seq H3K27ac — Primary Cell Subset

| Cell Type | Experiment | Reps | Species | Notes |
|---|---|---|---|---|
| T-helper 17 cell | ENCSR819NCZ | 1 | Human | Untreated |
| Memory B cell | ENCSR886LDA | 1 | Human | Untreated |
| Activated naive CD4+ T | ENCSR549RQM | 1 | Human | EXCLUDED (IL-2 + anti-CD3/CD28) |
| Activated naive CD8+ T | ENCSR033WMA | 1 | Human | EXCLUDED (anti-CD3/CD28 + IL-15) |
| Activated naive B cell | ENCSR709KME | 1 | Human | EXCLUDED (anti-CD40/IgM/IL-4) |
| Effector memory CD8+ T | ENCSR917BHP | 1 | Human | Untreated |
| Treg | ENCSR322MTA | 1 | Human | Untreated |
| Naive CD8+ T | ENCSR797POI | 1 | Human | Untreated |
| NK cell | ENCSR067KOO | 1 | Human | Untreated |

All Mint-ChIP experiments are **single replicate** and **human only**. No mouse Mint-ChIP H3K27ac.

### Roadmap Epigenomics H3K27ac (Human Only)

63 reference epigenomes with H3K27ac narrowPeak files. **Human only — no mouse equivalent.**

Relevant samples:

| EID | Cell Type | Notes |
|---|---|---|
| E029 | CD14+ monocytes (peripheral blood) | Primary cell |
| E032 | B cells (peripheral blood) | Primary cell |
| E034 | T cells (peripheral blood) | Primary cell |
| E046 | NK cells (peripheral blood) | Primary cell |
| E047 | CD8+ naive T cells (peripheral blood) | Primary cell |
| E048 | CD8+ memory T cells (peripheral blood) | Primary cell |
| E038 | CD4+ helper T cells (peripheral blood) | Primary cell |
| E039-E045 | T helper subsets (Th1/Th2/Th17/Treg/memory) | Primary cells |
| E050 | T effector/memory enriched | Primary cell |
| E062 | PBMC | Mixed — not usable |

**No liver, hepatocyte, endothelial, fibroblast, macrophage, smooth muscle, or
stromal samples in Roadmap H3K27ac.**

### Tier 2: Lara-Astiaso et al. 2014 (GSE60103) — Mouse Hematopoietic H3K27ac

Lara-Astiaso D, Ciabrelli F, Vasseur P et al. "Chromatin state dynamics during
blood formation." Science 345:943-949, 2014.

16 mouse hematopoietic cell types with H3K27ac ChIP-seq:

| GSM | Cell Type | Our 35-type mapping | Rank |
|---|---|---|---|
| GSM1441269 | LT-HSC | HSC (rank 32) | 32 |
| GSM1441270 | ST-HSC | HSC | 32 |
| GSM1441271 | MPP | Hematopoietic precursor (rank 33) | 33 |
| GSM1441272 | CMP | Hematopoietic precursor | 33 |
| GSM1441273 | GMP | Not direct match | — |
| GSM1441274 | MEP | Not in our 35 | — |
| GSM1441275 | EryA | Not in our 35 | — |
| GSM1441276 | EryB | Not in our 35 | — |
| GSM1441277 | GN (granulocyte/neutrophil) | Granulocyte (9) / Neutrophil (28) | 9/28 |
| GSM1441278 | Mono (monocyte) | Monocyte (rank 14) | 14 |
| GSM1441279 | MF (macrophage) | Macrophage (rank 20) | 20 |
| GSM1441280 | B cell | B cell (rank 19) | 19 |
| GSM1441281 | CD4 T cell | CD4+ T cell (rank 22) | 22 |
| GSM1441282 | CD8 T cell | CD8+ T cell (rank 1) | 1 |
| GSM1441283 | NK cell | NK cell (rank 13) | 13 |
| GSM1441284 | CLP | Not in our 35 | — |

**Critical limitation: 1 biological replicate per cell type.** No replication within
the Lara-Astiaso dataset. Published 2014 — genome assembly likely mm9 (would require
liftOver to mm10 or re-alignment to match ENCODE mm10 data).

Peak files: raw data available on GEO (GSE60103_RAW.tar, 7.9 GB). Processed peak
files may require peak calling from raw FASTQ/BAM — needs verification.

---

## Section 1b: Combined Coverage Summary (All Sources)

| Cell Type | Rank | Human H3K27ac | Mouse H3K27ac | Matched? |
|---|---|---|---|---|
| CD8+ T cell | 1 | ENCODE ENCSR835OJV (2r) | Lara-Astiaso CD8 (1r) | **YES** |
| Non-classical monocyte | 2 | None | None | NO |
| Endothelial cell | 3 | ENCODE ENCSR000ALB HUVEC (3r) | None | NO |
| Hepatocyte | 4 | Tissue: liver (4 exps, 2r each) | Tissue: liver ENCSR000CDH (2r) | **PROXY** |
| Smooth muscle cell | 5 | In vitro only (H9-derived) | None | NO |
| Pancreatic ductal cell | 6 | None | None | NO |
| Bladder urothelial cell | 7 | None | None | NO |
| Fibroblast | 8 | ENCODE ENCSR000APN (2r) | ENCODE ENCSR000CDI (2r, **embryonic**) | **FLAGGED** |
| Granulocyte | 9 | ENCODE neutrophil (3r) | Lara-Astiaso GN (1r) | **YES** (=neutrophil) |
| Adventitial cell | 10 | None | None | NO |
| Mature NK T cell | 11 | None | None | NO |
| Intermediate monocyte | 12 | None | None | NO |
| Natural killer cell | 13 | ENCODE ENCSR391EQV (1r) | Lara-Astiaso NK (1r) | **YES** (both 1r) |
| Monocyte | 14 | ENCODE ENCSR000ASJ CD14+ (2r) | Lara-Astiaso Mono (1r) | **YES** |
| Myeloid dendritic cell | 15 | None primary | None | NO |
| Enterocyte | 16 | None | None | NO |
| Large intestine goblet cell | 17 | None | None | NO |
| Luminal epi mammary | 18 | ENCODE ENCSR000ALW (2r) | None | NO |
| B cell | 19 | ENCODE ENCSR000AUP (2r) | Lara-Astiaso B (1r) | **YES** |
| Macrophage | 20 | **None** | ENCODE (2r) + Lara-Astiaso (1r) | **NO** (human absent) |
| Classical monocyte | 21 | ENCODE CD14+ (same as monocyte) | Lara-Astiaso Mono (same) | DUPLICATE of rank 14 |
| CD4+ T cell | 22 | ENCODE ENCSR120WKZ naive (2r) | Lara-Astiaso CD4 (1r) | **YES** |
| Mesenchymal stem cell | 23 | In vitro only (H1-derived) | None | NO |
| Plasma cell | 24 | None | None | NO |
| MSC adipose | 25 | None | None | NO |
| Myeloid leukocyte | 26 | None | None | NO |
| Fibroblast cardiac | 27 | None (dermal fibroblast is different) | None | NO |
| Neutrophil | 28 | ENCODE ENCSR267YXV (3r) | Lara-Astiaso GN (1r) | **YES** |
| T cell (generic) | 29 | ENCODE ENCSR222QLW (1r) | Ambiguous mapping | NO (ambiguous) |
| Basal cell | 30 | None | None | NO |
| Pancreatic acinar cell | 31 | None | None | NO |
| HSC | 32 | CMP CD34+ ENCSR891KSP (2r) | Lara-Astiaso LT/ST-HSC (1r) | **PARTIAL** (CMP ≠ HSC) |
| Hematopoietic precursor | 33 | CMP CD34+ (same as above) | Lara-Astiaso MPP/CMP (1r) | **PARTIAL** |
| Epithelial cell | 34 | ENCODE prostate epi (2r) | None | NO |
| Stromal cell | 35 | None | None | NO |

---

## Section 2: Matched Pair List

### Clean primary cell matched pairs (n=6)

| # | Cell Type | Rank | Human Source | H.Reps | Mouse Source | M.Reps | Flags |
|---|---|---|---|---|---|---|---|
| 1 | CD8+ T cell | 1 | ENCSR835OJV | 2 | Lara-Astiaso GSM1441282 | 1 | Mouse single rep |
| 2 | NK cell | 13 | ENCSR391EQV | 1 | Lara-Astiaso GSM1441283 | 1 | Both single rep |
| 3 | Monocyte (CD14+) | 14 | ENCSR000ASJ | 2 | Lara-Astiaso GSM1441278 | 1 | Mouse single rep |
| 4 | B cell | 19 | ENCSR000AUP | 2 | Lara-Astiaso GSM1441280 | 1 | Mouse single rep |
| 5 | CD4+ T cell | 22 | ENCSR120WKZ (naive) | 2 | Lara-Astiaso GSM1441281 | 1 | Mouse single rep |
| 6 | Neutrophil | 28 | ENCSR267YXV | 3 | Lara-Astiaso GSM1441277 (GN) | 1 | Mouse single rep; mouse term is "granulocyte" |

### Flagged matched pairs (additional n=2)

| # | Cell Type | Rank | Human Source | Mouse Source | Issue |
|---|---|---|---|---|---|
| 7 | Hepatocyte | 4 | Tissue: ENCSR864OOO liver (2r) | Tissue: ENCSR000CDH liver (2r) | Both bulk tissue, not primary cell. Hepatocytes 60-80% of liver but signal includes stellate, Kupffer, endothelial. |
| 8 | Fibroblast | 8 | ENCSR000APN dermis (2r) | ENCSR000CDI embryonic (2r) | Mouse is embryonic (E13.5) vs human adult. H3K27ac landscape differs between embryonic and adult fibroblasts. |

### Deduplication notes

- **Granulocyte (rank 9) and Neutrophil (rank 28)** map to the same Lara-Astiaso
  sample (GSM1441277 "GN"). These are counted as ONE data point (neutrophil, rank 28).
  Cannot use the same mouse sample for two different rigidity ranks.
- **Monocyte (rank 14) and Classical monocyte (rank 21)** map to the same ENCODE
  and Lara-Astiaso samples. Counted as ONE data point (monocyte, rank 14).
- **HSC (rank 32)** excluded: human CMP CD34+ is not equivalent to mouse LT-HSC.
  Different progenitor stages with distinct H3K27ac profiles.

### Asymmetric cell types (data in one species only)

| Cell Type | Rank | Has Human | Has Mouse | Blocker |
|---|---|---|---|---|
| Endothelial | 3 | ENCODE HUVEC (3r) | None | No mouse endothelial H3K27ac |
| Macrophage | 20 | None | ENCODE (2r) + Lara-Astiaso | No human macrophage H3K27ac |
| Mammary epi | 18 | ENCODE (2r) | None | No mouse mammary H3K27ac |
| Epithelial | 34 | ENCODE prostate (2r) | None | No mouse epithelial H3K27ac |

Human macrophage H3K27ac is completely absent from ENCODE, ENCODE in vitro,
Mint-ChIP, and Roadmap Epigenomics. Macrophage is the only cell type with rich
mouse data (two ENCODE experiments + Lara-Astiaso) but zero human data.

---

## Section 3: Range Coverage Assessment

### Rigidity ranks of matched pairs

**n=6 clean primary cell pairs:**

```
Rigid half (ranks 1-17):
  Rank  1: CD8+ T cell      ████████████████████████████████████ (most rigid)
  Rank 13: NK cell           ████████████████████████
  Rank 14: Monocyte          ███████████████████████

Flexible half (ranks 18-35):
  Rank 19: B cell            ███████████████████
  Rank 22: CD4+ T cell       ████████████████
  Rank 28: Neutrophil        █████████  (80th percentile)

Coverage: ranks 1-28 (80% of full range)
Split: 3 rigid / 3 flexible (even)
```

**n=8 with flagged proxies:**

```
Rigid half (ranks 1-17):
  Rank  1: CD8+ T cell       (most rigid)
  Rank  4: Hepatocyte        ** tissue proxy **
  Rank  8: Fibroblast        ** embryonic mouse **
  Rank 13: NK cell
  Rank 14: Monocyte

Flexible half (ranks 18-35):
  Rank 19: B cell
  Rank 22: CD4+ T cell
  Rank 28: Neutrophil

Coverage: ranks 1-28 (80% of full range)
Split: 5 rigid / 3 flexible (slightly rigid-heavy)
```

### Range assessment

**NOT range-restricted.** Both n=6 and n=8 matched sets span ranks 1-28, covering
both the rigid half (CD8+ T at rank 1) and the flexible half (neutrophil at rank 28).
The most flexible quarter (ranks 27-35) has one entry (neutrophil at 28); the most
flexible types (HSC rank 32, hematopoietic precursor rank 33, epithelial rank 34,
stromal rank 35) are absent due to no matched H3K27ac data.

**Comparison to T3-E ATAC-seq feasibility (which was BLOCKED):**
The ATAC-seq feasibility check found n=6-7 but range was restricted to the rigid half
only (Calderon+ImmGen immune cells clustered in ranks 1-14). The H3K27ac situation is
structurally better: the matched set spans both halves with an even 3/3 split at n=6.
This is the critical difference — a correlation spanning ranks 1-28 is interpretable
regardless of sample size concerns.

### Missing range segments

| Range segment | Types present | Types absent |
|---|---|---|
| Ranks 1-5 (most rigid) | CD8+T (1), Hepatocyte* (4) | Non-cl mono (2), Endothelial (3), SMC (5) |
| Ranks 6-12 | Fibroblast* (8) | Pancreatic ductal (6), Bladder (7), Granulocyte (9), Adventitial (10), NKT (11), Interm. mono (12) |
| Ranks 13-17 | NK (13), Monocyte (14) | Myeloid DC (15), Enterocyte (16), Goblet (17) |
| Ranks 18-24 | B cell (19), CD4+T (22) | Mammary epi (18), Macrophage (20), Class. mono (21), MSC (23), Plasma (24) |
| Ranks 25-35 (most flexible) | Neutrophil (28) | 9 types absent |

*asterisk = flagged proxy

---

## Section 4: Power Assessment

### Statistical thresholds (Spearman, two-tailed, alpha=0.05)

| n | |rho| required | Scenario |
|---|---|---|
| 6 | 0.829 | Clean primary cell pairs |
| 7 | 0.786 | If macrophage human source found |
| 8 | 0.738 | Including tissue/embryonic proxies |
| 9 | 0.683 | If one additional type identified |
| 10 | 0.648 | Expansion target |
| 12 | 0.591 | Optimistic expansion |
| 15 | 0.521 | Theoretical maximum |
| 20 | 0.447 | Not achievable with current data |
| 35 | 0.334 | Full coverage (not achievable) |

### Realistic assessment

**At n=6 (clean):** Need |rho|>=0.829 for significance. The pre-registered positive
threshold (rho>=0.50) falls BELOW the significance boundary. A moderate positive
effect (rho=0.50-0.80) would be directionally positive but statistically NS — the
same underpowered situation as T3-C (n=5, rho=0.600, p=0.285).

**At n=8 (with proxies):** Need |rho|>=0.738. Still requires strong effect. The
tissue/embryonic flags add uncertainty — are tissue-level H3K27ac peaks meaningful
proxies for cell-type-specific enhancers? For liver (hepatocyte-dominated), arguably
yes; for embryonic vs adult fibroblast, much more questionable.

**The null closure threshold (rho<0.35) IS testable:** At any n>=6, a clearly null
result (rho near zero, like the eight previous nulls with rho ranging -0.06 to 0.29)
would confidently trigger closure. The underpowerment only affects the POSITIVE
detection, not the null detection.

### Comparison to previous T3-E analyses

| Analysis | n | rho | p | Outcome |
|---|---|---|---|---|
| T3-E Step 2 (promoter phastCons) | 35 | -0.058 | 0.740 | 8th null |
| T3-E ATAC feasibility | 6-7 | — | — | BLOCKED (range-restricted) |
| T3-E H3K27ac feasibility | 6-8 | — | — | **This report** |

The H3K27ac analysis is structurally better than ATAC because:
1. Range spans both halves (not range-restricted)
2. Lara-Astiaso provides a unified mouse source (not multi-study patchwork)
3. H3K27ac specifically marks active enhancers (cleaner signal than ATAC open chromatin)

But structurally worse because:
1. All mouse immune data is single replicate (Lara-Astiaso)
2. Peak files from Lara-Astiaso may need reprocessing (2014 paper, likely mm9)
3. Source heterogeneity (ENCODE human + Lara-Astiaso mouse)

---

## Section 5: Enhancer Definition Feasibility

### Peak file format assessment

All ENCODE experiments have **narrowPeak** files available (both replicated peaks and
pseudoreplicated peaks in BED6+4 format, GRCh38 for human, mm10 for mouse). This is
the optimal format for the enhancer analysis:

1. **Distance-to-TSS filtering:** Can directly intersect H3K27ac narrowPeak BED files
   with our existing `identity_gene_tss_hg38.bed` (from T3-E step 2) to identify
   distal peaks (>2kb from any TSS). Standard bedtools operation.

2. **No peak calling required** for ENCODE experiments — peaks are pre-called by the
   ENCODE uniform processing pipeline.

3. **Lara-Astiaso peaks:** Raw data on GEO (FASTQ). Peak calls may need to be
   generated from BAM files using MACS2. Genome assembly likely mm9 — requires
   liftOver to mm10 for compatibility with ENCODE mouse tissue data, OR re-alignment
   to mm10.

### Enhancer identification pipeline

For each matched pair:
1. Take H3K27ac narrowPeak file
2. Remove peaks within 2kb of any annotated TSS (using GENCODE annotation)
3. Remaining distal peaks = active enhancer candidates
4. Intersect enhancer peaks with identity gene loci (e.g., within 500kb of identity
   gene TSS, or within TAD boundaries if available)
5. Compute enhancer conservation metric per cell type:
   - Option A: Fraction of human enhancer peaks that liftOver to mouse enhancer peaks
     at orthologous loci
   - Option B: Jaccard similarity of enhancer peaks at orthologous identity gene loci
6. Correlate enhancer conservation score with Procrustes rigidity rank

**Technical feasibility: HIGH.** All necessary tools exist (bedtools, liftOver,
MACS2). The identity gene TSS BED file already exists. The analysis is
computationally straightforward once peak files are obtained.

---

## Section 6: Blockers and Go/No-Go Recommendation

### 6.1 Assessment against decision criteria

| Criterion | Status | Details |
|---|---|---|
| n >= 10 AND full range | **NOT MET** | n=6-8, range spans both halves but n<10 |
| n >= 10 BUT range-restricted | N/A | Range is NOT restricted |
| n < 10 regardless of range | **TRIGGERED** | n=6 clean, n=8 with proxies |

### 6.2 Go/No-Go Recommendation

**UNDERPOWERED — advisor decision required before proceeding.**

The H3K27ac enhancer feasibility is structurally superior to the ATAC-seq feasibility
(which was BLOCKED by range restriction), but falls in the same underpowered territory
as T3-C tissue-stratified rigidity (n=5, rho=0.600, p=0.285).

**Arguments FOR proceeding:**

1. **Range coverage is good.** Unlike ATAC (range-restricted to rigid half), the H3K27ac
   matched set spans ranks 1-28 with an even 3/3 rigid/flexible split. A correlation
   is interpretable even at small n.

2. **Null detection is powered.** If the enhancer hypothesis is null (like the previous
   8 nulls), a near-zero rho will be clearly detectable. Only positive detection is
   underpowered.

3. **This is the last computationally testable hypothesis.** If we don't test it, the
   mechanistic story stops at "8 nulls, chromatin is the surviving hypothesis but
   untested." Testing it — even underpowered — provides more information than not
   testing it.

4. **A directional positive (rho=0.50-0.80, NS) combined with 8 nulls is scientifically
   meaningful** even if not individually significant. The pattern of 8 null correlations
   and one positive trend at the enhancer level tells a coherent mechanistic story.

5. **Computational cost is modest.** Peak file downloads + bedtools intersections +
   one Spearman correlation. Less than one day of compute.

**Arguments AGAINST proceeding:**

1. **Single replicate mouse data.** All Lara-Astiaso samples are unreplicated. Peak
   calls from single replicates have higher false positive/negative rates. IDR
   (irreproducible discovery rate) cannot be computed without replicates.

2. **Source heterogeneity.** ENCODE (human, uniform pipeline) + Lara-Astiaso (mouse,
   2014, different pipeline, likely mm9). Protocol and computational differences could
   create systematic bias.

3. **Underpowered positive is hard to interpret.** If rho=0.60 p=0.20 (like T3-C), we
   have another "suggestive but NS" result. Does this advance the story or just add
   ambiguity?

4. **Tissue proxy quality.** Liver H3K27ac as hepatocyte proxy mixes signal from
   multiple cell types. For distal enhancers (which are more cell-type-specific than
   promoters), this contamination is worse than for promoter-level analysis.

### 6.3 Potential expansion paths to increase n

| Path | Types added | New n | Feasibility |
|---|---|---|---|
| Include hepatocyte (tissue proxy) | +1 | 7 | Moderate — liver is hepatocyte-dominated |
| Include fibroblast (embryonic mouse) | +1 | 8 | Low — embryonic H3K27ac landscape differs |
| Find human macrophage H3K27ac (BLUEPRINT, IHEC, literature) | +1 | 7-9 | Unknown — no ENCODE source exists |
| Map Lara-Astiaso GMP to myeloid leukocyte (rank 26) | +1 | 7-9 | Tenuous — GMP is progenitor, not mature |
| Use Lara-Astiaso LT-HSC for HSC (rank 32) + ENCODE CMP CD34+ for human | +1 | 7-9 | Low — CMP and LT-HSC are different stages |
| Search for non-ENCODE mouse endothelial H3K27ac | +1 | 7-9 | Unknown — requires literature search |

**Best-case realistic n: 8** (adding hepatocyte tissue proxy and fibroblast with
embryonic flag). This requires |rho|>=0.738, still stringent but more accessible
than |rho|>=0.829.

### 6.4 Cell-type mapping decisions requiring advisor confirmation

1. **Monocyte mapping:** Human CD14+ monocyte (ENCODE) mapped to Lara-Astiaso "Mono"
   (mouse bone marrow monocyte). CD14+ selects classical monocytes. Reasonable match.

2. **Granulocyte/Neutrophil deduplication:** Lara-Astiaso "GN" sample is used for
   neutrophil (rank 28). Cannot also use for granulocyte (rank 9) — same data point.

3. **Hepatocyte tissue proxy:** Both human and mouse liver tissue H3K27ac used as
   hepatocyte proxy. Acceptable? Liver is 60-80% hepatocytes but enhancers are more
   cell-type-specific than promoters.

4. **Fibroblast embryonic mouse:** ENCODE mouse embryonic fibroblast (E13.5) paired
   with human adult dermal fibroblast. H3K27ac landscape changes substantially during
   development. Include or exclude?

5. **Naive CD4+ T as CD4+ T proxy:** Human ENCSR120WKZ is naive CD4+ T, mapped to
   our generic "CD4-positive, alpha-beta T cell" (rank 22). Mouse Lara-Astiaso CD4 is
   total CD4+ T. Subtype mismatch but both are unstimulated CD4+ T.

### 6.5 Failed queries and data gaps

No API failures. All ENCODE queries returned valid results (49 human primary, 6 mouse
primary, 296 human tissue, 91 mouse tissue, 173 Mint-ChIP, 27 in vitro).

The fundamental data gap is the same as ATAC-seq: **ENCODE's historical focus on human
immune cells and mouse developmental biology creates systematic human/mouse asymmetry
for cell-type-resolved epigenomic data.** This is a structural limitation of the
available public data, not a query failure.

### 6.6 Summary table

```
DATA AVAILABILITY SUMMARY
=========================

                    ENCODE primary    ENCODE tissue    Lara-Astiaso    COMBINED
                    Human  Mouse      Human  Mouse     (Mouse only)    Both?
CD8+ T cell         2r     0          -      -         1r              YES
NK cell             1r     0          -      -         1r              YES (both 1r)
Monocyte (CD14+)    2r     0          -      -         1r              YES
B cell              2r     0          -      -         1r              YES
CD4+ T cell         2r     0          -      -         1r              YES
Neutrophil          3r     0          -      -         1r (GN)         YES
Hepatocyte          0      0          4exps  1exp      -               PROXY (tissue)
Fibroblast          2r     2r(emb)    -      -         -               FLAGGED (embryonic)
Macrophage          0      2r+2r      -      -         1r              NO (human absent)
Endothelial         3r     0          -      0         -               NO (mouse absent)
All others          -      -          -      -         -               NO

Clean matched pairs: n = 6 (ranks 1, 13, 14, 19, 22, 28)
With proxies/flags:  n = 8 (add ranks 4, 8)
Range coverage:      Both halves, ranks 1-28 (80%)
Range-restricted:    NO
Power at n=6:        |ρ| ≥ 0.829 required
Power at n=8:        |ρ| ≥ 0.738 required
```

### 6.7 Final recommendation

**UNDERPOWERED but NOT BLOCKED.** Unlike the ATAC-seq feasibility (which was blocked
by range restriction), the H3K27ac enhancer analysis is structurally feasible: range
spans both halves, peak files are in the correct format, and the analysis pipeline is
straightforward. The limitation is statistical power (n=6-8, needing |rho|>=0.738-0.829).

The decision tree:
- If the enhancer correlation is clearly null (rho near zero): **9th mechanistic null,
  computational ceiling reached.** Write the paper with 9 nulls pointing to a mechanism
  below cis-regulatory elements — possibly 3D genome organization or trans-regulatory
  networks.
- If the enhancer correlation is strongly positive (rho>=0.80): **First positive
  mechanistic result.** Statistically significant even at n=6. Changes the paper
  substantially.
- If the enhancer correlation is moderately positive (rho=0.50-0.75): **Suggestive
  positive, NS.** Same situation as T3-C. Report as directional evidence alongside
  the 8 nulls. The contrast (8 nulls + 1 enhancer trend) tells a coherent story even
  without individual significance.

**Advisor decision required:** Proceed with n=6-8 knowing the underpowerment, or stop
the mechanistic pursuit at 8 nulls and write the paper as-is?

---

## Appendix: Raw Data Checkpoints

All ENCODE API query results saved to:
```
output/validation/t3e_enhancer_feasibility/encode_raw/
├── human_primary_cell_h3k27ac.json   (49 experiments)
├── mouse_primary_cell_h3k27ac.json   (6 experiments)
├── human_tissue_h3k27ac.json         (296 experiments)
├── mouse_tissue_h3k27ac.json         (91 experiments)
└── mint_chip_h3k27ac.json            (173 experiments)
```

Lara-Astiaso et al. 2014 metadata verified at GEO GSE60103 (16 H3K27ac samples,
GSM1441269-GSM1441284).

Roadmap Epigenomics H3K27ac narrowPeak file list verified at:
https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/
(63 epigenomes with H3K27ac).
