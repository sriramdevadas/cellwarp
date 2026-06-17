# Analysis Plan: Conserved-Contribution Gene-Set Gate
## Is the cross-species–conserved cell-type geometry carried by a coherent, expression-independent gene set?

*This is an internal analysis plan documenting the intended approach; it was
not registered with an external preregistration repository (consistent with the
manuscript's disclosure). It was frozen on 2026-06-05, before the checks below
were run, and is reproduced here as the pre-specification of record for the
conserved-contribution analysis (main Figure 2). Two pre-specifications are
included: Part I (the gene-set gate) and Part II (the donor-stability gate),
both frozen on the same date. The pre-specified positive-control transcription-
factor list (Part I, Section 7) was fixed before any check was computed.*

---

# Part I — Conserved-Contribution Gene-Set Gate

Frozen 2026-06-05, before running gate checks 1–4 and the secondary attribution.
Preliminary confirmations (data present; primary obs/null = 0.522; per-gene
Pearson r reproduced against the existing gene-conservation table; expression
Spearman ≈ 0.26) are descriptive and do not change any threshold below.

### Section 1 — Question

The 35-type human↔mouse Procrustes alignment is significantly conserved
(obs/null = 0.522, p = 1e-4). Is that conserved geometry carried by a
**coherent, expression-independent** set of genes whose **cell-type-specificity
pattern is conserved across species** — a "conserved-contribution" set that
(a) is more than just highly-expressed genes, (b) is biologically coherent, and
(c) recovers known identity biology? Meeting all pre-specified criteria in
Section 6 indicates such a set exists; failing to meet them is reported as a
negative result.

### Section 2 — Core quantity (non-circular headline)

For each ortholog g, the conservation score **C_g = Pearson r** between g's
expression profile across the 35 matched cell-type centroids in human vs mouse.
Pearson is per-gene mean-centred and scale-invariant, so it captures the
cell-type-specificity **pattern**, not magnitude (equivalent to z-scoring each
profile then correlating). It is computed directly on centroids, **independent
of the joint-PCA/Procrustes axes**. Genes with no defined r (constant profile in
either species) are excluded from ranking (~1,019). Robustness variant: Spearman
C_g.

Covariate "expression level" = mean of the gene's 35 human centroid values
(`mean_expression` in the existing table). It is **not** defined by loadings on
the between-type axes (that would be the circular `procrustes_contribution`
metric, used only for an explicit sanity contrast).

### Section 3 — Conserved-set definition (frozen)

- **Primary:** top quartile of valid C_g = **conserved set**; bottom quartile =
  **divergent set** (~3,985 genes each). Quartiles are symmetric and
  threshold-light.
- **Strict (sharpness):** C_g > 0.8 (~1,386 genes) for the positive control /
  enrichment.
- All conserved-set properties are evaluated against an **expression-matched
  background**, never against all genes.

### Section 4 — Expression-matched background (frozen procedure)

Bin all valid orthologs into 20 equal-frequency bins of `mean_expression`. An
expression-matched random set redraws, from each bin, the same count of genes
the target set has in that bin. Generate N = 1,000 matched sets (seeded). Every
conserved-set statistic gets a one-sided empirical p vs this matched null. This
removes the expression confound by construction.

### Section 5 — Gate checks & thresholds (frozen)

**Check 1 — Structure.** Hartigan's dip test on C_g (p < 0.05 ⇒ multimodal)
plus Gaussian-mixture BIC (1- vs 2-component). Structured/separable strengthens
the result; a smooth unimodal continuum means "thresholds arbitrary" — this
weakens but does **not** kill the result (per spec). Reported, not a hard gate.

**Check 2 — Expression-independence.** ρ = Spearman(C_g, mean_expression).
Frozen bands:
- |ρ| < 0.3 → conservation is NOT primarily an expression proxy → PASS check 2.
- 0.3 ≤ |ρ| < 0.5 → partial confound; the conserved set MUST still clear check 3
  vs the expression-matched background to pass overall.
- |ρ| ≥ 0.5 → conservation ≈ expression level → FAIL (signal is an expression
  artifact).
Also report decile-of-expression means of C_g and the partial structure.

**Check 3 — Coherence + positive control** (everything vs expression-matched
background):
- **3a TF positive control.** Pre-specified canonical lineage-defining master
  TFs (Section 7), intersected with the 16,959 orthologs that have a valid C_g.
  Statistic: median conservation percentile of the list. One-sided empirical p
  vs expression-matched draws of equal size (matching is essential — master TFs
  are typically low-abundance). PASS(3a) iff p < 0.05 AND > 50% of the testable
  TFs sit above the global median C_g.
- **3b Coherent enrichment.** Two pre-specified hypotheses, each vs
  expression-matched background, with BH-FDR across the tested families:
  - H3b-i: the conserved set is enriched for sequence-specific **TF / DNA-binding
    regulators** (local 296-TF CollecTRI list; GO "DNA-binding transcription
    factor activity" if Enrichr is reachable).
  - H3b-ii: the conserved set is enriched for **CellMarker 2.0 canonical identity
    markers** (reusing the Figure-S6 machinery).
  PASS(3b) iff ≥ 1 hypothesis is significant vs the expression-matched background
  after FDR, in the correct direction (ideally both). Descriptor: which identity
  programs / TF families dominate (coherence, not a smear).

**Check 4 — Stability (membership robustness).** Documented fallback: raw
Tabula Sapiens / Tabula Muris Senis single cells are absent from the deposit, so
no donor/cell resampling is possible at this stage (donor-level stability is
addressed separately in Part II). Frozen substitutes:
- **Type-jackknife:** leave-k-types-out (k = 1, k = 5) from the 35; recompute C_g
  and conserved-set membership; report median Jaccard of conserved-set membership
  plus Spearman of C_g vs full.
- **Cell-count restriction:** restrict to types with min(human, mouse) cell count
  above a threshold (a proxy for a cell-count cap, using
  `cell_type_inventory_passing.csv`); recompute; report membership overlap.
- **Estimator:** Pearson- vs Spearman-defined conserved-set Jaccard.
PASS(4) iff median conserved-set Jaccard ≥ 0.60 (leave-1-out) and ≥ 0.50
(leave-5-out), with high C_g rank-correlation. Labelled **type-level** stability
(not donor-level).

### Section 6 — Pre-specified decision criteria (frozen)

- The conserved-contribution set is supported iff: check 2 is not an expression
  proxy (|ρ| < 0.5 AND properties survive the matched background) AND 3a
  (identity-TF recovery) AND 3b (coherent expression-matched-significant
  enrichment) AND 4 (membership reasonably stable).
- It is reported as a negative result iff: conservation ≈ expression level, OR
  the conserved set is incoherent / fails the TF positive control, OR membership
  is unstable so the attribution is uninformative.
- An ambiguous outcome (e.g., a smooth continuum in check 1 but an
  expression-independent, coherent, stable set) is reported honestly as such.

### Section 7 — Pre-specified positive-control TF list (canonical lineage-defining / master regulators)

Mapped to the lineages of the 35 types; the final test is restricted to those
present in the 16,959 orthologs with a valid C_g.

- **T / NKT / thymocyte:** TCF7, LEF1, GATA3, TBX21, RUNX3, FOXP3, ETS1, BCL11B,
  ZEB1, TOX
- **NK:** EOMES, TBX21, ID2, NFIL3
- **B / plasma:** PAX5, EBF1, POU2AF1, SPIB, IRF4, PRDM1, XBP1, BACH2
- **Macrophage / monocyte / DC / granulocyte / neutrophil / microglia /
  myeloid:** SPI1, CEBPA, CEBPB, CEBPE, CEBPD, MAFB, IRF8, KLF4, BATF3
- **HSC / hematopoietic precursor:** GATA1, GATA2, TAL1, RUNX1, GFI1B, MEIS1,
  HLF, KLF1
- **Hepatocyte:** HNF4A, FOXA1, FOXA2, ONECUT1, NR1H4, HNF1A
- **Endothelial / vein endothelial:** ERG, FLI1, SOX17, SOX18, KLF2, FOXF1
- **Smooth muscle / cardiac myocyte / fibroblast / stromal / mesenchymal /
  adventitial:** MYOCD, SRF, GATA4, NKX2-5, TBX5, MEF2C, TCF21, PRRX1, TWIST2
- **Epithelial / basal / urothelial / luminal mammary / enterocyte / goblet /
  acinar / ductal:** TP63, GRHL2, ELF3, EHF, KLF5, CDX2, HNF1B, SPDEF, PTF1A,
  RBPJL, FOXA3, GATA6, SOX9, ASCL2

### Section 8 — Enrichment hypotheses (pre-specified)

- H1: the conserved set is enriched for sequence-specific TF / DNA-binding
  regulators (vs the expression-matched background).
- H2: the conserved set is enriched for CellMarker 2.0 canonical identity markers
  (vs the expression-matched background).
- (Descriptive) GO-BP convergence on identity / developmental programs, BH-FDR;
  descriptive only.

### Section 9 — Anti-circularity guards (explicit)

- The headline is the cross-species `pearson_r`, a comparison — never the
  within-data loading metric.
- The conserved set is always evaluated vs an **expression-matched** background —
  never vs all genes.
- "Conserved genes are conserved vs all genes" is baked in by 1:1-ortholog
  selection; the framing is always relative to the ortholog background only.
- No enrichment is reported without the pre-specified hypotheses plus
  multiple-testing correction.
- The secondary attribution is read as a sanity link, never as causal proof.

### Section 10 — Secondary geometry attribution (caveated, NOT the test)

Rebuild the joint-PCA + Procrustes obs/null using only (i) the conserved set,
(ii) the divergent set, and (iii) N expression-matched random sets of size equal
to the conserved set. **Caveat:** a low obs/null for the conserved set is partly
expected by selection (genes were chosen for cross-species pattern
conservation). The informative contrasts are conserved vs divergent, and
conserved vs expression-matched-random (does conservation beat expression-matched
genes at carrying the geometry?). This is a sanity link only.

---

# Part II — Donor-Stability Gate for C and the master-TF finding

Frozen 2026-06-05, before running the donor-split / ceiling / null / cap /
cross-protocol tests. This builds on Part I; **C and all controls are inherited
frozen** (Sections 2–9 above). The ONLY new operation is donor resampling.
Donor-structure facts (below) are descriptive and change no threshold.

### Section 1 — Frozen definitions (inherited, not redefined)

- **C** = per-gene cross-species Pearson r of the expression profile across the
  35 matched cell-type centroids (human vs mouse); computed on centroids,
  independent of PCA/Procrustes.
- **Master-TF positive control** = the curated lineage-defining TF list in
  Part I, Section 7.
- **Background** = expression-matched (20 equal-frequency bins on mean
  expression); the hardened control additionally matches specificity (Tau) via
  joint expression × Tau bins.
- **Conserved set** = top-quartile C. **obs/null** = Procrustes observed /
  permutation-null median (deposit anchor = 0.522).
- Pipeline: CELLxGENE Census 2025-11-08, the TS/TMS collections,
  `is_primary_data == True and disease == 'normal'`, `normalize_total(1e4)` +
  `log1p`, mean centroid (matching the deposit).

### Section 2 — Donor structure (reported, frozen context)

Human Tabula Sapiens: 804,639 cells / **24 donors** (10x 3'v3 718k/22d; 10x 5'v2
55k/5d; Smart-seq2 29.7k/17d; Smart-seq3 1.6k/2d). Mouse Tabula Muris Senis:
187,553 cells / **44 donors** (10x 3'v2 + Smart-seq2, 2 datasets). All 35 types
are donor-powered (human ≥ 15 donors for the abundant types, ≥ 4 for all but a
few; mouse all 35 types ≥ 15 donors). Donor-split is therefore feasible in
**both** species.

### Section 3 — Validity gate (must pass before any resampling)

Recompute centroids + C from the freshly-pulled Census cells (full pulled
population, per-type cap 10,000) and confirm reproduction of the deposit:
**Spearman/Pearson corr(C_new, C_deposit) ≥ 0.95 AND obs/null within
0.50–0.55** (anchor 0.522). If this is not reproduced, the pull/processing is
wrong; stop and fix before interpreting any stability result.

### Section 4 — Donor-split scheme (frozen)

- Per species, randomly partition the donor SET into two disjoint halves (human
  12/12, mouse 22/22), balanced by total cell count where possible. A donor is
  in the same half for all of its types. **100 random splits** (seeded).
- For each split and half, centroid_t = (Σ cells of type t from donors in the
  half) / count; compute C for half-A and half-B (each a full 35-type
  human-vs-mouse pair).
- Types with no cells in a half are dropped from that split's C (expected: none,
  given donor counts); report if any.

### Section 5 — Primary statistic (what the study builds on)

Per donor-half, does the **master-TF finding replicate independently**:
(i) **recovery** — median C-percentile of the master TFs; and
(ii) **enrichment** — conserved-set TF over-representation vs the
expression(+specificity)-matched background (empirical p).
Track (i) and (ii) in **each half across the 100 splits** (median + 2.5–97.5%
spread, and the fraction of halves that individually clear: median TF percentile
> 0.75 AND matched-background p < 0.05).

### Section 6 — Supporting / ceiling / null / power (frozen)

- **Cross-half C correlation:** Spearman(C_A, C_B) across 100 splits
  (distribution).
- **Conserved-membership agreement:** Jaccard of top-quartile (conserved)
  membership between halves; same for the bottom quartile.
- **Ceiling (cell-bootstrap test-retest within the SAME donors):** split each
  donor's cells in half (both halves keep identical donor composition),
  recompute C twice, correlate. 20 reps. The ceiling isolates cell-sampling
  noise; the **gap (ceiling − donor-split)** is the donor-specific instability —
  the standard localization of donor variance.
- **Null:** shuffle gene labels → chance cross-half C correlation (~0).
- **Power / cell-count cap:** rerun donor-split at per-type caps **{500, 2000,
  10000 (all pulled)}**. If cross-half C and the TF finding improve with cells,
  the result was power-limited; if flat and below the ceiling, it is genuine
  donor instability. Flag donor-underpowered types (few donors, or
  one-donor-dominated).
- **Cross-protocol (a within-atlas cross-"site" analogue):** compute C from
  10x-only centroids and from Smart-seq2-only centroids (per species), restricted
  to types powered in both protocols; report Spearman(C_10x, C_SS) and whether
  the master-TF finding holds in each. A donor-split pass plus a cross-protocol
  pass is much stronger than donor-split alone.

### Section 7 — Pre-specified stability criteria (frozen)

- **STABLE** iff: the master-TF finding replicates in **each donor half**
  (recovery median TF-percentile > 0.75 AND matched-background enrichment
  p < 0.05 in the large majority, ≥ 90%, of half-instances across splits), AND
  per-gene cross-half C Spearman beats the shuffle null and sits **close to the
  cell-bootstrap ceiling** (donor sampling adds little: ceiling − donor-split
  small, e.g. ≤ ~0.15), AND this holds at the cell-count cap (not a power
  artifact). If the cross-protocol check runs, it should also hold.
- **UNSTABLE** iff: the TF finding fails to replicate across halves, OR cross-half
  C collapses toward the null / far below the ceiling — i.e., C inherits a
  per-type ranking's donor-fragility.
- **AMBIGUOUS:** mixed (e.g., the TF finding is robust but cross-protocol fails,
  or there is a power-limited tail) → report honestly.

### Section 8 — Traps / do-not (frozen)

- Do not redefine C or the controls — this is a stability test of the existing
  result.
- Do not skip the validity reproduction (the data source switches from deposit
  centroids to Census cells).
- A high per-gene C correlation alone is NOT a pass (C is high-data and will
  correlate); the decision-relevant criterion is the master-TF headline
  replicating per donor-half.
- Do not conflate underpowered (small halves, noisy centroids) with unstable —
  that is what the ceiling and the cell-count cap separate.
- Donor-split is necessary, not sufficient: the full study must still add a
  cross-dataset (independent atlas pair) C replication. The within-atlas
  cross-protocol check is a down-payment on that, not a substitute.

---

# Registration Metadata

- **Date frozen:** 2026-06-05 (both Part I and Part II)
- **Registered by:** Sriram Devadas
- **Analysis to be executed by:** Sriram Devadas
- **Data source:** existing CellWarp 35-type human/mouse centroids
  (`output/phase2/scaled_35types/`) for Part I; CELLxGENE Census 2025-11-08
  (Tabula Sapiens / Tabula Muris Senis) for the Part II donor resampling.
- **Producers:** `analysis/conserved_contribution/` (gate and robustness) and
  `analysis/conserved_contribution/donor_stability/` (donor resampling).
- **Pre-specification commitment:** the positive-control TF list (Part I,
  Section 7) and all thresholds were fixed before any check was computed; the
  repository commit recording this document serves as the timestamp.
