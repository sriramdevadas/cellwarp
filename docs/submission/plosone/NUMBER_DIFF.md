# NUMBER_DIFF — PCOMPBIOL submission → PLOS ONE draft

Every headline/validated number in the manuscript, checked against the code for
this re-assembly. The underlying analyses are unchanged: `reproduce/validate.py`
stays **30/30 green** and the values below trace to tracked artifacts. This file
records (A) the `[VERIFY]` values confirmed against code, (B) numbers newly added
in this draft (the basal-ganglia three-pair replication, Fig 2), (C) numbers
**corrected or flagged** where the frozen draft disagreed with the code, and (D)
presentational renumbering. The `MANUSCRIPT_MD5` gate is decoupled during the
rewrite and re-pinned at the end (independent of this number check).

Rounding note: several abstract/Results values are intentionally reported to
2 significant figures (0.522→"0.52", 0.668→"0.67", 0.346→"0.35"). That rounding
is deliberate and is **not** a changed number.

---

## A. `[VERIFY]` values — confirmed against code (unchanged from submission)

| Manuscript value | Code value | Source artifact | validate.py? |
|---|---|---|---|
| Global config obs/null **0.52** (35 types) | 0.5222 | `output/phase2/scaled_35types/procrustes_results_35.json` | ✔ |
| Global config p **< 10⁻⁶** (1M perms) | permutation floor at 10⁶ | same | ✔ |
| Lineage-stratified obs/null **0.67** | 0.6683 | `output/validation/lineage_stratified/lineage_stratified_results.json` (`stratified_null`) | (Table 1 T02) |
| Lineage-stratified p **< 10⁻⁴** | 9.999e-05 | same | — |
| Mouse-lemur obs/null **0.35** | 0.346 | `analysis/mouse_lemur/procrustes_results.json` | (Table 1 T13) |
| Mouse-lemur p **< 10⁻⁴** | 9.999e-05 | same | — |
| Independent-PCA obs/null 0.473 | 0.473 | `analysis/independent_pca_sensitivity/independent_pca_results.json` | ✔ |
| Primary Layer-2 S **0.483 → 0.230** at k=5 | pre 0.4826, post 0.2295 | `output/mechanistic/ellipsoid_alignment/summary_stats.json` (`35type.mean_alignment.k=5`) | — |
| Layer-2 post-rotation p **< 10⁻⁴** (S=0.230) | 1e-4 | `output/mechanistic/ellipsoid_alignment/permutation_results.json` (T30) | — |
| Sun2023 obs/null 0.554, p<10⁻⁴ | 0.554 | Table 1 T19 | — |
| PanSci obs/null 0.552, p<10⁻⁴ | 0.552 | Table 1 T20 | — |
| CellHint obs/null **0.448**, p<10⁻⁴ | 0.448 | `output/validation/cellhint_replication/cellhint_replication.json` | ✔ |
| Pan-Census obs/null **0.811** (n=22) | 0.811 | `analysis/census_replication/replication_results.json` | (Table 1 T22) |
| Andrews obs/null 0.797, p=0.116 (n=6) | 0.797 | Table 1 T23 | — |
| MCA×HCA obs/null 1.003, p=0.542 (n=17) | 1.003 | Table 1 T24 | — |
| Within-human HCA×TS **0.728, p=0.003** (n=6) | 0.7281 / 0.00300 | `output/validation/hca_centroid_comparison/hca_centroid_comparison.json` | — |
| Macaque obs/null 0.810, raw p=0.0043; corrected p≈0.22 (NS) | 0.810 / 0.0043 / 0.2236 | Table 1 T12; `output/macaque_pipeline/` | — |
| Bootstrap test–retest ρ **0.99** | 0.994 | `analysis/simulation_study/simulation_results.json` (T65) | — |
| Cross-atlas ranking ρ **≈0.15** | Sun 0.146, PanSci 0.194 | Table 1 T25/T26 | — |
| Within/cross inversion ρ (qualitative) | −0.410 | Table 1 T59 | — |
| Simulation recovery ceiling ρ **≈0.42** | 0.42 | `analysis/simulation_study/…` (T64) | — |
| Matched-scale reversal 0.607 (H-H) vs 0.317 (H-M) | 0.607 / 0.317 | `output/phase2/negative_control_v2/negctrl_v2_results.json` (T16/T17) | — |
| Donor-split delta **+0.159**, 100/100 | median +0.1588, pct_positive 1.0 | `analysis/donor_split/donor_split_shared_pca_results.json` | — |
| Donor-split within/cross obs/null 0.375 / 0.527 | 0.3751 / 0.5273 | same | — |
| Ten mechanistic nulls, none significant | T42–T51 all NS | `output/mechanistic/…` | (T42 hk ✔) |
| Conserved-quartile obs/null **0.384** | 0.384 | `analysis/conserved_contribution/gate_results.json` (`secondary.conserved.ratio`) | ✔ |
| Divergent-set 0.709; expr-matched-random 0.525±0.012 | 0.709 / 0.525 | same (`secondary`) | ✔ |
| Master-TF median percentile **0.94** vs **0.54** (expr) / **0.76** (joint) | 0.94 / 0.54 / 0.76 | `analysis/conserved_contribution/gate_results.json`, Methods | ✔ (0.94) |

**All A-values reproduce.** None changed vs the PCOMPBIOL submission.

---

## B. New numbers in this draft — basal-ganglia (HMBA) three-pair Layer-2 replication (Fig 2)

Not present in the PCOMPBIOL submission; added as the paper's positive-result
figure. Source: `~/cellwarp_plans/cellwarp-bg-test/analysis/bg/results/layer2_results_{pair}.json`.
Compression = post-rotation S / pre-rotation S at k=5 (< 1 ⇒ centroid-optimal
alignment compresses the covariance axis). Both weightings shown (W0 unscaled;
W2 = per-gene-standardized Scheme B). Every pre/post permutation p = 9.999e-05
(≈10⁻⁴), and S_post > null at every k tested.

| Pair | n types | k5 compression W0 | k5 compression W2 | W2 rank-1 = canonical identity marker |
|---|---|---|---|---|
| Human–Macaque | 55 | 0.705 | 0.628 | **18 / 55** |
| Human–Marmoset | 52 | 0.627 | 0.594 | 7 / 52 |
| Macaque–Marmoset | 52 | 0.709 | 0.704 | 5 / 52 |

- **All three pairs compress at k=1, 3, 5, both weightings** (post < pre everywhere;
  verified across all 18 cells) — supports "compresses … at every subspace
  dimension tested, each beyond its permutation null (all p ≤ 10⁻⁴)".
- **W0 unscaled → 0 canonical markers** (55/52/52 all "other"); **W2 standardized →
  markers surface** — supports "under per-gene standardization the axis resolves to
  cell-identity markers rather than to housekeeping structure". Main text states
  this **without a fraction** (per CC-WRITE): the classifier is a conservative
  keyword matcher and the 18/55 is a known undercount; the marmoset arms (7, 5) are
  lower still. Fraction stays out of the main text.
- n types **52–55** ⇒ supports "52 to 55 matched cell types".

---

## C. Corrected / flagged — frozen prose vs code

### C1. FLAG (factual — draft yields to code): mouse-lemur divergence time
- **Frozen §1 draft:** "human and mouse lemur, separated by **roughly 90 million years**".
- **Code:** `analysis/mouse_lemur/procrustes_results.json` → `"divergence_mya": 75`;
  PCOMPBIOL submission body also uses **~75 Mya** for mouse lemur.
- **Why the draft is wrong:** ~90 Myr is the human–**mouse** (primate–rodent /
  euarchontoglires) split. Mouse lemur (*Microcebus murinus*) is a strepsirrhine
  **primate**; the human–lemur divergence is ~74–75 Myr (TimeTree), i.e. *shallower*
  than human–mouse. Writing "90 Myr" for the mouse-lemur pair conflates it with the
  rodent split.
- **The `[VERIFY]` note explicitly asked to confirm and make consistent** →
  **resolved to 75 Myr** in the assembled manuscript (`§1` mouse-lemur sentence).
  Human–mouse stays ~90 Myr; Discussion Block 3 "primate–rodent split" (~90 Myr,
  human–mouse) is correct and unchanged. **This is the single place the frozen
  prose was edited on a factual basis — surfaced here, not silently changed.**

### C2. FLAG (minor — awaiting coordinator): cells-per-type range in §2
- **Frozen §2 draft:** "with **thousands to over 100,000** cells per type".
- **Code:** per-type cell floor is **~108** (e.g. `SN-VTR CALB1 Dopa` nH 108 / nM 110;
  `STR Cholinergic GABA` 159; `Monocyte` 242) up to ~145,491 (`Oligo OPALIN`), in
  `layer2_results*.json` `per_type_n`.
- The upper bound ("over 100,000") is correct; the floor ("thousands") understates
  it — several matched types have a few hundred cells. **Suggested wording: "from a
  few hundred to over 100,000 cells per type."** Not yet applied to the frozen
  sentence pending coordinator confirmation (it is not a computed statistic, so no
  validate.py impact either way).

---

## D. Presentational renumbering (no value change)

- **Figures 7 → 5.** New Fig 1 (config conserved) = old Fig 1 panels minus
  bootstrap/LOOCV (→ SI); Fig 2 (two layers + **new BG panel**) = old Fig 4 +
  BG; Fig 3 (config robust) = old Fig 3/replication panels; Fig 4 (per-type not
  resolvable) = old Fig 6 + Fig S3 + simulation; Fig 5 (conserved identity genes)
  = old Fig 2/Fig 7 enrichment panels. Old Fig 2 (conserved-contribution) content
  moves into new Fig 5; old Fig 5 (macaque) demoted to SI.
- **Results sections re-ordered** to findings-first (config → two-layer+replication
  → robustness → per-type limit → identity genes). No statistic changed; only the
  narrative order and figure callouts.
- **Title/Abstract reframed** (measurement-vs-conservation framing). No numbers in
  the title; abstract numbers are the rounded A-values above.
- **"CellWarp" now first appears in Methods** (removed from title/abstract/intro),
  per CC-WRITE ground rule 4.

---

*Bottom line: no validated or headline statistic changed value. One factual
correction (mouse-lemur 90→75 Myr, restoring the submission's value and matching
the code). One minor wording flag (cells-per-type floor) awaiting confirmation.
All new numbers are additive (BG three-pair Layer-2 replication).*
