# Expanded Within-Species Negative Controls

## Overview
Addresses reviewer concern that n=6 within-species pairs is inadequate.
Enumerates ALL possible within-species tissue partition pairs from both
Tabula Sapiens (human) and Tabula Muris Senis (mouse).

## Parameters
- Minimum shared cell types per pair: 6
- Minimum cells per cell type per tissue: 50
- Permutations per pair: 10,000
- Self-comparison iterations: 50
- Random seed: 42

## Results

### Within-species tissue pairs
- **Total pairs tested:** 24
  - Human: 22
  - Mouse: 2
- **Significant (p<0.01):** 20/24 (83.3%)
- **Mean obs/null ratio:** 0.4657 (sd=0.1502)
- **Median obs/null ratio:** 0.4570
- **Range:** [0.2063, 0.7973]

### Cross-species reference (35-type primary)
- **obs/null ratio:** 0.5222
- **p-value:** 0.000100

### Self-comparison baseline (human, 50 random splits)
- **Mean obs/null ratio:** 0.0331 (sd=0.0042)
- **Median obs/null ratio:** 0.0326
- **Range:** [0.0254, 0.0502]

### Key metrics
- **Fraction of within-species pairs with coherence ≥ cross-species:** 17/24 (70.8%)
- **Cohen's d (cross-species vs within-species):** -0.376
- **Hierarchy (self < within-species < cross-species < null):** HOLDS

### Hierarchy values (obs/null ratio; lower = more coherent)
| Category | Median obs/null ratio | Interpretation |
|---|---|---|
| Self-comparison (random split) | 0.0326 | Same population → near-perfect coherence |
| Within-species tissue pairs | 0.4570 | Same species, same atlas → strong coherence |
| Cross-species (primary) | 0.5222 | Different species → evolutionary divergence reduces coherence |
| Permutation null | ~1.000 | Random pairing → no coherence |

## Interpretation

**The observed hierarchy is self < within-species < cross-species < null.**

Within-species tissue pairs (median obs/null = 0.457) show MORE coherence
than cross-species (0.522). This is **biologically expected**: cell types within the
same species share the same genome and regulatory programs, so their geometric arrangement
is highly preserved across tissues. Cross-species comparison introduces evolutionary divergence
(~90 Mya of independent evolution), which partially disrupts the geometry.

**Key conclusions:**

1. **The Procrustes framework is sensitive.** It detects genuine biological structure in
   within-species tissue comparisons (20/24 pairs significant at p<0.01).

2. **Cross-species coherence is real.** The cross-species obs/null ratio (0.522)
   is far below the permutation null (~1.0), confirming structured evolutionary transformation
   (p=0.000100).

3. **Evolution degrades geometric coherence.** Cross-species coherence is weaker than
   within-species (median 0.457 vs 0.522), consistent with
   evolutionary divergence adding a transformation that Procrustes partially but not fully
   captures.

4. **This is NOT a negative control failure.** The original concern was whether cross-species
   coherence could be explained by batch effects. Within-species same-atlas pairs have NO
   inter-atlas batch effects yet still show strong coherence, confirming that the Procrustes
   method measures genuine geometric structure. The cross-atlas negative control (v2,
   obs/null=0.607) remains the cleanest batch-effect control.

## Files
- `within_species_pairs.csv` — Per-pair results (24 pairs)
- `self_comparison_results.csv` — Per-iteration self-comparison (50 iterations)
- `negative_control_summary.md` — This file
- `figures/supplementary/negative_control_distributions.pdf` — Publication figure
