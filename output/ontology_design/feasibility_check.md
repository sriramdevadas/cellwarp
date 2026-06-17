# Feasibility Check — 100-Type Ontology Expansion

**Date:** 2026-03-16
**Purpose:** Assess feasibility without downloading expression data.

---

## 1. Gene Space Impact

**Question:** Will adding 65 new cell types change the shared ortholog space?

**Answer: NO.** The ortholog space is determined by species-level gene mapping, not by
cell type selection. The current shared ortholog space:
- Human-mouse 1:1 orthologs: 17,187 genes
- After intersection with both atlases: 16,959 genes
- Human-macaque 1:1 orthologs: 19,123 genes (DECISION-098)
- Three-way intersection (human-mouse-macaque): 15,028 genes

Adding new cell types queries the same gene space — PCA is performed on the full
ortholog gene matrix. New cell types add rows (centroids) to the same columns (genes).
The only scenario where gene space could change is if a supplementary atlas uses a
different gene annotation version, but this is handled at the data loading stage (Ensembl
gene ID harmonization), not at the cell type selection stage.

**Confirmation: PASS.** No gene space impact from expansion.

---

## 2. Computational Scaling

**Question:** How does compute time scale from n=35 to n=100?

### Current n=35 benchmarks (from Phase 2 runs)
| Step | Scaling | Est. n=35 time | Est. n=100 time |
|---|---|---|---|
| Centroid computation | O(cells × genes) per type | ~2 min | ~6 min (linear with total cells) |
| PCA | O(n × d²) where d=genes | ~1 min | ~1 min (d unchanged, n triples) |
| Procrustes alignment | O(n × k) where k=PCs | <1 sec | <1 sec (still tiny matrix) |
| Permutation test (10,000 perms) | O(perms × n × k) | ~30 sec | ~90 sec (linear in n) |
| LOOCV (leave-one-type-out) | O(n × Procrustes) | ~35 × <1 sec | ~100 × <1 sec |
| Residual decomposition | O(n × k) | <1 sec | <1 sec |
| Bootstrap (100 subsamples) | O(samples × pipeline) | ~5 min | ~15 min |

### Total estimated compute time
- **n=35 full pipeline:** ~10 minutes
- **n=100 full pipeline:** ~25 minutes
- **Scaling factor:** ~2.5× (approximately linear in n for most steps, sub-linear for
  PCA since gene dimension dominates)

### Data download scaling
The real scaling bottleneck is **data download**, not computation:
- Current 35 types: ~96,600 cells total (52,434 human + 44,175 mouse)
- Estimated 100 types: ~300,000-500,000 cells total (depends on new type sizes)
- Download time: ~30-60 min from Census (network-bound)
- Memory: ~8-16 GB for full AnnData objects (within M-series MacBook capacity)

**Assessment: FEASIBLE.** Computational scaling is linear and manageable. The bottleneck
is data download and QC, not Procrustes computation.

---

## 3. Tabula Sapiens Gaps

**Question:** Which of the 65 new types are absent from Tabula Sapiens?

### Types requiring supplementary human atlas

| Category | Count | Types | Recommended Atlas |
|---|---|---|---|
| **Neural (TS has no brain)** | **12** | oligodendrocyte, astrocyte, microglial cell, OPC, excitatory neuron, inhibitory neuron, Schwann cell, ependymal cell, Purkinje neuron, medium spiny neuron, dopaminergic neuron, cerebellar granule cell | **Allen Brain Cell Atlas** (Siletti et al. 2023, Science). 3M+ human brain cells. 10x Chromium. In CELLxGENE Census. |
| **Endocrine (TS has no thyroid/adrenal)** | **5** | thyroid follicular cell, adrenocortical cell, chromaffin cell, parathyroid cell + pancreatic beta/alpha (low count) | **Human Pancreas Analysis Program** (HPAP) for islet cells. Census endocrine datasets for others. |
| **Gastric (TS has no stomach)** | **2** | gastric chief cell, parietal cell | Census gastric datasets (limited availability). |
| **Vascular/neural overlap** | **1** | brain endothelial cell | Allen Brain (endothelial captured in brain dissociations). |

**Total: 20 of 65 new types require a supplementary human atlas.**

The Allen Brain Cell Atlas alone would resolve 13 of these 20 gaps (all neural + brain
endothelial). This is the single most impactful supplementary dataset.

### Types with low TS counts (require Census pooling)
| Type | TS Count | Census Pooling Feasible? |
|---|---|---|
| plasmacytoid dendritic cell | 99 | YES (blood datasets) |
| Langerhans cell | 76 | UNCERTAIN (skin datasets limited) |
| enteroendocrine cell | 25 | UNCERTAIN (GI datasets limited) |
| pancreatic delta cell | 9 | NO (too rare even pooled) |

---

## 4. Tabula Muris Senis Gaps

**Question:** Which of the 65 new types are absent or under 500 in TMS?

### Types requiring Census mouse pooling or supplementary mouse atlas

| Type | TMS Count | Census Mouse Feasible? | Risk Level |
|---|---|---|---|
| regulatory T cell | 170 | YES (Immunological Genome datasets) | LOW |
| mast cell | 76 | LIKELY (skin/lung datasets) | LOW |
| basophil | 325 | LIKELY (bone marrow datasets) | LOW |
| plasmacytoid DC | 122 | LIKELY (spleen datasets) | MEDIUM |
| erythrocyte | 356 | LIKELY (blood datasets) | LOW |
| pericyte | 99 | LIKELY (vascular/heart datasets) | MEDIUM |
| club cell | 37 | LIKELY (lung datasets) | MEDIUM |
| type I pneumocyte | 4 | **HIGH RISK** (fragile cells, poor capture) | HIGH |
| lymphatic endothelial | 140 | UNCERTAIN (vascular datasets) | MEDIUM |
| hepatic stellate cell | 39 | **HIGH RISK** (rare in scRNA-seq) | HIGH |

### Technology compatibility risk

**Critical lesson from MCA failure (DECISION-104):** The CellWarp pipeline requires
10x Chromium or comparable droplet technology. Supplementary mouse atlases MUST be
pre-screened for:
1. **Technology:** 10x Chromium 3' v2 or v3 only. Reject microwell-seq, sci-RNA-seq,
   Smart-seq2 (different gene detection characteristics).
2. **Gene detection:** Median genes/cell ≥2,000.
3. **Cell count:** ≥500 cells of the target type after QC filtering.

The Tabula Muris Senis FACS arm uses Smart-seq2 (plate-based, high-sensitivity) which
we have successfully used. Census pooling of 10x datasets from other labs introduces
batch effects that must be evaluated via the same negative control framework used for
Sun2023 and PanSci replications.

### Types likely NOT achievable even with Census pooling
| Type | Reason |
|---|---|
| Type I pneumocyte | Fragile, lysed during dissociation. Only 4 in TMS. |
| Hepatic stellate cell | Activated during dissociation, loses markers. 39 in TMS. |
| Pancreatic delta cell | Rare endocrine type. 463 in TMS but 9 in TS. |
| Parathyroid cell | No parathyroid tissue in standard atlases. |
| Enteroendocrine cell | 161 in TMS, 25 in TS. Rare even pooled. |

**Assessment: 5 types at HIGH RISK of not reaching 500-cell threshold.** These are
candidates for deferral or replacement. The remaining 60 new types are LIKELY achievable
with Census pooling + Allen Brain supplementation.

---

## 5. Summary Assessment

| Dimension | Status | Notes |
|---|---|---|
| Gene space | **PASS** | No impact — ortholog space is species-level |
| Compute scaling | **PASS** | ~2.5× current time, ~25 min total |
| Memory | **PASS** | ~16 GB peak, within M-series capacity |
| Allen Brain dependency | **REQUIRED** | Resolves 13/20 TS gaps (all neural + brain vascular) |
| Census pooling dependency | **REQUIRED** | Resolves ~30 types with cross-dataset pooling |
| High-risk types | **5 types** | AT1, hepatic stellate, delta cell, parathyroid, enteroendocrine |
| Realistic achievable ceiling | **~93-95 types** | 100 minus 5-7 high-risk types |

### Recommendation
Proceed with 100-type design but designate 5 high-risk types as PROVISIONAL. If
empirical Census queries confirm <500 cells for any provisional type, replace from
the candidate pool (105 candidates enumerated in candidate_pool.csv). The design
is robust to losing up to 7 types without compromising organ system coverage.
