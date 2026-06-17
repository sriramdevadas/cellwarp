CellWarp Aging Axis — Pre-Registration
Date: 2026-03-16
Author: [you]
Status: Committed before expression matrix analysis begins

Primary Hypothesis
The primary geometric signature of aging in mouse somatic cells is a directional centroid shift in transcriptomic space — not symmetric cluster expansion. Biological aging represents configuration drift in a complex system: the epigenetic controls maintaining a highly specialized cellular state degrade over time, causing 18-month cells to drift away from their 3-month transcriptomic coordinates toward a more generic, less specialized state. This drift will produce a measurable, consistent displacement vector for a given cell type — moving the center of gravity rather than inflating the noise.

Primary Outcome Measure
For each candidate cell type, compute Δcentroid = centroid(18m) − centroid(3m) in 1:1 ortholog gene space.
A positive result requires:

Δcentroid magnitude exceeds inter-donor variation in both age groups (young-young and old-old pairwise variation as bracketing estimates)
Result holds in ≥2 tissues for that cell type (cross-tissue consistency)


What Would Change My Mind
Two falsification conditions, either sufficient:

Symmetric scatter: 18-month cells show massive shape distortion without a unified directional shift — Δcentroid magnitude is small relative to cluster expansion. This would support H2 (diffuse activation) over H1.
Swamped by donor noise: Δcentroid magnitude between 3m and 18m is indistinguishable from old-old inter-donor variation. Identity shift cannot be separated from individual biological noise at n=4 old donors.


Gene Space
1:1 human-mouse orthologs (same space as core CellWarp analysis). Rationale: enables direct geometric comparison of the aging displacement vector against the cross-species evolutionary transformation. Primary downstream question: does a cell age in the same transcriptomic direction it evolved?

Candidate Cell Types
In priority order: satellite cell, endothelial cell, basal cell of epidermis, cardiac fibroblast, B cell. All LOW or MEDIUM compositional risk. HIGH-risk types (macrophage, myeloid, T cell subtypes) excluded from primary analysis.

Secondary Measurements
Regardless of H1 outcome, record:

Shape distortion magnitude per cell type (H2 signal)
LOOCV stability in old vs young (H3 signal)
Direction of Δcentroid relative to cross-species Procrustes transformation (the evolutionary alignment question)

The last point is exploratory and not a success/failure criterion — it's the question this pre-registration is designed to make answerable.
