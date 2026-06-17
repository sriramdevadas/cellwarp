# Three-Way Cell Type Intersection Report

**Date:** 2026-04-05
**Status:** NO-GO (8 types < 10-type threshold)

## Intersection (8 types)

These cell types are present in ALL THREE species-pair analyses
(human-mouse 35 types, human-macaque 20 types, human-mouse lemur 15 types):

1. B cell
2. CD4-positive, alpha-beta T cell
3. CD8-positive, alpha-beta T cell
4. T cell
5. endothelial cell
6. fibroblast
7. macrophage
8. natural killer cell

## Bottleneck

The macaque analysis (20 types) is the limiting factor. The following types
are present in BOTH human-mouse and human-mouse lemur but ABSENT from
human-macaque:

| Type | In Mouse | In Lemur | Macaque Issue |
|------|:--------:|:--------:|---------------|
| monocyte | YES | YES | Only subtypes: classical, intermediate, non-classical |
| neutrophil | YES | YES | Only parent: granulocyte |
| mature NK T cell | YES | YES | Absent |
| mesenchymal stem cell | YES | YES | Absent |
| plasma cell | YES | YES | Absent |
| enterocyte of epithelium of large intestine | YES | YES | Absent |
| pancreatic acinar cell | YES | YES | Absent |

## Potential Rescue via Controlled Aliasing

If the user permits aliasing at matching annotation granularity:
- **monocyte ↔ classical monocyte**: Could be justified (classical monocytes are
  the dominant subset). Would add 1 type → 9 types.
- **neutrophil ↔ granulocyte**: Neutrophils are ~60-70% of granulocytes.
  Biologically close but not identical. Would add 1 type → 10 types.

With BOTH aliases: 10 types (meets threshold).

**Risk:** These aliases introduce systematic bias — the macaque "granulocyte"
centroid includes eosinophils/basophils that are excluded from the mouse/lemur
"neutrophil" centroid. The "classical monocyte" centroid excludes intermediate
and non-classical monocytes that are included in the lemur "monocyte" centroid.
This is a judgment call for the user.

## Source Data

| Analysis | n_types | n_genes | obs/null | p-value |
|----------|---------|---------|----------|---------|
| Human-mouse | 35 | 16,959 | 0.522 | <10^-6 |
| Human-macaque | 20 | 13,927 | 0.841 | 0.0002 |
| Human-mouse lemur | 15 | 13,796 | 0.346 | 0.0001 |
