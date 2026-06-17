# Mantel Test: Cross-Species Distance Matrix Correlation

## Results

| Metric   | r      | p-value |
|----------|--------|---------|
| Pearson  | 0.7890 | 0.000100 |
| Spearman | 0.7368 | 0.000100 |

- Cell types: 35
- Permutations: 10,000

## Interpretation

The Mantel test evaluated whether the 35 x 35 pairwise Euclidean distance matrix among cell-type centroids in joint PCA space is correlated between human and mouse. The Pearson correlation was 0.7890 (p < 0.001, 10,000 permutations) and the Spearman rank correlation was 0.7368 (p = 0.000100), indicating a strong and statistically significant preservation of pairwise distance structure across species. This result supports the assumption underlying the Procrustes analysis: the geometric relationships among cell types in transcriptomic space are conserved between human and mouse, even under a weaker test that does not assume a single rigid-body transformation. Because the Mantel test only requires monotonic distance preservation (especially via the Spearman variant) rather than exact linear correspondence, this positive result provides independent evidence that the cross-species cell-type configuration is not an artefact of the Procrustes fitting procedure but reflects genuine biological conservation of the transcriptomic geometry.
