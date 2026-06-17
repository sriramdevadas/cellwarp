#!/usr/bin/env Rscript
# Phase 3 blinded dry-run — fgsea worker
#
# Usage:
#   Rscript scripts/phase3_fgsea_worker.R <rnk_csv> <gmt_path> <nperm> <seed> <min_size> <max_size> <out_csv>
#
# - rnk_csv: CSV with columns (gene, score) — per-axis GSEA ranking vector
# - gmt_path: GMT catalog file
# - nperm: permutations (10000 per §3.5)
# - seed: fgsea/sampling seed (per-axis deterministic; derived in caller)
# - min_size / max_size: size filter per §3.5 (20/500 for GO-BP & Reactome; 1/Inf for Hallmark)
# - out_csv: output CSV with columns pathway, pval, padj, ES, NES, size
#
# Runs fgsea against a pre-ranked gene list and writes results to a CSV.
# §3.5 protocol: 10,000 permutations, BH FDR within-axis (per run). No
# cross-axis aggregation here — caller stitches axis-level CSVs.

suppressPackageStartupMessages({
  library(fgsea)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 7)
rnk_csv  <- args[[1]]
gmt_path <- args[[2]]
nperm    <- as.integer(args[[3]])
seed     <- as.integer(args[[4]])
min_size <- as.integer(args[[5]])
max_size <- as.integer(args[[6]])
out_csv  <- args[[7]]

set.seed(seed)

rnk <- fread(rnk_csv)
stopifnot(all(c("gene", "score") %in% names(rnk)))
# Deduplicate by gene — if a gene appears multiple times, keep the max-abs score
# (shouldn't happen under §3.4 scoring since gene scores are per-gene, but defensive).
rnk <- rnk[, .(score = score[which.max(abs(score))]), by = gene]
ranks <- setNames(rnk$score, rnk$gene)
# fgsea requires non-NA numeric; order is not required (fgsea sorts internally).
ranks <- ranks[!is.na(ranks)]

pathways <- gmtPathways(gmt_path)

res <- fgsea(
  pathways = pathways,
  stats    = ranks,
  minSize  = min_size,
  maxSize  = max_size,
  nPermSimple = nperm
)

# Drop leadingEdge (list column — we don't need it for decision tree) for CSV.
res[, leadingEdge := NULL]
fwrite(res, out_csv)
