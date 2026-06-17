#!/usr/bin/env Rscript
# Extract RIRA cell-level metadata from RIRA.All.Metadata.rds to a CSV.
#
# Biology: RIRA (Robust Integrated Rhesus Atlas, Mahyari et al. 2025) is a
# 47-donor M. mulatta multi-tissue immune atlas. The metadata RDS carries
# per-cell annotations (cell type, donor, tissue) that our Python pipeline
# needs to compute per-type centroids. We extract once to CSV so downstream
# reproducibility does not require an R runtime.
#
# Inputs:
#   data/macaque/rira/RIRA.All.Metadata.rds   (55 MB; ~596K rows)
# Outputs:
#   data/macaque/rira/rira_metadata.csv       (~50 MB; per-cell labels)
#   data/macaque/rira/rira_metadata_summary.txt (column inventory)

suppressPackageStartupMessages({
  library(utils)
})

# Resolve repo root portably: $CELLWARP_ROOT, else this script's location, else cwd.
.cw_root <- Sys.getenv("CELLWARP_ROOT", unset = "")
if (.cw_root == "") {
  .m <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  .cw_root <- if (length(.m)) normalizePath(file.path(dirname(sub("^--file=", "", .m)), "..", "..")) else getwd()
}
IN_RDS  <- file.path(.cw_root, "data/macaque/rira/RIRA.All.Metadata.rds")
OUT_CSV <- file.path(.cw_root, "data/macaque/rira/rira_metadata.csv")
OUT_SUM <- file.path(.cw_root, "data/macaque/rira/rira_metadata_summary.txt")

cat("Loading RDS:", IN_RDS, "\n")
md <- readRDS(IN_RDS)
cat("Class:", class(md), "\n")
cat("Dim:", paste(dim(md), collapse=" x "), "\n")
cat("Columns:\n"); print(colnames(md))

# Write a short summary report
sink(OUT_SUM)
cat("RIRA Metadata Summary\n")
cat("=====================\n\n")
cat("Source:", IN_RDS, "\n")
cat("Class:", class(md)[1], "\n")
cat("Rows (cells):", nrow(md), "\n")
cat("Columns:", ncol(md), "\n\n")
cat("Column inventory:\n")
for (c in colnames(md)) {
  vals <- md[[c]]
  cls <- class(vals)[1]
  n_unique <- length(unique(vals))
  na_count <- sum(is.na(vals))
  cat(sprintf("  %-30s class=%-10s unique=%-8d NAs=%d\n",
              c, cls, n_unique, na_count))
  if (n_unique <= 60 && cls %in% c("character","factor","logical","integer")) {
    tab <- sort(table(vals), decreasing=TRUE)
    for (nm in names(tab)[1:min(60,length(tab))]) {
      cat(sprintf("      %-40s %d\n", nm, tab[[nm]]))
    }
  } else if (cls %in% c("character","factor") && n_unique <= 200) {
    tab <- sort(table(vals), decreasing=TRUE)
    cat("      (top 15 of", n_unique, "values)\n")
    for (nm in names(tab)[1:15]) {
      cat(sprintf("      %-40s %d\n", nm, tab[[nm]]))
    }
  }
}
sink()
cat("\nSummary written to:", OUT_SUM, "\n")

# Write full metadata as CSV
if (is.data.frame(md)) {
  df <- md
} else {
  df <- as.data.frame(md)
}

# Preserve row names (cell barcodes) as first column
if (!is.null(rownames(df))) {
  df <- cbind(cell_id=rownames(df), df)
}

write.csv(df, OUT_CSV, row.names=FALSE)
cat("Full metadata CSV written to:", OUT_CSV, "\n")
cat("File size:", format(file.info(OUT_CSV)$size/1e6, digits=3), "MB\n")
