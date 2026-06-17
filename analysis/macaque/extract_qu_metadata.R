#!/usr/bin/env Rscript
# Extract Qu per-cell metadata from MacFas.meta.data.rds (14.33 GB, Zenodo
# 10.5281/zenodo.5881495) to a lightweight CSV.
#
# The 14.3 GB size suggests this may be a Seurat object carrying metadata
# + reductions + assays, not just a data.frame. First we inspect class and
# structure to determine how to extract per-cell labels. We then write only
# the metadata columns (~20 cols × 174K cells ≈ <100 MB CSV).

suppressPackageStartupMessages({
  library(Matrix)
  # Seurat is needed only if the RDS turns out to be a Seurat object
})

# The Zenodo RDS carries counts + metadata in one object; uncompressed form
# can exceed R's default 24 GB vector heap limit. Raise it (macOS will
# extend into swap if physical RAM is insufficient).
mem.maxVSize(200 * 1024)  # MB → 200 GB

# Resolve repo root portably: $CELLWARP_ROOT, else this script's location, else cwd.
.cw_root <- Sys.getenv("CELLWARP_ROOT", unset = "")
if (.cw_root == "") {
  .m <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  .cw_root <- if (length(.m)) normalizePath(file.path(dirname(sub("^--file=", "", .m)), "..", "..")) else getwd()
}
IN_RDS  <- file.path(.cw_root, "data/macaque/qu_2022/MacFas.meta.data.rds")
OUT_CSV <- file.path(.cw_root, "data/macaque/qu_2022/qu_metadata.csv")
OUT_SUM <- file.path(.cw_root, "data/macaque/qu_2022/qu_metadata_summary.txt")

cat("Loading RDS:", IN_RDS, "\n")
t0 <- Sys.time()
obj <- readRDS(IN_RDS)
t1 <- Sys.time()
cat(sprintf("Loaded in %.1f minutes\n", as.numeric(difftime(t1, t0, units="mins"))))
cat("Class:", class(obj), "\n")

md <- NULL
if (is.data.frame(obj)) {
  cat("Object is a data.frame.\n")
  md <- obj
} else if (isS4(obj)) {
  cat("Object is S4. Slot names:\n"); print(slotNames(obj))
  # Direct @meta.data access avoids loading full Seurat package for class methods.
  if ("meta.data" %in% slotNames(obj)) {
    md <- obj@meta.data
    cat("Extracted @meta.data:", nrow(md), "rows x", ncol(md), "cols\n")
  }
} else if (is.list(obj)) {
  cat("Object is a list. Names:\n"); print(names(obj))
  if ("meta.data" %in% names(obj)) md <- obj[["meta.data"]]
}

if (is.null(md)) {
  stop("Could not find metadata in the loaded object.")
}

cat("\nMetadata dim:", dim(md), "\n")
cat("Columns:\n"); print(colnames(md))

# Write summary report
sink(OUT_SUM)
cat("Qu Metadata Summary\n==================\n\n")
cat("Source:", IN_RDS, "\n")
cat("Rows (cells):", nrow(md), "\n")
cat("Columns:", ncol(md), "\n\n")
cat("Column inventory:\n")
for (c in colnames(md)) {
  vals <- md[[c]]
  cls <- class(vals)[1]
  n_unique <- length(unique(vals))
  na_count <- sum(is.na(vals))
  cat(sprintf("  %-30s class=%-10s unique=%-8d NAs=%d\n", c, cls, n_unique, na_count))
  if (n_unique <= 80 && cls %in% c("character","factor","ordered","logical","integer")) {
    tab <- sort(table(vals, useNA="ifany"), decreasing=TRUE)
    for (nm in names(tab)[1:min(80,length(tab))]) {
      cat(sprintf("      %-50s %d\n", nm, tab[[nm]]))
    }
  } else if (cls %in% c("character","factor","ordered") && n_unique <= 300) {
    tab <- sort(table(vals, useNA="ifany"), decreasing=TRUE)
    cat(sprintf("      (top 20 of %d values)\n", n_unique))
    for (nm in names(tab)[1:20]) {
      cat(sprintf("      %-50s %d\n", nm, tab[[nm]]))
    }
  }
}
sink()

# Preserve row names as cell_id
if (!is.null(rownames(md))) {
  md <- cbind(cell_id=rownames(md), md)
}
write.csv(md, OUT_CSV, row.names=FALSE)
cat(sprintf("\nCSV written: %s (%.1f MB)\n", OUT_CSV, file.info(OUT_CSV)$size/1e6))
cat(sprintf("Summary written: %s\n", OUT_SUM))
