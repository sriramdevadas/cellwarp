#!/usr/bin/env Rscript
# Extract MSigDB Hallmark, GO-BP, and Reactome catalogs from msigdbr
# as GMT files for fgsea. Pinned version: MSigDB 2026.1.Hs via msigdbr 26.1.0.
#
# Biology: gene set catalogs define grouped gene memberships for pre-ranked GSEA.
# Hallmark (50 curated sets) is the §4 decision-tree primary; GO-BP and Reactome
# (size-filtered 20-500) are §6.5 corroboration catalogs.
#
# Math: GMT format is one tab-separated line per set: set_id \t description \t gene1 \t gene2 ...
# Size filter applies to gene-set cardinality; fgsea rejects sets outside [minSize, maxSize].

suppressPackageStartupMessages({
  library(msigdbr)
})

out_dir <- "data/phase3/catalogs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

db_version <- "2026.1.Hs"
species <- "Homo sapiens"

write_gmt <- function(df, key_col, path) {
  # GMT: tab-separated, header per set: set_id, description, gene1, gene2, ...
  con <- file(path, "w")
  on.exit(close(con))
  split_by <- split(df$gene_symbol, df[[key_col]])
  # carry a descriptor per set (set name); we index by term_id for corroboration
  desc_lookup <- tapply(df$gs_name, df[[key_col]], function(x) x[1])
  for (sid in names(split_by)) {
    genes <- unique(split_by[[sid]])
    # drop empty / NA genes
    genes <- genes[!is.na(genes) & nchar(genes) > 0]
    line <- paste(c(sid, as.character(desc_lookup[[sid]]), genes), collapse = "\t")
    writeLines(line, con)
  }
}

cat("msigdbr version:", as.character(packageVersion("msigdbr")), "\n")
cat("MSigDB db_version:", db_version, "\n\n")

# --- Hallmark: key by gs_name (H collection has 50 sets identified by name) ---
cat("== Hallmark ==\n")
h <- msigdbr(species = species, collection = "H")
cat("rows:", nrow(h), "unique sets:", length(unique(h$gs_name)), "\n")
# Use gs_name as the identifier (these are "HALLMARK_*" names)
h$term_id <- h$gs_name
hallmark_path <- file.path(out_dir, sprintf("msigdb_hallmark_v%s.gmt", db_version))
write_gmt(h, "term_id", hallmark_path)
cat("wrote:", hallmark_path, "\n\n")

# --- GO-BP: key by gs_exact_source (GO ID), description = gs_name ---
cat("== GO-BP ==\n")
gobp <- msigdbr(species = species, collection = "C5", subcollection = "GO:BP")
cat("rows:", nrow(gobp), "unique sets by gs_name:", length(unique(gobp$gs_name)), "\n")
cat("unique by gs_exact_source:", length(unique(gobp$gs_exact_source)), "\n")
# gs_exact_source is the GO ID
gobp$term_id <- gobp$gs_exact_source
# drop any rows where term_id is NA/empty
gobp <- gobp[!is.na(gobp$term_id) & nchar(gobp$term_id) > 0, ]
gobp_path <- file.path(out_dir, sprintf("msigdb_c5_gobp_v%s.gmt", db_version))
write_gmt(gobp, "term_id", gobp_path)
cat("wrote:", gobp_path, "\n\n")

# --- Reactome: key by gs_exact_source (R-HSA ID) ---
cat("== Reactome ==\n")
rxm <- msigdbr(species = species, collection = "C2", subcollection = "CP:REACTOME")
cat("rows:", nrow(rxm), "unique sets by gs_name:", length(unique(rxm$gs_name)), "\n")
cat("unique by gs_exact_source:", length(unique(rxm$gs_exact_source)), "\n")
rxm$term_id <- rxm$gs_exact_source
rxm <- rxm[!is.na(rxm$term_id) & nchar(rxm$term_id) > 0, ]
reactome_path <- file.path(out_dir, sprintf("msigdb_c2_reactome_v%s.gmt", db_version))
write_gmt(rxm, "term_id", reactome_path)
cat("wrote:", reactome_path, "\n\n")

cat("Catalog extraction complete.\n")
