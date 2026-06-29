#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
ascii_lib <- Sys.getenv("PHASE1_ASCII_R_LIB", unset = "")
if (nzchar(ascii_lib) && dir.exists(ascii_lib)) .libPaths(c(ascii_lib, .libPaths()))

root <- Sys.getenv("PHASE1_AUDIT_ROOT", unset = "")
if (!nzchar(root)) root <- getwd()
root <- chartr("\\", "/", root)

show_obj <- function(rel) {
  cat("\nFILE ", rel, "\n", sep = "")
  x <- readRDS(file.path(root, rel))
  print(class(x))
  if (is.matrix(x) || is.data.frame(x)) print(dim(x))
  if (is.list(x)) {
    print(names(x))
    print(utils::head(lapply(x, function(y) {
      if (is.matrix(y) || is.data.frame(y)) dim(y) else if (is.vector(y)) length(y) else class(y)
    }), 20))
  }
}

for (rel in c(
  "phase1_runtime/02_processed_data/tcga_expr_hgnc_log2cpm.rds",
  "phase1_runtime/02_processed_data/tcga_expr_clin_matched.rds",
  "phase1_runtime/04_prft_score/tcga_prft_score.rds",
  "phase1_runtime/03_gene_sets/prft_gene_sets_all.rds",
  "phase1_runtime/03_gene_sets/prft_gene_sets_main.rds",
  "phase1_runtime/03_gene_sets/prft_gene_sets_supplementary.rds"
)) {
  show_obj(rel)
}
