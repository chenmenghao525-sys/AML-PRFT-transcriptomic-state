#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

options(stringsAsFactors = FALSE)

dir.create("04_prft_score", recursive = TRUE, showWarnings = FALSE)
dir.create("14_tables", recursive = TRUE, showWarnings = FALSE)
dir.create("16_logs", recursive = TRUE, showWarnings = FALSE)

save_session_info <- function(path) {
  writeLines(capture.output(sessionInfo()), con = path)
}

zscore_safe <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) {
    return(rep(0, length(x)))
  }
  as.numeric(scale(x))
}

custom_rank_ssgsea_like <- function(expr_mat, gene_sets) {
  expr_mat <- as.matrix(expr_mat)
  if (!is.numeric(expr_mat)) {
    storage.mode(expr_mat) <- "numeric"
  }

  sample_ranks <- apply(expr_mat, 2, function(x) rank(x, ties.method = "average", na.last = "keep"))
  if (is.null(dim(sample_ranks))) {
    sample_ranks <- matrix(sample_ranks, ncol = 1)
    rownames(sample_ranks) <- rownames(expr_mat)
    colnames(sample_ranks) <- colnames(expr_mat)
  }

  n_genes <- nrow(expr_mat)
  score_mat <- matrix(NA_real_, nrow = length(gene_sets), ncol = ncol(expr_mat))
  rownames(score_mat) <- names(gene_sets)
  colnames(score_mat) <- colnames(expr_mat)

  for (i in seq_along(gene_sets)) {
    gs_name <- names(gene_sets)[i]
    available_genes <- intersect(unique(gene_sets[[i]]), rownames(expr_mat))

    if (length(available_genes) == 0) {
      warning("Gene set ", gs_name, " has 0 available genes. Returning NA scores.")
      next
    }

    mean_ranks <- colMeans(sample_ranks[available_genes, , drop = FALSE], na.rm = TRUE)
    normalized_scores <- (mean_ranks - 1) / max(n_genes - 1, 1)
    score_mat[i, ] <- zscore_safe(normalized_scores)
  }

  score_mat
}

obj <- readRDS("02_processed_data/tcga_expr_clin_matched.rds")
gene_sets <- readRDS("03_gene_sets/prft_gene_sets_all.rds")

expr <- obj$expr
if (is.null(expr) || !is.matrix(expr)) {
  expr <- as.matrix(expr)
}

if (is.null(rownames(expr)) || is.null(colnames(expr))) {
  stop("Expression matrix must have gene symbols as rownames and sample IDs as colnames.")
}

if (is.null(names(gene_sets)) || any(names(gene_sets) == "")) {
  stop("Gene sets must be a named list.")
}

coverage_dt <- rbindlist(lapply(names(gene_sets), function(nm) {
  gs <- unique(as.character(gene_sets[[nm]]))
  available <- intersect(gs, rownames(expr))
  missing <- setdiff(gs, rownames(expr))
  if (length(available) < 5) {
    warning("Gene set ", nm, " has fewer than 5 available genes in the expression matrix.")
  }
  data.table(
    gene_set_name = nm,
    total_genes = length(gs),
    available_genes = length(available),
    missing_genes = length(missing),
    available_gene_list = paste(available, collapse = ";"),
    missing_gene_list = paste(missing, collapse = ";")
  )
}))

fwrite(coverage_dt, "14_tables/tcga_gene_set_coverage.csv")

scoring_method <- NA_character_
scores <- NULL

gsva_available <- requireNamespace("GSVA", quietly = TRUE)

if (gsva_available) {
  scores <- tryCatch(
    {
      GSVA::gsva(
        expr = expr,
        gset.idx.list = gene_sets,
        method = "ssgsea",
        kcdf = "Gaussian",
        abs.ranking = TRUE
      )
    },
    error = function(e1) {
      message("GSVA::gsva(method = 'ssgsea') failed: ", conditionMessage(e1))
      tryCatch(
        {
          param <- GSVA::ssgseaParam(exprData = expr, geneSets = gene_sets)
          GSVA::gsva(param)
        },
        error = function(e2) {
          message("GSVA::ssgseaParam workflow also failed: ", conditionMessage(e2))
          NULL
        }
      )
    }
  )

  if (!is.null(scores)) {
    scoring_method <- "GSVA_ssgsea"
  }
}

if (is.null(scores)) {
  message("Using fallback scoring method: custom_rank_ssgsea_like")
  scores <- custom_rank_ssgsea_like(expr_mat = expr, gene_sets = gene_sets)
  scoring_method <- "custom_rank_ssgsea_like"
}

scores <- as.matrix(scores)
saveRDS(scores, "04_prft_score/tcga_ssgsea_scores.rds")
fwrite(as.data.table(scores, keep.rownames = "gene_set_name"), "04_prft_score/tcga_ssgsea_scores.csv")

method_dt <- data.table(
  scoring_method = scoring_method,
  expression_input = "02_processed_data/tcga_expr_clin_matched.rds::obj$expr",
  score_matrix_rows = nrow(scores),
  score_matrix_cols = ncol(scores)
)
fwrite(method_dt, "14_tables/tcga_ssgsea_method_used.csv")

save_session_info("16_logs/sessionInfo_08_calculate_ssgsea_scores.txt")

