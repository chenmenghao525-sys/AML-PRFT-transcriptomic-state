#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

ascii_lib <- Sys.getenv("PHASE1_ASCII_R_LIB", unset = "")
if (nzchar(ascii_lib) && dir.exists(ascii_lib)) {
  .libPaths(c(ascii_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

root <- Sys.getenv("PHASE1_AUDIT_ROOT", unset = "")
if (!nzchar(root)) root <- getwd()
root <- chartr("\\", "/", root)

input_obj <- readRDS(file.path(root, "phase1_runtime", "07_signature", "cross_platform_lasso_input_matrix_tcga.rds"))
six <- fread(file.path(root, "03_results_tables", "phase1_six_gene_coefficients.csv"))

expr <- input_obj$expr
clin <- as.data.table(input_obj$clin)
expr <- expr[six$gene_symbol[six$gene_symbol %in% rownames(expr)], , drop = FALSE]
clin <- clin[match(colnames(expr), sample_id)]
keep <- !is.na(clin$OS_time) & !is.na(clin$OS_status) & clin$OS_time > 0
expr <- expr[, keep, drop = FALSE]
clin <- clin[keep]

coef <- six$coefficient
names(coef) <- six$gene_symbol
score <- as.numeric(crossprod(coef[rownames(expr)], expr))

cat("n=", nrow(clin), " events=", sum(clin$OS_status == 1), "\n", sep = "")
cat("time_range=", paste(range(clin$OS_time), collapse = "-"), "\n", sep = "")
cat("score_range=", paste(range(score), collapse = "-"), "\n", sep = "")
print(table(clin$OS_status))
cat("timeROC_available=", requireNamespace("timeROC", quietly = TRUE), "\n", sep = "")
cat("survivalROC_available=", requireNamespace("survivalROC", quietly = TRUE), "\n", sep = "")

for (tp in c(365, 1095, 1825)) {
  cat("timepoint=", tp,
      " events_before=", sum(clin$OS_time <= tp & clin$OS_status == 1),
      " controls_beyond=", sum(clin$OS_time > tp), "\n", sep = "")
  tr <- tryCatch(
    timeROC::timeROC(T = clin$OS_time, delta = clin$OS_status, marker = score, cause = 1, times = tp, iid = FALSE),
    error = function(e) e
  )
  if (inherits(tr, "error")) {
    cat("timeROC_error=", conditionMessage(tr), "\n", sep = "")
  } else {
    cat("timeROC_AUC=", paste(tr$AUC, collapse = ","), "\n", sep = "")
  }
  sr <- tryCatch(
    survivalROC::survivalROC(
      Stime = clin$OS_time,
      status = clin$OS_status,
      marker = score,
      predict.time = tp,
      method = "NNE",
      span = 0.25 * (length(score)^(-0.20))
    ),
    error = function(e) e
  )
  if (inherits(sr, "error")) {
    cat("survivalROC_error=", conditionMessage(sr), "\n", sep = "")
  } else {
    cat("survivalROC_AUC=", sr$AUC, "\n", sep = "")
  }
}
