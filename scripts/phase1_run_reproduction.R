#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

runtime_dir <- "."
audit_root <- ".."

log_dir <- file.path("..", "05_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, "phase1_reproduction_log.txt")

append_log <- function(...) {
  msg <- paste0(...)
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

write_block <- function(text_lines) {
  if (length(text_lines) == 0) {
    return(invisible(NULL))
  }
  cat(paste(text_lines, collapse = "\n"), "\n", file = log_file, append = TRUE)
}

ascii_lib <- Sys.getenv("PHASE1_ASCII_R_LIB", unset = "")
if (!nzchar(ascii_lib)) {
  default_ascii_lib <- "C:/Users/ROBIN-~1/AppData/Local/Temp/phase1_R_libs"
  if (dir.exists(default_ascii_lib)) {
    ascii_lib <- default_ascii_lib
  }
}
if (nzchar(ascii_lib) && dir.exists(ascii_lib)) {
  .libPaths(c(ascii_lib, .libPaths()))
}

required_pkgs <- c(
  "data.table", "edgeR", "dplyr", "GSVA", "ggplot2", "limma",
  "pheatmap", "WGCNA", "survival", "glmnet", "survivalROC"
)
pkg_versions <- vapply(
  required_pkgs,
  function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      as.character(utils::packageVersion(pkg))
    } else {
      NA_character_
    }
  },
  character(1)
)

dir.create(file.path("..", "03_results_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("..", "04_figures", "phase1_KM_ROC"), recursive = TRUE, showWarnings = FALSE)

cat("", file = log_file)
append_log("Phase 1 reproduction started.")
append_log("Audit root: ", audit_root)
append_log("Runtime directory: ", runtime_dir)
append_log("ASCII R library path: ", if (nzchar(ascii_lib)) ascii_lib else "not set")
append_log("Fixed random seed policy: set.seed(1234)")
append_log("Phase 0 report: ", file.path(audit_root, "project_audit_report.md"))
append_log("Phase 0 inventory: ", file.path(audit_root, "05_logs", "file_inventory_phase0.csv"))
append_log("Required package versions:")
write_block(paste0("  - ", names(pkg_versions), ": ", ifelse(is.na(pkg_versions), "NOT_AVAILABLE", pkg_versions)))
append_log("sessionInfo():")
write_block(capture.output(sessionInfo()))

scripts_to_run <- c(
  "04_preprocess_tcga_laml_xena.R",
  "07_prepare_gene_sets.R",
  "08_calculate_ssgsea_scores.R",
  "09_calculate_prft_score_and_group.R",
  "10_deg_prft_high_low_tcga.R",
  "11_wgcna_prft_modules_tcga.R",
  "12_candidate_gene_selection_and_univariate_cox.R",
  "13_lasso_cox_prft_signature_tcga.R",
  "14_external_geo_dataset_feasibility_check.R",
  "15_cross_platform_candidate_coverage_check.R",
  "16_rebuild_cross_platform_lasso_signature_tcga.R",
  "17_external_validation_cross_platform_signature_geo.R",
  "phase1_finalize_outputs.R"
)

rscript_bin <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript_bin)) {
  rscript_bin <- "Rscript"
}

run_script <- function(script_name) {
  script_path <- file.path("15_scripts", script_name)
  if (!file.exists(script_path)) {
    stop("Missing script in runtime: ", script_path)
  }

  append_log("Running script: ", script_name)
  step_out <- tempfile(pattern = sub("\\.R$", "_", script_name), fileext = "_out.log")
  step_err <- tempfile(pattern = sub("\\.R$", "_", script_name), fileext = "_err.log")
  env_vars <- character(0)
  if (nzchar(ascii_lib) && dir.exists(ascii_lib)) {
    env_vars <- c(env_vars, paste0("PHASE1_ASCII_R_LIB=", ascii_lib))
  }

  status <- system2(
    command = rscript_bin,
    args = c("--vanilla", script_path),
    stdout = step_out,
    stderr = step_err,
    env = env_vars
  )

  step_lines <- c(
    tryCatch(readLines(step_out, warn = FALSE), error = function(e) character(0)),
    tryCatch(readLines(step_err, warn = FALSE), error = function(e) character(0))
  )
  if (length(step_lines) > 0) {
    write_block(paste0("  ", step_lines))
  }
  unlink(c(step_out, step_err))

  if (!identical(status, 0L)) {
    append_log("Script failed: ", script_name, " ; exit status = ", status)
    stop("Phase 1 reproduction stopped at ", script_name)
  }

  append_log("Completed script: ", script_name)
}

for (script_name in scripts_to_run) {
  run_script(script_name)
}

append_log("Phase 1 reproduction completed successfully.")
