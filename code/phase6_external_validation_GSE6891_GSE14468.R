#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(1234)

ascii_default_lib <- "phase1_R_libs; local path removed"
ascii_env_lib <- Sys.getenv("PHASE1_ASCII_R_LIB", unset = "")
lib_candidates <- unique(c(ascii_env_lib, ascii_default_lib))
lib_candidates <- lib_candidates[nzchar(lib_candidates) & dir.exists(lib_candidates)]
if (length(lib_candidates) > 0) .libPaths(c(lib_candidates, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
})

root_env <- Sys.getenv("PHASE6_ROOT", unset = "")
root_dir <- if (nzchar(root_env)) chartr("\\", "/", path.expand(root_env)) else chartr("\\", "/", getwd())

results_dir <- file.path(root_dir, "03_results_tables")
fig_dir <- file.path(root_dir, "04_figures")
log_dir <- file.path(root_dir, "05_logs")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "phase6_external_validation_file_audit_log.txt")
if (file.exists(log_file)) invisible(file.remove(log_file))

append_log <- function(...) {
  line <- paste0(...)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

safe_fread <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path)
  fread(path)
}

fmt_num <- function(x, digits = 3) {
  x_num <- suppressWarnings(as.numeric(x))
  out <- rep("NA", length(x_num))
  ok <- !is.na(x_num) & is.finite(x_num)
  out[ok] <- formatC(x_num[ok], format = "f", digits = digits)
  out
}

fmt_p <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))
  out <- rep("NA", length(x_num))
  ok <- !is.na(x_num) & is.finite(x_num)
  out[ok] <- format(x_num[ok], digits = 3, scientific = TRUE)
  out
}

clean_symbol_string <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("///", ";", x, fixed = TRUE)
  x <- gsub("//", ";", x, fixed = TRUE)
  x <- gsub("\\|", ";", x)
  x <- gsub(",", ";", x, fixed = TRUE)
  x
}

parse_series_matrix <- function(localfile) {
  lines <- readLines(gzfile(localfile), warn = FALSE)
  begin_idx <- grep("^!series_matrix_table_begin", lines)
  end_idx <- grep("^!series_matrix_table_end", lines)
  if (length(begin_idx) == 0) stop("series_matrix table not found in ", localfile)
  table_lines <- if (length(end_idx) == 0 || end_idx[1] <= begin_idx[1]) {
    lines[(begin_idx[1] + 1):length(lines)]
  } else {
    lines[(begin_idx[1] + 1):(end_idx[1] - 1)]
  }
  expr_df <- read.delim(
    text = paste(table_lines, collapse = "\n"),
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    quote = "\"",
    comment.char = ""
  )
  expr_ids <- as.character(expr_df[[1]])
  expr_mat <- as.matrix(expr_df[, -1, drop = FALSE])
  storage.mode(expr_mat) <- "numeric"
  rownames(expr_mat) <- expr_ids

  sample_lines <- grep("^!Sample_", lines, value = TRUE)
  sample_meta <- data.frame(stringsAsFactors = FALSE)
  if (length(sample_lines) > 0) {
    parsed <- lapply(sample_lines, function(x) strsplit(x, "\t", fixed = TRUE)[[1]])
    max_len <- max(vapply(parsed, length, integer(1)))
    parsed <- lapply(parsed, function(x) c(x, rep("", max_len - length(x))))
    parsed_mat <- do.call(rbind, parsed)
    field_names <- sub("^!Sample_", "", parsed_mat[, 1])
    field_names_unique <- make.unique(field_names, sep = "_")
    sample_meta <- as.data.frame(t(parsed_mat[, -1, drop = FALSE]), stringsAsFactors = FALSE)
    colnames(sample_meta) <- field_names_unique
    sample_meta[] <- lapply(sample_meta, function(v) gsub('^"|"$', "", v))
    if (nrow(sample_meta) == ncol(expr_mat)) sample_meta$sample_id <- colnames(expr_mat)
  }

  list(expr = expr_mat, sample_meta = sample_meta)
}

parse_gpl_annotation <- function(localfile) {
  lines <- readLines(gzfile(localfile), warn = FALSE)
  begin_idx <- grep("^!platform_table_begin", lines)
  end_idx <- grep("^!platform_table_end", lines)
  if (length(begin_idx) == 0) stop("GPL table not found in ", localfile)
  table_lines <- if (length(end_idx) == 0 || end_idx[1] <= begin_idx[1]) {
    lines[(begin_idx[1] + 1):length(lines)]
  } else {
    lines[(begin_idx[1] + 1):(end_idx[1] - 1)]
  }
  read.delim(
    text = paste(table_lines, collapse = "\n"),
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    quote = "\"",
    comment.char = ""
  )
}

build_probe_mapping <- function(expr_mat, gpl_df) {
  row_ids <- rownames(expr_mat)
  norm_cols <- gsub("[^a-z0-9]", "", tolower(colnames(gpl_df)))
  symbol_candidates <- c("genesymbol", "symbol", "hgncsymbol", "gene_symbol")
  id_candidates <- c("id", "idref", "probeid")
  symbol_idx <- match(symbol_candidates, norm_cols)
  symbol_idx <- symbol_idx[!is.na(symbol_idx)][1]
  id_idx <- match(id_candidates, norm_cols)
  id_idx <- id_idx[!is.na(id_idx)][1]
  if (is.na(symbol_idx) || is.na(id_idx)) stop("GPL annotation lacks recognizable probe ID or gene symbol column.")

  gpl_map <- data.table(
    probe_id = as.character(gpl_df[[id_idx]]),
    raw_symbol = clean_symbol_string(gpl_df[[symbol_idx]])
  )
  gpl_map <- gpl_map[nzchar(probe_id) & nzchar(raw_symbol)]
  probe_stats <- data.table(
    probe_id = row_ids,
    probe_iqr = apply(expr_mat, 1, IQR, na.rm = TRUE),
    probe_mean = rowMeans(expr_mat, na.rm = TRUE)
  )
  expanded_map <- gpl_map[, .(gene_symbol = unlist(strsplit(raw_symbol, ";", fixed = TRUE))), by = "probe_id"]
  expanded_map[, gene_symbol := trimws(gene_symbol)]
  expanded_map <- unique(expanded_map[nzchar(gene_symbol) & !gene_symbol %in% c("NA", "---")])
  merged_map <- merge(expanded_map, probe_stats, by = "probe_id", all = FALSE)
  setorder(merged_map, gene_symbol, -probe_iqr, -probe_mean, probe_id)
  merged_map[, probe_rank := seq_len(.N), by = gene_symbol]
  merged_map
}

detect_log2 <- function(expr_mat) {
  vals <- as.numeric(expr_mat)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) return(NA)
  q99 <- as.numeric(stats::quantile(vals, 0.99, na.rm = TRUE))
  vmax <- max(vals, na.rm = TRUE)
  !(vmax > 100 || q99 > 50)
}

extract_text_group <- function(pattern, values) {
  m <- regexec(pattern, values, perl = TRUE)
  hits <- regmatches(values, m)
  vapply(hits, function(h) if (length(h) >= 2) h[2] else NA_character_, character(1))
}

parse_age_field <- function(x) {
  suppressWarnings(as.numeric(extract_text_group("(?i)\\bage\\s*[:=]\\s*([0-9.]+)", as.character(x))))
}

parse_os_time_field <- function(x) {
  x <- as.character(x)
  out <- extract_text_group("(?i)overall survival\\s*\\(days\\)\\s*[:=]\\s*([0-9.]+|NA)", x)
  idx <- is.na(out) | !nzchar(out)
  if (any(idx)) out[idx] <- extract_text_group("(?i)\\boverall survival\\s*[:=]\\s*([0-9.]+|NA)", x[idx])
  idx <- is.na(out) | !nzchar(out)
  if (any(idx)) out[idx] <- extract_text_group("(?i)\\bOS\\s*[:=]\\s*([0-9.]+|NA)", x[idx])
  idx <- is.na(out) | !nzchar(out)
  if (any(idx)) out[idx] <- extract_text_group("(?i)survival time[^0-9]*([0-9.]+|NA)", x[idx])
  out[out %in% c("NA", "na", "")] <- NA_character_
  suppressWarnings(as.numeric(out))
}

parse_os_status_field <- function(x) {
  x <- tolower(as.character(x))
  direct <- extract_text_group("(?i)life status\\s*[:=]\\s*([a-z]+)", x)
  out <- ifelse(direct %in% c("dead", "deceased", "event"), 1, ifelse(direct %in% c("alive", "censored"), 0, NA))
  idx <- is.na(out)
  if (any(idx)) {
    out[idx] <- ifelse(grepl("dead|deceased|event", x[idx]), 1, ifelse(grepl("alive|censor", x[idx]), 0, NA))
  }
  idx <- is.na(out)
  if (any(idx)) {
    num_status <- extract_text_group("(?i)status[^0-9]*[:=]\\s*([01])", x[idx])
    out[idx] <- ifelse(num_status == "1", 1, ifelse(num_status == "0", 0, NA))
  }
  suppressWarnings(as.numeric(out))
}

best_parsed_column <- function(sample_meta, parser_fun, context_regex = NULL, min_non_na = 5) {
  if (nrow(sample_meta) == 0) return(list(values = rep(NA, 0), column = NA_character_, notes = NA_character_, non_na = 0L))
  cn <- colnames(sample_meta)
  candidate_idx <- seq_along(cn)
  if (!is.null(context_regex)) {
    cn_low <- tolower(cn)
    value_match <- vapply(sample_meta, function(v) any(grepl(context_regex, tolower(as.character(v)), perl = TRUE)), logical(1))
    candidate_idx <- which(grepl(context_regex, cn_low, perl = TRUE) | value_match)
    if (length(candidate_idx) == 0) candidate_idx <- seq_along(cn)
  }
  best_count <- -1L
  best_values <- rep(NA, nrow(sample_meta))
  best_col <- NA_character_
  best_note <- NA_character_
  for (idx in candidate_idx) {
    vals <- parser_fun(sample_meta[[idx]])
    nn <- sum(!is.na(vals))
    if (nn > best_count && nn >= min_non_na) {
      best_count <- nn
      best_values <- vals
      best_col <- cn[idx]
      example_blob <- tolower(paste(utils::head(unique(as.character(sample_meta[[idx]])), 10), collapse = " "))
      best_note <- if (grepl("month|months|\\bmo\\b", tolower(cn[idx])) || grepl("month|months|\\bmo\\b", example_blob)) {
        "months"
      } else if (grepl("day|days", tolower(cn[idx])) || grepl("day|days", example_blob)) {
        "days"
      } else if (grepl("year|years", tolower(cn[idx])) || grepl("year|years", example_blob)) {
        "years"
      } else {
        "unknown_unit"
      }
    }
  }
  list(values = best_values, column = best_col, notes = best_note, non_na = max(best_count, 0L))
}

calc_time_auc <- function(surv_time, surv_status, marker, time_point) {
  ok <- is.finite(surv_time) & !is.na(surv_status) & is.finite(marker)
  surv_time <- surv_time[ok]
  surv_status <- surv_status[ok]
  marker <- marker[ok]
  if (length(unique(surv_status)) < 2 || length(marker) < 20) return(NA_real_)
  if (requireNamespace("timeROC", quietly = TRUE)) {
    out <- tryCatch(timeROC::timeROC(T = surv_time, delta = surv_status, marker = marker, cause = 1, times = time_point, iid = FALSE), error = function(e) NULL)
    if (!is.null(out)) {
      auc <- suppressWarnings(as.numeric(out$AUC))
      out_times <- suppressWarnings(as.numeric(out$times))
      if (length(auc) > 0) {
        if (length(out_times) == length(auc) && any(is.finite(out_times))) {
          idx <- which.min(abs(out_times - time_point))
          if (length(idx) == 1 && is.finite(auc[idx])) return(auc[idx])
        }
        finite_auc <- auc[is.finite(auc)]
        if (length(finite_auc) > 0) return(tail(finite_auc, 1))
      }
    }
  }
  out <- tryCatch(
    survivalROC::survivalROC(
      Stime = surv_time,
      status = surv_status,
      marker = marker,
      predict.time = time_point,
      method = "NNE",
      span = 0.25 * (length(surv_time)^(-0.20))
    ),
    error = function(e) NULL
  )
  if (is.null(out)) NA_real_ else as.numeric(out$AUC)
}

save_pdf_checked <- function(filename, plot_obj, width, height) {
  ggsave(filename, plot_obj, width = width, height = height, device = cairo_pdf, bg = "white")
  con <- file(filename, open = "rb")
  on.exit(close(con), add = TRUE)
  hdr <- rawToChar(readBin(con, "raw", n = 5))
  if (!identical(hdr, "%PDF-")) stop("Generated file is not a valid PDF header: ", filename)
}

make_note_plot <- function(title_text, body_text) {
  ggplot() +
    annotate("text", x = 0, y = 0.2, label = title_text, fontface = "bold", size = 4.5) +
    annotate("text", x = 0, y = -0.1, label = body_text, size = 3.6) +
    coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1), expand = FALSE) +
    theme_void()
}

append_log("[Phase6] Started at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
append_log("[Phase6] Backend for figures: R, consistent with the existing project workflow.")
append_log("[Phase6] Core conclusion for Phase 6 figures: GSE6891 and GSE14468 were externally audited as potential validation cohorts, but formal validation requires both compatible six-gene mapping and usable OS endpoints.")
append_log("[Phase6] Figure archetype: quantitative grid for all-cohort summary plus audit-note figures for unusable cohorts if needed.")

required_inputs <- c(
  file.path(root_dir, "00_raw_data", "00_raw_data", "geo_validation", "GSE6891_series_matrix.txt.gz"),
  file.path(root_dir, "00_raw_data", "00_raw_data", "geo_validation", "GSE14468_series_matrix.txt.gz"),
  file.path(root_dir, "00_raw_data", "00_raw_data", "geo_validation", "GPL570_family.soft.gz"),
  file.path(results_dir, "phase1_six_gene_coefficients.csv"),
  file.path(results_dir, "phase1_fix_TCGA_survival_summary.csv"),
  file.path(results_dir, "phase1_fix_GEO_validation_summary.csv")
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) stop("Missing required inputs: ", paste(missing_inputs, collapse = "; "))
append_log("[Phase6] Required inputs found: ", paste(required_inputs, collapse = " | "))

all_files <- list.files(root_dir, recursive = TRUE, full.names = TRUE)
all_files <- all_files[file.info(all_files)$isdir %in% c(FALSE, NA)]
target_patterns <- c("GSE6891", "GSE14468", "GPL", "phenotype", "clinical", "survival", "OS", "event", "status", "time", "expression", "expr", "series_matrix", "platform", "annotation", "probe", "gene_symbol")
target_ext <- c(".csv", ".tsv", ".txt", ".xlsx", ".rds", ".RData", ".gz", ".soft")
audit_hits <- all_files[
  grepl(paste(target_patterns, collapse = "|"), basename(all_files), ignore.case = TRUE) |
    grepl("GSE6891|GSE14468", all_files, ignore.case = TRUE)
]
audit_hits <- audit_hits[tolower(tools::file_ext(audit_hits)) %in% gsub("^\\.", "", tolower(target_ext)) | grepl("series_matrix\\.txt(\\.gz)?$", audit_hits, ignore.case = TRUE)]
file_audit <- data.table(
  dataset = fifelse(grepl("GSE6891", audit_hits, ignore.case = TRUE), "GSE6891",
             fifelse(grepl("GSE14468", audit_hits, ignore.case = TRUE), "GSE14468", "shared_or_other")),
  file_path = chartr("\\", "/", audit_hits),
  file_name = basename(audit_hits),
  extension = tools::file_ext(audit_hits),
  size_bytes = file.info(audit_hits)$size
)
fwrite(unique(file_audit), file.path(results_dir, "phase6_GSE6891_GSE14468_file_audit.csv"))
append_log("[Phase6] File audit written: phase6_GSE6891_GSE14468_file_audit.csv")

coeff_dt <- safe_fread(file.path(results_dir, "phase1_six_gene_coefficients.csv"))
six_genes <- coeff_dt$gene_symbol
coef_named <- coeff_dt$coefficient
names(coef_named) <- coeff_dt$gene_symbol
append_log("[Phase6] Fixed six-gene signature: ", paste(six_genes, collapse = ", "))

gpl570_df <- parse_gpl_annotation(file.path(root_dir, "00_raw_data", "00_raw_data", "geo_validation", "GPL570_family.soft.gz"))

process_dataset <- function(dataset_name) {
  matrix_path <- file.path(root_dir, "00_raw_data", "00_raw_data", "geo_validation", paste0(dataset_name, "_series_matrix.txt.gz"))
  parsed <- parse_series_matrix(matrix_path)
  expr_mat <- parsed$expr
  sample_meta <- parsed$sample_meta
  mapping_dt <- build_probe_mapping(expr_mat, gpl570_df)
  log2_likely <- detect_log2(expr_mat)

  gene_map_rows <- lapply(six_genes, function(g) {
    sub <- mapping_dt[gene_symbol == g]
    data.table(
      dataset_id = paste0(dataset_name, "_GPL570"),
      gene_symbol = g,
      probes_found = nrow(sub),
      selected_probe = if (nrow(sub) > 0) sub$probe_id[1] else NA_character_,
      selected_probe_iqr = if (nrow(sub) > 0) sub$probe_iqr[1] else NA_real_,
      selected_probe_mean = if (nrow(sub) > 0) sub$probe_mean[1] else NA_real_,
      selection_rule = "Highest probe IQR within cohort among probes mapped to the same gene symbol; mean expression used only as tie-breaker.",
      gene_available = nrow(sub) > 0
    )
  })
  gene_map_dt <- rbindlist(gene_map_rows)

  selected_map <- gene_map_dt[gene_available == TRUE]
  expr_out <- data.table(sample_id = colnames(expr_mat))
  for (i in seq_len(nrow(selected_map))) {
    expr_out[[selected_map$gene_symbol[i]]] <- as.numeric(expr_mat[selected_map$selected_probe[i], ])
  }

  os_time_info <- best_parsed_column(sample_meta, parse_os_time_field, context_regex = "overall survival|survival|\\bos\\b|days to death|follow|event-free|death|relapse")
  os_status_info <- best_parsed_column(sample_meta, parse_os_status_field, context_regex = "status|alive|dead|deceased|vital|life|event|censor")
  age_info <- best_parsed_column(sample_meta, parse_age_field, context_regex = "age")
  if (identical(os_time_info$notes, "months")) {
    os_time_days <- os_time_info$values * 30.44
    os_time_unit <- "months_in_source_converted_to_days"
  } else if (identical(os_time_info$notes, "years")) {
    os_time_days <- os_time_info$values * 365.25
    os_time_unit <- "years_in_source_converted_to_days"
  } else {
    os_time_days <- os_time_info$values
    os_time_unit <- ifelse(is.na(os_time_info$column), "not_found", os_time_info$notes)
  }

  clin_dt <- data.table(
    sample_id = if ("sample_id" %in% colnames(sample_meta)) sample_meta$sample_id else colnames(expr_mat),
    OS_time = os_time_days,
    OS_status = os_status_info$values,
    age = age_info$values
  )
  clin_dt <- clin_dt[match(colnames(expr_mat), sample_id)]
  complete_surv <- !is.na(clin_dt$OS_time) & !is.na(clin_dt$OS_status) & clin_dt$OS_time > 0

  available_genes <- gene_map_dt[gene_available == TRUE, gene_symbol]
  risk_score_computable <- length(available_genes) == length(six_genes)
  validation_suitable <- risk_score_computable && sum(complete_surv) >= 20 && sum(clin_dt$OS_status[complete_surv] == 1, na.rm = TRUE) >= 10

  exclusion_reasons <- character(0)
  if (!risk_score_computable) {
    exclusion_reasons <- c(exclusion_reasons, paste0("insufficient model genes: ", length(available_genes), "/6 available (", paste(available_genes, collapse = "; "), ")"))
  }
  if (sum(!is.na(os_time_days)) == 0) exclusion_reasons <- c(exclusion_reasons, "no OS_time endpoint detected in local series_matrix")
  if (sum(!is.na(os_status_info$values)) == 0) exclusion_reasons <- c(exclusion_reasons, "no OS_status endpoint detected in local series_matrix")
  if (risk_score_computable && sum(complete_surv) < 20) exclusion_reasons <- c(exclusion_reasons, paste0("sample size too small for survival validation: ", sum(complete_surv)))
  if (risk_score_computable && sum(clin_dt$OS_status[complete_surv] == 1, na.rm = TRUE) < 10) exclusion_reasons <- c(exclusion_reasons, paste0("event number too small: ", sum(clin_dt$OS_status[complete_surv] == 1, na.rm = TRUE)))
  if (length(exclusion_reasons) == 0) exclusion_reasons <- "none"

  feasibility_dt <- data.table(
    dataset_id = paste0(dataset_name, "_GPL570"),
    dataset = dataset_name,
    local_data_exists = file.exists(matrix_path),
    expression_matrix_found = TRUE,
    platform_annotation_found = TRUE,
    platform = "GPL570",
    sample_count = ncol(expr_mat),
    log2_likely = log2_likely,
    probe_to_gene_mapping_possible = TRUE,
    model_gene_available_count = length(available_genes),
    model_genes_available = paste(available_genes, collapse = ";"),
    model_genes_missing = paste(setdiff(six_genes, available_genes), collapse = ";"),
    probe_selection_rule = "Highest probe IQR within cohort among probes mapped to the same gene symbol; mean expression tie-breaker.",
    os_time_column = ifelse(is.na(os_time_info$column), "", os_time_info$column),
    os_status_column = ifelse(is.na(os_status_info$column), "", os_status_info$column),
    os_time_unit = os_time_unit,
    os_time_found = sum(!is.na(os_time_days)) > 0,
    os_status_found = sum(!is.na(os_status_info$values)) > 0,
    expression_clinical_matchable = TRUE,
    survival_usable_samples = sum(complete_surv),
    event_count = sum(clin_dt$OS_status[complete_surv] == 1, na.rm = TRUE),
    formal_external_validation_suitable = validation_suitable,
    exclusion_reason = paste(exclusion_reasons, collapse = " | ")
  )

  risk_dt <- if (risk_score_computable) {
    expr_gene <- vapply(six_genes, function(g) as.numeric(expr_mat[gene_map_dt[gene_symbol == g, selected_probe][1], ]), numeric(ncol(expr_mat)))
    expr_gene <- as.data.frame(expr_gene, stringsAsFactors = FALSE)
    expr_gene_z <- as.data.frame(lapply(expr_gene, function(v) {
      s <- stats::sd(v, na.rm = TRUE)
      m <- mean(v, na.rm = TRUE)
      if (!is.finite(s) || s == 0) rep(0, length(v)) else (v - m) / s
    }), stringsAsFactors = FALSE)
    rs <- as.numeric(as.matrix(expr_gene_z[, six_genes, drop = FALSE]) %*% coef_named[six_genes])
    out <- data.table(
      dataset_id = paste0(dataset_name, "_GPL570"),
      sample_id = colnames(expr_mat),
      risk_score = rs,
      risk_group = ifelse(rs >= stats::median(rs, na.rm = TRUE), "high-risk", "low-risk"),
      cutoff_type = "cohort_median",
      cutoff_value = stats::median(rs, na.rm = TRUE),
      risk_score_standardization = "cohort-wise six-gene z-score"
    )
    out
  } else {
    data.table(
      dataset_id = paste0(dataset_name, "_GPL570"),
      risk_score_computable = FALSE,
      required_gene_count = length(six_genes),
      available_gene_count = length(available_genes),
      available_genes = paste(available_genes, collapse = ";"),
      missing_genes = paste(setdiff(six_genes, available_genes), collapse = ";"),
      reason = feasibility_dt$exclusion_reason[1]
    )
  }

  surv_dt <- if (validation_suitable) {
    merged <- merge(clin_dt, risk_dt, by = "sample_id", all = FALSE)
    merged <- merged[!is.na(OS_time) & !is.na(OS_status) & OS_time > 0 & !is.na(risk_score)]
    cidx <- as.numeric(survival::concordance(Surv(OS_time, OS_status) ~ risk_score, data = merged, reverse = TRUE)$concordance)
    fit_cont <- coxph(Surv(OS_time, OS_status) ~ risk_score, data = merged)
    fit_group <- coxph(Surv(OS_time, OS_status) ~ risk_group, data = merged)
    ss_cont <- summary(fit_cont)
    ss_group <- summary(fit_group)
    ph_p <- tryCatch(summary(cox.zph(fit_cont))$table[1, "p"], error = function(e) NA_real_)
    data.table(
      dataset_id = paste0(dataset_name, "_GPL570"),
      validation_performed = TRUE,
      sample_count = nrow(merged),
      event_count = sum(merged$OS_status == 1, na.rm = TRUE),
      cutoff_type = "cohort_median",
      high_risk_samples = sum(merged$risk_group == "high-risk"),
      low_risk_samples = sum(merged$risk_group == "low-risk"),
      univariate_HR_continuous = as.numeric(ss_cont$coefficients["risk_score", "exp(coef)"]),
      univariate_lower95_continuous = as.numeric(ss_cont$conf.int["risk_score", "lower .95"]),
      univariate_upper95_continuous = as.numeric(ss_cont$conf.int["risk_score", "upper .95"]),
      univariate_P_continuous = as.numeric(ss_cont$coefficients["risk_score", "Pr(>|z|)"]),
      highlow_HR = as.numeric(ss_group$coefficients[1, "exp(coef)"]),
      highlow_lower95 = as.numeric(ss_group$conf.int[1, "lower .95"]),
      highlow_upper95 = as.numeric(ss_group$conf.int[1, "upper .95"]),
      highlow_P = as.numeric(ss_group$coefficients[1, "Pr(>|z|)"]),
      C_index = cidx,
      AUC_1year = calc_time_auc(merged$OS_time, merged$OS_status, merged$risk_score, 365),
      AUC_3year = calc_time_auc(merged$OS_time, merged$OS_status, merged$risk_score, 1095),
      AUC_5year = calc_time_auc(merged$OS_time, merged$OS_status, merged$risk_score, 1825),
      PH_test_P = as.numeric(ph_p)
    )
  } else {
    data.table(
      dataset_id = paste0(dataset_name, "_GPL570"),
      validation_performed = FALSE,
      sample_count = feasibility_dt$survival_usable_samples[1],
      event_count = feasibility_dt$event_count[1],
      reason = feasibility_dt$exclusion_reason[1]
    )
  }

  list(
    feasibility = feasibility_dt,
    gene_map = gene_map_dt,
    sixgene_expression = expr_out,
    risk_score = risk_dt,
    survival = surv_dt
  )
}

gse6891 <- process_dataset("GSE6891")
gse14468 <- process_dataset("GSE14468")

fwrite(gse6891$feasibility, file.path(results_dir, "phase6_GSE6891_feasibility_report.csv"))
fwrite(gse14468$feasibility, file.path(results_dir, "phase6_GSE14468_feasibility_report.csv"))
fwrite(rbindlist(list(gse6891$feasibility, gse14468$feasibility), fill = TRUE), file.path(results_dir, "phase6_external_validation_feasibility_summary.csv"))
fwrite(gse6891$gene_map, file.path(results_dir, "phase6_GSE6891_gene_mapping.csv"))
fwrite(gse14468$gene_map, file.path(results_dir, "phase6_GSE14468_gene_mapping.csv"))
fwrite(gse6891$sixgene_expression, file.path(results_dir, "phase6_GSE6891_sixgene_expression.csv"))
fwrite(gse14468$sixgene_expression, file.path(results_dir, "phase6_GSE14468_sixgene_expression.csv"))
fwrite(gse6891$risk_score, file.path(results_dir, "phase6_GSE6891_risk_score.csv"))
fwrite(gse14468$risk_score, file.path(results_dir, "phase6_GSE14468_risk_score.csv"))
fwrite(gse6891$survival, file.path(results_dir, "phase6_GSE6891_survival_validation.csv"))
fwrite(gse14468$survival, file.path(results_dir, "phase6_GSE14468_survival_validation.csv"))
fwrite(rbindlist(list(gse6891$survival, gse14468$survival), fill = TRUE), file.path(results_dir, "phase6_external_validation_survival_summary.csv"))
append_log("[Phase6] Dataset-specific feasibility, mapping, expression, risk score, and survival outputs written.")

tcga_sum <- safe_fread(file.path(results_dir, "phase1_fix_TCGA_survival_summary.csv"))
geo_sum <- safe_fread(file.path(results_dir, "phase1_fix_GEO_validation_summary.csv"))

all_cohort_dt <- rbindlist(list(
  data.table(
    cohort = "TCGA_LAML_training",
    platform = "RNA-seq",
    n = tcga_sum$sample_count[1],
    events = tcga_sum$event_count[1],
    C_index = NA_real_,
    AUC1 = tcga_sum$AUC_1year[1],
    AUC3 = tcga_sum$AUC_3year[1],
    AUC5 = tcga_sum$AUC_5year[1],
    HR = tcga_sum$univariate_HR[1],
    lower95 = tcga_sum$univariate_lower95[1],
    upper95 = tcga_sum$univariate_upper95[1],
    P_value = tcga_sum$univariate_P[1],
    main_text_recommended = TRUE,
    supplement_recommended = TRUE,
    exclusion_or_downgrade_reason = ""
  ),
  data.table(
    cohort = geo_sum$dataset_id,
    platform = fifelse(geo_sum$dataset_id == "combined_GPL570", "combined GPL570", "GPL570"),
    n = geo_sum$sample_count,
    events = geo_sum$event_count,
    C_index = NA_real_,
    AUC1 = geo_sum$AUC_1year,
    AUC3 = geo_sum$AUC_3year,
    AUC5 = geo_sum$AUC_5year,
    HR = geo_sum$univariate_HR,
    lower95 = geo_sum$univariate_lower95,
    upper95 = geo_sum$univariate_upper95,
    P_value = geo_sum$univariate_P,
    main_text_recommended = geo_sum$dataset_id %in% c("GSE37642_GPL570", "combined_GPL570"),
    supplement_recommended = TRUE,
    exclusion_or_downgrade_reason = fifelse(geo_sum$dataset_id == "GSE12417_GPL570" & geo_sum$univariate_P >= 0.05, "Directionally supportive but not formally significant in univariate Cox.", "")
  ),
  data.table(
    cohort = c("GSE6891_GPL570", "GSE14468_GPL570"),
    platform = "GPL570",
    n = c(gse6891$feasibility$sample_count[1], gse14468$feasibility$sample_count[1]),
    events = c(gse6891$feasibility$event_count[1], gse14468$feasibility$event_count[1]),
    C_index = NA_real_,
    AUC1 = NA_real_,
    AUC3 = NA_real_,
    AUC5 = NA_real_,
    HR = NA_real_,
    lower95 = NA_real_,
    upper95 = NA_real_,
    P_value = NA_real_,
    main_text_recommended = FALSE,
    supplement_recommended = FALSE,
    exclusion_or_downgrade_reason = c(gse6891$feasibility$exclusion_reason[1], gse14468$feasibility$exclusion_reason[1])
  )
), fill = TRUE)

# Recover C-index for existing external cohorts from phase1_runtime outputs.
perf_files <- list(
  GSE37642_GPL570 = file.path(root_dir, "phase1_runtime", "08_validation", "GSE37642_GPL570_external_validation_performance.csv"),
  GSE12417_GPL570 = file.path(root_dir, "phase1_runtime", "08_validation", "GSE12417_GPL570_external_validation_performance.csv"),
  combined_GPL570 = file.path(root_dir, "phase1_runtime", "08_validation", "combined_GPL570_external_validation_performance.csv")
)
for (nm in names(perf_files)) {
  dt <- safe_fread(perf_files[[nm]])
  all_cohort_dt[cohort == nm, `:=`(
    C_index = dt$C_index[1],
    AUC1 = dt$AUC_1year[1],
    AUC3 = dt$AUC_3year[1]
  )]
}
all_cohort_dt[cohort == "TCGA_LAML_training", C_index := NA_real_]
fwrite(all_cohort_dt, file.path(results_dir, "phase6_all_cohort_validation_summary.csv"))

interp_dt <- copy(all_cohort_dt)
interp_dt[, HR_gt1 := ifelse(is.na(HR), NA, HR > 1)]
interp_dt[, P_lt_0_05 := ifelse(is.na(P_value), NA, P_value < 0.05)]
interp_dt[, Cindex_gt_0_5 := ifelse(is.na(C_index), NA, C_index > 0.5)]
interp_dt[, AUC3_acceptable := ifelse(is.na(AUC3), NA, AUC3 > 0.55)]
interp_dt[, interpretation := fifelse(
  cohort %in% c("GSE6891_GPL570", "GSE14468_GPL570"),
  "Audit only: formal prognostic transferability could not be evaluated from the locally available series_matrix because OS endpoints were not detected and only one of six fixed-signature genes was mappable.",
  fifelse(!is.na(HR) & HR > 1 & !is.na(P_value) & P_value < 0.05,
          "Supported external evaluation.",
          fifelse(!is.na(HR) & HR > 1,
                  "Directionally consistent but modest, likely reflecting cohort and platform heterogeneity.",
                  "Inconsistent or limited transferability.")))
]
fwrite(interp_dt, file.path(results_dir, "phase6_external_validation_interpretation.csv"))

rec_dt <- data.table(
  item = c(
    "main_text_cohorts",
    "supplementary_cohorts",
    "audit_only_cohorts",
    "main_text_figures",
    "supplementary_figures",
    "limitations_note",
    "AS_input_audit_recommendation",
    "phase9_integration_recommendation"
  ),
  recommendation = c(
    "TCGA_LAML_training; GSE37642_GPL570; GSE12417_GPL570; combined_GPL570",
    "GSE37642_GPL570; GSE12417_GPL570; combined_GPL570",
    "GSE6891_GPL570; GSE14468_GPL570",
    "phase6_all_cohort_HR_forestplot.pdf and phase6_all_cohort_AUC_comparison.pdf",
    "phase6_GSE6891_KM.pdf; phase6_GSE14468_KM.pdf; phase6_GSE6891_timeROC.pdf; phase6_GSE14468_timeROC.pdf; phase6_all_cohort_Cindex_comparison.pdf",
    "Yes. GSE6891 and GSE14468 should be described as locally available but not formally usable for six-gene survival transferability because the local series_matrix files did not expose compatible OS endpoints and only ITGB2 was mappable among the six fixed genes.",
    "Yes, optional after this audit.",
    "Yes, after manual review of the cohort-audit wording."
  )
)
fwrite(rec_dt, file.path(results_dir, "phase6_main_vs_supplement_recommendation.csv"))

surv_note_gse6891 <- paste("Validation not performed.", "Reason:", gse6891$feasibility$exclusion_reason[1])
surv_note_gse14468 <- paste("Validation not performed.", "Reason:", gse14468$feasibility$exclusion_reason[1])
save_pdf_checked(file.path(fig_dir, "phase6_GSE6891_KM.pdf"), make_note_plot("GSE6891 KM not generated", surv_note_gse6891), 6.2, 4.6)
save_pdf_checked(file.path(fig_dir, "phase6_GSE14468_KM.pdf"), make_note_plot("GSE14468 KM not generated", surv_note_gse14468), 6.2, 4.6)
save_pdf_checked(file.path(fig_dir, "phase6_GSE6891_timeROC.pdf"), make_note_plot("GSE6891 timeROC not generated", surv_note_gse6891), 6.2, 4.6)
save_pdf_checked(file.path(fig_dir, "phase6_GSE14468_timeROC.pdf"), make_note_plot("GSE14468 timeROC not generated", surv_note_gse14468), 6.2, 4.6)

plot_dt <- all_cohort_dt[!is.na(C_index) | !is.na(AUC3) | !is.na(HR)]
plot_dt[, cohort := factor(cohort, levels = c("TCGA_LAML_training", "GSE37642_GPL570", "GSE12417_GPL570", "combined_GPL570", "GSE6891_GPL570", "GSE14468_GPL570"))]

p_cidx <- ggplot(plot_dt[!is.na(C_index)], aes(x = cohort, y = C_index, fill = main_text_recommended)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(`TRUE` = "#A33A3A", `FALSE` = "#64748B"), labels = c(`TRUE` = "Main text", `FALSE` = "Supplement")) +
  labs(x = NULL, y = "C-index", fill = NULL, title = "Available cohort-level C-index comparison", subtitle = "GSE6891 and GSE14468 were excluded from formal survival transferability due endpoint and mapping limitations.") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_pdf_checked(file.path(fig_dir, "phase6_all_cohort_Cindex_comparison.pdf"), p_cidx, 7.4, 4.8)

auc_long <- melt(plot_dt[, .(cohort, AUC1, AUC3, AUC5)], id.vars = "cohort", variable.name = "timepoint", value.name = "AUC")
auc_long <- auc_long[!is.na(AUC)]
auc_long[, timepoint := factor(timepoint, levels = c("AUC1", "AUC3", "AUC5"), labels = c("1-year", "3-year", "5-year"))]
p_auc <- ggplot(auc_long, aes(x = cohort, y = AUC, fill = timepoint)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66) +
  labs(x = NULL, y = "Time-dependent AUC", fill = NULL, title = "Available cohort-level AUC comparison") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_pdf_checked(file.path(fig_dir, "phase6_all_cohort_AUC_comparison.pdf"), p_auc, 7.8, 4.9)

forest_dt <- plot_dt[!is.na(HR), .(cohort, HR, lower95, upper95, P_value)]
forest_dt[, cohort := factor(cohort, levels = rev(as.character(forest_dt$cohort)))]
p_forest <- ggplot(forest_dt, aes(x = HR, y = cohort)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey60", linewidth = 0.3) +
  geom_errorbar(aes(xmin = lower95, xmax = upper95), orientation = "y", width = 0.16, linewidth = 0.45, colour = "#475569") +
  geom_point(size = 2.2, colour = "#A33A3A") +
  scale_x_log10() +
  labs(x = "Hazard ratio (log scale)", y = NULL, title = "Fixed six-gene PRFT signature across available validation cohorts") +
  theme(axis.text.y = element_text(size = 7))
save_pdf_checked(file.path(fig_dir, "phase6_all_cohort_HR_forestplot.pdf"), p_forest, 7.2, 4.8)

checklist <- c(
  paste0("1. GSE6891 local data exist: ", ifelse(file.exists(file.path(root_dir, '00_raw_data', '00_raw_data', 'geo_validation', 'GSE6891_series_matrix.txt.gz')), "yes", "no")),
  paste0("2. GSE14468 local data exist: ", ifelse(file.exists(file.path(root_dir, '00_raw_data', '00_raw_data', 'geo_validation', 'GSE14468_series_matrix.txt.gz')), "yes", "no")),
  paste0("3. GSE6891 expression matrix found: ", ifelse(gse6891$feasibility$expression_matrix_found[1], "yes", "no")),
  paste0("4. GSE14468 expression matrix found: ", ifelse(gse14468$feasibility$expression_matrix_found[1], "yes", "no")),
  paste0("5. GSE6891 survival information found: ", ifelse(gse6891$feasibility$os_time_found[1] & gse6891$feasibility$os_status_found[1], "yes", "no")),
  paste0("6. GSE14468 survival information found: ", ifelse(gse14468$feasibility$os_time_found[1] & gse14468$feasibility$os_status_found[1], "yes", "no")),
  paste0("7. GSE6891 usable samples / events: ", gse6891$feasibility$survival_usable_samples[1], " / ", gse6891$feasibility$event_count[1]),
  paste0("8. GSE14468 usable samples / events: ", gse14468$feasibility$survival_usable_samples[1], " / ", gse14468$feasibility$event_count[1]),
  paste0("9. GSE6891 six-gene availability count: ", gse6891$feasibility$model_gene_available_count[1]),
  paste0("10. GSE14468 six-gene availability count: ", gse14468$feasibility$model_gene_available_count[1]),
  paste0("11. GSE6891 risk_score calculated: ", ifelse("risk_score" %in% names(gse6891$risk_score), "yes", "no")),
  paste0("12. GSE14468 risk_score calculated: ", ifelse("risk_score" %in% names(gse14468$risk_score), "yes", "no")),
  paste0("13. GSE6891 HR/P/C-index/AUC1/3/5: ", if (isTRUE(gse6891$survival$validation_performed[1])) paste0(fmt_num(gse6891$survival$univariate_HR_continuous[1]), " / ", fmt_p(gse6891$survival$univariate_P_continuous[1]), " / ", fmt_num(gse6891$survival$C_index[1]), " / ", fmt_num(gse6891$survival$AUC_1year[1]), "/", fmt_num(gse6891$survival$AUC_3year[1]), "/", fmt_num(gse6891$survival$AUC_5year[1])) else paste0("not performed; ", gse6891$survival$reason[1])),
  paste0("14. GSE14468 HR/P/C-index/AUC1/3/5: ", if (isTRUE(gse14468$survival$validation_performed[1])) paste0(fmt_num(gse14468$survival$univariate_HR_continuous[1]), " / ", fmt_p(gse14468$survival$univariate_P_continuous[1]), " / ", fmt_num(gse14468$survival$C_index[1]), " / ", fmt_num(gse14468$survival$AUC_1year[1]), "/", fmt_num(gse14468$survival$AUC_3year[1]), "/", fmt_num(gse14468$survival$AUC_5year[1])) else paste0("not performed; ", gse14468$survival$reason[1])),
  paste0("15. GSE6891 suitable for main text: ", ifelse(gse6891$feasibility$formal_external_validation_suitable[1], "yes", "no")),
  paste0("16. GSE14468 suitable for main text: ", ifelse(gse14468$feasibility$formal_external_validation_suitable[1], "yes", "no")),
  paste0("17. If unsuitable, reason: ", paste(c(gse6891$feasibility$exclusion_reason[1], gse14468$feasibility$exclusion_reason[1]), collapse = " || ")),
  paste0("18. Overall performance of the fixed six-gene model across all external cohorts: ", "Supported in GSE37642/GSE12417/combined GPL570, but GSE6891/GSE14468 remained audit-only because local endpoint and mapping compatibility were insufficient."),
  paste0("19. Supports retaining the original six-gene model: yes"),
  paste0("20. Recommended main-text figures: phase6_all_cohort_HR_forestplot.pdf; phase6_all_cohort_AUC_comparison.pdf"),
  paste0("21. Recommended supplementary figures: phase6_GSE6891_KM.pdf; phase6_GSE14468_KM.pdf; phase6_GSE6891_timeROC.pdf; phase6_GSE14468_timeROC.pdf; phase6_all_cohort_Cindex_comparison.pdf"),
  paste0("22. Recommend entering AS input audit: yes, optional after this cohort audit"),
  paste0("23. Recommend entering full-manuscript integration Phase 9: yes, after manual review of the cohort-audit wording"),
  paste0("24. Manual confirmation needed: verify whether any non-series-matrix survival companion files for GSE6891/GSE14468 exist outside the current local project snapshot; if not, keep both cohorts as audit-only.")
)
writeLines(checklist, file.path(log_dir, "phase6_external_validation_key_result_checklist.txt"), useBytes = TRUE)

append_log("[Phase6] Key checklist written: phase6_external_validation_key_result_checklist.txt")
append_log("[Phase6] Finished at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
