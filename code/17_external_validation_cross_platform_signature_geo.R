#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
  library(survivalROC)
})
source("15_scripts/plot_label_utils.R")

dir.create("08_validation", recursive = TRUE, showWarnings = FALSE)
dir.create("13_figures", recursive = TRUE, showWarnings = FALSE)
dir.create("14_tables", recursive = TRUE, showWarnings = FALSE)
dir.create("16_logs", recursive = TRUE, showWarnings = FALSE)

save_session_info <- function(path) {
  writeLines(capture.output(sessionInfo()), con = path)
}

required_files <- c(
  "07_signature/final_cross_platform_prft_signature_coefficients.csv",
  "14_tables/tcga_cross_platform_lasso_signature_summary.csv",
  "08_validation/geo_feasibility_summary.csv",
  "08_validation/geo_download_manifest.csv",
  "08_validation/cross_platform_gene_universe_by_dataset.csv",
  "08_validation/cross_platform_rebuild_recommended_strategy.csv",
  "00_raw_data/geo_validation/GPL570_family.soft.gz",
  "00_raw_data/geo_validation/GSE37642-GPL570_series_matrix.txt.gz",
  "00_raw_data/geo_validation/GSE12417-GPL570_series_matrix.txt.gz"
)
missing_required <- required_files[!file.exists(required_files)]
if (length(missing_required) > 0) {
  stop("Missing required files: ", paste(missing_required, collapse = "; "))
}

signature_coef <- fread("07_signature/final_cross_platform_prft_signature_coefficients.csv")
tcga_summary <- fread("14_tables/tcga_cross_platform_lasso_signature_summary.csv")
geo_feas <- fread("08_validation/geo_feasibility_summary.csv")
geo_manifest <- fread("08_validation/geo_download_manifest.csv")
cross_platform_universe <- fread("08_validation/cross_platform_gene_universe_by_dataset.csv")
strategy_dt <- fread("08_validation/cross_platform_rebuild_recommended_strategy.csv")

if (nrow(signature_coef) == 0) {
  stop("final_cross_platform_prft_signature_coefficients.csv is empty.")
}
if (!all(c("gene_symbol", "coefficient") %in% colnames(signature_coef))) {
  stop("Signature coefficient file must contain gene_symbol and coefficient columns.")
}
if (!"median_risk_score_cutoff" %in% colnames(tcga_summary)) {
  stop("tcga_cross_platform_lasso_signature_summary.csv lacks median_risk_score_cutoff.")
}

recommended_scenario <- strategy_dt$recommended_scenario[1]
if (!identical(recommended_scenario, "Scenario_B_GPL570_two_cohorts")) {
  stop("Expected recommended_scenario to be Scenario_B_GPL570_two_cohorts, got: ", recommended_scenario)
}

signature_genes <- signature_coef$gene_symbol
coef_named <- signature_coef$coefficient
names(coef_named) <- signature_coef$gene_symbol
tcga_cutoff <- as.numeric(tcga_summary$median_risk_score_cutoff[1])

target_datasets <- data.table(
  dataset_id = c("GSE37642_GPL570", "GSE12417_GPL570"),
  gse = c("GSE37642", "GSE12417"),
  matrix_file = c("GSE37642-GPL570_series_matrix.txt.gz", "GSE12417-GPL570_series_matrix.txt.gz")
)

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
  if (length(begin_idx) == 0) {
    stop("series_matrix table not found in ", localfile)
  }
  if (length(end_idx) == 0 || end_idx[1] <= begin_idx[1]) {
    table_lines <- lines[(begin_idx[1] + 1):length(lines)]
  } else {
    table_lines <- lines[(begin_idx[1] + 1):(end_idx[1] - 1)]
  }

  expr_df <- read.delim(
    text = paste(table_lines, collapse = "\n"),
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    quote = "\"",
    comment.char = ""
  )
  if (ncol(expr_df) < 2) {
    stop("Expression matrix has fewer than 2 columns in ", localfile)
  }

  id_col <- colnames(expr_df)[1]
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
    if (nrow(sample_meta) == ncol(expr_mat)) {
      sample_meta$sample_id <- colnames(expr_mat)
    }
  }

  list(expr = expr_mat, sample_meta = sample_meta, id_col = id_col)
}

parse_gpl_annotation <- function(localfile) {
  lines <- readLines(gzfile(localfile), warn = FALSE)
  begin_idx <- grep("^!platform_table_begin", lines)
  end_idx <- grep("^!platform_table_end", lines)
  if (length(begin_idx) == 0) {
    stop("GPL table not found in ", localfile)
  }
  if (length(end_idx) == 0 || end_idx[1] <= begin_idx[1]) {
    table_lines <- lines[(begin_idx[1] + 1):length(lines)]
  } else {
    table_lines <- lines[(begin_idx[1] + 1):(end_idx[1] - 1)]
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

  if (is.na(symbol_idx) || is.na(id_idx)) {
    stop("GPL annotation lacks recognizable probe ID or gene symbol column.")
  }

  gpl_map <- data.table(
    probe_id = as.character(gpl_df[[id_idx]]),
    raw_symbol = clean_symbol_string(gpl_df[[symbol_idx]])
  )
  gpl_map <- gpl_map[nzchar(probe_id) & nzchar(raw_symbol)]

  probe_iqr <- apply(expr_mat, 1, IQR, na.rm = TRUE)
  probe_stats <- data.table(
    probe_id = row_ids,
    probe_iqr = probe_iqr
  )

  expanded_map <- gpl_map[, .(gene_symbol = unlist(strsplit(raw_symbol, ";", fixed = TRUE))), by = "probe_id"]
  expanded_map[, gene_symbol := trimws(gene_symbol)]
  expanded_map <- unique(expanded_map[nzchar(gene_symbol) & !gene_symbol %in% c("NA", "---")])

  merged_map <- merge(expanded_map, probe_stats, by = "probe_id", all = FALSE)
  if (nrow(merged_map) == 0) {
    stop("No overlap between expression probe IDs and GPL annotation.")
  }

  setorder(merged_map, gene_symbol, -probe_iqr)
  merged_map[, probe_rank := seq_len(.N), by = gene_symbol]
  merged_map
}

extract_text_group <- function(pattern, values) {
  m <- regexec(pattern, values, perl = TRUE)
  hits <- regmatches(values, m)
  out <- vapply(hits, function(h) {
    if (length(h) >= 2) h[2] else NA_character_
  }, character(1))
  out
}

parse_age_field <- function(x) {
  x <- as.character(x)
  out <- extract_text_group("(?i)\\bage\\s*[:=]\\s*([0-9.]+)", x)
  suppressWarnings(as.numeric(out))
}

parse_os_time_field <- function(x) {
  x <- as.character(x)
  out <- extract_text_group("(?i)overall survival\\s*\\(days\\)\\s*[:=]\\s*([0-9.]+|NA)", x)
  missing_idx <- is.na(out) | !nzchar(out)
  if (any(missing_idx)) {
    out2 <- extract_text_group("(?i)\\bOS\\s*[:=]\\s*([0-9.]+|NA)", x[missing_idx])
    out[missing_idx] <- out2
  }
  missing_idx <- is.na(out) | !nzchar(out)
  if (any(missing_idx)) {
    out3 <- extract_text_group("(?i)survival time[^0-9]*([0-9.]+|NA)", x[missing_idx])
    out[missing_idx] <- out3
  }
  out[out %in% c("NA", "na", "")] <- NA_character_
  suppressWarnings(as.numeric(out))
}

parse_os_status_field <- function(x) {
  x <- tolower(as.character(x))

  direct <- extract_text_group("(?i)life status\\s*[:=]\\s*([a-z]+)", x)
  out <- ifelse(direct %in% c("dead", "deceased", "event"), 1,
                ifelse(direct %in% c("alive", "censored"), 0, NA))

  idx <- is.na(out)
  if (any(idx)) {
    num_status <- extract_text_group("(?i)status[^0-9]*[:=]\\s*([01])", x[idx])
    out[idx] <- ifelse(num_status == "1", 1,
                       ifelse(num_status == "0", 0, NA))
  }

  idx <- is.na(out)
  if (any(idx)) {
    vals <- trimws(gsub("^.*?:\\s*", "", x[idx]))
    out[idx] <- ifelse(vals %in% c("dead", "deceased", "event", "1"), 1,
                       ifelse(vals %in% c("alive", "censored", "0"), 0, NA))
  }
  suppressWarnings(as.numeric(out))
}

parse_sex_field <- function(x) {
  x <- tolower(as.character(x))
  out <- extract_text_group("(?i)\\b(?:sex|gender)\\s*[:=]\\s*([a-z]+)", x)
  out <- ifelse(out %in% c("m", "male"), "MALE",
                ifelse(out %in% c("f", "female"), "FEMALE", NA))
  out
}

best_parsed_column <- function(sample_meta, parser_fun, context_regex = NULL, min_non_na = 5) {
  if (nrow(sample_meta) == 0) {
    return(list(values = rep(NA, 0), column = NA_character_, notes = NA_character_))
  }
  cn <- colnames(sample_meta)
  candidate_idx <- seq_along(cn)
  if (!is.null(context_regex)) {
    cn_low <- tolower(cn)
    value_match <- vapply(sample_meta, function(v) {
      vals <- tolower(as.character(v))
      any(grepl(context_regex, vals, perl = TRUE))
    }, logical(1))
    candidate_idx <- which(grepl(context_regex, cn_low, perl = TRUE) | value_match)
    if (length(candidate_idx) == 0) {
      candidate_idx <- seq_along(cn)
    }
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
      if (grepl("month|months|\\bmo\\b", tolower(cn[idx])) || grepl("month|months|\\bmo\\b", example_blob)) {
        best_note <- "months"
      } else if (grepl("day|days", tolower(cn[idx])) || grepl("day|days", example_blob)) {
        best_note <- "days"
      } else {
        best_note <- "unknown_unit"
      }
    }
  }

  list(values = best_values, column = best_col, notes = best_note)
}

zscore_vector <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    rep(0, length(x))
  } else {
    (x - m) / s
  }
}

calc_time_roc <- function(surv_time, surv_status, marker, time_point) {
  if (requireNamespace("timeROC", quietly = TRUE)) {
    out <- tryCatch(
      timeROC::timeROC(
        T = surv_time,
        delta = surv_status,
        marker = marker,
        cause = 1,
        times = time_point,
        iid = FALSE
      ),
      error = function(e) NULL
    )
    if (!is.null(out)) {
      auc <- as.numeric(out$AUC[1])
      return(list(method = "timeROC", auc = auc, fp = NULL, tp = NULL))
    }
  }

  n <- length(stats::na.omit(surv_time))
  span_value <- 0.25 * (n^(-0.20))
  out <- tryCatch(
    survivalROC::survivalROC(
      Stime = surv_time,
      status = surv_status,
      marker = marker,
      predict.time = time_point,
      method = "NNE",
      span = span_value
    ),
    error = function(e) NULL
  )
  if (is.null(out)) {
    return(list(method = "failed", auc = NA_real_, fp = NULL, tp = NULL))
  }
  list(method = "survivalROC", auc = as.numeric(out$AUC), fp = out$FP, tp = out$TP)
}

plot_to_pdf_png <- function(base_path, width, height, plot_fun) {
  pdf_path <- paste0(base_path, ".pdf")
  png_path <- paste0(base_path, ".png")

  grDevices::pdf(pdf_path, width = width, height = height)
  plot_fun()
  grDevices::dev.off()

  grDevices::png(png_path, width = width, height = height, units = "in", res = 300)
  plot_fun()
  grDevices::dev.off()
}

choose_multiv_formula <- function(dt, include_dataset = FALSE) {
  covars <- character(0)
  if (include_dataset && "dataset_id" %in% colnames(dt) && data.table::uniqueN(dt$dataset_id) >= 2) {
    covars <- c(covars, "dataset_id")
  }

  if ("age" %in% colnames(dt) && sum(!is.na(dt$age)) >= max(20, floor(0.5 * nrow(dt)))) {
    covars <- c(covars, "age")
  }
  if ("sex" %in% colnames(dt)) {
    sex_non_na <- na.omit(as.character(dt$sex))
    if (length(sex_non_na) >= max(20, floor(0.5 * nrow(dt))) && data.table::uniqueN(sex_non_na) >= 2) {
      covars <- c(covars, "sex")
    }
  }

  if (length(covars) == 0) {
    rhs <- "risk_score"
  } else {
    rhs <- paste(c("risk_score", covars), collapse = " + ")
  }

  list(
    formula = as.formula(paste("Surv(OS_time, OS_status) ~", rhs)),
    label = rhs
  )
}

evaluate_dataset <- function(dt, dataset_id, cutoff_value, cutoff_label, include_dataset = FALSE) {
  dt <- copy(dt)
  dt <- dt[!is.na(OS_time) & !is.na(OS_status) & OS_time > 0 & !is.na(risk_score)]
  if (nrow(dt) == 0) {
    return(list(
      performance = data.table(
        dataset_id = dataset_id,
        sample_count = 0,
        event_count = 0,
        all_6_signature_genes_available = TRUE,
        risk_score_standardization = "cohort-wise gene-level z-score",
        cutoff_type = cutoff_label,
        cutoff_value = cutoff_value,
        high_risk_samples = 0,
        low_risk_samples = 0,
        logrank_P.Value = NA_real_,
        univariate_cox_HR = NA_real_,
        univariate_cox_lower95 = NA_real_,
        univariate_cox_upper95 = NA_real_,
        univariate_cox_P.Value = NA_real_,
        multivariate_model_used = "not_fitted",
        multivariate_cox_HR = NA_real_,
        multivariate_cox_lower95 = NA_real_,
        multivariate_cox_upper95 = NA_real_,
        multivariate_cox_P.Value = NA_real_,
        C_index = NA_real_,
        AUC_1year = NA_real_,
        AUC_2year = NA_real_,
        AUC_3year = NA_real_,
        ROC_method_used = NA_character_
      ),
      risk_dt = dt,
      fit_km = NULL,
      roc_list = list()
    ))
  }

  dt[, risk_group := ifelse(risk_score >= cutoff_value, "high_risk", "low_risk")]
  dt[, risk_group := factor(risk_group, levels = c("low_risk", "high_risk"))]

  surv_obj <- Surv(dt$OS_time, dt$OS_status)
  fit_km <- NULL
  logrank_p <- NA_real_
  if (data.table::uniqueN(dt$risk_group) == 2) {
    fit_km <- survfit(surv_obj ~ risk_group, data = dt)
    logrank <- survdiff(surv_obj ~ risk_group, data = dt)
    logrank_p <- 1 - pchisq(logrank$chisq, df = 1)
  }

  uni_fit <- coxph(surv_obj ~ risk_score, data = dt)
  uni_sum <- summary(uni_fit)
  uni_hr <- as.numeric(uni_sum$coefficients["risk_score", "exp(coef)"])
  uni_p <- as.numeric(uni_sum$coefficients["risk_score", "Pr(>|z|)"])
  uni_low <- as.numeric(uni_sum$conf.int["risk_score", "lower .95"])
  uni_up <- as.numeric(uni_sum$conf.int["risk_score", "upper .95"])

  multiv_choice <- choose_multiv_formula(dt, include_dataset = include_dataset)
  multiv_fit <- coxph(multiv_choice$formula, data = dt)
  multiv_sum <- summary(multiv_fit)
  if ("risk_score" %in% rownames(multiv_sum$coefficients)) {
    multiv_hr <- as.numeric(multiv_sum$coefficients["risk_score", "exp(coef)"])
    multiv_p <- as.numeric(multiv_sum$coefficients["risk_score", "Pr(>|z|)"])
    multiv_low <- as.numeric(multiv_sum$conf.int["risk_score", "lower .95"])
    multiv_up <- as.numeric(multiv_sum$conf.int["risk_score", "upper .95"])
  } else {
    multiv_hr <- NA_real_
    multiv_p <- NA_real_
    multiv_low <- NA_real_
    multiv_up <- NA_real_
  }

  c_index_obj <- survival::concordance(surv_obj ~ risk_score, data = dt, reverse = TRUE)
  c_index <- as.numeric(c_index_obj$concordance)

  roc_1 <- calc_time_roc(dt$OS_time, dt$OS_status, dt$risk_score, 365)
  roc_2 <- calc_time_roc(dt$OS_time, dt$OS_status, dt$risk_score, 730)
  roc_3 <- calc_time_roc(dt$OS_time, dt$OS_status, dt$risk_score, 1095)
  roc_method_used <- unique(c(roc_1$method, roc_2$method, roc_3$method))
  roc_method_used <- roc_method_used[roc_method_used != "failed"][1]

  performance <- data.table(
    dataset_id = dataset_id,
    sample_count = nrow(dt),
    event_count = sum(dt$OS_status == 1, na.rm = TRUE),
    all_6_signature_genes_available = TRUE,
    risk_score_standardization = "cohort-wise gene-level z-score",
    cutoff_type = cutoff_label,
    cutoff_value = cutoff_value,
    high_risk_samples = sum(dt$risk_group == "high_risk"),
    low_risk_samples = sum(dt$risk_group == "low_risk"),
    logrank_P.Value = logrank_p,
    univariate_cox_HR = uni_hr,
    univariate_cox_lower95 = uni_low,
    univariate_cox_upper95 = uni_up,
    univariate_cox_P.Value = uni_p,
    multivariate_model_used = multiv_choice$label,
    multivariate_cox_HR = multiv_hr,
    multivariate_cox_lower95 = multiv_low,
    multivariate_cox_upper95 = multiv_up,
    multivariate_cox_P.Value = multiv_p,
    C_index = c_index,
    AUC_1year = roc_1$auc,
    AUC_2year = roc_2$auc,
    AUC_3year = roc_3$auc,
    ROC_method_used = roc_method_used
  )

  list(
    performance = performance,
    risk_dt = dt,
    fit_km = fit_km,
    roc_list = list(`365` = roc_1, `730` = roc_2, `1095` = roc_3)
  )
}

make_km_plot <- function(fit_km, dt, logrank_p, title_text) {
  cols <- c(low_risk = "#225EA8", high_risk = "#D7301F")
  plot(fit_km, col = cols, lwd = 2, xlab = "Days", ylab = "Overall survival probability", main = pretty_label(title_text))
  legend("bottomleft", legend = pretty_label(c("low_risk", "high_risk")), col = cols, lwd = 2, bty = "n")
  text(
    x = max(dt$OS_time, na.rm = TRUE) * 0.6,
    y = 0.15,
    labels = paste0("Log-rank P = ", format(logrank_p, digits = 3, scientific = TRUE))
  )
}

make_risk_distribution_plot <- function(dt, title_text) {
  plot_dt <- copy(dt)
  setorder(plot_dt, risk_score)
  plot_dt[, sample_rank := seq_len(.N)]
  p <- ggplot(plot_dt, aes(x = sample_rank, y = risk_score, color = risk_group)) +
    geom_point(size = 2) +
    geom_hline(yintercept = unique(plot_dt$cutoff_value_used), linetype = 2) +
    scale_color_manual(
      values = c(low_risk = "#225EA8", high_risk = "#D7301F"),
      labels = pretty_label(c("low_risk", "high_risk"))
    ) +
    labs(x = "Samples ranked by risk score", y = "Risk score", color = NULL, title = pretty_label(title_text)) +
    theme_bw(base_size = 12)
  print(p)
}

make_time_roc_plot <- function(roc_list, title_text) {
  plot(0:1, 0:1, type = "n", xlab = "False positive rate", ylab = "True positive rate", main = pretty_label(title_text))
  abline(0, 1, lty = 2, col = "grey60")
  legend_items <- character(0)
  legend_cols <- character(0)

  if (!is.null(roc_list[["365"]]$fp)) {
    lines(roc_list[["365"]]$fp, roc_list[["365"]]$tp, col = "#1B9E77", lwd = 2)
    legend_items <- c(legend_items, paste0("1-year AUC = ", sprintf("%.3f", roc_list[["365"]]$auc)))
    legend_cols <- c(legend_cols, "#1B9E77")
  }
  if (!is.null(roc_list[["730"]]$fp)) {
    lines(roc_list[["730"]]$fp, roc_list[["730"]]$tp, col = "#D95F02", lwd = 2)
    legend_items <- c(legend_items, paste0("2-year AUC = ", sprintf("%.3f", roc_list[["730"]]$auc)))
    legend_cols <- c(legend_cols, "#D95F02")
  }
  if (!is.null(roc_list[["1095"]]$fp)) {
    lines(roc_list[["1095"]]$fp, roc_list[["1095"]]$tp, col = "#7570B3", lwd = 2)
    legend_items <- c(legend_items, paste0("3-year AUC = ", sprintf("%.3f", roc_list[["1095"]]$auc)))
    legend_cols <- c(legend_cols, "#7570B3")
  }
  if (length(legend_items) > 0) {
    legend("bottomright", legend = legend_items, col = legend_cols, lwd = 2, bty = "n")
  }
}

gpl570_df <- parse_gpl_annotation("00_raw_data/geo_validation/GPL570_family.soft.gz")

coverage_rows <- list()
main_rows <- list()
sensitivity_rows <- list()
summary_rows <- list()
risk_data_list <- list()
km_sensitivity_plots <- list()

for (i in seq_len(nrow(target_datasets))) {
  dataset_id_i <- target_datasets$dataset_id[i]
  gse_i <- target_datasets$gse[i]
  matrix_path <- file.path("00_raw_data/geo_validation", target_datasets$matrix_file[i])
  parsed <- parse_series_matrix(matrix_path)
  expr_mat <- parsed$expr
  sample_meta <- parsed$sample_meta

  mapping_dt <- build_probe_mapping(expr_mat, gpl570_df)
  sig_map <- mapping_dt[gene_symbol %in% signature_genes & probe_rank == 1]
  sig_map <- sig_map[match(signature_genes, gene_symbol)]

  available_vec <- signature_genes %in% sig_map$gene_symbol
  coverage_rows[[length(coverage_rows) + 1]] <- data.table(
    dataset_id = dataset_id_i,
    signature_gene = signature_genes,
    available = available_vec,
    probe_id_or_feature_id = sapply(signature_genes, function(g) {
      x <- sig_map[gene_symbol == g]$probe_id
      if (length(x) == 0) NA_character_ else x[1]
    }),
    mapping_source = "GPL570_annotation_IQR_max_probe"
  )

  missing_sig <- signature_genes[!available_vec]
  if (length(missing_sig) > 0) {
    summary_rows[[length(summary_rows) + 1]] <- data.table(
      dataset_id = dataset_id_i,
      validation_success = FALSE,
      notes = paste("Missing signature genes:", paste(missing_sig, collapse = ";"))
    )
    next
  }

  probe_ids <- sig_map$probe_id
  expr_sig <- expr_mat[probe_ids, , drop = FALSE]
  rownames(expr_sig) <- sig_map$gene_symbol
  expr_sig <- expr_sig[signature_genes, , drop = FALSE]

  expr_sig_z <- t(apply(expr_sig, 1, zscore_vector))
  risk_score <- rowSums(t(expr_sig_z) * coef_named[signature_genes])

  age_info <- best_parsed_column(sample_meta, parse_age_field, context_regex = "age")
  if (sum(!is.na(age_info$values)) > 0 && stats::median(age_info$values, na.rm = TRUE) > 120) {
    age_info$values <- age_info$values / 365.25
  }
  sex_info <- best_parsed_column(sample_meta, parse_sex_field, context_regex = "sex|gender|male|female")
  os_time_info <- best_parsed_column(sample_meta, parse_os_time_field, context_regex = "overall survival|survival|\\bos\\b|days to death|follow")
  os_status_info <- best_parsed_column(sample_meta, parse_os_status_field, context_regex = "status|alive|dead|deceased|vital|life")

  if (identical(os_time_info$notes, "months")) {
    os_time_info$values <- os_time_info$values * 30.44
  }

  clin_dt <- data.table(
    dataset_id = dataset_id_i,
    sample_id = if ("sample_id" %in% colnames(sample_meta)) as.character(sample_meta$sample_id) else colnames(expr_sig),
    patient_id = if ("geo_accession" %in% colnames(sample_meta)) as.character(sample_meta$geo_accession) else if ("sample_id" %in% colnames(sample_meta)) as.character(sample_meta$sample_id) else colnames(expr_sig),
    OS_time = os_time_info$values,
    OS_status = os_status_info$values,
    age = age_info$values,
    sex = sex_info$values,
    risk_score = as.numeric(risk_score)
  )
  clin_dt <- clin_dt[!is.na(OS_time) & !is.na(OS_status) & OS_time > 0]

  if (nrow(clin_dt) == 0) {
    summary_rows[[length(summary_rows) + 1]] <- data.table(
      dataset_id = dataset_id_i,
      validation_success = FALSE,
      notes = "No samples with complete OS_time and OS_status after filtering."
    )
    next
  }

  main_eval <- evaluate_dataset(clin_dt, dataset_id_i, tcga_cutoff, "TCGA_median_cutoff", include_dataset = FALSE)
  main_eval$risk_dt[, cutoff_value_used := tcga_cutoff]
  cohort_cutoff <- stats::median(main_eval$risk_dt$risk_score, na.rm = TRUE)
  sens_eval <- evaluate_dataset(clin_dt, dataset_id_i, cohort_cutoff, "cohort_median_sensitivity", include_dataset = FALSE)

  risk_out <- copy(main_eval$risk_dt)
  risk_out[, risk_group_tcga_cutoff := as.character(risk_group)]
  risk_out[, risk_group_cohort_median := as.character(ifelse(risk_score >= cohort_cutoff, "high_risk", "low_risk"))]
  fwrite(risk_out, file.path("08_validation", paste0(dataset_id_i, "_external_validation_risk_score.csv")))
  fwrite(main_eval$performance, file.path("08_validation", paste0(dataset_id_i, "_external_validation_performance.csv")))

  plot_to_pdf_png(
    file.path("13_figures", paste0("Figure7_", dataset_id_i, "_KM_TCGAcutoff")),
    7, 6,
    function() make_km_plot(main_eval$fit_km, main_eval$risk_dt, main_eval$performance$logrank_P.Value[1], paste0(dataset_id_i, " KM (TCGA cutoff)"))
  )
  plot_to_pdf_png(
    file.path("13_figures", paste0("Figure7_", dataset_id_i, "_timeROC")),
    7, 6,
    function() make_time_roc_plot(main_eval$roc_list, paste0(dataset_id_i, " time-dependent ROC"))
  )
  plot_to_pdf_png(
    file.path("13_figures", paste0("Figure7_", dataset_id_i, "_risk_score_distribution")),
    8, 6,
    function() make_risk_distribution_plot(main_eval$risk_dt, paste0(dataset_id_i, " risk score distribution"))
  )

  main_rows[[length(main_rows) + 1]] <- main_eval$performance
  sensitivity_rows[[length(sensitivity_rows) + 1]] <- data.table(
    dataset_id = dataset_id_i,
    cutoff_type = "cohort_median_sensitivity",
    cutoff_value = cohort_cutoff,
    high_risk_samples = sens_eval$performance$high_risk_samples[1],
    low_risk_samples = sens_eval$performance$low_risk_samples[1],
    logrank_P.Value = sens_eval$performance$logrank_P.Value[1]
  )
  summary_rows[[length(summary_rows) + 1]] <- data.table(
    dataset_id = dataset_id_i,
    validation_success = TRUE,
    notes = paste(
      paste0("OS_time_column=", os_time_info$column),
      paste0("OS_status_column=", os_status_info$column),
      paste0("age_column=", age_info$column),
      paste0("sex_column=", sex_info$column),
      sep = "; "
    )
  )
  risk_data_list[[dataset_id_i]] <- risk_out
  km_sensitivity_plots[[dataset_id_i]] <- list(
    fit = sens_eval$fit_km,
    dt = sens_eval$risk_dt,
    p = sens_eval$performance$logrank_P.Value[1]
  )
}

coverage_dt <- rbindlist(coverage_rows, fill = TRUE)
fwrite(coverage_dt, "08_validation/external_validation_signature_gene_coverage.csv")

main_dt <- if (length(main_rows) > 0) rbindlist(main_rows, fill = TRUE) else data.table()
sens_dt <- if (length(sensitivity_rows) > 0) rbindlist(sensitivity_rows, fill = TRUE) else data.table()
summary_dt <- if (length(summary_rows) > 0) rbindlist(summary_rows, fill = TRUE) else data.table()

combined_success <- all(target_datasets$dataset_id %in% names(risk_data_list))
if (combined_success) {
  combined_dt <- rbindlist(risk_data_list, fill = TRUE)
  combined_dt[, dataset_id_factor := factor(dataset_id)]
  combined_eval <- evaluate_dataset(combined_dt, "combined_GPL570", tcga_cutoff, "TCGA_median_cutoff", include_dataset = TRUE)
  combined_eval$risk_dt[, cutoff_value_used := tcga_cutoff]
  combined_cohort_cutoff <- stats::median(combined_eval$risk_dt$risk_score, na.rm = TRUE)
  combined_sens_eval <- evaluate_dataset(combined_dt, "combined_GPL570", combined_cohort_cutoff, "cohort_median_sensitivity", include_dataset = TRUE)

  combined_risk_out <- copy(combined_eval$risk_dt)
  combined_risk_out[, risk_group_tcga_cutoff := as.character(risk_group)]
  combined_risk_out[, risk_group_cohort_median := as.character(ifelse(risk_score >= combined_cohort_cutoff, "high_risk", "low_risk"))]
  fwrite(combined_risk_out, "08_validation/combined_GPL570_external_validation_risk_score.csv")
  fwrite(combined_eval$performance, "08_validation/combined_GPL570_external_validation_performance.csv")

  plot_to_pdf_png(
    file.path("13_figures", "Figure7_combined_GPL570_KM_TCGAcutoff"),
    7, 6,
    function() make_km_plot(combined_eval$fit_km, combined_eval$risk_dt, combined_eval$performance$logrank_P.Value[1], "Combined GPL570 KM (TCGA cutoff)")
  )
  plot_to_pdf_png(
    file.path("13_figures", "Figure7_combined_GPL570_timeROC"),
    7, 6,
    function() make_time_roc_plot(combined_eval$roc_list, "Combined GPL570 time-dependent ROC")
  )
  plot_to_pdf_png(
    file.path("13_figures", "Figure7_combined_GPL570_risk_score_distribution"),
    8, 6,
    function() make_risk_distribution_plot(combined_eval$risk_dt, "Combined GPL570 risk score distribution")
  )

  main_dt <- rbindlist(list(main_dt, combined_eval$performance), fill = TRUE)
  sens_dt <- rbindlist(list(
    sens_dt,
    data.table(
      dataset_id = "combined_GPL570",
      cutoff_type = "cohort_median_sensitivity",
      cutoff_value = combined_cohort_cutoff,
      high_risk_samples = combined_sens_eval$performance$high_risk_samples[1],
      low_risk_samples = combined_sens_eval$performance$low_risk_samples[1],
      logrank_P.Value = combined_sens_eval$performance$logrank_P.Value[1]
    )
  ), fill = TRUE)
  summary_dt <- rbindlist(list(
    summary_dt,
    data.table(
      dataset_id = "combined_GPL570",
      validation_success = TRUE,
      notes = "Combined analysis with dataset_id adjustment."
    )
  ), fill = TRUE)
  km_sensitivity_plots[["combined_GPL570"]] <- list(
    fit = combined_sens_eval$fit_km,
    dt = combined_sens_eval$risk_dt,
    p = combined_sens_eval$performance$logrank_P.Value[1]
  )
}

fwrite(main_dt, "14_tables/external_validation_main_results.csv")
fwrite(summary_dt, "08_validation/external_validation_summary.csv")
fwrite(sens_dt, "08_validation/external_validation_sensitivity_cohort_median.csv")

manuscript_dt <- copy(main_dt)
fwrite(manuscript_dt, "14_tables/external_validation_summary_for_manuscript.csv")

if ("GSE37642_GPL570" %in% main_dt$dataset_id) {
  fwrite(main_dt[dataset_id == "GSE37642_GPL570"], "08_validation/GSE37642_GPL570_external_validation_performance.csv")
}
if ("GSE12417_GPL570" %in% main_dt$dataset_id) {
  fwrite(main_dt[dataset_id == "GSE12417_GPL570"], "08_validation/GSE12417_GPL570_external_validation_performance.csv")
}
if ("combined_GPL570" %in% main_dt$dataset_id) {
  fwrite(main_dt[dataset_id == "combined_GPL570"], "08_validation/combined_GPL570_external_validation_performance.csv")
}

if (length(km_sensitivity_plots) > 0) {
  plot_to_pdf_png(
    file.path("13_figures", "FigureS_external_validation_KM_cohort_median_sensitivity"),
    width = 6 * length(km_sensitivity_plots),
    height = 5,
    function() {
      op <- par(mfrow = c(1, length(km_sensitivity_plots)))
      on.exit(par(op), add = TRUE)
      for (nm in names(km_sensitivity_plots)) {
        item <- km_sensitivity_plots[[nm]]
        make_km_plot(item$fit, item$dt, item$p, paste0(nm, " cohort median sensitivity"))
      }
    }
  )
}

save_session_info("16_logs/sessionInfo_17_external_validation_cross_platform_signature_geo.txt")
