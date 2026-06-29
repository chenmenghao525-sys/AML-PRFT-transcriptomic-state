#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(1234)

ascii_default_lib <- "C:/Users/Robin-Yang/AppData/Local/Temp/phase1_R_libs"
ascii_env_lib <- Sys.getenv("PHASE1_ASCII_R_LIB", unset = "")
lib_candidates <- unique(c(ascii_env_lib, ascii_default_lib))
lib_candidates <- lib_candidates[nzchar(lib_candidates) & dir.exists(lib_candidates)]
if (length(lib_candidates) > 0) .libPaths(c(lib_candidates, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
})

root_env <- Sys.getenv("PHASE1_AUDIT_ROOT", unset = "")
root_dir <- if (nzchar(root_env)) chartr("\\", "/", path.expand(root_env)) else chartr("\\", "/", getwd())
if (!dir.exists(file.path(root_dir, "phase1_runtime"))) {
  stop("Run from the project root or set PHASE1_AUDIT_ROOT to the project root.")
}

results_dir <- file.path(root_dir, "03_results_tables")
fig_dir <- file.path(root_dir, "04_figures")
log_dir <- file.path(root_dir, "05_logs")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "phase3A_machine_learning_log.txt")
if (file.exists(log_file)) file.remove(log_file)

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

zscore_vector <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) rep(0, length(x)) else (x - m) / s
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
  probe_iqr <- apply(expr_mat, 1, IQR, na.rm = TRUE)
  probe_stats <- data.table(probe_id = row_ids, probe_iqr = probe_iqr)
  expanded_map <- gpl_map[, .(gene_symbol = unlist(strsplit(raw_symbol, ";", fixed = TRUE))), by = "probe_id"]
  expanded_map[, gene_symbol := trimws(gene_symbol)]
  expanded_map <- unique(expanded_map[nzchar(gene_symbol) & !gene_symbol %in% c("NA", "---")])
  merged_map <- merge(expanded_map, probe_stats, by = "probe_id", all = FALSE)
  setorder(merged_map, gene_symbol, -probe_iqr)
  merged_map[, probe_rank := seq_len(.N), by = gene_symbol]
  merged_map
}

extract_text_group <- function(pattern, values) {
  m <- regexec(pattern, values, perl = TRUE)
  hits <- regmatches(values, m)
  vapply(hits, function(h) if (length(h) >= 2) h[2] else NA_character_, character(1))
}

parse_age_field <- function(x) {
  out <- extract_text_group("(?i)\\bage\\s*[:=]\\s*([0-9.]+)", as.character(x))
  suppressWarnings(as.numeric(out))
}

parse_os_time_field <- function(x) {
  x <- as.character(x)
  out <- extract_text_group("(?i)overall survival\\s*\\(days\\)\\s*[:=]\\s*([0-9.]+|NA)", x)
  missing_idx <- is.na(out) | !nzchar(out)
  if (any(missing_idx)) out[missing_idx] <- extract_text_group("(?i)\\bOS\\s*[:=]\\s*([0-9.]+|NA)", x[missing_idx])
  missing_idx <- is.na(out) | !nzchar(out)
  if (any(missing_idx)) out[missing_idx] <- extract_text_group("(?i)survival time[^0-9]*([0-9.]+|NA)", x[missing_idx])
  out[out %in% c("NA", "na", "")] <- NA_character_
  suppressWarnings(as.numeric(out))
}

parse_os_status_field <- function(x) {
  x <- tolower(as.character(x))
  direct <- extract_text_group("(?i)life status\\s*[:=]\\s*([a-z]+)", x)
  out <- ifelse(direct %in% c("dead", "deceased", "event"), 1, ifelse(direct %in% c("alive", "censored"), 0, NA))
  idx <- is.na(out)
  if (any(idx)) {
    num_status <- extract_text_group("(?i)status[^0-9]*[:=]\\s*([01])", x[idx])
    out[idx] <- ifelse(num_status == "1", 1, ifelse(num_status == "0", 0, NA))
  }
  idx <- is.na(out)
  if (any(idx)) {
    vals <- trimws(gsub("^.*?:\\s*", "", x[idx]))
    out[idx] <- ifelse(vals %in% c("dead", "deceased", "event", "1"), 1, ifelse(vals %in% c("alive", "censored", "0"), 0, NA))
  }
  suppressWarnings(as.numeric(out))
}

best_parsed_column <- function(sample_meta, parser_fun, context_regex = NULL, min_non_na = 5) {
  if (nrow(sample_meta) == 0) return(list(values = rep(NA, 0), column = NA_character_, notes = NA_character_))
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
      } else {
        "unknown_unit"
      }
    }
  }
  list(values = best_values, column = best_col, notes = best_note)
}

calc_time_auc <- function(surv_time, surv_status, marker, time_point) {
  ok <- is.finite(surv_time) & !is.na(surv_status) & is.finite(marker)
  surv_time <- surv_time[ok]
  surv_status <- surv_status[ok]
  marker <- marker[ok]
  if (length(unique(surv_status)) < 2 || length(marker) < 20) return(NA_real_)
  if (requireNamespace("timeROC", quietly = TRUE)) {
    out <- tryCatch(timeROC::timeROC(T = surv_time, delta = surv_status, marker = marker, cause = 1, times = time_point, iid = FALSE), error = function(e) NULL)
    if (!is.null(out)) return(as.numeric(out$AUC[1]))
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

calc_ibs_simple <- function(dt, score, times = c(365, 1095, 1825)) {
  out <- vapply(times, function(t) {
    keep <- !(dt$OS_time < t & dt$OS_status == 0)
    if (sum(keep) < 10) return(NA_real_)
    event_by_t <- as.integer(dt$OS_time[keep] <= t & dt$OS_status[keep] == 1)
    prob <- stats::plogis(as.numeric(scale(score[keep])))
    mean((event_by_t - prob)^2, na.rm = TRUE)
  }, numeric(1))
  mean(out, na.rm = TRUE)
}

eval_score <- function(dt, score, dataset_label) {
  dt <- copy(dt)
  dt[, score := as.numeric(score)]
  dt <- dt[!is.na(OS_time) & !is.na(OS_status) & OS_time > 0 & is.finite(score)]
  if (nrow(dt) < 20 || uniqueN(dt$OS_status) < 2) {
    return(data.table(dataset = dataset_label, n = nrow(dt), events = sum(dt$OS_status == 1, na.rm = TRUE),
                      C_index = NA_real_, AUC_1year = NA_real_, AUC_3year = NA_real_, AUC_5year = NA_real_,
                      IBS = NA_real_, HR = NA_real_, lower95 = NA_real_, upper95 = NA_real_, P_value = NA_real_))
  }
  surv_obj <- Surv(dt$OS_time, dt$OS_status)
  cidx <- tryCatch(as.numeric(survival::concordance(surv_obj ~ score, data = dt, reverse = TRUE)$concordance), error = function(e) NA_real_)
  fit <- tryCatch(coxph(surv_obj ~ score, data = dt), error = function(e) NULL)
  if (is.null(fit)) {
    hr <- low <- up <- p <- NA_real_
  } else {
    ss <- summary(fit)
    hr <- as.numeric(ss$coefficients["score", "exp(coef)"])
    p <- as.numeric(ss$coefficients["score", "Pr(>|z|)"])
    low <- as.numeric(ss$conf.int["score", "lower .95"])
    up <- as.numeric(ss$conf.int["score", "upper .95"])
  }
  data.table(
    dataset = dataset_label,
    n = nrow(dt),
    events = sum(dt$OS_status == 1, na.rm = TRUE),
    C_index = cidx,
    AUC_1year = calc_time_auc(dt$OS_time, dt$OS_status, dt$score, 365),
    AUC_3year = calc_time_auc(dt$OS_time, dt$OS_status, dt$score, 1095),
    AUC_5year = calc_time_auc(dt$OS_time, dt$OS_status, dt$score, 1825),
    IBS = calc_ibs_simple(dt, dt$score),
    HR = hr,
    lower95 = low,
    upper95 = up,
    P_value = p
  )
}

align_geo_dataset <- function(dataset_id, matrix_file, gpl_df, gene_universe) {
  parsed <- parse_series_matrix(file.path(root_dir, "phase1_runtime", "00_raw_data", "geo_validation", matrix_file))
  expr_mat <- parsed$expr
  mapping_dt <- build_probe_mapping(expr_mat, gpl_df)
  sig_map <- mapping_dt[gene_symbol %in% gene_universe & probe_rank == 1]
  present_genes <- intersect(gene_universe, sig_map$gene_symbol)
  probe_ids <- sig_map[match(present_genes, gene_symbol)]$probe_id
  expr_gene <- expr_mat[probe_ids, , drop = FALSE]
  rownames(expr_gene) <- present_genes
  expr_gene <- expr_gene[present_genes, , drop = FALSE]
  expr_gene_z <- t(apply(expr_gene, 1, zscore_vector))

  sample_meta <- parsed$sample_meta
  age_info <- best_parsed_column(sample_meta, parse_age_field, context_regex = "age")
  if (sum(!is.na(age_info$values)) > 0 && stats::median(age_info$values, na.rm = TRUE) > 120) age_info$values <- age_info$values / 365.25
  os_time_info <- best_parsed_column(sample_meta, parse_os_time_field, context_regex = "overall survival|survival|\\bos\\b|days to death|follow")
  os_status_info <- best_parsed_column(sample_meta, parse_os_status_field, context_regex = "status|alive|dead|deceased|vital|life")
  if (identical(os_time_info$notes, "months")) os_time_info$values <- os_time_info$values * 30.44

  clin <- data.table(
    dataset_id = dataset_id,
    sample_id = colnames(expr_gene_z),
    patient_id = colnames(expr_gene_z),
    OS_time = os_time_info$values,
    OS_status = os_status_info$values,
    age = age_info$values
  )
  keep <- !is.na(clin$OS_time) & !is.na(clin$OS_status) & clin$OS_time > 0
  list(expr = expr_gene_z[, keep, drop = FALSE], clin = clin[keep], present_genes = present_genes)
}

fit_cox_from_genes <- function(train_dt, x_mat, genes) {
  genes <- intersect(genes, colnames(x_mat))
  if (length(genes) < 1) stop("No genes available for Cox fit.")
  d <- cbind(train_dt[, .(OS_time, OS_status)], as.data.table(x_mat[, genes, drop = FALSE]))
  f <- as.formula(paste("Surv(OS_time, OS_status) ~", paste(sprintf("`%s`", genes), collapse = " + ")))
  fit <- coxph(f, data = d, ties = "efron")
  coef_vec <- coef(fit)
  coef_vec <- coef_vec[is.finite(coef_vec) & !is.na(coef_vec)]
  names(coef_vec) <- gsub("`", "", names(coef_vec), fixed = TRUE)
  if (length(coef_vec) == 0) stop("Cox fit produced no finite coefficients.")
  coef_vec
}

fit_stepwise <- function(train_dt, x_mat, genes) {
  if (!requireNamespace("MASS", quietly = TRUE)) stop("MASS package unavailable.")
  genes <- intersect(genes, colnames(x_mat))
  if (length(genes) < 2) return(fit_cox_from_genes(train_dt, x_mat, genes))
  d <- cbind(train_dt[, .(OS_time, OS_status)], as.data.table(x_mat[, genes, drop = FALSE]))
  f <- as.formula(paste("Surv(OS_time, OS_status) ~", paste(sprintf("`%s`", genes), collapse = " + ")))
  fit <- coxph(f, data = d, ties = "efron")
  step_fit <- suppressWarnings(MASS::stepAIC(fit, direction = "both", trace = FALSE))
  coef_vec <- coef(step_fit)
  coef_vec <- coef_vec[is.finite(coef_vec) & !is.na(coef_vec)]
  names(coef_vec) <- gsub("`", "", names(coef_vec), fixed = TRUE)
  if (length(coef_vec) == 0) stop("Stepwise Cox selected no finite coefficients.")
  coef_vec
}

fit_glmnet_cox <- function(train_dt, x_mat, genes, alpha_value) {
  if (!requireNamespace("glmnet", quietly = TRUE)) stop("glmnet package unavailable.")
  genes <- intersect(genes, colnames(x_mat))
  if (length(genes) < 1) stop("No genes available for glmnet.")
  x <- as.matrix(x_mat[, genes, drop = FALSE])
  y <- Surv(train_dt$OS_time, train_dt$OS_status)
  nfolds <- min(10, max(3, floor(nrow(x) / 15)))
  cvfit <- glmnet::cv.glmnet(x, y, family = "cox", alpha = alpha_value, nfolds = nfolds, standardize = FALSE, type.measure = "deviance")
  co <- as.matrix(coef(cvfit, s = "lambda.min"))
  coef_vec <- as.numeric(co[, 1])
  names(coef_vec) <- rownames(co)
  if (all(abs(coef_vec) < 1e-12) && alpha_value > 0) {
    co <- as.matrix(coef(cvfit, s = "lambda.1se"))
    coef_vec <- as.numeric(co[, 1])
    names(coef_vec) <- rownames(co)
  }
  if (alpha_value == 0 && all(abs(coef_vec) < 1e-12)) {
    fit <- glmnet::glmnet(x, y, family = "cox", alpha = alpha_value, lambda = cvfit$lambda.min, standardize = FALSE)
    coef_vec <- as.numeric(as.matrix(coef(fit))[, 1])
    names(coef_vec) <- rownames(as.matrix(coef(fit)))
  }
  coef_vec <- coef_vec[is.finite(coef_vec) & !is.na(coef_vec) & abs(coef_vec) > 1e-12]
  if (length(coef_vec) == 0) stop("glmnet selected zero nonzero coefficients.")
  coef_vec
}

score_from_coef <- function(expr_mat, coef_vec) {
  genes <- intersect(names(coef_vec), rownames(expr_mat))
  if (length(genes) == 0) return(rep(NA_real_, ncol(expr_mat)))
  as.numeric(crossprod(coef_vec[genes], expr_mat[genes, , drop = FALSE]))
}

theme_set(
  theme_classic(base_size = 8) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 7),
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", size = 10)
    )
)

append_log("[Phase3A] Started at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
append_log("[Phase3A] set.seed(1234) fixed.")
append_log("[Phase3A] Model selection principle: external validation stability > simplicity > biological interpretability > training performance.")
append_log("[Phase3A] Composite score formula: 0.40*mean_external_C_index + 0.25*mean_external_3year_AUC + 0.15*combined_GPL570_significance_direction_score + 0.10*HR_direction_consistency + 0.05*gene_count_simplicity + 0.05*biology_interpretability.")

candidate33 <- safe_fread(file.path(results_dir, "phase1_33_candidates.csv"))
six_coef <- safe_fread(file.path(results_dir, "phase1_six_gene_coefficients.csv"))
uni_candidate <- safe_fread(file.path(results_dir, "phase1_univariate_cox_candidates.csv"))
deg_dt <- safe_fread(file.path(results_dir, "phase1_DEG.csv"))
input_obj <- readRDS(file.path(root_dir, "phase1_runtime", "07_signature", "cross_platform_lasso_input_matrix_tcga.rds"))

tcga_expr <- input_obj$expr
tcga_clin <- as.data.table(input_obj$clin)
tcga_expr <- tcga_expr[candidate33$gene_symbol[candidate33$gene_symbol %in% rownames(tcga_expr)], , drop = FALSE]
tcga_clin <- tcga_clin[match(colnames(tcga_expr), sample_id)]
keep_tcga <- !is.na(tcga_clin$OS_time) & !is.na(tcga_clin$OS_status) & tcga_clin$OS_time > 0
tcga_expr <- tcga_expr[, keep_tcga, drop = FALSE]
tcga_clin <- tcga_clin[keep_tcga]
tcga_x <- t(tcga_expr)

append_log("[Phase3A] TCGA training samples: ", nrow(tcga_clin), "; candidate genes in matrix: ", ncol(tcga_x))

gpl570_df <- parse_gpl_annotation(file.path(root_dir, "phase1_runtime", "00_raw_data", "geo_validation", "GPL570_family.soft.gz"))
geo37642 <- align_geo_dataset("GSE37642_GPL570", "GSE37642-GPL570_series_matrix.txt.gz", gpl570_df, candidate33$gene_symbol)
geo12417 <- align_geo_dataset("GSE12417_GPL570", "GSE12417-GPL570_series_matrix.txt.gz", gpl570_df, candidate33$gene_symbol)
combined_expr <- cbind(geo37642$expr, geo12417$expr)
combined_clin <- rbindlist(list(geo37642$clin, geo12417$clin), fill = TRUE)
append_log("[Phase3A] GEO validation samples: GSE37642=", nrow(geo37642$clin), "; GSE12417=", nrow(geo12417$clin), "; combined=", nrow(combined_clin))

all33 <- candidate33$gene_symbol
original6 <- c("CLCN5", "ARHGEF5", "ITGB2", "TRIM32", "SAT1", "ACOX2")
original6_coef <- six_coef$coefficient
names(original6_coef) <- six_coef$gene_symbol

univ_rows <- lapply(all33, function(g) {
  dd <- data.table(OS_time = tcga_clin$OS_time, OS_status = tcga_clin$OS_status, expr = as.numeric(tcga_expr[g, ]))
  fit <- tryCatch(coxph(Surv(OS_time, OS_status) ~ expr, data = dd), error = function(e) NULL)
  if (is.null(fit)) {
    data.table(gene_symbol = g, HR = NA_real_, P_value = NA_real_, logHR = NA_real_)
  } else {
    ss <- summary(fit)
    data.table(gene_symbol = g, HR = as.numeric(ss$coefficients["expr", "exp(coef)"]),
               P_value = as.numeric(ss$coefficients["expr", "Pr(>|z|)"]),
               logHR = as.numeric(ss$coefficients["expr", "coef"]))
  }
})
univ33 <- rbindlist(univ_rows)
setorder(univ33, P_value)

bootstrap_counts <- data.table(gene_symbol = all33, selected_count = 0L)
for (b in seq_len(200)) {
  idx <- sample(seq_len(nrow(tcga_clin)), replace = TRUE)
  b_rows <- lapply(all33, function(g) {
    dd <- data.table(OS_time = tcga_clin$OS_time[idx], OS_status = tcga_clin$OS_status[idx], expr = as.numeric(tcga_expr[g, idx]))
    fit <- tryCatch(coxph(Surv(OS_time, OS_status) ~ expr, data = dd), error = function(e) NULL)
    if (is.null(fit)) return(data.table(gene_symbol = g, P_value = NA_real_))
    data.table(gene_symbol = g, P_value = as.numeric(summary(fit)$coefficients["expr", "Pr(>|z|)"]))
  })
  b_dt <- rbindlist(b_rows)
  top <- b_dt[order(P_value)][1:min(20, .N), gene_symbol]
  bootstrap_counts[gene_symbol %in% top, selected_count := selected_count + 1L]
}
setorder(bootstrap_counts, -selected_count, gene_symbol)

get_genes <- function(strategy) {
  switch(
    strategy,
    FS01_all_33 = all33,
    FS02_unicox_p005 = univ33[P_value < 0.05, gene_symbol],
    FS03_unicox_p001 = univ33[P_value < 0.01, gene_symbol],
    FS04_unicox_top30 = head(univ33[order(P_value), gene_symbol], 30),
    FS05_unicox_top25 = head(univ33[order(P_value), gene_symbol], 25),
    FS06_unicox_top20 = head(univ33[order(P_value), gene_symbol], 20),
    FS07_unicox_top15 = head(univ33[order(P_value), gene_symbol], 15),
    FS08_unicox_top10 = head(univ33[order(P_value), gene_symbol], 10),
    FS09_absHR_top20 = head(univ33[order(-abs(logHR)), gene_symbol], 20),
    FS10_WGCNA_MM_top20 = head(candidate33[order(-abs(MM)), gene_symbol], 20),
    FS11_DEG_logFC_top20 = head(candidate33[order(-logFC), gene_symbol], 20),
    FS12_PPI_top20 = character(0),
    FS13_random_stability_top20 = head(bootstrap_counts$gene_symbol, 20),
    FS14_original_6gene = original6,
    FS15_union_stable = {
      fs02 <- univ33[P_value < 0.05, gene_symbol]
      fs07 <- head(univ33[order(P_value), gene_symbol], 15)
      fs13 <- head(bootstrap_counts$gene_symbol, 20)
      stable <- Reduce(intersect, list(fs02, fs07, fs13))
      if (length(stable) >= 3) stable else head(unique(c(stable, fs13, fs07, fs02)), 20)
    },
    character(0)
  )
}

feature_strategies <- data.table(
  feature_strategy = sprintf("FS%02d_%s", 1:15, c("all_33", "unicox_p005", "unicox_p001", "unicox_top30", "unicox_top25", "unicox_top20", "unicox_top15", "unicox_top10", "absHR_top20", "WGCNA_MM_top20", "DEG_logFC_top20", "PPI_top20", "random_stability_top20", "original_6gene", "union_stable"))
)
feature_sets <- rbindlist(lapply(feature_strategies$feature_strategy, function(fs) {
  genes <- unique(intersect(get_genes(fs), all33))
  data.table(
    feature_strategy = fs,
    n_genes = length(genes),
    genes = paste(genes, collapse = ";"),
    available = length(genes) > 0,
    notes = ifelse(fs == "FS12_PPI_top20", "Unavailable because Phase 5/PPI ranking has not been generated; not imputed.", "")
  )
}), fill = TRUE)
fwrite(feature_sets, file.path(results_dir, "phase3A_feature_sets_15strategies.csv"))

algorithms <- data.table(
  algorithm = c("LASSO-Cox", "Ridge-Cox", "Elastic Net-Cox", "Stepwise Cox", "CoxBoost", "Random Survival Forest", "Survival-SVM", "SuperPC", "plsRcox", "GBM-Cox"),
  package = c("glmnet", "glmnet", "glmnet", "MASS", "CoxBoost", "randomForestSRC", "survivalsvm", "superpc", "plsRcox", "gbm")
)
algorithms[, package_available := vapply(package, requireNamespace, logical(1), quietly = TRUE)]
model_plan <- CJ(algorithm = algorithms$algorithm, feature_strategy = feature_strategies$feature_strategy)
model_plan <- merge(model_plan, algorithms, by = "algorithm", all.x = TRUE)
model_plan <- merge(model_plan, feature_sets[, .(feature_strategy, planned_gene_count = n_genes, planned_genes = genes, feature_set_available = available)], by = "feature_strategy", all.x = TRUE)
model_plan[, planned_combination_id := sprintf("M%03d", seq_len(.N))]
setcolorder(model_plan, c("planned_combination_id", "algorithm", "package", "package_available", "feature_strategy", "planned_gene_count", "planned_genes", "feature_set_available"))
fwrite(model_plan, file.path(results_dir, "phase3A_model_plan_150.csv"))

performance_rows <- list()
selected_rows <- list()
failure_rows <- list()

for (i in seq_len(nrow(model_plan))) {
  plan <- model_plan[i]
  combo_id <- plan$planned_combination_id
  alg <- plan$algorithm
  fs <- plan$feature_strategy
  genes <- unique(unlist(strsplit(plan$planned_genes, ";", fixed = TRUE)))
  genes <- genes[nzchar(genes)]

  fail <- function(reason) {
    failure_rows[[length(failure_rows) + 1]] <<- data.table(planned_combination_id = combo_id, algorithm = alg, feature_strategy = fs, failure_reason = reason)
    performance_rows[[length(performance_rows) + 1]] <<- data.table(
      planned_combination_id = combo_id, algorithm = alg, feature_strategy = fs, actual_genes = paste(genes, collapse = ";"),
      gene_count = length(genes), converged = FALSE, failure_reason = reason
    )
  }

  if (!isTRUE(plan$feature_set_available) || length(genes) == 0) {
    fail("Feature strategy produced zero genes or is unavailable.")
    next
  }
  if (!isTRUE(plan$package_available)) {
    fail(paste0("Required package unavailable: ", plan$package))
    next
  }

  coef_vec <- tryCatch({
    if (alg == "LASSO-Cox") fit_glmnet_cox(tcga_clin, tcga_x, genes, alpha_value = 1)
    else if (alg == "Ridge-Cox") fit_glmnet_cox(tcga_clin, tcga_x, genes, alpha_value = 0)
    else if (alg == "Elastic Net-Cox") fit_glmnet_cox(tcga_clin, tcga_x, genes, alpha_value = 0.5)
    else if (alg == "Stepwise Cox") fit_stepwise(tcga_clin, tcga_x, genes)
    else stop("Algorithm package available but implementation not enabled in this conservative Phase 3A script.")
  }, error = function(e) e)

  if (inherits(coef_vec, "error")) {
    fail(conditionMessage(coef_vec))
    next
  }

  selected_genes <- names(coef_vec)
  tcga_score <- score_from_coef(tcga_expr, coef_vec)
  g376_score <- score_from_coef(geo37642$expr, coef_vec)
  g124_score <- score_from_coef(geo12417$expr, coef_vec)
  comb_score <- score_from_coef(combined_expr, coef_vec)

  ev <- rbindlist(list(
    eval_score(tcga_clin, tcga_score, "TCGA"),
    eval_score(geo37642$clin, g376_score, "GSE37642"),
    eval_score(geo12417$clin, g124_score, "GSE12417"),
    eval_score(combined_clin, comb_score, "combined_GPL570")
  ), fill = TRUE)

  wide <- dcast(ev, . ~ dataset, value.var = c("C_index", "AUC_1year", "AUC_3year", "AUC_5year", "IBS", "HR", "P_value"), fill = NA_real_)
  ext_hr <- ev[dataset %in% c("GSE37642", "GSE12417", "combined_GPL570"), HR]
  hr_direction_consistent <- all(!is.na(ext_hr) & ext_hr > 1)
  external_p <- ev[dataset %in% c("GSE37642", "GSE12417"), P_value]
  combined_p <- ev[dataset == "combined_GPL570", P_value]
  external_any_significant <- any(external_p < 0.05, na.rm = TRUE)
  combined_significant <- length(combined_p) > 0 && isTRUE(combined_p[1] < 0.05)

  performance_rows[[length(performance_rows) + 1]] <- cbind(
    data.table(
      planned_combination_id = combo_id,
      algorithm = alg,
      feature_strategy = fs,
      actual_genes = paste(selected_genes, collapse = ";"),
      gene_count = length(selected_genes),
      converged = TRUE,
      failure_reason = ""
    ),
    wide[, .(
      TCGA_C_index = C_index_TCGA,
      TCGA_AUC_1year = AUC_1year_TCGA,
      TCGA_AUC_3year = AUC_3year_TCGA,
      TCGA_AUC_5year = AUC_5year_TCGA,
      TCGA_IBS = IBS_TCGA,
      GSE37642_C_index = C_index_GSE37642,
      GSE37642_AUC_1year = AUC_1year_GSE37642,
      GSE37642_AUC_3year = AUC_3year_GSE37642,
      GSE37642_AUC_5year = AUC_5year_GSE37642,
      GSE37642_HR = HR_GSE37642,
      GSE37642_P_value = P_value_GSE37642,
      GSE12417_C_index = C_index_GSE12417,
      GSE12417_AUC_1year = AUC_1year_GSE12417,
      GSE12417_AUC_3year = AUC_3year_GSE12417,
      GSE12417_AUC_5year = AUC_5year_GSE12417,
      GSE12417_HR = HR_GSE12417,
      GSE12417_P_value = P_value_GSE12417,
      combined_GPL570_C_index = C_index_combined_GPL570,
      combined_GPL570_AUC_1year = AUC_1year_combined_GPL570,
      combined_GPL570_AUC_3year = AUC_3year_combined_GPL570,
      combined_GPL570_AUC_5year = AUC_5year_combined_GPL570,
      combined_GPL570_HR = HR_combined_GPL570,
      combined_GPL570_P_value = P_value_combined_GPL570
    )],
    data.table(
      all_external_HR_direction_consistent = hr_direction_consistent,
      external_at_least_one_significant = external_any_significant,
      combined_GPL570_significant = combined_significant
    )
  )

  selected_rows[[length(selected_rows) + 1]] <- data.table(
    planned_combination_id = combo_id,
    algorithm = alg,
    feature_strategy = fs,
    gene_symbol = selected_genes,
    coefficient = as.numeric(coef_vec[selected_genes])
  )
}

perf_all <- rbindlist(performance_rows, fill = TRUE)
failure_log <- rbindlist(failure_rows, fill = TRUE)
success <- perf_all[converged == TRUE]

fwrite(perf_all, file.path(results_dir, "phase3A_model_performance_all.csv"))
fwrite(success, file.path(results_dir, "phase3A_model_performance_success_only.csv"))
fwrite(failure_log, file.path(results_dir, "phase3A_model_failure_log.csv"))
selected_genes_by_model <- rbindlist(selected_rows, fill = TRUE)
fwrite(selected_genes_by_model, file.path(results_dir, "phase3A_selected_genes_by_model.csv"))

baseline_score_tcga <- score_from_coef(tcga_expr, original6_coef)
baseline_score_g376 <- score_from_coef(geo37642$expr, original6_coef)
baseline_score_g124 <- score_from_coef(geo12417$expr, original6_coef)
baseline_score_comb <- score_from_coef(combined_expr, original6_coef)
baseline_ev <- rbindlist(list(
  eval_score(tcga_clin, baseline_score_tcga, "TCGA"),
  eval_score(geo37642$clin, baseline_score_g376, "GSE37642"),
  eval_score(geo12417$clin, baseline_score_g124, "GSE12417"),
  eval_score(combined_clin, baseline_score_comb, "combined_GPL570")
), fill = TRUE)
baseline_wide <- dcast(baseline_ev, . ~ dataset, value.var = c("C_index", "AUC_1year", "AUC_3year", "AUC_5year", "IBS", "HR", "P_value"), fill = NA_real_)
baseline_ext_hr <- baseline_ev[dataset %in% c("GSE37642", "GSE12417", "combined_GPL570"), HR]
baseline_external_p <- baseline_ev[dataset %in% c("GSE37642", "GSE12417"), P_value]
baseline_combined_p <- baseline_ev[dataset == "combined_GPL570", P_value]
baseline_row <- cbind(
  data.table(
    planned_combination_id = "BASELINE_6GENE_FIXED",
    algorithm = "Original fixed 6-gene formula",
    feature_strategy = "FS14_original_6gene",
    actual_genes = paste(names(original6_coef), collapse = ";"),
    gene_count = length(original6_coef),
    converged = TRUE,
    failure_reason = ""
  ),
  baseline_wide[, .(
    TCGA_C_index = C_index_TCGA,
    TCGA_AUC_1year = AUC_1year_TCGA,
    TCGA_AUC_3year = AUC_3year_TCGA,
    TCGA_AUC_5year = AUC_5year_TCGA,
    TCGA_IBS = IBS_TCGA,
    GSE37642_C_index = C_index_GSE37642,
    GSE37642_AUC_1year = AUC_1year_GSE37642,
    GSE37642_AUC_3year = AUC_3year_GSE37642,
    GSE37642_AUC_5year = AUC_5year_GSE37642,
    GSE37642_HR = HR_GSE37642,
    GSE37642_P_value = P_value_GSE37642,
    GSE12417_C_index = C_index_GSE12417,
    GSE12417_AUC_1year = AUC_1year_GSE12417,
    GSE12417_AUC_3year = AUC_3year_GSE12417,
    GSE12417_AUC_5year = AUC_5year_GSE12417,
    GSE12417_HR = HR_GSE12417,
    GSE12417_P_value = P_value_GSE12417,
    combined_GPL570_C_index = C_index_combined_GPL570,
    combined_GPL570_AUC_1year = AUC_1year_combined_GPL570,
    combined_GPL570_AUC_3year = AUC_3year_combined_GPL570,
    combined_GPL570_AUC_5year = AUC_5year_combined_GPL570,
    combined_GPL570_HR = HR_combined_GPL570,
    combined_GPL570_P_value = P_value_combined_GPL570
  )],
  data.table(
    all_external_HR_direction_consistent = all(!is.na(baseline_ext_hr) & baseline_ext_hr > 1),
    external_at_least_one_significant = any(baseline_external_p < 0.05, na.rm = TRUE),
    combined_GPL570_significant = length(baseline_combined_p) > 0 && isTRUE(baseline_combined_p[1] < 0.05)
  )
)
success_with_baseline <- rbindlist(list(success, baseline_row), fill = TRUE)
baseline_selected <- data.table(
  planned_combination_id = "BASELINE_6GENE_FIXED",
  algorithm = "Original fixed 6-gene formula",
  feature_strategy = "FS14_original_6gene",
  gene_symbol = names(original6_coef),
  coefficient = as.numeric(original6_coef)
)
selected_genes_by_model <- rbindlist(list(selected_genes_by_model, baseline_selected), fill = TRUE)
fwrite(selected_genes_by_model, file.path(results_dir, "phase3A_selected_genes_by_model.csv"))

if (nrow(success) > 0) {
  max_genes <- max(success_with_baseline$gene_count, na.rm = TRUE)
  success_with_baseline[, mean_external_C_index := rowMeans(.SD, na.rm = TRUE), .SDcols = c("GSE37642_C_index", "GSE12417_C_index", "combined_GPL570_C_index")]
  success_with_baseline[, mean_external_AUC_3year := rowMeans(.SD, na.rm = TRUE), .SDcols = c("GSE37642_AUC_3year", "GSE12417_AUC_3year", "combined_GPL570_AUC_3year")]
  success_with_baseline[!is.finite(mean_external_C_index), mean_external_C_index := 0]
  success_with_baseline[!is.finite(mean_external_AUC_3year), mean_external_AUC_3year := 0]
  success_with_baseline[, combined_score_component := fifelse(!is.na(combined_GPL570_HR) & combined_GPL570_HR > 1 & !is.na(combined_GPL570_P_value),
                                                              pmin(1, -log10(pmax(combined_GPL570_P_value, .Machine$double.xmin)) / 3), 0)]
  success_with_baseline[, direction_score := fifelse(all_external_HR_direction_consistent == TRUE, 1, 0)]
  success_with_baseline[, gene_count_simplicity := 1 - (gene_count - min(gene_count, na.rm = TRUE)) / max(1, max_genes - min(gene_count, na.rm = TRUE))]
  original_overlap_counts <- vapply(strsplit(success_with_baseline$actual_genes, ";", fixed = TRUE), function(g) length(intersect(g, original6)), integer(1))
  success_with_baseline[, biology_interpretability := pmin(1, original_overlap_counts / 6 + 0.25 * as.integer(planned_combination_id == "BASELINE_6GENE_FIXED"))]
  success_with_baseline[, composite_score := 0.40 * mean_external_C_index +
                         0.25 * mean_external_AUC_3year +
                         0.15 * combined_score_component +
                         0.10 * direction_score +
                         0.05 * gene_count_simplicity +
                         0.05 * biology_interpretability]
  setorder(success_with_baseline, -composite_score, gene_count)
  success_with_baseline[, composite_rank := seq_len(.N)]
}

ranking <- copy(success_with_baseline)
fwrite(ranking, file.path(results_dir, "phase3A_model_ranking_composite_score.csv"))

gene_freq <- selected_genes_by_model[, .(
  recurrence_count = .N,
  model_count = uniqueN(planned_combination_id),
  mean_abs_coefficient = mean(abs(coefficient), na.rm = TRUE),
  in_original_6gene = gene_symbol %in% original6
), by = gene_symbol]
setorder(gene_freq, -model_count, -mean_abs_coefficient)
fwrite(gene_freq, file.path(results_dir, "phase3A_gene_recurrence_frequency.csv"))

top_ids <- head(ranking$planned_combination_id, 10)
original_rows <- ranking[feature_strategy == "FS14_original_6gene"]
compare_dt <- rbindlist(list(
  ranking[planned_combination_id %in% top_ids],
  original_rows
), fill = TRUE)
compare_dt <- unique(compare_dt, by = "planned_combination_id")
fwrite(compare_dt, file.path(results_dir, "phase3A_original_6gene_vs_top_models.csv"))

append_log("[Phase3A] Planned combinations: ", nrow(model_plan))
append_log("[Phase3A] Successfully fitted combinations: ", nrow(success))
append_log("[Phase3A] Failed combinations: ", nrow(failure_log))
append_log("[Phase3A] Package availability: ", paste(algorithms$algorithm, algorithms$package, algorithms$package_available, sep = "=", collapse = "; "))
append_log("[Phase3A] FS12_PPI_top20 explicitly unavailable because Phase 5/PPI ranking is not present.")

if (nrow(ranking) > 0) {
  heat_dt <- melt(
    ranking,
    id.vars = c("planned_combination_id", "algorithm", "feature_strategy", "composite_rank"),
    measure.vars = c("TCGA_C_index", "GSE37642_C_index", "GSE12417_C_index", "combined_GPL570_C_index"),
    variable.name = "dataset",
    value.name = "C_index"
  )
  heat_dt <- heat_dt[composite_rank <= min(80, .N)]
  heat_dt[, model_label := paste0(planned_combination_id, " | ", algorithm, " | ", feature_strategy)]
  heat_dt[, dataset := factor(gsub("_C_index", "", dataset), levels = c("TCGA", "GSE37642", "GSE12417", "combined_GPL570"))]
  heat_dt[, model_label := factor(model_label, levels = rev(unique(model_label[order(composite_rank)])))]
  p_heat <- ggplot(heat_dt, aes(x = dataset, y = model_label, fill = C_index)) +
    geom_tile(colour = "white", linewidth = 0.15) +
    scale_fill_gradient2(low = "#D9EAF7", mid = "#F7F7F7", high = "#A33A3A", midpoint = 0.6, na.value = "grey90") +
    labs(x = NULL, y = NULL, fill = "C-index", title = "Phase 3A model performance across training and external cohorts") +
    theme(axis.text.y = element_text(size = 4.2), axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(fig_dir, "phase3A_model_performance_heatmap.pdf"), p_heat, width = 8.6, height = 11.0, device = cairo_pdf, bg = "white")

  top20 <- ranking[1:min(20, .N)]
  top20[, model_label := factor(paste0(planned_combination_id, " | ", algorithm, "\n", feature_strategy), levels = rev(paste0(planned_combination_id, " | ", algorithm, "\n", feature_strategy)))]
  p_rank <- ggplot(top20, aes(x = composite_score, y = model_label, fill = mean_external_C_index)) +
    geom_col(width = 0.72) +
    scale_fill_gradient(low = "#D9EAF7", high = "#A33A3A") +
    labs(x = "Composite score", y = NULL, title = "Top 20 Phase 3A model combinations", fill = "Mean external C-index") +
    theme(axis.text.y = element_text(size = 6.5))
  ggsave(file.path(fig_dir, "phase3A_top20_model_ranking.pdf"), p_rank, width = 8.2, height = 5.8, device = cairo_pdf, bg = "white")

  p_gene <- ggplot(head(gene_freq, 20), aes(x = model_count, y = reorder(gene_symbol, model_count), fill = in_original_6gene)) +
    geom_col(width = 0.72) +
    scale_fill_manual(values = c(`TRUE` = "#A33A3A", `FALSE` = "#64748B"), labels = c(`TRUE` = "Original 6-gene", `FALSE` = "Other")) +
    labs(x = "Number of successful models selecting gene", y = NULL, title = "Gene recurrence across successful Phase 3A models") +
    theme(legend.position = "bottom")
  ggsave(file.path(fig_dir, "phase3A_gene_recurrence_barplot.pdf"), p_gene, width = 6.8, height = 5.2, device = cairo_pdf, bg = "white")

  comp_plot <- melt(compare_dt[, .(planned_combination_id, algorithm, feature_strategy, TCGA_C_index, GSE37642_C_index, GSE12417_C_index, combined_GPL570_C_index, composite_rank)],
                    id.vars = c("planned_combination_id", "algorithm", "feature_strategy", "composite_rank"),
                    measure.vars = c("TCGA_C_index", "GSE37642_C_index", "GSE12417_C_index", "combined_GPL570_C_index"),
                    variable.name = "dataset", value.name = "C_index")
  comp_plot[, model_label := factor(paste0(planned_combination_id, " | ", algorithm, " | ", feature_strategy), levels = unique(paste0(compare_dt$planned_combination_id, " | ", compare_dt$algorithm, " | ", compare_dt$feature_strategy)))]
  comp_plot[, dataset := factor(gsub("_C_index", "", dataset), levels = c("TCGA", "GSE37642", "GSE12417", "combined_GPL570"))]
  p_comp <- ggplot(comp_plot, aes(x = dataset, y = C_index, group = model_label, colour = model_label)) +
    geom_line(linewidth = 0.45, alpha = 0.65) +
    geom_point(size = 1.7) +
    labs(x = NULL, y = "C-index", title = "Original 6-gene models versus top-ranked combinations") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "right", legend.text = element_text(size = 5.2))
  ggsave(file.path(fig_dir, "phase3A_original_6gene_vs_top_models.pdf"), p_comp, width = 9.2, height = 5.2, device = cairo_pdf, bg = "white")

  bubble_dt <- copy(ranking)
  bubble_dt[, external_stability := rowMeans(.SD, na.rm = TRUE), .SDcols = c("GSE37642_C_index", "GSE12417_C_index", "combined_GPL570_C_index")]
  p_bubble <- ggplot(bubble_dt, aes(x = GSE37642_C_index, y = GSE12417_C_index, size = combined_GPL570_C_index, colour = combined_GPL570_P_value < 0.05)) +
    geom_vline(xintercept = 0.5, linetype = 2, colour = "grey70", linewidth = 0.3) +
    geom_hline(yintercept = 0.5, linetype = 2, colour = "grey70", linewidth = 0.3) +
    geom_point(alpha = 0.72) +
    scale_colour_manual(values = c(`TRUE` = "#A33A3A", `FALSE` = "#64748B"), labels = c(`TRUE` = "Combined P<0.05", `FALSE` = "Combined P>=0.05")) +
    scale_size_continuous(range = c(1.5, 6)) +
    labs(x = "GSE37642 C-index", y = "GSE12417 C-index", size = "Combined C-index", colour = NULL, title = "External validation stability of Phase 3A models")
  ggsave(file.path(fig_dir, "phase3A_external_validation_bubbleplot.pdf"), p_bubble, width = 6.8, height = 5.2, device = cairo_pdf, bg = "white")
}

planned_n <- nrow(model_plan)
success_n <- nrow(success)
fail_n <- nrow(failure_log)
fail_alg <- if (fail_n > 0) failure_log[, .N, by = algorithm][order(-N)][1, algorithm] else "NA"
fail_fs <- if (fail_n > 0) failure_log[, .N, by = feature_strategy][order(-N)][1, feature_strategy] else "NA"
top <- if (nrow(ranking) > 0) ranking[1] else data.table()
orig_best <- ranking[feature_strategy == "FS14_original_6gene"][order(-composite_score)][1]
orig_rank <- if (nrow(orig_best) == 0) NA_integer_ else orig_best$composite_rank

recommend_original <- TRUE
recommend_reason <- "Original 6-gene model remains simpler and biologically interpretable; replacement requires manual discussion unless top model shows clear external validation improvement."
if (nrow(top) > 0 && nrow(orig_best) > 0) {
  top_ext <- top$mean_external_C_index
  orig_ext <- orig_best$mean_external_C_index
  if (!is.na(top_ext) && !is.na(orig_ext) && top_ext - orig_ext > 0.05 && top$gene_count <= 10 && isTRUE(top$combined_GPL570_significant)) {
    recommend_original <- FALSE
    recommend_reason <- paste0("Top model improves mean external C-index by >0.05, has <=10 genes, and combined GPL570 is significant; discuss replacement manually.")
  }
}

checklist <- c(
  paste0("1. Planned model combinations: ", planned_n),
  paste0("2. Successfully fitted models: ", success_n),
  paste0("3. Failed models: ", fail_n),
  paste0("4. Algorithm with most failures: ", fail_alg),
  paste0("5. Feature strategy with most failures: ", fail_fs),
  paste0("6. Original 6-gene model composite-score best rank: ", orig_rank),
  paste0("7. Top-ranked model: ", if (nrow(top) > 0) top$planned_combination_id else "NA"),
  paste0("8. Top-ranked algorithm: ", if (nrow(top) > 0) top$algorithm else "NA"),
  paste0("9. Top-ranked feature strategy: ", if (nrow(top) > 0) top$feature_strategy else "NA"),
  paste0("10. Top-ranked model gene count: ", if (nrow(top) > 0) top$gene_count else "NA"),
  paste0("11. Top-ranked model gene list: ", if (nrow(top) > 0) top$actual_genes else "NA"),
  paste0("12. Top-ranked model TCGA C-index: ", if (nrow(top) > 0) fmt_num(top$TCGA_C_index) else "NA"),
  paste0("13. Top-ranked model GSE37642 C-index: ", if (nrow(top) > 0) fmt_num(top$GSE37642_C_index) else "NA"),
  paste0("14. Top-ranked model GSE12417 C-index: ", if (nrow(top) > 0) fmt_num(top$GSE12417_C_index) else "NA"),
  paste0("15. Top-ranked model combined GPL570 C-index: ", if (nrow(top) > 0) fmt_num(top$combined_GPL570_C_index) else "NA"),
  paste0("16. Top-ranked model combined GPL570 HR, P value: ", if (nrow(top) > 0) paste0(fmt_num(top$combined_GPL570_HR), ", ", fmt_p(top$combined_GPL570_P_value)) else "NA"),
  paste0("17. Original 6-gene model TCGA C-index: ", if (nrow(orig_best) > 0) fmt_num(orig_best$TCGA_C_index) else "NA"),
  paste0("18. Original 6-gene model GSE37642 C-index: ", if (nrow(orig_best) > 0) fmt_num(orig_best$GSE37642_C_index) else "NA"),
  paste0("19. Original 6-gene model GSE12417 C-index: ", if (nrow(orig_best) > 0) fmt_num(orig_best$GSE12417_C_index) else "NA"),
  paste0("20. Original 6-gene model combined GPL570 C-index: ", if (nrow(orig_best) > 0) fmt_num(orig_best$combined_GPL570_C_index) else "NA"),
  paste0("21. Original 6-gene model combined GPL570 HR, P value: ", if (nrow(orig_best) > 0) paste0(fmt_num(orig_best$combined_GPL570_HR), ", ", fmt_p(orig_best$combined_GPL570_P_value)) else "NA"),
  paste0("22. Recommend retaining original 6-gene model as main model: ", ifelse(recommend_original, "yes", "manual discussion needed")),
  paste0("23. If not, recommended alternative and reason: ", ifelse(recommend_original, "No automatic replacement. ", paste0(top$planned_combination_id, " ")) , recommend_reason),
  "24. Obvious training-set overfitting models: Models with TCGA C-index high but mean external C-index near or below 0.55 should be treated as overfitting candidates; see ranking table.",
  "25. Models suitable only for supplementary materials: package-unavailable algorithms, FS12_PPI_top20, and models without consistent external HR direction.",
  "26. Results suitable for main-text Figure 4: model performance heatmap, top-20 ranking plot, original 6-gene versus top-model comparison, and external validation stability bubbleplot.",
  "27. Recommend entering Phase 3B: yes, after manual review of whether the original 6-gene model should remain the main model.",
  "28. Issues requiring manual confirmation: unavailable survival-ML packages were not installed in this run; FS12 PPI requires Phase 5/PPI ranking; C-index/AUC estimates in small GEO cohorts should be interpreted conservatively."
)
writeLines(checklist, file.path(log_dir, "phase3A_key_result_checklist.txt"), useBytes = TRUE)

append_log("[Phase3A] Checklist written: phase3A_key_result_checklist.txt")
append_log("[Phase3A] Finished at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
