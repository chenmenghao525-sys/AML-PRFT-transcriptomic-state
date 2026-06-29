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

root_env <- Sys.getenv("PHASE5_ROOT", unset = "")
root_dir <- if (nzchar(root_env)) chartr("\\", "/", path.expand(root_env)) else chartr("\\", "/", getwd())
if (!dir.exists(file.path(root_dir, "phase1_runtime"))) {
  stop("Run from the project root or set PHASE5_ROOT to the project root.")
}

results_dir <- file.path(root_dir, "03_results_tables")
fig_dir <- file.path(root_dir, "04_figures")
log_dir <- file.path(root_dir, "05_logs")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "phase3A_PPI_mini_update_log.txt")
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

calc_ibs_simple <- function(dt, score, times = c(365, 1095, 1825)) {
  if (requireNamespace("pec", quietly = TRUE)) {
    dt_pec <- data.frame(OS_time = dt$OS_time, OS_status = dt$OS_status, score = as.numeric(score))
    dt_pec <- dt_pec[is.finite(dt_pec$OS_time) & !is.na(dt_pec$OS_status) & is.finite(dt_pec$score), , drop = FALSE]
    if (nrow(dt_pec) >= 20 && length(unique(dt_pec$OS_status)) >= 2) {
      pec_val <- tryCatch({
        fit <- coxph(Surv(OS_time, OS_status) ~ score, data = dt_pec, x = TRUE, y = TRUE)
        max_time <- min(max(times), max(dt_pec$OS_time, na.rm = TRUE))
        times_use <- times[times <= max_time]
        if (length(times_use) == 0) times_use <- max_time
        pe <- pec::pec(list(score_model = fit), formula = Surv(OS_time, OS_status) ~ 1, data = dt_pec, times = times_use, exact = FALSE, cens.model = "marginal", splitMethod = "none", verbose = FALSE)
        ib <- tryCatch(pec::crps(pe), error = function(e) NULL)
        if (is.null(ib)) NA_real_ else as.numeric(ib[1, "score_model"])
      }, error = function(e) NA_real_)
      if (is.finite(pec_val)) return(pec_val)
    }
  }
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

make_train_data <- function(train_dt, x_mat, genes) {
  genes <- intersect(genes, colnames(x_mat))
  if (length(genes) < 1) stop("No genes available for model fitting.")
  d <- cbind(data.frame(OS_time = train_dt$OS_time, OS_status = train_dt$OS_status), as.data.frame(x_mat[, genes, drop = FALSE]))
  colnames(d) <- make.names(colnames(d))
  list(data = d, genes = genes, safe_genes = make.names(genes))
}

expr_to_newdata <- function(expr_mat, genes) {
  genes <- intersect(genes, rownames(expr_mat))
  nd <- as.data.frame(t(expr_mat[genes, , drop = FALSE]))
  colnames(nd) <- make.names(colnames(nd))
  nd
}

fit_coxboost_model <- function(train_dt, x_mat, genes) {
  if (!requireNamespace("CoxBoost", quietly = TRUE)) stop("CoxBoost package unavailable.")
  genes <- intersect(genes, colnames(x_mat))
  if (length(genes) < 1) stop("No genes available for CoxBoost.")
  x <- as.matrix(x_mat[, genes, drop = FALSE])
  maxstep <- min(100, max(20, length(genes) * 5))
  cv <- tryCatch(CoxBoost::cv.CoxBoost(time = train_dt$OS_time, status = train_dt$OS_status, x = x, maxstepno = maxstep, K = 5, type = "verweij"), error = function(e) NULL)
  stepno <- if (!is.null(cv) && !is.null(cv$optimal.step)) cv$optimal.step else min(50, maxstep)
  fit <- CoxBoost::CoxBoost(time = train_dt$OS_time, status = train_dt$OS_status, x = x, stepno = stepno)
  list(type = "coxboost", fit = fit, genes = genes, stepno = stepno)
}

fit_rsf_model <- function(train_dt, x_mat, genes) {
  genes <- intersect(genes, colnames(x_mat))
  if (length(genes) < 1) stop("No genes available for Random Survival Forest.")
  td <- make_train_data(train_dt, x_mat, genes)
  if (requireNamespace("randomForestSRC", quietly = TRUE)) {
    f <- as.formula(paste("Surv(OS_time, OS_status) ~", paste(td$safe_genes, collapse = " + ")))
    fit <- randomForestSRC::rfsrc(f, data = td$data, ntree = 300, nodesize = max(5, floor(nrow(td$data) / 20)), seed = 1234, forest = TRUE)
    return(list(type = "rfsrc", fit = fit, genes = td$genes))
  }
  if (requireNamespace("ranger", quietly = TRUE)) {
    f <- as.formula(paste("Surv(OS_time, OS_status) ~", paste(td$safe_genes, collapse = " + ")))
    fit <- ranger::ranger(f, data = td$data, num.trees = 300, mtry = max(1, floor(sqrt(length(td$genes)))), seed = 1234, write.forest = TRUE)
    return(list(type = "ranger_survival", fit = fit, genes = td$genes))
  }
  stop("Neither randomForestSRC nor ranger is available.")
}

fit_survivalsvm_model <- function(train_dt, x_mat, genes) {
  if (!requireNamespace("survivalsvm", quietly = TRUE)) stop("survivalsvm package unavailable.")
  td <- make_train_data(train_dt, x_mat, genes)
  f <- as.formula(paste("Surv(OS_time, OS_status) ~", paste(td$safe_genes, collapse = " + ")))
  fit <- survivalsvm::survivalsvm(f, data = td$data, type = "regression", gamma.mu = 1, opt.meth = "quadprog")
  list(type = "survivalsvm", fit = fit, genes = td$genes)
}

fit_superpc_model <- function(train_dt, x_mat, genes) {
  if (!requireNamespace("superpc", quietly = TRUE)) stop("superpc package unavailable.")
  genes <- intersect(genes, colnames(x_mat))
  if (length(genes) < 2) stop("SuperPC requires at least two genes.")
  sdata <- list(x = t(as.matrix(x_mat[, genes, drop = FALSE])), y = train_dt$OS_time, censoring.status = train_dt$OS_status, featurenames = genes)
  fit <- superpc::superpc.train(sdata, type = "survival")
  list(type = "superpc", fit = fit, train_data = sdata, genes = genes, threshold = 0, n.components = 1)
}

fit_gbm_cox_model <- function(train_dt, x_mat, genes) {
  if (!requireNamespace("gbm", quietly = TRUE)) stop("gbm package unavailable.")
  td <- make_train_data(train_dt, x_mat, genes)
  f <- as.formula(paste("Surv(OS_time, OS_status) ~", paste(td$safe_genes, collapse = " + ")))
  fit <- gbm::gbm(f, data = td$data, distribution = "coxph", n.trees = 500, interaction.depth = 2, shrinkage = 0.03, bag.fraction = 0.7, train.fraction = 1.0, verbose = FALSE)
  list(type = "gbm_cox", fit = fit, genes = td$genes, n.trees = 500)
}

fit_model_by_algorithm <- function(alg, train_dt, x_mat, genes) {
  if (alg == "LASSO-Cox") {
    coef_vec <- fit_glmnet_cox(train_dt, x_mat, genes, alpha_value = 1)
    return(list(type = "linear", coef = coef_vec, genes = names(coef_vec)))
  }
  if (alg == "Ridge-Cox") {
    coef_vec <- fit_glmnet_cox(train_dt, x_mat, genes, alpha_value = 0)
    return(list(type = "linear", coef = coef_vec, genes = names(coef_vec)))
  }
  if (alg == "Elastic Net-Cox") {
    coef_vec <- fit_glmnet_cox(train_dt, x_mat, genes, alpha_value = 0.5)
    return(list(type = "linear", coef = coef_vec, genes = names(coef_vec)))
  }
  if (alg == "Stepwise Cox") {
    coef_vec <- fit_stepwise(train_dt, x_mat, genes)
    return(list(type = "linear", coef = coef_vec, genes = names(coef_vec)))
  }
  if (alg == "CoxBoost") return(fit_coxboost_model(train_dt, x_mat, genes))
  if (alg == "Random Survival Forest") return(fit_rsf_model(train_dt, x_mat, genes))
  if (alg == "Survival-SVM") return(fit_survivalsvm_model(train_dt, x_mat, genes))
  if (alg == "SuperPC") return(fit_superpc_model(train_dt, x_mat, genes))
  if (alg == "GBM-Cox") return(fit_gbm_cox_model(train_dt, x_mat, genes))
  stop("Unknown or unsupported algorithm: ", alg)
}

predict_model_score <- function(model, expr_mat) {
  genes <- intersect(model$genes, rownames(expr_mat))
  if (length(genes) == 0) return(rep(NA_real_, ncol(expr_mat)))
  if (model$type == "linear") return(score_from_coef(expr_mat, model$coef))
  if (model$type == "coxboost") {
    x <- as.matrix(t(expr_mat[model$genes, , drop = FALSE]))
    sc <- predict(model$fit, newdata = x, type = "lp", at.step = model$stepno)
    return(as.numeric(sc))
  }
  if (model$type == "rfsrc") {
    nd <- expr_to_newdata(expr_mat, model$genes)
    pr <- predict(model$fit, newdata = nd)
    return(as.numeric(pr$predicted))
  }
  if (model$type == "ranger_survival") {
    nd <- expr_to_newdata(expr_mat, model$genes)
    pr <- predict(model$fit, data = nd)
    if (!is.null(pr$chf)) return(as.numeric(rowSums(pr$chf, na.rm = TRUE)))
    return(as.numeric(rowMeans(1 - pr$survival, na.rm = TRUE)))
  }
  if (model$type == "survivalsvm") {
    nd <- expr_to_newdata(expr_mat, model$genes)
    pr <- predict(model$fit, nd)
    return(as.numeric(pr$predicted))
  }
  if (model$type == "superpc") {
    xnew <- list(x = as.matrix(expr_mat[model$genes, , drop = FALSE]), y = rep(1, ncol(expr_mat)), censoring.status = rep(1, ncol(expr_mat)), featurenames = model$genes)
    pr <- superpc::superpc.predict(model$fit, model$train_data, xnew, threshold = model$threshold, n.components = model$n.components, prediction.type = "continuous")
    return(as.numeric(pr$v.pred))
  }
  if (model$type == "gbm_cox") {
    nd <- expr_to_newdata(expr_mat, model$genes)
    return(as.numeric(predict(model$fit, newdata = nd, n.trees = model$n.trees, type = "link")))
  }
  rep(NA_real_, ncol(expr_mat))
}

theme_set(
  theme_classic(base_size = 6.5) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 6.2),
      legend.position = "bottom",
      legend.title = element_text(size = 6.2),
      legend.text = element_text(size = 5.8),
      plot.title = element_text(face = "bold", size = 7),
      panel.grid = element_blank()
    )
)

save_pdf_checked <- function(filename, plot_obj, width, height) {
  ggsave(filename, plot_obj, width = width, height = height, device = cairo_pdf, bg = "white")
  con <- file(filename, open = "rb")
  on.exit(close(con), add = TRUE)
  hdr <- rawToChar(readBin(con, "raw", n = 5))
  if (!identical(hdr, "%PDF-")) stop("Generated file is not a valid PDF header: ", filename)
}

append_log("[Phase3A-PPI-mini] Started at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
append_log("[Phase3A-PPI-mini] set.seed(1234) fixed.")
append_log("[Phase3A-PPI-mini] This update only adds FS12_PPI_top20 as an additional network-derived feature strategy and does not rerun all Phase 3A combinations.")
append_log("[Phase3A-PPI-mini] Composite score formula preserved from Phase 3A-fix: 0.35*mean_external_C_index + 0.25*mean_external_3year_AUC + 0.15*combined_GPL570_significance_direction_score + 0.10*HR_direction_consistency + 0.05*gene_count_simplicity + 0.05*biology_interpretability + 0.05*anti_overfitting_score.")
append_log("[Phase3A-PPI-mini] Missing metric handling preserved: NA metrics are excluded by weight renormalization rather than replaced with zero.")
append_log("[Phase3A-PPI-mini] Figure contract: quantitative grid plus comparison panels; core conclusion is whether FS12_PPI_top20 provides complementary prognostic information relative to the retained original six-gene signature.")

required_inputs <- c(
  file.path(results_dir, "phase5_FS12_PPI_top20.csv"),
  file.path(results_dir, "phase5_consensus_hub_genes.csv"),
  file.path(results_dir, "phase5_PPI_hub_overlap_with_PRFT_modules.csv"),
  file.path(results_dir, "phase3A_fix_model_performance_success_only.csv"),
  file.path(results_dir, "phase3A_fix_model_ranking_composite_score.csv"),
  file.path(results_dir, "phase1_six_gene_coefficients.csv"),
  file.path(root_dir, "01_processed_data", "02_processed_data", "tcga_expr_clin_matched.rds"),
  file.path(root_dir, "phase1_runtime", "00_raw_data", "geo_validation", "GPL570_family.soft.gz"),
  file.path(root_dir, "phase1_runtime", "00_raw_data", "geo_validation", "GSE37642-GPL570_series_matrix.txt.gz"),
  file.path(root_dir, "phase1_runtime", "00_raw_data", "geo_validation", "GSE12417-GPL570_series_matrix.txt.gz")
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) stop("Missing required inputs: ", paste(missing_inputs, collapse = "; "))
append_log("[Phase3A-PPI-mini] Input files loaded from: ", paste(required_inputs, collapse = " | "))

ppi_top20 <- safe_fread(file.path(results_dir, "phase5_FS12_PPI_top20.csv"))
consensus_hubs <- safe_fread(file.path(results_dir, "phase5_consensus_hub_genes.csv"))
hub_overlap <- safe_fread(file.path(results_dir, "phase5_PPI_hub_overlap_with_PRFT_modules.csv"))
phase3A_success <- safe_fread(file.path(results_dir, "phase3A_fix_model_performance_success_only.csv"))
phase3A_ranking <- safe_fread(file.path(results_dir, "phase3A_fix_model_ranking_composite_score.csv"))
six_coef <- safe_fread(file.path(results_dir, "phase1_six_gene_coefficients.csv"))

fs12_genes <- unique(ppi_top20$gene_symbol)
append_log("[Phase3A-PPI-mini] FS12_PPI_top20 genes: ", paste(fs12_genes, collapse = ", "))
append_log("[Phase3A-PPI-mini] Consensus hub top10 from Phase 5: ", paste(head(consensus_hubs$gene_symbol, 10), collapse = ", "))
append_log("[Phase3A-PPI-mini] Read prior Phase3A-fix ranking successfully: ", ifelse(nrow(phase3A_ranking) > 0, "yes", "no"))
append_log("[Phase3A-PPI-mini] Read prior Phase3A-fix successful-model table successfully: ", ifelse(nrow(phase3A_success) > 0, "yes", "no"))

input_obj <- readRDS(file.path(root_dir, "01_processed_data", "02_processed_data", "tcga_expr_clin_matched.rds"))
tcga_expr <- input_obj$expr
tcga_clin <- as.data.table(input_obj$clin)
keep_tcga <- !is.na(tcga_clin$OS_time) & !is.na(tcga_clin$OS_status) & tcga_clin$OS_time > 0
tcga_expr <- tcga_expr[, keep_tcga, drop = FALSE]
tcga_clin <- tcga_clin[keep_tcga]
tcga_x <- t(tcga_expr)
append_log("[Phase3A-PPI-mini] TCGA source: 01_processed_data/02_processed_data/tcga_expr_clin_matched.rds")
append_log("[Phase3A-PPI-mini] TCGA survival-eligible samples used for FS12 modeling: ", nrow(tcga_clin))

gpl570_df <- parse_gpl_annotation(file.path(root_dir, "phase1_runtime", "00_raw_data", "geo_validation", "GPL570_family.soft.gz"))
geo37642 <- align_geo_dataset("GSE37642", "GSE37642-GPL570_series_matrix.txt.gz", gpl570_df, fs12_genes)
geo12417 <- align_geo_dataset("GSE12417", "GSE12417-GPL570_series_matrix.txt.gz", gpl570_df, fs12_genes)
combined_expr <- cbind(geo37642$expr, geo12417$expr)
combined_clin <- rbindlist(list(geo37642$clin, geo12417$clin), fill = TRUE)
append_log("[Phase3A-PPI-mini] GEO validation sources: GSE37642-GPL570, GSE12417-GPL570, combined GPL570.")
append_log("[Phase3A-PPI-mini] GEO survival-eligible samples: GSE37642=", nrow(geo37642$clin), "; GSE12417=", nrow(geo12417$clin), "; combined=", nrow(combined_clin))

extra_geo_files <- c(
  GSE6891 = file.path(root_dir, "phase1_runtime", "00_raw_data", "geo_validation", "GSE6891_series_matrix.txt.gz"),
  GSE14468 = file.path(root_dir, "phase1_runtime", "00_raw_data", "geo_validation", "GSE14468_series_matrix.txt.gz")
)

availability_rows <- list()
make_availability_row <- function(dataset_id, gene_universe, expr_genes, modeling_threshold = 10, notes = "") {
  present <- intersect(gene_universe, expr_genes)
  missing <- setdiff(gene_universe, expr_genes)
  data.table(
    dataset_id = dataset_id,
    total_FS12_genes = length(gene_universe),
    available_genes = length(present),
    missing_genes = length(missing),
    available_ratio = length(present) / max(1, length(gene_universe)),
    available_gene_list = paste(present, collapse = ";"),
    missing_gene_list = paste(missing, collapse = ";"),
    suitable_for_formal_modeling = length(present) >= modeling_threshold,
    notes = notes
  )
}
availability_rows[[length(availability_rows) + 1]] <- make_availability_row("TCGA-LAML", fs12_genes, rownames(tcga_expr), notes = "Training cohort; same OS_time > 0 filtering as Phase 3A-fix for modeling.")
availability_rows[[length(availability_rows) + 1]] <- make_availability_row("GSE37642_GPL570", fs12_genes, rownames(geo37642$expr), notes = "External validation cohort.")
availability_rows[[length(availability_rows) + 1]] <- make_availability_row("GSE12417_GPL570", fs12_genes, rownames(geo12417$expr), notes = "External validation cohort.")
availability_rows[[length(availability_rows) + 1]] <- make_availability_row("combined_GPL570", fs12_genes, rownames(combined_expr), notes = "Combined GPL570 external validation cohort.")
for (nm in names(extra_geo_files)) {
  if (file.exists(extra_geo_files[[nm]])) {
    parsed_extra <- tryCatch(parse_series_matrix(extra_geo_files[[nm]]), error = function(e) NULL)
    if (is.null(parsed_extra)) {
      availability_rows[[length(availability_rows) + 1]] <- data.table(
        dataset_id = nm,
        total_FS12_genes = length(fs12_genes),
        available_genes = NA_integer_,
        missing_genes = NA_integer_,
        available_ratio = NA_real_,
        available_gene_list = "",
        missing_gene_list = "",
        suitable_for_formal_modeling = FALSE,
        notes = "Local series matrix exists but parsing/annotation was not completed in this mini update; feasibility only."
      )
    } else {
      availability_rows[[length(availability_rows) + 1]] <- make_availability_row(
        nm,
        fs12_genes,
        rownames(parsed_extra$expr),
        notes = "Local series matrix exists; feasibility/availability check only, not used for formal modeling."
      )
    }
  } else {
    availability_rows[[length(availability_rows) + 1]] <- data.table(
      dataset_id = nm,
      total_FS12_genes = length(fs12_genes),
      available_genes = NA_integer_,
      missing_genes = NA_integer_,
      available_ratio = NA_real_,
      available_gene_list = "",
      missing_gene_list = "",
      suitable_for_formal_modeling = FALSE,
      notes = "Local dataset not found."
    )
  }
}
availability_dt <- rbindlist(availability_rows, fill = TRUE)
fwrite(availability_dt, file.path(results_dir, "phase3A_PPI_FS12_gene_availability.csv"))
append_log("[Phase3A-PPI-mini] Gene availability written: phase3A_PPI_FS12_gene_availability.csv")

algorithms <- data.table(
  algorithm = c("LASSO-Cox", "Elastic Net-Cox", "Ridge-Cox", "CoxBoost", "Random Survival Forest", "Survival-SVM", "GBM-Cox", "SuperPC"),
  package = c("glmnet", "glmnet", "glmnet", "CoxBoost", "randomForestSRC_or_ranger", "survivalsvm", "gbm", "superpc")
)
algorithms[, package_available := fifelse(
  package == "randomForestSRC_or_ranger",
  requireNamespace("randomForestSRC", quietly = TRUE) || requireNamespace("ranger", quietly = TRUE),
  vapply(package, requireNamespace, logical(1), quietly = TRUE)
)]
algorithms[, package_note := fifelse(
  algorithm == "Random Survival Forest" & requireNamespace("randomForestSRC", quietly = TRUE),
  "randomForestSRC used",
  fifelse(algorithm == "Random Survival Forest" & !requireNamespace("randomForestSRC", quietly = TRUE) & requireNamespace("ranger", quietly = TRUE),
          "ranger survival forest used as fallback",
          fifelse(package_available, "primary package available", paste0("package unavailable: ", package)))
)]

performance_rows <- list()
failure_rows <- list()
selected_rows <- list()
ev_rows <- list()

tcga_present <- intersect(fs12_genes, rownames(tcga_expr))
if (length(tcga_present) < 2) stop("FS12_PPI_top20 has fewer than two genes in TCGA; cannot run formal survival modeling.")
append_log("[Phase3A-PPI-mini] TCGA FS12-available genes: ", length(tcga_present), " / ", length(fs12_genes), " -> ", paste(tcga_present, collapse = ", "))

for (i in seq_len(nrow(algorithms))) {
  alg <- algorithms$algorithm[i]
  pkg_note <- algorithms$package_note[i]
  model_id <- sprintf("FS12_%02d", i)

  fail <- function(reason) {
    failure_rows[[length(failure_rows) + 1]] <<- data.table(
      planned_combination_id = model_id,
      algorithm = alg,
      feature_strategy = "FS12_PPI_top20_new",
      failure_reason = reason
    )
  }

  if (!isTRUE(algorithms$package_available[i])) {
    fail(paste0("Required package unavailable: ", algorithms$package[i]))
    next
  }

  model_obj <- tryCatch(fit_model_by_algorithm(alg, tcga_clin, tcga_x, tcga_present), error = function(e) e)
  if (inherits(model_obj, "error")) {
    fail(conditionMessage(model_obj))
    next
  }

  selected_genes <- model_obj$genes
  score_list <- tryCatch({
    list(
      tcga = predict_model_score(model_obj, tcga_expr),
      g376 = predict_model_score(model_obj, geo37642$expr),
      g124 = predict_model_score(model_obj, geo12417$expr),
      comb = predict_model_score(model_obj, combined_expr)
    )
  }, error = function(e) e)
  if (inherits(score_list, "error")) {
    fail(paste0("Prediction failed: ", conditionMessage(score_list)))
    next
  }

  ev <- rbindlist(list(
    eval_score(tcga_clin, score_list$tcga, "TCGA"),
    eval_score(geo37642$clin, score_list$g376, "GSE37642"),
    eval_score(geo12417$clin, score_list$g124, "GSE12417"),
    eval_score(combined_clin, score_list$comb, "combined_GPL570")
  ), fill = TRUE)

  wide <- dcast(ev, . ~ dataset, value.var = c("C_index", "AUC_1year", "AUC_3year", "AUC_5year", "IBS", "HR", "lower95", "upper95", "P_value"), fill = NA_real_)
  ext_hr <- ev[dataset %in% c("GSE37642", "GSE12417", "combined_GPL570"), HR]
  hr_direction_consistent <- length(ext_hr) > 0 && all(!is.na(ext_hr) & ext_hr > 1)
  external_p <- ev[dataset %in% c("GSE37642", "GSE12417"), P_value]
  combined_p <- ev[dataset == "combined_GPL570", P_value]
  external_any_significant <- any(external_p < 0.05, na.rm = TRUE)
  combined_significant <- length(combined_p) > 0 && isTRUE(combined_p[1] < 0.05)

  performance_rows[[length(performance_rows) + 1]] <- cbind(
    data.table(
      planned_combination_id = model_id,
      algorithm = alg,
      package_note = pkg_note,
      feature_strategy = "FS12_PPI_top20_new",
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
      TCGA_HR = HR_TCGA,
      TCGA_lower95 = lower95_TCGA,
      TCGA_upper95 = upper95_TCGA,
      TCGA_P_value = P_value_TCGA,
      GSE37642_C_index = C_index_GSE37642,
      GSE37642_AUC_1year = AUC_1year_GSE37642,
      GSE37642_AUC_3year = AUC_3year_GSE37642,
      GSE37642_AUC_5year = AUC_5year_GSE37642,
      GSE37642_IBS = IBS_GSE37642,
      GSE37642_HR = HR_GSE37642,
      GSE37642_lower95 = lower95_GSE37642,
      GSE37642_upper95 = upper95_GSE37642,
      GSE37642_P_value = P_value_GSE37642,
      GSE12417_C_index = C_index_GSE12417,
      GSE12417_AUC_1year = AUC_1year_GSE12417,
      GSE12417_AUC_3year = AUC_3year_GSE12417,
      GSE12417_AUC_5year = AUC_5year_GSE12417,
      GSE12417_IBS = IBS_GSE12417,
      GSE12417_HR = HR_GSE12417,
      GSE12417_lower95 = lower95_GSE12417,
      GSE12417_upper95 = upper95_GSE12417,
      GSE12417_P_value = P_value_GSE12417,
      combined_GPL570_C_index = C_index_combined_GPL570,
      combined_GPL570_AUC_1year = AUC_1year_combined_GPL570,
      combined_GPL570_AUC_3year = AUC_3year_combined_GPL570,
      combined_GPL570_AUC_5year = AUC_5year_combined_GPL570,
      combined_GPL570_IBS = IBS_combined_GPL570,
      combined_GPL570_HR = HR_combined_GPL570,
      combined_GPL570_lower95 = lower95_combined_GPL570,
      combined_GPL570_upper95 = upper95_combined_GPL570,
      combined_GPL570_P_value = P_value_combined_GPL570
    )],
    data.table(
      all_external_HR_direction_consistent = hr_direction_consistent,
      external_at_least_one_significant = external_any_significant,
      combined_GPL570_significant = combined_significant
    )
  )

  ev_out <- copy(ev)
  ev_out[, planned_combination_id := model_id]
  ev_out[, algorithm := alg]
  ev_out[, feature_strategy := "FS12_PPI_top20_new"]
  ev_rows[[length(ev_rows) + 1]] <- ev_out

  selected_rows[[length(selected_rows) + 1]] <- data.table(
    planned_combination_id = model_id,
    algorithm = alg,
    feature_strategy = "FS12_PPI_top20_new",
    gene_symbol = selected_genes,
    coefficient = if (model_obj$type == "linear") as.numeric(model_obj$coef[selected_genes]) else NA_real_,
    model_type = model_obj$type
  )
}

perf_list <- performance_rows
perf_dt <- if (length(perf_list) > 0) rbindlist(perf_list, fill = TRUE) else data.table()
failure_dt <- if (length(failure_rows) > 0) rbindlist(failure_rows, fill = TRUE) else data.table(planned_combination_id = character(), algorithm = character(), feature_strategy = character(), failure_reason = character())
ev_dt <- if (length(ev_rows) > 0) rbindlist(ev_rows, fill = TRUE) else data.table(dataset = character(), n = integer(), events = integer(), C_index = numeric(), AUC_1year = numeric(), AUC_3year = numeric(), AUC_5year = numeric(), IBS = numeric(), HR = numeric(), lower95 = numeric(), upper95 = numeric(), P_value = numeric(), planned_combination_id = character(), algorithm = character(), feature_strategy = character())
selected_dt <- if (length(selected_rows) > 0) rbindlist(selected_rows, fill = TRUE) else data.table()

fwrite(perf_dt, file.path(results_dir, "phase3A_PPI_mini_model_performance.csv"))
fwrite(ev_dt, file.path(results_dir, "phase3A_PPI_mini_external_validation.csv"))
fwrite(failure_dt, file.path(results_dir, "phase3A_PPI_mini_algorithm_failure_log.csv"))
append_log("[Phase3A-PPI-mini] FS12 mini-model performance written.")
append_log("[Phase3A-PPI-mini] Successful FS12 algorithms: ", ifelse(nrow(perf_dt) == 0, "none", paste(perf_dt$algorithm, collapse = "; ")))
append_log("[Phase3A-PPI-mini] Failed FS12 algorithms: ", ifelse(nrow(failure_dt) == 0, "none", paste(paste0(failure_dt$algorithm, " [", failure_dt$failure_reason, "]"), collapse = "; ")))

original6 <- c("CLCN5", "ARHGEF5", "ITGB2", "TRIM32", "SAT1", "ACOX2")
if (nrow(perf_dt) > 0) {
  max_genes_all <- max(c(phase3A_ranking$gene_count, perf_dt$gene_count), na.rm = TRUE)
  min_genes_all <- min(c(phase3A_ranking$gene_count, perf_dt$gene_count), na.rm = TRUE)
  perf_dt[, mean_external_C_index := rowMeans(.SD, na.rm = TRUE), .SDcols = c("GSE37642_C_index", "GSE12417_C_index", "combined_GPL570_C_index")]
  perf_dt[, mean_external_AUC_3year := rowMeans(.SD, na.rm = TRUE), .SDcols = c("GSE37642_AUC_3year", "GSE12417_AUC_3year", "combined_GPL570_AUC_3year")]
  perf_dt[, combined_score_component := fifelse(!is.na(combined_GPL570_HR) & combined_GPL570_HR > 1 & !is.na(combined_GPL570_P_value),
                                                pmin(1, -log10(pmax(combined_GPL570_P_value, .Machine$double.xmin)) / 3), 0)]
  perf_dt[, direction_score := fifelse(all_external_HR_direction_consistent == TRUE, 1, 0)]
  perf_dt[, gene_count_simplicity := 1 - (gene_count - min_genes_all) / max(1, max_genes_all - min_genes_all)]
  original_overlap_counts <- vapply(strsplit(perf_dt$actual_genes, ";", fixed = TRUE), function(g) length(intersect(g, original6)), integer(1))
  perf_dt[, biology_interpretability := pmin(1, original_overlap_counts / 6)]
  perf_dt[, external_overfit_gap := TCGA_C_index - mean_external_C_index]
  perf_dt[, anti_overfitting_score := pmax(0, pmin(1, 1 - pmax(0, external_overfit_gap) / 0.20))]
  component_matrix <- as.matrix(perf_dt[, .(
    mean_external_C_index,
    mean_external_AUC_3year,
    combined_score_component,
    direction_score,
    gene_count_simplicity,
    biology_interpretability,
    anti_overfitting_score
  )])
  component_weights <- c(0.35, 0.25, 0.15, 0.10, 0.05, 0.05, 0.05)
  available_matrix <- is.finite(component_matrix)
  weighted_values <- sweep(component_matrix, 2, component_weights, `*`)
  weighted_values[!available_matrix] <- 0
  available_weights <- sweep(available_matrix, 2, component_weights, `*`)
  denom <- rowSums(available_weights)
  score <- rowSums(weighted_values) / denom
  score[denom <= 0] <- NA_real_
  perf_dt[, composite_score := score]
  perf_dt[, composite_available_weight := denom]
  perf_dt[, overfitting_flag := is.finite(external_overfit_gap) & TCGA_C_index >= 0.70 & external_overfit_gap >= 0.12]
}

integrated_ranking <- rbindlist(list(phase3A_ranking, perf_dt), fill = TRUE)
setorder(integrated_ranking, -composite_score, gene_count)
integrated_ranking[, composite_rank := seq_len(.N)]
fwrite(integrated_ranking, file.path(results_dir, "phase3A_PPI_mini_integrated_model_ranking.csv"))
append_log("[Phase3A-PPI-mini] Integrated ranking written with FS12 models appended to prior Phase 3A-fix ranking.")

baseline_row <- integrated_ranking[planned_combination_id == "BASELINE_6GENE_FIXED"][1]
phase3a_top_row <- integrated_ranking[planned_combination_id == "M002"][1]
fs12_best_row <- if (nrow(perf_dt) > 0) integrated_ranking[feature_strategy == "FS12_PPI_top20_new"][order(composite_rank)][1] else data.table()

compare_dt <- rbindlist(list(
  if (nrow(baseline_row) > 0) baseline_row else NULL,
  if (nrow(phase3a_top_row) > 0) phase3a_top_row else NULL,
  if (nrow(fs12_best_row) > 0) fs12_best_row else NULL
), fill = TRUE)
if (ncol(compare_dt) > 0) {
  compare_dt[, model_label := fifelse(
    planned_combination_id == "BASELINE_6GENE_FIXED", "Original fixed 6-gene",
    fifelse(planned_combination_id == "M002", "Phase3A-fix top model", "Best FS12_PPI_top20 model")
  )]
  compare_dt[, comparison_note := fifelse(
    planned_combination_id == "BASELINE_6GENE_FIXED", "Retained baseline main model",
    fifelse(planned_combination_id == "M002", "Prior top-ranked model from Phase 3A-fix", "Additional network-derived feature strategy")
  )]
}
fwrite(compare_dt, file.path(results_dir, "phase3A_PPI_vs_original_6gene_comparison.csv"))

recommendation_text <- "FS12_PPI_top20 was evaluated as an additional network-derived feature strategy."
replacement_recommendation <- "No"
main_text_priority <- "Original fixed 6-gene signature remains the main model."
fs12_position <- "Supplement"
if (nrow(fs12_best_row) > 0 && nrow(baseline_row) > 0) {
  ext_delta <- fs12_best_row$mean_external_C_index - baseline_row$mean_external_C_index
  if (is.finite(ext_delta) && ext_delta > 0.05 && isTRUE(fs12_best_row$combined_GPL570_significant)) {
    recommendation_text <- "FS12_PPI_top20 provided additional prognostic information, but the original six-gene signature was retained for interpretability and cross-platform stability."
    fs12_position <- "Main-text optional comparison plus supplement"
  } else if (is.finite(ext_delta) && abs(ext_delta) <= 0.03) {
    recommendation_text <- "FS12_PPI_top20 showed comparable but not replacement-level performance."
    fs12_position <- "Supplement with brief main-text mention"
  } else {
    recommendation_text <- "Network-derived features provided complementary but not superior prognostic performance."
    fs12_position <- "Supplement"
  }
}

recommendation_dt <- data.table(
  topic = c(
    "FS12_main_text_suitability",
    "FS12_supplement_suitability",
    "replace_original_6gene",
    "retain_original_6gene_main_model",
    "methods_sentence_recommendation",
    "GSE6891_GSE14468_followup",
    "AS_input_audit_followup"
  ),
  recommendation = c(
    ifelse(fs12_position == "Main-text optional comparison plus supplement", "Conditional brief comparison", "No as standalone main-text replacement"),
    "Yes",
    replacement_recommendation,
    "Yes",
    "Use: network-derived feature strategy was additionally evaluated.",
    "Optional, only as extra validation after manual review of local data suitability.",
    "Optional after manual review; not triggered automatically by this mini update."
  ),
  rationale = c(
    recommendation_text,
    "FS12 is best framed as an added network-derived feature strategy rather than a model replacement.",
    "Do not automatically replace the original fixed six-gene signature.",
    "The original six-gene signature remains simpler and cross-platform stable.",
    "This wording stays within the requested writing boundary.",
    "Local presence was checked, but formal modeling was not forced in this update.",
    "The user explicitly requested stopping after this mini update."
  )
)
fwrite(recommendation_dt, file.path(results_dir, "phase3A_PPI_mini_main_vs_supplement_recommendation.csv"))

metric_summary <- function(row, prefix, include_ibs = FALSE) {
  if (nrow(row) == 0) return("NA")
  auc <- paste0(fmt_num(row[[paste0(prefix, "_AUC_1year")]]), "/", fmt_num(row[[paste0(prefix, "_AUC_3year")]]), "/", fmt_num(row[[paste0(prefix, "_AUC_5year")]]))
  base <- paste0("C-index=", fmt_num(row[[paste0(prefix, "_C_index")]]), ", AUC1/3/5=", auc)
  if (include_ibs) base <- paste0(base, ", IBS=", fmt_num(row[[paste0(prefix, "_IBS")]]))
  base
}
external_summary <- function(row, prefix) {
  if (nrow(row) == 0) return("NA")
  paste0(metric_summary(row, prefix, include_ibs = TRUE), ", HR=", fmt_num(row[[paste0(prefix, "_HR")]]), " (", fmt_num(row[[paste0(prefix, "_lower95")]]), "-", fmt_num(row[[paste0(prefix, "_upper95")]]), "), P=", fmt_p(row[[paste0(prefix, "_P_value")]]))
}

if (nrow(integrated_ranking) > 0) {
  heat_dt <- melt(
    integrated_ranking[feature_strategy == "FS12_PPI_top20_new"],
    id.vars = c("planned_combination_id", "algorithm", "feature_strategy", "composite_rank"),
    measure.vars = c("TCGA_C_index", "GSE37642_C_index", "GSE12417_C_index", "combined_GPL570_C_index"),
    variable.name = "dataset",
    value.name = "C_index"
  )
  if (nrow(heat_dt) > 0) {
    heat_dt[, dataset := factor(gsub("_C_index", "", dataset), levels = c("TCGA", "GSE37642", "GSE12417", "combined_GPL570"))]
    heat_dt[, model_label := factor(paste0(planned_combination_id, " | ", algorithm), levels = rev(unique(paste0(planned_combination_id, " | ", algorithm)[order(composite_rank)])))]
    p_heat <- ggplot(heat_dt, aes(x = dataset, y = model_label, fill = C_index)) +
      geom_tile(colour = "white", linewidth = 0.2) +
      scale_fill_gradient2(low = "#D9EAF7", mid = "#F7F7F7", high = "#A33A3A", midpoint = 0.6, na.value = "grey90") +
      labs(x = NULL, y = NULL, fill = "C-index", title = "FS12_PPI_top20 model performance across training and external cohorts") +
      theme(axis.text.y = element_text(size = 5.6), axis.text.x = element_text(angle = 30, hjust = 1))
    save_pdf_checked(file.path(fig_dir, "phase3A_PPI_mini_model_performance_heatmap.pdf"), p_heat, width = 6.8, height = 4.8)
  }

  if (nrow(compare_dt) > 0) {
    comp_long <- melt(
      compare_dt[, .(model_label, planned_combination_id, algorithm, feature_strategy, TCGA_C_index, GSE37642_C_index, GSE12417_C_index, combined_GPL570_C_index)],
      id.vars = c("model_label", "planned_combination_id", "algorithm", "feature_strategy"),
      variable.name = "dataset",
      value.name = "C_index"
    )
    comp_long[, dataset := factor(gsub("_C_index", "", dataset), levels = c("TCGA", "GSE37642", "GSE12417", "combined_GPL570"))]
    p_comp <- ggplot(comp_long, aes(x = dataset, y = C_index, group = model_label, colour = model_label)) +
      geom_line(linewidth = 0.5) +
      geom_point(size = 2) +
      labs(x = NULL, y = "C-index", colour = NULL, title = "Original six-gene model versus Phase3A-fix top model and best FS12 model") +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
    save_pdf_checked(file.path(fig_dir, "phase3A_PPI_vs_original_6gene_comparison.pdf"), p_comp, width = 7.2, height = 4.8)

    ext_bar <- comp_long[dataset != "TCGA"]
    p_ext <- ggplot(ext_bar, aes(x = dataset, y = C_index, fill = model_label)) +
      geom_col(position = position_dodge(width = 0.72), width = 0.66) +
      labs(x = NULL, y = "External C-index", fill = NULL, title = "External validation comparison of baseline, prior top model, and best FS12 model") +
      theme(axis.text.x = element_text(angle = 25, hjust = 1))
    save_pdf_checked(file.path(fig_dir, "phase3A_PPI_external_validation_comparison.pdf"), p_ext, width = 7.0, height = 4.6)
  }

  if (nrow(fs12_best_row) > 0) {
    best_score <- ev_dt[planned_combination_id == fs12_best_row$planned_combination_id & dataset == "TCGA"]
    score_vec <- NULL
    if (nrow(best_score) > 0) {
      score_model <- perf_dt[planned_combination_id == fs12_best_row$planned_combination_id][1]
      score_vec <- rep(NA_real_, nrow(tcga_clin))
      rownames_placeholder <- NULL
      if (nrow(score_model) > 0 && score_model$algorithm %in% perf_dt$algorithm) {
        model_obj <- tryCatch(fit_model_by_algorithm(score_model$algorithm, tcga_clin, tcga_x, unlist(strsplit(score_model$actual_genes, ";", fixed = TRUE))), error = function(e) NULL)
        if (!is.null(model_obj)) score_vec <- predict_model_score(model_obj, tcga_expr)
      }
    }
    if (!is.null(score_vec) && any(is.finite(score_vec))) {
      km_dt <- data.table(OS_time = tcga_clin$OS_time, OS_status = tcga_clin$OS_status, score = as.numeric(score_vec))
      km_dt <- km_dt[is.finite(score) & OS_time > 0 & !is.na(OS_status)]
      if (nrow(km_dt) >= 20) {
        km_dt[, risk_group := ifelse(score >= median(score, na.rm = TRUE), "High", "Low")]
        fit <- survfit(Surv(OS_time, OS_status) ~ risk_group, data = km_dt)
        if (requireNamespace("survminer", quietly = TRUE)) {
          p_km <- survminer::ggsurvplot(
            fit,
            data = km_dt,
            risk.table = FALSE,
            pval = TRUE,
            palette = c("Low" = "#4C78A8", "High" = "#A33A3A"),
            legend.title = NULL,
            legend.labs = c("High", "Low")
          )$plot + labs(title = "Best FS12_PPI_top20 model in TCGA-LAML")
          save_pdf_checked(file.path(fig_dir, "phase3A_PPI_best_model_KM.pdf"), p_km, width = 5.2, height = 4.6)
        } else {
          append_log("[Phase3A-PPI-mini] survminer unavailable, so KM figure was not produced.")
        }
      }
    }
  }
}

checklist <- c(
  paste0("1. FS12_PPI_top20 read successfully: ", ifelse(length(fs12_genes) > 0, "yes", "no")),
  paste0("2. FS12_PPI_top20 gene list: ", paste(fs12_genes, collapse = ", ")),
  paste0("3. TCGA FS12 available gene count: ", availability_dt[dataset_id == "TCGA-LAML", available_genes]),
  paste0("4. GSE37642 FS12 available gene count: ", availability_dt[dataset_id == "GSE37642_GPL570", available_genes]),
  paste0("5. GSE12417 FS12 available gene count: ", availability_dt[dataset_id == "GSE12417_GPL570", available_genes]),
  paste0("6. combined GPL570 FS12 available gene count: ", availability_dt[dataset_id == "combined_GPL570", available_genes]),
  paste0("7. GSE6891/GSE14468 local data present: ", paste(availability_dt[dataset_id %in% c("GSE6891", "GSE14468"), paste0(dataset_id, "=", ifelse(grepl('not found', notes), 'no', 'yes'))], collapse = "; ")),
  paste0("8. Number of successfully rerun algorithms: ", nrow(perf_dt)),
  paste0("9. Failed algorithms and reasons: ", ifelse(nrow(failure_dt) == 0, "none", paste(paste0(failure_dt$algorithm, ": ", failure_dt$failure_reason), collapse = " || "))),
  paste0("10. Best FS12_PPI model: ", if (nrow(fs12_best_row) > 0) paste0(fs12_best_row$planned_combination_id, " / ", fs12_best_row$algorithm) else "NA"),
  paste0("11. Best FS12_PPI model TCGA C-index: ", metric_summary(fs12_best_row, "TCGA", include_ibs = TRUE)),
  paste0("12. Best FS12_PPI model TCGA 1/3/5-year AUC: ", if (nrow(fs12_best_row) > 0) paste0(fmt_num(fs12_best_row$TCGA_AUC_1year), "/", fmt_num(fs12_best_row$TCGA_AUC_3year), "/", fmt_num(fs12_best_row$TCGA_AUC_5year)) else "NA"),
  paste0("13. Best FS12_PPI external validation performance: ", paste(c(
    paste0("GSE37642[", external_summary(fs12_best_row, "GSE37642"), "]"),
    paste0("GSE12417[", external_summary(fs12_best_row, "GSE12417"), "]"),
    paste0("combined[", external_summary(fs12_best_row, "combined_GPL570"), "]")
  ), collapse = " ; ")),
  paste0("14. Best FS12_PPI model HR/P: ", external_summary(fs12_best_row, "combined_GPL570")),
  paste0("15. Original fixed 6-gene model performance: ", paste(c(
    metric_summary(baseline_row, "TCGA", include_ibs = TRUE),
    external_summary(baseline_row, "GSE37642"),
    external_summary(baseline_row, "GSE12417"),
    external_summary(baseline_row, "combined_GPL570")
  ), collapse = " ; ")),
  paste0("16. Phase3A-fix top model performance: ", paste(c(
    metric_summary(phase3a_top_row, "TCGA", include_ibs = TRUE),
    external_summary(phase3a_top_row, "GSE37642"),
    external_summary(phase3a_top_row, "GSE12417"),
    external_summary(phase3a_top_row, "combined_GPL570")
  ), collapse = " ; ")),
  paste0("17. Is FS12_PPI superior to the original 6-gene model: ", ifelse(nrow(fs12_best_row) == 0, "no evaluable model", ifelse(fs12_best_row$mean_external_C_index > baseline_row$mean_external_C_index + 0.05, "not automatically acted on; only additional information", "no or comparable only"))),
  paste0("18. Recommend replacing the original 6-gene main model: ", replacement_recommendation),
  paste0("19. Recommended main-text figure(s): ", "Keep original six-gene main-model panels; optional brief comparison can use phase3A_PPI_vs_original_6gene_comparison.pdf only if space allows."),
  paste0("20. Recommended supplementary figure(s): ", paste(c("phase3A_PPI_mini_model_performance_heatmap.pdf", "phase3A_PPI_vs_original_6gene_comparison.pdf", "phase3A_PPI_external_validation_comparison.pdf", if (file.exists(file.path(fig_dir, "phase3A_PPI_best_model_KM.pdf"))) "phase3A_PPI_best_model_KM.pdf" else NULL), collapse = "; ")),
  paste0("21. Recommend entering GSE6891/GSE14468 supplemental validation: ", "optional manual follow-up only"),
  paste0("22. Recommend entering AS input audit: ", "optional manual follow-up only"),
  paste0("23. Issues requiring manual confirmation: ", "FS12 genes are network-derived and should be described as complementary; original six-gene signature should remain retained unless the user explicitly decides otherwise after review.")
)
writeLines(checklist, file.path(log_dir, "phase3A_PPI_mini_key_result_checklist.txt"), useBytes = TRUE)

append_log("[Phase3A-PPI-mini] Checklist written: phase3A_PPI_mini_key_result_checklist.txt")
append_log("[Phase3A-PPI-mini] Recommendation summary: ", recommendation_text)
append_log("[Phase3A-PPI-mini] Finished at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
