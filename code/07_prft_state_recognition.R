#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(1234)

ascii_lib <- Sys.getenv("PHASE1_ASCII_R_LIB", unset = "")
if (nzchar(ascii_lib) && dir.exists(ascii_lib)) .libPaths(c(ascii_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pROC)
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

log_file <- file.path(log_dir, "phase3B_PRFT_state_ML_log.txt")
if (file.exists(log_file)) file.remove(log_file)

append_log <- function(...) {
  line <- paste0(...)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

fmt_num <- function(x, digits = 3) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep("NA", length(x))
  ok <- is.finite(x) & !is.na(x)
  out[ok] <- formatC(x[ok], format = "f", digits = digits)
  out
}

fmt_p <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep("NA", length(x))
  ok <- is.finite(x) & !is.na(x)
  out[ok] <- format(x[ok], digits = 3, scientific = TRUE)
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

zscore_vector <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) rep(0, length(x)) else (x - m) / s
}

safe_auc <- function(y, prob) {
  ok <- !is.na(y) & is.finite(prob)
  y <- as.integer(y[ok])
  prob <- as.numeric(prob[ok])
  if (length(unique(y)) < 2 || length(prob) < 5) return(NA_real_)
  as.numeric(pROC::auc(response = y, predictor = prob, levels = c(0, 1), direction = "<", quiet = TRUE))
}

average_precision <- function(y, prob) {
  ok <- !is.na(y) & is.finite(prob)
  y <- as.integer(y[ok])
  prob <- as.numeric(prob[ok])
  if (length(unique(y)) < 2 || sum(y == 1) == 0) return(NA_real_)
  ord <- order(prob, decreasing = TRUE)
  y <- y[ord]
  tp <- cumsum(y == 1)
  precision <- tp / seq_along(y)
  sum(precision[y == 1]) / sum(y == 1)
}

calc_metrics <- function(y, prob, threshold = 0.5) {
  ok <- !is.na(y) & is.finite(prob)
  y <- as.integer(y[ok])
  prob <- as.numeric(prob[ok])
  if (length(unique(y)) < 2 || length(y) < 5) {
    return(data.table(AUROC = NA_real_, AUPRC = NA_real_, accuracy = NA_real_, balanced_accuracy = NA_real_,
                      F1 = NA_real_, MCC = NA_real_, sensitivity = NA_real_, specificity = NA_real_))
  }
  pred <- as.integer(prob >= threshold)
  tp <- sum(pred == 1 & y == 1)
  tn <- sum(pred == 0 & y == 0)
  fp <- sum(pred == 1 & y == 0)
  fn <- sum(pred == 0 & y == 1)
  sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  f1 <- if (is.finite(precision + sens) && (precision + sens) > 0) 2 * precision * sens / (precision + sens) else NA_real_
  denom <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- if (denom > 0) (tp * tn - fp * fn) / denom else NA_real_
  data.table(
    AUROC = safe_auc(y, prob),
    AUPRC = average_precision(y, prob),
    accuracy = mean(pred == y),
    balanced_accuracy = mean(c(sens, spec), na.rm = TRUE),
    F1 = f1,
    MCC = mcc,
    sensitivity = sens,
    specificity = spec
  )
}

scale_fit <- function(x) {
  center <- colMeans(x, na.rm = TRUE)
  scale <- apply(x, 2, stats::sd, na.rm = TRUE)
  scale[!is.finite(scale) | scale == 0] <- 1
  list(center = center, scale = scale)
}

scale_apply <- function(x, sf) {
  x <- sweep(x, 2, sf$center[colnames(x)], "-")
  sweep(x, 2, sf$scale[colnames(x)], "/")
}

make_cv_splits <- function(y, k = 5, repeats = 10) {
  y <- as.integer(y)
  splits <- list()
  for (r in seq_len(repeats)) {
    fold_id <- rep(NA_integer_, length(y))
    for (cls in sort(unique(y))) {
      idx <- sample(which(y == cls))
      fold_id[idx] <- rep(seq_len(k), length.out = length(idx))
    }
    for (fold in seq_len(k)) {
      test <- which(fold_id == fold)
      train <- setdiff(seq_along(y), test)
      splits[[length(splits) + 1L]] <- list(repeat_id = r, fold = fold, train = train, test = test)
    }
  }
  splits
}

train_predict <- function(algorithm, x_train, y_train, x_test) {
  y_train <- as.integer(y_train)
  if (algorithm == "Elastic Net Logistic Regression") {
    if (!requireNamespace("glmnet", quietly = TRUE)) stop("Required package unavailable: glmnet")
    fit <- glmnet::cv.glmnet(as.matrix(x_train), y_train, family = "binomial", alpha = 0.5, type.measure = "auc", nfolds = 5)
    return(as.numeric(predict(fit, newx = as.matrix(x_test), s = "lambda.1se", type = "response")))
  }
  if (algorithm == "Random Forest") {
    if (!requireNamespace("ranger", quietly = TRUE)) stop("Required package unavailable: ranger")
    d <- data.frame(y = factor(ifelse(y_train == 1, "high", "low"), levels = c("low", "high")), x_train, check.names = FALSE)
    fit <- ranger::ranger(y ~ ., data = d, probability = TRUE, num.trees = 500, importance = "permutation", seed = 1234)
    pr <- predict(fit, data = data.frame(x_test, check.names = FALSE))$predictions
    return(as.numeric(pr[, "high"]))
  }
  if (algorithm == "XGBoost") {
    if (!requireNamespace("xgboost", quietly = TRUE)) stop("Required package unavailable: xgboost")
  }
  if (algorithm == "CatBoost") {
    if (!requireNamespace("catboost", quietly = TRUE)) stop("Required package unavailable: catboost")
  }
  if (algorithm == "LightGBM") {
    if (!requireNamespace("lightgbm", quietly = TRUE)) stop("Required package unavailable: lightgbm")
  }
  if (algorithm == "SVM radial") {
    if (!requireNamespace("e1071", quietly = TRUE)) stop("Required package unavailable: e1071")
    fit <- e1071::svm(
      x = as.matrix(x_train),
      y = factor(ifelse(y_train == 1, "high", "low"), levels = c("low", "high")),
      kernel = "radial",
      probability = TRUE,
      scale = FALSE
    )
    pred <- predict(fit, as.matrix(x_test), probability = TRUE)
    probs <- attr(pred, "probabilities")
    return(as.numeric(probs[, "high"]))
  }
  if (algorithm == "GBM") {
    if (!requireNamespace("gbm", quietly = TRUE)) stop("Required package unavailable: gbm")
    d <- data.frame(y = y_train, x_train, check.names = FALSE)
    fit <- gbm::gbm(y ~ ., data = d, distribution = "bernoulli", n.trees = 300, interaction.depth = 2,
                    shrinkage = 0.03, bag.fraction = 0.8, train.fraction = 1, verbose = FALSE)
    return(as.numeric(predict(fit, newdata = data.frame(x_test, check.names = FALSE), n.trees = 300, type = "response")))
  }
  if (algorithm == "kNN") {
    if (!requireNamespace("class", quietly = TRUE)) stop("Required package unavailable: class")
    k <- max(3, round(sqrt(nrow(x_train))))
    pred <- class::knn(train = as.matrix(x_train), test = as.matrix(x_test),
                       cl = factor(ifelse(y_train == 1, "high", "low"), levels = c("low", "high")),
                       k = k, prob = TRUE)
    p_win <- attr(pred, "prob")
    return(ifelse(pred == "high", p_win, 1 - p_win))
  }
  if (algorithm == "Naive Bayes") {
    classes <- c(0, 1)
    prior <- vapply(classes, function(cls) mean(y_train == cls), numeric(1))
    means <- sapply(classes, function(cls) colMeans(x_train[y_train == cls, , drop = FALSE], na.rm = TRUE))
    sds <- sapply(classes, function(cls) apply(x_train[y_train == cls, , drop = FALSE], 2, stats::sd, na.rm = TRUE))
    sds[!is.finite(sds) | sds == 0] <- 1
    logprob <- sapply(seq_along(classes), function(i) {
      rowSums(stats::dnorm(as.matrix(x_test), mean = matrix(means[, i], nrow(x_test), ncol(x_test), byrow = TRUE),
                           sd = matrix(sds[, i], nrow(x_test), ncol(x_test), byrow = TRUE), log = TRUE)) + log(prior[i])
    })
    logprob <- logprob - apply(logprob, 1, max)
    prob <- exp(logprob)
    prob[, 2] / rowSums(prob)
  } else if (algorithm == "Logistic Regression baseline") {
    d <- data.frame(y = y_train, x_train, check.names = FALSE)
    fit <- stats::glm(y ~ ., data = d, family = stats::binomial())
    as.numeric(stats::predict(fit, newdata = data.frame(x_test, check.names = FALSE), type = "response"))
  } else {
    stop("Unknown algorithm: ", algorithm)
  }
}

fit_full_model <- function(algorithm, x_train, y_train) {
  y_train <- as.integer(y_train)
  if (algorithm == "Elastic Net Logistic Regression") {
    fit <- glmnet::cv.glmnet(as.matrix(x_train), y_train, family = "binomial", alpha = 0.5, type.measure = "auc", nfolds = 5)
    return(list(algorithm = algorithm, fit = fit))
  }
  if (algorithm == "Random Forest") {
    d <- data.frame(y = factor(ifelse(y_train == 1, "high", "low"), levels = c("low", "high")), x_train, check.names = FALSE)
    fit <- ranger::ranger(y ~ ., data = d, probability = TRUE, num.trees = 500, importance = "permutation", seed = 1234)
    return(list(algorithm = algorithm, fit = fit))
  }
  if (algorithm == "SVM radial") {
    fit <- e1071::svm(x = as.matrix(x_train), y = factor(ifelse(y_train == 1, "high", "low"), levels = c("low", "high")),
                      kernel = "radial", probability = TRUE, scale = FALSE)
    return(list(algorithm = algorithm, fit = fit))
  }
  if (algorithm == "GBM") {
    d <- data.frame(y = y_train, x_train, check.names = FALSE)
    fit <- gbm::gbm(y ~ ., data = d, distribution = "bernoulli", n.trees = 300, interaction.depth = 2,
                    shrinkage = 0.03, bag.fraction = 0.8, train.fraction = 1, verbose = FALSE)
    return(list(algorithm = algorithm, fit = fit))
  }
  if (algorithm == "kNN") {
    return(list(algorithm = algorithm, fit = list(x = as.matrix(x_train), y = y_train, k = max(3, round(sqrt(nrow(x_train)))))))
  }
  if (algorithm == "Naive Bayes") {
    classes <- c(0, 1)
    prior <- vapply(classes, function(cls) mean(y_train == cls), numeric(1))
    means <- sapply(classes, function(cls) colMeans(x_train[y_train == cls, , drop = FALSE], na.rm = TRUE))
    sds <- sapply(classes, function(cls) apply(x_train[y_train == cls, , drop = FALSE], 2, stats::sd, na.rm = TRUE))
    sds[!is.finite(sds) | sds == 0] <- 1
    return(list(algorithm = algorithm, fit = list(prior = prior, means = means, sds = sds)))
  }
  if (algorithm == "Logistic Regression baseline") {
    d <- data.frame(y = y_train, x_train, check.names = FALSE)
    fit <- stats::glm(y ~ ., data = d, family = stats::binomial())
    return(list(algorithm = algorithm, fit = fit))
  }
  stop("No full-model fit for unavailable algorithm: ", algorithm)
}

predict_full_model <- function(model, x_test) {
  algorithm <- model$algorithm
  if (algorithm == "Elastic Net Logistic Regression") {
    return(as.numeric(predict(model$fit, newx = as.matrix(x_test), s = "lambda.1se", type = "response")))
  }
  if (algorithm == "Random Forest") {
    pr <- predict(model$fit, data = data.frame(x_test, check.names = FALSE))$predictions
    return(as.numeric(pr[, "high"]))
  }
  if (algorithm == "SVM radial") {
    pred <- predict(model$fit, as.matrix(x_test), probability = TRUE)
    return(as.numeric(attr(pred, "probabilities")[, "high"]))
  }
  if (algorithm == "GBM") {
    return(as.numeric(predict(model$fit, newdata = data.frame(x_test, check.names = FALSE), n.trees = 300, type = "response")))
  }
  if (algorithm == "kNN") {
    pred <- class::knn(train = model$fit$x, test = as.matrix(x_test),
                       cl = factor(ifelse(model$fit$y == 1, "high", "low"), levels = c("low", "high")),
                       k = model$fit$k, prob = TRUE)
    p_win <- attr(pred, "prob")
    return(ifelse(pred == "high", p_win, 1 - p_win))
  }
  if (algorithm == "Naive Bayes") {
    classes <- c(0, 1)
    logprob <- sapply(seq_along(classes), function(i) {
      rowSums(stats::dnorm(as.matrix(x_test), mean = matrix(model$fit$means[, i], nrow(x_test), ncol(x_test), byrow = TRUE),
                           sd = matrix(model$fit$sds[, i], nrow(x_test), ncol(x_test), byrow = TRUE), log = TRUE)) + log(model$fit$prior[i])
    })
    logprob <- logprob - apply(logprob, 1, max)
    prob <- exp(logprob)
    return(prob[, 2] / rowSums(prob))
  }
  if (algorithm == "Logistic Regression baseline") {
    return(as.numeric(stats::predict(model$fit, newdata = data.frame(x_test, check.names = FALSE), type = "response")))
  }
  rep(NA_real_, nrow(x_test))
}

permutation_importance <- function(model, x, y, nrep = 5) {
  base_prob <- predict_full_model(model, x)
  base_auc <- safe_auc(y, base_prob)
  rbindlist(lapply(colnames(x), function(g) {
    drops <- replicate(nrep, {
      xp <- x
      xp[, g] <- sample(xp[, g])
      p <- predict_full_model(model, xp)
      base_auc - safe_auc(y, p)
    })
    data.table(gene_safe = g, importance = mean(drops, na.rm = TRUE), importance_sd = stats::sd(drops, na.rm = TRUE))
  }))
}

model_specific_importance <- function(model, x, y) {
  alg <- model$algorithm
  if (alg == "Elastic Net Logistic Regression") {
    co <- as.matrix(stats::coef(model$fit, s = "lambda.1se"))
    dt <- data.table(gene_safe = rownames(co), importance = abs(as.numeric(co[, 1])), importance_type = "abs_coefficient")
    return(dt[gene_safe != "(Intercept)"])
  }
  if (alg == "Random Forest") {
    imp <- ranger::importance(model$fit)
    return(data.table(gene_safe = names(imp), importance = as.numeric(imp), importance_type = "ranger_permutation"))
  }
  if (alg == "GBM") {
    sm <- suppressWarnings(gbm::summary.gbm(model$fit, plotit = FALSE))
    return(data.table(gene_safe = as.character(sm$var), importance = as.numeric(sm$rel.inf), importance_type = "gbm_relative_influence"))
  }
  if (alg == "Logistic Regression baseline") {
    co <- stats::coef(model$fit)
    dt <- data.table(gene_safe = names(co), importance = abs(as.numeric(co)), importance_type = "abs_coefficient")
    return(dt[gene_safe != "(Intercept)"])
  }
  pi <- permutation_importance(model, x, y, nrep = 3)
  pi[, importance_type := "permutation_auc_drop"]
  pi
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
  expr_df <- read.delim(text = paste(table_lines, collapse = "\n"), sep = "\t", header = TRUE,
                        check.names = FALSE, quote = "\"", comment.char = "")
  expr_ids <- as.character(expr_df[[1]])
  expr_mat <- as.matrix(expr_df[, -1, drop = FALSE])
  storage.mode(expr_mat) <- "numeric"
  rownames(expr_mat) <- expr_ids
  list(expr = expr_mat)
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
  read.delim(text = paste(table_lines, collapse = "\n"), sep = "\t", header = TRUE,
             check.names = FALSE, quote = "\"", comment.char = "")
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
  gpl_map <- data.table(probe_id = as.character(gpl_df[[id_idx]]), raw_symbol = clean_symbol_string(gpl_df[[symbol_idx]]))
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

align_geo_expr <- function(dataset_id, matrix_file, gpl_df, gene_universe) {
  parsed <- parse_series_matrix(file.path(root_dir, "phase1_runtime", "00_raw_data", "geo_validation", matrix_file))
  expr_mat <- parsed$expr
  mapping_dt <- build_probe_mapping(expr_mat, gpl_df)
  sig_map <- mapping_dt[gene_symbol %in% gene_universe & probe_rank == 1]
  present_genes <- intersect(gene_universe, sig_map$gene_symbol)
  if (length(present_genes) < 5) stop(dataset_id, " has too few mapped genes.")
  probe_ids <- sig_map[match(present_genes, gene_symbol)]$probe_id
  expr_gene <- expr_mat[probe_ids, , drop = FALSE]
  rownames(expr_gene) <- present_genes
  expr_gene <- expr_gene[present_genes, , drop = FALSE]
  expr_gene_z <- t(apply(expr_gene, 1, zscore_vector))
  rownames(expr_gene_z) <- rownames(expr_gene)
  colnames(expr_gene_z) <- colnames(expr_gene)
  expr_gene_z
}

score_geo_prft <- function(expr_gene_z, gene_sets) {
  sample_ids <- colnames(expr_gene_z)
  prot <- intersect(gene_sets$Proteostasis_core, rownames(expr_gene_z))
  ferro <- intersect(gene_sets$Ferroptosis_tolerance_set, rownames(expr_gene_z))
  if (length(prot) < 5 || length(ferro) < 5) stop("Insufficient gene-set coverage for GEO PRFT state scoring.")
  prot_score <- colMeans(expr_gene_z[prot, , drop = FALSE], na.rm = TRUE)
  ferro_score <- colMeans(expr_gene_z[ferro, , drop = FALSE], na.rm = TRUE)
  prft <- (zscore_vector(prot_score) + zscore_vector(ferro_score)) / 2
  data.table(sample_id = sample_ids, PRFT_score = prft,
             y = as.integer(prft >= stats::median(prft, na.rm = TRUE)))
}

make_safe_matrix <- function(expr_mat, samples, genes) {
  genes <- intersect(genes, rownames(expr_mat))
  x <- t(expr_mat[genes, samples, drop = FALSE])
  safe <- make.names(genes, unique = TRUE)
  colnames(x) <- safe
  list(x = as.matrix(x), genes = genes, safe = safe)
}

append_log("[Phase3B] Started at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
append_log("[Phase3B] set.seed(1234) fixed.")
append_log("[Phase3B] Task scope: PRFT-high state recognition and interpretability only; survival outcome is not used as label.")
append_log("[Phase3B] Figure contract: quantitative grid; evidence chain includes model performance, ROC/PR curves, feature importance/SHAP fallback, and Phase3A-Phase3B gene overlap.")

target_packages <- c("glmnet", "ranger", "randomForest", "xgboost", "catboost", "lightgbm", "e1071", "gbm", "class",
                     "naivebayes", "klaR", "pROC", "PRROC", "caret", "vip", "iml", "fastshap", "SHAPforxgboost",
                     "data.table", "ggplot2", "pheatmap", "ggrepel")
pkg_dt <- rbindlist(lapply(target_packages, function(pkg) {
  available <- requireNamespace(pkg, quietly = TRUE)
  data.table(package = pkg, available = available,
             version = if (available) as.character(utils::packageVersion(pkg)) else NA_character_)
}))
append_log("[Phase3B] Package availability: ", paste(pkg_dt[, paste(package, available, version, sep = "=")], collapse = "; "))

expr_all <- readRDS(file.path(root_dir, "phase1_runtime", "02_processed_data", "tcga_expr_hgnc_log2cpm.rds"))
prft <- as.data.table(readRDS(file.path(root_dir, "phase1_runtime", "04_prft_score", "tcga_prft_score.rds")))
candidate33 <- fread(file.path(results_dir, "phase1_33_candidates.csv"))
candidate715 <- fread(file.path(results_dir, "phase1_715_candidates.csv"))
phase3a_freq <- fread(file.path(results_dir, "phase3A_fix_gene_recurrence_frequency.csv"))
gene_sets_all <- readRDS(file.path(root_dir, "phase1_runtime", "03_gene_sets", "prft_gene_sets_all.rds"))

prft <- prft[sample_id %in% colnames(expr_all) & !is.na(PRFT_score)]
setorder(prft, sample_id)
samples <- prft$sample_id
expr_all <- expr_all[, samples, drop = FALSE]
prft[, y := as.integer(PRFT_score >= stats::median(PRFT_score, na.rm = TRUE))]
prft[, PRFT_group_cv := ifelse(y == 1, "PRFT_high", "PRFT_low")]

original6 <- c("CLCN5", "ITGB2", "ARHGEF5", "TRIM32", "SAT1", "ACOX2")
phase3a_top20 <- head(phase3a_freq[order(-recurrence_count, gene_symbol), gene_symbol], 20)
phase3a_recurrent <- phase3a_freq[recurrence_count > 0, gene_symbol]
core_axis_genes <- unique(c(
  gene_sets_all$Proteostasis_core,
  gene_sets_all$Ferroptosis_tolerance_set,
  gene_sets_all$SLC7A11_GPX4_GSH_axis,
  gene_sets_all$JAK2_STAT5_PDL1_set
))
core_axis_genes <- setdiff(core_axis_genes, unique(c(gene_sets_all$SUMOylation_set, gene_sets_all$NEDDylation_set)))

feature_set_static <- list(
  "FS-A_original_6gene" = intersect(original6, rownames(expr_all)),
  "FS-B_33_cross_platform_candidates" = intersect(candidate33$gene_symbol, rownames(expr_all)),
  "FS-C_phase3A_recurrence_top20" = intersect(phase3a_top20, rownames(expr_all)),
  "FS-D_715_PRFT_correlation_top50" = intersect(head(candidate715[order(-abs(GS_PRFT)), gene_symbol], 50), rownames(expr_all)),
  "FS-E_33_intersect_phase3A_recurrent" = intersect(intersect(candidate33$gene_symbol, phase3a_recurrent), rownames(expr_all)),
  "FS-F_core_axis_explanatory_genes" = intersect(core_axis_genes, rownames(expr_all))
)

feature_set_dt <- rbindlist(lapply(names(feature_set_static), function(fs) {
  genes <- feature_set_static[[fs]]
  data.table(feature_set = fs, n_genes = length(genes), genes = paste(genes, collapse = ";"),
             notes = ifelse(fs == "FS-F_core_axis_explanatory_genes",
                            "Contains gene axes used to define PRFT biology; explanatory and label-proximal, not independent causal proof.",
                            ifelse(fs == "FS-D_715_PRFT_correlation_top50",
                                   "For cross-validation, top50 PRFT-correlation genes are reselected within each training fold.",
                                   "")))
}))
append_log("[Phase3B] Feature sets: ", paste(feature_set_dt[, paste0(feature_set, " n=", n_genes)], collapse = "; "))
append_log("[Phase3B] Repeated CV setting: 5-fold repeated 10 times due computational burden from multiple classifiers and feature sets.")
append_log("[Phase3B] Circularity guard: FS-D is reselected inside training folds; FS-F is treated as explanatory because it overlaps PRFT score-defining axes.")

select_features <- function(fs, train_idx = NULL) {
  if (fs != "FS-D_715_PRFT_correlation_top50" || is.null(train_idx)) return(feature_set_static[[fs]])
  train_samples <- samples[train_idx]
  genes <- intersect(candidate715$gene_symbol, rownames(expr_all))
  cors <- vapply(genes, function(g) {
    suppressWarnings(stats::cor(as.numeric(expr_all[g, train_samples]), prft$PRFT_score[train_idx], method = "spearman", use = "complete.obs"))
  }, numeric(1))
  head(genes[order(-abs(cors), genes)], 50)
}

algorithms <- c(
  "Elastic Net Logistic Regression",
  "Random Forest",
  "XGBoost",
  "CatBoost",
  "LightGBM",
  "SVM radial",
  "GBM",
  "kNN",
  "Naive Bayes",
  "Logistic Regression baseline"
)

failure_log <- data.table()
cv_predictions <- list()
performance_rows <- list()
splits <- make_cv_splits(prft$y, k = 5, repeats = 10)

for (fs in names(feature_set_static)) {
  for (alg in algorithms) {
    combo <- paste(alg, fs, sep = " | ")
    append_log("[Phase3B] Running ", combo)
    pred_rows <- list()
    combo_failed <- FALSE
    fail_reason <- NA_character_
    for (sp in splits) {
      genes <- select_features(fs, sp$train)
      genes <- intersect(genes, rownames(expr_all))
      if (length(genes) < 2) {
        combo_failed <- TRUE
        fail_reason <- "Feature strategy produced fewer than 2 available genes."
        break
      }
      mat <- make_safe_matrix(expr_all, samples, genes)
      x <- mat$x
      sf <- scale_fit(x[sp$train, , drop = FALSE])
      x_train <- scale_apply(x[sp$train, , drop = FALSE], sf)
      x_test <- scale_apply(x[sp$test, , drop = FALSE], sf)
      prob <- tryCatch(train_predict(alg, x_train, prft$y[sp$train], x_test), error = function(e) e)
      if (inherits(prob, "error")) {
        combo_failed <- TRUE
        fail_reason <- conditionMessage(prob)
        break
      }
      pred_rows[[length(pred_rows) + 1L]] <- data.table(
        algorithm = alg,
        feature_set = fs,
        repeat_id = sp$repeat_id,
        fold = sp$fold,
        sample_id = samples[sp$test],
        y = prft$y[sp$test],
        prob = pmin(pmax(as.numeric(prob), 0), 1),
        n_genes = length(genes),
        genes = paste(genes, collapse = ";")
      )
    }
    if (combo_failed) {
      failure_log <- rbind(failure_log, data.table(algorithm = alg, feature_set = fs, failure_reason = fail_reason), fill = TRUE)
      next
    }
    pred_dt <- rbindlist(pred_rows, fill = TRUE)
    cv_predictions[[combo]] <- pred_dt
    agg <- pred_dt[, .(prob = mean(prob, na.rm = TRUE), y = y[1]), by = sample_id]
    met <- calc_metrics(agg$y, agg$prob)
    performance_rows[[length(performance_rows) + 1L]] <- cbind(
      data.table(algorithm = alg, feature_set = fs, n_genes_median = stats::median(pred_dt$n_genes, na.rm = TRUE),
                 converged = TRUE, failure_reason = NA_character_),
      met
    )
  }
}

performance <- rbindlist(performance_rows, fill = TRUE)
if (nrow(failure_log) > 0) {
  failed_perf <- failure_log[, .(algorithm, feature_set, n_genes_median = NA_real_, converged = FALSE,
                                 failure_reason, AUROC = NA_real_, AUPRC = NA_real_, accuracy = NA_real_,
                                 balanced_accuracy = NA_real_, F1 = NA_real_, MCC = NA_real_,
                                 sensitivity = NA_real_, specificity = NA_real_)]
  performance_all <- rbindlist(list(performance, failed_perf), fill = TRUE)
} else {
  performance_all <- copy(performance)
}
setorder(performance, -AUROC, -AUPRC, -MCC)
setorder(performance_all, algorithm, feature_set)

fwrite(performance_all, file.path(results_dir, "phase3B_classification_model_performance.csv"))
fwrite(failure_log, file.path(results_dir, "phase3B_model_failure_log.csv"))

best <- if (nrow(performance) > 0) performance[1] else data.table()
append_log("[Phase3B] Successful model-feature combinations: ", nrow(performance))
append_log("[Phase3B] Failed model-feature combinations: ", nrow(failure_log))
append_log("[Phase3B] Best CV model: ", if (nrow(best) > 0) paste(best$algorithm, best$feature_set, sep = " + ") else "NA")

safe_train_full <- function(alg, fs) {
  genes <- select_features(fs, NULL)
  genes <- intersect(genes, rownames(expr_all))
  mat <- make_safe_matrix(expr_all, samples, genes)
  sf <- scale_fit(mat$x)
  x_scaled <- scale_apply(mat$x, sf)
  model <- fit_full_model(alg, x_scaled, prft$y)
  list(model = model, x_scaled = x_scaled, genes = mat$genes, safe = mat$safe, scaler = sf)
}

importance_all <- list()
full_models <- list()
for (i in seq_len(nrow(performance))) {
  row <- performance[i]
  fm <- tryCatch(safe_train_full(row$algorithm, row$feature_set), error = function(e) e)
  if (inherits(fm, "error")) {
    failure_log <- rbind(failure_log, data.table(algorithm = row$algorithm, feature_set = row$feature_set,
                                                 failure_reason = paste("Full model importance failed:", conditionMessage(fm))), fill = TRUE)
    next
  }
  key <- paste(row$algorithm, row$feature_set, sep = " | ")
  full_models[[key]] <- fm
  gene_map <- data.table(gene_safe = fm$safe, gene_symbol = fm$genes)
  imp <- tryCatch(model_specific_importance(fm$model, fm$x_scaled, prft$y), error = function(e) e)
  if (inherits(imp, "error")) {
    pi <- tryCatch(permutation_importance(fm$model, fm$x_scaled, prft$y, nrep = 3), error = function(e) e)
    if (inherits(pi, "error")) next
    imp <- pi[, importance_type := "permutation_auc_drop"]
  }
  if (!"importance_sd" %in% names(imp)) imp[, importance_sd := NA_real_]
  if (!"importance_type" %in% names(imp)) imp[, importance_type := "model_specific_importance"]
  imp <- merge(imp, gene_map[, .(gene_safe, gene_symbol)], by = "gene_safe", all.x = TRUE)
  imp[, `:=`(algorithm = row$algorithm, feature_set = row$feature_set)]
  importance_all[[length(importance_all) + 1L]] <- imp[, .(algorithm, feature_set, gene_symbol, importance, importance_sd, importance_type)]
}
importance_dt <- rbindlist(importance_all, fill = TRUE)
importance_dt[!is.finite(importance), importance := 0]
setorder(importance_dt, algorithm, feature_set, -importance)
fwrite(importance_dt, file.path(results_dir, "phase3B_feature_importance_all_models.csv"))

tree_perf <- performance[algorithm %in% c("Random Forest", "GBM")]
shap_success <- FALSE
shap_note <- "Exact SHAP not run because xgboost/catboost/fastshap/SHAPforxgboost were unavailable; permutation importance from the best available tree model is reported as an interpretable fallback."
shap_source <- if (nrow(tree_perf) > 0) tree_perf[1] else best
shap_key <- if (nrow(shap_source) > 0) paste(shap_source$algorithm, shap_source$feature_set, sep = " | ") else NA_character_
shap_dt <- data.table()
if (!is.na(shap_key) && shap_key %in% names(full_models)) {
  fm <- full_models[[shap_key]]
  pi <- permutation_importance(fm$model, fm$x_scaled, prft$y, nrep = 10)
  gene_map <- data.table(gene_safe = fm$safe, gene_symbol = fm$genes)
  shap_dt <- merge(pi, gene_map, by = "gene_safe", all.x = TRUE)
  shap_dt[, `:=`(algorithm = shap_source$algorithm, feature_set = shap_source$feature_set,
                 rank = frank(-importance, ties.method = "first"),
                 explanation_method = "permutation_importance_fallback")]
  setorder(shap_dt, rank)
}
fwrite(shap_dt[, .(rank, gene_symbol, importance, importance_sd, algorithm, feature_set, explanation_method)],
       file.path(results_dir, "phase3B_SHAP_top_features.csv"))

top_imp <- if (nrow(shap_dt) > 0) head(shap_dt$gene_symbol, 20) else character(0)
phase3a_overlap <- merge(
  data.table(gene_symbol = unique(c(top_imp, importance_dt[, head(gene_symbol[order(-importance)], 50)], phase3a_freq$gene_symbol))),
  phase3a_freq,
  by = "gene_symbol",
  all.x = TRUE
)
phase3b_gene_summary <- importance_dt[, .(
  phase3B_model_count = uniqueN(paste(algorithm, feature_set)),
  phase3B_mean_importance = mean(importance, na.rm = TRUE),
  phase3B_max_importance = max(importance, na.rm = TRUE)
), by = gene_symbol]
phase3a_overlap <- merge(phase3a_overlap, phase3b_gene_summary, by = "gene_symbol", all.x = TRUE)
phase3a_overlap[, in_original_6gene := gene_symbol %in% original6]
setorder(phase3a_overlap, -phase3B_max_importance, -recurrence_count)
fwrite(phase3a_overlap, file.path(results_dir, "phase3B_overlap_with_phase3A_recurrent_genes.csv"))

original6_rank <- importance_dt[gene_symbol %in% original6, .(
  best_rank = min(frank(-importance, ties.method = "first"), na.rm = TRUE),
  mean_importance = mean(importance, na.rm = TRUE),
  max_importance = max(importance, na.rm = TRUE),
  model_count = uniqueN(paste(algorithm, feature_set))
), by = gene_symbol]
shap_rank <- if (nrow(shap_dt) > 0) shap_dt[gene_symbol %in% original6, .(gene_symbol, shap_or_fallback_rank = rank, shap_or_fallback_importance = importance)] else data.table(gene_symbol = original6)
original6_rank <- merge(data.table(gene_symbol = original6), original6_rank, by = "gene_symbol", all.x = TRUE)
original6_rank <- merge(original6_rank, shap_rank, by = "gene_symbol", all.x = TRUE)
setorder(original6_rank, shap_or_fallback_rank, -max_importance)
fwrite(original6_rank, file.path(results_dir, "phase3B_original_6gene_importance_rank.csv"))

external_rows <- list()
if (nrow(best) > 0) {
  external_note <- tryCatch({
    gpl570_df <- parse_gpl_annotation(file.path(root_dir, "phase1_runtime", "00_raw_data", "geo_validation", "GPL570_family.soft.gz"))
    geo_gene_universe <- unique(c(unlist(feature_set_static), gene_sets_all$Proteostasis_core, gene_sets_all$Ferroptosis_tolerance_set))
    geo37642 <- align_geo_expr("GSE37642", "GSE37642-GPL570_series_matrix.txt.gz", gpl570_df, geo_gene_universe)
    geo12417 <- align_geo_expr("GSE12417", "GSE12417-GPL570_series_matrix.txt.gz", gpl570_df, geo_gene_universe)
    common_geo_genes <- intersect(rownames(geo37642), rownames(geo12417))
    geo_list <- list(
      GSE37642 = geo37642,
      GSE12417 = geo12417,
      combined_GPL570 = cbind(geo37642[common_geo_genes, , drop = FALSE], geo12417[common_geo_genes, , drop = FALSE])
    )
    top_for_external <- head(performance, 5)
    for (i in seq_len(nrow(top_for_external))) {
      row <- top_for_external[i]
      for (dataset_name in names(geo_list)) {
        geo_expr <- geo_list[[dataset_name]]
        label <- score_geo_prft(geo_expr, gene_sets_all)
        genes <- intersect(select_features(row$feature_set, NULL), intersect(rownames(expr_all), rownames(geo_expr)))
        if (length(genes) < 2) {
          external_rows[[length(external_rows) + 1L]] <- data.table(
            dataset = dataset_name, algorithm = row$algorithm, feature_set = row$feature_set,
            n = nrow(label), n_genes = length(genes), AUROC = NA_real_, AUPRC = NA_real_,
            note = "Too few common genes for external PRFT-state validation."
          )
          next
        }
        train_mat <- make_safe_matrix(expr_all, samples, genes)
        sf <- scale_fit(train_mat$x)
        x_train <- scale_apply(train_mat$x, sf)
        model <- fit_full_model(row$algorithm, x_train, prft$y)
        geo_samples <- intersect(label$sample_id, colnames(geo_expr))
        geo_x0 <- t(geo_expr[genes, geo_samples, drop = FALSE])
        colnames(geo_x0) <- train_mat$safe
        geo_x <- scale_apply(as.matrix(geo_x0), sf)
        prob <- predict_full_model(model, geo_x)
        label_sub <- label[match(geo_samples, sample_id)]
        met <- calc_metrics(label_sub$y, prob)
        external_rows[[length(external_rows) + 1L]] <- cbind(
          data.table(dataset = dataset_name, algorithm = row$algorithm, feature_set = row$feature_set,
                     n = nrow(label_sub), n_genes = length(genes), note = "GEO label defined by within-cohort median approximate PRFT score using available gene-set expression."),
          met[, .(AUROC, AUPRC, accuracy, balanced_accuracy, F1, MCC)]
        )
      }
    }
    "External GEO PRFT-state validation completed using approximate within-cohort PRFT score labels."
  }, error = function(e) {
    append_log("[Phase3B] External GEO PRFT-state validation skipped/failed: ", conditionMessage(e))
    paste("External GEO PRFT-state validation failed:", conditionMessage(e))
  })
  append_log("[Phase3B] ", external_note)
}
external_dt <- rbindlist(external_rows, fill = TRUE)
if (nrow(external_dt) == 0) {
  external_dt <- data.table(dataset = character(), algorithm = character(), feature_set = character(), n = integer(),
                            n_genes = integer(), AUROC = numeric(), AUPRC = numeric(), note = character())
}
fwrite(external_dt, file.path(results_dir, "phase3B_external_state_validation_if_available.csv"))

theme_set(
  theme_classic(base_size = 8) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 7),
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", size = 9),
      panel.grid = element_blank()
    )
)

perf_plot_dt <- copy(performance)
perf_plot_dt[, label := sprintf("%.2f", AUROC)]
p_heat <- ggplot(perf_plot_dt, aes(x = feature_set, y = algorithm, fill = AUROC)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = label), size = 2.0) +
  scale_fill_gradient(low = "#edf2f7", high = "#c0392b", na.value = "grey90", limits = c(0.5, 1)) +
  labs(x = NULL, y = NULL, fill = "AUROC", title = "PRFT-high state classification performance") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 6), axis.text.y = element_text(size = 6))
ggsave(file.path(fig_dir, "phase3B_model_performance_heatmap.pdf"), p_heat, width = 9.5, height = 4.8, device = cairo_pdf, bg = "white")

curve_data <- function(y, prob, model_label) {
  th <- sort(unique(c(0, 1, prob)), decreasing = TRUE)
  rbindlist(lapply(th, function(t) {
    pred <- as.integer(prob >= t)
    tp <- sum(pred == 1 & y == 1)
    fp <- sum(pred == 1 & y == 0)
    tn <- sum(pred == 0 & y == 0)
    fn <- sum(pred == 0 & y == 1)
    data.table(
      model = model_label,
      threshold = t,
      TPR = if ((tp + fn) > 0) tp / (tp + fn) else NA_real_,
      FPR = if ((fp + tn) > 0) fp / (fp + tn) else NA_real_,
      recall = if ((tp + fn) > 0) tp / (tp + fn) else NA_real_,
      precision = if ((tp + fp) > 0) tp / (tp + fp) else 1
    )
  }))
}

top5 <- head(performance, 5)
curves <- list()
for (i in seq_len(nrow(top5))) {
  key <- paste(top5$algorithm[i], top5$feature_set[i], sep = " | ")
  pred_dt <- cv_predictions[[key]]
  agg <- pred_dt[, .(prob = mean(prob, na.rm = TRUE), y = y[1]), by = sample_id]
  curves[[length(curves) + 1L]] <- curve_data(agg$y, agg$prob, paste0(top5$algorithm[i], "\n", top5$feature_set[i]))
}
curve_dt <- rbindlist(curves, fill = TRUE)
p_roc <- ggplot(curve_dt, aes(x = FPR, y = TPR, colour = model)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey75", linewidth = 0.35) +
  geom_path(linewidth = 0.65, alpha = 0.95) +
  coord_equal() +
  labs(x = "False positive rate", y = "True positive rate", title = "Cross-validated AUROC curves for top PRFT-state models")
ggsave(file.path(fig_dir, "phase3B_AUROC_curves_top_models.pdf"), p_roc, width = 6.5, height = 5.2, device = cairo_pdf, bg = "white")

p_pr <- ggplot(curve_dt, aes(x = recall, y = precision, colour = model)) +
  geom_path(linewidth = 0.65, alpha = 0.95) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Recall", y = "Precision", title = "Cross-validated AUPRC curves for top PRFT-state models")
ggsave(file.path(fig_dir, "phase3B_AUPRC_curves_top_models.pdf"), p_pr, width = 6.5, height = 5.2, device = cairo_pdf, bg = "white")

shap_plot_dt <- head(shap_dt[order(-importance)], 30)
if (nrow(shap_plot_dt) == 0) shap_plot_dt <- data.table(gene_symbol = "NA", importance = 0)
p_shap <- ggplot(shap_plot_dt, aes(x = reorder(gene_symbol, importance), y = importance, fill = gene_symbol %in% original6)) +
  geom_col(width = 0.72) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#c0392b", "FALSE" = "#4c78a8")) +
  labs(x = NULL, y = "Permutation AUROC drop", title = "Interpretable contribution to PRFT-high state recognition")
ggsave(file.path(fig_dir, "phase3B_SHAP_summary_plot.pdf"), p_shap, width = 6.2, height = 5.6, device = cairo_pdf, bg = "white")

imp_bar <- importance_dt[, .(importance = max(importance, na.rm = TRUE)), by = gene_symbol][order(-importance)][1:min(.N, 30)]
p_imp <- ggplot(imp_bar, aes(x = reorder(gene_symbol, importance), y = importance, fill = gene_symbol %in% original6)) +
  geom_col(width = 0.72) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#c0392b", "FALSE" = "#5f6f7f")) +
  labs(x = NULL, y = "Maximum model importance", title = "Top PRFT-state recognition features across models")
ggsave(file.path(fig_dir, "phase3B_feature_importance_barplot.pdf"), p_imp, width = 6.2, height = 5.6, device = cairo_pdf, bg = "white")

overlap_plot <- phase3a_overlap[!is.na(recurrence_count) & !is.na(phase3B_max_importance)][1:min(.N, 40)]
p_overlap <- ggplot(overlap_plot, aes(x = recurrence_count, y = phase3B_max_importance, colour = in_original_6gene, label = gene_symbol)) +
  geom_point(size = 2.0, alpha = 0.85) +
  scale_colour_manual(values = c("TRUE" = "#c0392b", "FALSE" = "#4c78a8")) +
  labs(x = "Phase 3A recurrence count", y = "Phase 3B maximum importance", title = "Concordance between survival-ML recurrence and PRFT-state interpretability")
if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_overlap <- p_overlap +
    ggrepel::geom_text_repel(data = overlap_plot[in_original_6gene == TRUE | frank(-phase3B_max_importance) <= 8],
                             size = 2.2, max.overlaps = 30, show.legend = FALSE)
} else {
  p_overlap <- p_overlap +
    geom_text(data = overlap_plot[in_original_6gene == TRUE | frank(-phase3B_max_importance) <= 6],
              size = 2.0, hjust = -0.05, vjust = 0.4, show.legend = FALSE)
  append_log("[Phase3B] ggrepel unavailable; overlap plot used geom_text fallback.")
}
ggsave(file.path(fig_dir, "phase3B_overlap_phase3A_phase3B_genes.pdf"), p_overlap, width = 6.4, height = 4.8, device = cairo_pdf, bg = "white")

success_alg <- sort(unique(performance$algorithm))
failed_alg <- if (nrow(failure_log) > 0) {
  failure_log[, .(reason = paste(unique(failure_reason), collapse = " | ")), by = algorithm][order(algorithm), paste0(algorithm, ": ", reason)]
} else character(0)
external_done <- nrow(external_dt) > 0 && any(is.finite(external_dt$AUROC))
external_summary <- if (external_done) {
  paste(external_dt[is.finite(AUROC), paste0(dataset, " ", algorithm, "+", feature_set, " AUROC=", fmt_num(AUROC), ", AUPRC=", fmt_num(AUPRC))], collapse = " | ")
} else "not completed"
top10_shap <- if (nrow(shap_dt) > 0) head(shap_dt$gene_symbol, 10) else character(0)
original6_in_top10 <- intersect(original6, top10_shap)
orig6_forward <- sum(original6_rank$shap_or_fallback_rank <= 20, na.rm = TRUE) >= 4
phase3a_phase3b_consistent <- length(intersect(head(phase3a_freq$gene_symbol, 20), head(shap_dt$gene_symbol, 20))) >= 3

checklist <- c(
  paste0("1. TCGA samples used for PRFT-high/low classification: ", nrow(prft)),
  paste0("2. PRFT-high samples: ", sum(prft$y == 1)),
  paste0("3. PRFT-low samples: ", sum(prft$y == 0)),
  paste0("4. Successful classification algorithms: ", paste(success_alg, collapse = "; ")),
  paste0("5. Failed algorithms and reasons: ", ifelse(length(failed_alg) == 0, "none", paste(failed_alg, collapse = " || "))),
  paste0("6. Best model: ", if (nrow(best) > 0) best$algorithm else "NA"),
  paste0("7. Best model feature set: ", if (nrow(best) > 0) best$feature_set else "NA"),
  paste0("8. Best model TCGA CV AUROC: ", if (nrow(best) > 0) fmt_num(best$AUROC) else "NA"),
  paste0("9. Best model TCGA CV AUPRC: ", if (nrow(best) > 0) fmt_num(best$AUPRC) else "NA"),
  paste0("10. Best model F1: ", if (nrow(best) > 0) fmt_num(best$F1) else "NA"),
  paste0("11. Best model MCC: ", if (nrow(best) > 0) fmt_num(best$MCC) else "NA"),
  paste0("12. External GEO state-recognition validation completed: ", ifelse(external_done, "yes", "no")),
  paste0("13. External validation AUROC/AUPRC if available: ", external_summary),
  paste0("14. SHAP completed: ", ifelse(shap_success, "yes", "no; permutation-importance fallback used")),
  paste0("15. Top 10 SHAP/fallback genes: ", paste(top10_shap, collapse = "; ")),
  paste0("16. Original 6 genes in SHAP/fallback top 10: ", ifelse(length(original6_in_top10) == 0, "none", paste(original6_in_top10, collapse = "; "))),
  paste0("17. Original 6 genes overall importance ranks are forward: ", ifelse(orig6_forward, "yes", "partial/no")),
  paste0("18. Phase3A high-frequency genes and Phase3B important genes consistent: ", ifelse(phase3a_phase3b_consistent, "yes", "partial/no")),
  "19. Recommend placing Phase3B in main Figure 4 or Figure 5: yes, as an interpretability/state-recognition panel with circularity caveat.",
  "20. Recommend entering Phase 3C: yes, after manual review of SHAP fallback and external GEO label approximation.",
  paste0("21. Issues requiring manual confirmation: ", shap_note, " GEO PRFT labels use approximate within-cohort gene-set scoring; FS-F is label-proximal and should not be described as independent causal evidence.")
)
writeLines(checklist, file.path(log_dir, "phase3B_key_result_checklist.txt"), useBytes = TRUE)

append_log("[Phase3B] SHAP note: ", shap_note)
append_log("[Phase3B] Checklist written: phase3B_key_result_checklist.txt")
append_log("[Phase3B] Finished at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
