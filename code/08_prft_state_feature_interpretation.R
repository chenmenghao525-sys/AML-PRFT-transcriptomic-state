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

log_file <- file.path(log_dir, "phase3B_fix_log.txt")
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

calc_metrics <- function(y, prob) {
  ok <- !is.na(y) & is.finite(prob)
  y <- as.integer(y[ok])
  prob <- as.numeric(prob[ok])
  if (length(unique(y)) < 2 || length(y) < 5) {
    return(data.table(AUROC = NA_real_, AUPRC = NA_real_, accuracy = NA_real_,
                      balanced_accuracy = NA_real_, F1 = NA_real_, MCC = NA_real_,
                      sensitivity = NA_real_, specificity = NA_real_))
  }
  pred <- as.integer(prob >= 0.5)
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
  data.table(AUROC = safe_auc(y, prob), AUPRC = average_precision(y, prob),
             accuracy = mean(pred == y), balanced_accuracy = mean(c(sens, spec), na.rm = TRUE),
             F1 = f1, MCC = mcc, sensitivity = sens, specificity = spec)
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
      splits[[length(splits) + 1L]] <- list(repeat_id = r, fold = fold,
                                            train = setdiff(seq_along(y), test), test = test)
    }
  }
  splits
}

make_safe_matrix <- function(expr_mat, samples, genes) {
  genes <- intersect(genes, rownames(expr_mat))
  x <- t(expr_mat[genes, samples, drop = FALSE])
  safe <- make.names(genes, unique = TRUE)
  colnames(x) <- safe
  list(x = as.matrix(x), genes = genes, safe = safe)
}

train_model <- function(algorithm, x_train, y_train) {
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
  if (algorithm == "GBM") {
    d <- data.frame(y = y_train, x_train, check.names = FALSE)
    fit <- gbm::gbm(y ~ ., data = d, distribution = "bernoulli", n.trees = 300, interaction.depth = 2,
                    shrinkage = 0.03, bag.fraction = 0.8, train.fraction = 1, verbose = FALSE)
    return(list(algorithm = algorithm, fit = fit))
  }
  if (algorithm == "XGBoost") {
    dtrain <- xgboost::xgb.DMatrix(data = as.matrix(x_train), label = y_train)
    fit <- xgboost::xgb.train(
      params = list(objective = "binary:logistic", eval_metric = "auc", max_depth = 2,
                    eta = 0.05, subsample = 0.85, colsample_bytree = 0.85, nthread = 1),
      data = dtrain,
      nrounds = 120,
      verbose = 0
    )
    return(list(algorithm = algorithm, fit = fit))
  }
  stop("Unknown algorithm: ", algorithm)
}

predict_model <- function(model, x_test) {
  if (model$algorithm == "Elastic Net Logistic Regression") {
    return(as.numeric(predict(model$fit, newx = as.matrix(x_test), s = "lambda.1se", type = "response")))
  }
  if (model$algorithm == "Random Forest") {
    pr <- predict(model$fit, data = data.frame(x_test, check.names = FALSE))$predictions
    return(as.numeric(pr[, "high"]))
  }
  if (model$algorithm == "GBM") {
    return(as.numeric(predict(model$fit, newdata = data.frame(x_test, check.names = FALSE), n.trees = 300, type = "response")))
  }
  if (model$algorithm == "XGBoost") {
    return(as.numeric(predict(model$fit, newdata = as.matrix(x_test))))
  }
  rep(NA_real_, nrow(x_test))
}

permutation_importance <- function(model, x, y, nrep = 8) {
  base_prob <- predict_model(model, x)
  base_auc <- safe_auc(y, base_prob)
  rbindlist(lapply(colnames(x), function(g) {
    drops <- replicate(nrep, {
      xp <- x
      xp[, g] <- sample(xp[, g])
      base_auc - safe_auc(y, predict_model(model, xp))
    })
    data.table(gene_safe = g, importance = mean(drops, na.rm = TRUE), importance_sd = stats::sd(drops, na.rm = TRUE),
               explanation_method = "permutation_importance_fallback")
  }))
}

xgb_shap_importance <- function(model, x, genes) {
  contrib <- predict(model$fit, newdata = as.matrix(x), predcontrib = TRUE)
  contrib <- as.matrix(contrib)
  contrib <- contrib[, colnames(contrib) != "BIAS", drop = FALSE]
  dt <- data.table(gene_safe = colnames(contrib), mean_abs_shap = colMeans(abs(contrib), na.rm = TRUE))
  dt <- merge(dt, data.table(gene_safe = colnames(x), gene_symbol = genes), by = "gene_safe", all.x = TRUE)
  dt <- dt[!is.na(gene_symbol) & nzchar(gene_symbol)]
  setorder(dt, -mean_abs_shap)
  dt[, rank := seq_len(.N)]
  dt[, explanation_method := "xgboost_predcontrib_SHAP"]
  dt
}

parse_series_matrix <- function(localfile) {
  lines <- readLines(gzfile(localfile), warn = FALSE)
  begin_idx <- grep("^!series_matrix_table_begin", lines)
  end_idx <- grep("^!series_matrix_table_end", lines)
  table_lines <- if (length(end_idx) == 0 || end_idx[1] <= begin_idx[1]) lines[(begin_idx[1] + 1):length(lines)] else lines[(begin_idx[1] + 1):(end_idx[1] - 1)]
  expr_df <- read.delim(text = paste(table_lines, collapse = "\n"), sep = "\t", header = TRUE,
                        check.names = FALSE, quote = "\"", comment.char = "")
  expr_ids <- as.character(expr_df[[1]])
  expr_mat <- as.matrix(expr_df[, -1, drop = FALSE])
  storage.mode(expr_mat) <- "numeric"
  rownames(expr_mat) <- expr_ids
  expr_mat
}

parse_gpl_annotation <- function(localfile) {
  lines <- readLines(gzfile(localfile), warn = FALSE)
  begin_idx <- grep("^!platform_table_begin", lines)
  end_idx <- grep("^!platform_table_end", lines)
  table_lines <- if (length(end_idx) == 0 || end_idx[1] <= begin_idx[1]) lines[(begin_idx[1] + 1):length(lines)] else lines[(begin_idx[1] + 1):(end_idx[1] - 1)]
  read.delim(text = paste(table_lines, collapse = "\n"), sep = "\t", header = TRUE,
             check.names = FALSE, quote = "\"", comment.char = "")
}

build_probe_mapping <- function(expr_mat, gpl_df) {
  norm_cols <- gsub("[^a-z0-9]", "", tolower(colnames(gpl_df)))
  symbol_idx <- match(c("genesymbol", "symbol", "hgncsymbol", "gene_symbol"), norm_cols)
  symbol_idx <- symbol_idx[!is.na(symbol_idx)][1]
  id_idx <- match(c("id", "idref", "probeid"), norm_cols)
  id_idx <- id_idx[!is.na(id_idx)][1]
  if (is.na(symbol_idx) || is.na(id_idx)) stop("GPL annotation lacks recognizable probe ID or gene symbol column.")
  gpl_map <- data.table(probe_id = as.character(gpl_df[[id_idx]]), raw_symbol = clean_symbol_string(gpl_df[[symbol_idx]]))
  gpl_map <- gpl_map[nzchar(probe_id) & nzchar(raw_symbol)]
  probe_iqr <- apply(expr_mat, 1, IQR, na.rm = TRUE)
  expanded <- gpl_map[, .(gene_symbol = unlist(strsplit(raw_symbol, ";", fixed = TRUE))), by = probe_id]
  expanded[, gene_symbol := trimws(gene_symbol)]
  expanded <- unique(expanded[nzchar(gene_symbol) & !gene_symbol %in% c("NA", "---")])
  merged <- merge(expanded, data.table(probe_id = rownames(expr_mat), probe_iqr = probe_iqr), by = "probe_id")
  setorder(merged, gene_symbol, -probe_iqr)
  merged[, probe_rank := seq_len(.N), by = gene_symbol]
  merged
}

align_geo_expr <- function(matrix_file, gpl_df, gene_universe) {
  expr_mat <- parse_series_matrix(file.path(root_dir, "phase1_runtime", "00_raw_data", "geo_validation", matrix_file))
  mapping <- build_probe_mapping(expr_mat, gpl_df)
  sig_map <- mapping[gene_symbol %in% gene_universe & probe_rank == 1]
  genes <- intersect(gene_universe, sig_map$gene_symbol)
  probes <- sig_map[match(genes, gene_symbol)]$probe_id
  expr_gene <- expr_mat[probes, , drop = FALSE]
  rownames(expr_gene) <- genes
  expr_gene <- expr_gene[genes, , drop = FALSE]
  expr_gene_z <- t(apply(expr_gene, 1, zscore_vector))
  rownames(expr_gene_z) <- rownames(expr_gene)
  colnames(expr_gene_z) <- colnames(expr_gene)
  expr_gene_z
}

score_geo_prft <- function(expr_gene_z, gene_sets) {
  prot <- intersect(gene_sets$Proteostasis_core, rownames(expr_gene_z))
  ferro <- intersect(gene_sets$Ferroptosis_tolerance_set, rownames(expr_gene_z))
  if (length(prot) < 5 || length(ferro) < 5) stop("Insufficient gene-set coverage for GEO PRFT scoring.")
  prot_score <- colMeans(expr_gene_z[prot, , drop = FALSE], na.rm = TRUE)
  ferro_score <- colMeans(expr_gene_z[ferro, , drop = FALSE], na.rm = TRUE)
  prft_score <- (zscore_vector(prot_score) + zscore_vector(ferro_score)) / 2
  data.table(sample_id = colnames(expr_gene_z), PRFT_score = prft_score,
             y = as.integer(prft_score >= stats::median(prft_score, na.rm = TRUE)))
}

append_log("[Phase3B-fix] Started at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
append_log("[Phase3B-fix] Scope: PRFT-high state recognition and model interpretation; no survival outcome labels and no causal language.")
append_log("[Phase3B-fix] Figure contract: quantitative grid separating independent-like feature sets from label-proximal FS-F.")

target_packages <- c("xgboost", "fastshap", "DALEX", "iml", "vip", "pROC", "PRROC", "caret", "glmnet", "ranger", "gbm")
pkg_dt <- rbindlist(lapply(target_packages, function(pkg) {
  available <- requireNamespace(pkg, quietly = TRUE)
  data.table(package = pkg, available = available,
             version = if (available) as.character(utils::packageVersion(pkg)) else NA_character_,
             note = ifelse(available, "available after Phase3B-fix package audit/install attempt",
                           ifelse(pkg == "fastshap", "not available for current R 4.5 binary repository during install attempt",
                                  "not available after package audit/install attempt")))
}))
fwrite(pkg_dt, file.path(results_dir, "phase3B_fix_package_availability.csv"))
append_log("[Phase3B-fix] Package availability: ", paste(pkg_dt[, paste(package, available, version, sep = "=")], collapse = "; "))

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

original6 <- c("CLCN5", "ITGB2", "ARHGEF5", "TRIM32", "SAT1", "ACOX2")
phase3a_top20 <- head(phase3a_freq[order(-recurrence_count, gene_symbol), gene_symbol], 20)
phase3a_recurrent <- phase3a_freq[recurrence_count > 0, gene_symbol]
core_axis_genes <- setdiff(unique(c(gene_sets_all$Proteostasis_core, gene_sets_all$Ferroptosis_tolerance_set,
                                    gene_sets_all$SLC7A11_GPX4_GSH_axis, gene_sets_all$JAK2_STAT5_PDL1_set)),
                           unique(c(gene_sets_all$SUMOylation_set, gene_sets_all$NEDDylation_set)))
feature_sets <- list(
  "FS-A_original_6gene" = intersect(original6, rownames(expr_all)),
  "FS-B_33_cross_platform_candidates" = intersect(candidate33$gene_symbol, rownames(expr_all)),
  "FS-C_phase3A_recurrence_top20" = intersect(phase3a_top20, rownames(expr_all)),
  "FS-E_33_intersect_phase3A_recurrent" = intersect(intersect(candidate33$gene_symbol, phase3a_recurrent), rownames(expr_all)),
  "FS-F_core_axis_explanatory_genes" = intersect(core_axis_genes, rownames(expr_all))
)
feature_class <- data.table(
  feature_set = names(feature_sets),
  interpretation_class = c("independent-like", "independent-like", "independent-like", "independent-like", "label-proximal explanatory"),
  n_genes = vapply(feature_sets, length, integer(1))
)
append_log("[Phase3B-fix] Feature-set classes: ", paste(feature_class[, paste(feature_set, interpretation_class, n_genes, sep = "=")], collapse = "; "))
append_log("[Phase3B-fix] FS-D_715_PRFT_correlation_top50 is excluded from focused retraining because it is PRFT-correlation selected and treated as intermediate/cautious evidence.")

algorithms <- c("Elastic Net Logistic Regression", "Random Forest", "GBM", "XGBoost")
splits <- make_cv_splits(prft$y, k = 5, repeats = 10)
performance_rows <- list()
prediction_cache <- list()
failure_rows <- list()

for (fs in names(feature_sets)) {
  genes <- feature_sets[[fs]]
  mat <- make_safe_matrix(expr_all, samples, genes)
  for (alg in algorithms) {
    if (alg == "XGBoost" && !requireNamespace("xgboost", quietly = TRUE)) {
      failure_rows[[length(failure_rows) + 1L]] <- data.table(algorithm = alg, feature_set = fs, failure_reason = "Required package unavailable: xgboost")
      next
    }
    append_log("[Phase3B-fix] Running ", alg, " | ", fs)
    pred_rows <- list()
    failed <- FALSE
    fail_reason <- NA_character_
    for (sp in splits) {
      sf <- scale_fit(mat$x[sp$train, , drop = FALSE])
      x_train <- scale_apply(mat$x[sp$train, , drop = FALSE], sf)
      x_test <- scale_apply(mat$x[sp$test, , drop = FALSE], sf)
      model <- tryCatch(train_model(alg, x_train, prft$y[sp$train]), error = function(e) e)
      if (inherits(model, "error")) {
        failed <- TRUE
        fail_reason <- conditionMessage(model)
        break
      }
      prob <- tryCatch(predict_model(model, x_test), error = function(e) e)
      if (inherits(prob, "error")) {
        failed <- TRUE
        fail_reason <- conditionMessage(prob)
        break
      }
      pred_rows[[length(pred_rows) + 1L]] <- data.table(sample_id = samples[sp$test], y = prft$y[sp$test],
                                                        prob = pmin(pmax(prob, 0), 1),
                                                        repeat_id = sp$repeat_id, fold = sp$fold)
    }
    if (failed) {
      failure_rows[[length(failure_rows) + 1L]] <- data.table(algorithm = alg, feature_set = fs, failure_reason = fail_reason)
      next
    }
    pred_dt <- rbindlist(pred_rows)
    prediction_cache[[paste(alg, fs, sep = " | ")]] <- pred_dt
    agg <- pred_dt[, .(prob = mean(prob, na.rm = TRUE), y = y[1]), by = sample_id]
    met <- calc_metrics(agg$y, agg$prob)
    performance_rows[[length(performance_rows) + 1L]] <- cbind(
      data.table(algorithm = alg, feature_set = fs, interpretation_class = feature_class[feature_set == fs, interpretation_class],
                 n_genes = length(genes), converged = TRUE),
      met
    )
  }
}

perf <- rbindlist(performance_rows, fill = TRUE)
failures <- rbindlist(failure_rows, fill = TRUE)
setorder(perf, feature_set, -AUROC, -AUPRC)
fwrite(perf, file.path(results_dir, "phase3B_fix_focused_model_performance.csv"))
fwrite(failures, file.path(results_dir, "phase3B_fix_model_failure_log.csv"))

best_by_fs <- perf[, .SD[order(-AUROC, -AUPRC)][1], by = feature_set]

train_full <- function(alg, fs) {
  genes <- feature_sets[[fs]]
  mat <- make_safe_matrix(expr_all, samples, genes)
  sf <- scale_fit(mat$x)
  x_scaled <- scale_apply(mat$x, sf)
  model <- train_model(alg, x_scaled, prft$y)
  list(model = model, x_scaled = x_scaled, genes = mat$genes, safe = mat$safe, scaler = sf)
}

full_models <- list()
importance_rows <- list()
for (i in seq_len(nrow(best_by_fs))) {
  row <- best_by_fs[i]
  fm <- train_full(row$algorithm, row$feature_set)
  full_models[[row$feature_set]] <- fm
  if (row$algorithm == "XGBoost") {
    imp <- xgb_shap_importance(fm$model, fm$x_scaled, fm$genes)
  } else {
    pi <- permutation_importance(fm$model, fm$x_scaled, prft$y, nrep = 10)
    imp <- merge(pi, data.table(gene_safe = fm$safe, gene_symbol = fm$genes), by = "gene_safe", all.x = TRUE)
    setorder(imp, -importance)
    imp[, rank := seq_len(.N)]
  }
  imp[, `:=`(algorithm = row$algorithm, feature_set = row$feature_set,
             interpretation_class = row$interpretation_class)]
  if (!"mean_abs_shap" %in% names(imp)) imp[, mean_abs_shap := NA_real_]
  if (!"importance" %in% names(imp)) imp[, importance := mean_abs_shap]
  if (!"importance_sd" %in% names(imp)) imp[, importance_sd := NA_real_]
  importance_rows[[length(importance_rows) + 1L]] <- imp[, .(feature_set, interpretation_class, algorithm, rank, gene_symbol,
                                                             importance, mean_abs_shap, importance_sd, explanation_method)]
}
xgb_shap_rows <- list()
if (requireNamespace("xgboost", quietly = TRUE)) {
  for (fs in names(feature_sets)) {
    genes <- feature_sets[[fs]]
    mat <- make_safe_matrix(expr_all, samples, genes)
    sf <- scale_fit(mat$x)
    x_scaled <- scale_apply(mat$x, sf)
    xgb_model <- train_model("XGBoost", x_scaled, prft$y)
    shap_imp <- xgb_shap_importance(xgb_model, x_scaled, mat$genes)
    shap_imp[, `:=`(
      feature_set = fs,
      interpretation_class = feature_class[feature_set == fs, interpretation_class],
      algorithm = "XGBoost",
      importance = mean_abs_shap,
      importance_sd = NA_real_,
      explanation_method = "xgboost_predcontrib_SHAP"
    )]
    xgb_shap_rows[[length(xgb_shap_rows) + 1L]] <- shap_imp[, .(feature_set, interpretation_class, algorithm, rank, gene_symbol,
                                                               importance, mean_abs_shap, importance_sd, explanation_method)]
  }
}
importance_dt <- rbindlist(importance_rows, fill = TRUE)
if (length(xgb_shap_rows) > 0) {
  xgb_shap_dt <- rbindlist(xgb_shap_rows, fill = TRUE)
  importance_dt <- rbindlist(list(importance_dt, xgb_shap_dt), fill = TRUE)
} else {
  xgb_shap_dt <- data.table()
}
setorder(importance_dt, feature_set, rank)
fwrite(importance_dt, file.path(results_dir, "phase3B_fix_SHAP_or_importance_top_features.csv"))

orig6_imp <- importance_dt[gene_symbol %in% original6]
orig6_imp <- merge(data.table(gene_symbol = original6), orig6_imp, by = "gene_symbol", all.x = TRUE)
fwrite(orig6_imp, file.path(results_dir, "phase3B_fix_original_6gene_importance.csv"))

gpl <- parse_gpl_annotation(file.path(root_dir, "phase1_runtime", "00_raw_data", "geo_validation", "GPL570_family.soft.gz"))
geo_universe <- unique(c(unlist(feature_sets), gene_sets_all$Proteostasis_core, gene_sets_all$Ferroptosis_tolerance_set))
geo37642 <- align_geo_expr("GSE37642-GPL570_series_matrix.txt.gz", gpl, geo_universe)
geo12417 <- align_geo_expr("GSE12417-GPL570_series_matrix.txt.gz", gpl, geo_universe)
common_geo <- intersect(rownames(geo37642), rownames(geo12417))
geo_list <- list(
  GSE37642 = geo37642,
  GSE12417 = geo12417,
  combined_GPL570 = cbind(geo37642[common_geo, , drop = FALSE], geo12417[common_geo, , drop = FALSE])
)

external_rows <- list()
for (i in seq_len(nrow(best_by_fs))) {
  row <- best_by_fs[i]
  for (dataset_name in names(geo_list)) {
    geo_expr <- geo_list[[dataset_name]]
    label <- score_geo_prft(geo_expr, gene_sets_all)
    genes <- intersect(feature_sets[[row$feature_set]], intersect(rownames(expr_all), rownames(geo_expr)))
    if (length(genes) < 2) {
      external_rows[[length(external_rows) + 1L]] <- data.table(dataset = dataset_name, algorithm = row$algorithm,
                                                                feature_set = row$feature_set, n = nrow(label),
                                                                n_genes = length(genes), AUROC = NA_real_,
                                                                AUPRC = NA_real_, note = "Too few common genes.")
      next
    }
    train_mat <- make_safe_matrix(expr_all, samples, genes)
    sf <- scale_fit(train_mat$x)
    x_train <- scale_apply(train_mat$x, sf)
    model <- train_model(row$algorithm, x_train, prft$y)
    geo_samples <- intersect(label$sample_id, colnames(geo_expr))
    geo_x0 <- t(geo_expr[genes, geo_samples, drop = FALSE])
    colnames(geo_x0) <- train_mat$safe
    geo_x <- scale_apply(as.matrix(geo_x0), sf)
    prob <- predict_model(model, geo_x)
    label_sub <- label[match(geo_samples, sample_id)]
    met <- calc_metrics(label_sub$y, prob)
    external_rows[[length(external_rows) + 1L]] <- cbind(
      data.table(dataset = dataset_name, algorithm = row$algorithm, feature_set = row$feature_set,
                 interpretation_class = row$interpretation_class, n = nrow(label_sub), n_genes = length(genes),
                 note = "GEO label defined by within-cohort median approximate PRFT score using available gene-set expression."),
      met[, .(AUROC, AUPRC, accuracy, balanced_accuracy, F1, MCC)]
    )
  }
}
external_dt <- rbindlist(external_rows, fill = TRUE)
fwrite(external_dt, file.path(results_dir, "phase3B_fix_external_validation_by_feature_set.csv"))

gene_overlap <- merge(
  importance_dt[, .(phase3B_best_rank = min(rank, na.rm = TRUE),
                    phase3B_max_importance = max(importance, na.rm = TRUE),
                    phase3B_feature_sets = paste(unique(feature_set), collapse = ";")),
                by = gene_symbol],
  phase3a_freq,
  by = "gene_symbol",
  all = TRUE
)
gene_overlap[, `:=`(in_original_6gene = gene_symbol %in% original6,
                    in_core_axis = gene_symbol %in% core_axis_genes)]
setorder(gene_overlap, phase3B_best_rank, -recurrence_count)
fwrite(gene_overlap, file.path(results_dir, "phase3B_fix_phase3A_phase3B_overlap.csv"))

summary_dt <- merge(best_by_fs[, .(feature_set, interpretation_class, best_algorithm = algorithm,
                                   CV_AUROC = AUROC, CV_AUPRC = AUPRC, CV_F1 = F1, CV_MCC = MCC)],
                    external_dt[, .(external_mean_AUROC = mean(AUROC, na.rm = TRUE),
                                     external_mean_AUPRC = mean(AUPRC, na.rm = TRUE)),
                                by = feature_set],
                    by = "feature_set", all.x = TRUE)
summary_dt[, interpretation_boundary := fifelse(
  interpretation_class == "label-proximal explanatory",
  "Label-proximal/explanatory only; do not describe as independent validation or causal evidence.",
  "Independent-like feature set for PRFT-high state recognition; still not a survival/prognostic model."
)]
fwrite(summary_dt, file.path(results_dir, "phase3B_fix_independent_vs_label_proximal_summary.csv"))

theme_set(
  theme_classic(base_size = 8) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 7),
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", size = 9)
    )
)

p_heat <- ggplot(perf, aes(x = feature_set, y = algorithm, fill = AUROC)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = fmt_num(AUROC)), size = 2.2) +
  scale_fill_gradient(low = "#edf2f7", high = "#b23b3b", limits = c(0.5, 1), na.value = "grey90") +
  labs(x = NULL, y = NULL, fill = "AUROC", title = "Focused PRFT-high state recognition models") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 6.5), axis.text.y = element_text(size = 6.5))
ggsave(file.path(fig_dir, "phase3B_fix_model_performance_heatmap.pdf"), p_heat, width = 7.8, height = 4.3, device = cairo_pdf, bg = "white")

p_ext <- ggplot(external_dt, aes(x = feature_set, y = AUROC, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.65) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.35, colour = "grey50") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "External AUROC", title = "External PRFT-state recognition by feature set") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 6.5))
ggsave(file.path(fig_dir, "phase3B_fix_external_validation_barplot.pdf"), p_ext, width = 7.8, height = 4.4, device = cairo_pdf, bg = "white")

plot_importance <- function(fs, out_name, title) {
  dt <- importance_dt[feature_set == fs & explanation_method == "xgboost_predcontrib_SHAP"][order(rank)]
  if (nrow(dt) == 0) dt <- importance_dt[feature_set == fs][order(rank)]
  dt <- dt[1:min(.N, 25)]
  if (nrow(dt) == 0) dt <- data.table(gene_symbol = "NA", importance = 0)
  p <- ggplot(dt, aes(x = reorder(gene_symbol, importance), y = importance, fill = gene_symbol %in% original6)) +
    geom_col(width = 0.72) +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#b23b3b", "FALSE" = "#4c78a8")) +
    labs(x = NULL, y = unique(dt$explanation_method)[1], title = title)
  ggsave(file.path(fig_dir, out_name), p, width = 6.3, height = 4.8, device = cairo_pdf, bg = "white")
}
plot_importance("FS-A_original_6gene", "phase3B_fix_SHAP_or_importance_FS_A_6gene.pdf", "Original 6-gene state-recognition importance")
plot_importance("FS-B_33_cross_platform_candidates", "phase3B_fix_SHAP_or_importance_FS_B_or_E_33gene.pdf", "33-gene candidate state-recognition importance")
plot_importance("FS-F_core_axis_explanatory_genes", "phase3B_fix_SHAP_or_importance_FS_F_core_axis.pdf", "Core-axis label-proximal explanatory importance")

overlap_plot <- gene_overlap[is.finite(phase3B_max_importance) & is.finite(recurrence_count)]
overlap_plot <- overlap_plot[order(phase3B_best_rank)][1:min(.N, 40)]
p_overlap <- ggplot(overlap_plot, aes(x = recurrence_count, y = phase3B_max_importance,
                                      colour = in_original_6gene | in_core_axis, label = gene_symbol)) +
  geom_point(size = 2, alpha = 0.85) +
  ggrepel::geom_text_repel(data = overlap_plot[in_original_6gene == TRUE | in_core_axis == TRUE | phase3B_best_rank <= 10],
                           size = 2.2, max.overlaps = 35, show.legend = FALSE) +
  scale_colour_manual(values = c("TRUE" = "#b23b3b", "FALSE" = "#4c78a8")) +
  labs(x = "Phase 3A survival-ML recurrence count", y = "Phase 3B state-recognition importance",
       title = "Complementarity of survival-ML genes and PRFT-state features")
ggsave(file.path(fig_dir, "phase3B_fix_phase3A_phase3B_overlap.pdf"), p_overlap, width = 6.4, height = 4.8, device = cairo_pdf, bg = "white")

xgb_success <- requireNamespace("xgboost", quietly = TRUE) && any(perf$algorithm == "XGBoost")
explainer_success <- any(pkg_dt[package %in% c("fastshap", "DALEX", "iml") & available == TRUE, available])
true_shap_success <- xgb_success && any(importance_dt$explanation_method == "xgboost_predcontrib_SHAP")
explanation_method <- if (true_shap_success) "xgboost predcontrib SHAP" else "permutation importance fallback"

best_line <- function(fs) {
  row <- best_by_fs[feature_set == fs]
  if (nrow(row) == 0) return("NA")
  paste0(row$algorithm, " AUROC=", fmt_num(row$AUROC), ", AUPRC=", fmt_num(row$AUPRC))
}
external_line <- function(fs) {
  dt <- external_dt[feature_set == fs]
  if (nrow(dt) == 0) return("NA")
  paste(dt[, paste0(dataset, " AUROC=", fmt_num(AUROC), ", AUPRC=", fmt_num(AUPRC))], collapse = " | ")
}

importance_for_text <- if (exists("xgb_shap_dt") && nrow(xgb_shap_dt) > 0) xgb_shap_dt else importance_dt
independent_orig <- importance_for_text[gene_symbol %in% original6 & feature_set %in% c("FS-A_original_6gene", "FS-B_33_cross_platform_candidates", "FS-E_33_intersect_phase3A_recurrent")]
orig_importance_text <- if (nrow(independent_orig) > 0) {
  paste(independent_orig[order(feature_set, rank), paste0(feature_set, ":", gene_symbol, "(rank ", rank, ")")], collapse = "; ")
} else "not detected"
core_axis_text <- paste(importance_for_text[feature_set == "FS-F_core_axis_explanatory_genes"][order(rank)][1:min(.N, 10),
                                      paste0(gene_symbol, "(rank ", rank, ")")], collapse = "; ")

checklist <- c(
  paste0("1. XGBoost successful: ", ifelse(xgb_success, "yes", "no")),
  paste0("2. At least one of fastshap/DALEX/iml successful: ", ifelse(explainer_success, "yes", "no")),
  paste0("3. True SHAP successful: ", ifelse(true_shap_success, "yes", "no")),
  paste0("4. If no true SHAP, explanation method used: ", explanation_method),
  paste0("5. FS-A best model and AUROC/AUPRC: ", best_line("FS-A_original_6gene")),
  paste0("6. FS-B best model and AUROC/AUPRC: ", best_line("FS-B_33_cross_platform_candidates")),
  paste0("7. FS-C best model and AUROC/AUPRC: ", best_line("FS-C_phase3A_recurrence_top20")),
  paste0("8. FS-E best model and AUROC/AUPRC: ", best_line("FS-E_33_intersect_phase3A_recurrent")),
  paste0("9. FS-F best model and AUROC/AUPRC: ", best_line("FS-F_core_axis_explanatory_genes")),
  paste0("10. FS-A external validation: ", external_line("FS-A_original_6gene")),
  paste0("11. FS-B external validation: ", external_line("FS-B_33_cross_platform_candidates")),
  paste0("12. FS-C external validation: ", external_line("FS-C_phase3A_recurrence_top20")),
  paste0("13. FS-E external validation: ", external_line("FS-E_33_intersect_phase3A_recurrent")),
  paste0("14. FS-F external validation: ", external_line("FS-F_core_axis_explanatory_genes")),
  paste0("15. Original 6-gene importance in independent feature sets: ", orig_importance_text),
  paste0("16. Core-axis gene importance in FS-F: ", core_axis_text),
  "17. FS-F confirmed as label-proximal: yes; explanatory only and not independent validation or causal evidence.",
  "18. Phase3B suitable for main text: yes, if framed as PRFT-high state recognition and model interpretation with circularity caveat.",
  "19. Suggested panel: Figure 4/5 supplemental ML interpretability panel; keep FS-F as explanatory subpanel and independent-like FS-A/FS-B/FS-C/FS-E in main comparison.",
  "20. Recommend entering Phase3C: yes, after manual review of whether FS-F should be main-text or supplementary.",
  "21. Issues requiring manual confirmation: fastshap unavailable for R 4.5 binary repository; FS-F is label-proximal; avoid causal terms such as driver/proved/mediated/activated/reversed."
)
writeLines(checklist, file.path(log_dir, "phase3B_fix_key_result_checklist.txt"), useBytes = TRUE)

append_log("[Phase3B-fix] XGBoost successful: ", ifelse(xgb_success, "yes", "no"))
append_log("[Phase3B-fix] True SHAP successful: ", ifelse(true_shap_success, "yes", "no"))
append_log("[Phase3B-fix] Checklist written: phase3B_fix_key_result_checklist.txt")
append_log("[Phase3B-fix] Finished at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
