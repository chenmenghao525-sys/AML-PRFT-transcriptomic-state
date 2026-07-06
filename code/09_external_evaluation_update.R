#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(1234)

ascii_lib <- Sys.getenv("PHASE1_ASCII_R_LIB", unset = "")
if (nzchar(ascii_lib) && dir.exists(ascii_lib)) .libPaths(c(ascii_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
})

root <- Sys.getenv("PHASE1_AUDIT_ROOT", unset = "")
if (!nzchar(root)) root <- getwd()
root <- chartr("\\", "/", root)
results_dir <- file.path(root, "03_results_tables")
log_dir <- file.path(root, "05_logs")
log_file <- file.path(log_dir, "phase3B_PRFT_state_ML_log.txt")

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
                      balanced_accuracy = NA_real_, F1 = NA_real_, MCC = NA_real_))
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
             F1 = f1, MCC = mcc)
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
  if (algorithm == "Logistic Regression baseline") {
    d <- data.frame(y = y_train, x_train, check.names = FALSE)
    fit <- stats::glm(y ~ ., data = d, family = stats::binomial())
    return(list(algorithm = algorithm, fit = fit))
  }
  stop("External validation updater supports top parametric/tree models only: ", algorithm)
}

predict_model <- function(model, x_test) {
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
  if (algorithm == "Logistic Regression baseline") {
    return(as.numeric(stats::predict(model$fit, newdata = data.frame(x_test, check.names = FALSE), type = "response")))
  }
  rep(NA_real_, nrow(x_test))
}

parse_series_matrix <- function(localfile) {
  lines <- readLines(gzfile(localfile), warn = FALSE)
  begin_idx <- grep("^!series_matrix_table_begin", lines)
  end_idx <- grep("^!series_matrix_table_end", lines)
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
  expr_mat
}

parse_gpl_annotation <- function(localfile) {
  lines <- readLines(gzfile(localfile), warn = FALSE)
  begin_idx <- grep("^!platform_table_begin", lines)
  end_idx <- grep("^!platform_table_end", lines)
  table_lines <- if (length(end_idx) == 0 || end_idx[1] <= begin_idx[1]) {
    lines[(begin_idx[1] + 1):length(lines)]
  } else {
    lines[(begin_idx[1] + 1):(end_idx[1] - 1)]
  }
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
  expr_mat <- parse_series_matrix(file.path(root, "phase1_runtime", "00_raw_data", "geo_validation", matrix_file))
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

append_log("[Phase3B-external] Updating external GEO PRFT-state validation without rerunning CV.")

expr_all <- readRDS(file.path(root, "phase1_runtime", "02_processed_data", "tcga_expr_hgnc_log2cpm.rds"))
prft <- as.data.table(readRDS(file.path(root, "phase1_runtime", "04_prft_score", "tcga_prft_score.rds")))
candidate33 <- fread(file.path(results_dir, "phase1_33_candidates.csv"))
candidate715 <- fread(file.path(results_dir, "phase1_715_candidates.csv"))
phase3a_freq <- fread(file.path(results_dir, "phase3A_fix_gene_recurrence_frequency.csv"))
gene_sets_all <- readRDS(file.path(root, "phase1_runtime", "03_gene_sets", "prft_gene_sets_all.rds"))
performance <- fread(file.path(results_dir, "phase3B_classification_model_performance.csv"))

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
  "FS-D_715_PRFT_correlation_top50" = intersect(head(candidate715[order(-abs(GS_PRFT)), gene_symbol], 50), rownames(expr_all)),
  "FS-E_33_intersect_phase3A_recurrent" = intersect(intersect(candidate33$gene_symbol, phase3a_recurrent), rownames(expr_all)),
  "FS-F_core_axis_explanatory_genes" = intersect(core_axis_genes, rownames(expr_all))
)

gpl <- parse_gpl_annotation(file.path(root, "phase1_runtime", "00_raw_data", "geo_validation", "GPL570_family.soft.gz"))
geo_universe <- unique(c(unlist(feature_sets), gene_sets_all$Proteostasis_core, gene_sets_all$Ferroptosis_tolerance_set))
geo37642 <- align_geo_expr("GSE37642-GPL570_series_matrix.txt.gz", gpl, geo_universe)
geo12417 <- align_geo_expr("GSE12417-GPL570_series_matrix.txt.gz", gpl, geo_universe)
common <- intersect(rownames(geo37642), rownames(geo12417))
geo_list <- list(
  GSE37642 = geo37642,
  GSE12417 = geo12417,
  combined_GPL570 = cbind(geo37642[common, , drop = FALSE], geo12417[common, , drop = FALSE])
)

top_for_external <- performance[converged == TRUE & algorithm %in% c("Elastic Net Logistic Regression", "Random Forest", "SVM radial", "GBM", "Logistic Regression baseline")]
top_for_external <- top_for_external[order(-AUROC, -AUPRC)][1:min(.N, 5)]
external_rows <- list()

for (i in seq_len(nrow(top_for_external))) {
  row <- top_for_external[i]
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
    model <- tryCatch(train_model(row$algorithm, x_train, prft$y), error = function(e) e)
    if (inherits(model, "error")) {
      external_rows[[length(external_rows) + 1L]] <- data.table(dataset = dataset_name, algorithm = row$algorithm,
                                                                feature_set = row$feature_set, n = nrow(label),
                                                                n_genes = length(genes), AUROC = NA_real_,
                                                                AUPRC = NA_real_,
                                                                note = paste("Training failed:", conditionMessage(model)))
      next
    }
    geo_samples <- intersect(label$sample_id, colnames(geo_expr))
    geo_x0 <- t(geo_expr[genes, geo_samples, drop = FALSE])
    colnames(geo_x0) <- train_mat$safe
    geo_x <- scale_apply(as.matrix(geo_x0), sf)
    prob <- tryCatch(predict_model(model, geo_x), error = function(e) e)
    if (inherits(prob, "error")) {
      external_rows[[length(external_rows) + 1L]] <- data.table(dataset = dataset_name, algorithm = row$algorithm,
                                                                feature_set = row$feature_set, n = length(geo_samples),
                                                                n_genes = length(genes), AUROC = NA_real_,
                                                                AUPRC = NA_real_,
                                                                note = paste("Prediction failed:", conditionMessage(prob)))
      next
    }
    label_sub <- label[match(geo_samples, sample_id)]
    met <- calc_metrics(label_sub$y, prob)
    external_rows[[length(external_rows) + 1L]] <- cbind(
      data.table(dataset = dataset_name, algorithm = row$algorithm, feature_set = row$feature_set,
                 n = nrow(label_sub), n_genes = length(genes),
                 note = "GEO label defined by within-cohort median approximate PRFT score using available gene-set expression."),
      met
    )
  }
}

external_dt <- rbindlist(external_rows, fill = TRUE)
fwrite(external_dt, file.path(results_dir, "phase3B_external_state_validation_if_available.csv"))

external_done <- nrow(external_dt) > 0 && any(is.finite(external_dt$AUROC))
external_summary <- if (external_done) {
  paste(external_dt[is.finite(AUROC), paste0(dataset, " ", algorithm, "+", feature_set,
                                             " AUROC=", fmt_num(AUROC), ", AUPRC=", fmt_num(AUPRC))],
        collapse = " | ")
} else {
  "not completed"
}

checklist_path <- file.path(log_dir, "phase3B_key_result_checklist.txt")
checklist <- readLines(checklist_path, warn = FALSE)
checklist[grepl("^12\\.", checklist)] <- paste0("12. External GEO state-recognition validation completed: ", ifelse(external_done, "yes", "no"))
checklist[grepl("^13\\.", checklist)] <- paste0("13. External validation AUROC/AUPRC if available: ", external_summary)
writeLines(checklist, checklist_path, useBytes = TRUE)

append_log("[Phase3B-external] External GEO PRFT-state validation updated: ", ifelse(external_done, "completed", "not completed"))
append_log("[Phase3B-external] ", external_summary)
