#!/usr/bin/env Rscript

# Archived old model, not used for main text.
# This script generated the superseded training-only 9-gene model.
# Human Genomics main-text reporting should use the final cross-platform
# six-gene workflow and outputs instead.

options(stringsAsFactors = FALSE)

local_lib <- normalizePath("17_tmp/R_libs", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
})
source("15_scripts/plot_label_utils.R")

if (!requireNamespace("glmnet", quietly = TRUE)) {
  stop("glmnet is not installed. Please install glmnet before running this script.")
}

suppressPackageStartupMessages({
  library(glmnet)
})

dir.create("07_signature", recursive = TRUE, showWarnings = FALSE)
dir.create("13_figures", recursive = TRUE, showWarnings = FALSE)
dir.create("14_tables", recursive = TRUE, showWarnings = FALSE)
dir.create("16_logs", recursive = TRUE, showWarnings = FALSE)

save_session_info <- function(path) {
  writeLines(capture.output(sessionInfo()), con = path)
}

zscore_cols <- function(mat) {
  out <- scale(mat)
  out[, apply(out, 2, function(x) all(is.na(x)))] <- 0
  out[is.na(out)] <- 0
  out
}

extract_nonzero_coefs <- function(fit_obj, s_value) {
  cf <- as.matrix(coef(fit_obj, s = s_value))
  cf_dt <- data.table(
    gene_symbol = rownames(cf),
    coefficient = as.numeric(cf[, 1])
  )
  cf_dt <- cf_dt[coefficient != 0]
  cf_dt
}

make_km_df <- function(survfit_obj) {
  s <- summary(survfit_obj)
  data.frame(
    time = s$time,
    surv = s$surv,
    strata = sub("^risk_group=", "", s$strata),
    stringsAsFactors = FALSE
  )
}

set.seed(20260622)

obj <- readRDS("02_processed_data/tcga_expr_clin_matched.rds")
prft_df <- as.data.table(readRDS("04_prft_score/tcga_prft_score.rds"))
main_candidates_csv <- fread("07_signature/main_candidate_genes_for_lasso.csv")
main_candidates_txt <- readLines("07_signature/main_candidate_genes_for_lasso.txt", warn = FALSE)
uni_integrated <- fread("07_signature/univariate_cox_integrated_results.csv")

expr <- as.matrix(obj$expr)
storage.mode(expr) <- "numeric"

candidate_genes_input <- unique(c(main_candidates_csv$gene_symbol, main_candidates_txt))
candidate_genes_found <- intersect(candidate_genes_input, rownames(expr))

samples_available_before_survival_filter <- ncol(expr)

clin_cols <- intersect(
  c("sample_id", "patient_id", "OS_time", "OS_status", "PRFT_score", "PRFT_group", "age", "sex"),
  colnames(prft_df)
)
clin_dt <- unique(prft_df[, ..clin_cols])
clin_dt <- clin_dt[!is.na(OS_time) & !is.na(OS_status)]
excluded_nonpositive_time <- clin_dt[OS_time <= 0, sample_id]
clin_dt <- clin_dt[OS_time > 0]

common_samples <- intersect(colnames(expr), clin_dt$sample_id)
clin_dt <- clin_dt[match(common_samples, sample_id)]
expr_sub <- expr[candidate_genes_found, common_samples, drop = FALSE]

if (!identical(colnames(expr_sub), clin_dt$sample_id)) {
  stop("Sample order mismatch after survival filtering.")
}

x <- t(expr_sub)
x_z <- zscore_cols(x)
rownames(x_z) <- rownames(x)
colnames(x_z) <- colnames(x)

y <- survival::Surv(clin_dt$OS_time, clin_dt$OS_status)

saveRDS(x_z, "07_signature/lasso_input_matrix_tcga.rds")
fwrite(clin_dt, "07_signature/lasso_input_clinical_tcga.csv")

cvfit <- cv.glmnet(
  x = x_z,
  y = y,
  family = "cox",
  alpha = 1,
  nfolds = 10,
  standardize = FALSE
)

fit <- glmnet(
  x = x_z,
  y = y,
  family = "cox",
  alpha = 1,
  standardize = FALSE
)

saveRDS(cvfit, "07_signature/lasso_cv_fit_tcga.rds")
saveRDS(fit, "07_signature/lasso_fit_tcga.rds")

lambda_min <- cvfit$lambda.min
lambda_1se <- cvfit$lambda.1se

coef_min <- extract_nonzero_coefs(cvfit, lambda_min)
coef_1se <- extract_nonzero_coefs(cvfit, lambda_1se)

fwrite(coef_1se, "07_signature/lasso_signature_coefficients_lambda_1se.csv")
fwrite(coef_min, "07_signature/lasso_signature_coefficients_lambda_min.csv")

selection_counts <- setNames(rep(0L, length(candidate_genes_found)), candidate_genes_found)
for (i in 1:100) {
  set.seed(20260622 + i)
  rep_cvfit <- cv.glmnet(
    x = x_z,
    y = y,
    family = "cox",
    alpha = 1,
    nfolds = 10,
    standardize = FALSE
  )
  rep_genes <- extract_nonzero_coefs(rep_cvfit, rep_cvfit$lambda.1se)$gene_symbol
  if (length(rep_genes) == 0) {
    rep_genes <- extract_nonzero_coefs(rep_cvfit, rep_cvfit$lambda.min)$gene_symbol
  }
  if (length(rep_genes) > 0) {
    selection_counts[rep_genes] <- selection_counts[rep_genes] + 1L
  }
}

freq_dt <- data.table(
  gene_symbol = names(selection_counts),
  selection_count = as.integer(selection_counts),
  selection_frequency = as.numeric(selection_counts) / 100
)
fwrite(freq_dt, "07_signature/lasso_repeated_selection_frequency.csv")

lambda_rule_used <- "lambda.1se"
lambda_used <- lambda_1se
chosen_coef <- copy(coef_1se)
if (nrow(chosen_coef) == 0) {
  lambda_rule_used <- "lambda.min because lambda.1se selected 0 genes"
  lambda_used <- lambda_min
  chosen_coef <- copy(coef_min)
}

final_coef <- merge(chosen_coef, freq_dt, by = "gene_symbol", all.x = TRUE, sort = FALSE)
final_coef <- merge(
  final_coef,
  uni_integrated[, .(gene_symbol, univariate_HR = HR, univariate_P.Value = P.Value, logFC, adj.P.Val, module_color, GS_PRFT, MM)],
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)

if (nrow(final_coef) > 10) {
  setorder(final_coef, -selection_frequency, univariate_P.Value, -univariate_HR, -GS_PRFT, -MM)
  final_coef <- final_coef[1:10]
  lambda_rule_used <- paste0("stability-filtered ", ifelse(lambda_used == lambda_1se, "lambda.1se", "lambda.min"), " model")
}

if (nrow(final_coef) < 3) {
  lambda_rule_used <- paste0(lambda_rule_used, "; final model has fewer than 3 genes")
}

final_coef[, abs_coefficient := abs(coefficient)]
setorder(final_coef, -abs_coefficient)
final_coef[, abs_coefficient := NULL]
fwrite(final_coef, "07_signature/final_prft_signature_coefficients.csv")
writeLines(final_coef$gene_symbol, "07_signature/final_prft_signature_genes.txt")

final_genes <- final_coef$gene_symbol
final_betas <- final_coef$coefficient
names(final_betas) <- final_genes

risk_score <- as.numeric(x_z[, final_genes, drop = FALSE] %*% final_betas)
median_cutoff <- stats::median(risk_score, na.rm = TRUE)
risk_group <- ifelse(risk_score >= median_cutoff, "high_risk", "low_risk")

risk_dt <- data.table(
  sample_id = clin_dt$sample_id,
  patient_id = if ("patient_id" %in% colnames(clin_dt)) clin_dt$patient_id else NA_character_,
  OS_time = clin_dt$OS_time,
  OS_status = clin_dt$OS_status,
  risk_score = risk_score,
  risk_group = risk_group,
  PRFT_score = if ("PRFT_score" %in% colnames(clin_dt)) clin_dt$PRFT_score else NA_real_,
  PRFT_group = if ("PRFT_group" %in% colnames(clin_dt)) clin_dt$PRFT_group else NA_character_,
  age = if ("age" %in% colnames(clin_dt)) clin_dt$age else NA_real_,
  sex = if ("sex" %in% colnames(clin_dt)) clin_dt$sex else NA_character_
)
fwrite(risk_dt, "07_signature/tcga_risk_score_by_sample.csv")

km_fit <- survfit(Surv(OS_time, OS_status) ~ risk_group, data = risk_dt)
km_diff <- survdiff(Surv(OS_time, OS_status) ~ risk_group, data = risk_dt)
logrank_p <- 1 - pchisq(km_diff$chisq, df = length(km_diff$n) - 1)

km_df <- make_km_df(km_fit)
p_km <- ggplot(km_df, aes(x = time, y = surv, color = strata)) +
  geom_step(linewidth = 1) +
  annotate("text", x = max(risk_dt$OS_time) * 0.65, y = 0.15, label = paste0("Log-rank P = ", signif(logrank_p, 3))) +
  scale_color_manual(
    values = c(high_risk = "#D62728", low_risk = "#1F77B4"),
    labels = pretty_label(c("high_risk", "low_risk"))
  ) +
  labs(x = "Time (days)", y = "Overall survival probability", color = NULL, title = "TCGA-LAML high vs low risk") +
  theme_bw(base_size = 11)
ggsave("13_figures/Figure5_TCGA_KM_high_low_risk.pdf", p_km, width = 7, height = 5.5)
ggsave("13_figures/Figure5_TCGA_KM_high_low_risk.png", p_km, width = 7, height = 5.5, dpi = 300)

uni_risk_fit <- coxph(Surv(OS_time, OS_status) ~ risk_score, data = risk_dt)
uni_risk_sum <- summary(uni_risk_fit)

multivar_formula_used <- NA_character_
multi_fit <- NULL
multivar_reason <- "risk_score + age + sex"
if ("age" %in% colnames(risk_dt) && sum(!is.na(risk_dt$age)) / nrow(risk_dt) >= 0.7) {
  if ("sex" %in% colnames(risk_dt) && sum(!is.na(risk_dt$sex)) / nrow(risk_dt) >= 0.7 && length(unique(na.omit(risk_dt$sex))) >= 2) {
    multi_dt <- risk_dt[complete.cases(risk_dt[, .(OS_time, OS_status, risk_score, age, sex)])]
    multi_dt$sex <- factor(multi_dt$sex)
    multi_fit <- coxph(Surv(OS_time, OS_status) ~ risk_score + age + sex, data = multi_dt)
    multivar_formula_used <- "risk_score + age + sex"
  } else {
    multi_dt <- risk_dt[complete.cases(risk_dt[, .(OS_time, OS_status, risk_score, age)])]
    multi_fit <- coxph(Surv(OS_time, OS_status) ~ risk_score + age, data = multi_dt)
    multivar_formula_used <- "risk_score + age"
    multivar_reason <- "sex missing or not modelable; used risk_score + age"
  }
} else {
  multi_dt <- risk_dt[complete.cases(risk_dt[, .(OS_time, OS_status, risk_score)])]
  multi_fit <- coxph(Surv(OS_time, OS_status) ~ risk_score, data = multi_dt)
  multivar_formula_used <- "risk_score only"
  multivar_reason <- "age missing severely; used risk_score only"
}
multi_sum <- summary(multi_fit)

c_index <- as.numeric(uni_risk_sum$concordance[1])

roc_method_used <- "not_run"
auc_1 <- NA_real_
auc_2 <- NA_real_
auc_3 <- NA_real_
roc_plot_created <- FALSE

if (requireNamespace("timeROC", quietly = TRUE)) {
  roc_obj <- tryCatch(
    timeROC::timeROC(
      T = risk_dt$OS_time,
      delta = risk_dt$OS_status,
      marker = risk_dt$risk_score,
      cause = 1,
      weighting = "marginal",
      times = c(365, 730, 1095),
      iid = FALSE
    ),
    error = function(e) e
  )
  if (!inherits(roc_obj, "error")) {
    roc_method_used <- "timeROC"
    auc_1 <- roc_obj$AUC[1]
    auc_2 <- roc_obj$AUC[2]
    auc_3 <- roc_obj$AUC[3]
    grDevices::pdf("13_figures/Figure5_TCGA_timeROC.pdf", width = 6, height = 6)
    plot(roc_obj, time = 365, col = "#1F77B4", title = FALSE)
    plot(roc_obj, time = 730, add = TRUE, col = "#D62728")
    plot(roc_obj, time = 1095, add = TRUE, col = "#2CA02C")
    legend("bottomright", legend = c(
      paste0("1-year AUC = ", round(auc_1, 3)),
      paste0("2-year AUC = ", round(auc_2, 3)),
      paste0("3-year AUC = ", round(auc_3, 3))
    ), col = c("#1F77B4", "#D62728", "#2CA02C"), lwd = 2, bty = "n")
    title(main = "TCGA-LAML time-dependent ROC")
    dev.off()
    grDevices::png("13_figures/Figure5_TCGA_timeROC.png", width = 1800, height = 1800, res = 300)
    plot(roc_obj, time = 365, col = "#1F77B4", title = FALSE)
    plot(roc_obj, time = 730, add = TRUE, col = "#D62728")
    plot(roc_obj, time = 1095, add = TRUE, col = "#2CA02C")
    legend("bottomright", legend = c(
      paste0("1-year AUC = ", round(auc_1, 3)),
      paste0("2-year AUC = ", round(auc_2, 3)),
      paste0("3-year AUC = ", round(auc_3, 3))
    ), col = c("#1F77B4", "#D62728", "#2CA02C"), lwd = 2, bty = "n")
    title(main = "TCGA-LAML time-dependent ROC")
    dev.off()
    roc_plot_created <- TRUE
  } else {
    roc_method_used <- paste0("timeROC_failed: ", conditionMessage(roc_obj))
  }
} else if (requireNamespace("survivalROC", quietly = TRUE)) {
  roc_method_used <- "survivalROC"
  roc1 <- tryCatch(survivalROC::survivalROC(Stime = risk_dt$OS_time, status = risk_dt$OS_status, marker = risk_dt$risk_score, predict.time = 365, method = "KM"), error = function(e) e)
  roc2 <- tryCatch(survivalROC::survivalROC(Stime = risk_dt$OS_time, status = risk_dt$OS_status, marker = risk_dt$risk_score, predict.time = 730, method = "KM"), error = function(e) e)
  roc3 <- tryCatch(survivalROC::survivalROC(Stime = risk_dt$OS_time, status = risk_dt$OS_status, marker = risk_dt$risk_score, predict.time = 1095, method = "KM"), error = function(e) e)
  if (!inherits(roc1, "error") && !inherits(roc2, "error") && !inherits(roc3, "error")) {
    auc_1 <- roc1$AUC
    auc_2 <- roc2$AUC
    auc_3 <- roc3$AUC
    grDevices::pdf("13_figures/Figure5_TCGA_timeROC.pdf", width = 6, height = 6)
    plot(roc1$FP, roc1$TP, type = "l", col = "#1F77B4", xlab = "False positive rate", ylab = "True positive rate", main = "TCGA-LAML time-dependent ROC")
    lines(roc2$FP, roc2$TP, col = "#D62728")
    lines(roc3$FP, roc3$TP, col = "#2CA02C")
    abline(0, 1, lty = 2, col = "grey50")
    legend("bottomright", legend = c(
      paste0("1-year AUC = ", round(auc_1, 3)),
      paste0("2-year AUC = ", round(auc_2, 3)),
      paste0("3-year AUC = ", round(auc_3, 3))
    ), col = c("#1F77B4", "#D62728", "#2CA02C"), lwd = 2, bty = "n")
    dev.off()
    grDevices::png("13_figures/Figure5_TCGA_timeROC.png", width = 1800, height = 1800, res = 300)
    plot(roc1$FP, roc1$TP, type = "l", col = "#1F77B4", xlab = "False positive rate", ylab = "True positive rate", main = "TCGA-LAML time-dependent ROC")
    lines(roc2$FP, roc2$TP, col = "#D62728")
    lines(roc3$FP, roc3$TP, col = "#2CA02C")
    abline(0, 1, lty = 2, col = "grey50")
    legend("bottomright", legend = c(
      paste0("1-year AUC = ", round(auc_1, 3)),
      paste0("2-year AUC = ", round(auc_2, 3)),
      paste0("3-year AUC = ", round(auc_3, 3))
    ), col = c("#1F77B4", "#D62728", "#2CA02C"), lwd = 2, bty = "n")
    dev.off()
    roc_plot_created <- TRUE
  } else {
    roc_method_used <- "survivalROC_failed"
  }
}

grDevices::pdf("13_figures/Figure5_LASSO_cv_curve.pdf", width = 7, height = 5.5)
plot(cvfit)
dev.off()
grDevices::png("13_figures/Figure5_LASSO_cv_curve.png", width = 2100, height = 1650, res = 300)
plot(cvfit)
dev.off()

grDevices::pdf("13_figures/Figure5_LASSO_coefficient_path.pdf", width = 7, height = 5.5)
plot(fit, xvar = "lambda", label = FALSE)
dev.off()
grDevices::png("13_figures/Figure5_LASSO_coefficient_path.png", width = 2100, height = 1650, res = 300)
plot(fit, xvar = "lambda", label = FALSE)
dev.off()

freq_plot_dt <- copy(freq_dt)
setorder(freq_plot_dt, -selection_frequency, -selection_count)
freq_plot_dt <- freq_plot_dt[1:min(20, .N)]
freq_plot_dt[, gene_symbol := factor(gene_symbol, levels = rev(gene_symbol))]
p_freq <- ggplot(freq_plot_dt, aes(x = gene_symbol, y = selection_frequency)) +
  geom_col(fill = "#4C78A8") +
  coord_flip() +
  labs(x = NULL, y = "Selection frequency", title = "Top repeated LASSO selection frequencies") +
  theme_bw(base_size = 11)
ggsave("13_figures/Figure5_LASSO_selection_frequency.pdf", p_freq, width = 7, height = 6)
ggsave("13_figures/Figure5_LASSO_selection_frequency.png", p_freq, width = 7, height = 6, dpi = 300)

risk_plot_dt <- copy(risk_dt)
setorder(risk_plot_dt, risk_score)
risk_plot_dt[, sample_index := .I]
p_risk <- ggplot(risk_plot_dt, aes(x = sample_index, y = risk_score, color = risk_group)) +
  geom_point(size = 1.8) +
  geom_hline(yintercept = median_cutoff, linetype = "dashed", color = "grey40") +
  scale_color_manual(
    values = c(high_risk = "#D62728", low_risk = "#1F77B4"),
    labels = pretty_label(c("high_risk", "low_risk"))
  ) +
  labs(x = "Sample rank", y = "Risk score", color = NULL, title = "TCGA-LAML risk score distribution") +
  theme_bw(base_size = 11)
ggsave("13_figures/Figure5_TCGA_risk_score_distribution.pdf", p_risk, width = 7, height = 5.5)
ggsave("13_figures/Figure5_TCGA_risk_score_distribution.png", p_risk, width = 7, height = 5.5, dpi = 300)

uni_risk_row <- data.table(
  model = "univariate",
  formula_used = "Surv(OS_time, OS_status) ~ risk_score",
  HR_risk_score = unname(uni_risk_sum$coefficients[1, "exp(coef)"]),
  lower95_risk_score = unname(uni_risk_sum$conf.int[1, "lower .95"]),
  upper95_risk_score = unname(uni_risk_sum$conf.int[1, "upper .95"]),
  P.Value_risk_score = unname(uni_risk_sum$coefficients[1, "Pr(>|z|)"])
)
multi_coef_row_idx <- which(rownames(multi_sum$coefficients) == "risk_score")[1]
multi_risk_row <- data.table(
  model = "multivariate",
  formula_used = multivar_formula_used,
  HR_risk_score = unname(multi_sum$coefficients[multi_coef_row_idx, "exp(coef)"]),
  lower95_risk_score = unname(multi_sum$conf.int[multi_coef_row_idx, "lower .95"]),
  upper95_risk_score = unname(multi_sum$conf.int[multi_coef_row_idx, "upper .95"]),
  P.Value_risk_score = unname(multi_sum$coefficients[multi_coef_row_idx, "Pr(>|z|)"])
)
cox_risk_dt <- rbindlist(list(uni_risk_row, multi_risk_row), fill = TRUE)
fwrite(cox_risk_dt, "07_signature/tcga_univariate_multivariate_cox_risk.csv")

perf_dt <- data.table(
  metric = c(
    "samples_used_for_lasso",
    "events_used_for_lasso",
    "lambda_min",
    "lambda_1se",
    "lambda_used",
    "logrank_P.Value",
    "univariate_cox_HR_risk_score",
    "univariate_cox_P.Value_risk_score",
    "multivariate_cox_HR_risk_score",
    "multivariate_cox_P.Value_risk_score",
    "C_index",
    "AUC_1year",
    "AUC_2year",
    "AUC_3year",
    "ROC_method_used",
    "multivariate_reason"
  ),
  value = c(
    nrow(risk_dt),
    sum(risk_dt$OS_status == 1, na.rm = TRUE),
    lambda_min,
    lambda_1se,
    lambda_used,
    logrank_p,
    uni_risk_row$HR_risk_score,
    uni_risk_row$P.Value_risk_score,
    multi_risk_row$HR_risk_score,
    multi_risk_row$P.Value_risk_score,
    c_index,
    auc_1,
    auc_2,
    auc_3,
    roc_method_used,
    multivar_reason
  )
)
fwrite(perf_dt, "07_signature/tcga_lasso_model_performance.csv")

summary_dt <- data.table(
  candidate_genes_input = length(candidate_genes_input),
  candidate_genes_found_in_expr = length(candidate_genes_found),
  samples_available_before_survival_filter = samples_available_before_survival_filter,
  samples_used_for_lasso = nrow(risk_dt),
  events_used_for_lasso = sum(risk_dt$OS_status == 1, na.rm = TRUE),
  lambda_min = lambda_min,
  lambda_1se = lambda_1se,
  lambda_used = lambda_used,
  lambda_rule_used = lambda_rule_used,
  genes_selected_lambda_min = nrow(coef_min),
  genes_selected_lambda_1se = nrow(coef_1se),
  final_signature_genes_count = nrow(final_coef),
  final_signature_genes = paste(final_coef$gene_symbol, collapse = ";"),
  median_risk_score_cutoff = median_cutoff,
  high_risk_samples = sum(risk_dt$risk_group == "high_risk", na.rm = TRUE),
  low_risk_samples = sum(risk_dt$risk_group == "low_risk", na.rm = TRUE),
  logrank_P.Value = logrank_p,
  univariate_cox_HR_risk_score = uni_risk_row$HR_risk_score,
  univariate_cox_P.Value_risk_score = uni_risk_row$P.Value_risk_score,
  multivariate_cox_HR_risk_score = multi_risk_row$HR_risk_score,
  multivariate_cox_P.Value_risk_score = multi_risk_row$P.Value_risk_score,
  C_index = c_index,
  AUC_1year = auc_1,
  AUC_2year = auc_2,
  AUC_3year = auc_3,
  ROC_method_used = roc_method_used
)
fwrite(summary_dt, "14_tables/tcga_lasso_signature_summary.csv")

save_session_info("16_logs/sessionInfo_13_lasso_cox_prft_signature_tcga.txt")
