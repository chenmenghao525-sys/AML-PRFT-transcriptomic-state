#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(glmnet)
  library(ggplot2)
  library(survivalROC)
})
source("15_scripts/plot_label_utils.R")

dir.create("07_signature", recursive = TRUE, showWarnings = FALSE)
dir.create("13_figures", recursive = TRUE, showWarnings = FALSE)
dir.create("14_tables", recursive = TRUE, showWarnings = FALSE)
dir.create("16_logs", recursive = TRUE, showWarnings = FALSE)

save_session_info <- function(path) {
  writeLines(capture.output(sessionInfo()), con = path)
}

required_files <- c(
  "02_processed_data/tcga_expr_clin_matched.rds",
  "04_prft_score/tcga_prft_score.rds",
  "08_validation/cross_platform_rebuild_candidate_genes.csv",
  "08_validation/cross_platform_rebuild_recommended_strategy.csv",
  "07_signature/univariate_cox_integrated_results.csv",
  "05_deg/deg_prft_high_vs_low_tcga_all.csv",
  "06_wgcna/wgcna_gene_level_statistics.csv"
)
missing_required <- required_files[!file.exists(required_files)]
if (length(missing_required) > 0) {
  stop("Missing required files: ", paste(missing_required, collapse = "; "))
}

if (!requireNamespace("glmnet", quietly = TRUE)) {
  stop("glmnet is required but not installed.")
}
if (!requireNamespace("survivalROC", quietly = TRUE)) {
  stop("survivalROC is required but not installed.")
}

set.seed(20260622)

find_first_col <- function(dt, candidates) {
  nms <- colnames(dt)
  idx <- match(tolower(candidates), tolower(nms))
  idx <- idx[!is.na(idx)][1]
  if (is.na(idx)) NA_character_ else nms[idx]
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

zscore_rows <- function(mat) {
  z <- t(apply(mat, 1, zscore_vector))
  rownames(z) <- rownames(mat)
  colnames(z) <- colnames(mat)
  z
}

order_by_priority <- function(dt, p_col, hr_col, gs_col, mm_col) {
  dt[order(get(p_col), -get(hr_col), -get(gs_col), -get(mm_col))]
}

nonzero_genes_from_coef <- function(coef_obj) {
  cf <- as.matrix(coef_obj)
  genes <- rownames(cf)[abs(cf[, 1]) > 0]
  genes[!is.na(genes) & nzchar(genes)]
}

extract_coef_table <- function(fit, s_value) {
  cf <- as.matrix(coef(fit, s = s_value))
  out <- data.table(
    gene_symbol = rownames(cf),
    coefficient = as.numeric(cf[, 1])
  )
  out <- out[abs(coefficient) > 0]
  out <- out[order(-abs(coefficient))]
  out
}

choose_lambda_rule <- function(cv_fit) {
  beta_mat <- as.matrix(cv_fit$glmnet.fit$beta)
  lambda_df <- data.table(
    lambda = cv_fit$glmnet.fit$lambda,
    cvm = cv_fit$cvm,
    cvsd = cv_fit$cvsd,
    gene_count = colSums(beta_mat != 0)
  )

  lambda_min <- cv_fit$lambda.min
  lambda_1se <- cv_fit$lambda.1se

  count_min <- lambda_df[which.min(abs(lambda - lambda_min))]$gene_count
  count_1se <- lambda_df[which.min(abs(lambda - lambda_1se))]$gene_count

  if (count_1se >= 3 && count_1se <= 10) {
    return(list(lambda_used = lambda_1se, lambda_rule_used = "lambda.1se", lambda_table = lambda_df))
  }

  if (count_min >= 3 && count_min <= 10) {
    return(list(lambda_used = lambda_min, lambda_rule_used = "lambda.min", lambda_table = lambda_df))
  }

  candidate_df <- lambda_df[gene_count >= 3 & gene_count <= 10]
  if (nrow(candidate_df) > 0) {
    setorder(candidate_df, cvm, gene_count, -lambda)
    return(list(
      lambda_used = candidate_df$lambda[1],
      lambda_rule_used = "searched_path_3_to_10_genes",
      lambda_table = lambda_df
    ))
  }

  if (count_1se == 0) {
    return(list(lambda_used = lambda_min, lambda_rule_used = "lambda.min_fallback_after_lambda.1se_zero_genes", lambda_table = lambda_df))
  }

  list(lambda_used = lambda_min, lambda_rule_used = "lambda.min_fallback_no_3_to_10_gene_model", lambda_table = lambda_df)
}

calc_time_roc <- function(time_point, surv_time, surv_status, marker) {
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
  out
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

obj <- readRDS("02_processed_data/tcga_expr_clin_matched.rds")
prft <- readRDS("04_prft_score/tcga_prft_score.rds")
rebuild_candidates <- fread("08_validation/cross_platform_rebuild_candidate_genes.csv")
strategy_dt <- fread("08_validation/cross_platform_rebuild_recommended_strategy.csv")
cox_dt <- fread("07_signature/univariate_cox_integrated_results.csv")
deg_dt <- fread("05_deg/deg_prft_high_vs_low_tcga_all.csv")
wgcna_dt <- fread("06_wgcna/wgcna_gene_level_statistics.csv")

expr <- obj$expr
clin <- as.data.table(obj$clin)
prft_dt <- as.data.table(prft)

recommended_scenario_value <- strategy_dt$recommended_scenario[1]
if (is.na(recommended_scenario_value) || !nzchar(recommended_scenario_value)) {
  stop("recommended_scenario is missing in cross_platform_rebuild_recommended_strategy.csv")
}

scenario_dt <- rebuild_candidates[recommended_scenario == recommended_scenario_value & covered_in_selected_scenario == TRUE]
if (nrow(scenario_dt) == 0) {
  stop("No genes found for recommended_scenario: ", recommended_scenario_value)
}

col_hr <- find_first_col(scenario_dt, c("HR"))
col_cox_p <- find_first_col(scenario_dt, c("Cox_P.Value", "P.Value", "cox_P.Value"))
col_cox_fdr <- find_first_col(scenario_dt, c("Cox_FDR", "FDR", "cox_FDR"))
col_gs <- find_first_col(scenario_dt, c("GS_PRFT"))
col_mm <- find_first_col(scenario_dt, c("MM", "module_membership", "kME"))
col_adjp <- find_first_col(scenario_dt, c("adj.P.Val"))
col_logfc <- find_first_col(scenario_dt, c("logFC"))

if (any(is.na(c(col_hr, col_cox_p, col_gs, col_mm, col_adjp, col_logfc)))) {
  stop("Failed to identify one or more required columns in cross_platform_rebuild_candidate_genes.csv")
}

cross_platform_pool_all <- unique(scenario_dt$gene_symbol)
cross_platform_pool_risk <- unique(scenario_dt[get(col_hr) > 1 & get(col_cox_p) < 0.05]$gene_symbol)
cross_platform_pool_strict <- unique(
  scenario_dt[
    get(col_hr) > 1 &
      get(col_cox_p) < 0.05 &
      get(col_logfc) > 0.5 &
      get(col_adjp) < 0.05 &
      get(col_gs) > 0.30 &
      get(col_mm) > 0.50
  ]$gene_symbol
)

fwrite(data.table(gene_symbol = sort(cross_platform_pool_all)), "07_signature/cross_platform_lasso_input_genes_all.csv")
fwrite(data.table(gene_symbol = sort(cross_platform_pool_strict)), "07_signature/cross_platform_lasso_input_genes_strict.csv")

lasso_input_rule_used <- if (length(cross_platform_pool_strict) >= 20 && length(cross_platform_pool_strict) <= 80) {
  "strict_pool_used"
} else {
  "risk_pool_used"
}

lasso_input_genes <- if (lasso_input_rule_used == "strict_pool_used") {
  cross_platform_pool_strict
} else {
  cross_platform_pool_risk
}

input_dt <- scenario_dt[match(lasso_input_genes, gene_symbol)]
if (nrow(input_dt) == 0) {
  stop("No LASSO input genes available after applying input rule.")
}

if (lasso_input_rule_used == "risk_pool_used" && nrow(input_dt) > 80) {
  input_dt <- order_by_priority(input_dt, col_cox_p, col_hr, col_gs, col_mm)
  input_dt <- input_dt[1:80]
  lasso_input_rule_used <- "risk_pool_used_top80"
}

lasso_input_genes <- unique(input_dt$gene_symbol)

expr_genes_found <- intersect(lasso_input_genes, rownames(expr))
input_used_dt <- input_dt[gene_symbol %in% expr_genes_found]
input_used_dt <- order_by_priority(input_used_dt, col_cox_p, col_hr, col_gs, col_mm)
fwrite(input_used_dt, "07_signature/cross_platform_lasso_input_genes_used.csv")

if (length(expr_genes_found) == 0) {
  stop("No LASSO input genes found in TCGA expression matrix.")
}

clin_before <- copy(clin)
samples_available_before_survival_filter <- nrow(clin_before)

clin_model <- clin_before[!is.na(OS_time) & !is.na(OS_status) & OS_time > 0]
excluded_samples_survival <- samples_available_before_survival_filter - nrow(clin_model)

if (nrow(clin_model) == 0) {
  stop("No samples available after survival filtering.")
}

events_used_for_lasso <- sum(clin_model$OS_status == 1, na.rm = TRUE)
samples_used_for_lasso <- nrow(clin_model)

sample_ids <- clin_model$sample_id
expr_sub <- expr[expr_genes_found, sample_ids, drop = FALSE]
expr_z <- zscore_rows(expr_sub)

clin_model_out <- copy(clin_model)
fwrite(clin_model_out, "07_signature/cross_platform_lasso_input_clinical_tcga.csv")
saveRDS(list(expr = expr_z, clin = clin_model_out), "07_signature/cross_platform_lasso_input_matrix_tcga.rds")

x <- t(expr_z)
x <- as.matrix(x)
storage.mode(x) <- "numeric"
y <- survival::Surv(time = clin_model$OS_time, event = clin_model$OS_status)

cv_fit <- cv.glmnet(
  x = x,
  y = y,
  family = "cox",
  nfolds = 10,
  alpha = 1,
  standardize = FALSE,
  cox.ties = "breslow"
)
glmnet_fit <- cv_fit$glmnet.fit

saveRDS(cv_fit, "07_signature/cross_platform_lasso_cv_fit_tcga.rds")
saveRDS(glmnet_fit, "07_signature/cross_platform_lasso_fit_tcga.rds")

lambda_choice <- choose_lambda_rule(cv_fit)
lambda_used <- lambda_choice$lambda_used
lambda_rule_used <- lambda_choice$lambda_rule_used

coef_lambda_min_dt <- extract_coef_table(glmnet_fit, cv_fit$lambda.min)
coef_lambda_1se_dt <- extract_coef_table(glmnet_fit, cv_fit$lambda.1se)
fwrite(coef_lambda_min_dt, "07_signature/cross_platform_signature_coefficients_lambda_min.csv")
fwrite(coef_lambda_1se_dt, "07_signature/cross_platform_signature_coefficients_lambda_1se.csv")

repeated_selection <- vector("list", 100)
selection_count <- integer(length(expr_genes_found))
names(selection_count) <- expr_genes_found

for (i in seq_len(100)) {
  set.seed(20260622 + i)
  cv_rep <- cv.glmnet(
    x = x,
    y = y,
    family = "cox",
    nfolds = 10,
    alpha = 1,
    standardize = FALSE,
    cox.ties = "breslow"
  )
  choice_rep <- choose_lambda_rule(cv_rep)
  genes_rep <- nonzero_genes_from_coef(coef(cv_rep$glmnet.fit, s = choice_rep$lambda_used))
  if (length(genes_rep) > 0) {
    selection_count[genes_rep] <- selection_count[genes_rep] + 1L
  }
}

freq_dt <- data.table(
  gene_symbol = names(selection_count),
  selection_count = as.integer(selection_count),
  selection_frequency = as.numeric(selection_count) / 100
)
setorder(freq_dt, -selection_frequency, gene_symbol)
fwrite(freq_dt, "07_signature/cross_platform_lasso_repeated_selection_frequency.csv")

final_coef_dt <- extract_coef_table(glmnet_fit, lambda_used)
final_lambda_gene_count <- nrow(final_coef_dt)

final_model_rule <- lambda_rule_used
if (final_lambda_gene_count > 10) {
  freq_sub <- freq_dt[gene_symbol %in% final_coef_dt$gene_symbol]
  setorder(freq_sub, -selection_frequency, gene_symbol)
  keep_genes <- head(freq_sub$gene_symbol, 10)
  x_top <- x[, keep_genes, drop = FALSE]
  refit_df <- data.frame(
    OS_time = clin_model$OS_time,
    OS_status = clin_model$OS_status,
    x_top,
    check.names = FALSE
  )
  cox_formula <- as.formula(
    paste("Surv(OS_time, OS_status) ~", paste(sprintf("`%s`", keep_genes), collapse = " + "))
  )
  refit <- coxph(cox_formula, data = refit_df)
  refit_coef <- stats::coef(refit)
  final_coef_dt <- data.table(
    gene_symbol = names(refit_coef),
    coefficient = as.numeric(refit_coef)
  )
  final_coef_dt <- final_coef_dt[order(-abs(coefficient))]
  final_model_rule <- paste0(lambda_rule_used, "_stability_filtered_top10")
}

fwrite(final_coef_dt, "07_signature/final_cross_platform_prft_signature_coefficients.csv")
writeLines(final_coef_dt$gene_symbol, con = "07_signature/final_cross_platform_prft_signature_genes.txt")

if (nrow(final_coef_dt) == 0) {
  stop("Final cross-platform signature contains zero genes.")
}

coef_named <- final_coef_dt$coefficient
names(coef_named) <- final_coef_dt$gene_symbol
final_gene_mat <- expr_z[final_coef_dt$gene_symbol, , drop = FALSE]
risk_score <- rowSums(t(final_gene_mat) * coef_named)

risk_cutoff <- stats::median(risk_score, na.rm = TRUE)
risk_group <- ifelse(risk_score >= risk_cutoff, "high_risk", "low_risk")

risk_dt <- data.table(
  sample_id = sample_ids,
  patient_id = clin_model$patient_id,
  OS_time = clin_model$OS_time,
  OS_status = clin_model$OS_status,
  age = clin_model$age,
  sex = clin_model$sex,
  FAB = clin_model$FAB,
  WBC = clin_model$WBC,
  risk_score = as.numeric(risk_score),
  risk_group = risk_group
)
fwrite(risk_dt, "07_signature/tcga_cross_platform_risk_score_by_sample.csv")

surv_obj <- Surv(risk_dt$OS_time, risk_dt$OS_status)
fit_km <- survfit(surv_obj ~ risk_group, data = risk_dt)
logrank <- survdiff(surv_obj ~ risk_group, data = risk_dt)
logrank_p <- 1 - pchisq(logrank$chisq, df = length(logrank$n) - 1)

uni_fit <- coxph(surv_obj ~ risk_score, data = risk_dt)
uni_sum <- summary(uni_fit)

sex_non_missing <- sum(!is.na(risk_dt$sex))
sex_levels <- unique(na.omit(as.character(risk_dt$sex)))
age_missing_rate <- mean(is.na(risk_dt$age))
multiv_formula_label <- "risk_score + age + sex"

if (sex_non_missing == 0 || length(sex_levels) < 2) {
  multiv_formula <- surv_obj ~ risk_score + age
  multiv_formula_label <- "risk_score + age"
} else if (age_missing_rate > 0.30) {
  multiv_formula <- surv_obj ~ risk_score
  multiv_formula_label <- "risk_score_only_age_missing_severe"
} else {
  multiv_formula <- surv_obj ~ risk_score + age + sex
}

multiv_data <- copy(risk_dt)
multiv_fit <- coxph(multiv_formula, data = multiv_data)
multiv_sum <- summary(multiv_fit)

if ("risk_score" %in% rownames(multiv_sum$coefficients)) {
  multiv_hr <- as.numeric(multiv_sum$coefficients["risk_score", "exp(coef)"])
  multiv_p <- as.numeric(multiv_sum$coefficients["risk_score", "Pr(>|z|)"])
} else {
  multiv_hr <- NA_real_
  multiv_p <- NA_real_
}

c_index_obj <- survival::concordance(surv_obj ~ risk_score, data = risk_dt, reverse = TRUE)
c_index <- as.numeric(c_index_obj$concordance)

roc_1 <- calc_time_roc(365, risk_dt$OS_time, risk_dt$OS_status, risk_dt$risk_score)
roc_2 <- calc_time_roc(730, risk_dt$OS_time, risk_dt$OS_status, risk_dt$risk_score)
roc_3 <- calc_time_roc(1095, risk_dt$OS_time, risk_dt$OS_status, risk_dt$risk_score)

auc_1 <- if (!is.null(roc_1)) roc_1$AUC else NA_real_
auc_2 <- if (!is.null(roc_2)) roc_2$AUC else NA_real_
auc_3 <- if (!is.null(roc_3)) roc_3$AUC else NA_real_
roc_method_used <- "survivalROC"

perf_dt <- data.table(
  recommended_scenario = recommended_scenario_value,
  lasso_input_rule_used = lasso_input_rule_used,
  lambda_rule_used = final_model_rule,
  lambda_min = cv_fit$lambda.min,
  lambda_1se = cv_fit$lambda.1se,
  lambda_used = lambda_used,
  samples_used_for_lasso = samples_used_for_lasso,
  events_used_for_lasso = events_used_for_lasso,
  final_signature_genes_count = nrow(final_coef_dt),
  final_signature_genes = paste(final_coef_dt$gene_symbol, collapse = ";"),
  median_risk_score_cutoff = risk_cutoff,
  high_risk_samples = sum(risk_dt$risk_group == "high_risk"),
  low_risk_samples = sum(risk_dt$risk_group == "low_risk"),
  logrank_P.Value = logrank_p,
  univariate_cox_HR_risk_score = as.numeric(uni_sum$coefficients["risk_score", "exp(coef)"]),
  univariate_cox_P.Value_risk_score = as.numeric(uni_sum$coefficients["risk_score", "Pr(>|z|)"]),
  multivariate_cox_HR_risk_score = multiv_hr,
  multivariate_cox_P.Value_risk_score = multiv_p,
  multivariate_model_used = multiv_formula_label,
  C_index = c_index,
  AUC_1year = auc_1,
  AUC_2year = auc_2,
  AUC_3year = auc_3,
  ROC_method_used = roc_method_used
)
fwrite(perf_dt, "07_signature/tcga_cross_platform_lasso_model_performance.csv")

cox_risk_dt <- data.table(
  model = c("univariate", "multivariate"),
  formula = c("Surv(OS_time, OS_status) ~ risk_score", multiv_formula_label),
  HR_risk_score = c(
    as.numeric(uni_sum$coefficients["risk_score", "exp(coef)"]),
    multiv_hr
  ),
  P.Value_risk_score = c(
    as.numeric(uni_sum$coefficients["risk_score", "Pr(>|z|)"]),
    multiv_p
  )
)
fwrite(cox_risk_dt, "07_signature/tcga_cross_platform_univariate_multivariate_cox_risk.csv")

gse37642_cov <- all(vapply(final_coef_dt$gene_symbol, function(g) {
  any(rebuild_candidates$gene_symbol == g & grepl("GSE37642_GPL570", rebuild_candidates$covered_in_datasets, fixed = TRUE))
}, logical(1)))
gse12417_cov <- all(vapply(final_coef_dt$gene_symbol, function(g) {
  any(rebuild_candidates$gene_symbol == g & grepl("GSE12417_GPL570", rebuild_candidates$covered_in_datasets, fixed = TRUE))
}, logical(1)))

summary_dt <- data.table(
  recommended_scenario = recommended_scenario_value,
  datasets_for_future_validation = strategy_dt$datasets_for_future_validation[1],
  cross_platform_pool_all_genes = length(cross_platform_pool_all),
  cross_platform_pool_risk_genes = length(cross_platform_pool_risk),
  cross_platform_pool_strict_genes = length(cross_platform_pool_strict),
  lasso_input_rule_used = lasso_input_rule_used,
  lasso_input_genes_count = length(lasso_input_genes),
  candidate_genes_found_in_expr = length(expr_genes_found),
  samples_available_before_survival_filter = samples_available_before_survival_filter,
  samples_used_for_lasso = samples_used_for_lasso,
  events_used_for_lasso = events_used_for_lasso,
  lambda_min = cv_fit$lambda.min,
  lambda_1se = cv_fit$lambda.1se,
  lambda_used = lambda_used,
  lambda_rule_used = final_model_rule,
  genes_selected_lambda_min = nrow(coef_lambda_min_dt),
  genes_selected_lambda_1se = nrow(coef_lambda_1se_dt),
  final_signature_genes_count = nrow(final_coef_dt),
  final_signature_genes = paste(final_coef_dt$gene_symbol, collapse = ";"),
  all_final_genes_covered_in_GSE37642_GPL570 = gse37642_cov,
  all_final_genes_covered_in_GSE12417_GPL570 = gse12417_cov,
  standardization_for_validation = "cohort-wise gene-level z-score",
  median_risk_score_cutoff = risk_cutoff,
  high_risk_samples = sum(risk_dt$risk_group == "high_risk"),
  low_risk_samples = sum(risk_dt$risk_group == "low_risk"),
  logrank_P.Value = logrank_p,
  univariate_cox_HR_risk_score = as.numeric(uni_sum$coefficients["risk_score", "exp(coef)"]),
  univariate_cox_P.Value_risk_score = as.numeric(uni_sum$coefficients["risk_score", "Pr(>|z|)"]),
  multivariate_cox_HR_risk_score = multiv_hr,
  multivariate_cox_P.Value_risk_score = multiv_p,
  C_index = c_index,
  AUC_1year = auc_1,
  AUC_2year = auc_2,
  AUC_3year = auc_3,
  ROC_method_used = roc_method_used
)
fwrite(summary_dt, "14_tables/tcga_cross_platform_lasso_signature_summary.csv")

plot_to_pdf_png("13_figures/Figure6_cross_platform_LASSO_cv_curve", 7, 6, function() {
  plot(cv_fit)
  title(main = "Cross-platform compatible LASSO-Cox CV curve")
})

plot_to_pdf_png("13_figures/Figure6_cross_platform_LASSO_coefficient_path", 7, 6, function() {
  plot(glmnet_fit, xvar = "lambda", label = FALSE)
  title(main = "Cross-platform compatible LASSO coefficient path")
})

freq_plot_dt <- copy(freq_dt)
freq_plot_dt[, gene_symbol := factor(gene_symbol, levels = rev(head(gene_symbol, 30)))]
plot_to_pdf_png("13_figures/Figure6_cross_platform_LASSO_selection_frequency", 8, 7, function() {
  top_dt <- copy(freq_dt[1:min(30, .N)])
  top_dt[, gene_symbol := factor(gene_symbol, levels = rev(gene_symbol))]
  p <- ggplot(top_dt, aes(x = gene_symbol, y = selection_frequency)) +
    geom_col(fill = "#2C7FB8") +
    coord_flip() +
    labs(x = NULL, y = "Selection frequency", title = "Repeated LASSO selection frequency") +
    theme_bw(base_size = 12)
  print(p)
})

plot_to_pdf_png("13_figures/Figure6_TCGA_cross_platform_risk_score_distribution", 8, 6, function() {
  plot_dt <- copy(risk_dt)
  plot_dt[, sample_rank := seq_len(.N)]
  setorder(plot_dt, risk_score)
  plot_dt[, sample_rank := seq_len(.N)]
  p <- ggplot(plot_dt, aes(x = sample_rank, y = risk_score, color = risk_group)) +
    geom_point(size = 2) +
    geom_hline(yintercept = risk_cutoff, linetype = 2) +
    scale_color_manual(
      values = c(high_risk = "#D7301F", low_risk = "#225EA8"),
      labels = pretty_label(c("high_risk", "low_risk"))
    ) +
    labs(x = "Samples ranked by risk score", y = "Risk score", color = NULL, title = "TCGA cross-platform risk score distribution") +
    theme_bw(base_size = 12)
  print(p)
})

plot_to_pdf_png("13_figures/Figure6_TCGA_cross_platform_KM_high_low_risk", 7, 6, function() {
  km_cols <- c("high_risk" = "#D7301F", "low_risk" = "#225EA8")
  plot(fit_km, col = km_cols[names(km_cols)], lwd = 2, xlab = "Days", ylab = "Overall survival probability", main = "TCGA cross-platform signature")
  legend("bottomleft", legend = pretty_label(c("high_risk", "low_risk")), col = km_cols, lwd = 2, bty = "n")
  text(
    x = max(risk_dt$OS_time, na.rm = TRUE) * 0.6,
    y = 0.2,
    labels = paste0("Log-rank P = ", format(logrank_p, digits = 3, scientific = TRUE))
  )
})

plot_to_pdf_png("13_figures/Figure6_TCGA_cross_platform_timeROC", 7, 6, function() {
  plot(0:1, 0:1, type = "n", xlab = "False positive rate", ylab = "True positive rate", main = "TCGA cross-platform time-dependent ROC")
  abline(0, 1, lty = 2, col = "grey60")
  legend_items <- character(0)
  legend_cols <- character(0)

  if (!is.null(roc_1)) {
    lines(roc_1$FP, roc_1$TP, col = "#1B9E77", lwd = 2)
    legend_items <- c(legend_items, paste0("1-year AUC = ", sprintf("%.3f", auc_1)))
    legend_cols <- c(legend_cols, "#1B9E77")
  }
  if (!is.null(roc_2)) {
    lines(roc_2$FP, roc_2$TP, col = "#D95F02", lwd = 2)
    legend_items <- c(legend_items, paste0("2-year AUC = ", sprintf("%.3f", auc_2)))
    legend_cols <- c(legend_cols, "#D95F02")
  }
  if (!is.null(roc_3)) {
    lines(roc_3$FP, roc_3$TP, col = "#7570B3", lwd = 2)
    legend_items <- c(legend_items, paste0("3-year AUC = ", sprintf("%.3f", auc_3)))
    legend_cols <- c(legend_cols, "#7570B3")
  }
  if (length(legend_items) > 0) {
    legend("bottomright", legend = legend_items, col = legend_cols, lwd = 2, bty = "n")
  }
})

save_session_info("16_logs/sessionInfo_16_rebuild_cross_platform_lasso_signature_tcga.txt")
