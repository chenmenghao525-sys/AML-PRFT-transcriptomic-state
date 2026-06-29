#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
})
source("15_scripts/plot_label_utils.R")

dir.create("07_signature", recursive = TRUE, showWarnings = FALSE)
dir.create("13_figures", recursive = TRUE, showWarnings = FALSE)
dir.create("14_tables", recursive = TRUE, showWarnings = FALSE)
dir.create("16_logs", recursive = TRUE, showWarnings = FALSE)

save_session_info <- function(path) {
  writeLines(capture.output(sessionInfo()), con = path)
}

detect_mm_col <- function(dt) {
  mm_candidates <- c("MM", "MM_PRFT_module", "module_membership", "kME")
  mm_found <- intersect(mm_candidates, colnames(dt))
  if (length(mm_found) == 0) {
    stop("No MM-like column found in candidate pool. Checked: ", paste(mm_candidates, collapse = ", "))
  }
  mm_found[1]
}

run_univariate_cox <- function(expr_mat, clin_df, gene_list, candidate_set_label) {
  results <- vector("list", length(gene_list))

  surv_obj <- survival::Surv(clin_df$OS_time, clin_df$OS_status)

  for (i in seq_along(gene_list)) {
    gene <- gene_list[i]
    if (!gene %in% rownames(expr_mat)) {
      results[[i]] <- data.frame(
        gene_symbol = gene,
        HR = NA_real_,
        lower95 = NA_real_,
        upper95 = NA_real_,
        coef = NA_real_,
        P.Value = NA_real_,
        FDR = NA_real_,
        cox_direction = NA_character_,
        mean_expression = NA_real_,
        missing_rate = 1,
        candidate_set = candidate_set_label,
        stringsAsFactors = FALSE
      )
      next
    }

    gene_expr <- as.numeric(expr_mat[gene, clin_df$sample_id])
    missing_rate <- mean(is.na(gene_expr))
    mean_expression <- mean(gene_expr, na.rm = TRUE)

    fit_df <- data.frame(OS_time = clin_df$OS_time, OS_status = clin_df$OS_status, gene_expression = gene_expr)
    fit_df <- fit_df[stats::complete.cases(fit_df), , drop = FALSE]

    if (nrow(fit_df) < 20 || stats::sd(fit_df$gene_expression, na.rm = TRUE) == 0) {
      results[[i]] <- data.frame(
        gene_symbol = gene,
        HR = NA_real_,
        lower95 = NA_real_,
        upper95 = NA_real_,
        coef = NA_real_,
        P.Value = NA_real_,
        FDR = NA_real_,
        cox_direction = NA_character_,
        mean_expression = mean_expression,
        missing_rate = missing_rate,
        candidate_set = candidate_set_label,
        stringsAsFactors = FALSE
      )
      next
    }

    fit <- tryCatch(
      survival::coxph(survival::Surv(OS_time, OS_status) ~ gene_expression, data = fit_df),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      results[[i]] <- data.frame(
        gene_symbol = gene,
        HR = NA_real_,
        lower95 = NA_real_,
        upper95 = NA_real_,
        coef = NA_real_,
        P.Value = NA_real_,
        FDR = NA_real_,
        cox_direction = NA_character_,
        mean_expression = mean_expression,
        missing_rate = missing_rate,
        candidate_set = candidate_set_label,
        stringsAsFactors = FALSE
      )
      next
    }

    fit_sum <- summary(fit)
    coef_val <- unname(fit_sum$coefficients[1, "coef"])
    hr_val <- unname(fit_sum$coefficients[1, "exp(coef)"])
    p_val <- unname(fit_sum$coefficients[1, "Pr(>|z|)"])
    lower95 <- unname(fit_sum$conf.int[1, "lower .95"])
    upper95 <- unname(fit_sum$conf.int[1, "upper .95"])

    results[[i]] <- data.frame(
      gene_symbol = gene,
      HR = hr_val,
      lower95 = lower95,
      upper95 = upper95,
      coef = coef_val,
      P.Value = p_val,
      FDR = NA_real_,
      cox_direction = ifelse(is.na(hr_val), NA_character_, ifelse(hr_val > 1, "risk", "protective")),
      mean_expression = mean_expression,
      missing_rate = missing_rate,
      candidate_set = candidate_set_label,
      stringsAsFactors = FALSE
    )
  }

  result_dt <- as.data.table(rbindlist(results, fill = TRUE))
  result_dt[, FDR := p.adjust(P.Value, method = "BH")]
  result_dt
}

obj <- readRDS("02_processed_data/tcga_expr_clin_matched.rds")
prft_df <- as.data.table(readRDS("04_prft_score/tcga_prft_score.rds"))
candidate_pool <- fread("06_wgcna/wgcna_deg_up_intersect_prft_module_genes.csv")
wgcna_gene_level <- fread("06_wgcna/wgcna_gene_level_statistics.csv")
deg_all <- fread("05_deg/deg_prft_high_vs_low_tcga_all.csv")

expr <- as.matrix(obj$expr)
storage.mode(expr) <- "numeric"

required_prft_cols <- c("sample_id", "OS_time", "OS_status")
if (!all(required_prft_cols %in% colnames(prft_df))) {
  stop("tcga_prft_score.rds must contain sample_id, OS_time, OS_status.")
}

candidate_pool_all_genes <- unique(candidate_pool$gene_symbol)
if (length(candidate_pool_all_genes) != 715) {
  warning("Initial candidate_pool_all gene count is ", length(candidate_pool_all_genes), ", expected 715.")
}
writeLines(candidate_pool_all_genes, "07_signature/candidate_pool_all_genes.txt")

mm_col <- detect_mm_col(candidate_pool)
candidate_pool[, MM_used := get(mm_col)]

strict_candidate <- copy(candidate_pool)[
  logFC > 0.5 &
    adj.P.Val < 0.05 &
    GS_PRFT > 0.30 &
    MM_used > 0.50
]
strict_rule_used <- "logFC > 0.5; adj.P.Val < 0.05; GS_PRFT > 0.30; MM > 0.50"

if (nrow(strict_candidate) < 20) {
  strict_candidate <- copy(candidate_pool)[
    logFC > 0.5 &
      adj.P.Val < 0.05 &
      GS_PRFT > 0.25 &
      MM_used > 0.40
  ]
  strict_rule_used <- "Relaxed strict filter: logFC > 0.5; adj.P.Val < 0.05; GS_PRFT > 0.25; MM > 0.40"
}

if (nrow(strict_candidate) < 20) {
  strict_candidate <- copy(candidate_pool)
  strict_rule_used <- "Strict filter yielded <20 genes; fallback to candidate_pool_all as main Cox input"
}

candidate_pool_strict_genes <- unique(strict_candidate$gene_symbol)
writeLines(candidate_pool_strict_genes, "07_signature/candidate_pool_strict_genes.txt")

cox_clin <- unique(prft_df[, .(sample_id, OS_time, OS_status)])
cox_clin <- cox_clin[!is.na(OS_time) & !is.na(OS_status)]
removed_nonpositive_ids <- cox_clin[OS_time <= 0, sample_id]
cox_clin <- cox_clin[OS_time > 0]

common_samples <- intersect(colnames(expr), cox_clin$sample_id)
expr_cox <- expr[, common_samples, drop = FALSE]
cox_clin <- cox_clin[match(common_samples, sample_id)]

if (!identical(colnames(expr_cox), cox_clin$sample_id)) {
  stop("Expression matrix and Cox clinical data sample IDs are not aligned.")
}

cox_all_dt <- run_univariate_cox(expr_cox, cox_clin, candidate_pool_all_genes, "candidate_pool_all")
cox_strict_dt <- run_univariate_cox(expr_cox, cox_clin, candidate_pool_strict_genes, "candidate_pool_strict")

fwrite(cox_all_dt, "07_signature/univariate_cox_candidate_pool_all.csv")
fwrite(cox_strict_dt, "07_signature/univariate_cox_candidate_pool_strict.csv")

candidate_pool_info <- unique(candidate_pool[, .(
  gene_symbol, module_color, GS_PRFT, GS_PRFT_P, MM = MM_used, MM_P, logFC, adj.P.Val
)])

cox_integrated <- merge(
  cox_all_dt,
  candidate_pool_info,
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)
fwrite(cox_integrated, "07_signature/univariate_cox_integrated_results.csv")

main_candidates <- copy(cox_integrated)[
  HR > 1 &
    P.Value < 0.05 &
    logFC > 0.5 &
    adj.P.Val < 0.05 &
    GS_PRFT > 0.30 &
    MM > 0.50
]
filtering_rule_used <- "HR > 1; P.Value < 0.05; logFC > 0.5; adj.P.Val < 0.05; GS_PRFT > 0.30; MM > 0.50"

if (nrow(main_candidates) < 10) {
  main_candidates <- copy(cox_integrated)[
    HR > 1 &
      P.Value < 0.05 &
      logFC > 0.5 &
      adj.P.Val < 0.05 &
      GS_PRFT > 0.25 &
      MM > 0.40
  ]
  filtering_rule_used <- "Relaxed main candidate rule: HR > 1; P.Value < 0.05; logFC > 0.5; adj.P.Val < 0.05; GS_PRFT > 0.25; MM > 0.40"
}

if (nrow(main_candidates) < 10) {
  main_candidates <- copy(cox_integrated)[HR > 1 & P.Value < 0.05]
  filtering_rule_used <- "Fallback main candidate rule: HR > 1 and P.Value < 0.05"
}

setorder(main_candidates, P.Value, -HR, -GS_PRFT, -MM, -logFC)
cox_preselected_genes <- unique(main_candidates$gene_symbol)

if (nrow(main_candidates) > 80) {
  main_candidates <- main_candidates[1:80]
  filtering_rule_used <- paste0(filtering_rule_used, "; trimmed to top 80 by P.Value, HR, GS_PRFT, MM, logFC")
}

fwrite(main_candidates, "07_signature/main_candidate_genes_for_lasso.csv")
writeLines(main_candidates$gene_symbol, "07_signature/main_candidate_genes_for_lasso.txt")
writeLines(cox_preselected_genes, "07_signature/cox_preselected_genes.txt")

cox_hr1_p005 <- cox_integrated[HR > 1 & P.Value < 0.05]
cox_hr1_fdr025 <- cox_integrated[HR > 1 & FDR < 0.25]
cox_hr1_fdr010 <- cox_integrated[HR > 1 & FDR < 0.10]

summary_dt <- data.table(
  samples_with_complete_survival = nrow(cox_clin),
  candidate_pool_all_genes = length(candidate_pool_all_genes),
  candidate_pool_strict_genes = length(candidate_pool_strict_genes),
  univariate_cox_all_success = sum(!is.na(cox_all_dt$P.Value)),
  univariate_cox_strict_success = sum(!is.na(cox_strict_dt$P.Value)),
  cox_P_less_0.05_HR_greater_1 = nrow(cox_hr1_p005),
  cox_FDR_less_0.25_HR_greater_1 = nrow(cox_hr1_fdr025),
  cox_FDR_less_0.10_HR_greater_1 = nrow(cox_hr1_fdr010),
  main_candidate_genes_count = length(unique(main_candidates$gene_symbol)),
  lasso_input_genes_count = length(unique(main_candidates$gene_symbol)),
  filtering_rule_used = filtering_rule_used
)
fwrite(summary_dt, "14_tables/tcga_univariate_cox_summary.csv")

forest_dt <- copy(cox_integrated)[HR > 1 & P.Value < 0.05]
setorder(forest_dt, P.Value)
forest_dt <- forest_dt[1:min(20, .N)]
if (nrow(forest_dt) > 0) {
  forest_dt[, gene_symbol := factor(gene_symbol, levels = rev(gene_symbol))]
  p_forest <- ggplot(forest_dt, aes(x = gene_symbol, y = HR, ymin = lower95, ymax = upper95)) +
    geom_pointrange(color = "#1F77B4") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
    coord_flip() +
    labs(x = NULL, y = "Hazard ratio (95% CI)", title = "Top 20 univariate Cox candidate genes") +
    theme_bw(base_size = 11)
  ggsave("13_figures/Figure4_univariate_cox_forest_top20.pdf", p_forest, width = 7, height = 6.5)
  ggsave("13_figures/Figure4_univariate_cox_forest_top20.png", p_forest, width = 7, height = 6.5, dpi = 300)
}

filter_counts_dt <- data.table(
  step = c(
    "candidate_pool_all",
    "candidate_pool_strict",
    "HR>1 & P<0.05",
    "HR>1 & FDR<0.25",
    "main_candidate_genes",
    "lasso_input_genes"
  ),
  count = c(
    length(candidate_pool_all_genes),
    length(candidate_pool_strict_genes),
    nrow(cox_hr1_p005),
    nrow(cox_hr1_fdr025),
    length(unique(main_candidates$gene_symbol)),
    length(unique(main_candidates$gene_symbol))
  )
)
filter_counts_dt[, step_label := pretty_label(step)]

p_counts <- ggplot(filter_counts_dt, aes(x = factor(step_label, levels = step_label), y = count)) +
  geom_col(fill = "#4C78A8") +
  geom_text(aes(label = count), vjust = -0.3, size = 3.2) +
  labs(x = NULL, y = "Gene count", title = "Candidate gene filtering counts") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("13_figures/Figure4_candidate_filtering_counts.pdf", p_counts, width = 7.5, height = 5)
ggsave("13_figures/Figure4_candidate_filtering_counts.png", p_counts, width = 7.5, height = 5, dpi = 300)

save_session_info("16_logs/sessionInfo_12_candidate_gene_selection_and_univariate_cox.txt")
