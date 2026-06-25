#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
source("15_scripts/plot_label_utils.R")

dir.create("04_prft_score", recursive = TRUE, showWarnings = FALSE)
dir.create("13_figures", recursive = TRUE, showWarnings = FALSE)
dir.create("14_tables", recursive = TRUE, showWarnings = FALSE)
dir.create("16_logs", recursive = TRUE, showWarnings = FALSE)

save_session_info <- function(path) {
  writeLines(capture.output(sessionInfo()), con = path)
}

zscore_safe <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) {
    return(rep(0, length(x)))
  }
  as.numeric(scale(x))
}

scores <- readRDS("04_prft_score/tcga_ssgsea_scores.rds")
obj <- readRDS("02_processed_data/tcga_expr_clin_matched.rds")

required_sets <- c("Proteostasis_core", "Ferroptosis_tolerance_set")
missing_sets <- setdiff(required_sets, rownames(scores))
if (length(missing_sets) > 0) {
  stop("Missing required ssGSEA rows: ", paste(missing_sets, collapse = ", "))
}

proteostasis_score <- as.numeric(scores["Proteostasis_core", ])
ferroptosis_tolerance_score <- as.numeric(scores["Ferroptosis_tolerance_set", ])
sample_ids <- colnames(scores)

z_proteostasis <- zscore_safe(proteostasis_score)
z_ferroptosis_tolerance <- zscore_safe(ferroptosis_tolerance_score)
prft_score <- (z_proteostasis + z_ferroptosis_tolerance) / 2
cutoff <- stats::median(prft_score, na.rm = TRUE)
prft_group <- ifelse(prft_score >= cutoff, "PRFT_high", "PRFT_low")

prft_df <- data.frame(
  sample_id = sample_ids,
  Proteostasis_core_score = proteostasis_score,
  Ferroptosis_tolerance_set_score = ferroptosis_tolerance_score,
  z_Proteostasis_core = z_proteostasis,
  z_Ferroptosis_tolerance_set = z_ferroptosis_tolerance,
  PRFT_score = prft_score,
  PRFT_group = prft_group,
  stringsAsFactors = FALSE
)

clin <- obj$clin
merge_cols <- intersect(c("sample_id", "patient_id", "OS_time", "OS_status", "age", "sex", "FAB", "WBC"), colnames(clin))
prft_df <- merge(prft_df, clin[, merge_cols, drop = FALSE], by = "sample_id", all.x = TRUE, sort = FALSE)
prft_df <- prft_df[match(sample_ids, prft_df$sample_id), , drop = FALSE]

fwrite(prft_df, "04_prft_score/tcga_prft_score.csv")
fwrite(prft_df[, c("sample_id", "PRFT_group")], "04_prft_score/tcga_prft_group.csv")
saveRDS(prft_df, "04_prft_score/tcga_prft_score.rds")

group_summary <- data.table(
  metric = c("total_samples", "PRFT_high", "PRFT_low", "median_cutoff"),
  value = c(
    nrow(prft_df),
    sum(prft_df$PRFT_group == "PRFT_high", na.rm = TRUE),
    sum(prft_df$PRFT_group == "PRFT_low", na.rm = TRUE),
    cutoff
  )
)
fwrite(group_summary, "14_tables/tcga_prft_group_summary.csv")

plot_df <- prft_df

p1 <- ggplot(plot_df, aes(x = PRFT_score)) +
  geom_histogram(bins = 30, fill = "#4C78A8", color = "white") +
  geom_vline(xintercept = cutoff, linetype = "dashed", color = "#D62728", linewidth = 0.8) +
  labs(
    title = "PRFT score distribution",
    x = "PRFT score",
    y = "Sample count"
  ) +
  theme_bw(base_size = 11)

ggsave("13_figures/Figure1_PRFT_score_distribution.pdf", p1, width = 6, height = 4.5)
ggsave("13_figures/Figure1_PRFT_score_distribution.png", p1, width = 6, height = 4.5, dpi = 300)

cor_val <- suppressWarnings(stats::cor(
  plot_df$Proteostasis_core_score,
  plot_df$Ferroptosis_tolerance_set_score,
  method = "spearman",
  use = "complete.obs"
))

p2 <- ggplot(plot_df, aes(x = Proteostasis_core_score, y = Ferroptosis_tolerance_set_score)) +
  geom_point(color = "#4C78A8", alpha = 0.8, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "#D62728", linewidth = 0.8) +
  annotate(
    "text",
    x = Inf, y = Inf,
    hjust = 1.1, vjust = 1.5,
    label = paste0("Spearman rho = ", round(cor_val, 3))
  ) +
  labs(
    title = "PRFT subscore correlation",
    x = pretty_label("Proteostasis_core_score"),
    y = pretty_label("Ferroptosis_tolerance_set_score")
  ) +
  theme_bw(base_size = 11)

ggsave("13_figures/Figure1_PRFT_subscore_correlation.pdf", p2, width = 6, height = 4.5)
ggsave("13_figures/Figure1_PRFT_subscore_correlation.png", p2, width = 6, height = 4.5, dpi = 300)

save_session_info("16_logs/sessionInfo_09_calculate_prft_score_and_group.txt")
