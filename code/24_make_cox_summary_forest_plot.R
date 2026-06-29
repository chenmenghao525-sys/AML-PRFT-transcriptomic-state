#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(survival)
  library(grid)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
figure_dir <- file.path(project_root, "13_figures")
main_figure_dir <- file.path(figure_dir, "main_figures")
table_dir <- file.path(project_root, "14_tables")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(main_figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

ext_file <- file.path(project_root, "14_tables", "external_validation_main_results.csv")
tcga_perf_file <- file.path(project_root, "07_signature", "tcga_cross_platform_lasso_model_performance.csv")
tcga_summary_file <- file.path(project_root, "14_tables", "tcga_cross_platform_lasso_signature_summary.csv")
tcga_risk_file <- file.path(project_root, "07_signature", "tcga_cross_platform_risk_score_by_sample.csv")
tcga_cox_reference_file <- file.path(project_root, "07_signature", "tcga_cross_platform_univariate_multivariate_cox_risk.csv")

required_files <- c(ext_file, tcga_perf_file, tcga_summary_file, tcga_risk_file, tcga_cox_reference_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required result files: ", paste(missing_files, collapse = "; "))
}

ext_dt <- fread(ext_file)
tcga_perf_dt <- fread(tcga_perf_file)
tcga_summary_dt <- fread(tcga_summary_file)
tcga_risk_dt <- fread(tcga_risk_file)
tcga_ref_dt <- fread(tcga_cox_reference_file)

fmt_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return(formatC(p, format = "e", digits = 2))
  formatC(p, format = "f", digits = 3)
}

fmt_hr <- function(x) {
  out <- rep(NA_character_, length(x))
  ok <- !is.na(x)
  out[ok] <- formatC(x[ok], format = "f", digits = 2)
  out
}

cohort_map <- c(
  "GSE37642_GPL570" = "GSE37642 GPL570",
  "GSE12417_GPL570" = "GSE12417 GPL570",
  "combined_GPL570" = "Combined GPL570 validation cohort"
)

extract_risk_row <- function(sum_obj, coef_name = "risk_score") {
  idx <- which(rownames(sum_obj$coefficients) == coef_name)[1]
  if (is.na(idx)) {
    stop("Could not find coefficient row for ", coef_name)
  }
  list(
    HR = as.numeric(sum_obj$coefficients[idx, "exp(coef)"]),
    lower = as.numeric(sum_obj$conf.int[idx, "lower .95"]),
    upper = as.numeric(sum_obj$conf.int[idx, "upper .95"]),
    p = as.numeric(sum_obj$coefficients[idx, "Pr(>|z|)"])
  )
}

calc_tcga_rows <- function(tcga_risk_dt, tcga_ref_dt) {
  dt <- copy(tcga_risk_dt)
  keep_cols <- intersect(c("sample_id", "patient_id", "OS_time", "OS_status", "age", "sex", "risk_score"), names(dt))
  dt <- dt[, ..keep_cols]
  dt <- dt[!is.na(OS_time) & !is.na(OS_status) & !is.na(risk_score) & OS_time > 0]
  if (!nrow(dt)) {
    stop("No TCGA samples available for Cox regression after filtering.")
  }

  surv_obj <- Surv(dt$OS_time, dt$OS_status)
  uni_fit <- coxph(surv_obj ~ risk_score, data = dt)
  uni_sum <- summary(uni_fit)
  uni_row <- extract_risk_row(uni_sum)

  sex_ok <- "sex" %in% names(dt) &&
    sum(!is.na(dt$sex)) > 0 &&
    length(unique(na.omit(as.character(dt$sex)))) >= 2
  age_ok <- "age" %in% names(dt) && sum(!is.na(dt$age)) / nrow(dt) >= 0.70

  if (age_ok && sex_ok) {
    multi_dt <- dt[complete.cases(dt[, .(OS_time, OS_status, risk_score, age, sex)])]
    multi_dt[, sex := factor(sex)]
    multiv_formula <- Surv(OS_time, OS_status) ~ risk_score + age + sex
    covariates <- "risk score + age + sex"
    note <- "Multivariate Cox adjusted for age and sex."
  } else if (age_ok) {
    multi_dt <- dt[complete.cases(dt[, .(OS_time, OS_status, risk_score, age)])]
    multiv_formula <- Surv(OS_time, OS_status) ~ risk_score + age
    covariates <- "risk score + age"
    note <- "Sex was unavailable or not modelable; multivariate Cox adjusted for age only."
  } else {
    stop("Age availability is insufficient for the requested multivariate Cox model.")
  }

  if (!nrow(multi_dt)) {
    stop("No TCGA samples available for multivariate Cox regression after filtering.")
  }

  multi_fit <- coxph(multiv_formula, data = multi_dt)
  multi_sum <- summary(multi_fit)
  multi_row <- extract_risk_row(multi_sum)

  ref_uni_hr <- as.numeric(tcga_ref_dt[model == "univariate", HR_risk_score][1])
  ref_uni_p <- as.numeric(tcga_ref_dt[model == "univariate", P.Value_risk_score][1])
  ref_multi_hr <- as.numeric(tcga_ref_dt[model == "multivariate", HR_risk_score][1])
  ref_multi_p <- as.numeric(tcga_ref_dt[model == "multivariate", P.Value_risk_score][1])

  if (any(is.na(c(ref_uni_hr, ref_uni_p, ref_multi_hr, ref_multi_p)))) {
    stop("Reference TCGA cross-platform Cox file is incomplete.")
  }

  if (abs(uni_row$HR - ref_uni_hr) > 1e-6 || abs(uni_row$p - ref_uni_p) > 1e-12) {
    stop("Recomputed TCGA univariate Cox result does not match the fixed six-gene reference file closely enough.")
  }
  if (abs(multi_row$HR - ref_multi_hr) > 1e-6 || abs(multi_row$p - ref_multi_p) > 1e-12) {
    stop("Recomputed TCGA multivariate Cox result does not match the fixed six-gene reference file closely enough.")
  }

  rbind(
    data.table(
      cohort = "TCGA training cohort",
      analysis_type = "Univariate Cox",
      sample_count = nrow(dt),
      event_count = sum(dt$OS_status == 1, na.rm = TRUE),
      HR = uni_row$HR,
      lower_95_CI = uni_row$lower,
      upper_95_CI = uni_row$upper,
      P_value = uni_row$p,
      covariates = "risk score",
      source_model = "cross-platform six-gene model",
      source_file = basename(tcga_risk_file),
      note = "TCGA Cox was recomputed from the fixed six-gene risk score by sample."
    ),
    data.table(
      cohort = "TCGA training cohort",
      analysis_type = "Multivariate Cox",
      sample_count = nrow(multi_dt),
      event_count = sum(multi_dt$OS_status == 1, na.rm = TRUE),
      HR = multi_row$HR,
      lower_95_CI = multi_row$lower,
      upper_95_CI = multi_row$upper,
      P_value = multi_row$p,
      covariates = covariates,
      source_model = "cross-platform six-gene model",
      source_file = basename(tcga_risk_file),
      note = note
    )
  )
}

tcga_rows <- calc_tcga_rows(tcga_risk_dt, tcga_ref_dt)

ext_rows <- rbindlist(lapply(seq_len(nrow(ext_dt)), function(i) {
  row <- ext_dt[i]
  uni <- data.table(
    cohort = cohort_map[[row$dataset_id]],
    analysis_type = "Univariate Cox",
    sample_count = as.integer(row$sample_count),
    event_count = as.integer(row$event_count),
    HR = as.numeric(row$univariate_cox_HR),
    lower_95_CI = as.numeric(row$univariate_cox_lower95),
    upper_95_CI = as.numeric(row$univariate_cox_upper95),
    P_value = as.numeric(row$univariate_cox_P.Value),
    covariates = "risk score",
    source_model = "cross-platform six-gene model",
    source_file = basename(ext_file),
    note = "External validation used the fixed cross-platform six-gene model."
  )
  mv <- data.table(
    cohort = cohort_map[[row$dataset_id]],
    analysis_type = "Multivariate Cox",
    sample_count = as.integer(row$sample_count),
    event_count = as.integer(row$event_count),
    HR = as.numeric(row$multivariate_cox_HR),
    lower_95_CI = as.numeric(row$multivariate_cox_lower95),
    upper_95_CI = as.numeric(row$multivariate_cox_upper95),
    P_value = as.numeric(row$multivariate_cox_P.Value),
    covariates = gsub("_", " ", as.character(row$multivariate_model_used), fixed = TRUE),
    source_model = "cross-platform six-gene model",
    source_file = basename(ext_file),
    note = "External validation used the fixed cross-platform six-gene model."
  )
  rbind(uni, mv)
}), fill = TRUE)

plot_dt <- rbind(tcga_rows, ext_rows, fill = TRUE)
plot_dt <- plot_dt[!is.na(cohort) & !is.na(HR)]
plot_dt[, CI_available := !is.na(lower_95_CI) & !is.na(upper_95_CI)]
plot_dt[, cohort := factor(cohort, levels = c(
  "TCGA training cohort",
  "GSE37642 GPL570",
  "GSE12417 GPL570",
  "Combined GPL570 validation cohort"
))]
plot_dt[, analysis_type := factor(analysis_type, levels = c("Univariate Cox", "Multivariate Cox"))]
plot_dt <- plot_dt[order(cohort, analysis_type)]

out_dt <- copy(plot_dt)[, .(
  cohort,
  analysis_type,
  sample_count,
  event_count,
  HR,
  lower_95_CI,
  upper_95_CI,
  P_value,
  covariates,
  source_model,
  source_file,
  note
)]
fwrite(out_dt, file.path(table_dir, "Figure4F_cox_summary_table.csv"))

plot_dt[, hr_label := fmt_hr(HR)]
plot_dt[, ci_label := paste0("(", fmt_hr(lower_95_CI), ", ", fmt_hr(upper_95_CI), ")")]
plot_dt[, p_label := vapply(P_value, fmt_p, character(1))]
plot_dt[, sample_event_label := paste0(sample_count, " / ", event_count)]
plot_dt[, row_id := .I]

xmax_ci <- max(plot_dt$upper_95_CI, na.rm = TRUE)
x_text1 <- xmax_ci * 1.18
x_text2 <- xmax_ci * 1.52
x_text3 <- xmax_ci * 1.98
x_text4 <- xmax_ci * 2.38

p <- ggplot(plot_dt, aes(x = HR, y = row_id, color = analysis_type)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#808B96", linewidth = 0.5) +
  geom_errorbar(
    aes(y = row_id, xmin = lower_95_CI, xmax = upper_95_CI),
    orientation = "y",
    width = 0.16,
    linewidth = 0.7,
    na.rm = TRUE
  ) +
  geom_point(size = 2.8) +
  scale_color_manual(values = c("Univariate Cox" = "#315B7D", "Multivariate Cox" = "#B6423C")) +
  scale_y_continuous(
    breaks = plot_dt$row_id,
    labels = paste0(as.character(plot_dt$cohort), " | ", as.character(plot_dt$analysis_type)),
    trans = "reverse"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.45))) +
  annotate("text", x = x_text1, y = 0.45, label = "HR", fontface = "bold", hjust = 0, size = 3.7, color = "#1F2D3A") +
  annotate("text", x = x_text2, y = 0.45, label = "95% CI", fontface = "bold", hjust = 0, size = 3.7, color = "#1F2D3A") +
  annotate("text", x = x_text3, y = 0.45, label = "P value", fontface = "bold", hjust = 0, size = 3.7, color = "#1F2D3A") +
  annotate("text", x = x_text4, y = 0.45, label = "Samples / events", fontface = "bold", hjust = 0, size = 3.7, color = "#1F2D3A") +
  geom_text(aes(x = x_text1, label = hr_label), hjust = 0, size = 3.5, color = "#1F2D3A", show.legend = FALSE) +
  geom_text(aes(x = x_text2, label = ci_label), hjust = 0, size = 3.5, color = "#1F2D3A", show.legend = FALSE) +
  geom_text(aes(x = x_text3, label = p_label), hjust = 0, size = 3.5, color = "#1F2D3A", show.legend = FALSE) +
  geom_text(aes(x = x_text4, label = sample_event_label), hjust = 0, size = 3.5, color = "#1F2D3A", show.legend = FALSE) +
  labs(
    title = "Cox summary of the PRFT-related risk score across training and validation cohorts",
    subtitle = "Training and external validation cohorts are shown without re-estimating the fixed six-gene risk score",
    x = "Hazard ratio",
    y = NULL,
    color = "Analysis"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", color = "#1F2D3A"),
    plot.subtitle = element_text(hjust = 0.5, color = "#4E5D6C"),
    axis.text.y = element_text(size = 9.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(16, 180, 16, 16)
  ) +
  coord_cartesian(clip = "off")

out_base <- file.path(figure_dir, "Figure4F_cox_summary_forest_plot")
ggsave(paste0(out_base, ".pdf"), plot = p, width = 13.6, height = 5.8, units = "in")
ggsave(paste0(out_base, ".png"), plot = p, width = 13.6, height = 5.8, units = "in", dpi = 300)

file.copy(paste0(out_base, ".pdf"), file.path(main_figure_dir, "Figure4F_cox_summary_forest_plot.pdf"), overwrite = TRUE)
file.copy(paste0(out_base, ".png"), file.path(main_figure_dir, "Figure4F_cox_summary_forest_plot.png"), overwrite = TRUE)

message("Figure 4F Cox summary forest plot generated successfully.")
message("Files used: ", paste(basename(required_files), collapse = "; "))
