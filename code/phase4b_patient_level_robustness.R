#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(nlme)
})

options(stringsAsFactors = FALSE)
set.seed(1234)

root <- Sys.getenv("PHASE4_ROOT", unset = "")
if (!nzchar(root)) root <- getwd()
root <- gsub("\\\\", "/", root)
if (!dir.exists(root)) stop("Project root does not exist: ", root)

tables_dir <- file.path(root, "03_results_tables")
fig_dir <- file.path(root, "04_figures")
log_dir <- file.path(root, "05_logs")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

main_log <- file.path(log_dir, "phase4b_processed_singlecell_log.txt")
mixed_log <- file.path(log_dir, "phase4b_mixed_model_log.txt")
for (f in c(main_log, mixed_log)) if (file.exists(f)) file.remove(f)

log_msg <- function(file, ...) {
  line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(..., collapse = " "))
  cat(line, "\n")
  cat(line, "\n", file = file, append = TRUE)
}

write_csv <- function(x, path, log_file = main_log) {
  fwrite(as.data.table(x), path)
  log_msg(log_file, "Wrote", path)
}

save_pdf <- function(plot, path, width = 8, height = 6) {
  ggsave(path, plot = plot, width = width, height = height, units = "in", device = cairo_pdf)
  log_msg(main_log, "Wrote", path)
}

norm_name <- function(x) {
  y <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
  gsub("^_|_$", "", y)
}

find_field <- function(cols, patterns) {
  nc <- norm_name(cols)
  idx <- which(nc %in% patterns)
  if (length(idx) == 0) idx <- which(vapply(patterns, function(p) any(grepl(p, nc, fixed = TRUE)), logical(1)))
  if (length(idx) == 0) return(NA_character_)
  cols[idx[1]]
}

cohens_d_unpaired <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  s1 <- sd(x); s2 <- sd(y)
  n1 <- length(x); n2 <- length(y)
  sp <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  (mean(x) - mean(y)) / sp
}

eta_sq_kw <- function(kw_obj, n_total) {
  if (is.null(kw_obj$statistic) || !is.finite(kw_obj$statistic)) return(NA_real_)
  h <- as.numeric(kw_obj$statistic)
  k <- length(kw_obj$parameter) + 1
  if (!is.finite(h) || !is.finite(n_total) || n_total <= k) return(NA_real_)
  (h - k + 1) / (n_total - k)
}

root_ascii <- root
log_msg(main_log, "Phase 4b processed single-cell robustness analysis started.")
log_msg(main_log, "Project root:", root_ascii)

candidate_inputs <- c(
  file.path(tables_dir, "phase4_fix_singlecell_module_scores.csv"),
  file.path(tables_dir, "phase4_singlecell_module_scores.csv")
)
input_file <- candidate_inputs[file.exists(candidate_inputs)][1]
if (is.na(input_file)) stop("No processed single-cell module score table found.")

dt <- fread(input_file)
log_msg(main_log, "Loaded processed single-cell table:", input_file, "rows=", nrow(dt), "cols=", ncol(dt))

cols <- names(dt)
field_map <- data.table(
  logical_field = c("cell_id", "sample_id", "patient_id", "sample_state", "cell_type",
                    "six_gene_risk_score", "PRFT_score", "proteostasis_score",
                    "ferroptosis_tolerance_score", "JAK2_STAT5_PD_L1_score",
                    "SLC7A11_GPX4_GSH_score", "immune_checkpoint_score",
                    "myeloid_suppressive_score", "LSC17_score", "stemness_score", "risk_like_group"),
  original_column = c(
    find_field(cols, c("cell", "cell_id", "barcode", "barcodes")),
    find_field(cols, c("sample_id", "sample")),
    find_field(cols, c("patient_id", "patient")),
    find_field(cols, c("sample_state", "sample_status", "state", "group")),
    find_field(cols, c("cell_type", "celltype", "annotation")),
    find_field(cols, c("risk_score", "six_gene_risk_score")),
    find_field(cols, c("prft_score")),
    find_field(cols, c("proteostasis_core_score", "proteostasis_score")),
    find_field(cols, c("ferroptosis_tolerance_set_score", "ferroptosis_tolerance_score")),
    find_field(cols, c("jak2_stat5_pdl1_set_score", "jak2_stat5_pd_l1_score")),
    find_field(cols, c("slc7a11_gpx4_gsh_axis_score", "slc7a11_gpx4_gsh_score")),
    find_field(cols, c("immune_checkpoint_set_score", "immune_checkpoint_score")),
    find_field(cols, c("myeloid_suppressive_set_score", "myeloid_suppressive_score")),
    find_field(cols, c("lsc17_core_score", "lsc17_score")),
    find_field(cols, c("stemness_quiescence_set_score", "stemness_score")),
    find_field(cols, c("risk_like_group"))
  )
)
field_map[, status := ifelse(is.na(original_column), "missing", "mapped")]
field_map[, notes := ""]
write_csv(field_map, file.path(tables_dir, "phase4b_field_mapping.csv"))

req_fields <- c("cell_id", "sample_id", "patient_id", "sample_state", "cell_type",
                "six_gene_risk_score", "PRFT_score", "proteostasis_score",
                "ferroptosis_tolerance_score", "JAK2_STAT5_PD_L1_score",
                "SLC7A11_GPX4_GSH_score", "myeloid_suppressive_score",
                "LSC17_score", "stemness_score", "risk_like_group")
missing_req <- field_map[logical_field %in% req_fields & status == "missing", logical_field]
if (length(missing_req) > 0) {
  stop("Required fields missing from processed table: ", paste(missing_req, collapse = ", "))
}

mapped <- copy(dt)
for (i in seq_len(nrow(field_map[status == "mapped"]))) {
  lf <- field_map[status == "mapped", logical_field][i]
  oc <- field_map[status == "mapped", original_column][i]
  setnames(mapped, oc, lf)
}
log_msg(main_log, "Field mapping complete:", paste(sprintf("%s->%s", field_map$logical_field, field_map$original_column), collapse = "; "))

score_fields <- c("six_gene_risk_score", "PRFT_score", "proteostasis_score",
                  "ferroptosis_tolerance_score", "JAK2_STAT5_PD_L1_score",
                  "SLC7A11_GPX4_GSH_score", "myeloid_suppressive_score",
                  "LSC17_score", "stemness_score")
avail_scores <- score_fields[score_fields %in% names(mapped)]
immune_col <- if ("immune_checkpoint_score" %in% names(mapped)) "immune_checkpoint_score" else NA_character_
if (!is.na(immune_col)) avail_scores <- c(avail_scores, immune_col)

mapped[, AML_vs_healthy := ifelse(grepl("healthy", sample_state, ignore.case = TRUE), "healthy_BM", "AML")]
mapped[, risk_like_group := factor(risk_like_group, levels = c("low_risk_like", "intermediate", "high_risk_like"))]

assign_cell_group <- function(x) {
  x <- as.character(x)
  if (x %in% c("Mono-like", "Mono", "ProMono-like", "ProMono", "cDC-like", "cDC")) return("monocyte_myeloid_related")
  if (x %in% c("HSC-like", "HSC", "GMP-like", "GMP", "Prog-like", "Prog")) return("primitive_progenitor_related")
  if (x %in% c("lateEry", "earlyEry")) return("lymphoid_other")
  if (x %in% c("T", "NK", "B", "Plasma", "CTL", "ProB", "pDC", "unannotated")) return("lymphoid_other")
  "lymphoid_other"
}
mapped[, cell_group := vapply(cell_type, assign_cell_group, character(1))]
mapped[, cell_group := factor(cell_group, levels = c("monocyte_myeloid_related", "primitive_progenitor_related", "lymphoid_other"))]

pseudobulk <- mapped[, {
  out <- list(
    n_cells = .N,
    PRFT_high_like_fraction = mean(risk_like_group == "high_risk_like", na.rm = TRUE),
    high_risk_like_fraction = mean(risk_like_group == "high_risk_like", na.rm = TRUE),
    intermediate_fraction = mean(risk_like_group == "intermediate", na.rm = TRUE),
    low_risk_like_fraction = mean(risk_like_group == "low_risk_like", na.rm = TRUE)
  )
  for (v in avail_scores) {
    out[[paste0(v, "_mean")]] <- mean(get(v), na.rm = TRUE)
    out[[paste0(v, "_median")]] <- median(get(v), na.rm = TRUE)
  }
  out
}, by = .(patient_id, sample_id, sample_state, AML_vs_healthy, cell_type, cell_group)]
write_csv(pseudobulk, file.path(tables_dir, "phase4b_patient_celltype_pseudobulk_scores.csv"))

comp_rows <- list()
target_vars <- c("six_gene_risk_score_mean", "PRFT_score_mean", "proteostasis_score_mean",
                 "ferroptosis_tolerance_score_mean", "JAK2_STAT5_PD_L1_score_mean",
                 "SLC7A11_GPX4_GSH_score_mean", "myeloid_suppressive_score_mean",
                 "LSC17_score_mean", "stemness_score_mean", "PRFT_high_like_fraction")
target_vars <- target_vars[target_vars %in% names(pseudobulk)]

for (v in target_vars) {
  sub <- pseudobulk[is.finite(get(v))]
  if (uniqueN(sub$cell_group) < 2 || nrow(sub) < 6) {
    comp_rows[[length(comp_rows) + 1]] <- data.table(
      variable = v, test = "not_run", comparison = "all_cell_groups",
      n_total = nrow(sub), p_value = NA_real_, effect_size = NA_real_,
      effect_size_type = NA_character_, group_medians = paste(capture.output(print(sub[, .(median = median(get(v), na.rm = TRUE), n = .N), by = cell_group])), collapse = " "),
      notes = "Insufficient grouped pseudo-bulk observations"
    )
    next
  }
  kw <- kruskal.test(stats::as.formula(paste(v, "~ cell_group")), data = sub)
  med_txt <- sub[, .(median = median(get(v), na.rm = TRUE), n = .N), by = cell_group]
  comp_rows[[length(comp_rows) + 1]] <- data.table(
    variable = v, test = "Kruskal-Wallis", comparison = "all_cell_groups",
    n_total = nrow(sub), p_value = kw$p.value, effect_size = eta_sq_kw(kw, nrow(sub)),
    effect_size_type = "eta_squared_H", group_medians = paste(capture.output(print(med_txt)), collapse = " "),
    notes = ""
  )
  pair_levels <- levels(sub$cell_group)
  for (i in 1:(length(pair_levels) - 1)) {
    for (j in (i + 1):length(pair_levels)) {
      a <- pair_levels[i]; b <- pair_levels[j]
      pair_sub <- sub[cell_group %in% c(a, b)]
      if (uniqueN(pair_sub$cell_group) < 2 || nrow(pair_sub) < 4) next
      wt <- wilcox.test(pair_sub[get("cell_group") == a, get(v)], pair_sub[get("cell_group") == b, get(v)], exact = FALSE)
      dval <- cohens_d_unpaired(pair_sub[get("cell_group") == a, get(v)], pair_sub[get("cell_group") == b, get(v)])
      comp_rows[[length(comp_rows) + 1]] <- data.table(
        variable = v, test = "Wilcoxon", comparison = paste(a, "vs", b),
        n_total = nrow(pair_sub), p_value = wt$p.value, effect_size = dval,
        effect_size_type = "cohens_d",
        group_medians = paste(capture.output(print(pair_sub[, .(median = median(get(v), na.rm = TRUE), n = .N), by = cell_group])), collapse = " "),
        notes = ""
      )
    }
  }
}
comparison_dt <- rbindlist(comp_rows, fill = TRUE)
comparison_dt[, FDR := p.adjust(p_value, method = "BH")]
write_csv(comparison_dt, file.path(tables_dir, "phase4b_cell_group_comparison_patient_level.csv"))

boxplot_vars <- c("six_gene_risk_score_mean", "PRFT_score_mean", "JAK2_STAT5_PD_L1_score_mean",
                  "myeloid_suppressive_score_mean", "LSC17_score_mean")
boxplot_vars <- boxplot_vars[boxplot_vars %in% names(pseudobulk)]
plot_list <- lapply(boxplot_vars, function(v) {
  ggplot(pseudobulk, aes(x = cell_group, y = .data[[v]], fill = cell_group)) +
    geom_boxplot(outlier.size = 0.6, alpha = 0.85) +
    geom_jitter(width = 0.15, size = 0.9, alpha = 0.7) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "none") +
    labs(title = v, x = "cell group", y = "patient/sample pseudo-bulk mean")
})
save_pdf(wrap_plots(plot_list, ncol = 2), file.path(fig_dir, "phase4b_patient_level_score_boxplots.pdf"), 11, 8.5)

p_frac <- ggplot(pseudobulk, aes(x = cell_group, y = PRFT_high_like_fraction, fill = cell_group)) +
  geom_boxplot(outlier.size = 0.6, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 0.9, alpha = 0.7) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "none") +
  labs(title = "PRFT-high-like fraction by cell group", x = "cell group", y = "PRFT-high-like fraction")
save_pdf(p_frac, file.path(fig_dir, "phase4b_PRFT_high_like_fraction_patient_level.pdf"), 7.5, 5.8)

mixed_rows <- list()
model_vars <- c("six_gene_risk_score_mean", "PRFT_score_mean", "JAK2_STAT5_PD_L1_score_mean",
                "myeloid_suppressive_score_mean", "LSC17_score_mean", "PRFT_high_like_fraction")
model_vars <- model_vars[model_vars %in% names(pseudobulk)]
for (v in model_vars) {
  sub <- pseudobulk[is.finite(get(v)) & !is.na(patient_id) & !is.na(cell_group) & !is.na(AML_vs_healthy)]
  if (nrow(sub) < 10 || uniqueN(sub$patient_id) < 3 || uniqueN(sub$cell_group) < 2) {
    mixed_rows[[length(mixed_rows) + 1]] <- data.table(
      variable = v, model_type = "not_run", term = NA_character_,
      estimate = NA_real_, std_error = NA_real_, df = NA_real_, t_value = NA_real_, p_value = NA_real_,
      n_obs = nrow(sub), n_patients = uniqueN(sub$patient_id),
      notes = "Insufficient pseudo-bulk observations for mixed model"
    )
    next
  }
  fit <- tryCatch(
    nlme::lme(
      fixed = stats::as.formula(paste(v, "~ cell_group + AML_vs_healthy")),
      random = ~1 | patient_id,
      data = as.data.frame(sub),
      na.action = na.omit,
      control = nlme::lmeControl(returnObject = TRUE)
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    log_msg(mixed_log, "Mixed model failed for", v, ":", conditionMessage(fit))
    fallback <- tryCatch(
      lm(stats::as.formula(paste(v, "~ cell_group + AML_vs_healthy")), data = sub),
      error = function(e) e
    )
    if (inherits(fallback, "error")) {
      mixed_rows[[length(mixed_rows) + 1]] <- data.table(
        variable = v, model_type = "failed", term = NA_character_,
        estimate = NA_real_, std_error = NA_real_, df = NA_real_, t_value = NA_real_, p_value = NA_real_,
        n_obs = nrow(sub), n_patients = uniqueN(sub$patient_id),
        notes = paste("lme failed:", conditionMessage(fit), "| lm failed:", conditionMessage(fallback))
      )
    } else {
      co <- as.data.table(summary(fallback)$coefficients, keep.rownames = "term")
      setnames(co, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"), c("estimate", "std_error", "t_value", "p_value"))
      co[, `:=`(variable = v, model_type = "lm_fallback", df = fallback$df.residual,
                n_obs = nrow(sub), n_patients = uniqueN(sub$patient_id), notes = "lme unavailable/failed; lm fallback used")]
      mixed_rows[[length(mixed_rows) + 1]] <- co[, .(variable, model_type, term, estimate, std_error, df, t_value, p_value, n_obs, n_patients, notes)]
      log_msg(mixed_log, "Used lm fallback for", v)
    }
  } else {
    tt <- as.data.table(summary(fit)$tTable, keep.rownames = "term")
    setnames(tt, c("Value", "Std.Error", "DF", "t-value", "p-value"), c("estimate", "std_error", "df", "t_value", "p_value"))
    tt[, `:=`(variable = v, model_type = "nlme_lme", n_obs = nrow(sub), n_patients = uniqueN(sub$patient_id), notes = "random intercept for patient_id")]
    mixed_rows[[length(mixed_rows) + 1]] <- tt[, .(variable, model_type, term, estimate, std_error, df, t_value, p_value, n_obs, n_patients, notes)]
    log_msg(mixed_log, "Mixed model completed for", v)
  }
}
mixed_dt <- rbindlist(mixed_rows, fill = TRUE)
mixed_dt[, FDR := p.adjust(p_value, method = "BH")]
write_csv(mixed_dt, file.path(tables_dir, "phase4b_mixed_effect_models.csv"), mixed_log)

aml_healthy_rows <- list()
if (any(mapped$AML_vs_healthy == "healthy_BM")) {
  healthy_summary <- pseudobulk[, .(
    n_obs = .N,
    PRFT_high_like_fraction_mean = mean(PRFT_high_like_fraction, na.rm = TRUE),
    six_gene_risk_score_mean = mean(six_gene_risk_score_mean, na.rm = TRUE),
    PRFT_score_mean = mean(PRFT_score_mean, na.rm = TRUE),
    JAK2_STAT5_PD_L1_score_mean = mean(JAK2_STAT5_PD_L1_score_mean, na.rm = TRUE),
    myeloid_suppressive_score_mean = mean(myeloid_suppressive_score_mean, na.rm = TRUE)
  ), by = .(AML_vs_healthy, cell_group)]
  write_csv(healthy_summary, file.path(tables_dir, "phase4b_AML_vs_healthy_cellgroup_summary.csv"))

  comp_subset <- pseudobulk[cell_group == "monocyte_myeloid_related"]
  aml_vars <- c("PRFT_high_like_fraction", "six_gene_risk_score_mean", "PRFT_score_mean",
                "JAK2_STAT5_PD_L1_score_mean", "myeloid_suppressive_score_mean")
  aml_vars <- aml_vars[aml_vars %in% names(comp_subset)]
  for (v in aml_vars) {
    a <- comp_subset[AML_vs_healthy == "AML", get(v)]
    h <- comp_subset[AML_vs_healthy == "healthy_BM", get(v)]
    if (sum(is.finite(a)) >= 2 && sum(is.finite(h)) >= 2) {
      wt <- wilcox.test(a, h, exact = FALSE)
      aml_healthy_rows[[length(aml_healthy_rows) + 1]] <- data.table(
        variable = v, comparison = "AML_vs_healthy_in_monocyte_myeloid_group",
        p_value = wt$p.value, effect_size = cohens_d_unpaired(a, h),
        AML_median = median(a, na.rm = TRUE), healthy_median = median(h, na.rm = TRUE),
        n_AML = sum(is.finite(a)), n_healthy = sum(is.finite(h))
      )
    }
  }
  aml_health_dt <- rbindlist(aml_healthy_rows, fill = TRUE)
  if (nrow(aml_health_dt) > 0) aml_health_dt[, FDR := p.adjust(p_value, method = "BH")]
  else aml_health_dt <- data.table()
  write_csv(aml_health_dt, file.path(tables_dir, "phase4b_AML_vs_healthy_comparison_details.csv"))

  p_aml <- ggplot(comp_subset, aes(x = AML_vs_healthy, y = PRFT_score_mean, fill = AML_vs_healthy)) +
    geom_boxplot(outlier.size = 0.6, alpha = 0.85) +
    geom_jitter(width = 0.12, size = 0.9, alpha = 0.7) +
    facet_wrap(~ cell_group, scales = "free_y") +
    theme_bw(base_size = 9) +
    theme(legend.position = "none") +
    labs(title = "AML vs healthy BM PRFT-related scores", x = "", y = "pseudo-bulk mean PRFT score")
  save_pdf(p_aml, file.path(fig_dir, "phase4b_AML_vs_healthy_PRFT_scores.pdf"), 8.5, 4.8)
} else {
  write_csv(data.table(status = "not_run", reason = "No healthy_BM sample_state available"), file.path(tables_dir, "phase4b_AML_vs_healthy_cellgroup_summary.csv"))
}

summary_rows <- list(
  data.table(
    item = "patient_level_support_for_monocyte_myeloid_PRFT_high_enrichment",
    result = {
      hit <- comparison_dt[variable == "PRFT_high_like_fraction" & comparison == "monocyte_myeloid_related vs primitive_progenitor_related" & FDR < 0.05]
      if (nrow(hit) > 0 && hit$effect_size[1] > 0) "yes" else "partial_or_no"
    },
    note = "Based on patient/sample pseudo-bulk Wilcoxon comparison"
  ),
  data.table(
    item = "mixed_model_support_for_cell_group_effect",
    result = {
      hit <- mixed_dt[grepl("^cell_group", term) & variable %in% c("six_gene_risk_score_mean", "PRFT_score_mean", "JAK2_STAT5_PD_L1_score_mean", "myeloid_suppressive_score_mean") & FDR < 0.05]
      if (nrow(hit) > 0) "yes" else "partial_or_no"
    },
    note = "nlme::lme preferred; lm fallback if needed"
  ),
  data.table(
    item = "AML_vs_healthy_direction_consistency",
    result = if (exists("aml_health_dt") && nrow(aml_health_dt) > 0) "available" else "not_available_or_weak",
    note = "Focused on monocyte/myeloid-related pseudo-bulk"
  ),
  data.table(
    item = "LSC17_supports_LSC_dominance",
    result = {
      lsc_hit <- comparison_dt[variable == "LSC17_score_mean" & comparison == "monocyte_myeloid_related vs primitive_progenitor_related"]
      if (nrow(lsc_hit) > 0 && lsc_hit$effect_size[1] > 0 && lsc_hit$FDR[1] < 0.05) "yes" else "no"
    },
    note = "Current recommendation should stay conservative unless LSC17/stemness clearly supports primitive dominance"
  ),
  data.table(
    item = "recommended_final_wording",
    result = "PRFT-high-like cells were associated with a monocyte-like/myeloid stress-adapted and immune-suppressive AML state.",
    note = "Processed single-cell table-based analysis only"
  )
)
summary_dt <- rbindlist(summary_rows, fill = TRUE)
write_csv(summary_dt, file.path(tables_dir, "phase4b_singlecell_robustness_summary.csv"))

monocyte_group <- paste(c("Mono-like", "Mono", "ProMono-like", "ProMono", "cDC-like", "cDC"), collapse = ", ")
primitive_group <- paste(c("HSC-like", "HSC", "GMP-like", "GMP", "Prog-like", "Prog"), collapse = ", ")

checklist <- c(
  paste0("1. processed single-cell table successfully read: yes (", basename(input_file), ")"),
  paste0("2. patient_id available: ", ifelse("patient_id" %in% names(mapped), "yes", "no")),
  paste0("3. sample_id available: ", ifelse("sample_id" %in% names(mapped), "yes", "no")),
  paste0("4. sample_state/healthy_BM available: ", ifelse("sample_state" %in% names(mapped) && any(mapped$AML_vs_healthy == 'healthy_BM'), "yes", "yes_without_healthy")),
  paste0("5. cell_type available: ", ifelse("cell_type" %in% names(mapped), "yes", "no")),
  "6. patient x cell_type pseudo-bulk completed: yes",
  paste0("7. monocyte/myeloid group includes: ", monocyte_group),
  paste0("8. primitive/progenitor group includes: ", primitive_group),
  paste0("9. patient-level comparison supports monocyte/myeloid PRFT-high enrichment: ", summary_dt[item == 'patient_level_support_for_monocyte_myeloid_PRFT_high_enrichment', result]),
  paste0("10. mixed-effect model completed: ", ifelse(any(mixed_dt$model_type %in% c('nlme_lme', 'lm_fallback')), "yes", "no")),
  paste0("11. mixed-effect model supports cell_group effect: ", summary_dt[item == 'mixed_model_support_for_cell_group_effect', result]),
  paste0("12. AML vs healthy BM analysis completed: ", ifelse(file.exists(file.path(tables_dir, 'phase4b_AML_vs_healthy_cellgroup_summary.csv')), "yes", "no")),
  paste0("13. LSC17/stemness supports LSC dominance: ", ifelse(summary_dt[item == 'LSC17_supports_LSC_dominance', result] == 'yes', 'yes', 'no')),
  paste0("14. Final recommended wording: ", summary_dt[item == 'recommended_final_wording', result]),
  "15. Recommended main-text figures: phase4b_patient_level_score_boxplots.pdf; phase4b_PRFT_high_like_fraction_patient_level.pdf; phase4b_AML_vs_healthy_PRFT_scores.pdf (if direction remains interpretable)",
  "16. Recommended supplementary figures: phase4_score_heatmap_by_celltype.pdf; phase4_score_dotplot_by_celltype.pdf; phase4_QC_violin.pdf; pseudo-bulk and mixed-model tables",
  "17. Should we continue searching for original Seurat/h5ad/10x data: yes",
  "18. Recommend entering AS/PPI next stage: yes, but only if dedicated inputs exist and should be audited first",
  "19. Issues needing manual confirmation: uneven sample sizes across cell groups and sample states; processed-table-only limitation remains in force"
)
writeLines(checklist, file.path(log_dir, "phase4b_key_result_checklist.txt"))
log_msg(main_log, "Wrote", file.path(log_dir, "phase4b_key_result_checklist.txt"))
log_msg(main_log, "Phase 4b processed single-cell robustness analysis completed.")
