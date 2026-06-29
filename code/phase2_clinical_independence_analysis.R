#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(1234)

ascii_default_lib <- "phase1_R_libs; local path removed"
ascii_env_lib <- Sys.getenv("PHASE1_ASCII_R_LIB", unset = "")
lib_candidates <- unique(c(ascii_env_lib, ascii_default_lib))
lib_candidates <- lib_candidates[nzchar(lib_candidates) & dir.exists(lib_candidates)]
if (length(lib_candidates) > 0) {
  .libPaths(c(lib_candidates, .libPaths()))
}

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
})

root_env <- Sys.getenv("PHASE1_AUDIT_ROOT", unset = "")
if (nzchar(root_env)) {
  root_dir <- chartr("\\", "/", path.expand(root_env))
} else {
  root_dir <- chartr("\\", "/", getwd())
}
if (!dir.exists(file.path(root_dir, "phase1_runtime"))) {
  stop("Run from the project root or set PHASE1_AUDIT_ROOT to the project root.")
}

results_dir <- file.path(root_dir, "03_results_tables")
fig_dir <- file.path(root_dir, "04_figures")
log_dir <- file.path(root_dir, "05_logs")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "phase2_clinical_independence_log.txt")
if (file.exists(log_file)) file.remove(log_file)

append_log <- function(...) {
  line <- paste0(...)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

safe_fread <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path)
  fread(path)
}

blank_to_na <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == "" | toupper(x) %in% c("NA", "N/A", "NULL", "NOT AVAILABLE", "UNKNOWN")] <- NA_character_
  x
}

num_safely <- function(x) {
  suppressWarnings(as.numeric(blank_to_na(x)))
}

fmt_num <- function(x, digits = 3) {
  if (length(x) == 0 || is.na(x)) return("NA")
  formatC(as.numeric(x), digits = digits, format = "f")
}

fmt_p <- function(x) {
  if (length(x) == 0 || is.na(x)) return("NA")
  format(as.numeric(x), digits = 3, scientific = TRUE)
}

count_non_missing <- function(x) sum(!is.na(x))
missing_rate <- function(x) mean(is.na(x))

parse_positive_negative <- function(txt, positive_pattern, negative_pattern) {
  txt <- blank_to_na(txt)
  out <- rep(NA_character_, length(txt))
  out[!is.na(txt) & grepl(positive_pattern, txt, ignore.case = TRUE)] <- "mutated_or_positive"
  out[!is.na(txt) & is.na(out) & grepl(negative_pattern, txt, ignore.case = TRUE)] <- "wildtype_or_negative"
  factor(out, levels = c("wildtype_or_negative", "mutated_or_positive"))
}

parse_idh1 <- function(txt) {
  txt <- blank_to_na(txt)
  out <- rep(NA_character_, length(txt))
  pos <- "IDH1 R132 Positive|IDH1 R140 Positive|IDH1 R172 Positive"
  neg <- "IDH1 R132 Negative|IDH1 R140 Negative|IDH1 R172 Negative"
  out[!is.na(txt) & grepl(pos, txt, ignore.case = TRUE)] <- "mutated_or_positive"
  out[!is.na(txt) & is.na(out) & grepl(neg, txt, ignore.case = TRUE)] <- "wildtype_or_negative"
  factor(out, levels = c("wildtype_or_negative", "mutated_or_positive"))
}

make_factor <- function(x) {
  x <- blank_to_na(x)
  factor(x)
}

theme_set(
  theme_classic(base_size = 8) +
    theme(
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.title = element_blank()
    )
)

append_log("[Phase2] Started at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
append_log("[Phase2] Figure contract: clinical independence and utility; R backend only; outputs are PDF plus source tables.")
append_log("[Phase2] No Phase 1-fix files will be modified; no six-gene formula/model rebuilding will be performed.")

risk_path <- file.path(root_dir, "phase1_runtime", "07_signature", "tcga_cross_platform_risk_score_by_sample.csv")
raw_clin_path <- file.path(root_dir, "00_raw_data", "00_raw_data", "tcga_laml_raw_xena", "TCGA.LAML.sampleMap_LAML_clinicalMatrix.tsv")
coef_path <- file.path(results_dir, "phase1_six_gene_coefficients.csv")
if (!file.exists(coef_path)) {
  coef_path <- file.path(root_dir, "phase1_runtime", "07_signature", "final_cross_platform_prft_signature_coefficients.csv")
}

phase1_fix_surv <- safe_fread(file.path(results_dir, "phase1_fix_TCGA_survival_summary.csv"))
phase1_sample_flow <- safe_fread(file.path(results_dir, "phase1_fix_sample_flow.csv"))
coef_dt <- safe_fread(coef_path)
risk_dt <- safe_fread(risk_path)
raw_clin <- safe_fread(raw_clin_path)

risk_dt <- risk_dt[!is.na(OS_time) & !is.na(OS_status) & OS_time > 0 & !is.na(risk_score)]
raw_keep <- raw_clin[, .(
  sample_id = sampleID,
  age_raw = num_safely(age_at_initial_pathologic_diagnosis),
  sex_raw = blank_to_na(gender),
  WBC_raw = num_safely(lab_procedure_leukocyte_result_unspecified_value),
  PB_blast_percent = num_safely(lab_procedure_blast_cell_outcome_percentage_value),
  BM_blast_percent = num_safely(lab_procedure_bone_marrow_blast_cell_outcome_percent_value),
  FAB_subtype = blank_to_na(leukemia_french_american_british_morphology_code),
  cytogenetic_risk = blank_to_na(acute_myeloid_leukemia_calgb_cytogenetics_risk_category),
  cytogenetic_abnormality = blank_to_na(cytogenetic_abnormality),
  molecular_result = blank_to_na(molecular_analysis_abnormality_testing_result)
)]

clin <- merge(risk_dt, raw_keep, by = "sample_id", all.x = TRUE, sort = FALSE)
clin[, age := fifelse(is.na(age), age_raw, age)]
clin[, sex := fifelse(is.na(blank_to_na(sex)), sex_raw, blank_to_na(sex))]
clin[, WBC := fifelse(is.na(num_safely(WBC)), WBC_raw, num_safely(WBC))]
clin[, FAB := fifelse(is.na(blank_to_na(FAB)), FAB_subtype, blank_to_na(FAB))]
clin[, risk_group := factor(risk_group, levels = c("low_risk", "high_risk"))]
clin[, sex := factor(toupper(sex))]
clin[, WBC_log10 := ifelse(!is.na(WBC) & WBC >= 0, log10(WBC + 1), NA_real_)]
clin[, FAB_subtype := factor(FAB)]
clin[, cytogenetic_risk := factor(cytogenetic_risk)]
clin[, FLT3_status := parse_positive_negative(molecular_result, "FLT3 Mutation Positive", "FLT3 Mutation Negative")]
clin[, NPM1_status := parse_positive_negative(molecular_result, "NPMc Positive", "NPMc Negative")]
clin[, IDH1_status := parse_idh1(molecular_result)]
clin[, IDH2_status := factor(rep(NA_character_, .N), levels = c("wildtype_or_negative", "mutated_or_positive"))]
clin[, DNMT3A_status := factor(rep(NA_character_, .N), levels = c("wildtype_or_negative", "mutated_or_positive"))]
clin[, TP53_status := factor(rep(NA_character_, .N), levels = c("wildtype_or_negative", "mutated_or_positive"))]
clin[, RUNX1_status := factor(rep(NA_character_, .N), levels = c("wildtype_or_negative", "mutated_or_positive"))]
clin[, CEBPA_status := factor(rep(NA_character_, .N), levels = c("wildtype_or_negative", "mutated_or_positive"))]
clin[, TET2_status := factor(rep(NA_character_, .N), levels = c("wildtype_or_negative", "mutated_or_positive"))]
clin[, ASXL1_status := factor(rep(NA_character_, .N), levels = c("wildtype_or_negative", "mutated_or_positive"))]
clin[, ELN_risk := factor(rep(NA_character_, .N))]

append_log("[Phase2] TCGA risk-score analysis set: ", nrow(clin), " samples.")
append_log("[Phase2] Local clinical matrix fields used: age, gender, WBC, PB blast, BM blast, FAB, CALGB cytogenetic risk, cytogenetic abnormality, molecular testing string.")
append_log("[Phase2] Molecular field supports FLT3, NPMc proxy for NPM1, and IDH1. DNMT3A/TP53/RUNX1/IDH2/CEBPA/TET2/ASXL1 and ELN risk were not found as directly codable local fields.")
append_log("[Phase2] Table 1 uses Wilcoxon/t-test for continuous variables and Chi-square/Fisher tests for categorical variables; simulated Fisher is used only when exact Fisher fails for larger sparse tables.")

total_n <- nrow(clin)
variable_defs <- data.table(
  variable = c(
    "risk_score", "age", "sex", "WBC", "WBC_log10", "BM_blast_percent", "PB_blast_percent",
    "FAB_subtype", "cytogenetic_risk", "cytogenetic_abnormality", "ELN_risk",
    "FLT3_status", "NPM1_status", "DNMT3A_status", "TP53_status", "RUNX1_status",
    "IDH1_status", "IDH2_status", "CEBPA_status", "TET2_status", "ASXL1_status"
  ),
  variable_type = c(
    "continuous", "continuous", "categorical", "continuous", "continuous", "continuous", "continuous",
    "categorical", "categorical", "categorical", "categorical",
    "categorical", "categorical", "categorical", "categorical", "categorical",
    "categorical", "categorical", "categorical", "categorical", "categorical"
  )
)

available_detail <- function(var, type) {
  x <- clin[[var]]
  n_non <- count_non_missing(x)
  n_miss <- total_n - n_non
  miss_prop <- n_miss / total_n
  levels_nonmissing <- if (type == "categorical") length(unique(na.omit(as.character(x)))) else NA_integer_
  recommended <- n_non >= 120
  reason <- "Adequate availability for main model."
  if (var %in% c("risk_score", "age", "sex")) {
    recommended <- n_non >= 140
    reason <- "Core adjustment variable."
  }
  if (var == "WBC_log10") {
    recommended <- n_non >= 120
    reason <- "WBC was log10-transformed for modeling because of expected skewness."
  }
  if (var == "WBC") {
    recommended <- FALSE
    reason <- "Raw WBC retained for Table 1 and plots; WBC_log10 is used in Cox models."
  }
  if (var %in% c("BM_blast_percent", "PB_blast_percent", "FAB_subtype", "cytogenetic_risk", "FLT3_status", "NPM1_status", "IDH1_status")) {
    if (n_non >= 100 && (type == "continuous" || levels_nonmissing >= 2)) {
      recommended <- TRUE
      reason <- "Available for main or sensitivity modeling if complete-case size remains adequate."
    } else {
      recommended <- FALSE
      reason <- "Insufficient non-missing observations or non-informative levels for main model."
    }
  }
  if (var %in% c("ELN_risk", "DNMT3A_status", "TP53_status", "RUNX1_status", "IDH2_status", "CEBPA_status", "TET2_status", "ASXL1_status")) {
    recommended <- FALSE
    reason <- "No directly codable local TCGA field found; not forced into any multivariable model."
  }
  data.table(
    variable = var,
    variable_type = type,
    total_samples = total_n,
    non_missing_samples = n_non,
    missing_samples = n_miss,
    missing_proportion = miss_prop,
    observed_levels = ifelse(is.na(levels_nonmissing), NA_character_, as.character(levels_nonmissing)),
    recommend_for_main_model = ifelse(recommended, "yes", "no"),
    reason_not_in_main_model = ifelse(recommended, "", reason)
  )
}

missingness_dt <- rbindlist(Map(available_detail, variable_defs$variable, variable_defs$variable_type))
fwrite(missingness_dt, file.path(results_dir, "phase2_clinical_variable_missingness.csv"))

continuous_vars <- c("age", "WBC", "BM_blast_percent", "PB_blast_percent", "risk_score")
categorical_vars <- c("sex", "FAB_subtype", "cytogenetic_risk", "FLT3_status", "NPM1_status", "IDH1_status")

format_cont <- function(x) {
  if (sum(!is.na(x)) == 0) return("NA")
  paste0(fmt_num(median(x, na.rm = TRUE)), " [", fmt_num(quantile(x, 0.25, na.rm = TRUE)), ", ", fmt_num(quantile(x, 0.75, na.rm = TRUE)), "]")
}

table1_cont <- function(var) {
  dt <- clin[!is.na(get(var)) & !is.na(risk_group)]
  low <- dt[risk_group == "low_risk", get(var)]
  high <- dt[risk_group == "high_risk", get(var)]
  if (length(low) < 3 || length(high) < 3) {
    p <- NA_real_
    test <- "not_tested"
  } else {
    shapiro_ok <- FALSE
    if (length(low) <= 5000 && length(high) <= 5000) {
      shapiro_ok <- tryCatch(shapiro.test(low)$p.value > 0.05 && shapiro.test(high)$p.value > 0.05, error = function(e) FALSE)
    }
    if (shapiro_ok) {
      p <- tryCatch(t.test(get(var) ~ risk_group, data = dt)$p.value, error = function(e) NA_real_)
      test <- "t-test"
    } else {
      p <- tryCatch(wilcox.test(get(var) ~ risk_group, data = dt)$p.value, error = function(e) NA_real_)
      test <- "Wilcoxon rank-sum"
    }
  }
  data.table(
    variable = var,
    level = "median [IQR]",
    overall = format_cont(dt[[var]]),
    low_risk = format_cont(low),
    high_risk = format_cont(high),
    p_value = p,
    test = test
  )
}

table1_cat <- function(var) {
  dt <- clin[!is.na(get(var)) & !is.na(risk_group)]
  if (nrow(dt) == 0 || length(unique(dt[[var]])) < 2) {
    return(data.table(variable = var, level = NA_character_, overall = NA_character_, low_risk = NA_character_, high_risk = NA_character_, p_value = NA_real_, test = "not_tested"))
  }
  tab <- table(dt[[var]], dt$risk_group)
  expected <- suppressWarnings(chisq.test(tab)$expected)
  if (any(expected < 5)) {
    p <- tryCatch(
      fisher.test(tab)$p.value,
      error = function(e) fisher.test(tab, simulate.p.value = TRUE, B = 10000)$p.value
    )
    test <- ifelse(all(dim(tab) == c(2, 2)), "Fisher exact", "Fisher exact or simulated Fisher for larger table")
  } else {
    p <- chisq.test(tab, correct = FALSE)$p.value
    test <- "Chi-square"
  }
  levels_here <- rownames(tab)
  out <- rbindlist(lapply(levels_here, function(lev) {
    overall_n <- sum(dt[[var]] == lev, na.rm = TRUE)
    low_n <- sum(dt[[var]] == lev & dt$risk_group == "low_risk", na.rm = TRUE)
    high_n <- sum(dt[[var]] == lev & dt$risk_group == "high_risk", na.rm = TRUE)
    data.table(
      variable = var,
      level = lev,
      overall = paste0(overall_n, " (", fmt_num(100 * overall_n / nrow(dt), 1), "%)"),
      low_risk = paste0(low_n, " (", fmt_num(100 * low_n / sum(dt$risk_group == "low_risk"), 1), "%)"),
      high_risk = paste0(high_n, " (", fmt_num(100 * high_n / sum(dt$risk_group == "high_risk"), 1), "%)"),
      p_value = p,
      test = test
    )
  }))
  out
}

table1_dt <- rbindlist(c(lapply(continuous_vars, table1_cont), lapply(categorical_vars, table1_cat)), fill = TRUE)
fwrite(table1_dt, file.path(results_dir, "phase2_clinical_characteristics_by_risk.csv"))

cox_extract <- function(fit, model_name, n, events, ph_p, ph_pass, aic, cindex) {
  fit_sum <- summary(fit)
  coef_dt <- as.data.table(fit_sum$coefficients, keep.rownames = "term")
  ci_dt <- as.data.table(fit_sum$conf.int, keep.rownames = "term")
  out <- merge(coef_dt, ci_dt[, .(term, lower95 = `lower .95`, upper95 = `upper .95`)], by = "term", all.x = TRUE)
  data.table(
    model = model_name,
    term = out$term,
    n = n,
    events = events,
    HR = out$`exp(coef)`,
    lower95 = out$lower95,
    upper95 = out$upper95,
    P_value = out$`Pr(>|z|)`,
    C_index = cindex,
    AIC = aic,
    PH_global_p = ph_p,
    PH_assumption_pass = ph_pass
  )
}

run_cox_model <- function(model_name, vars, data = clin) {
  needed <- c("OS_time", "OS_status", vars)
  dt <- data[complete.cases(data[, ..needed])]
  for (v in vars) {
    if (is.factor(dt[[v]]) || is.character(dt[[v]])) {
      dt[, (v) := factor(get(v))]
      if (length(unique(dt[[v]])) < 2) {
        append_log("[Phase2] ", model_name, " skipped because ", v, " has <2 levels after complete-case filtering.")
        return(NULL)
      }
    }
  }
  if (nrow(dt) < 40 || sum(dt$OS_status == 1, na.rm = TRUE) < 15) {
    append_log("[Phase2] ", model_name, " skipped for insufficient complete-case size/events.")
    return(NULL)
  }
  form <- as.formula(paste("Surv(OS_time, OS_status) ~", paste(vars, collapse = " + ")))
  fit <- tryCatch(
    withCallingHandlers(
      coxph(form, data = dt, x = TRUE, y = TRUE),
      warning = function(w) {
        append_log("[Phase2] ", model_name, " warning: ", conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    append_log("[Phase2] ", model_name, " failed: ", conditionMessage(fit))
    return(NULL)
  }
  zph <- tryCatch(cox.zph(fit), error = function(e) NULL)
  ph_p <- if (is.null(zph)) NA_real_ else as.numeric(zph$table["GLOBAL", "p"])
  ph_pass <- ifelse(is.na(ph_p), "not_tested", ifelse(ph_p >= 0.05, "yes", "no"))
  cindex <- as.numeric(summary(fit)$concordance[1])
  out <- cox_extract(
    fit = fit,
    model_name = model_name,
    n = nrow(dt),
    events = sum(dt$OS_status == 1, na.rm = TRUE),
    ph_p = ph_p,
    ph_pass = ph_pass,
    aic = AIC(fit),
    cindex = cindex
  )
  list(model = model_name, vars = vars, data = dt, fit = fit, table = out, ph = zph)
}

univ_vars <- c("risk_score", "age", "sex", "WBC_log10", "BM_blast_percent", "PB_blast_percent", "FAB_subtype", "cytogenetic_risk", "FLT3_status", "NPM1_status", "IDH1_status")
univ_models <- Filter(Negate(is.null), lapply(univ_vars, function(v) run_cox_model(paste0("univariate_", v), v)))
univ_dt <- rbindlist(lapply(univ_models, `[[`, "table"), fill = TRUE)
fwrite(univ_dt, file.path(results_dir, "phase2_univariate_clinical_cox.csv"))

model_specs <- list()
model_specs[["Model 1"]] <- c("risk_score", "age", "sex")
model_specs[["Model 2"]] <- c("risk_score", "age", "sex", "WBC_log10")
if (count_non_missing(clin$cytogenetic_risk) >= 100 && length(unique(na.omit(clin$cytogenetic_risk))) >= 2) {
  model_specs[["Model 3"]] <- c("risk_score", "age", "sex", "WBC_log10", "cytogenetic_risk")
  model3_note <- "cytogenetic_risk used because it is clinically stronger than FAB and had adequate availability."
} else if (count_non_missing(clin$FAB_subtype) >= 100 && length(unique(na.omit(clin$FAB_subtype))) >= 2) {
  model_specs[["Model 3"]] <- c("risk_score", "age", "sex", "WBC_log10", "FAB_subtype")
  model3_note <- "FAB_subtype used because cytogenetic_risk was not adequate."
} else {
  model_specs[["Model 3"]] <- NULL
  model3_note <- "Neither cytogenetic_risk nor FAB_subtype was adequate for Model 3."
}
mut_vars <- c()
for (v in c("FLT3_status", "NPM1_status", "IDH1_status")) {
  if (count_non_missing(clin[[v]]) >= 100 && length(unique(na.omit(clin[[v]]))) >= 2) mut_vars <- c(mut_vars, v)
}
if (length(mut_vars) > 0) {
  model_specs[["Model 4"]] <- c("risk_score", "age", "sex", "WBC_log10", mut_vars)
} else {
  model_specs[["Model 4"]] <- NULL
}
if (count_non_missing(clin$ELN_risk) >= 50 && length(unique(na.omit(clin$ELN_risk))) >= 2) {
  model_specs[["Model 5"]] <- c("risk_score", "age", "sex", "ELN_risk")
} else {
  model_specs[["Model 5"]] <- NULL
}

append_log("[Phase2] Model 3 rule: ", model3_note)
append_log("[Phase2] Model 4 mutation variables used: ", ifelse(length(mut_vars) > 0, paste(mut_vars, collapse = ", "), "none"))
append_log("[Phase2] Model 5 ELN risk unavailable locally; model skipped unless a codable ELN_risk field is present.")

multi_models <- list()
for (nm in names(model_specs)) {
  vars <- model_specs[[nm]]
  if (is.null(vars)) {
    append_log("[Phase2] ", nm, " not run because required variable availability was inadequate.")
  } else {
    multi_models[[nm]] <- run_cox_model(nm, vars)
  }
}
multi_models <- Filter(Negate(is.null), multi_models)
multi_dt <- rbindlist(lapply(multi_models, `[[`, "table"), fill = TRUE)
fwrite(multi_dt, file.path(results_dir, "phase2_multivariable_cox_models.csv"))

model_perf <- unique(multi_dt[, .(model, n, events, C_index, AIC, PH_global_p, PH_assumption_pass)])
fwrite(model_perf, file.path(results_dir, "phase2_model_cindex_AIC_PHtest.csv"))

risk_rows <- multi_dt[term == "risk_score"]
best_model_name <- if ("Model 2" %in% names(multi_models)) {
  "Model 2"
} else {
  "Model 1"
}
best_model <- multi_models[[best_model_name]]
nomogram_vars <- data.table(
  selected_model = best_model_name,
  variable = best_model$vars,
  role = ifelse(best_model$vars == "risk_score", "six-gene PRFT risk score", "clinical covariate"),
  rationale = "Selected by conservative availability and clinical interpretability; no variables were imputed."
)
fwrite(nomogram_vars, file.path(results_dir, "phase2_nomogram_model_variables.csv"))

association_rows <- list()
for (v in c("age", "WBC", "BM_blast_percent", "PB_blast_percent")) {
  dt <- clin[!is.na(risk_score) & !is.na(get(v))]
  if (nrow(dt) >= 10) {
    ct <- suppressWarnings(cor.test(dt$risk_score, dt[[v]], method = "spearman"))
    association_rows[[v]] <- data.table(variable = v, variable_type = "continuous", n = nrow(dt), statistic = unname(ct$estimate), p_value = ct$p.value, test = "Spearman correlation")
  }
}
for (v in categorical_vars[-1]) {
  dt <- clin[!is.na(risk_score) & !is.na(get(v))]
  if (nrow(dt) >= 10 && length(unique(dt[[v]])) >= 2) {
    p <- tryCatch(kruskal.test(risk_score ~ get(v), data = dt)$p.value, error = function(e) NA_real_)
    association_rows[[v]] <- data.table(variable = v, variable_type = "categorical", n = nrow(dt), statistic = NA_real_, p_value = p, test = "Kruskal-Wallis")
  }
}
association_dt <- rbindlist(association_rows, fill = TRUE)
fwrite(association_dt, file.path(results_dir, "phase2_risk_score_clinical_association.csv"))

plot_forest <- function(dt, out_path, title_text, terms_filter = NULL) {
  pdt <- copy(dt)
  if (!is.null(terms_filter)) pdt <- pdt[term %in% terms_filter]
  pdt <- pdt[!is.na(HR) & !is.na(lower95) & !is.na(upper95)]
  pdt[, label := paste(model, term, sep = ": ")]
  pdt[, label := factor(label, levels = rev(label))]
  p <- ggplot(pdt, aes(x = HR, y = label)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey55", linewidth = 0.35) +
    geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0.16, linewidth = 0.35, colour = "#3B4A5A") +
    geom_point(size = 1.7, colour = "#B34A4A") +
    scale_x_log10() +
    labs(x = "Hazard ratio (log scale)", y = NULL, title = title_text) +
    theme(axis.text.y = element_text(size = 6.5))
  ggsave(out_path, p, width = 7.2, height = max(3.2, 0.28 * nrow(pdt) + 1.2), device = cairo_pdf, bg = "white")
}

plot_forest(multi_dt[model == "Model 1"], file.path(fig_dir, "phase2_forestplot_model1.pdf"), "Model 1: PRFT risk score adjusted for age and sex")
plot_forest(multi_dt, file.path(fig_dir, "phase2_forestplot_all_models.pdf"), "Risk score stability across clinical adjustment models")

nomogram_success <- FALSE
if (requireNamespace("rms", quietly = TRUE)) {
  nomogram_success <- tryCatch({
    nomo_cols <- c("OS_time", "OS_status", best_model$vars)
    nomo_dt <- best_model$data[, ..nomo_cols]
    dd <- rms::datadist(nomo_dt)
    old_opt <- options(datadist = "dd")
    on.exit(options(old_opt), add = TRUE)
    form <- as.formula(paste("Surv(OS_time, OS_status) ~", paste(best_model$vars, collapse = " + ")))
    rms_fit <- rms::cph(form, data = nomo_dt, x = TRUE, y = TRUE, surv = TRUE, time.inc = 1095)
    surv_fun <- rms::Survival(rms_fit)
    nom <- rms::nomogram(
      rms_fit,
      fun = list(
        function(x) surv_fun(365, x),
        function(x) surv_fun(1095, x),
        function(x) surv_fun(1825, x)
      ),
      funlabel = c("1-year survival", "3-year survival", "5-year survival")
    )
    pdf(file.path(fig_dir, "phase2_nomogram.pdf"), width = 9, height = 6)
    plot(nom, xfrac = 0.42)
    dev.off()
    TRUE
  }, error = function(e) {
    append_log("[Phase2] rms nomogram failed: ", conditionMessage(e))
    FALSE
  })
}
if (!nomogram_success) {
  pdf(file.path(fig_dir, "phase2_nomogram.pdf"), width = 7, height = 4.5)
  plot.new()
  text(0.5, 0.6, "Nomogram not generated", cex = 1.2)
  text(0.5, 0.45, "rms package unavailable or model failed; see log.", cex = 0.8)
  dev.off()
}

predict_risk_at <- function(fit, newdata, time_days) {
  bh <- basehaz(fit, centered = FALSE)
  h0 <- approx(bh$time, bh$hazard, xout = time_days, rule = 2)$y
  lp <- predict(fit, newdata = newdata, type = "lp")
  1 - exp(-h0 * exp(lp))
}

plot_calibration <- function(model_obj, time_days, out_path) {
  dt <- copy(model_obj$data)
  dt[, pred_risk := predict_risk_at(model_obj$fit, dt, time_days)]
  dt <- dt[!is.na(pred_risk)]
  dt[, group := cut(pred_risk, breaks = unique(quantile(pred_risk, probs = seq(0, 1, length.out = 5), na.rm = TRUE)), include.lowest = TRUE)]
  cal_dt <- dt[, {
    sf <- survfit(Surv(OS_time, OS_status) ~ 1, data = .SD)
    ss <- summary(sf, times = time_days, extend = TRUE)
    observed_risk <- ifelse(length(ss$surv) == 0, NA_real_, 1 - ss$surv[1])
    data.table(predicted_risk = mean(pred_risk, na.rm = TRUE), observed_risk = observed_risk, n = .N)
  }, by = group]
  p <- ggplot(cal_dt, aes(x = predicted_risk, y = observed_risk, size = n)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55", linewidth = 0.35) +
    geom_point(colour = "#3B6EA8") +
    geom_line(colour = "#3B6EA8", linewidth = 0.35) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "Mean predicted risk", y = "Observed Kaplan-Meier risk", title = paste0(time_days / 365, "-year calibration"))
  ggsave(out_path, p, width = 4.5, height = 4.2, device = cairo_pdf, bg = "white")
}

calibration_success <- TRUE
for (tp in c(365, 1095, 1825)) {
  out_name <- switch(as.character(tp), "365" = "phase2_calibration_1year.pdf", "1095" = "phase2_calibration_3year.pdf", "1825" = "phase2_calibration_5year.pdf")
  ok <- tryCatch({ plot_calibration(best_model, tp, file.path(fig_dir, out_name)); TRUE }, error = function(e) { append_log("[Phase2] Calibration failed at ", tp, " days: ", conditionMessage(e)); FALSE })
  calibration_success <- calibration_success && ok
}

dca_success <- TRUE
tryCatch({
  dca_time <- 1095
  dca_vars <- unique(c(best_model$vars, "OS_time", "OS_status"))
  dca_dt <- best_model$data[complete.cases(best_model$data[, ..dca_vars])]
  dca_dt <- dca_dt[!(OS_time < dca_time & OS_status == 0)]
  dca_dt[, event_by_time := as.integer(OS_time <= dca_time & OS_status == 1)]
  risk_fit <- coxph(Surv(OS_time, OS_status) ~ risk_score, data = dca_dt)
  clinical_vars <- setdiff(best_model$vars, "risk_score")
  clinical_fit <- coxph(as.formula(paste("Surv(OS_time, OS_status) ~", paste(clinical_vars, collapse = " + "))), data = dca_dt)
  combined_fit <- best_model$fit
  pred_dt <- data.table(
    event = dca_dt$event_by_time,
    risk_score_alone = predict_risk_at(risk_fit, dca_dt, dca_time),
    clinical_model_alone = predict_risk_at(clinical_fit, dca_dt, dca_time),
    combined_model = predict_risk_at(combined_fit, dca_dt, dca_time)
  )
  thresholds <- seq(0.05, 0.80, by = 0.01)
  nb_one <- function(prob, event, pt) {
    treat <- prob >= pt
    tp <- sum(treat & event == 1)
    fp <- sum(treat & event == 0)
    n <- length(event)
    tp / n - fp / n * pt / (1 - pt)
  }
  dca_dt_long <- rbindlist(lapply(thresholds, function(pt) {
    prevalence <- mean(pred_dt$event)
    data.table(
      threshold = pt,
      model = c("risk_score_alone", "clinical_model_alone", "combined_model", "treat_all", "treat_none"),
      net_benefit = c(
        nb_one(pred_dt$risk_score_alone, pred_dt$event, pt),
        nb_one(pred_dt$clinical_model_alone, pred_dt$event, pt),
        nb_one(pred_dt$combined_model, pred_dt$event, pt),
        prevalence - (1 - prevalence) * pt / (1 - pt),
        0
      )
    )
  }))
  p <- ggplot(dca_dt_long, aes(x = threshold, y = net_benefit, colour = model)) +
    geom_line(linewidth = 0.45) +
    scale_colour_manual(values = c(risk_score_alone = "#B34A4A", clinical_model_alone = "#4F6F52", combined_model = "#3B6EA8", treat_all = "grey45", treat_none = "black")) +
    labs(x = "Threshold probability", y = "Net benefit", title = "Decision curve analysis at 3 years")
  ggsave(file.path(fig_dir, "phase2_DCA.pdf"), p, width = 6.2, height = 4.6, device = cairo_pdf, bg = "white")
}, error = function(e) {
  dca_success <<- FALSE
  append_log("[Phase2] DCA failed: ", conditionMessage(e))
  pdf(file.path(fig_dir, "phase2_DCA.pdf"), width = 7, height = 4.5)
  plot.new()
  text(0.5, 0.6, "DCA not generated", cex = 1.2)
  text(0.5, 0.45, conditionMessage(e), cex = 0.75)
  dev.off()
})

pdf(file.path(fig_dir, "phase2_risk_score_clinical_boxplots.pdf"), width = 7, height = 5)
for (v in c("FAB_subtype", "cytogenetic_risk", "FLT3_status", "NPM1_status", "IDH1_status")) {
  dt <- clin[!is.na(risk_score) & !is.na(get(v))]
  if (nrow(dt) >= 10 && length(unique(dt[[v]])) >= 2) {
    p <- ggplot(dt, aes(x = get(v), y = risk_score, fill = get(v))) +
      geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.3, show.legend = FALSE) +
      geom_jitter(width = 0.12, size = 0.8, alpha = 0.55, show.legend = FALSE) +
      labs(x = v, y = "Risk score", title = paste("Risk score by", v)) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1))
    print(p)
  }
}
for (v in c("age", "WBC", "BM_blast_percent", "PB_blast_percent")) {
  dt <- clin[!is.na(risk_score) & !is.na(get(v))]
  if (nrow(dt) >= 10) {
    p <- ggplot(dt, aes(x = get(v), y = risk_score)) +
      geom_point(size = 1.1, alpha = 0.65, colour = "#3B6EA8") +
      geom_smooth(method = "lm", se = TRUE, linewidth = 0.4, colour = "#B34A4A") +
      labs(x = v, y = "Risk score", title = paste("Risk score and", v))
    print(p)
  }
}
dev.off()

main_model_names <- intersect(c("Model 1", "Model 2", "Model 3", "Model 4"), names(multi_models))
stable_sig <- all(risk_rows[model %in% main_model_names, P_value < 0.05], na.rm = TRUE)

get_risk_line <- function(model_name) {
  row <- risk_rows[model == model_name]
  if (nrow(row) == 0) return("not run")
  paste0("HR ", fmt_num(row$HR[1]), " (95% CI ", fmt_num(row$lower95[1]), "-", fmt_num(row$upper95[1]), "), P=", fmt_p(row$P_value[1]), "; n=", row$n[1], ", events=", row$events[1])
}

checklist_lines <- c(
  paste0("1. TCGA进入Phase 2的样本数：", nrow(clin)),
  paste0("2. risk_score样本数：", count_non_missing(clin$risk_score)),
  paste0("3. 有完整OS信息样本数：", sum(!is.na(clin$OS_time) & clin$OS_time > 0 & !is.na(clin$OS_status))),
  paste0("4. age可用样本数：", count_non_missing(clin$age)),
  paste0("5. sex可用样本数：", count_non_missing(clin$sex)),
  paste0("6. WBC可用样本数：", count_non_missing(clin$WBC)),
  paste0("7. FAB可用样本数：", count_non_missing(clin$FAB_subtype)),
  paste0("8. cytogenetic risk可用样本数：", count_non_missing(clin$cytogenetic_risk)),
  paste0("9. ELN risk是否可用：", ifelse(count_non_missing(clin$ELN_risk) > 0, "是", "否")),
  paste0("10. FLT3可用样本数：", count_non_missing(clin$FLT3_status)),
  paste0("11. NPM1可用样本数：", count_non_missing(clin$NPM1_status), "（使用NPMc字段作为NPM1 proxy）"),
  paste0("12. DNMT3A可用样本数：", count_non_missing(clin$DNMT3A_status)),
  paste0("13. TP53可用样本数：", count_non_missing(clin$TP53_status)),
  paste0("14. RUNX1可用样本数：", count_non_missing(clin$RUNX1_status)),
  paste0("15. IDH1/2可用样本数：", count_non_missing(clin$IDH1_status), "（IDH1可解析；IDH2未在本地字段中发现）"),
  paste0("16. Model 1 risk_score HR、95%CI、P值：", get_risk_line("Model 1")),
  paste0("17. Model 2 risk_score HR、95%CI、P值：", get_risk_line("Model 2")),
  paste0("18. Model 3 risk_score HR、95%CI、P值：", get_risk_line("Model 3")),
  paste0("19. Model 4 risk_score HR、95%CI、P值：", get_risk_line("Model 4")),
  paste0("20. Model 5 risk_score HR、95%CI、P值：", get_risk_line("Model 5")),
  paste0("21. risk_score是否在主要多因素模型中稳定显著：", ifelse(stable_sig, "是", "否")),
  paste0("22. 最推荐写入主文的多因素模型：", best_model_name, "（", paste(best_model$vars, collapse = " + "), "）"),
  paste0("23. 是否成功生成nomogram：", ifelse(nomogram_success, "是", "否")),
  paste0("24. 是否成功生成calibration：", ifelse(calibration_success, "是", "否")),
  paste0("25. 是否成功生成DCA：", ifelse(dca_success, "是", "否")),
  paste0("26. 是否建议进入Phase 3：", ifelse(nrow(multi_dt) > 0 && "Model 1" %in% names(multi_models), "是", "否")),
  paste0("27. 需要人工确认的问题：ELN、DNMT3A、TP53、RUNX1、IDH2、CEBPA、TET2、ASXL1未在本地Xena临床矩阵中直接获得；NPM1使用NPMc字段作为proxy；WBC使用log10(WBC+1)进入Cox。")
)
checklist_text <- paste(enc2utf8(checklist_lines), collapse = "\n")
writeBin(charToRaw(paste0(checklist_text, "\n")), file.path(log_dir, "phase2_key_result_checklist.txt"))

append_log("[Phase2] Missingness table written: phase2_clinical_variable_missingness.csv")
append_log("[Phase2] Table 1 written: phase2_clinical_characteristics_by_risk.csv")
append_log("[Phase2] Multivariable models run: ", paste(names(multi_models), collapse = ", "))
append_log("[Phase2] Best model for nomogram/calibration/DCA: ", best_model_name, " with variables: ", paste(best_model$vars, collapse = ", "))
append_log("[Phase2] PH-test summary:")
for (i in seq_len(nrow(model_perf))) {
  append_log("  - ", model_perf$model[i], ": global PH p=", fmt_p(model_perf$PH_global_p[i]), "; pass=", model_perf$PH_assumption_pass[i])
}
append_log("[Phase2] Nomogram success: ", nomogram_success)
append_log("[Phase2] Calibration success: ", calibration_success)
append_log("[Phase2] DCA success: ", dca_success)
append_log("[Phase2] Finished at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
