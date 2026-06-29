#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(1234)

ascii_default_lib <- "C:/Users/Robin-Yang/AppData/Local/Temp/phase1_R_libs"
ascii_env_lib <- Sys.getenv("PHASE1_ASCII_R_LIB", unset = "")
lib_candidates <- unique(c(ascii_env_lib, ascii_default_lib))
lib_candidates <- lib_candidates[nzchar(lib_candidates) & dir.exists(lib_candidates)]
if (length(lib_candidates) > 0) .libPaths(c(lib_candidates, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
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

audit_log <- file.path(log_dir, "phase2b_mutation_audit_log.txt")
if (file.exists(audit_log)) file.remove(audit_log)

append_log <- function(...) {
  line <- paste0(...)
  cat(line, "\n")
  cat(line, "\n", file = audit_log, append = TRUE)
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

num_safely <- function(x) suppressWarnings(as.numeric(blank_to_na(x)))

fmt_num <- function(x, digits = 3) {
  if (length(x) == 0) return(character(0))
  x_num <- suppressWarnings(as.numeric(x))
  out <- rep("NA", length(x_num))
  ok <- !is.na(x_num)
  out[ok] <- formatC(x_num[ok], format = "f", digits = digits)
  out
}

fmt_p <- function(x) {
  if (length(x) == 0) return(character(0))
  x_num <- suppressWarnings(as.numeric(x))
  out <- rep("NA", length(x_num))
  ok <- !is.na(x_num)
  out[ok] <- format(x_num[ok], digits = 3, scientific = TRUE)
  out
}

theme_set(
  theme_classic(base_size = 8) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank()
    )
)

append_log("[Phase2b] Started at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
append_log("[Phase2b] No Phase 2 original statistical result will be modified.")
append_log("[Phase2b] Figure contract: R backend, quantitative clinical grid, PDF/PNG export, clear HR/CI/P labels.")

phase2_models <- safe_fread(file.path(results_dir, "phase2_multivariable_cox_models.csv"))
phase2_perf <- safe_fread(file.path(results_dir, "phase2_model_cindex_AIC_PHtest.csv"))
risk_dt <- safe_fread(file.path(root_dir, "phase1_runtime", "07_signature", "tcga_cross_platform_risk_score_by_sample.csv"))
raw_clin <- safe_fread(file.path(root_dir, "00_raw_data", "00_raw_data", "tcga_laml_raw_xena", "TCGA.LAML.sampleMap_LAML_clinicalMatrix.tsv"))
risk_dt <- risk_dt[!is.na(OS_time) & !is.na(OS_status) & OS_time > 0 & !is.na(risk_score)]

gene_targets <- c("ELN", "DNMT3A", "TP53", "RUNX1", "IDH2", "CEBPA", "TET2", "ASXL1", "FLT3", "NPM1", "NPM", "IDH1")
file_keywords <- "(mutation|maf|clinical|cbio|cBioPortal|gdc|tcga|laml|xena|molecular|genomic|wustl|hiseq)"
candidate_ext <- "\\.(csv|tsv|txt|maf|md|R|r|json|xml|soft|gz)$"
all_files <- list.files(root_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
all_files <- all_files[!grepl("phase1_runtime/17_tmp/R_libs|/R_libs/|\\\\R_libs\\\\", all_files)]
rel_path <- sub(paste0("^", gsub("([\\^\\$\\.\\|\\?\\*\\+\\(\\)\\[\\{\\\\])", "\\\\\\1", root_dir), "/?"), "", chartr("\\", "/", all_files))
candidate_by_name <- grepl(file_keywords, basename(all_files), ignore.case = TRUE) | grepl(file_keywords, rel_path, ignore.case = TRUE)
candidates <- all_files[candidate_by_name]

scan_file <- function(path) {
  path_norm <- chartr("\\", "/", path)
  info <- file.info(path)
  suffix <- tolower(tools::file_ext(path))
  name_hit_genes <- paste(gene_targets[grepl(paste(gene_targets, collapse = "|"), basename(path), ignore.case = TRUE)], collapse = ";")
  content_hit_genes <- character(0)
  content_scanned <- FALSE
  content_error <- NA_character_
  if (!is.na(info$size) && info$size <= 2e7 && grepl(candidate_ext, path, ignore.case = TRUE)) {
    content_scanned <- TRUE
    lines <- tryCatch(readLines(path, warn = FALSE, n = 4000), error = function(e) e)
    if (inherits(lines, "error")) {
      content_error <- conditionMessage(lines)
    } else {
      hay <- paste(lines, collapse = "\n")
      content_hit_genes <- gene_targets[grepl(paste0("\\b", gene_targets, "\\b"), hay, ignore.case = TRUE)]
    }
  }
  parse_status <- "not_parsed"
  useful_fields <- ""
  if (grepl("clinicalMatrix\\.tsv$", path_norm, ignore.case = TRUE)) {
    parse_status <- "parsed_for_phase2b"
    useful_fields <- "FLT3; NPMc_as_NPM1_proxy; IDH1; cytogenetic_risk; FAB; WBC; blasts"
  } else if (length(content_hit_genes) > 0 && grepl("mutation|maf|clinical|molecular", path_norm, ignore.case = TRUE)) {
    parse_status <- "candidate_reviewed_not_used"
    useful_fields <- paste(content_hit_genes, collapse = ";")
  }
  data.table(
    relative_path = sub(paste0("^", root_dir, "/?"), "", path_norm),
    file_name = basename(path),
    suffix = suffix,
    size_bytes = info$size,
    name_keyword_hit = candidate_by_name[match(path, all_files)],
    content_scanned = content_scanned,
    name_gene_hits = name_hit_genes,
    content_gene_hits = paste(content_hit_genes, collapse = ";"),
    parse_status = parse_status,
    useful_fields = useful_fields,
    content_scan_error = content_error
  )
}

mutation_audit <- rbindlist(lapply(candidates, scan_file), fill = TRUE)
mutation_audit[, parsed_rank := as.integer(parse_status == "parsed_for_phase2b")]
mutation_audit[, content_gene_rank := as.integer(nzchar(content_gene_hits))]
setorder(mutation_audit, -parsed_rank, -content_gene_rank, relative_path)
mutation_audit[, c("parsed_rank", "content_gene_rank") := NULL]
fwrite(mutation_audit, file.path(results_dir, "phase2b_mutation_file_audit.csv"))

parse_positive_negative <- function(txt, positive_pattern, negative_pattern) {
  txt <- blank_to_na(txt)
  out <- rep(NA_character_, length(txt))
  out[!is.na(txt) & grepl(positive_pattern, txt, ignore.case = TRUE)] <- "mutated_or_positive"
  out[!is.na(txt) & is.na(out) & grepl(negative_pattern, txt, ignore.case = TRUE)] <- "wildtype_or_negative"
  out
}

parse_idh1 <- function(txt) {
  txt <- blank_to_na(txt)
  out <- rep(NA_character_, length(txt))
  pos <- "IDH1 R132 Positive|IDH1 R140 Positive|IDH1 R172 Positive"
  neg <- "IDH1 R132 Negative|IDH1 R140 Negative|IDH1 R172 Negative"
  out[!is.na(txt) & grepl(pos, txt, ignore.case = TRUE)] <- "mutated_or_positive"
  out[!is.na(txt) & is.na(out) & grepl(neg, txt, ignore.case = TRUE)] <- "wildtype_or_negative"
  out
}

clin_mut <- raw_clin[, .(
  sample_id = sampleID,
  molecular_result = blank_to_na(molecular_analysis_abnormality_testing_result),
  cytogenetic_risk = blank_to_na(acute_myeloid_leukemia_calgb_cytogenetics_risk_category)
)]
clin_mut[, FLT3_status := parse_positive_negative(molecular_result, "FLT3 Mutation Positive", "FLT3 Mutation Negative")]
clin_mut[, NPM1_status := parse_positive_negative(molecular_result, "NPMc Positive", "NPMc Negative")]
clin_mut[, IDH1_status := parse_idh1(molecular_result)]
for (v in c("DNMT3A_status", "TP53_status", "RUNX1_status", "IDH2_status", "CEBPA_status", "TET2_status", "ASXL1_status", "ELN_risk")) {
  clin_mut[, (v) := NA_character_]
}
supp_mut <- merge(risk_dt[, .(sample_id, patient_id, risk_score, risk_group, OS_time, OS_status)], clin_mut, by = "sample_id", all.x = TRUE, sort = FALSE)
fwrite(supp_mut, file.path(results_dir, "phase2b_supplemental_mutation_matrix.csv"))

extra_mut_candidates_reviewed <- mutation_audit[
  parse_status != "parsed_for_phase2b" &
    grepl("maf|mutation|cbio|gdc|wustl", relative_path, ignore.case = TRUE) &
    nzchar(content_gene_hits)
]
extra_mut_files <- mutation_audit[
  parse_status == "parsed_for_phase2b" &
    !grepl("clinicalMatrix\\.tsv$", relative_path, ignore.case = TRUE)
]
append_log("[Phase2b] Candidate mutation/clinical files audited: ", nrow(mutation_audit))
append_log("[Phase2b] Additional local parse-ready mutation files beyond Xena clinical matrix: ", nrow(extra_mut_files))
append_log("[Phase2b] Candidate mutation-like files with target gene hits reviewed but not parsed: ", nrow(extra_mut_candidates_reviewed))
append_log("[Phase2b] Supplemental mutation matrix source: Xena clinicalMatrix molecular_analysis_abnormality_testing_result.")
append_log("[Phase2b] Parsed variables: FLT3, NPMc-as-NPM1 proxy, IDH1. ELN/DNMT3A/TP53/RUNX1/IDH2/CEBPA/TET2/ASXL1 remain locally unavailable.")

term_label <- function(term) {
  map <- c(
    risk_score = "PRFT risk score",
    age = "Age",
    sexMALE = "Male sex",
    WBC_log10 = "WBC, log10(x + 1)",
    "cytogenetic_riskIntermediate/Normal" = "Cytogenetic risk: intermediate/normal",
    "cytogenetic_riskPoor" = "Cytogenetic risk: poor",
    "FLT3_statusmutated_or_positive" = "FLT3 mutated/positive",
    "NPM1_statusmutated_or_positive" = "NPM1/NPMc positive",
    "IDH1_statusmutated_or_positive" = "IDH1 mutated/positive"
  )
  out <- unname(map[term])
  out[is.na(out)] <- term[is.na(out)]
  out
}

make_forest_data <- function(dt) {
  out <- copy(dt)
  out <- out[!is.na(HR) & !is.na(lower95) & !is.na(upper95)]
  out[, term_label := term_label(term)]
  out[, hr_text := paste0(fmt_num(HR), " (", fmt_num(lower95), "-", fmt_num(upper95), ")")]
  out[, p_text := ifelse(P_value < 0.001, "P<0.001", paste0("P=", fmt_p(P_value)))]
  out[, text_label := paste0(hr_text, "; ", p_text)]
  out
}

save_plot_pdf_png <- function(plot, base_path, width, height, dpi = 450) {
  ggsave(paste0(base_path, ".pdf"), plot, width = width, height = height, device = cairo_pdf, bg = "white")
  png(paste0(base_path, ".png"), width = width * dpi, height = height * dpi, res = dpi, type = "cairo", bg = "white")
  print(plot)
  dev.off()
}

plot_model2_forest <- function() {
  pdt <- make_forest_data(phase2_models[model == "Model 2"])
  pdt[, term_label := factor(term_label, levels = rev(c("PRFT risk score", "Age", "Male sex", "WBC, log10(x + 1)")))]
  x_max <- max(pdt$upper95, na.rm = TRUE) * 2.2
  ggplot(pdt, aes(y = term_label, x = HR)) +
    geom_vline(xintercept = 1, linetype = 2, linewidth = 0.35, colour = "grey55") +
    geom_errorbar(aes(xmin = lower95, xmax = upper95), orientation = "y", width = 0.16, linewidth = 0.45, colour = "#334155") +
    geom_point(size = 2.2, colour = "#A33A3A") +
    geom_text(aes(x = x_max, label = text_label), hjust = 1, size = 3.0, colour = "#111827") +
    scale_x_log10(limits = c(0.45, x_max), breaks = c(0.5, 1, 2, 4, 8)) +
    labs(x = "Hazard ratio (log scale)", y = NULL, title = "Model 2: PRFT risk score adjusted for age, sex and WBC") +
    theme(plot.title = element_text(face = "bold", size = 10), axis.text.y = element_text(size = 8.5), plot.margin = margin(8, 18, 8, 8))
}

plot_all_forest <- function() {
  pdt <- make_forest_data(phase2_models)
  pdt[, term_label := factor(term_label, levels = rev(unique(term_label)))]
  x_max <- max(pdt$upper95, na.rm = TRUE) * 3.0
  ggplot(pdt, aes(y = term_label, x = HR, colour = model)) +
    geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3, colour = "grey55") +
    geom_errorbar(aes(xmin = lower95, xmax = upper95), orientation = "y", width = 0.14, linewidth = 0.35) +
    geom_point(size = 1.7) +
    geom_text(aes(x = x_max, label = text_label), hjust = 1, size = 2.3, colour = "#111827") +
    facet_grid(model ~ ., scales = "free_y", space = "free_y") +
    scale_x_log10(limits = c(0.25, x_max), breaks = c(0.25, 0.5, 1, 2, 4, 8, 16)) +
    scale_colour_manual(values = c("Model 1" = "#64748B", "Model 2" = "#A33A3A", "Model 3" = "#3B6EA8", "Model 4" = "#4F6F52")) +
    labs(x = "Hazard ratio (log scale)", y = NULL, title = "Clinical adjustment models 1-4") +
    theme(plot.title = element_text(face = "bold", size = 10), axis.text.y = element_text(size = 6.7), legend.position = "none", plot.margin = margin(8, 18, 8, 8))
}

save_plot_pdf_png(plot_model2_forest(), file.path(fig_dir, "phase2b_forestplot_model2_publication"), 7.6, 3.8)
save_plot_pdf_png(plot_all_forest(), file.path(fig_dir, "phase2b_forestplot_all_models_publication"), 8.4, 8.8)
append_log("[Phase2b] Forest plots exported as PDF and PNG.")

clin <- merge(
  risk_dt,
  raw_clin[, .(
    sample_id = sampleID,
    age_raw = num_safely(age_at_initial_pathologic_diagnosis),
    sex_raw = blank_to_na(gender),
    WBC_raw = num_safely(lab_procedure_leukocyte_result_unspecified_value)
  )],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)
clin[, age := fifelse(is.na(age), age_raw, age)]
clin[, sex := factor(toupper(fifelse(is.na(blank_to_na(sex)), sex_raw, blank_to_na(sex))), levels = c("FEMALE", "MALE"), labels = c("Female", "Male"))]
clin[, WBC := fifelse(is.na(num_safely(WBC)), WBC_raw, num_safely(WBC))]
clin[, WBC_log10 := ifelse(!is.na(WBC) & WBC >= 0, log10(WBC + 1), NA_real_)]
model2_dt <- clin[complete.cases(clin[, .(OS_time, OS_status, risk_score, age, sex, WBC_log10)])]
model2_fit <- coxph(Surv(OS_time, OS_status) ~ risk_score + age + sex + WBC_log10, data = model2_dt, x = TRUE, y = TRUE)

nomogram_suitable <- FALSE
if (requireNamespace("rms", quietly = TRUE)) {
  nomogram_suitable <- tryCatch({
    nomo_dt <- copy(model2_dt[, .(
      OS_time,
      OS_status,
      PRFT_score = risk_score,
      Age = age,
      Sex = sex,
      log10_WBC = WBC_log10
    )])
    dd <- rms::datadist(nomo_dt)
    old_opt <- options(datadist = "dd")
    on.exit(options(old_opt), add = TRUE)
    rms_fit <- rms::cph(Surv(OS_time, OS_status) ~ PRFT_score + Age + Sex + log10_WBC, data = nomo_dt, x = TRUE, y = TRUE, surv = TRUE, time.inc = 1095)
    surv_fun <- rms::Survival(rms_fit)
    nom <- rms::nomogram(
      rms_fit,
      fun = list(
        function(x) surv_fun(365, x),
        function(x) surv_fun(1095, x),
        function(x) surv_fun(1825, x)
      ),
      funlabel = c("1-year survival", "3-year survival", "5-year survival"),
      lp = FALSE
    )
    pdf(file.path(fig_dir, "phase2b_nomogram_model2_publication.pdf"), width = 11.5, height = 7.2, family = "sans")
    plot(nom, xfrac = 0.33, cex.var = 0.82, cex.axis = 0.62, lmgp = 0.20)
    dev.off()
    TRUE
  }, error = function(e) {
    append_log("[Phase2b] Nomogram generation failed: ", conditionMessage(e))
    FALSE
  })
}
if (!nomogram_suitable) {
  pdf(file.path(fig_dir, "phase2b_nomogram_model2_publication.pdf"), width = 8, height = 4.5)
  plot.new()
  text(0.5, 0.58, "Model 2 nomogram not generated", cex = 1.1)
  text(0.5, 0.44, "rms unavailable or plotting failed; see Phase 2b log.", cex = 0.8)
  dev.off()
}

basehaz_surv <- function(fit, newdata, time_days) {
  bh <- basehaz(fit, centered = FALSE)
  h0 <- approx(bh$time, bh$hazard, xout = time_days, rule = 2)$y
  lp <- predict(fit, newdata = newdata, type = "lp")
  exp(-h0 * exp(lp))
}

calibration_one <- function(time_days) {
  dt <- copy(model2_dt)
  dt[, predicted_survival := basehaz_surv(model2_fit, dt, time_days)]
  probs <- unique(quantile(dt$predicted_survival, probs = seq(0, 1, length.out = 5), na.rm = TRUE))
  dt[, group := cut(predicted_survival, breaks = probs, include.lowest = TRUE)]
  out <- dt[, {
    sf <- survfit(Surv(OS_time, OS_status) ~ 1, data = .SD)
    ss <- summary(sf, times = time_days, extend = TRUE)
    obs <- ifelse(length(ss$surv) == 0, NA_real_, ss$surv[1])
    data.table(
      n = .N,
      mean_predicted_survival = mean(predicted_survival, na.rm = TRUE),
      observed_survival = obs,
      absolute_error = abs(mean(predicted_survival, na.rm = TRUE) - obs)
    )
  }, by = group]
  out[, time_year := time_days / 365]
  out
}
cal_dt <- rbindlist(lapply(c(365, 1095, 1825), calibration_one), fill = TRUE)
setcolorder(cal_dt, c("time_year", "group", "n", "mean_predicted_survival", "observed_survival", "absolute_error"))
fwrite(cal_dt, file.path(results_dir, "phase2b_calibration_metrics.csv"))
cal_plot <- ggplot(cal_dt, aes(x = mean_predicted_survival, y = observed_survival, group = 1)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.35, colour = "grey55") +
  geom_line(linewidth = 0.45, colour = "#3B6EA8") +
  geom_point(size = 1.9, colour = "#A33A3A") +
  facet_wrap(~ time_year, labeller = as_labeller(c(`1` = "1-year", `3` = "3-year", `5` = "5-year"))) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Mean predicted survival probability", y = "Observed survival probability", title = "Model 2 calibration")
ggsave(file.path(fig_dir, "phase2b_calibration_model2_1_3_5year.pdf"), cal_plot, width = 7.4, height = 3.2, device = cairo_pdf, bg = "white")

calibration_suitable <- all(cal_dt$absolute_error <= 0.25, na.rm = TRUE)
append_log("[Phase2b] Calibration max absolute error: ", fmt_num(max(cal_dt$absolute_error, na.rm = TRUE)), ". Suitable for main text = ", calibration_suitable)

dca_time <- 1095
dca_dt <- copy(model2_dt)
dca_dt <- dca_dt[!(OS_time < dca_time & OS_status == 0)]
dca_dt[, event_by_time := as.integer(OS_time <= dca_time & OS_status == 1)]
risk_fit <- coxph(Surv(OS_time, OS_status) ~ risk_score, data = dca_dt)
clinical_fit <- coxph(Surv(OS_time, OS_status) ~ age + sex + WBC_log10, data = dca_dt)
combined_fit <- coxph(Surv(OS_time, OS_status) ~ risk_score + age + sex + WBC_log10, data = dca_dt)
pred_dt <- data.table(
  event = dca_dt$event_by_time,
  risk_score_alone = 1 - basehaz_surv(risk_fit, dca_dt, dca_time),
  clinical_model_alone = 1 - basehaz_surv(clinical_fit, dca_dt, dca_time),
  combined_model = 1 - basehaz_surv(combined_fit, dca_dt, dca_time)
)
thresholds <- seq(0.05, 0.80, by = 0.01)
nb_one <- function(prob, event, pt) {
  treat <- prob >= pt
  tp <- sum(treat & event == 1)
  fp <- sum(treat & event == 0)
  n <- length(event)
  tp / n - fp / n * pt / (1 - pt)
}
prevalence <- mean(pred_dt$event)
dca_long <- rbindlist(lapply(thresholds, function(pt) {
  data.table(
    threshold = pt,
    model = c("risk_score alone", "clinical model alone", "combined model", "treat all", "treat none"),
    net_benefit = c(
      nb_one(pred_dt$risk_score_alone, pred_dt$event, pt),
      nb_one(pred_dt$clinical_model_alone, pred_dt$event, pt),
      nb_one(pred_dt$combined_model, pred_dt$event, pt),
      prevalence - (1 - prevalence) * pt / (1 - pt),
      0
    ),
    endpoint = "3-year mortality"
  )
}))
fwrite(dca_long, file.path(results_dir, "phase2b_DCA_net_benefit.csv"))
dca_plot <- ggplot(dca_long, aes(x = threshold, y = net_benefit, colour = model)) +
  geom_line(linewidth = 0.55) +
  scale_colour_manual(values = c("risk_score alone" = "#A33A3A", "clinical model alone" = "#4F6F52", "combined model" = "#3B6EA8", "treat all" = "grey45", "treat none" = "black")) +
  labs(x = "Threshold probability", y = "Net benefit", title = "Decision curve analysis at 3 years")
ggsave(file.path(fig_dir, "phase2b_DCA_3year_publication.pdf"), dca_plot, width = 6.6, height = 4.5, device = cairo_pdf, bg = "white")

dca_suitable <- any(dca_long[model == "combined model" & threshold >= 0.10 & threshold <= 0.50, net_benefit] >
  dca_long[model == "clinical model alone" & threshold >= 0.10 & threshold <= 0.50, net_benefit], na.rm = TRUE)
append_log("[Phase2b] DCA suitable for main text = ", dca_suitable, " based on net benefit over clinical model in threshold range 0.10-0.50.")

png_files <- c(
  file.path(fig_dir, "phase2b_forestplot_model2_publication.png"),
  file.path(fig_dir, "phase2b_forestplot_all_models_publication.png")
)
append_log("[Phase2b] PNG forest outputs: ", paste(basename(png_files), collapse = "; "))

has_extra <- nrow(extra_mut_files) > 0
supp_counts <- function(var) sum(!is.na(supp_mut[[var]]))
checklist <- c(
  paste0("1. 是否找到额外突变文件：", ifelse(has_extra, "是", "否")),
  paste0("2. 是否补充DNMT3A：", ifelse(supp_counts("DNMT3A_status") > 0, "是", "否")),
  paste0("3. 是否补充TP53：", ifelse(supp_counts("TP53_status") > 0, "是", "否")),
  paste0("4. 是否补充RUNX1：", ifelse(supp_counts("RUNX1_status") > 0, "是", "否")),
  paste0("5. 是否补充IDH2：", ifelse(supp_counts("IDH2_status") > 0, "是", "否")),
  paste0("6. 是否补充CEBPA：", ifelse(supp_counts("CEBPA_status") > 0, "是", "否")),
  paste0("7. 是否补充TET2：", ifelse(supp_counts("TET2_status") > 0, "是", "否")),
  paste0("8. 是否补充ASXL1：", ifelse(supp_counts("ASXL1_status") > 0, "是", "否")),
  paste0("9. 是否补充ELN risk：", ifelse(supp_counts("ELN_risk") > 0, "是", "否")),
  paste0("10. Model 2森林图是否重画成功：", ifelse(file.exists(file.path(fig_dir, "phase2b_forestplot_model2_publication.pdf")), "是", "否")),
  paste0("11. Nomogram是否适合主文：", ifelse(nomogram_suitable, "是（Model 2，宽画布并简化Sex标签）", "否，建议补充材料或方法说明")),
  paste0("12. Calibration是否适合主文：", ifelse(calibration_suitable, "是", "否，曲线偏离较大时建议仅补充材料")),
  paste0("13. DCA是否适合主文：", ifelse(dca_suitable, "是", "否，建议补充材料")),
  "14. 建议主文保留哪些Phase 2图：phase2b_forestplot_model2_publication.pdf；可视版DCA如版面允许可入主文。",
  "15. 建议补充材料保留哪些Phase 2图：phase2b_forestplot_all_models_publication.pdf；phase2b_nomogram_model2_publication.pdf；phase2b_calibration_model2_1_3_5year.pdf；phase2b_DCA_3year_publication.pdf；phase2b_supplemental_mutation_matrix.csv。",
  "16. 是否建议进入Phase 3：是",
  "17. 需要人工确认的问题：本地未找到额外可解析MAF/cBioPortal/GDC突变文件；NPM1仍使用NPMc proxy；IDH2及DNMT3A/TP53/RUNX1/CEBPA/TET2/ASXL1/ELN仍缺失；nomogram需人工打开PDF确认无标签重叠。"
)
checklist_con <- file(file.path(log_dir, "phase2b_key_result_checklist.txt"), open = "w", encoding = "UTF-8")
writeLines(enc2utf8(checklist), checklist_con, useBytes = TRUE)
close(checklist_con)

append_log("[Phase2b] Checklist written: phase2b_key_result_checklist.txt")
append_log("[Phase2b] Finished at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
