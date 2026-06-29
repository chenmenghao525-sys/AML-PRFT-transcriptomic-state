#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})
source("15_scripts/plot_label_utils.R")

dir.create("11_drug", recursive = TRUE, showWarnings = FALSE)
dir.create("13_figures", recursive = TRUE, showWarnings = FALSE)
dir.create("14_tables", recursive = TRUE, showWarnings = FALSE)
dir.create("16_logs", recursive = TRUE, showWarnings = FALSE)

# Formula A is locked and must not be manually changed.
coef_candidates <- c(
  "metadata/LOCKED_PRFT_six_gene_formula_A_coefficients.csv",
  "../metadata/LOCKED_PRFT_six_gene_formula_A_coefficients.csv"
)
coef_path <- coef_candidates[file.exists(coef_candidates)][1]
if (is.na(coef_path) || length(coef_path) == 0) {
  stop("locked Formula A coefficient metadata not found. Checked: ", paste(coef_candidates, collapse = "; "))
}
coef_dt <- data.table::fread(coef_path)
expected_signature_genes <- c("CLCN5", "ARHGEF5", "TRIM32", "ITGB2", "SAT1", "ACOX2")
missing_signature_genes <- setdiff(expected_signature_genes, coef_dt$gene_symbol)
if (length(missing_signature_genes) > 0) {
  stop("Locked Formula A coefficient metadata missing genes: ", paste(missing_signature_genes, collapse = ", "))
}
signature_coef <- stats::setNames(coef_dt$coefficient, coef_dt$gene_symbol)[expected_signature_genes]

save_session_info <- function(path) {
  writeLines(capture.output(sessionInfo()), con = path)
}

get_existing_path <- function(candidates, label) {
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit) || length(hit) == 0) {
    stop(label, " not found. Checked: ", paste(candidates, collapse = "; "))
  }
  normalizePath(hit, winslash = "/", mustWork = TRUE)
}

zscore_vector <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(rep(0, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / s
}

zscore_rows <- function(mat) {
  out <- t(apply(mat, 1, zscore_vector))
  rownames(out) <- rownames(mat)
  colnames(out) <- colnames(mat)
  out
}

deduplicate_by_iqr <- function(mat, symbols) {
  symbols <- toupper(trimws(symbols))
  keep <- !is.na(symbols) & nzchar(symbols)
  mat <- mat[keep, , drop = FALSE]
  symbols <- symbols[keep]

  if (!anyDuplicated(symbols)) {
    rownames(mat) <- symbols
    return(list(
      matrix = mat,
      duplicate_records = data.table(
        gene_symbol = character(),
        rows_before = integer(),
        row_kept_index = integer()
      )
    ))
  }

  iqr_vec <- apply(mat, 1, stats::IQR, na.rm = TRUE)
  split_idx <- split(seq_along(symbols), symbols)

  keep_idx <- integer(length(split_idx))
  dup_dt <- vector("list", length(split_idx))

  ii <- 0L
  for (nm in names(split_idx)) {
    ii <- ii + 1L
    idx <- split_idx[[nm]]
    chosen <- idx[which.max(iqr_vec[idx])]
    keep_idx[ii] <- chosen
    dup_dt[[ii]] <- data.table(
      gene_symbol = nm,
      rows_before = length(idx),
      row_kept_index = chosen
    )
  }

  keep_idx <- sort(keep_idx)
  out <- mat[keep_idx, , drop = FALSE]
  rownames(out) <- symbols[keep_idx]
  list(matrix = out, duplicate_records = rbindlist(dup_dt))
}

standardize_drug_category <- function(drug_name) {
  x <- tolower(drug_name)
  cats <- character()

  if (grepl("cytarabine|daunorubicin|idarubicin|azacitidine|decitabine|venetoclax|gilteritinib|midostaurin|quizartinib|sorafenib", x)) {
    cats <- c(cats, "AML_standard_related")
  }
  if (grepl("ruxolitinib|tofacitinib|fedratinib|jak inhibitor|\\bjak\\b", x)) {
    cats <- c(cats, "JAK_STAT_related")
  }
  if (grepl("venetoclax|navitoclax|abt-199|abt-263|bcl", x)) {
    cats <- c(cats, "BCL2_apoptosis_related")
  }
  if (grepl("bortezomib|carfilzomib|ixazomib|hsp90|geldanamycin|tanespimycin|auy922|mg132|proteasome", x)) {
    cats <- c(cats, "Proteostasis_stress_related")
  }
  if (grepl("rapamycin|everolimus|temsirolimus|pi3k|akt inhibitor|mtor", x)) {
    cats <- c(cats, "PI3K_AKT_mTOR_related")
  }
  if (grepl("erastin|rsl3|sorafenib|sulfasalazine|fin56|ferroptosis|glutathione|gpx4", x)) {
    cats <- c(cats, "Oxidative_stress_ferroptosis_adjacent")
  }

  if (length(cats) == 0) {
    return("Uncategorized")
  }
  paste(unique(cats), collapse = ";")
}

interpret_metric_direction <- function(metric_name) {
  metric_name <- tolower(metric_name)
  if (metric_name %in% c("auc", "area_under_curve", "normalized_auc")) {
    return("Higher AUC indicates lower ex vivo sensitivity / greater relative resistance.")
  }
  if (metric_name %in% c("dss", "drug_sensitivity_score")) {
    return("Higher DSS indicates greater ex vivo sensitivity.")
  }
  if (metric_name %in% c("ic50", "ec50")) {
    return("Higher IC50/EC50 indicates lower ex vivo sensitivity / greater relative resistance.")
  }
  "Metric direction needs confirmation."
}

interpret_association <- function(metric_name, rho) {
  if (!is.finite(rho)) {
    return(NA_character_)
  }
  metric_name <- tolower(metric_name)
  if (metric_name %in% c("auc", "area_under_curve", "normalized_auc", "ic50", "ec50")) {
    if (rho > 0) {
      return("Higher risk associated with greater ex vivo resistance")
    }
    return("Higher risk associated with greater ex vivo sensitivity")
  }
  if (metric_name %in% c("dss", "drug_sensitivity_score")) {
    if (rho > 0) {
      return("Higher risk associated with greater ex vivo sensitivity")
    }
    return("Higher risk associated with greater ex vivo resistance")
  }
  NA_character_
}

pick_metric_column <- function(dt) {
  candidate_order <- c(
    "auc", "AUC", "normalized_auc", "normalizedAUC", "area_under_curve",
    "dss", "DSS", "drug_sensitivity_score", "ic50", "IC50", "ec50", "EC50"
  )
  keep <- candidate_order[candidate_order %in% names(dt)]
  if (length(keep) == 0) {
    stop("No supported drug response metric found in BeatAML response table.")
  }
  keep[1]
}

open_graphics_device <- function(file_path, width = 9, height = 6, type = c("pdf", "png")) {
  type <- match.arg(type)
  if (type == "pdf") {
    grDevices::pdf(file_path, width = width, height = height)
  } else {
    grDevices::png(file_path, width = width, height = height, units = "in", res = 300)
  }
}

plot_placeholder <- function(file_path, title_text, body_text, type = c("pdf", "png")) {
  type <- match.arg(type)
  open_graphics_device(file_path, width = 8, height = 5, type = type)
  graphics::plot.new()
  graphics::title(main = title_text)
  graphics::text(0.5, 0.55, body_text, cex = 1)
  grDevices::dev.off()
}

expr_candidates <- c(
  "00_raw_data/beataml/beataml_expr.txt",
  "00_raw_data/BeatAML/beataml_expr.txt",
  "FQD_bioinformatics_work/data/beataml_expr.txt"
)

drug_candidates <- c(
  "00_raw_data/beataml/beataml_auc.txt",
  "00_raw_data/BeatAML/beataml_auc.txt",
  "FQD_bioinformatics_work/data/beataml_auc.txt"
)

expr_file <- get_existing_path(expr_candidates, "BeatAML expression file")
drug_file <- get_existing_path(drug_candidates, "BeatAML drug response file")

manifest_dt <- data.table(
  component = c("expression", "drug_response"),
  file_path = c(expr_file, drug_file),
  file_size_bytes = c(file.info(expr_file)$size, file.info(drug_file)$size)
)
fwrite(manifest_dt, "11_drug/beataml_download_manifest.csv")

message("Reading BeatAML expression: ", expr_file)
expr_dt <- fread(expr_file)

required_expr_cols <- c("display_label")
missing_expr_cols <- setdiff(required_expr_cols, names(expr_dt))
if (length(missing_expr_cols) > 0) {
  stop("BeatAML expression file is missing required columns: ", paste(missing_expr_cols, collapse = ", "))
}

metadata_cols <- intersect(c("stable_id", "display_label", "description", "biotype"), names(expr_dt))
sample_cols <- setdiff(names(expr_dt), metadata_cols)

if (length(sample_cols) == 0) {
  stop("No sample columns detected in BeatAML expression file.")
}

expr_dt <- expr_dt[!is.na(display_label) & trimws(display_label) != ""]
if ("biotype" %in% names(expr_dt)) {
  expr_dt <- expr_dt[is.na(biotype) | biotype == "protein_coding"]
}

expr_mat <- as.matrix(expr_dt[, ..sample_cols])
storage.mode(expr_mat) <- "numeric"
gene_symbols <- toupper(trimws(expr_dt$display_label))

dedup_res <- deduplicate_by_iqr(expr_mat, gene_symbols)
expr_mat <- dedup_res$matrix

coverage_dt <- data.table(
  gene_symbol = names(signature_coef),
  available = names(signature_coef) %in% rownames(expr_mat),
  rows_matched = as.integer(names(signature_coef) %in% rownames(expr_mat))
)
coverage_dt[, note := ifelse(available, "present", "missing")]
fwrite(coverage_dt, "11_drug/beataml_expression_signature_gene_coverage.csv")

if (!all(coverage_dt$available)) {
  stop(
    "Not all six signature genes are available in BeatAML expression matrix. Missing: ",
    paste(coverage_dt$gene_symbol[!coverage_dt$available], collapse = ", ")
  )
}

expr_type <- if (all(abs(expr_mat - round(expr_mat)) < 1e-8, na.rm = TRUE) && max(expr_mat, na.rm = TRUE) > 100) {
  "raw_counts_like"
} else {
  "normalized_expression_like"
}

signature_expr <- expr_mat[names(signature_coef), , drop = FALSE]
signature_z <- zscore_rows(signature_expr)
risk_score <- colSums(signature_z * signature_coef[rownames(signature_z)])
risk_median <- stats::median(risk_score, na.rm = TRUE)
risk_group <- ifelse(risk_score >= risk_median, "high_risk", "low_risk")

risk_dt <- data.table(
  sample_id = colnames(signature_z),
  risk_score = as.numeric(risk_score),
  risk_group = risk_group
)
for (g in rownames(signature_z)) {
  risk_dt[[paste0("z_", g)]] <- as.numeric(signature_z[g, risk_dt$sample_id])
}
setorder(risk_dt, -risk_score)
fwrite(risk_dt, "11_drug/beataml_risk_score_by_sample.csv")

message("Reading BeatAML drug response: ", drug_file)
drug_dt <- fread(drug_file)

sample_col <- c("sample_id", "dbgap_rnaseq_sample", "rnaseq_sample", "sample")[c("sample_id", "dbgap_rnaseq_sample", "rnaseq_sample", "sample") %in% names(drug_dt)][1]
drug_name_col <- c("drug_name", "inhibitor", "compound", "drug")[c("drug_name", "inhibitor", "compound", "drug") %in% names(drug_dt)][1]
metric_col <- pick_metric_column(drug_dt)

if (is.na(sample_col) || is.na(drug_name_col)) {
  stop("Unable to identify sample ID or drug name column in BeatAML drug response file.")
}

if ("paper_inclusion" %in% names(drug_dt)) {
  drug_dt <- drug_dt[isTRUE(paper_inclusion) | paper_inclusion == TRUE]
}
if ("type" %in% names(drug_dt)) {
  drug_dt <- drug_dt[tolower(type) == "single-agent"]
}

drug_clean <- data.table(
  sample_id = trimws(as.character(drug_dt[[sample_col]])),
  drug_name = trimws(as.character(drug_dt[[drug_name_col]])),
  response_value = suppressWarnings(as.numeric(drug_dt[[metric_col]])),
  metric_name = metric_col
)

if ("curve_type" %in% names(drug_dt)) {
  drug_clean[, curve_type := as.character(drug_dt$curve_type)]
}
if ("status" %in% names(drug_dt)) {
  drug_clean[, status := as.character(drug_dt$status)]
}

drug_clean <- drug_clean[
  !is.na(sample_id) & sample_id != "" &
    !is.na(drug_name) & drug_name != "" &
    is.finite(response_value)
]

if (nrow(drug_clean) == 0) {
  stop("No usable BeatAML drug response rows remained after filtering.")
}

metric_direction_note <- interpret_metric_direction(metric_col)
drug_clean[, metric_direction_note := metric_direction_note]

drug_clean <- drug_clean[
  ,
  .(
    response_value = mean(response_value, na.rm = TRUE),
    metric_name = unique(metric_name)[1],
    metric_direction_note = unique(metric_direction_note)[1],
    curve_type = if ("curve_type" %in% names(.SD)) unique(curve_type)[1] else NA_character_,
    status = if ("status" %in% names(.SD)) unique(status)[1] else NA_character_
  ),
  by = .(sample_id, drug_name)
]

fwrite(drug_clean, "11_drug/beataml_drug_response_clean.csv")

matched_dt <- merge(drug_clean, risk_dt, by = "sample_id", all.x = FALSE, all.y = FALSE)

per_drug_summary <- matched_dt[
  ,
  .(
    matched_samples = uniqueN(sample_id),
    high_risk_samples = uniqueN(sample_id[risk_group == "high_risk"]),
    low_risk_samples = uniqueN(sample_id[risk_group == "low_risk"]),
    response_median = stats::median(response_value, na.rm = TRUE)
  ),
  by = .(drug_name, metric_name, metric_direction_note)
]
per_drug_summary[, eligible_for_testing := matched_samples >= 20]

fwrite(per_drug_summary, "11_drug/beataml_drug_matched_sample_summary.csv")

test_drugs <- per_drug_summary[eligible_for_testing == TRUE, drug_name]

correlation_res <- matched_dt[drug_name %in% test_drugs, {
  ok <- is.finite(risk_score) & is.finite(response_value)
  n_ok <- sum(ok)
  if (n_ok < 20) {
    .(
      matched_samples = n_ok,
      spearman_rho = NA_real_,
      P.Value = NA_real_
    )
  } else {
    ct <- suppressWarnings(stats::cor.test(risk_score[ok], response_value[ok], method = "spearman", exact = FALSE))
    .(
      matched_samples = n_ok,
      spearman_rho = unname(ct$estimate),
      P.Value = ct$p.value
    )
  }
}, by = .(drug_name, metric_name, metric_direction_note)]

correlation_res[, FDR := stats::p.adjust(P.Value, method = "BH")]
correlation_res[, association_interpretation := vapply(
  spearman_rho,
  function(x) interpret_association(metric_name[1], x),
  character(1)
), by = .(drug_name, metric_name, metric_direction_note)]

diff_res <- matched_dt[drug_name %in% test_drugs, {
  dt <- .SD[is.finite(response_value) & !is.na(risk_group)]
  med_high <- stats::median(dt$response_value[dt$risk_group == "high_risk"], na.rm = TRUE)
  med_low <- stats::median(dt$response_value[dt$risk_group == "low_risk"], na.rm = TRUE)
  n_high <- uniqueN(dt$sample_id[dt$risk_group == "high_risk"])
  n_low <- uniqueN(dt$sample_id[dt$risk_group == "low_risk"])
  pval <- if (n_high >= 3 && n_low >= 3) {
    suppressWarnings(stats::wilcox.test(response_value ~ risk_group, data = dt)$p.value)
  } else {
    NA_real_
  }
  .(
    matched_samples = uniqueN(dt$sample_id),
    high_risk_samples = n_high,
    low_risk_samples = n_low,
    median_high = med_high,
    median_low = med_low,
    difference = med_high - med_low,
    P.Value = pval
  )
}, by = .(drug_name, metric_name, metric_direction_note)]

diff_res[, FDR := stats::p.adjust(P.Value, method = "BH")]
diff_res[, group_difference_interpretation := fifelse(
  tolower(metric_name) %in% c("auc", "area_under_curve", "normalized_auc", "ic50", "ec50") & difference > 0,
  "Higher high-risk metric suggests greater ex vivo resistance",
  fifelse(
    tolower(metric_name) %in% c("auc", "area_under_curve", "normalized_auc", "ic50", "ec50") & difference < 0,
    "Lower high-risk metric suggests greater ex vivo sensitivity",
    fifelse(
      tolower(metric_name) %in% c("dss", "drug_sensitivity_score") & difference > 0,
      "Higher high-risk metric suggests greater ex vivo sensitivity",
      fifelse(
        tolower(metric_name) %in% c("dss", "drug_sensitivity_score") & difference < 0,
        "Lower high-risk metric suggests greater ex vivo resistance",
        NA_character_
      )
    )
  )
)]

correlation_res[, abs_spearman_rho := abs(spearman_rho)]
diff_res[, abs_difference := abs(difference)]
setorder(correlation_res, FDR, P.Value, -abs_spearman_rho)
setorder(diff_res, FDR, P.Value, -abs_difference)

fwrite(correlation_res, "11_drug/beataml_risk_score_drug_correlation_all.csv")
fwrite(diff_res, "11_drug/beataml_high_low_risk_drug_response_difference_all.csv")

annot_dt <- merge(
  correlation_res,
  diff_res[, .(drug_name, median_high, median_low, difference, wilcox_P.Value = P.Value, wilcox_FDR = FDR, group_difference_interpretation)],
  by = "drug_name",
  all = TRUE
)
annot_dt[, drug_category := vapply(drug_name, standardize_drug_category, character(1))]
annot_dt[, is_key_drug := grepl("venetoclax|ruxolitinib|sorafenib|bortezomib", drug_name, ignore.case = TRUE)]
annot_dt[, category_count := lengths(strsplit(drug_category, ";", fixed = TRUE))]

fwrite(annot_dt, "11_drug/beataml_drug_category_annotated_results.csv")

top_positive <- correlation_res[is.finite(spearman_rho)][order(-spearman_rho, FDR, P.Value)][1:min(10, .N)]
top_negative <- correlation_res[is.finite(spearman_rho)][order(spearman_rho, FDR, P.Value)][1:min(10, .N)]
top_fdr <- correlation_res[is.finite(FDR)][order(FDR, P.Value, -abs(spearman_rho))][1:min(20, .N)]
key_drugs <- annot_dt[is_key_drug == TRUE][order(FDR, P.Value, -abs(spearman_rho))]

top_hits <- rbindlist(list(
  if (nrow(top_fdr) > 0) data.table(hit_group = "top_by_FDR", top_fdr) else NULL,
  if (nrow(top_positive) > 0) data.table(hit_group = "top_positive_correlation", top_positive) else NULL,
  if (nrow(top_negative) > 0) data.table(hit_group = "top_negative_correlation", top_negative) else NULL,
  if (nrow(key_drugs) > 0) data.table(hit_group = "key_drugs", key_drugs) else NULL
), fill = TRUE)
top_hits <- unique(top_hits, by = c("hit_group", "drug_name"))
fwrite(top_hits, "11_drug/beataml_drug_sensitivity_top_hits.csv")

summary_dt <- data.table(
  beataml_expression_samples = ncol(expr_mat),
  beataml_drug_response_samples = uniqueN(drug_clean$sample_id),
  matched_expression_drug_samples = uniqueN(matched_dt$sample_id),
  signature_genes_available = sum(coverage_dt$available),
  drug_response_metric_used = metric_col,
  metric_direction_note = metric_direction_note,
  drugs_tested_total = uniqueN(drug_clean$drug_name),
  drugs_with_matched_samples_ge_20 = length(test_drugs),
  significant_drugs_FDR_less_0.05 = sum(correlation_res$FDR < 0.05, na.rm = TRUE),
  significant_drugs_P_less_0.05 = sum(correlation_res$P.Value < 0.05, na.rm = TRUE),
  significant_drugs_FDR_less_0.05_names = paste(correlation_res[FDR < 0.05, drug_name], collapse = ";"),
  significant_drugs_P_less_0.05_names = paste(correlation_res[P.Value < 0.05, drug_name], collapse = ";"),
  top_positive_correlated_drugs = paste(top_positive$drug_name, collapse = ";"),
  top_negative_correlated_drugs = paste(top_negative$drug_name, collapse = ";"),
  venetoclax_available = any(grepl("venetoclax", drug_clean$drug_name, ignore.case = TRUE)),
  ruxolitinib_available = any(grepl("ruxolitinib", drug_clean$drug_name, ignore.case = TRUE)),
  sorafenib_available = any(grepl("sorafenib", drug_clean$drug_name, ignore.case = TRUE)),
  bortezomib_available = any(grepl("bortezomib", drug_clean$drug_name, ignore.case = TRUE)),
  proteasome_or_HSP90_drugs_available = any(grepl("bortezomib|carfilzomib|ixazomib|hsp90|geldanamycin|tanespimycin|auy922|mg132|proteasome", drug_clean$drug_name, ignore.case = TRUE)),
  ferroptosis_related_drugs_available = any(grepl("erastin|rsl3|sorafenib|sulfasalazine|fin56|ferroptosis|glutathione|gpx4", drug_clean$drug_name, ignore.case = TRUE)),
  exploratory_interpretation_note = "BeatAML drug response analysis was exploratory and based on ex vivo drug sensitivity associations rather than clinical treatment outcomes.",
  expression_data_type = expr_type,
  expression_source = expr_file,
  drug_response_source = drug_file
)
fwrite(summary_dt, "14_tables/beataml_drug_sensitivity_summary.csv")

plot_bubble <- function() {
  plot_dt <- unique(rbindlist(list(
    annot_dt[drug_category != "Uncategorized"],
    top_positive,
    top_negative
  ), fill = TRUE), by = "drug_name")

  plot_dt <- plot_dt[is.finite(spearman_rho)]
  if (nrow(plot_dt) == 0) {
    plot_placeholder(
      "13_figures/Figure10_BeatAML_risk_score_drug_correlation_bubble.pdf",
      "BeatAML risk score vs drug response",
      "No eligible drug correlation results were available.",
      type = "pdf"
    )
    plot_placeholder(
      "13_figures/Figure10_BeatAML_risk_score_drug_correlation_bubble.png",
      "BeatAML risk score vs drug response",
      "No eligible drug correlation results were available.",
      type = "png"
    )
    return(invisible(NULL))
  }

  plot_dt[, neg_log10_fdr := -log10(pmax(FDR, 1e-300))]
  plot_dt[, drug_name := factor(drug_name, levels = plot_dt[order(spearman_rho), drug_name])]

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ggplot2 <- asNamespace("ggplot2")
    p <- ggplot2$ggplot(
      plot_dt,
      ggplot2$aes(x = spearman_rho, y = drug_name, size = matched_samples, color = neg_log10_fdr)
    ) +
      ggplot2$geom_point(alpha = 0.85) +
      ggplot2$geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
      ggplot2$scale_color_gradient(low = "#2C7BB6", high = "#D7191C") +
      ggplot2$labs(
        x = "Spearman rho: risk score vs drug response",
        y = NULL,
        size = "Matched samples",
        color = "-log10(FDR)",
        title = "BeatAML risk score and ex vivo drug response"
      ) +
      ggplot2$theme_bw(base_size = 11)

    ggplot2$ggsave("13_figures/Figure10_BeatAML_risk_score_drug_correlation_bubble.pdf", p, width = 9, height = 7)
    ggplot2$ggsave("13_figures/Figure10_BeatAML_risk_score_drug_correlation_bubble.png", p, width = 9, height = 7, dpi = 300)
  } else {
    open_graphics_device("13_figures/Figure10_BeatAML_risk_score_drug_correlation_bubble.pdf", width = 9, height = 7, type = "pdf")
    graphics::par(mar = c(5, 10, 3, 1))
    graphics::plot(
      plot_dt$spearman_rho,
      seq_len(nrow(plot_dt)),
      pch = 19,
      cex = 0.6 + plot_dt$matched_samples / max(plot_dt$matched_samples) * 1.6,
      col = grDevices::colorRampPalette(c("#2C7BB6", "#D7191C"))(100)[cut(plot_dt$neg_log10_fdr, 100)],
      yaxt = "n",
      xlab = "Spearman rho: risk score vs drug response",
      ylab = "",
      main = "BeatAML risk score and ex vivo drug response"
    )
    graphics::axis(2, at = seq_len(nrow(plot_dt)), labels = plot_dt$drug_name, las = 2)
    graphics::abline(v = 0, lty = 2, col = "grey50")
    grDevices::dev.off()

    open_graphics_device("13_figures/Figure10_BeatAML_risk_score_drug_correlation_bubble.png", width = 9, height = 7, type = "png")
    graphics::par(mar = c(5, 10, 3, 1))
    graphics::plot(
      plot_dt$spearman_rho,
      seq_len(nrow(plot_dt)),
      pch = 19,
      cex = 0.6 + plot_dt$matched_samples / max(plot_dt$matched_samples) * 1.6,
      col = grDevices::colorRampPalette(c("#2C7BB6", "#D7191C"))(100)[cut(plot_dt$neg_log10_fdr, 100)],
      yaxt = "n",
      xlab = "Spearman rho: risk score vs drug response",
      ylab = "",
      main = "BeatAML risk score and ex vivo drug response"
    )
    graphics::axis(2, at = seq_len(nrow(plot_dt)), labels = plot_dt$drug_name, las = 2)
    graphics::abline(v = 0, lty = 2, col = "grey50")
    grDevices::dev.off()
  }
}

plot_boxplot <- function() {
  top_drugs <- diff_res[is.finite(P.Value)][order(FDR, P.Value, -abs(difference))][1:min(6, .N), drug_name]
  plot_dt <- matched_dt[drug_name %in% top_drugs]

  if (length(top_drugs) == 0 || nrow(plot_dt) == 0) {
    plot_placeholder(
      "13_figures/Figure10_BeatAML_high_low_risk_top_drugs_boxplot.pdf",
      "BeatAML high- vs low-risk drug response",
      "No eligible grouped drug comparison results were available.",
      type = "pdf"
    )
    plot_placeholder(
      "13_figures/Figure10_BeatAML_high_low_risk_top_drugs_boxplot.png",
      "BeatAML high- vs low-risk drug response",
      "No eligible grouped drug comparison results were available.",
      type = "png"
    )
    return(invisible(NULL))
  }

  plot_dt[, drug_name := factor(drug_name, levels = top_drugs)]
  plot_dt[, risk_group := pretty_factor(risk_group, levels = c("low_risk", "high_risk"))]

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ggplot2 <- asNamespace("ggplot2")
    p <- ggplot2$ggplot(
      plot_dt,
      ggplot2$aes(x = risk_group, y = response_value, fill = risk_group)
    ) +
      ggplot2$geom_boxplot(outlier.shape = NA, alpha = 0.85) +
      ggplot2$geom_jitter(width = 0.12, alpha = 0.35, size = 0.8) +
      ggplot2$facet_wrap(~drug_name, scales = "free_y") +
      ggplot2$scale_fill_manual(values = c("high risk" = "#D55E00", "low risk" = "#2E86AB")) +
      ggplot2$labs(
        x = NULL,
        y = paste0(metric_col, " response"),
        title = "BeatAML top drug response differences by risk group"
      ) +
      ggplot2$theme_bw(base_size = 11) +
      ggplot2$theme(legend.position = "none")

    ggplot2$ggsave("13_figures/Figure10_BeatAML_high_low_risk_top_drugs_boxplot.pdf", p, width = 11, height = 7)
    ggplot2$ggsave("13_figures/Figure10_BeatAML_high_low_risk_top_drugs_boxplot.png", p, width = 11, height = 7, dpi = 300)
  } else {
    plot_placeholder(
      "13_figures/Figure10_BeatAML_high_low_risk_top_drugs_boxplot.pdf",
      "BeatAML high- vs low-risk drug response",
      "ggplot2 not installed; boxplot fallback not implemented.",
      type = "pdf"
    )
    plot_placeholder(
      "13_figures/Figure10_BeatAML_high_low_risk_top_drugs_boxplot.png",
      "BeatAML high- vs low-risk drug response",
      "ggplot2 not installed; boxplot fallback not implemented.",
      type = "png"
    )
  }
}

plot_category_heatmap <- function() {
  plot_dt <- annot_dt[drug_category != "Uncategorized" & is.finite(spearman_rho)]
  if (nrow(plot_dt) == 0) {
    plot_placeholder(
      "13_figures/Figure10_BeatAML_drug_category_heatmap.pdf",
      "BeatAML category heatmap",
      "No categorized drugs were available for plotting.",
      type = "pdf"
    )
    plot_placeholder(
      "13_figures/Figure10_BeatAML_drug_category_heatmap.png",
      "BeatAML category heatmap",
      "No categorized drugs were available for plotting.",
      type = "png"
    )
    return(invisible(NULL))
  }

  plot_dt <- plot_dt[
    ,
    {
      sp <- strsplit(drug_category, ";", fixed = TRUE)[[1]]
      data.table(category = sp)
    },
    by = .(drug_name, spearman_rho)
  ]
  plot_dt[, category := pretty_factor(category, levels = unique(category))]
  plot_dt <- plot_dt[, .SD[which.max(abs(spearman_rho))], by = .(category, drug_name)]
  plot_dt <- plot_dt[
    ,
    .SD[order(-abs(spearman_rho))[1:min(3, .N)]],
    by = category
  ]

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ggplot2 <- asNamespace("ggplot2")
    p <- ggplot2$ggplot(
      plot_dt,
      ggplot2$aes(x = drug_name, y = category, fill = spearman_rho)
    ) +
      ggplot2$geom_tile(color = "white") +
      ggplot2$scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0) +
      ggplot2$labs(
        x = NULL,
        y = NULL,
        fill = "Spearman rho",
        title = "BeatAML risk score associations across drug categories"
      ) +
      ggplot2$theme_bw(base_size = 11) +
      ggplot2$theme(axis.text.x = ggplot2$element_text(angle = 45, hjust = 1))

    ggplot2$ggsave("13_figures/Figure10_BeatAML_drug_category_heatmap.pdf", p, width = 11, height = 5.5)
    ggplot2$ggsave("13_figures/Figure10_BeatAML_drug_category_heatmap.png", p, width = 11, height = 5.5, dpi = 300)
  } else {
    plot_placeholder(
      "13_figures/Figure10_BeatAML_drug_category_heatmap.pdf",
      "BeatAML category heatmap",
      "ggplot2 not installed; heatmap fallback not implemented.",
      type = "pdf"
    )
    plot_placeholder(
      "13_figures/Figure10_BeatAML_drug_category_heatmap.png",
      "BeatAML category heatmap",
      "ggplot2 not installed; heatmap fallback not implemented.",
      type = "png"
    )
  }
}

plot_selected_drugs <- function() {
  selected_pattern <- "venetoclax|ruxolitinib|sorafenib|bortezomib"
  plot_dt <- matched_dt[grepl(selected_pattern, drug_name, ignore.case = TRUE)]
  if (nrow(plot_dt) == 0) {
    plot_placeholder(
      "13_figures/Figure10_BeatAML_selected_drug_scatter.pdf",
      "BeatAML selected drug scatter",
      "Venetoclax, ruxolitinib, sorafenib, and bortezomib were not available.",
      type = "pdf"
    )
    plot_placeholder(
      "13_figures/Figure10_BeatAML_selected_drug_scatter.png",
      "BeatAML selected drug scatter",
      "Venetoclax, ruxolitinib, sorafenib, and bortezomib were not available.",
      type = "png"
    )
    return(invisible(NULL))
  }
  plot_dt[, risk_group := pretty_factor(risk_group, levels = c("low_risk", "high_risk"))]

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ggplot2 <- asNamespace("ggplot2")
    p <- ggplot2$ggplot(
      plot_dt,
      ggplot2$aes(x = risk_score, y = response_value)
    ) +
      ggplot2$geom_point(alpha = 0.65, color = "#2E86AB") +
      ggplot2$geom_smooth(method = "lm", se = FALSE, color = "#D55E00", linewidth = 0.7) +
      ggplot2$facet_wrap(~drug_name, scales = "free_y") +
      ggplot2$labs(
        x = "Fixed 6-gene risk score",
        y = paste0(metric_col, " response"),
        title = "BeatAML selected drugs: risk score vs response"
      ) +
      ggplot2$theme_bw(base_size = 11)

    ggplot2$ggsave("13_figures/Figure10_BeatAML_selected_drug_scatter.pdf", p, width = 11, height = 5.5)
    ggplot2$ggsave("13_figures/Figure10_BeatAML_selected_drug_scatter.png", p, width = 11, height = 5.5, dpi = 300)
  } else {
    plot_placeholder(
      "13_figures/Figure10_BeatAML_selected_drug_scatter.pdf",
      "BeatAML selected drug scatter",
      "ggplot2 not installed; scatter fallback not implemented.",
      type = "pdf"
    )
    plot_placeholder(
      "13_figures/Figure10_BeatAML_selected_drug_scatter.png",
      "BeatAML selected drug scatter",
      "ggplot2 not installed; scatter fallback not implemented.",
      type = "png"
    )
  }
}

plot_bubble()
plot_boxplot()
plot_category_heatmap()
plot_selected_drugs()

save_session_info("16_logs/sessionInfo_20_beataml_drug_sensitivity_analysis.txt")

message("BeatAML drug sensitivity analysis completed.")
