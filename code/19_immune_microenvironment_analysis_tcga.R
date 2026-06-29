#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
source("15_scripts/plot_label_utils.R")

dir.create("10_immune", recursive = TRUE, showWarnings = FALSE)
dir.create("13_figures", recursive = TRUE, showWarnings = FALSE)
dir.create("14_tables", recursive = TRUE, showWarnings = FALSE)
dir.create("16_logs", recursive = TRUE, showWarnings = FALSE)

save_session_info <- function(path) {
  writeLines(capture.output(sessionInfo()), con = path)
}

find_input_file <- function(default_path, filename = basename(default_path)) {
  if (file.exists(default_path)) {
    return(normalizePath(default_path, winslash = "/", mustWork = TRUE))
  }
  hits <- list.files(".", pattern = paste0("^", gsub("\\.", "\\\\.", filename), "$"), recursive = TRUE, full.names = TRUE)
  if (length(hits) > 0) {
    return(normalizePath(hits[1], winslash = "/", mustWork = TRUE))
  }
  stop("Required file not found: ", default_path)
}

required_paths <- list(
  expr_clin = find_input_file("02_processed_data/tcga_expr_clin_matched.rds"),
  risk = find_input_file("07_signature/tcga_cross_platform_risk_score_by_sample.csv"),
  ssgsea = find_input_file("04_prft_score/tcga_ssgsea_scores.rds"),
  signature_key = find_input_file("09_enrichment/signature_key_gene_difference_high_vs_low_risk.csv"),
  gene_sets = find_input_file("03_gene_sets/prft_gene_sets_all.rds")
)

rank_score_matrix <- function(expr_mat, gene_sets) {
  sample_ids <- colnames(expr_mat)
  gene_ids <- rownames(expr_mat)

  rank_mat <- apply(expr_mat, 2, function(v) rank(v, ties.method = "average", na.last = "keep"))
  if (is.null(dim(rank_mat))) {
    rank_mat <- matrix(rank_mat, nrow = nrow(expr_mat), ncol = ncol(expr_mat))
    rownames(rank_mat) <- gene_ids
    colnames(rank_mat) <- sample_ids
  } else {
    rownames(rank_mat) <- gene_ids
    colnames(rank_mat) <- sample_ids
  }

  score_mat <- matrix(NA_real_, nrow = length(gene_sets), ncol = ncol(expr_mat), dimnames = list(names(gene_sets), sample_ids))
  coverage_rows <- vector("list", length(gene_sets))

  for (i in seq_along(gene_sets)) {
    gs_name <- names(gene_sets)[i]
    genes_total <- unique(gene_sets[[i]])
    genes_available <- intersect(genes_total, gene_ids)
    genes_missing <- setdiff(genes_total, gene_ids)

    if (length(genes_available) > 0) {
      mean_rank <- colMeans(rank_mat[genes_available, , drop = FALSE], na.rm = TRUE)
      score01 <- mean_rank / nrow(expr_mat)
      score_z <- {
        s <- sd(score01, na.rm = TRUE)
        m <- mean(score01, na.rm = TRUE)
        if (!is.finite(s) || s == 0) rep(0, length(score01)) else (score01 - m) / s
      }
      score_mat[gs_name, ] <- score_z
    }

    coverage_rows[[i]] <- data.table(
      immune_signature = gs_name,
      total_genes = length(genes_total),
      available_genes = length(genes_available),
      missing_genes = length(genes_missing),
      available_gene_list = paste(genes_available, collapse = ";"),
      missing_gene_list = paste(genes_missing, collapse = ";")
    )
  }

  list(score_mat = score_mat, coverage_dt = rbindlist(coverage_rows))
}

safe_cor_test <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(list(estimate = NA_real_, p.value = NA_real_))
  }
  out <- suppressWarnings(tryCatch(
    cor.test(x, y, method = "spearman", exact = FALSE),
    error = function(e) NULL
  ))
  if (is.null(out)) {
    return(list(estimate = NA_real_, p.value = NA_real_))
  }
  list(estimate = unname(out$estimate), p.value = out$p.value)
}

safe_wilcox <- function(x, group) {
  ok <- is.finite(x) & !is.na(group)
  x <- x[ok]
  group <- as.character(group[ok])
  if (length(unique(group)) != 2) {
    return(NULL)
  }
  suppressWarnings(tryCatch(
    wilcox.test(x ~ group),
    error = function(e) NULL
  ))
}

zscore_vector <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) rep(0, length(x)) else (x - m) / s
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

expr_obj <- readRDS(required_paths$expr_clin)
risk_dt <- as.data.table(read.csv(required_paths$risk, check.names = FALSE))
ssgsea_mat <- readRDS(required_paths$ssgsea)
signature_key_dt <- as.data.table(read.csv(required_paths$signature_key, check.names = FALSE))
gene_sets_all <- readRDS(required_paths$gene_sets)

expr_mat <- expr_obj$expr
if (!all(c("sample_id", "risk_score", "risk_group") %in% colnames(risk_dt))) {
  stop("Risk score file must contain sample_id, risk_score, and risk_group.")
}

risk_dt <- unique(risk_dt[, .(sample_id, patient_id, OS_time, OS_status, age, sex, FAB, WBC, risk_score, risk_group)], by = "sample_id")
sample_ids_use <- intersect(risk_dt$sample_id, colnames(expr_mat))
risk_dt <- risk_dt[sample_id %in% sample_ids_use]
setorder(risk_dt, sample_id)
expr_use <- expr_mat[, risk_dt$sample_id, drop = FALSE]

immune_gene_sets <- list(
  Myeloid_suppressive_extended = c("S100A8", "S100A9", "IL10", "TGFB1", "ARG1", "FCGR3A", "LILRB1", "LILRB2", "LILRB4", "CD14", "CD163", "MSR1", "MRC1", "CSF1R", "TYROBP", "ITGAM", "AIF1"),
  Monocyte_macrophage_like = c("CD14", "LST1", "FCGR3A", "CSF1R", "TYROBP", "AIF1", "MS4A7", "C1QA", "C1QB", "C1QC", "LILRB4", "SPI1"),
  M2_macrophage_like = c("CD163", "MRC1", "MSR1", "IL10", "TGFB1", "CCL18", "MAF", "VSIG4"),
  Neutrophil_inflammatory_like = c("S100A8", "S100A9", "CXCR2", "FCGR3B", "CSF3R", "LCN2", "MPO", "ELANE", "CEACAM8"),
  T_cell_exhaustion_extended = c("PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX", "ENTPD1", "CXCL13", "EOMES", "BATF"),
  Immune_checkpoint_extended = c("CD274", "PDCD1LG2", "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "IDO1", "LGALS9", "VSIR", "VTCN1", "CD80", "CD86"),
  IFN_gamma_response = c("IFNG", "STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11", "GBP1", "GBP5", "IDO1", "HLA-DRA", "HLA-DRB1"),
  Antigen_presentation = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "HLA-DRA", "HLA-DRB1", "CD74", "CIITA"),
  Cytotoxic_T_NK = c("CD8A", "CD8B", "GZMB", "GZMA", "PRF1", "NKG7", "GNLY", "KLRD1", "KLRK1"),
  PD1_PDL1_axis = c("CD274", "PDCD1", "PDCD1LG2", "JAK2", "STAT1", "IRF1", "IFNGR1", "IFNGR2")
)

immune_score_res <- rank_score_matrix(expr_use, immune_gene_sets)
immune_score_mat <- immune_score_res$score_mat
immune_cov_dt <- immune_score_res$coverage_dt

saveRDS(immune_score_mat, "10_immune/tcga_custom_immune_signature_scores.rds")
immune_score_dt <- as.data.table(immune_score_mat, keep.rownames = "immune_signature")
fwrite(immune_score_dt, "10_immune/tcga_custom_immune_signature_scores.csv")
fwrite(immune_cov_dt, "10_immune/tcga_immune_signature_gene_coverage.csv")

## A. risk_score vs immune signatures
immune_corr_rows <- lapply(rownames(immune_score_mat), function(sig) {
  ct <- safe_cor_test(risk_dt$risk_score, as.numeric(immune_score_mat[sig, risk_dt$sample_id]))
  cov_row <- immune_cov_dt[immune_signature == sig]
  data.table(
    immune_signature = sig,
    rho = ct$estimate,
    P.Value = ct$p.value,
    available_genes = cov_row$available_genes[1],
    missing_genes = cov_row$missing_genes[1]
  )
})
immune_corr_dt <- rbindlist(immune_corr_rows)
immune_corr_dt[, FDR := p.adjust(P.Value, method = "BH")]
immune_corr_dt[, direction := fifelse(is.na(rho), NA_character_, fifelse(rho > 0, "positive", fifelse(rho < 0, "negative", "zero")))]
immune_corr_dt <- immune_corr_dt[order(FDR, -abs(rho))]
fwrite(immune_corr_dt, "10_immune/risk_score_immune_signature_correlation.csv")

## B. high vs low immune signature differences
immune_diff_rows <- lapply(rownames(immune_score_mat), function(sig) {
  vals <- as.numeric(immune_score_mat[sig, risk_dt$sample_id])
  wt <- safe_wilcox(vals, risk_dt$risk_group)
  med_high <- stats::median(vals[risk_dt$risk_group == "high_risk"], na.rm = TRUE)
  med_low <- stats::median(vals[risk_dt$risk_group == "low_risk"], na.rm = TRUE)
  data.table(
    immune_signature = sig,
    median_high = med_high,
    median_low = med_low,
    difference = med_high - med_low,
    P.Value = if (is.null(wt)) NA_real_ else wt$p.value
  )
})
immune_diff_dt <- rbindlist(immune_diff_rows)
immune_diff_dt[, FDR := p.adjust(P.Value, method = "BH")]
immune_diff_dt[, direction := fifelse(is.na(difference), NA_character_, fifelse(difference > 0, "high_risk_higher", fifelse(difference < 0, "low_risk_higher", "no_difference")))]
immune_diff_dt <- immune_diff_dt[order(FDR, -abs(difference))]
fwrite(immune_diff_dt, "10_immune/high_low_risk_immune_signature_difference.csv")

compare_gene_group <- function(genes, out_path) {
  rows <- lapply(genes, function(g) {
    if (!(g %in% rownames(expr_use))) {
      return(data.table(
        gene_symbol = g,
        available_in_expr = FALSE,
        median_high = NA_real_,
        median_low = NA_real_,
        difference = NA_real_,
        P.Value = NA_real_
      ))
    }
    vals <- as.numeric(expr_use[g, risk_dt$sample_id])
    wt <- safe_wilcox(vals, risk_dt$risk_group)
    med_high <- stats::median(vals[risk_dt$risk_group == "high_risk"], na.rm = TRUE)
    med_low <- stats::median(vals[risk_dt$risk_group == "low_risk"], na.rm = TRUE)
    data.table(
      gene_symbol = g,
      available_in_expr = TRUE,
      median_high = med_high,
      median_low = med_low,
      difference = med_high - med_low,
      P.Value = if (is.null(wt)) NA_real_ else wt$p.value
    )
  })
  dt <- rbindlist(rows)
  dt[, FDR := p.adjust(P.Value, method = "BH")]
  dt[, direction := fifelse(is.na(difference), NA_character_, fifelse(difference > 0, "high_risk_higher", fifelse(difference < 0, "low_risk_higher", "no_difference")))]
  dt <- dt[order(FDR, -abs(difference))]
  fwrite(dt, out_path)
  dt
}

## C. checkpoint genes
checkpoint_genes <- c("CD274", "PDCD1LG2", "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "IDO1", "LGALS9", "VSIR", "VTCN1", "CD80", "CD86")
checkpoint_diff_dt <- compare_gene_group(checkpoint_genes, "10_immune/high_low_risk_checkpoint_gene_difference.csv")

## D. myeloid suppressive genes
myeloid_genes <- c("S100A8", "S100A9", "IL10", "TGFB1", "ARG1", "FCGR3A", "LILRB1", "LILRB2", "LILRB4", "CD14", "CD163", "MSR1", "MRC1", "CSF1R", "TYROBP", "ITGAM")
myeloid_diff_dt <- compare_gene_group(myeloid_genes, "10_immune/high_low_risk_myeloid_suppressive_gene_difference.csv")

## E. risk_score vs key immune genes
key_immune_genes <- c("CD274", "JAK2", "STAT1", "IRF1", "S100A8", "S100A9", "IL10", "TGFB1", "FCGR3A", "LILRB4", "CTLA4", "PDCD1", "LAG3", "HAVCR2", "TIGIT")
key_immune_corr_rows <- lapply(key_immune_genes, function(g) {
  if (!(g %in% rownames(expr_use))) {
    return(data.table(
      gene_symbol = g,
      available_in_expr = FALSE,
      rho = NA_real_,
      P.Value = NA_real_
    ))
  }
  ct <- safe_cor_test(risk_dt$risk_score, as.numeric(expr_use[g, risk_dt$sample_id]))
  data.table(
    gene_symbol = g,
    available_in_expr = TRUE,
    rho = ct$estimate,
    P.Value = ct$p.value
  )
})
key_immune_corr_dt <- rbindlist(key_immune_corr_rows)
key_immune_corr_dt[, FDR := p.adjust(P.Value, method = "BH")]
key_immune_corr_dt <- key_immune_corr_dt[order(FDR, -abs(rho))]
fwrite(key_immune_corr_dt, "10_immune/risk_score_key_immune_gene_correlation.csv")

## Summary
get_sig_value <- function(sig_name, col_name) {
  x <- immune_corr_dt[immune_signature == sig_name][[col_name]]
  if (length(x) == 0) NA else x[1]
}

summary_dt <- data.table(
  risk_score_myeloid_suppressive_extended_rho = get_sig_value("Myeloid_suppressive_extended", "rho"),
  risk_score_myeloid_suppressive_extended_P = get_sig_value("Myeloid_suppressive_extended", "P.Value"),
  risk_score_checkpoint_extended_rho = get_sig_value("Immune_checkpoint_extended", "rho"),
  risk_score_checkpoint_extended_P = get_sig_value("Immune_checkpoint_extended", "P.Value"),
  risk_score_PD1_PDL1_axis_rho = get_sig_value("PD1_PDL1_axis", "rho"),
  risk_score_PD1_PDL1_axis_P = get_sig_value("PD1_PDL1_axis", "P.Value"),
  risk_score_T_cell_exhaustion_extended_rho = get_sig_value("T_cell_exhaustion_extended", "rho"),
  risk_score_T_cell_exhaustion_extended_P = get_sig_value("T_cell_exhaustion_extended", "P.Value"),
  high_risk_significantly_higher_immune_signatures_FDR_less_0.05 = paste(immune_diff_dt[FDR < 0.05 & direction == "high_risk_higher"]$immune_signature, collapse = ";"),
  high_risk_significantly_higher_checkpoint_genes_FDR_less_0.05 = paste(checkpoint_diff_dt[FDR < 0.05 & direction == "high_risk_higher"]$gene_symbol, collapse = ";"),
  high_risk_significantly_higher_myeloid_suppressive_genes_FDR_less_0.05 = paste(myeloid_diff_dt[FDR < 0.05 & direction == "high_risk_higher"]$gene_symbol, collapse = ";"),
  bulk_RNAseq_interpretation_note = "Immune-related scores were interpreted as bulk transcriptional signatures rather than direct immune cell fractions."
)
fwrite(summary_dt, "14_tables/tcga_immune_microenvironment_summary.csv")

## Figures
plot_to_pdf_png(
  "13_figures/Figure9_risk_score_immune_signature_correlation",
  8, 5.5,
  function() {
    pdt <- copy(immune_corr_dt)
    pdt[, immune_signature := factor(pretty_label(immune_signature), levels = rev(pretty_label(immune_signature)))]
    p <- ggplot(pdt, aes(x = rho, y = immune_signature, color = FDR < 0.05)) +
      geom_segment(aes(x = 0, xend = rho, y = immune_signature, yend = immune_signature), linewidth = 0.8, color = "grey70") +
      geom_point(size = 3) +
      theme_bw(base_size = 12) +
      labs(x = "Spearman rho with risk score", y = NULL, title = "Risk score correlations with immune signatures")
    print(p)
  }
)
immune_plot_dt <- melt(
  as.data.table(t(immune_score_mat), keep.rownames = "sample_id")[risk_dt, on = "sample_id"],
  id.vars = c("sample_id", "patient_id", "OS_time", "OS_status", "age", "sex", "FAB", "WBC", "risk_score", "risk_group"),
  measure.vars = rownames(immune_score_mat),
  variable.name = "immune_signature",
  value.name = "score"
)
immune_plot_dt[, risk_group := pretty_factor(risk_group, levels = c("low_risk", "high_risk"))]
immune_plot_dt[, immune_signature := pretty_factor(immune_signature, levels = rownames(immune_score_mat))]
plot_to_pdf_png(
  "13_figures/Figure9_high_low_risk_immune_signature_scores",
  10, 8,
  function() {
    p <- ggplot(immune_plot_dt, aes(x = risk_group, y = score, fill = risk_group)) +
      geom_boxplot(outlier.size = 0.5) +
      facet_wrap(~ immune_signature, scales = "free_y", ncol = 3) +
      theme_bw(base_size = 11) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1)) +
      labs(x = NULL, y = "Immune signature score", title = "Immune signatures in high- vs low-risk groups")
    print(p)
  }
)

make_heatmap_dt <- function(genes) {
  available <- intersect(genes, rownames(expr_use))
  if (length(available) == 0) {
    return(data.table())
  }
  expr_sub <- expr_use[available, risk_dt$sample_id, drop = FALSE]
  expr_z <- t(apply(expr_sub, 1, zscore_vector))
  dt <- as.data.table(as.table(expr_z))
  colnames(dt) <- c("gene_symbol", "sample_id", "z_expr")
  dt <- merge(dt, risk_dt[, .(sample_id, risk_group, risk_score)], by = "sample_id", all.x = TRUE)
  sample_order <- risk_dt[order(risk_group, risk_score)]$sample_id
  dt[, sample_id := factor(sample_id, levels = sample_order)]
  dt
}

checkpoint_heat_dt <- make_heatmap_dt(checkpoint_genes)
plot_to_pdf_png(
  "13_figures/Figure9_checkpoint_gene_expression_high_vs_low_risk",
  11, 6,
  function() {
    p <- ggplot(checkpoint_heat_dt, aes(x = sample_id, y = gene_symbol, fill = z_expr)) +
      geom_tile() +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      labs(x = "TCGA samples ordered by risk group", y = NULL, fill = "z", title = "Checkpoint gene expression in high- vs low-risk groups")
    print(p)
  }
)

myeloid_heat_dt <- make_heatmap_dt(myeloid_genes)
plot_to_pdf_png(
  "13_figures/Figure9_myeloid_suppressive_gene_expression_high_vs_low_risk",
  11, 6,
  function() {
    p <- ggplot(myeloid_heat_dt, aes(x = sample_id, y = gene_symbol, fill = z_expr)) +
      geom_tile() +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      labs(x = "TCGA samples ordered by risk group", y = NULL, fill = "z", title = "Myeloid suppressive gene expression in high- vs low-risk groups")
    print(p)
  }
)

pd1_pdl1_axis_genes <- c("CD274", "PDCD1", "PDCD1LG2", "JAK2", "STAT1", "IRF1")
pd1_box_rows <- lapply(pd1_pdl1_axis_genes, function(g) {
  if (!(g %in% rownames(expr_use))) return(NULL)
  data.table(
    sample_id = risk_dt$sample_id,
    risk_group = risk_dt$risk_group,
    gene_symbol = g,
    expression = as.numeric(expr_use[g, risk_dt$sample_id])
  )
})
pd1_box_dt <- rbindlist(pd1_box_rows, fill = TRUE)
pd1_box_dt[, risk_group := pretty_factor(risk_group, levels = c("low_risk", "high_risk"))]
plot_to_pdf_png(
  "13_figures/Figure9_PD1_PDL1_axis_key_genes",
  9, 6,
  function() {
    p <- ggplot(pd1_box_dt, aes(x = risk_group, y = expression, fill = risk_group)) +
      geom_boxplot(outlier.size = 0.5) +
      facet_wrap(~ gene_symbol, scales = "free_y", ncol = 3) +
      theme_bw(base_size = 11) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1)) +
      labs(x = NULL, y = "Expression", title = "PD-1/PD-L1 axis key genes in high- vs low-risk groups")
    print(p)
  }
)

save_session_info("16_logs/sessionInfo_19_immune_microenvironment_analysis_tcga.txt")
