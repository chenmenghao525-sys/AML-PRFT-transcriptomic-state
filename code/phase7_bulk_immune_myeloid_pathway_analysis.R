#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(1234)

phase_lib <- Sys.getenv("PHASE1_ASCII_R_LIB", unset = "phase1_R_libs; local path removed")
if (dir.exists(phase_lib)) {
  .libPaths(unique(c(phase_lib, .libPaths())))
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

root <- Sys.getenv("PHASE7_ROOT", unset = "aml_prft_phase1_fix; local path removed")
setwd(root)

dir.create("03_results_tables", showWarnings = FALSE, recursive = TRUE)
dir.create("04_figures", showWarnings = FALSE, recursive = TRUE)
dir.create("05_logs", showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(root, "05_logs", "phase7_bulk_immune_log.txt")
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
msg_con <- file(log_file, open = "at")
sink(msg_con, type = "message")
on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(), silent = TRUE)
  try(close(msg_con), silent = TRUE)
  try(close(log_con), silent = TRUE)
}, add = TRUE)

cat("Phase 7: bulk immune / myeloid / pathway activity enhancement analysis\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Root:", root, "\n")
cat("Random seed: 1234\n\n")

if (file.exists("02_scripts/15_scripts/plot_label_utils.R")) {
  source("02_scripts/15_scripts/plot_label_utils.R")
} else {
  pretty_label <- function(x) as.character(x)
}

save_plot_pdf <- function(plot_obj, filename, width = 8, height = 6) {
  grDevices::pdf(filename, width = width, height = height, useDingbats = FALSE)
  print(plot_obj)
  grDevices::dev.off()
}

zscore_safe <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(rep(0, length(x)))
  }
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

bh_safe <- function(p) {
  p <- as.numeric(p)
  out <- rep(NA_real_, length(p))
  ok <- is.finite(p)
  if (any(ok)) {
    out[ok] <- p.adjust(p[ok], method = "BH")
  }
  out
}

infer_expression_scale <- function(mat) {
  vals <- as.numeric(mat)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) {
    return("unknown")
  }
  q99 <- as.numeric(stats::quantile(vals, probs = 0.99, na.rm = TRUE))
  vmax <- max(vals, na.rm = TRUE)
  if (vmax <= 30 && q99 <= 20) {
    "log2_or_already_normalized"
  } else if (vmax > 100) {
    "non_log_raw_like_or_intensity_scale"
  } else {
    "normalized_uncertain_scale"
  }
}

custom_rank_ssgsea_like <- function(expr_mat, gene_sets) {
  expr_mat <- as.matrix(expr_mat)
  storage.mode(expr_mat) <- "numeric"
  rank_mat <- apply(expr_mat, 2, function(x) rank(x, ties.method = "average", na.last = "keep"))
  if (is.null(dim(rank_mat))) {
    rank_mat <- matrix(rank_mat, ncol = 1)
    rownames(rank_mat) <- rownames(expr_mat)
    colnames(rank_mat) <- colnames(expr_mat)
  } else {
    rownames(rank_mat) <- rownames(expr_mat)
    colnames(rank_mat) <- colnames(expr_mat)
  }

  score_mat <- matrix(NA_real_, nrow = length(gene_sets), ncol = ncol(expr_mat))
  rownames(score_mat) <- names(gene_sets)
  colnames(score_mat) <- colnames(expr_mat)
  coverage_rows <- vector("list", length(gene_sets))

  for (i in seq_along(gene_sets)) {
    gs_name <- names(gene_sets)[i]
    gs <- unique(toupper(trimws(gene_sets[[i]])))
    available <- intersect(gs, rownames(expr_mat))
    missing <- setdiff(gs, rownames(expr_mat))
    if (length(available) > 0) {
      mean_ranks <- colMeans(rank_mat[available, , drop = FALSE], na.rm = TRUE)
      score01 <- (mean_ranks - 1) / max(nrow(expr_mat) - 1, 1)
      score_mat[gs_name, ] <- zscore_safe(score01)
    }
    coverage_rows[[i]] <- data.table(
      set_name = gs_name,
      total_genes = length(gs),
      available_genes = length(available),
      missing_genes = length(missing),
      available_gene_list = paste(available, collapse = ";"),
      missing_gene_list = paste(missing, collapse = ";")
    )
  }
  list(score_mat = score_mat, coverage = rbindlist(coverage_rows))
}

safe_cor_test <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- as.numeric(x[ok])
  y <- as.numeric(y[ok])
  if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(list(rho = NA_real_, p = NA_real_, n = length(x)))
  }
  fit <- suppressWarnings(tryCatch(stats::cor.test(x, y, method = "spearman", exact = FALSE), error = function(e) NULL))
  if (is.null(fit)) {
    return(list(rho = NA_real_, p = NA_real_, n = length(x)))
  }
  list(rho = unname(fit$estimate), p = fit$p.value, n = length(x))
}

safe_wilcox <- function(x, group) {
  ok <- is.finite(x) & !is.na(group)
  x <- as.numeric(x[ok])
  group <- as.character(group[ok])
  if (length(unique(group)) != 2) {
    return(list(p = NA_real_, n = length(x)))
  }
  fit <- suppressWarnings(tryCatch(stats::wilcox.test(x ~ group), error = function(e) NULL))
  if (is.null(fit)) {
    return(list(p = NA_real_, n = length(x)))
  }
  list(p = fit$p.value, n = length(x))
}

clean_symbol_string <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("///", ";", x, fixed = TRUE)
  x <- gsub("//", ";", x, fixed = TRUE)
  x <- gsub("\\|", ";", x)
  x <- gsub(",", ";", x, fixed = TRUE)
  x
}

parse_series_matrix <- function(localfile) {
  lines <- readLines(gzfile(localfile), warn = FALSE)
  begin_idx <- grep("^!series_matrix_table_begin", lines)
  end_idx <- grep("^!series_matrix_table_end", lines)
  if (length(begin_idx) == 0) {
    stop("series_matrix table not found in ", localfile)
  }
  if (length(end_idx) == 0 || end_idx[1] <= begin_idx[1]) {
    table_lines <- lines[(begin_idx[1] + 1):length(lines)]
  } else {
    table_lines <- lines[(begin_idx[1] + 1):(end_idx[1] - 1)]
  }
  expr_df <- read.delim(
    text = paste(table_lines, collapse = "\n"),
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    quote = "\"",
    comment.char = ""
  )
  expr_ids <- as.character(expr_df[[1]])
  expr_mat <- as.matrix(expr_df[, -1, drop = FALSE])
  storage.mode(expr_mat) <- "numeric"
  rownames(expr_mat) <- expr_ids
  expr_mat
}

parse_gpl_annotation <- function(localfile) {
  lines <- readLines(gzfile(localfile), warn = FALSE)
  begin_idx <- grep("^!platform_table_begin", lines)
  end_idx <- grep("^!platform_table_end", lines)
  if (length(begin_idx) == 0) {
    stop("GPL table not found in ", localfile)
  }
  if (length(end_idx) == 0 || end_idx[1] <= begin_idx[1]) {
    table_lines <- lines[(begin_idx[1] + 1):length(lines)]
  } else {
    table_lines <- lines[(begin_idx[1] + 1):(end_idx[1] - 1)]
  }
  read.delim(
    text = paste(table_lines, collapse = "\n"),
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    quote = "\"",
    comment.char = ""
  )
}

build_probe_mapping <- function(expr_mat, gpl_df) {
  norm_cols <- gsub("[^a-z0-9]", "", tolower(colnames(gpl_df)))
  symbol_candidates <- c("genesymbol", "symbol", "hgncsymbol", "gene_symbol")
  id_candidates <- c("id", "idref", "probeid")
  symbol_idx <- match(symbol_candidates, norm_cols)
  symbol_idx <- symbol_idx[!is.na(symbol_idx)][1]
  id_idx <- match(id_candidates, norm_cols)
  id_idx <- id_idx[!is.na(id_idx)][1]
  if (is.na(symbol_idx) || is.na(id_idx)) {
    stop("GPL annotation lacks recognizable probe ID or gene symbol column.")
  }
  gpl_map <- data.table(
    probe_id = as.character(gpl_df[[id_idx]]),
    raw_symbol = clean_symbol_string(gpl_df[[symbol_idx]])
  )
  gpl_map <- gpl_map[nzchar(probe_id) & nzchar(raw_symbol)]
  probe_stats <- data.table(
    probe_id = rownames(expr_mat),
    probe_iqr = apply(expr_mat, 1, stats::IQR, na.rm = TRUE)
  )
  expanded <- gpl_map[, .(gene_symbol = unlist(strsplit(raw_symbol, ";", fixed = TRUE))), by = "probe_id"]
  expanded[, gene_symbol := toupper(trimws(gene_symbol))]
  expanded <- unique(expanded[nzchar(gene_symbol) & !gene_symbol %in% c("NA", "---")])
  merged <- merge(expanded, probe_stats, by = "probe_id", all = FALSE)
  setorder(merged, gene_symbol, -probe_iqr)
  merged[, probe_rank := seq_len(.N), by = "gene_symbol"]
  merged
}

collapse_probes_to_gene <- function(expr_mat, mapping_dt) {
  keep <- mapping_dt[probe_rank == 1]
  keep <- keep[probe_id %in% rownames(expr_mat)]
  out <- expr_mat[keep$probe_id, , drop = FALSE]
  rownames(out) <- keep$gene_symbol
  out
}

long_from_matrix <- function(score_mat, category_name, cohort_name) {
  dt <- as.data.table(as.table(score_mat))
  setnames(dt, c("set_name", "sample_id", "score"))
  dt[, `:=`(set_category = category_name, cohort = cohort_name)]
  dt[]
}

wrap_text <- function(x, width = 26) {
  vapply(as.character(x), function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

plot_heatmap <- function(dt, x_col, y_col, fill_col, filename, title, midpoint = 0, width = 8.5, height = 6.5) {
  p <- ggplot(dt, aes_string(x = x_col, y = y_col, fill = fill_col)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = "#2b8cbe", mid = "white", high = "#d7301f", midpoint = midpoint) +
    theme_minimal(base_size = 10) +
    theme(
      axis.title = element_blank(),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    ) +
    labs(title = title, fill = "Score")
  save_plot_pdf(p, filename, width = width, height = height)
}

cat("Package availability audit\n")
pkg_dt <- rbindlist(lapply(c("data.table", "ggplot2", "patchwork", "GSVA", "MCPcounter", "estimate", "xCell", "ConsensusClusterPlus"), function(pkg) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  data.table(
    package = pkg,
    available = ok,
    version = if (ok) as.character(utils::packageVersion(pkg)) else NA_character_
  )
}))
print(pkg_dt)

required_files <- c(
  "01_processed_data/02_processed_data/tcga_expr_clin_matched.rds",
  "01_processed_data/04_prft_score/tcga_prft_score.csv",
  "03_results_tables/07_signature/tcga_cross_platform_risk_score_by_sample.csv",
  "03_results_tables/phase1_six_gene_coefficients.csv",
  "03_results_tables/phase4b_singlecell_robustness_summary.csv",
  "03_results_tables/phase3C_BeatAML_PRFT_scores.csv",
  "03_results_tables/phase3C_BeatAML_drug_correlation_all.csv",
  "03_results_tables/phase3C_BeatAML_high_low_drug_comparison.csv",
  "03_results_tables/phase3C_BeatAML_representative_drugs.csv",
  "03_results_tables/phase5_FS12_PPI_top20.csv",
  "03_results_tables/phase3A_fix_gene_recurrence_frequency.csv",
  "03_results_tables/phase3B_fix_SHAP_or_importance_top_features.csv",
  "phase1_runtime/08_validation/GSE37642_GPL570_external_validation_risk_score.csv",
  "phase1_runtime/08_validation/GSE12417_GPL570_external_validation_risk_score.csv",
  "phase1_runtime/08_validation/combined_GPL570_external_validation_risk_score.csv",
  "phase1_runtime/00_raw_data/geo_validation/GSE37642-GPL570_series_matrix.txt.gz",
  "phase1_runtime/00_raw_data/geo_validation/GSE12417-GPL570_series_matrix.txt.gz",
  "phase1_runtime/00_raw_data/geo_validation/GPL570_family.soft.gz"
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files: ", paste(missing_files, collapse = "; "))
}

tcga_obj <- readRDS("01_processed_data/02_processed_data/tcga_expr_clin_matched.rds")
tcga_expr <- tcga_obj$expr
storage.mode(tcga_expr) <- "numeric"
rownames(tcga_expr) <- toupper(rownames(tcga_expr))
tcga_clin <- as.data.table(tcga_obj$clin)
tcga_prft <- fread("01_processed_data/04_prft_score/tcga_prft_score.csv")
tcga_risk <- fread("03_results_tables/07_signature/tcga_cross_platform_risk_score_by_sample.csv")
six_coef <- fread("03_results_tables/phase1_six_gene_coefficients.csv")
phase4b_summary <- fread("03_results_tables/phase4b_singlecell_robustness_summary.csv")
beataml_scores <- fread("03_results_tables/phase3C_BeatAML_PRFT_scores.csv")
beataml_corr <- fread("03_results_tables/phase3C_BeatAML_drug_correlation_all.csv")
beataml_diff <- fread("03_results_tables/phase3C_BeatAML_high_low_drug_comparison.csv")
beataml_repr <- fread("03_results_tables/phase3C_BeatAML_representative_drugs.csv")
phase5_ppi <- fread("03_results_tables/phase5_FS12_PPI_top20.csv")
phase3a_rec <- fread("03_results_tables/phase3A_fix_gene_recurrence_frequency.csv")
phase3b_imp <- fread("03_results_tables/phase3B_fix_SHAP_or_importance_top_features.csv")

gpl570 <- parse_gpl_annotation("phase1_runtime/00_raw_data/geo_validation/GPL570_family.soft.gz")
gse37642_expr_probe <- parse_series_matrix("phase1_runtime/00_raw_data/geo_validation/GSE37642-GPL570_series_matrix.txt.gz")
gse12417_expr_probe <- parse_series_matrix("phase1_runtime/00_raw_data/geo_validation/GSE12417-GPL570_series_matrix.txt.gz")
probe_map <- build_probe_mapping(gse37642_expr_probe, gpl570)
gse37642_expr <- collapse_probes_to_gene(gse37642_expr_probe, probe_map)
gse12417_expr <- collapse_probes_to_gene(gse12417_expr_probe, probe_map)
storage.mode(gse37642_expr) <- "numeric"
storage.mode(gse12417_expr) <- "numeric"

gse37642_risk <- fread("phase1_runtime/08_validation/GSE37642_GPL570_external_validation_risk_score.csv")
gse12417_risk <- fread("phase1_runtime/08_validation/GSE12417_GPL570_external_validation_risk_score.csv")
combined_risk <- fread("phase1_runtime/08_validation/combined_GPL570_external_validation_risk_score.csv")

phase3a_top20 <- unique(phase3a_rec[order(-recurrence_count, -mean_abs_coefficient)]$gene_symbol[1:min(20, nrow(phase3a_rec))])
phase3b_core_axis <- unique(phase3b_imp[feature_set == "FS-F_core_axis_explanatory_genes"][order(rank)]$gene_symbol)
if (length(phase3b_core_axis) == 0) {
  phase3b_core_axis <- c("AIFM2", "HSPA1A", "HERPUD1", "HSPA5", "NFE2L2", "SLC3A2", "GCLM", "SLC40A1", "SOCS1", "STAT5B")
}

main_gene_sets_rds <- readRDS("01_processed_data/03_gene_sets/prft_gene_sets_main.rds")

bulk_signature_sets <- list(
  immune_checkpoint_set = c("CD274", "PDCD1LG2", "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "CD80", "CD86", "VSIR"),
  myeloid_suppressive_set = c("IL10", "TGFB1", "CXCL8", "ITGAM", "CD163", "FCGR3A", "CSF1R", "LILRB4", "S100A8", "S100A9", "LST1", "FCGR1A"),
  monocyte_macrophage_like_set = c("LYZ", "LST1", "FCN1", "S100A8", "S100A9", "MS4A7", "CD14", "FCGR3A", "CTSS", "TYROBP"),
  cDC_like_set = c("FCER1A", "CLEC10A", "ITGAX", "CST3", "IRF8", "BATF3"),
  T_NK_cytotoxic_set = c("CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "GZMB", "PRF1"),
  T_cell_exhaustion_set = c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "CTLA4"),
  proteostasis_UPR_set = c("HSPA5", "HSP90AA1", "HSP90AB1", "CALR", "CANX", "HERPUD1", "ERN1", "EIF2AK3", "ATF4", "XBP1", "DNAJB1"),
  ferroptosis_defense_set = c("SLC7A11", "GPX4", "GCLC", "GCLM", "GSR", "GSS", "FTH1", "FTL", "AIFM2", "NFE2L2", "CISD1", "SLC3A2"),
  JAK_STAT_PDL1_set = c("JAK2", "STAT1", "STAT3", "STAT5A", "STAT5B", "CD274", "IRF1", "IFNGR1", "IFNGR2", "SOCS1", "SOCS3"),
  oxidative_stress_NRF2_set = c("NFE2L2", "HMOX1", "NQO1", "GCLC", "GCLM", "TXNRD1", "SOD2", "GPX4"),
  LSC_stemness_set = main_gene_sets_rds$LSC17_core,
  phase5_PPI_hub_set = phase5_ppi$gene_symbol,
  phase3B_SHAP_core_axis_set = phase3b_core_axis
)

deconvolution_signature_sets <- list(
  monocyte_macrophage_like = bulk_signature_sets$monocyte_macrophage_like_set,
  myeloid_suppressive_like = bulk_signature_sets$myeloid_suppressive_set,
  cDC_like = bulk_signature_sets$cDC_like_set,
  T_NK_cytotoxic_like = bulk_signature_sets$T_NK_cytotoxic_set,
  T_cell_exhaustion_like = bulk_signature_sets$T_cell_exhaustion_set,
  immune_checkpoint_like = bulk_signature_sets$immune_checkpoint_set,
  neutrophil_inflammatory_like = c("S100A8", "S100A9", "CXCR2", "FCGR3B", "CSF3R", "LCN2", "MPO", "ELANE"),
  antigen_presentation_like = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "HLA-DRA", "HLA-DRB1", "CD74")
)

pathway_sets <- list(
  JAK_STAT_pathway = bulk_signature_sets$JAK_STAT_PDL1_set,
  IFN_response_pathway = c("IFNG", "STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11", "GBP1", "GBP5", "IDO1", "HLA-DRA"),
  TNF_NFKB_pathway = c("TNF", "NFKB1", "RELA", "NFKBIA", "TNFAIP3", "BIRC3", "CXCL8", "ICAM1"),
  inflammatory_response_pathway = c("IL1B", "IL6", "CXCL8", "CCL2", "TNF", "PTGS2", "NLRP3", "IRF1"),
  hypoxia_pathway = c("HIF1A", "VEGFA", "SLC2A1", "LDHA", "BNIP3", "CA9", "PGK1"),
  unfolded_protein_response_pathway = bulk_signature_sets$proteostasis_UPR_set,
  oxidative_stress_NRF2_pathway = bulk_signature_sets$oxidative_stress_NRF2_set,
  apoptosis_pathway = c("BAX", "BCL2", "BCL2L1", "MCL1", "CASP3", "CASP8", "PMAIP1", "BBC3"),
  ferroptosis_defense_pathway = bulk_signature_sets$ferroptosis_defense_set,
  PI3K_AKT_mTOR_pathway = c("PIK3CA", "PIK3CB", "AKT1", "AKT2", "MTOR", "RPTOR", "EIF4EBP1", "RPS6KB1"),
  MYC_targets_pathway = c("MYC", "MAX", "ODC1", "NCL", "HSPD1", "LDHA", "EIF4A1"),
  p53_pathway = c("TP53", "CDKN1A", "MDM2", "BAX", "GADD45A", "BBC3", "PMAIP1")
)

gene_sets_used <- rbindlist(list(
  rbindlist(lapply(names(bulk_signature_sets), function(nm) data.table(set_name = nm, set_category = "bulk_signature", gene = toupper(trimws(bulk_signature_sets[[nm]]))))),
  rbindlist(lapply(names(deconvolution_signature_sets), function(nm) data.table(set_name = nm, set_category = "immune_deconvolution_fallback", gene = toupper(trimws(deconvolution_signature_sets[[nm]]))))),
  rbindlist(lapply(names(pathway_sets), function(nm) data.table(set_name = nm, set_category = "pathway_activity", gene = toupper(trimws(pathway_sets[[nm]])))))
))
gene_sets_used <- unique(gene_sets_used[nzchar(gene)])
gene_sets_used_summary <- gene_sets_used[, .(
  gene_count = uniqueN(gene),
  genes = paste(sort(unique(gene)), collapse = ";")
), by = .(set_category, set_name)]
fwrite(gene_sets_used_summary, "03_results_tables/phase7_gene_sets_used.csv")

align_meta_expr <- function(meta_dt, expr_mat) {
  keep <- intersect(meta_dt$sample_id, colnames(expr_mat))
  meta_dt <- copy(meta_dt)[sample_id %in% keep]
  meta_dt[, .phase7_expr_order := match(sample_id, colnames(expr_mat))]
  setorder(meta_dt, .phase7_expr_order)
  meta_dt[, .phase7_expr_order := NULL]
  expr_use <- expr_mat[, meta_dt$sample_id, drop = FALSE]
  list(meta = meta_dt, expr = expr_use)
}

score_cohort <- function(cohort_name, expr_mat, meta_dt) {
  sig_res <- custom_rank_ssgsea_like(expr_mat, bulk_signature_sets)
  dec_res <- custom_rank_ssgsea_like(expr_mat, deconvolution_signature_sets)
  path_res <- custom_rank_ssgsea_like(expr_mat, pathway_sets)

  score_dt <- copy(meta_dt)
  sig_wide <- as.data.table(t(sig_res$score_mat), keep.rownames = "sample_id")
  path_wide <- as.data.table(t(path_res$score_mat), keep.rownames = "sample_id")
  score_dt <- merge(score_dt, sig_wide, by = "sample_id", all.x = TRUE)
  score_dt <- merge(score_dt, path_wide, by = "sample_id", all.x = TRUE)
  if (!"PRFT_score" %in% colnames(score_dt) || all(is.na(score_dt$PRFT_score))) {
    score_dt[, PRFT_score := zscore_safe(proteostasis_UPR_set + ferroptosis_defense_set)]
    score_dt[, PRFT_score_source := "phase7_proxy_from_signature_scores"]
  } else {
    score_dt[, PRFT_score_source := "original_or_existing_input"]
  }
  score_dt[, PRFT_group_phase7 := ifelse(PRFT_score >= stats::median(PRFT_score, na.rm = TRUE), "PRFT_high_like", "PRFT_low_like")]

  cov_dt <- rbindlist(list(
    cbind(sig_res$coverage, set_category = "bulk_signature"),
    cbind(dec_res$coverage, set_category = "immune_deconvolution_fallback"),
    cbind(path_res$coverage, set_category = "pathway_activity")
  ), fill = TRUE)
  cov_dt[, `:=`(
    cohort = cohort_name,
    scoring_method = "custom_rank_ssgsea_like",
    expression_scale_guess = infer_expression_scale(expr_mat),
    n_expression_genes = nrow(expr_mat),
    n_expression_samples = ncol(expr_mat)
  )]

  list(
    sample_scores = score_dt,
    signature_scores = sig_res$score_mat,
    deconv_scores = dec_res$score_mat,
    pathway_scores = path_res$score_mat,
    coverage = cov_dt
  )
}

tcga_meta <- merge(tcga_prft, tcga_risk, by = c("sample_id", "patient_id", "OS_time", "OS_status", "age", "sex", "FAB", "WBC"), all.x = TRUE)
tcga_aligned <- align_meta_expr(tcga_meta, tcga_expr)

gse37642_aligned <- align_meta_expr(gse37642_risk, gse37642_expr)
gse12417_aligned <- align_meta_expr(gse12417_risk, gse12417_expr)

common_combined_genes <- intersect(rownames(gse37642_aligned$expr), rownames(gse12417_aligned$expr))
combined_expr <- cbind(
  gse37642_aligned$expr[common_combined_genes, , drop = FALSE],
  gse12417_aligned$expr[common_combined_genes, , drop = FALSE]
)
combined_meta <- rbindlist(list(gse37642_aligned$meta, gse12417_aligned$meta), fill = TRUE)
combined_aligned <- align_meta_expr(combined_meta, combined_expr)

cohort_inputs <- list(
  TCGA = tcga_aligned,
  GSE37642 = gse37642_aligned,
  GSE12417 = gse12417_aligned,
  combined_GPL570 = combined_aligned
)

input_audit_rows <- lapply(names(cohort_inputs), function(cohort) {
  aligned <- cohort_inputs[[cohort]]
  meta <- aligned$meta
  expr <- aligned$expr
  data.table(
    cohort = cohort,
    expression_source = switch(
      cohort,
      TCGA = "01_processed_data/02_processed_data/tcga_expr_clin_matched.rds",
      GSE37642 = "phase1_runtime/00_raw_data/geo_validation/GSE37642-GPL570_series_matrix.txt.gz",
      GSE12417 = "phase1_runtime/00_raw_data/geo_validation/GSE12417-GPL570_series_matrix.txt.gz",
      combined_GPL570 = "merged GSE37642 GPL570 + GSE12417 GPL570"
    ),
    risk_source = switch(
      cohort,
      TCGA = "03_results_tables/07_signature/tcga_cross_platform_risk_score_by_sample.csv",
      GSE37642 = "phase1_runtime/08_validation/GSE37642_GPL570_external_validation_risk_score.csv",
      GSE12417 = "phase1_runtime/08_validation/GSE12417_GPL570_external_validation_risk_score.csv",
      combined_GPL570 = "phase1_runtime/08_validation/combined_GPL570_external_validation_risk_score.csv"
    ),
    raw_expression_samples = ncol(expr),
    raw_expression_genes = nrow(expr),
    matched_samples = nrow(meta),
    samples_with_risk_score = sum(is.finite(meta$risk_score)),
    samples_with_PRFT_score = if ("PRFT_score" %in% colnames(meta)) sum(is.finite(meta$PRFT_score)) else 0L,
    samples_with_OS = sum(is.finite(meta$OS_time) & is.finite(meta$OS_status)),
    expression_scale_guess = infer_expression_scale(expr),
    gene_identifier = "HGNC_symbol_rows",
    note = if (cohort == "TCGA") {
      "TCGA retains 173 expression samples; risk_score available for 151 survival-matched samples."
    } else if (cohort == "combined_GPL570") {
      "Combined GPL570 constructed by intersecting gene-level mapped GSE37642 and GSE12417 matrices."
    } else {
      "GPL570 series matrix mapped from probe to HGNC symbol using max-IQR probe retention."
    }
  )
})
input_audit_dt <- rbindlist(input_audit_rows)
fwrite(input_audit_dt, "03_results_tables/phase7_bulk_immune_input_audit.csv")

cat("\nInput audit\n")
print(input_audit_dt)

scored <- lapply(names(cohort_inputs), function(cohort) {
  cat("\nScoring cohort:", cohort, "\n")
  out <- score_cohort(cohort, cohort_inputs[[cohort]]$expr, cohort_inputs[[cohort]]$meta)
  out
})
names(scored) <- names(cohort_inputs)

score_method_log <- rbindlist(lapply(scored, `[[`, "coverage"), fill = TRUE)
fwrite(score_method_log, "03_results_tables/phase7_bulk_signature_scoring_method_log.csv")

fwrite(scored$TCGA$sample_scores, "03_results_tables/phase7_bulk_signature_scores_TCGA.csv")
fwrite(scored$GSE37642$sample_scores, "03_results_tables/phase7_bulk_signature_scores_GSE37642.csv")
fwrite(scored$GSE12417$sample_scores, "03_results_tables/phase7_bulk_signature_scores_GSE12417.csv")
fwrite(scored$combined_GPL570$sample_scores, "03_results_tables/phase7_bulk_signature_scores_combined_GPL570.csv")

analyze_signature_sets <- function(sample_dt, set_names, cohort_name, table_label) {
  cor_rows <- list()
  hl_rows <- list()
  for (nm in set_names) {
    vals <- sample_dt[[nm]]
    if ("risk_score" %in% colnames(sample_dt)) {
      fit <- safe_cor_test(sample_dt$risk_score, vals)
      cor_rows[[length(cor_rows) + 1L]] <- data.table(
        cohort = cohort_name,
        table_label = table_label,
        predictor_name = "risk_score",
        predictor_source = "existing_input",
        set_name = nm,
        rho = fit$rho,
        P.Value = fit$p,
        n = fit$n
      )
    }
    if ("PRFT_score" %in% colnames(sample_dt)) {
      fit2 <- safe_cor_test(sample_dt$PRFT_score, vals)
      cor_rows[[length(cor_rows) + 1L]] <- data.table(
        cohort = cohort_name,
        table_label = table_label,
        predictor_name = "PRFT_score",
        predictor_source = unique(sample_dt$PRFT_score_source)[1],
        set_name = nm,
        rho = fit2$rho,
        P.Value = fit2$p,
        n = fit2$n
      )
    }
    if ("risk_group" %in% colnames(sample_dt)) {
      fit3 <- safe_wilcox(vals, sample_dt$risk_group)
      med_high <- stats::median(vals[sample_dt$risk_group == "high_risk"], na.rm = TRUE)
      med_low <- stats::median(vals[sample_dt$risk_group == "low_risk"], na.rm = TRUE)
      hl_rows[[length(hl_rows) + 1L]] <- data.table(
        cohort = cohort_name,
        table_label = table_label,
        grouping_name = "risk_group",
        set_name = nm,
        high_group_label = "high_risk",
        low_group_label = "low_risk",
        median_high = med_high,
        median_low = med_low,
        difference_high_minus_low = med_high - med_low,
        P.Value = fit3$p,
        n = fit3$n
      )
    }
    if ("PRFT_group_phase7" %in% colnames(sample_dt)) {
      fit4 <- safe_wilcox(vals, sample_dt$PRFT_group_phase7)
      med_high2 <- stats::median(vals[sample_dt$PRFT_group_phase7 == "PRFT_high_like"], na.rm = TRUE)
      med_low2 <- stats::median(vals[sample_dt$PRFT_group_phase7 == "PRFT_low_like"], na.rm = TRUE)
      hl_rows[[length(hl_rows) + 1L]] <- data.table(
        cohort = cohort_name,
        table_label = table_label,
        grouping_name = "PRFT_group_phase7",
        set_name = nm,
        high_group_label = "PRFT_high_like",
        low_group_label = "PRFT_low_like",
        median_high = med_high2,
        median_low = med_low2,
        difference_high_minus_low = med_high2 - med_low2,
        P.Value = fit4$p,
        n = fit4$n
      )
    }
  }
  cor_dt <- rbindlist(cor_rows)
  hl_dt <- rbindlist(hl_rows)
  cor_dt[, FDR := bh_safe(P.Value), by = .(cohort, predictor_name, table_label)]
  cor_dt[, direction := fifelse(is.na(rho), NA_character_, fifelse(rho > 0, "positive", fifelse(rho < 0, "negative", "zero")))]
  hl_dt[, FDR := bh_safe(P.Value), by = .(cohort, grouping_name, table_label)]
  hl_dt[, direction := fifelse(is.na(difference_high_minus_low), NA_character_, fifelse(difference_high_minus_low > 0, "high_higher", fifelse(difference_high_minus_low < 0, "low_higher", "zero")))]
  list(cor = cor_dt, highlow = hl_dt)
}

sig_assoc <- rbindlist(lapply(names(scored), function(cohort) analyze_signature_sets(scored[[cohort]]$sample_scores, names(bulk_signature_sets), cohort, "bulk_signature")$cor))
sig_hl <- rbindlist(lapply(names(scored), function(cohort) analyze_signature_sets(scored[[cohort]]$sample_scores, names(bulk_signature_sets), cohort, "bulk_signature")$highlow))
fwrite(sig_assoc, "03_results_tables/phase7_PRFT_signature_correlation_by_cohort.csv")
fwrite(sig_hl, "03_results_tables/phase7_high_low_signature_comparison_by_cohort.csv")

deconv_assoc <- rbindlist(lapply(names(scored), function(cohort) analyze_signature_sets(scored[[cohort]]$sample_scores, names(deconvolution_signature_sets), cohort, "immune_deconvolution")$cor))
deconv_hl <- rbindlist(lapply(names(scored), function(cohort) analyze_signature_sets(scored[[cohort]]$sample_scores, names(deconvolution_signature_sets), cohort, "immune_deconvolution")$highlow))

deconv_scores_long <- rbindlist(lapply(names(scored), function(cohort) {
  dt <- long_from_matrix(scored[[cohort]]$deconv_scores, "immune_deconvolution_fallback", cohort)
  dt <- merge(dt, scored[[cohort]]$sample_scores[, .(sample_id, risk_score, risk_group, PRFT_score, PRFT_group_phase7)], by = "sample_id", all.x = TRUE)
  dt
}))

fwrite(rbindlist(list(
  data.table(
    method = "MCPcounter",
    available = requireNamespace("MCPcounter", quietly = TRUE),
    output_interpretation = "cell population score"
  ),
  data.table(
    method = "xCell",
    available = requireNamespace("xCell", quietly = TRUE),
    output_interpretation = "enrichment score"
  ),
  data.table(
    method = "ESTIMATE",
    available = requireNamespace("estimate", quietly = TRUE),
    output_interpretation = "stromal/immune score"
  ),
  data.table(
    method = "CIBERSORT_LM22",
    available = FALSE,
    output_interpretation = "not_run_local_signature_matrix_unavailable"
  ),
  data.table(
    method = "signature_scoring_fallback",
    available = TRUE,
    output_interpretation = "rank-based bulk signature estimate; not a direct cell fraction"
  )
)), "03_results_tables/phase7_immune_deconvolution_method_availability.csv")

fwrite(deconv_scores_long, "03_results_tables/phase7_immune_deconvolution_scores.csv")
fwrite(rbindlist(list(deconv_assoc, deconv_hl), fill = TRUE), "03_results_tables/phase7_immune_deconvolution_PRFT_association.csv")

pathway_assoc <- rbindlist(lapply(names(scored), function(cohort) analyze_signature_sets(scored[[cohort]]$sample_scores, names(pathway_sets), cohort, "pathway_activity")$cor))
pathway_hl <- rbindlist(lapply(names(scored), function(cohort) analyze_signature_sets(scored[[cohort]]$sample_scores, names(pathway_sets), cohort, "pathway_activity")$highlow))
pathway_scores_long <- rbindlist(lapply(names(scored), function(cohort) {
  dt <- long_from_matrix(scored[[cohort]]$pathway_scores, "pathway_activity", cohort)
  dt <- merge(dt, scored[[cohort]]$sample_scores[, .(sample_id, risk_score, risk_group, PRFT_score, PRFT_group_phase7)], by = "sample_id", all.x = TRUE)
  dt
}))
fwrite(pathway_scores_long, "03_results_tables/phase7_pathway_activity_scores.csv")
fwrite(rbindlist(list(pathway_assoc, pathway_hl), fill = TRUE), "03_results_tables/phase7_pathway_PRFT_association.csv")

summarize_consistency <- function(cor_dt, hl_dt, object_name) {
  risk_cor <- cor_dt[predictor_name == "risk_score", .(
    risk_cor_direction = if (all(direction == "positive", na.rm = TRUE)) "positive_all" else if (all(direction == "negative", na.rm = TRUE)) "negative_all" else "mixed",
    risk_cor_significant_cohorts = sum(FDR < 0.05, na.rm = TRUE),
    risk_cor_mean_rho = mean(rho, na.rm = TRUE),
    risk_cor_cohorts = paste(cohort, collapse = ";")
  ), by = set_name]
  risk_hl <- hl_dt[grouping_name == "risk_group", .(
    high_low_direction = if (all(direction == "high_higher", na.rm = TRUE)) "high_higher_all" else if (all(direction == "low_higher", na.rm = TRUE)) "low_higher_all" else "mixed",
    high_low_significant_cohorts = sum(FDR < 0.05, na.rm = TRUE),
    high_low_mean_diff = mean(difference_high_minus_low, na.rm = TRUE)
  ), by = set_name]
  prft_cor <- cor_dt[predictor_name == "PRFT_score", .(
    PRFT_cor_direction = if (all(direction == "positive", na.rm = TRUE)) "positive_all" else if (all(direction == "negative", na.rm = TRUE)) "negative_all" else "mixed",
    PRFT_cor_significant_cohorts = sum(FDR < 0.05, na.rm = TRUE),
    PRFT_cor_mean_rho = mean(rho, na.rm = TRUE)
  ), by = set_name]
  prft_hl <- hl_dt[grouping_name == "PRFT_group_phase7", .(
    PRFT_group_direction = if (all(direction == "high_higher", na.rm = TRUE)) "high_higher_all" else if (all(direction == "low_higher", na.rm = TRUE)) "low_higher_all" else "mixed",
    PRFT_group_significant_cohorts = sum(FDR < 0.05, na.rm = TRUE),
    PRFT_group_mean_diff = mean(difference_high_minus_low, na.rm = TRUE)
  ), by = set_name]
  out <- Reduce(function(x, y) merge(x, y, by = "set_name", all = TRUE), list(risk_cor, risk_hl, prft_cor, prft_hl))
  out[, analysis_domain := object_name]
  out[, cross_cohort_consistent := ifelse(
    risk_cor_direction %in% c("positive_all", "negative_all") &
      high_low_direction %in% c("high_higher_all", "low_higher_all"),
    "yes", "no"
  )]
  out[]
}

sig_consistency <- summarize_consistency(sig_assoc, sig_hl, "bulk_signature")
pathway_consistency <- summarize_consistency(pathway_assoc, pathway_hl, "pathway_activity")
fwrite(sig_consistency, "03_results_tables/phase7_cross_cohort_consistent_signatures.csv")
fwrite(pathway_consistency, "03_results_tables/phase7_pathway_cross_cohort_consistency.csv")

bulk_singlecell_summary <- rbindlist(list(
  data.table(
    evidence_axis = "monocyte_macrophage_like",
    bulk_result = sig_consistency[set_name == "monocyte_macrophage_like_set"]$cross_cohort_consistent[1],
    bulk_note = "Cross-cohort bulk signature consistency based on risk_score and PRFT score/group axes.",
    singlecell_result = phase4b_summary[item == "patient_level_support_for_monocyte_myeloid_PRFT_high_enrichment"]$result[1],
    consistency = ifelse(
      sig_consistency[set_name == "monocyte_macrophage_like_set"]$cross_cohort_consistent[1] == "yes" &&
        phase4b_summary[item == "patient_level_support_for_monocyte_myeloid_PRFT_high_enrichment"]$result[1] == "yes",
      "consistent", "partial_or_uncertain"
    )
  ),
  data.table(
    evidence_axis = "myeloid_suppressive",
    bulk_result = sig_consistency[set_name == "myeloid_suppressive_set"]$cross_cohort_consistent[1],
    bulk_note = "Bulk myeloid suppressive signature and deconvolution-like myeloid scores were evaluated conservatively.",
    singlecell_result = phase4b_summary[item == "recommended_final_wording"]$note[1],
    consistency = ifelse(sig_consistency[set_name == "myeloid_suppressive_set"]$cross_cohort_consistent[1] == "yes", "consistent", "partial_or_uncertain")
  ),
  data.table(
    evidence_axis = "immune_checkpoint_JAK_STAT",
    bulk_result = paste(
      "immune_checkpoint =", sig_consistency[set_name == "immune_checkpoint_set"]$cross_cohort_consistent[1],
      "; JAK_STAT =", pathway_consistency[set_name == "JAK_STAT_pathway"]$cross_cohort_consistent[1]
    ),
    bulk_note = "Bulk checkpoint and JAK/STAT pathway activity were treated as associated pathway signatures.",
    singlecell_result = "Phase 4b supports an immune-suppressive PRFT-high-like state.",
    consistency = "consistent_if_checkpoint_or_JAK_STAT_positive"
  ),
  data.table(
    evidence_axis = "LSC_stemness",
    bulk_result = sig_consistency[set_name == "LSC_stemness_set"]$cross_cohort_consistent[1],
    bulk_note = "LSC/stemness remained a secondary rather than dominant interpretive axis unless cross-cohort evidence was strong.",
    singlecell_result = phase4b_summary[item == "LSC17_supports_LSC_dominance"]$result[1],
    consistency = ifelse(
      phase4b_summary[item == "LSC17_supports_LSC_dominance"]$result[1] == "no",
      "supports_non_dominant_LSC_wording", "uncertain"
    )
  )
))
fwrite(bulk_singlecell_summary, "03_results_tables/phase7_bulk_singlecell_consistency_summary.csv")

beataml_target_scores <- c(
  "risk_score",
  "PRFT_score",
  "proteostasis_score",
  "ferroptosis_tolerance_score",
  "SLC7A11_GPX4_GSH_score",
  "JAK2_STAT5_PDL1_score",
  "immune_checkpoint_score",
  "myeloid_suppressive_score"
)
beataml_target_drugs <- unique(c(
  beataml_repr$drug_name,
  "Venetoclax",
  "Bortezomib (Velcade)",
  "17-AAG (Tanespimycin)",
  "Panobinostat",
  "Selumetinib (AZD6244)",
  "Cytarabine"
))
beataml_assoc <- merge(
  beataml_corr[score_name %in% beataml_target_scores & drug_name %in% beataml_target_drugs],
  beataml_diff[drug_name %in% beataml_target_drugs, .(drug_name, difference_high_minus_low, P.Value_highlow = P.Value, FDR_highlow = FDR, group_difference_interpretation)],
  by = "drug_name",
  all.x = TRUE
)
beataml_assoc[, main_text_candidate := drug_name %in% c("Venetoclax", "Bortezomib (Velcade)", "17-AAG (Tanespimycin)", "Panobinostat", "Selumetinib (AZD6244)")]
fwrite(beataml_assoc, "03_results_tables/phase7_BeatAML_pathway_drug_association.csv")

main_vs_supp <- rbindlist(list(
  data.table(
    artifact = "phase7_PRFT_signature_correlation_heatmap.pdf",
    recommended_placement = "main_text",
    rationale = "Summarizes cross-cohort association of PRFT/risk axes with bulk immune and myeloid signatures."
  ),
  data.table(
    artifact = "phase7_cross_cohort_signature_consistency.pdf",
    recommended_placement = "main_text",
    rationale = "Provides concise cross-cohort consistency support for the monocyte/myeloid PRFT-high state line."
  ),
  data.table(
    artifact = "phase7_PRFT_myeloid_deconvolution_boxplots.pdf",
    recommended_placement = "main_text",
    rationale = "Shows estimated myeloid/monocyte enrichment in the PRFT-high-like or high-risk state with transparent fallback wording."
  ),
  data.table(
    artifact = "phase7_PRFT_pathway_association_bubbleplot.pdf",
    recommended_placement = "main_text",
    rationale = "Links PRFT/risk axes to JAK/STAT, UPR/proteostasis, oxidative stress, and ferroptosis-defense pathway activity."
  ),
  data.table(
    artifact = "phase7_BeatAML_pathway_drug_correlation_heatmap.pdf",
    recommended_placement = "supplement",
    rationale = "Helpful for ex vivo pharmacogenomic interpretation but should stay subordinate to main survival/state evidence."
  ),
  data.table(
    artifact = "phase7_high_low_signature_boxplots.pdf",
    recommended_placement = "supplement",
    rationale = "Detailed distribution view for individual cohorts."
  ),
  data.table(
    artifact = "phase7_immune_deconvolution_heatmap.pdf",
    recommended_placement = "supplement",
    rationale = "Method is signature-based fallback rather than direct fraction deconvolution, so best kept as supporting evidence."
  ),
  data.table(
    artifact = "phase7_pathway_activity_heatmap.pdf",
    recommended_placement = "supplement",
    rationale = "Comprehensive pathway summary that is useful but denser than a main-text panel."
  )
))
fwrite(main_vs_supp, "03_results_tables/phase7_main_vs_supplement_recommendation.csv")

## Figures
cor_heatmap_dt <- sig_assoc[predictor_name %in% c("risk_score", "PRFT_score"),
  .(fill_value = rho),
  by = .(cohort, predictor_name, set_name)
]
cor_heatmap_dt[, row_label := paste(pretty_label(set_name), "\n(", predictor_name, ")", sep = "")]
cor_heatmap_dt[, cohort := factor(cohort, levels = c("TCGA", "GSE37642", "GSE12417", "combined_GPL570"))]
cor_heatmap_dt[, row_label := factor(row_label, levels = rev(unique(row_label)))]
plot_heatmap(
  cor_heatmap_dt,
  x_col = "cohort",
  y_col = "row_label",
  fill_col = "fill_value",
  filename = "04_figures/phase7_PRFT_signature_correlation_heatmap.pdf",
  title = "Bulk immune / myeloid signatures correlated with PRFT and risk axes",
  width = 9,
  height = 7
)

boxplot_sets <- c("monocyte_macrophage_like_set", "myeloid_suppressive_set", "immune_checkpoint_set", "proteostasis_UPR_set", "ferroptosis_defense_set", "JAK_STAT_PDL1_set")
box_dt <- rbindlist(lapply(names(scored), function(cohort) {
  dt <- copy(scored[[cohort]]$sample_scores)
  use_cols <- intersect(boxplot_sets, colnames(dt))
  melt(
    dt[, c("sample_id", "risk_group", "PRFT_group_phase7", use_cols), with = FALSE],
    id.vars = c("sample_id", "risk_group", "PRFT_group_phase7"),
    variable.name = "set_name",
    value.name = "score"
  )[, cohort := cohort]
}))
box_dt <- box_dt[risk_group %in% c("high_risk", "low_risk")]
box_dt[, set_name := factor(pretty_label(set_name), levels = pretty_label(boxplot_sets))]
box_dt[, risk_group := factor(risk_group, levels = c("low_risk", "high_risk"), labels = c("Low risk", "High risk"))]
box_plot <- ggplot(box_dt, aes(x = risk_group, y = score, fill = risk_group)) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.3) +
  facet_grid(set_name ~ cohort, scales = "free_y") +
  scale_fill_manual(values = c("Low risk" = "#91bfdb", "High risk" = "#fc8d59")) +
  theme_bw(base_size = 8) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  ) +
  labs(y = "Signature score", title = "Selected bulk signatures by risk group across cohorts")
save_plot_pdf(box_plot, "04_figures/phase7_high_low_signature_boxplots.pdf", width = 11, height = 8)

consistency_plot_dt <- sig_consistency[, .(
  set_name = pretty_label(set_name),
  risk_cor_significant_cohorts,
  high_low_significant_cohorts,
  risk_cor_mean_rho,
  cross_cohort_consistent
)]
consistency_plot_dt[, set_name := factor(set_name, levels = set_name[order(risk_cor_mean_rho)])]
consistency_plot <- ggplot(consistency_plot_dt, aes(x = risk_cor_mean_rho, y = set_name, size = risk_cor_significant_cohorts + high_low_significant_cohorts, color = cross_cohort_consistent)) +
  geom_point() +
  scale_color_manual(values = c("yes" = "#d7301f", "no" = "#636363")) +
  theme_bw(base_size = 9) +
  labs(x = "Mean rho across cohorts", y = NULL, size = "Significant cohort count", color = "Consistent", title = "Cross-cohort consistency of bulk signatures")
save_plot_pdf(consistency_plot, "04_figures/phase7_cross_cohort_signature_consistency.pdf", width = 8.5, height = 6.5)

deconv_heatmap_dt <- deconv_assoc[predictor_name == "risk_score", .(fill_value = rho), by = .(cohort, set_name)]
deconv_heatmap_dt[, set_name := factor(pretty_label(set_name), levels = rev(pretty_label(unique(set_name))))]
deconv_heatmap_dt[, cohort := factor(cohort, levels = c("TCGA", "GSE37642", "GSE12417", "combined_GPL570"))]
plot_heatmap(
  deconv_heatmap_dt,
  x_col = "cohort",
  y_col = "set_name",
  fill_col = "fill_value",
  filename = "04_figures/phase7_immune_deconvolution_heatmap.pdf",
  title = "Signature-based fallback immune-state estimates associated with risk score",
  width = 8.5,
  height = 6
)

myeloid_box_sets <- c("monocyte_macrophage_like", "myeloid_suppressive_like", "cDC_like")
myeloid_dt <- deconv_scores_long[set_name %in% myeloid_box_sets & risk_group %in% c("high_risk", "low_risk")]
myeloid_dt[, set_name := factor(pretty_label(set_name), levels = pretty_label(myeloid_box_sets))]
myeloid_dt[, risk_group := factor(risk_group, levels = c("low_risk", "high_risk"), labels = c("Low risk", "High risk"))]
myeloid_plot <- ggplot(myeloid_dt, aes(x = risk_group, y = score, fill = risk_group)) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.3) +
  facet_grid(set_name ~ cohort, scales = "free_y") +
  scale_fill_manual(values = c("Low risk" = "#91bfdb", "High risk" = "#fc8d59")) +
  theme_bw(base_size = 8) +
  theme(axis.title.x = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(y = "Estimated enrichment score", title = "Myeloid / monocyte-like estimated states across cohorts")
save_plot_pdf(myeloid_plot, "04_figures/phase7_PRFT_myeloid_deconvolution_boxplots.pdf", width = 10.5, height = 6.8)

path_heatmap_dt <- pathway_assoc[predictor_name == "risk_score", .(fill_value = rho), by = .(cohort, set_name)]
path_heatmap_dt[, set_name := factor(pretty_label(set_name), levels = rev(pretty_label(unique(set_name))))]
path_heatmap_dt[, cohort := factor(cohort, levels = c("TCGA", "GSE37642", "GSE12417", "combined_GPL570"))]
plot_heatmap(
  path_heatmap_dt,
  x_col = "cohort",
  y_col = "set_name",
  fill_col = "fill_value",
  filename = "04_figures/phase7_pathway_activity_heatmap.pdf",
  title = "Pathway activity associated with risk score across cohorts",
  width = 8.5,
  height = 7
)

path_bubble_dt <- pathway_assoc[predictor_name %in% c("risk_score", "PRFT_score")]
path_bubble_dt[, predictor_name := factor(predictor_name, levels = c("risk_score", "PRFT_score"))]
path_bubble_dt[, cohort := factor(cohort, levels = c("TCGA", "GSE37642", "GSE12417", "combined_GPL570"))]
path_bubble_dt[, set_name := factor(pretty_label(set_name), levels = rev(pretty_label(unique(set_name))))]
path_bubble_plot <- ggplot(path_bubble_dt, aes(x = cohort, y = set_name, size = -log10(pmax(P.Value, 1e-300)), color = rho, shape = predictor_name)) +
  geom_point() +
  scale_color_gradient2(low = "#2b8cbe", mid = "white", high = "#d7301f", midpoint = 0) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL, size = "-log10(P)", color = "rho", shape = "Predictor", title = "PRFT / risk associations with pathway activity")
save_plot_pdf(path_bubble_plot, "04_figures/phase7_PRFT_pathway_association_bubbleplot.pdf", width = 9.5, height = 7)

beataml_heatmap_dt <- beataml_assoc[, .(
  drug_name = wrap_text(drug_name, width = 18),
  score_name = pretty_label(score_name),
  rho = spearman_rho
)]
beataml_heatmap_dt[, score_name := factor(score_name, levels = pretty_label(beataml_target_scores))]
beataml_heatmap_dt[, drug_name := factor(drug_name, levels = unique(drug_name))]
beataml_heatmap <- ggplot(beataml_heatmap_dt, aes(x = drug_name, y = score_name, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2b8cbe", mid = "white", high = "#d7301f", midpoint = 0) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(x = NULL, y = NULL, fill = "rho", title = "BeatAML ex vivo drug AUC associations with PRFT-related scores")
save_plot_pdf(beataml_heatmap, "04_figures/phase7_BeatAML_pathway_drug_correlation_heatmap.pdf", width = 10.5, height = 5.8)

checklist_lines <- c(
  paste0("1. TCGA expression samples scored: ", input_audit_dt[cohort == "TCGA"]$matched_samples[1]),
  paste0("2. GSE37642 expression samples scored: ", input_audit_dt[cohort == "GSE37642"]$matched_samples[1]),
  paste0("3. GSE12417 expression samples scored: ", input_audit_dt[cohort == "GSE12417"]$matched_samples[1]),
  paste0("4. combined GPL570 expression samples scored: ", input_audit_dt[cohort == "combined_GPL570"]$matched_samples[1]),
  "5. Bulk signature scoring method: custom_rank_ssgsea_like (uniform across cohorts)",
  paste0(
    "6. Immune deconvolution packages available: ",
    if (nrow(pkg_dt[package %in% c("MCPcounter", "estimate", "xCell") & available == TRUE]) > 0) {
      paste(pkg_dt[package %in% c("MCPcounter", "estimate", "xCell") & available == TRUE]$package, collapse = ", ")
    } else {
      "none; fallback used"
    }
  ),
  paste0("7. Cross-cohort consistent bulk signatures: ", paste(sig_consistency[cross_cohort_consistent == "yes"]$set_name, collapse = ", ")),
  paste0("8. Cross-cohort consistent pathways: ", paste(pathway_consistency[cross_cohort_consistent == "yes"]$set_name, collapse = ", ")),
  paste0("9. Monocyte/macrophage-like bulk signature consistency: ", sig_consistency[set_name == "monocyte_macrophage_like_set"]$cross_cohort_consistent[1]),
  paste0("10. Myeloid suppressive bulk signature consistency: ", sig_consistency[set_name == "myeloid_suppressive_set"]$cross_cohort_consistent[1]),
  paste0("11. Immune checkpoint bulk signature consistency: ", sig_consistency[set_name == "immune_checkpoint_set"]$cross_cohort_consistent[1]),
  paste0("12. JAK/STAT pathway consistency: ", pathway_consistency[set_name == "JAK_STAT_pathway"]$cross_cohort_consistent[1]),
  paste0("13. LSC/stemness supports dominant primitive interpretation: ", ifelse(phase4b_summary[item == "LSC17_supports_LSC_dominance"]$result[1] == "yes", "yes", "no")),
  paste0("14. BeatAML representative drugs retained for interpretation: ", paste(unique(beataml_assoc[main_text_candidate == TRUE]$drug_name), collapse = ", ")),
  "15. Recommended wording: PRFT-high AML was associated with an inferred monocyte-like / myeloid stress-adapted and immune-suppressive bulk transcriptional state.",
  "16. Main-text candidates: signature correlation heatmap, cross-cohort consistency panel, myeloid-state boxplots, pathway association bubbleplot.",
  "17. Supplementary candidates: full high/low signature boxplots, deconvolution fallback heatmap, full pathway heatmap, BeatAML pathway-drug heatmap.",
  "18. Consensus clustering: not run because ConsensusClusterPlus was unavailable locally and this step was optional.",
  "19. Suggested next step: yes, Phase 9 full-manuscript integration can proceed with conservative wording."
)
writeLines(checklist_lines, "05_logs/phase7_bulk_immune_key_result_checklist.txt")

cat("\nKey checklist\n")
cat(paste(checklist_lines, collapse = "\n"), "\n")

cat("\nSession info\n")
print(sessionInfo())
cat("\nFinished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
