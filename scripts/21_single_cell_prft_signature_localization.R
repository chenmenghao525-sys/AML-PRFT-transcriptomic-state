#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
  library(irlba)
  library(patchwork)
  library(matrixStats)
  library(ggrepel)
})

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
plot_only_mode <- "--plot-only" %in% args

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

dir.create(file.path(project_root, "12_single_cell"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(project_root, "13_figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(project_root, "14_tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(project_root, "16_logs"), showWarnings = FALSE, recursive = TRUE)

label_file <- file.path(project_root, "15_scripts", "plot_label_utils.R")
if (file.exists(label_file)) {
  source(label_file)
} else {
  pretty_label <- function(x) {
    x <- as.character(x)
    gsub("_", " ", x, fixed = TRUE)
  }
}

sc_pretty_label <- function(x) {
  x <- as.character(x)
  label_map <- c(
    "risk_score" = "risk score",
    "PRFT_score" = "PRFT score",
    "high_risk" = "high risk",
    "low_risk" = "low risk",
    "high_risk_like" = "high-risk-like cells",
    "low_risk_like" = "low-risk-like cells",
    "Proteostasis_core" = "proteostasis score",
    "Proteostasis_core_score" = "proteostasis score",
    "Ferroptosis_tolerance_set" = "ferroptosis-tolerance score",
    "Ferroptosis_tolerance_set_score" = "ferroptosis-tolerance score",
    "SLC7A11_GPX4_GSH_axis" = "SLC7A11/GPX4-GSH axis score",
    "SLC7A11_GPX4_GSH_axis_score" = "SLC7A11/GPX4-GSH axis score",
    "JAK2_STAT5_PDL1_set" = "JAK2/STAT5/PD-L1 score",
    "JAK2_STAT5_PDL1_set_score" = "JAK2/STAT5/PD-L1 score",
    "Myeloid_suppressive_set" = "myeloid suppressive score",
    "Myeloid_suppressive_set_score" = "myeloid suppressive score",
    "Immune_checkpoint_set" = "immune checkpoint score",
    "Immune_checkpoint_set_score" = "immune checkpoint score",
    "Stemness_quiescence_set" = "stemness/quiescence score",
    "Stemness_quiescence_set_score" = "stemness/quiescence score",
    "LSC17_core" = "LSC17 score",
    "LSC17_core_score" = "LSC17 score",
    "Relapse_resistance_set" = "relapse-resistance score",
    "Relapse_resistance_set_score" = "relapse-resistance score"
  )
  out <- ifelse(x %in% names(label_map), label_map[x], x)
  out <- gsub("_", " ", out, fixed = TRUE)
  out
}

save_plot_dual <- function(plot_obj, file_base, width = 8, height = 6) {
  ggsave(paste0(file_base, ".pdf"), plot = plot_obj, width = width, height = height, units = "in")
  ggsave(paste0(file_base, ".png"), plot = plot_obj, width = width, height = height, units = "in", dpi = 300)
}

get_embedding_spec <- function(embedding_method) {
  if (identical(embedding_method, "PCA_fallback")) {
    return(list(
      x_col = "PC1",
      y_col = "PC2",
      x_lab = "PC1",
      y_lab = "PC2",
      title_prefix = "PCA projection",
      file_prefix_celltypes = "Figure11_scRNA_PCA_projection_celltypes",
      file_prefix_risk = "Figure11_scRNA_PCA_projection_risk_score",
      file_prefix_scores = "Figure11_scRNA_PCA_projection_PRFT_related_scores"
    ))
  }
  list(
    x_col = "UMAP_1",
    y_col = "UMAP_2",
    x_lab = "UMAP 1",
    y_lab = "UMAP 2",
    title_prefix = "UMAP",
    file_prefix_celltypes = "Figure11_scRNA_UMAP_celltypes",
    file_prefix_risk = "Figure11_scRNA_UMAP_risk_score",
    file_prefix_scores = "Figure11_scRNA_UMAP_PRFT_related_scores"
  )
}

ensure_embedding_coordinates <- function(cell_meta, embedding_method) {
  embed_spec <- get_embedding_spec(embedding_method)
  if (all(c(embed_spec$x_col, embed_spec$y_col) %in% names(cell_meta))) {
    return(cell_meta)
  }
  if (identical(embedding_method, "PCA_fallback")) {
    numeric_cols <- intersect(c(
      "risk_score",
      "PRFT_score",
      "Proteostasis_core_score",
      "Ferroptosis_tolerance_set_score",
      "SLC7A11_GPX4_GSH_axis_score",
      "JAK2_STAT5_PDL1_set_score",
      "Myeloid_suppressive_set_score",
      "Immune_checkpoint_set_score",
      "Stemness_quiescence_set_score",
      "LSC17_core_score",
      "Relapse_resistance_set_score",
      "primitive_annotation_score",
      "nFeature_RNA",
      "nCount_RNA",
      grep("^expr_", names(cell_meta), value = TRUE)
    ), names(cell_meta))
    if (length(numeric_cols) < 2) {
      stop("PCA_fallback plotting requires at least two numeric columns to reconstruct PC1/PC2.")
    }
    embed_mat <- as.matrix(cell_meta[, ..numeric_cols])
    storage.mode(embed_mat) <- "double"
    for (j in seq_len(ncol(embed_mat))) {
      col_j <- embed_mat[, j]
      if (all(!is.finite(col_j))) {
        embed_mat[, j] <- 0
      } else {
        col_j[!is.finite(col_j)] <- mean(col_j[is.finite(col_j)], na.rm = TRUE)
        embed_mat[, j] <- col_j
      }
    }
    pc_fit <- stats::prcomp(embed_mat, center = TRUE, scale. = TRUE)
    cell_meta[[embed_spec$x_col]] <- pc_fit$x[, 1]
    cell_meta[[embed_spec$y_col]] <- pc_fit$x[, 2]
    return(cell_meta)
  }
  stop("Embedding columns for plotting are missing and could not be reconstructed.")
}

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(..., collapse = " "), "\n")
}

zscore_vec <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(rep(0, length(x)))
  }
  (x - m) / s
}

generate_figure11_plots <- function(project_root, cell_meta, sig_expr_summary, embedding_method) {
  embed_spec <- get_embedding_spec(embedding_method)
  cell_meta <- ensure_embedding_coordinates(cell_meta, embedding_method)
  plot_dt <- as.data.table(cell_meta)
  plot_dt$cell_type <- factor(plot_dt$cell_type)

  p_celltypes <- ggplot(
    plot_dt,
    aes_string(x = embed_spec$x_col, y = embed_spec$y_col, color = "cell_type")
  ) +
    geom_point(size = 0.22, alpha = 0.7) +
    labs(
      title = paste(embed_spec$title_prefix, "by annotated cell type"),
      x = embed_spec$x_lab,
      y = embed_spec$y_lab,
      color = "cell type"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right"
    )
  save_plot_dual(p_celltypes, file.path(project_root, "13_figures", embed_spec$file_prefix_celltypes), 9.8, 7.2)

  p_risk <- ggplot(
    plot_dt,
    aes_string(x = embed_spec$x_col, y = embed_spec$y_col, color = "risk_score")
  ) +
    geom_point(size = 0.22, alpha = 0.75) +
    scale_color_gradient2(low = "#315B7D", mid = "white", high = "#B6423C") +
    labs(
      title = paste(embed_spec$title_prefix, "colored by risk score"),
      x = embed_spec$x_lab,
      y = embed_spec$y_lab,
      color = "risk score"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5))
  save_plot_dual(p_risk, file.path(project_root, "13_figures", embed_spec$file_prefix_risk), 8.2, 6.4)

  embed_score_panels <- c(
    "Proteostasis_core_score",
    "Ferroptosis_tolerance_set_score",
    "Myeloid_suppressive_set_score",
    "JAK2_STAT5_PDL1_set_score"
  )
  embed_score_panels <- intersect(embed_score_panels, names(plot_dt))
  embed_plots <- lapply(embed_score_panels, function(sc) {
    ggplot(
      plot_dt,
      aes_string(x = embed_spec$x_col, y = embed_spec$y_col, color = sc)
    ) +
      geom_point(size = 0.2, alpha = 0.75) +
      scale_color_gradient2(low = "#315B7D", mid = "white", high = "#B6423C") +
      labs(
        title = sc_pretty_label(sc),
        x = embed_spec$x_lab,
        y = embed_spec$y_lab,
        color = sc_pretty_label(sc)
      ) +
      theme_bw(base_size = 10) +
      theme(plot.title = element_text(size = 10, hjust = 0.5))
  })
  p_embed_combo <- wrap_plots(embed_plots, ncol = 2) +
    plot_annotation(title = paste(embed_spec$title_prefix, "of PRFT-related scores"))
  save_plot_dual(p_embed_combo, file.path(project_root, "13_figures", embed_spec$file_prefix_scores), 10.5, 8.5)

  celltype_order <- cell_meta[, .(median_risk = median(risk_score, na.rm = TRUE), n_cells = .N), by = .(cell_type)][n_cells >= 30][order(median_risk), cell_type]
  violin_dt <- plot_dt[cell_type %in% celltype_order]
  violin_dt$cell_type <- factor(violin_dt$cell_type, levels = celltype_order)
  p_violin <- ggplot(violin_dt, aes(x = cell_type, y = risk_score, fill = cell_type)) +
    geom_violin(scale = "width", trim = TRUE, color = NA, alpha = 0.85) +
    geom_boxplot(width = 0.12, outlier.size = 0.12, fill = "white") +
    labs(
      title = "risk score across annotated cell types",
      x = "cell type",
      y = "risk score"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  save_plot_dual(p_violin, file.path(project_root, "13_figures", "Figure11_scRNA_celltype_risk_score_violin"), 10.5, 6.2)

  heatmap_scores <- intersect(c(
    "risk_score",
    "PRFT_score",
    "Proteostasis_core_score",
    "Ferroptosis_tolerance_set_score",
    "SLC7A11_GPX4_GSH_axis_score",
    "Myeloid_suppressive_set_score",
    "Immune_checkpoint_set_score",
    "Stemness_quiescence_set_score",
    "LSC17_core_score",
    "Relapse_resistance_set_score",
    "JAK2_STAT5_PDL1_set_score"
  ), names(cell_meta))
  heatmap_mean <- cell_meta[, lapply(.SD, mean, na.rm = TRUE), by = .(cell_type), .SDcols = heatmap_scores]
  heatmap_long <- melt(as.data.table(heatmap_mean), id.vars = "cell_type", variable.name = "score_name", value.name = "mean_score")
  heatmap_long[, scaled_score := zscore_vec(mean_score), by = score_name]
  p_heat <- ggplot(heatmap_long, aes(x = sc_pretty_label(score_name), y = cell_type, fill = scaled_score)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "#315B7D", mid = "white", high = "#B6423C") +
    labs(
      title = "Cell-type summary of PRFT-related scores",
      x = "",
      y = "cell type",
      fill = "scaled mean"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  save_plot_dual(p_heat, file.path(project_root, "13_figures", "Figure11_scRNA_celltype_PRFT_score_heatmap"), 9.2, 7.8)

  dot_dt <- sig_expr_summary[cell_type %in% celltype_order]
  p_dot <- ggplot(dot_dt, aes(x = gene_symbol, y = factor(cell_type, levels = celltype_order), size = pct_detected, color = mean_expr)) +
    geom_point() +
    scale_color_gradient2(low = "#315B7D", mid = "white", high = "#B6423C") +
    labs(
      title = "6-gene signature localization across annotated cell types",
      x = "signature gene",
      y = "cell type",
      size = "fraction detected",
      color = "mean expression"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5))
  save_plot_dual(p_dot, file.path(project_root, "13_figures", "Figure11_scRNA_signature_gene_dotplot"), 8.8, 6.4)

  compare_scores <- intersect(c(
    "PRFT_score",
    "Proteostasis_core_score",
    "Ferroptosis_tolerance_set_score",
    "SLC7A11_GPX4_GSH_axis_score",
    "JAK2_STAT5_PDL1_set_score",
    "Myeloid_suppressive_set_score",
    "Immune_checkpoint_set_score",
    "Stemness_quiescence_set_score",
    "LSC17_core_score",
    "Relapse_resistance_set_score"
  ), names(cell_meta))
  high_low_dt <- cell_meta[risk_like_group %in% c("high_risk_like", "low_risk_like")]
  comp_plot_dt <- melt(
    high_low_dt[, c("risk_like_group", compare_scores), with = FALSE],
    id.vars = "risk_like_group",
    variable.name = "score_name",
    value.name = "score_value"
  )
  comp_plot_dt$risk_like_group <- factor(
    comp_plot_dt$risk_like_group,
    levels = c("low_risk_like", "high_risk_like"),
    labels = c("low-risk-like cells", "high-risk-like cells")
  )
  comp_plot_dt$score_name <- factor(comp_plot_dt$score_name, levels = compare_scores, labels = sc_pretty_label(compare_scores))
  p_comp <- ggplot(comp_plot_dt, aes(x = risk_like_group, y = score_value, fill = risk_like_group)) +
    geom_violin(scale = "width", trim = TRUE, color = NA, alpha = 0.85) +
    geom_boxplot(width = 0.12, outlier.size = 0.1, fill = "white") +
    facet_wrap(~ score_name, scales = "free_y", ncol = 3) +
    labs(
      title = "Program differences between high-risk-like and low-risk-like cells",
      x = "",
      y = "score",
      fill = "group"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
  save_plot_dual(p_comp, file.path(project_root, "13_figures", "Figure11_scRNA_high_low_risk_like_comparison"), 11.5, 8.8)
}

safe_mean <- function(x) {
  if (!length(x)) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

extract_sample_label <- function(path, suffix) {
  x <- basename(path)
  x <- sub("^GSM[0-9]+_", "", x)
  sub(suffix, "", x, fixed = TRUE)
}

parse_sample_state <- function(sample_id) {
  if (grepl("^BM", sample_id)) {
    return(data.frame(
      patient = sample_id,
      day = NA_real_,
      sample_state = "healthy_BM",
      stringsAsFactors = FALSE
    ))
  }
  patient <- sub("-D[0-9]+$", "", sample_id)
  day <- suppressWarnings(as.numeric(sub("^.*-D([0-9]+)$", "\\1", sample_id)))
  sample_state <- ifelse(!is.na(day) & day > 0, "treatment_persistent", "diagnostic_D0")
  data.frame(patient = patient, day = day, sample_state = sample_state, stringsAsFactors = FALSE)
}

detect_extracted_dir <- function(project_root) {
  candidates <- c(
    file.path(project_root, "12_single_cell", "GSE116256_extracted"),
    file.path(project_root, "00_raw_data", "GSE116256", "extracted"),
    file.path(project_root, "outputs", "aml_mrd_virtual_ko_drug_reversal", "runs", "strong_q2_v2_20260619", "data", "raw", "GSE116256", "extracted"),
    file.path(project_root, "outputs", "aml_mrd_virtual_ko_drug_reversal", "CellDeathDisease_submission_package_20260619", "bioinformatics_data", "raw", "GSE116256_RAW")
  )
  hits <- candidates[dir.exists(candidates)]
  if (length(hits)) {
    return(normalizePath(hits[[1]], winslash = "/", mustWork = TRUE))
  }
  NA_character_
}

extract_tar_if_needed <- function(project_root) {
  tar_candidates <- c(
    file.path(project_root, "FQD_bioinformatics_work", "data", "GSE116256_RAW.tar"),
    file.path(project_root, "outputs", "aml_mrd_virtual_ko_drug_reversal", "runs", "strong_q2_v2_20260619", "data", "raw", "GSE116256", "GSE116256_RAW.tar")
  )
  tar_file <- tar_candidates[file.exists(tar_candidates)][1]
  if (is.na(tar_file) || !nzchar(tar_file)) {
    return(NA_character_)
  }
  extract_dir <- file.path(project_root, "12_single_cell", "GSE116256_extracted")
  dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
  if (!length(list.files(extract_dir, pattern = "\\.dem\\.txt\\.gz$"))) {
    message_ts("Extracting local GSE116256 tar archive to", extract_dir)
    utils::untar(tar_file, exdir = extract_dir)
  }
  if (length(list.files(extract_dir, pattern = "\\.dem\\.txt\\.gz$"))) {
    return(normalizePath(extract_dir, winslash = "/", mustWork = TRUE))
  }
  NA_character_
}

rank_score_from_logexpr <- function(logexpr, gene_sets_upper) {
  if (!nrow(logexpr) || !ncol(logexpr)) {
    return(list(raw = data.frame(), coverage = data.frame()))
  }
  ranks <- matrixStats::colRanks(
    logexpr,
    ties.method = "average",
    preserveShape = TRUE
  )
  ranks <- (ranks - 1) / max(1, nrow(logexpr) - 1)
  storage.mode(ranks) <- "double"

  raw_scores <- list()
  coverage_rows <- list()
  for (nm in names(gene_sets_upper)) {
    genes_all <- unique(gene_sets_upper[[nm]])
    genes_avail <- intersect(genes_all, rownames(logexpr))
    genes_missing <- setdiff(genes_all, genes_avail)
    coverage_rows[[length(coverage_rows) + 1]] <- data.frame(
      feature_type = "gene_set",
      feature_name = nm,
      total_genes = length(genes_all),
      available_genes = length(genes_avail),
      missing_genes = length(genes_missing),
      available_gene_list = paste(genes_avail, collapse = ";"),
      missing_gene_list = paste(genes_missing, collapse = ";"),
      stringsAsFactors = FALSE
    )
    if (!length(genes_avail)) {
      raw_scores[[nm]] <- rep(NA_real_, ncol(logexpr))
    } else {
      raw_scores[[nm]] <- colMeans(ranks[genes_avail, , drop = FALSE], na.rm = TRUE)
    }
  }
  list(
    raw = as.data.frame(raw_scores, stringsAsFactors = FALSE),
    coverage = do.call(rbind, coverage_rows)
  )
}

read_sample_data <- function(dem_file, anno_file) {
  dem <- fread(dem_file)
  anno <- fread(anno_file)
  names(dem)[1] <- "Gene"
  list(dem = dem, anno = anno)
}

scale_rows <- function(mat) {
  mat <- as.matrix(mat)
  row_means <- rowMeans(mat, na.rm = TRUE)
  row_sds <- apply(mat, 1, stats::sd, na.rm = TRUE)
  row_sds[!is.finite(row_sds) | row_sds == 0] <- 1
  out <- (mat - row_means) / row_sds
  out[!is.finite(out)] <- 0
  out
}

if (plot_only_mode) {
  message_ts("Running in plot-only mode using cached single-cell outputs.")
  cell_scores_path <- file.path(project_root, "12_single_cell", "sc_cell_scores.csv")
  sig_expr_path <- file.path(project_root, "12_single_cell", "sc_signature_gene_expression_by_celltype.csv")
  summary_path <- file.path(project_root, "14_tables", "single_cell_prft_signature_summary.csv")
  if (!file.exists(cell_scores_path) || !file.exists(sig_expr_path) || !file.exists(summary_path)) {
    stop("Plot-only mode requires existing sc_cell_scores.csv, sc_signature_gene_expression_by_celltype.csv, and single_cell_prft_signature_summary.csv.")
  }
  cell_meta <- fread(cell_scores_path)
  sig_expr_summary <- fread(sig_expr_path)
  summary_tab_cached <- fread(summary_path)
  embedding_method <- if ("embedding_method" %in% names(summary_tab_cached)) summary_tab_cached$embedding_method[1] else "PCA_fallback"
  cell_meta <- as.data.table(cell_meta)
  if (!"cell_type" %in% names(cell_meta) && "CellType" %in% names(cell_meta)) {
    cell_meta[, cell_type := CellType]
  }
  generate_figure11_plots(project_root, cell_meta, sig_expr_summary, embedding_method)
  writeLines(capture.output(sessionInfo()), file.path(project_root, "16_logs", "sessionInfo_21_single_cell_prft_signature_localization.txt"))
  message_ts("Plot-only Figure11 regeneration completed.")
  message_ts("Embedding method:", embedding_method)
  quit(save = "no", status = 0)
}

signature_coef_path <- file.path(project_root, "00_LOCKED_FORMULA", "LOCKED_PRFT_six_gene_formula_A_coefficients.csv")
gene_set_path <- file.path(project_root, "03_gene_sets", "prft_gene_sets_all.rds")
if (!file.exists(signature_coef_path)) {
  stop("Missing signature coefficient file: ", signature_coef_path)
}
if (!file.exists(gene_set_path)) {
  stop("Missing gene set file: ", gene_set_path)
}

signature_coef <- fread(signature_coef_path)
setnames(signature_coef, tolower(names(signature_coef)))
if (!all(c("gene_symbol", "coefficient") %in% names(signature_coef))) {
  stop("Signature coefficient file must contain gene_symbol and coefficient columns.")
}
signature_coef$gene_symbol <- toupper(signature_coef$gene_symbol)

gene_sets <- readRDS(gene_set_path)
required_sets <- c(
  "Proteostasis_core",
  "Ferroptosis_tolerance_set",
  "SLC7A11_GPX4_GSH_axis",
  "JAK2_STAT5_PDL1_set",
  "Myeloid_suppressive_set",
  "Immune_checkpoint_set",
  "Stemness_quiescence_set",
  "LSC17_core",
  "Relapse_resistance_set"
)
missing_sets <- setdiff(required_sets, names(gene_sets))
if (length(missing_sets)) {
  stop("Missing required gene sets: ", paste(missing_sets, collapse = ", "))
}

gene_sets_use <- gene_sets[required_sets]
gene_sets_upper <- lapply(gene_sets_use, function(x) unique(toupper(x)))
signature_genes <- unique(signature_coef$gene_symbol)

extracted_dir <- detect_extracted_dir(project_root)
if (is.na(extracted_dir)) {
  extracted_dir <- extract_tar_if_needed(project_root)
}
if (is.na(extracted_dir)) {
  stop(
    "Could not locate extracted GSE116256 data. ",
    "Please provide a local directory containing matched .dem.txt.gz and .anno.txt.gz files."
  )
}
message_ts("Using local GSE116256 directory:", extracted_dir)

dem_files <- sort(list.files(extracted_dir, pattern = "\\.dem\\.txt\\.gz$", full.names = TRUE))
anno_files <- sort(list.files(extracted_dir, pattern = "\\.anno\\.txt\\.gz$", full.names = TRUE))
if (!length(dem_files) || !length(anno_files)) {
  stop("No dem/anno files found in: ", extracted_dir)
}

dem_samples <- vapply(dem_files, extract_sample_label, character(1), suffix = ".dem.txt.gz")
anno_samples <- vapply(anno_files, extract_sample_label, character(1), suffix = ".anno.txt.gz")
common_samples <- intersect(dem_samples, anno_samples)
common_samples <- common_samples[grepl("^(AML|BM)", common_samples)]
if (!length(common_samples)) {
  stop("No matched GSE116256 dem/anno sample pairs found.")
}

dem_files <- dem_files[match(common_samples, dem_samples)]
anno_files <- anno_files[match(common_samples, anno_samples)]
pair_manifest <- data.frame(
  sample_id = common_samples,
  dem_file = dem_files,
  anno_file = anno_files,
  stringsAsFactors = FALSE
)

message_ts("Matched GSE116256 sample pairs:", nrow(pair_manifest))

## Pass 1: estimate global gene variance for lightweight embedding.
global_sum <- NULL
global_sumsq <- NULL
global_n <- 0L
raw_cells_total <- 0L
qc_rows <- list()
gene_order_global <- NULL

for (i in seq_len(nrow(pair_manifest))) {
  sample_id <- pair_manifest$sample_id[i]
  message_ts("Pass 1 variable-gene scan:", sample_id, "(", i, "/", nrow(pair_manifest), ")")
  dat <- read_sample_data(pair_manifest$dem_file[i], pair_manifest$anno_file[i])
  dem <- dat$dem
  anno <- dat$anno
  genes <- toupper(dem$Gene)
  expr <- as.matrix(dem[, -1, with = FALSE])
  storage.mode(expr) <- "double"
  if (is.null(gene_order_global)) {
    gene_order_global <- genes
    global_sum <- setNames(numeric(length(genes)), genes)
    global_sumsq <- setNames(numeric(length(genes)), genes)
  }
  lib_size <- colSums(expr, na.rm = TRUE)
  mt_idx <- grepl("^MT-", genes)
  percent_mt <- if (any(mt_idx)) colSums(expr[mt_idx, , drop = FALSE], na.rm = TRUE) / pmax(lib_size, 1) * 100 else rep(0, ncol(expr))
  nfeature <- if ("NumberOfGenes" %in% names(anno)) anno$NumberOfGenes else colSums(expr > 0, na.rm = TRUE)
  if (!"Cell" %in% names(anno)) {
    stop("Annotation file missing Cell column: ", pair_manifest$anno_file[i])
  }
  cells_common <- intersect(colnames(expr), anno$Cell)
  expr <- expr[, cells_common, drop = FALSE]
  anno <- anno[match(cells_common, anno$Cell), ]
  lib_size <- lib_size[match(cells_common, colnames(dem)[-1])]
  percent_mt <- percent_mt[match(cells_common, colnames(dem)[-1])]
  nfeature <- nfeature[match(cells_common, anno$Cell)]
  keep <- !is.na(nfeature) & nfeature > 200 & nfeature < 6000 & !is.na(percent_mt) & percent_mt < 20
  raw_cells_total <- raw_cells_total + length(cells_common)
  qc_rows[[length(qc_rows) + 1]] <- data.frame(
    dataset = "GSE116256",
    sample_id = sample_id,
    cells_before_qc = length(cells_common),
    cells_after_qc = sum(keep),
    median_nFeature_RNA = stats::median(nfeature, na.rm = TRUE),
    median_nCount_RNA = stats::median(lib_size, na.rm = TRUE),
    median_percent_mt = stats::median(percent_mt, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  if (!any(keep)) {
    next
  }
  expr_keep <- expr[, keep, drop = FALSE]
  lib_keep <- colSums(expr_keep, na.rm = TRUE)
  logexpr <- log1p(t(t(expr_keep) / pmax(lib_keep, 1)) * 10000)
  rownames(logexpr) <- genes
  rs <- rowSums(logexpr, na.rm = TRUE)
  rs2 <- rowSums(logexpr^2, na.rm = TRUE)
  global_sum[names(rs)] <- global_sum[names(rs)] + rs
  global_sumsq[names(rs2)] <- global_sumsq[names(rs2)] + rs2
  global_n <- global_n + ncol(logexpr)
  rm(dat, dem, anno, expr, expr_keep, logexpr)
  gc(verbose = FALSE)
}

if (global_n == 0L) {
  stop("No cells passed QC in GSE116256.")
}

global_mean <- global_sum / global_n
global_var <- (global_sumsq / global_n) - (global_mean^2)
global_var[!is.finite(global_var)] <- 0
global_var <- sort(global_var, decreasing = TRUE)
top_variable_genes <- names(global_var)[seq_len(min(1200, length(global_var)))]
required_embedding_genes <- unique(c(unlist(gene_sets_upper, use.names = FALSE), signature_genes))
selected_features <- unique(c(top_variable_genes, required_embedding_genes))
message_ts("Selected embedding features:", length(selected_features))

## Pass 2: cell-level scoring and feature extraction.
cell_meta_list <- list()
feature_mat_list <- list()
coverage_list <- list()

for (i in seq_len(nrow(pair_manifest))) {
  sample_id <- pair_manifest$sample_id[i]
  message_ts("Pass 2 scoring:", sample_id, "(", i, "/", nrow(pair_manifest), ")")
  dat <- read_sample_data(pair_manifest$dem_file[i], pair_manifest$anno_file[i])
  dem <- dat$dem
  anno <- dat$anno
  genes <- toupper(dem$Gene)
  expr <- as.matrix(dem[, -1, with = FALSE])
  storage.mode(expr) <- "double"
  rownames(expr) <- genes
  if (anyDuplicated(rownames(expr))) {
    expr <- rowsum(expr, group = rownames(expr), reorder = FALSE)
  }
  if (!"Cell" %in% names(anno)) {
    stop("Annotation file missing Cell column: ", pair_manifest$anno_file[i])
  }
  cells_common <- intersect(colnames(expr), anno$Cell)
  expr <- expr[, cells_common, drop = FALSE]
  anno <- anno[match(cells_common, anno$Cell), ]
  lib_size <- colSums(expr, na.rm = TRUE)
  mt_idx <- grepl("^MT-", rownames(expr))
  percent_mt <- if (any(mt_idx)) colSums(expr[mt_idx, , drop = FALSE], na.rm = TRUE) / pmax(lib_size, 1) * 100 else rep(0, ncol(expr))
  nfeature <- if ("NumberOfGenes" %in% names(anno)) anno$NumberOfGenes else colSums(expr > 0, na.rm = TRUE)
  keep <- !is.na(nfeature) & nfeature > 200 & nfeature < 6000 & !is.na(percent_mt) & percent_mt < 20
  if (!any(keep)) {
    next
  }
  expr <- expr[, keep, drop = FALSE]
  anno <- anno[keep, ]
  lib_size <- lib_size[keep]
  percent_mt <- percent_mt[keep]
  nfeature <- nfeature[keep]
  logexpr <- log1p(t(t(expr) / pmax(lib_size, 1)) * 10000)

  score_res <- rank_score_from_logexpr(logexpr, gene_sets_upper)
  coverage_df <- score_res$coverage
  if (nrow(coverage_df)) {
    coverage_df$sample_id <- sample_id
    coverage_list[[length(coverage_list) + 1]] <- coverage_df
  }

  state_df <- parse_sample_state(sample_id)
  meta <- data.frame(
    Cell = colnames(logexpr),
    sample_id = sample_id,
    patient_id = state_df$patient,
    day = state_df$day,
    sample_state = state_df$sample_state,
    PredictionRefined = if ("PredictionRefined" %in% names(anno)) anno$PredictionRefined else NA_character_,
    PredictionRF2 = if ("PredictionRF2" %in% names(anno)) anno$PredictionRF2 else NA_character_,
    CellType = if ("CellType" %in% names(anno)) anno$CellType else NA_character_,
    nFeature_RNA = nfeature,
    nCount_RNA = lib_size,
    percent_mt = percent_mt,
    TranscriptomeUMIs = if ("TranscriptomeUMIs" %in% names(anno)) anno$TranscriptomeUMIs else NA_real_,
    Score_HSC = if ("Score_HSC" %in% names(anno)) anno$Score_HSC else NA_real_,
    Score_Prog = if ("Score_Prog" %in% names(anno)) anno$Score_Prog else NA_real_,
    Score_GMP = if ("Score_GMP" %in% names(anno)) anno$Score_GMP else NA_real_,
    stringsAsFactors = FALSE
  )
  meta$primitive_annotation_score <- pmax(meta$Score_HSC, meta$Score_Prog, meta$Score_GMP, na.rm = TRUE)
  meta$primitive_annotation_score[!is.finite(meta$primitive_annotation_score)] <- NA_real_

  score_df <- score_res$raw
  score_df <- setNames(score_df, paste0(names(score_df), "_raw"))
  meta <- cbind(meta, score_df)

  for (g in signature_genes) {
    meta[[paste0("expr_", g)]] <- if (g %in% rownames(logexpr)) logexpr[g, ] else NA_real_
  }

  feature_keep <- intersect(selected_features, rownames(logexpr))
  feat_mat <- logexpr[feature_keep, , drop = FALSE]
  if (length(feature_keep) < length(selected_features)) {
    missing_feat <- setdiff(selected_features, feature_keep)
    fill_mat <- matrix(0, nrow = length(missing_feat), ncol = ncol(logexpr), dimnames = list(missing_feat, colnames(logexpr)))
    feat_mat <- rbind(feat_mat, fill_mat)
  }
  feat_mat <- feat_mat[selected_features, , drop = FALSE]

  cell_meta_list[[length(cell_meta_list) + 1]] <- meta
  feature_mat_list[[length(feature_mat_list) + 1]] <- feat_mat

  rm(dat, dem, anno, expr, logexpr, meta, feat_mat, score_res)
  gc(verbose = FALSE)
}

cell_meta <- rbindlist(cell_meta_list, fill = TRUE)
feature_mat <- do.call(cbind, feature_mat_list)
feature_mat <- feature_mat[, cell_meta$Cell, drop = FALSE]

if (!nrow(cell_meta) || !ncol(feature_mat)) {
  stop("No single-cell data available after QC.")
}

coverage_combined <- rbindlist(coverage_list, fill = TRUE)
coverage_summary <- coverage_combined[, .(
  total_genes = max(total_genes, na.rm = TRUE),
  available_genes = max(available_genes, na.rm = TRUE),
  missing_genes = min(missing_genes, na.rm = TRUE),
  available_gene_list = available_gene_list[which.max(available_genes)],
  missing_gene_list = missing_gene_list[which.min(missing_genes)]
), by = .(feature_type, feature_name)]

signature_coverage <- lapply(signature_genes, function(g) {
  present <- g %in% rownames(feature_mat)
  data.frame(
    feature_type = "signature_gene",
    feature_name = g,
    total_genes = 1L,
    available_genes = as.integer(present),
    missing_genes = as.integer(!present),
    available_gene_list = if (present) g else "",
    missing_gene_list = if (present) "" else g,
    stringsAsFactors = FALSE
  )
})
coverage_summary <- rbind(coverage_summary, do.call(rbind, signature_coverage))
coverage_summary <- as.data.table(coverage_summary)
fwrite(coverage_summary, file.path(project_root, "12_single_cell", "sc_signature_gene_coverage.csv"))

score_cols_raw <- grep("_raw$", names(cell_meta), value = TRUE)
for (col in score_cols_raw) {
  new_col <- sub("_raw$", "_score", col)
  cell_meta[[new_col]] <- zscore_vec(cell_meta[[col]])
}

for (g in signature_genes) {
  expr_col <- paste0("expr_", g)
  z_col <- paste0("z_", g)
  cell_meta[[z_col]] <- zscore_vec(cell_meta[[expr_col]])
}

risk_components <- sapply(seq_len(nrow(signature_coef)), function(i) {
  gene <- signature_coef$gene_symbol[i]
  coef <- signature_coef$coefficient[i]
  z_col <- paste0("z_", gene)
  coef * cell_meta[[z_col]]
})
if (is.null(dim(risk_components))) {
  risk_components <- matrix(risk_components, ncol = 1)
}
cell_meta$risk_score <- rowSums(risk_components, na.rm = TRUE)

if (!all(c("Proteostasis_core_score", "Ferroptosis_tolerance_set_score") %in% names(cell_meta))) {
  stop("Missing proteostasis or ferroptosis-tolerance scores needed for PRFT score.")
}
cell_meta$PRFT_score <- (
  zscore_vec(cell_meta$Proteostasis_core_score) +
    zscore_vec(cell_meta$Ferroptosis_tolerance_set_score)
) / 2

q75 <- stats::quantile(cell_meta$risk_score, probs = 0.75, na.rm = TRUE)
q25 <- stats::quantile(cell_meta$risk_score, probs = 0.25, na.rm = TRUE)
cell_meta$risk_like_group <- "intermediate"
cell_meta$risk_like_group[cell_meta$risk_score >= q75] <- "high_risk_like"
cell_meta$risk_like_group[cell_meta$risk_score <= q25] <- "low_risk_like"

cell_meta$cell_type <- cell_meta$CellType
cell_meta$cell_type[is.na(cell_meta$cell_type) | cell_meta$cell_type == ""] <- "unannotated"
cell_meta$PredictionRefined[is.na(cell_meta$PredictionRefined) | cell_meta$PredictionRefined == ""] <- "unassigned"

## Embedding and lightweight clustering.
embedding_method <- "PCA_fallback"
feature_mat_scaled <- scale_rows(feature_mat)
n_pcs <- min(20L, nrow(feature_mat_scaled) - 1L, ncol(feature_mat_scaled) - 1L)
if (is.na(n_pcs) || n_pcs < 2L) {
  stop("Too few cells or features available for embedding.")
}
pc_fit <- irlba::prcomp_irlba(t(feature_mat_scaled), n = n_pcs, center = TRUE, scale. = FALSE)
pc_df <- as.data.frame(pc_fit$x[, seq_len(n_pcs), drop = FALSE])

if (requireNamespace("uwot", quietly = TRUE)) {
  umap_mat <- uwot::umap(
    pc_df[, 1:min(20, ncol(pc_df)), drop = FALSE],
    n_neighbors = 30,
    min_dist = 0.3,
    metric = "cosine",
    verbose = FALSE,
    ret_model = FALSE
  )
  colnames(umap_mat) <- c("UMAP_1", "UMAP_2")
  embedding_method <- "uwot_umap"
} else {
  umap_mat <- as.matrix(pc_df[, 1:2, drop = FALSE])
  colnames(umap_mat) <- c("PC1", "PC2")
}
cell_meta[, colnames(umap_mat) := as.data.table(umap_mat)]

set.seed(20260623)
pc_for_cluster <- pc_df[, 1:min(10, ncol(pc_df)), drop = FALSE]
cluster_k <- min(12, max(8, length(unique(cell_meta$cell_type))))
cluster_fit <- stats::kmeans(pc_for_cluster, centers = cluster_k, iter.max = 50, nstart = 10)
cell_meta$cluster <- paste0("Cluster ", cluster_fit$cluster)

score_keep <- c(
  "risk_score",
  "PRFT_score",
  "Proteostasis_core_score",
  "Ferroptosis_tolerance_set_score",
  "SLC7A11_GPX4_GSH_axis_score",
  "JAK2_STAT5_PDL1_set_score",
  "Myeloid_suppressive_set_score",
  "Immune_checkpoint_set_score",
  "Stemness_quiescence_set_score",
  "LSC17_core_score",
  "Relapse_resistance_set_score"
)
score_keep <- intersect(score_keep, names(cell_meta))

cell_scores_out <- cell_meta[, c(
  "Cell", "sample_id", "patient_id", "day", "sample_state", "PredictionRefined", "PredictionRF2",
  "cell_type", "cluster", "nFeature_RNA", "nCount_RNA", "percent_mt", "primitive_annotation_score",
  intersect(c("UMAP_1", "UMAP_2", "PC1", "PC2"), names(cell_meta)),
  score_keep, "risk_like_group", paste0("expr_", signature_genes)
), with = FALSE]
fwrite(cell_scores_out, file.path(project_root, "12_single_cell", "sc_cell_scores.csv"))

qc_summary <- rbindlist(qc_rows, fill = TRUE)
qc_summary <- rbind(
  qc_summary,
  data.frame(
    dataset = "GSE116256",
    sample_id = "ALL",
    cells_before_qc = sum(qc_summary$cells_before_qc, na.rm = TRUE),
    cells_after_qc = sum(qc_summary$cells_after_qc, na.rm = TRUE),
    median_nFeature_RNA = median(cell_meta$nFeature_RNA, na.rm = TRUE),
    median_nCount_RNA = median(cell_meta$nCount_RNA, na.rm = TRUE),
    median_percent_mt = median(cell_meta$percent_mt, na.rm = TRUE),
    stringsAsFactors = FALSE
  ),
  fill = TRUE
)
fwrite(qc_summary, file.path(project_root, "12_single_cell", "sc_qc_summary.csv"))

annotation_summary <- cell_meta[, .(
  n_cells = .N
), by = .(PredictionRefined, cell_type, sample_state)]
setorder(annotation_summary, -n_cells)
fwrite(annotation_summary, file.path(project_root, "12_single_cell", "sc_cell_annotation_summary.csv"))

cluster_summary_long <- rbindlist(lapply(score_keep, function(sc) {
  cell_meta[, .(
    n_cells = .N,
    mean_score = mean(get(sc), na.rm = TRUE),
    median_score = median(get(sc), na.rm = TRUE)
  ), by = .(cluster)][, score_name := sc]
}), fill = TRUE)
fwrite(cluster_summary_long, file.path(project_root, "12_single_cell", "sc_cluster_score_summary.csv"))

celltype_summary <- rbindlist(lapply(score_keep, function(sc) {
  cell_meta[, .(
    n_cells = .N,
    mean_score = mean(get(sc), na.rm = TRUE),
    median_score = median(get(sc), na.rm = TRUE)
  ), by = .(cell_type)][, score_name := sc]
}), fill = TRUE)
setorder(celltype_summary, score_name, -mean_score)
fwrite(celltype_summary, file.path(project_root, "12_single_cell", "sc_celltype_score_summary.csv"))

sig_expr_summary <- rbindlist(lapply(signature_genes, function(g) {
  expr_col <- paste0("expr_", g)
  cell_meta[, .(
    n_cells = .N,
    mean_expr = mean(get(expr_col), na.rm = TRUE),
    median_expr = median(get(expr_col), na.rm = TRUE),
    pct_detected = mean(get(expr_col) > 0, na.rm = TRUE)
  ), by = .(cell_type)][, gene_symbol := g]
}), fill = TRUE)
setorder(sig_expr_summary, gene_symbol, -mean_expr)
fwrite(sig_expr_summary, file.path(project_root, "12_single_cell", "sc_signature_gene_expression_by_celltype.csv"))

compare_scores <- intersect(score_keep, c(
  "PRFT_score",
  "Proteostasis_core_score",
  "Ferroptosis_tolerance_set_score",
  "SLC7A11_GPX4_GSH_axis_score",
  "JAK2_STAT5_PDL1_set_score",
  "Myeloid_suppressive_set_score",
  "Immune_checkpoint_set_score",
  "Stemness_quiescence_set_score",
  "LSC17_core_score",
  "Relapse_resistance_set_score"
))

high_low_dt <- cell_meta[risk_like_group %in% c("high_risk_like", "low_risk_like")]
high_low_comp <- rbindlist(lapply(compare_scores, function(sc) {
  x <- high_low_dt[risk_like_group == "high_risk_like", get(sc)]
  y <- high_low_dt[risk_like_group == "low_risk_like", get(sc)]
  wt <- tryCatch(wilcox.test(x, y), error = function(e) NULL)
  data.frame(
    score_name = sc,
    high_risk_like_mean = mean(x, na.rm = TRUE),
    low_risk_like_mean = mean(y, na.rm = TRUE),
    delta = mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE),
    P.Value = if (is.null(wt)) NA_real_ else wt$p.value,
    high_risk_like_cells = sum(high_low_dt$risk_like_group == "high_risk_like"),
    low_risk_like_cells = sum(high_low_dt$risk_like_group == "low_risk_like"),
    stringsAsFactors = FALSE
  )
}), fill = TRUE)
high_low_comp$FDR <- p.adjust(high_low_comp$P.Value, method = "BH")
setorder(high_low_comp, FDR, -delta)
fwrite(high_low_comp, file.path(project_root, "12_single_cell", "sc_high_low_risk_like_cell_comparison.csv"))

## Summary table.
top_group_for_score <- function(score_name, group_var = "cell_type") {
  tmp <- cell_meta[, .(mean_score = mean(get(score_name), na.rm = TRUE)), by = c(group_var)]
  tmp <- tmp[order(-mean_score)]
  if (!nrow(tmp)) return(NA_character_)
  paste0(tmp[[group_var]][1], " (mean=", sprintf("%.3f", tmp$mean_score[1]), ")")
}

high_programs <- high_low_comp[FDR < 0.05 & delta > 0, sc_pretty_label(score_name)]
if (!length(high_programs)) {
  high_programs <- "none at FDR < 0.05"
}

summary_tab <- data.frame(
  dataset_used = "GSE116256",
  cells_before_qc = raw_cells_total,
  cells_after_qc = nrow(cell_meta),
  clusters_detected = length(unique(cell_meta$cluster)),
  cell_types_detected = length(unique(cell_meta$cell_type)),
  signature_genes_available = sum(coverage_summary$feature_type == "signature_gene" & coverage_summary$available_genes > 0),
  risk_score_highest_celltype_or_cluster = top_group_for_score("risk_score", "cell_type"),
  proteostasis_highest_celltype_or_cluster = top_group_for_score("Proteostasis_core_score", "cell_type"),
  ferroptosis_tolerance_highest_celltype_or_cluster = top_group_for_score("Ferroptosis_tolerance_set_score", "cell_type"),
  myeloid_suppressive_highest_celltype_or_cluster = top_group_for_score("Myeloid_suppressive_set_score", "cell_type"),
  JAK2_STAT5_PDL1_highest_celltype_or_cluster = top_group_for_score("JAK2_STAT5_PDL1_set_score", "cell_type"),
  risk_score_cor_PRFT_score = suppressWarnings(cor(cell_meta$risk_score, cell_meta$PRFT_score, method = "spearman", use = "pairwise.complete.obs")),
  risk_score_cor_ferroptosis_tolerance_score = suppressWarnings(cor(cell_meta$risk_score, cell_meta$Ferroptosis_tolerance_set_score, method = "spearman", use = "pairwise.complete.obs")),
  risk_score_cor_myeloid_suppressive_score = suppressWarnings(cor(cell_meta$risk_score, cell_meta$Myeloid_suppressive_set_score, method = "spearman", use = "pairwise.complete.obs")),
  high_risk_like_cells_count = sum(cell_meta$risk_like_group == "high_risk_like"),
  low_risk_like_cells_count = sum(cell_meta$risk_like_group == "low_risk_like"),
  high_risk_like_enriched_programs = paste(high_programs, collapse = "; "),
  interpretation_note = "Single-cell analysis was used to localize PRFT-related transcriptional states and did not establish causal mechanisms.",
  embedding_method = embedding_method,
  stringsAsFactors = FALSE
)
fwrite(summary_tab, file.path(project_root, "14_tables", "single_cell_prft_signature_summary.csv"))

generate_figure11_plots(project_root, cell_meta, sig_expr_summary, embedding_method)

writeLines(capture.output(sessionInfo()), file.path(project_root, "16_logs", "sessionInfo_21_single_cell_prft_signature_localization.txt"))

message_ts("Single-cell analysis completed.")
message_ts("Cells after QC:", nrow(cell_meta))
message_ts("Embedding method:", embedding_method)
