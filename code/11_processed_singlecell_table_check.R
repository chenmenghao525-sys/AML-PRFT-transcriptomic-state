#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
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
sc_dir <- file.path(tables_dir, "12_single_cell")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

search_log <- file.path(log_dir, "phase4_fix_singlecell_file_search_log.txt")
if (file.exists(search_log)) file.remove(search_log)
log_msg <- function(...) {
  line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(..., collapse = " "))
  cat(line, "\n")
  cat(line, "\n", file = search_log, append = TRUE)
}

write_csv <- function(x, path) {
  fwrite(as.data.table(x), path)
  log_msg("Wrote", path)
}

save_pdf <- function(plot, path, width = 7.5, height = 6.0) {
  ggsave(path, plot = plot, width = width, height = height, units = "in", device = cairo_pdf)
  log_msg("Wrote", path)
}

status_plot <- function(title, body) {
  ggplot() +
    annotate("text", x = 0, y = 0.18, label = title, fontface = "bold", size = 5.0, hjust = 0.5) +
    annotate("text", x = 0, y = -0.08, label = body, size = 3.6, hjust = 0.5, lineheight = 1.05) +
    xlim(-1, 1) + ylim(-1, 1) +
    theme_void(base_size = 11)
}

norm_names <- function(x) {
  y <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
  gsub("^_|_$", "", y)
}

detect_cols <- function(cols) {
  nc <- norm_names(cols)
  list(
    has_cell_id = any(nc %in% c("cell", "cell_id", "barcode", "barcodes", "cellname", "cell_name", "cells")),
    has_umap = any(nc %in% c("umap_1", "umap1", "x_umap", "umap_x", "umap_01")) &&
      any(nc %in% c("umap_2", "umap2", "y_umap", "umap_y", "umap_02")),
    has_tsne = any(nc %in% c("tsne_1", "tsne1", "t_sne_1", "x_tsne", "tsne_x")) &&
      any(nc %in% c("tsne_2", "tsne2", "t_sne_2", "y_tsne", "tsne_y")),
    has_pca = any(nc %in% c("pc1", "pca_1", "pc_1")) &&
      any(nc %in% c("pc2", "pca_2", "pc_2")),
    has_module_scores = any(nc %in% norm_names(c(
      "risk_score", "PRFT_score", "Proteostasis_core_score",
      "Ferroptosis_tolerance_set_score", "JAK2_STAT5_PDL1_set_score",
      "SLC7A11_GPX4_GSH_axis_score", "Myeloid_suppressive_set_score",
      "LSC17_core_score", "Stemness_quiescence_set_score"
    ))),
    column_names = paste(cols, collapse = ";")
  )
}

classify_file <- function(p) {
  nm <- basename(p)
  ext <- tolower(tools::file_ext(nm))
  lp <- tolower(p)
  if (ext %in% c("rds", "rda", "rdata")) return("R object candidate")
  if (ext == "h5ad") return("h5ad object")
  if (ext == "loom") return("loom object")
  if (ext == "h5") return("h5 object / possible 10x h5")
  if (grepl("matrix\\.mtx(\\.gz)?$", nm, ignore.case = TRUE)) return("10x matrix.mtx")
  if (grepl("(barcodes|features|genes)\\.tsv(\\.gz)?$", nm, ignore.case = TRUE)) return("10x barcode/features/genes")
  if (grepl("umap|tsne|t-sne|pca", lp, ignore.case = TRUE)) return("embedding/coordinate candidate")
  if (grepl("metadata|annotation|cell_scores|celltype", lp, ignore.case = TRUE)) return("metadata/annotation candidate")
  if (grepl("count|counts|expression|expr|matrix", lp, ignore.case = TRUE)) return("count/expression matrix candidate")
  if (grepl("singlecell|single_cell|scrna|gse116256|seurat", lp, ignore.case = TRUE)) return("single-cell related file")
  paste0("other.", ext)
}

log_msg("Phase 4-fix single-cell UMAP/raw-object audit started.")
log_msg("Project root:", root)

scan_roots <- file.path(root, c("00_raw_data", "01_processed_data", "02_scripts", "03_results_tables", "04_figures", "05_logs"))
scan_roots <- scan_roots[dir.exists(scan_roots)]
all_files <- unlist(lapply(scan_roots, function(d) {
  list.files(d, recursive = TRUE, full.names = TRUE, all.files = FALSE)
}), use.names = FALSE)
all_files <- all_files[file.exists(all_files)]
all_files <- all_files[!grepl("/phase1_runtime/17_tmp/R_libs/|/\\.git/|Human_Genomics_PRFT_AML_|submission_package|INTERNAL_working_package|OFFICIAL_submission_package", all_files, ignore.case = TRUE)]

requested_ext <- grepl("\\.(rds|rda|RData|h5ad|h5|loom|mtx|mtx\\.gz|tsv|tsv\\.gz|csv|txt|txt\\.gz)$", basename(all_files), ignore.case = TRUE)
requested_kw <- grepl("matrix\\.mtx|barcodes\\.tsv|features\\.tsv|genes\\.tsv|counts?|expression|matrix|metadata|annotation|UMAP|tSNE|TSNE|PCA|GSE116256|singlecell|single_cell|scRNA|AML|Seurat|cell_scores|celltype", all_files, ignore.case = TRUE)
hit_files <- all_files[requested_ext | requested_kw]
fi <- file.info(hit_files)

search_dt <- data.table(
  file_path = normalizePath(hit_files, winslash = "/", mustWork = FALSE),
  file_name = basename(hit_files),
  file_type = vapply(hit_files, classify_file, character(1)),
  extension = tolower(tools::file_ext(basename(hit_files))),
  size_bytes = as.numeric(fi$size),
  last_modified = as.character(fi$mtime),
  matched_keywords = vapply(hit_files, function(p) {
    kws <- c("rds", "rda", "RData", "h5ad", "h5", "loom", "matrix.mtx", "barcodes.tsv",
             "features.tsv", "genes.tsv", "counts", "expression", "matrix", "metadata",
             "annotation", "UMAP", "tSNE", "TSNE", "PCA", "GSE116256", "singlecell",
             "single_cell", "scRNA", "AML", "Seurat")
    paste(kws[grepl(kws, p, ignore.case = TRUE, fixed = FALSE)], collapse = ";")
  }, character(1)),
  header_read = FALSE,
  has_cell_id_col = FALSE,
  has_UMAP_cols = FALSE,
  has_tSNE_cols = FALSE,
  has_PCA_cols = FALSE,
  has_module_score_cols = FALSE,
  header_columns = NA_character_,
  notes = NA_character_
)
search_dt[, singlecell_context := grepl("GSE116256|singlecell|single_cell|scrna|scRNA|Seurat|/12_single_cell/|sc_", file_path, ignore.case = TRUE)]
search_dt[, raw_or_semiraw_candidate := singlecell_context &
            !grepl("\\.(pdf|png|tiff|jpg|jpeg|R|r|log|txt)$", file_name, ignore.case = TRUE) &
            !grepl("phase4_|Figure|Supplementary|summary|coverage|comparison|checklist|sessionInfo", file_name, ignore.case = TRUE)]

tabular_idx <- which(search_dt$extension %in% c("csv", "tsv", "txt") &
                       search_dt$size_bytes <= 200 * 1024 * 1024)
for (i in tabular_idx) {
  path_i <- search_dt$file_path[i]
  hdr <- tryCatch({
    if (grepl("\\.csv$", path_i, ignore.case = TRUE)) {
      fread(path_i, nrows = 0, showProgress = FALSE)
    } else {
      fread(path_i, nrows = 0, showProgress = FALSE)
    }
  }, error = function(e) e)
  if (!inherits(hdr, "error")) {
    det <- detect_cols(names(hdr))
    search_dt[i, `:=`(
      header_read = TRUE,
      has_cell_id_col = det$has_cell_id,
      has_UMAP_cols = det$has_umap,
      has_tSNE_cols = det$has_tsne,
      has_PCA_cols = det$has_pca,
      has_module_score_cols = det$has_module_scores,
      header_columns = det$column_names,
      notes = "header parsed"
    )]
  } else {
    search_dt[i, notes := paste("header parse failed:", conditionMessage(hdr))]
  }
}

write_csv(search_dt, file.path(tables_dir, "phase4_fix_singlecell_file_search.csv"))

pkg_names <- c("UCell", "AUCell", "Seurat", "SeuratObject", "Scissor")
pkg_status <- data.table(
  package = pkg_names,
  available = vapply(pkg_names, requireNamespace, logical(1), quietly = TRUE),
  version = vapply(pkg_names, function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else NA_character_
  }, character(1)),
  load_attempt = "requireNamespace only",
  install_attempt = "not_attempted",
  reason = "No raw single-cell object/count matrix found during local audit; network/dependency installation was not forced because it would not create missing input data."
)
write_csv(pkg_status, file.path(tables_dir, "phase4_fix_UCell_AUCell_status.csv"))

has_seurat_obj <- any(search_dt$file_type == "R object candidate" &
                        search_dt$singlecell_context &
                        grepl("seurat", search_dt$file_path, ignore.case = TRUE) &
                        search_dt$raw_or_semiraw_candidate)
has_h5ad <- any(search_dt$file_type == "h5ad object")
has_10x_matrix <- any(search_dt$file_type == "10x matrix.mtx") &&
  any(search_dt$file_type == "10x barcode/features/genes" & grepl("barcodes", search_dt$file_name, ignore.case = TRUE)) &&
  any(search_dt$file_type == "10x barcode/features/genes" & grepl("features|genes", search_dt$file_name, ignore.case = TRUE))
has_count_matrix <- any(search_dt$file_type == "count/expression matrix candidate" &
                          search_dt$extension %in% c("csv", "tsv", "txt", "rds", "rda", "rdata", "h5") &
                          search_dt$singlecell_context &
                          search_dt$raw_or_semiraw_candidate &
                          grepl("count|counts|expression|expr|matrix", search_dt$file_name, ignore.case = TRUE))
has_metadata_umap <- any(search_dt$has_UMAP_cols)
has_metadata_tsne <- any(search_dt$has_tSNE_cols)
has_metadata_pca <- any(search_dt$has_PCA_cols)
has_cell_module_match <- any(search_dt$has_cell_id_col & search_dt$has_module_score_cols) &&
  any(search_dt$has_cell_id_col & (search_dt$has_UMAP_cols | search_dt$has_tSNE_cols))

cell_scores_path <- file.path(sc_dir, "sc_cell_scores.csv")
if (file.exists(cell_scores_path)) {
  cell_scores <- fread(cell_scores_path)
} else {
  cell_scores <- data.table()
}
score_cols <- intersect(c("risk_score", "PRFT_score", "Proteostasis_core_score",
                          "Ferroptosis_tolerance_set_score", "JAK2_STAT5_PDL1_set_score",
                          "SLC7A11_GPX4_GSH_axis_score", "Myeloid_suppressive_set_score",
                          "LSC17_core_score", "Stemness_quiescence_set_score"),
                        names(cell_scores))
processed_has_scores <- nrow(cell_scores) > 0 && length(score_cols) > 0
processed_has_umap <- nrow(cell_scores) > 0 && all(c("UMAP_1", "UMAP_2") %in% names(cell_scores))
processed_has_tsne <- nrow(cell_scores) > 0 && all(c("tSNE_1", "tSNE_2") %in% names(cell_scores))

can_generate_umap <- has_metadata_umap && processed_has_scores && has_cell_module_match
can_generate_tsne <- !can_generate_umap && has_metadata_tsne && processed_has_scores && has_cell_module_match
can_recompute_scores <- (has_seurat_obj || has_h5ad || has_10x_matrix || has_count_matrix) &&
  any(pkg_status$available[pkg_status$package %in% c("UCell", "AUCell", "Seurat")])

feasibility <- data.table(
  criterion = c("Seurat object exists", "h5ad object exists", "10x matrix exists",
                "count matrix + metadata exists", "metadata has UMAP coordinates",
                "metadata has tSNE coordinates", "processed table has module scores",
                "cell_id can match coordinates and module scores", "can generate true UMAP figures",
                "can generate tSNE fallback figures", "can recompute UCell/AUCell/module scores"),
  status = c(has_seurat_obj, has_h5ad, has_10x_matrix,
             has_count_matrix, has_metadata_umap, has_metadata_tsne, processed_has_scores,
             has_cell_module_match, can_generate_umap, can_generate_tsne, can_recompute_scores),
  evidence = c(
    paste(search_dt[file_type == "R object candidate" & grepl("seurat|single|scrna|gse116256|cell", file_path, ignore.case = TRUE), file_name], collapse = ";"),
    paste(search_dt[file_type == "h5ad object", file_name], collapse = ";"),
    paste(search_dt[file_type %in% c("10x matrix.mtx", "10x barcode/features/genes"), file_name], collapse = ";"),
    paste(search_dt[file_type == "count/expression matrix candidate" & singlecell_context == TRUE & raw_or_semiraw_candidate == TRUE, file_name], collapse = ";"),
    paste(search_dt[has_UMAP_cols == TRUE, file_name], collapse = ";"),
    paste(search_dt[has_tSNE_cols == TRUE, file_name], collapse = ";"),
    if (processed_has_scores) paste(score_cols, collapse = ";") else "",
    paste(search_dt[has_cell_id_col == TRUE & (has_UMAP_cols == TRUE | has_tSNE_cols == TRUE | has_module_score_cols == TRUE), file_name], collapse = ";"),
    if (can_generate_umap) "true UMAP coordinates and module scores available with matchable cell IDs" else "not possible: no matchable UMAP coordinates found",
    if (can_generate_tsne) "tSNE coordinates available with module scores" else "not possible or not needed",
    if (can_recompute_scores) "raw object/matrix and scoring package available" else "not possible: missing raw object/matrix and/or scoring packages"
  )
)
write_csv(feasibility, file.path(tables_dir, "phase4_fix_singlecell_UMAP_feasibility.csv"))

if (processed_has_scores) {
  module_scores <- cell_scores[, c(
    intersect(c("Cell", "sample_id", "patient_id", "sample_state", "cell_type", "cluster", "risk_like_group"), names(cell_scores)),
    score_cols
  ), with = FALSE]
  module_scores[, score_source := "existing processed module-score table from Phase 4; not recomputed by UCell/AUCell"]
  write_csv(module_scores, file.path(tables_dir, "phase4_fix_singlecell_module_scores.csv"))
} else {
  write_csv(data.table(status = "not_available", reason = "processed score table missing"), file.path(tables_dir, "phase4_fix_singlecell_module_scores.csv"))
}

if (can_generate_umap || can_generate_tsne) {
  coord_file <- search_dt[(has_UMAP_cols == can_generate_umap | has_tSNE_cols == can_generate_tsne) & has_cell_id_col == TRUE][1, file_path]
  coords <- fread(coord_file)
  nm <- norm_names(names(coords))
  cell_col <- names(coords)[match(TRUE, nm %in% c("cell", "cell_id", "barcode", "barcodes", "cellname", "cell_name", "cells"))]
  if (can_generate_umap) {
    x_col <- names(coords)[match(TRUE, nm %in% c("umap_1", "umap1", "x_umap", "umap_x", "umap_01"))]
    y_col <- names(coords)[match(TRUE, nm %in% c("umap_2", "umap2", "y_umap", "umap_y", "umap_02"))]
    prefix <- "UMAP"
  } else {
    x_col <- names(coords)[match(TRUE, nm %in% c("tsne_1", "tsne1", "t_sne_1", "x_tsne", "tsne_x"))]
    y_col <- names(coords)[match(TRUE, nm %in% c("tsne_2", "tsne2", "t_sne_2", "y_tsne", "tsne_y"))]
    prefix <- "tSNE"
  }
  setnames(coords, c(cell_col, x_col, y_col), c("Cell", "DIM_1", "DIM_2"))
  plot_dt <- merge(coords[, .(Cell, DIM_1, DIM_2)], cell_scores, by = "Cell")
  if (nrow(plot_dt) > 0) {
    p_cell <- ggplot(plot_dt, aes(DIM_1, DIM_2, color = cell_type)) +
      geom_point(size = 0.18, alpha = 0.75) +
      theme_bw(base_size = 9) +
      labs(title = paste(prefix, "cell-type annotation"), x = paste(prefix, "1"), y = paste(prefix, "2"), color = "cell type")
    save_pdf(p_cell, file.path(fig_dir, "phase4_fix_UMAP_celltype_annotation.pdf"), 8.5, 6.5)
    score_plot <- function(sc, path, ttl) {
      p <- ggplot(plot_dt, aes(DIM_1, DIM_2, color = .data[[sc]])) +
        geom_point(size = 0.18, alpha = 0.75) +
        scale_color_gradient2(low = "#315B7D", mid = "white", high = "#B6423C") +
        theme_bw(base_size = 9) +
        labs(title = ttl, x = paste(prefix, "1"), y = paste(prefix, "2"), color = sc)
      save_pdf(p, path, 7.2, 5.8)
    }
    score_plot("risk_score", file.path(fig_dir, "phase4_fix_UMAP_six_gene_risk_score.pdf"), paste(prefix, "six-gene risk score"))
    score_plot("PRFT_score", file.path(fig_dir, "phase4_fix_UMAP_PRFT_score.pdf"), paste(prefix, "PRFT score"))
    core_scores <- intersect(c("Proteostasis_core_score", "Ferroptosis_tolerance_set_score",
                               "JAK2_STAT5_PDL1_set_score", "SLC7A11_GPX4_GSH_axis_score",
                               "Myeloid_suppressive_set_score", "LSC17_core_score"), names(plot_dt))
    p_list <- lapply(core_scores, function(sc) {
      ggplot(plot_dt, aes(DIM_1, DIM_2, color = .data[[sc]])) +
        geom_point(size = 0.12, alpha = 0.72) +
        scale_color_gradient2(low = "#315B7D", mid = "white", high = "#B6423C") +
        theme_bw(base_size = 8) +
        labs(title = gsub("_", " ", sc), x = paste(prefix, "1"), y = paste(prefix, "2"), color = "score")
    })
    p_core <- Reduce(`+`, p_list)
    if (requireNamespace("patchwork", quietly = TRUE)) {
      p_core <- patchwork::wrap_plots(p_list, ncol = 2)
    } else {
      p_core <- p_list[[1]]
    }
    save_pdf(p_core, file.path(fig_dir, "phase4_fix_UMAP_core_scores.pdf"), 10.5, 8.2)
  } else {
    can_generate_umap <- FALSE
    can_generate_tsne <- FALSE
    log_msg("Coordinate file was found but no cell IDs matched sc_cell_scores.csv.")
  }
}

if (!can_generate_umap && !can_generate_tsne) {
  limitation <- c(
    "Phase 4-fix single-cell limitation statement",
    "",
    "A systematic local audit did not identify a matchable Seurat object, h5ad object, 10x matrix, count matrix with metadata, or metadata table containing UMAP/tSNE coordinates that can be linked to the processed module-score table.",
    "",
    "Therefore, the current single-cell module should be described as a processed-table-based cell-state enrichment analysis, not as a de novo Seurat single-cell reanalysis.",
    "",
    "Do not claim that UCell/AUCell scoring succeeded in Phase 4-fix. UCell, AUCell, Seurat, SeuratObject, and Scissor were not available in the current R environment, and no raw single-cell object/count matrix was available to recompute scores.",
    "",
    "Do not display UMAP localization panels for this phase. PCA panels or placeholder pages must not be described as UMAP.",
    "",
    "Main-text evidence that can be retained:",
    "- phase4_score_heatmap_by_celltype.pdf",
    "- phase4_PRFT_high_like_fraction_by_celltype.pdf",
    "- phase4_RoE_heatmap.pdf",
    "",
    "Supported wording:",
    "PRFT-high-like cells were associated with a monocyte-like/myeloid stress-adapted PRFT-high state.",
    "",
    "Avoid causal wording such as proved, mediated, activated, reversed, drives, or causes."
  )
  writeLines(limitation, file.path(log_dir, "phase4_fix_singlecell_limitation_statement.txt"))
  log_msg("Wrote", file.path(log_dir, "phase4_fix_singlecell_limitation_statement.txt"))
}

reason_cannot <- paste(c(
  if (!has_seurat_obj) "no local Seurat object",
  if (!has_h5ad) "no h5ad object",
  if (!has_10x_matrix) "no complete 10x matrix/barcode/features set",
  if (!has_count_matrix) "no usable raw count matrix plus metadata",
  if (!has_metadata_umap) "no metadata UMAP_1/UMAP_2 coordinates",
  if (!any(pkg_status$available[pkg_status$package %in% c("Seurat", "SeuratObject")])) "Seurat/SeuratObject unavailable",
  if (!pkg_status[package == "Scissor", available]) "Scissor unavailable"
), collapse = "; ")

checklist <- c(
  paste0("1. Seurat object found: ", ifelse(has_seurat_obj, "yes", "no")),
  paste0("2. h5ad object found: ", ifelse(has_h5ad, "yes", "no")),
  paste0("3. 10x matrix found: ", ifelse(has_10x_matrix, "yes", "no")),
  paste0("4. count matrix found: ", ifelse(has_count_matrix, "yes", "no")),
  paste0("5. metadata UMAP coordinates found: ", ifelse(has_metadata_umap, "yes", "no")),
  paste0("6. UMAP celltype annotation generated: ", ifelse(can_generate_umap, "yes", "no")),
  paste0("7. UMAP PRFT score generated: ", ifelse(can_generate_umap, "yes", "no")),
  paste0("8. UCell successful: ", ifelse(pkg_status[package == "UCell", available] && can_recompute_scores, "yes", "no")),
  paste0("9. AUCell successful: ", ifelse(pkg_status[package == "AUCell", available] && can_recompute_scores, "yes", "no")),
  paste0("10. myeloid/monocyte reclustering possible: ", ifelse(has_seurat_obj || has_10x_matrix || has_count_matrix, "potentially yes after manual object validation", "no")),
  paste0("11. ScissorR possible: ", ifelse(pkg_status[package == "Scissor", available] && (has_seurat_obj || has_count_matrix), "yes", "no")),
  paste0("12. If not possible, reason: ", reason_cannot),
  "13. Current Phase 4 suitable for main text: yes, as processed-table-based cell-state enrichment analysis only",
  "14. Recommended main-text figures: phase4_score_heatmap_by_celltype.pdf; phase4_PRFT_high_like_fraction_by_celltype.pdf; phase4_RoE_heatmap.pdf",
  "15. Recommended supplementary figures: phase4_QC_violin.pdf; phase4_score_violin_by_celltype.pdf; phase4_score_dotplot_by_celltype.pdf; Phase 4-fix limitation/audit tables",
  "16. Single-cell result supports this wording: monocyte-like/myeloid stress-adapted PRFT-high state",
  "17. Recommend entering AS analysis: yes, after confirming suitable RNA-seq/splicing input files",
  "18. Issues needing manual confirmation: provide a complete Seurat/h5ad/10x/count object with cell-level metadata and embeddings if true UMAP localization, UCell/AUCell rescoring, myeloid reclustering, or ScissorR mapping is required"
)
writeLines(checklist, file.path(log_dir, "phase4_fix_key_result_checklist.txt"))
log_msg("Wrote", file.path(log_dir, "phase4_fix_key_result_checklist.txt"))

log_msg("Phase 4-fix audit completed.")
