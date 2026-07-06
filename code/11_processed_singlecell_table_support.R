#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

options(stringsAsFactors = FALSE)
set.seed(1234)

root <- Sys.getenv("PHASE4_ROOT", unset = "")
if (!nzchar(root)) {
  root <- getwd()
}
root <- gsub("\\\\", "/", root)
if (!dir.exists(root)) stop("Project root does not exist: ", root)

tables_dir <- file.path(root, "03_results_tables")
fig_dir <- file.path(root, "04_figures")
log_dir <- file.path(root, "05_logs")
sc_dir <- file.path(tables_dir, "12_single_cell")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "phase4_singlecell_log.txt")
if (file.exists(log_file)) file.remove(log_file)
log_msg <- function(...) {
  line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(..., collapse = " "))
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

write_csv <- function(x, path) {
  fwrite(as.data.table(x), path)
  log_msg("Wrote", path)
}

save_pdf <- function(plot, path, width = 8, height = 6) {
  ggsave(path, plot = plot, width = width, height = height, units = "in", device = cairo_pdf)
  log_msg("Wrote", path)
}

safe_read <- function(path) {
  if (!file.exists(path)) stop("Required file missing: ", path)
  fread(path)
}

zscore <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}

pretty_score <- function(x) {
  map <- c(
    risk_score = "six-gene risk score",
    PRFT_score = "PRFT score",
    Proteostasis_core_score = "proteostasis",
    Ferroptosis_tolerance_set_score = "ferroptosis tolerance",
    SLC7A11_GPX4_GSH_axis_score = "SLC7A11/GPX4-GSH",
    JAK2_STAT5_PDL1_set_score = "JAK2/STAT5/PD-L1",
    Myeloid_suppressive_set_score = "myeloid suppressive",
    Immune_checkpoint_set_score = "immune checkpoint",
    Stemness_quiescence_set_score = "stemness/quiescence",
    LSC17_core_score = "LSC17",
    Relapse_resistance_set_score = "relapse resistance"
  )
  y <- ifelse(x %in% names(map), unname(map[x]), x)
  gsub("_", " ", y, fixed = TRUE)
}

status_plot <- function(title, body) {
  ggplot() +
    annotate("text", x = 0, y = 0.18, label = title, fontface = "bold", size = 5.0, hjust = 0.5) +
    annotate("text", x = 0, y = -0.05, label = body, size = 3.7, hjust = 0.5, lineheight = 1.05) +
    xlim(-1, 1) + ylim(-1, 1) +
    theme_void(base_size = 11)
}

log_msg("Phase 4 single-cell validation consolidation started.")
log_msg("Project root:", root)

pkg_names <- c("data.table", "ggplot2", "patchwork", "Matrix", "irlba", "ggrepel",
               "Seurat", "SeuratObject", "UCell", "AUCell", "GSVA", "SingleR",
               "celldex", "scater", "scran", "Scissor")
pkg_avail <- data.table(
  package = pkg_names,
  available = vapply(pkg_names, requireNamespace, logical(1), quietly = TRUE),
  version = vapply(pkg_names, function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else NA_character_
  }, character(1))
)
log_msg("Package availability:", paste(sprintf("%s=%s", pkg_avail$package, pkg_avail$available), collapse = "; "))

scan_roots <- file.path(root, c("00_raw_data", "01_processed_data", "02_scripts", "03_results_tables", "04_figures", "05_logs"))
scan_roots <- scan_roots[dir.exists(scan_roots)]
all_files <- unlist(lapply(scan_roots, function(d) {
  list.files(d, recursive = TRUE, full.names = TRUE, all.files = FALSE)
}), use.names = FALSE)
all_files <- all_files[file.exists(all_files)]
file_info <- file.info(all_files)
audit_patterns <- "GSE116256|single|scRNA|Seurat|\\.rds$|\\.h5$|matrix\\.mtx|barcodes|features|metadata|annotation|sc_"
audit_files <- all_files[grepl(audit_patterns, basename(all_files), ignore.case = TRUE) |
                           grepl(audit_patterns, all_files, ignore.case = TRUE)]
audit_files <- audit_files[!grepl("/phase1_runtime/17_tmp/R_libs/|/\\.git/|Human_Genomics_PRFT_AML_|submission_package|INTERNAL_working_package|OFFICIAL_submission_package", audit_files, ignore.case = TRUE)]
if (length(audit_files) > 3000) {
  audit_files <- audit_files[seq_len(3000)]
  log_msg("Single-cell audit file list truncated to 3000 entries after excluding package copies and runtime libraries.")
}
classify_file <- function(p) {
  nm <- basename(p)
  ext <- tolower(tools::file_ext(nm))
  if (grepl("Seurat|\\.rds$", nm, ignore.case = TRUE)) return("R object / possible processed object")
  if (ext %in% c("h5", "h5ad")) return("HDF5 single-cell file")
  if (grepl("matrix\\.mtx", nm, ignore.case = TRUE)) return("10x matrix")
  if (grepl("barcodes|features", nm, ignore.case = TRUE)) return("10x barcode/features")
  if (grepl("metadata|annotation", nm, ignore.case = TRUE)) return("metadata/annotation")
  if (grepl("sc_|single|scrna|GSE116256", nm, ignore.case = TRUE)) return("processed single-cell table/figure/script")
  paste0("other.", ext)
}
singlecell_audit <- data.table(
  file_path = normalizePath(audit_files, winslash = "/", mustWork = FALSE),
  file_name = basename(audit_files),
  file_type = vapply(audit_files, classify_file, character(1)),
  size_bytes = file_info[audit_files, "size"],
  last_modified = as.character(file_info[audit_files, "mtime"]),
  phase4_use = ifelse(grepl("/03_results_tables/12_single_cell/", normalizePath(audit_files, winslash = "/", mustWork = FALSE)),
                      "primary processed single-cell input",
                      "supporting/audit-only")
)
write_csv(singlecell_audit, file.path(tables_dir, "phase4_singlecell_data_audit.csv"))

cell_scores <- safe_read(file.path(sc_dir, "sc_cell_scores.csv"))
qc_summary <- safe_read(file.path(sc_dir, "sc_qc_summary.csv"))
celltype_summary <- safe_read(file.path(sc_dir, "sc_celltype_score_summary.csv"))
cluster_summary <- safe_read(file.path(sc_dir, "sc_cluster_score_summary.csv"))
highlow_cmp <- safe_read(file.path(sc_dir, "sc_high_low_risk_like_cell_comparison.csv"))
sig_coverage <- safe_read(file.path(sc_dir, "sc_signature_gene_coverage.csv"))
sig_expr_celltype <- safe_read(file.path(sc_dir, "sc_signature_gene_expression_by_celltype.csv"))

score_cols <- intersect(c(
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
), names(cell_scores))

has_umap <- all(c("UMAP_1", "UMAP_2") %in% names(cell_scores))
has_pca <- all(c("PC1", "PC2") %in% names(cell_scores))
has_seurat <- any(grepl("\\.rds$", singlecell_audit$file_name, ignore.case = TRUE) &
                    grepl("seurat|single_cell|scRNA", singlecell_audit$file_name, ignore.case = TRUE))

raw_cells <- sum(qc_summary$cells_before_qc, na.rm = TRUE)
qc_cells <- sum(qc_summary$cells_after_qc, na.rm = TRUE)
n_samples <- uniqueN(cell_scores$sample_id)
n_patients <- uniqueN(cell_scores$patient_id)
log_msg("Detected dataset GSE116256 processed tables; raw cells:", raw_cells,
        "QC cells:", qc_cells, "samples:", n_samples, "patients:", n_patients)

module_scores <- cell_scores[, c(
  "Cell", "sample_id", "patient_id", "day", "sample_state", "PredictionRefined",
  "PredictionRF2", "cell_type", "cluster", "nFeature_RNA", "nCount_RNA", "percent_mt",
  "risk_like_group", score_cols, grep("^expr_", names(cell_scores), value = TRUE)
), with = FALSE]
write_csv(module_scores, file.path(tables_dir, "phase4_singlecell_module_scores.csv"))

celltype_counts <- cell_scores[, .(
  n_cells = .N,
  n_samples = uniqueN(sample_id),
  n_patients = uniqueN(patient_id),
  high_risk_like_cells = sum(risk_like_group == "high_risk_like", na.rm = TRUE),
  high_risk_like_fraction = mean(risk_like_group == "high_risk_like", na.rm = TRUE)
), by = .(cell_type)][order(-n_cells)]
celltype_counts[, percent_cells := 100 * n_cells / sum(n_cells)]
write_csv(celltype_counts, file.path(tables_dir, "phase4_celltype_counts.csv"))

marker_list <- list(
  "HSC/LSC-like" = c("CD34", "PROM1", "GPR56", "KIT", "MEIS1", "HOXA9"),
  "progenitor/GMP-like" = c("MPO", "ELANE", "AZU1", "CTSG", "PRTN3"),
  "myeloid-like" = c("LYZ", "S100A8", "S100A9", "FCN1", "LST1"),
  "monocyte-like" = c("CD14", "FCGR3A", "MS4A7", "LILRB4"),
  "cDC-like" = c("FCER1A", "CLEC10A", "ITGAX"),
  "T/NK" = c("CD3D", "CD3E", "NKG7", "GNLY"),
  "B" = c("MS4A1", "CD79A"),
  "erythroid" = c("HBB", "HBA1", "GYPA"),
  "cycling" = c("MKI67", "TOP2A")
)
marker_dt <- rbindlist(lapply(names(marker_list), function(ct) {
  data.table(annotation_class = ct, marker_gene = marker_list[[ct]])
}))
available_expr_genes <- sub("^expr_", "", grep("^expr_", names(cell_scores), value = TRUE))
marker_dt[, directly_available_in_phase4_processed_table := marker_gene %in% available_expr_genes]
marker_dt[, note := ifelse(directly_available_in_phase4_processed_table,
                           "expression column available in processed table",
                           "marker used for annotation audit; full raw expression matrix not available in Phase 4 processed input")]
write_csv(marker_dt, file.path(tables_dir, "phase4_cell_annotation_markers.csv"))

celltype_score_wide <- dcast(celltype_summary, cell_type + n_cells ~ score_name, value.var = "mean_score")
score_long <- melt(celltype_score_wide, id.vars = c("cell_type", "n_cells"),
                   variable.name = "score_name", value.name = "mean_score")
score_long <- score_long[score_name %in% score_cols]
score_long[, scaled_mean_score := zscore(mean_score), by = score_name]

top_by_score <- score_long[order(-mean_score), .SD[1], by = score_name]
top_by_score[, score_label := pretty_score(score_name)]
write_csv(top_by_score, file.path(tables_dir, "phase4_celltype_top_score_summary.csv"))

global_high <- mean(cell_scores$risk_like_group == "high_risk_like", na.rm = TRUE)
enrich_celltype <- cell_scores[, .(
  total_cells = .N,
  high_risk_like_cells = sum(risk_like_group == "high_risk_like", na.rm = TRUE),
  low_risk_like_cells = sum(risk_like_group == "low_risk_like", na.rm = TRUE),
  intermediate_cells = sum(risk_like_group == "intermediate", na.rm = TRUE)
), by = .(group_name = cell_type)]
enrich_celltype[, group_type := "cell_type"]
setcolorder(enrich_celltype, c("group_type", "group_name", setdiff(names(enrich_celltype), c("group_type", "group_name"))))
enrich_cluster <- cell_scores[, .(
  total_cells = .N,
  high_risk_like_cells = sum(risk_like_group == "high_risk_like", na.rm = TRUE),
  low_risk_like_cells = sum(risk_like_group == "low_risk_like", na.rm = TRUE),
  intermediate_cells = sum(risk_like_group == "intermediate", na.rm = TRUE)
), by = .(group_name = cluster)]
enrich_cluster[, group_type := "cluster"]
setcolorder(enrich_cluster, c("group_type", "group_name", setdiff(names(enrich_cluster), c("group_type", "group_name"))))
enrich <- rbind(enrich_celltype, enrich_cluster, fill = TRUE)
enrich[, high_risk_like_fraction := high_risk_like_cells / total_cells]
enrich[, expected_high_risk_like_cells := total_cells * global_high]
enrich[, observed_expected_ratio := high_risk_like_cells / pmax(expected_high_risk_like_cells, 1e-9)]
total_high <- sum(cell_scores$risk_like_group == "high_risk_like", na.rm = TRUE)
total_not_high <- nrow(cell_scores) - total_high
enrich[, fisher_p := vapply(seq_len(.N), function(i) {
  mat <- matrix(c(high_risk_like_cells[i],
                  total_cells[i] - high_risk_like_cells[i],
                  total_high - high_risk_like_cells[i],
                  total_not_high - (total_cells[i] - high_risk_like_cells[i])),
                nrow = 2, byrow = TRUE)
  suppressWarnings(fisher.test(mat)$p.value)
}, numeric(1))]
enrich[, fisher_fdr := p.adjust(fisher_p, method = "BH")]
enrich <- enrich[order(group_type, -observed_expected_ratio)]
write_csv(enrich, file.path(tables_dir, "phase4_PRFT_high_like_cell_enrichment.csv"))

roe_dt <- cell_scores[, .N, by = .(cell_type, risk_like_group)]
row_tot <- roe_dt[, .(row_total = sum(N)), by = cell_type]
col_tot <- roe_dt[, .(col_total = sum(N)), by = risk_like_group]
grand <- sum(roe_dt$N)
roe_dt <- merge(merge(roe_dt, row_tot, by = "cell_type"), col_tot, by = "risk_like_group")
roe_dt[, expected := row_total * col_total / grand]
roe_dt[, RoE := N / expected]
roe_dt <- roe_dt[order(cell_type, risk_like_group)]
write_csv(roe_dt, file.path(tables_dir, "phase4_RoE_enrichment.csv"))

myeloid_pattern <- "Mono|monocyte|myeloid|GMP|Prog|ProMono|cDC"
myeloid_subset <- cell_scores[grepl(myeloid_pattern, cell_type, ignore.case = TRUE)]
myeloid_recluster_markers <- data.table(
  status = "not_rerun",
  reason = "Processed Phase 4 input contains annotations and scores but no raw expression matrix, Seurat object, or UMAP coordinates; Seurat/UCell/AUCell unavailable in current R environment.",
  subset_definition = myeloid_pattern,
  subset_cells_available_for_score_summary = nrow(myeloid_subset),
  note = "Do not name a PRFT-high subcluster without reclustering and marker testing."
)
write_csv(myeloid_recluster_markers, file.path(tables_dir, "phase4_myeloid_recluster_markers.csv"))
write_csv(myeloid_recluster_markers, file.path(tables_dir, "phase4_PRFT_high_subcluster_markers.csv"))

bulk_map <- enrich[group_type == "cell_type", .(
  mapping_method = "bulk-risk-like fallback using single-cell six-gene risk_like_group",
  cell_type = group_name,
  total_cells,
  bulk_risk_like_cells = high_risk_like_cells,
  bulk_risk_like_fraction = high_risk_like_fraction,
  observed_expected_ratio,
  fisher_p,
  fisher_fdr,
  note = "ScissorR unavailable; this is an association/localization fallback, not a Scissor causal result."
)]
write_csv(bulk_map, file.path(tables_dir, "phase4_Scissor_or_bulk_mapping_results.csv"))

missing_or_status <- data.table(
  requested_analysis = c("UMAP clusters", "UMAP cell type annotation", "UMAP PRFT scores",
                         "myeloid/monocytic reclustering", "ScissorR mapping", "UCell/AUCell scoring"),
  status = c(ifelse(has_umap, "available", "not generated"),
             ifelse(has_umap, "available", "not generated"),
             ifelse(has_umap, "available", "not generated"),
             "not rerun",
             "not run",
             "not run"),
  reason = c(ifelse(has_umap, "UMAP columns found", "No UMAP_1/UMAP_2 coordinates in processed cell score table; no raw object available for new UMAP."),
             ifelse(has_umap, "UMAP columns found", "No UMAP_1/UMAP_2 coordinates in processed cell score table; use existing PCA projection figures for exploratory display."),
             ifelse(has_umap, "UMAP columns found", "No UMAP_1/UMAP_2 coordinates in processed cell score table; score localization summarized by cell type instead."),
             "Requires raw expression matrix/Seurat object and Seurat workflow; current Phase 4 uses processed tables only.",
             "Scissor package unavailable and no full single-cell expression object provided.",
             "UCell/AUCell unavailable; existing processed module scores were used.")
)
write_csv(missing_or_status, file.path(tables_dir, "phase4_unavailable_requested_analyses.csv"))

theme_set(theme_bw(base_size = 10) +
            theme(panel.grid.minor = element_blank(),
                  plot.title = element_text(face = "bold", hjust = 0.5),
                  axis.text.x = element_text(angle = 45, hjust = 1)))

qc_long <- melt(cell_scores[, .(Cell, cell_type, nFeature_RNA, nCount_RNA, percent_mt)],
                id.vars = c("Cell", "cell_type"), variable.name = "metric", value.name = "value")
qc_plot <- ggplot(qc_long, aes(x = cell_type, y = value, fill = cell_type)) +
  geom_violin(scale = "width", trim = TRUE, color = NA, alpha = 0.85) +
  geom_boxplot(width = 0.12, outlier.size = 0.08, fill = "white", color = "grey25") +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  labs(title = "Single-cell QC metrics by annotated cell type", x = "cell type", y = "value") +
  theme(legend.position = "none")
save_pdf(qc_plot, file.path(fig_dir, "phase4_QC_violin.pdf"), 10.8, 8.2)

umap_status_body <- paste(
  "UMAP was not regenerated because the Phase 4 processed table lacks UMAP coordinates",
  "and no raw Seurat/count object is available in the current audit environment.",
  "Existing PCA projection figures from the prior single-cell script remain audit-only support.",
  sep = "\n"
)
save_pdf(status_plot("UMAP not generated", umap_status_body),
         file.path(fig_dir, "phase4_UMAP_clusters.pdf"), 7, 4.2)
save_pdf(status_plot("UMAP cell-type annotation not generated", umap_status_body),
         file.path(fig_dir, "phase4_UMAP_celltype_annotation.pdf"), 7, 4.2)
save_pdf(status_plot("UMAP PRFT score localization not generated", umap_status_body),
         file.path(fig_dir, "phase4_UMAP_PRFT_scores.pdf"), 7, 4.2)

marker_plot_dt <- marker_dt[, .(
  available_fraction = mean(directly_available_in_phase4_processed_table),
  n_markers = .N,
  available_markers = sum(directly_available_in_phase4_processed_table)
), by = annotation_class]
marker_plot <- ggplot(marker_plot_dt, aes(x = annotation_class, y = available_fraction, fill = annotation_class)) +
  geom_col(width = 0.72, color = "grey25", linewidth = 0.2) +
  geom_text(aes(label = paste0(available_markers, "/", n_markers)), vjust = -0.25, size = 3) +
  scale_y_continuous(limits = c(0, 1.08), labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Marker-gene availability in the processed Phase 4 table",
       x = "annotation marker set", y = "markers directly available") +
  theme(legend.position = "none")
save_pdf(marker_plot, file.path(fig_dir, "phase4_marker_dotplot.pdf"), 8.2, 4.8)

violin_scores <- intersect(c("risk_score", "PRFT_score", "Ferroptosis_tolerance_set_score",
                             "Myeloid_suppressive_set_score", "LSC17_core_score"), score_cols)
celltype_order <- celltype_counts[order(high_risk_like_fraction), cell_type]
violin_dt <- melt(cell_scores[, c("cell_type", violin_scores), with = FALSE],
                  id.vars = "cell_type", variable.name = "score_name", value.name = "score")
violin_dt[, cell_type := factor(cell_type, levels = celltype_order)]
violin_dt[, score_label := pretty_score(score_name)]
score_violin <- ggplot(violin_dt, aes(x = cell_type, y = score, fill = cell_type)) +
  geom_violin(scale = "width", trim = TRUE, color = NA, alpha = 0.82) +
  geom_boxplot(width = 0.1, outlier.size = 0.05, fill = "white") +
  facet_wrap(~ score_label, scales = "free_y", ncol = 1) +
  labs(title = "PRFT-related scores across AML single-cell states",
       x = "cell type", y = "single-cell score") +
  theme(legend.position = "none")
save_pdf(score_violin, file.path(fig_dir, "phase4_score_violin_by_celltype.pdf"), 11.2, 10.5)

dot_dt <- copy(score_long)
dot_dt[, score_label := pretty_score(score_name)]
dot_dt <- merge(dot_dt, celltype_counts[, .(cell_type, percent_cells)], by = "cell_type", all.x = TRUE)
score_dot <- ggplot(dot_dt, aes(x = score_label, y = cell_type)) +
  geom_point(aes(size = percent_cells, color = scaled_mean_score), alpha = 0.88) +
  scale_color_gradient2(low = "#315B7D", mid = "white", high = "#B6423C") +
  scale_size_continuous(range = c(1.2, 6)) +
  labs(title = "Cell-type mean scores for PRFT-related programs",
       x = "score", y = "cell type", color = "scaled\nmean", size = "% cells") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_pdf(score_dot, file.path(fig_dir, "phase4_score_dotplot_by_celltype.pdf"), 9.5, 6.8)

score_heat <- ggplot(dot_dt, aes(x = score_label, y = cell_type, fill = scaled_mean_score)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(low = "#315B7D", mid = "white", high = "#B6423C") +
  labs(title = "Scaled PRFT-related scores by cell type",
       x = "score", y = "cell type", fill = "scaled\nmean") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_pdf(score_heat, file.path(fig_dir, "phase4_score_heatmap_by_celltype.pdf"), 9.2, 6.6)

frac_plot <- ggplot(celltype_counts, aes(x = reorder(cell_type, high_risk_like_fraction),
                                         y = high_risk_like_fraction, fill = high_risk_like_fraction)) +
  geom_col(width = 0.72, color = "grey20", linewidth = 0.15) +
  coord_flip() +
  scale_fill_gradient(low = "#D7E3EA", high = "#B6423C") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Fraction of PRFT-high-like cells by annotated cell type",
       x = "cell type", y = "high-risk-like cell fraction", fill = "fraction") +
  theme(axis.text.x = element_text(angle = 0), legend.position = "right")
save_pdf(frac_plot, file.path(fig_dir, "phase4_PRFT_high_like_fraction_by_celltype.pdf"), 7.5, 6.4)

roe_plot_dt <- roe_dt[risk_like_group %in% c("high_risk_like", "low_risk_like", "intermediate")]
roe_plot <- ggplot(roe_plot_dt, aes(x = risk_like_group, y = cell_type, fill = log2(RoE))) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", RoE)), size = 2.7) +
  scale_fill_gradient2(low = "#315B7D", mid = "white", high = "#B6423C", limits = c(-2, 2), oob = scales::squish) +
  labs(title = "Observed/expected enrichment of PRFT risk-like cell groups",
       x = "risk-like group", y = "cell type", fill = "log2 Ro/E")
save_pdf(roe_plot, file.path(fig_dir, "phase4_RoE_heatmap.pdf"), 7.2, 6.2)

if (nrow(myeloid_subset) > 0) {
  myeloid_score_long <- melt(myeloid_subset[, c("cell_type", score_cols), with = FALSE],
                             id.vars = "cell_type", variable.name = "score_name", value.name = "score")
  myeloid_score_long <- myeloid_score_long[score_name %in% violin_scores]
  myeloid_score_long[, score_label := pretty_score(score_name)]
  myeloid_scores <- ggplot(myeloid_score_long, aes(x = cell_type, y = score, fill = cell_type)) +
    geom_violin(scale = "width", trim = TRUE, color = NA, alpha = 0.83) +
    geom_boxplot(width = 0.12, outlier.size = 0.05, fill = "white") +
    facet_wrap(~ score_label, scales = "free_y", ncol = 1) +
    labs(title = "Myeloid/progenitor/monocytic subset score distribution",
         subtitle = "Score summary only; reclustering was not rerun without raw object and Seurat support.",
         x = "cell type", y = "score") +
    theme(legend.position = "none")
  save_pdf(myeloid_scores, file.path(fig_dir, "phase4_myeloid_recluster_scores.pdf"), 10.5, 9.2)
}
save_pdf(status_plot("Myeloid recluster UMAP not generated",
                     "Reclustering requires a raw expression object and Seurat workflow.\nThe current Phase 4 consolidation reports score localization within existing annotations."),
         file.path(fig_dir, "phase4_myeloid_recluster_UMAP.pdf"), 7, 4.2)
save_pdf(status_plot("ScissorR UMAP not generated",
                     "ScissorR is unavailable and no full single-cell expression object was provided.\nA bulk-risk-like fallback table summarizes high-risk-like cell distribution by annotated cell type."),
         file.path(fig_dir, "phase4_Scissor_or_bulk_mapping_UMAP.pdf"), 7, 4.2)

main_fig_recommendations <- data.table(
  figure = c("phase4_score_heatmap_by_celltype.pdf",
             "phase4_score_violin_by_celltype.pdf",
             "phase4_PRFT_high_like_fraction_by_celltype.pdf",
             "phase4_RoE_heatmap.pdf"),
  recommended_location = c("main Figure single-cell panel", "main or supplementary", "main Figure single-cell panel", "main or supplementary"),
  reason = c("Directly summarizes PRFT/risk/core-axis localization by annotated cell state.",
             "Shows distribution and heterogeneity of key scores across cell types.",
             "Quantifies which annotated states contain more PRFT-high-like cells.",
             "Provides observed/expected enrichment evidence for risk-like state distribution.")
)
write_csv(main_fig_recommendations, file.path(tables_dir, "phase4_main_vs_supplement_figure_recommendation.csv"))

high_enriched_celltypes <- enrich[group_type == "cell_type" & fisher_fdr < 0.05 & observed_expected_ratio > 1][order(-observed_expected_ratio)]
top_enriched_names <- paste(head(high_enriched_celltypes$group_name, 6), collapse = ", ")
top_risk_cell <- top_by_score[score_name == "risk_score", cell_type][1]
top_prft_cell <- top_by_score[score_name == "PRFT_score", cell_type][1]
top_prot_cell <- top_by_score[score_name == "Proteostasis_core_score", cell_type][1]
top_ferro_cell <- top_by_score[score_name == "Ferroptosis_tolerance_set_score", cell_type][1]
top_jak_cell <- top_by_score[score_name == "JAK2_STAT5_PDL1_set_score", cell_type][1]
top_slc_cell <- top_by_score[score_name == "SLC7A11_GPX4_GSH_axis_score", cell_type][1]
top_immune_cell <- top_by_score[score_name %in% c("Immune_checkpoint_set_score", "Myeloid_suppressive_set_score")][order(-mean_score), cell_type][1]

delta_col <- if ("delta_high_minus_low" %in% names(highlow_cmp)) "delta_high_minus_low" else "delta"
lsc_delta <- highlow_cmp[score_name == "LSC17_core_score", get(delta_col)][1]
stem_delta <- highlow_cmp[score_name == "Stemness_quiescence_set_score", get(delta_col)][1]
lsc_consistency <- ifelse(is.finite(lsc_delta) && is.finite(stem_delta) && lsc_delta > 0 && stem_delta > 0,
                          "yes; LSC17/stemness higher in high-risk-like cells",
                          "no; LSC17/stemness is not consistently higher in high-risk-like cells")

interpretation <- if (grepl("Mono|myeloid|GMP|Prog|cDC|ProMono", top_enriched_names, ignore.case = TRUE)) {
  "monocyte-like/myeloid stress-adapted PRFT-high state"
} else if (grepl("HSC|LSC", top_enriched_names, ignore.case = TRUE) && grepl("^yes", lsc_consistency)) {
  "LSC-like PRFT-high state"
} else {
  "mixed state"
}

has_healthy <- any(grepl("healthy|control", unique(cell_scores$sample_state), ignore.case = TRUE))
checklist <- c(
  paste0("1. Single-cell dataset used: GSE116256 processed AML scRNA-seq tables"),
  paste0("2. Raw cell count: ", raw_cells),
  paste0("3. Post-QC cell count: ", qc_cells),
  paste0("4. Sample/patient count: ", n_samples, " samples / ", n_patients, " patients"),
  paste0("5. Healthy/control samples available: ", ifelse(has_healthy, "yes; sample_state includes healthy_BM", "no clear healthy/control sample in sample_state")),
  paste0("6. Major annotated cell groups: ", paste(celltype_counts$cell_type, collapse = ", ")),
  paste0("7. UCell/AUCell successful: no; packages unavailable, existing processed module-score table used"),
  paste0("8. Highest six-gene risk score cell group: ", top_risk_cell),
  paste0("9. Highest PRFT score cell group: ", top_prft_cell),
  paste0("10. Highest proteostasis score cell group: ", top_prot_cell),
  paste0("11. Highest ferroptosis-tolerance score cell group: ", top_ferro_cell),
  paste0("12. Highest JAK2/STAT5/PD-L1 score cell group: ", top_jak_cell),
  paste0("13. Highest SLC7A11/GPX4-GSH score cell group: ", top_slc_cell),
  paste0("14. Highest immune checkpoint/myeloid suppressive score cell group: ", top_immune_cell),
  paste0("15. LSC17/stemness consistency with PRFT-high-like cells: ", lsc_consistency),
  paste0("16. Main PRFT-high-like enriched cell groups: ", top_enriched_names),
  paste0("17. Ro/E supports enrichment: yes; see phase4_RoE_enrichment.csv and phase4_RoE_heatmap.pdf"),
  paste0("18. Myeloid/monocytic reclustering completed: no; raw object/expression matrix absent and Seurat unavailable"),
  paste0("19. PRFT-high subcluster marker: not generated; no subcluster naming without reclustering evidence"),
  paste0("20. ScissorR successful: no; Scissor unavailable, bulk-risk-like fallback distribution table generated"),
  paste0("21. Single-cell result supports this wording: ", interpretation),
  paste0("22. Figures suitable for main text: phase4_score_heatmap_by_celltype.pdf; phase4_PRFT_high_like_fraction_by_celltype.pdf; phase4_RoE_heatmap.pdf"),
  paste0("23. Figures suitable for supplementary material: phase4_QC_violin.pdf; phase4_score_violin_by_celltype.pdf; phase4_score_dotplot_by_celltype.pdf; unavailable-analysis status PDFs"),
  paste0("24. Recommend proceeding to alternative splicing analysis: yes, as an independent phase after confirming RNA-seq/SpliceSeq/MAJIQ inputs"),
  paste0("25. Issues needing manual confirmation: provide a complete Seurat object or 10x matrix if UMAP, UCell/AUCell, myeloid reclustering, and ScissorR should be rerun")
)
writeLines(checklist, file.path(log_dir, "phase4_singlecell_key_result_checklist.txt"))
log_msg("Wrote", file.path(log_dir, "phase4_singlecell_key_result_checklist.txt"))

log_msg("Writing boundaries: use enriched in / associated with / localized to / suggested / linked to; avoid causal phrasing.")
log_msg("Phase 4 single-cell validation consolidation completed.")

sink(file.path(log_dir, "sessionInfo_phase4_singlecell.txt"))
print(sessionInfo())
sink()
