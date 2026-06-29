#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

local_lib <- normalizePath("17_tmp/R_libs", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
source("15_scripts/plot_label_utils.R")

if (!requireNamespace("WGCNA", quietly = TRUE)) {
  stop("WGCNA is not installed. Please install WGCNA before running this script.")
}

suppressPackageStartupMessages({
  library(WGCNA)
})

dir.create("06_wgcna", recursive = TRUE, showWarnings = FALSE)
dir.create("13_figures", recursive = TRUE, showWarnings = FALSE)
dir.create("14_tables", recursive = TRUE, showWarnings = FALSE)
dir.create("16_logs", recursive = TRUE, showWarnings = FALSE)

allowWGCNAThreads()

save_session_info <- function(path) {
  writeLines(capture.output(sessionInfo()), con = path)
}

coerce_numeric <- function(x) {
  out <- suppressWarnings(as.numeric(x))
  out
}

pick_trait_cols <- function(df, requested) {
  keep <- intersect(requested, colnames(df))
  if (length(keep) == 0) {
    return(data.frame(row.names = rownames(df)))
  }
  out <- df[, keep, drop = FALSE]
  for (nm in colnames(out)) {
    out[[nm]] <- coerce_numeric(out[[nm]])
  }
  out
}

obj <- readRDS("02_processed_data/tcga_expr_clin_matched.rds")
prft_df <- as.data.table(readRDS("04_prft_score/tcga_prft_score.rds"))
deg_up <- fread("05_deg/deg_prft_high_upregulated_tcga.csv")
deg_all <- fread("05_deg/deg_prft_high_vs_low_tcga_all.csv")
ssgsea_scores <- readRDS("04_prft_score/tcga_ssgsea_scores.rds")

expr <- as.matrix(obj$expr)
storage.mode(expr) <- "numeric"
if (is.null(rownames(expr)) || is.null(colnames(expr))) {
  stop("Expression matrix must have gene symbols as rownames and sample IDs as colnames.")
}

if (!all(c("sample_id", "PRFT_score", "z_Proteostasis_core", "z_Ferroptosis_tolerance_set",
           "Proteostasis_core_score", "Ferroptosis_tolerance_set_score",
           "OS_status", "OS_time", "age") %in% colnames(prft_df))) {
  stop("tcga_prft_score.rds is missing required columns for WGCNA traits.")
}

ssgsea_scores <- as.matrix(ssgsea_scores)
score_dt <- data.table(
  sample_id = colnames(ssgsea_scores)
)
score_rows_to_extract <- c(
  "SLC7A11_GPX4_GSH_axis",
  "SUMOylation_set",
  "NEDDylation_set",
  "JAK2_STAT5_PDL1_set",
  "LSC17_core",
  "Stemness_quiescence_set",
  "Relapse_resistance_set",
  "Immune_checkpoint_set",
  "T_cell_exhaustion_set",
  "Myeloid_suppressive_set"
)
for (nm in score_rows_to_extract) {
  col_nm <- paste0(nm, "_score")
  if (nm %in% rownames(ssgsea_scores)) {
    score_dt[[col_nm]] <- as.numeric(ssgsea_scores[nm, ])
  }
}

trait_base <- unique(prft_df[, .(
  sample_id,
  PRFT_score,
  z_Proteostasis_core,
  z_Ferroptosis_tolerance_set,
  Proteostasis_core_score,
  Ferroptosis_tolerance_set_score,
  OS_status,
  OS_time,
  age
)])

trait_dt <- merge(trait_base, score_dt, by = "sample_id", all.x = TRUE, sort = FALSE)

common_samples <- intersect(colnames(expr), trait_dt$sample_id)
if (length(common_samples) == 0) {
  stop("No overlapping samples between expression matrix and trait data.")
}

expr <- expr[, common_samples, drop = FALSE]
trait_dt <- trait_dt[match(common_samples, trait_dt$sample_id)]
if (!identical(colnames(expr), trait_dt$sample_id)) {
  stop("Sample order mismatch between expression matrix and trait data after matching.")
}

samples_input <- ncol(expr)
genes_input <- nrow(expr)

mad_values <- apply(expr, 1, mad, na.rm = TRUE)
mad_values[is.na(mad_values)] <- -Inf
if (nrow(expr) > 8000) {
  selected_genes <- names(sort(mad_values, decreasing = TRUE))[1:8000]
} else {
  selected_genes <- rownames(expr)
}

expr_top <- expr[selected_genes, , drop = FALSE]
datExpr0 <- t(expr_top)
saveRDS(datExpr0, "06_wgcna/wgcna_input_expr_topMAD8000.rds")

gsg <- goodSamplesGenes(datExpr0, verbose = 3)
removed_samples_good <- character(0)
removed_genes_good <- character(0)
if (!gsg$allOK) {
  if (sum(!gsg$goodSamples) > 0) {
    removed_samples_good <- rownames(datExpr0)[!gsg$goodSamples]
  }
  if (sum(!gsg$goodGenes) > 0) {
    removed_genes_good <- colnames(datExpr0)[!gsg$goodGenes]
  }
  datExpr <- datExpr0[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  trait_dt <- trait_dt[match(rownames(datExpr), trait_dt$sample_id)]
} else {
  datExpr <- datExpr0
}

sampleTree <- hclust(dist(datExpr), method = "average")
grDevices::pdf("13_figures/Figure3_WGCNA_sample_clustering.pdf", width = 10, height = 6)
par(cex = 0.7)
plot(sampleTree, main = "Sample clustering for WGCNA", sub = "", xlab = "", cex.lab = 1.2)
dev.off()
grDevices::png("13_figures/Figure3_WGCNA_sample_clustering.png", width = 3000, height = 1800, res = 300)
par(cex = 0.7)
plot(sampleTree, main = "Sample clustering for WGCNA", sub = "", xlab = "", cex.lab = 1.2)
dev.off()

samples_removed <- removed_samples_good
sample_removal_reason <- if (length(samples_removed) > 0) "Removed by goodSamplesGenes" else "No sample removed; no obvious outlier removed."

trait_df <- as.data.frame(trait_dt)
rownames(trait_df) <- trait_df$sample_id
trait_df$sample_id <- NULL

requested_traits <- c(
  "PRFT_score",
  "z_Proteostasis_core",
  "z_Ferroptosis_tolerance_set",
  "Proteostasis_core_score",
  "Ferroptosis_tolerance_set_score",
  "SLC7A11_GPX4_GSH_axis_score",
  "SUMOylation_set_score",
  "NEDDylation_set_score",
  "JAK2_STAT5_PDL1_set_score",
  "LSC17_core_score",
  "Stemness_quiescence_set_score",
  "Relapse_resistance_set_score",
  "Immune_checkpoint_set_score",
  "T_cell_exhaustion_set_score",
  "Myeloid_suppressive_set_score",
  "OS_status",
  "OS_time",
  "age"
)
trait_numeric <- pick_trait_cols(trait_df, requested_traits)
trait_numeric <- trait_numeric[rownames(datExpr), , drop = FALSE]

trait_out <- data.table(sample_id = rownames(trait_numeric), trait_numeric, keep.rownames = FALSE)
fwrite(trait_out, "06_wgcna/wgcna_sample_trait_data.csv")

powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 5)
sft_df <- as.data.table(sft$fitIndices)
fwrite(sft_df, "06_wgcna/wgcna_soft_threshold_table.csv")

sft_r2_col <- if ("SFT.R.sq" %in% colnames(sft_df)) "SFT.R.sq" else "SFT.R.sq"
candidate_power_df <- sft_df[!is.na(get(sft_r2_col))]
power_ge_085 <- candidate_power_df[get(sft_r2_col) >= 0.85]
if (nrow(power_ge_085) > 0) {
  soft_power <- power_ge_085[order(Power)][1, Power]
  soft_power_reason <- "First power with scale-free topology fit R^2 >= 0.85"
} else {
  fallback_row <- candidate_power_df[order(-get(sft_r2_col), -mean.k., Power)][1]
  soft_power <- fallback_row$Power
  soft_power_reason <- "No power reached R^2 >= 0.85; selected highest R^2 with reasonable mean connectivity"
}

grDevices::pdf("13_figures/Figure3_WGCNA_soft_threshold.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
plot(sft_df$Power, sft_df$SFT.R.sq, type = "b", xlab = "Soft threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", main = "Scale independence")
abline(h = 0.85, col = "red", lty = 2)
abline(v = soft_power, col = "blue", lty = 2)
plot(sft_df$Power, sft_df$mean.k., type = "b", xlab = "Soft threshold (power)", ylab = "Mean connectivity", main = "Mean connectivity")
abline(v = soft_power, col = "blue", lty = 2)
dev.off()
grDevices::png("13_figures/Figure3_WGCNA_soft_threshold.png", width = 3000, height = 1500, res = 300)
par(mfrow = c(1, 2))
plot(sft_df$Power, sft_df$SFT.R.sq, type = "b", xlab = "Soft threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", main = "Scale independence")
abline(h = 0.85, col = "red", lty = 2)
abline(v = soft_power, col = "blue", lty = 2)
plot(sft_df$Power, sft_df$mean.k., type = "b", xlab = "Soft threshold (power)", ylab = "Mean connectivity", main = "Mean connectivity")
abline(v = soft_power, col = "blue", lty = 2)
dev.off()

net <- blockwiseModules(
  datExpr,
  power = soft_power,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  numericLabels = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)

module_colors <- net$colors
module_levels <- sort(unique(module_colors))
module_assignments <- data.table(
  gene_symbol = colnames(datExpr),
  module_color = module_colors
)
fwrite(module_assignments, "06_wgcna/wgcna_module_assignments.csv")

grDevices::pdf("13_figures/Figure3_WGCNA_gene_dendrogram_modules.pdf", width = 12, height = 6)
plotDendroAndColors(
  net$dendrograms[[1]],
  net$colors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene dendrogram and module colors"
)
dev.off()
grDevices::png("13_figures/Figure3_WGCNA_gene_dendrogram_modules.png", width = 3600, height = 1800, res = 300)
plotDendroAndColors(
  net$dendrograms[[1]],
  net$colors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene dendrogram and module colors"
)
dev.off()

MEs0 <- moduleEigengenes(datExpr, colors = module_colors)$eigengenes
MEs <- orderMEs(MEs0)

module_trait_cor <- cor(MEs, trait_numeric, use = "pairwise.complete.obs", method = "pearson")
module_trait_p <- corPvalueStudent(module_trait_cor, nrow(datExpr))

module_trait_cor_dt <- as.data.table(module_trait_cor, keep.rownames = "module")
module_trait_p_dt <- as.data.table(module_trait_p, keep.rownames = "module")
fwrite(module_trait_cor_dt, "06_wgcna/wgcna_module_trait_cor.csv")
fwrite(module_trait_p_dt, "06_wgcna/wgcna_module_trait_pvalue.csv")

text_matrix <- paste(signif(module_trait_cor, 2), "\n(", signif(module_trait_p, 2), ")", sep = "")
dim(text_matrix) <- dim(module_trait_cor)

grDevices::pdf("13_figures/Figure3_WGCNA_module_trait_heatmap.pdf", width = 12, height = 8)
labeledHeatmap(
  Matrix = module_trait_cor,
  xLabels = pretty_label(colnames(trait_numeric)),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = text_matrix,
  setStdMargins = FALSE,
  cex.text = 0.6,
  zlim = c(-1, 1),
  main = "Module-trait relationships"
)
dev.off()
grDevices::png("13_figures/Figure3_WGCNA_module_trait_heatmap.png", width = 3600, height = 2400, res = 300)
labeledHeatmap(
  Matrix = module_trait_cor,
  xLabels = pretty_label(colnames(trait_numeric)),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = text_matrix,
  setStdMargins = FALSE,
  cex.text = 0.6,
  zlim = c(-1, 1),
  main = "Module-trait relationships"
)
dev.off()

prft_module_stats <- data.table(
  module = rownames(module_trait_cor),
  module_color = sub("^ME", "", rownames(module_trait_cor)),
  cor_with_PRFT_score = module_trait_cor[, "PRFT_score"],
  P.Value = module_trait_p[, "PRFT_score"]
)
prft_module_stats[, FDR := p.adjust(P.Value, method = "BH")]
module_sizes <- table(module_colors)
prft_module_stats[, module_size := as.integer(module_sizes[module_color])]
prft_module_stats[, selection_type := fifelse(cor_with_PRFT_score > 0.30 & FDR < 0.05, "FDR_significant",
                                              fifelse(cor_with_PRFT_score > 0.30 & P.Value < 0.05, "exploratory_rawP", "not_selected"))]
prft_module_stats <- prft_module_stats[order(-cor_with_PRFT_score, P.Value)]
fwrite(prft_module_stats, "06_wgcna/wgcna_prft_module_summary.csv")

selected_modules_fdr <- prft_module_stats[selection_type == "FDR_significant", module_color]
if (length(selected_modules_fdr) > 0) {
  prft_positive_modules <- selected_modules_fdr
  module_selection_rule <- "cor_with_PRFT_score > 0.30 and FDR < 0.05"
} else {
  prft_positive_modules <- prft_module_stats[selection_type == "exploratory_rawP", module_color]
  module_selection_rule <- "No FDR-significant module; exploratory modules selected by cor_with_PRFT_score > 0.30 and raw P < 0.05"
}

prft_positive_module_genes <- module_assignments[module_color %in% prft_positive_modules]
fwrite(prft_positive_module_genes, "06_wgcna/wgcna_prft_positive_module_genes.csv")
writeLines(prft_positive_module_genes$gene_symbol, "06_wgcna/wgcna_prft_positive_module_genes.txt")

if (length(prft_positive_modules) > 0) {
  barplot_dt <- prft_module_stats[module_color %in% prft_positive_modules]
  barplot_dt$module_color_factor <- factor(barplot_dt$module_color, levels = barplot_dt$module_color)
  p_bar <- ggplot(barplot_dt, aes(x = module_color_factor, y = cor_with_PRFT_score, fill = module_color_factor)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(label = paste0("n=", module_size)), vjust = -0.4, size = 3) +
    scale_fill_manual(values = setNames(barplot_dt$module_color, barplot_dt$module_color)) +
    labs(x = "Module", y = "Correlation with PRFT score", title = "PRFT-positive WGCNA modules") +
    theme_bw(base_size = 11)
  ggsave("13_figures/Figure3_WGCNA_PRFT_positive_modules_barplot.pdf", p_bar, width = 7, height = 5)
  ggsave("13_figures/Figure3_WGCNA_PRFT_positive_modules_barplot.png", p_bar, width = 7, height = 5, dpi = 300)
}

gene_module_membership <- as.data.frame(cor(datExpr, MEs, use = "pairwise.complete.obs"))
gene_module_membership_p <- as.data.frame(corPvalueStudent(as.matrix(gene_module_membership), nrow(datExpr)))
colnames(gene_module_membership) <- names(MEs)
colnames(gene_module_membership_p) <- paste0(names(MEs), "_P")

gs_prft <- cor(datExpr, trait_numeric$PRFT_score, use = "pairwise.complete.obs")
gs_prft_p <- corPvalueStudent(as.matrix(gs_prft), nrow(datExpr))

gene_level_dt <- data.table(
  gene_symbol = colnames(datExpr),
  module_color = module_colors,
  GS_PRFT = as.numeric(gs_prft),
  GS_PRFT_P = as.numeric(gs_prft_p)
)

assigned_me_names <- paste0("ME", module_colors)
gene_level_dt[, MM := NA_real_]
gene_level_dt[, MM_P := NA_real_]
for (i in seq_len(nrow(gene_level_dt))) {
  me_name <- assigned_me_names[i]
  if (me_name %in% colnames(gene_module_membership)) {
    gene_level_dt$MM[i] <- gene_module_membership[i, me_name]
    gene_level_dt$MM_P[i] <- gene_module_membership_p[i, paste0(me_name, "_P")]
  }
}
fwrite(gene_level_dt, "06_wgcna/wgcna_gene_level_statistics.csv")

deg_up_genes <- unique(deg_up$gene_symbol)
candidate_pool <- merge(
  prft_positive_module_genes,
  gene_level_dt,
  by = c("gene_symbol", "module_color"),
  all.x = TRUE,
  sort = FALSE
)
candidate_pool <- merge(
  candidate_pool,
  deg_all[, .(gene_symbol, logFC, P.Value, adj.P.Val)],
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)
candidate_pool <- merge(
  candidate_pool,
  prft_module_stats[, .(module_color, module_cor_with_PRFT = cor_with_PRFT_score, module_P = P.Value, module_FDR = FDR, selection_type)],
  by = "module_color",
  all.x = TRUE,
  sort = FALSE
)
candidate_pool <- candidate_pool[gene_symbol %in% deg_up_genes]
candidate_pool[, module_priority := fifelse(selection_type == "FDR_significant", 1L,
                                            fifelse(selection_type == "exploratory_rawP", 0L, -1L))]
candidate_pool[, abs_GS_PRFT := abs(GS_PRFT)]
candidate_pool[, abs_MM := abs(MM)]
setorder(candidate_pool, -module_priority, -abs_GS_PRFT, -abs_MM, adj.P.Val, -logFC)
candidate_pool[, c("module_priority", "abs_GS_PRFT", "abs_MM") := NULL]
fwrite(candidate_pool, "06_wgcna/wgcna_deg_up_intersect_prft_module_genes.csv")
writeLines(candidate_pool$gene_symbol, "06_wgcna/wgcna_candidate_pool_genes.txt")

summary_dt <- data.table(
  samples_input = samples_input,
  samples_used = nrow(datExpr),
  samples_removed = length(samples_removed),
  genes_input = genes_input,
  genes_used_for_wgcna = ncol(datExpr),
  soft_power_selected = soft_power,
  network_type = "signed",
  min_module_size = 30,
  merge_cut_height = 0.25,
  total_modules_detected = length(setdiff(unique(module_colors), "grey")),
  PRFT_positive_modules_count = length(prft_positive_modules),
  PRFT_positive_module_genes_count = nrow(prft_positive_module_genes),
  PRFT_high_upregulated_DEGs_count = length(deg_up_genes),
  candidate_pool_genes_count = nrow(candidate_pool),
  candidate_pool_definition = "PRFT_high_upregulated_DEGs intersected with genes from PRFT_score positively correlated WGCNA modules.",
  sample_removal_reason = sample_removal_reason,
  soft_power_reason = soft_power_reason,
  module_selection_rule = module_selection_rule
)
fwrite(summary_dt, "14_tables/tcga_wgcna_summary.csv")

save_session_info("16_logs/sessionInfo_11_wgcna_prft_modules_tcga.txt")
