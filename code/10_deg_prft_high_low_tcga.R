#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
source("15_scripts/plot_label_utils.R")

dir.create("05_deg", recursive = TRUE, showWarnings = FALSE)
dir.create("13_figures", recursive = TRUE, showWarnings = FALSE)
dir.create("14_tables", recursive = TRUE, showWarnings = FALSE)
dir.create("16_logs", recursive = TRUE, showWarnings = FALSE)

save_session_info <- function(path) {
  writeLines(capture.output(sessionInfo()), con = path)
}

zscore_rows <- function(mat) {
  mat <- as.matrix(mat)
  t(apply(mat, 1, function(x) {
    s <- stats::sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) {
      rep(0, length(x))
    } else {
      as.numeric(scale(x))
    }
  }))
}

obj <- readRDS("02_processed_data/tcga_expr_clin_matched.rds")
prft_score_df <- readRDS("04_prft_score/tcga_prft_score.rds")
prft_group_df <- fread("04_prft_score/tcga_prft_group.csv")

expr <- obj$expr
if (is.null(expr) || is.null(rownames(expr)) || is.null(colnames(expr))) {
  stop("obj$expr must be a matrix-like object with gene symbols as rownames and sample IDs as colnames.")
}
expr <- as.matrix(expr)
storage.mode(expr) <- "numeric"

required_group_cols <- c("sample_id", "PRFT_group")
if (!all(required_group_cols %in% colnames(prft_group_df))) {
  stop("PRFT group file must contain columns: ", paste(required_group_cols, collapse = ", "))
}

if (!all(required_group_cols %in% colnames(prft_score_df))) {
  stop("PRFT score RDS must contain columns: ", paste(required_group_cols, collapse = ", "))
}

prft_group_df <- unique(prft_group_df[, .(sample_id, PRFT_group)])
prft_score_df <- unique(as.data.table(prft_score_df)[, .(sample_id, PRFT_group)])

group_check <- merge(
  prft_group_df,
  prft_score_df,
  by = "sample_id",
  suffixes = c(".csv", ".rds"),
  all = FALSE
)

if (nrow(group_check) == 0) {
  stop("No overlapping sample IDs between PRFT group CSV and PRFT score RDS.")
}

if (any(group_check$PRFT_group.csv != group_check$PRFT_group.rds)) {
  stop("Inconsistent PRFT_group assignments between tcga_prft_group.csv and tcga_prft_score.rds.")
}

common_samples <- intersect(colnames(expr), prft_group_df$sample_id)
if (length(common_samples) == 0) {
  stop("No overlapping sample IDs between expression matrix and PRFT group file.")
}

prft_group_use <- prft_group_df[match(common_samples, prft_group_df$sample_id)]
expr_use <- expr[, common_samples, drop = FALSE]

if (!identical(colnames(expr_use), prft_group_use$sample_id)) {
  stop("Expression matrix column order and PRFT group sample_id order are not aligned after matching.")
}

group_factor <- factor(prft_group_use$PRFT_group, levels = c("PRFT_low", "PRFT_high"))
if (any(is.na(group_factor))) {
  stop("PRFT_group contains unexpected values. Expected PRFT_low / PRFT_high.")
}

design <- model.matrix(~ 0 + group_factor)
colnames(design) <- c("PRFT_low", "PRFT_high")

fit <- limma::lmFit(expr_use, design)
contrast_matrix <- limma::makeContrasts(PRFT_high - PRFT_low, levels = design)
fit2 <- limma::contrasts.fit(fit, contrast_matrix)
fit2 <- limma::eBayes(fit2)

deg_all <- limma::topTable(fit2, number = Inf, sort.by = "P", adjust.method = "BH")
deg_all$gene_symbol <- rownames(deg_all)
deg_all <- deg_all[, c("gene_symbol", setdiff(colnames(deg_all), "gene_symbol"))]

deg_all_dt <- as.data.table(deg_all)
deg_all_dt[, regulation := fifelse(adj.P.Val < 0.05 & logFC > 0.5, "Upregulated_in_PRFT_high",
                                   fifelse(adj.P.Val < 0.05 & logFC < -0.5, "Downregulated_in_PRFT_high",
                                           "Not_significant"))]
deg_all_dt[, ranking_metric := sign(logFC) * -log10(P.Value)]
deg_all_dt[, neg_log10_adjP := -log10(pmax(adj.P.Val, .Machine$double.xmin))]

deg_up <- deg_all_dt[adj.P.Val < 0.05 & logFC > 0.5]
deg_down <- deg_all_dt[adj.P.Val < 0.05 & logFC < -0.5]

fwrite(deg_all_dt, "05_deg/deg_prft_high_vs_low_tcga_all.csv")
fwrite(deg_up, "05_deg/deg_prft_high_upregulated_tcga.csv")
fwrite(deg_down, "05_deg/deg_prft_high_downregulated_tcga.csv")
writeLines(deg_up$gene_symbol, con = "05_deg/deg_prft_high_up_genes.txt")
writeLines(deg_down$gene_symbol, con = "05_deg/deg_prft_high_down_genes.txt")

ranked_gene_list <- deg_all_dt[, .(gene_symbol, logFC, P.Value, adj.P.Val, ranking_metric)]
fwrite(ranked_gene_list, "05_deg/deg_ranked_gene_list_tcga.csv")

volcano_dt <- copy(deg_all_dt)
top_up_labels <- deg_up[order(adj.P.Val, -logFC)][1:min(10, nrow(deg_up)), gene_symbol]
top_down_labels <- deg_down[order(adj.P.Val, logFC)][1:min(10, nrow(deg_down)), gene_symbol]
volcano_dt[, label := ifelse(gene_symbol %in% c(top_up_labels, top_down_labels), gene_symbol, "")]
fwrite(volcano_dt, "05_deg/deg_volcano_input_tcga.csv")

summary_dt <- data.table(
  total_genes_tested = nrow(deg_all_dt),
  upregulated_genes = nrow(deg_up),
  downregulated_genes = nrow(deg_down),
  significant_genes_total = nrow(deg_up) + nrow(deg_down),
  threshold_adjP = 0.05,
  threshold_logFC = 0.5,
  PRFT_high_samples = sum(group_factor == "PRFT_high"),
  PRFT_low_samples = sum(group_factor == "PRFT_low")
)
fwrite(summary_dt, "14_tables/tcga_deg_summary.csv")

volcano_plot <- ggplot(volcano_dt, aes(x = logFC, y = neg_log10_adjP, color = regulation)) +
  geom_point(alpha = 0.75, size = 1.3) +
  scale_color_manual(
    values = c(
      Upregulated_in_PRFT_high = "#D62728",
      Downregulated_in_PRFT_high = "#1F77B4",
      Not_significant = "grey70"
    ),
    breaks = c("Upregulated_in_PRFT_high", "Downregulated_in_PRFT_high", "Not_significant"),
    labels = pretty_label(c("Upregulated_in_PRFT_high", "Downregulated_in_PRFT_high", "Not_significant"))
  ) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_text(
    data = volcano_dt[label != ""],
    aes(label = label),
    size = 2.7,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  labs(
    title = "Differential expression: PRFT high vs PRFT low",
    x = "logFC (PRFT high vs PRFT low)",
    y = "-log10(adj.P.Val)",
    color = NULL
  ) +
  theme_bw(base_size = 11)

ggsave("13_figures/Figure2_PRFT_DEG_volcano.pdf", volcano_plot, width = 7, height = 5.5)
ggsave("13_figures/Figure2_PRFT_DEG_volcano.png", volcano_plot, width = 7, height = 5.5, dpi = 300)

sig_deg <- deg_all_dt[adj.P.Val < 0.05 & abs(logFC) > 0.5]
if (nrow(sig_deg) > 0) {
  heatmap_gene_n <- min(50, nrow(sig_deg))
  top_heatmap_genes <- sig_deg[order(adj.P.Val, -abs(logFC))][1:heatmap_gene_n, gene_symbol]
  heatmap_mat <- expr_use[top_heatmap_genes, , drop = FALSE]
  heatmap_mat_z <- zscore_rows(heatmap_mat)
  rownames(heatmap_mat_z) <- top_heatmap_genes

  sample_order <- order(group_factor, colnames(expr_use))
  heatmap_mat_z <- heatmap_mat_z[, sample_order, drop = FALSE]
  heatmap_groups <- as.character(group_factor[sample_order])
  names(heatmap_groups) <- colnames(heatmap_mat_z)

  saveRDS(heatmap_mat_z, "05_deg/deg_heatmap_matrix_top50.rds")

  if (requireNamespace("pheatmap", quietly = TRUE)) {
    ann_col <- data.frame(PRFT_group = pretty_label(heatmap_groups), row.names = names(heatmap_groups))
    ann_colors <- list(PRFT_group = c("PRFT low" = "#1F77B4", "PRFT high" = "#D62728"))

    pheatmap::pheatmap(
      heatmap_mat_z,
      annotation_col = ann_col,
      annotation_colors = ann_colors,
      show_colnames = FALSE,
      scale = "none",
      clustering_method = "complete",
      filename = "13_figures/Figure2_PRFT_DEG_heatmap_top50.pdf",
      width = 8,
      height = 10
    )
    pheatmap::pheatmap(
      heatmap_mat_z,
      annotation_col = ann_col,
      annotation_colors = ann_colors,
      show_colnames = FALSE,
      scale = "none",
      clustering_method = "complete",
      filename = "13_figures/Figure2_PRFT_DEG_heatmap_top50.png",
      width = 8,
      height = 10
    )
  } else {
    col_side_colors <- ifelse(heatmap_groups == "PRFT_high", "#D62728", "#1F77B4")

    grDevices::pdf("13_figures/Figure2_PRFT_DEG_heatmap_top50.pdf", width = 8, height = 10)
    stats::heatmap(
      heatmap_mat_z,
      Colv = NA,
      scale = "none",
      col = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
      ColSideColors = col_side_colors,
      labCol = FALSE,
      margins = c(5, 8),
      main = "Top DEGs: PRFT high vs PRFT low"
    )
    grDevices::dev.off()

    grDevices::png("13_figures/Figure2_PRFT_DEG_heatmap_top50.png", width = 2400, height = 3000, res = 300)
    stats::heatmap(
      heatmap_mat_z,
      Colv = NA,
      scale = "none",
      col = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
      ColSideColors = col_side_colors,
      labCol = FALSE,
      margins = c(5, 8),
      main = "Top DEGs: PRFT high vs PRFT low"
    )
    grDevices::dev.off()
  }
} else {
  saveRDS(matrix(numeric(0), nrow = 0, ncol = 0), "05_deg/deg_heatmap_matrix_top50.rds")
  warning("No significant DEGs detected under adj.P.Val < 0.05 and |logFC| > 0.5. Heatmap not generated.")
}

save_session_info("16_logs/sessionInfo_10_deg_prft_high_low_tcga.txt")
