#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
figure_dir <- file.path(project_root, "13_figures")
main_figure_dir <- file.path(figure_dir, "main_figures")
table_dir <- file.path(project_root, "14_tables")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(main_figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

coef_file <- file.path(project_root, "00_LOCKED_FORMULA", "LOCKED_PRFT_six_gene_formula_A_coefficients.csv")
if (!file.exists(coef_file)) {
  stop("Missing required coefficient file: ", coef_file)
}

coef_dt <- fread(coef_file)
colnames(coef_dt) <- tolower(colnames(coef_dt))
if (!all(c("gene_symbol", "coefficient") %in% colnames(coef_dt))) {
  stop("Coefficient file must contain gene_symbol and coefficient columns.")
}

coef_dt <- coef_dt[, .(gene_symbol = as.character(gene_symbol), coefficient = as.numeric(coefficient))]
coef_dt <- coef_dt[order(-coefficient)]
coef_dt[, coefficient_round := sprintf("%.8f", coefficient)]

out_table <- coef_dt[, .(
  gene_symbol,
  coefficient = coefficient_round
)]
fwrite(out_table, file.path(table_dir, "Figure3D_signature_coefficients_formula_table.csv"))

save_plot_dual <- function(plot_obj, file_base, width = 11, height = 7.8) {
  ggsave(paste0(file_base, ".pdf"), plot = plot_obj, width = width, height = height, units = "in")
  ggsave(paste0(file_base, ".png"), plot = plot_obj, width = width, height = height, units = "in", dpi = 300)
}

coef_plot_dt <- copy(coef_dt)
coef_plot_dt[, gene_symbol := factor(gene_symbol, levels = rev(gene_symbol))]

p_bar <- ggplot(coef_plot_dt, aes(x = gene_symbol, y = coefficient)) +
  geom_col(fill = "#5B8DB8", width = 0.68) +
  geom_text(aes(label = coefficient_round), hjust = -0.08, size = 3.7, color = "#1F2D3A") +
  coord_flip(clip = "off") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(
    title = "Cross-platform six-gene PRFT-related signature",
    subtitle = "Gene coefficients derived from the fixed cross-platform model",
    x = NULL,
    y = "Gene coefficient"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", color = "#1F2D3A"),
    plot.subtitle = element_text(hjust = 0.5, color = "#4E5D6C"),
    axis.text.y = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(12, 20, 12, 12)
  )

formula_coef_dt <- coef_dt[order(match(gene_symbol, c("CLCN5", "ARHGEF5", "TRIM32", "ITGB2", "SAT1", "ACOX2")))]
formula_lines <- c(
  "risk score =",
  paste0(
    sprintf("%.10f", formula_coef_dt$coefficient),
    " x z(",
    formula_coef_dt$gene_symbol,
    ")",
    c(rep(" +", nrow(formula_coef_dt) - 1), "")
  )
)

formula_df <- data.frame(
  x = 0.5,
  y = seq(0.88, 0.34, length.out = length(formula_lines)),
  label = formula_lines,
  stringsAsFactors = FALSE
)

p_formula <- ggplot() +
  annotate("rect", xmin = 0.06, xmax = 0.94, ymin = 0.14, ymax = 0.95, fill = "#F7F7F7", color = "#BCC7D1", linewidth = 0.5) +
  geom_text(
    data = formula_df,
    aes(x = x, y = y, label = label),
    family = "mono",
    size = c(5.2, rep(4.2, length(formula_lines) - 1)),
    fontface = c("bold", rep("plain", length(formula_lines) - 1)),
    color = "#1F2D3A",
    lineheight = 1.05
  ) +
  annotate(
    "text",
    x = 0.5,
    y = 0.08,
    label = "Gene-wise z-score was applied within each cohort before using the fixed coefficients.",
    size = 3.7,
    color = "#4E5D6C"
  ) +
  annotate(
    "text",
    x = 0.5,
    y = 0.02,
    label = "This panel summarizes the fixed PRFT-related signature and does not re-estimate the model.",
    size = 3.4,
    color = "#4E5D6C"
  ) +
  labs(
    title = "Fixed risk score formula",
    subtitle = "Formula applied unchanged across validation cohorts",
    x = NULL,
    y = NULL
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", color = "#1F2D3A"),
    plot.subtitle = element_text(hjust = 0.5, color = "#4E5D6C"),
    plot.margin = margin(12, 16, 14, 16)
  )

p_combined <- p_bar + p_formula + plot_layout(widths = c(1.1, 1))

out_base <- file.path(figure_dir, "Figure3D_cross_platform_signature_coefficients_formula")
save_plot_dual(p_combined, out_base, width = 13, height = 7.8)

file.copy(paste0(out_base, ".pdf"), file.path(main_figure_dir, "Figure3D_cross_platform_signature_coefficients_formula.pdf"), overwrite = TRUE)
file.copy(paste0(out_base, ".png"), file.path(main_figure_dir, "Figure3D_cross_platform_signature_coefficients_formula.png"), overwrite = TRUE)

message("Figure 3D coefficient and formula figure generated successfully.")
message("Coefficient file used: ", coef_file)
