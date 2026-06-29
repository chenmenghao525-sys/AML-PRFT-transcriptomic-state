#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(grid)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
figure_dir <- file.path(project_root, "13_figures")
main_figure_dir <- file.path(figure_dir, "main_figures")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(main_figure_dir, showWarnings = FALSE, recursive = TRUE)

save_plot_dual <- function(plot_obj, file_base, width = 12, height = 7) {
  ggsave(paste0(file_base, ".pdf"), plot = plot_obj, width = width, height = height, units = "in", bg = "white")
  ggsave(paste0(file_base, ".png"), plot = plot_obj, width = width, height = height, units = "in", dpi = 300, bg = "white")
}

theme_pub_void <- function() {
  theme_void(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16, color = "#1F2D3A"),
      plot.subtitle = element_text(hjust = 0.5, size = 10.5, color = "#51606F"),
      plot.margin = margin(16, 18, 16, 18),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

draw_box <- function(plot_obj, x, y, label, width, height, fill, border = "#44515E",
                     text_color = "#1F2D3A", text_size = 3.75, fontface = "plain",
                     lineheight = 1.02) {
  rect_df <- data.frame(
    xmin = x - width / 2,
    xmax = x + width / 2,
    ymin = y - height / 2,
    ymax = y + height / 2
  )
  text_df <- data.frame(x = x, y = y, label = label)

  plot_obj +
    geom_rect(
      data = rect_df,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = fill,
      color = border,
      linewidth = 0.45,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = text_df,
      aes(x = x, y = y, label = label),
      color = text_color,
      size = text_size,
      fontface = fontface,
      lineheight = lineheight,
      inherit.aes = FALSE
    )
}

draw_down_arrow <- function(plot_obj, x, y_start, y_end, color = "#4E5D6C", linewidth = 0.55) {
  plot_obj +
    geom_segment(
      data = data.frame(x = x, y = y_start, xend = x, yend = y_end),
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.14, "inches"), type = "closed"),
      linewidth = linewidth,
      color = color,
      lineend = "round",
      inherit.aes = FALSE
    )
}

draw_line <- function(plot_obj, x, y_start, y_end, color = "#8C3D2E", linewidth = 0.5, linetype = "dashed") {
  plot_obj +
    geom_segment(
      data = data.frame(x = x, y = y_start, xend = x, yend = y_end),
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = linewidth,
      color = color,
      linetype = linetype,
      lineend = "round",
      inherit.aes = FALSE
    )
}

# Figure 1A -------------------------------------------------------------------

workflow_base <- ggplot() +
  coord_cartesian(xlim = c(0.35, 4.35), ylim = c(0.2, 4.85), clip = "off") +
  labs(
    title = "Study workflow",
    subtitle = "Computational framework for PRFT program discovery, signature construction, and downstream validation",
    x = NULL,
    y = NULL
  ) +
  theme_pub_void()

col_x <- c(1.00, 2.35, 3.70)
row_y <- c(4.2, 3.2, 2.2, 1.2)
box_h <- 0.36

workflow_labels <- list(
  list(x = col_x[1], y = row_y[1], label = "TCGA-LAML training cohort", width = 0.98, fill = "#DCEAF4"),
  list(x = col_x[1], y = row_y[2], label = "PRFT-related gene set\npreparation", width = 1.02, fill = "#DCEAF4"),
  list(x = col_x[1], y = row_y[3], label = "PRFT score calculation", width = 0.92, fill = "#DCEAF4"),
  list(x = col_x[1], y = row_y[4], label = "Differential expression\nanalysis", width = 0.98, fill = "#DCEAF4"),
  list(x = col_x[2], y = row_y[1], label = "WGCNA", width = 0.74, fill = "#E6F0E3"),
  list(x = col_x[2], y = row_y[2], label = "Candidate gene selection", width = 0.96, fill = "#E6F0E3"),
  list(x = col_x[2], y = row_y[3], label = "Cross-platform\nLASSO-Cox signature construction", width = 1.18, fill = "#E6F0E3"),
  list(x = col_x[2], y = row_y[4], label = "TCGA internal evaluation", width = 0.96, fill = "#E6F0E3"),
  list(x = col_x[3], y = row_y[1], label = "GEO external validation", width = 0.92, fill = "#F6E4DE"),
  list(x = col_x[3], y = row_y[2], label = "Biological\ncharacterization", width = 0.88, fill = "#F6E4DE"),
  list(x = col_x[3], y = row_y[3], label = "Immune microenvironment\nanalysis", width = 1.04, fill = "#F6E4DE"),
  list(x = col_x[3], y = row_y[4], label = "BeatAML drug sensitivity and\nsingle-cell localization analyses", width = 1.26, fill = "#F6E4DE")
)

p_workflow <- workflow_base +
  annotate("text", x = col_x[1], y = 4.68, label = "Discovery", size = 4.8, fontface = "bold", color = "#315B7D") +
  annotate("text", x = col_x[2], y = 4.68, label = "Signature construction", size = 4.8, fontface = "bold", color = "#4B6A3D") +
  annotate("text", x = col_x[3], y = 4.68, label = "Validation and interpretation", size = 4.8, fontface = "bold", color = "#8C3D2E") +
  annotate(
    "text",
    x = 2.35,
    y = 0.42,
    label = "Workflow schematic based on completed analyses only",
    size = 3.7,
    color = "#566573"
  )

for (box in workflow_labels) {
  p_workflow <- draw_box(
    p_workflow,
    x = box$x,
    y = box$y,
    label = box$label,
    width = box$width,
    height = box_h,
    fill = box$fill,
    text_size = 3.7
  )
}

for (x in col_x) {
  for (idx in 1:3) {
    p_workflow <- draw_down_arrow(
      p_workflow,
      x = x,
      y_start = row_y[idx] - box_h / 2 - 0.05,
      y_end = row_y[idx + 1] + box_h / 2 + 0.05
    )
  }
}

# Figure 1B -------------------------------------------------------------------

definition_base <- ggplot() +
  coord_cartesian(xlim = c(0.15, 4.95), ylim = c(0.0, 4.65), clip = "off") +
  labs(
    title = "PRFT conceptual definition",
    subtitle = "PRFT score is defined from the proteostasis core and ferroptosis tolerance set",
    x = NULL,
    y = NULL
  ) +
  theme_pub_void() +
  theme(plot.margin = margin(20, 20, 18, 20))

center_x <- 2.5
center_y <- 2.8
center_h <- 0.62
center_w <- 1.24
top_y <- 4.0
top_h <- 0.36
top_w <- 1.05
bottom_y <- 1.1
bottom_h <- 0.42

p_definition <- definition_base +
  annotate(
    "text",
    x = 2.5,
    y = 0.42,
    label = "Conceptual transcriptomic framework for interpretation only",
    size = 3.7,
    color = "#566573"
  ) +
  annotate(
    "text",
    x = 2.5,
    y = 0.12,
    label = "Auxiliary signatures support biological characterization and do not imply a proven causal mechanism",
    size = 3.25,
    color = "#566573"
  )

p_definition <- draw_box(
  p_definition,
  x = 1.5,
  y = top_y,
  label = "Proteostasis core",
  width = top_w,
  height = top_h,
  fill = "#DCEAF4",
  text_size = 3.8
)
p_definition <- draw_box(
  p_definition,
  x = 3.5,
  y = top_y,
  label = "Ferroptosis tolerance set",
  width = 1.18,
  height = top_h,
  fill = "#DCEAF4",
  text_size = 3.8
)
p_definition <- draw_box(
  p_definition,
  x = center_x,
  y = center_y,
  label = "PRFT score\nMean of z-scored core signatures",
  width = center_w,
  height = center_h,
  fill = "#F3E6B3",
  border = "#8A6D1D",
  text_color = "#6B520F",
  text_size = 4.0,
  fontface = "bold",
  lineheight = 1.05
)

bottom_boxes <- list(
  list(x = 0.9, label = "SLC7A11/GPX4-GSH axis", width = 1.08),
  list(x = 2.0, label = "JAK2/STAT5/PD-L1", width = 0.98),
  list(x = 3.1, label = "Immune checkpoint and\nmyeloid suppressive features", width = 1.18),
  list(x = 4.2, label = "Stemness and\nLSC/MRD-like related features", width = 1.12)
)

for (box in bottom_boxes) {
  p_definition <- draw_box(
    p_definition,
    x = box$x,
    y = bottom_y,
    label = box$label,
    width = box$width,
    height = bottom_h,
    fill = "#F6E4DE",
    text_size = 3.45
  )
}

p_definition <- p_definition +
  geom_segment(
    data = data.frame(
      x = c(1.5, 3.5),
      y = c(top_y - top_h / 2 - 0.04, top_y - top_h / 2 - 0.04),
      xend = c(center_x - 0.36, center_x + 0.36),
      yend = c(center_y + center_h / 2 + 0.05, center_y + center_h / 2 + 0.05)
    ),
    aes(x = x, y = y, xend = xend, yend = yend),
    arrow = arrow(length = unit(0.14, "inches"), type = "closed"),
    linewidth = 0.55,
    color = "#4E5D6C",
    lineend = "round",
    inherit.aes = FALSE
  )

bottom_start_x <- c(1.98, 2.28, 2.72, 3.02)
bottom_end_x <- vapply(bottom_boxes, function(box) box$x, numeric(1))
for (i in seq_along(bottom_start_x)) {
  p_definition <- p_definition +
    geom_segment(
      data = data.frame(
        x = bottom_start_x[i],
        y = center_y - center_h / 2 - 0.05,
        xend = bottom_end_x[i],
        yend = bottom_y + bottom_h / 2 + 0.07
      ),
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.5,
      color = "#8C3D2E",
      linetype = "dashed",
      lineend = "round",
      inherit.aes = FALSE
    )
}

# Combined --------------------------------------------------------------------

p_workflow_tagged <- p_workflow + labs(tag = "A") +
  theme(
    plot.tag = element_text(face = "bold", size = 16, color = "#1F2D3A"),
    plot.tag.position = c(0.01, 0.99),
    plot.margin = margin(18, 18, 42, 18)
  )

p_definition_tagged <- p_definition + labs(tag = "B") +
  theme(
    plot.tag = element_text(face = "bold", size = 16, color = "#1F2D3A"),
    plot.tag.position = c(0.01, 0.99),
    plot.margin = margin(44, 18, 18, 18)
  )

p_combined <- p_workflow_tagged / p_definition_tagged +
  plot_layout(heights = c(1, 0.92))

save_plot_dual(p_workflow, file.path(figure_dir, "Figure1A_workflow_schematic"), width = 12, height = 7.6)
save_plot_dual(p_definition, file.path(figure_dir, "Figure1B_PRFT_definition_schematic"), width = 12, height = 6.8)
save_plot_dual(p_combined, file.path(figure_dir, "Figure1_workflow_and_PRFT_definition_combined"), width = 12, height = 13.6)

file.copy(
  c(
    file.path(figure_dir, "Figure1A_workflow_schematic.pdf"),
    file.path(figure_dir, "Figure1A_workflow_schematic.png"),
    file.path(figure_dir, "Figure1B_PRFT_definition_schematic.pdf"),
    file.path(figure_dir, "Figure1B_PRFT_definition_schematic.png"),
    file.path(figure_dir, "Figure1_workflow_and_PRFT_definition_combined.pdf"),
    file.path(figure_dir, "Figure1_workflow_and_PRFT_definition_combined.png")
  ),
  main_figure_dir,
  overwrite = TRUE
)

message("Figure 1A and Figure 1B schematic files generated successfully.")
