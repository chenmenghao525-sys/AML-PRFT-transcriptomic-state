#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(grid)
  library(svglite)
  library(ragg)
})

options(stringsAsFactors = FALSE)

out_dir <- "Human_Genomics_PRFT_AML_submission_package/graphical_abstract"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("16_logs", recursive = TRUE, showWarnings = FALSE)

font_family <- "Arial"

cols <- list(
  ink = "#1F2933",
  muted = "#5E6A75",
  line = "#6B7280",
  blue = "#315B9A",
  blue_light = "#EAF1FA",
  green = "#2E8B57",
  green_light = "#EAF7EF",
  amber = "#C9881D",
  amber_light = "#FFF6E4",
  orange = "#C75A28",
  orange_light = "#FCEEE8",
  purple = "#6554A6",
  purple_light = "#F0EEFA",
  teal = "#2C8C86",
  teal_light = "#E8F5F3",
  grey_light = "#F6F8FA"
)

text_gp <- function(size = 8, col = cols$ink, face = "plain") {
  gpar(fontfamily = font_family, fontsize = size, col = col, fontface = face)
}

fill_gp <- function(fill = "white", col = cols$line, lwd = 0.8) {
  gpar(fill = fill, col = col, lwd = lwd)
}

draw_box <- function(x, y, w, h, fill, col, r = unit(0.06, "snpc"), lwd = 0.9) {
  grid.roundrect(
    x = unit(x, "native"), y = unit(y, "native"),
    width = unit(w, "native"), height = unit(h, "native"),
    r = r,
    gp = fill_gp(fill, col, lwd)
  )
}

draw_label <- function(label, x, y, size = 8, col = cols$ink, face = "plain", just = "centre", lineheight = 0.92) {
  grid.text(
    label,
    x = unit(x, "native"), y = unit(y, "native"),
    just = just,
    gp = gpar(fontfamily = font_family, fontsize = size, col = col, fontface = face, lineheight = lineheight)
  )
}

draw_arrow <- function(x0, y0, x1, y1, col = cols$line, lwd = 1.2) {
  grid.segments(
    unit(x0, "native"), unit(y0, "native"),
    unit(x1, "native"), unit(y1, "native"),
    arrow = arrow(length = unit(0.12, "inches"), type = "closed"),
    gp = gpar(col = col, lwd = lwd, lineend = "round")
  )
}

draw_dashed <- function(x0, y0, x1, y1, col = cols$line, lwd = 0.85) {
  grid.segments(
    unit(x0, "native"), unit(y0, "native"),
    unit(x1, "native"), unit(y1, "native"),
    gp = gpar(col = col, lwd = lwd, lty = "33", lineend = "round")
  )
}

draw_chip <- function(x, y, w, h, label, fill = "white", col = cols$line, size = 7.2, face = "plain") {
  draw_box(x, y, w, h, fill, col, r = unit(0.08, "snpc"), lwd = 0.7)
  draw_label(label, x, y, size = size, col = cols$ink, face = face)
}

draw_database_icon <- function(x, y, col = cols$blue) {
  for (i in 0:2) {
    grid.roundrect(
      unit(x, "native"), unit(y - i * 0.55, "native"),
      width = unit(3.6, "native"), height = unit(0.55, "native"),
      r = unit(0.035, "snpc"),
      gp = gpar(fill = ifelse(i == 0, cols$blue_light, "white"), col = col, lwd = 0.85)
    )
  }
}

draw_network_icon <- function(x, y, col = cols$blue) {
  pts <- data.frame(dx = c(-1.1, 0, 1.1, -0.4, 0.7), dy = c(0.25, 0.8, 0.25, -0.65, -0.65))
  edges <- list(c(1, 2), c(2, 3), c(1, 4), c(3, 5), c(4, 5), c(2, 5))
  for (e in edges) {
    grid.segments(
      unit(x + pts$dx[e[1]], "native"), unit(y + pts$dy[e[1]], "native"),
      unit(x + pts$dx[e[2]], "native"), unit(y + pts$dy[e[2]], "native"),
      gp = gpar(col = col, lwd = 0.7)
    )
  }
  for (i in seq_len(nrow(pts))) {
    grid.circle(unit(x + pts$dx[i], "native"), unit(y + pts$dy[i], "native"), r = unit(0.15, "native"),
                gp = gpar(fill = "white", col = col, lwd = 0.8))
  }
}

draw_heat_strip_icon <- function(x, y) {
  pal <- c("#2C7FB8", "#FFFFFF", "#D7301F", "#FFFFFF", "#2C7FB8", "#D7301F")
  for (i in seq_along(pal)) {
    grid.rect(
      unit(x - 1.5 + (i - 0.5) * 0.5, "native"), unit(y, "native"),
      width = unit(0.46, "native"), height = unit(1.8, "native"),
      gp = gpar(fill = pal[i], col = NA)
    )
  }
  grid.rect(unit(x, "native"), unit(y, "native"), width = unit(3.1, "native"), height = unit(1.9, "native"),
            gp = gpar(fill = NA, col = cols$green, lwd = 0.7))
}

draw_signature_icon <- function(x, y) {
  draw_chip(x, y + 0.5, 3.7, 0.62, "LASSO-Cox", fill = cols$green_light, col = cols$green, size = 5.8, face = "bold")
  grid.lines(unit(c(x - 1.3, x - 0.4, x + 0.3, x + 1.2), "native"), unit(c(y - 0.45, y - 0.1, y - 0.35, y + 0.05), "native"),
             gp = gpar(col = cols$green, lwd = 1.0))
}

draw_survival_icon <- function(x, y, col = cols$orange) {
  grid.lines(unit(c(x - 1.2, x - 0.4, x - 0.4, x + 0.35, x + 0.35, x + 1.2), "native"),
             unit(c(y + 0.65, y + 0.65, y + 0.25, y + 0.25, y - 0.15, y - 0.15), "native"),
             gp = gpar(col = col, lwd = 1.0))
  grid.lines(unit(c(x - 1.2, x + 1.2), "native"), unit(c(y - 0.55, y - 0.55), "native"), gp = gpar(col = cols$muted, lwd = 0.55))
}

draw_immune_icon <- function(x, y, col = cols$teal) {
  grid.circle(unit(x - 0.45, "native"), unit(y + 0.22, "native"), r = unit(0.45, "native"),
              gp = gpar(fill = cols$teal_light, col = col, lwd = 0.8))
  grid.circle(unit(x + 0.55, "native"), unit(y - 0.18, "native"), r = unit(0.37, "native"),
              gp = gpar(fill = "white", col = col, lwd = 0.8))
  draw_label("PD-L1", x + 0.1, y - 0.9, size = 5.4, col = col, face = "bold")
}

draw_drug_icon <- function(x, y, col = cols$purple) {
  for (i in 0:2) for (j in 0:1) {
    grid.circle(unit(x - 0.85 + i * 0.85, "native"), unit(y + 0.35 - j * 0.75, "native"), r = unit(0.2, "native"),
                gp = gpar(fill = cols$purple_light, col = col, lwd = 0.7))
  }
  grid.lines(unit(c(x - 1.2, x - 0.2, x + 1.2), "native"), unit(c(y - 1.0, y - 0.55, y - 0.1), "native"),
             gp = gpar(col = col, lwd = 0.9))
}

draw_single_cell_icon <- function(x, y, col = cols$orange) {
  set.seed(42)
  px <- c(-0.85, -0.45, -0.15, 0.25, 0.55, 0.9, -0.2, 0.35)
  py <- c(0.55, 0.1, 0.65, 0.28, -0.25, 0.28, -0.5, -0.62)
  for (i in seq_along(px)) {
    grid.circle(unit(x + px[i], "native"), unit(y + py[i], "native"), r = unit(0.13, "native"),
                gp = gpar(fill = ifelse(i %% 2 == 0, cols$orange_light, "white"), col = col, lwd = 0.7))
  }
}

draw_cell_hero <- function(x, y) {
  grid.circle(unit(x, "native"), unit(y, "native"), r = unit(5.3, "native"),
              gp = gpar(fill = cols$amber_light, col = cols$amber, lwd = 1.2))
  grid.circle(unit(x - 0.8, "native"), unit(y - 0.2, "native"), r = unit(1.55, "native"),
              gp = gpar(fill = "white", col = cols$amber, lwd = 0.8))
  draw_network_icon(x - 2.4, y + 1.55, col = cols$blue)
  grid.polygon(
    x = unit(c(x + 1.25, x + 2.55, x + 2.25, x + 1.9, x + 1.55), "native"),
    y = unit(c(y + 1.8, y + 1.8, y + 0.55, y + 0.05, y + 0.55), "native"),
    gp = gpar(fill = cols$orange_light, col = cols$orange, lwd = 0.9)
  )
  draw_label("Fe2+", x + 1.9, y + 1.15, size = 5.5, col = cols$orange, face = "bold")
  draw_label("PRFT\nstate", x, y - 0.2, size = 10.5, col = cols$ink, face = "bold")
  draw_chip(x - 3.55, y - 4.15, 6.2, 1.05, "Proteostasis", fill = "white", col = cols$blue, size = 5.2, face = "bold")
  draw_chip(x + 3.65, y - 4.15, 7.0, 1.05, "Ferroptosis tolerance", fill = "white", col = cols$orange, size = 5.0, face = "bold")
}

draw_graphical_abstract <- function() {
  grid.newpage()
  pushViewport(viewport(xscale = c(0, 100), yscale = c(0, 32), clip = "off"))
  grid.rect(gp = gpar(fill = "white", col = NA))

  draw_label(
    "Proteostasis-associated ferroptosis tolerance defines a myeloid-suppressive AML state",
    50, 30.2, size = 12.8, col = cols$ink, face = "bold"
  )
  draw_label(
    "Cross-platform PRFT-related transcriptomic signature integrating TCGA, GEO, BeatAML and single-cell RNA-seq",
    50, 28.35, size = 7.1, col = cols$muted
  )

  # Module anchors.
  draw_box(13, 17.1, 21.2, 17.0, cols$blue_light, cols$blue, lwd = 1.0)
  draw_label("Program discovery", 13, 24.55, size = 8.8, col = cols$blue, face = "bold")
  draw_database_icon(6.4, 21.1, col = cols$blue)
  draw_label("TCGA-LAML\ntraining cohort", 13.3, 21.0, size = 6.7)
  draw_chip(13, 17.5, 12.7, 1.55, "PRFT score", fill = "white", col = cols$blue, size = 7.2, face = "bold")
  draw_chip(9.2, 14.9, 8.3, 1.35, "Proteostasis\ncore", fill = "white", col = cols$blue, size = 5.8)
  draw_label("+", 13, 14.9, size = 9.2, col = cols$blue, face = "bold")
  draw_chip(16.8, 14.9, 8.8, 1.35, "Ferroptosis\ntolerance set", fill = "white", col = cols$blue, size = 5.8)
  draw_label("GEO validation\nBeatAML\nGSE116256 scRNA-seq", 13, 10.55, size = 6.1, col = cols$muted)

  draw_box(36, 17.1, 19.5, 17.0, cols$green_light, cols$green, lwd = 1.0)
  draw_label("Signature construction", 36, 24.55, size = 8.8, col = cols$green, face = "bold")
  draw_heat_strip_icon(30.5, 21.1)
  draw_label("PRFT high vs low\ntranscriptomic contrast", 38.1, 21.1, size = 6.4)
  draw_signature_icon(30.5, 17.2)
  draw_label("WGCNA modules\nCross-platform filtering", 38.1, 17.1, size = 6.4)
  draw_chip(36, 12.2, 14.1, 1.55, "Final six-gene signature", fill = "white", col = cols$green, size = 6.8, face = "bold")

  draw_cell_hero(57.5, 17.5)
  draw_label("Six-gene PRFT signature", 57.5, 8.3, size = 8.2, col = cols$purple, face = "bold")
  gene_y <- c(6.6, 5.1)
  gene_x <- c(48.7, 57.5, 66.3)
  genes <- matrix(c("CLCN5", "ARHGEF5", "ITGB2", "TRIM32", "SAT1", "ACOX2"), nrow = 2, byrow = TRUE)
  for (r in 1:2) for (c in 1:3) {
    draw_chip(gene_x[c], gene_y[r], 7.2, 1.25, genes[r, c], fill = cols$purple_light, col = cols$purple, size = 6.3, face = "bold")
  }

  draw_box(85, 17.1, 23.0, 17.0, cols$orange_light, cols$orange, lwd = 1.0)
  draw_label("Validation and interpretation", 85, 24.55, size = 8.6, col = cols$orange, face = "bold")

  draw_survival_icon(77.3, 21.3, col = cols$orange)
  draw_label("Reproducible survival\nstratification", 87.2, 21.25, size = 6.2)
  draw_immune_icon(77.3, 17.45, col = cols$teal)
  draw_label("Immune checkpoint and\nmyeloid suppressive features", 88.2, 17.45, size = 5.8)
  draw_drug_icon(77.3, 13.45, col = cols$purple)
  draw_label("BeatAML ex vivo\ndrug-sensitivity associations", 88.4, 13.45, size = 5.8)
  draw_single_cell_icon(77.3, 9.55, col = cols$orange)
  draw_label("Single-cell localization\nto myeloid-like AML states", 88.1, 9.9, size = 5.6)

  draw_arrow(23.6, 17.1, 26.6, 17.1, col = cols$line)
  draw_arrow(45.9, 17.1, 50.8, 17.1, col = cols$line)
  draw_arrow(64.4, 17.1, 73.0, 17.1, col = cols$line)

  draw_box(50, 2.3, 75.0, 2.5, "white", "#D9DEE5", lwd = 0.8)
  draw_label(
    "A cross-platform PRFT-related signature stratifies prognosis and is associated with ferroptosis-tolerance and immune-related transcriptional states in AML.",
    50, 2.3, size = 6.6, col = cols$ink, face = "bold"
  )

  popViewport()
}

export_graphical_abstract <- function() {
  base <- file.path(out_dir, "graphical_abstract_PRFT_AML")

  grDevices::cairo_pdf(paste0(base, ".pdf"), width = 9.2, height = 3.0, family = font_family)
  draw_graphical_abstract()
  dev.off()

  svglite::svglite(paste0(base, ".svg"), width = 9.2, height = 3.0, system_fonts = list(sans = font_family))
  draw_graphical_abstract()
  dev.off()

  ragg::agg_png(paste0(base, ".png"), width = 920, height = 300, units = "px", res = 100, background = "white")
  draw_graphical_abstract()
  dev.off()

  ragg::agg_png(paste0(base, "_highres.png"), width = 2760, height = 900, units = "px", res = 300, background = "white")
  draw_graphical_abstract()
  dev.off()
}

export_graphical_abstract()

readme <- c(
  "Graphical abstract: PRFT-related AML transcriptional state",
  "",
  "This graphical abstract was fully redrawn as a clean horizontal concept figure rather than a reused workflow panel.",
  "Left module: public AML resources and PRFT score definition from the proteostasis core and ferroptosis tolerance set.",
  "Middle module: transcriptomic signature construction using PRFT group contrast, WGCNA modules and cross-platform filtering.",
  "Central module: six-gene PRFT-related signature (CLCN5, ARHGEF5, ITGB2, TRIM32, SAT1 and ACOX2).",
  "Right module: reproducible survival stratification, immune-associated features, BeatAML ex vivo drug-sensitivity associations and single-cell localization.",
  "The figure intentionally avoids causal wording, treatment-response claims, low-confidence embedding labels and archived-model information."
)
writeLines(readme, file.path(out_dir, "graphical_abstract_README.txt"))

checklist <- data.frame(
  item = c(
    "layout",
    "standard_png",
    "highres_png",
    "pdf",
    "svg",
    "forbidden_terms",
    "underscore_variable_names",
    "old_model_information",
    "statistics_added",
    "human_genomics_style"
  ),
  status = c(
    "redrawn horizontal concept layout",
    "920 x 300 px",
    "2760 x 900 px",
    "generated",
    "generated",
    "not present",
    "not present in display text",
    "not present",
    "none",
    "clean white-background graphical abstract"
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(checklist, file.path(out_dir, "graphical_abstract_checklist.csv"), row.names = FALSE)
utils::write.csv(checklist, "Human_Genomics_PRFT_AML_submission_package/graphical_abstract_checklist.csv", row.names = FALSE)

writeLines(capture.output(sessionInfo()), "16_logs/sessionInfo_29_redraw_graphical_abstract_prft_aml.txt")
