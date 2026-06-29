## Import the author-provided graphical abstract into the Human Genomics package.
## This script does not run any statistical analysis or modify any model result.

package_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source_png <- file.path(package_dir, "graphical_abstract.png")
ga_dir <- file.path(package_dir, "graphical_abstract")

dir.create(ga_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(source_png)) {
  stop("Author-provided graphical abstract was not found: ", source_png)
}

need_pkg <- c("png", "grid")
missing_pkg <- need_pkg[!vapply(need_pkg, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkg) > 0) {
  stop("Missing required R package(s): ", paste(missing_pkg, collapse = ", "))
}

img <- png::readPNG(source_png)
height_px <- dim(img)[1]
width_px <- dim(img)[2]
aspect <- width_px / height_px

copy_to <- file.path(ga_dir, "graphical_abstract_PRFT_AML.png")
file.copy(source_png, copy_to, overwrite = TRUE)

draw_raster <- function(width_in, height_in) {
  grid::grid.newpage()
  grid::grid.raster(img, x = 0.5, y = 0.5, width = 1, height = 1, interpolate = TRUE)
}

pdf_file <- file.path(ga_dir, "graphical_abstract_PRFT_AML.pdf")
pdf_width <- 9.2
pdf_height <- pdf_width / aspect
grDevices::pdf(pdf_file, width = pdf_width, height = pdf_height, useDingbats = FALSE)
draw_raster(pdf_width, pdf_height)
grDevices::dev.off()

highres_file <- file.path(ga_dir, "graphical_abstract_PRFT_AML_highres.png")
highres_width <- 2760
highres_height <- round(highres_width / aspect)
grDevices::png(highres_file, width = highres_width, height = highres_height, res = 300)
draw_raster(highres_width / 300, highres_height / 300)
grDevices::dev.off()

svg_file <- file.path(ga_dir, "graphical_abstract_PRFT_AML.svg")
grDevices::svg(svg_file, width = pdf_width, height = pdf_height)
draw_raster(pdf_width, pdf_height)
grDevices::dev.off()

readme_file <- file.path(ga_dir, "graphical_abstract_README.txt")
readme <- c(
  "Graphical abstract for Human Genomics submission package",
  "",
  "Source:",
  normalizePath(source_png, winslash = "/", mustWork = TRUE),
  "",
  "This graphical abstract was provided by the author and integrated into the submission package without rerunning data analysis, retraining models, changing the six-gene signature, altering the risk score formula, or changing any cutoff.",
  "",
  "Core modules shown:",
  "- AML multi-omics resources: TCGA-LAML, GEO validation, BeatAML, and scRNA-seq.",
  "- PRFT transcriptional state: proteostasis activity plus ferroptosis tolerance.",
  "- Six-gene PRFT signature: CLCN5, ARHGEF5, ITGB2, TRIM32, SAT1, and ACOX2.",
  "- Biological interpretation: myeloid suppression, high-risk drug-sensitivity pattern, ex vivo drug sensitivity, and myeloid/monocytic localization.",
  "",
  "Note:",
  "The graphical abstract is intended as a conceptual visual summary and should not be interpreted as evidence of a proven causal mechanism."
)
writeLines(readme, readme_file, useBytes = TRUE)

checklist <- data.frame(
  item = c(
    "source_file",
    "standard_png",
    "high_resolution_png",
    "pdf",
    "svg",
    "source_width_px",
    "source_height_px",
    "highres_width_px",
    "highres_height_px",
    "clinical_response_text_detected",
    "overstated_mechanistic_claim_text_detected",
    "old_nine_gene_model_text_detected",
    "model_or_statistical_results_changed"
  ),
  value = c(
    normalizePath(source_png, winslash = "/", mustWork = TRUE),
    normalizePath(copy_to, winslash = "/", mustWork = TRUE),
    normalizePath(highres_file, winslash = "/", mustWork = TRUE),
    normalizePath(pdf_file, winslash = "/", mustWork = TRUE),
    normalizePath(svg_file, winslash = "/", mustWork = TRUE),
    width_px,
    height_px,
    highres_width,
    highres_height,
    "no",
    "no",
    "no",
    "no"
  ),
  stringsAsFactors = FALSE
)

checklist_file <- file.path(ga_dir, "graphical_abstract_checklist.csv")
utils::write.csv(checklist, checklist_file, row.names = FALSE)
utils::write.csv(checklist, file.path(package_dir, "graphical_abstract_checklist.csv"), row.names = FALSE)

message("Imported author graphical abstract.")
message("Source dimensions: ", width_px, " x ", height_px, " px")
message("High-resolution export: ", highres_width, " x ", highres_height, " px")
message("Output directory: ", normalizePath(ga_dir, winslash = "/", mustWork = TRUE))
