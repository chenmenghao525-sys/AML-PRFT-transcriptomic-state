#!/usr/bin/env Rscript

# TCGA-LAML preprocessing script
# - Reads raw SummarizedExperiment from 01_download_tcga_laml.R
# - Extracts raw count matrix
# - Converts Ensembl IDs to HGNC gene symbols
# - Removes genes without symbols
# - Resolves duplicated gene symbols by keeping the row with the highest mean raw count
#   This is more conservative than averaging raw counts across duplicated symbols and preserves
#   count-based assumptions before TMM normalization.
# - Filters low-expression genes
# - Performs TMM normalization and converts to log2 CPM
# - Harmonizes sample_id / patient_id with clinical data
# - Derives OS_time and OS_status using compatible clinical field lookup

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(S4Vectors)
  library(data.table)
  library(edgeR)
  library(dplyr)
  library(stringr)
  library(AnnotationDbi)
})

options(stringsAsFactors = FALSE)

input_se <- file.path("00_raw_data", "tcga_laml_raw", "tcga_laml_se.rds")
input_clin <- file.path("01_metadata", "clinical_tcga_raw.csv")

output_expr <- file.path("02_processed_data", "tcga_expr_hgnc_log2cpm.rds")
output_clin <- file.path("02_processed_data", "tcga_clin_clean.rds")
output_matched <- file.path("02_processed_data", "tcga_expr_clin_matched.rds")
output_summary <- file.path("14_tables", "tcga_sample_summary.csv")
output_session <- file.path("16_logs", "sessionInfo_04_preprocess_tcga_laml.txt")
output_gene_annot <- file.path("01_metadata", "tcga_gene_annotation_hgnc.csv")
output_sample_map <- file.path("01_metadata", "tcga_sample_mapping.csv")
output_duplicate_samples <- file.path("14_tables", "tcga_duplicate_samples.csv")

dir.create("02_processed_data", recursive = TRUE, showWarnings = FALSE)
dir.create("01_metadata", recursive = TRUE, showWarnings = FALSE)
dir.create("14_tables", recursive = TRUE, showWarnings = FALSE)
dir.create("16_logs", recursive = TRUE, showWarnings = FALSE)

save_session_info <- function(path) {
  writeLines(capture.output(sessionInfo()), con = path)
}

pick_first_existing <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) {
    return(rep(NA, nrow(df)))
  }
  df[[hit[1]]]
}

pick_existing_names <- function(df, candidates) {
  candidates[candidates %in% colnames(df)]
}

pick_preferred_sample_per_patient <- function(sample_ids) {
  sample_type_code <- substr(sample_ids, 14, 15)
  sample_type_rank <- dplyr::case_when(
    sample_type_code == "01" ~ 1L,
    sample_type_code == "03" ~ 2L,
    sample_type_code == "06" ~ 3L,
    TRUE ~ 9L
  )

  data.frame(
    sample_id = sample_ids,
    patient_id = substr(sample_ids, 1, 12),
    sample_type_code = sample_type_code,
    sample_type_rank = sample_type_rank,
    stringsAsFactors = FALSE
  ) %>%
    arrange(.data$patient_id, .data$sample_type_rank, .data$sample_id)
}

coalesce_columns <- function(df, candidates) {
  out <- rep(NA, nrow(df))
  hits <- candidates[candidates %in% colnames(df)]
  if (length(hits) == 0) {
    return(out)
  }
  for (nm in hits) {
    idx <- is.na(out) | out == ""
    out[idx] <- df[[nm]][idx]
  }
  out
}

as_numeric_safely <- function(x) {
  suppressWarnings(as.numeric(x))
}

if (!file.exists(input_se)) {
  stop("Missing input file: ", input_se)
}

if (!file.exists(input_clin)) {
  stop("Missing input file: ", input_clin)
}

message("Reading SummarizedExperiment")
se <- readRDS(input_se)
clinical_raw <- fread(input_clin, data.table = FALSE)

if (!inherits(se, "SummarizedExperiment")) {
  stop("Input object is not a SummarizedExperiment.")
}

available_assays <- assayNames(se)
message("Available assays: ", paste(available_assays, collapse = ", "))

preferred_assays <- c("unstranded", "htseq_counts", "counts")
assay_to_use <- preferred_assays[preferred_assays %in% available_assays][1]
if (is.na(assay_to_use) || length(assay_to_use) == 0) {
  assay_to_use <- available_assays[1]
  message("Preferred count assay not found; using first assay: ", assay_to_use)
}

expr_raw <- assay(se, assay_to_use)
if (!is.matrix(expr_raw) && !inherits(expr_raw, "Matrix")) {
  expr_raw <- as.matrix(expr_raw)
}

mode(expr_raw) <- "numeric"

sample_ids_full <- colnames(expr_raw)
sample_ids <- substr(sample_ids_full, 1, 16)
patient_ids <- substr(sample_ids_full, 1, 12)
colnames(expr_raw) <- sample_ids

row_df <- as.data.frame(rowData(se))
rownames(row_df) <- rownames(se)

ensembl_id <- rownames(expr_raw)
ensembl_id_clean <- sub("\\..*$", "", ensembl_id)

symbol_candidates <- c("gene_name", "external_gene_name", "symbol", "hgnc_symbol")
symbol_from_rowdata_cols <- symbol_candidates[symbol_candidates %in% colnames(row_df)]

if (length(symbol_from_rowdata_cols) > 0) {
  gene_symbol <- row_df[[symbol_from_rowdata_cols[1]]]
  message("Using gene symbol column from rowData: ", symbol_from_rowdata_cols[1])
} else {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("No gene symbol column found in rowData and org.Hs.eg.db is not installed.")
  }
  message("Gene symbol not found in rowData; mapping via org.Hs.eg.db")
  gene_symbol <- mapIds(
    x = getNamespace("org.Hs.eg.db"),
    keys = ensembl_id_clean,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
}

gene_annot <- data.frame(
  ensembl_id = ensembl_id,
  ensembl_id_clean = ensembl_id_clean,
  gene_symbol = unname(gene_symbol),
  stringsAsFactors = FALSE
)

keep_symbol <- !is.na(gene_annot$gene_symbol) & gene_annot$gene_symbol != ""
expr_symbol <- expr_raw[keep_symbol, , drop = FALSE]
gene_annot <- gene_annot[keep_symbol, , drop = FALSE]

row_mean_raw <- rowMeans(expr_symbol, na.rm = TRUE)
gene_annot$row_mean_raw <- row_mean_raw

gene_annot <- gene_annot %>%
  arrange(.data$gene_symbol, desc(.data$row_mean_raw))

expr_symbol <- expr_symbol[match(gene_annot$ensembl_id, rownames(expr_symbol)), , drop = FALSE]

dedup_index <- !duplicated(gene_annot$gene_symbol)
expr_dedup <- expr_symbol[dedup_index, , drop = FALSE]
gene_annot_dedup <- gene_annot[dedup_index, , drop = FALSE]
rownames(expr_dedup) <- gene_annot_dedup$gene_symbol

message("Genes after symbol mapping and deduplication: ", nrow(expr_dedup))
fwrite(as.data.table(gene_annot_dedup), file = output_gene_annot)

sample_map_all <- pick_preferred_sample_per_patient(sample_ids) %>%
  mutate(keep_sample = !duplicated(.data$patient_id))

duplicate_sample_map <- sample_map_all %>%
  group_by(.data$patient_id) %>%
  filter(n() > 1) %>%
  ungroup()

if (nrow(duplicate_sample_map) > 0) {
  message("Detected duplicated patient_id entries; keeping one sample per patient.")
  fwrite(as.data.table(duplicate_sample_map), file = output_duplicate_samples)
} else {
  fwrite(as.data.table(duplicate_sample_map), file = output_duplicate_samples)
}

sample_map_kept <- sample_map_all %>%
  filter(.data$keep_sample) %>%
  select(sample_id, patient_id, sample_type_code, sample_type_rank)

fwrite(as.data.table(sample_map_all), file = output_sample_map)

expr_dedup <- expr_dedup[, sample_map_kept$sample_id, drop = FALSE]
sample_ids <- colnames(expr_dedup)
patient_ids <- substr(sample_ids, 1, 12)

group_factor <- factor(rep("all", ncol(expr_dedup)))
dge <- DGEList(counts = expr_dedup)

keep_expr <- filterByExpr(dge, group = group_factor)
dge <- dge[keep_expr, , keep.lib.sizes = FALSE]

message("Genes after low-expression filtering: ", nrow(dge))

dge <- calcNormFactors(dge, method = "TMM")
expr_log2cpm <- cpm(dge, log = TRUE, prior.count = 1)

expr_range <- range(expr_log2cpm, na.rm = TRUE)
message("log2 CPM range: ", paste(round(expr_range, 3), collapse = " to "))

clinical_raw <- as.data.frame(clinical_raw)

if ("submitter_id" %in% colnames(clinical_raw)) {
  clinical_raw$patient_id <- substr(clinical_raw$submitter_id, 1, 12)
} else if ("case_submitter_id" %in% colnames(clinical_raw)) {
  clinical_raw$patient_id <- substr(clinical_raw$case_submitter_id, 1, 12)
} else {
  clinical_raw$patient_id <- substr(coalesce_columns(clinical_raw, c("bcr_patient_barcode", "case_id", "submitter_id")), 1, 12)
}

clinical_raw$patient_id <- toupper(clinical_raw$patient_id)
patient_ids <- toupper(patient_ids)
colnames(expr_log2cpm) <- sample_ids

os_status_raw <- toupper(as.character(coalesce_columns(
  clinical_raw,
  c("vital_status", "overall_survival_status", "os_status")
)))

days_to_death <- as_numeric_safely(coalesce_columns(
  clinical_raw,
  c("days_to_death", "death_days_to", "days_to_death_demographic")
))
days_to_last_follow_up <- as_numeric_safely(coalesce_columns(
  clinical_raw,
  c("days_to_last_follow_up", "days_to_last_followup", "days_to_last_known_alive", "last_contact_days_to")
))

os_time <- ifelse(!is.na(days_to_death), days_to_death, days_to_last_follow_up)
os_status <- ifelse(os_status_raw %in% c("DEAD", "1"), 1L,
                    ifelse(os_status_raw %in% c("ALIVE", "0"), 0L, NA_integer_))

age_years <- as_numeric_safely(coalesce_columns(
  clinical_raw,
  c("age_at_index", "age_at_diagnosis", "age_at_initial_pathologic_diagnosis")
))

if (all(is.na(age_years)) && "days_to_birth" %in% colnames(clinical_raw)) {
  age_years <- abs(as_numeric_safely(clinical_raw$days_to_birth)) / 365.25
}

if (!all(is.na(age_years)) && stats::median(age_years, na.rm = TRUE) > 120) {
  message("Detected age field likely stored in days; converting to years.")
  age_years <- age_years / 365.25
}

sex <- coalesce_columns(clinical_raw, c("gender", "sex"))
fab <- coalesce_columns(
  clinical_raw,
  c("french_american_british_classification", "fab_classification", "morphology", "diagnosis")
)
wbc <- as_numeric_safely(coalesce_columns(
  clinical_raw,
  c("white_blood_cell_count", "wbc", "wbc_at_diagnosis")
))

mutation_patterns <- c(
  "FLT3", "NPM1", "TP53", "DNMT3A", "IDH1", "IDH2", "RUNX1",
  "TET2", "CEBPA", "KIT", "WT1", "NRAS", "KRAS", "ASXL1"
)
mutation_cols <- grep(
  paste(mutation_patterns, collapse = "|"),
  colnames(clinical_raw),
  ignore.case = TRUE,
  value = TRUE
)

clin_clean <- data.frame(
  patient_id = clinical_raw$patient_id,
  sample_id = NA_character_,
  OS_time = os_time,
  OS_status = os_status,
  age = age_years,
  sex = sex,
  FAB = fab,
  WBC = wbc,
  stringsAsFactors = FALSE
)

if (length(mutation_cols) > 0) {
  clin_clean <- cbind(clin_clean, clinical_raw[, mutation_cols, drop = FALSE])
} else {
  message("No mutation-related columns detected in raw clinical file; mutation fields will not be added.")
}

clin_clean <- clin_clean %>%
  filter(!is.na(patient_id) & patient_id != "") %>%
  distinct(patient_id, .keep_all = TRUE)

sample_map <- sample_map_kept %>%
  select(sample_id, patient_id)

clin_matched <- sample_map %>%
  left_join(clin_clean, by = "patient_id") %>%
  mutate(sample_id = .data$sample_id.x) %>%
  select(-sample_id.x, -sample_id.y)

matched_sample_ids <- intersect(colnames(expr_log2cpm), clin_matched$sample_id)
expr_matched <- expr_log2cpm[, matched_sample_ids, drop = FALSE]
clin_matched <- clin_matched %>%
  filter(sample_id %in% matched_sample_ids) %>%
  arrange(match(sample_id, matched_sample_ids))

if (!identical(colnames(expr_matched), clin_matched$sample_id)) {
  stop("Expression matrix columns and clinical sample_id could not be aligned.")
}

samples_with_os_time <- sum(!is.na(clin_matched$OS_time))
samples_with_os_status <- sum(!is.na(clin_matched$OS_status))
samples_with_complete_survival <- sum(!is.na(clin_matched$OS_time) & !is.na(clin_matched$OS_status))

sample_summary <- data.frame(
  metric = c(
    "raw_genes",
    "genes_with_symbol",
    "genes_after_dedup",
    "genes_after_low_expression_filter",
    "raw_samples",
    "samples_after_patient_dedup",
    "matched_samples",
    "matched_patients",
    "samples_with_OS_time",
    "samples_with_OS_status",
    "samples_with_complete_survival"
  ),
  value = c(
    nrow(expr_raw),
    sum(keep_symbol),
    nrow(expr_dedup),
    nrow(expr_log2cpm),
    ncol(expr_raw),
    nrow(sample_map_kept),
    ncol(expr_matched),
    length(unique(clin_matched$patient_id)),
    samples_with_os_time,
    samples_with_os_status,
    samples_with_complete_survival
  ),
  stringsAsFactors = FALSE
)

saveRDS(expr_log2cpm, file = output_expr)
saveRDS(clin_matched, file = output_clin)
saveRDS(
  list(
    expr = expr_matched,
    clin = clin_matched
  ),
  file = output_matched
)
fwrite(as.data.table(sample_summary), file = output_summary)

save_session_info(output_session)
message("Saved processed expression, cleaned clinical data, matched object, summary table, and sessionInfo.")
