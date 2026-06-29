#!/usr/bin/env Rscript

# Preprocess TCGA-LAML expression and clinical/survival files downloaded from UCSC Xena / GDC Xena Hub.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(dplyr)
})

options(stringsAsFactors = FALSE)

raw_dir <- file.path("00_raw_data", "tcga_laml_raw_xena")
manifest_file <- file.path(raw_dir, "tcga_laml_xena_manifest.csv")

output_expr <- file.path("02_processed_data", "tcga_expr_hgnc_log2cpm.rds")
output_clin <- file.path("02_processed_data", "tcga_clin_clean.rds")
output_matched <- file.path("02_processed_data", "tcga_expr_clin_matched.rds")
output_summary <- file.path("14_tables", "tcga_sample_summary.csv")
output_gene_annot <- file.path("01_metadata", "tcga_gene_annotation_hgnc.csv")
output_sample_map <- file.path("01_metadata", "tcga_sample_mapping.csv")
output_duplicate_samples <- file.path("14_tables", "tcga_duplicate_samples.csv")
output_session <- file.path("16_logs", "sessionInfo_04_preprocess_tcga_laml_xena.txt")

dir.create("01_metadata", recursive = TRUE, showWarnings = FALSE)
dir.create("02_processed_data", recursive = TRUE, showWarnings = FALSE)
dir.create("14_tables", recursive = TRUE, showWarnings = FALSE)
dir.create("16_logs", recursive = TRUE, showWarnings = FALSE)

save_session_info <- function(path) {
  writeLines(capture.output(sessionInfo()), con = path)
}

coalesce_columns <- function(df, candidates) {
  hits <- candidates[candidates %in% colnames(df)]
  out <- rep(NA, nrow(df))
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

derive_patient_id <- function(sample_id, patient_raw) {
  patient_norm <- toupper(trimws(as.character(patient_raw)))
  sample_norm <- toupper(trimws(as.character(sample_id)))

  out <- ifelse(
    !is.na(patient_norm) & grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", patient_norm),
    substr(patient_norm, 1, 12),
    NA_character_
  )

  fallback_idx <- is.na(out) | out == ""
  out[fallback_idx] <- ifelse(
    !is.na(sample_norm[fallback_idx]) & grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", sample_norm[fallback_idx]),
    substr(sample_norm[fallback_idx], 1, 12),
    out[fallback_idx]
  )

  out
}

is_sample_like <- function(x) {
  grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", x, ignore.case = TRUE)
}

parse_gene_symbol <- function(ids) {
  ids <- as.character(ids)
  ids <- sub("\\..*$", "", ids)
  out <- ids

  has_pipe <- grepl("\\|", ids)
  if (any(has_pipe)) {
    split_ids <- strsplit(ids[has_pipe], "\\|", fixed = FALSE)
    parsed <- vapply(split_ids, function(tokens) {
      tokens <- tokens[tokens != "" & !is.na(tokens)]
      non_ensg <- tokens[!grepl("^ENSG", tokens, ignore.case = TRUE)]
      if (length(non_ensg) > 0) {
        non_ensg[1]
      } else {
        tokens[1]
      }
    }, character(1))
    out[has_pipe] <- parsed
  }

  out <- sub("\\..*$", "", out)
  out <- trimws(out)
  out[out == "" | grepl("^ENSG", out, ignore.case = TRUE)] <- NA_character_
  out
}

build_sample_mapping <- function(sample_ids) {
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
    arrange(.data$patient_id, .data$sample_type_rank, .data$sample_id) %>%
    mutate(keep_sample = !duplicated(.data$patient_id))
}

read_table_auto <- function(path) {
  fread(path, data.table = FALSE)
}

if (!file.exists(manifest_file)) {
  stop("Missing manifest file: ", manifest_file)
}

manifest <- fread(manifest_file, data.table = FALSE)
expr_row <- manifest[manifest$role == "expression" & manifest$success %in% c(TRUE, "TRUE"), , drop = FALSE]
if (nrow(expr_row) == 0) {
  expr_files <- list.files(raw_dir, pattern = "HiSeqV2|htseq|fpkm", full.names = TRUE)
  if (length(expr_files) == 0) {
    stop("No successful expression download recorded in manifest and no expression file found in raw_dir.")
  }
  expr_row <- data.frame(
    role = "expression",
    source = "TCGA Xena Hub",
    data_type = "normalized",
    dest = expr_files[1],
    stringsAsFactors = FALSE
  )
}

pheno_row <- manifest[manifest$role == "phenotype" & manifest$success %in% c(TRUE, "TRUE"), , drop = FALSE]
surv_row <- manifest[manifest$role == "survival" & manifest$success %in% c(TRUE, "TRUE"), , drop = FALSE]

if (nrow(pheno_row) == 0) {
  pheno_files <- list.files(raw_dir, pattern = "clinicalMatrix|phenotype", full.names = TRUE)
  if (length(pheno_files) > 0) {
    pheno_row <- data.frame(
      role = "phenotype",
      source = "TCGA Xena Hub",
      data_type = "clinicalMatrix",
      dest = pheno_files[1],
      stringsAsFactors = FALSE
    )
  }
}

if (nrow(surv_row) == 0) {
  surv_files <- list.files(raw_dir, pattern = "survival", full.names = TRUE)
  if (length(surv_files) > 0) {
    surv_row <- data.frame(
      role = "survival",
      source = "Xena",
      data_type = "survival",
      dest = surv_files[1],
      stringsAsFactors = FALSE
    )
  }
}

expr_path <- expr_row$dest[1]
if (!file.exists(expr_path)) {
  stop("Expression file recorded in manifest does not exist: ", expr_path)
}

expression_data_type <- expr_row$data_type[1]
expression_source <- expr_row$source[1]
is_raw_counts <- identical(expression_data_type, "HTSeq counts")

message("Reading expression file: ", expr_path)
expr_df <- read_table_auto(expr_path)
if (ncol(expr_df) < 2) {
  stop("Expression table has fewer than 2 columns.")
}

first_col <- colnames(expr_df)[1]
first_values <- as.character(expr_df[[1]])
sample_like_first_col <- mean(is_sample_like(first_values), na.rm = TRUE)
sample_like_colnames <- mean(is_sample_like(colnames(expr_df)[-1]), na.rm = TRUE)

if (sample_like_first_col > 0.5) {
  message("Detected sample-by-gene orientation; transposing to gene-by-sample.")
  expr_mat0 <- as.matrix(expr_df[, -1, drop = FALSE])
  rownames(expr_mat0) <- expr_df[[1]]
  expr_mat <- t(expr_mat0)
  storage.mode(expr_mat) <- "numeric"
  gene_ids <- rownames(expr_mat)
  sample_ids_raw <- colnames(expr_mat)
} else if (sample_like_colnames > 0.5) {
  message("Detected gene-by-sample orientation.")
  expr_mat <- as.matrix(expr_df[, -1, drop = FALSE])
  rownames(expr_mat) <- expr_df[[1]]
  storage.mode(expr_mat) <- "numeric"
  gene_ids <- rownames(expr_mat)
  sample_ids_raw <- colnames(expr_mat)
} else {
  stop(
    paste(
      "Unable to determine expression matrix orientation.",
      "First 5 values in first column:",
      paste(utils::head(first_values, 5), collapse = ", "),
      "First 5 column names:",
      paste(utils::head(colnames(expr_df), 5), collapse = ", ")
    )
  )
}

sample_ids <- substr(sample_ids_raw, 1, 16)
colnames(expr_mat) <- sample_ids

gene_symbols <- parse_gene_symbol(gene_ids)
if (all(is.na(gene_symbols))) {
  stop(
    paste(
      "Unable to parse gene symbols from expression matrix gene IDs.",
      "Examples:",
      paste(utils::head(gene_ids, 10), collapse = ", ")
    )
  )
}

gene_annot <- data.frame(
  original_gene_id = gene_ids,
  gene_symbol = gene_symbols,
  stringsAsFactors = FALSE
)

keep_symbol <- !is.na(gene_annot$gene_symbol) & gene_annot$gene_symbol != ""
expr_mat <- expr_mat[keep_symbol, , drop = FALSE]
gene_annot <- gene_annot[keep_symbol, , drop = FALSE]
gene_annot$row_mean <- rowMeans(expr_mat, na.rm = TRUE)

if (is_raw_counts) {
  gene_annot <- gene_annot %>%
    arrange(.data$gene_symbol, desc(.data$row_mean))
  expr_mat <- expr_mat[match(gene_annot$original_gene_id, rownames(expr_mat)), , drop = FALSE]
  dedup_index <- !duplicated(gene_annot$gene_symbol)
  expr_dedup <- expr_mat[dedup_index, , drop = FALSE]
  gene_annot_dedup <- gene_annot[dedup_index, , drop = FALSE]
  rownames(expr_dedup) <- gene_annot_dedup$gene_symbol
} else {
  gene_annot$gene_symbol <- as.character(gene_annot$gene_symbol)
  expr_dt <- as.data.table(expr_mat)
  expr_dt[, gene_symbol := gene_annot$gene_symbol]
  expr_dedup_dt <- expr_dt[, lapply(.SD, mean, na.rm = TRUE), by = gene_symbol]
  expr_dedup <- as.matrix(expr_dedup_dt[, -1, drop = FALSE])
  rownames(expr_dedup) <- expr_dedup_dt$gene_symbol
  storage.mode(expr_dedup) <- "numeric"
  gene_annot_dedup <- gene_annot %>%
    group_by(.data$gene_symbol) %>%
    summarise(
      original_gene_id = paste(utils::head(.data$original_gene_id, 3), collapse = ";"),
      row_mean = mean(.data$row_mean, na.rm = TRUE),
      .groups = "drop"
    )
}

fwrite(as.data.table(gene_annot_dedup), file = output_gene_annot)

sample_map_all <- build_sample_mapping(colnames(expr_dedup))
duplicate_sample_map <- sample_map_all %>%
  group_by(.data$patient_id) %>%
  filter(n() > 1) %>%
  ungroup()
fwrite(as.data.table(duplicate_sample_map), file = output_duplicate_samples)
fwrite(as.data.table(sample_map_all), file = output_sample_map)

sample_map_kept <- sample_map_all %>%
  filter(.data$keep_sample) %>%
  select(sample_id, patient_id, sample_type_code, sample_type_rank)

expr_dedup <- expr_dedup[, sample_map_kept$sample_id, drop = FALSE]

if (is_raw_counts) {
  dge <- DGEList(counts = expr_dedup)
  keep_expr <- filterByExpr(dge, group = factor(rep("all", ncol(expr_dedup))))
  dge <- dge[keep_expr, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  expr_final <- cpm(dge, log = TRUE, prior.count = 1)
} else {
  expr_work <- expr_dedup
  if (max(expr_work, na.rm = TRUE) > 50 || stats::quantile(expr_work, 0.99, na.rm = TRUE) > 20) {
    expr_work <- log2(expr_work + 1)
  }
  keep_expr <- rowSums(expr_work > 0, na.rm = TRUE) >= max(2, ceiling(0.05 * ncol(expr_work)))
  expr_final <- expr_work[keep_expr, , drop = FALSE]
}

phenotype_df <- NULL
survival_df <- NULL
if (nrow(pheno_row) > 0 && file.exists(pheno_row$dest[1])) {
  phenotype_df <- read_table_auto(pheno_row$dest[1])
}
if (nrow(surv_row) > 0 && file.exists(surv_row$dest[1])) {
  survival_df <- read_table_auto(surv_row$dest[1])
}

if (is.null(phenotype_df) && is.null(survival_df)) {
  stop("Neither phenotype nor survival file is available for clinical preprocessing.")
}

prepare_clinical_table <- function(df, source_name) {
  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }

  cn <- colnames(df)
  sample_col_candidates <- c("sample", "SampleID", "sampleID", "submitter_id.samples", "sampleID")
  patient_col_candidates <- c("_PATIENT", "patient", "submitter_id", "case_submitter_id", "patient_id")

  sample_raw <- coalesce_columns(df, sample_col_candidates)
  patient_raw <- coalesce_columns(df, patient_col_candidates)

  sample_id <- ifelse(
    !is.na(sample_raw) & sample_raw != "",
    substr(toupper(as.character(sample_raw)), 1, 16),
    NA_character_
  )
  patient_id <- derive_patient_id(sample_id, patient_raw)

  os_time <- as_numeric_safely(coalesce_columns(df, c("OS.time", "OS_time", "_OS.time", "_OS_time")))
  os_status_raw <- coalesce_columns(df, c("OS", "_OS", "vital_status", "overall_survival_status", "os_status"))

  if (all(is.na(os_time))) {
    days_to_death <- as_numeric_safely(coalesce_columns(df, c("days_to_death", "death_days_to", "days_to_death.demographic")))
    days_to_last_follow_up <- as_numeric_safely(coalesce_columns(df, c("days_to_last_follow_up", "days_to_last_followup", "days_to_last_known_alive")))
    os_time <- ifelse(!is.na(days_to_death), days_to_death, days_to_last_follow_up)
  }

  os_status_chr <- toupper(as.character(os_status_raw))
  os_status <- ifelse(
    os_status_chr %in% c("DEAD", "DECEASED", "1"), 1L,
    ifelse(os_status_chr %in% c("ALIVE", "LIVING", "0"), 0L, NA_integer_)
  )

  age <- as_numeric_safely(coalesce_columns(df, c(
    "age_at_index", "age_at_diagnosis", "age_at_initial_pathologic_diagnosis",
    "age_at_initial_pathologic_diagnosis.demographic", "_AGE"
  )))
  if (all(is.na(age)) && "days_to_birth" %in% cn) {
    age <- abs(as_numeric_safely(df$days_to_birth)) / 365.25
  }
  if (!all(is.na(age)) && stats::median(age, na.rm = TRUE) > 120) {
    age <- age / 365.25
  }

  sex <- coalesce_columns(df, c("gender", "sex", "gender.demographic", "_SEX"))
  fab <- coalesce_columns(df, c("french_american_british_classification", "fab_classification", "morphology", "diagnosis"))
  wbc <- as_numeric_safely(coalesce_columns(df, c("white_blood_cell_count", "wbc", "wbc_at_diagnosis")))

  out <- data.frame(
    source = source_name,
    patient_id = patient_id,
    sample_id = sample_id,
    OS_time = os_time,
    OS_status = os_status,
    age = age,
    sex = sex,
    FAB = fab,
    WBC = wbc,
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$patient_id) & out$patient_id != "", , drop = FALSE]
  out
}

survival_clean <- prepare_clinical_table(survival_df, "survival")
phenotype_clean <- prepare_clinical_table(phenotype_df, "phenotype")

sample_map <- sample_map_kept %>%
  select(sample_id, patient_id)

clin_join <- sample_map
if (!is.null(survival_clean) && nrow(survival_clean) > 0) {
  survival_patient <- survival_clean %>%
    arrange(.data$patient_id) %>%
    distinct(.data$patient_id, .keep_all = TRUE) %>%
    select(-sample_id, -source)
  clin_join <- clin_join %>% left_join(survival_patient, by = "patient_id")
}

if (!is.null(phenotype_clean) && nrow(phenotype_clean) > 0) {
  phenotype_patient <- phenotype_clean %>%
    arrange(.data$patient_id) %>%
    distinct(.data$patient_id, .keep_all = TRUE) %>%
    select(-sample_id, -source)

  if (!all(c("OS_time", "OS_status", "age", "sex", "FAB", "WBC") %in% colnames(clin_join))) {
    clin_join <- clin_join %>% left_join(phenotype_patient, by = "patient_id")
  } else {
    for (nm in c("OS_time", "OS_status", "age", "sex", "FAB", "WBC")) {
      if (!nm %in% colnames(phenotype_patient)) {
        next
      }
      idx <- is.na(clin_join[[nm]]) | clin_join[[nm]] == ""
      clin_join[[nm]][idx] <- phenotype_patient[[nm]][match(clin_join$patient_id[idx], phenotype_patient$patient_id)]
    }
  }
}

clin_clean <- clin_join %>%
  mutate(sample_id = .data$sample_id) %>%
  select(patient_id, sample_id, OS_time, OS_status, age, sex, FAB, WBC)

matched_sample_ids <- intersect(colnames(expr_final), clin_clean$sample_id)
expr_matched <- expr_final[, matched_sample_ids, drop = FALSE]
clin_matched <- clin_clean %>%
  filter(.data$sample_id %in% matched_sample_ids) %>%
  arrange(match(.data$sample_id, matched_sample_ids))

if (!identical(colnames(expr_matched), clin_matched$sample_id)) {
  stop(
    paste(
      "Expression matrix columns and clin_clean$sample_id could not be aligned.",
      "Expression sample IDs:",
      paste(utils::head(colnames(expr_final), 5), collapse = ", "),
      "Clinical sample IDs:",
      paste(utils::head(clin_clean$sample_id, 5), collapse = ", "),
      "Clinical patient IDs:",
      paste(utils::head(clin_clean$patient_id, 5), collapse = ", ")
    )
  )
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
    "samples_with_complete_survival",
    "expression_data_type",
    "expression_source"
  ),
  value = c(
    length(gene_ids),
    sum(keep_symbol),
    nrow(expr_dedup),
    nrow(expr_final),
    length(sample_ids_raw),
    nrow(sample_map_kept),
    ncol(expr_matched),
    length(unique(clin_matched$patient_id)),
    samples_with_os_time,
    samples_with_os_status,
    samples_with_complete_survival,
    expression_data_type,
    expression_source
  ),
  stringsAsFactors = FALSE
)

saveRDS(expr_final, file = output_expr)
saveRDS(clin_matched, file = output_clin)
saveRDS(list(expr = expr_matched, clin = clin_matched), file = output_matched)
fwrite(as.data.table(sample_summary), file = output_summary)
save_session_info(output_session)

message("Preprocessing completed successfully.")
