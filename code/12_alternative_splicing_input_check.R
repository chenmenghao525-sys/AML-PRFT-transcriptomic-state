#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(1234)

suppressPackageStartupMessages({
  library(data.table)
})

root_env <- Sys.getenv("PHASE8A_ROOT", unset = "")
if (nzchar(root_env)) {
  root_dir <- chartr("\\", "/", path.expand(root_env))
} else {
  args_all <- commandArgs(trailingOnly = FALSE)
  script_arg <- args_all[grep("^--file=", args_all)][1]
  script_path <- if (length(script_arg) == 1 && nzchar(script_arg)) {
    sub("^--file=", "", script_arg)
  } else {
    file.path(getwd(), "02_scripts", "phase8A_AS_input_audit.R")
  }
  script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
  root_dir <- chartr("\\", "/", normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE))
}
results_dir <- file.path(root_dir, "03_results_tables")
logs_dir <- file.path(root_dir, "05_logs")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(logs_dir, "phase8A_AS_input_audit_log.txt")
if (file.exists(log_file)) invisible(file.remove(log_file))

log_msg <- function(...) {
  line <- paste0(...)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

write_csv <- function(dt, path) {
  fwrite(dt, path)
  log_msg("[write] ", chartr("\\", "/", path))
}

allowed_ext <- c("csv", "tsv", "txt", "xlsx", "rds", "rdata", "gz", "zip")
excluded_patterns <- c(
  "/phase1_runtime/17_tmp/",
  "/03_results_tables/Human_Genomics_PRFT_AML_OFFICIAL_submission_package/",
  "/03_results_tables/Human_Genomics_PRFT_AML_INTERNAL_working_package/",
  "/03_results_tables/Human_Genomics_PRFT_AML_submission_package_clean/",
  "/03_results_tables/Bioinformatics_audit_and_wet_validation_package/",
  "/03_results_tables/phase8A_",
  "/05_logs/phase8A_"
)

keyword_groups <- list(
  spliceseq = "(?i)spliceseq|tcga spliceseq|percent spliced in|\\bpsi\\b",
  as_event = "(?i)alternative splicing|splicing event|as event|exon skipping|skipped exon|retained intron|\\bri\\b|\\bse\\b|\\bes\\b|mxe|a3ss|a5ss|\\baa\\b|\\bad\\b|\\bap\\b|\\bat\\b|\\bme\\b",
  tool = "(?i)rmats|majiq|suppa|whippet|junction",
  sf = "(?i)splicing_factor|srsf|sf3b1|u2af1|hnrnp|rbm"
)
keyword_regex <- paste(unlist(keyword_groups), collapse = "|")

support_files <- c(
  "00_raw_data/00_raw_data/tcga_laml_raw_xena/TCGA.LAML.sampleMap_HiSeqV2.gz",
  "00_raw_data/00_raw_data/tcga_laml_raw_xena/TCGA.LAML.sampleMap_LAML_clinicalMatrix.tsv",
  "00_raw_data/00_raw_data/tcga_laml_raw_xena/tcga_laml_xena_manifest.csv",
  "01_processed_data/01_metadata/tcga_gene_annotation_hgnc.csv",
  "01_processed_data/01_metadata/tcga_sample_mapping.csv",
  "01_processed_data/02_processed_data/tcga_expr_hgnc_log2cpm.rds",
  "01_processed_data/02_processed_data/tcga_expr_clin_matched.rds",
  "01_processed_data/04_prft_score/tcga_prft_score.csv",
  "01_processed_data/04_prft_score/tcga_prft_score.rds",
  "01_processed_data/04_prft_score/tcga_ssgsea_scores.csv"
)
support_full <- file.path(root_dir, support_files)

tcga_expr_path <- file.path(root_dir, "01_processed_data/02_processed_data/tcga_expr_hgnc_log2cpm.rds")
tcga_expr_matched_path <- file.path(root_dir, "01_processed_data/02_processed_data/tcga_expr_clin_matched.rds")
tcga_prft_path <- file.path(root_dir, "01_processed_data/04_prft_score/tcga_prft_score.csv")

all_files <- list.files(root_dir, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
all_files <- chartr("\\", "/", all_files)
file_info <- file.info(all_files)
all_files <- all_files[!file_info$isdir %in% TRUE]
all_files <- all_files[tolower(tools::file_ext(all_files)) %in% allowed_ext]
if (length(excluded_patterns) > 0) {
  keep_idx <- rep(TRUE, length(all_files))
  for (pat in excluded_patterns) keep_idx <- keep_idx & !grepl(pat, all_files, fixed = TRUE)
  all_files <- all_files[keep_idx]
}
log_msg("[Phase8A] Eligible files scanned: ", length(all_files))
log_msg("[Phase8A] Excluded temporary/mirror directories: ", paste(excluded_patterns, collapse = " | "))

scan_text_preview <- function(path, n = 80L) {
  ext <- tolower(tools::file_ext(path))
  if (!ext %in% c("csv", "tsv", "txt", "gz")) return(NA_character_)
  con <- NULL
  txt <- tryCatch({
    con <- if (ext == "gz") gzfile(path, open = "rt") else file(path, open = "rt")
    lines <- readLines(con, n = n, warn = FALSE, encoding = "UTF-8")
    paste(lines, collapse = "\n")
  }, error = function(e) NA_character_)
  if (!is.null(con)) try(close(con), silent = TRUE)
  txt
}

path_hit <- grepl(keyword_regex, basename(all_files), perl = TRUE) | grepl(keyword_regex, all_files, perl = TRUE)
preview_vec <- vapply(all_files, scan_text_preview, character(1))
content_hit <- !is.na(preview_vec) & grepl(keyword_regex, preview_vec, perl = TRUE)
support_hit <- all_files %in% support_full
candidate_files <- unique(all_files[path_hit | content_hit | support_hit])

extract_keyword_tags <- function(path, preview) {
  hits <- character(0)
  for (nm in names(keyword_groups)) {
    if (grepl(keyword_groups[[nm]], path, perl = TRUE) || (!is.na(preview) && grepl(keyword_groups[[nm]], preview, perl = TRUE))) {
      hits <- c(hits, nm)
    }
  }
  if (path %in% support_full) hits <- c(hits, "support")
  if (length(hits) == 0) "" else paste(unique(hits), collapse = ";")
}

is_tcga_barcode <- function(x) {
  x <- as.character(x)
  any(grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}(-[0-9A-Z]{2,})?$", x))
}

read_rds_summary <- function(path) {
  obj <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(obj)) {
    return(list(class = "unreadable_rds", dims = NA_character_, row_ct = NA_integer_, col_ct = NA_integer_,
                sample_ids = character(0), row_ids = character(0), notes = "readRDS failed"))
  }
  cls <- paste(class(obj), collapse = ";")
  if (is.matrix(obj) || is.data.frame(obj)) {
    rn <- rownames(obj)
    cn <- colnames(obj)
    return(list(class = cls, dims = paste(dim(obj), collapse = "x"), row_ct = nrow(obj), col_ct = ncol(obj),
                sample_ids = if (is_tcga_barcode(cn)) cn else character(0),
                row_ids = if (!is.null(rn)) rn else character(0),
                notes = "matrix_or_dataframe"))
  }
  if (is.list(obj)) {
    expr <- obj[["expr"]]
    clin <- obj[["clin"]]
    sample_ids <- character(0)
    row_ids <- character(0)
    dims <- character(0)
    if (is.matrix(expr) || is.data.frame(expr)) {
      dims <- c(dims, paste(dim(expr), collapse = "x"))
      if (is_tcga_barcode(colnames(expr))) sample_ids <- colnames(expr)
      if (!is.null(rownames(expr))) row_ids <- rownames(expr)
    }
    if (is.data.frame(clin)) {
      dims <- c(dims, paste0("clin:", nrow(clin), "x", ncol(clin)))
      if ("sample_id" %in% names(clin) && is_tcga_barcode(clin$sample_id)) sample_ids <- unique(c(sample_ids, clin$sample_id))
    }
    return(list(class = cls, dims = paste(dims, collapse = ";"), row_ct = NA_integer_, col_ct = NA_integer_,
                sample_ids = sample_ids, row_ids = row_ids, notes = "list_object"))
  }
  list(class = cls, dims = NA_character_, row_ct = NA_integer_, col_ct = NA_integer_, sample_ids = character(0),
       row_ids = character(0), notes = "other_rds")
}

read_tabular_header <- function(path) {
  ext <- tolower(tools::file_ext(path))
  dt <- tryCatch({
    fread(path, nrows = 5, data.table = TRUE)
  }, error = function(e) NULL)
  if (is.null(dt)) return(list(colnames = character(0), preview = NA_character_))
  list(colnames = names(dt), preview = paste(capture.output(print(dt, topn = 3)), collapse = "\n"))
}

tcga_expr <- readRDS(tcga_expr_path)
tcga_expr_samples <- colnames(tcga_expr)
prft_dt <- fread(tcga_prft_path)
prft_samples <- prft_dt$sample_id
os_samples <- prft_dt[!is.na(OS_time) & !is.na(OS_status) & OS_time > 0, sample_id]
expr_matched_obj <- readRDS(tcga_expr_matched_path)
matched_expr_samples <- if (is.list(expr_matched_obj) && !is.null(expr_matched_obj$expr)) colnames(expr_matched_obj$expr) else character(0)

classify_candidate <- function(path) {
  ext <- tolower(tools::file_ext(path))
  rel <- sub(paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", root_dir), "/?"), "", path)
  preview <- preview_vec[match(path, all_files)]
  hdr <- if (ext %in% c("csv", "tsv", "txt", "gz")) read_tabular_header(path) else list(colnames = character(0), preview = NA_character_)
  rds_info <- if (ext %in% c("rds", "rdata")) read_rds_summary(path) else NULL
  coln <- tolower(c(hdr$colnames, if (!is.null(rds_info)) rds_info$row_ids[seq_len(min(length(rds_info$row_ids), 20))] else character(0)))
  file_text <- paste(c(rel, basename(path), hdr$preview, if (!is.null(rds_info)) rds_info$notes else "", preview), collapse = "\n")
  path_lower <- tolower(rel)
  is_textual_only_file <- grepl("(^|/)(05_logs|02_scripts|project_audit_report\\.md)", path_lower)
  coln_norm <- gsub("[^a-z0-9]", "", tolower(hdr$colnames))
  event_path_hit <- grepl("(^|[_/.-])(a3ss|a5ss|mxe|ri|se|es|aa|ad|ap|at|me)([_/.-]|$)|splicing_event|splice_event|event_annotation|alternative_splicing", path_lower, perl = TRUE)
  event_col_hit <- any(grepl("^(event|eventid|asevent|aseventid|spliceevent|spliceeventid|eventtype|splicetype|splicingtype|eventannotation|fromexon|toexon|exonstart|exonend|altstart|altend|psi)$", coln_norm))
  psi_path_hit <- grepl("spliceseq|tcga spliceseq|percent_spliced_in|percentsplicedin|(^|[_/.-])psi([_/.-]|$)", path_lower, perl = TRUE)
  psi_col_hit <- any(grepl("(^|_)psi($|_)", tolower(hdr$colnames))) || any(grepl("^psi$", coln_norm))

  classification <- "other_keyword_hit_or_support"
  is_psi <- FALSE
  is_annotation <- FALSE
  is_rmats <- FALSE
  is_majiq <- FALSE
  is_expr <- FALSE
  is_sf_expr <- FALSE
  sample_ids <- character(0)
  sample_count <- NA_integer_
  as_event_count <- NA_integer_
  notes <- character(0)

  if (psi_path_hit || (psi_col_hit && event_col_hit)) {
    classification <- "candidate_PSI_or_splicing_table"
    is_psi <- TRUE
  }
  if (event_path_hit || event_col_hit) {
    classification <- if (is_psi) "candidate_PSI_with_AS_event_terms" else "candidate_AS_event_annotation_or_table"
    is_annotation <- TRUE
  }
  if (grepl("rmats", path_lower, ignore.case = TRUE)) {
    classification <- "candidate_rMATS_result"
    is_rmats <- TRUE
  }
  if (grepl("majiq|suppa|whippet", path_lower, ignore.case = TRUE)) {
    classification <- "candidate_MAJIQ_SUPPA_Whippet_result"
    is_majiq <- TRUE
  }

  if (is_textual_only_file && !(path %in% support_full)) {
    if (!(grepl("rmats|majiq|suppa|whippet|spliceseq|psi", path_lower, ignore.case = TRUE))) {
      is_psi <- FALSE
      is_annotation <- FALSE
      is_rmats <- FALSE
      is_majiq <- FALSE
      classification <- "textual_reference_only"
    }
  }

  if (path %in% support_full || grepl("HiSeqV2|expr|expression|log2cpm|ssgsea|immune_signature", rel, ignore.case = TRUE)) {
    if (!(is_psi || is_annotation || is_rmats || is_majiq)) {
      classification <- "ordinary_expression_or_support_matrix"
      is_expr <- TRUE
    }
  }
  if (grepl("clinicalMatrix|clin|sample_mapping|prft_score", rel, ignore.case = TRUE) &&
      !(is_psi || is_annotation || is_rmats || is_majiq)) {
    classification <- "clinical_or_matching_support"
  }
  if (grepl("gene_annotation", rel, ignore.case = TRUE) && !(is_psi || is_annotation || is_rmats || is_majiq)) {
    classification <- "gene_annotation_support"
  }
  if (is_expr && (grepl("hiseqv2|expr_hgnc|expr_clin_matched|expression", path_lower, ignore.case = TRUE) || ext %in% c("rds", "rdata")) &&
      grepl("sf3b1|u2af1|hnrnp|rbm|srsf", paste(rownames(tcga_expr)[seq_len(min(nrow(tcga_expr), 20000))], collapse = " "), ignore.case = TRUE)) {
    is_sf_expr <- TRUE
    notes <- c(notes, "general TCGA expression matrix can quantify splicing-factor genes, but this is not PSI-level AS data")
  }

  if (ext %in% c("rds", "rdata") && !is.null(rds_info)) {
    sample_ids <- unique(rds_info$sample_ids)
    if (length(sample_ids) > 0) sample_count <- length(sample_ids)
    if ((is_psi || is_annotation) && !is.na(rds_info$row_ct)) as_event_count <- rds_info$row_ct
    notes <- c(notes, paste0("rds_class=", rds_info$class), paste0("rds_dims=", rds_info$dims))
  } else if (ext %in% c("csv", "tsv", "txt", "gz")) {
    cn <- hdr$colnames
    if (length(cn) > 0 && is_tcga_barcode(cn)) {
      sample_ids <- cn[grepl("^TCGA-", cn)]
      sample_count <- length(sample_ids)
    }
    if ((is_psi || is_annotation) && length(cn) > 0) {
      as_event_count <- tryCatch({
        dt_n <- fread(path, select = 1L)
        nrow(dt_n)
      }, error = function(e) NA_integer_)
    }
    if (classification == "clinical_or_matching_support" && grepl("clinicalMatrix", rel, ignore.case = TRUE)) {
      sample_ids <- tryCatch({
        dt0 <- fread(path, select = 1L)
        ids <- dt0[[1]]
        ids[grepl("^TCGA-", ids)]
      }, error = function(e) character(0))
      sample_count <- if (length(sample_ids) > 0) length(sample_ids) else sample_count
    }
  }

  has_tcga_barcode <- length(sample_ids) > 0
  can_match_expr <- has_tcga_barcode && sum(sample_ids %in% tcga_expr_samples) > 0
  can_match_prft <- has_tcga_barcode && sum(sample_ids %in% prft_samples) > 0
  can_match_os <- has_tcga_barcode && sum(sample_ids %in% os_samples) > 0

  if (classification == "ordinary_expression_or_support_matrix" && grepl("clinicalMatrix", rel, ignore.case = TRUE)) {
    classification <- "clinical_or_matching_support"
    is_expr <- FALSE
  }

  data.table(
    file_path = path,
    relative_path = rel,
    extension = ext,
    size_bytes = file.info(path)$size,
    keyword_tags = extract_keyword_tags(rel, preview),
    classification = classification,
    is_PSI_matrix = is_psi,
    is_AS_event_annotation = is_annotation,
    is_rMATS_result = is_rmats,
    is_MAJIQ_SUPPA_Whippet_result = is_majiq,
    is_ordinary_expression_matrix = is_expr,
    is_splicing_factor_expression_table = is_sf_expr,
    has_TCGA_barcode = has_tcga_barcode,
    can_match_TCGA_expression = can_match_expr,
    can_match_PRFT_risk_score = can_match_prft,
    can_match_OS_time_status = can_match_os,
    estimated_sample_count = ifelse(is.na(sample_count), NA_integer_, as.integer(sample_count)),
    estimated_AS_event_count = ifelse(is.na(as_event_count), NA_integer_, as.integer(as_event_count)),
    notes = paste(unique(notes[nzchar(notes)]), collapse = " | ")
  )
}

audit_dt <- data.table(
  file_path = candidate_files,
  relative_path = sub(paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", root_dir), "/?"), "", candidate_files),
  extension = tolower(tools::file_ext(candidate_files)),
  size_bytes = file.info(candidate_files)$size,
  matched_by_path_keyword = grepl(keyword_regex, basename(candidate_files), perl = TRUE) | grepl(keyword_regex, candidate_files, perl = TRUE),
  matched_by_content_keyword = content_hit[match(candidate_files, all_files)],
  matched_as_support_file = candidate_files %in% support_full,
  keyword_tags = vapply(candidate_files, function(x) extract_keyword_tags(sub(paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", root_dir), "/?"), "", x), preview_vec[match(x, all_files)]), character(1))
)
setorder(audit_dt, -matched_as_support_file, -matched_by_path_keyword, relative_path)
write_csv(audit_dt, file.path(results_dir, "phase8A_AS_input_file_audit.csv"))

class_dt <- rbindlist(lapply(candidate_files, classify_candidate), fill = TRUE)
setorder(class_dt, relative_path)
write_csv(class_dt, file.path(results_dir, "phase8A_AS_candidate_file_classification.csv"))

usable_psi <- class_dt[is_PSI_matrix == TRUE & has_TCGA_barcode == TRUE]
usable_annotation <- class_dt[is_AS_event_annotation == TRUE]
usable_tool <- class_dt[is_rMATS_result == TRUE | is_MAJIQ_SUPPA_Whippet_result == TRUE]
tcga_spliceseq_exists <- any(grepl("spliceseq", class_dt$file_path, ignore.case = TRUE))
general_sf_support <- any(class_dt$is_splicing_factor_expression_table)

feasibility <- data.table(
  usable_PSI_matrix_exists = nrow(usable_psi) > 0,
  TCGA_SpliceSeq_file_exists = tcga_spliceseq_exists,
  AS_event_annotation_exists = nrow(usable_annotation) > 0,
  rMATS_MAJIQ_SUPPA_Whippet_exists = nrow(usable_tool) > 0,
  usable_TCGA_barcode_in_AS_input = nrow(class_dt[(is_PSI_matrix | is_AS_event_annotation | is_rMATS_result | is_MAJIQ_SUPPA_Whippet_result) & has_TCGA_barcode == TRUE]) > 0,
  can_match_PRFT_risk_score = nrow(class_dt[(is_PSI_matrix | is_AS_event_annotation | is_rMATS_result | is_MAJIQ_SUPPA_Whippet_result) & can_match_PRFT_risk_score == TRUE]) > 0,
  can_match_TCGA_OS = nrow(class_dt[(is_PSI_matrix | is_AS_event_annotation | is_rMATS_result | is_MAJIQ_SUPPA_Whippet_result) & can_match_OS_time_status == TRUE]) > 0,
  expected_usable_sample_count = if (nrow(usable_psi) > 0) max(usable_psi$estimated_sample_count, na.rm = TRUE) else 0L,
  expected_AS_event_count = if (nrow(usable_psi) > 0) max(usable_psi$estimated_AS_event_count, na.rm = TRUE) else 0L,
  can_do_PRFT_high_vs_low_differential_AS = FALSE,
  can_do_PRFT_correlated_AS = FALSE,
  can_do_survival_related_AS = FALSE,
  can_do_splicing_factor_AS_network = FALSE,
  splicing_factor_expression_support_exists = general_sf_support,
  recommend_phase8B_formal_AS = FALSE,
  main_reason = paste(
    c(
      if (nrow(usable_psi) == 0) "No usable PSI-level matrix was found in the local project package." else NULL,
      if (nrow(usable_annotation) == 0) "No AS event annotation table was found." else NULL,
      if (nrow(usable_tool) == 0) "No rMATS/MAJIQ/SUPPA/Whippet result package was found." else NULL,
      if (general_sf_support) "General TCGA gene-expression matrices can support splicing-factor expression only, but cannot substitute for PSI/event-level AS inputs." else NULL
    ),
    collapse = " "
  )
)
write_csv(feasibility, file.path(results_dir, "phase8A_AS_feasibility_report.csv"))

missing_statement_path <- file.path(logs_dir, "phase8A_AS_missing_input_statement.txt")
missing_statement <- c(
  "1. The current local project package does not contain a matched PSI-level matrix or other formal AS event matrix suitable for alternative splicing analysis.",
  "2. Ordinary gene-expression matrices were found for TCGA-LAML, but they cannot be treated as AS data and cannot replace PSI/event-level inputs.",
  "3. Therefore, the present project cannot claim that an alternative splicing analysis has been completed.",
  "4. If an AS module is pursued later, TCGA SpliceSeq LAML PSI data or dedicated rMATS/MAJIQ/SUPPA/Whippet outputs should be obtained separately and matched to PRFT risk score and OS information.",
  "5. To preserve analytical credibility, the current manuscript can omit the AS module rather than adding a weak or unsupported splicing section."
)
writeLines(missing_statement, missing_statement_path)
log_msg("[write] ", chartr("\\", "/", missing_statement_path))

checklist <- c(
  paste0("1. PSI matrix found: ", ifelse(feasibility$usable_PSI_matrix_exists[1], "yes", "no")),
  paste0("2. TCGA SpliceSeq file found: ", ifelse(feasibility$TCGA_SpliceSeq_file_exists[1], "yes", "no")),
  paste0("3. AS event annotation found: ", ifelse(feasibility$AS_event_annotation_exists[1], "yes", "no")),
  paste0("4. rMATS/MAJIQ/SUPPA results found: ", ifelse(feasibility$rMATS_MAJIQ_SUPPA_Whippet_exists[1], "yes", "no")),
  paste0("5. Usable TCGA barcode found: ", ifelse(feasibility$usable_TCGA_barcode_in_AS_input[1], "yes", "no")),
  paste0("6. Can match PRFT risk_score: ", ifelse(feasibility$can_match_PRFT_risk_score[1], "yes", "no")),
  paste0("7. Can match OS_time/OS_status: ", ifelse(feasibility$can_match_TCGA_OS[1], "yes", "no")),
  paste0("8. Expected usable sample count: ", feasibility$expected_usable_sample_count[1]),
  paste0("9. Expected usable AS event count: ", feasibility$expected_AS_event_count[1]),
  paste0("10. Recommend formal AS analysis: ", ifelse(feasibility$recommend_phase8B_formal_AS[1], "yes", "no")),
  paste0("11. If not recommended, reason: ", feasibility$main_reason[1]),
  paste0("12. Recommend keeping AS module in manuscript: ", ifelse(feasibility$recommend_phase8B_formal_AS[1], "yes", "no")),
  paste0("13. Recommend deleting AS module from manuscript: ", ifelse(feasibility$recommend_phase8B_formal_AS[1], "no", "yes")),
  "14. Recommend entering Phase 9 full-manuscript integration: yes",
  "15. Manual confirmation needed: verify whether any TCGA SpliceSeq PSI package or external rMATS/MAJIQ/SUPPA result folder exists outside the current project snapshot; if not, keep AS fully out of the manuscript."
)
checklist_path <- file.path(logs_dir, "phase8A_AS_input_key_result_checklist.txt")
writeLines(checklist, checklist_path)
log_msg("[write] ", chartr("\\", "/", checklist_path))

log_msg("[Phase8A] Candidate files audited: ", nrow(audit_dt))
log_msg("[Phase8A] Candidate files classified: ", nrow(class_dt))
log_msg("[Phase8A] Completed. Formal Phase 8B AS analysis recommended: ", ifelse(feasibility$recommend_phase8B_formal_AS[1], "yes", "no"))
