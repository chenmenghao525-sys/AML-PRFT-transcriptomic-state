#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
options(timeout = 1200)
options(download.file.method = "libcurl")

suppressPackageStartupMessages({
  library(data.table)
})

dir.create("00_raw_data/geo_validation", recursive = TRUE, showWarnings = FALSE)
dir.create("08_validation", recursive = TRUE, showWarnings = FALSE)
dir.create("16_logs", recursive = TRUE, showWarnings = FALSE)

save_session_info <- function(path) {
  writeLines(capture.output(sessionInfo()), con = path)
}

signature_genes <- c(
  "ANKRD5", "SAT1", "CLCN5", "ITGB2", "CPNE8",
  "C7orf50", "METTL7B", "ARHGEF35", "ACOX2"
)

get_bucket <- function(acc) {
  prefix <- gsub("[0-9]", "", acc)
  digits <- gsub("[^0-9]", "", acc)
  bucket_digits <- substr(digits, 1, max(1, nchar(digits) - 3))
  paste0(prefix, bucket_digits, "nnn")
}

get_gpl_bucket <- function(gpl) {
  prefix <- gsub("[0-9]", "", gpl)
  digits <- gsub("[^0-9]", "", gpl)
  if (nchar(digits) <= 3) {
    return(paste0(prefix, "nnn"))
  }
  bucket_digits <- substr(digits, 1, nchar(digits) - 3)
  paste0(prefix, bucket_digits, "nnn")
}

download_to_file <- function(url, destfile) {
  out <- tryCatch({
    if (file.exists(destfile) && file.info(destfile)$size > 0) {
      return(TRUE)
    }
    if (file.exists(destfile) && file.info(destfile)$size == 0) {
      unlink(destfile)
    }
    utils::download.file(url, destfile = destfile, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) {
    structure(FALSE, error_message = conditionMessage(e))
  })
  out
}

read_url_text <- function(url) {
  tryCatch(
    paste(readLines(url, warn = FALSE), collapse = "\n"),
    error = function(e) NA_character_
  )
}

extract_matrix_files <- function(listing_text) {
  if (is.na(listing_text) || !nzchar(listing_text)) {
    return(character(0))
  }
  hits <- gregexpr("GSE[0-9A-Za-z_-]+(_|-)?[A-Za-z0-9_-]*series_matrix[^\"'<> ]*txt\\.gz", listing_text, perl = TRUE)
  files <- regmatches(listing_text, hits)[[1]]
  unique(files)
}

parse_series_matrix <- function(localfile) {
  lines <- readLines(gzfile(localfile), warn = FALSE)

  begin_idx <- grep("^!series_matrix_table_begin", lines)
  end_idx <- grep("^!series_matrix_table_end", lines)
  if (length(begin_idx) == 0) {
    stop("series_matrix table not found.")
  }
  if (length(end_idx) == 0 || end_idx[1] <= begin_idx[1]) {
    table_lines <- lines[(begin_idx[1] + 1):length(lines)]
  } else {
    table_lines <- lines[(begin_idx[1] + 1):(end_idx[1] - 1)]
  }

  expr_df <- read.delim(
    text = paste(table_lines, collapse = "\n"),
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    quote = "\"",
    comment.char = ""
  )
  if (ncol(expr_df) < 2) {
    stop("Expression matrix has fewer than 2 columns.")
  }

  id_col <- colnames(expr_df)[1]
  expr_ids <- as.character(expr_df[[1]])
  expr_mat <- as.matrix(expr_df[, -1, drop = FALSE])
  storage.mode(expr_mat) <- "numeric"
  rownames(expr_mat) <- expr_ids

  sample_lines <- grep("^!Sample_", lines, value = TRUE)
  sample_meta <- data.frame(stringsAsFactors = FALSE)
  if (length(sample_lines) > 0) {
    parsed <- lapply(sample_lines, function(x) strsplit(x, "\t", fixed = TRUE)[[1]])
    max_len <- max(vapply(parsed, length, integer(1)))
    parsed <- lapply(parsed, function(x) c(x, rep("", max_len - length(x))))
    parsed_mat <- do.call(rbind, parsed)
    field_names <- sub("^!Sample_", "", parsed_mat[, 1])
    field_names_unique <- make.unique(field_names, sep = "_")
    sample_meta <- as.data.frame(t(parsed_mat[, -1, drop = FALSE]), stringsAsFactors = FALSE)
    colnames(sample_meta) <- field_names_unique
    sample_meta[] <- lapply(sample_meta, function(v) gsub('^"|"$', "", v))
    if (nrow(sample_meta) == ncol(expr_mat)) {
      sample_meta$sample_id <- colnames(expr_mat)
    }
  }

  series_platform <- grep("^!Series_platform_id", lines, value = TRUE)
  sample_platform <- if (nrow(sample_meta) > 0 && "platform_id" %in% colnames(sample_meta)) unique(sample_meta$platform_id) else character(0)
  platform_ids <- unique(c(
    gsub('^!Series_platform_id\\s+|^!Series_platform_id\\t', "", series_platform),
    sample_platform
  ))
  platform_ids <- gsub('^"|"$', "", platform_ids)
  platform_ids <- platform_ids[nzchar(platform_ids)]

  list(
    expr = expr_mat,
    sample_meta = sample_meta,
    platform_ids = platform_ids,
    id_col = id_col
  )
}

download_gpl_annotation <- function(gpl) {
  bucket <- get_gpl_bucket(gpl)
  url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/platforms/%s/%s/soft/%s_family.soft.gz", bucket, gpl, gpl)
  dest <- file.path("00_raw_data/geo_validation", paste0(gpl, "_family.soft.gz"))
  ok <- download_to_file(url, dest)
  list(ok = ok, url = url, dest = dest)
}

parse_gpl_annotation <- function(localfile) {
  lines <- readLines(gzfile(localfile), warn = FALSE)
  begin_idx <- grep("^!platform_table_begin", lines)
  end_idx <- grep("^!platform_table_end", lines)
  if (length(begin_idx) == 0) {
    stop("GPL table not found.")
  }
  if (length(end_idx) == 0 || end_idx[1] <= begin_idx[1]) {
    table_lines <- lines[(begin_idx[1] + 1):length(lines)]
  } else {
    table_lines <- lines[(begin_idx[1] + 1):(end_idx[1] - 1)]
  }
  gpl_df <- read.delim(
    text = paste(table_lines, collapse = "\n"),
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    quote = "\"",
    comment.char = ""
  )
  gpl_df
}

clean_symbol_string <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("///", ";", x, fixed = TRUE)
  x <- gsub("//", ";", x, fixed = TRUE)
  x <- gsub("\\|", ";", x)
  x <- gsub(",", ";", x, fixed = TRUE)
  x <- gsub(";", ";", x, fixed = TRUE)
  x
}

map_expr_to_symbol <- function(expr_mat, gpl_df = NULL) {
  row_ids <- rownames(expr_mat)
  direct_hit_count <- sum(signature_genes %in% row_ids)

  if (direct_hit_count >= 5) {
    expr_dt <- data.table(probe_id = row_ids, expr_mat, keep.rownames = FALSE)
    expr_dt[, gene_symbol := probe_id]
    return(list(
      success = TRUE,
      expr_symbol = expr_dt,
      symbol_source = "expression_rowname_direct",
      symbol_column = "rowname"
    ))
  }

  if (is.null(gpl_df)) {
    return(list(success = FALSE, reason = "No GPL annotation available and rownames do not appear to be gene symbols."))
  }

  norm_cols <- gsub("[^a-z0-9]", "", tolower(colnames(gpl_df)))
  symbol_candidates <- c("genesymbol", "symbol", "hgncsymbol", "gene_symbol")
  id_candidates <- c("id", "idref", "probeid")

  symbol_idx <- match(symbol_candidates, norm_cols)
  symbol_idx <- symbol_idx[!is.na(symbol_idx)][1]
  id_idx <- match(id_candidates, norm_cols)
  id_idx <- id_idx[!is.na(id_idx)][1]

  if (is.na(symbol_idx) || is.na(id_idx)) {
    return(list(success = FALSE, reason = "GPL annotation lacks recognizable probe ID or gene symbol column."))
  }

  gpl_map <- data.table(
    probe_id = as.character(gpl_df[[id_idx]]),
    raw_symbol = clean_symbol_string(gpl_df[[symbol_idx]])
  )
  gpl_map <- gpl_map[nzchar(probe_id) & nzchar(raw_symbol)]

  probe_iqr <- apply(expr_mat, 1, IQR, na.rm = TRUE)
  probe_stats <- data.table(probe_id = row_ids, probe_iqr = probe_iqr)
  expanded_map <- gpl_map[, .(gene_symbol = unlist(strsplit(raw_symbol, ";", fixed = TRUE))), by = "probe_id"]
  expanded_map[, gene_symbol := trimws(gene_symbol)]
  expanded_map <- unique(expanded_map[nzchar(gene_symbol) & gene_symbol != "---"])

  merged_map <- merge(expanded_map, probe_stats, by = "probe_id", all.x = FALSE, all.y = FALSE)
  if (nrow(merged_map) == 0) {
    return(list(success = FALSE, reason = "No overlap between expression probe IDs and GPL annotation."))
  }
  setorder(merged_map, gene_symbol, -probe_iqr)
  dedup_map <- merged_map[, .SD[1], by = gene_symbol]
  expr_sub <- expr_mat[dedup_map$probe_id, , drop = FALSE]
  expr_symbol <- data.table(gene_symbol = dedup_map$gene_symbol, expr_sub, keep.rownames = FALSE)

  list(
    success = TRUE,
    expr_symbol = expr_symbol,
    symbol_source = "GPL_annotation",
    symbol_column = colnames(gpl_df)[symbol_idx]
  )
}

detect_log2 <- function(expr_mat) {
  vals <- as.numeric(expr_mat)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) return(NA)
  q99 <- stats::quantile(vals, 0.99, na.rm = TRUE)
  vmax <- max(vals, na.rm = TRUE)
  if (vmax > 100 || q99 > 50) FALSE else TRUE
}

detect_clinical_candidates <- function(sample_meta) {
  if (nrow(sample_meta) == 0) return(data.table())
  out <- list()
  for (nm in colnames(sample_meta)) {
    vals <- as.character(sample_meta[[nm]])
    vals <- vals[!is.na(vals) & nzchar(vals)]
    if (length(vals) == 0) next
    nm_low <- tolower(nm)
    example <- paste(utils::head(unique(vals), 3), collapse = " | ")
    value_blob <- tolower(paste(utils::head(unique(vals), 20), collapse = " | "))
    if (grepl("overall survival|survival|os|event|status|dead|death|deceased|vital|follow|efs|rfs|relapse", nm_low) ||
        grepl("overall survival|survival|os|event|status|dead|death|deceased|vital|follow|efs|rfs|relapse", value_blob)) {
      out[[length(out) + 1]] <- data.table(
        candidate_column = nm,
        example_values = example
      )
    }
  }
  if (length(out) == 0) data.table() else rbindlist(out, fill = TRUE)
}

parse_numeric_field <- function(x) {
  x <- as.character(x)
  extract_first_group <- function(pattern, values) {
    m <- regexec(pattern, values, perl = TRUE)
    hits <- regmatches(values, m)
    out <- vapply(hits, function(h) {
      if (length(h) >= 2) h[2] else NA_character_
    }, character(1))
    suppressWarnings(as.numeric(out))
  }

  os_eq <- extract_first_group("(?i)\\bOS\\s*=\\s*([0-9.]+)", x)
  if (sum(!is.na(os_eq)) > 0) {
    return(os_eq)
  }

  overall_survival <- extract_first_group("(?i)overall survival[^0-9]*([0-9.]+)", x)
  if (sum(!is.na(overall_survival)) > 0) {
    return(overall_survival)
  }

  x <- gsub("^.*?:\\s*", "", x)
  x <- gsub("[^0-9.\\-]", "", x)
  suppressWarnings(as.numeric(x))
}

extract_os_info <- function(sample_meta) {
  if (nrow(sample_meta) == 0) {
    return(list(
      has_os_time = FALSE, has_os_status = FALSE, os_time = NULL, os_status = NULL,
      survival_samples = 0, event_count = 0, notes = "No sample metadata available."
    ))
  }

  cn <- colnames(sample_meta)
  cn_low <- tolower(cn)

  time_idx_name <- which(
    grepl("overall.survival|overall survival|os.time|survival.time|days_to_death|days_to_last_follow_up", cn_low) &
      !grepl("efs|rfs|relapse", cn_low)
  )
  time_idx_value <- which(vapply(sample_meta, function(v) {
    vals <- tolower(as.character(v))
    any(grepl("overall survival|survival time|os.time|os time|\\bos\\b|days to death|last follow", vals)) &&
      any(grepl("[0-9]", vals))
  }, logical(1)))
  time_idx <- unique(c(time_idx_name, time_idx_value))

  status_idx_name <- which(
    grepl("os.status|overall.survival.status|vital.status|life.status", cn_low) &
      !grepl("efs|rfs|relapse", cn_low)
  )
  status_idx_value <- which(vapply(sample_meta, function(v) {
    vals_raw <- tolower(as.character(v))
    vals <- trimws(gsub("^.*?:\\s*", "", vals_raw))
    mapped <- ifelse(vals %in% c("1", "dead", "deceased", "event"), 1,
                     ifelse(vals %in% c("0", "alive", "censored"), 0, NA))
    has_context <- any(grepl("life status|vital status|overall survival status|os status|alive|dead|deceased|censored", vals_raw))
    has_context && sum(!is.na(mapped)) >= max(10, floor(0.3 * length(mapped)))
  }, logical(1)))
  status_idx <- unique(c(status_idx_name, status_idx_value))

  choose_best_time <- function(idxs) {
    for (idx in idxs) {
      vals <- parse_numeric_field(sample_meta[[idx]])
      if (sum(!is.na(vals)) >= max(10, floor(0.3 * length(vals)))) return(idx)
    }
    NA_integer_
  }
  choose_best_status <- function(idxs) {
    for (idx in idxs) {
      vals <- tolower(as.character(sample_meta[[idx]]))
      vals2 <- trimws(gsub("^.*?:\\s*", "", vals))
      mapped <- ifelse(vals2 %in% c("1", "dead", "deceased", "event"), 1,
                       ifelse(vals2 %in% c("0", "alive", "censored"), 0, NA))
      if (sum(!is.na(mapped)) >= max(10, floor(0.3 * length(mapped)))) return(idx)
    }
    NA_integer_
  }

  time_col_idx <- choose_best_time(time_idx)
  status_col_idx <- choose_best_status(status_idx)

  os_time <- NULL
  os_status <- NULL
  notes <- character(0)

  if (!is.na(time_col_idx)) {
    os_time <- parse_numeric_field(sample_meta[[time_col_idx]])
    col_name <- cn[time_col_idx]
    col_blob <- tolower(paste(utils::head(unique(as.character(sample_meta[[time_col_idx]])), 10), collapse = " "))
    if (grepl("month|months|mo", tolower(col_name)) || grepl("month|months| mo", col_blob)) {
      os_time <- os_time * 30.44
      notes <- c(notes, "OS_time converted from months to days.")
    } else if (!(grepl("day|days", tolower(col_name)) || grepl("day|days", col_blob))) {
      notes <- c(notes, "OS_time unit uncertain.")
    }
  }

  if (!is.na(status_col_idx)) {
    vals <- tolower(as.character(sample_meta[[status_col_idx]]))
    vals <- trimws(gsub("^.*?:\\s*", "", vals))
    os_status <- ifelse(vals %in% c("1", "dead", "deceased", "event"), 1,
                        ifelse(vals %in% c("0", "alive", "censored"), 0, NA))
  }

  complete_n <- if (!is.null(os_time) && !is.null(os_status)) sum(!is.na(os_time) & !is.na(os_status)) else 0
  event_n <- if (!is.null(os_status)) sum(os_status == 1, na.rm = TRUE) else 0

  list(
    has_os_time = !is.null(os_time),
    has_os_status = !is.null(os_status),
    os_time = os_time,
    os_status = os_status,
    survival_samples = complete_n,
    event_count = event_n,
    notes = paste(unique(notes), collapse = " ")
  )
}

gse_list <- c("GSE37642", "GSE6891", "GSE12417", "GSE14468")

manifest_list <- list()
feas_list <- list()
coverage_list <- list()
clinical_candidate_list <- list()

write_current_outputs <- function() {
  manifest_dt <- if (length(manifest_list) > 0) rbindlist(manifest_list, fill = TRUE) else data.table()
  feas_dt <- if (length(feas_list) > 0) rbindlist(feas_list, fill = TRUE) else data.table()
  coverage_dt <- if (length(coverage_list) > 0) rbindlist(coverage_list, fill = TRUE) else data.table()
  clinical_dt <- if (length(clinical_candidate_list) > 0) rbindlist(clinical_candidate_list, fill = TRUE) else data.table()
  usable_dt <- feas_dt[usable_for_validation == TRUE]

  fwrite(manifest_dt, "08_validation/geo_download_manifest.csv")
  fwrite(feas_dt, "08_validation/geo_feasibility_summary.csv")
  fwrite(coverage_dt, "08_validation/geo_signature_gene_coverage.csv")
  fwrite(clinical_dt, "08_validation/geo_clinical_field_candidates.csv")
  fwrite(usable_dt, "08_validation/geo_usable_datasets.csv")
}

for (gse in gse_list) {
  bucket <- get_bucket(gse)
  matrix_url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/", bucket, gse)
  listing_text <- read_url_text(matrix_url)
  matrix_files <- extract_matrix_files(listing_text)
  if (length(matrix_files) == 0) {
    local_cached <- list.files(
      "00_raw_data/geo_validation",
      pattern = paste0("^", gse, ".*series_matrix.*txt\\.gz$"),
      full.names = FALSE
    )
    matrix_files <- unique(local_cached)
  }

  if (length(matrix_files) == 0) {
    manifest_list[[length(manifest_list) + 1]] <- data.table(
      gse = gse,
      file_type = "series_matrix",
      file_name = NA_character_,
      url = matrix_url,
      download_status = "failed_listing_or_no_matrix",
      local_path = NA_character_,
      file_size = NA_real_,
      notes = "Could not list matrix directory or no series_matrix file found."
    )
    feas_list[[length(feas_list) + 1]] <- data.table(
      gse = gse,
      matrix_file = NA_character_,
      platform = NA_character_,
      series_matrix_downloaded = FALSE,
      expression_parsed = FALSE,
      sample_count = NA_integer_,
      log2_likely = NA,
      gene_symbol_available = FALSE,
      gene_symbol_source = NA_character_,
      gene_symbol_column = NA_character_,
      signature_genes_covered = 0,
      missing_signature_genes = paste(signature_genes, collapse = ";"),
      has_os_time = FALSE,
      has_os_status = FALSE,
      survival_samples = 0,
      event_count = 0,
      usable_for_validation = FALSE,
        signature_availability = "not_available",
        notes = "Series matrix unavailable."
      )
      write_current_outputs()
      next
    }

  for (mf in matrix_files) {
    matrix_file_url <- paste0(matrix_url, mf)
    matrix_dest <- file.path("00_raw_data/geo_validation", mf)
    ok <- download_to_file(matrix_file_url, matrix_dest)
    manifest_list[[length(manifest_list) + 1]] <- data.table(
      gse = gse,
      file_type = "series_matrix",
      file_name = mf,
      url = matrix_file_url,
      download_status = if (isTRUE(ok)) "downloaded" else "failed",
      local_path = if (isTRUE(ok)) normalizePath(matrix_dest, winslash = "/", mustWork = FALSE) else NA_character_,
      file_size = if (isTRUE(ok) && file.exists(matrix_dest)) file.info(matrix_dest)$size else NA_real_,
      notes = if (isTRUE(ok)) "" else attr(ok, "error_message")
    )

    if (!isTRUE(ok)) {
      write_current_outputs()
      next
    }

    parsed <- tryCatch(parse_series_matrix(matrix_dest), error = function(e) e)
    if (inherits(parsed, "error")) {
      feas_list[[length(feas_list) + 1]] <- data.table(
        gse = gse,
        matrix_file = mf,
        platform = NA_character_,
        series_matrix_downloaded = TRUE,
        expression_parsed = FALSE,
        sample_count = NA_integer_,
        log2_likely = NA,
        gene_symbol_available = FALSE,
        gene_symbol_source = NA_character_,
        gene_symbol_column = NA_character_,
        signature_genes_covered = 0,
        missing_signature_genes = paste(signature_genes, collapse = ";"),
        has_os_time = FALSE,
        has_os_status = FALSE,
        survival_samples = 0,
        event_count = 0,
        usable_for_validation = FALSE,
        signature_availability = "not_available",
        notes = conditionMessage(parsed)
      )
      write_current_outputs()
      next
    }

    sample_meta <- parsed$sample_meta
    platforms <- if (length(parsed$platform_ids) > 0) paste(parsed$platform_ids, collapse = ";") else NA_character_

    clinical_candidates <- detect_clinical_candidates(sample_meta)
    if (nrow(clinical_candidates) > 0) {
      clinical_candidates[, `:=`(gse = gse, matrix_file = mf, platform = platforms)]
      clinical_candidate_list[[length(clinical_candidate_list) + 1]] <- clinical_candidates[, .(gse, matrix_file, platform, candidate_column, example_values)]
    }

    log2_flag <- detect_log2(parsed$expr)

    gpl_df <- NULL
    file_platform_match <- regmatches(mf, regexpr("GPL[0-9]+", mf))
    gpl_to_use <- if (length(file_platform_match) > 0 && nzchar(file_platform_match[1])) {
      file_platform_match[1]
    } else if (length(parsed$platform_ids) > 0) {
      parsed$platform_ids[1]
    } else {
      NA_character_
    }

    if (!any(signature_genes %in% rownames(parsed$expr)) && !is.na(gpl_to_use)) {
      gpl_try <- download_gpl_annotation(gpl_to_use)
      manifest_list[[length(manifest_list) + 1]] <- data.table(
        gse = gse,
        file_type = "GPL_annotation",
        file_name = basename(gpl_try$dest),
        url = gpl_try$url,
        download_status = if (isTRUE(gpl_try$ok)) "downloaded" else "failed",
        local_path = if (isTRUE(gpl_try$ok)) normalizePath(gpl_try$dest, winslash = "/", mustWork = FALSE) else NA_character_,
        file_size = if (isTRUE(gpl_try$ok) && file.exists(gpl_try$dest)) file.info(gpl_try$dest)$size else NA_real_,
        notes = if (isTRUE(gpl_try$ok)) "" else attr(gpl_try$ok, "error_message")
      )
      if (isTRUE(gpl_try$ok)) {
        gpl_df <- tryCatch(parse_gpl_annotation(gpl_try$dest), error = function(e) NULL)
      }
    }

    mapped <- map_expr_to_symbol(parsed$expr, gpl_df)
    if (!isTRUE(mapped$success)) {
      feas_list[[length(feas_list) + 1]] <- data.table(
        gse = gse,
        matrix_file = mf,
        platform = platforms,
        series_matrix_downloaded = TRUE,
        expression_parsed = TRUE,
        sample_count = ncol(parsed$expr),
        log2_likely = log2_flag,
        gene_symbol_available = FALSE,
        gene_symbol_source = NA_character_,
        gene_symbol_column = NA_character_,
        signature_genes_covered = 0,
        missing_signature_genes = paste(signature_genes, collapse = ";"),
        has_os_time = FALSE,
        has_os_status = FALSE,
        survival_samples = 0,
        event_count = 0,
        usable_for_validation = FALSE,
        signature_availability = "not_available",
        notes = mapped$reason
      )
      write_current_outputs()
      next
    }

    expr_symbol_dt <- mapped$expr_symbol
    gene_symbols_present <- unique(expr_symbol_dt$gene_symbol)
    covered <- signature_genes[signature_genes %in% gene_symbols_present]
    missing <- setdiff(signature_genes, covered)
    coverage_n <- length(covered)

    os_info <- extract_os_info(sample_meta)
    usable <- coverage_n >= 7 && os_info$survival_samples >= 80 && os_info$event_count >= 30
    signature_availability <- if (coverage_n == 9) {
      "full_signature_available"
    } else if (coverage_n %in% c(7, 8)) {
      "partial_signature_available"
    } else {
      "not_recommended_for_main_validation"
    }

    feas_list[[length(feas_list) + 1]] <- data.table(
      gse = gse,
      matrix_file = mf,
      platform = platforms,
      series_matrix_downloaded = TRUE,
      expression_parsed = TRUE,
      sample_count = ncol(parsed$expr),
      log2_likely = log2_flag,
      gene_symbol_available = TRUE,
      gene_symbol_source = mapped$symbol_source,
      gene_symbol_column = mapped$symbol_column,
      signature_genes_covered = coverage_n,
      missing_signature_genes = paste(missing, collapse = ";"),
      has_os_time = os_info$has_os_time,
      has_os_status = os_info$has_os_status,
      survival_samples = os_info$survival_samples,
      event_count = os_info$event_count,
      usable_for_validation = usable,
      signature_availability = signature_availability,
      notes = os_info$notes
    )

    for (sg in signature_genes) {
      coverage_list[[length(coverage_list) + 1]] <- data.table(
        gse = gse,
        matrix_file = mf,
        platform = platforms,
        signature_gene = sg,
        covered = sg %in% covered
      )
    }
    write_current_outputs()
  }
}
write_current_outputs()

save_session_info("16_logs/sessionInfo_14_external_geo_dataset_feasibility_check.txt")
