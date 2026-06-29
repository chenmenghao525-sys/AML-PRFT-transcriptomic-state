#!/usr/bin/env Rscript

# Download TCGA-LAML expression and clinical/survival files from UCSC Xena / GDC Xena Hub
# without using TCGAbiolinks.

suppressPackageStartupMessages({
  library(data.table)
})

options(stringsAsFactors = FALSE)
options(timeout = 1200)

raw_dir <- file.path("00_raw_data", "tcga_laml_raw_xena")
log_dir <- "16_logs"
log_file <- file.path(log_dir, "download_tcga_laml_xena_log.txt")
session_file <- file.path(log_dir, "sessionInfo_01_download_tcga_laml_xena.txt")
manifest_file <- file.path(raw_dir, "tcga_laml_xena_manifest.csv")

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_message <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

save_session_info <- function(path) {
  writeLines(capture.output(sessionInfo()), con = path)
}

probe_url <- function(url) {
  con <- NULL
  ok <- FALSE
  err <- NA_character_
  tryCatch({
    con <- url(url, open = "rb", encoding = "bytes")
    readBin(con, what = "raw", n = 64)
    ok <- TRUE
  }, error = function(e) {
    err <<- conditionMessage(e)
  }, finally = {
    if (!is.null(con)) {
      try(close(con), silent = TRUE)
    }
  })
  list(ok = ok, error = err)
}

download_with_fallback <- function(candidates, role, required = TRUE) {
  errors <- character(0)
  for (i in seq_len(nrow(candidates))) {
    row <- candidates[i, ]
    url <- row$url
    dest <- file.path(raw_dir, row$dest_name)

    log_message("Testing [", role, "] ", url)
    probe <- probe_url(url)
    if (!isTRUE(probe$ok)) {
      msg <- paste0("Probe failed: ", url, " | ", probe$error)
      log_message(msg)
      errors <- c(errors, msg)
      next
    }

    log_message("Downloading [", role, "] ", url)
    dl_ok <- tryCatch({
      download.file(url, destfile = dest, mode = "wb", quiet = FALSE, method = "libcurl")
      TRUE
    }, error = function(e) {
      msg <- paste0("Download failed: ", url, " | ", conditionMessage(e))
      log_message(msg)
      errors <<- c(errors, msg)
      FALSE
    })

    if (!dl_ok || !file.exists(dest)) {
      next
    }

    file_info <- file.info(dest)
    if (is.na(file_info$size) || file_info$size <= 0) {
      msg <- paste0("Downloaded file is empty: ", dest)
      log_message(msg)
      errors <- c(errors, msg)
      next
    }

    log_message(
      "Downloaded [", role, "] to ", normalizePath(dest, winslash = "/"),
      " | size=", format(file_info$size, big.mark = ",")
    )

    return(list(
      success = TRUE,
      role = role,
      source = row$source,
      data_type = row$data_type,
      url = url,
      dest = dest,
      size = unname(file_info$size),
      errors = paste(errors, collapse = " || ")
    ))
  }

  if (required) {
    stop(
      paste(
        c(paste0("All candidates failed for role: ", role), errors),
        collapse = "\n"
      )
    )
  }

  list(
    success = FALSE,
    role = role,
    source = NA_character_,
    data_type = NA_character_,
    url = NA_character_,
    dest = NA_character_,
    size = NA_real_,
    errors = paste(errors, collapse = " || ")
  )
}

expression_candidates <- rbindlist(list(
  data.table(
    source = "GDC Xena Hub",
    data_type = "HTSeq counts",
    url = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LAML.htseq_counts.tsv.gz",
    dest_name = "TCGA-LAML.htseq_counts.tsv.gz"
  ),
  data.table(
    source = "GDC Xena Hub",
    data_type = "HTSeq counts",
    url = "https://gdc.xenahubs.net/download/TCGA-LAML.htseq_counts.tsv.gz",
    dest_name = "TCGA-LAML.htseq_counts.tsv.gz"
  ),
  data.table(
    source = "GDC Xena Hub",
    data_type = "FPKM-UQ",
    url = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LAML.htseq_fpkm-uq.tsv.gz",
    dest_name = "TCGA-LAML.htseq_fpkm-uq.tsv.gz"
  ),
  data.table(
    source = "GDC Xena Hub",
    data_type = "FPKM-UQ",
    url = "https://gdc.xenahubs.net/download/TCGA-LAML.htseq_fpkm-uq.tsv.gz",
    dest_name = "TCGA-LAML.htseq_fpkm-uq.tsv.gz"
  ),
  data.table(
    source = "GDC Xena Hub",
    data_type = "FPKM",
    url = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LAML.htseq_fpkm.tsv.gz",
    dest_name = "TCGA-LAML.htseq_fpkm.tsv.gz"
  ),
  data.table(
    source = "GDC Xena Hub",
    data_type = "FPKM",
    url = "https://gdc.xenahubs.net/download/TCGA-LAML.htseq_fpkm.tsv.gz",
    dest_name = "TCGA-LAML.htseq_fpkm.tsv.gz"
  ),
  data.table(
    source = "TCGA Xena Hub",
    data_type = "normalized",
    url = "https://tcga.xenahubs.net/download/TCGA.LAML.sampleMap/HiSeqV2.gz",
    dest_name = "TCGA.LAML.sampleMap_HiSeqV2.gz"
  )
), fill = TRUE)

phenotype_candidates <- rbindlist(list(
  data.table(
    source = "TCGA Xena Hub",
    data_type = "clinicalMatrix",
    url = "https://tcga.xenahubs.net/download/TCGA.LAML.sampleMap/LAML_clinicalMatrix",
    dest_name = "TCGA.LAML.sampleMap_LAML_clinicalMatrix.tsv"
  ),
  data.table(
    source = "GDC Xena Hub",
    data_type = "phenotype",
    url = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LAML.GDC_phenotype.tsv.gz",
    dest_name = "TCGA-LAML.GDC_phenotype.tsv.gz"
  ),
  data.table(
    source = "GDC Xena Hub",
    data_type = "phenotype",
    url = "https://gdc.xenahubs.net/download/TCGA-LAML.GDC_phenotype.tsv.gz",
    dest_name = "TCGA-LAML.GDC_phenotype.tsv.gz"
  )
), fill = TRUE)

survival_candidates <- rbindlist(list(
  data.table(
    source = "GDC Xena Hub",
    data_type = "survival",
    url = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LAML.survival.tsv",
    dest_name = "TCGA-LAML.survival.tsv"
  ),
  data.table(
    source = "GDC Xena Hub",
    data_type = "survival",
    url = "https://gdc.xenahubs.net/download/TCGA-LAML.survival.tsv",
    dest_name = "TCGA-LAML.survival.tsv"
  )
), fill = TRUE)

if (file.exists(log_file)) {
  file.remove(log_file)
}
log_message("Starting TCGA-LAML Xena download workflow")

expr_res <- download_with_fallback(expression_candidates, role = "expression", required = TRUE)
if (!identical(expr_res$data_type, "HTSeq counts")) {
  log_message("Expression fallback activated: using ", expr_res$data_type, " instead of HTSeq counts")
}

pheno_res <- download_with_fallback(phenotype_candidates, role = "phenotype", required = FALSE)
surv_res <- download_with_fallback(survival_candidates, role = "survival", required = FALSE)

manifest <- rbindlist(list(
  as.data.table(expr_res),
  as.data.table(pheno_res),
  as.data.table(surv_res)
), fill = TRUE)

manifest[, download_time := format(Sys.time(), "%Y-%m-%d %H:%M:%S")]
fwrite(manifest, file = manifest_file)

save_session_info(session_file)
log_message("Saved manifest to ", normalizePath(manifest_file, winslash = "/"))
log_message("Saved sessionInfo to ", normalizePath(session_file, winslash = "/"))
log_message("Completed TCGA-LAML Xena download workflow")
