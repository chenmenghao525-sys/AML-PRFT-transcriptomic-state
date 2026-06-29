#!/usr/bin/env Rscript

# TCGA-LAML download script
# - Downloads RNA-seq expression data from GDC using TCGAbiolinks
# - Prefers STAR-counts and falls back to HTSeq-counts if needed
# - Downloads raw clinical data
# - Saves raw objects only; no downstream preprocessing is performed here

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(data.table)
})

options(stringsAsFactors = FALSE)

project_id <- "TCGA-LAML"

raw_dir <- file.path("00_raw_data", "tcga_laml_raw")
meta_dir <- "01_metadata"
log_dir <- "16_logs"

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

save_session_info <- function(path) {
  writeLines(capture.output(sessionInfo()), con = path)
}

query_and_prepare <- function(project_id) {
  workflow_candidates <- c("STAR - Counts", "HTSeq - Counts")
  error_log <- list()

  for (wf in workflow_candidates) {
    message("Trying workflow.type = ", wf)

    query <- tryCatch(
      GDCquery(
        project = project_id,
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        workflow.type = wf
      ),
      error = function(e) {
        message("Query failed for ", wf, ": ", conditionMessage(e))
        return(NULL)
      }
    )

    if (is.null(query)) {
      next
    }

    if (nrow(getResults(query)) == 0) {
      message("No records returned for ", wf)
      next
    }

    message("Downloading expression files for ", wf)
    download_ok <- tryCatch({
      GDCdownload(
        query,
        method = "api",
        files.per.chunk = 20
      )
      TRUE
    }, error = function(e) {
      msg <- paste0("Download failed for ", wf, ": ", conditionMessage(e))
      message(msg)
      error_log[[length(error_log) + 1]] <<- msg
      FALSE
    })

    if (!download_ok) {
      next
    }

    message("Preparing SummarizedExperiment for ", wf)
    se <- tryCatch(
      GDCprepare(query),
      error = function(e) {
        msg <- paste0("GDCprepare failed for ", wf, ": ", conditionMessage(e))
        message(msg)
        error_log[[length(error_log) + 1]] <<- msg
        return(NULL)
      }
    )

    if (is.null(se)) {
      next
    }

    metadata(se)$download_workflow_type <- wf
    return(se)
  }

  stop(
    paste(
      c("No usable TCGA-LAML RNA-seq workflow was available.", unlist(error_log)),
      collapse = "\n"
    )
  )
}

message("Starting TCGA-LAML download")

tcga_se <- query_and_prepare(project_id)
saveRDS(tcga_se, file = file.path(raw_dir, "tcga_laml_se.rds"))
message("Saved raw SummarizedExperiment to ", file.path(raw_dir, "tcga_laml_se.rds"))

message("Downloading raw clinical data")
clinical_raw <- tryCatch(
  GDCquery_clinic(project = project_id, type = "clinical"),
  error = function(e) {
    stop("Failed to download TCGA clinical data: ", conditionMessage(e))
  }
)

fwrite(
  as.data.table(clinical_raw),
  file = file.path(meta_dir, "clinical_tcga_raw.csv"),
  na = ""
)
message("Saved raw clinical data to ", file.path(meta_dir, "clinical_tcga_raw.csv"))

save_session_info(file.path(log_dir, "sessionInfo_01_download_tcga_laml.txt"))
message("Saved sessionInfo")
