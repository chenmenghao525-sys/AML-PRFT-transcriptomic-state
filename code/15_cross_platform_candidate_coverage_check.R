#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

dir.create("08_validation", recursive = TRUE, showWarnings = FALSE)
dir.create("14_tables", recursive = TRUE, showWarnings = FALSE)
dir.create("16_logs", recursive = TRUE, showWarnings = FALSE)

save_session_info <- function(path) {
  writeLines(capture.output(sessionInfo()), con = path)
}

required_files <- c(
  "06_wgcna/wgcna_deg_up_intersect_prft_module_genes.csv",
  "07_signature/main_candidate_genes_for_lasso.csv",
  "07_signature/final_prft_signature_coefficients.csv",
  "07_signature/univariate_cox_integrated_results.csv",
  "05_deg/deg_prft_high_vs_low_tcga_all.csv",
  "06_wgcna/wgcna_gene_level_statistics.csv",
  "08_validation/geo_feasibility_summary.csv",
  "08_validation/geo_signature_gene_coverage.csv",
  "08_validation/geo_download_manifest.csv"
)

missing_required <- required_files[!file.exists(required_files)]
if (length(missing_required) > 0) {
  stop("Missing required input files: ", paste(missing_required, collapse = "; "))
}

get_dataset_id <- function(gse, gpl) {
  paste(gse, gpl, sep = "_")
}

extract_gpl_from_text <- function(x) {
  hit <- regmatches(x, regexpr("GPL[0-9]+", x))
  if (length(hit) == 0 || is.na(hit[1]) || !nzchar(hit[1])) {
    NA_character_
  } else {
    hit[1]
  }
}

clean_symbol_string <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("///", ";", x, fixed = TRUE)
  x <- gsub("//", ";", x, fixed = TRUE)
  x <- gsub("\\|", ";", x)
  x <- gsub(",", ";", x, fixed = TRUE)
  x
}

parse_series_matrix <- function(localfile) {
  lines <- readLines(gzfile(localfile), warn = FALSE)
  begin_idx <- grep("^!series_matrix_table_begin", lines)
  end_idx <- grep("^!series_matrix_table_end", lines)

  if (length(begin_idx) == 0) {
    stop("series_matrix table not found in ", localfile)
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
    stop("Expression matrix has fewer than 2 columns in ", localfile)
  }

  expr_ids <- as.character(expr_df[[1]])
  expr_mat <- as.matrix(expr_df[, -1, drop = FALSE])
  storage.mode(expr_mat) <- "numeric"
  rownames(expr_mat) <- expr_ids
  expr_mat
}

parse_gpl_annotation <- function(localfile) {
  lines <- readLines(gzfile(localfile), warn = FALSE)
  begin_idx <- grep("^!platform_table_begin", lines)
  end_idx <- grep("^!platform_table_end", lines)

  if (length(begin_idx) == 0) {
    stop("GPL table not found in ", localfile)
  }

  if (length(end_idx) == 0 || end_idx[1] <= begin_idx[1]) {
    table_lines <- lines[(begin_idx[1] + 1):length(lines)]
  } else {
    table_lines <- lines[(begin_idx[1] + 1):(end_idx[1] - 1)]
  }

  read.delim(
    text = paste(table_lines, collapse = "\n"),
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    quote = "\"",
    comment.char = ""
  )
}

map_expr_to_symbol <- function(expr_mat, gpl_df = NULL) {
  row_ids <- rownames(expr_mat)

  if (is.null(gpl_df)) {
    gene_symbols <- unique(trimws(row_ids))
    gene_symbols <- gene_symbols[nzchar(gene_symbols) & !gene_symbols %in% c("NA", "---")]
    return(sort(unique(gene_symbols)))
  }

  norm_cols <- gsub("[^a-z0-9]", "", tolower(colnames(gpl_df)))
  symbol_candidates <- c("genesymbol", "symbol", "hgncsymbol", "gene_symbol")
  id_candidates <- c("id", "idref", "probeid")

  symbol_idx <- match(symbol_candidates, norm_cols)
  symbol_idx <- symbol_idx[!is.na(symbol_idx)][1]
  id_idx <- match(id_candidates, norm_cols)
  id_idx <- id_idx[!is.na(id_idx)][1]

  if (is.na(symbol_idx) || is.na(id_idx)) {
    stop("GPL annotation lacks recognizable probe ID or gene symbol column.")
  }

  gpl_map <- data.table(
    probe_id = as.character(gpl_df[[id_idx]]),
    raw_symbol = clean_symbol_string(gpl_df[[symbol_idx]])
  )
  gpl_map <- gpl_map[nzchar(probe_id) & nzchar(raw_symbol)]

  probe_stats <- data.table(
    probe_id = row_ids,
    probe_iqr = apply(expr_mat, 1, IQR, na.rm = TRUE)
  )

  expanded_map <- gpl_map[, .(gene_symbol = unlist(strsplit(raw_symbol, ";", fixed = TRUE))), by = "probe_id"]
  expanded_map[, gene_symbol := trimws(gene_symbol)]
  expanded_map <- unique(expanded_map[nzchar(gene_symbol) & !gene_symbol %in% c("NA", "---")])

  merged_map <- merge(expanded_map, probe_stats, by = "probe_id", all.x = FALSE, all.y = FALSE)
  if (nrow(merged_map) == 0) {
    stop("No overlap between expression probe IDs and GPL annotation.")
  }

  setorder(merged_map, gene_symbol, -probe_iqr)
  dedup_map <- merged_map[, .SD[1], by = gene_symbol]
  sort(unique(dedup_map$gene_symbol))
}

feas_dt <- fread("08_validation/geo_feasibility_summary.csv")
manifest_dt <- fread("08_validation/geo_download_manifest.csv")

target_datasets <- data.table(
  dataset_id = c(
    "GSE37642_GPL570", "GSE37642_GPL96",
    "GSE12417_GPL570", "GSE12417_GPL96",
    "GSE6891_GPL570", "GSE14468_GPL570"
  ),
  gse = c("GSE37642", "GSE37642", "GSE12417", "GSE12417", "GSE6891", "GSE14468"),
  gpl = c("GPL570", "GPL96", "GPL570", "GPL96", "GPL570", "GPL570")
)

resolve_dataset_row <- function(gse_value, gpl_value) {
  row <- feas_dt[
    gse == gse_value &
      grepl(gpl_value, matrix_file, fixed = TRUE)
  ]

  if (nrow(row) == 0 && gse_value %in% c("GSE6891", "GSE14468")) {
    row <- feas_dt[gse == gse_value]
  }

  if (nrow(row) == 0) {
    stop("No feasibility row found for ", gse_value, " / ", gpl_value)
  }

  row[1]
}

universe_rows <- list()
dataset_gene_sets <- list()

for (i in seq_len(nrow(target_datasets))) {
  gse_i <- target_datasets$gse[i]
  gpl_i <- target_datasets$gpl[i]
  dataset_id_i <- target_datasets$dataset_id[i]

  feas_row <- resolve_dataset_row(gse_i, gpl_i)
  matrix_file <- feas_row$matrix_file[[1]]
  matrix_path <- file.path("00_raw_data/geo_validation", matrix_file)
  if (!file.exists(matrix_path)) {
    stop("Missing series matrix file: ", matrix_path)
  }

  expr_mat <- parse_series_matrix(matrix_path)

  gpl_file <- paste0(gpl_i, "_family.soft.gz")
  gpl_path <- file.path("00_raw_data/geo_validation", gpl_file)
  if (!file.exists(gpl_path)) {
    stop("Missing GPL annotation file: ", gpl_path)
  }

  gpl_df <- parse_gpl_annotation(gpl_path)
  gene_universe <- map_expr_to_symbol(expr_mat, gpl_df)
  dataset_gene_sets[[dataset_id_i]] <- gene_universe

  universe_rows[[length(universe_rows) + 1]] <- data.table(
    dataset_id = dataset_id_i,
    GSE = gse_i,
    GPL = gpl_i,
    sample_count = as.integer(feas_row$sample_count[[1]]),
    gene_universe_size = length(gene_universe),
    gene_symbol = gene_universe
  )
}

gene_universe_dt <- rbindlist(universe_rows, fill = TRUE)
fwrite(gene_universe_dt, "08_validation/cross_platform_gene_universe_by_dataset.csv")

candidate_pool_715_dt <- fread("06_wgcna/wgcna_deg_up_intersect_prft_module_genes.csv")
lasso_input_67_dt <- fread("07_signature/main_candidate_genes_for_lasso.csv")
final_sig_9_dt <- fread("07_signature/final_prft_signature_coefficients.csv")
cox_dt <- fread("07_signature/univariate_cox_integrated_results.csv")
deg_all_dt <- fread("05_deg/deg_prft_high_vs_low_tcga_all.csv")
wgcna_gene_dt <- fread("06_wgcna/wgcna_gene_level_statistics.csv")

safe_unique_genes <- function(x) {
  sort(unique(trimws(as.character(x[nzchar(as.character(x)) & !is.na(x)]))))
}

gene_sets <- list(
  current_9_gene_signature = safe_unique_genes(final_sig_9_dt$gene_symbol),
  lasso_input_67_genes = safe_unique_genes(lasso_input_67_dt$gene_symbol),
  candidate_pool_715_genes = safe_unique_genes(candidate_pool_715_dt$gene_symbol),
  PRFT_high_upregulated_DEGs = safe_unique_genes(
    deg_all_dt[adj.P.Val < 0.05 & logFC > 0.5]$gene_symbol
  ),
  candidate_pool_HR_gt1_CoxP_lt0.05 = safe_unique_genes(
    cox_dt[HR > 1 & P.Value < 0.05]$gene_symbol
  ),
  candidate_pool_HR_gt1_CoxFDR_lt0.25 = safe_unique_genes(
    cox_dt[HR > 1 & FDR < 0.25]$gene_symbol
  ),
  strict_candidate_pool = safe_unique_genes(
    candidate_pool_715_dt[
      logFC > 0.5 &
        adj.P.Val < 0.05 &
        GS_PRFT > 0.30 &
        MM > 0.50
    ]$gene_symbol
  )
)

strict_with_cox <- safe_unique_genes(
  merge(
    candidate_pool_715_dt[
      logFC > 0.5 &
        adj.P.Val < 0.05 &
        GS_PRFT > 0.30 &
        MM > 0.50,
      .(gene_symbol, logFC, adj.P.Val, module_color, GS_PRFT, MM)
    ],
    cox_dt[HR > 1 & P.Value < 0.05, .(gene_symbol)],
    by = "gene_symbol"
  )$gene_symbol
)
gene_sets[["strict_candidate_pool_with_cox"]] <- strict_with_cox

coverage_rows <- list()

for (i in seq_len(nrow(target_datasets))) {
  dataset_id_i <- target_datasets$dataset_id[i]
  gse_i <- target_datasets$gse[i]
  gpl_i <- target_datasets$gpl[i]
  feas_row <- resolve_dataset_row(gse_i, gpl_i)
  universe_i <- dataset_gene_sets[[dataset_id_i]]

  for (gene_set_name in names(gene_sets)) {
    input_genes <- gene_sets[[gene_set_name]]
    covered_genes <- sort(intersect(input_genes, universe_i))
    missing_genes <- sort(setdiff(input_genes, universe_i))

    coverage_rows[[length(coverage_rows) + 1]] <- data.table(
      dataset_id = dataset_id_i,
      GSE = gse_i,
      GPL = gpl_i,
      sample_count = as.integer(feas_row$sample_count[[1]]),
      OS_available = isTRUE(feas_row$has_os_time[[1]]) && isTRUE(feas_row$has_os_status[[1]]) && as.integer(feas_row$survival_samples[[1]]) > 0,
      OS_complete_samples = as.integer(feas_row$survival_samples[[1]]),
      events = as.integer(feas_row$event_count[[1]]),
      gene_set_name = gene_set_name,
      input_gene_count = length(input_genes),
      covered_gene_count = length(covered_genes),
      missing_gene_count = length(missing_genes),
      covered_gene_list = paste(covered_genes, collapse = ";"),
      missing_gene_list = paste(missing_genes, collapse = ";")
    )
  }
}

coverage_by_dataset_dt <- rbindlist(coverage_rows, fill = TRUE)
fwrite(coverage_by_dataset_dt, "08_validation/cross_platform_candidate_coverage_by_dataset.csv")

scenario_defs <- list(
  Scenario_A_GPL96_two_cohorts = c("GSE37642_GPL96", "GSE12417_GPL96"),
  Scenario_B_GPL570_two_cohorts = c("GSE37642_GPL570", "GSE12417_GPL570"),
  Scenario_C_GSE37642_GPL96_only = c("GSE37642_GPL96"),
  Scenario_D_GSE37642_GPL570_only = c("GSE37642_GPL570"),
  Scenario_E_GSE12417_GPL96_only = c("GSE12417_GPL96"),
  Scenario_F_GSE12417_GPL570_only = c("GSE12417_GPL570")
)

scenario_rows <- list()
scenario_stats <- list()

for (scenario_name in names(scenario_defs)) {
  ds_ids <- scenario_defs[[scenario_name]]
  scenario_universe <- Reduce(intersect, dataset_gene_sets[ds_ids])

  stat_row <- list(
    scenario = scenario_name,
    datasets_included = paste(ds_ids, collapse = ";"),
    candidate_pool_covered_count = length(intersect(gene_sets$candidate_pool_715_genes, scenario_universe)),
    risk_count = length(intersect(gene_sets$candidate_pool_HR_gt1_CoxP_lt0.05, scenario_universe)),
    strict_count = length(intersect(gene_sets$strict_candidate_pool_with_cox, scenario_universe))
  )
  scenario_stats[[scenario_name]] <- stat_row

  for (gene_set_name in names(gene_sets)) {
    covered_genes <- sort(intersect(gene_sets[[gene_set_name]], scenario_universe))
    scenario_rows[[length(scenario_rows) + 1]] <- data.table(
      scenario = scenario_name,
      datasets_included = paste(ds_ids, collapse = ";"),
      gene_set_name = gene_set_name,
      input_gene_count = length(gene_sets[[gene_set_name]]),
      scenario_common_covered_gene_count = length(covered_genes),
      scenario_common_covered_gene_list = paste(covered_genes, collapse = ";")
    )
  }
}

coverage_by_scenario_dt <- rbindlist(scenario_rows, fill = TRUE)
fwrite(coverage_by_scenario_dt, "08_validation/cross_platform_candidate_coverage_by_scenario.csv")

priority_order <- c(
  "Scenario_A_GPL96_two_cohorts",
  "Scenario_B_GPL570_two_cohorts",
  "Scenario_C_GSE37642_GPL96_only",
  "Scenario_D_GSE37642_GPL570_only"
)

recommended_scenario <- "none"
recommended_reason <- "No scenario met the predefined rebuild criteria."
rebuild_recommended <- FALSE

for (scenario_name in priority_order) {
  st <- scenario_stats[[scenario_name]]
  if (!is.null(st) && st$candidate_pool_covered_count >= 100 && st$risk_count >= 20) {
    recommended_scenario <- scenario_name
    recommended_reason <- "First priority scenario meeting candidate_pool_715 >= 100 and HR > 1 with Cox P < 0.05 >= 20."
    rebuild_recommended <- TRUE
    break
  }
}

if (!rebuild_recommended) {
  best_priority <- priority_order[which.max(vapply(priority_order, function(x) scenario_stats[[x]]$risk_count, numeric(1)))]
  best_stats <- scenario_stats[[best_priority]]
  recommended_reason <- paste0(
    "No scenario met the predefined rebuild criteria. Best priority scenario by available risk genes was ",
    best_priority,
    " (candidate_pool_covered_count=",
    best_stats$candidate_pool_covered_count,
    ", risk_count=",
    best_stats$risk_count,
    ", strict_count=",
    best_stats$strict_count,
    ")."
  )
}

scenario_for_output <- if (rebuild_recommended) recommended_scenario else priority_order[which.max(vapply(priority_order, function(x) scenario_stats[[x]]$risk_count, numeric(1)))]
datasets_for_output <- scenario_defs[[scenario_for_output]]
scenario_output_universe <- Reduce(intersect, dataset_gene_sets[datasets_for_output])

candidate_pool_annotated <- merge(
  candidate_pool_715_dt[, .(gene_symbol, module_color, GS_PRFT, MM, logFC, adj.P.Val)],
  cox_dt[, .(gene_symbol, HR, Cox_P.Value = P.Value, Cox_FDR = FDR)],
  by = "gene_symbol",
  all.x = TRUE
)
candidate_pool_annotated <- unique(candidate_pool_annotated, by = "gene_symbol")
candidate_pool_annotated[, covered_in_selected_scenario := gene_symbol %in% scenario_output_universe]
candidate_pool_annotated[, covered_in_datasets := paste(datasets_for_output[gene_symbol %in% dataset_gene_sets[[datasets_for_output[1]]]], collapse = ";")]

covered_flags <- lapply(datasets_for_output, function(ds) candidate_pool_annotated$gene_symbol %in% dataset_gene_sets[[ds]])
covered_matrix <- do.call(cbind, covered_flags)
if (is.null(dim(covered_matrix))) {
  covered_matrix <- matrix(covered_matrix, ncol = 1)
}
candidate_pool_annotated[, covered_in_datasets := apply(covered_matrix, 1, function(x) paste(datasets_for_output[as.logical(x)], collapse = ";"))]

rebuild_candidate_dt <- candidate_pool_annotated[covered_in_selected_scenario == TRUE]
rebuild_candidate_dt[, recommended_scenario := scenario_for_output]
rebuild_candidate_dt[, risk_gene_flag := !is.na(Cox_P.Value) & Cox_P.Value < 0.05 & !is.na(HR) & HR > 1]
rebuild_candidate_dt[, strict_gene_flag := !is.na(logFC) & logFC > 0.5 & !is.na(adj.P.Val) & adj.P.Val < 0.05 & !is.na(GS_PRFT) & GS_PRFT > 0.30 & !is.na(MM) & MM > 0.50]
setcolorder(
  rebuild_candidate_dt,
  c("gene_symbol", "logFC", "adj.P.Val", "HR", "Cox_P.Value", "Cox_FDR", "module_color", "GS_PRFT", "MM", "covered_in_datasets", "recommended_scenario")
)

setorder(
  rebuild_candidate_dt,
  -risk_gene_flag,
  -strict_gene_flag,
  Cox_P.Value,
  -GS_PRFT,
  -MM
)
rebuild_candidate_dt[, c("risk_gene_flag", "strict_gene_flag") := NULL]

fwrite(rebuild_candidate_dt, "08_validation/cross_platform_rebuild_candidate_genes.csv")
writeLines(rebuild_candidate_dt$gene_symbol, con = "08_validation/cross_platform_rebuild_candidate_genes.txt")

strategy_dt <- data.table(
  recommended_scenario = recommended_scenario,
  reason = recommended_reason,
  candidate_pool_covered_count = scenario_stats[[scenario_for_output]]$candidate_pool_covered_count,
  risk_cox_P_less_0.05_HR_greater_1_count = scenario_stats[[scenario_for_output]]$risk_count,
  strict_candidate_count = scenario_stats[[scenario_for_output]]$strict_count,
  datasets_for_future_validation = paste(datasets_for_output, collapse = ";"),
  whether_rebuild_lasso_is_recommended = rebuild_recommended
)
fwrite(strategy_dt, "08_validation/cross_platform_rebuild_recommended_strategy.csv")

summary_dt <- data.table(
  metric = c(
    "GSE37642_GPL96_candidate_pool_715_covered",
    "GSE12417_GPL96_candidate_pool_715_covered",
    "Scenario_A_GPL96_two_cohorts_candidate_pool_715_common",
    "GSE37642_GPL570_candidate_pool_715_covered",
    "GSE12417_GPL570_candidate_pool_715_covered",
    "Scenario_B_GPL570_two_cohorts_candidate_pool_715_common",
    "recommended_scenario_risk_gene_count",
    "recommended_scenario_strict_gene_count"
  ),
  value = c(
    coverage_by_dataset_dt[dataset_id == "GSE37642_GPL96" & gene_set_name == "candidate_pool_715_genes"]$covered_gene_count[1],
    coverage_by_dataset_dt[dataset_id == "GSE12417_GPL96" & gene_set_name == "candidate_pool_715_genes"]$covered_gene_count[1],
    coverage_by_scenario_dt[scenario == "Scenario_A_GPL96_two_cohorts" & gene_set_name == "candidate_pool_715_genes"]$scenario_common_covered_gene_count[1],
    coverage_by_dataset_dt[dataset_id == "GSE37642_GPL570" & gene_set_name == "candidate_pool_715_genes"]$covered_gene_count[1],
    coverage_by_dataset_dt[dataset_id == "GSE12417_GPL570" & gene_set_name == "candidate_pool_715_genes"]$covered_gene_count[1],
    coverage_by_scenario_dt[scenario == "Scenario_B_GPL570_two_cohorts" & gene_set_name == "candidate_pool_715_genes"]$scenario_common_covered_gene_count[1],
    scenario_stats[[scenario_for_output]]$risk_count,
    scenario_stats[[scenario_for_output]]$strict_count
  )
)
fwrite(summary_dt, "14_tables/cross_platform_candidate_coverage_summary.csv")

save_session_info("16_logs/sessionInfo_15_cross_platform_candidate_coverage_check.txt")
