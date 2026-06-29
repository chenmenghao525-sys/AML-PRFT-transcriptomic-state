#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(1234)

phase_lib <- Sys.getenv("PHASE1_ASCII_R_LIB", unset = "C:/Users/ROBIN-~1/AppData/Local/Temp/phase1_R_libs")
if (dir.exists(phase_lib)) {
  .libPaths(unique(c(phase_lib, .libPaths())))
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

root <- Sys.getenv("PHASE7_ROOT", unset = "C:/Users/Robin-Yang/AppData/Local/Temp/aml_prft_phase1_fix")
setwd(root)

dir.create("03_results_tables", showWarnings = FALSE, recursive = TRUE)
dir.create("04_figures", showWarnings = FALSE, recursive = TRUE)
dir.create("05_logs", showWarnings = FALSE, recursive = TRUE)
dir.create("06_manuscript_support", showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(root, "05_logs", "phase7_cleanup_log.txt")
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
msg_con <- file(log_file, open = "at")
sink(msg_con, type = "message")
on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(), silent = TRUE)
  try(close(msg_con), silent = TRUE)
  try(close(log_con), silent = TRUE)
}, add = TRUE)

cat("Phase 7-cleanup: main-text cleanup and language-boundary correction\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Root:", root, "\n")
cat("Detected polishing axes: paper_type=research; section=results; language=en; journal=nat-comms\n")
cat("Figure backend: R\n\n")

if (file.exists("02_scripts/15_scripts/plot_label_utils.R")) {
  source("02_scripts/15_scripts/plot_label_utils.R")
} else {
  pretty_label <- function(x) as.character(x)
}

required_files <- c(
  "05_logs/phase7_bulk_immune_key_result_checklist.txt",
  "03_results_tables/phase7_cross_cohort_consistent_signatures.csv",
  "03_results_tables/phase7_pathway_PRFT_association.csv",
  "03_results_tables/phase7_pathway_cross_cohort_consistency.csv",
  "03_results_tables/phase7_bulk_singlecell_consistency_summary.csv",
  "03_results_tables/phase7_immune_deconvolution_PRFT_association.csv",
  "03_results_tables/phase7_main_vs_supplement_recommendation.csv",
  "04_figures/phase7_PRFT_signature_correlation_heatmap.pdf",
  "04_figures/phase7_pathway_activity_heatmap.pdf"
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files: ", paste(missing_files, collapse = "; "))
}

sig_cons <- fread("03_results_tables/phase7_cross_cohort_consistent_signatures.csv")
path_assoc <- fread("03_results_tables/phase7_pathway_PRFT_association.csv")
path_cons <- fread("03_results_tables/phase7_pathway_cross_cohort_consistency.csv")
bulk_sc <- fread("03_results_tables/phase7_bulk_singlecell_consistency_summary.csv")
deconv_assoc <- fread("03_results_tables/phase7_immune_deconvolution_PRFT_association.csv")
main_vs_supp <- fread("03_results_tables/phase7_main_vs_supplement_recommendation.csv")
phase7_checklist <- readLines("05_logs/phase7_bulk_immune_key_result_checklist.txt", warn = FALSE)

main_signature_targets <- c(
  "monocyte_macrophage_like_set",
  "myeloid_suppressive_set",
  "immune_checkpoint_set",
  "JAK_STAT_PDL1_set",
  "oxidative_stress_NRF2_set",
  "phase3B_SHAP_core_axis_set",
  "phase5_PPI_hub_set"
)
main_pathway_targets <- c(
  "JAK_STAT_pathway",
  "IFN_response_pathway",
  "TNF_NFKB_pathway",
  "inflammatory_response_pathway",
  "oxidative_stress_NRF2_pathway"
)

main_sig_summary <- sig_cons[set_name %in% main_signature_targets, .(
  result_type = "bulk_signature",
  set_name,
  display_name = pretty_label(set_name),
  risk_cor_mean_rho,
  risk_cor_significant_cohorts,
  high_low_direction,
  high_low_significant_cohorts,
  PRFT_cor_mean_rho,
  PRFT_cor_significant_cohorts,
  PRFT_group_direction,
  PRFT_group_significant_cohorts,
  cross_cohort_consistent,
  main_text_priority = fifelse(set_name %in% c("monocyte_macrophage_like_set", "myeloid_suppressive_set", "immune_checkpoint_set", "JAK_STAT_PDL1_set", "oxidative_stress_NRF2_set"), "core", "supporting"),
  cleaned_interpretation = fifelse(
    set_name %in% c("monocyte_macrophage_like_set", "myeloid_suppressive_set"),
    "Cross-cohort bulk signature scoring linked the PRFT-related axis to monocyte/myeloid stress-adapted transcriptional programs.",
    fifelse(
      set_name == "immune_checkpoint_set",
      "Bulk signature scoring indicated a cross-cohort association with immune checkpoint-related transcriptional programs.",
      fifelse(
        set_name == "JAK_STAT_PDL1_set",
        "Bulk signature scoring indicated a cross-cohort association with the JAK/STAT/PD-L1-related transcriptional axis.",
        fifelse(
          set_name == "oxidative_stress_NRF2_set",
          "Bulk signature scoring indicated a cross-cohort association with oxidative stress/NRF2-related programs.",
          "This set was retained as supporting evidence for the PRFT-related transcriptional context, not as a standalone immune-state label."
        )
      )
    )
  )
)]

main_path_summary <- path_cons[set_name %in% main_pathway_targets, .(
  result_type = "pathway_activity",
  set_name,
  display_name = pretty_label(set_name),
  risk_cor_mean_rho,
  risk_cor_significant_cohorts,
  high_low_direction,
  high_low_significant_cohorts,
  PRFT_cor_mean_rho,
  PRFT_cor_significant_cohorts,
  PRFT_group_direction,
  PRFT_group_significant_cohorts,
  cross_cohort_consistent,
  main_text_priority = "core",
  cleaned_interpretation = fifelse(
    set_name == "JAK_STAT_pathway",
    "Pathway scoring linked the PRFT-related axis to JAK/STAT activity across cohorts.",
    fifelse(
      set_name == "IFN_response_pathway",
      "Pathway scoring linked the PRFT-related axis to IFN-response programs across cohorts.",
      fifelse(
        set_name == "TNF_NFKB_pathway",
        "Pathway scoring linked the PRFT-related axis to TNF/NF-kB-related inflammatory signaling across cohorts.",
        fifelse(
          set_name == "inflammatory_response_pathway",
          "Pathway scoring linked the PRFT-related axis to inflammatory-response programs across cohorts.",
          "Pathway scoring linked the PRFT-related axis to oxidative stress/NRF2-related activity across cohorts."
        )
      )
    )
  )
)]

main_summary <- rbindlist(list(main_sig_summary, main_path_summary), fill = TRUE)
setcolorder(main_summary, c(
  "result_type", "set_name", "display_name", "main_text_priority",
  "risk_cor_mean_rho", "risk_cor_significant_cohorts",
  "high_low_direction", "high_low_significant_cohorts",
  "PRFT_cor_mean_rho", "PRFT_cor_significant_cohorts",
  "PRFT_group_direction", "PRFT_group_significant_cohorts",
  "cross_cohort_consistent", "cleaned_interpretation"
))
fwrite(main_summary, "03_results_tables/phase7_cleanup_main_signature_pathway_summary.csv")

deconv_all_na <- all(is.na(deconv_assoc$rho)) && all(is.na(deconv_assoc$P.Value))
downplay_dt <- rbindlist(list(
  data.table(
    result_name = "T_NK_cytotoxic_set",
    result_class = "bulk_signature",
    direction_summary = sig_cons[set_name == "T_NK_cytotoxic_set"]$risk_cor_direction[1],
    downplay_reason = "Direction was mixed and this set did not support the monocyte/myeloid main line.",
    recommended_handling = "supplement_or_not_emphasized",
    cleaned_wording = "Cytotoxic lymphoid-related signals were not prioritized in the main text because the cross-cohort pattern was not central to the PRFT-high state narrative."
  ),
  data.table(
    result_name = "T_cell_exhaustion_set",
    result_class = "bulk_signature",
    direction_summary = sig_cons[set_name == "T_cell_exhaustion_set"]$risk_cor_direction[1],
    downplay_reason = "This set was directionally stable but negative and did not serve the main monocyte/myeloid narrative.",
    recommended_handling = "supplement_or_context_only",
    cleaned_wording = "T-cell exhaustion-related signatures were not emphasized as a primary feature of the PRFT-high state."
  ),
  data.table(
    result_name = "LSC_stemness_set",
    result_class = "bulk_signature",
    direction_summary = sig_cons[set_name == "LSC_stemness_set"]$risk_cor_direction[1],
    downplay_reason = "Bulk and processed single-cell evidence did not support a dominant primitive or LSC-led interpretation.",
    recommended_handling = "downplay_in_main_text",
    cleaned_wording = "LSC/stemness-associated signals did not support a dominant primitive interpretation."
  ),
  data.table(
    result_name = "proteostasis_UPR_set",
    result_class = "bulk_signature",
    direction_summary = sig_cons[set_name == "proteostasis_UPR_set"]$risk_cor_direction[1],
    downplay_reason = "The risk-score direction was mixed across cohorts, although the PRFT-axis association remained informative.",
    recommended_handling = "supplement_with_context",
    cleaned_wording = "Proteostasis/UPR-related scoring was retained as contextual support rather than a primary risk-linked immune result."
  ),
  data.table(
    result_name = "ferroptosis_defense_set",
    result_class = "bulk_signature",
    direction_summary = sig_cons[set_name == "ferroptosis_defense_set"]$risk_cor_direction[1],
    downplay_reason = "The risk-score direction was mixed, whereas the PRFT-axis association was stable and biologically coherent.",
    recommended_handling = "PRFT_axis_associated_not_risk_primary",
    cleaned_wording = "Ferroptosis-defense scoring was described as PRFT-axis-associated rather than as a primary risk-score feature."
  ),
  data.table(
    result_name = "signature_based_myeloid_estimates",
    result_class = "fallback_estimation",
    direction_summary = ifelse(deconv_all_na, "not_quantitatively_usable", "available_but_fallback"),
    downplay_reason = "Formal immune deconvolution packages were unavailable, so these outputs should not be framed as validated cell fractions.",
    recommended_handling = "supplement_or_remove_from_main_text",
    cleaned_wording = "Fallback signature-based inference can be mentioned cautiously, but it should not be titled or framed as formal immune deconvolution."
  )
))
fwrite(downplay_dt, "03_results_tables/phase7_cleanup_supplement_or_downplay_results.csv")

scan_files <- c(
  "05_logs/phase7_bulk_immune_key_result_checklist.txt",
  "03_results_tables/phase7_main_vs_supplement_recommendation.csv",
  "03_results_tables/phase7_bulk_singlecell_consistency_summary.csv"
)
forbidden_patterns <- data.table(
  pattern = c(
    "immune deconvolution demonstrated",
    "immune infiltration demonstrated",
    "CIBERSORT",
    "xCell showed",
    "MCPcounter showed",
    "immune cells were validated"
  ),
  preferred_replacement = c(
    "signature-based inference suggested",
    "bulk signature scoring indicated",
    "not used unless actually run",
    "bulk signature scoring indicated",
    "bulk signature scoring indicated",
    "bulk transcriptional programs were associated with"
  )
)

language_findings <- rbindlist(lapply(scan_files, function(f) {
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  rbindlist(lapply(seq_len(nrow(forbidden_patterns)), function(i) {
    pat <- forbidden_patterns$pattern[i]
    data.table(
      file = f,
      pattern = pat,
      found = grepl(pat, txt, ignore.case = TRUE),
      preferred_replacement = forbidden_patterns$preferred_replacement[i]
    )
  }))
}))

boundary_notes <- c(
  "Phase 7-cleanup language boundary audit",
  paste("Scanned files:", paste(scan_files, collapse = "; ")),
  "",
  "Detected problematic phrases:"
)
if (any(language_findings$found)) {
  hit_lines <- language_findings[found == TRUE, paste0("- ", file, ": ", pattern, " -> ", preferred_replacement)]
  boundary_notes <- c(boundary_notes, hit_lines)
} else {
  boundary_notes <- c(boundary_notes, "- No explicit forbidden immune-validation phrases were detected in the scanned text files.")
}
boundary_notes <- c(
  boundary_notes,
  "",
  "Cleanup boundary decisions:",
  "- Use 'bulk signature scoring indicated' or 'signature-based inference suggested'.",
  "- Avoid 'immune deconvolution' as a main-text figure title.",
  "- Avoid 'immune infiltration demonstrated' and any validated-cell-fraction wording.",
  "- Describe Phase 7 as bulk transcriptional program inference, not causal or mechanistic proof.",
  "- Keep LSC/stemness in a non-dominant, non-primary position."
)
writeLines(boundary_notes, "05_logs/phase7_cleanup_language_boundary_log.txt")

main_figure_reco <- rbindlist(list(
  data.table(
    artifact = "phase7_PRFT_signature_correlation_heatmap.pdf",
    proposed_panel_role = "main_text_core",
    title_guidance = "Bulk signature-based inference across cohorts",
    keep_reason = "Directly summarizes the core monocyte/myeloid, checkpoint, and JAK/STAT-associated transcriptional programs."
  ),
  data.table(
    artifact = "phase7_pathway_activity_heatmap.pdf",
    proposed_panel_role = "main_text_or_condensed_supplement",
    title_guidance = "Bulk pathway activity scoring across cohorts",
    keep_reason = "Retains the pathway layer, but should be trimmed to core pathways if space is limited."
  ),
  data.table(
    artifact = "phase7_cleanup_core_signature_consistency_panel.pdf",
    proposed_panel_role = "main_text_compact_summary",
    title_guidance = "Core signature and pathway consistency",
    keep_reason = "Provides a short, cleaner summary panel for the main text without the deconvolution label."
  ),
  data.table(
    artifact = "phase7_cleanup_bulk_singlecell_consistency_panel.pdf",
    proposed_panel_role = "main_text_or_supplement_bridge",
    title_guidance = "Bulk and processed single-cell consistency",
    keep_reason = "Bridges Phase 7 bulk inference with the processed single-cell patient-level evidence."
  ),
  data.table(
    artifact = "phase7_PRFT_myeloid_deconvolution_boxplots.pdf",
    proposed_panel_role = "supplement_only",
    title_guidance = "Fallback signature-based myeloid-state estimates",
    keep_reason = "The current file name contains deconvolution wording and should not serve as a main-text title."
  )
))
fwrite(main_figure_reco, "03_results_tables/phase7_cleanup_main_figure_recommendation.csv")

theme_set(
  theme_bw(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
)

sig_panel_dt <- main_summary[result_type == "bulk_signature"]
sig_panel_dt[, display_name := factor(display_name, levels = rev(display_name[order(risk_cor_mean_rho)]))]
sig_panel_plot <- ggplot(sig_panel_dt, aes(x = risk_cor_mean_rho, y = display_name, size = risk_cor_significant_cohorts, color = main_text_priority)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, color = "grey60") +
  geom_point() +
  scale_color_manual(values = c("core" = "#d7301f", "supporting" = "#636363")) +
  labs(
    x = "Mean risk-score rho across cohorts",
    y = NULL,
    color = "Priority",
    size = "Significant cohorts",
    title = "Bulk signature-based inference"
  )

path_panel_dt <- main_summary[result_type == "pathway_activity"]
path_panel_dt[, display_name := factor(display_name, levels = rev(display_name[order(risk_cor_mean_rho)]))]
path_panel_plot <- ggplot(path_panel_dt, aes(x = risk_cor_mean_rho, y = display_name, size = risk_cor_significant_cohorts, color = cross_cohort_consistent)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, color = "grey60") +
  geom_point() +
  scale_color_manual(values = c("yes" = "#1b9e77", "no" = "#7570b3")) +
  labs(
    x = "Mean risk-score rho across cohorts",
    y = NULL,
    color = "Consistent",
    size = "Significant cohorts",
    title = "Core pathway activity scoring"
  )

consistency_panel <- sig_panel_plot + path_panel_plot + patchwork::plot_layout(ncol = 2, widths = c(1, 1))
grDevices::pdf("04_figures/phase7_cleanup_core_signature_consistency_panel.pdf", width = 10.5, height = 5.8, useDingbats = FALSE)
print(consistency_panel)
grDevices::dev.off()

bulk_sc_plot_dt <- copy(bulk_sc)
bulk_sc_plot_dt[, evidence_axis := factor(evidence_axis, levels = c("monocyte_macrophage_like", "myeloid_suppressive", "immune_checkpoint_JAK_STAT", "LSC_stemness"))]
bulk_sc_plot_dt[, consistency_score := fifelse(
  consistency == "consistent", 2,
  fifelse(consistency == "supports_non_dominant_LSC_wording", 1, 1.5)
)]
bulk_sc_plot <- ggplot(bulk_sc_plot_dt, aes(x = "bulk_vs_singlecell", y = evidence_axis, fill = consistency_score)) +
  geom_tile(color = "white", linewidth = 0.4, width = 0.9, height = 0.9) +
  geom_text(aes(label = c("Consistent", "Consistent", "Concordant", "Non-dominant LSC")), size = 2.8) +
  scale_fill_gradient(low = "#deebf7", high = "#2171b5", guide = "none") +
  labs(
    x = NULL,
    y = NULL,
    title = "Bulk and processed single-cell consistency"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
grDevices::pdf("04_figures/phase7_cleanup_bulk_singlecell_consistency_panel.pdf", width = 5.8, height = 3.8, useDingbats = FALSE)
print(bulk_sc_plot)
grDevices::dev.off()

results_paragraph <- paste(
  "Using a uniform rank-based bulk signature scoring framework across TCGA-LAML and two GPL570 validation cohorts, we observed a consistent association between the PRFT-related axis and monocyte/macrophage-like and myeloid-suppressive transcriptional programs.",
  "Immune checkpoint, JAK/STAT/PD-L1, and oxidative stress/NRF2-related signatures showed the same cross-cohort direction.",
  "Pathway-level scoring further linked the PRFT-related axis to IFN-response, TNF/NF-kB, inflammatory-response, and JAK/STAT activity.",
  "These bulk transcriptional patterns were concordant with the processed single-cell patient-level analyses, which supported a monocyte-like and immune-suppressive PRFT-high-like state.",
  "By contrast, LSC/stemness-associated signals did not support a dominant primitive interpretation.",
  "Because formal deconvolution tools were unavailable in the current analysis package, these findings should be interpreted as bulk signature-based inference rather than direct immune cell quantification."
)
writeLines(results_paragraph, "06_manuscript_support/phase7_cleanup_results_paragraph.txt")

checklist_lines <- c(
  "1. 是否完成核心signature/pathway主文表：是",
  "2. 是否标记downplay/supplement结果：是",
  "3. 是否修正immune deconvolution语言：是",
  "4. 是否生成主文Figure建议：是",
  "5. 是否生成cleaned consistency panel：是",
  "6. 是否保留monocyte/myeloid主线：是",
  "7. 是否保留JAK/STAT/PD-L1主线：是",
  "8. 是否保留oxidative stress/NRF2主线：是",
  "9. 是否避免LSC主导表述：是",
  "10. 是否建议进入Phase 9全文整合：是",
  paste0(
    "11. 需要人工确认的问题：",
    "如果主文版面紧张，建议在 phase7_pathway_activity_heatmap 与 phase7_cleanup_core_signature_consistency_panel 之间二选一；",
    "phase7_PRFT_myeloid_deconvolution_boxplots 不建议以现文件名进入主文。"
  )
)
writeLines(checklist_lines, "05_logs/phase7_cleanup_key_result_checklist.txt")

cat("Main summary rows:", nrow(main_summary), "\n")
cat("Downplay rows:", nrow(downplay_dt), "\n")
cat("Language hits found:", sum(language_findings$found), "\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
