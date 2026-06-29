pretty_label <- function(x) {
  x <- as.character(x)
  label_map <- c(
    "PRFT_score" = "PRFT score",
    "risk_score" = "risk score",
    "high_risk" = "high risk",
    "low_risk" = "low risk",
    "PRFT_high" = "PRFT high",
    "PRFT_low" = "PRFT low",
    "PRFT_high_enriched" = "PRFT-high enriched",
    "PRFT_low_enriched" = "PRFT-low enriched",
    "Proteostasis_core" = "proteostasis signature",
    "Proteostasis_core_score" = "proteostasis signature score",
    "z_Proteostasis_core" = "z-scored proteostasis signature",
    "Ferroptosis_tolerance_set" = "ferroptosis-tolerance signature",
    "Ferroptosis_tolerance_set_score" = "ferroptosis-tolerance signature score",
    "z_Ferroptosis_tolerance_set" = "z-scored ferroptosis-tolerance signature",
    "SLC7A11_GPX4_GSH_axis" = "SLC7A11/GPX4-GSH axis",
    "SLC7A11_GPX4_GSH_axis_score" = "SLC7A11/GPX4-GSH axis score",
    "Myeloid_suppressive_set" = "myeloid suppressive signature",
    "Myeloid_suppressive_set_score" = "myeloid suppressive signature score",
    "Myeloid_suppressive_extended" = "myeloid suppressive signature",
    "Immune_checkpoint_set" = "immune checkpoint signature",
    "Immune_checkpoint_set_score" = "immune checkpoint signature score",
    "Immune_checkpoint_extended" = "immune checkpoint signature",
    "PD1_PDL1_axis" = "PD-1/PD-L1 axis",
    "JAK2_STAT5_PDL1_set" = "JAK2/STAT5/PD-L1 signature",
    "JAK2_STAT5_PDL1_set_score" = "JAK2/STAT5/PD-L1 signature score",
    "T_cell_exhaustion_set" = "T-cell exhaustion signature",
    "T_cell_exhaustion_set_score" = "T-cell exhaustion signature score",
    "T_cell_exhaustion_extended" = "T-cell exhaustion signature",
    "Monocyte_macrophage_like" = "monocyte/macrophage-like signature",
    "M2_macrophage_like" = "M2 macrophage-like signature",
    "Neutrophil_inflammatory_like" = "neutrophil inflammatory-like signature",
    "IFN_gamma_response" = "IFN-gamma response",
    "Cytotoxic_T_NK" = "cytotoxic T/NK signature",
    "Antigen_presentation" = "antigen presentation",
    "SUMOylation_set" = "SUMOylation signature",
    "SUMOylation_set_score" = "SUMOylation signature score",
    "NEDDylation_set" = "NEDDylation signature",
    "NEDDylation_set_score" = "NEDDylation signature score",
    "Ferroptosis_driver_set" = "ferroptosis driver signature",
    "Ferroptosis_driver_set_score" = "ferroptosis driver signature score",
    "Stemness_quiescence_set" = "stemness/quiescence signature",
    "Stemness_quiescence_set_score" = "stemness/quiescence signature score",
    "Relapse_resistance_set" = "relapse-resistance signature",
    "Relapse_resistance_set_score" = "relapse-resistance signature score",
    "LSC17_core" = "LSC17 signature",
    "LSC17_core_score" = "LSC17 signature score",
    "Upregulated_in_PRFT_high" = "Upregulated in PRFT high",
    "Downregulated_in_PRFT_high" = "Downregulated in PRFT high",
    "Not_significant" = "Not significant",
    "candidate_pool_all" = "candidate pool: all",
    "candidate_pool_strict" = "candidate pool: strict",
    "main_candidate_genes" = "main candidate genes",
    "lasso_input_genes" = "LASSO input genes",
    "PRFT_low" = "PRFT low",
    "PRFT_high" = "PRFT high",
    "AML_standard_related" = "AML-related therapy",
    "JAK_STAT_related" = "JAK/STAT-related therapy",
    "BCL2_apoptosis_related" = "BCL2/apoptosis-related therapy",
    "Proteostasis_stress_related" = "proteostasis/stress-related therapy",
    "PI3K_AKT_mTOR_related" = "PI3K/AKT/mTOR-related therapy",
    "Oxidative_stress_ferroptosis_adjacent" = "oxidative stress/ferroptosis-adjacent therapy",
    "combined_GPL570" = "combined GPL570",
    "GSE37642_GPL570" = "GSE37642 GPL570",
    "GSE12417_GPL570" = "GSE12417 GPL570",
    "high_risk_enriched" = "high-risk enriched",
    "low_risk_enriched" = "low-risk enriched",
    "high_risk_higher" = "higher in high risk",
    "low_risk_higher" = "higher in low risk",
    "no_difference" = "no difference",
    "FDR_significant" = "FDR significant",
    "exploratory_rawP" = "exploratory raw P",
    "not_selected" = "not selected"
  )
  out <- ifelse(x %in% names(label_map), label_map[x], x)
  out <- gsub("_", " ", out, fixed = TRUE)
  out
}

pretty_factor <- function(x, levels = NULL) {
  x_chr <- as.character(x)
  if (is.null(levels)) {
    levels <- unique(x_chr)
  }
  factor(pretty_label(x_chr), levels = pretty_label(levels))
}

pretty_labeller <- function(variable, value) {
  pretty_label(value)
}
