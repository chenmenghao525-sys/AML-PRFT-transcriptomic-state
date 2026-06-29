#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

options(stringsAsFactors = FALSE)

dir.create("03_gene_sets", recursive = TRUE, showWarnings = FALSE)
dir.create("14_tables", recursive = TRUE, showWarnings = FALSE)
dir.create("16_logs", recursive = TRUE, showWarnings = FALSE)

save_session_info <- function(path) {
  writeLines(capture.output(sessionInfo()), con = path)
}

write_gmt <- function(gene_sets, file) {
  con <- file(file, open = "wt")
  on.exit(close(con), add = TRUE)
  for (nm in names(gene_sets)) {
    genes <- unique(gene_sets[[nm]])
    line <- paste(c(nm, "PRFT_gene_set", genes), collapse = "\t")
    writeLines(line, con = con)
  }
}

prft_gene_sets_main <- list(
  Proteostasis_core = c(
    "HSPA1A","HSPA5","HSP90AA1","HSP90AB1","DNAJB1",
    "CANX","CALR","P4HB",
    "XBP1","ATF4","ERN1","EIF2AK3",
    "SEL1L","SYVN1","DERL1","HERPUD1",
    "PSMB5","PSMC4","UBE2G2","VCP",
    "SQSTM1","LAMP2"
  ),
  Ferroptosis_tolerance_set = c(
    "SLC7A11","SLC3A2","GPX4","GCLC","GCLM","GSS",
    "AIFM2","GCH1",
    "NFE2L2","TXNRD1",
    "FTH1","FTL","SLC40A1",
    "DHODH","CISD1"
  ),
  SLC7A11_GPX4_GSH_axis = c(
    "SLC7A11","SLC3A2","GPX4","GCLC","GCLM",
    "GSS","GSR","GLS","TXNRD1"
  ),
  LSC17_core = c(
    "DNMT3B","ZBTB46","NYNRIN","ARHGAP22","LAPTM4B",
    "MMRN1","DPYSL3","KIAA0125","CDK6","CPXM1",
    "SOCS2","SMIM24","EMP1","NGFRAP1","AKR1C3",
    "ADGRG1","CD34"
  ),
  Stemness_quiescence_set = c(
    "CD34","PROM1","MEIS1","BMI1","SOCS2",
    "GATA2","KLF4","FOXO3","HLF","MECOM"
  ),
  Relapse_resistance_set = c(
    "AKR1C3","LAPTM4B","EMP1","SOCS2","CDK6",
    "BCL2","MCL1","ABCB1","GATA2","CXCR4"
  ),
  Immune_checkpoint_set = c(
    "CD274","PDCD1LG2","VSIR","VTCN1","IDO1","LGALS9"
  ),
  T_cell_exhaustion_set = c(
    "HAVCR2","LAG3","TIGIT","CTLA4","PDCD1","TOX","ENTPD1"
  ),
  Myeloid_suppressive_set = c(
    "S100A8","S100A9","IL10","TGFB1","ARG1","CXCL8","FCGR3A"
  )
)

prft_gene_sets_supplementary <- list(
  SUMOylation_set = c(
    "SUMO1","SUMO2","SUMO3","SAE1","UBA2","UBE2I",
    "PIAS1","PIAS3","PIAS4","RANBP2",
    "SENP1","SENP2","SENP3","SENP5","SENP6","SENP7"
  ),
  NEDDylation_set = c(
    "NEDD8","NAE1","UBA3","UBE2M","UBE2F",
    "RBX1","RBX2","CUL1","CUL2","CUL3","CUL4A","CUL4B","CUL5",
    "DCUN1D1","DCUN1D3"
  ),
  Ferroptosis_driver_set = c(
    "ACSL4","LPCAT3","ALOX5","ALOX12","ALOX15",
    "NCOA4","HMOX1","TFRC","SAT1","POR","PEBP1"
  ),
  JAK2_STAT5_PDL1_set = c(
    "JAK2","STAT5A","STAT5B","CD274","JAK1",
    "STAT1","IRF1","IFNGR1","IFNGR2","SOCS1","SOCS3",
    "BCL2L1","PIM1"
  )
)

prft_gene_sets_all <- c(prft_gene_sets_main, prft_gene_sets_supplementary)

saveRDS(prft_gene_sets_main, "03_gene_sets/prft_gene_sets_main.rds")
saveRDS(prft_gene_sets_supplementary, "03_gene_sets/prft_gene_sets_supplementary.rds")
saveRDS(prft_gene_sets_all, "03_gene_sets/prft_gene_sets_all.rds")
write_gmt(prft_gene_sets_all, "03_gene_sets/prft_gene_sets.gmt")

table_s1 <- rbindlist(lapply(names(prft_gene_sets_all), function(nm) {
  data.table(
    gene_set_name = nm,
    category = if (nm %in% names(prft_gene_sets_main)) "main" else "supplementary",
    gene_count = length(unique(prft_gene_sets_all[[nm]])),
    genes = paste(unique(prft_gene_sets_all[[nm]]), collapse = ";")
  )
}))

fwrite(table_s1, "14_tables/TableS1_prft_gene_sets.csv")
save_session_info("16_logs/sessionInfo_07_prepare_gene_sets.txt")

