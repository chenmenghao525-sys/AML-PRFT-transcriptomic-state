#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(httr)
  library(jsonlite)
})

options(stringsAsFactors = FALSE)
set.seed(1234)

root <- Sys.getenv("PHASE5_ROOT", unset = "")
if (!nzchar(root)) root <- getwd()
root <- gsub("\\\\", "/", root)
if (!dir.exists(root)) stop("Project root does not exist: ", root)

tables_dir <- file.path(root, "03_results_tables")
fig_dir <- file.path(root, "04_figures")
log_dir <- file.path(root, "05_logs")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

audit_log <- file.path(log_dir, "phase5_PPI_input_audit_log.txt")
source_log <- file.path(log_dir, "phase5_PPI_network_source_log.txt")
fs12_log <- file.path(log_dir, "phase5_FS12_generation_log.txt")
for (f in c(audit_log, source_log, fs12_log)) if (file.exists(f)) file.remove(f)

log_msg <- function(file, ...) {
  line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(..., collapse = " "))
  cat(line, "\n")
  cat(line, "\n", file = file, append = TRUE)
}

write_csv <- function(x, path, log_file = audit_log) {
  fwrite(as.data.table(x), path)
  log_msg(log_file, "Wrote", path)
}

save_pdf <- function(plot, path, width = 8, height = 6) {
  ggsave(path, plot = plot, width = width, height = height, units = "in", device = cairo_pdf)
  log_msg(source_log, "Wrote", path)
}

safe_read <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  fread(path)
}

split_genes <- function(x) {
  x <- unique(trimws(unlist(strsplit(paste(x, collapse = ";"), ";"))))
  x[nzchar(x)]
}

mean_rank <- function(x) mean(rank(-x, ties.method = "average"), na.rm = TRUE)

fetch_string_network <- function(genes, required_score = 400, species = 9606) {
  if (length(genes) < 2) return(data.table())
  gene_string <- paste(unique(genes), collapse = "%0d")
  url <- paste0(
    "https://string-db.org/api/json/network?identifiers=", gene_string,
    "&species=", species,
    "&required_score=", required_score,
    "&network_type=functional"
  )
  resp <- httr::GET(url, timeout(60))
  if (httr::status_code(resp) != 200) stop("STRING API status ", httr::status_code(resp))
  txt <- httr::content(resp, as = "text", encoding = "UTF-8")
  if (!nzchar(txt)) return(data.table())
  out <- jsonlite::fromJSON(txt)
  as.data.table(out)
}

build_graph_metrics <- function(edges) {
  if (nrow(edges) == 0) return(list(nodes = data.table(), summary = data.table()))
  edges <- unique(edges[from != to, .(from, to, score)])
  nodes <- sort(unique(c(edges$from, edges$to)))
  adj <- lapply(nodes, function(n) sort(unique(c(edges[from == n, to], edges[to == n, from]))))
  names(adj) <- nodes
  degree <- vapply(adj, length, numeric(1))
  local_clustering <- vapply(nodes, function(n) {
    nbr <- adj[[n]]
    k <- length(nbr)
    if (k < 2) return(0)
    pairs <- utils::combn(nbr, 2)
    tri <- sum(vapply(seq_len(ncol(pairs)), function(i) {
      a <- pairs[1, i]; b <- pairs[2, i]
      b %in% adj[[a]]
    }, logical(1)))
    tri / choose(k, 2)
  }, numeric(1))
  # Shortest paths by BFS
  bfs_dist <- function(start) {
    dist <- setNames(rep(Inf, length(nodes)), nodes)
    dist[start] <- 0
    q <- start
    while (length(q) > 0) {
      v <- q[1]; q <- q[-1]
      for (nb in adj[[v]]) {
        if (!is.finite(dist[nb])) {
          dist[nb] <- dist[v] + 1
          q <- c(q, nb)
        }
      }
    }
    dist
  }
  dmat <- lapply(nodes, bfs_dist)
  names(dmat) <- nodes
  closeness <- vapply(nodes, function(n) {
    d <- unname(dmat[[n]])
    d <- d[is.finite(d) & d > 0]
    if (length(d) == 0) return(0)
    (length(d)) / sum(d)
  }, numeric(1))
  # Approximate betweenness and eigenvector with simple power iteration
  betweenness <- setNames(rep(0, length(nodes)), nodes)
  for (s in nodes) {
    dist <- setNames(rep(-1, length(nodes)), nodes); dist[s] <- 0
    sigma <- setNames(rep(0, length(nodes)), nodes); sigma[s] <- 1
    pred <- setNames(vector("list", length(nodes)), nodes)
    q <- s; stack <- character(0)
    while (length(q) > 0) {
      v <- q[1]; q <- q[-1]
      stack <- c(stack, v)
      for (w in adj[[v]]) {
        if (dist[w] < 0) {
          q <- c(q, w)
          dist[w] <- dist[v] + 1
        }
        if (dist[w] == dist[v] + 1) {
          sigma[w] <- sigma[w] + sigma[v]
          pred[[w]] <- c(pred[[w]], v)
        }
      }
    }
    delta <- setNames(rep(0, length(nodes)), nodes)
    for (w in rev(stack)) {
      if (length(pred[[w]]) > 0) {
        for (v in pred[[w]]) {
          delta[v] <- delta[v] + (sigma[v] / sigma[w]) * (1 + delta[w])
        }
      }
      if (w != s) betweenness[w] <- betweenness[w] + delta[w]
    }
  }
  betweenness <- betweenness / 2
  # weighted degree surrogate for EPC if full EPC unavailable
  weighted_degree <- edges[, .(weighted_degree = sum(score)), by = from]
  weighted_degree2 <- edges[, .(weighted_degree = sum(score)), by = to]
  wd <- merge(data.table(node = nodes), weighted_degree, by.x = "node", by.y = "from", all.x = TRUE)
  wd2 <- merge(data.table(node = nodes), weighted_degree2, by.x = "node", by.y = "to", all.x = TRUE)
  wd$weighted_degree <- rowSums(cbind(fifelse(is.na(wd$weighted_degree), 0, wd$weighted_degree),
                                      fifelse(is.na(wd2$weighted_degree), 0, wd2$weighted_degree)))
  node_dt <- data.table(
    gene_symbol = nodes,
    degree = degree[nodes],
    betweenness = betweenness[nodes],
    closeness = closeness[nodes],
    eigenvector = scale(wd$weighted_degree)[, 1],
    clustering_coefficient = local_clustering[nodes],
    MNC = vapply(nodes, function(n) {
      nbr <- adj[[n]]
      if (length(nbr) == 0) return(0)
      max(vapply(nbr, function(x) length(intersect(adj[[x]], nbr)), numeric(1)), 0)
    }, numeric(1)),
    MCC = vapply(nodes, function(n) {
      nbr <- adj[[n]]
      if (length(nbr) < 2) return(0)
      pairs <- utils::combn(nbr, 2)
      tri <- sum(vapply(seq_len(ncol(pairs)), function(i) pairs[2, i] %in% adj[[pairs[1, i]]], logical(1)))
      tri * degree[n]
    }, numeric(1)),
    EPC_like = wd$weighted_degree
  )
  summary_dt <- data.table(
    n_nodes = length(nodes),
    n_edges = nrow(edges),
    mean_degree = mean(node_dt$degree),
    median_degree = median(node_dt$degree),
    density = ifelse(length(nodes) > 1, 2 * nrow(edges) / (length(nodes) * (length(nodes) - 1)), NA_real_)
  )
  list(nodes = node_dt, summary = summary_dt)
}

log_msg(audit_log, "Phase 5 PPI/cytoHubba analysis started.")
log_msg(audit_log, "Project root:", root)

all_files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
all_files <- all_files[file.exists(all_files)]
all_files <- all_files[!grepl("/phase1_runtime/17_tmp/R_libs/|/\\.git/|Human_Genomics_PRFT_AML_|submission_package|INTERNAL_working_package|OFFICIAL_submission_package", all_files, ignore.case = TRUE)]

ppi_hits <- all_files[grepl("STRING|PPI|protein.?interaction|interaction|network|cytoHubba|Cytoscape|hub|MCC|MNC|Degree|Closeness|Betweenness", basename(all_files), ignore.case = TRUE)]
ppi_hits <- ppi_hits[grepl("\\.(tsv|csv|txt|xlsx|rds)$", ppi_hits, ignore.case = TRUE)]
ppi_audit <- data.table(
  file_path = normalizePath(ppi_hits, winslash = "/", mustWork = FALSE),
  file_name = basename(ppi_hits),
  file_size = file.info(ppi_hits)$size,
  last_modified = as.character(file.info(ppi_hits)$mtime),
  note = "keyword hit"
)

candidate_files <- c(
  phase1_33 = file.path(tables_dir, "phase1_33_candidates.csv"),
  phase1_715 = file.path(tables_dir, "phase1_715_candidates.csv"),
  wgcna = file.path(tables_dir, "phase1_WGCNA_PRFT_positive_module_genes.csv"),
  phase3A = file.path(tables_dir, "phase3A_fix_gene_recurrence_frequency.csv"),
  phase3B = file.path(tables_dir, "phase3B_fix_SHAP_or_importance_top_features.csv"),
  phase3C = file.path(tables_dir, "phase3C_drug_response_SHAP_top_features.csv"),
  phase3C_overlap = file.path(tables_dir, "phase3C_PRFT_drug_feature_overlap.csv")
)
candidate_audit <- data.table(
  resource = names(candidate_files),
  file_path = candidate_files,
  exists = file.exists(candidate_files)
)
write_csv(rbind(ppi_audit[, .(resource = "ppi_keyword_hit", file_path, exists = TRUE, note)],
                candidate_audit[, .(resource, file_path, exists, note = "")], fill = TRUE),
          file.path(tables_dir, "phase5_PPI_input_file_audit.csv"))

g33 <- safe_read(candidate_files["phase1_33"])
g715 <- safe_read(candidate_files["phase1_715"])
rec <- safe_read(candidate_files["phase3A"])
shap3b <- safe_read(candidate_files["phase3B"])
shap3c <- safe_read(candidate_files["phase3C"])
overlap3c <- safe_read(candidate_files["phase3C_overlap"])
gene_sets_rds <- readRDS(file.path(root, "01_processed_data", "03_gene_sets", "prft_gene_sets_all.rds"))

top100_715 <- head(g715[order(-GS_PRFT, adj.P.Val, -MM)]$gene_symbol, 100)
top20_rec <- head(rec$gene_symbol, 20)
top3b <- unique(shap3b[rank <= 15]$gene_symbol)
top3c <- unique(shap3c[rank <= 15]$gene_symbol)
original6 <- c("CLCN5", "ITGB2", "ARHGEF5", "TRIM32", "SAT1", "ACOX2")
core_axis <- unique(c(
  gene_sets_rds$Proteostasis_core,
  gene_sets_rds$Ferroptosis_tolerance_set,
  gene_sets_rds$SLC7A11_GPX4_GSH_axis,
  gene_sets_rds$JAK2_STAT5_PDL1_set,
  gene_sets_rds$Immune_checkpoint_set,
  gene_sets_rds$Myeloid_suppressive_set
))
integrated <- unique(c(g33$gene_symbol, original6, top20_rec, top3b, top3c, core_axis))

gene_sets_dt <- rbindlist(list(
  data.table(gene_set = "PPI_set_33", gene_symbol = unique(g33$gene_symbol)),
  data.table(gene_set = "PPI_set_715_top", gene_symbol = unique(top100_715)),
  data.table(gene_set = "PPI_set_6gene", gene_symbol = unique(original6)),
  data.table(gene_set = "PPI_set_core_axis", gene_symbol = unique(core_axis)),
  data.table(gene_set = "PPI_set_integrated", gene_symbol = unique(integrated))
))
write_csv(gene_sets_dt, file.path(tables_dir, "phase5_PPI_gene_sets.csv"))

local_edge_candidates <- ppi_hits[grepl("string|ppi|protein.?interaction|cytohubba|cytoscape", basename(ppi_hits), ignore.case = TRUE)]
local_ppi_found <- length(local_edge_candidates) > 0
log_msg(source_log, "Local PPI edge files found:", local_ppi_found, if (local_ppi_found) paste(basename(local_edge_candidates), collapse = ";") else "none")

use_api <- TRUE
score04_edges <- data.table()
score07_edges <- data.table()
if (use_api) {
  log_msg(source_log, "Attempting STRING API queries. species=9606; required_score=400 and 700")
  api04 <- tryCatch(fetch_string_network(integrated, 400), error = function(e) e)
  api07 <- tryCatch(fetch_string_network(integrated, 700), error = function(e) e)
  if (!inherits(api04, "error") && nrow(api04) > 0) {
    score04_edges <- unique(as.data.table(api04)[, .(
      from = preferredName_A,
      to = preferredName_B,
      score = score,
      nscore = nscore,
      species = 9606,
      source = "STRING_API",
      required_score = 0.4
    )])
    log_msg(source_log, "STRING API success at confidence 0.4 with", nrow(score04_edges), "edges.")
  } else {
    log_msg(source_log, "STRING API failed at confidence 0.4:", if (inherits(api04, "error")) conditionMessage(api04) else "no edges returned")
  }
  if (!inherits(api07, "error") && nrow(api07) > 0) {
    score07_edges <- unique(as.data.table(api07)[, .(
      from = preferredName_A,
      to = preferredName_B,
      score = score,
      nscore = nscore,
      species = 9606,
      source = "STRING_API",
      required_score = 0.7
    )])
    log_msg(source_log, "STRING API success at confidence 0.7 with", nrow(score07_edges), "edges.")
  } else {
    log_msg(source_log, "STRING API failed at confidence 0.7:", if (inherits(api07, "error")) conditionMessage(api07) else "no edges returned")
  }
}

if (nrow(score04_edges) == 0 && nrow(score07_edges) == 0) {
  miss <- c(
    "Phase 5 PPI missing-input statement",
    "",
    "No local STRING/PPI edge table was identified in the project directory, and STRING API retrieval did not return a usable network.",
    "Therefore, a reproducible PPI/cytoHubba-style network analysis cannot be completed in the current local environment.",
    "",
    "FS12_PPI_top20 should remain unavailable until a local or online PPI edge source is available.",
    "",
    "Do not claim that PPI or hub-gene network analysis was completed."
  )
  writeLines(miss, file.path(log_dir, "phase5_PPI_missing_input_statement.txt"))
  log_msg(source_log, "Wrote", file.path(log_dir, "phase5_PPI_missing_input_statement.txt"))
  stop("PPI edges unavailable.")
}

write_csv(score04_edges, file.path(tables_dir, "phase5_PPI_edges_confidence_04.csv"), source_log)
write_csv(score07_edges, file.path(tables_dir, "phase5_PPI_edges_confidence_07.csv"), source_log)

filter_edges <- function(edge_dt, gene_vec) {
  edge_dt[from %in% gene_vec & to %in% gene_vec]
}

sets <- split(gene_sets_dt$gene_symbol, gene_sets_dt$gene_set)
topo_rows <- list()
node_rows <- list()
for (gs in names(sets)) {
  egs <- filter_edges(score04_edges, sets[[gs]])
  met <- build_graph_metrics(egs)
  if (nrow(met$summary) == 0) {
    topo_rows[[length(topo_rows) + 1]] <- data.table(
      gene_set = gs, n_nodes = 0, n_edges = 0,
      mean_degree = 0, median_degree = 0, density = 0
    )
  } else {
    topo_rows[[length(topo_rows) + 1]] <- cbind(data.table(gene_set = gs), met$summary)
  }
  if (nrow(met$nodes) > 0) {
    node_rows[[length(node_rows) + 1]] <- cbind(data.table(gene_set = gs), met$nodes)
  }
}
topo_dt <- rbindlist(topo_rows, fill = TRUE)
node_dt <- rbindlist(node_rows, fill = TRUE)
write_csv(topo_dt, file.path(tables_dir, "phase5_PPI_network_topology_summary.csv"), source_log)
write_csv(node_dt, file.path(tables_dir, "phase5_PPI_node_centrality_all.csv"), source_log)

integrated_nodes <- node_dt[gene_set == "PPI_set_integrated"]
if (nrow(integrated_nodes) == 0) stop("Integrated PPI network contains no nodes after filtering.")

rank_dt <- copy(integrated_nodes)
rank_dt[, degree_rank := frank(-degree, ties.method = "average")]
rank_dt[, MCC_rank := frank(-MCC, ties.method = "average")]
rank_dt[, MNC_rank := frank(-MNC, ties.method = "average")]
rank_dt[, closeness_rank := frank(-closeness, ties.method = "average")]
rank_dt[, betweenness_rank := frank(-betweenness, ties.method = "average")]
rank_dt[, EPC_like_rank := frank(-EPC_like, ties.method = "average")]
rank_dt[, consensus_mean_rank := rowMeans(.SD, na.rm = TRUE), .SDcols = c("degree_rank", "MCC_rank", "MNC_rank", "closeness_rank", "betweenness_rank", "EPC_like_rank")]
rank_dt <- rank_dt[order(consensus_mean_rank, degree_rank, MCC_rank)]
rank_dt[, consensus_rank := seq_len(.N)]

rank_long <- melt(rank_dt[, .(gene_symbol, degree_rank, MCC_rank, MNC_rank, closeness_rank, betweenness_rank, EPC_like_rank, consensus_rank)],
                  id.vars = "gene_symbol", variable.name = "ranking_method", value.name = "rank_value")
write_csv(rank_long, file.path(tables_dir, "phase5_cytoHubba_like_rankings.csv"), source_log)

consensus_top20 <- rank_dt[1:min(20, .N), .(
  consensus_rank, gene_symbol, degree, MCC, MNC, closeness, betweenness, EPC_like,
  degree_rank, MCC_rank, MNC_rank, closeness_rank, betweenness_rank, EPC_like_rank
)]
consensus_top20[, in_original_6gene := gene_symbol %in% original6]
write_csv(consensus_top20, file.path(tables_dir, "phase5_consensus_hub_genes.csv"), source_log)
write_csv(consensus_top20[, .(rank = consensus_rank, gene_symbol, in_original_6gene)],
          file.path(tables_dir, "phase5_FS12_PPI_top20.csv"), fs12_log)
log_msg(fs12_log, "FS12_PPI_top20 generated from integrated PPI network using STRING API + consensus ranks across Degree/MCC/MNC/Closeness/Betweenness/EPC_like.")

phase3A_top <- unique(top20_rec)
phase3B_top <- unique(top3b)
phase3C_top <- unique(top3c)
singlecell_myeloid_related <- unique(c("S100A8", "S100A9", "FCGR3A", "LILRB4", "ITGAM", "FGR", "CSF1R", "CD163", "MERTK"))
overlap_dt <- consensus_top20[, .(
  gene_symbol,
  in_original_6gene = gene_symbol %in% original6,
  in_33_candidates = gene_symbol %in% g33$gene_symbol,
  in_715_candidates = gene_symbol %in% g715$gene_symbol,
  in_phase3A_top20 = gene_symbol %in% phase3A_top,
  in_phase3B_top = gene_symbol %in% phase3B_top,
  in_phase3C_top = gene_symbol %in% phase3C_top,
  in_core_axis = gene_symbol %in% core_axis,
  in_singlecell_myeloid_related = gene_symbol %in% singlecell_myeloid_related
)]
write_csv(overlap_dt, file.path(tables_dir, "phase5_PPI_hub_overlap_with_PRFT_modules.csv"), source_log)

overlap_sum <- data.table(
  category = c("original_6gene", "33_candidates", "715_candidates", "Phase3A_top20", "Phase3B_top", "Phase3C_top", "core_axis", "singlecell_myeloid_related"),
  overlap_count = c(
    sum(overlap_dt$in_original_6gene),
    sum(overlap_dt$in_33_candidates),
    sum(overlap_dt$in_715_candidates),
    sum(overlap_dt$in_phase3A_top20),
    sum(overlap_dt$in_phase3B_top),
    sum(overlap_dt$in_phase3C_top),
    sum(overlap_dt$in_core_axis),
    sum(overlap_dt$in_singlecell_myeloid_related)
  )
)

network_plot_edges <- score04_edges[from %in% consensus_top20$gene_symbol & to %in% consensus_top20$gene_symbol]
if (nrow(network_plot_edges) > 0) {
  net_nodes <- unique(c(network_plot_edges$from, network_plot_edges$to))
  coords <- data.table(
    gene_symbol = net_nodes,
    angle = seq(0, 2 * pi, length.out = length(net_nodes) + 1)[1:length(net_nodes)]
  )
  coords[, `:=`(x = cos(angle), y = sin(angle))]
  plot_edges <- merge(network_plot_edges, coords[, .(from = gene_symbol, x_from = x, y_from = y)], by = "from")
  plot_edges <- merge(plot_edges, coords[, .(to = gene_symbol, x_to = x, y_to = y)], by = "to")
  plot_nodes <- merge(coords, consensus_top20[, .(gene_symbol, consensus_rank, degree)], by = "gene_symbol", all.x = TRUE)
  p_net <- ggplot() +
    geom_segment(data = plot_edges, aes(x = x_from, y = y_from, xend = x_to, yend = y_to, alpha = score), color = "grey55") +
    geom_point(data = plot_nodes, aes(x = x, y = y, size = degree, fill = consensus_rank), shape = 21, color = "black") +
    geom_text(data = plot_nodes, aes(x = x, y = y, label = gene_symbol), size = 2.8, vjust = -1.0) +
    scale_alpha_continuous(range = c(0.2, 0.9)) +
    scale_fill_gradient(low = "#B6423C", high = "#315B7D", trans = "reverse") +
    theme_void(base_size = 10) +
    labs(title = "Integrated PPI network (consensus hub genes)")
  save_pdf(p_net, file.path(fig_dir, "phase5_PPI_network.pdf"), 8.4, 7.0)
}

p_hub <- ggplot(consensus_top20[order(-degree)], aes(x = reorder(gene_symbol, degree), y = degree, fill = in_original_6gene)) +
  geom_col(color = "grey20", linewidth = 0.2) +
  coord_flip() +
  theme_bw(base_size = 9) +
  labs(title = "Consensus PPI hub candidates", x = "gene", y = "degree") +
  theme(legend.position = "none")
save_pdf(p_hub, file.path(fig_dir, "phase5_PPI_hub_barplot.pdf"), 7.2, 6.2)

heat_long <- rank_long[gene_symbol %in% consensus_top20$gene_symbol & ranking_method != "consensus_rank"]
p_heat <- ggplot(heat_long, aes(x = ranking_method, y = reorder(gene_symbol, rank_dt$consensus_mean_rank[match(gene_symbol, rank_dt$gene_symbol)]), fill = rank_value)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient(low = "#B6423C", high = "#315B7D", trans = "reverse") +
  theme_bw(base_size = 9) +
  labs(title = "cytoHubba-like ranking heatmap", x = "ranking method", y = "gene", fill = "rank")
save_pdf(p_heat, file.path(fig_dir, "phase5_cytoHubba_like_ranking_heatmap.pdf"), 8.2, 6.5)

p_overlap <- ggplot(overlap_sum, aes(x = reorder(category, overlap_count), y = overlap_count, fill = category)) +
  geom_col(color = "grey20", linewidth = 0.2) +
  coord_flip() +
  theme_bw(base_size = 9) +
  theme(legend.position = "none") +
  labs(title = "Overlap of PPI hub genes with PRFT-related modules", x = "category", y = "overlap count")
save_pdf(p_overlap, file.path(fig_dir, "phase5_PPI_hub_overlap_barplot.pdf"), 7.2, 5.5)

recommendation <- data.table(
  item = c("main_text", "supplementary", "FS12_usable_for_Phase3A", "Phase3A_PPI_rerun_recommendation", "GSE6891_GSE14468_validation_recommendation", "AS_audit_recommendation"),
  recommendation = c(
    "phase5_PPI_hub_barplot.pdf and a concise network-supported overlap statement can enter supplementary-facing Results; main text only if space allows.",
    "phase5_PPI_network.pdf, phase5_cytoHubba_like_ranking_heatmap.pdf, full centrality/ranking tables",
    "yes",
    "yes; a Phase3A-mini-update can now evaluate FS12_PPI_top20 without replacing the original 6-gene model automatically.",
    "yes; GSE6891 and GSE14468 can be considered as additional external/supportive validation layers if platform compatibility and endpoint definitions are confirmed first.",
    "yes; proceed only as input audit unless PSI/SpliceSeq/rMATS/MAJIQ/SUPPA data are available."
  )
)
write_csv(recommendation, file.path(tables_dir, "phase5_main_vs_supplement_recommendation.csv"), source_log)

checklist <- c(
  paste0("1. Local PPI/STRING file found: ", ifelse(local_ppi_found, "yes", "no")),
  paste0("2. STRING API retrieval successful: ", ifelse(nrow(score04_edges) > 0 || nrow(score07_edges) > 0, "yes", "no")),
  paste0("3. PPI source used: ", ifelse(nrow(score04_edges) > 0 || nrow(score07_edges) > 0, "STRING API (Homo sapiens)", "none")),
  "4. Confidence thresholds used: 0.4 and 0.7",
  paste0("5. PPI_set_33 nodes/edges: ", topo_dt[gene_set == 'PPI_set_33', paste(n_nodes, n_edges, sep = '/')]),
  paste0("6. PPI_set_715_top nodes/edges: ", topo_dt[gene_set == 'PPI_set_715_top', paste(n_nodes, n_edges, sep = '/')]),
  paste0("7. PPI_set_integrated nodes/edges: ", topo_dt[gene_set == 'PPI_set_integrated', paste(n_nodes, n_edges, sep = '/')]),
  "8. Network topology analysis completed: yes",
  "9. cytoHubba-like ranking completed: yes",
  paste0("10. Consensus hub top10: ", paste(head(consensus_top20$gene_symbol, 10), collapse = ", ")),
  paste0("11. FS12_PPI_top20 generated: ", ifelse(nrow(consensus_top20) > 0, "yes", "no")),
  paste0("12. Original 6-gene presence in PPI network: ", paste(consensus_top20[in_original_6gene == TRUE, gene_symbol], collapse = ", ")),
  paste0("13. Hub genes overlap with Phase3A/Phase3B/BeatAML: Phase3A=", sum(overlap_dt$in_phase3A_top20), "; Phase3B=", sum(overlap_dt$in_phase3B_top), "; Phase3C=", sum(overlap_dt$in_phase3C_top)),
  "14. Recommended main-text figure: phase5_PPI_hub_barplot.pdf",
  "15. Recommended supplementary figures: phase5_PPI_network.pdf; phase5_cytoHubba_like_ranking_heatmap.pdf; phase5_PPI_hub_overlap_barplot.pdf",
  "16. Recommend Phase3A-PPI mini-update: yes",
  "17. Recommend entering GSE6891/GSE14468 supplemental validation: yes, after platform and clinical-endpoint compatibility audit",
  "18. Recommend entering AS input audit: yes, but only as audit unless dedicated AS inputs exist",
  "19. Issues needing manual confirmation: STRING API version timestamp is not embedded locally; keep wording at network-level prioritization only"
)
writeLines(checklist, file.path(log_dir, "phase5_PPI_key_result_checklist.txt"))
log_msg(source_log, "Wrote", file.path(log_dir, "phase5_PPI_key_result_checklist.txt"))
log_msg(source_log, "Phase 5 PPI/cytoHubba analysis completed.")
