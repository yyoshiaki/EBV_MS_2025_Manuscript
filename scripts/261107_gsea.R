# Production source: EBV_scRNAseq/scripts/260327_nebula_reactome_go_gsea.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Hs.eg.db)
})
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 4)
input_dir <- args[[1]]
output_root <- args[[2]]
config_path <- args[[3]]
selected_configs <- strsplit(args[[4]], ",")[[1]]
table_dir <- file.path(output_root, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
th_mean <- 0.5
config_df <- read.delim(config_path, check.names = FALSE)
config_df <- subset(config_df, !grepl("^Bcell_", name))
config_df <- subset(config_df, name %in% selected_configs)
if (nrow(config_df) == 0) {
  stop("No selected configs matched the nebula config table.")
}

reactome_results <- list()
go_results <- list()

set.seed(1)

for (i in seq_len(nrow(config_df))) {
  cfg_name <- config_df$name[[i]]
  deg_path <- file.path(table_dir, paste0("deg_", cfg_name, ".tsv"))
  gene_path <- file.path(input_dir, paste0(cfg_name, "_genes.tsv"))

  if (!file.exists(deg_path) || !file.exists(gene_path)) {
    stop(paste("Missing matched DEG or gene metadata input for", cfg_name))
  }

  deg_df <- read.delim(deg_path, check.names = FALSE)
  deg_df$gene_symbol <- if ("gene" %in% names(deg_df)) deg_df$gene else deg_df$names

  gene_df <- read.delim(gene_path, check.names = FALSE)
  deg_df <- left_join(deg_df, gene_df, by = c("gene_symbol" = "gene"), suffix = c("", "_input"))
  if ("mean_input" %in% names(deg_df)) {
    deg_df$mean_post <- deg_df$mean_input
  } else if ("mean" %in% names(deg_df)) {
    deg_df$mean_post <- deg_df$mean
  } else {
    deg_df$mean_post <- NA_real_
  }

  rank_df <- deg_df %>%
    filter(!is.na(gene_symbol), gene_symbol != "", !is.na(logfoldchanges), mean_post > th_mean) %>%
    distinct(gene_symbol, .keep_all = TRUE)

  if (nrow(rank_df) == 0) next

  entrez_df <- bitr(
    unique(rank_df$gene_symbol),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )

  rank_df <- inner_join(rank_df, entrez_df, by = c("gene_symbol" = "SYMBOL")) %>%
    group_by(ENTREZID) %>%
    arrange(desc(abs(logfoldchanges)), .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    arrange(desc(logfoldchanges))

  gene_list <- rank_df$logfoldchanges
  names(gene_list) <- rank_df$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)

  reactome_df <- tibble()
  reactome_obj <- gsePathway(
      geneList = gene_list,
      organism = "human",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      verbose = FALSE,
      seed = TRUE
  )
  if (!is.null(reactome_obj)) {
    reactome_df <- as.data.frame(reactome_obj) %>%
      filter(p.adjust < 0.05) %>%
      mutate(Cluster = cfg_name)
  }

  go_df <- tibble()
  go_obj <- gseGO(
      geneList = gene_list,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      verbose = FALSE,
      seed = TRUE
  )
  if (!is.null(go_obj)) {
    go_df <- as.data.frame(go_obj) %>%
      filter(p.adjust < 0.05) %>%
      mutate(Cluster = cfg_name)
  }

  reactome_results[[cfg_name]] <- reactome_df
  go_results[[cfg_name]] <- go_df

}

write.table(bind_rows(reactome_results), file.path(table_dir, "reactome_gsea.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bind_rows(go_results), file.path(table_dir, "go_bp_gsea.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
