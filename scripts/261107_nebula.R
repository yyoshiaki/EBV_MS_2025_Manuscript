# Production source: EBV_scRNAseq/scripts/260903_run_pbmc_ebv_b_nebula_corrected_offset.R

suppressPackageStartupMessages({library(Matrix); library(nebula)})
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3)
input_dir <- args[[1]]
output_root <- args[[2]]
config_path <- args[[3]]
fit_dir <- file.path(output_root, "nebula", "fits")
table_dir <- file.path(output_root, "tables")
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
selected_configs <- if (length(args) >= 4 && nzchar(args[[4]])) strsplit(args[[4]], ",")[[1]] else NULL
min_cells_per_gene_runtime <- if (length(args) >= 5) as.integer(args[[5]]) else 50L
ncore_runtime <- if (length(args) >= 6) as.integer(args[[6]]) else 4L
config_df <- read.delim(config_path, check.names = FALSE)
if (!is.null(selected_configs)) {
  config_df <- subset(config_df, name %in% selected_configs)
}
if (nrow(config_df) == 0) stop("No matching config selected.")

th_mean <- 0.5
th_lfc <- 0.25
th_padj <- 0.05

for (cfg_name in config_df$name) {
  cell_metadata <- read.delim(file.path(input_dir, paste0(cfg_name, "_cell_metadata.tsv")), check.names = FALSE)
  gene_table <- read.delim(file.path(input_dir, paste0(cfg_name, "_genes.tsv")), check.names = FALSE)
  count_matrix <- readMM(file.path(input_dir, paste0(cfg_name, "_counts.mtx")))
  count_matrix <- as(count_matrix, "dgCMatrix")
  rownames(count_matrix) <- gene_table$gene
  colnames(count_matrix) <- cell_metadata$cell_key

  cell_metadata$EBV_detected <- factor(cell_metadata$EBV_detected, levels = c("False", "True"))
  cell_metadata$sample <- factor(cell_metadata$sample)
  n_true <- sum(cell_metadata$EBV_detected == "True")
  n_false <- sum(cell_metadata$EBV_detected == "False")
  if (n_true < 5 || n_false < 5) {
    next
  }
  keep_gene_idx <- which(gene_table$n_cells_by_counts >= min_cells_per_gene_runtime)
  gene_table <- gene_table[keep_gene_idx, , drop = FALSE]
  count_matrix <- count_matrix[keep_gene_idx, , drop = FALSE]
  id <- factor(cell_metadata$donor_demux)
  pred <- model.matrix(~ EBV_detected + sample, data = cell_metadata)
  offset <- cell_metadata$library_size_counts

  fit <- nebula(
    count = count_matrix,
    id = id,
    pred = pred,
    offset = offset,
    method = "LN",
    model = "NBGMM",
    cpc = 0.005,
    ncore = ncore_runtime
  )

  saveRDS(fit, file.path(fit_dir, paste0(cfg_name, "_nebula_fit.rds")))

  summary_df <- as.data.frame(fit$summary)
  summary_df$names <- rownames(summary_df)
  merge_key <- if ("gene" %in% names(summary_df)) "gene" else "names"
  effect_col <- grep("^logFC_.*EBV_detected.*True$", names(summary_df), value = TRUE)
  if (length(effect_col) == 0) {
    effect_col <- grep("^logFC_.*EBV_detected", names(summary_df), value = TRUE)
  }
  if (length(effect_col) != 1) {
    stop(paste("Could not identify a unique EBV effect column for", cfg_name))
  }
  se_col <- sub("^logFC_", "se_", effect_col)
  p_col <- sub("^logFC_", "p_", effect_col)
  if (!se_col %in% names(summary_df) || !p_col %in% names(summary_df)) {
    stop(paste("Missing se/p columns for", cfg_name))
  }

  result_df <- merge(summary_df, gene_table, by.x = merge_key, by.y = "gene", all.x = TRUE, sort = FALSE)
  if (!"gene" %in% names(result_df) && merge_key %in% names(result_df)) {
    result_df$gene <- result_df[[merge_key]]
  }
  result_df$logfoldchanges <- result_df[[effect_col]]
  result_df$se <- result_df[[se_col]]
  result_df$pvals <- result_df[[p_col]]
  result_df$pvals_adj <- p.adjust(result_df$pvals, method = "BH")
  result_df$strict_deg <- with(result_df, pvals_adj < th_padj & abs(logfoldchanges) > th_lfc & mean > th_mean)
  result_df <- result_df[order(result_df$pvals, -abs(result_df$logfoldchanges)), ]

  write.table(result_df, file.path(table_dir, paste0("deg_", cfg_name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

}
