# Production source: EBV_scRNAseq/scripts/260326_run_scoper_define_clones.R

suppressPackageStartupMessages({
  library(data.table)
  library(scoper)
  library(shazam)
})

# Inputs are donor-split productive IGH AIRR files and the eligible donor list.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
donor_manifest <- fread(args[[1]])
ebv_positive_donors <- fread(args[[2]])$donor
out_dir <- args[[3]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fallback_threshold <- 0.01
nproc <- 1
threshold_rows <- list()

for (donor_id in ebv_positive_donors) {
  donor_files <- donor_manifest[donor == donor_id][["output_airr"]]
  donor_dt <- rbindlist(lapply(donor_files, fread), use.names = TRUE, fill = TRUE)
  if ("locus" %in% colnames(donor_dt)) {
    donor_dt[, locus := as.character(locus)]
  } else {
    donor_dt[, locus := "IGH"]
  }
  if ("clone_id" %in% colnames(donor_dt)) {
    setnames(donor_dt, "clone_id", "source_clone_id")
  }
  setorder(donor_dt, lineage_obs_name, -consensus_count)
  donor_cell_dt <- unique(donor_dt, by = "lineage_obs_name")

  threshold_source <- "fixed_fallback"
  dist_dt <- distToNearest(
    donor_cell_dt,
    cellIdColumn = "lineage_obs_name",
    locusColumn = "locus",
    onlyHeavy = TRUE,
    nproc = nproc
  )
  dist_values <- dist_dt[["dist_nearest"]]
  threshold_value <- fallback_threshold
  threshold_fit <- NULL
  gmm_fit <- NULL
  finite_dist <- dist_values[is.finite(dist_values)]
  positive_dist <- finite_dist[finite_dist > 0]
  fit_dist <- finite_dist
  if (length(fit_dist) > 10000) {
    set.seed(1)
    fit_dist <- sample(fit_dist, 10000)
  }
  density_threshold <- NA_real_
  gmm_threshold <- NA_real_
  hist_threshold <- NA_real_
  min_auto_threshold <- fallback_threshold
  if (length(positive_dist) >= 20) {
    min_auto_threshold <- max(
      fallback_threshold,
      as.numeric(stats::quantile(positive_dist, probs = 0.05, na.rm = TRUE, names = FALSE))
    )
  }
  max_auto_threshold <- if (length(positive_dist) >= 20) {
    as.numeric(stats::quantile(positive_dist, probs = 0.98, na.rm = TRUE, names = FALSE))
  } else {
    Inf
  }

  if (length(fit_dist) >= 20 && length(unique(fit_dist)) >= 5) {
    threshold_fit <- tryCatch(findThreshold(fit_dist, method = "density"), error = function(e) NULL)
    if (!is.null(threshold_fit) && is.finite(threshold_fit@threshold)) {
      density_threshold <- threshold_fit@threshold
    }
    if (is.finite(density_threshold) &&
        density_threshold >= min_auto_threshold &&
        density_threshold <= max_auto_threshold) {
      threshold_value <- density_threshold
        threshold_source <- "density"
    }

    if (threshold_source == "fixed_fallback") {
      gmm_fit <- tryCatch(
        findThreshold(fit_dist, method = "gmm", cutoff = "optimal"),
        error = function(e) NULL
      )
      if (!is.null(gmm_fit) && is.finite(gmm_fit@threshold)) {
        gmm_threshold <- gmm_fit@threshold
      }
      if (is.finite(gmm_threshold) &&
          gmm_threshold >= min_auto_threshold &&
          gmm_threshold <= max_auto_threshold) {
        threshold_value <- gmm_threshold
        threshold_source <- "gmm"
      }
    }

    if (threshold_source == "fixed_fallback") {
      hist_obj <- hist(fit_dist, breaks = "FD", plot = FALSE)
      if (length(hist_obj$counts) >= 5) {
        valley_idx <- which(diff(sign(diff(hist_obj$counts))) == 2) + 1
        peak_idx <- which.max(hist_obj$counts)
        valley_idx <- valley_idx[valley_idx > peak_idx]
        valley_mid <- hist_obj$mids[valley_idx]
        valley_mid <- valley_mid[valley_mid >= min_auto_threshold & valley_mid <= max_auto_threshold]
        if (length(valley_mid) > 0) {
          hist_threshold <- valley_mid[[1]]
          threshold_value <- hist_threshold
          threshold_source <- "histogram_valley"
        }
      }
    }
  }

  clone_result <- hierarchicalClones(
    as.data.frame(donor_cell_dt),
    threshold = threshold_value,
    cell_id = NULL,
    locus = "locus",
    only_heavy = TRUE,
    split_light = FALSE,
    nproc = nproc
  )
  clone_dt <- as.data.table(clone_result@db)
  clone_dt[, scoper_clone_id_raw := as.character(clone_id)]
  clone_dt[, scoper_clone_id := sprintf("%s:%s", donor_id, scoper_clone_id_raw)]

  fwrite(clone_dt, file.path(out_dir, paste0(donor_id, "_clone_assignments.tsv")), sep = "\t")
  threshold_rows[[donor_id]] <- data.table(
    donor = donor_id,
    n_cells = donor_cell_dt[, uniqueN(lineage_obs_name)],
    threshold = threshold_value,
    threshold_source = threshold_source
  )
}
fwrite(rbindlist(threshold_rows), file.path(out_dir, "clone_thresholds.tsv"), sep = "\t")
