# Production source: EBV_scRNAseq/scripts/260326_run_scoper_define_clones.R
# Prepared 2026-10-15 EDT: Explicit donor-split AIRR manifest/eligible-donor file/output root; omit diagnostic histogram drawing; retain threshold selection and per-donor clone calling.
# Statistical settings retained; source hash in docs/261015_provenance.tsv.

suppressPackageStartupMessages({
  library(data.table)
  library(scoper)
  library(shazam)
})

args <- list(
  version = "260326_scoperclone",
  input_manifest = "",
  positive_donors = "",
  output_root = "",
  fallback_threshold = 0.01,
  nproc = 1,
  donors = ""
)
cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) > 0) {
  for (i in seq(1, length(cli_args), by = 2)) {
    key <- gsub("^--", "", cli_args[i])
    if (!key %in% names(args) || i == length(cli_args)) {
      next
    }
    value <- cli_args[i + 1]
    if (key == "fallback_threshold") {
      args[[key]] <- as.numeric(value)
    } else if (key == "nproc") {
      args[[key]] <- as.integer(value)
    } else {
      args[[key]] <- value
    }
  }
}

stopifnot(nzchar(args$output_root), nzchar(args$input_manifest), nzchar(args$positive_donors))
out_dir <- args$output_root
tables_dir <- file.path(out_dir, "tables")
scoper_dir <- file.path(out_dir, "scoper_clones")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(scoper_dir, recursive = TRUE, showWarnings = FALSE)

manifest_path <- args$input_manifest
if (manifest_path == "") {
  manifest_path <- file.path(tables_dir, sprintf("%s_donor_split_airr_manifest.tsv", args$version))
}
donor_manifest <- fread(manifest_path)
donor_summary_path <- args$positive_donors
ebv_positive_donors <- fread(donor_summary_path)$donor
if (args$donors != "") {
  requested_donors <- unlist(strsplit(args$donors, ",", fixed = TRUE))
  requested_donors <- trimws(requested_donors)
  ebv_positive_donors <- intersect(ebv_positive_donors, requested_donors)
}

assignment_rows <- list()
clone_summary_rows <- list()
threshold_rows <- list()

for (donor_id in ebv_positive_donors) {
  message(sprintf("[scoper] start donor=%s", donor_id))
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
    nproc = args$nproc
  )
  dist_values <- dist_dt[["dist_nearest"]]
  threshold_value <- args$fallback_threshold
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
  min_auto_threshold <- args$fallback_threshold
  if (length(positive_dist) >= 20) {
    min_auto_threshold <- max(
      args$fallback_threshold,
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
    nproc = args$nproc
  )
  clone_dt <- as.data.table(clone_result@db)
  clone_dt[, scoper_clone_id_raw := as.character(clone_id)]
  clone_dt[, scoper_clone_id := sprintf("%s:%s", donor_id, scoper_clone_id_raw)]

  clone_size_dt <- clone_dt[
    ,
    .(
      scoper_clone_size = uniqueN(lineage_obs_name),
      n_samples = uniqueN(lineage_source_sample),
      n_ebv_cells = sum(as.character(lineage_EBV_detected) %chin% c("TRUE", "True", "T", "true", "1")),
      cluster_L2_levels = paste(sort(unique(as.character(lineage_cluster_L2))), collapse = ";")
    ),
    by = .(scoper_clone_id_raw, scoper_clone_id)
  ]
  clone_size_dt[, donor := donor_id]

  assignment_dt <- merge(
    clone_dt[
      ,
      .(
        donor = donor_id,
        lineage_obs_name,
        cell_id,
        lineage_source_sample,
        lineage_cluster_L2,
        lineage_EBV_detected,
        scoper_clone_id_raw,
        scoper_clone_id
      )
    ],
    clone_size_dt,
    by = c("donor", "scoper_clone_id_raw", "scoper_clone_id"),
    all.x = TRUE,
    sort = FALSE
  )
  assignment_dt[, scoper_threshold := threshold_value]
  assignment_dt[, scoper_threshold_source := threshold_source]

  clone_summary_dt <- assignment_dt[
    ,
    .(
      scoper_clone_size = unique(scoper_clone_size),
      n_samples = unique(n_samples),
      n_ebv_cells = unique(n_ebv_cells),
      cluster_L2_levels = unique(cluster_L2_levels),
      scoper_threshold = unique(scoper_threshold),
      scoper_threshold_source = unique(scoper_threshold_source)
    ),
    by = .(donor, scoper_clone_id_raw, scoper_clone_id)
  ]
  setorder(clone_summary_dt, donor, -scoper_clone_size, -n_ebv_cells, scoper_clone_id_raw)

  fwrite(assignment_dt, file.path(scoper_dir, sprintf("%s_scoper_clone_assignments.tsv", donor_id)), sep = "\t")
  fwrite(clone_summary_dt, file.path(scoper_dir, sprintf("%s_scoper_clone_summary.tsv", donor_id)), sep = "\t")

  threshold_rows[[donor_id]] <- data.table(
    donor = donor_id,
    n_cells = donor_cell_dt[, uniqueN(lineage_obs_name)],
    threshold = threshold_value,
    threshold_source = threshold_source
  )
  assignment_rows[[donor_id]] <- assignment_dt
  clone_summary_rows[[donor_id]] <- clone_summary_dt
  message(sprintf("[scoper] done donor=%s n_cells=%s n_clones=%s threshold=%s source=%s",
                  donor_id,
                  donor_cell_dt[, uniqueN(lineage_obs_name)],
                  clone_summary_dt[, .N],
                  signif(threshold_value, 4),
                  threshold_source))
}

assignment_table <- rbindlist(assignment_rows, use.names = TRUE, fill = TRUE)
clone_summary_table <- rbindlist(clone_summary_rows, use.names = TRUE, fill = TRUE)
threshold_table <- rbindlist(threshold_rows, use.names = TRUE, fill = TRUE)

fwrite(assignment_table, file.path(tables_dir, sprintf("%s_scoper_clone_assignments.tsv", args$version)), sep = "\t")
fwrite(clone_summary_table, file.path(tables_dir, sprintf("%s_scoper_clone_summary.tsv", args$version)), sep = "\t")
fwrite(threshold_table, file.path(tables_dir, sprintf("%s_scoper_thresholds.tsv", args$version)), sep = "\t")
