# Production source: EBV_GWAS/script/260311_EBV_metaanalysis.R
# Prepared 2026-10-15 EDT: CSV/output CLI and equivalent base-R table preprocessing; omit forest/funnel drawing; keep metabin OR/DL/Hartung-Knapp settings.
# Statistical settings retained; source hash in docs/261015_provenance.tsv.

library(meta)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
df <- read.csv(args[[1]], check.names = FALSE)
df <- df[!is.na(df$MS_pos), ]
df_meta <- transform(df, MS_total = MS_pos + MS_neg, Ctrl_total = Ctrl_pos + Ctrl_neg)
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
m1 <- metabin(
  event.e = MS_pos,
  n.e     = MS_total,
  event.c = Ctrl_pos,
  n.c     = Ctrl_total,
  studlab = Study,
  data    = df_meta,
  sm      = "OR",
  method.tau = "DL",
  hakn = TRUE
)

write.csv(data.frame(Study = df_meta$Study, logOR = m1$TE, se = m1$seTE), file.path(out_dir, "study_effects.csv"), row.names = FALSE)
write.csv(data.frame(model = c("fixed", "random"), logOR = c(m1$TE.common, m1$TE.random), se = c(m1$seTE.common, m1$seTE.random), lower_logOR = c(m1$lower.common, m1$lower.random), upper_logOR = c(m1$upper.common, m1$upper.random)), file.path(out_dir, "pooled_effects.csv"), row.names = FALSE)
saveRDS(m1, file.path(out_dir, "meta_fit.rds"))
capture.output(sessionInfo(), file = file.path(out_dir, "R_sessionInfo.txt"))
