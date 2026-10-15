# Existing environments

Updated 2026-10-15 EDT. Use an existing environment; do not create a `.venv` or install into the notes vault. These version records are not a fully tested lockfile.

| Component | Version/evidence |
| --- | --- |
| Single-cell Python runtime queried on production host | anndata 0.10.9; numpy 1.26.4; pandas 2.0.3; scipy 1.13.0; scanpy 1.10.2 |
| NEBULA runtime from frozen SLE validation sessionInfo | R 4.4.1; Matrix 1.7-0; nebula 1.5.6; Bioconductor 3.19 container family |
| Core GSEA requirements | dplyr, tibble, clusterProfiler, ReactomePA, org.Hs.eg.db; exact annotation/package versions not yet frozen for release |
| Core BCR requirements | data.table, scoper, shazam; exact production versions not yet frozen for release |
| Meta-analysis local smoke runtime | R 4.4.1; meta 8.1-0; use identical version for a frozen rerun; original production meta version unverified |

The retained preparation script imports only anndata, numpy, pandas and scipy. No Scanpy analysis is invoked by this export. NEBULA/GSEA/BCR packages belong in an R environment compatible with the recorded production stack; automatic installation scripts are intentionally omitted. `hakn=TRUE` in the original meta code is retained and may emit a deprecation warning in current meta versions. Record `sessionInfo()` and versions alongside each authorized rerun. The core fit/meta entry points record sessionInfo automatically.
