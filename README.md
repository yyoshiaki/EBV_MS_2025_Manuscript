# EBV-MS analysis code

Analysis code accompanying our [medRxiv paper](https://www.medrxiv.org/content/10.64898/2025.12.11.25342083v2).

- **Single-cell main:** [`261108_single_cell_main.py`](scripts/261108_single_cell_main.py) — 10x counts, QC, integration, B-cell selection and reference mapping.
- **EBER-SCOPE/PCR:** [`261108_run_eber_scope.sh`](scripts/261108_run_eber_scope.sh) → [`261108_single_cell_ebv_pcr.py`](scripts/261108_single_cell_ebv_pcr.py) — targeted detection and union with VIRTUS3 calls.
- **Final labels:** [`261108_single_cell_labels.py`](scripts/261108_single_cell_labels.py) → NEBULA, GSEA and scDRS; BCR clone analysis is also included.
- **GWAS:** `aou_*` scripts cover phenotypes, QC and REGENIE; MAGMA and EBV-DNA study meta-analysis are included separately.

Related tools: [scRefMap](https://github.com/yyoshiaki/screfmapping/tree/main/python), [EBER-SCOPE](https://github.com/yyoshiaki/EBER-SCOPE), [VIRTUS3](https://github.com/yyoshiaki/VIRTUS3).
Arguments are at the top of each script; sources are in [`docs/261107_provenance.tsv`](docs/261107_provenance.tsv).
Participant data, UKB GWAS, cross-cohort GWAS meta-analysis, HLA and PRS code are not included.

Updated November 8, 2026 (EST).

[GitHub](https://github.com/yyoshiaki/EBV_MS_2025_Manuscript) · [Archived code (DOI)](https://doi.org/10.5281/zenodo.24070191)
