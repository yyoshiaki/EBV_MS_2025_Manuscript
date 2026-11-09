#!/usr/bin/env bash
# Sources: EBV_scRNAseq/scripts/260818_scdrs_pcr_union_EBV_bt.sbatch and ..._MS.sbatch
# Usage: bash scripts/261107_run_scdrs.sh PREPARED_DIR GENE_SET.gs OUTPUT_DIR
set -euo pipefail
if [[ $# -ne 3 ]]; then
    echo "Usage: $0 PREPARED_DIR GENE_SET.gs OUTPUT_DIR" >&2
    exit 2
fi
mkdir -p "$3"
scdrs compute-score \
    --h5ad-file "$1/res.magic.scdrs.noraw.h5ad" \
    --h5ad-species human \
    --gs-file "$2" \
    --gs-species human \
    --cov-file "$1/cov.tsv" \
    --flag-filter-data True \
    --flag-raw-count False \
    --flag-return-ctrl-raw-score False \
    --flag-return-ctrl-norm-score True \
    --out-folder "$3"
