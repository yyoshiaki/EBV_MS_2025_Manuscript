#!/usr/bin/env bash
# Source: hc865/EBV/magma.sh
# Preserve MAGMA SNP/P columns and variant-specific N; no build conversion is performed.
# Usage: bash scripts/261107_run_magma.sh LD_BFILE SUMSTATS GENE_ANNOT OUTPUT_PREFIX
set -euo pipefail
if [[ $# -ne 4 ]]; then
    echo "Usage: $0 LD_BFILE SUMSTATS GENE_ANNOT OUTPUT_PREFIX" >&2
    exit 2
fi
magma --bfile "$1" \
    --pval "$2" use='SNP,P' ncol='N' \
    --gene-annot "$3" \
    --out "$4"
