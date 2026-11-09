#!/usr/bin/env bash
# EBER1 detection in targeted 10x 5-prime v3 paired-end reads.
# Source: EBER-SCOPE v0.1.0; study anchor and alignment-span thresholds.
# Usage: bash SCRIPT R1.fastq.gz R2.fastq.gz SAMPLE WHITELIST.txt.gz OUTPUT_PREFIX
# SAMPLE must match the GEX sample column used for PCR integration.
set -euo pipefail
mkdir -p "$(dirname "$5")"
eber-scope detect \
    --r1 "$1" --r2 "$2" --sample "$3" \
    --barcode-whitelist "$4" --chemistry 5p-v3 \
    --anchor-mismatches 1 --threshold 0.60 \
    --output-prefix "$5"
