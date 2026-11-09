#!/usr/bin/env bash
# Source: 250825_regenie_step1.ipynb, cell 11 (zero-based); REGENIE 3.2.4 source image.
# Usage: bash SCRIPT QC_ROOT INPUT_DIR OUTPUT_DIR PREFIX
# Use absolute directory paths. Run once in a dedicated output directory.
# Source PLINK PSAM columns before rewrite are IID and SEX, without FID.
# Keep this directory in place: Step 2 uses prediction lists and LOCO files here.
set -euo pipefail
pgen_pruned="$1"
input_dir="$2"
output_dir="$3"
prefix="$4"
pheno_bt_file="$input_dir/phenotype.bt.txt"
pheno_ct_file="$input_dir/phenotype.ct.txt"
cov_file="$input_dir/covariates.txt"
mkdir -p "$output_dir"
cd "$output_dir"
ancestries=(eur afr amr eas sas mid)

for anc in "${ancestries[@]}"; do

    ls "${pgen_pruned}"/*/chr*_"${anc}"_filtered_LDprune.pvar \
      | sed 's/\.pvar$//' \
      > "pmerge-list_${anc}.txt"

    plink2 \
      --pmerge-list pmerge-list_${anc}.txt \
      --make-pgen \
      --out merged_LDprune.allchr_${anc}

    cp "merged_LDprune.allchr_${anc}.psam" "merged_LDprune.allchr_${anc}.psam.orig"
    awk 'NR==1 { print "#FID IID SEX"; next } { print 0, $1, $2 }' \
          "merged_LDprune.allchr_${anc}.psam.orig" \
          > "merged_LDprune.allchr_${anc}.psam"

    out_prefix="${prefix}_${anc}"

    regenie \
    --step 1 \
    --pgen     "merged_LDprune.allchr_${anc}" \
    --phenoFile "${pheno_bt_file}" \
    --phenoCol  "EBVpos" \
    --covarFile "${cov_file}" \
    --catCovarList "site_id" \
    --bsize 1000 \
    --lowmem \
    --ref-first \
    --bt \
    --out "${out_prefix}_bt_step1"

    regenie \
    --step   1 \
    --pgen     "merged_LDprune.allchr_${anc}" \
    --phenoFile "${pheno_ct_file}" \
    --phenoCol  "log10_norm_EBV_reads" \
    --covarFile "${cov_file}" \
    --catCovarList "site_id" \
    --bsize    1000 \
    --ref-first \
    --force-step1 \
    --out     "${out_prefix}_ct_step1"

done
