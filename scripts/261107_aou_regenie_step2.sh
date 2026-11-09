#!/usr/bin/env bash
# Source: 250825_regenie_step2.ipynb, cell 15 (zero-based); REGENIE 3.2.4 source image.
# CHROM argument: numeric 1-22, without the chr prefix.
# Usage: bash SCRIPT PGEN_PREFIX CHROM INPUT_DIR SNPLIST.gz STEP1_DIR OUTPUT_DIR PREFIX
# Use absolute paths; one chromosome per output directory. Run with Step 1 cwd
# so its relative prediction/LOCO paths resolve.
# MAC >=40 is global, before ancestry keep lists; no later ancestry-specific MAC filter.
# As in the source, .bi.pvar CHROM is expected numeric; PSAM has IID/SEX without FID.
set -euo pipefail
pgen_file="${1}.pgen"
psam_file="${1}.psam"
pvar_file="${1}.pvar"
chrom="$2"
input_dir="$3"
snplist="$4"
regenie_step1="$5"
output_dir="$6"
prefix="$7"
pheno_bt_file="$input_dir/phenotype.bt.txt"
pheno_ct_file="$input_dir/phenotype.ct.txt"
cov_file="$input_dir/covariates.txt"
mkdir -p "$output_dir"
cd "$regenie_step1"
mkdir -p "${output_dir}/plink_maf_hwe_filtered"

tmp_prefix="${output_dir}/plink_maf_hwe_filtered/chr${chrom}_filtered"

plink2 \
      --pgen    "${pgen_file}" \
      --psam    "${psam_file}" \
      --pvar    "${pvar_file}" \
      --make-pgen \
      --max-alleles 2 \
      --min-alleles 2 \
      --mac 40 \
      --out "${tmp_prefix}.bi"

mv "${tmp_prefix}.bi.pvar" "${tmp_prefix}.bi.pvar.orig"

awk 'BEGIN { OFS = "\t" }
     /^##/ { print; next }
     /^#CHROM/ { print; next }
     {
       chr = $1
       pos = $2
       ref = $4
       alt = $5
       new_id = "chr" chr "-" pos "-" ref "-" alt
       printf("%s\t%s\t%s", chr, pos, new_id)
       for (i = 4; i <= NF; i++) {
         printf("\t%s", $i)
       }
       printf("\n")
     }
    ' "${tmp_prefix}.bi.pvar.orig" > "${tmp_prefix}.bi.pvar"

zcat "${snplist}" > "${output_dir}/snplist.ids"

plink2 \
  --pfile "${tmp_prefix}.bi" \
  --extract "${output_dir}/snplist.ids" \
  --make-pgen \
  --out "${tmp_prefix}"

cp "${tmp_prefix}.psam" "${tmp_prefix}.psam.orig"
awk 'NR==1 { print "#FID IID SEX"; next } { print 0, $1, $2 }' \
  "${tmp_prefix}.psam.orig" \
  > "${tmp_prefix}.psam"

ancestries=(eur afr amr eas sas mid)

for anc in "${ancestries[@]}"; do

    keep_file="${input_dir}/fid_iid_filtered_ancestry_pred_${anc}.txt"

    out_prefix="${output_dir}/${prefix}_chr${chrom}_${anc}"

    regenie \
    --step 2 \
    --pgen     "${tmp_prefix}" \
    --phenoFile "${pheno_bt_file}" \
    --phenoCol  "EBVpos" \
    --covarFile "${cov_file}" \
    --catCovarList "site_id" \
    --keep "${keep_file}" \
    --pred     "${regenie_step1}/${prefix}_${anc}_bt_step1_pred.list" \
    --bsize 400 \
    --verbose \
    --ref-first \
    --bt \
    --firth --approx \
    --pThresh 0.01 \
    --out "${out_prefix}_bt_step2"

    regenie \
    --step   2 \
    --pgen     "${tmp_prefix}" \
    --phenoFile "${pheno_ct_file}" \
    --phenoCol  "log10_norm_EBV_reads" \
    --covarFile "${cov_file}" \
    --catCovarList "site_id" \
    --keep "${keep_file}" \
    --pred     "${regenie_step1}/${prefix}_${anc}_ct_step1_pred.list" \
    --bsize    400 \
    --verbose  \
    --ref-first \
    --out     "${out_prefix}_ct_step2"

done
