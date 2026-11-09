#!/usr/bin/env bash
# Source: 250825_variantQC_LDpruning.ipynb, cell 9 (zero-based).
# CHROM argument: numeric 1-22, without the chr prefix.
# Usage: bash SCRIPT PGEN_PREFIX CHROM INPUT_DIR QC_ROOT
# One chromosome per call; source pvar CHROM values include the chr prefix.
# LD pruning uses a window of 1000 variants.
set -euo pipefail
pgen_file="${1}.pgen"
psam_file="${1}.psam"
pvar_file="${1}.pvar"
chrom="$2"
input_dir="$3"
qc_dir="$4/$chrom"
mkdir -p "$qc_dir"
pvar_file_rename="${qc_dir}/renamed.pvar"
awk '  \
  /^##/ { print; next }  \
  /^#CHROM/ { print; next }  \
  { \
    split($5, A, ","); \
    $3 = $1 "-" $2 "-" $4 "-" A[1]; \
    print; \
  }' OFS="\t" "${pvar_file}" > "${pvar_file_rename}"

pvar_file="${pvar_file_rename}"

ancestries=(eur afr amr eas sas mid)

for anc in "${ancestries[@]}"; do

  keep_file="${input_dir}/fid_iid_filtered_ancestry_pred_${anc}.txt"
  if [[ "${anc}" == "afr" ]]; then
    ld_opt="--indep-pairwise 1000 100 0.15"
  else
    ld_opt="--indep-pairwise 1000 100 0.3"
  fi

  filt_pref="${qc_dir}/chr${chrom}_${anc}_filtered"
  plink2 \
    --pgen    "${pgen_file}" \
    --psam    "${psam_file}" \
    --pvar    "${pvar_file_rename}" \
    --keep "${keep_file}" \
    --maf     0.01 \
    --max-maf 0.99 \
    --hwe     1e-6 keep-fewhet \
    --max-alleles 2 \
    --min-alleles 2 \
    --geno    0.05 \
    --write-snplist allow-dups \
    --freq \
    --hardy \
    --missing \
    --out     "${filt_pref}"

  filt_pref_LDprune="${qc_dir}/chr${chrom}_${anc}_filtered_LDprune"
plink2 \
    --pgen    "${pgen_file}" \
    --psam    "${psam_file}" \
    --pvar    "${pvar_file_rename}" \
    --keep "${keep_file}" \
    --maf     0.01 \
    --max-maf 0.99 \
    --hwe     1e-6 keep-fewhet \
    --geno 0.05 \
    --max-alleles 2 \
    --min-alleles 2 \
    ${ld_opt} \
    --out     "${filt_pref_LDprune}"

    plink2 \
  --pgen "${pgen_file}" \
  --psam "${psam_file}" \
  --pvar "${pvar_file_rename}" \
   --keep "${keep_file}" \
  --extract "${filt_pref_LDprune}.prune.in" \
  --make-pgen \
  --out "${filt_pref_LDprune}"

done
