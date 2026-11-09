# Combine VIRTUS3 calls with EBER-SCOPE targeted-PCR calls.
# Sources: 260416_build_pcr_augmented.py; EBER-SCOPE v0.1.0 merge_calls.py.
# Usage: python SCRIPT BCELLS.h5ad TARGETED_MANIFEST.csv OUTPUT.h5ad
# Manifest: sample, cell_summary (EBER-SCOPE *.cell_summary.tsv).
# One summary per GEX sample; sample names must match adata.obs['sample'].
# A qualifying read pair is sufficient; unique UMI count does not set positivity.
import sys
from pathlib import Path
import anndata as ad
import pandas as pd

input_h5ad, manifest, output = sys.argv[1:]
adata = ad.read_h5ad(input_h5ad)
targeted = []
for _, row in pd.read_csv(manifest, dtype={'sample': str}).iterrows():
    table = pd.read_csv(row['cell_summary'], sep='\t')
    table['sample'] = row['sample']
    table['barcode16'] = table['barcode'].astype(str).str.replace(r'-1$', '', regex=True)
    targeted.append(table[['sample', 'barcode16', 'n_qualifying_read_pairs',
                           'n_unique_umis', 'max_summed_alignment_span_fraction']])
targeted = pd.concat(targeted, ignore_index=True)

existing = adata.obs[['sample']].copy()
existing['sample'] = existing['sample'].astype(str)
existing['barcode16'] = adata.obs_names.to_series().str.split('-').str[0].values
existing['cell_id'] = adata.obs_names
combined = existing.merge(targeted, on=['sample', 'barcode16'], how='left',
                          sort=False, validate='one_to_one').set_index('cell_id')
combined = combined.reindex(adata.obs_names)
for column in ['n_qualifying_read_pairs', 'n_unique_umis']:
    adata.obs[column] = combined[column].fillna(0).astype(int).to_numpy()
adata.obs['max_summed_alignment_span_fraction'] = (
    combined['max_summed_alignment_span_fraction'].fillna(0.0).to_numpy())
adata.obs['EBV_detected_existing'] = adata.obs['EBV_detected'].astype(bool)
adata.obs['EBER1_pcr_v2_hit'] = adata.obs['n_qualifying_read_pairs'].gt(0)
adata.obs['EBV_detected_union_pcr_v2'] = (
    adata.obs['EBV_detected_existing'] | adata.obs['EBER1_pcr_v2_hit'])
adata.obs['EBV_detected_rescued_by_pcr_v2'] = (
    adata.obs['EBV_detected_union_pcr_v2'] & ~adata.obs['EBV_detected_existing'])
adata.obs['EBV_detection_source_pcr_v2'] = 'negative'
adata.obs.loc[adata.obs['EBV_detected_existing'], 'EBV_detection_source_pcr_v2'] = 'existing_only'
adata.obs.loc[adata.obs['EBER1_pcr_v2_hit'], 'EBV_detection_source_pcr_v2'] = 'pcr_only'
adata.obs.loc[adata.obs['EBV_detected_existing'] & adata.obs['EBER1_pcr_v2_hit'],
              'EBV_detection_source_pcr_v2'] = 'overlap'
Path(output).parent.mkdir(parents=True, exist_ok=True)
adata.write_h5ad(output)
