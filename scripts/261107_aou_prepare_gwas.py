# Source: 250825_prep_regenie_vis_metadata_rmHIVid_binary.ipynb, selected cells.
# Usage: python SCRIPT METADATA.tsv RELATED.tsv FLAGGED.tsv OUTPUT_DIR
# Metadata contains the original joined fields, PC1-PC10 and boolean HIV /
# Immunodeficiency flags serialized as True/False. Immunodeficiency: ICD10 D80-D89; HIV: phecode ID_057.1.
# Related and flagged inputs contain sample_id and s, respectively.
import sys
from pathlib import Path
import numpy as np
import pandas as pd

metadata, related, flagged, output = sys.argv[1:]
out = Path(output)
out.mkdir(parents=True, exist_ok=True)
df_ebv = pd.read_csv(metadata, sep='\t')
df_rel = pd.read_csv(related, sep='\t')
df_flagged = pd.read_csv(flagged, sep='\t')

df_ebv['norm_EBV_reads'] = df_ebv['EBV_reads'] / df_ebv['mean_coverage']
df_ebv['log10_norm_EBV_reads'] = np.log10(df_ebv['norm_EBV_reads'] + 1)
df_ebv['EBVpos'] = (df_ebv['EBV_reads'] > 0) * 1
df_ebv['age_at_collection'] = (
    pd.to_datetime(df_ebv['biosample_collection_date'].str.split(' ').str.get(0))
    - pd.to_datetime(df_ebv['date_of_birth'])
).dt.days / 365.2
df_ebv['FID'] = 0
df_ebv['IID'] = df_ebv['person_id']
df_ebv['sex_at_birth_male'] = (df_ebv['sex_at_birth'] == 'Male') * 1
df_ebv_bld = df_ebv[df_ebv['sample_source'] == 'WHOLE BLOOD'].copy()

# As in the source, exclusions are enforced through keep lists, not these files.
covariates = ['FID', 'IID', 'sex_at_birth_male', 'age_at_collection', 'site_id']
covariates += [f'PC{i}' for i in range(1, 11)]
df_ebv_bld[covariates].to_csv(out / 'covariates.txt', sep=' ', index=False)
df_ebv_bld[['FID', 'IID', 'EBVpos']].to_csv(out / 'phenotype.bt.txt', sep=' ', index=False)
df_ebv_bld[['FID', 'IID', 'log10_norm_EBV_reads']].to_csv(
    out / 'phenotype.ct.txt', sep=' ', index=False)

df_exclude = df_ebv.loc[
    (~df_ebv['sex_at_birth'].isin(['Male', 'Female']))
    | (df_ebv['sample_source'] != 'WHOLE BLOOD')
    | (df_ebv['age_at_collection'].isna()), ['FID', 'IID']]
df_filtered = df_ebv_bld[~df_ebv_bld.IID.isin(df_exclude.IID)]
df_filtered = df_filtered[(df_filtered['Immunodeficiency'] == False)
                          & (df_filtered['HIV'] == False)]
df_filtered = df_filtered[~df_filtered.IID.isin(df_rel['sample_id'])]
df_filtered = df_filtered[~df_filtered.IID.isin(df_flagged['s'])]
for anc in df_ebv['ancestry_pred'].unique():
    df_filtered.loc[df_filtered['ancestry_pred'] == anc, ['FID', 'IID']].to_csv(
        out / f'fid_iid_filtered_ancestry_pred_{anc}.txt',
        sep=' ', header=False, index=False)
