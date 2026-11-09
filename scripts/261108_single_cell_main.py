# PBMC gene-expression QC, B-cell selection and reference mapping.
# Source: EBV_scRNAseq/scripts/260528_rebuild_bcell_h5mu.py
# Usage: python SCRIPT MANIFEST.csv CELL_METADATA.tsv CELLTYPIST.pkl REFERENCE.h5ad OUTPUT.h5ad
# Manifest: sample, obs_suffix, path_h5 (10x H5), path_virus (VIRTUS3 H5AD).
# Metadata: cell IDs in the first column, donor_demux and pre-QC donor/library annotations;
# provide donor_demux for every cell, including non-pooled libraries;
# exclude sample/obs_suffix columns, which are supplied by the manifest.
# Cell IDs use the 10x barcode plus '-' plus obs_suffix, as in the source build.
# Reference: the study B-cell Symphony object, with cluster_L2 labels.
import re
import sys
from pathlib import Path

import anndata as ad
import celltypist
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import symphonypy as sp

manifest, metadata, celltypist_model, reference, output = sys.argv[1:]
np.random.seed(0)
df_build = pd.read_csv(manifest, dtype={'sample': str, 'obs_suffix': str})
gex = {}
for _, row in df_build.iterrows():
    cells = sc.read_10x_h5(row['path_h5'], gex_only=False)
    cells.var_names_make_unique()
    cells = cells[:, cells.var['feature_types'] == 'Gene Expression'].copy()
    cells.obs['sample'] = row['sample']
    cells.obs['obs_suffix'] = row['obs_suffix']
    gex[row['obs_suffix']] = cells
adata = sc.concat(gex, join='outer', index_unique='-')
adata.obs = adata.obs.join(pd.read_csv(metadata, sep='\t', index_col=0), how='left')
if adata.obs['donor_demux'].isna().any():
    raise ValueError('CELL_METADATA must provide donor_demux for every GEX cell, including non-pooled libraries.')

# Gene-expression QC and integration
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt", "ribo", "hb"],
    percent_top=None,
    log1p=False,
    inplace=True,
)
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e5)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=5000, batch_key="sample")
patterns_to_exclude = [
    "^IGKV", "^IGLV", "^IGHV", "^IGLC",
    "^TRAV", "^TRAJ",
    "^TRBV", "^TRBD", "^TRBJ",
    "^TRGV", "^TRGJ",
    "^TRDV", "^TRDD", "^TRDJ",
    "^MT-", "^RPL", "^RPS",
]
adata.var["highly_variable"] = (
    (~adata.var_names.str.contains("|".join(patterns_to_exclude), flags=re.IGNORECASE))
    & adata.var["highly_variable"]
)

sc.tl.pca(adata)
sce.pp.harmony_integrate(adata, "sample", max_iter_harmony=50)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50, use_rep="X_pca_harmony")
sc.tl.umap(adata)
sc.tl.leiden(adata)

# Viral counts and EBV detection
dict_adata_virus = {}
for _, row in df_build.iterrows():
    adata_virus_single = sc.read(row["path_virus"])
    # Normalize only the terminal 10x barcode token here.
    # The sample suffix is added once by sc.concat(index_unique="-"), matching GEX obs_names.
    adata_virus_single.obs_names = (
        adata_virus_single.obs_names.to_series().astype(str).str.replace(r"-(\d+)$", "-1", regex=True)
    ).values
    dict_adata_virus[row["obs_suffix"]] = adata_virus_single

adata_virus = sc.concat(dict_adata_virus, index_unique="-")
expr = pd.DataFrame(
    adata_virus.X.toarray() if hasattr(adata_virus.X, "toarray") else adata_virus.X,
    index=adata_virus.obs_names,
    columns=adata_virus.var_names,
)
expr_sum = expr.groupby(level=0).sum()
adata_virus_agg = ad.AnnData(
    X=expr_sum.values,
    obs=pd.DataFrame(index=expr_sum.index),
    var=adata_virus.var.copy(),
)
sc.pp.calculate_qc_metrics(adata_virus_agg, percent_top=None, log1p=False, inplace=True)
adata_virus_agg.obs.columns = adata_virus_agg.obs.columns + "_EBV"

virus_df = sc.get.obs_df(
    adata_virus_agg,
    keys=list(adata_virus_agg.obs.columns) + list(adata_virus_agg.var.index),
)
adata.obs = pd.merge(
    adata.obs,
    virus_df,
    left_index=True,
    right_index=True,
    how="left",
)
virus_fill_cols = list(adata_virus_agg.obs.columns) + list(adata_virus_agg.var.index)
adata.obs[virus_fill_cols] = adata.obs[virus_fill_cols].fillna(0)
adata.obs["total_counts_EBV_log1p"] = np.log1p(adata.obs["total_counts_EBV"])
adata.obs["EBV_detected"] = (adata.obs["total_counts_EBV"] > 0).astype(bool)


# CellTypist B-cell selection
bcell_labels = [
    "Age-associated B cells",
    "B cells",
    "Cycling B cells",
    "Follicular B cells",
    "Germinal center B cells",
    "Large pre-B cells",
    "Memory B cells",
    "Naive B cells",
    "Plasma cells",
    "Plasmablasts",
    "Pre-pro-B cells",
    "Pro-B cells",
    "Proliferative germinal center B cells",
    "Small pre-B cells",
    "Transitional B cells",
]

adata_ct = adata.copy()
adata_ct.X = adata_ct.layers["counts"].copy()
sc.pp.normalize_total(adata_ct, target_sum=1e4)
sc.pp.log1p(adata_ct)
adata_ct.raw = adata_ct.copy()
predictions = celltypist.annotate(
    adata_ct,
    model=str(celltypist_model),
    majority_voting=True,
)
df_pred = predictions.predicted_labels.copy()
df_pred.columns = [
    "predicted_labels_Immune_All_Low",
    "over_clustering_Immune_All_Low",
    "majority_voting_Immune_All_Low",
]
adata.obs = pd.merge(
    adata.obs,
    df_pred,
    how="left",
    left_index=True,
    right_index=True,
)
adata = adata[adata.obs["majority_voting_Immune_All_Low"].isin(bcell_labels)].copy()

# B-cell reference mapping
adata_ref = sc.read(reference)
# Preserve the local Harmony UMAP before Symphony ingest overwrites X_umap.
adata.obsm["X_umap_local"] = adata.obsm["X_umap"].copy()
sp.tl.map_embedding(adata, adata_ref, key="sample")
sp.tl.transfer_labels_kNN(adata, adata_ref, "cluster_L2")
sp.tl.ingest(adata, adata_ref)
adata.obsm["X_umap_refmap"] = adata.obsm["X_umap"].copy()
adata.obsm["X_umap"] = adata.obsm["X_umap_local"].copy()

adata.obs.columns = adata.obs.columns.str.replace('/', '_', regex=False)
for column in adata.obs.columns:
    if adata.obs[column].dtype == 'object' or pd.api.types.is_string_dtype(adata.obs[column]):
        adata.obs[column] = adata.obs[column].astype(str)

Path(output).parent.mkdir(parents=True, exist_ok=True)
adata.write_h5ad(output)
