# Study labels: retain the submitted label for shared cells and use rebuilt labels for new cells.
# Source: EBV_scRNAseq/scripts/260819_build_hybrid_cluster_map.py
# Usage: python SCRIPT SUBMITTED_OBS.csv.gz REBUILT.h5ad OUTPUT.tsv.gz
# Submitted table: cell IDs as index and cluster_L2. Rebuilt H5AD: PCR integration output.
# The output follows the H5AD cell order for NEBULA/scDRS preparation.
import sys
from pathlib import Path
import pandas as pd
import anndata as ad

submitted_path, rebuilt_path, output = sys.argv[1:]
submitted = pd.read_csv(submitted_path, index_col=0)
adata = ad.read_h5ad(rebuilt_path, backed='r')
rebuild = adata.obs[['cluster_L2']].copy()
adata.file.close()
submitted['cluster_L2'] = submitted['cluster_L2'].astype(str)
rebuild['cluster_L2'] = rebuild['cluster_L2'].astype(str)
hybrid = rebuild[['cluster_L2']].rename(columns={'cluster_L2': 'cluster_L2_rebuild'})
hybrid['cluster_L2_submitted'] = submitted['cluster_L2'].reindex(hybrid.index)
hybrid['cluster_L2_hybrid'] = hybrid['cluster_L2_submitted'].fillna(hybrid['cluster_L2_rebuild'])
hybrid['cluster_label_source'] = '260528_new_only'
hybrid.loc[hybrid['cluster_L2_submitted'].notna(), 'cluster_label_source'] = '260408_submitted'
Path(output).parent.mkdir(parents=True, exist_ok=True)
hybrid.to_csv(output, sep='\t')
