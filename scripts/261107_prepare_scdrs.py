# Source: EBV_scRNAseq/scripts/260818_prepare_magic_scdrs_pcr_union.py
import argparse
from pathlib import Path

import numpy as np
import scanpy as sc
import scanpy.external as sce


parser = argparse.ArgumentParser(description="Prepare frozen PCR-union counts for MAGIC/scDRS")
parser.add_argument("--input-h5ad", required=True)
parser.add_argument("--output-dir", required=True)
args = parser.parse_args()
input_h5ad = Path(args.input_h5ad)
output_dir = Path(args.output_dir)
output_dir.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(input_h5ad)
assigned_qc = ~adata.obs["donor_demux"].astype(str).isin(
    ["doublet", "unassigned", "pool", "nan", ""]
)
adata = adata[assigned_qc].copy()
adata.obs["EBV_pos"] = (
    adata.obs["EBV_detected_union_pcr_v2"].astype(str).eq("True")
)

cov = adata.obs[["n_genes_by_counts"]].copy()
cov["const"] = 1
cov.to_csv(output_dir / "cov.tsv", sep="\t")
adata.obs[["donor_demux", "cluster_L2", "EBV_pos"]].to_csv(
    output_dir / "assigned_union_obs.tsv", sep="\t"
)

adata.X = adata.layers["counts"].astype(np.float32)
adata.layers.clear()
adata.raw = None
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sce.pp.magic(adata, name_list="all_genes", knn=5)
adata.raw = None
adata.write_h5ad(output_dir / "res.magic.scdrs.noraw.h5ad")

