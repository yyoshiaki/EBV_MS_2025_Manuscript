# Production source: EBV_scRNAseq/scripts/260819_prepare_pbmc_ebv_b_nebula_hybrid.py
# Prepared 2026-10-15 EDT: Explicit input/config/output arguments; boolean status check; no model/filter changes.
# Statistical settings retained; source hash in docs/261015_provenance.tsv.

# Reference:
# - scripts/260327_prepare_pbmc_ebv_b_nebula_inputs.py
# - scripts/260416_build_pcr_augmented.py

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmwrite

parser = argparse.ArgumentParser(description="Prepare PCR-union EBV+/- inputs for Nebula.")
parser.add_argument("--analysis-name", required=True)
parser.add_argument("--input-h5ad", required=True)
parser.add_argument("--cluster-map", required=True)
parser.add_argument("--config", required=True)
parser.add_argument("--output-root", required=True)
parser.add_argument("--configs", default="")
parser.add_argument("--exclude-donors", default="")
args = parser.parse_args()

input_h5ad = Path(args.input_h5ad)
hybrid_map_path = Path(args.cluster_map)
output_root = Path(args.output_root) / args.analysis_name
output_dir = output_root / "nebula" / "inputs"
config_tsv = Path(args.config)
output_dir.mkdir(parents=True, exist_ok=True)

config_df = pd.read_csv(config_tsv, sep="\t")
selected_configs = {x for x in args.configs.split(",") if x}
if selected_configs - set(config_df["name"]):
    raise ValueError("Unknown requested config name.")
if selected_configs:
    config_df = config_df.loc[config_df["name"].isin(selected_configs)].copy()
configs = []
for _, row in config_df.iterrows():
    configs.append(
        {
            "name": row["name"],
            "clusters": str(row["clusters"]).split("|"),
            "min_cells_per_gene": int(row["min_cells_per_gene"]),
            "min_cells_per_donor": int(row["min_cells_per_donor"]),
            "n_top": int(row["n_top"]),
        }
    )

exclude_donor_demux = {"doublet", "unassigned", "pool"}
exclude_donor_demux.update(x for x in args.exclude_donors.split(",") if x)

adata = ad.read_h5ad(input_h5ad, backed="r")
if not pd.api.types.is_bool_dtype(adata.obs["EBV_detected_union_pcr_v2"]):
    raise ValueError("EBV union status must be boolean, not text.")
hybrid_map = pd.read_csv(hybrid_map_path, sep="\t", index_col=0)["cluster_L2_hybrid"]
assert adata.obs_names.equals(hybrid_map.index)
cluster_l2 = hybrid_map.astype(str)
var_names = adata.var_names.to_numpy().astype(str)
obs_names = adata.obs_names.to_numpy().astype(str)

summary_rows = []

for cfg in configs:
    initial_mask = (
        cluster_l2.isin(cfg["clusters"])
        & ~adata.obs["donor_demux"].astype(str).isin(exclude_donor_demux)
    )
    initial_idx = np.flatnonzero(initial_mask.to_numpy())
    initial_meta = adata.obs.iloc[initial_idx][
        ["sample", "donor_demux", "cluster_L2", "EBV_detected_union_pcr_v2"]
    ].copy()
    initial_meta["cluster_L2"] = cluster_l2.iloc[initial_idx].to_numpy()
    initial_meta = initial_meta.rename(columns={"EBV_detected_union_pcr_v2": "EBV_detected"})
    donor_counts = initial_meta["donor_demux"].astype(str).value_counts()
    keep_donors = donor_counts[donor_counts >= cfg["min_cells_per_donor"]].index
    donor_keep_mask = initial_meta["donor_demux"].astype(str).isin(keep_donors).to_numpy()
    kept_meta = initial_meta.loc[donor_keep_mask].copy()
    kept_idx = initial_idx[donor_keep_mask]

    counts = adata.layers["counts"][kept_idx, :].tocsr()
    normalized = adata.X[kept_idx, :].tocsr()
    gene_detected = np.asarray((counts > 0).sum(axis=0)).ravel()
    keep_genes = gene_detected >= cfg["min_cells_per_gene"]

    count_matrix = counts[:, keep_genes].transpose().tocsr()
    mean_log_expr = np.asarray(normalized[:, keep_genes].mean(axis=0)).ravel()
    kept_gene_names = var_names[keep_genes]
    library_size = np.asarray(counts.sum(axis=1)).ravel()

    kept_meta = kept_meta.reset_index(drop=True)
    kept_meta["cell_key"] = obs_names[kept_idx]
    kept_meta["sample"] = kept_meta["sample"].astype(str)
    kept_meta["donor_demux"] = kept_meta["donor_demux"].astype(str)
    kept_meta["cluster_L2"] = kept_meta["cluster_L2"].astype(str)
    kept_meta["EBV_detected"] = kept_meta["EBV_detected"].astype(bool).map({False: "False", True: "True"})
    kept_meta["library_size_counts"] = library_size

    order = np.lexsort((kept_meta["cell_key"], kept_meta["sample"], kept_meta["donor_demux"]))
    kept_meta = kept_meta.iloc[order].reset_index(drop=True)
    count_matrix = count_matrix[:, order]

    gene_table = pd.DataFrame(
        {
            "gene": kept_gene_names,
            "mean": mean_log_expr,
            "n_cells_by_counts": gene_detected[keep_genes],
        }
    )

    mmwrite(output_dir / f"{cfg['name']}_counts.mtx", count_matrix)
    kept_meta.to_csv(output_dir / f"{cfg['name']}_cell_metadata.tsv", sep="\t", index=False)
    gene_table.to_csv(output_dir / f"{cfg['name']}_genes.tsv", sep="\t", index=False)

    summary_rows.append(
        {
            "name": cfg["name"],
            "analysis_name": args.analysis_name,
            "excluded_donors": ",".join(sorted(exclude_donor_demux)),
            "n_cells": int(kept_meta.shape[0]),
            "n_genes": int(gene_table.shape[0]),
            "n_donors": int(kept_meta["donor_demux"].nunique()),
            "n_samples": int(kept_meta["sample"].nunique()),
            "n_ebv_true": int((kept_meta["EBV_detected"] == "True").sum()),
            "n_ebv_false": int((kept_meta["EBV_detected"] == "False").sum()),
            "min_cells_per_donor": int(kept_meta["donor_demux"].value_counts().min()),
        }
    )

pd.DataFrame(summary_rows).to_csv(output_dir / "nebula_input_summary.tsv", sep="\t", index=False)
