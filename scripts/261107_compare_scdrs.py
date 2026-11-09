# Source: EBV_scRNAseq/scripts/260819_render_scdrs_hybrid.py
import argparse
from pathlib import Path

import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


parser = argparse.ArgumentParser(description="Compare frozen scDRS scores by PCR-union status")
parser.add_argument("--prepared-dir", required=True)
parser.add_argument("--hybrid-map", required=True)
parser.add_argument("--ebv-score", required=True)
parser.add_argument("--ms-score", required=True)
parser.add_argument("--output-dir", required=True)
args = parser.parse_args()
data_dir = Path(args.prepared_dir)
hybrid_map_path = Path(args.hybrid_map)
table_dir = Path(args.output_dir)
table_dir.mkdir(parents=True, exist_ok=True)

obs = pd.read_csv(data_dir / "assigned_union_obs.tsv", sep="\t", index_col=0)
hybrid_map = pd.read_csv(hybrid_map_path, sep="\t", index_col=0)["cluster_L2_hybrid"]
if not obs.index.isin(hybrid_map.index).all():
    raise ValueError("Hybrid labels must cover every assigned cell.")
obs["cluster_L2"] = hybrid_map.reindex(obs.index).astype(str)
obs["EBV_pos"] = obs["EBV_pos"].astype(str).eq("True")
cluster_order = [
    "Bnaive SOX4", "Bnaive", "Bnaive MX1", "Bnaive preASC", "Busm",
    "Bmem int", "Bsm FCER2", "Bsm", "ABC", "ASC",
]
score_specs = [
    ("EBV_bt", Path(args.ebv_score)),
    ("MS", Path(args.ms_score)),
]

for score_name, score_path in score_specs:
    score = pd.read_csv(score_path, sep="\t", index_col=0)[["norm_score"]]
    plot_df = obs.join(score, how="inner").dropna(subset=["norm_score"])
    score_col = f"{score_name}_magic"
    plot_df = plot_df.rename(columns={"norm_score": score_col})
    positive_counts = plot_df.loc[plot_df["EBV_pos"], "cluster_L2"].value_counts()
    keep_clusters = [x for x in cluster_order if positive_counts.get(x, 0) > 30]
    plot_df = plot_df.loc[plot_df["cluster_L2"].isin(keep_clusters)].copy()
    if not keep_clusters:
        raise ValueError(f"No cluster has more than 30 positive cells for {score_name}.")
    stats_rows = []
    for cluster in keep_clusters:
        sub = plot_df.loc[plot_df["cluster_L2"] == cluster]
        positive = sub.loc[sub["EBV_pos"], score_col]
        negative = sub.loc[sub["EBV_pos"].eq(False), score_col]
        if negative.empty:
            raise ValueError(f"No negative cells in {cluster} for {score_name}.")
        statistic, pvalue = mannwhitneyu(
            positive, negative, alternative="two-sided", method="asymptotic"
        )
        stats_rows.append({
            "score": score_name,
            "cluster_L2": cluster,
            "n_positive": len(positive),
            "n_negative": len(negative),
            "U": statistic,
            "pvalue": pvalue,
        })
    stats_df = pd.DataFrame(stats_rows)
    stats_df["p_adj"] = multipletests(stats_df["pvalue"], method="fdr_bh")[1]
    stats_df.to_csv(
        table_dir / f"polygenic_{score_name}_mannwhitney_by_cluster.tsv",
        sep="\t",
        index=False,
    )
