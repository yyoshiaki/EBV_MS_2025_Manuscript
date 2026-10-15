# Minimal input contracts

Updated 2026-10-15 EDT. Real inputs are not distributed in this repository.

## Differential expression

- H5AD: cells × genes, unique gene symbols and cell identifiers; sparse nonnegative integer `layers["counts"]`; sparse normalized log-expression `X`. `X` must preserve the frozen production preprocessing; the preparation script does not normalize it.
- `obs`: `sample`, `donor_demux`, `cluster_L2`, boolean `EBV_detected_union_pcr_v2`. Strings such as `"False"` are rejected because Python would otherwise treat nonempty strings as true. Donor/sample identifiers must be pseudonymous and stay in the private work area.
- Cluster-map TSV: cell-index first column and `cluster_L2_hybrid`; index must equal H5AD obs names in the same order. This is a frozen upstream annotation dependency.
- Config TSV: `name`, `clusters` (pipe-delimited), `n_top` (legacy field; no effect on core model), `min_cells_per_gene`, `min_cells_per_donor`. The supplied table is a production configuration, not participant data.
- Prepared outputs: for each config, genes × cells `<name>_counts.mtx`; `<name>_cell_metadata.tsv` containing sample, donor, cluster, boolean-as-text EBV state, cell key and library-size counts; `<name>_genes.tsv` containing gene, mean log expression and number of detected cells. Rows/columns and metadata order are aligned. Library size is computed before the preparation gene filter.
- In fit calls, effective gene-detection threshold is the maximum of the config threshold and runtime threshold. Invalid donor labels `doublet`, `unassigned`, `pool` and user-specified exclusions are removed before minimum-donor-cell filtering. Each EBV group must have ≥5 cells to fit.

## BCR

`--input_manifest` is a TSV with `donor` and `output_airr` (path to each donor's AIRR TSV). `--positive_donors` is a TSV with column `donor`, limited to donors eligible under the paper rule. AIRR rows must already represent productive heavy chains; this upstream restriction is not reconstructed by this core. Required AIRR columns include `sequence_id`, `sequence_alignment`, `germline_alignment`, `v_call`, `j_call`, `junction`, `junction_length`, `consensus_count`, `lineage_obs_name`, `cell_id`, `lineage_source_sample`, `lineage_cluster_L2`, `lineage_EBV_detected`; `locus` is IGH. Verify exact package contracts against the production schema before a real-data run. No real receptor sequences are included.

## Study meta-analysis

CSV columns: `Study`, `MS_pos`, `MS_neg`, `Ctrl_pos`, `Ctrl_neg`, with nonnegative integer counts and positive group totals. Pos/neg are event/non-event counts, not total sample sizes. Source behavior drops rows missing `MS_pos` and does not apply a separate Exclude column. Resolve study inclusion in the frozen input before execution. The example CSV is entirely synthetic.
