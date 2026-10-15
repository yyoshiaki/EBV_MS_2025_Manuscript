# EBV and multiple sclerosis analysis core

Code supporting *Genetic Determinants of Epstein–Barr Virus Persistence and Multiple Sclerosis Converge in Infected Memory B Cells* (Science Immunology manuscript aei0178).

**Private preparation branch; incomplete paper-wide reproduction. Updated 2026-10-15 EDT.** Five production-derived computational scripts are retained. Drawing, layout, exploratory notebooks, old versions and upstream tool implementations are excluded. No participant data or controlled-access files are included. No release, DOI or public archive has been created.

## Core and execution order

| Step | Script | Purpose |
| --- | --- | --- |
| 1 | `scripts/261015_prepare_nebula.py` | Frozen annotated counts and EBV-union status → ordered subset/pool inputs |
| 2 | `scripts/261015_fit_nebula.R` | Donor-aware EBV+/− differential expression with corrected library-size scaling |
| 3 | `scripts/261015_enrich_pathways.R` | Reactome/GO BP GSEA from subset effect estimates |
| Independent | `scripts/261015_define_bcr_clones.R` | Per-donor productive IGH clone calling |
| Independent | `scripts/261015_meta_ebv_dna.R` | Study-level MS versus control EBV DNA positivity meta-analysis |

Run from the repository root with an existing environment from `environments/README.md`. All inputs and generated outputs belong in an external private work directory. For example:

```bash
python scripts/261015_prepare_nebula.py   --input-h5ad "$INPUT_H5AD" --cluster-map "$CLUSTER_MAP"   --config config/260318_pbmc_bcell_nebula_configs.tsv   --output-root "$WORK_DIR/prepared" --analysis-name primary   --configs Bnaive,Bnaive_MX1,Busm,Bmem_int,Bsm,ABC,ASC,Bcell_memoryB
Rscript scripts/261015_fit_nebula.R   "$WORK_DIR/prepared/primary/nebula/inputs" "$WORK_DIR/fitted/primary"   config/260318_pbmc_bcell_nebula_configs.tsv   Bnaive,Bnaive_MX1,Busm,Bmem_int,Bsm,ABC,ASC,Bcell_memoryB 50 4
Rscript scripts/261015_enrich_pathways.R   "$WORK_DIR/prepared/primary/nebula/inputs" "$WORK_DIR/fitted/primary"   config/260318_pbmc_bcell_nebula_configs.tsv   Bnaive,Bnaive_MX1,Busm,Bmem_int,Bsm,ABC,ASC
```

Donor exclusion sensitivity reuses steps 1–2 with `--exclude-donors "$EXCLUDED_DONOR"` and a distinct `--analysis-name sensitivity`, prepared output path and fit output directory. Do not reuse the primary gene metadata for a sensitivity GSEA. The extracted GSEA preserves the production exclusion of `Bcell_` pooled configurations; its scope is subset enrichment, not pooled-memory enrichment.

BCR inputs must already be split by donor and restricted to productive IGH:

```bash
Rscript scripts/261015_define_bcr_clones.R   --input_manifest "$AIRR_MANIFEST" --positive_donors "$ELIGIBLE_DONORS"   --output_root "$WORK_DIR/bcr" --nproc 1
Rscript scripts/261015_meta_ebv_dna.R "$STUDY_COUNTS" "$WORK_DIR/meta"
```

## Statistical settings

NEBULA retains `~ EBV_detected + sample`, donor random effect, raw positive library-size scaling factors (no pre-log), `LN`, `NBGMM`, `cpc=0.005`, config-specific preparation filters, and runtime detection threshold 50. BH adjustment is within each tested subset/pool. The `deg_*.tsv` files contain **all tested genes**, with `strict_deg` marking the original strict predicate; strict DEGs require `pvals_adj < 0.05`, `abs(logfoldchanges) > 0.25` and `mean > 0.5`. Fold changes are natural-log values. Memory pooling is defined in the frozen config table.

GSEA retains mean >0.5, ranking by log fold change across eligible genes, SYMBOL→ENTREZ mapping with maximum absolute effect for duplicate Entrez IDs, human Reactome/GO BP, gene-set sizes 10–500, seed 1 and FDR <0.05. Genuine enrichment exceptions now stop execution. Meta-analysis retains OR, DerSimonian–Laird heterogeneity and Hartung–Knapp intervals. No new statistical method is introduced.

## Inputs, provenance and limits

See [input schema](docs/261015_input_schema.md), [coverage and remaining gates](docs/261015_coverage.md), [source hashes and adaptations](docs/261015_provenance.tsv), and [validation](validation/README.md). These scripts start from frozen analysis inputs; they do not reconstruct those inputs from raw reads. Missing GWAS/HLA/PRS paths and upstream annotation are explicit completion gates. The original scripts remain unchanged on the production host.

Upstream software is referenced in [dependencies and citation](docs/261015_dependencies.md), without vendoring. License selection and the final archival citation remain pending author confirmation. This branch must remain private until those gates and paper-wide coverage are resolved and the author approves publication.
