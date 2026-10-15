# Bounded validation

Verified 2026-10-15 17:07 EDT. This is software validation with synthetic inputs, not a real-data paper rerun.

| Check | Observed result |
| --- | --- |
| Python syntax / all four R scripts parse | PASS |
| Production numerical-block comparison | NEBULA invocation, Entrez/ranking reduction and hierarchical clone invocation retained verbatim |
| Actual preparation entry point, 240 synthetic cells / 6 fictitious donors | PASS: matrix/metadata ordering, raw library sizes, mean expression, EBV balance and donor exclusion |
| Invalid flags and typo config | PASS: text EBV flags rejected, including `python -O`; unknown config rejected |
| Actual NEBULA 1.5.6 entry point | PASS: 19 genes fitted; detected-in-49 gene excluded, 50/51 retained; engineered positive effect recovered; finite effects and BH adjustment checked |
| Model invocation | LN/NBGMM/cpc=0.005; 240-cell run reports no HL fallback. A preliminary 120-cell run auto-switched to HL because cells per donor <30, so the larger run was used to resolve that limitation. Requested ncore=2 was reduced by the runtime to available cores. |
| Meta-analysis synthetic three-study table | PASS: study and pooled log OR equal log(2.25) to 1e-12; exported random-model interval equals fitted HK interval |
| Repository privacy/payload scan | PASS for checked file types, internal host paths, credential patterns, legacy participant labels and long nucleotide strings |

The final Slurm synthetic job was `13414682`, COMPLETED, exit 0, elapsed 1m09s. It used isolated scratch and no participant data. Matrix emits a deprecation warning for the inherited sparse coercion; meta emits a deprecation warning for inherited `hakn`. These were preserved to avoid changing the production calculation during extraction.

To repeat preparation tests in an existing Python environment containing anndata/numpy/pandas/scipy:

```bash
python validation/261015_smoke_prepare.py --work-dir "$WORK_DIR/synthetic"
Rscript scripts/261015_fit_nebula.R \
  "$WORK_DIR/synthetic/prepared/primary/nebula/inputs" \
  "$WORK_DIR/synthetic/fitted/primary" \
  config/260318_pbmc_bcell_nebula_configs.tsv Bcell_memoryB 50 2
Rscript scripts/261015_meta_ebv_dna.R examples/261015_synthetic_meta.csv "$WORK_DIR/synthetic/meta"
python validation/261015_audit_repository.py
```

GSEA and BCR clone entry points were parsed and source-compared, but were **not executed end to end**. Their annotation/package versions and real input provenance remain gates. The repository scan is bounded and does not certify de-identification or publication permission.

Independent review found six items: four accepted (fail on missing GSEA inputs, export fitted meta confidence bounds, reject typo configs, retain flag checks under optimization); one rejected by direct R 4.4.1 execution (`True`/`False` remain character, not logical); one addressed by the actual NEBULA boundary test already in progress. Review logs are kept outside the repository.
