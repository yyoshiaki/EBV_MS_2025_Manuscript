# Computational coverage and release gates

Updated 2026-10-15 EDT. Mapping uses the 14 September 2026 manuscript and supplementary materials, not old figure numbers in production filenames.

| Main result | Retained core | Remaining dependency or gap |
| --- | --- | --- |
| Fig. 1 genetic architecture; Fig. 2 HLA | None yet | Final cohort QC/REGENIE/MR-MEGA/METAL, HLA/conditional association and PRS code has not been verified/exported. Older AoU-only launchers were excluded because they do not establish the final cross-cohort analysis. |
| Fig. 3C, study meta-analysis | `261015_meta_ebv_dna.R` | Frozen study-count input and its inclusion decision; adjusted cohort qPCR regression is not included. |
| Fig. 4 B-cell clonality | `261015_define_bcr_clones.R` | Upstream productive/donor-split AIRR preparation and downstream repertoire/network comparisons remain outside the retained core. |
| Fig. 5A–C, subset and pooled-memory programs | preparation → NEBULA → subset GSEA | Frozen counts/normalization, PCR-union calls and hybrid label mapping; selected gene panels/figure composition excluded. Pooled-memory GSEA is not supplied by this subset GSEA script. |
| Fig. 5D, scDRS | External tool dependency only | Final MAGIC/scDRS inputs, covariates, gene scores and cohort-specific test invocation still need verification. |
| Supplementary sensitivity | Same preparation/NEBULA entry points with explicit exclusion | Author supplies the frozen exclusion privately; exact numerical equivalence to final sensitivity is not claimed. External SLE/IM-HLH and MS–HC computations are not included in this minimal pass. |

These are real production-derived calculations, but the repository is **not paper-wide reproducibility-complete**. No large recomputation was run, and synthetic tests do not validate real-data biological conclusions. The privacy checks cover repository files, not unavailable production datasets.

Before public release: obtain and verify the missing core genetic paths; settle remaining upstream input provenance and permitted data access; compare against frozen paper aggregate outputs; confirm software/database versions and code ownership/license; confirm final manuscript citation; obtain author publication approval. Public release/Zenodo/DOI are intentionally not performed.
