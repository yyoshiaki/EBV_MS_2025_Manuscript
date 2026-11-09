# Sources: 250825_variantQC_LDpruning.ipynb cell 21 and
# 250825_regenie_step2.ipynb cell 11 (zero-based).
# Usage: python SCRIPT QC_ROOT OUTPUT.snplist.gz
import sys
from pathlib import Path
import pandas as pd

qc = Path(sys.argv[1])
for anc in ['eur', 'afr', 'amr', 'eas', 'sas', 'mid']:
    snps = pd.concat([pd.read_csv(qc / str(i) / f'chr{i}_{anc}_filtered.snplist',
                                 header=None) for i in range(1, 23)]).set_index(0)
    snps[anc] = 1
    if anc == 'eur':
        merged = snps.copy()
    else:
        merged = pd.merge(merged, snps, left_index=True, right_index=True, how='outer')
merged = merged.fillna(0).astype(int)
merged[merged.sum(axis=1) >= 4][[]].to_csv(sys.argv[2], header=False)
