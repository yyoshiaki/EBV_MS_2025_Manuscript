"""Run the actual preparation entry point using entirely synthetic cells."""
import argparse
from pathlib import Path
import subprocess
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

parser = argparse.ArgumentParser()
parser.add_argument('--work-dir', required=True)
args = parser.parse_args()
work = Path(args.work_dir).resolve()
work.mkdir(parents=True, exist_ok=True)
repo = Path(__file__).resolve().parents[1]
rng = np.random.default_rng(1015)
n = 240
positive = np.arange(n) % 2 == 0
counts = rng.poisson(3, size=(n, 20)) + 1
counts[:, 0] = 0
counts[:49, 0] = 2
counts[:, 1] = 0
counts[:50, 1] = 2
counts[:, 2] = 0
counts[:51, 2] = 2
counts[positive, 3] *= 4
counts[positive, 4] *= 2
obs = pd.DataFrame({
    'sample': np.where(np.arange(n) < 120, 'SYNTH_SAMPLE_A', 'SYNTH_SAMPLE_B'),
    'donor_demux': [f'SYNTH_DONOR_{i // 40}' for i in range(n)],
    'cluster_L2': np.resize(['Busm', 'Bmem int', 'Bsm', 'Bsm FCER2', 'ABC'], n),
    'EBV_detected_union_pcr_v2': positive,
}, index=[f'SYNTH_CELL_{i:03d}' for i in range(n)])
genes = ['SYNTH_49', 'SYNTH_50', 'SYNTH_51', 'HLA-A', 'CD86'] + [f'SYNTH_GENE_{i}' for i in range(15)]
normalized = np.log1p(counts / counts.sum(axis=1)[:, None] * 1e5)
a = ad.AnnData(X=sparse.csr_matrix(normalized), obs=obs, var=pd.DataFrame(index=genes))
a.layers['counts'] = sparse.csr_matrix(counts)
a.write_h5ad(work / 'synthetic.h5ad')
pd.DataFrame({'cluster_L2_hybrid': obs.cluster_L2}, index=obs.index).to_csv(work / 'clusters.tsv', sep='\t')
cmd = [sys.executable, str(repo / 'scripts/261015_prepare_nebula.py'),
       '--input-h5ad', str(work / 'synthetic.h5ad'), '--cluster-map', str(work / 'clusters.tsv'),
       '--config', str(repo / 'config/260318_pbmc_bcell_nebula_configs.tsv'),
       '--output-root', str(work / 'prepared'), '--analysis-name', 'primary', '--configs', 'Bcell_memoryB']
subprocess.run(cmd, check=True)
p = work / 'prepared/primary/nebula/inputs'
m = pd.read_csv(p / 'Bcell_memoryB_cell_metadata.tsv', sep='\t', dtype={'EBV_detected': str})
g = pd.read_csv(p / 'Bcell_memoryB_genes.tsv', sep='\t')
x = mmread(p / 'Bcell_memoryB_counts.mtx').toarray()
expected_indices = obs.index.get_indexer(m.cell_key)
np.testing.assert_array_equal(x, counts[expected_indices].T)
np.testing.assert_array_equal(m.library_size_counts, counts.sum(axis=1)[expected_indices])
np.testing.assert_allclose(g['mean'], normalized.mean(axis=0))
assert g.n_cells_by_counts.tolist()[:3] == [49, 50, 51]
assert m.EBV_detected.eq('True').sum() == 120 and m.EBV_detected.eq('False').sum() == 120
assert m.donor_demux.nunique() == 6
assert m.cell_key.tolist() == m.sort_values(['donor_demux', 'sample', 'cell_key']).cell_key.tolist()
subprocess.run(cmd + ['--exclude-donors', 'SYNTH_DONOR_0', '--analysis-name', 'sensitivity'], check=True)
s = pd.read_csv(work / 'prepared/sensitivity/nebula/inputs/Bcell_memoryB_cell_metadata.tsv', sep='\t')
assert len(s) == 200 and not s.donor_demux.eq('SYNTH_DONOR_0').any()
# Demonstrate that text False cannot silently become True.
a.obs['EBV_detected_union_pcr_v2'] = a.obs['EBV_detected_union_pcr_v2'].astype(str)
a.write_h5ad(work / 'invalid_text_flags.h5ad')
bad = cmd.copy()
bad[bad.index('--input-h5ad') + 1] = str(work / 'invalid_text_flags.h5ad')
for invalid_command in [bad, [sys.executable, '-O', *bad[1:]]]:
    result = subprocess.run(invalid_command, capture_output=True, text=True)
    assert result.returncode != 0 and 'must be boolean' in result.stderr
result = subprocess.run(cmd + ['--configs', 'UNKNOWN_CONFIG'], capture_output=True, text=True)
assert result.returncode != 0 and 'Unknown requested config' in result.stderr
print('PASS: real entry point, matrix alignment, raw offsets, gene means, 49/50/51 boundary, status balance, donor exclusion, text-flag rejection.')
