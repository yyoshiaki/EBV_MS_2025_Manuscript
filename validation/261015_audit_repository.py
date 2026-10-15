"""Fail if the review candidate contains data payloads or obvious private material."""
from pathlib import Path
import re

repo = Path(__file__).resolve().parents[1]
patterns = {
    'private host/path': r'/Users/|/home/|/vast/|/gpfs/',
    'credential': r'-----BEGIN .*PRIVATE KEY|\bsk-(?:proj-)?[A-Za-z0-9_-]{20,}|\bgh[pousr]_[A-Za-z0-9]{20,}',
    'legacy participant label': r'\b(?:MS|DRE|HC)[0-9]{3,}\b',
    'long receptor/nucleotide sequence': r'\b[ACGTN]{40,}\b',
}
failures = []
files = []
for p in sorted(repo.rglob('*')):
    if not p.is_file() or '.git' in p.relative_to(repo).parts or '__pycache__' in p.parts:
        continue
    rel = p.relative_to(repo)
    if str(rel).split('/')[0] in {'work', 'outputs', 'results', 'data'}:
        failures.append(f'{rel}: generated/private-data directory')
        continue
    if p.suffix not in {'.md', '.py', '.R', '.tsv', '.csv'} and p.name != '.gitignore':
        failures.append(f'{rel}: unexpected file type')
        continue
    # Pattern definitions necessarily contain the literals they are checking.
    if p.resolve() != Path(__file__).resolve():
        text = p.read_text()
        for label, pattern in patterns.items():
            if re.search(pattern, text):
                failures.append(f'{rel}: {label}')
    files.append(rel)
assert not failures, '\n'.join(failures)
print(f'PASS: {len(files)} text files; no prohibited payload types, internal host paths, tested credential patterns, legacy participant IDs or long nucleotide sequences.')
print('This bounded scan is not proof of complete de-identification or permission to publish.')
