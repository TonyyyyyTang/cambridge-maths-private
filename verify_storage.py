#!/usr/bin/env python3
"""Fail on missing/corrupt/duplicate records; list unregistered stored files."""
import csv
import hashlib
from pathlib import Path
import sys

here=Path(__file__).resolve().parent
rows=list(csv.DictReader((here/'manifest.csv').open(newline='')))
seen=set(); errors=[]
for row in rows:
    path=row['local_path']
    if path in seen: errors.append('DUPLICATE '+path)
    seen.add(path)
    if Path(path).is_absolute() or '..' in Path(path).parts or not path.startswith(('raw/','my-handwritten/')):
        errors.append('UNSAFE '+path); continue
    file=here/path
    if not file.is_file(): errors.append('MISSING '+path)
    elif hashlib.sha256(file.read_bytes()).hexdigest()!=row['sha256']: errors.append('MISMATCH '+path)
if sys.argv[1:] == ['untracked']:
    for directory in ('raw','my-handwritten'):
        for file in (here/directory).rglob('*'):
            if file.is_file() and file.name!='.DS_Store' and file.relative_to(here).as_posix() not in seen:
                errors.append('UNTRACKED '+file.relative_to(here).as_posix())
print('\n'.join(errors) if errors else f'PASS: {len(rows)} archive records')
sys.exit(bool(errors))
