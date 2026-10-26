#!/usr/bin/env python3
"""Locate saved QC, optionally restore an attached cohort-only QC snapshot, and audit.

The snapshot is gzip/base64 CSV on stdin; it never enters Git. Run this with the
Python interpreter documented by upstream's TEDANA audit. No pip/conda changes.
"""
from __future__ import annotations

import argparse
import base64
import csv
import gzip
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import uuid
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[2]
FLAGS = [f'{t}_{kind}_run{r}_ok' for t in ('trust', 'ugr') for kind in ('beh', 'ev', 'mri') for r in (1, 2)]
# Full flags have the order "run1_full_ok", unlike the component flags above.
FLAGS += [f'{t}_run{r}_full_ok' for t in ('trust', 'ugr') for r in (1, 2)]


def read_csv(text):
    return list(csv.DictReader(io.StringIO(text.lstrip('\ufeff'))))


def cohort_view(rows, cohort):
    seen, view, included = set(), {}, set()
    for row in rows:
        sub = str(int(float(row['sub_id'])))
        if sub in seen:
            raise ValueError('Duplicate QC participant ID')
        seen.add(sub)
        if row.get('include_final', '').strip().upper() in ('TRUE','1'):
            included.add(sub)
        if sub not in cohort:
            continue
        if row.get('include_final', '').strip().upper() not in ('TRUE', '1'):
            raise ValueError('Manuscript participant not marked include_final')
        values = {}
        for key in FLAGS:
            value = row.get(key, '').strip().upper()
            if value not in ('TRUE', 'FALSE', '1', '0'):
                raise ValueError(f'Missing/unknown saved QC flag: {key}')
            values[key] = value in ('TRUE', '1')
        view[sub] = values
    if included != cohort:
        raise ValueError('Corrected include_final identities must match manuscript N=132 exactly')
    if set(view) != cohort:
        raise ValueError('QC table does not cover the exact manuscript cohort')
    return view


def candidates(roots):
    roots = [str(Path(r).resolve()) for r in roots if Path(r).is_dir()]
    if shutil.which('rg'):
        result = subprocess.run(['rg', '--files', '--hidden', '--no-ignore',
            '-g', 'master_subject_table_full*.csv', '-g', '!**/.git/**', *roots],
            text=True, capture_output=True)
        if result.returncode not in (0, 1):
            print('QC search notices:', result.stderr.strip(), flush=True)
        return sorted({Path(line) for line in result.stdout.splitlines()})
    found = []
    for root in roots:
        for directory, dirs, files in os.walk(root):
            dirs[:] = [d for d in dirs if d not in ('.git', 'node_modules', '.venv')]
            found.extend(Path(directory)/name for name in files if name.startswith('master_subject_table_full') and name.endswith('.csv'))
    return sorted(set(found))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--qc', type=Path, help='Use an explicit existing full QC CSV')
    parser.add_argument('--qc-snapshot-stdin', action='store_true', help='Read gzip/base64 cohort-only QC CSV from stdin')
    parser.add_argument('--upstream-root', type=Path, default=Path('/ZPOOL/data/projects/rf1-sra-linux2'))
    args = parser.parse_args()
    # Import only after parsing --help, then verify the actual chosen environment.
    import insula_ica as pipeline
    pipeline.software()
    print('Python:', sys.executable, sys.version.split()[0], flush=True)
    print('FSL:', pipeline.software()['fsl'], flush=True)
    cohort = set((ROOT / 'code/sublist_n132.txt').read_text().split())
    if len(cohort) != 132:
        raise ValueError('Expected the unchanged 132-participant manuscript list')
    snapshot_text, expected = None, None
    if args.qc_snapshot_stdin:
        snapshot_text = gzip.decompress(base64.b64decode(sys.stdin.read(), validate=False)).decode('utf-8')
        expected = cohort_view(read_csv(snapshot_text), cohort)
    paths = [args.qc] if args.qc else candidates([ROOT, args.upstream_root, Path.home() / 'Downloads'])
    matches = []
    for path in paths:
        try:
            rows = read_csv(path.read_text(encoding='utf-8-sig'))
            view = cohort_view(rows, cohort)
            if expected is not None and view != expected:
                print('Skipped QC with flags differing from the attachment:', path, flush=True)
                continue
            matches.append((path, rows, view))
        except (OSError, ValueError, KeyError) as exc:
            print(f'Skipped QC {path}: {exc}', flush=True)
    # Without an attached reference, require agreeing flags across candidate files.
    if len(matches) > 1 and any(m[2] != matches[0][2] for m in matches[1:]):
        raise ValueError('Candidate QC files disagree; specify --qc PATH')
    if not matches and (snapshot_text is None or args.qc):
        raise ValueError('No matching QC CSV found; supply --qc or the attached cohort snapshot')
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid.uuid4().hex[:6]
    out = ROOT / 'code/insula_ica/output' / f'linux2-audit-{stamp}'
    out.mkdir(parents=True)
    if matches:
        path, rows, view = matches[0]
        qc = path.resolve()
        print('Using matching existing QC:', qc, flush=True)
        print('Whole table rows:', len(rows), '; include_final TRUE:', sum(
            r.get('include_final', '').strip().upper() in ('TRUE', '1') for r in rows), flush=True)
    else:
        qc = out / 'inputs/master_subject_table_full.csv'
        qc.parent.mkdir()
        qc.write_text(snapshot_text)
        view = expected
        print('No matching linux2 copy found; restored attachment-derived N=132 QC snapshot:', qc, flush=True)
        print('Snapshot preserves original run QC flags; ages and noncohort rows are omitted.', flush=True)
    print('Analysis cohort: 132 (unchanged code/sublist_n132.txt)', flush=True)
    for task in ('trust', 'ugr'):
        count = sum(values[f'{task}_run{r}_full_ok'] for values in view.values() for r in (1, 2))
        print(f'Original QC-eligible {task} runs: {count} (not yet verified image availability)', flush=True)
    c = json.loads((ROOT / 'code/insula_ica/config.example.json').read_text())
    c.update(master_qc=str(qc), output_root=str(out), upstream_root=str(args.upstream_root.resolve()),
             derivatives_root=str(ROOT/'derivatives/insula_ica'/f'feat-residuals-{stamp}'))
    config_path = out / 'config.json'
    config_path.write_text(json.dumps(c, indent=2) + '\n')
    print('Config:', config_path, flush=True)
    result = subprocess.run([sys.executable, str(ROOT / 'code/insula_ica/insula_ica.py'), '--config', str(config_path), 'audit'])
    report = out / 'audit/report.md'
    if report.is_file():
        print('\n' + report.read_text(), flush=True)
        manifest = pipeline.read_table(out / 'audit/manifest.tsv')
        blocked = [r for r in manifest if r['status'] == 'blocked']
        for row in blocked[:10]:
            print('BLOCKED RUN:', row['key'], row['reason'], flush=True)
        print('\nSend back this output plus:', report, out / 'audit/manifest.tsv', flush=True)
    if result.returncode == 2 and report.exists():
        print('FEAT residual inventory completed with blockers; no FEAT rerun or production analysis was started.', flush=True)
        return 0
    return result.returncode


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (ValueError, KeyError, OSError, ImportError) as exc:
        print(f'BLOCKED: {exc}', file=sys.stderr)
        raise SystemExit(2)
