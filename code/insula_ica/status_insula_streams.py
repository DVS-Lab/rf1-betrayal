#!/usr/bin/env python3
"""Read-only, bounded status for baseline ICA, sm-0 FEAT and task-preserved ICA."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import subprocess

from status_group_ica import log_errors, short


ROOT = Path(__file__).resolve().parents[2]
STREAMS = (
    ('Baseline sm-5 ICA', 'linux2-paired-group-ica-*',
     ('run_group_ica.py', 'resume_group_ica.py'), 'paired-group-ica', 'rf1-paired-group-ica.*'),
    ('Unsmoothed sm-0 FEAT', 'linux2-sm-0-feat-*',
     ('run_unsmoothed_feat.py', 'recover_ugr_unsmoothed.py'), 'sm-0-feat', 'rf1-sm0-feat.*'),
    ('Task-preserved sensitivity', 'linux2-task-preserved-ica-*',
     ('task_preserved_ica.py',), 'task-preserved-ica', 'rf1-task-preserved-ica.*'),
)


def read(path):
    return json.loads(path.read_text())


def processes():
    try:
        output = subprocess.check_output(
            ['ps', '-u', str(os.getuid()), '-o', 'pid=,ppid=,stat=,etime=,args='],
            text=True, stderr=subprocess.DEVNULL, timeout=10)
    except (OSError, subprocess.SubprocessError):
        return None
    rows = []
    for line in output.splitlines():
        parts = line.split(None, 4)
        if len(parts) == 5 and not parts[2].startswith('Z'):
            rows.append(dict(pid=int(parts[0]), ppid=int(parts[1]), status=parts[2],
                             elapsed=parts[3], args=parts[4]))
    return rows


def pipeline(rows, scripts, root):
    runners = [p for p in rows if Path(p['args'].split()[0]).name.startswith('python')
               and any(re.search(r'(?:^|\s)(?:\S*/)?' + re.escape(s) + r'(?:\s|$)',
                                 p['args']) for s in scripts)]
    ids = {p['pid'] for p in runners}
    while True:
        expanded = ids | {p['pid'] for p in rows if p['ppid'] in ids}
        if expanded == ids:
            break
        ids = expanded
    active = [p for p in rows if p['pid'] in ids or (root and str(root) in p['args'])]
    counts = Counter()
    for item in active:
        argv = item['args'].split()
        exe = Path(argv[0]).name
        if exe in ('sh', 'bash', 'csh', 'tcsh', 'fslWish', 'wish', 'tclsh') and len(argv) > 1:
            exe = Path(argv[1]).name
        counts[exe] += 1
    return runners, active, {name: counts[name] for name in ('feat', 'film_gls', 'melodic')}


def launcher(lock, pattern):
    pointer = ROOT / 'code/insula_ica/output' / (lock + '-launch.lock') / 'log'
    if pointer.is_file():
        path = Path(pointer.read_text().strip())
        if path.is_file():
            return path
    paths = [p for p in Path('/tmp').glob(pattern)
             if p.is_file() and p.stat().st_uid == os.getuid()]
    return max(paths, key=lambda p: p.stat().st_mtime_ns) if paths else None


def fit_files(fits, errors):
    """Count saved success markers/maps; no image loading or large file hashing."""
    present = []
    for key, item in fits.items():
        if item.get('state') != 'complete':
            continue
        try:
            maps = Path(item['oriented_maps'])
            marker = read(maps.parent / 'complete.json')
            if marker['oriented_maps'] != str(maps) or not maps.is_file():
                raise ValueError('completion marker/map missing or inconsistent')
            present.append(key)
        except (OSError, ValueError, KeyError) as exc:
            errors.append(f'{key}: saved completion check: {exc}')
    return present


def brief_error(error):
    if 'No convergence after five deterministic starts' in error:
        return error.split(': No convergence')[0] + ': no convergence in 5 starts × 10,000 iterations'
    if 'returned non-zero exit status' in error:
        key = error.split(': Command ')[0]
        code = re.search(r'returned non-zero exit status (\d+)', error)
        return key + ': MELODIC exit ' + (code[1] if code else '?')
    return short(error)


def inspect(spec, rows):
    title, pattern, scripts, lock, log_pattern = spec
    records = sorted(p for p in (ROOT / 'code/insula_ica/records').glob(pattern)
                     if (p / 'state.json').is_file())
    record = records[-1] if records else None
    state = read(record / 'state.json') if record else {}
    root = Path(state['output_root']) if state.get('output_root') else None
    runners, active, counts = pipeline(rows or [], scripts, root)
    errors = []
    for collection in ('runs', 'fits', 'native_variance_normalization_pilot'):
        for key, item in state.get(collection, {}).items():
            if item.get('state') == 'failed':
                errors.append(f'{key}: {item.get("error", "failed")}')
    if state.get('error'):
        errors.append('Runner: ' + state['error'])
    publication = record / 'publication-error.txt' if record else None
    if publication and publication.is_file():
        errors.append('Git publication: ' + publication.read_text())
    log = launcher(lock, log_pattern)
    tail_errors = log_errors(log)
    run_done = [k for k, v in state.get('runs', {}).items() if v.get('state') == 'complete']
    fit_done = [k for k, v in state.get('fits', {}).items() if v.get('state') == 'complete']
    lines = [title, 'Record: ' + (record.name if record else 'none found'),
             'Runner: ' + ('unknown (process inspection unavailable)' if rows is None else
                          ', '.join(f'PID {p["pid"]}, elapsed {p["elapsed"]}' for p in runners)
                          if runners else 'not running') + '; state: ' + state.get('state', 'not recorded')]
    # Separate reported progress from inexpensive on-disk completion checks.
    disk_available = bool(root and root.is_dir())
    disk = {}
    if lock == 'sm-0-feat':
        plan = read(record / 'input-plan.json') if record and (record / 'input-plan.json').is_file() else []
        task_done = Counter(k.split('_task-')[1].split('_')[0] for k in run_done)
        for item in plan:
            if item.get('status') == 'blocked':
                errors.append(item['key'] + ': ' + item.get('error', 'blocked'))
        missing = [p for p in plan if p.get('status') == 'missing_source_inputs']
        if disk_available:
            disk['FEAT'] = 0
            for key in run_done:
                try:
                    marker = read(root / 'feat' / (key + '_sm-0.feat') / 'sm0_complete.json')
                    qc = read(record / 'qc' / (key + '.json'))
                    outputs = marker['outputs']
                    if marker['qc'] != qc or not outputs or not all(Path(p['path']).is_file() for p in outputs):
                        raise ValueError('QC or required output missing/inconsistent')
                    disk['FEAT'] += 1
                except (OSError, ValueError, KeyError) as exc:
                    errors.append(f'{key}: saved FEAT completion check: {exc}')
        lines += [f'FEAT complete: {len(run_done)}/490 (Trust {task_done["trust"]}/241; UGR {task_done["ugr"]}/249); ICA: not launched by this stream',
                  'Input audit: ' + '; '.join(f'{k} {v}' for k, v in state.get('input_status_counts', {}).items())]
        if missing:
            lines.append(f'Blocked by missing original inputs: {len(missing)} runs; exact paths in input-plan.json')
    else:
        if disk_available:
            disk['prepared markers'] = sum((root / 'run_arrays' / (k + '_complete.json')).is_file() for k in run_done)
            disk['ICA fits'] = len(fit_files(state.get('fits', {}), errors))
        expected = 24 if lock == 'paired-group-ica' else 6
        lines.append(f'Inputs prepared: {len(run_done)}/490; ICA complete: {len(fit_done)}/{expected}; new FEATs: 0 (reuses sm-5 inputs)')
        if lock == 'paired-group-ica':
            auto = sum(k.endswith('/dim-0') for k in fit_done)
            pilot = state.get('native_variance_normalization_pilot', {})
            lines.append(f'Fixed ICA: {len(fit_done)-auto}/18; automatic: {auto}/6; separate VN pilot: '
                         f'{sum(v.get("state")=="complete" for v in pilot.values())}/2')
            if disk_available:
                disk['VN fits'] = len(fit_files(pilot, errors))
    lines.append('Saved completion checks: ' + ('; '.join(f'{k} {v}' for k, v in disk.items())
                 if disk_available else 'output root unavailable here; counts above are recorded progress'))
    if rows is not None:
        lines.append(f'Active FSL processes: FEAT {counts["feat"]}; FILM {counts["film_gls"]}; MELODIC {counts["melodic"]}')
    lines.append(f'Recorded/check errors: {len(errors)}' + (' (first 5 below)' if errors else ''))
    lines.extend('  ' + brief_error(e) for e in errors[:5])
    lines.append('Launcher errors, bounded tail: ' + (' | '.join(tail_errors) if tail_errors else 'none found' if log else 'log unavailable'))
    return lines, dict(record=str(record) if record else None, state=state, output_checks=disk,
                       errors=errors, runners=runners, active_pipeline_processes=active,
                       launcher_log=str(log) if log else None, launcher_errors=tail_errors)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--brief', action='store_true', help='Compact output (also the default)')
    parser.parse_args()
    rows = processes()
    results = {}
    for spec in STREAMS:
        try:
            lines, snapshot = inspect(spec, rows)
        except (OSError, ValueError, KeyError) as exc:
            lines, snapshot = [spec[0], 'Status inspection error: ' + short(exc)], dict(error=str(exc))
        print('\n'.join(lines) + '\n')
        results[spec[0]] = snapshot
    saved = ROOT / 'code/insula_ica/output' / ('all-streams-status-' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
    saved.mkdir(parents=True)
    (saved / 'snapshot.json').write_text(json.dumps(results, indent=2) + '\n')
    print('Full snapshot (including all errors and log paths): ' + str(saved / 'snapshot.json'))
    print('Status checks only; no jobs started or restarted.')


if __name__ == '__main__':
    main()
