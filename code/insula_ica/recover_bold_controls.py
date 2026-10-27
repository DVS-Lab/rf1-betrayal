#!/usr/bin/env python3
"""Collect completed controls without launching any preprocessing or model fits.

Preserve the original failed collection record. Write a new recovery subfolder
with comparisons of the already completed output images.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import traceback
import uuid

import nibabel as nib
import numpy as np

from compare_bold_controls import find_bold
from rebuild_trust import ROOT, require, sha
from run_bold_controls import compare, report, save, verify_snapshot


def completed_case(item, record):
    log = record / (item['name'] + '.fmriprep-log.txt')
    require(log.is_file(), 'Missing preprocessing log: ' + str(log))
    require('fMRIPrep finished successfully!' in log.read_text(errors='replace'),
            'Preprocessing success marker absent: ' + str(log))
    outputs = {}
    for row in item['runs']:
        bold = find_bold(Path(item['case']), row)
        image = nib.load(bold)
        require(len(image.shape) == 4 and image.shape[-1] == row['nvol'], 'BOLD time count differs: ' + str(bold))
        if 'tr' in row:
            require(np.isclose(image.header.get_zooms()[3], row['tr'], atol=1e-4), 'BOLD TR differs: ' + str(bold))
        outputs[row['key']] = str(bold)
    return dict(name=item['name'], returncode=0, preprocessing_success_marker_verified=True,
                preprocessing_log_sha256=sha(log), bold_outputs=outputs)


def recover(record):
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid.uuid4().hex[:6]
    output = record / ('recovered-' + stamp)
    output.mkdir(exist_ok=False)
    status = dict(state='collecting_existing_outputs', cases={}, input_changes=[],
                  source_record=str(record), preprocessing_launched=False)
    save(output / 'status.json', status)
    try:
        plan = json.loads((record / 'plan.json').read_text())
        cases = plan['cases']
        rows = {r['key']: r for case in cases for r in case['runs']}
        snapshots = json.loads((record / 'input-sha256.json').read_text())
        print('Checking retained inputs; no preprocessing will run.', flush=True)
        status['input_changes'] = verify_snapshot(snapshots)
        status['input_changes_at_collection_start'] = list(status['input_changes'])
        # A later raw-input edit does not prevent reading completed BOLD outputs.
        # Retain and flag drift rather than discarding otherwise useful metrics.
        before = {}
        for case in cases:
            try:
                result = completed_case(case, record)
                for path in result['bold_outputs'].values():
                    stat = Path(path).stat()
                    before[path] = dict(size=stat.st_size, mtime_ns=stat.st_mtime_ns)
                status['cases'][case['name']] = result
                print(case['name'] + ': completed outputs found', flush=True)
            except Exception as exc:
                status['cases'][case['name']] = dict(name=case['name'], error=str(exc))
                print(case['name'] + ': ' + str(exc), flush=True)
        save(output / 'status.json', status)
        print('Comparing saved BOLD frames and full time series.', flush=True)
        results = compare(list(rows.values()), cases, status)
        save(output / 'results.json', results)
        identities = {}
        for path, previous in before.items():
            stat = Path(path).stat()
            require(previous == dict(size=stat.st_size, mtime_ns=stat.st_mtime_ns), 'BOLD changed during collection: ' + path)
            identities[path] = dict(**previous, sha256=sha(path))
        save(output / 'bold-output-sha256.json', identities)
        status['input_changes'] = verify_snapshot(snapshots)
        settings = json.loads((record / 'control-settings-checks.json').read_text())
        invalid_settings = any(c.get('saved_config_matches_controls') is False or c.get('repeat_seed_identity') is False for c in settings)
        errors = any('error' in case for case in status['cases'].values())
        incomplete = any(not m.get('compared') for row in results.values() for m in row['full_series'].values())
        status.update(state='complete_with_errors' if errors or incomplete or invalid_settings or status['input_changes'] else 'complete',
                      finished_utc=datetime.now(timezone.utc).isoformat())
        save(output / 'status.json', status)
        (output / 'report.md').write_text(report(results, status))
        save(output / 'collection-provenance.json', dict(
            source_record=str(record), original_record_untouched=True,
            source_file_sha256={name: sha(record / name) for name in
                ('plan.json', 'status.json', 'input-sha256.json', 'control-settings-checks.json')},
            recovery_script_sha256=sha(Path(__file__)),
            code_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()))
        (output / 'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in sorted(output.iterdir()) if p.is_file()))
        print('Recovery state:', status['state'], flush=True)
        print('Saved recovery:', output, flush=True)
        print('No fMRIPrep, FEAT, TEDANA, or Warpkit jobs were launched.', flush=True)
        return output
    except Exception as exc:
        status.update(state='failed', error=str(exc))
        save(output / 'status.json', status)
        (output / 'failure.txt').write_text(traceback.format_exc())
        print('Collection error:', exc, flush=True)
        print('Saved recovery:', output, flush=True)
        return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record', type=Path, required=True)
    args = parser.parse_args()
    record = args.record.resolve()
    require(record.is_relative_to((ROOT / 'code/insula_ica/records').resolve()), 'Expected an existing repository record')
    require((record / 'plan.json').is_file(), 'Missing completed control plan')
    recover(record)


if __name__ == '__main__':
    main()
