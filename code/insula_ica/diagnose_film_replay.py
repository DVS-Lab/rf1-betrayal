#!/usr/bin/env python3
"""Separate statistical-fit differences from the Trust preprocessing rebuild.

Replay one original UGR FILM command on its retained filtered data/full design,
changing only the results destination and resolving the current FILM binary.
Never run fMRIPrep, TEDANA, FEAT prestats, or the cohort reconstruction here.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import uuid

import nibabel as nib
import numpy as np

import insula_ica as core
import feat_residuals as residuals
from rebuild_trust import compare_fit, require, save_json


def agreement(a, b):
    valid = np.isfinite(a) & np.isfinite(b)
    a, b = np.asarray(a)[valid], np.asarray(b)[valid]
    require(len(a) > 0, 'No finite common voxels')
    denominator = float(np.linalg.norm(a))
    return dict(voxels=len(a),
                spatial_r=float(np.corrcoef(a, b)[0, 1]) if a.std() > 0 and b.std() > 0 else None,
                relative_rms_difference=float(np.linalg.norm(b - a) / denominator) if denominator else None,
                max_absolute_difference=float(np.max(np.abs(b - a))))


def diagnose(audit, subject, run, output, execute):
    output = output.resolve()
    with (audit / 'manifest.tsv').open() as stream:
        rows = [row for row in csv.DictReader(stream, delimiter='\t')
                if row['task'] == 'ugr' and row['subject'].removeprefix('sub-') == subject
                and int(row['run']) == run and row['status'] == 'included']
    require(len(rows) == 1, 'Expected one included UGR run with a verified replay route')
    row = rows[0]
    require(row['recovery_mode'] == 'replay_film', 'Diagnostic requires the retained-filtered-data route')
    original = Path(row['original_feat']).resolve()
    require(original.is_relative_to((core.ROOT / 'derivatives/fsl').resolve()), 'Unexpected historical FEAT root')
    require(row['prewhiten'] in ('0', '1'), 'Unknown original prewhitening setting')
    filtered = original / 'filtered_func_data.nii.gz'
    image = nib.load(filtered)
    require(len(image.shape) == 4, 'Retained filtered data are not 4D')
    nvol = image.shape[3]
    require(len(core.fsl_matrix(original / 'design.mat')) == nvol, 'Original data/design length mismatch')
    target = output / (row['key'] + '.feat')
    cmd, log_source = residuals.replay_command(original, target, row['prewhiten'] == '1')
    require(shutil.which(cmd[0]) is not None, 'film_gls is unavailable; source the normal FSL setup')
    require(output.is_relative_to((core.ROOT / 'derivatives/insula_ica').resolve()), 'Use an isolated diagnostic root')
    require(not output.exists(), 'Diagnostic output exists; preserve it and use a fresh destination')
    output.mkdir(parents=True)
    paths = [filtered, original / 'design.mat', original / 'design.con', original / 'design.fts',
             original / 'mask.nii.gz', original / 'mean_func.nii.gz', log_source,
             *sorted((original / 'stats').glob('*.nii.gz'))]
    before = residuals.files_identity(paths)
    version_file = Path(os.environ.get('FSLDIR', '/unavailable')) / 'etc/fslversion'
    plan = dict(key=row['key'], original_feat=str(original), argv=cmd, cwd=str(original),
                output=str(target), prewhiten=row['prewhiten'], nvol=nvol,
                original_inputs=before, film_binary_sha256=core.sha(cmd[0]),
                fsl_version=version_file.read_text().strip() if version_file.is_file() else None)
    save_json(output / 'plan.json', plan)
    print('Diagnostic output:', output, flush=True)
    print('Exact same-input replay:', shlex.join(cmd), flush=True)
    if not execute:
        print('DRY RUN; add --execute to run this one statistical fit.', flush=True)
        return
    target.mkdir()
    for name in ('design.mat', 'design.con', 'design.fts', 'mask.nii.gz', 'mean_func.nii.gz'):
        source = original / name
        if source.is_file():
            shutil.copy2(source, target / name)
    try:
        with (output / 'film.log').open('w') as stream:
            subprocess.run(cmd, cwd=original, stdout=stream, stderr=subprocess.STDOUT, check=True)
    finally:
        after = residuals.files_identity(paths)
        require(before == after, 'Historical inputs changed during diagnostic; inspect before continuing')
    issue = None
    try:
        compare_fit(dict(key=row['key'], original_feat=str(original), nvol=nvol), target)
    except ValueError as exc:
        issue = str(exc)
    comparison = json.loads((target / 'historical_comparison.json').read_text())
    mask_image = nib.load(original / 'mask.nii.gz')
    mask = mask_image.get_fdata() > 0
    metrics = {}
    task_contrast = int(row['contrast_index'])
    for name in (f'cope{task_contrast}.nii.gz', f'zstat{task_contrast}.nii.gz'):
        old_image, new_image = nib.load(original / 'stats' / name), nib.load(target / 'stats' / name)
        require(core.same_grid(old_image, new_image) and core.same_grid(old_image, mask_image), 'Statistic grid differs')
        metrics[name] = agreement(old_image.get_fdata()[mask], new_image.get_fdata()[mask])
    passed = issue is None
    result = dict(key=row['key'], same_input_film_passed=passed, comparison_error=issue,
                  original_inputs_unchanged=True, task_contrast_metrics=metrics,
                  maps_matched=sum(d['matched'] for d in comparison['comparisons']),
                  maps_compared=len(comparison['comparisons']), fsl_version=plan['fsl_version'])
    save_json(output / 'result.json', result)
    print(json.dumps(result, indent=2), flush=True)
    print('Interpretation:', 'This retained-data fit matches; investigate Trust BOLD/prestats or Trust-specific settings next.'
          if passed else 'Current replay does not reproduce this retained-data fit; investigate fitting/runtime or retained-input provenance first.', flush=True)
    print('This diagnostic does not unlock or launch the production rebuild.', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit-dir', type=Path, required=True)
    parser.add_argument('--subject', default='10369')
    parser.add_argument('--run', type=int, default=1, choices=(1, 2))
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid.uuid4().hex[:6]
    output = core.ROOT / 'derivatives/insula_ica/diagnostics' / ('ugr-film-' + stamp)
    diagnose(args.audit_dir.resolve(), args.subject.removeprefix('sub-'), args.run, output, args.execute)


if __name__ == '__main__':
    main()
