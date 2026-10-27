#!/usr/bin/env python3
"""Read retained provenance and intermediate checkpoints for the Trust pilots.

No preprocessing, model fitting, source edits, or production launch. Only a
fresh diagnostic JSON is written inside the existing isolated rebuild root.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import difflib
import html
import json
from pathlib import Path
import re
import tomllib
import uuid

import nibabel as nib
import numpy as np

import insula_ica as core
from diagnose_film_replay import agreement
from rebuild_trust import save_json


def matrix_match(left, right):
    result = dict(left=str(left), right=str(right), available=left.is_file() and right.is_file())
    if not result['available']:
        result['equal'] = None
        return result
    try:
        a, b = np.loadtxt(left, ndmin=2), np.loadtxt(right, ndmin=2)
    except (ValueError, OSError) as exc:
        result.update(equal=None, error=str(exc))
        return result
    result.update(left_shape=list(a.shape), right_shape=list(b.shape),
                  equal=bool(a.shape == b.shape and np.array_equal(a, b)))
    if a.shape == b.shape:
        result['max_absolute_difference'] = float(np.max(np.abs(a - b)))
    return result


def tedana_check(directory, prefix, saved):
    metrics = directory / (prefix + '_desc-tedana_metrics.tsv')
    mixing = directory / (prefix + '_desc-ICA_mixing.tsv')
    result = dict(metrics=str(metrics), mixing=str(mixing),
                  available=metrics.is_file() and mixing.is_file(), rejected_suffix_matches=None)
    description = directory / 'dataset_description.json'
    if description.is_file():
        metadata = json.loads(description.read_text())
        result['dataset_generated_by'] = metadata.get('GeneratedBy', [])
    if not result['available']:
        result['status'] = 'Component identities cannot be checked from these missing source files.'
        return result
    try:
        with metrics.open(encoding='utf-8-sig') as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        rejected = [int(re.fullmatch(r'ICA_(\d+)', row['Component'].strip()).group(1))
                    for row in rows if row['classification'].strip().lower() == 'rejected']
        result['rejected_components'] = rejected
        result['n_rejected'] = len(rejected)
        if not rejected:
            result['status'] = 'Source metrics list zero rejected components.'
            return result
        timecourses = np.loadtxt(mixing, delimiter='\t', skiprows=1, ndmin=2)[:, rejected]
        combined = np.loadtxt(saved, ndmin=2)
        match = combined.shape[0] == timecourses.shape[0] and combined.shape[1] >= len(rejected)
        match = match and np.allclose(combined[:, -len(rejected):], timecourses, atol=1e-6, rtol=1e-6)
        result['rejected_suffix_matches'] = bool(match)
        result['status'] = 'Checked positional rejected-component selection and suffix ordering from the archived generator.'
    except (ValueError, KeyError, IndexError, AttributeError) as exc:
        result['status'] = 'Source verification failed: ' + str(exc)
    return result


def provenance(directory, subject):
    result = dict(root=str(directory), dataset_generated_by=[], execution_configs=[],
                  subject_html_versions=[], html_software_versions={})
    description = directory / 'dataset_description.json'
    if description.is_file():
        data = json.loads(description.read_text())
        result['dataset_generated_by'] = data.get('GeneratedBy', [data['PipelineDescription']] if 'PipelineDescription' in data else [])
    candidates = set()
    for base in (directory / subject / 'log', directory / subject / 'ses-01/log'):
        candidates.update(base.glob('*/fmriprep.toml'))
        if (base / 'fmriprep.toml').is_file():
            candidates.add(base / 'fmriprep.toml')
    fields = dict(environment=('version', 'nipype_version', 'templateflow_version'),
                  execution=('task_id', 'bids_filters', 'output_spaces', 'bids_dir'),
                  workflow=('me_t2s_fit_method', 'run_reconall', 'bold2anat_init', 'bold2anat_dof',
                            'use_bbr', 'slice_time_ref', 'dummy_scans', 'fmap_bspline', 'fmap_demean',
                            'use_syn_sdc', 'force_syn', 'skull_strip_template', 'skull_strip_fixed_seed'),
                  seeds=('master', 'ants', 'numpy'), nipype=('nprocs', 'omp_nthreads'))
    for path in sorted(candidates):
        try:
            data = tomllib.loads(path.read_text())
            result['execution_configs'].append(dict(path=str(path),
                settings={section: {key: data.get(section, {})[key] for key in keys
                                    if key in data.get(section, {})} for section, keys in fields.items()}))
        except (ValueError, OSError) as exc:
            result['execution_configs'].append(dict(path=str(path), error=str(exc)))
    report = directory / (subject + '.html')
    if report.is_file():
        text = html.unescape(re.sub('<[^>]+>', ' ', report.read_text(errors='replace')))
        result['subject_html_versions'] = sorted(set(re.findall(r'fmriprep(?:\s+version)?[\s-]+v?(\d+\.\d+\.\d+)', text, re.I)))
        for software in ('ANTs', 'FSL', 'Nipype', 'FreeSurfer', 'tedana', 'TemplateFlow'):
            values = sorted(set(re.findall(r'\b' + software + r'\b[^\d]{0,40}(\d+(?:\.\d+)+)', text, re.I)))
            if values:
                result['html_software_versions'][software] = values
    result['scope_note'] = 'Dataset/subject metadata may cover multiple runs or be overwritten; inspect execution filters and dates before assigning it to historical Trust.'
    return result


def checkpoint(old, new, name, mask):
    left, right = old / name, new / name
    result = dict(available=left.is_file() and right.is_file())
    if not result['available']:
        return result
    a, b = nib.load(left), nib.load(right)
    result['same_grid'] = bool(core.same_grid(a, b))
    if not result['same_grid']:
        return result
    x, y = a.get_fdata().squeeze(), b.get_fdata().squeeze()
    if x.shape != mask.shape or y.shape != mask.shape:
        result['error'] = 'Not compatible 3D checkpoints'
        return result
    result.update(agreement(x[mask], y[mask]))
    finite = mask & np.isfinite(x) & np.isfinite(y)
    a, b = x[finite], y[finite]
    denominator = float(a @ a)
    if denominator > 0:
        scale = float((a @ b) / denominator)
        result['new_on_old_global_scale'] = scale
        result['relative_rms_after_global_scale'] = float(np.linalg.norm(b - scale*a) / np.linalg.norm(a))
    return result


def commands(feat, replacements):
    path = feat / 'report_log.html'
    if not path.is_file():
        return []
    text = html.unescape(re.sub('<[^>]*>', '\n', path.read_text(errors='replace')))
    found = []
    for line in text.splitlines():
        match = re.search(r'\b(?:fslroi|fslmaths|fslstats|bet2|susan|film_gls|feat_model|mcflirt|slicetimer)\s+', line)
        if not match:
            continue
        line = line[match.start():].strip()
        for value, label in sorted(replacements.items(), key=lambda item: -len(item[0])):
            if value:
                line = line.replace(value, label)
        found.append(line)
    return found


def inspect(root, upstream):
    plan = json.loads((root / 'plan-pilot.json').read_text())
    selected = {key for cmd in plan['commands'] for key in cmd['runs']}
    runs = [row for row in plan['all_runs'] if row['key'] in selected]
    result = dict(container=json.loads((root / 'container.json').read_text()), subjects={}, runs={})
    print('CURRENT CONTAINER:', result['container']['version'], result['container']['path'], flush=True)
    for subject in sorted({row['subject'] for row in runs}):
        pair = {label: provenance(directory, subject) for label, directory in
                (('historical', upstream / 'derivatives/fmriprep-24'), ('rebuilt', root / 'fmriprep-24'))}
        result['subjects'][subject] = pair
        print('\n=== VERSION/CONFIG EVIDENCE:', subject, '===')
        for label, info in pair.items():
            print(label, 'dataset GeneratedBy:', json.dumps(info['dataset_generated_by']))
            print(label, 'subject HTML versions:', info['subject_html_versions'])
            print(label, 'HTML tool versions:', json.dumps(info['html_software_versions']))
            print(label, 'execution TOMLs:', len(info['execution_configs']))
            for config in info['execution_configs'][-3:]:
                print(json.dumps(config))
        print('Metadata scope: dataset/subject evidence; not automatically proof for the original Trust run.')
    for row in runs:
        key, subject = row['key'], row['subject']
        old = Path(row['original_feat'])
        new = root / 'feat_reconstruction' / (key + '.feat')
        original_fsf = root / 'historical_fsf' / (key + '.fsf')
        original = core.fsl_settings(original_fsf)
        saved = Path(row['confounds'])
        prefix = f'{subject}_ses-01_task-trust_run-{row["run"]}'
        info = dict(original_feat=str(old), rebuilt_feat=str(new),
                    saved_matrix_snapshot_hash_matches=core.sha(saved) == row['confounds_sha256'],
                    saved_vs_original_feat=matrix_match(saved, old / 'confoundevs.txt'),
                    saved_vs_rebuilt_feat=matrix_match(saved, new / 'confoundevs.txt'),
                    tedana=tedana_check(upstream / 'derivatives/tedana-24' / subject / 'ses-01', prefix, saved))
        print('\n===', key, '===')
        print('Saved matrix matches frozen hash:', info['saved_matrix_snapshot_hash_matches'])
        for name in ('saved_vs_original_feat', 'saved_vs_rebuilt_feat'):
            print(name, json.dumps(info[name]))
        print('TEDANA component verification:', json.dumps(info['tedana']))
        if not (new / 'design.fsf').is_file():
            info['rebuilt_feat_available'] = False
            print('No rebuilt FEAT design yet; checkpoint comparison skipped.')
            result['runs'][key] = info
            continue
        revised = core.fsl_settings(new / 'design.fsf')
        ignored = {'feat_files(1)', 'confoundev_files(1)', 'fmri(outputdir)'}
        differences = {key: [original.get(key), revised.get(key)] for key in sorted(set(original) | set(revised))
                       if key not in ignored and original.get(key) != revised.get(key)}
        info['fsf_setting_differences_except_locations'] = differences
        print('FSF setting differences excluding the three relocated paths:', json.dumps(differences))
        info['full_model_checks'] = {}
        for name in ('design.mat', 'design.con'):
            if (old / name).is_file() and (new / name).is_file():
                a, b = core.fsl_matrix(old / name), core.fsl_matrix(new / name)
                values = dict(old_shape=list(a.shape), new_shape=list(b.shape),
                              match=bool(a.shape == b.shape and np.allclose(a, b, atol=1e-8, rtol=1e-8)))
                if a.shape == b.shape:
                    values['max_absolute_difference'] = float(np.max(np.abs(a-b)))
                info['full_model_checks'][name] = values
                print(name, json.dumps(values))
        print('Reference selection:', {key: original.get(key) for key in
              ('fmri(npts)', 'fmri(ndelete)', 'fmri(alternateReference_yn)', 'fmri(mc)', 'fmri(st)', 'fmri(regunwarp_yn)')})
        if (old / 'mask.nii.gz').is_file() and (new / 'mask.nii.gz').is_file():
            a, b = nib.load(old / 'mask.nii.gz'), nib.load(new / 'mask.nii.gz')
            if core.same_grid(a, b):
                common = (a.get_fdata() > 0) & (b.get_fdata() > 0)
                info['checkpoints'] = {name: checkpoint(old, new, name, common)
                                       for name in ('example_func.nii.gz', 'mean_func.nii.gz')}
                for name, values in info['checkpoints'].items():
                    print(name, json.dumps(values))
        replacements = {row['old_bold']: '<BOLD>', row['new_bold']: '<BOLD>',
                        original.get('confoundev_files(1)', ''): '<CONFOUNDS>', str(saved): '<CONFOUNDS>',
                        str(old): '<FEAT>', str(new): '<FEAT>'}
        before, after = commands(old, replacements), commands(new, replacements)
        info['normalized_old_commands'] = before
        info['normalized_new_commands'] = after
        diff = list(difflib.unified_diff(before, after, fromfile='historical commands', tofile='rebuilt commands', lineterm=''))
        info['normalized_command_diff'] = diff
        print('Logged commands available:', len(before), 'historical;', len(after), 'rebuilt')
        print('Logged preprocessing/fit command differences (first 24 lines; full comparison saved):')
        print('\n'.join(diff[:24]) if diff else 'No normalized command differences found (check source-log availability).')
        result['runs'][key] = info
    directory = root / 'diagnostics'
    directory.mkdir(exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid.uuid4().hex[:6]
    destination = directory / ('trust-preprocessing-evidence-' + stamp + '.json')
    save_json(destination, result)
    print('\nFull diagnostic:', destination, flush=True)
    print('No preprocessing or model fits were run; this check does not launch production.', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rebuild-root', type=Path, required=True)
    parser.add_argument('--upstream-root', type=Path, default=Path('/ZPOOL/data/projects/rf1-sra-linux2'))
    args = parser.parse_args()
    root = args.rebuild_root.resolve()
    core.require(root.is_relative_to((core.ROOT / 'derivatives/insula_ica').resolve()), 'Expected the existing isolated rebuild root')
    inspect(root, args.upstream_root.resolve())


if __name__ == '__main__':
    main()
