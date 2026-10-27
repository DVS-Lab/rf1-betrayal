#!/usr/bin/env python3
"""Collect focused BOLD provenance evidence without launching preprocessing.

Sources are limited to pilot subjects at the historical/rebuilt RF1 paths and
the exact inputs in the frozen plan. The output contains no images or nuisance
time series. Unknown historical inputs remain unknown, rather than being
inferred from the current BIDS files or a dataset-level version label.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import html
import json
import os
from pathlib import Path
import re
import subprocess
import tarfile
import tomllib
import uuid

import nibabel as nib
import numpy as np

from rebuild_trust import ROOT, require, sha
from diagnose_film_replay import agreement


MAX_TEXT = 4 * 1024 * 1024
METADATA_KEYS = (
    'GeneratedBy', 'PipelineDescription', 'Sources', 'RawSources', 'SpatialReference',
    'SkullStripped', 'TaskName', 'RepetitionTime', 'EchoTime', 'EchoTimes',
    'PhaseEncodingDirection', 'TotalReadoutTime', 'EffectiveEchoSpacing',
    'SliceTiming', 'MultibandAccelerationFactor', 'NumberOfVolumesDiscardedByScanner',
    'NumberOfVolumesDiscardedByUser', 'IntendedFor', 'B0FieldIdentifier',
    'B0FieldSource', 'Units', 'SoftwareVersions', 'PulseSequenceDetails',
)


def save(path, value):
    path.write_text(json.dumps(value, indent=2, default=str, allow_nan=False) + '\n')


def identity(path):
    path = Path(path)
    try:
        if not path.is_file():
            return dict(path=str(path), exists=False)
        stat = path.stat()
        return dict(path=str(path), exists=True, size=stat.st_size, mtime_ns=stat.st_mtime_ns)
    except OSError as exc:
        return dict(path=str(path), exists=False, error=str(exc))


def frozen_check(path, expected):
    """Size/mtime is a drift check, never a historical content-hash proof."""
    result = identity(path)
    if not result['exists']:
        result['matches_rebuild_snapshot'] = False
    elif 'sha256' in expected:
        try:
            result['sha256'] = sha(path)
            result['matches_rebuild_snapshot'] = result['sha256'] == expected['sha256']
        except OSError as exc:
            result.update(error=str(exc), matches_rebuild_snapshot=False)
        result['comparison'] = 'sha256'
    else:
        result['comparison'] = 'size_and_mtime_only'
        result['matches_rebuild_snapshot'] = all(result[k] == expected[k] for k in ('size', 'mtime_ns'))
    return result


def task_scope(config):
    execution = config.get('execution', {})
    task = execution.get('task_id')
    filters = execution.get('bids_filters', {})
    if not task and isinstance(filters, dict):
        task = filters.get('bold', {}).get('task') if isinstance(filters.get('bold'), dict) else None
    if not task:
        return 'task_not_recorded'
    tasks = task if isinstance(task, list) else [task]
    return 'trust_named' if 'trust' in tasks else 'other_task'


def settings(config):
    # Preserve all workflow keys, including ones not anticipated by an older
    # version of this diagnostic. Paths/resources are recorded separately.
    result = {f'workflow.{k}': v for k, v in config.get('workflow', {}).items()}
    for section in ('environment', 'seeds', 'nipype'):
        result.update({f'{section}.{k}': v for k, v in config.get(section, {}).items()})
    for key in ('output_spaces', 'output_layout', 'anat_derivatives', 'derivatives',
                'bids_filters', 'task_id', 'echo_idx', 'participant_label'):
        if key in config.get('execution', {}):
            result[f'execution.{key}'] = config['execution'][key]
    return result


def setting_differences(old, new):
    a, b = settings(old), settings(new)
    return {key: dict(historical_present=key in a, rebuilt_present=key in b,
                      historical=a.get(key), rebuilt=b.get(key))
            for key in sorted(set(a) | set(b)) if key not in a or key not in b or a[key] != b[key]}


def excerpt(text):
    text = html.unescape(re.sub('<[^>]+>', ' ', text))
    text = re.sub(r'[ \t]+', ' ', text)
    # Keep report methods as connected text; do not infer an acquisition date
    # or a causal version difference from a keyword match.
    return text.strip()[:80000]


class Collector:
    def __init__(self, output):
        self.output = output
        self.sources = []
        self.errors = []
        (output / 'sources').mkdir()

    def read(self, path, kind):
        info = identity(path)
        if not info['exists']:
            return None
        try:
            if info['size'] > MAX_TEXT:
                raise ValueError('Text source exceeds 4 MiB; not copied')
            data = path.read_bytes()
            info['sha256'] = hashlib.sha256(data).hexdigest()
            text = data.decode('utf-8', errors='replace')
            if kind == 'toml':
                value = tomllib.loads(text)
                payload, suffix = data, '.toml'
            elif kind == 'json':
                full = json.loads(text)
                value = {k: full[k] for k in METADATA_KEYS if k in full}
                payload = (json.dumps(value, indent=2) + '\n').encode()
                suffix = '.selected.json'
            else:
                value = excerpt(text)
                payload, suffix = (value + '\n').encode(), '.txt'
            name = hashlib.sha256(str(path).encode()).hexdigest()[:12] + '-' + path.name + suffix
            destination = self.output / 'sources' / name
            destination.write_bytes(payload)
            info.update(kind=kind, saved_as=str(destination.relative_to(self.output)),
                        saved_sha256=hashlib.sha256(payload).hexdigest())
            self.sources.append(info)
            return value
        except (OSError, ValueError, TypeError) as exc:
            self.errors.append(dict(path=str(path), error=str(exc)))
            return None

    def provenance(self, directory, subject, keys):
        result = dict(root=str(directory), configs=[], methods=[], metadata=[])
        for path in (directory / 'dataset_description.json', directory / (subject + '.html')):
            value = self.read(path, 'json' if path.suffix == '.json' else 'text')
            if value is not None:
                result['metadata'].append(dict(path=str(path), value=value))
                if path.suffix == '.html':
                    result['subject_report_fmriprep_versions'] = sorted(set(re.findall(
                        r'fmriprep(?:\s+version)?[\s-]+v?(\d+\.\d+\.\d+)', value, re.I)))
        base = directory / subject
        if not base.is_dir():
            result['subject_directory_available'] = False
            return result
        result['subject_directory_available'] = True
        # Traverse only this pilot subject, never other projects or subjects.
        # Read only explicit config/report/sidecar formats, never Nipype pickles.
        for parent, dirs, files in os.walk(base, followlinks=False, onerror=lambda e: self.errors.append(dict(error=str(e)))):
            dirs[:] = [d for d in dirs if not d.startswith('.')]
            for name in sorted(files):
                path = Path(parent) / name
                if name == 'fmriprep.toml':
                    config = self.read(path, 'toml')
                    if config is not None:
                        result['configs'].append(dict(path=str(path), task_scope=task_scope(config),
                                                      settings=settings(config), full=config))
                elif name.endswith(('.html', '.md', '.bib')) and any(x in name.lower() for x in ('summary', 'about', 'citation', 'methods')):
                    value = self.read(path, 'text')
                    if value is not None:
                        result['methods'].append(dict(path=str(path), text=value))
                elif name.endswith('.json') and (any(name.startswith(k + '_') for k in keys) or
                        ('T1w' in name and ('desc-preproc' in name or 'xfm' in name))):
                    value = self.read(path, 'json')
                    if value is not None:
                        result['metadata'].append(dict(path=str(path), value=value))
        for folder in (directory / 'logs', directory / 'log'):
            for name in ('CITATION.md', 'CITATION.html', 'CITATION.bib'):
                path = folder / name
                value = self.read(path, 'text')
                if value is not None:
                    result['methods'].append(dict(path=str(path), text=value, scope='dataset_level'))
        return result


def frame_comparison(reference, candidate, mask, index):
    a, b = nib.load(reference), nib.load(candidate)
    info = dict(path=str(candidate), shape=list(b.shape), frame_index=index)
    same = a.shape[:3] == b.shape[:3] and np.allclose(a.affine, b.affine, atol=1e-4)
    info['same_grid'] = bool(same)
    if not same or len(b.shape) != 4 or index >= b.shape[3]:
        info['compared'] = False
        return info
    x = a.get_fdata().squeeze()
    # Array proxy reads only the requested frame, not the entire 4D series.
    y = np.asarray(b.dataobj[..., index], dtype=float)
    finite = (mask > 0) & np.isfinite(x) & np.isfinite(y)
    info.update(agreement(x[finite], y[finite]))
    xa, ya = x[finite], y[finite]
    denom = float(xa @ xa)
    if denom:
        scale = float((xa @ ya) / denom)
        info.update(new_on_old_global_scale=scale,
                    relative_rms_after_global_scale=float(np.linalg.norm(ya-scale*xa) / np.sqrt(denom)))
    info['compared'] = True
    return info


def build_report(data):
    lines = ['# Trust BOLD preprocessing provenance', '',
             'No preprocessing, model fits, or production jobs were launched.', '',
             'The original nuisance fit is already preserved; this check does not require standalone TEDANA files.', '',
             'Historical dataset/subject metadata can have been overwritten. Config comparisons below are candidate evidence, not proof that a config generated the old Trust images.', '',
             '## Version and execution evidence', '']
    lines += ['Rebuild container record: `' + json.dumps(data['container']) + '`', '']
    for subject, pair in data['subjects'].items():
        lines += [f'### {subject}', '']
        for label in ('historical', 'rebuilt'):
            p = pair[label]
            lines += [f'{label}: {len(p["configs"])} execution configs; {len(p["methods"])} methods/summary reports.']
            lines += [f'Subject HTML version strings: `{p.get("subject_report_fmriprep_versions", [])}` (subject-level evidence only).']
            for item in p['metadata']:
                if isinstance(item['value'], dict) and ('GeneratedBy' in item['value'] or 'PipelineDescription' in item['value']):
                    lines += [f'- Dataset-level version evidence: `{json.dumps(item["value"])}`']
            for cfg in p['configs']:
                lines += [f'- `{cfg["path"]}` — {cfg["task_scope"]}; environment: `{json.dumps(cfg["full"].get("environment", {}), default=str)}`']
            for method in p['methods']:
                lines += ['', f'Report `{method["path"]}`:', '', method['text'][:5000], '']
        if not pair['config_comparisons']:
            lines += ['', 'No historical/rebuilt candidate config pair is available. Exact historical execution differences remain unknown.']
        for comp in pair['config_comparisons']:
            lines += ['', f'Config comparison: `{comp["historical"]}` vs `{comp["rebuilt"]}`', '',
                      '```json', json.dumps(comp['differences'], indent=2, default=str), '```']
    lines += ['', '## Frozen rebuild inputs', '',
              'These comparisons detect changes since the rebuild snapshot, not changes since the original historical analysis. Image size/mtime matches are not content-hash proof.', '']
    for key, row in data['runs'].items():
        checks = row['frozen_input_checks']
        mismatches = [x for x in checks if x.get('matches_rebuild_snapshot') is False]
        lines += [f'### {key}', '', f'{len(checks)} snapshot checks; {len(mismatches)} missing/changed.', '']
        for check in mismatches:
            lines += ['- `' + json.dumps(check) + '`']
        lines += ['Historical reference vs rebuilt outputs (same midpoint index; historical FEAT mask):', '']
        for item in row['rebuilt_frame_candidates']:
            if item.get('compared'):
                lines += [f'- `{item["path"]}`: r={item.get("spatial_r")}; relative RMS={item.get("relative_rms_difference")}; RMS after global scale={item.get("relative_rms_after_global_scale")}']
            else:
                lines += [f'- `{item["path"]}`: not compared ({item.get("error", "different grid or incompatible dimensions")}).']
        if not row['rebuilt_frame_candidates']:
            lines += ['- No comparable rebuilt-output candidates found.']
    lines += ['', '## Errors and limits', '', f'{len(data["errors"])} collection/read errors (see evidence.json).', '',
              'No config difference is automatically attributed as the cause. Check task/session coverage and timestamps before choosing a controlled rerun. Matching a midpoint frame alone cannot prove full-series equivalence.', '',
              'Full settings, report text, source identities and snapshot details are in evidence.json and sources/. No images, full nuisance matrices, or participant tables are packaged.', '']
    lines += ['- `' + json.dumps(e) + '`' for e in data['errors']]
    return '\n'.join(lines)


def collect(root, upstream, output):
    plan = json.loads((root / 'plan-pilot.json').read_text())
    selected = {key for cmd in plan['commands'] for key in cmd['runs']}
    runs = [r for r in plan['all_runs'] if r['key'] in selected]
    collector = Collector(output)
    try:
        commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = None
    data = dict(schema_version=1, created_utc=datetime.now(timezone.utc).isoformat(), code_commit=commit,
                root=str(root), upstream=str(upstream), container=json.loads((root / 'container.json').read_text()),
                subjects={}, runs={}, errors=collector.errors, sources=collector.sources)
    for subject in sorted({r['subject'] for r in runs}):
        print('Collecting configs/reports:', subject, flush=True)
        keys = [r['key'] for r in runs if r['subject'] == subject]
        pair = {label: collector.provenance(folder, subject, keys) for label, folder in
                (('historical', upstream / 'derivatives/fmriprep-24'), ('rebuilt', root / 'fmriprep-24'))}
        pair['config_comparisons'] = [dict(historical=a['path'], rebuilt=b['path'],
                                         differences=setting_differences(a['full'], b['full']))
            for a in pair['historical']['configs'] for b in pair['rebuilt']['configs']
            if a['task_scope'] != 'other_task' and b['task_scope'] != 'other_task']
        data['subjects'][subject] = pair
    for row in runs:
        key = row['key']
        print('Checking frozen inputs and available echo/combined frames:', key, flush=True)
        checks = [frozen_check(Path(p), expected) for p, expected in row.get('raw_images', {}).items()]
        checks += [frozen_check(Path(p), dict(sha256=h)) for p, h in row.get('raw_metadata_sha256', {}).items()]
        frozen = root / 'frozen_warpkit' / (row['subject'] + '.json')
        if frozen.is_file():
            warpkit = json.loads(frozen.read_text()).get(key, {})
            checks += [frozen_check(Path(p), expected) for p, expected in warpkit.get('files', {}).items()]
        raw_metadata = []
        for p in row.get('raw_metadata_sha256', {}):
            value = collector.read(Path(p), 'json')
            if value is not None:
                raw_metadata.append(dict(path=p, metadata=value))
        old = Path(row['original_feat'])
        results = []
        func = Path(row['new_bold']).parent
        try:
            mask_img = nib.load(old / 'mask.nii.gz')
            reference = old / 'example_func.nii.gz'
            ref_img = nib.load(reference)
            require(ref_img.shape[:3] == mask_img.shape[:3] and np.allclose(ref_img.affine, mask_img.affine),
                    'Historical reference/mask grid mismatch')
            mask = mask_img.get_fdata()
            # Include echo-specific/T1w outputs; grid checks prevent accidental
            # voxelwise comparison across spaces. Never resample in this check.
            for candidate in sorted(func.glob(key + '_*desc-preproc_bold.nii.gz')):
                try:
                    results.append(frame_comparison(reference, candidate, mask, int(row['nvol']) // 2))
                except (OSError, ValueError, nib.filebasedimages.ImageFileError) as exc:
                    results.append(dict(path=str(candidate), compared=False, error=str(exc)))
        except (OSError, ValueError, nib.filebasedimages.ImageFileError) as exc:
            collector.errors.append(dict(run=key, error=str(exc)))
        data['runs'][key] = dict(frozen_input_checks=checks, raw_metadata=raw_metadata,
                                rebuilt_frame_candidates=results)
    # Include the preceding full diagnostic, if present, so a truncated terminal
    # paste no longer loses the version section. No recomputation or TEDANA read.
    previous = sorted((root / 'diagnostics').glob('trust-preprocessing-evidence-*.json'))
    if previous:
        path = previous[-1]
        try:
            value = json.loads(path.read_text())
            save(output / 'previous-preprocessing-evidence.json', value)
            data['previous_diagnostic'] = dict(path=str(path), sha256=sha(path),
                                               saved_as='previous-preprocessing-evidence.json')
        except (OSError, ValueError) as exc:
            collector.errors.append(dict(path=str(path), error=str(exc)))
    save(output / 'evidence.json', data)
    (output / 'report.md').write_text(build_report(data))
    files = sorted(p for p in output.rglob('*') if p.is_file())
    (output / 'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.relative_to(output)}\n' for p in files))
    archive = output.with_suffix('.tar.gz')
    with tarfile.open(archive, 'w:gz') as bundle:
        bundle.add(output, arcname=output.name)
    print('\n' + (output / 'report.md').read_text(), flush=True)
    print('Saved record:', output, flush=True)
    print('Transfer archive:', archive, flush=True)
    print('No preprocessing or model fits were run; production remains unchanged.', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rebuild-root', type=Path, required=True)
    parser.add_argument('--upstream-root', type=Path, default=Path('/ZPOOL/data/projects/rf1-sra-linux2'))
    args = parser.parse_args()
    root = args.rebuild_root.resolve()
    require(root.is_relative_to((ROOT / 'derivatives/insula_ica').resolve()), 'Expected isolated rebuild root')
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid.uuid4().hex[:6]
    output = ROOT / 'code/insula_ica/records' / ('linux2-bold-provenance-' + stamp)
    output.mkdir(parents=True, exist_ok=False)
    print('Saving GitHub record to:', output, flush=True)
    collect(root, args.upstream_root.resolve(), output)


if __name__ == '__main__':
    main()
