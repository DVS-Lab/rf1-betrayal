#!/usr/bin/env python3
"""Collect sub-10834 registration reports from its recorded preprocessing root.

Read-only diagnostic: no rerun, transform edit, model fit, or exclusion. Only
small report assets, headers, configuration summaries and log excerpts are
written to a fresh Git record. No NIfTI data are copied.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import shutil
import uuid

import nibabel as nib

import insula_ica as core
from inspect_trust_preprocessing import provenance
from run_bold_controls import publish, save


SUBJECT = 'sub-10834'
SOURCE = 'linux2-trust-regeneration-20261027T162428Z-ef53c5'
REPORT_WORDS = ('coreg', 'bbregister', 'flirtbbr', 'norm', 'brainmask',
                'skullstrip', 'sdc', 'summary', 'validation', 'desc-rois')
LOG_PATTERN = re.compile(r'warn|error|fallback|reject|boundary|\bbbr\b|coreg|bold.?reg|'
                         r'anat.?norm|normalization|skull.?strip|single.?band|sbref', re.I)
LOG_NOISE = re.compile(r'getdefaultlocale|rcond|standardize|PySocks|dtype from|palette', re.I)


def selected_report(path):
    return (any(word in path.name.lower() for word in REPORT_WORDS)
            or path.name.endswith('_dseg.svg')
            or bool(re.search(r'_space-[^/]+_T1w\.svg$', path.name)))


def sbref_selection(rows, extra_logs=()):
    """Read actual selected references from retained report-generation log lines.

    A different run is a provenance flag, not by itself a registration failure.
    Missing log evidence remains unknown, never a pass.
    """
    evidence, logs = {}, {}
    pair = re.compile(r'Report - setting before \(([^)]+)\) and after \(([^)]+)\)')
    entities = re.compile(r'(sub-[^/_]+_ses-[^/_]+_task-[^/_]+_run-(\d+))')
    for row in rows:
        output = Path(row['new_bold']).parents[4]
        logs[output / (row['subject'] + '.fmriprep.log')] = None
    logs.update({Path(p): None for p in extra_logs})
    for log in logs:
        if not log.is_file():
            continue
        with log.open(errors='replace') as stream:
            for number, line in enumerate(stream, 1):
                fallback = re.search(r'No single-band-reference found for (\S+)', line)
                if fallback:
                    bold = entities.search(fallback[1])
                    if bold:
                        evidence[bold[1]] = dict(bold_run=bold[2], source_log=str(log), line_number=number,
                            status='bold_derived_reference_no_sbref')
                match = pair.search(line)
                if not match:
                    continue
                before, after = match.groups()
                if not before.endswith('_desc-coreg_boldref.nii.gz') or not after.endswith('_sbref.nii.gz'):
                    continue
                bold, reference = entities.search(Path(before).name), entities.search(Path(after).name)
                if bold and reference:
                    evidence[bold[1]] = dict(bold_run=bold[2], sbref_run=reference[2], selected_sbref=after,
                        source_log=str(log), line_number=number,
                        status='same_run' if bold[2] == reference[2] else 'different_run_reference')
    result=[]
    for row in rows:
        selected=evidence.get(row['key'])
        if selected is None:
            # Stock-filter reruns expose exactly one BOLD run. Their logs can
            # omit the SBRef comparison figure but explicitly list selected
            # references during workflow construction. Attribute that list only
            # after verifying the one-run view and all selected reference files.
            case=Path(row['new_bold']).parents[4]
            view=case/'bids-view'/row['subject']/'ses-01/func'
            expected={Path(p).name for p in row.get('raw_images',{})}
            actual={p.name for p in view.glob('*_bold.nii*')}
            log=case/(row['subject']+'.fmriprep.log')
            if expected and expected==actual and all(
                    (m:=entities.search(name)) and m[1]==row['key'] for name in actual) and log.is_file():
                choices=[]
                for number,line in enumerate(log.read_text(errors='replace').splitlines(),1):
                    if 'Using single-band reference file(s)' not in line:
                        continue
                    names=re.findall(r'sub-[^,\s/]+_sbref\.nii(?:\.gz)?',line)
                    if not names:
                        continue
                    choices.append((number,names))
                if choices and all(all((m:=entities.search(name)) and m[1]==row['key'] and
                                      (view/name).is_file() for name in names) for _,names in choices):
                    selected=dict(status='same_run',bold_run=entities.search(row['key'])[2],
                        sbref_run=entities.search(row['key'])[2],source_log=str(log),
                        line_number=choices[0][0],selected_sbrefs=choices[0][1],
                        evidence_type='explicit_selection_log_with_verified_single_run_input_view',
                        scoped_bold_files=sorted(actual))
        result.append(dict(key=row['key'],**(selected or dict(status='unknown_no_retained_log_evidence'))))
    return result


def log_excerpt(path, max_hits=160):
    lines = path.read_text(errors='replace').splitlines()
    hits = [i for i, line in enumerate(lines) if LOG_PATTERN.search(line) and not LOG_NOISE.search(line)]
    # Save context around the last relevant messages, plus the final log lines.
    selected = set(range(max(0, len(lines)-15), len(lines)))
    for i in hits[-max_hits:]:
        selected.update(range(max(0, i-1), min(len(lines), i+3)))
    text = '\n'.join(f'{i+1}: {lines[i]}' for i in sorted(selected))
    return text + '\n', dict(total_lines=len(lines), matching_lines=len(hits),
                             matching_lines_retained=min(len(hits), max_hits))


def collect(source, record):
    rows = [r for r in json.loads((source / 'input-registry.json').read_text())
            if r['subject'] == SUBJECT]
    plans = [p for p in json.loads((source / 'plan.json').read_text()) if p['subject'] == SUBJECT]
    core.require(len(rows) == 2 and len(plans) == 1, 'Expected two accepted sub-10834 runs and one recorded command')
    func = Path(rows[0]['new_bold']).parent
    subject_root = func.parent.parent
    core.require(subject_root.name == SUBJECT and subject_root.is_dir(), 'Recorded participant output is absent')
    # Restrict all report/header inspection to this exact recorded participant.
    derivatives = subject_root.parent
    assets = []
    figures = sorted(p for p in subject_root.rglob('*') if p.is_file() and p.suffix.lower() in ('.svg', '.png'))
    for path in figures:
        selected = selected_report(path)
        info = dict(path=str(path), size_bytes=path.stat().st_size, selected=selected)
        if selected:
            core.require(info['size_bytes'] <= 12*1024*1024, 'Oversize registration report; inspect separately: ' + str(path))
            relative = Path('figures') / path.relative_to(subject_root)
            target = record / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, target)
            info.update(saved=str(relative), sha256=core.sha(path))
        assets.append(info)
    headers, errors = [], []
    for path in sorted(subject_root.rglob('*.nii.gz')):
        # T1w, BOLD reference, masks, segmentation, and combined BOLD headers.
        if not any(s in path.name for s in ('T1w.nii.gz', '_boldref.nii.gz', '_mask.nii.gz', '_dseg.nii.gz', 'desc-preproc_bold.nii.gz')):
            continue
        try:
            image = nib.load(path)
            headers.append(dict(path=str(path), shape=list(image.shape), zooms=list(map(float, image.header.get_zooms())),
                affine=image.affine.tolist(), axis_codes=list(nib.aff2axcodes(image.affine)),
                qform_code=int(image.header['qform_code']), sform_code=int(image.header['sform_code'])))
        except Exception as exc:
            errors.append(dict(path=str(path), error=str(exc)))
    transforms = [dict(path=str(p), size_bytes=p.stat().st_size) for p in sorted(subject_root.rglob('*_xfm.*')) if p.is_file()]
    logs = []
    log = derivatives.parent / (SUBJECT + '.fmriprep.log')
    if log.is_file():
        text, counts = log_excerpt(log)
        (record / 'fmriprep-log-excerpt.txt').write_text(text)
        logs.append(dict(path=str(log), saved='fmriprep-log-excerpt.txt', **counts))
    details = provenance(derivatives, SUBJECT)
    save(record / 'diagnostic.json', dict(subject=SUBJECT, source_record=str(source),
        recorded_plan=plans[0], original_runs=rows, provenance=details,
        report_assets=assets, image_headers=headers, transforms=transforms, logs=logs, errors=errors,
        actions=dict(preprocessing_run=False, transforms_changed=False, exclusions=[])))
    copied = [a for a in assets if a.get('saved')]
    lines = ['# sub-10834 registration diagnostic', '',
        'Read-only collection from the exact completed preprocessing location. No source images or transforms were changed. '
        'No fMRIPrep, FEAT, TEDANA or Warpkit rerun was launched.', '',
        'Both regenerated Trust runs show gross displacement relative to the atlas insula. '
        'These reports are needed to distinguish anatomical normalization, BOLD-to-anatomical registration, '
        'brain masking or distortion-correction problems before selecting a correction.', '',
        'The complete fMRIPrep report remains on linux2: `' + str(derivatives / (SUBJECT + '.html')) + '`.', '',
        '## Registration and masking report assets', '']
    for asset in copied:
        lines += ['### ' + Path(asset['path']).name, '', '![' + Path(asset['path']).name + '](' + asset['saved'] + ')', '']
    if not copied:
        lines += ['No matching report assets were found. All available figure paths are listed in diagnostic.json.', '']
    lines += ['## Provenance and logs', '',
        'See diagnostic.json for image geometry, original command/work directory, registration settings and source figure paths. '
        'If present, fmriprep-log-excerpt.txt contains matching warning/registration lines and the end of the log.', '',
        f'Copied report assets: {len(copied)}. Header read errors: {len(errors)}. No NIfTI files were copied.', '']
    (record / 'README.md').write_text('\n'.join(lines))
    return len(copied), len(errors)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--publish', action='store_true')
    args = parser.parse_args()
    source = core.ROOT / 'code/insula_ica/records' / SOURCE
    core.require(source.is_dir(), 'Missing completed regeneration record')
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid.uuid4().hex[:6]
    record = core.ROOT / 'code/insula_ica/records' / ('linux2-registration-' + SUBJECT + '-' + stamp)
    record.mkdir(parents=True)
    assets, errors = collect(source, record)
    print(f'{SUBJECT}: {assets} registration/masking report assets collected; header errors: {errors}', flush=True)
    print('Record:', record, flush=True)
    if args.publish:
        publish(record, 'Collect sub-10834 registration reports before focused correction')
    raise SystemExit(bool(errors or not assets))


if __name__ == '__main__':
    main()
