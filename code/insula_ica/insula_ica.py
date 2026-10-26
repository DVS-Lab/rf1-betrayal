#!/usr/bin/env python3
"""Bounded insula ICA: audit -> prepare -> ICA -> regress -> summarize.

No existing FEAT or upstream derivatives are modified. See README before use.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict

import nibabel as nib
from nibabel.processing import resample_from_to
import numpy as np
from scipy import linalg, stats

ROOT = Path(__file__).resolve().parents[2]
ROIS = ('left', 'right', 'bilateral')
TASKS = {'trust': ('01', 10, 'rec-def'), 'ugr': ('3', 11, 'offer (un)fairness (pmod)')}
BASE_CONFOUNDS = [f'a_comp_cor_{i:02d}' for i in range(6)] + [
    'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z', 'framewise_displacement']


def require(ok, message):
    if not ok:
        raise ValueError(message)


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def read_table(path, delimiter='\t'):
    with Path(path).open(newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def write_table(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text('')
        return
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def identity(path, content=True):
    p = Path(path).resolve()
    s = p.stat()
    return {'path': str(p), 'size': s.st_size, 'mtime_ns': s.st_mtime_ns,
            'sha256': sha(p) if content else None}


def run(cmd, directory):
    """Argument arrays only; save complete command, stdout, stderr, and exit."""
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / 'commands.jsonl').open('a') as f:
        f.write(json.dumps({'argv': list(map(str, cmd)), 'shell_display': shlex.join(list(map(str, cmd)))}) + '\n')
    with (directory / 'commands.log').open('a') as f:
        f.write(shlex.join(list(map(str, cmd))) + '\n')
        f.flush()
        result = subprocess.run(list(map(str, cmd)), stdout=f, stderr=subprocess.STDOUT)
        f.write(f'\nexit={result.returncode}\n')
    require(result.returncode == 0, f'Command failed; see {directory / "commands.log"}')


def save_image(path, data, reference):
    header = reference.header.copy()
    header.set_data_dtype(np.asarray(data).dtype)
    image = nib.Nifti1Image(data, reference.affine, header)
    image.set_qform(reference.affine, code=1)
    image.set_sform(reference.affine, code=1)
    nib.save(image, path)


def same_grid(a, b):
    return a.shape[:3] == b.shape[:3] and np.allclose(a.affine, b.affine, atol=1e-4, rtol=0)


def load_numeric(path):
    a = np.loadtxt(path, ndmin=2)
    require(np.isfinite(a).all(), f'Nonfinite matrix: {path}')
    return a


def basis(x):
    """Orthonormal basis of the column space, robust to zero/redundant confounds."""
    x = np.asarray(x, dtype=float)
    scale = np.linalg.norm(x, axis=0)
    x = x[:, scale > 1e-10] / scale[scale > 1e-10]
    if x.shape[1] == 0:
        return np.zeros((x.shape[0], 0))
    u, s, _ = linalg.svd(x, full_matrices=False)
    return u[:, s > s[0] * 1e-10]


def residualize(q, y):
    return y - q @ (q.T @ y)


def fsl_settings(path):
    return {k: v.strip().strip('"') for k, v in re.findall(r'^set\s+(\S+)\s+(.+?)\s*$', Path(path).read_text(), re.M)}


def fsl_matrix(path):
    text = Path(path).read_text()
    require('/Matrix' in text, f'Not a VEST matrix: {path}')
    import io
    return np.loadtxt(io.StringIO(text.split('/Matrix', 1)[1]), ndmin=2)


def task_design(fsf_path, dest, task, nvol, tr, registry):
    """Render task-only EVs using original FSF and FSL's own feat_model."""
    d = fsl_settings(fsf_path)
    expected, contrast, name = TASKS[task]
    require(float(d['fmri(tr)']) == float(tr) or np.isclose(float(d['fmri(tr)']), tr, atol=1e-5), 'FSF TR mismatch')
    require(int(d['fmri(npts)']) == nvol and int(d['fmri(ndelete)']) == 0, 'FSF length/deletion mismatch')
    require(d['fmri(temphp_yn)'] == '0', 'Existing FSF high-pass is enabled: review filtering before proceeding')
    require(d[f'fmri(conname_real.{contrast})'] == name, 'Unexpected task contrast label')
    nev = int(d['fmri(evs_real)'])
    require(nev == int(d['fmri(evs_orig)']) == (10 if task == 'trust' else 11), 'Unexpected EV/derivative layout')
    titles = [d[f'fmri(evtitle{i})'] for i in range(1, nev + 1)]
    expected_titles = (['c_C', 'c_F', 'c_S', 'C_def', 'C_rec', 'F_def', 'F_rec', 'S_def', 'S_rec', 'miss']
                       if task == 'trust' else ['nonsocial_high_constant', 'nonsocial_high_pmod',
                       'nonsocial_low_constant', 'nonsocial_low_pmod', 'social_high_constant',
                       'social_high_pmod', 'social_low_constant', 'social_low_pmod', 'rt_constant', 'rt_pmod', 'miss'])
    require(titles == expected_titles, 'Unexpected EV ordering')
    events = []
    for i in range(1, nev + 1):
        require(d.get(f'fmri(deriv_yn{i})', '0') == '0', 'Unexpected temporal derivative')
        if d[f'fmri(shape{i})'] == '10':
            continue
        require(d[f'fmri(shape{i})'] == '3', 'Expected existing three-column EV')
        ev = Path(d[f'fmri(custom{i})'])
        require(ev.is_absolute() and ev.is_file(), f'Missing existing EV: {ev}')
        a = load_numeric(ev)
        require(a.shape[1] == 3 and (a[:, 1] >= 0).all(), f'Invalid EV timing: {ev}')
        require((a[:, 0] >= 0).all() and (a[:, 0] < nvol * tr).all() and
                (a[:, 0] + a[:, 1] <= nvol * tr + tr).all(), f'EV outside run: {ev}')
        registry.append(identity(ev))
        events.append({'ev': i, 'title': titles[i-1], 'path': str(ev), 'rows': len(a),
                       'amplitude_min': float(a[:, 2].min()), 'amplitude_max': float(a[:, 2].max())})
    source = Path(fsf_path).read_text()
    source = re.sub(r'^set fmri\(confoundevs\).*$', 'set fmri(confoundevs) 0', source, flags=re.M)
    # Write only inside analysis output; feat_model never runs FEAT or changes inputs.
    dest.mkdir(parents=True, exist_ok=True)
    stem = dest / 'task'
    stem.with_suffix('.fsf').write_text(source)
    run(['feat_model', stem], dest)
    x = fsl_matrix(stem.with_suffix('.mat'))
    require(x.shape == (nvol, nev) and np.isfinite(x).all(), 'Rendered task design shape mismatch')
    c = np.array([float(d.get(f'fmri(con_real{contrast}.{i})', 0)) for i in range(1, nev + 1)])
    expected_c = (np.array([0, 0, 0, -1, 1, -1, 1, -1, 1, 0]) if task == 'trust'
                  else np.array([0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0]))
    require(np.array_equal(c, expected_c), 'Unexpected contrast weights')
    np.save(dest / 'task.npy', x)
    write_json(dest / 'events.json', events)
    write_json(dest / 'contrast.json', {'original_label': name, 'original_weights': c.tolist(), 'model': expected})
    return x, c


def nuisance(combined_path, fmriprep_path, nvol, tr, cutoff):
    """Verify legacy matrix prefix, NSS and cosine terms; add drift only if absent."""
    rows = read_table(fmriprep_path)
    require(len(rows) == nvol, 'fMRIPrep confound length mismatch')
    names = list(rows[0])
    require(all(c in names for c in BASE_CONFOUNDS), 'Required motion/aCompCor columns missing')
    selected = BASE_CONFOUNDS + [c for c in names if c.startswith('cosine')] + [c for c in names if c.startswith('non_steady_state')]
    def number(s):
        return 0.0 if s in ('n/a', 'NaN', 'nan', '') else float(s)
    base = np.array([[number(row[c]) for c in selected] for row in rows])
    combined = load_numeric(combined_path)
    require(combined.shape[0] == nvol and combined.shape[1] >= base.shape[1], 'Combined confound shape mismatch')
    require(np.allclose(combined[:, :base.shape[1]], base, atol=1e-5, rtol=1e-5),
            'Legacy confounds differ from the verified generator layout; do not infer column identities')
    nss = [c for c in selected if c.startswith('non_steady_state')]
    flags = np.array([[number(row[c]) for c in nss] for row in rows]) if nss else np.zeros((nvol, 0))
    require(np.isin(flags, [0, 1]).all() and (flags.sum(axis=0) == 1).all() and (flags.sum(axis=1) <= 1).all(), 'Invalid NSS columns')
    nss_count = int(flags.sum())
    require(np.array_equal(np.flatnonzero(flags.sum(axis=1)), np.arange(nss_count)), 'Noncontiguous NSS volumes')
    drift = np.zeros((nvol, 0))
    if not any(c.startswith('cosine') for c in selected):
        k = int(np.floor(2 * nvol * tr / cutoff))
        drift = np.cos(np.pi * (np.arange(nvol)[:, None] + .5) * np.arange(1, k + 1) / nvol)
    n = np.column_stack([np.ones(nvol), combined, drift])
    q = basis(n)
    return q, {'nss': nss_count, 'nuisance_rank': q.shape[1], 'columns': combined.shape[1],
               'added_dct_columns': drift.shape[1], 'known_columns': selected,
               'tedana_columns': combined.shape[1] - base.shape[1],
               'drift_policy': 'existing cosines' if drift.shape[1] == 0 else f'joint DCT cutoff {cutoff}s'}


def anatomical_masks(atlas_path, xml_path, reference):
    labels = [x for x in ET.parse(xml_path).findall('.//label') if (x.text or '').strip().lower() == 'insular cortex']
    require(len(labels) == 1, 'Atlas must contain one Insular Cortex label')
    label = int(labels[0].attrib['index']) + 1
    atlas = nib.load(atlas_path)
    inds = np.argwhere(np.asanyarray(atlas.dataobj) == label)
    xyz = nib.affines.apply_affine(atlas.affine, inds)
    masks, counts = {}, {}
    for roi in ROIS:
        chosen = inds[xyz[:, 0] < 0] if roi == 'left' else inds[xyz[:, 0] > 0] if roi == 'right' else inds
        binary = np.zeros(atlas.shape, dtype=np.uint8)
        binary[tuple(chosen.T)] = 1
        image = resample_from_to(nib.Nifti1Image(binary, atlas.affine), (reference.shape[:3], reference.affine), order=0)
        masks[roi] = np.asanyarray(image.dataobj) > .5
        counts[roi] = {'native_voxels': len(chosen), 'grid_voxels': int(masks[roi].sum())}
        require(masks[roi].sum() > 7, f'Empty/tiny {roi} ROI')
    return masks, counts, label


def config(path):
    c = json.loads(Path(path).read_text())
    for key in ('master_qc', 'output_root'):
        c[key] = str((ROOT / c[key]).resolve()) if not Path(c[key]).is_absolute() else c[key]
    require(c['expected_tr'] > 0 and c['drift_cutoff_seconds'] > 0, 'TR/drift cutoff must be positive')
    require(0 < c['minimum_roi_coverage'] <= 1 and c['min_residual_dof'] > 0, 'Invalid coverage/dof gate')
    require(c['max_condition_number'] > 1, 'Invalid condition-number gate')
    return c


def software():
    needed = ('melodic', 'fsl_glm', 'feat_model')
    require(all(shutil.which(x) for x in needed), 'Source FSL environment: melodic, fsl_glm, feat_model required')
    fsldir = Path(os.environ['FSLDIR'])
    return {'fsl': (fsldir / 'etc/fslversion').read_text().strip(),
            'python': sys.version, 'numpy': np.__version__, 'scipy': __import__('scipy').__version__,
            'nibabel': nib.__version__, 'executables': {x: shutil.which(x) for x in needed}}


def audit(c):
    out = Path(c['output_root']) / 'audit'
    require(not out.exists(), f'Audit exists: {out}. Choose a fresh output_root to preserve provenance.')
    out.mkdir(parents=True)
    blockers, registry, rows = [], [], []
    sw = software()
    for exe in ('melodic', 'fsl_glm'):
        p = subprocess.run([exe, '--help'], capture_output=True, text=True)
        (out / f'{exe}_help.txt').write_text(p.stdout + p.stderr)
    help_text = (out / 'melodic_help.txt').read_text()
    require(all(flag in help_text for flag in ('--tr', '--nobet', '--Oall', '--disableMigp', '--vn', '--update_mask')), 'Installed MELODIC lacks required options; review version')
    cohort_path = ROOT / 'code/sublist_n132.txt'
    cohort = [x.removeprefix('sub-') for x in cohort_path.read_text().split()]
    require(len(cohort) == len(set(cohort)) == c['expected_cohort_size'], 'Cohort identity/size mismatch')
    registry.append(identity(cohort_path))
    qc_path = Path(c['master_qc'])
    qc = {}
    if qc_path.exists():
        registry.append(identity(qc_path))
        records = read_table(qc_path, ',')
        require(len({r['sub_id'] for r in records}) == len(records), 'Duplicate master QC subjects')
        qc = {str(int(float(r['sub_id']))): r for r in records}
    else:
        blockers.append(f'Missing saved manuscript run-aligned QC table: {qc_path}')
    review = c['provenance_review']
    for key in ('signal_description', 'evidence', 'ugr_coding_evidence'):
        if not review.get(key):
            blockers.append(f'Unresolved provenance_review.{key}')
    if review.get('tedana_noise_not_already_removed') is not True:
        blockers.append('Must verify selected BOLD has not already had the appended TEDANA regressors removed')
    if not isinstance(review.get('ugr_pmods_increase_with_fairness'), bool):
        blockers.append('UGR pmod direction must be verified from existing EV/event coding')
    upstream = Path(c['upstream_root'])
    reference, common = None, None
    for sub in cohort:
        for task, (model, _, _) in TASKS.items():
            for run_number in (1, 2):
                key = f'sub-{sub}_ses-01_task-{task}_run-{run_number}'
                prefix = f'{key}_part-mag'
                func = upstream / 'derivatives' / c['fmriprep_derivative'] / f'sub-{sub}/ses-01/func'
                bold = func / f'{prefix}_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz'
                conf = upstream / 'derivatives' / c['combined_confounds_derivative'] / f'sub-{sub}/{key}_desc-TedanaPlusConfounds.tsv'
                fconf = func / f'{prefix}_desc-confounds_timeseries.tsv'
                mask = func / f'{prefix}_space-MNI152NLin6Asym_desc-brain_mask.nii.gz'
                sidecar = bold.with_name(bold.name.removesuffix('.nii.gz') + '.json')
                fsf = ROOT / f'derivatives/fsl/sub-{sub}/ses-01/L1_task-{task}_model-{model}_type-act_run-{run_number}_sm-5.feat/design.fsf'
                row = {'key': key, 'subject': sub, 'task': task, 'run': run_number,
                       'bold': str(bold), 'confounds': str(conf), 'fmriprep_confounds': str(fconf),
                       'brain_mask': str(mask), 'bold_json': str(sidecar), 'source_fsf': str(fsf)}
                row.update(original_behavioral_ok=qc.get(sub, {}).get(f'{task}_beh_run{run_number}_ok', ''),
                           original_ev_ok=qc.get(sub, {}).get(f'{task}_ev_run{run_number}_ok', ''),
                           original_mriqc_ok=qc.get(sub, {}).get(f'{task}_mri_run{run_number}_ok', ''))
                rows.append(row)
                value = qc.get(sub, {}).get(f'{task}_run{run_number}_full_ok', '')
                if value.upper() in ('FALSE', '0'):
                    row.update(status='excluded_original_qc', reason='saved full_ok = false')
                    continue
                if value.upper() not in ('TRUE', '1'):
                    row.update(status='blocked', reason='missing or unknown original run QC')
                    continue
                try:
                    require(qc[sub].get('include_final', '').upper() in ('TRUE', '1'), 'Cohort member not include_final in master QC')
                    paths = [bold, conf, fconf, mask, sidecar, fsf]
                    require(all(p.is_file() for p in paths), 'Missing: ' + '; '.join(str(p) for p in paths if not p.is_file()))
                    im, brain = nib.load(bold), nib.load(mask)
                    require(len(im.shape) == 4 and im.shape[3] > 10 and same_grid(im, brain), 'Not aligned whole-brain 4D BOLD/mask')
                    require(im.header.get_xyzt_units()[1] == 'sec', 'Unverified NIfTI temporal units')
                    tr = float(im.header.get_zooms()[3])
                    require(np.isclose(tr, c['expected_tr'], atol=1e-5), f'TR mismatch {tr}')
                    qform, qcode = im.get_qform(coded=True)
                    sform, scode = im.get_sform(coded=True)
                    require(not(qcode and scode) or np.allclose(qform, sform, atol=1e-3), 'Conflicting qform/sform')
                    meta = json.loads(sidecar.read_text())
                    if 'RepetitionTime' in meta:
                        require(np.isclose(float(meta['RepetitionTime']), tr, atol=1e-5), 'BOLD sidecar TR mismatch')
                    nvol = im.shape[3]
                    source_settings = fsl_settings(fsf)
                    require(Path(source_settings['feat_files(1)']).resolve() == bold.resolve(), 'Published FSF BOLD source differs from configured derivative')
                    require(Path(source_settings['confoundev_files(1)']).resolve() == conf.resolve(), 'Published FSF confounds differ')
                    registry.extend(identity(p, content=p != bold) for p in paths)
                    q, ninfo = nuisance(conf, fconf, nvol, tr, c['drift_cutoff_seconds'])
                    x, contrast = task_design(fsf, out / key, task, nvol, tr, registry)
                    xr = residualize(q, x)
                    require(np.allclose(contrast @ np.linalg.pinv(xr, rcond=1e-10) @ xr, contrast, atol=1e-7),
                            'Task contrast not estimable after nuisance projection')
                    dof = nvol - q.shape[1] - basis(xr).shape[1] - 7
                    require(dof >= c['min_residual_dof'], f'Insufficient residual dof {dof}')
                    if reference is None:
                        reference = im
                        common = np.asanyarray(brain.dataobj) > 0
                    else:
                        require(same_grid(reference, im), 'Incompatible BOLD geometry across accepted runs')
                        common &= np.asanyarray(brain.dataobj) > 0
                    np.save(out / key / 'nuisance_q.npy', q)
                    row.update(status='included', reason='', tr=tr, nvol=nvol, duration=nvol * tr,
                               residual_dof_lower_bound=dof, **{k:v for k,v in ninfo.items() if k != 'known_columns'})
                    write_json(out / key / 'nuisance.json', ninfo)
                except (ValueError, KeyError, OSError) as exc:
                    row.update(status='blocked', reason=str(exc))
    included = [r for r in rows if r['status'] == 'included']
    blocked = [r for r in rows if r['status'] == 'blocked']
    if blocked:
        blockers.append(f'{len(blocked)} runs cannot be adjudicated/loaded; see manifest.tsv')
    if not included:
        blockers.append('No verified usable 4D runs; production prohibited')
    coverage = {}
    if reference is not None:
        fsldir = Path(os.environ['FSLDIR'])
        atlas = fsldir / 'data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr25-2mm.nii.gz'
        xml = fsldir / 'data/atlases/HarvardOxford-Cortical.xml'
        registry.extend(identity(p) for p in (atlas, xml))
        masks, counts, label = anatomical_masks(atlas, xml, reference)
        save_image(out / 'wholebrain_coverage.nii.gz', common.astype('uint8'), reference)
        for roi in ROIS:
            covered = masks[roi] & common
            lost = masks[roi] & ~common
            xyz = nib.affines.apply_affine(reference.affine, np.argwhere(lost))
            coverage[roi] = dict(counts[roi], final_voxels=int(covered.sum()), label=label,
                retained_fraction=float(covered.sum()/masks[roi].sum()),
                lost_mni_bounds=[xyz.min(axis=0).tolist(), xyz.max(axis=0).tolist()] if len(xyz) else [])
            if coverage[roi]['retained_fraction'] < c['minimum_roi_coverage']:
                blockers.append(f'{roi} ROI coverage below prespecified threshold')
            save_image(out / f'{roi}_mask.nii.gz', covered.astype('uint8'), reference)
            save_image(out / f'{roi}_lost_coverage.nii.gz', lost.astype('uint8'), reference)
    counts = Counter(r['task'] for r in included)
    contribution = []
    for sub in cohort:
        for task in TASKS:
            rr = [r for r in included if r['subject'] == sub and r['task'] == task]
            contribution.append({'subject': sub, 'task': task, 'runs': len(rr), 'seconds': sum(r['duration'] for r in rr)})
            if not rr:
                blockers.append(f'sub-{sub}: no verified {task} run')
    write_table(out / 'manifest.tsv', rows)
    write_table(out / 'contributions.tsv', contribution)
    # Deterministic pilot: span nuisance rank and NSS, retain every run for selected subjects.
    subscore = sorted({r['subject'] for r in included}, key=lambda s: (
        max(r['nuisance_rank'] for r in included if r['subject'] == s),
        max(r['nss'] for r in included if r['subject'] == s), s))
    pilot = [subscore[i] for i in sorted(set(np.linspace(0, len(subscore)-1, min(5, len(subscore))).astype(int)))] if subscore else []
    upstream_commit = subprocess.run(['git', '-C', str(upstream), 'rev-parse', 'HEAD'], capture_output=True, text=True)
    if upstream_commit.returncode:
        blockers.append('Cannot record upstream checkout commit')
    record = {'passed': not blockers, 'blockers': blockers, 'software': sw, 'config': c,
              'pipeline_sha256': sha(Path(__file__)),
              'analysis_commit': subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip(),
              'upstream_commit': upstream_commit.stdout.strip(), 'inputs': registry, 'coverage': coverage,
              'cohort_n': len(cohort), 'run_counts': dict(counts), 'pilot_subjects': pilot,
              'weighting': c['ica_weighting'], 'created_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}
    write_json(out / 'audit.json', record)
    report = [f'# Insula ICA audit\n\nStatus: {"PASS" if not blockers else "BLOCKED"}',
              f'Cohort: {len(cohort)}; verified runs: {dict(counts)}; blocked runs: {len(blocked)}.',
              'Paths are script-derived candidates until status=included; no alternate images are substituted.',
              '## Coverage\n\n' + json.dumps(coverage, indent=2),
              '## Blockers\n\n' + '\n'.join('- ' + b for b in blockers)]
    (out / 'report.md').write_text('\n\n'.join(report) + '\n')
    print(f'{out / "report.md"}: {"PASS" if not blockers else "BLOCKED"}')
    return 0 if not blockers else 2


def checked_audit(c):
    out = Path(c['output_root'])
    audit_path = out / 'audit/audit.json'
    a = json.loads(audit_path.read_text())
    require(a['passed'], 'Audit has blockers; inspect audit/report.md')
    require(a['pipeline_sha256'] == sha(Path(__file__)), 'Pipeline changed since audit; re-audit in a fresh output_root')
    ignored = {'network_templates', 'network_template_notes'}
    require({k:v for k,v in a['config'].items() if k not in ignored} == {k:v for k,v in c.items() if k not in ignored},
            'Analysis config changed since audit; use fresh output_root and re-audit')
    require(a['software'] == software(), 'Software differs from audited environment')
    for item in a['inputs']:
        now = identity(item['path'], content=item['sha256'] is not None)
        require(now == item, f'Input changed since audit: {item["path"]}')
    return out, a, [r for r in read_table(out / 'audit/manifest.tsv') if r['status'] == 'included']


def analysis_dir(out, scope):
    return out / scope


def scoped_rows(rows, a, scope):
    if scope == 'pilot':
        require(len(a['pilot_subjects']) >= 3, 'Pilot needs at least three subjects')
        return [r for r in rows if r['subject'] in a['pilot_subjects']]
    return rows


def load_bold(row):
    image = nib.load(row['bold'])
    return image, image.get_fdata(dtype=np.float32)


def prepare(c):
    out, a, rows = checked_audit(c)
    dest = out / 'prepared'
    require(not dest.exists(), 'Prepared directory exists; choose a fresh output_root for a new analysis')
    dest.mkdir()
    wb = nib.load(out / 'audit/wholebrain_coverage.nii.gz')
    mask = np.asanyarray(wb.dataobj) > 0
    diagnostics = []
    for row in rows:
        im, raw = load_bold(row)
        if not (dest / 'reference_mean.nii.gz').exists():
            save_image(dest / 'reference_mean.nii.gz', raw.mean(axis=3).astype('float32'), im)
        y = raw[mask].T.astype(float)
        require(np.isfinite(y).all() and (y.std(axis=0) > 1e-8).all(), f'Invalid/zero-variance BOLD inside fixed coverage: {row["key"]}')
        q = np.load(out / 'audit' / row['key'] / 'nuisance_q.npy')
        # Scale by run/voxel raw mean: retain an interpretable percent-signal unit.
        mean = y.mean(axis=0)
        require((mean > 0).all(), 'Nonpositive BOLD mean inside coverage')
        clean = residualize(q, 100 * y / mean)
        leak = float(np.linalg.norm(q.T @ clean) / max(np.linalg.norm(clean), 1e-12))
        require(leak < 1e-8, 'Joint projection orthogonality failure')
        data = np.zeros(im.shape, dtype=np.float32)
        data[mask] = clean.T
        path = dest / f'{row["key"]}_clean.nii.gz'
        save_image(path, data, im)
        diagnostics.append({'key': row['key'], 'projection_leakage': leak,
                            'residual_dof': int(row['nvol']) - q.shape[1], 'mean_clean_sd': float(clean.std(axis=0).mean())})
        del raw, y, clean, data
    write_table(dest / 'cleaning_qc.tsv', diagnostics)
    write_json(dest / 'complete.json', {'audit_sha256': sha(out / 'audit/audit.json'), 'n_runs': len(rows)})


def orient_maps(data, mask):
    signs = []
    for k in range(data.shape[3]):
        v = data[..., k][mask]
        sign = 1 if v[np.argmax(np.abs(v))] >= 0 else -1
        data[..., k] *= sign
        signs.append(sign)
    return data, signs


def ica(c, scope):
    out, a, rows = checked_audit(c)
    rows = scoped_rows(rows, a, scope)
    require((out / 'prepared/complete.json').exists(), 'Prepare must complete first')
    dest = analysis_dir(out, scope)
    if scope == 'full':
        require((out / 'pilot/regression_complete.json').exists(), 'Complete pilot ICA + regress first')
    require(not (dest / 'ica_complete.json').exists(), 'ICA already completed')
    dest.mkdir(parents=True, exist_ok=True)
    # One input session per participant/task, so unequal valid-run counts do
    # not create unequal numbers of subject-level PCA blocks. Within each input,
    # normalize each original run/voxel, and weight runs by 1/sqrt(T_run).
    # Disable variance normalization so these within-session weights persist.
    # Standard non-MIGP group MELODIC then whitens each participant/task input.
    require(c['ica_weighting'] == 'participant_task_sessions', 'Unsupported ICA balance policy')
    groups = defaultdict(list)
    for row in rows:
        groups[(row['subject'], row['task'])].append(row)
    weighting = []
    union = np.asanyarray(nib.load(out / 'audit/bilateral_mask.nii.gz').dataobj) > 0
    inputs = []
    discovery = dest / 'discovery'
    discovery.mkdir(exist_ok=True)
    for (subject, task), rr in sorted(groups.items()):
        rr = sorted(rr, key=lambda r: int(r['run']))
        total = sum(int(r['nvol'])-int(r['nss']) for r in rr)
        backing = discovery / 'concat.tmp.dat'
        data = np.memmap(backing, mode='w+', shape=(*union.shape, total), dtype='float32')
        data[:] = 0
        start = 0
        for row in rr:
            im = nib.load(out / 'prepared' / f'{row["key"]}_clean.nii.gz')
            nss = int(row['nss'])
            y = np.asanyarray(im.dataobj)[union, nss:].astype(float)
            y -= y.mean(axis=1, keepdims=True)
            sd = y.std(axis=1, keepdims=True)
            require((sd > 1e-8).all(), 'Zero steady-state variance inside ICA mask')
            n = y.shape[1]
            weight = 1 / np.sqrt(len(rr) * n)
            block = np.zeros((*union.shape, n), dtype=np.float32)
            block[union] = (y / sd * weight).astype('float32')
            data[..., start:start+n] = block
            weighting.append({'key': row['key'], 'participant_task_input': f'sub-{subject}_task-{task}',
                'original_steady_frames': n, 'ica_frames': n, 'within_session_weight': float(weight),
                'start_in_session': start, 'stop_in_session': start+n, 'nss_omitted_from_discovery': nss})
            start += n
        balanced = discovery / f'sub-{subject}_task-{task}.nii.gz'
        save_image(balanced, data, im)
        del data
        backing.unlink()
        inputs.append(str(balanced))
    listing = dest / 'ica_inputs.txt'
    listing.write_text('\n'.join(inputs) + '\n')
    write_table(dest / 'ica_contributions.tsv', weighting)
    for roi in ROIS:
        mask_path = out / f'audit/{roi}_mask.nii.gz'
        mask = np.asanyarray(nib.load(mask_path).dataobj) > 0
        for dim in (5, 3, 7):
            folder = dest / f'{roi}/dim-{dim}'
            require(not folder.exists(), f'Partial ICA directory exists: {folder}; preserve it and use a fresh output_root')
            folder.parent.mkdir(parents=True, exist_ok=True)
            run(['melodic', '-i', listing, '-o', folder, '--tr=' + str(c['expected_tr']), '--nobet',
                 '-a', 'concat', '-m', mask_path, '--bgimage=' + str(out / 'prepared/reference_mean.nii.gz'), '--report', '--Oall', '-d', str(dim),
                 '--disableMigp', '--vn', '--update_mask', '--maxit=2000', '--maxrestart=10'], dest)
            actual_mask = np.asanyarray(nib.load(folder / 'mask.nii.gz').dataobj) > 0
            require(np.array_equal(actual_mask, mask), 'MELODIC changed the fixed anatomical coverage mask')
            image = nib.load(folder / 'melodic_IC.nii.gz')
            data = image.get_fdata()
            require(data.shape[3] == dim and same_grid(image, nib.load(mask_path)), f'Unexpected ICA dimensionality/grid: requested {dim}, returned {data.shape}; inspect convergence' )
            data, signs = orient_maps(data, mask)
            save_image(folder / 'oriented_IC.nii.gz', data.astype('float32'), image)
            write_json(folder / 'signs.json', {'signs': signs, 'rule': 'largest absolute within-mask loading is positive; ties in array order'})
        primary = nib.load(dest / f'{roi}/dim-5/oriented_IC.nii.gz').get_fdata()[mask]
        stability = []
        for dim in (3, 7):
            other = nib.load(dest / f'{roi}/dim-{dim}/oriented_IC.nii.gz').get_fdata()[mask]
            for k in range(5):
                for j in range(dim):
                    stability.append({'primary_component': k+1, 'sensitivity_dim': dim, 'sensitivity_component': j+1,
                                      'spatial_r': float(np.corrcoef(primary[:, k], other[:, j])[0, 1])})
        write_table(dest / roi / 'stability.tsv', stability)
    write_json(dest / 'ica_complete.json', {'audit_sha256': sha(out / 'audit/audit.json'), 'scope': scope, 'n_runs': len(rows)})


def contrast_fit(x, y, contrast, nuisance_rank, min_dof, max_cond):
    """OLS point estimate with Newey-West temporal covariance; estimability checked."""
    pinv = np.linalg.pinv(x, rcond=1e-10)
    require(np.allclose(contrast @ pinv @ x, contrast, atol=1e-7), 'Task contrast not estimable after nuisance projection')
    rank = basis(x).shape[1]
    nonzero_s = linalg.svdvals(x)
    nonzero_s = nonzero_s[nonzero_s > nonzero_s[0] * 1e-10]
    condition = float(nonzero_s[0]/nonzero_s[-1])
    require(condition < max_cond, f'Ill-conditioned task design {condition}')
    dof = len(y) - nuisance_rank - rank
    require(dof >= min_dof, 'Task residual dof below gate')
    beta = pinv @ y
    resid = y - x @ beta
    weights = contrast @ pinv
    lag = max(1, int(np.floor(4 * (len(y)/100)**(2/9))))
    influence = weights[:, None] * resid
    var = np.sum(influence**2, axis=0)
    for j in range(1, lag+1):
        var += 2*(1-j/(lag+1))*np.sum(influence[j:]*influence[:-j], axis=0)
    var *= len(y)/dof
    require((var > 0).all(), 'Nonpositive HAC contrast variance')
    return contrast @ beta, var, {'dof': dof, 'hac_lag': lag, 'condition': condition}


def regress(c, scope):
    out, a, rows = checked_audit(c)
    rows = scoped_rows(rows, a, scope)
    dest = analysis_dir(out, scope)
    require((dest / 'ica_complete.json').exists(), 'ICA must complete first')
    require(not (dest / 'regression_complete.json').exists(), 'Regression already complete')
    wb_image = nib.load(out / 'audit/wholebrain_coverage.nii.gz')
    wb = np.asanyarray(wb_image.dataobj) > 0
    effects, qc = [], []
    for roi in ROIS:
        mask = np.asanyarray(nib.load(out / f'audit/{roi}_mask.nii.gz').dataobj) > 0
        maps = nib.load(dest / f'{roi}/dim-5/oriented_IC.nii.gz').get_fdata()[mask]
        spatial = np.column_stack([np.ones(mask.sum()), maps - maps.mean(axis=0)])
        require(np.linalg.matrix_rank(spatial) == 6, 'Rank deficient spatial ICA basis')
        map_pinv = np.linalg.pinv(spatial)
        for row in rows:
            rd = dest / roi / 'runs' / row['key']
            rd.mkdir(parents=True, exist_ok=True)
            im = nib.load(out / 'prepared' / f'{row["key"]}_clean.nii.gz')
            data = im.get_fdata(dtype=np.float32)
            q = np.load(out / 'audit' / row['key'] / 'nuisance_q.npy')
            x = np.load(out / 'audit' / row['key'] / 'task.npy')
            xr = residualize(q, x)
            # Stage 1: all component maps + spatial intercept simultaneously.
            tc = (map_pinv @ data[mask])[1:].T
            tc = residualize(q, tc)
            scale = tc.std(axis=0, ddof=1)
            require((scale > 1e-8).all(), 'Zero-variance component')
            tc /= scale
            np.savetxt(rd / 'timecourses.tsv', tc, delimiter='\t')
            write_json(rd / 'scales.json', {'raw_timecourse_sd': scale.tolist(), 'unit': 'one run-specific temporal SD'})
            # Stage 2: include nuisance basis explicitly along with projected
            # task EVs and ALL components; FWL avoids rank-deficient zero EVs.
            adjustment = basis(np.column_stack([q, xr]))
            partial_tc = residualize(adjustment, tc)
            require(np.linalg.matrix_rank(partial_tc) == 5, 'Task-adjusted component basis rank deficient')
            condition = float(np.linalg.cond(partial_tc))
            require(condition < c['max_condition_number'], 'Ill-conditioned stage-2 design')
            dof = len(tc) - adjustment.shape[1] - 5
            require(dof >= c['min_residual_dof'], 'Insufficient stage-2 dof')
            np.savetxt(rd / 'stage2_design.tsv', np.column_stack([tc, xr, q]), delimiter='\t')
            y = data[wb].T.astype(float)
            coefficients = np.linalg.pinv(partial_tc) @ residualize(adjustment, y)
            result = np.zeros((*wb.shape, 5), dtype=np.float32)
            result[wb] = coefficients.T
            save_image(rd / 'partial_connectivity.nii.gz', result, wb_image)
            # Same stage-1 TCs, original task columns projected with identical Q.
            original = json.loads((out / 'audit' / row['key'] / 'contrast.json').read_text())
            contrast = np.array(original['original_weights'])
            sign = -1 if row['task'] == 'trust' or c['provenance_review']['ugr_pmods_increase_with_fairness'] else 1
            effect, variance, info = contrast_fit(xr, tc, sign * contrast, q.shape[1], c['min_residual_dof'], c['max_condition_number'])
            for k in range(5):
                effects.append({'roi': roi, 'subject': row['subject'], 'task': row['task'], 'run': row['run'],
                                'component': k+1, 'effect': float(effect[k]), 'variance_hac': float(variance[k]),
                                'sign_multiplier': sign, **info})
            qc.append({'roi': roi, 'key': row['key'], 'stage2_dof': dof, 'stage2_condition': condition,
                       'nuisance_tc_leakage': float(np.linalg.norm(q.T @ tc)/np.linalg.norm(tc)),
                       'task_only_variance_fraction': float(np.linalg.norm(basis(xr).T @ tc)**2 / np.linalg.norm(tc)**2)})
    write_table(dest / 'run_effects.tsv', effects)
    write_table(dest / 'regression_qc.tsv', qc)
    # Equal run means prevent high nominal temporal precision from dominating
    # participant estimates. Group inference uses participants, never runs.
    grouped = defaultdict(list)
    for r in effects:
        grouped[(r['roi'], r['subject'], r['task'], r['component'])].append(r)
    participant = []
    for (roi, sub, task, k), rr in grouped.items():
        participant.append({'roi': roi, 'subject': sub, 'task': task, 'component': k, 'runs': len(rr),
                            'effect': float(np.mean([r['effect'] for r in rr])),
                            'se_hac_independent_runs': float(np.sqrt(sum(r['variance_hac'] for r in rr))/len(rr))})
    write_table(dest / 'participant_effects.tsv', participant)
    write_json(dest / 'regression_complete.json', {'scope': scope, 'audit_sha256': sha(out / 'audit/audit.json'),
                                                 'runs': len(rows), 'subjects': len({r['subject'] for r in rows})})


def bonferroni(p, family=30):
    return np.minimum(np.asarray(p) * family, 1)


def summarize(c, scope):
    out, a, rows = checked_audit(c)
    rows = scoped_rows(rows, a, scope)
    dest = analysis_dir(out, scope)
    require((dest / 'regression_complete.json').exists(), 'Complete regress first')
    templates = c['network_templates']
    require(templates and all(t.get('name') and t.get('citation') and Path(t['path']).is_file() for t in templates),
            'Supply independent reference-network maps and citations before final summary')
    wb_image = nib.load(out / 'audit/wholebrain_coverage.nii.gz')
    wb = np.asanyarray(wb_image.dataobj) > 0
    refmaps = []
    for t in templates:
        image = nib.load(t['path'])
        require(len(image.shape) == 3, 'Network templates must be 3D MNI maps')
        resampled = resample_from_to(image, (wb.shape, wb_image.affine), order=1)
        refmaps.append(resampled.get_fdata())
    write_json(dest / 'network_template_provenance.json', [dict(t, **identity(t['path'])) for t in templates])
    group = []
    participants = read_table(dest / 'participant_effects.tsv')
    for roi in ROIS:
        for task in TASKS:
            for k in range(1, 6):
                vals = np.array([float(r['effect']) for r in participants if r['roi'] == roi and r['task'] == task and int(r['component']) == k])
                require(len(vals) >= 3 and vals.std(ddof=1) > 0, 'Insufficient participant variation for group test')
                effect, se = vals.mean(), stats.sem(vals)
                t = effect / se
                p = float(2 * stats.t.sf(abs(t), len(vals)-1))
                delta = stats.t.ppf(.975, len(vals)-1)*se
                group.append({'roi': roi, 'task': task, 'component': k, 'n': len(vals), 'effect': float(effect),
                              'se': float(se), 'ci_low': float(effect-delta), 'ci_high': float(effect+delta),
                              't': float(t), 'p': p, 'p_bonferroni_30': min(30*p, 1)})
    write_table(dest / 'group_task_effects.tsv', group)
    fingerprints, anatomy = [], []
    for roi in ROIS:
        seed = np.asanyarray(nib.load(out / f'audit/{roi}_mask.nii.gz').dataobj) > 0
        outside = wb & ~seed
        maps = nib.load(dest / f'{roi}/dim-5/oriented_IC.nii.gz').get_fdata()
        for k in range(5):
            v = maps[..., k][seed]
            xyz = nib.affines.apply_affine(wb_image.affine, np.argwhere(seed))
            peak = int(np.argmax(np.abs(v)))
            anatomy.append({'roi': roi, 'component': k+1, 'peak_x': float(xyz[peak, 0]), 'peak_y': float(xyz[peak, 1]),
                            'peak_z': float(xyz[peak, 2]), 'positive_voxels': int((v>0).sum()), 'negative_voxels': int((v<0).sum())})
        # Save participant/task means; tasks have equal weight in pooled subject mean.
        group_sum = np.zeros((*wb.shape, 5), dtype=np.float64)
        participant_count = 0
        folder = dest / roi / 'summary'
        folder.mkdir(exist_ok=True)
        for sub in sorted({r['subject'] for r in rows}):
            tasks = []
            for task in TASKS:
                rr = [r for r in rows if r['subject'] == sub and r['task'] == task]
                require(rr, f'Missing participant/task {sub}/{task}')
                total = np.zeros((*wb.shape, 5), dtype=np.float32)
                for r in rr:
                    total += nib.load(dest / roi / 'runs' / r['key'] / 'partial_connectivity.nii.gz').get_fdata(dtype=np.float32)/len(rr)
                save_image(folder / f'sub-{sub}_task-{task}_mean.nii.gz', total, wb_image)
                tasks.append(total)
            pooled = np.mean(tasks, axis=0)
            save_image(folder / f'sub-{sub}_pooled_mean.nii.gz', pooled, wb_image)
            group_sum += pooled
            participant_count += 1
        mean = (group_sum/participant_count).astype('float32')
        save_image(folder / 'group_descriptive_mean.nii.gz', mean, wb_image)
        for k in range(5):
            for t, refmap in zip(templates, refmaps):
                valid = outside & np.isfinite(refmap)
                require(valid.sum() > 100 and refmap[valid].std() > 0, 'Insufficient reference template variance/coverage')
                fingerprints.append({'roi': roi, 'component': k+1, 'template': t['name'],
                                     'outside_insula_r': float(np.corrcoef(mean[..., k][valid], refmap[valid])[0, 1]),
                                     'n_voxels': int(valid.sum()), 'descriptive_only': True})
    write_table(dest / 'template_correspondence.tsv', fingerprints)
    write_table(dest / 'component_anatomy.tsv', anatomy)
    figures(dest, group, wb_image)
    summary = [f'# Insula ICA results ({scope})',
        f'Participants: {len({r["subject"] for r in rows})}; runs: {len(rows)}. All 15 primary components retained.',
        'Joint nuisance projection preceded ICA. Drift and confounds were modeled together; task EVs were retained for discovery and projected consistently in later models.',
        'Connectivity maps are descriptive partial regression coefficients, adjusted for all five component timecourses, nuisance terms and the original task design. They carry no voxelwise significance claim.',
        'Task tests use participant-level equal-run estimates, two-sided one-sample t tests, and Bonferroni correction across all 30 primary ROI × component × task tests.',
        'TG and UG coefficients have different regressor units. No TG–UG interaction, hemisphere difference, unique psychological encoding, independent replication, or validation of the original gPPI is established.',
        'The ICA basis was estimated on these task data. Effects, including corrected p values, are exploratory. Different component maps alone do not establish statistically different network fingerprints.',
        f'Corrected task tests below .05: {sum(r["p_bonferroni_30"] < .05 for r in group)} / 30. All estimates, including null and ambiguous results, are in group_task_effects.tsv.',
        'Network correspondence is an outside-ROI descriptive correlation with independent templates; weak or competing correlations should remain unlabeled. Inspect all ICA reports, coverage loss, component spectra, nuisance burden and 3/7 correspondence before manuscript interpretation.']
    (dest / 'results_summary.md').write_text('\n\n'.join(summary) + '\n')


def figures(dest, group, reference):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(3, 2, figsize=(11, 10), sharex=True)
    for i, roi in enumerate(ROIS):
        for j, task in enumerate(TASKS):
            rr = [r for r in group if r['roi'] == roi and r['task'] == task]
            ax = axes[i, j]
            vals = np.array([r['effect'] for r in rr])
            ci = np.array([[r['ci_low'], r['ci_high']] for r in rr])
            ax.errorbar(np.arange(1, 6), vals, yerr=np.array([vals-ci[:, 0], ci[:, 1]-vals]), fmt='o', capsize=3)
            ax.axhline(0, color='.5', lw=1)
            ax.set(title=f'{roi.capitalize()} — {task.upper()}', xlabel='Component', ylabel='Task effect (TC SD / EV unit)')
            ax.set_xticks(range(1, 6))
    fig.suptitle('Exploratory component task responses; 95% pointwise participant CIs\nTG: nonreciprocated > reciprocated; UG: increasing unfairness')
    fig.tight_layout(rect=(0, 0, 1, .94))
    fig.savefig(dest / 'task_responses.png', dpi=180)
    plt.close(fig)
    mean_path = dest.parent / 'prepared/reference_mean.nii.gz'
    bg_image = nib.as_closest_canonical(nib.load(mean_path))
    bg = bg_image.get_fdata()
    fig, axes = plt.subplots(15, 2, figsize=(11, 28))
    for i, roi in enumerate(ROIS):
        maps = nib.as_closest_canonical(nib.load(dest / f'{roi}/dim-5/oriented_IC.nii.gz')).get_fdata()
        conn = nib.as_closest_canonical(nib.load(dest / f'{roi}/summary/group_descriptive_mean.nii.gz')).get_fdata()
        for k in range(5):
            for j, data in enumerate((maps[..., k], conn[..., k])):
                z = np.unravel_index(np.argmax(np.abs(data)), data.shape)[2]
                nonzero = np.abs(data[np.abs(data)>0])
                scale = float(np.percentile(nonzero, 98)) if len(nonzero) else 1.
                threshold = float(np.percentile(nonzero, 70)) if len(nonzero) else 0.
                shown = np.ma.masked_where(np.abs(data[:, :, z].T) <= threshold, data[:, :, z].T)
                ax = axes[5*i+k, j]
                ax.imshow(bg[:, :, z].T, origin='lower', cmap='gray', vmin=0, vmax=np.percentile(bg, 99))
                artist = ax.imshow(shown, origin='lower', cmap='RdBu_r', vmin=-scale, vmax=scale, alpha=.85)
                z_mm = nib.affines.apply_affine(bg_image.affine, [0, 0, z])[2]
                ax.set_title(f'{roi} IC{k+1}: {"anatomy" if j==0 else "connectivity (descriptive)"}; z={z_mm:g} mm', fontsize=10)
                ax.text(.01, .05, 'L', transform=ax.transAxes, color='white')
                ax.text(.95, .05, 'R', transform=ax.transAxes, color='white')
                ax.axis('off')
                fig.colorbar(artist, ax=ax, fraction=.03, pad=.02, label='loading' if j==0 else '% signal / TC SD')
    fig.suptitle('All primary components; anatomy and descriptive partial connectivity\nDisplay: largest 30% of nonzero absolute loadings, peak axial slice; no significance threshold', y=1.002)
    fig.tight_layout()
    fig.savefig(dest / 'anatomy_connectivity.png', dpi=140, bbox_inches='tight')
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--scope', choices=('pilot', 'full'), default='pilot')
    p.add_argument('stage', choices=('audit', 'prepare', 'ica', 'regress', 'summarize'))
    args = p.parse_args()
    c = config(args.config)
    try:
        if args.stage == 'audit':
            return audit(c)
        if args.stage == 'prepare':
            prepare(c)
        else:
            globals()[args.stage](c, args.scope)
        return 0
    except (ValueError, OSError, KeyError) as exc:
        print(f'BLOCKED: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
