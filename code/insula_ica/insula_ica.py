#!/usr/bin/env python3
"""Bounded insula ICA: audit -> reconstruct -> verify -> prepare -> ICA -> regress -> summarize.

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
    require(c.get('input_mode') == 'feat_full_model_residuals', 'Legacy pre-ICA nuisance-cleaning configs are retired; use the revised example config')
    for key in ('master_qc', 'output_root', 'derivatives_root'):
        c[key] = str((ROOT / c[key]).resolve()) if not Path(c[key]).is_absolute() else c[key]
    require(c['expected_tr'] > 0, 'TR must be positive')
    require(0 < c['minimum_roi_coverage'] <= 1 and c['min_residual_dof'] > 0, 'Invalid coverage/dof gate')
    require(c['max_condition_number'] > 1, 'Invalid condition-number gate')
    return c


def software():
    needed = ('melodic', 'fsl_glm', 'feat_model', 'feat', 'film_gls')
    require(all(shutil.which(x) for x in needed), 'Source FSL environment: melodic, fsl_glm, feat_model, feat, film_gls required')
    fsldir = Path(os.environ['FSLDIR'])
    return {'fsl': (fsldir / 'etc/fslversion').read_text().strip(),
            'python': sys.version, 'numpy': np.__version__, 'scipy': __import__('scipy').__version__,
            'nibabel': nib.__version__, 'executables': {x: shutil.which(x) for x in needed}}


def analysis_dir(out, scope):
    return out / scope


def scoped_rows(rows, a, scope):
    if scope == 'pilot':
        require(len(a['pilot_subjects']) >= 3, 'Pilot needs at least three subjects')
        return [r for r in rows if r['subject'] in a['pilot_subjects']]
    return rows


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
            im = nib.load(row['residual'])
            nss = 0  # Keep every original residual time sample; no new censoring.
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
        'Discovery and connectivity used the same native full-task-plus-nuisance FEAT residuals. No additional nuisance regression or temporal filter was applied. Original task COPEs supplied a separate task-preserved spatial projection.',
        'Connectivity maps are descriptive partial/background coefficients from simultaneous regression of all five residual component timecourses into the same whole-brain residuals. The existing task+nuisance removal is inherited from FEAT; no EVs were appended again. Native prewhitened residuals, if explicitly approved, characterize innovation-scale connectivity rather than unwhitened BOLD covariance.',
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
            ax.set(title=f'{roi.capitalize()} — {task.upper()}', xlabel='Component', ylabel='Projected task contrast (FEAT units)')
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
                fig.colorbar(artist, ax=ax, fraction=.03, pad=.02, label='loading' if j==0 else 'native residual / TC SD')
    fig.suptitle('All primary components; anatomy and descriptive partial connectivity\nDisplay: largest 30% of nonzero absolute loadings, peak axial slice; no significance threshold', y=1.002)
    fig.tight_layout()
    fig.savefig(dest / 'anatomy_connectivity.png', dpi=140, bbox_inches='tight')
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--scope', choices=('pilot', 'full'), default='pilot')
    p.add_argument('--execute', action='store_true', help='Execute reconstruction; default reconstruct is dry-run only')
    p.add_argument('stage', choices=('audit', 'reconstruct', 'verify', 'prepare', 'ica', 'regress', 'summarize'))
    args = p.parse_args()
    c = config(args.config)
    try:
        if args.stage == 'audit':
            return audit(c)
        if args.stage == 'reconstruct':
            reconstruct(c, args.scope, args.execute)
        elif args.stage == 'verify':
            return verify(c, args.scope)
        elif args.stage == 'prepare':
            prepare(c)
        else:
            globals()[args.stage](c, args.scope)
        return 0
    except (ValueError, OSError, KeyError) as exc:
        print(f'BLOCKED: {exc}', file=sys.stderr)
        return 2


# Ensure the phase-0 module shares this exact module when invoked as a script.
sys.modules.setdefault('insula_ica', sys.modules[__name__])
from feat_residuals import audit, reconstruct, verify, checked_audit, prepare, regress

if __name__ == '__main__':
    raise SystemExit(main())
