#!/usr/bin/env python3
"""Stage-1 network timecourses from task-preserved BOLD, then joint temporal GLS.

Reuses completed ICA maps. Never reads an existing voxelwise task COPE, fits an
ICA, reruns FEAT preprocessing, or uses task-regressed residuals as task data.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import threading
import traceback
import uuid

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from scipy import stats

import explore_completed_ica as exploration
import run_group_ica as g
from run_bold_controls import save, publish

ROOT = g.core.ROOT
RECORDS = ROOT/'code/insula_ica/records'
POLICY = 'network-timecourse-joint-GLS-v1'


def vest(path, matrix):
    matrix = np.atleast_2d(matrix)
    with path.open('w') as stream:
        stream.write(f'/NumWaves {matrix.shape[1]}\n/NumPoints {matrix.shape[0]}\n'
                     f'/NumContrasts {matrix.shape[0]}\n/Matrix\n')
        np.savetxt(stream, matrix, fmt='%.15g')


def contrasts_for(row, settings, design, saved_contrasts):
    """Primary stored contrast plus real condition estimates, never inferred means."""
    ntask = int(settings['fmri(evs_real)'])
    c = saved_contrasts[int(row['contrast_index'])-1].copy()
    g.core.require(len(c) == design.shape[1] and not np.any(c[ntask:]),
                   'Task contrast unexpectedly includes nuisance columns')
    if row['task'] == 'trust':
        expected = np.array([0, 0, 0, -1, 1, -1, 1, -1, 1])
        g.core.require(np.array_equal(c[:9], expected) and not np.any(c[9:]), 'Unexpected Trust contrast')
        for k, title in enumerate(('C_def', 'C_rec', 'F_def', 'F_rec', 'S_def', 'S_rec'), 4):
            g.core.require(settings[f'fmri(evtitle{k})'] == title, 'Trust EV order changed')
        betrayal = np.zeros_like(c); betrayal[[3, 5, 7]] = 1/3
        reciprocity = np.zeros_like(c); reciprocity[[4, 6, 8]] = 1/3
        # Keep the original sum contrast's units. Condition bars average C/F/S;
        # their difference is exactly one third of the primary sum contrast.
        cs = np.stack([-c, betrayal, reciprocity])
        labels = ['betrayal_minus_reciprocity_sum', 'betrayal_mean_CFS', 'reciprocity_mean_CFS']
    else:
        expected = np.array([0, 1, 0, 1, 0, 1, 0, 1, 0])
        g.core.require(np.array_equal(c[:9], expected) and not np.any(c[9:]), 'Unexpected UGR contrast')
        cs = c[None, :]; labels = ['offer_pmod_as_stored']
    centered = design-design.mean(axis=0)
    g.core.require(np.allclose(cs @ np.linalg.pinv(centered) @ centered, cs, atol=1e-6, rtol=1e-6),
                   'Requested network contrast is not estimable in this run')
    return cs, labels


def task_timecourses(projector, bold):
    g.core.require(np.isfinite(bold).all(), 'Nonfinite task-preserved BOLD')
    # Same simultaneous spatial model for all five maps plus spatial intercept.
    # Temporal centering removes the baseline, not task/nuisance effects.
    tc = (projector @ bold).T
    tc -= tc.mean(axis=0)
    g.core.require(tc.shape[1] == 5 and np.isfinite(tc).all() and (tc.std(axis=0) > 1e-8).all(),
                   'Invalid stage-1 timecourses')
    return tc


def fit_timecourses(tc, design, contrasts, folder, tr, film):
    """Each component has its own FILM autocorrelation estimate, with no smoothing."""
    folder.mkdir()
    centered = design-design.mean(axis=0)
    vest(folder/'design.mat', centered); vest(folder/'design.con', contrasts)
    np.savetxt(folder/'stage1_timecourses.tsv', tc, delimiter='\t', fmt='%.12g')
    # Five independent artificial voxels carry five responses, not five
    # predictors. Without --sa there is no spatial autocorrelation pooling.
    image = nib.Nifti1Image(tc.T.reshape(5, 1, 1, len(tc)).astype('float32'), np.eye(4))
    image.header.set_xyzt_units('mm', 'sec'); image.header.set_zooms((1, 1, 1, tr))
    nib.save(image, folder/'network_timecourses.nii.gz')
    argv = [str(film), '--in='+str(folder/'network_timecourses.nii.gz'),
            '--rn='+str(folder/'stats'), '--pd='+str(folder/'design.mat'),
            '--con='+str(folder/'design.con'), '--thr=-1000000000']
    save(folder/'command.json', argv)
    with (folder/'film.log').open('w') as stream:
        subprocess.run(argv, stdout=stream, stderr=subprocess.STDOUT, check=True)
    estimates = []; variances = []
    for k in range(len(contrasts)):
        estimates.append(nib.load(folder/f'stats/cope{k+1}.nii.gz').get_fdata().reshape(5))
        variances.append(nib.load(folder/f'stats/varcope{k+1}.nii.gz').get_fdata().reshape(5))
    estimates = np.stack(estimates); variances = np.stack(variances)
    g.core.require(np.isfinite(estimates).all() and np.isfinite(variances).all() and (variances > 0).all(),
                   'Invalid component-level GLS results')
    # A direct temporal OLS fit is a diagnostic, not a projected voxelwise fit.
    ols = contrasts @ np.linalg.pinv(centered) @ tc
    return estimates, variances, ols


def summarize(rows, subjects):
    participant = []; group = []
    fields = ('mode', 'roi', 'task', 'component', 'contrast')
    keys = sorted({tuple(r[k] for k in fields) for r in rows})
    for key in keys:
        subset = [r for r in rows if tuple(r[k] for k in fields) == key]
        values = []
        for subject in subjects:
            rr = [r for r in subset if r['subject'] == subject]
            g.core.require(bool(rr), 'Missing participant/task network result')
            value = float(np.mean([r['effect'] for r in rr])); values.append(value)
            participant.append(dict(zip(fields, key), subject=subject, runs=len(rr), effect=value))
        values = np.array(values); n = len(values); se = float(stats.sem(values)); mean = float(values.mean())
        g.core.require(n == 132 and se > 0 and np.isfinite(values).all(), 'Incomplete or nonvarying group results')
        delta = float(stats.t.ppf(.975, n-1)*se)
        primary = key[-1] in ('betrayal_minus_reciprocity_sum', 'offer_pmod_as_stored')
        p = float(2*stats.t.sf(abs(mean/se), n-1))
        group.append(dict(zip(fields, key), n=n, effect=mean, se=se, ci_low=mean-delta, ci_high=mean+delta,
                          primary=primary, p=p if primary else '',
                          p_bonferroni_60=min(1., p*60) if primary else '',
                          note='Direct component-timecourse GLS; condition CIs descriptive; same-cohort exploratory ICA.'))
    g.core.require(sum(r['primary'] for r in group) == 60, 'All 60 primary tests required')
    return participant, group


def figures(group, folder):
    for mode in g.inputs.MODES:
        fig, axes = plt.subplots(3, 2, figsize=(12, 10))
        for i, roi in enumerate(g.core.ROIS):
            ax = axes[i, 0]
            for label, shift, color in [('betrayal_mean_CFS', -.18, '#d87535'), ('reciprocity_mean_CFS', .18, '#367bad')]:
                rr = sorted((r for r in group if (r['mode'], r['roi'], r['contrast']) == (mode, roi, label)), key=lambda r:r['component'])
                y = np.array([r['effect'] for r in rr]); lo = np.array([r['ci_low'] for r in rr]); hi = np.array([r['ci_high'] for r in rr])
                ax.bar(np.arange(1, 6)+shift, y, width=.34, color=color, label=label.split('_')[0].capitalize(),
                       yerr=[y-lo, hi-y], capsize=3)
            ax.axhline(0, color='.4', lw=1); ax.set(title=roi+' | Trust conditions (C/F/S mean)', xticks=range(1, 6), ylabel='Network timecourse / task-EV unit')
            ax.legend(frameon=False, fontsize=8)
            rr = sorted((r for r in group if (r['mode'], r['roi'], r['contrast']) == (mode, roi, 'betrayal_minus_reciprocity_sum')), key=lambda r:r['component'])
            ax = axes[i, 1]; y = np.array([r['effect'] for r in rr]); lo = np.array([r['ci_low'] for r in rr]); hi = np.array([r['ci_high'] for r in rr])
            ax.bar(range(1, 6), y, color='#6b658d', yerr=[y-lo, hi-y], capsize=3)
            ax.axhline(0, color='.4', lw=1); ax.set(title=roi+' | Betrayal - reciprocity (C/F/S sum)', xticks=range(1, 6), ylabel='Direct temporal GLS contrast')
            for r in rr:
                if r['p_bonferroni_60'] < .05: ax.annotate('*', (r['component'], r['ci_high']), xytext=(0, 4), textcoords='offset points', ha='center')
        fig.suptitle(f'{mode} ICA basis | Task-preserved BOLD -> stage-1 timecourses -> joint task+nuisance GLS\nN=132; 95% pointwise participant CIs; * Bonferroni across all 60 primary tests', fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, .94)); fig.savefig(folder/f'{mode}_trust_bars.png', dpi=160); plt.close(fig)


def run(args, record, output, state):
    source = json.loads((args.exploration_record/'source_snapshot.json').read_text())
    source_record = Path(source['source_record'])
    root, rows, completed = exploration.snapshot(source_record)
    reference = nib.load(Path(rows[0]['active_feat'])/'mean_func.nii.gz')
    primary = {}; identities = {}
    for mode in g.inputs.MODES:
        for roi in g.core.ROIS:
            item = completed[f'{mode}/{roi}/dim-5']; path = Path(item['oriented_maps'])
            image = nib.load(path); mask_image = nib.load(root/'masks'/(roi+'.nii.gz')); mask = mask_image.get_fdata() > 0
            g.core.require(g.core.same_grid(image, reference) and g.core.same_grid(mask_image, reference), 'Primary map/mask grid changed')
            primary[mode, roi] = dict(mask=mask, weights=exploration.weights(image.get_fdata()[mask]))
            identities[f'{mode}/{roi}'] = dict(map=g.core.identity(path), mask=g.core.identity(root/'masks'/(roi+'.nii.gz')))
    save(record/'map-provenance.json', identities); save(record/'accepted-inputs.json', rows)
    film = Path(os.environ.get('FSLDIR', '/usr/local/fsl'))/'bin/film_gls'
    g.core.require(film.is_file(), 'FSL film_gls is required')
    software = dict(film=g.core.identity(film), fsl=(film.parents[1]/'etc/fslversion').read_text().strip())
    save(record/'software.json', software)
    shared = dict(policy=POLICY, maps=identities, software=software)
    gate = g.inputs.MemoryGate(); lock = threading.Lock()

    def worker(row):
        key = row['key']; folder = output/key; marker = folder/'complete.json'
        try:
            design, inputs = g.inputs.inspect_source(row, reference)
            settings = g.core.fsl_settings(Path(row['active_feat'])/'design.fsf')
            g.core.require(settings['fmri(temphp_yn)'] == '0' and settings['fmri(ndelete)'] == '0', 'Unexpected filtering/deletion')
            cs, labels = contrasts_for(row, settings, design, g.core.fsl_matrix(Path(row['active_feat'])/'design.con'))
            # Existing residuals/PEs need not participate in the task estimator;
            # retain only data, model, nuisance and QC provenance in its signature.
            selected = [v for v in inputs if Path(v['path']).name in ('filtered_func_data.nii.gz', 'design.mat', 'design.con', 'confoundevs.txt', 'mask.nii.gz') or v['path'] == str(Path(row['confounds']).resolve())]
            signature = dict(shared, inputs=selected, fsf=g.core.identity(Path(row['active_feat'])/'design.fsf'), tr=row['tr'])
            if marker.is_file():
                old = json.loads(marker.read_text())
                g.core.require(old['signature'] == signature and all(g.core.identity(v['path']) == v for v in old['outputs']), 'Completed timecourse analysis changed')
                result = old['result']
            else:
                g.core.require(not folder.exists(), 'Partial timecourse fit exists; preserve it and use --retry-failed')
                folder.mkdir(); result = []
                with gate.reserve(4):
                    data = np.asarray(nib.load(row['filtered_data']).dataobj, dtype='float32')
                    g.core.require(data.shape == (*reference.shape, row['nvol']), 'Task-preserved 4D dimensions changed')
                    for (mode, roi), current in primary.items():
                        tc = task_timecourses(current['weights'], data[current['mask']])
                        estimates, variances, ols = fit_timecourses(tc, design, cs, folder/f'{mode}_{roi}', row['tr'], film)
                        for j, label in enumerate(labels):
                            for k in range(5):
                                result.append(dict(mode=mode, roi=roi, subject=row['subject'], task=row['task'], run=row['run'],
                                    component=k+1, contrast=label, effect=float(estimates[j, k]), variance=float(variances[j, k]),
                                    direct_temporal_OLS_diagnostic=float(ols[j, k])))
                    del data
                g.core.require(all(g.core.identity(v['path'], v['sha256'] is not None) == v for v in selected), 'Task inputs changed during fit')
                outputs = [g.core.identity(p) for p in folder.rglob('*') if p.is_file()]
                save(marker, dict(signature=signature, result=result, outputs=outputs))
            status = dict(state='complete', network_GLS_fits=6)
        except Exception as exc:
            status = dict(state='failed', error=str(exc)); result = None
            (record/'errors'/(key+'.txt')).write_text(traceback.format_exc())
        with lock:
            state['runs'][key] = status; save(record/'state.json', state)
            print(f'Network temporal GLS {len(state["runs"])}/490 | {key}: {status["state"]}', flush=True)
        return result

    pilots = [rows[0], next(r for r in rows if r['task'] == 'ugr' and not np.isclose(r['tr'], 1.710)),
              next(r for r in rows if np.isclose(r['tr'], 1.710))]
    state['state'] = 'running_network_timecourse_pilot'; save(record/'state.json', state)
    results = [worker(r) for r in pilots]
    g.core.require(all(r is not None for r in results), 'Pilot failed; no cohort expansion')
    state.update(state='network_timecourse_pilot_passed', pilot_passed=True); save(record/'state.json', state)
    if args.publish: publish(record, 'Publish direct network timecourse GLS pilot and frozen methods')
    if args.pilot_only: return
    remaining = [r for r in rows if r not in pilots]
    state['state'] = 'running_all_network_timecourse_GLS'; save(record/'state.json', state)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results.extend(f.result() for f in as_completed([pool.submit(worker, r) for r in remaining]))
    g.core.require(all(r is not None for r in results), 'Some temporal GLS fits failed; no incomplete-cohort results summarized')
    effects = [r for rr in results for r in rr]
    participants, group = summarize(effects, sorted({r['subject'] for r in rows}))
    g.core.write_table(record/'run_network_task_effects.tsv', effects)
    g.core.write_table(record/'participant_network_task_effects.tsv', participants)
    g.core.write_table(record/'group_network_task_effects.tsv', group)
    figures(group, record/'figures')
    state.update(state='network_timecourse_GLS_complete_review_pending', primary_tests=60)
    save(record/'state.json', state)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exploration-record', type=Path); p.add_argument('--resume-record', type=Path)
    p.add_argument('--retry-failed', action='store_true'); p.add_argument('--jobs', type=int, default=8)
    p.add_argument('--pilot-only', action='store_true'); p.add_argument('--publish', action='store_true')
    args = p.parse_args(); g.core.require(1 <= args.jobs <= 16, 'Use 1-16 RAM-gated 4D readers')
    lockdir = ROOT/'code/insula_ica/output'; lockdir.mkdir(exist_ok=True)
    with (lockdir/'network-timecourse.lock').open('a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if args.resume_record:
            record = args.resume_record; state = json.loads((record/'state.json').read_text())
            args.exploration_record = Path(state['exploration_record']); output = Path(state['output_root'])
            g.core.require(state['policy'] == POLICY and output.resolve().is_relative_to((ROOT/'derivatives/insula_ica').resolve())
                           and output.name.startswith('network-timecourses-'), 'Wrong resume record/root')
            if args.retry_failed:
                for key, value in list(state['runs'].items()):
                    if value['state'] == 'failed' and (output/key).is_dir():
                        # Preserve failed material instead of overwriting or deleting it.
                        (output/key).rename(output/(key+'.failed-'+uuid.uuid4().hex[:6]))
        else:
            if args.exploration_record is None:
                args.exploration_record = sorted(RECORDS.glob('linux2-completed-ICA-exploration-*'))[-1]
            stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:6]
            record = RECORDS/('linux2-network-timecourses-'+stamp); record.mkdir()
            output = ROOT/'derivatives/insula_ica'/('network-timecourses-'+stamp); output.mkdir(parents=True)
            state = dict(state='preparing', policy=POLICY, exploration_record=str(args.exploration_record),
                         output_root=str(output), runs={}, pilot_passed=False, jobs=args.jobs)
        state.pop('error', None)
        for name in ('errors', 'figures'): (record/name).mkdir(exist_ok=True)
        save(record/'state.json', state)
        (record/'README.md').write_text('Network task analysis uses stage-1 spatial-regression timecourses from task-preserved, '
            'unsmoothed filtered BOLD, then direct component-level FILM GLS with the complete saved task+nuisance design. '
            'Both map bases use unwhitened task-preserved BOLD here. No old task maps enter the estimator. '
            'No new ICA, fMRIPrep, FEAT preprocessing, nuisance filtering or PPI is run. '
            'See ../2026-10-29-network-timecourse-methods/README.md for methods.\n')
        try:
            run(args, record, output, state)
        except Exception as exc:
            state.update(state='network_timecourse_analysis_failed', error=str(exc)); save(record/'state.json', state)
            (record/'failure.txt').write_text(traceback.format_exc()); print('ERROR:', exc, flush=True)
        if args.publish:
            try: publish(record, 'Publish direct network timecourse task analysis state and results')
            except Exception as exc:
                (record/'publication-error.txt').write_text(str(exc)); print('Git publication failed:', exc, flush=True)
        print('Final state:', state['state'], '| record:', record, flush=True)
        raise SystemExit(0 if state['state'] in ('network_timecourse_pilot_passed', 'network_timecourse_GLS_complete_review_pending') else 1)


if __name__ == '__main__': main()
