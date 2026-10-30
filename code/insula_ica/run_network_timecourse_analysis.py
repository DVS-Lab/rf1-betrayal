#!/usr/bin/env python3
"""Stage-1 network timecourses from task-preserved BOLD, then standardized joint OLS.

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
POLICY = 'FEAT-input-network-timecourse-zscored-OLS-v2'
BAR_LABELS = ('ugr_nonsocial', 'ugr_social', 'trust_C_rec_minus_def', 'trust_F_rec_minus_def', 'trust_S_rec_minus_def')
BAR_TITLES = ('UGR nonsocial', 'UGR social', 'Trust C: rec - def', 'Trust F: rec - def', 'Trust S: rec - def')


def vest(path, matrix):
    matrix = np.atleast_2d(matrix)
    with path.open('w') as stream:
        stream.write(f'/NumWaves {matrix.shape[1]}\n/NumPoints {matrix.shape[0]}\n'
                     f'/NumContrasts {matrix.shape[0]}\n/Matrix\n')
        np.savetxt(stream, matrix, fmt='%.15g')


def contrasts_for(row, settings, design, saved_contrasts):
    """Weights apply AFTER z-scoring and OLS, as requested for the five bars."""
    ntask = int(settings['fmri(evs_real)'])
    c = saved_contrasts[int(row['contrast_index'])-1]
    g.core.require(len(c) == design.shape[1] and not np.any(c[ntask:]),
                   'Task contrast unexpectedly includes nuisance columns')
    if row['task'] == 'trust':
        expected = np.array([0, 0, 0, -1, 1, -1, 1, -1, 1])
        g.core.require(np.array_equal(c[:9], expected) and not np.any(c[9:]), 'Unexpected Trust contrast')
        for k, title in enumerate(('C_def', 'C_rec', 'F_def', 'F_rec', 'S_def', 'S_rec'), 4):
            g.core.require(settings[f'fmri(evtitle{k})'] == title, 'Trust EV order changed')
        cs = np.zeros((3, design.shape[1]))
        for k, (negative, positive) in enumerate(((3, 4), (5, 6), (7, 8))):
            cs[k, negative] = -1; cs[k, positive] = 1
        labels = list(BAR_LABELS[2:])
    else:
        titles = ('nonsocial_high_constant', 'nonsocial_high_pmod', 'nonsocial_low_constant',
                  'nonsocial_low_pmod', 'social_high_constant', 'social_high_pmod',
                  'social_low_constant', 'social_low_pmod')
        for k, title in enumerate(titles, 1):
            g.core.require(settings[f'fmri(evtitle{k})'] == title, 'UGR EV order changed')
        expected = np.array([0, 1, 0, 1, 0, 1, 0, 1, 0])
        g.core.require(np.array_equal(c[:9], expected) and not np.any(c[9:]), 'Unexpected UGR contrast')
        cs = np.zeros((6, design.shape[1]))
        cs[0, [1, 3]] = .5; cs[1, [5, 7]] = .5
        for k, index in enumerate((1, 3, 5, 7), 2): cs[k, index] = 1
        labels = list(BAR_LABELS[:2])+[titles[i] for i in (1, 3, 5, 7)]
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


def zscore_columns(values, allow_constant=False):
    values = np.asarray(values, dtype=float)
    g.core.require(values.ndim == 2 and len(values) > 2 and np.isfinite(values).all(), 'Invalid data/model matrix')
    means = values.mean(axis=0); centered = values-means
    sd = centered.std(axis=0, ddof=1)
    active = sd > 1e-12*np.maximum(1., np.max(np.abs(values), axis=0))
    g.core.require(allow_constant or active.all(), 'A network timecourse is constant')
    z = np.zeros_like(centered); z[:, active] = centered[:, active]/sd[active]
    return z, means, sd, active


def fit_timecourses(tc, design, contrasts, folder):
    """Z-score each response and every task/confound column; one simultaneous OLS."""
    folder.mkdir()
    zy, ymean, ysd, _ = zscore_columns(tc)
    zx, xmean, xsd, active = zscore_columns(design, allow_constant=True)
    inverse = np.linalg.pinv(zx, rcond=1e-10)
    singular = np.linalg.svd(zx, compute_uv=False)
    g.core.require(np.allclose(contrasts@inverse@zx, contrasts, atol=1e-6, rtol=1e-6),
                   'Requested standardized task contrast is not estimable')
    beta = inverse@zy; estimates = contrasts@beta
    g.core.require(np.isfinite(beta).all(), 'Nonfinite standardized OLS coefficients')
    np.savetxt(folder/'stage1_timecourses.tsv', tc, delimiter='\t', fmt='%.12g')
    np.save(folder/'standardized_timecourses.npy', zy)
    np.save(folder/'standardized_design.npy', zx)
    np.save(folder/'standardized_betas.npy', beta)
    save(folder/'standardization.json', dict(method='columnwise sample SD, ddof=1, within run',
        response_mean=ymean.tolist(), response_SD=ysd.tolist(),
        design_mean=xmean.tolist(), design_SD=xsd.tolist(),
        constant_design_columns_zeroed=np.flatnonzero(~active).tolist(),
        rank=int(np.count_nonzero(singular > singular[0]*1e-10)), columns=design.shape[1],
        pseudoinverse_relative_cutoff=1e-10,
        estimator='pinv(Z(full saved L1 task+nuisance design)) @ Z(stage1 timecourses)',
        temporal_filter=False, nuisance_removed_before_fit=False, whitening=False))
    return estimates, beta


def feat_input(row, settings, reference):
    """Exact 4D FEAT input from active design.fsf, not filtered data or residuals."""
    path = Path(settings['feat_files(1)'])
    g.core.require(path.is_absolute() and path.is_file(), 'Missing saved 4D input to FEAT: '+str(path))
    image = nib.load(path)
    g.core.require(g.core.same_grid(image, reference) and image.shape == (*reference.shape, row['nvol']),
                   'Original FEAT input grid/time count differs from saved model')
    g.core.require(np.isclose(image.header.get_zooms()[3], row['tr'], atol=1e-5), 'Original FEAT input TR differs')
    return path


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
        group.append(dict(zip(fields, key), n=n, effect=mean, se=se, ci_low=mean-delta, ci_high=mean+delta,
                          plotted=key[-1] in BAR_LABELS,
                          note='Standardized joint OLS; pointwise participant CI; exploratory task review.'))
    g.core.require(sum(r['plotted'] for r in group) == 150, 'All 150 planned bar estimates required')
    return participant, group


def figures(group, folder):
    """Exactly five bars per subregion; two UGR averages and three Trust contrasts."""
    colors = ('#68a9cf', '#186da2', '#e6a658', '#bf7340', '#88532e')
    def draw(ax, values, components):
        for j, (label, title, color) in enumerate(zip(BAR_LABELS, BAR_TITLES, colors)):
            rr = [next(r for r in values if r['component'] == k and r['contrast'] == label) for k in components]
            y = np.array([r['effect'] for r in rr]); lo = np.array([r['ci_low'] for r in rr]); hi = np.array([r['ci_high'] for r in rr])
            ax.bar(np.arange(len(components))+(j-2)*.155, y, width=.145, color=color, label=title,
                   yerr=[y-lo, hi-y], capsize=2)
        ax.axhline(0, color='.4', lw=1)
        ax.set(xticks=range(len(components)), xticklabels=['IC'+str(k) for k in components],
               ylabel='Standardized coefficient / contrast', xlabel='Insula component')
        ax.spines[['top', 'right']].set_visible(False)
    for mode in g.inputs.MODES:
        for roi in g.core.ROIS:
            values = [r for r in group if r['mode'] == mode and r['roi'] == roi]
            fig, ax = plt.subplots(figsize=(13, 5.4)); draw(ax, values, range(1, 6))
            ax.set_title(f'{mode} ICA basis | {roi} insula | standardized joint OLS', pad=16)
            ax.legend(ncol=3, loc='upper center', bbox_to_anchor=(.5, -.17), frameon=False, fontsize=9)
            fig.text(.5, .015, 'N=132; equal run means within participant; 95% pointwise participant CIs. UGR: high/low averaged. Trust: reciprocate minus defect.', ha='center', fontsize=8)
            fig.tight_layout(rect=(0, .07, 1, 1)); fig.savefig(folder/f'{mode}_{roi}_five_task_bars.png', dpi=160); plt.close(fig)
            if roi == 'bilateral':
                fig, ax = plt.subplots(figsize=(9.5, 5.4)); draw(ax, values, (3, 4))
                ax.set_title(f'{mode} | bilateral IC3 and IC4 | standardized joint OLS')
                ax.legend(ncol=2, loc='upper center', bbox_to_anchor=(.5, -.17), frameon=False, fontsize=9)
                fig.text(.5, .015, 'N=132; 95% pointwise participant CIs; no connectivity inference or between-task difference test.', ha='center', fontsize=8)
                fig.tight_layout(rect=(0, .07, 1, 1)); fig.savefig(folder/f'{mode}_bilateral_IC3_IC4_five_task_bars.png', dpi=160); plt.close(fig)


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
    software = dict(python=__import__('sys').version, numpy=np.__version__,
                    scipy=__import__('scipy').__version__, nibabel=nib.__version__)
    save(record/'software.json', software)
    shared = dict(policy=POLICY, maps=identities, software=software)
    gate = g.inputs.MemoryGate(); lock = threading.Lock()

    def worker(row):
        key = row['key']; folder = output/key; marker = folder/'complete.json'
        try:
            design, inputs = g.inputs.inspect_source(row, reference)
            settings = g.core.fsl_settings(Path(row['active_feat'])/'design.fsf')
            g.core.require(settings['fmri(temphp_yn)'] == '0' and settings['fmri(ndelete)'] == '0', 'Unexpected filtering/deletion')
            bold_path = feat_input(row, settings, reference)
            cs, labels = contrasts_for(row, settings, design, g.core.fsl_matrix(Path(row['active_feat'])/'design.con'))
            # Existing residuals/PEs need not participate in the task estimator;
            # retain only data, model, nuisance and QC provenance in its signature.
            selected = [v for v in inputs if Path(v['path']).name in ('filtered_func_data.nii.gz', 'design.mat', 'design.con', 'confoundevs.txt', 'mask.nii.gz') or v['path'] == str(Path(row['confounds']).resolve())]
            signature = dict(shared, inputs=selected, fsf=g.core.identity(Path(row['active_feat'])/'design.fsf'), tr=row['tr'], FEAT_input=g.core.identity(bold_path, False))
            if marker.is_file():
                old = json.loads(marker.read_text())
                g.core.require(old['signature'] == signature and all(g.core.identity(v['path']) == v for v in old['outputs']), 'Completed timecourse analysis changed')
                result = old['result']
            else:
                g.core.require(not folder.exists(), 'Partial timecourse fit exists; preserve it and use --retry-failed')
                folder.mkdir(); result = []
                with gate.reserve(4):
                    data = np.asarray(nib.load(bold_path).dataobj, dtype='float32')
                    g.core.require(data.shape == (*reference.shape, row['nvol']), 'Task-preserved 4D dimensions changed')
                    for (mode, roi), current in primary.items():
                        tc = task_timecourses(current['weights'], data[current['mask']])
                        estimates, _ = fit_timecourses(tc, design, cs, folder/f'{mode}_{roi}')
                        for j, label in enumerate(labels):
                            for k in range(5):
                                result.append(dict(mode=mode, roi=roi, subject=row['subject'], task=row['task'], run=row['run'],
                                    component=k+1, contrast=label, effect=float(estimates[j, k]),
                                    plotted=label in BAR_LABELS, source_FEAT_input=str(bold_path)))
                    del data
                g.core.require(g.core.identity(bold_path, False) == signature['FEAT_input'], 'FEAT input changed during fit')
                g.core.require(all(g.core.identity(v['path'], v['sha256'] is not None) == v for v in selected), 'Task inputs changed during fit')
                outputs = [g.core.identity(p) for p in folder.rglob('*') if p.is_file()]
                save(marker, dict(signature=signature, result=result, outputs=outputs))
            status = dict(state='complete', network_OLS_fits=6)
        except Exception as exc:
            status = dict(state='failed', error=str(exc)); result = None
            (record/'errors'/(key+'.txt')).write_text(traceback.format_exc())
        with lock:
            state['runs'][key] = status; save(record/'state.json', state)
            print(f'Network standardized OLS {len(state["runs"])}/490 | {key}: {status["state"]}', flush=True)
        return result

    pilots = [rows[0], next(r for r in rows if r['task'] == 'ugr' and not np.isclose(r['tr'], 1.710)),
              next(r for r in rows if np.isclose(r['tr'], 1.710))]
    state['state'] = 'running_network_timecourse_pilot'; save(record/'state.json', state)
    results = [worker(r) for r in pilots]
    g.core.require(all(r is not None for r in results), 'Pilot failed; no cohort expansion')
    state.update(state='network_timecourse_pilot_passed', pilot_passed=True); save(record/'state.json', state)
    if args.publish: publish(record, 'Publish direct standardized OLS task pilot and frozen methods')
    if args.pilot_only: return
    remaining = [r for r in rows if r not in pilots]
    state['state'] = 'running_all_standardized_network_OLS'; save(record/'state.json', state)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results.extend(f.result() for f in as_completed([pool.submit(worker, r) for r in remaining]))
    g.core.require(all(r is not None for r in results), 'Some standardized OLS fits failed; no incomplete-cohort results summarized')
    effects = [r for rr in results for r in rr]
    participants, group = summarize(effects, sorted({r['subject'] for r in rows}))
    g.core.write_table(record/'run_network_task_effects.tsv', effects)
    g.core.write_table(record/'participant_network_task_effects.tsv', participants)
    g.core.write_table(record/'group_network_task_effects.tsv', group)
    figures(group, record/'figures')
    save(record/'bar-layout.json', dict(labels=BAR_LABELS, titles=BAR_TITLES,
        standardization='Within-run columnwise z-scores (ddof=1) for every task/confound column and each component response',
        estimator='simultaneous OLS; no whitening; no prior task/nuisance removal',
        input='Exact 4D feat_files(1) from active saved design.fsf, restricted to the ICA ROI before spatial regression',
        group='Equal run means within each participant/task; N=132; pointwise 95% CIs',
        inference='Descriptive priority review; no randomise or multiple-testing significance labels.'))
    state.update(state='standardized_task_bars_complete_review_pending', bar_estimates=150)
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
            '4D input to FEAT, then joint OLS with z-scored responses and every column of the complete saved task+nuisance design. '
            'Five bars per component: UGR nonsocial/social pmods averaged over high/low; Trust C/F/S reciprocate minus defect. '
            'No connectivity, randomise, ICA, fMRIPrep, FEAT preprocessing, FILM, nuisance filtering or PPI is run. '
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
        raise SystemExit(0 if state['state'] in ('network_timecourse_pilot_passed', 'standardized_task_bars_complete_review_pending') else 1)


if __name__ == '__main__': main()
