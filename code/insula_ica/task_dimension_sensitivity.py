#!/usr/bin/env python3
"""Task OLS on existing 3/7-component maps plus ROI means; reuse completed 5D fits."""
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import fcntl
import json
from pathlib import Path
import shutil
import threading
import traceback
import uuid

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from scipy import stats

import run_network_timecourse_analysis as n
import pooled_network_tasks as pooled
from run_bold_controls import save, publish

POLICY = 'fixed-3-7-task-OLS-and-insula-means-reuse-pooled-5D-v1'
FIELDS = ('mode', 'roi', 'dimensionality', 'series', 'component', 'task', 'contrast')


def freeze_json(path, value):
    if path.is_file():
        n.g.core.require(json.loads(path.read_text()) == value, 'Frozen task-sensitivity source changed; use a fresh record')
    else: save(path, value)


def copy_frozen(source, target):
    if target.is_file(): n.g.core.require(n.g.core.sha(source) == n.g.core.sha(target), 'Published ICA map/mask source changed')
    else: shutil.copyfile(source, target)


def spatial_weights(maps):
    maps = np.asarray(maps, dtype=float)
    x = np.column_stack([np.ones(len(maps)), maps-maps.mean(axis=0)])
    n.g.core.require(maps.ndim == 2 and maps.shape[1] in (3, 7) and
                     np.isfinite(x).all() and np.linalg.matrix_rank(x) == maps.shape[1]+1,
                     'Nonestimable spatial component model')
    return np.linalg.pinv(x)[1:]


def column_layout(bases, masks):
    layout = []
    for (mode, roi, dim) in bases:
        for k in range(1, dim+1):
            layout.append(dict(mode=mode, roi=roi, dimensionality=dim, series='component', component=k))
    for roi in masks:
        layout.append(dict(mode='insula_mean', roi=roi, dimensionality=0, series='insula_mean', component=0))
    return layout


def extract_responses(data, bases, masks):
    """Only ROI voxels enter the spatial fits; ROI mean is a separate response."""
    masked = {roi:np.asarray(data[mask], dtype=float) for roi, mask in masks.items()}
    n.g.core.require(all(np.isfinite(y).all() for y in masked.values()), 'Nonfinite BOLD inside an ICA mask')
    arrays = [(item['weights']@masked[roi]).T for (mode, roi, dim), item in bases.items()]
    arrays.extend(masked[roi].mean(axis=0)[:, None] for roi in masks)
    result = np.column_stack(arrays)
    n.g.core.require(result.shape[1] == len(column_layout(bases, masks)) and np.isfinite(result).all(), 'Invalid response layout')
    return result


def preserve_partial(path):
    if path.exists(): path.rename(path.with_name(path.name+'.interrupted-'+uuid.uuid4().hex[:6]))


def cached_extract(folder, signature, loader, bases, masks):
    """An extraction checkpoint survives interrupted temporal fitting."""
    marker = folder/'extraction_complete.json'; target = folder/'stage1_responses.npy'
    if marker.is_file():
        old = json.loads(marker.read_text())
        n.g.core.require(old['signature'] == signature and n.g.core.identity(target) == old['output'],
                         'Cached extraction inputs/output changed; use a fresh record')
        return target, True
    preserve_partial(target)
    data = loader(); responses = extract_responses(data, bases, masks); del data
    temporary = target.with_suffix('.tmp')
    with temporary.open('wb') as stream: np.save(stream, responses)
    temporary.replace(target)
    save(marker, dict(signature=signature, output=n.g.core.identity(target)))
    return target, False


def run_cached_model(row, folder, signature, model, cs, labels, layout, bases, masks, loader):
    """One run: a single 4D read at most, then all component/mean OLS fits jointly."""
    folder.mkdir(parents=True, exist_ok=True)
    marker = folder/'complete.json'
    if marker.is_file():
        old = json.loads(marker.read_text())
        n.g.core.require(old['signature'] == signature and all(n.g.core.identity(p['path']) == p for p in old['outputs']),
                         'Completed task-sensitivity inputs/outputs changed; use a fresh record')
        return old['result'], dict(extraction_reused=True, fit_reused=True)
    target, extraction_reused = cached_extract(folder, signature, loader, bases, masks)
    fit = folder/'fit'; preserve_partial(fit)
    tc = np.load(target)
    estimates, _ = n.fit_timecourses(tc, model, cs, fit)
    effects = []
    for j, label in enumerate(labels):
        for k, item in enumerate(layout):
            effects.append(dict(item, subject=row['subject'], task=row['task'], run=int(row['run']),
                                contrast=label, effect=float(estimates[j, k])))
    outputs = [n.g.core.identity(path) for path in (target, folder/'extraction_complete.json')]
    outputs.extend(n.g.core.identity(path) for path in fit.iterdir() if path.is_file())
    save(marker, dict(signature=signature, result=effects, outputs=outputs))
    return effects, dict(extraction_reused=extraction_reused, fit_reused=False)


def verify_task_source(signature, marker):
    old = json.loads(marker.read_text())['signature']
    n.g.core.require(all(signature[k] == old[k] for k in ('fsf', 'FEAT_input', 'tr')),
                     'Task-preserved BOLD/FSF/TR differs from the completed five-component analysis')
    previous = {v['path']:v for v in old['inputs']}
    n.g.core.require(all(previous.get(v['path']) == v for v in signature['model_inputs']),
                     'Task design/nuisance/mask differs from the completed five-component analysis')


def aggregate(effects, subjects):
    grouped = defaultdict(lambda:defaultdict(dict))
    for row in effects:
        key = tuple(row[k] for k in FIELDS); subject = row['subject']; run = int(row['run'])
        n.g.core.require(run not in grouped[key][subject], 'Duplicate run estimate')
        grouped[key][subject][run] = float(row['effect'])
    participants = []; groups = []
    for key, by_subject in sorted(grouped.items()):
        n.g.core.require(set(by_subject) == set(subjects) and len(subjects) == 132, 'Incomplete manuscript cohort in task sensitivity')
        values = []
        for subject in subjects:
            rr = by_subject[subject]; value = float(np.mean([rr[k] for k in sorted(rr)])); values.append(value)
            participants.append(dict(zip(FIELDS, key), subject=subject, runs=len(rr), effect=value))
        values = np.array(values); mean = float(values.mean()); se = float(stats.sem(values)); delta = float(stats.t.ppf(.975, 131)*se)
        n.g.core.require(np.isfinite(values).all(), 'Nonfinite group effects')
        groups.append(dict(zip(FIELDS, key), n=132, effect=mean, se=se, ci_low=mean-delta, ci_high=mean+delta,
                           inference='Exploratory pointwise participant CI; no permutation test or significance selection'))
    return participants, groups


def verify_five_reuse(groups, source):
    original = n.g.core.read_table(source/'group_network_task_effects.tsv'); error = 0.
    selected = [r for r in groups if r['dimensionality'] == 5]
    n.g.core.require(len(selected) == len(original) == 150, 'Completed five-component results must all be reused')
    for row in selected:
        old = next(r for r in original if all(str(row[k]) == r[k] for k in ('mode', 'roi', 'component', 'task', 'contrast')))
        for field in ('effect', 'se', 'ci_low', 'ci_high'):
            error = max(error, abs(row[field]-float(old[field])))
    n.g.core.require(error < 1e-12, 'Reused five-component summaries differ')
    return dict(group_estimates_reused=150, max_absolute_summary_difference=error,
                source=n.g.core.identity(source/'group_network_task_effects.tsv'), new_five_component_fits=0)


def plot(groups, record):
    colors = ('#68a9cf', '#186da2', '#e6a658', '#bf7340', '#88532e')
    titles = ('UGR nonsocial', 'UGR social', 'Trust computer: rec - def', 'Trust friend: rec - def', 'Trust stranger: rec - def')
    for mode in n.g.inputs.MODES:
        for dim in (3, 5, 7):
            fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharey=True)
            for ax, roi in zip(axes, n.g.core.ROIS):
                values = [r for r in groups if r['roi'] == roi and
                          ((r['mode'] == mode and r['dimensionality'] == dim) or r['series'] == 'insula_mean')]
                for j, (label, title, color) in enumerate(zip(n.BAR_LABELS, titles, colors)):
                    rr = [next(r for r in values if r['component'] == k and r['contrast'] == label) for k in range(dim+1)]
                    y = np.array([r['effect'] for r in rr]); lo = np.array([r['ci_low'] for r in rr]); hi = np.array([r['ci_high'] for r in rr])
                    ax.bar(np.arange(dim+1)+(j-2)*.155, y, width=.145, color=color, label=title, yerr=[y-lo, hi-y], capsize=2)
                ax.axhline(0, color='.4', lw=.9); ax.axvline(.5, color='.65', lw=.8, linestyle=':')
                ax.set(xticks=range(dim+1), xticklabels=['Insula mean']+['IC'+str(k) for k in range(1, dim+1)],
                       ylabel='Standardized coefficient / contrast', title=roi.capitalize()+' insula')
                ax.spines[['top', 'right']].set_visible(False)
            fig.suptitle(f'{mode} discovery basis | {dim} components plus insula-mean reference', fontsize=14)
            handles, labels = axes[0].get_legend_handles_labels()
            fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(.5, .95), ncol=3, frameon=False, fontsize=9)
            fig.text(.5, .025, 'N=132; all 490 accepted runs. UGR pmods pooled before z-scoring; Trust rec-def. Pointwise 95% CIs.', ha='center', fontsize=8.5)
            fig.text(.5, .006, 'Mean is the same independent reference across dimensions/bases, not an extra confound. Five-component fits are reused unchanged.', ha='center', fontsize=8)
            fig.tight_layout(rect=(0, .06, 1, .89)); fig.savefig(record/'figures'/f'{mode}_dim-{dim}_task_bars_with_mean.png', dpi=160); plt.close(fig)
    fig, ax = plt.subplots(figsize=(11, 5))
    for j, (label, title, color) in enumerate(zip(n.BAR_LABELS, titles, colors)):
        rr = [next(r for r in groups if r['roi'] == roi and r['series'] == 'insula_mean' and r['contrast'] == label) for roi in n.g.core.ROIS]
        y = np.array([r['effect'] for r in rr]); lo = np.array([r['ci_low'] for r in rr]); hi = np.array([r['ci_high'] for r in rr])
        ax.bar(np.arange(3)+(j-2)*.155, y, width=.145, color=color, label=title, yerr=[y-lo, hi-y], capsize=3)
    ax.axhline(0, color='.4', lw=.9); ax.set(xticks=range(3), xticklabels=n.g.core.ROIS,
        ylabel='Standardized coefficient / contrast', title='Insula-mean task responses | same masks and full temporal model')
    ax.legend(ncol=3, loc='upper center', bbox_to_anchor=(.5, -.12), frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, .04, 1, 1)); fig.savefig(record/'figures/insula_mean_task_bars.png', dpi=160); plt.close(fig)


def load_maps(pooled_source, record):
    state = json.loads((pooled_source/'state.json').read_text()); base = Path(state['source_record'])
    n.g.core.require(state['state'] == 'pooled_five_bar_review_ready' and state['policy'] == pooled.POLICY and state['accepted_runs'] == 490,
                     'Expected completed pooled five-component task results')
    freeze_json(record/'source-task-provenance.json', dict(source_record=str(pooled_source), policy=state['policy'],
        run_estimates=n.g.core.identity(pooled_source/'run_network_task_effects.tsv'),
        group_estimates=n.g.core.identity(pooled_source/'group_network_task_effects.tsv')))
    base_state = json.loads((base/'state.json').read_text()); exploration = Path(base_state['exploration_record'])
    source_record = Path(json.loads((exploration/'source_snapshot.json').read_text())['source_record'])
    root, rows, fits = n.exploration.snapshot(source_record)
    expected = {(mode, roi, dim) for mode in n.g.inputs.MODES for roi in n.g.core.ROIS for dim in (3, 7)}
    n.g.core.require(all(f'{m}/{r}/dim-{d}' in fits for m, r, d in expected), 'All 12 completed 3/7 ICA solutions required; no ICA will be launched')
    base_rows = json.loads((base/'accepted-inputs.json').read_text())
    n.g.core.require(rows == base_rows, 'ICA and five-component task registries differ')
    reference = nib.load(Path(rows[0]['active_feat'])/'mean_func.nii.gz'); masks = {}; mask_ids = {}; bases = {}; map_ids = {}; comparisons = []
    old_maps = json.loads((base/'map-provenance.json').read_text())
    (record/'primary_maps').mkdir(exist_ok=True)
    for roi in n.g.core.ROIS:
        path = root/'masks'/(roi+'.nii.gz'); image = nib.load(path)
        n.g.core.require(n.g.core.same_grid(image, reference), 'ICA mask grid changed')
        masks[roi] = image.get_fdata() > 0; mask_ids[roi] = n.g.core.identity(path)
        copy_frozen(path, record/'primary_maps'/(roi+'_mask.nii.gz'))
    for mode in n.g.inputs.MODES:
        for roi in n.g.core.ROIS:
            parent = fits[f'{mode}/{roi}/dim-5']; parent_path = Path(parent['oriented_maps'])
            n.g.core.require(n.g.core.identity(parent_path) == old_maps[f'{mode}/{roi}']['map'] and
                             mask_ids[roi] == old_maps[f'{mode}/{roi}']['mask'], 'Five-component source/mask changed')
            parent_maps = nib.load(parent_path).get_fdata()[masks[roi]]
            for dim in (3, 7):
                item = fits[f'{mode}/{roi}/dim-{dim}']; path = Path(item['oriented_maps']); image = nib.load(path); data = image.get_fdata()
                n.g.core.require(image.shape == (*reference.shape, dim) and n.g.core.same_grid(image, reference) and
                                 np.isfinite(data).all() and np.count_nonzero(data[~masks[roi]]) == 0, 'Invalid restored ICA map geometry/support')
                bases[mode, roi, dim] = dict(weights=spatial_weights(data[masks[roi]]))
                map_ids[f'{mode}/{roi}/dim-{dim}'] = n.g.core.identity(path)
                copy_frozen(path, record/'primary_maps'/f'{mode}_{roi}_dim-{dim}.nii.gz')
                comparisons.append(dict(mode=mode, roi=roi, sensitivity_dimension=dim,
                    **n.g.compare_maps(parent_maps, data[masks[roi]], labels=('primary', 'sensitivity'))))
    freeze_json(record/'map-provenance.json', dict(maps=map_ids, masks=mask_ids, ICA_source=str(source_record)))
    freeze_json(record/'spatial_correspondences.json', comparisons)
    freeze_json(record/'accepted-inputs.json', rows); freeze_json(record/'response_layout.json', column_layout(bases, masks))
    return rows, reference, masks, bases, dict(maps=map_ids, masks=mask_ids), base, comparisons


def aligned_comparisons(groups, comparisons):
    result = []
    for item in comparisons:
        for match in item['matches']:
            sign = match['sensitivity_sign_for_primary_comparison']
            for label in n.BAR_LABELS:
                parent = next(r for r in groups if r['mode'] == item['mode'] and r['roi'] == item['roi'] and
                              r['dimensionality'] == 5 and r['component'] == match['primary_component'] and r['contrast'] == label)
                child = next(r for r in groups if r['mode'] == item['mode'] and r['roi'] == item['roi'] and
                             r['dimensionality'] == item['sensitivity_dimension'] and r['component'] == match['sensitivity_component'] and r['contrast'] == label)
                lo, hi = (child['ci_low'], child['ci_high']) if sign == 1 else (-child['ci_high'], -child['ci_low'])
                result.append(dict(mode=item['mode'], roi=item['roi'], primary_component=match['primary_component'],
                    sensitivity_dimension=item['sensitivity_dimension'], sensitivity_component=match['sensitivity_component'],
                    spatial_r=match['spatial_r'], sign_for_comparison=sign, contrast=label, primary_effect=parent['effect'],
                    aligned_sensitivity_effect=sign*child['effect'], aligned_sensitivity_ci_low=lo, aligned_sensitivity_ci_high=hi,
                    note='Spatial matching fixed without task effects; descriptive paired-basis comparison, not independent replication'))
    return result


def execute(args, record, output, state):
    rows, reference, masks, bases, spatial_signature, base, comparisons = load_maps(args.source_record, record)
    layout = column_layout(bases, masks)
    n.g.core.require(len(layout) == 63 and len(rows) == 490, 'Unexpected sensitivity/mean response count')
    software = dict(numpy=np.__version__, scipy=__import__('scipy').__version__, nibabel=nib.__version__,
        implementation_sha256={Path(module.__file__).name:n.g.core.sha(module.__file__)
            for module in (__import__(__name__), n, pooled, n.g.inputs, n.g.core)})
    shared = dict(policy=POLICY, spatial=spatial_signature, software=software, layout=layout)
    save(record/'software.json', software); gate = n.g.inputs.MemoryGate(); lock = threading.Lock()
    base_output = Path(json.loads((base/'state.json').read_text())['output_root'])

    def worker(row):
        key = row['key']; folder = output/'runs'/key
        try:
            design, identities = n.g.inputs.inspect_source(row, reference)
            feat = Path(row['active_feat']); settings = n.g.core.fsl_settings(feat/'design.fsf')
            n.g.core.require(settings['fmri(temphp_yn)'] == '0' and settings['fmri(ndelete)'] == '0', 'Unexpected temporal filtering/deletion')
            bold = n.feat_input(row, settings, reference)
            saved = n.g.core.fsl_matrix(feat/'design.con')
            original_cs, labels = n.contrasts_for(row, settings, design, saved)
            if row['task'] == 'ugr': model, cs, _ = pooled.pool_design(design, settings); labels = list(n.BAR_LABELS[:2])
            else: model, cs = design, original_cs
            selected = [v for v in identities if Path(v['path']).name in ('design.mat', 'design.con', 'confoundevs.txt', 'mask.nii.gz')
                        or v['path'] == str(Path(row['confounds']).resolve())]
            signature = dict(shared, model_inputs=selected, fsf=n.g.core.identity(feat/'design.fsf'),
                             FEAT_input=n.g.core.identity(bold, False), tr=row['tr'])
            verify_task_source(signature, base_output/key/'complete.json')
            def loader():
                image = nib.load(bold); data = np.asarray(image.dataobj, dtype='float32')
                n.g.core.require(data.shape == (*reference.shape, row['nvol']), 'BOLD time/grid changed')
                return data
            with gate.reserve(4):
                result, cache = run_cached_model(row, folder, signature, model, cs, labels, layout, bases, masks, loader)
            n.g.core.require(n.g.core.identity(bold, False) == signature['FEAT_input'] and
                             all(n.g.core.identity(v['path'], v['sha256'] is not None) == v for v in selected), 'Run inputs changed during analysis')
            status = dict(state='complete', **cache, new_component_responses=60, mean_responses=3)
            prior_error = record/'errors'/(key+'.txt')
            if prior_error.exists(): prior_error.rename(record/'attempts'/(key+'-'+uuid.uuid4().hex[:6]+'.txt'))
        except Exception as exc:
            result = None; status = dict(state='failed', error=str(exc)); (record/'errors'/(key+'.txt')).write_text(traceback.format_exc())
        with lock:
            state['runs'][key] = status; save(record/'state.json', state)
            print(f'Task sensitivity {sum(v["state"] == "complete" for v in state["runs"].values())}/490 | {key}: {status["state"]}', flush=True)
        return key, result

    pilots = [rows[0], next(r for r in rows if r['key'] == 'sub-11031_ses-01_task-ugr_run-1'),
              next(r for r in rows if r['key'] == 'sub-10817_ses-01_task-ugr_run-1')]
    results = dict(worker(row) for row in pilots)
    n.g.core.require(all(v is not None for v in results.values()), 'Task-sensitivity pilot failed; no cohort expansion')
    state.update(state='pilot_passed_expanding_task_sensitivity', pilot_passed=True); save(record/'state.json', state)
    if args.publish:
        try: publish(record, 'Publish 3/7-component and insula-mean task pilot; reuse finished 5D results')
        except Exception as exc:
            (record/'pilot-publication-error.txt').write_text(str(exc)); print('Pilot Git publication failed; continuing analysis:', exc, flush=True)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for future in as_completed([pool.submit(worker, row) for row in rows if row not in pilots]):
            key, value = future.result(); results[key] = value
    n.g.core.require(len(results) == 490 and all(v is not None for v in results.values()), 'Some sensitivity runs failed; no incomplete-cohort summary')
    effects = [r for key in sorted(results) for r in results[key]]
    old = n.g.core.read_table(args.source_record/'run_network_task_effects.tsv')
    effects.extend(dict(r, component=int(r['component']), run=int(r['run']), effect=float(r['effect']),
                        dimensionality=5, series='component') for r in old)
    subjects = sorted({r['subject'] for r in rows}); participant, groups = aggregate(effects, subjects)
    n.g.core.require(len(groups) == 465 and sum(r['series'] == 'insula_mean' for r in groups) == 15,
                     'Expected 450 component and 15 unique insula-mean group estimates')
    save(record/'reused-five-component-verification.json', verify_five_reuse(groups, args.source_record))
    flags = [r for r in effects if abs(r['effect']) > 5]
    save(record/'coefficient_review_flags.json', dict(flags=flags, review_threshold=5., exclusions_applied=False))
    for name, values in [('run_task_effects.tsv', effects), ('participant_task_effects.tsv', participant), ('group_task_effects.tsv', groups)]:
        n.g.core.write_table(record/name, values)
    n.g.core.write_table(record/'aligned_task_comparisons.tsv', aligned_comparisons(groups, comparisons))
    plot(groups, record)
    (record/'report.md').write_text('# Task responses across ICA dimensions and insula means\n\n'
        'All 132 manuscript participants and 490 accepted runs are included. '
        'The 3- and 7-component task fits use existing fixed ICA maps; the completed pooled 5-component results are reused unchanged '
        '(see reused-five-component-verification.json). The insula mean is the actual masked input-BOLD average, fitted as a separate response.\n\n'
        'Each series has five bars: pooled nonsocial/social UGR fairness slopes and computer/friend/stranger Trust reciprocate-minus-defect contrasts. '
        'All responses and complete task+nuisance design columns are z-scored within run before simultaneous OLS. '
        'Runs are averaged within participant/task; bars show equal-participant means with pointwise 95% CIs. '
        'These are descriptive estimates, with no permutation or multiple-comparison inference.\n\n'
        '\n'.join(f'- [{path.stem}](figures/{path.name})' for path in sorted((record/'figures').glob('*.png')))+
        '\n\nThe mean reference is shared across bases/dimensions and is not an additional nuisance regressor. '
        'Component numbers are local to each decomposition. Spatial correspondences and sign-aligned descriptive comparisons are saved separately; '
        'matching uses only spatial maps, and every component remains in its own plot.\n')
    state.update(state='task_dimensions_and_mean_complete_review_pending', group_estimates=465,
                 component_bar_estimates=450, mean_bar_estimates=15, five_component_fits_reused=True,
                 coefficient_review_flag_count=len(flags), excluded_participants=0, excluded_runs=0)
    save(record/'state.json', state)


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--source-record', type=Path)
    p.add_argument('--resume-record', type=Path); p.add_argument('--new-record', action='store_true')
    p.add_argument('--jobs', type=int, default=8); p.add_argument('--publish', action='store_true'); args = p.parse_args()
    n.g.core.require(1 <= args.jobs <= 16, 'Use 1-16 RAM-gated 4D readers')
    if args.source_record is None: args.source_record = sorted(n.RECORDS.glob('linux2-pooled-network-tasks-*'))[-1]
    if args.resume_record is None and not args.new_record:
        for candidate in sorted(n.RECORDS.glob('linux2-task-dimensions-*'), reverse=True):
            if not (candidate/'state.json').is_file(): continue
            prior = json.loads((candidate/'state.json').read_text())
            if prior.get('policy') == POLICY and Path(prior['source_record']).resolve() == args.source_record.resolve():
                args.resume_record = candidate; break
    lockdir = n.ROOT/'code/insula_ica/output'; lockdir.mkdir(exist_ok=True)
    with (lockdir/'task-dimension-sensitivity.lock').open('a') as handle:
        try: fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError: print('Task-dimension sensitivity is already running; no second runner started.'); return
        if args.resume_record:
            record = args.resume_record; state = json.loads((record/'state.json').read_text()); args.source_record = Path(state['source_record'])
            output = Path(state['output_root'])
            n.g.core.require(state['policy'] == POLICY and record.resolve().is_relative_to(n.RECORDS.resolve()) and
                             output.resolve().is_relative_to((n.ROOT/'derivatives/insula_ica').resolve()) and
                             output.name.startswith('task-dimensions-'), 'Wrong resume record/root')
        else:
            stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:6]
            record = n.RECORDS/('linux2-task-dimensions-'+stamp); record.mkdir()
            output = n.ROOT/'derivatives/insula_ica'/('task-dimensions-'+stamp); output.mkdir(parents=True)
            state = dict(state='preparing_task_sensitivity', policy=POLICY, source_record=str(args.source_record),
                         output_root=str(output), runs={}, pilot_passed=False)
        for name in ('errors', 'figures', 'attempts'): (record/name).mkdir(exist_ok=True)
        for prior in [record/'failure.txt', *record.glob('*publication-error.txt')]:
            if prior.exists(): prior.rename(record/'attempts'/(prior.stem+'-'+uuid.uuid4().hex[:6]+'.txt'))
        state.pop('error', None); state['jobs'] = args.jobs
        state['last_started_UTC'] = datetime.now(timezone.utc).isoformat(); save(record/'state.json', state)
        (record/'README.md').write_text('Task sensitivity on existing 3/7-component ICA maps; completed pooled 5-component results reused unchanged. '
            'Includes one independently fitted insula-mean response per ROI, shared across dimensions and discovery bases. '
            'Exact masked 4D input to FEAT, columnwise z-scoring, pooled high/low UGR pmods, same full task+nuisance OLS; no extra cleaning/filtering/whitening. '
            'No ICA, FEAT, fMRIPrep, connectivity or randomise jobs. Extraction/fit checkpoints support automatic resume. '
            'All 132 participants/490 accepted runs retained. See ../../README.md for the current workflow.\n')
        try: execute(args, record, output, state)
        except Exception as exc:
            state.update(state='task_dimension_sensitivity_failed', error=str(exc)); save(record/'state.json', state)
            (record/'failure.txt').write_text(traceback.format_exc()); print('ERROR:', exc, flush=True)
        if args.publish:
            try: publish(record, 'Publish 3/5/7 task sensitivity and insula-mean reference without rerunning ICA')
            except Exception as exc: (record/'publication-error.txt').write_text(str(exc)); print('Git publication failed:', exc, flush=True)
        print('Final state:', state['state'], '| record:', record, flush=True)
        raise SystemExit(0 if state['state'] == 'task_dimensions_and_mean_complete_review_pending' else 1)


if __name__ == '__main__': main()
