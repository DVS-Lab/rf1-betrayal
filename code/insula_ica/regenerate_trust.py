#!/usr/bin/env python3
"""Prepare new Trust analysis inputs with fMRIPrep 24.1.1 and the old FEAT model.

The user accepted new preprocessing on 2026-10-27. Historical statistic equality
is recorded but is no longer required. Model/confound identity remains required.
No ICA or participant exclusions are performed by this preprocessing driver.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
import threading
import traceback
import uuid

import nibabel as nib
import numpy as np

import insula_ica as core
import rebuild_trust as rebuild
from compare_bold_controls import frame_metrics
from run_bold_controls import publish, save


SEED = 20261027
LOCK = threading.Lock()


def production_command(runtime, image, upstream, output, work, subject):
    argv = rebuild.command(runtime, image, upstream, output, work, subject, 14, 28000)
    # Record a fixed seed while retaining normal multi-thread throughput. Exact
    # repeatability is not claimed for eight OpenMP threads.
    argv += ['--random-seed', str(SEED), '--skull-strip-fixed-seed', '--omp-nthreads', '8', '--notrack']
    return argv


def exact_confounds(row, target):
    preserved = Path(row['confounds'])
    core.require(core.sha(preserved) == row['confounds_sha256'], 'Preserved nuisance bytes changed: ' + row['key'])
    original = Path(row['original_feat']) / 'confoundevs.txt'
    old = np.loadtxt(original if original.is_file() else preserved, ndmin=2)
    new = np.loadtxt(target / 'confoundevs.txt', ndmin=2)
    saved = np.loadtxt(preserved, ndmin=2)
    core.require(old.shape == saved.shape == new.shape and np.array_equal(old, saved) and np.array_equal(saved, new),
                 'Actual original/rebuilt nuisance values differ: ' + row['key'])
    return dict(shape=list(saved.shape), exact_numeric_equality=True, saved_sha256=core.sha(preserved))


def reuse_control(subject, selected, output, containers_sha256):
    """Adopt only a prespecified completed 24.1.1 A arm; never choose by fit."""
    record = core.ROOT / 'code/insula_ica/records/linux2-bold-controls-20261027T113051Z-b1efa8'
    recovered = record / 'recovered-20261027T153848Z-31fc23'
    if not (recovered / 'status.json').is_file():
        return None
    status = json.loads((recovered / 'status.json').read_text())
    name = subject + '-v24-a'
    case_status = status.get('cases', {}).get(name)
    if not case_status or case_status.get('returncode') != 0:
        return None
    inventory = json.loads((record / 'containers.json').read_text())
    core.require(any(x.get('package_versions', {}).get('fmriprep') == '24.1.1' and
                     x.get('sha256') == containers_sha256 for x in inventory), 'Reusable control image differs')
    outputs = case_status['bold_outputs']
    identities = json.loads((recovered / 'bold-output-sha256.json').read_text())
    for row in selected:
        source = Path(outputs[row['key']])
        core.require(source.is_file() and core.sha(source) == identities[str(source)]['sha256'],
                     'Completed control BOLD changed: ' + str(source))
    example = Path(outputs[selected[0]['key']])
    source_subject = example.parents[2].resolve()
    core.require(source_subject.name == subject and source_subject.is_relative_to((core.ROOT / 'derivatives/insula_ica').resolve()),
                 'Unexpected control subject directory')
    destination = output / 'fmriprep-24' / subject
    core.require(not destination.exists() and not destination.is_symlink(), 'Existing incomplete subject output; preserve and inspect: ' + str(destination))
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.symlink_to(source_subject, target_is_directory=True)
    html = source_subject.parent / (subject + '.html')
    if html.is_file() and not (destination.parent / html.name).exists():
        (destination.parent / html.name).symlink_to(html)
    rebuild.verify_outputs(selected)
    return dict(source='completed_v24_a_control', source_record=str(recovered), subject_directory=str(source_subject),
                note='Reused controlled output; seed/thread settings remain in its original execution config.')


def qc_run(row, output):
    bold = nib.load(row['new_bold'])
    target = output / 'feat_reconstruction' / (row['key'] + '.feat')
    residual, mask = nib.load(target / 'stats/res4d.nii.gz'), nib.load(target / 'mask.nii.gz')
    core.require(core.same_grid(bold, residual) and core.same_grid(bold, mask), 'BOLD/residual/mask grid differs')
    core.require(residual.shape == bold.shape, 'Residual time count differs')
    active = mask.get_fdata() > 0
    core.require(active.any(), 'Empty reconstructed FEAT mask')
    values = np.asanyarray(residual.dataobj)[active]
    core.require(np.isfinite(values).all(), 'Nonfinite residual data')
    variable = np.std(values, axis=1) > 0
    core.require(variable.any(), 'Constant residual data')
    active[active] = variable
    del values
    fsldir = Path(os.environ['FSLDIR'])
    roi, counts, label = core.anatomical_masks(
        fsldir / 'data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr25-2mm.nii.gz',
        fsldir / 'data/atlases/HarvardOxford-Cortical.xml', mask)
    coverage = {name: dict(total_voxels=int(m.sum()), retained_voxels=int((m & active).sum()),
                           retained_fraction=float((m & active).sum()/m.sum())) for name, m in roi.items()}
    frame = frame_metrics(Path(row['original_feat']) / 'example_func.nii.gz', row['new_bold'],
                          Path(row['original_feat']) / 'mask.nii.gz', int(row['nvol']) // 2)
    warnings = [name + '_insula_coverage_below_0.8' for name, info in coverage.items() if info['retained_fraction'] < .8]
    if frame.get('spatial_r') is None or frame['spatial_r'] < .95:
        warnings.append('historical_reference_spatial_r_below_0.95_review_alignment')
    report = json.loads((target / 'historical_comparison.json').read_text())
    return dict(key=row['key'], confounds=exact_confounds(row, target), coverage=coverage,
                atlas_label=label, atlas_counts=counts, historical_frame=frame,
                historical_statistics_matched=sum(x['matched'] for x in report['comparisons']),
                historical_statistics_compared=len(report['comparisons']),
                finite_nonconstant_residuals=True, constant_mask_voxels=int((~variable).sum()), review_flags=warnings,
                review_note='Review flags request visual inspection; they do not exclude runs or stop cohort preprocessing.',
                fmriprep_report=str(output / 'fmriprep-24' / (row['subject'] + '.html')),
                residual=str(target / 'stats/res4d.nii.gz'), filtered_data=str(target / 'filtered_func_data.nii.gz'),
                residual_scale='native FILM voxelwise-prewhitened; scale decision remains separate from this input preparation')


def run(args, record, state):
    output = args.output_root.resolve()
    core.require(output.is_relative_to((core.ROOT / 'derivatives/insula_ica').resolve()), 'Expected isolated output root')
    runtime = shutil.which('apptainer') or shutil.which('singularity')
    core.require(runtime is not None and shutil.which('feat'), 'Apptainer/Singularity and the normal FSL setup are required')
    core.require('FSLDIR' in os.environ, 'FSLDIR is not set')
    output.mkdir(parents=True, exist_ok=True)
    lock = output / '.execution.lock'
    lock.mkdir(exist_ok=False)
    rebuild.save_json(lock / 'owner.json', dict(pid=os.getpid(), hostname=os.uname().nodename, driver='analysis_regeneration'))
    try:
        upstream = args.upstream_root.resolve()
        rows = rebuild.load_runs(args.audit_dir.resolve(), upstream, output)
        subjects = sorted({r['subject'] for r in rows})
        core.require(len(subjects) == 132 and len(rows) == 241, 'Expected N132 / 241 accepted Trust runs')
        rebuild.freeze_json(output / 'frozen_runs.json', [{k: v for k, v in r.items() if k != 'warpkit'} for r in rows])
        image = Path('/ZPOOL/data/tools/fmriprep-24.1.1.simg')
        version = subprocess.check_output([runtime, 'run', '--cleanenv', str(image), '--version'], text=True, stderr=subprocess.STDOUT)
        rebuild.record_container(output, image, version)
        container = json.loads((output / 'container.json').read_text())
        decision = dict(date='2026-10-27', policy='new_analysis_preprocessing_original_joint_FEAT_model',
            user_accepted_nonidentical_historical_BOLD=True, historical_statistic_equality_required=False,
            full_design_contrast_and_nuisance_identity_required=True, fmriprep_version='24.1.1',
            seed_for_new_runs=SEED, omp_threads_for_new_runs=8, concurrent_subjects=8,
            completed_outputs_reused_with_original_settings_recorded=True, auto_exclusions=False,
            reason='User accepted regenerated 24.1.1 BOLD with preserved original confounds for the new insula analysis.')
        save(record / 'decision.json', decision)
        rebuild.freeze_json(output / 'analysis_regeneration_policy.json', decision)
        save(record / 'container.json', container)
        save(record / 'input-registry.json', rows)
        plan = []
        filters = output / 'filters'; filters.mkdir(exist_ok=True)
        for subject in subjects:
            selected = [r for r in rows if r['subject'] == subject]
            core.require(all(r['warpkit']['ready'] for r in selected), 'Missing usable Warpkit BIDS inputs: ' + subject)
            save(filters / (subject + '.json'), rebuild.filter_for(selected))
            work = Path('/ZPOOL/data/scratch') / os.environ.get('USER', 'tug87422') / 'rf1-betrayal-trust24-production' / output.name / subject
            argv = production_command(runtime, image, upstream, output, work, subject)
            plan.append(dict(subject=subject, work=str(work), argv=argv, keys=[r['key'] for r in selected]))
        save(record / 'plan.json', plan)
        subprocess.run([runtime, 'exec', '--cleanenv', '-B', f'{output}:/rebuild:ro', str(image),
            'python', '-c', rebuild.FILTER_PREFLIGHT, *[f'/rebuild/filters/{s}.json' for s in subjects]],
            stdout=subprocess.DEVNULL, check=True)
        state.update(state='running', cohort_n=132, accepted_runs=241)
        state['code_commit'] = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=core.ROOT, text=True).strip()
        save(record / 'state.json', state)
        env = dict(os.environ)
        for prefix in ('APPTAINERENV_', 'SINGULARITYENV_'):
            env.update({prefix + 'TEMPLATEFLOW_HOME': '/opt/templateflow', prefix + 'MPLCONFIGDIR': '/opt/mplconfigdir'})

        def launch(item):
            subject = item['subject']
            selected = [r for r in rows if r['subject'] == subject]
            result = dict(subject=subject, runs=item['keys'])
            try:
                marker = output / (subject + '.complete.json')
                if marker.is_file():
                    rebuild.verify_outputs(selected)
                    result['preprocessing'] = 'reused_existing_completed_output'
                else:
                    adopted = reuse_control(subject, selected, output, container['sha256'])
                    if adopted:
                        result['preprocessing'] = adopted
                    else:
                        Path(item['work']).mkdir(parents=True, exist_ok=True)
                        log = output / (subject + '.fmriprep.log')
                        print('Running fMRIPrep:', subject, flush=True)
                        with log.open('a') as stream:
                            stream.write('\n=== New analysis regeneration ' + datetime.now(timezone.utc).isoformat() + ' ===\n')
                            stream.flush()
                            subprocess.run(item['argv'], env=env, stdout=stream, stderr=subprocess.STDOUT, check=True)
                        rebuild.verify_outputs(selected)
                        result['preprocessing'] = 'regenerated_24.1.1_fixed_seed_omp8'
                    save(marker, dict(runs=item['keys'], stage='fmriprep_only', analysis_policy=str(output / 'analysis_regeneration_policy.json'),
                                      provenance=result['preprocessing']))
                summaries = []
                for row in selected:
                    rebuild.fit_run(row, output, strict_historical=False)
                    qc = qc_run(row, output)
                    save(record / 'qc' / (row['key'] + '.json'), qc)
                    summaries.append(qc)
                result.update(state='complete', review_flags={q['key']: q['review_flags'] for q in summaries if q['review_flags']})
                result['preprocessing_marker'] = json.loads(marker.read_text())
            except Exception as exc:
                result.update(state='failed', error=str(exc))
                (record / 'errors' / (subject + '.txt')).write_text(traceback.format_exc())
            with LOCK:
                state['subjects'][subject] = result
                save(record / 'state.json', state)
            print(subject + ': ' + result['state'], flush=True)

        with ThreadPoolExecutor(max_workers=8) as pool:
            futures = [pool.submit(launch, item) for item in plan]
            for future in as_completed(futures): future.result()
        failures = [x for x in state['subjects'].values() if x['state'] == 'failed']
        state.update(state='complete_with_errors' if failures else 'inputs_prepared_qc_review_pending',
                     failed_subjects=len(failures), completed_subjects=132-len(failures),
                     finished_utc=datetime.now(timezone.utc).isoformat())
        save(record / 'state.json', state)
        (record / 'report.md').write_text('# Trust analysis input preparation\n\n' +
            f'Completed subjects: {132-len(failures)}/132; failed: {len(failures)}.\n\n' +
            'Version 24.1.1; original complete task+nuisance FEAT model preserved. Historical statistic differences are informational.\n\n' +
            'QC flags and report paths are in qc/. No runs were silently excluded. Registration/coverage review and residual-scale validation precede ICA.\n\n' +
            '\n'.join(f'- {x["subject"]}: {x["error"]}' for x in failures) + '\n')
    finally:
        (lock / 'owner.json').unlink(missing_ok=True)
        lock.rmdir()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit-dir', type=Path, required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--upstream-root', type=Path, default=Path('/ZPOOL/data/projects/rf1-sra-linux2'))
    parser.add_argument('--publish', action='store_true')
    args = parser.parse_args()
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid.uuid4().hex[:6]
    record = core.ROOT / 'code/insula_ica/records' / ('linux2-trust-regeneration-' + stamp)
    record.mkdir(parents=True); (record / 'qc').mkdir(); (record / 'errors').mkdir()
    state = dict(state='preflight', started_utc=datetime.now(timezone.utc).isoformat(), subjects={})
    print('Record:', record, flush=True)
    save(record / 'state.json', state)
    try:
        run(args, record, state)
    except Exception as exc:
        state.update(state='failed', error=str(exc))
        save(record / 'state.json', state)
        (record / 'failure.txt').write_text(traceback.format_exc())
        print('ERROR:', exc, flush=True)
    if args.publish:
        try:
            publish(record, 'Record Trust analysis regeneration with preserved nuisance models')
        except Exception as exc:
            (record / 'publication-error.txt').write_text(str(exc))
            print('Git publication failed; record saved:', exc, flush=True)
    print('Final state:', state['state'], '| record:', record, flush=True)
    raise SystemExit(0 if state['state'] == 'inputs_prepared_qc_review_pending' else 1)


if __name__ == '__main__':
    main()
