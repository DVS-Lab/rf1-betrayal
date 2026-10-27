#!/usr/bin/env python3
"""Independent fMRIPrep repeat/version controls; never run FEAT or production."""
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

from compare_bold_controls import find_bold, frame_metrics, series_metrics
from diagnose_bold_provenance import Collector
from rebuild_trust import ROOT, filter_for, require, sha, warpkit_inputs


IMAGES = ('fmriprep-24.1.1.simg', 'fmriprep-25.2.4.simg', 'fmriprep-25.2.5.simg',
          'fmriprep-warpkit-andrew.sif', 'fmriprep-warpkit.sif')
ARMS = {'v24-a': '24.1.1', 'v24-b': '24.1.1', 'v25-a': '25.2.5'}
SEED = 20261027
LOCK = threading.Lock()
WORK_BASE = Path('/ZPOOL/data/scratch') / os.environ.get('USER', 'tug87422') / 'rf1-betrayal-bold-controls'


def save(path, value):
    # Unique temporary path also supports concurrent case completion.
    tmp = path.with_name(path.name + '.' + uuid.uuid4().hex + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, default=str, allow_nan=False) + '\n')
    tmp.replace(path)


def now():
    return datetime.now(timezone.utc).isoformat()


def publish(record, message):
    """Commit this record only, preserving unrelated staged/local changes."""
    require(subprocess.check_output(['git', 'branch', '--show-current'], cwd=ROOT, text=True).strip() == 'insula-ica',
            'Publishing requires branch insula-ica')
    subprocess.run(['git', 'fetch', 'origin', 'insula-ica'], cwd=ROOT, check=True, timeout=120)
    subprocess.run(['git', 'merge', '--ff-only', 'origin/insula-ica'], cwd=ROOT, check=True, timeout=120)
    relative = str(record.relative_to(ROOT))
    subprocess.run(['git', 'add', '--', relative], cwd=ROOT, check=True)
    subprocess.run(['git', 'commit', '--only', '-m', message, '--', relative], cwd=ROOT, check=True)
    subprocess.run(['git', 'push', 'origin', 'HEAD:insula-ica'], cwd=ROOT, check=True, timeout=120)


def select_runs(plan):
    pilots = {c['subject'] for c in plan['commands']}
    rows = plan['all_runs']
    existing = 'sub-10402'
    require(existing in pilots, 'Expected existing sub-10402 pilot')
    candidates = sorted({r['subject'] for r in rows if r['subject'] not in pilots and r['warpkit']['ready']})
    # Prefer one additional participant with both accepted runs, chosen before
    # inspecting any new version results. This normally selects sub-10478.
    additional = next((s for s in candidates if sum(r['subject'] == s for r in rows) == 2), None)
    require(additional is not None, 'No additional eligible two-run subject')
    return [r for r in rows if r['subject'] in (existing, additional)]


def command(runtime, image, upstream, case, work, subject, nprocs, mem, controlled=True):
    argv = [runtime, 'run', '--cleanenv',
        '-B', f'{upstream}/bids:/bids:ro', '-B', f'{case}:/control', '-B', f'{work}:/scratch',
        '-B', '/ZPOOL/data/tools/templateflow:/opt/templateflow:ro',
        '-B', f'{case}/mplconfig:/opt/mplconfigdir', '-B', '/ZPOOL/data/tools/licenses:/opts:ro',
        str(image), '/bids', '/control/fmriprep', 'participant',
        '--participant-label', subject.removeprefix('sub-'), '--task-id', 'trust',
        '--bids-filter-file', '/control/filter.json', '--stop-on-first-crash', '--skip-bids-validation',
        '--nthreads', str(nprocs), '--omp-nthreads', '1' if controlled else '8', '--mem-mb', str(mem),
        '--me-t2s-fit-method', 'curvefit',
        '--bold2anat-init', 't1w', '--bold2anat-dof', '6', '--slice-time-ref', '0.5',
        '--me-output-echos', '--output-spaces', 'T1w', 'MNI152NLin6Asym',
        '--fs-no-reconall', '--fs-license-file', '/opts/fs_license.txt', '--notrack', '-w', '/scratch']
    if controlled:
        argv += ['--random-seed', str(SEED), '--skull-strip-fixed-seed']
    return argv


def image_inventory(runtime, tools, record):
    result = []
    snippet = ('import importlib.metadata as m,json; '
               'print(json.dumps({x:m.version(x) for x in '
               '["fmriprep","nipype","templateflow","numpy","nibabel"]}))')
    for name in IMAGES:
        path = tools / name
        info = dict(path=str(path), exists=path.is_file())
        print('Inspecting container:', name, flush=True)
        if path.is_file():
            stat = path.stat()
            info.update(size=stat.st_size, mtime_ns=stat.st_mtime_ns, sha256=sha(path))
            for label, argv in (
                ('labels', [runtime, 'inspect', '--json', str(path)]),
                ('packages', [runtime, 'exec', '--cleanenv', str(path), 'python', '-c', snippet]),
            ):
                try:
                    proc = subprocess.run(argv, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=180)
                    info[label] = dict(returncode=proc.returncode, output=proc.stdout)
                    if label == 'packages' and proc.returncode == 0:
                        info['package_versions'] = next(json.loads(line) for line in reversed(proc.stdout.splitlines()) if line.startswith('{'))
                except (OSError, ValueError, StopIteration, subprocess.TimeoutExpired) as exc:
                    info[label] = dict(error=str(exc))
        result.append(info)
    save(record / 'containers.json', result)
    return result


def input_snapshot(runs, upstream):
    paths = set()
    for row in runs:
        paths.update(Path(row['original_feat']) / name for name in ('example_func.nii.gz', 'mask.nii.gz'))
        if Path(row['new_bold']).is_file():
            paths.add(Path(row['new_bold']))
        paths.update(Path(p) for p in row['raw_images'])
        paths.update(Path(p) for p in row['raw_metadata_sha256'])
        for p, expected in row['raw_images'].items():
            stat = Path(p).stat()
            require(stat.st_size == expected['size'] and stat.st_mtime_ns == expected['mtime_ns'],
                    'Raw image changed since recovery snapshot: ' + p)
        for p, expected in row['raw_metadata_sha256'].items():
            require(sha(p) == expected, 'Raw metadata changed since recovery snapshot: ' + p)
        warpkit = warpkit_inputs(upstream, row['subject'], int(row['run']), [Path(p) for p in row['raw_images']])
        require(warpkit['ready'], 'Warpkit inputs not ready: ' + row['key'])
        paths.update(Path(p) for p in warpkit['files'])
        anatomy = upstream / 'bids' / row['subject'] / 'ses-01/anat'
        require(any(anatomy.glob('*_T1w.nii*')), 'Missing T1w for ' + row['subject'])
        paths.update(p for p in anatomy.iterdir() if p.is_file() and
                     p.name.endswith(('_T1w.nii.gz', '_T2w.nii.gz', '_T1w.json', '_T2w.json')))
    paths.update(p for p in (upstream / 'bids').glob('*.json') if p.is_file())
    # Fingerprint the shared templates available to both versions. The bind is
    # read-only in every control, so no case can change/download shared assets.
    for template in ('MNI152NLin6Asym', 'MNI152NLin2009cAsym', 'OASIS30ANTs'):
        folder = Path('/ZPOOL/data/tools/templateflow') / ('tpl-' + template)
        require(folder.is_dir(), 'Missing template cache: ' + str(folder))
        paths.update(p for p in folder.rglob('*') if p.is_file() and
                     p.name.endswith(('.nii.gz', '.nii', '.json', '.h5', '.txt')))
    paths.add(Path('/ZPOOL/data/tools/licenses/fs_license.txt'))
    return {str(p): dict(size=p.stat().st_size, sha256=sha(p)) for p in sorted(paths)}


def verify_snapshot(snapshot):
    return [p for p, original in snapshot.items() if not Path(p).is_file() or
            Path(p).stat().st_size != original['size'] or sha(p) != original['sha256']]


def launch_case(item, env, record, status):
    print('Starting:', item['name'], flush=True)
    began = now()
    try:
        env = dict(env)
        if item['controlled']:
            for prefix in ('APPTAINERENV_', 'SINGULARITYENV_'):
                env.update({prefix + key: '1' for key in ('ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS',
                    'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')})
        with (record / (item['name'] + '.fmriprep-log.txt')).open('w') as log:
            proc = subprocess.run(item['argv'], env=env, stdout=log, stderr=subprocess.STDOUT)
        result = dict(name=item['name'], started_utc=began, finished_utc=now(), returncode=proc.returncode)
        require(proc.returncode == 0, f'fMRIPrep exited {proc.returncode}; see {item["name"]}.fmriprep-log.txt')
        result['bold_outputs'] = {r['key']: str(find_bold(Path(item['case']), r)) for r in item['runs']}
    except Exception as exc:
        result = dict(name=item['name'], started_utc=began, finished_utc=now(), error=str(exc))
    with LOCK:
        status['cases'][item['name']] = result
        save(record / 'status.json', status)
    print('Finished:', item['name'], 'ERROR' if 'error' in result else 'OK', flush=True)
    return result


def compare(runs, cases, status):
    by_subject = {(c['subject'], c['arm']): c for c in cases}
    result = {}
    for row in runs:
        old = Path(row['original_feat'])
        item = dict(historical_frames={}, full_series={})
        for arm in ARMS:
            case = by_subject[row['subject'], arm]
            state = status['cases'].get(case['name'], {})
            if 'error' in state or 'bold_outputs' not in state:
                item['historical_frames'][arm] = dict(error=state.get('error', 'Case not completed'))
                continue
            bold = Path(state['bold_outputs'][row['key']])
            try:
                item['historical_frames'][arm] = frame_metrics(old / 'example_func.nii.gz', bold,
                    old / 'mask.nii.gz', int(row['nvol']) // 2)
            except Exception as exc:
                item['historical_frames'][arm] = dict(error=str(exc))
        for label, left, right in (('v24_repeat', 'v24-a', 'v24-b'), ('v25_vs_v24', 'v24-a', 'v25-a')):
            a, b = (status['cases'].get(by_subject[row['subject'], arm]['name'], {}) for arm in (left, right))
            if 'bold_outputs' not in a or 'bold_outputs' not in b:
                item['full_series'][label] = dict(error='One or both cases did not complete')
                continue
            try:
                item['full_series'][label] = series_metrics(a['bold_outputs'][row['key']], b['bold_outputs'][row['key']], old / 'mask.nii.gz')
            except Exception as exc:
                item['full_series'][label] = dict(error=str(exc))
        ordinary = by_subject.get((row['subject'], 'v24-default-a'))
        ordinary_b = by_subject.get((row['subject'], 'v24-default-b'))
        if ordinary and ordinary_b:
            state = status['cases'].get(ordinary['name'], {})
            state_b = status['cases'].get(ordinary_b['name'], {})
            try:
                require('bold_outputs' in state and 'bold_outputs' in state_b, 'Default repeat did not complete')
                bold = state['bold_outputs'][row['key']]
                item['historical_frames']['v24-default-a'] = frame_metrics(old / 'example_func.nii.gz', bold,
                    old / 'mask.nii.gz', int(row['nvol']) // 2)
                item['historical_frames']['v24-default-b'] = frame_metrics(old / 'example_func.nii.gz',
                    state_b['bold_outputs'][row['key']], old / 'mask.nii.gz', int(row['nvol']) // 2)
                item['full_series']['v24_default_repeat'] = series_metrics(bold,
                    state_b['bold_outputs'][row['key']], old / 'mask.nii.gz')
                item['full_series']['v24_default_vs_previous_pilot'] = series_metrics(
                    row['new_bold'], bold, old / 'mask.nii.gz')
            except Exception as exc:
                item['full_series']['v24_default_repeat'] = dict(error=str(exc))
                item['full_series']['v24_default_vs_previous_pilot'] = dict(error=str(exc))
        result[row['key']] = item
    return result


def report(results, status):
    lines = ['# Independent Trust fMRIPrep controls', '',
        'Two 24.1.1 runs use separate fresh output/work directories, the same seed, fixed skull-stripping seed and one OpenMP thread. The 25.2.5 case uses the same requested settings. No FEAT, standalone TEDANA, Warpkit fit, or production expansion is run.', '',
        'These controls measure repeatability under explicit controls. They do not measure the full variation of the earlier unseeded, eight-OpenMP-thread configuration or establish the historical version.', '',
        'Two fresh unseeded 24.1.1 runs with eight OpenMP threads repeat the sub-10402 pilot settings. They are compared with each other and with the retained pilot BOLD. This is one empirical default-setting pair, not an estimate of the full distribution of possible reruns.', '',
        f'Status: {status["state"]}; {len(status.get("input_changes", []))} input changes detected.', '',
        '## Historical midpoint frames', '',
        '| Run | 24 A relative RMS | 24 B relative RMS | 25 A relative RMS |',
        '|---|---:|---:|---:|']
    for key, item in results.items():
        values = [item['historical_frames'].get(arm, {}).get('relative_rms_difference') for arm in ARMS]
        lines.append('| ' + key + ' | ' + ' | '.join(f'{100*v:.5f}%' if v is not None else 'unavailable / different grid' for v in values) + ' |')
    lines += ['', '## All time points within the historical FEAT mask', '',
        '| Run / comparison | Relative RMS | Voxel-demeaned RMS | Median voxel temporal r | Exact equality |',
        '|---|---:|---:|---:|---|']
    for key, item in results.items():
        for label, info in item['full_series'].items():
            if not info.get('compared'):
                lines.append(f'| {key} / {label} | unavailable | unavailable | unavailable | unavailable |')
                continue
            rms, centered = info.get('relative_rms_difference'), info.get('voxel_demeaned_relative_rms')
            median = (info.get('voxel_temporal_correlation_percentiles') or {}).get('median')
            lines.append(f'| {key} / {label} | {rms} | {centered} | {median} | {info["exact_numeric_equality"]} |')
    lines += ['', 'A smaller historical-frame difference means closer agreement for that frame, not proof of the correct historical pipeline or better scientific validity. Grid mismatches are reported without resampling. Numerical equality is reported as evidence, not used to unlock production.', '',
              'The original 24.1.1 pilot used uncontrolled seeds, no fixed skull-strip seed and eight OpenMP threads; strict historical map equality is not an established acceptance standard for that configuration.', '',
              'Case exit status and errors are in status.json. Full metrics are in results.json; runtime labels/packages are in containers.json; configs/methods and known versionless outputs are in provenance.json and sources/.', '']
    for key, item in status['cases'].items():
        if 'error' in item:
            lines.append(f'- {key}: {item["error"]}')
    return '\n'.join(lines) + '\n'


def run(args, record, status):
    root, upstream = args.rebuild_root.resolve(), args.upstream_root.resolve()
    plan = json.loads((root / 'plan-pilot.json').read_text())
    runs = select_runs(plan)
    for row in runs:
        for name in ('example_func.nii.gz', 'mask.nii.gz'):
            require((Path(row['original_feat']) / name).is_file(),
                    'Missing historical comparison input: ' + row['key'] + '/' + name)
    runtime = shutil.which('apptainer') or shutil.which('singularity')
    require(runtime is not None, 'Apptainer/Singularity is not available')
    version_images = {version: args.tools_root / f'fmriprep-{version}.simg' for version in set(ARMS.values())}
    inventories = image_inventory(runtime, args.tools_root, record)
    for version, image in version_images.items():
        info = next(x for x in inventories if x['path'] == str(image))
        require(info.get('package_versions', {}).get('fmriprep') == version,
                'Required image missing or runtime version unverified: ' + str(image))
    original_container = json.loads((root / 'container.json').read_text())
    repeat_container = next(x for x in inventories if x['path'] == str(version_images['24.1.1']))
    require(original_container['sha256'] == repeat_container['sha256'], '24.1.1 image differs from existing pilot image')
    for row in runs:
        if row['subject'] == 'sub-10402':
            require(Path(row['new_bold']).is_file(), 'Retained pilot BOLD is missing for the default repeat: ' + row['key'])
    print('Hashing selected inputs/templates for independent controls.', flush=True)
    snapshot = input_snapshot(runs, upstream)
    snapshot.update({x['path']: dict(size=x['size'], sha256=x['sha256']) for x in inventories
                     if x['path'] in {str(p) for p in version_images.values()}})
    save(record / 'input-sha256.json', snapshot)
    session = record.name
    experiment = ROOT / 'derivatives/insula_ica' / session
    experiment.mkdir(parents=True, exist_ok=False)
    cases = []
    for subject in sorted({r['subject'] for r in runs}):
        selected = [r for r in runs if r['subject'] == subject]
        arms = dict(ARMS)
        if subject == 'sub-10402':
            arms.update({'v24-default-a': '24.1.1', 'v24-default-b': '24.1.1'})
        for arm, version in arms.items():
            name = subject + '-' + arm
            case = experiment / name
            case.mkdir()
            (case / 'mplconfig').mkdir()
            save(case / 'filter.json', filter_for(selected))
            work = WORK_BASE / session / name
            require(not work.exists(), 'Control work directory already exists: ' + str(work))
            work.mkdir(parents=True)
            controlled = not arm.startswith('v24-default')
            argv = command(runtime, version_images[version], upstream, case, work, subject, args.nprocs, args.mem_mb, controlled)
            cases.append(dict(name=name, subject=subject, arm=arm, version=version,
                              controlled=controlled, case=str(case), work=str(work), argv=argv, runs=selected))
    save(record / 'plan.json', dict(created_utc=now(), seed=SEED, omp_nthreads=1,
        max_concurrent_subject_jobs=args.jobs, cases=cases, selection_note='Existing sub-10402 plus first additional eligible two-run subject; no outcome-based selection.',
        code_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()))
    if not args.execute:
        status['state'] = 'planned_only'
        save(record / 'status.json', status)
        return
    env = dict(os.environ)
    for prefix in ('APPTAINERENV_', 'SINGULARITYENV_'):
        env.update({prefix + 'TEMPLATEFLOW_HOME': '/opt/templateflow',
                    prefix + 'MPLCONFIGDIR': '/opt/mplconfigdir'})
    status['state'] = 'running'
    save(record / 'status.json', status)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(launch_case, c, env, record, status) for c in cases]
        for future in as_completed(futures):
            future.result()
    print('Comparing historical frames and complete control time series.', flush=True)
    results = compare(runs, cases, status)
    save(record / 'results.json', results)
    collector = Collector(record)
    provenance = {}
    for subject in sorted({r['subject'] for r in runs}):
        keys = [r['key'] for r in runs if r['subject'] == subject]
        arms = [c['arm'] for c in cases if c['subject'] == subject]
        provenance[subject] = {arm: collector.provenance(experiment / (subject + '-' + arm) / 'fmriprep', subject, keys)
                               for arm in arms}
        if subject == 'sub-10402':
            provenance[subject]['previous_pilot'] = collector.provenance(root / 'fmriprep-24', subject, keys)
        # The April 2026 Linux2 wrapper explicitly used this versionless root.
        # Read only the exact same subject/run paths, never a project-wide scan.
        current = upstream / 'derivatives/fmriprep'
        provenance[subject]['existing_versionless_outputs'] = collector.provenance(current, subject, keys)
        for row in (r for r in runs if r['subject'] == subject):
            candidate = current / subject / 'ses-01/func' / Path(row['new_bold']).name
            if candidate.is_file():
                try:
                    results[row['key']]['existing_versionless_frame'] = frame_metrics(
                        Path(row['original_feat']) / 'example_func.nii.gz', candidate,
                        Path(row['original_feat']) / 'mask.nii.gz', int(row['nvol']) // 2)
                except Exception as exc:
                    results[row['key']]['existing_versionless_frame'] = dict(error=str(exc))
    save(record / 'results.json', results)
    save(record / 'provenance.json', dict(subjects=provenance, sources=collector.sources, errors=collector.errors))
    config_checks = []
    for subject, subject_info in provenance.items():
        saved_seeds = {}
        for arm in ARMS:
            configs = subject_info[arm]['configs']
            config = configs[0]['full'] if len(configs) == 1 else {}
            saved_seeds[arm] = config.get('seeds', {})
            matches = (len(configs) == 1 and config.get('seeds', {}).get('master') == SEED and
                       config.get('nipype', {}).get('omp_nthreads') == 1 and
                       config.get('workflow', {}).get('skull_strip_fixed_seed') is True)
            config_checks.append(dict(subject=subject, arm=arm, saved_config_matches_controls=matches,
                                      actual_seeds=saved_seeds[arm]))
        for arm in ('v24-default-a', 'v24-default-b'):
            if arm not in subject_info:
                continue
            configs = subject_info[arm]['configs']
            config = configs[0]['full'] if len(configs) == 1 else {}
            matches = (len(configs) == 1 and config.get('nipype', {}).get('omp_nthreads') == 8 and
                       config.get('workflow', {}).get('skull_strip_fixed_seed') is False)
            config_checks.append(dict(subject=subject, arm=arm, saved_config_matches_controls=matches,
                                      actual_seeds=config.get('seeds', {})))
        config_checks.append(dict(subject=subject, repeat_seed_identity=saved_seeds['v24-a'] == saved_seeds['v24-b']))
    save(record / 'control-settings-checks.json', config_checks)
    invalid_config = any(c.get('saved_config_matches_controls') is False or c.get('repeat_seed_identity') is False for c in config_checks)
    status['input_changes'] = verify_snapshot(snapshot)
    failures = any('error' in x for x in status['cases'].values())
    incomplete = any(not x.get('compared') for r in results.values() for x in r['full_series'].values())
    status['state'] = 'complete_with_errors' if failures or incomplete or invalid_config or status['input_changes'] else 'complete'
    status['finished_utc'] = now()
    save(record / 'status.json', status)
    (record / 'report.md').write_text(report(results, status))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rebuild-root', type=Path, required=True)
    parser.add_argument('--upstream-root', type=Path, default=Path('/ZPOOL/data/projects/rf1-sra-linux2'))
    parser.add_argument('--tools-root', type=Path, default=Path('/ZPOOL/data/tools'))
    parser.add_argument('--jobs', type=int, default=8)
    parser.add_argument('--nprocs', type=int, default=14)
    parser.add_argument('--mem-mb', type=int, default=28000)
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--publish', action='store_true', help='Commit/push only the experiment record on completion, including errors')
    args = parser.parse_args()
    require(0 < args.jobs <= 8 and args.nprocs > 0 and args.mem_mb > 0, 'Invalid resource caps')
    require(args.rebuild_root.resolve().is_relative_to((ROOT / 'derivatives/insula_ica').resolve()), 'Unexpected rebuild root')
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid.uuid4().hex[:6]
    record = ROOT / 'code/insula_ica/records' / ('linux2-bold-controls-' + stamp)
    record.mkdir(parents=True, exist_ok=False)
    status = dict(state='preflight', started_utc=now(), cases={})
    save(record / 'status.json', status)
    print('Record:', record, flush=True)
    exit_code = 0
    try:
        run(args, record, status)
        exit_code = int(status['state'] == 'complete_with_errors')
    except Exception as exc:
        status.update(state='failed', error=str(exc), finished_utc=now())
        save(record / 'status.json', status)
        (record / 'failure.txt').write_text(traceback.format_exc())
        print('ERROR:', exc, flush=True)
        exit_code = 1
    if args.publish:
        try:
            publish(record, 'Record independent Trust fMRIPrep repeat/version controls')
            print('Record pushed to GitHub.', flush=True)
        except Exception as exc:
            (record / 'publication-error.txt').write_text(str(exc) + '\n')
            print('Git publication failed; results remain saved:', exc, flush=True)
            exit_code = 1
    print('Final state:', status['state'], '| record:', record, flush=True)
    raise SystemExit(exit_code)


if __name__ == '__main__':
    main()
