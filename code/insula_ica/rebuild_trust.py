#!/usr/bin/env python3
"""Regenerate missing Trust BOLD with archived cluster fMRIPrep settings.

The saved FEAT confound matrix, not a new TEDANA decomposition, remains the
authority for reconstructing the historical task+nuisance model. Outputs and
working directories are isolated from both upstream and historical FEATs.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess

import nibabel as nib
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
VERSION = '24.1.1'
HISTORICAL_SOURCE = ('DVS-Lab/rf1-sra-linux2@a19c46ef:code/fmriprep24-hpc.sh; '
                     'code/ses-01.json; code/tedana-hpc.sh; code/warpkit-hpc.sh')

# Exercise the installed container's own JSON converters before starting jobs.
# Extract just these functions: importing all of config initializes unrelated
# runtime state, whereas these converters require only PyBIDS Query enums.
FILTER_PREFLIGHT = r'''
import ast
import importlib.util
import json
from pathlib import Path
import sys
from bids.layout import Query

root = Path(importlib.util.find_spec('fmriprep').origin).parent
def load_function(path, name, namespace):
    tree = ast.parse(path.read_text())
    matches = [node for node in ast.walk(tree)
               if isinstance(node, ast.FunctionDef) and node.name == name]
    if len(matches) != 1:
        raise RuntimeError('Cannot identify installed filter converter: ' + name)
    exec(compile(ast.Module(body=matches, type_ignores=[]), str(path), 'exec'), namespace)
    return namespace[name]

parser_namespace = {}
load_function(root / 'cli/parser.py', '_process_value', parser_namespace)
json_hook = load_function(root / 'cli/parser.py', '_filter_pybids_none_any', parser_namespace)
convert = load_function(root / 'config.py', '_process_value', {'Query': Query})
for name in sys.argv[1:]:
    filters = json.loads(Path(name).read_text(), object_hook=json_hook)
    for query in filters.values():
        for value in query.values():
            convert(value)
    print('Container filter preflight passed:', name, flush=True)
'''


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def setting(text, key):
    lines = re.findall(r'^set\s+' + re.escape(key) + r'\s+(.+)$', text, re.M)
    require(len(lines) == 1, 'Missing/ambiguous FSF setting: ' + key)
    tokens = shlex.split(lines[0])
    require(len(tokens) == 1, 'Ambiguous FSF value: ' + key)
    return tokens[0]


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def preserve(source, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        require(sha(source) == sha(destination), 'Preserved input changed: ' + str(source))
    else:
        shutil.copy2(source, destination)


def freeze_json(path, value):
    if path.is_file():
        require(json.loads(path.read_text()) == value,
                'Rebuild inputs/settings changed; use a fresh output root: ' + str(path))
    else:
        save_json(path, value)


def parsed_version(output):
    # --version writes a standalone application/version line. Release notices,
    # proxy warnings, and update-check failures are observations, not identity.
    matches = re.findall(r'^\s*(?:fmriprep(?:\s+version)?[\s-]+v?)?'
                         r'(\d+\.\d+\.\d+(?:[A-Za-z0-9.+-]*))\s*$', output, re.I | re.M)
    require(len(set(matches)) == 1, 'Cannot identify a unique fMRIPrep version line: ' + output)
    return matches[0]


def record_container(output, image, version_output):
    version = parsed_version(version_output)
    require(version == VERSION, 'Container is not fMRIPrep ' + VERSION + ': ' + version)
    identity = dict(schema_version=2, path=str(image.resolve()), sha256=sha(image), version=version)
    record = output / 'container.json'
    observation = dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),
                       **identity, version_command_output=version_output)
    with (output / 'container_observations.jsonl').open('a') as stream:
        stream.write(json.dumps(observation) + '\n')
    if record.is_file():
        previous = json.loads(record.read_text())
        require(previous.get('sha256') == identity['sha256'],
                'Container image SHA256 changed; original container.json preserved. '
                'Use the original image or a fresh output root. Current observation: ' +
                str(output / 'container_observations.jsonl'))
        require(parsed_version(previous['version']) == version,
                'Recorded container version differs; original container.json preserved')
        if previous.get('schema_version') != 2:
            preserve(record, output / 'container.legacy.json')
            print('Normalized legacy container record; SHA256 and parsed version match.', flush=True)
    # Only after verifying both scientific identifiers may a legacy record or
    # relocated byte-identical image update the stable record. Keep raw notices
    # in the append-only observations and preserve the original legacy bytes.
    save_json(record, identity)


def subject_jobs(pilot, jobs):
    selected = (2 if pilot else 8) if jobs is None else jobs
    require(selected > 0, 'Invalid concurrent subject count')
    return selected


def warpkit_inputs(upstream, sub, run, raw_inputs):
    """The archived Warpkit script copies its first fieldmap volume into BIDS.

    Presence of raw BOLD echoes does not establish distortion-correction inputs.
    Check only the run-specific products and their actual BOLD associations.
    """
    stem = f'{sub}_ses-01_acq-trust_run-{run}'
    directory = upstream / 'bids' / sub / 'ses-01/fmap'
    paths = [directory / (stem + suffix) for suffix in
             ('_fieldmap.nii.gz', '_magnitude.nii.gz', '_fieldmap.json', '_magnitude.json')]
    problems = ['Missing ' + str(p) for p in paths if not p.is_file()]
    if not problems:
        try:
            metadata = json.loads(paths[2].read_text())
            require(metadata.get('Units') == 'Hz', 'Fieldmap Units must be Hz')
            links = metadata.get('IntendedFor', [])
            if isinstance(links, str):
                links = [links]
            normalized = set()
            for link in links:
                value = link.split(':')[-1].lstrip('/')
                if value.startswith(sub + '/'):
                    value = value[len(sub) + 1:]
                normalized.add(value)
            identifiers = metadata.get('B0FieldIdentifier', [])
            if isinstance(identifiers, str):
                identifiers = [identifiers]
            for raw in (p for p in raw_inputs if p.name.endswith('.nii.gz')):
                rel = raw.relative_to(upstream / 'bids' / sub).as_posix()
                source = json.loads(Path(str(raw)[:-7] + '.json').read_text()).get('B0FieldSource', [])
                if isinstance(source, str):
                    source = [source]
                require(rel in normalized or bool(set(source) & set(identifiers)),
                        'Fieldmap is not associated with ' + rel)
            fieldmap, magnitude = nib.load(paths[0]), nib.load(paths[1])
            require(len(fieldmap.shape) == 3 or (len(fieldmap.shape) == 4 and fieldmap.shape[3] == 1),
                    'BIDS fieldmap must be one volume')
            require(fieldmap.shape[:3] == magnitude.shape[:3] and
                    np.allclose(fieldmap.affine, magnitude.affine, atol=1e-4),
                    'BIDS fieldmap/magnitude geometry differs')
        except (ValueError, KeyError, OSError) as exc:
            problems.append(str(exc))
    retained_4d = (upstream / 'derivatives/warpkit' / sub / 'ses-01' /
                   f'{sub}_ses-01_task-trust_run-{run}_fieldmaps.nii')
    return dict(ready=not problems, problems=problems,
                files={str(p): dict(size=p.stat().st_size, mtime_ns=p.stat().st_mtime_ns,
                                     **({'sha256': sha(p)} if p.suffix == '.json' else {}))
                       for p in paths if p.is_file()},
                retained_4d_fieldmaps=str(retained_4d), retained_4d_exists=retained_4d.is_file())


def load_runs(audit, upstream, output):
    with (audit / 'manifest.tsv').open() as f:
        rows = [r for r in csv.DictReader(f, delimiter='\t')
                if r['task'] == 'trust' and r['status'] != 'excluded_original_qc']
    require(len(rows) == 241, 'Expected all 241 QC-eligible Trust runs')
    require(len({r['key'] for r in rows}) == len(rows), 'Duplicate Trust run')
    result = []
    for row in rows:
        sub = 'sub-' + row['subject'].removeprefix('sub-')
        run = int(row['run'])
        old = Path(row['original_feat'])
        require(old.resolve().is_relative_to((ROOT / 'derivatives/fsl').resolve()),
                'Historical FEAT is outside this checkout')
        fsf = old / 'design.fsf'
        if not fsf.is_file():
            fsf = old.parent / f'L1_{sub}_task-trust_model-01_type-act_run-{run}.fsf'
        text = fsf.read_text()
        bold = Path(setting(text, 'feat_files(1)'))
        expected = upstream / 'derivatives/fmriprep-24' / sub / 'ses-01/func'
        require(bold.parent.resolve() == expected.resolve(), 'Unexpected saved Trust BOLD root')
        conf = Path(setting(text, 'confoundev_files(1)'))
        matrix = conf if conf.is_file() else old / 'confoundevs.txt'
        require(matrix.is_file(), 'Missing original combined nuisance matrix: ' + row['key'])
        nvol = int(setting(text, 'fmri(npts)'))
        require(int(setting(text, 'fmri(ndelete)')) == 0, 'Unexpected volume deletion')
        nuisance = np.loadtxt(matrix, ndmin=2)
        require(nuisance.shape[0] == nvol and np.isfinite(nuisance).all(),
                'Invalid original nuisance matrix: ' + row['key'])
        if conf.is_file() and (old / 'confoundevs.txt').is_file():
            require(np.allclose(nuisance, np.loadtxt(old / 'confoundevs.txt', ndmin=2),
                                atol=1e-6, rtol=1e-6), 'Original nuisance copies disagree')
        prefix = f'{sub}_ses-01_task-trust_run-{run}'
        raw = upstream / 'bids' / sub / 'ses-01/func'
        raw_inputs = []
        for echo in range(1, 5):
            stem = f'{prefix}_echo-{echo}_part-mag_bold'
            for suffix in ('.nii.gz', '.json'):
                path = raw / (stem + suffix)
                if not path.is_file():
                    path = raw / (f'{prefix}_echo-{echo}_bold' + suffix)
                require(path.is_file(), 'Missing raw input: ' + str(path))
                raw_inputs.append(path)
            image = nib.load(raw_inputs[-2])
            require(len(image.shape) == 4 and image.shape[3] == nvol,
                    'Raw/FEAT volume count differs: ' + str(raw_inputs[-2]))
            metadata = json.loads(raw_inputs[-1].read_text())
            require(np.isclose(metadata['RepetitionTime'], float(setting(text, 'fmri(tr)'))),
                    'Raw/FEAT TR differs')
            require(float(metadata['EchoTime']) > 0, 'Missing/invalid echo time')
        require(any((upstream / 'bids' / sub / 'ses-01/anat').glob('*_T1w.nii*')),
                'Missing session-01 T1w: ' + sub)
        saved = output / 'original_confounds' / sub / (prefix + '_desc-TedanaPlusConfounds.tsv')
        preserve(matrix, saved)
        original = output / 'historical_fsf' / (row['key'] + '.fsf')
        preserve(fsf, original)
        result.append(dict(key=row['key'], subject=sub, run=run, nvol=nvol,
                           tr=float(setting(text, 'fmri(tr)')), old_bold=str(bold),
                           new_bold=str(output / 'fmriprep-24' / sub / 'ses-01/func' / bold.name),
                           original_feat=str(old), original_fsf=str(fsf), fsf_sha256=sha(fsf),
                           confounds=str(saved), confounds_sha256=sha(matrix),
                           nuisance_columns=nuisance.shape[1],
                           warpkit=warpkit_inputs(upstream, sub, run, raw_inputs),
                           raw_images={str(p): dict(size=p.stat().st_size, mtime_ns=p.stat().st_mtime_ns)
                                       for p in raw_inputs if p.name.endswith('.nii.gz')},
                           raw_metadata_sha256={str(p): sha(p) for p in raw_inputs if p.suffix == '.json'},
                           tedana_action='reuse_original_combined_FEAT_confounds'))
    return result


def filter_for(runs):
    # Keep the archived anatomical and fieldmap queries. Restrict only functional
    # inputs to session 01, Trust and the originally accepted runs.
    return dict(t1w=dict(datatype='anat', session='01', suffix='T1w'),
                t2w=dict(datatype='anat', session='01', suffix='T2w'),
                fmap=dict(datatype='fmap'),
                **{k: dict(datatype='func', session='01', suffix=k, part=[None, 'mag'],
                           # fMRIPrep 24.1.1's config deserializer applies a
                           # string membership check to each non-Query value.
                           # Integer runs crash it before preprocessing starts.
                           task='trust', run=[str(x) for x in sorted({int(r['run']) for r in runs})])
                   for k in ('sbref', 'bold')})


def command(runtime, image, upstream, output, work, subject, nprocs, mem):
    return [runtime, 'run', '--cleanenv',
            '-B', f'{upstream}/bids:/bids:ro',
            '-B', f'{output}:/rebuild', '-B', f'{work}:/scratch',
            '-B', '/ZPOOL/data/tools/templateflow:/opt/templateflow',
            '-B', '/ZPOOL/data/tools/mplconfigdir:/opt/mplconfigdir',
            '-B', '/ZPOOL/data/tools/licenses:/opts:ro',
            str(image), '/bids', '/rebuild/fmriprep-24', 'participant',
            '--participant-label', subject.removeprefix('sub-'), '--task-id', 'trust',
            '--bids-filter-file', f'/rebuild/filters/{subject}.json',
            '--stop-on-first-crash', '--skip-bids-validation',
            '--nthreads', str(nprocs), '--mem-mb', str(mem), '--me-output-echos',
            '--output-spaces', 'T1w', 'MNI152NLin6Asym', '--fs-no-reconall',
            '--fs-license-file', '/opts/fs_license.txt', '-w', '/scratch']


def recovery_fsf(text, bold, confounds, target):
    for key, value in [('feat_files(1)', bold), ('confoundev_files(1)', confounds),
                       ('fmri(outputdir)', target)]:
        require(not any(c in str(value) for c in ('"', '\n', '\r', '$', '[', ']')), 'Unsupported Tcl path')
        text, n = re.subn(r'^set\s+' + re.escape(key) + r'\s+.*$',
                          lambda _: f'set {key} "{value}"', text, flags=re.M)
        require(n == 1, 'Missing/ambiguous FSF setting: ' + key)
    return text


def compare_fit(run, active):
    # Reconstruction of a historical model is distinct from merely obtaining a
    # successful preprocessing exit status. Never accept new task maps silently.
    import insula_ica as core
    old = Path(run['original_feat'])
    differences = []
    for name in ('design.mat', 'design.con'):
        left, right = core.fsl_matrix(old / name), core.fsl_matrix(active / name)
        require(left.shape == right.shape and np.allclose(left, right, atol=1e-8, rtol=1e-8),
                'Regenerated full design/contrasts differ: ' + run['key'])
    for source in sorted((old / 'stats').glob('*.nii.gz')):
        if not re.fullmatch(r'(cope|zstat)\d+\.nii\.gz', source.name):
            continue
        left, right = nib.load(source), nib.load(active / 'stats' / source.name)
        require(left.shape == right.shape and np.allclose(left.affine, right.affine, atol=1e-5),
                'Regenerated statistic grid differs: ' + run['key'])
        a, b = left.get_fdata(), right.get_fdata()
        matched = bool(np.allclose(a, b, atol=1e-4, rtol=1e-4))
        differences.append(dict(map=source.name, matched=matched, max_abs_difference=float(np.max(np.abs(a-b)))))
    require(differences, 'No historical statistic maps to compare')
    save_json(active / 'historical_comparison.json', differences)
    require(all(d['matched'] for d in differences),
            'Historical COPE/ZSTAT comparison failed; inspect ' + str(active / 'historical_comparison.json') +
            '. New preprocessing is not established as equivalent; do not expand the rerun.')
    residual = nib.load(active / 'stats/res4d.nii.gz')
    require(len(residual.shape) == 4 and residual.shape[-1] == run['nvol'], 'Residual time dimension differs')


def fit_run(run, output):
    target = output / 'feat_reconstruction' / (run['key'] + '.feat')
    fsf = output / 'feat_reconstruction' / (run['key'] + '.fsf')
    fsf.parent.mkdir(exist_ok=True)
    text = recovery_fsf(Path(run['original_fsf']).read_text(), run['new_bold'], run['confounds'], target)
    if fsf.is_file():
        require(fsf.read_text() == text, 'Reconstruction FSF changed')
    else:
        fsf.write_text(text)
    if target.exists():
        require((target / 'stats/res4d.nii.gz').is_file(),
                'Partial FEAT exists; preserve it and inspect rather than relaunch: ' + str(target))
    else:
        with fsf.with_suffix('.log').open('w') as log:
            subprocess.run(['feat', str(fsf)], stdout=log, stderr=subprocess.STDOUT, check=True)
    require((target / 'stats/res4d.nii.gz').is_file(),
            'FEAT has not produced residuals; check the log/scheduler: ' + str(target))
    compare_fit(run, target)


def verify_outputs(runs):
    for r in runs:
        path = Path(r['new_bold'])
        require(path.is_file(), 'Missing regenerated BOLD: ' + str(path))
        new = nib.load(path)
        old = nib.load(Path(r['original_feat']) / 'mean_func.nii.gz')
        require(new.shape == (*old.shape, r['nvol']) and
                np.allclose(new.affine, old.affine, atol=1e-4),
                'Regenerated BOLD has different historical grid/time dimensions: ' + str(path))
        require(np.isclose(new.header.get_zooms()[3], r['tr'], atol=1e-4), 'Regenerated TR differs')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--audit-dir', type=Path, required=True)
    p.add_argument('--upstream-root', type=Path, default=Path('/ZPOOL/data/projects/rf1-sra-linux2'))
    p.add_argument('--output-root', type=Path)
    p.add_argument('--image', type=Path)
    p.add_argument('--fetch-image', action='store_true', help='Pull nipreps/fmriprep:24.1.1 if absent')
    p.add_argument('--pilot', action='store_true', help='One two-run participant and one single-run participant')
    p.add_argument('--jobs', type=int, help='Concurrent subjects (default: 2 for pilot, 8 for production)')
    p.add_argument('--nprocs', type=int, default=14)
    p.add_argument('--mem-mb', type=int, default=28000)
    p.add_argument('--execute', action='store_true')
    p.add_argument('--fit', action='store_true', help='Reconstruct FEAT with original combined confounds and compare every COPE/ZSTAT')
    args = p.parse_args()
    args.jobs = subject_jobs(args.pilot, args.jobs)
    require(args.jobs > 0 and args.nprocs > 0 and args.mem_mb > 0, 'Invalid resource allocation')
    audit = args.audit_dir.resolve()
    upstream = args.upstream_root.resolve()
    output = (args.output_root or ROOT / 'derivatives/insula_ica' /
              ('trust-rebuild-24-' + audit.parent.name)).resolve()
    require(output.is_relative_to((ROOT / 'derivatives/insula_ica').resolve()) and
            output != (ROOT / 'derivatives/insula_ica').resolve(), 'Use a dedicated isolated reconstruction root')
    output.mkdir(parents=True, exist_ok=True)
    runs = load_runs(audit, upstream, output)
    save_json(output / 'warpkit_inventory.json', {r['key']: r['warpkit'] for r in runs})
    print('Warpkit BIDS fieldmaps + magnitude + metadata/links ready:',
          str(sum(r['warpkit']['ready'] for r in runs)) + '/241', flush=True)
    for run in [r for r in runs if not r['warpkit']['ready']][:5]:
        print(run['key'], '; '.join(run['warpkit']['problems']), flush=True)
    subjects = sorted({r['subject'] for r in runs})
    if args.pilot:
        ready_subjects = [s for s in subjects if all(r['warpkit']['ready'] for r in runs if r['subject'] == s)]
        subjects = []
        for count in (2, 1):
            candidates = [s for s in ready_subjects if sum(r['subject'] == s for r in runs) == count]
            require(candidates, f'No Warpkit-ready {count}-run participant for the pilot; see {output / "warpkit_inventory.json"}')
            subjects.append(candidates[0])
    require(all(r['warpkit']['ready'] for r in runs if r['subject'] in subjects),
            'Selected Trust runs lack usable Warpkit BIDS products; see warpkit_inventory.json. Do not run without SDC.')
    freeze_json(output / 'frozen_runs.json', [{k: v for k, v in r.items() if k != 'warpkit'} for r in runs])
    runtime = shutil.which('apptainer') or shutil.which('singularity') or 'apptainer'
    image = args.image or Path('/ZPOOL/data/tools/fmriprep-24.1.1.simg')
    if not image.is_file() and args.image is None:
        image = output / 'containers/fmriprep-24.1.1.sif'
    filters = output / 'filters'
    filters.mkdir(exist_ok=True)
    plan = []
    for subject in subjects:
        selected = [r for r in runs if r['subject'] == subject]
        (output / 'frozen_warpkit').mkdir(exist_ok=True)
        freeze_json(output / 'frozen_warpkit' / (subject + '.json'), {r['key']: r['warpkit'] for r in selected})
        save_json(filters / (subject + '.json'), filter_for(selected))
        work = Path('/ZPOOL/data/scratch') / os.environ.get('USER', 'tug87422') / 'rf1-betrayal-trust24' / output.name / subject
        cmd = command(runtime, image, upstream, output, work, subject, args.nprocs, args.mem_mb)
        plan.append(dict(subject=subject, argv=cmd, work=str(work), runs=[r['key'] for r in selected]))
    save_json(output / ('plan-pilot.json' if args.pilot else 'plan-full.json'),
              dict(historical_source=HISTORICAL_SOURCE, fmriprep_version=VERSION,
                   note='Historical settings with Trust/run filtering and resource caps; equivalence is not assumed.',
                   jobs=args.jobs, all_runs=runs, commands=plan))
    print('Output:', output, flush=True)
    print(f'Original combined confounds preserved: {len(runs)}/241; no new TEDANA fit needed for historical FEAT recovery.', flush=True)
    print(f'{len(subjects)} participants; {sum(len(x["runs"]) for x in plan)} Trust runs selected.', flush=True)
    for x in plan[:2]:
        print(shlex.join(x['argv']), flush=True)
    if not args.execute:
        print('DRY RUN; add --execute to launch. Missing image can be pulled with --fetch-image.')
        return
    if not args.pilot:
        require((output / 'pilot-verified.json').is_file(),
                'First run --pilot --fit --execute and pass historical comparisons before the full rebuild')
        require(json.loads((output / 'pilot-verified.json').read_text())['frozen_runs_sha256'] ==
                sha(output / 'frozen_runs.json'), 'Pilot used different inputs')
    if args.fit:
        require(shutil.which('feat'), 'FSL feat is unavailable; source the normal FSL setup')
    require(shutil.which(runtime), 'Apptainer/Singularity is unavailable')
    for path in (upstream / 'bids', Path('/ZPOOL/data/tools/templateflow'),
                 Path('/ZPOOL/data/tools/mplconfigdir'), Path('/ZPOOL/data/tools/licenses/fs_license.txt')):
        require(path.exists(), 'Missing runtime input: ' + str(path))
    if not image.is_file() and args.fetch_image:
        image.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run([runtime, 'pull', str(image), 'docker://nipreps/fmriprep:24.1.1'], check=True)
    require(image.is_file(), 'Missing 24.1.1 image; use --image PATH or --fetch-image')
    version = subprocess.check_output([runtime, 'run', '--cleanenv', str(image), '--version'], text=True, stderr=subprocess.STDOUT)
    record_container(output, image, version)
    env = dict(os.environ, APPTAINERENV_TEMPLATEFLOW_HOME='/opt/templateflow',
               SINGULARITYENV_TEMPLATEFLOW_HOME='/opt/templateflow',
               APPTAINERENV_MPLCONFIGDIR='/opt/mplconfigdir', SINGULARITYENV_MPLCONFIGDIR='/opt/mplconfigdir')
    subprocess.run([runtime, 'exec', '--cleanenv', '-B', f'{output}:/rebuild:ro',
                    str(image), 'python', '-c', FILTER_PREFLIGHT,
                    *[f'/rebuild/filters/{subject}.json' for subject in subjects]],
                   env=env, check=True)
    def launch(item):
        selected = [r for r in runs if r['subject'] == item['subject']]
        done = output / (item['subject'] + '.complete.json')
        if done.is_file():
            verify_outputs(selected)
            print('Already complete:', item['subject'], flush=True)
            if args.fit:
                for run in selected:
                    fit_run(run, output)
            return
        Path(item['work']).mkdir(parents=True, exist_ok=True)
        log = output / (item['subject'] + '.fmriprep.log')
        print('Running:', item['subject'], 'log:', log, flush=True)
        with log.open('a') as f:
            f.write('\n=== fMRIPrep attempt started ' + datetime.now(timezone.utc).isoformat() + ' ===\n')
            f.flush()
            subprocess.run(item['argv'], env=env, stdout=f, stderr=subprocess.STDOUT, check=True)
        verify_outputs(selected)
        save_json(done, dict(runs=item['runs'], stage='fmriprep_only',
                            historical_FEAT_equivalence='not_yet_tested'))
        print('Completed fMRIPrep:', item['subject'], flush=True)
        if args.fit:
            for run in selected:
                fit_run(run, output)
    lock = output / '.execution.lock'
    require(not lock.exists(), 'Rebuild execution lock exists; inspect the running job before restarting: ' + str(lock))
    lock.mkdir()
    save_json(lock / 'owner.json', dict(pid=os.getpid(), hostname=os.uname().nodename))
    try:
        with ThreadPoolExecutor(max_workers=args.jobs) as executor:
            futures = [executor.submit(launch, item) for item in plan]
            try:
                for future in as_completed(futures):
                    future.result()
            except BaseException:
                for future in futures:
                    future.cancel()
                raise
        if args.pilot and args.fit:
            save_json(output / 'pilot-verified.json', dict(subjects=subjects, maps='all historical COPEs/ZSTATs',
                                                         atol=1e-4, rtol=1e-4, frozen_runs_sha256=sha(output / 'frozen_runs.json')))
    finally:
        (lock / 'owner.json').unlink()
        lock.rmdir()
    print('Selected Trust preprocessing complete. Next: isolated FEAT reconstruction and historical comparisons.', flush=True)


if __name__ == '__main__':
    main()
