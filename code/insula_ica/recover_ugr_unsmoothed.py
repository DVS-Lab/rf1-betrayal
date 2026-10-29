#!/usr/bin/env python3
"""Pilot then expand UGR 24.1.1 with cached anatomy; resume only missing sm-0 FEATs."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import threading
import traceback
from types import SimpleNamespace
import uuid

import nibabel as nib
import numpy as np

import rerun_trust_references as stock
import run_unsmoothed_feat as sm0
from compare_bold_controls import frame_metrics
from diagnose_trust_registration import sbref_selection
from run_bold_controls import publish, save


core = stock.core
ROOT = core.ROOT
RECORDS = ROOT / 'code/insula_ica/records'
WORK = Path('/ZPOOL/data/scratch')/os.environ.get('USER','tug87422')/'rf1-ugr-anatomy-reuse'


def available_gib():
    mem = Path('/proc/meminfo')
    available = next(int(x.split()[1])/1024**2 for x in mem.read_text().splitlines()
                     if x.startswith('MemAvailable:')) if mem.is_file() else 32
    # Respect a container/job cgroup allocation when one exists.
    for base in (Path('/sys/fs/cgroup'),):
        try:
            maximum = (base/'memory.max').read_text().strip()
            if maximum != 'max':
                available = min(available, (int(maximum)-int((base/'memory.current').read_text()))/1024**3)
        except (OSError, ValueError):
            pass
    return available


class Budget:
    def __init__(self, cpus=None):
        self.cpus = cpus or (len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else os.cpu_count() or 1)
        self.cpu = 0
        self.mem = 0
        self.condition = threading.Condition()

    @contextmanager
    def reserve(self, cpu, mem):
        cpu = min(cpu, self.cpus)
        with self.condition:
            while self.cpu+cpu > self.cpus or available_gib()-self.mem < mem+12:
                self.condition.wait(timeout=5)
            self.cpu += cpu
            self.mem += mem
        try:
            yield
        finally:
            with self.condition:
                self.cpu -= cpu
                self.mem -= mem
                self.condition.notify_all()


def pilots(rows):
    counts = {s: sum(r['subject']==s for r in rows) for s in sorted({r['subject'] for r in rows})}
    single = next(s for s, n in counts.items() if n == 1)
    double = next(s for s, n in counts.items() if n == 2)
    return [single, double]


def assert_no_active_outputs(roots):
    """A killed submitter can leave children alive; never compete with them."""
    proc=Path('/proc')
    if not proc.is_dir():return
    roots=[p.resolve() for p in roots]
    for entry in proc.glob('[0-9]*'):
        try:
            if int(entry.name)==os.getpid() or entry.stat().st_uid!=os.getuid():continue
            cmd=(entry/'cmdline').read_bytes().decode(errors='replace').replace('\0',' ')
            cwd=(entry/'cwd').resolve(strict=True)
            core.require(not any(cwd.is_relative_to(root) or str(root)+'/' in cmd for root in roots),
                         'An existing process still uses the recovery outputs; do not duplicate it: PID '+entry.name)
        except (OSError, RuntimeError):
            continue


def raw_inputs(row, raw):
    subject = 'sub-'+row['subject'].removeprefix('sub-')
    folder = raw/subject/'ses-01/func'
    paths = []
    for echo in range(1, 5):
        image = folder/f'{row["key"]}_echo-{echo}_part-mag_bold.nii.gz'
        sidecar = image.with_name(image.name[:-7]+'.json')
        core.require(image.is_file() and sidecar.is_file(), 'Missing raw UGR echo/JSON: '+str(image))
        header = nib.load(image)
        core.require(len(header.shape)==4 and header.shape[3]==row['nvol'], 'Raw UGR time count mismatch')
        metadata = json.loads(sidecar.read_text())
        core.require(np.isclose(metadata['RepetitionTime'],row['tr']), 'Raw UGR TR mismatch')
        paths.append(image)
    fmap = stock.rebuild.warpkit_inputs(raw.parent,subject,int(row['run']),paths,task='ugr')
    core.require(fmap['ready'],'Missing associated UGR Warpkit BIDS products: '+str(fmap['problems']))
    return dict(images=[str(p) for p in paths], warpkit=fmap)


def build_command(runtime, case, work, subject, anatomy):
    argv = stock.command(runtime,case,work,subject,anatomy)
    argv[argv.index('--task-id')+1] = 'ugr'
    argv[argv.index('--nthreads')+1] = '8'
    argv[argv.index('--omp-nthreads')+1] = '4'
    argv[argv.index('--mem-mb')+1] = '12000'
    return argv


def bold_qc(row, case):
    old = Path(row['active_feat'])
    settings = core.fsl_settings(old/'design.fsf')
    name = Path(settings['feat_files(1)']).name
    path = case/'fmriprep-24'/('sub-'+row['subject'].removeprefix('sub-'))/'ses-01/func'/name
    image = nib.load(path)
    reference = nib.load(old/'mean_func.nii.gz')
    core.require(core.same_grid(image,reference) and image.shape==(*reference.shape,row['nvol']), 'Regenerated UGR grid/time count differs')
    core.require(np.isclose(image.header.get_zooms()[3],row['tr']), 'Regenerated UGR TR differs')
    maskpath = path.with_name(path.name.replace('_desc-preproc_bold.nii.gz','_desc-brain_mask.nii.gz'))
    mask = nib.load(maskpath)
    core.require(core.same_grid(mask,reference), 'Regenerated UGR mask grid differs')
    active = mask.get_fdata()>0
    core.require(active.any(),'Empty UGR brain mask')
    values = np.asarray(image.dataobj,dtype='float32')
    core.require(np.isfinite(values).all(), 'Nonfinite UGR BOLD')
    active &= values.std(axis=3)>1e-8
    del values
    atlas = Path(os.environ['FSLDIR'])/'data/atlases'
    masks,_,_ = core.anatomical_masks(atlas/'HarvardOxford/HarvardOxford-cort-maxprob-thr25-2mm.nii.gz',
                                     atlas/'HarvardOxford-Cortical.xml',reference)
    coverage = {r:float((m&active).sum()/m.sum()) for r,m in masks.items()}
    frame = frame_metrics(old/'example_func.nii.gz',path,old/'mask.nii.gz',row['nvol']//2)
    flags = [r+'_coverage_below_0.8' for r,v in coverage.items() if v<.8]
    if frame.get('spatial_r') is None or frame['spatial_r']<.95:
        flags.append('historical_reference_spatial_r_below_0.95_review_alignment')
    core.require(core.sha(row['confounds'])==row['confounds_sha256'],'Saved nuisance matrix changed')
    return dict(key=row['key'],path=str(path),identity=core.identity(path,False),
                coverage=coverage,historical_frame=frame,review_flags=flags,
                exact_saved_confounds_sha256=row['confounds_sha256'],finite_nonconstant_bold=True,
                historical_equality_required=False)


def run(args, record, state, root):
    previous = json.loads((args.sm0_record/'state.json').read_text())
    feat_root = Path(previous['output_root'])
    core.require(feat_root.resolve().is_relative_to((ROOT/'derivatives/insula_ica').resolve()),'Invalid sm-0 root')
    assert_no_active_outputs([root,feat_root])
    source = Path(previous['source_record'])
    rows = json.loads((source/'accepted-inputs.json').read_text())
    ugr = [r for r in rows if r['task']=='ugr']
    core.require(len(ugr)==249 and len({r['subject'] for r in ugr})==132,'Expected 249 UGR runs / 132 subjects')
    subjects = sorted({r['subject'] for r in ugr})
    pilot_subjects = pilots(ugr)
    # Anatomy source is the already generated same-version Trust cohort, not
    # smoothed FEAT data or a new-version anatomical substitution.
    trust = json.loads((RECORDS/'linux2-trust-regeneration-20261027T162428Z-ef53c5/input-registry.json').read_text())
    anchor = {r['subject'].removeprefix('sub-'):r for r in trust}
    raw = stock.RAW
    actual_filter = raw.parent/'code/fmriprep_config.json'
    core.require(json.loads(actual_filter.read_text())==json.loads(stock.STOCK.read_text()), 'Stock BIDS filter queries changed; no filter was edited')
    shutil.copyfile(actual_filter,record/'fmriprep_config.json')
    shutil.copyfile(raw.parent/'code/fmriprep.sh',record/'original-fmriprep.sh')
    runtime = shutil.which('apptainer') or shutil.which('singularity')
    core.require(runtime and shutil.which('feat') and os.environ.get('FSLDIR'),'Container/FSL setup missing')
    container = json.loads((RECORDS/'linux2-trust-regeneration-20261027T162428Z-ef53c5/container.json').read_text())
    core.require(stock.rebuild.sha(stock.IMAGE)==container['sha256'],'24.1.1 container differs from the verified Trust image')
    save(record/'container.json',container)
    save(record/'accepted-ugr-inputs.json',ugr)
    budget = Budget()
    state.update(state='auditing_UGR_inputs',pilot_subjects=pilot_subjects,
                 fmriprep_worker_slots=args.fmriprep_jobs,feat_worker_slots=args.feat_jobs,
                 usable_cpus=budget.cpus,initial_available_GiB=available_gib(),sm0_output_root=str(feat_root))
    save(record/'state.json',state)
    audit = {}
    for row in ugr:
        try:
            audit[row['key']] = dict(state='ready',**raw_inputs(row,raw))
        except Exception as exc:
            audit[row['key']] = dict(state='blocked',error=str(exc))
    save(record/'raw-warpkit-audit.json',audit)
    print(f'UGR raw/Warpkit ready: {sum(x["state"]=="ready" for x in audit.values())}/249',flush=True)
    core.require(all(audit[r['key']]['state']=='ready' for r in ugr if r['subject'] in pilot_subjects), 'UGR pilot raw/Warpkit inputs incomplete; see raw-warpkit-audit.json')
    env = dict(os.environ)
    for prefix in ('APPTAINERENV_','SINGULARITYENV_'):
        env.update({prefix+'TEMPLATEFLOW_HOME':'/opt/templateflow',prefix+'MPLCONFIGDIR':'/opt/mplconfigdir'})
    lock = threading.Lock()
    overrides = {}
    override_path = record/'regenerated-ugr-bold-overrides.json'
    save(override_path,overrides)
    feat_record = record/'ugr-feat';(feat_record/'qc').mkdir(parents=True,exist_ok=True)
    for name in ('fsf','logs','feat'):(feat_root/name).mkdir(parents=True,exist_ok=True)
    for name in ('cases','qc','errors'):(record/name).mkdir(exist_ok=True)
    save(record/'policy.json',dict(anatomy_reuse_required=True,stock_filter_unchanged=True,
        same_saved_task_and_TEDANA_nuisance_model=True,preprocessing_version='24.1.1',
        pilot_subjects=pilot_subjects,automatic_expansion_requires_all_pilot_BOLD_and_FEAT_QC=True,
        pilot_minimum_insula_coverage=.8,pilot_minimum_historical_frame_r=.90,
        no_new_TEDANA_fit=True,no_new_Warpkit_fit=True,existing_Trust_sm0_retained=True,
        fmriprep_threads=8,omp_threads=4,fmriprep_memory_MB=12000,
        cpu_and_available_memory_gate=True,headroom_GiB=12,
        outputs_root=str(root),sm0_output_root=str(feat_root),
        note='New UGR preprocessing is for the unsmoothed sensitivity. Original sm-5 and ICA inputs are not replaced.'))
    (record/'README.md').write_text('# Automated UGR anatomy reuse and unsmoothed FEAT recovery\n\n'
        'A single-run and a two-run participant form the pilot. Each accepted run has an isolated read-only BIDS view, '
        'unchanged stock query JSON, anatomy-only derivative cache, and separate work/output tree. '
        'The pilot must pass installed-container anatomy-graph checks, actual normalization-reuse log checks, '
        'BOLD grid/time/TR/finiteness, correct-run reference evidence, minimum 80% left/right/bilateral insula coverage, '
        'historical reference spatial correlation >=0.90, and the original joint FEAT design/contrast/nuisance identity checks '
        'at zero smoothing before automatic expansion. These pilot thresholds gate expansion, not participant exclusion. '
        'Full-cohort QC flags remain available for scientific review. Historical voxelwise equality is not required.\n\n'
        'Up to 32 subject workers run fMRIPrep (8 threads / 4 OpenMP / 12,000 MB each), and up to 64 workers run FEAT. '
        'Both share a CPU and available-memory gate with 12 GiB headroom; queue limits are not promises of simultaneous execution. '
        'No anatomy normalization, TEDANA fitting, Warpkit fitting, original sm-5 FEAT, or ICA is requested. '
        'Regenerated UGR BOLD overrides are explicit and identity checked. Existing completed Trust sm-0 fits are reused. '
        'The combined sm-0 record is published after recovery and validates all 490 accepted FEATs. '
        'See state.json, raw-warpkit-audit.json, cases/, qc/, ugr-feat/qc/, and errors/.\n')

    def checkpoint(title):
        save(record/'state.json',state)
        if args.publish:
            try:
                publish(record,title)
                (record/'publication-error.txt').unlink(missing_ok=True)
            except Exception as exc:
                (record/'publication-error.txt').write_text(str(exc))
                print('Git publication failed (processing continues):',exc,flush=True)

    def run_case(row):
        subject = row['subject']
        selected = [row]
        sub = 'sub-'+subject.removeprefix('sub-')
        # One-run views retain unambiguous actual SBRef selection evidence and
        # independent work/output trees; the JSON queries stay unchanged.
        case = root/sub/('run-'+str(row['run']))
        result = dict(state='failed',runs=[r['key'] for r in selected])
        try:
            core.require(all(audit[r['key']]['state']=='ready' for r in selected),'UGR raw/Warpkit audit blocked; see raw-warpkit-audit.json')
            anatomy_row = dict(anchor[subject.removeprefix('sub-')])
            anatomy_row['raw_images'] = [p for r in selected for p in audit[r['key']]['images']]
            work = WORK/root.name/sub/('run-'+str(row['run']))
            plan_path = case/'execution-plan.json'
            if not case.exists():
                anatomy, identities = stock.stage_inputs(anatomy_row,case,stock=actual_filter,task='ugr')
                work.mkdir(parents=True,exist_ok=False)
                argv = build_command(runtime,case,work,sub,anatomy)
                save(plan_path,dict(argv=argv,anatomical_inputs=identities,raw_inputs={r['key']:audit[r['key']] for r in selected}))
            else:
                saved = json.loads(plan_path.read_text())
                argv, identities = saved['argv'], saved['anatomical_inputs']
            core.require((case/'fmriprep_config.json').read_bytes()==actual_filter.read_bytes(),'Staged stock filter differs')
            core.require(identities==[core.identity(x['path'],False) for x in identities],'Cached anatomical inputs changed')
            log = case/(sub+'.fmriprep.log')
            done = log.is_file() and 'fMRIPrep finished successfully!' in log.read_text(errors='replace')
            if not done:
                with budget.reserve(8,12):
                    pre = argv[:argv.index(str(stock.IMAGE))]+[str(stock.IMAGE),'python','-c',stock.ANATOMY_PREFLIGHT,sub,'/rebuild/anatomy-preflight.json']
                    pre[1] = 'exec'
                    with (case/'anatomy-preflight.log').open('w') as stream:
                        subprocess.run(pre,env=env,stdout=stream,stderr=subprocess.STDOUT,check=True)
                    if log.exists():
                        shutil.copyfile(log,case/('prior-'+record.name+'.log'))
                    print('UGR fMRIPrep:',sub,'| cached anatomy preflight passed',flush=True)
                    stock.execute(argv,env,log)
            text = log.read_text(errors='replace')
            preflight = json.loads((case/'anatomy-preflight.json').read_text())
            core.require(preflight['anatomical_registration_nodes_absent'] and
                         'ANAT Stage 4: Found pre-computed registrations' in text and
                         'fMRIPrep finished successfully!' in text,'Missing actual anatomy reuse/success evidence')
            core.require(identities==[core.identity(x['path'],False) for x in identities],'Anatomy changed during processing')
            run_qc = {}
            for row in selected:
                qc = bold_qc(row,case)
                # sbref_selection reads this row's own case/run log, not the old Trust log.
                used = sbref_selection([dict(row,subject=sub,new_bold=qc['path'],raw_images=audit[row['key']]['images'])])[0]
                core.require(used['status'] in ('same_run','bold_derived_reference_no_sbref'),
                             'UGR run reference not established: '+used['status'])
                qc['reference_selection'] = used
                save(record/'qc'/(row['key']+'.json'),qc)
                run_qc[row['key']] = qc
            with lock:
                for key,qc in run_qc.items():
                    overrides[key] = dict(path=qc['path'],identity=qc['identity'],source_record=str(record),
                                          anatomical_reuse_verified=True,qc=str(record/'qc'/(key+'.json')))
                save(override_path,overrides)
            result = dict(state='complete',runs=list(run_qc),anatomical_reuse_verified=True,
                          reused_fmriprep=done,review_flags={k:q['review_flags'] for k,q in run_qc.items() if q['review_flags']})
            evidence=[line.strip() for line in text.splitlines() if
                      'ANAT Stage 4: Found pre-computed registrations' in line or 'fMRIPrep finished successfully!' in line]
            save(record/'cases'/(row['key']+'.json'),dict(result,argv=argv,anatomical_inputs=identities,
                                                       anatomy_preflight=preflight,actual_run_success_evidence=evidence))
        except Exception as exc:
            result['error'] = str(exc)
            (record/'errors'/(row['key']+'.txt')).write_text(traceback.format_exc())
            for name in ('anatomy-preflight.log',sub+'.fmriprep.log'):
                path = case/name
                if path.is_file():
                    with path.open('rb') as stream:
                        stream.seek(max(0,path.stat().st_size-32768))
                        tail = stream.read().decode(errors='replace')
                    (record/'errors'/(row['key']+'-'+name+'.tail.txt')).write_text(tail)
        with lock:
            state['runs'][row['key']] = result
            save(record/'state.json',state)
        print('UGR '+row['key']+': '+result['state']+(' | '+result['error'] if 'error' in result else ''),flush=True)
        return selected if result['state']=='complete' else []

    def subject_worker(subject):
        selected = [r for r in ugr if r['subject']==subject]
        ready = [row for r in selected for row in run_case(r)]
        with lock:
            state['subjects']['sub-'+subject.removeprefix('sub-')] = dict(
                state='complete' if len(ready)==len(selected) else 'failed',
                completed_runs=len(ready),expected_runs=len(selected))
            save(record/'state.json',state)
        return ready

    def feat_worker(row):
        try:
            with budget.reserve(1,3):
                item = sm0.plan(row,feat_root,overrides[row['key']])
                core.require(item['status']=='ready','UGR sm-0 task/model inputs missing')
                result = sm0.fit(item,feat_root,feat_record)
                qc = json.loads((feat_record/'qc'/(row['key']+'.json')).read_text())
                core.require(min(qc['coverage'].values())>=.8,'UGR sm-0 insula coverage below 0.8; visual review needed')
        except Exception as exc:
            result = dict(state='failed',error=str(exc))
            (record/'errors'/(row['key']+'-feat.txt')).write_text(traceback.format_exc())
        with lock:
            state['feat_runs'][row['key']] = result
            save(record/'state.json',state)
        print('UGR sm-0 '+row['key']+': '+result['state'],flush=True)

    state['state'] = 'running_UGR_pilot'
    save(record/'state.json',state)
    with ThreadPoolExecutor(max_workers=2) as pool:
        ready = [row for done in pool.map(subject_worker,pilot_subjects) for row in done]
    expected_pilot = [r for r in ugr if r['subject'] in pilot_subjects]
    core.require(len(ready)==len(expected_pilot),'UGR pilot fMRIPrep failed; expansion stopped')
    for row in expected_pilot:
        qc = json.loads((record/'qc'/(row['key']+'.json')).read_text())
        core.require(min(qc['coverage'].values())>=.8 and (qc['historical_frame'].get('spatial_r') or 0)>=.90,
                     'UGR pilot coverage/alignment QC needs review; expansion stopped')
    with ThreadPoolExecutor(max_workers=min(3,args.feat_jobs)) as pool:
        list(pool.map(feat_worker,ready))
    core.require(all(state['feat_runs'][r['key']]['state']=='complete' for r in ready),'UGR pilot sm-0 model/QC failed; expansion stopped')
    state.update(state='UGR_pilot_passed_expanding',pilot_passed=True)
    checkpoint('UGR anatomy-reuse pilot passed; automatically expand missing unsmoothed inputs')
    # Functional processing and new FEAT jobs overlap; both share one CPU/RAM gate.
    with ThreadPoolExecutor(max_workers=args.feat_jobs) as feat_pool:
        feat_futures = []
        with ThreadPoolExecutor(max_workers=args.fmriprep_jobs) as pool:
            futures = [pool.submit(subject_worker,s) for s in subjects if s not in pilot_subjects]
            for future in as_completed(futures):
                feat_futures.extend(feat_pool.submit(feat_worker,row) for row in future.result())
        state['state'] = 'finishing_missing_UGR_FEATs'
        save(record/'state.json',state)
        for future in as_completed(feat_futures):
            future.result()
    checkpoint('Record UGR anatomy reuse and unsmoothed FEAT recovery results')
    # Reuse all old Trust successes and all new UGR completion markers, and publish
    # one normal sm-0 record covering the entire accepted 490-run registry.
    combined = RECORDS/('linux2-sm-0-feat-'+record.name.removeprefix('linux2-ugr-anatomy-reuse-')+'-ugr-resume')
    combined.mkdir(exist_ok=False)
    for name in ('qc','errors'):(combined/name).mkdir()
    all_state = dict(state='preflight',runs={},output_root=str(feat_root),source_record=str(source),
                     code_commit=state['code_commit'],ugr_recovery_record=str(record))
    sm0.run(SimpleNamespace(source_record=source,jobs=args.feat_jobs,execute=True,bold_overrides=override_path),
            feat_root,combined,all_state)
    state['combined_sm0_record'] = str(combined)
    if args.publish:
        try:
            publish(combined,'Resume missing UGR sm-0 FEATs and retain all completed Trust fits')
        except Exception as exc:
            (combined/'publication-error.txt').write_text(str(exc))
            print('Git publication failed:',exc,flush=True)
    state.update(state='UGR_and_sm0_complete_QC_review_pending' if all_state['state']=='sm0_FEAT_complete_QC_pending'
                 and len(overrides)==249 and all(v['state']=='complete' for v in state['feat_runs'].values())
                 else 'UGR_recovery_incomplete_review_errors',finished_utc=datetime.now(timezone.utc).isoformat())
    checkpoint('Finish automated UGR pilot, anatomy reuse and missing sm-0 FEAT recovery')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sm0-record',type=Path)
    parser.add_argument('--resume-record',type=Path)
    parser.add_argument('--fmriprep-jobs',type=int,default=32)
    parser.add_argument('--feat-jobs',type=int,default=64)
    parser.add_argument('--publish',action='store_true')
    args = parser.parse_args()
    core.require(1<=args.fmriprep_jobs<=32 and 1<=args.feat_jobs<=64,'Worker limits: 32 fMRIPrep / 64 FEAT')
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:6]
    if args.resume_record is None and args.sm0_record is None:
        existing=sorted(p for p in RECORDS.glob('linux2-ugr-anatomy-reuse-*') if (p/'state.json').is_file())
        if existing:args.resume_record=existing[-1]
    if args.resume_record:
        prior = json.loads((args.resume_record/'state.json').read_text())
        root = Path(prior['output_root'])
        args.sm0_record = Path(prior['sm0_source_record'])
        print('Reusing prior UGR output root:',root,flush=True)
    else:
        root = ROOT/'derivatives/insula_ica'/('ugr-anatomy-reuse-'+stamp)
    core.require(root.resolve().is_relative_to((ROOT/'derivatives/insula_ica').resolve()) and root.name.startswith('ugr-anatomy-reuse-'),'Invalid UGR root')
    if args.sm0_record is None:
        choices = sorted(p for p in RECORDS.glob('linux2-sm-0-feat-*') if (p/'state.json').is_file())
        core.require(bool(choices),'No existing sm-0 FEAT record found')
        args.sm0_record = choices[-1]
    record = RECORDS/('linux2-ugr-anatomy-reuse-'+stamp)
    record.mkdir()
    root.mkdir(parents=True,exist_ok=True)
    state = dict(state='preflight',subjects={},runs={},feat_runs={},output_root=str(root),sm0_source_record=str(args.sm0_record),
                 started_utc=datetime.now(timezone.utc).isoformat(),
                 code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip())
    save(record/'state.json',state)
    output = ROOT/'code/insula_ica/output';output.mkdir(exist_ok=True)
    with (output/'ugr-anatomy-reuse.lock').open('a') as lock, (output/'sm-0-feat.lock').open('a') as feat_lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        fcntl.flock(feat_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        try:
            run(args,record,state,root)
        except Exception as exc:
            state.update(state='failed',error=str(exc))
            (record/'failure.txt').write_text(traceback.format_exc())
            save(record/'state.json',state)
            print('ERROR:',exc,flush=True)
            if args.publish:
                try:publish(record,'Record UGR anatomy-reuse pilot or recovery failure without overwriting prior fits')
                except Exception as publication_error:
                    (record/'publication-error.txt').write_text(str(publication_error))
                    print('Git publication failed:',publication_error,flush=True)
    print('Final state:',state['state'],'| record:',record,flush=True)
    raise SystemExit(0 if state['state']=='UGR_and_sm0_complete_QC_review_pending' else 1)


if __name__ == '__main__':
    main()
