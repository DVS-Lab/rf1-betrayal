#!/usr/bin/env python3
"""Hand off the stock-filter queue, retaining active FEAT and completed outputs."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import argparse
import fcntl
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import threading
import time
import traceback
import uuid

import nibabel as nib
import numpy as np

import rerun_trust_references as original
from diagnose_trust_registration import sbref_selection
from regenerate_trust import exact_confounds, qc_run
from export_trust_qc import montage
from run_bold_controls import save, publish


ROOT=original.core.ROOT


def validate_worker_limits(jobs,fmriprep_jobs):
    original.core.require(1<=fmriprep_jobs<=jobs<=64,'Require 1 <= fmriprep jobs <= total jobs <= 64')


def processes():
    """Linux process identity includes start ticks, so reused PIDs are distinct."""
    result={}
    for entry in Path('/proc').iterdir():
        if not entry.name.isdigit():
            continue
        try:
            if entry.stat().st_uid!=os.getuid():
                continue
            fields=(entry/'stat').read_text().rsplit(')',1)[1].split()
            args=[s.decode(errors='replace') for s in (entry/'cmdline').read_bytes().split(b'\0') if s]
            result[int(entry.name)]=dict(pid=int(entry.name),ppid=int(fields[1]),
                start=fields[19],state=fields[0],args=args,cwd=str((entry/'cwd').resolve()))
        except (OSError,ValueError,IndexError):
            continue
    return result


def alive(identity,table):
    current=table.get(identity['pid'])
    return bool(current and current['start']==identity['start'] and current['state'] not in ('Z','X'))


def runner(table):
    candidates=[p for p in table.values() if p['cwd']==str(ROOT) and p['state'] not in ('Z','X') and
                any(Path(a).name=='rerun_trust_references.py' for a in p['args']) and
                Path(p['args'][0]).name.startswith('python')]
    original.core.require(len(candidates)<=1,'More than one original runner; no signals sent')
    return candidates[0] if candidates else None


def feat_children(table,parent,root):
    children=[p for p in table.values() if p['ppid']==parent['pid'] and p['state'] not in ('Z','X')]
    for child in children:
        fsfs=[Path(a) for a in child['args'] if a.endswith('.fsf')]
        if not any(Path(a).name=='feat' for a in child['args']) or len(fsfs)!=1 or not fsfs[0].is_relative_to(root):
            return None
    return children


def stopped_threads(pid):
    states=[p.read_text().rsplit(')',1)[1].split()[0] for p in Path(f'/proc/{pid}/task').glob('*/stat')]
    return bool(states) and all(s in ('T','t','Z','X') for s in states)


def handoff(root,record):
    """Pause only the original Python scheduler at a FEAT-only boundary.

    FEAT logs go directly to files, so those children do not need a Python pipe
    reader. Never stop a scheduler with an active container/preflight child.
    Only the scheduler receives signals; its FEAT children continue unchanged.
    """
    announced=0
    while True:
        table=processes();parent=runner(table)
        if parent is None:
            save(record/'handoff.json',dict(original_runner_absent=True,inherited_feat=[]))
            return []
        stopped=False
        try:
            original.core.require(alive(parent,processes()),'Original process identity changed')
            os.kill(parent['pid'],signal.SIGSTOP);stopped=True
            for _ in range(100):
                table=processes()
                if not alive(parent,table):
                    stopped=False
                    break
                if table[parent['pid']]['state'] in ('T','t') and stopped_threads(parent['pid']):
                    break
                time.sleep(.02)
            if not stopped:
                continue
            original.core.require(table[parent['pid']]['state'] in ('T','t') and stopped_threads(parent['pid']),
                                  'Scheduler threads did not all stop')
            children=feat_children(table,parent,root)
            if children is None:
                os.kill(parent['pid'],signal.SIGCONT);stopped=False
                if time.monotonic()-announced>60:
                    print('Waiting for current fMRIPrep/preflight jobs to reach FEAT; no compute job stopped.',flush=True)
                    announced=time.monotonic()
                time.sleep(5)
                continue
            details=dict(original_runner=parent,inherited_feat=children,
                policy='Original Python scheduler only is terminated; all active FEAT processes continue.',
                time_utc=datetime.now(timezone.utc).isoformat())
            save(record/'handoff.json',details)
            # SIGTERM is pending before continuation. Default Python SIGTERM
            # exits the scheduler without dispatching additional queued work.
            os.kill(parent['pid'],signal.SIGTERM)
            os.kill(parent['pid'],signal.SIGCONT);stopped=False
            for _ in range(100):
                if not alive(parent,processes()):
                    print(f'Handoff complete; {len(children)} active FEAT jobs retained.',flush=True)
                    return children
                time.sleep(.1)
            raise ValueError('Original scheduler did not exit; replacement processing has not started')
        finally:
            if stopped and alive(parent,processes()):
                os.kill(parent['pid'],signal.SIGCONT)


def case_processes(case,inherited):
    table=processes()
    live=[p for p in inherited if alive(p,table) and any(str(case)+'/' in a for a in p['args'])]
    for p in table.values():
        if p['pid']==os.getpid() or p['state'] in ('Z','X'):
            continue
        if Path(p['cwd']).is_relative_to(case) or any(a.startswith(str(case)+'/') for a in p['args']):
            live.append(p)
    return live


class MemoryBudget:
    """Bound fMRIPrep count and reserve RAM for new stages; existing jobs drain."""
    def __init__(self,fmriprep_jobs):
        self.limit=fmriprep_jobs;self.fmriprep=0;self.reserved=0;self.lock=threading.Condition();self.last_notice={}

    @staticmethod
    def available():
        lines=Path('/proc/meminfo').read_text().splitlines()
        return next(int(line.split()[1])/1024**2 for line in lines if line.startswith('MemAvailable:'))

    def acquire(self,kind):
        amount=28 if kind=='fmriprep' else 4
        with self.lock:
            while self.available()-self.reserved<amount+12 or (kind=='fmriprep' and self.fmriprep>=self.limit):
                if time.monotonic()-self.last_notice.get(kind,0)>60:
                    print(f'Waiting for {kind} capacity: available {self.available():.1f} GiB, '
                          f'reservations {self.reserved} GiB, fMRIPrep {self.fmriprep}/{self.limit}.',flush=True)
                    self.last_notice[kind]=time.monotonic()
                self.lock.wait(timeout=10)
            self.reserved+=amount
            self.fmriprep+=kind=='fmriprep'

    def release(self,kind):
        with self.lock:
            self.reserved-=28 if kind=='fmriprep' else 4
            self.fmriprep-=kind=='fmriprep'
            self.lock.notify_all()


def reuse_qc(row,case,source,record):
    """Reuse completed numerical QC; recheck small model/confound/header assets."""
    old=source/'qc'/(row['key']+'.json');png=old.with_suffix('.png')
    target=case/'feat_reconstruction'/(row['key']+'.feat')
    report=target/'historical_comparison.json'
    if not old.is_file() or not png.is_file() or not report.is_file():
        return None
    qc=json.loads(old.read_text());comparison=json.loads(report.read_text())
    assets=[Path(row['new_bold']),target/'stats/res4d.nii.gz',target/'mask.nii.gz',
            target/'design.mat',target/'design.con',target/'confoundevs.txt']
    if not all(p.is_file() and p.stat().st_mtime_ns<=old.stat().st_mtime_ns for p in assets):
        return None
    original.core.require(qc.get('finite_nonconstant_residuals') is True,'Prior QC did not validate residuals')
    original.core.require(Path(qc['residual']).resolve()==assets[1].resolve(),'Prior QC refers to a different fit')
    original.core.require(not comparison['missing_original'] and not comparison['missing_reconstructed'] and
        len(comparison['comparisons'])==2*comparison['ncontrasts'] and
        all(d['grid_matches'] for d in comparison['comparisons']),
        'Prior FEAT comparison is incomplete')
    original.rebuild.verify_outputs([row])
    residual=nib.load(assets[1]);bold=nib.load(assets[0])
    original.core.require(original.core.same_grid(bold,residual) and bold.shape==residual.shape,'Residual header changed')
    for name in ('design.mat','design.con'):
        a=original.core.fsl_matrix(Path(row['original_feat'])/name);b=original.core.fsl_matrix(target/name)
        original.core.require(a.shape==b.shape and np.allclose(a,b,atol=1e-8,rtol=1e-8),'Full model changed')
    exact_confounds(row,target)
    shutil.copyfile(old,record/'qc'/old.name);shutil.copyfile(png,record/'qc'/png.name)
    return qc


def run(args,record,state):
    source=args.record.resolve();old=json.loads((source/'state.json').read_text())
    root=Path(old['output_root']).resolve()
    original.core.require(root.is_relative_to((ROOT/'derivatives/insula_ica').resolve()) and
        root.name.startswith('trust-stock-filter-'),'Unexpected output root')
    # This lock is shared by all replacement submitters for this exact root.
    with (root/'.resume.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        runtime=shutil.which('apptainer') or shutil.which('singularity')
        original.core.require(runtime and shutil.which('feat') and os.environ.get('FSLDIR'),'Normal container/FSL environment required')
        source_registry=ROOT/'code/insula_ica/records'/original.SOURCE
        expected=json.loads((source_registry/'container.json').read_text())
        original.core.require(original.rebuild.sha(original.IMAGE)==expected['sha256'],'Container digest changed')
        stock=source/'fmriprep_config.json'
        original.core.require(json.loads(stock.read_text())==json.loads(original.STOCK.read_text()),'Original stock filter changed')
        rows=json.loads((source_registry/'input-registry.json').read_text())
        audit=json.loads((source/'reference-selection-audit.json').read_text())
        repair=ROOT/'code/insula_ica/records/linux2-registration-fix-sub-10834-20261028T151154Z-acdf91/review-decision.json'
        repaired={r['key'] for r in json.loads(repair.read_text())['selected_inputs']}
        selected=original.choose_runs(rows,audit,repaired)
        save(record/'reference-selection-audit.json',audit)
        shutil.copyfile(stock,record/stock.name)
        shutil.copyfile(source/'original-fmriprep.sh',record/'original-fmriprep.sh')
        state.update(state='handoff_pending',output_root=str(root),selected_runs=len(selected),
                     source_record=str(source),jobs=args.jobs,fmriprep_jobs=args.fmriprep_jobs)
        save(record/'state.json',state)
        inherited=handoff(root,record)
        # Capture the old record only after its scheduler has exited.
        save(record/'previous-state.json',json.loads((source/'state.json').read_text()))
        state.update(state='running',inherited_processes=inherited);save(record/'state.json',state)
        for name in ('cases','qc','errors'):(record/name).mkdir()
        env=dict(os.environ)
        for prefix in ('APPTAINERENV_','SINGULARITYENV_'):
            env.update({prefix+'TEMPLATEFLOW_HOME':'/opt/templateflow',prefix+'MPLCONFIGDIR':'/opt/mplconfigdir'})
        budget=MemoryBudget(args.fmriprep_jobs);state_lock=threading.Lock()

        def launch(original_row):
            subject=original_row['subject'];key=original_row['key'];case=root/subject
            result=dict(subject=subject,key=key,state='failed')
            try:
                while case_processes(case,inherited):
                    time.sleep(5)
                original.core.require(original.core.sha(original_row['confounds'])==original_row['confounds_sha256'],'Original confounds changed')
                if not case.exists():
                    source_anat,identities=original.stage_inputs(original_row,case,stock=stock)
                    work=original.WORK/root.name/subject;work.mkdir(parents=True,exist_ok=False)
                    argv=original.command(runtime,case,work,subject,source_anat)
                    plan=dict(argv=argv,anatomical_inputs=identities,original_row=original_row)
                else:
                    plan_path=source/'cases'/(subject+'.json')
                    original.core.require(plan_path.is_file(),'Existing case has no recorded plan')
                    plan=json.loads(plan_path.read_text());argv=plan['argv'];identities=plan['anatomical_inputs']
                save(record/'cases'/(subject+'.json'),dict(plan,state='resuming'))
                original.core.require((case/'fmriprep_config.json').read_bytes()==stock.read_bytes(),'Case filter changed')
                original.core.require(identities==[original.core.identity(x['path'],content=False) for x in identities],'Cached anatomy changed')
                log=case/(subject+'.fmriprep.log')
                done=log.is_file() and 'fMRIPrep finished successfully!' in log.read_text(errors='replace')
                if not done:
                    budget.acquire('fmriprep')
                    try:
                        raw_paths=[Path(p) for p in original_row['raw_images']]
                        original.core.require(original.rebuild.warpkit_inputs(original.RAW.parent,subject,int(original_row['run']),raw_paths)['ready'],'Existing Warpkit inputs missing')
                        preflight=argv[:argv.index(str(original.IMAGE))]+[str(original.IMAGE),'python','-c',original.ANATOMY_PREFLIGHT,subject,'/rebuild/anatomy-preflight.json']
                        preflight[1]='exec'
                        with (case/'anatomy-preflight.log').open('w') as stream:
                            subprocess.run(preflight,env=env,stdout=stream,stderr=subprocess.STDOUT,check=True)
                        if log.exists():
                            shutil.copyfile(log,case/('prior-fmriprep-'+record.name+'.log'))
                        print('fMRIPrep unfinished/new only:',subject,flush=True)
                        original.execute(argv,env,log)
                    finally:
                        budget.release('fmriprep')
                text=log.read_text(errors='replace')
                original.core.require('fMRIPrep finished successfully!' in text and
                    'ANAT Stage 4: Found pre-computed registrations' in text,'Missing preprocessing/anatomical reuse evidence')
                preflight=json.loads((case/'anatomy-preflight.json').read_text())
                original.core.require(preflight['anatomical_registration_nodes_absent'],'Anatomical reuse not established')
                row=dict(original_row,new_bold=str(case/'fmriprep-24'/subject/'ses-01/func'/Path(original_row['new_bold']).name))
                used=sbref_selection([row])[0]
                original.core.require(used['status'] in ('same_run','bold_derived_reference_no_sbref'),
                    'Actual reference selection remains unverified: '+used['status'])
                qc=reuse_qc(row,case,source,record)
                reused=qc is not None
                if qc is None:
                    budget.acquire('feat_qc')
                    try:
                        original.rebuild.verify_outputs([row])
                        original.rebuild.fit_run(row,case,strict_historical=False)
                        qc=qc_run(row,case);save(record/'qc'/(key+'.json'),qc)
                        montage(qc,record/'qc'/(key+'.png'),Path(os.environ['FSLDIR']))
                    finally:
                        budget.release('feat_qc')
                result.update(state='complete_qc_review_pending',row=row,reference_selection=used,
                    anatomical_reuse_verified=True,fmriprep_reused=done,numerical_qc_reused=reused,
                    review_flags=qc['review_flags'])
                save(record/'cases'/(subject+'.json'),dict(argv=argv,anatomical_inputs=identities,
                    anatomy_preflight=preflight,qc_source_record=str(source) if reused else str(record),**result))
            except Exception as exc:
                result['error']=str(exc);(record/'errors'/(subject+'.txt')).write_text(traceback.format_exc())
            with state_lock:
                state['subjects'][subject]=result;save(record/'state.json',state)
            print(subject+': '+result['state']+((' | '+result['error']) if 'error' in result else
                (' | existing QC retained' if result['numerical_qc_reused'] else '')),flush=True)

        # Count inherited FEAT jobs within the new worker limit immediately.
        selected.sort(key=lambda r:(not any(any(str(root/r['subject'])+'/' in a for a in p['args']) for p in inherited),r['subject']))
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            for future in as_completed([pool.submit(launch,r) for r in selected]):future.result()
        failed=sum(s['state']=='failed' for s in state['subjects'].values())
        state.update(state='complete_with_errors' if failed else 'functional_reruns_complete_qc_pending',
            completed=len(selected)-failed,failed=failed,finished_utc=datetime.now(timezone.utc).isoformat())
        save(record/'state.json',state)
        save(record/'input-registry.json',[s['row'] for s in state['subjects'].values() if s['state']=='complete_qc_review_pending'])
        (record/'README.md').write_text('# Trust queue resumed with existing outputs\n\n'+
            f'Workers: {args.jobs}; new fMRIPrep limit: {args.fmriprep_jobs}. Completed: {len(selected)-failed}/{len(selected)}; failed: {failed}.\n\n'+
            'Existing FEAT children were allowed to finish. Successful fMRIPrep runs were not repeated. '
            'Completed FEAT/QC outputs were retained with model, confound and header checks. '
            'Reference selection was verified from actual log evidence. Unknown or mismatched selection remains a failure. '
            'New stages reserve 28 GiB (fMRIPrep) or 4 GiB (FEAT/QC), plus 12 GiB available-memory headroom. '
            'No ICA or exclusions were performed.\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record',required=True,type=Path)
    parser.add_argument('--jobs',type=int,default=64)
    parser.add_argument('--fmriprep-jobs',type=int,default=2)
    parser.add_argument('--publish',action='store_true')
    args=parser.parse_args()
    original.core.require(Path('/proc/self/stat').is_file(),'Linux process inspection is required')
    validate_worker_limits(args.jobs,args.fmriprep_jobs)
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:6]+'-resume'
    record=ROOT/'code/insula_ica/records'/('linux2-trust-stock-filter-'+stamp);record.mkdir()
    state=dict(state='preflight',subjects={},started_utc=datetime.now(timezone.utc).isoformat(),
        code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip())
    try:run(args,record,state)
    except Exception as exc:
        state.update(state='failed',error=str(exc));(record/'failure.txt').write_text(traceback.format_exc());print('ERROR:',exc,flush=True)
    save(record/'state.json',state)
    if args.publish:
        try:publish(record,'Record resumed Trust stock-filter queue with retained completed outputs')
        except Exception as exc:print('Git publication failed:',exc,flush=True)
    print('Final state:',state['state'],'| record:',record,flush=True)
    raise SystemExit(0 if state['state']=='functional_reruns_complete_qc_pending' else 1)


if __name__=='__main__':main()
