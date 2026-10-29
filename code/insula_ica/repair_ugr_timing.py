#!/usr/bin/env python3
"""Correct the authorized 10817 UGR TR, refit models and finish sm-0 recovery."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import traceback
import uuid

import nibabel as nib
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import run_timing as timing
import recover_ugr_unsmoothed as recovery
import run_unsmoothed_feat as sm0
from diagnose_trust_registration import sbref_selection
from run_bold_controls import publish, save


core=recovery.core
ROOT=core.ROOT
RECORDS=ROOT/'code/insula_ica/records'
KEY='sub-10817_ses-01_task-ugr_run-1'
EVIDENCE=RECORDS/'linux2-ugr-tr-diagnostic-20261029T113559Z-9253e5/evidence.json'
CORRECTION=dict(key=KEY,historical_tr_seconds=1.615,corrected_tr_seconds=1.710,
    authorization='user_2026-10-29',basis='All four current raw JSONs and NIfTI headers agree; user authorized corrected reanalysis.',
    historical_acquisition_provenance_resolved=False,event_times_rescaled=False)


def latest(pattern):
    choices=sorted(p for p in RECORDS.glob(pattern) if (p/'state.json').is_file())
    core.require(bool(choices),'No record found: '+pattern)
    return choices[-1]


def audit_timing(rows):
    results=[]
    for row in rows:
        sub='sub-'+row['subject'].removeprefix('sub-')
        try:
            actual=timing.acquisition_timing(recovery.stock.RAW/sub/'ses-01/func',row['key'])
            fsf=Path(row['active_feat'])/'design.fsf'
            if not fsf.is_file():fsf=Path(row['source_fsf'])
            modeled=timing.fsf_tr(fsf.read_text())
            status='matches' if np.isclose(actual['tr_seconds'],modeled,rtol=0,atol=1e-5) else 'model_TR_mismatch'
            results.append(dict(actual,modeled_tr_seconds=modeled,status=status))
        except Exception as exc:
            results.append(dict(key=row['key'],status='unresolved',error=str(exc)))
    return results


def event_audit(row):
    settings=core.fsl_settings(Path(row['active_feat'])/'design.fsf');events=[]
    for ev in range(1,int(settings['fmri(evs_orig)'])+1):
        if settings.get(f'fmri(shape{ev})')!='3':continue
        path=Path(settings[f'fmri(custom{ev})']);values=core.load_numeric(path)
        core.require(values.shape[1]==3 and np.isfinite(values).all() and (values[:,1]>=0).all(),'Invalid seconds-based task events')
        end=float(np.max(values[:,0]+values[:,1]))
        core.require(end<=row['nvol']*CORRECTION['corrected_tr_seconds']+1.,'Task EV extends beyond corrected scan duration')
        events.append(dict(ev=ev,path=str(path),sha256=core.sha(path),max_event_end_seconds=end))
    return dict(events=events,onsets_durations_retained_in_seconds=True,scanner_zero_offset_independently_verified=False,
        historical_scan_duration_seconds=row['nvol']*CORRECTION['historical_tr_seconds'],
        corrected_scan_duration_seconds=row['nvol']*CORRECTION['corrected_tr_seconds'])


def controlled_fit(row,bold,tr,root,label):
    folder=root/'sm5-controls';folder.mkdir(exist_ok=True)
    target=folder/(label+'.feat');fsf=folder/(label+'.fsf');marker=folder/(label+'.complete.json')
    original=Path(row['active_feat']);text=recovery.stock.rebuild.recovery_fsf((original/'design.fsf').read_text(),bold,row['confounds'],target)
    text=timing.stamp(text,tr)
    signature=dict(bold=core.identity(bold,False),fsf=text,confounds_sha256=row['confounds_sha256'])
    if marker.is_file():
        saved=json.loads(marker.read_text())
        core.require(saved['signature']==signature and all(core.identity(p['path'],p['sha256'] is not None)==p for p in saved['outputs']),
                     'Completed timing-control output/input changed')
        return saved['result']
    core.require(not target.exists(),'Partial timing-control FEAT exists; preserve and inspect: '+str(target))
    fsf.write_text(text);base=fsf.with_suffix('')
    with (folder/(label+'.model.log')).open('w') as stream:
        subprocess.run(['feat_model',str(base),row['confounds']],stdout=stream,stderr=subprocess.STDOUT,check=True)
    with (folder/(label+'.feat.log')).open('w') as stream:
        subprocess.run(['feat',str(fsf)],stdout=stream,stderr=subprocess.STDOUT,check=True)
    for name,expected in [('design.mat',base.with_suffix('.mat')),('design.con',base.with_suffix('.con'))]:
        a,b=core.fsl_matrix(expected),core.fsl_matrix(target/name)
        core.require(a.shape==b.shape and np.allclose(a,b,rtol=1e-8,atol=1e-8),'Timing-control model mismatch: '+name)
    a,b=core.fsl_matrix(original/'design.con'),core.fsl_matrix(target/'design.con')
    core.require(a.shape==b.shape and np.allclose(a,b,rtol=1e-8,atol=1e-8),'Timing repair changed contrast definitions')
    a,b=core.load_numeric(row['confounds']),core.load_numeric(target/'confoundevs.txt')
    core.require(a.shape==b.shape and np.array_equal(a,b),'Timing repair changed nuisance matrix')
    design=core.fsl_matrix(target/'design.mat');old=core.fsl_matrix(original/'design.mat')
    changed=not np.allclose(design,old,rtol=1e-8,atol=1e-8)
    core.require(changed==bool(np.isclose(tr,CORRECTION['corrected_tr_seconds'])),'Expected task-design timing change not established')
    data=nib.load(target/'stats/res4d.nii.gz');mask=nib.load(target/'mask.nii.gz').get_fdata()>0
    core.require(mask.any() and np.isclose(timing.image_tr(target/'stats/res4d.nii.gz')[0],tr,atol=1e-5), 'Residual TR/mask mismatch')
    values=np.asarray(data.dataobj,dtype='float32')[mask]
    core.require(np.isfinite(values).all() and (values.std(axis=1)>0).any(),'Invalid control residuals')
    core.require(float((target/'stats/dof').read_text())>=37,'Insufficient control degrees of freedom')
    paths=[target/name for name in ['filtered_func_data.nii.gz','stats/res4d.nii.gz','stats/threshac1.nii.gz',
        'stats/sigmasquareds.nii.gz','design.mat','design.con','confoundevs.txt','mask.nii.gz','mean_func.nii.gz']]
    paths.extend(target/f'stats/pe{k+1}.nii.gz' for k in range(design.shape[1]))
    contrasts=core.fsl_matrix(target/'design.con').shape[0]
    paths.extend(target/f'stats/{kind}{k}.nii.gz' for kind in ('cope','zstat') for k in range(1,contrasts+1))
    result=dict(state='complete',feat=str(target),tr_seconds=tr,smoothing_mm=5,task_design_changed=changed,
                contrast_definitions_preserved=True,exact_original_nuisance_matrix=True,
                role='corrected_analysis_candidate' if changed else 'wrong_TR_diagnostic_control_only')
    save(marker,dict(signature=signature,result=result,outputs=[core.identity(p,not p.name.endswith('.nii.gz')) for p in paths]))
    return result


def compare_controls(row,control,corrected):
    left,right=Path(control['feat']),Path(corrected['feat'])
    mask=(nib.load(left/'mask.nii.gz').get_fdata()>0)&(nib.load(right/'mask.nii.gz').get_fdata()>0)
    n=core.fsl_matrix(left/'design.con').shape[0];settings=core.fsl_settings(left/'design.fsf');comparisons=[]
    for index in range(1,n+1):
        for kind in ('cope','zstat'):
            a=nib.load(left/f'stats/{kind}{index}.nii.gz').get_fdata()[mask]
            b=nib.load(right/f'stats/{kind}{index}.nii.gz').get_fdata()[mask]
            valid=np.isfinite(a)&np.isfinite(b);a,b=a[valid],b[valid]
            comparisons.append(dict(map=f'{kind}{index}',contrast_name=settings.get(f'fmri(conname_real.{index})'),
                spatial_r=float(np.corrcoef(a,b)[0,1]) if a.std()>0 and b.std()>0 else None,
                relative_RMS_difference=float(np.linalg.norm(b-a)/np.linalg.norm(a)) if np.linalg.norm(a)>0 else None,
                maximum_absolute_difference=float(np.max(np.abs(b-a)))))
    return dict(key=row['key'],comparison='Same regenerated BOLD, smoothing, task-event files and nuisance values; TR/model sampling changes only.',
                historical_TR=control['tr_seconds'],corrected_TR=corrected['tr_seconds'],comparisons=comparisons,
                main_UGR_contrast_index=11,group_statistics_recomputed=False)


def run(record,root,state):
    source=sm0.SOURCE;rows=json.loads((source/'accepted-inputs.json').read_text())
    row=next(r for r in rows if r['key']==KEY);corrected_row=dict(row,tr=CORRECTION['corrected_tr_seconds'])
    evidence=json.loads(EVIDENCE.read_text())
    measured=timing.acquisition_timing(recovery.stock.RAW/'sub-10817/ses-01/func',KEY)
    core.require(np.isclose(measured['tr_seconds'],CORRECTION['corrected_tr_seconds'],rtol=0,atol=1e-5),'Approved acquisition TR no longer matches raw data')
    for item in evidence['echoes']:
        core.require(core.sha(item['json']['path'])==item['json']['sha256'],'Approved timing sidecar changed')
        core.require(core.identity(item['image']['path'],False)==item['image']['identity'],'Approved raw image changed')
    save(record/'correction-policy.json',CORRECTION);save(record/'event-timing-audit.json',event_audit(row))
    state['state']='auditing_cohort_and_Trust_timing';save(record/'state.json',state)
    audit=audit_timing(rows);save(record/'cohort-timing-audit.json',audit)
    trust=[x for x in audit if x['key'].startswith('sub-10817_') and '_task-trust_' in x['key']]
    save(record/'sub-10817-Trust-timing.json',trust)
    for item in trust:print('Trust timing:',item['key'],item.get('tr_seconds'),item.get('modeled_tr_seconds'),item['status'],flush=True)
    base=latest('linux2-ugr-anatomy-reuse-*');prior=json.loads((base/'state.json').read_text())
    combined_source=Path(prior['combined_sm0_record']);previous=json.loads((combined_source/'state.json').read_text())
    feat_root=Path(previous['output_root']);recovery.assert_no_active_outputs([root,feat_root])
    core.require(sum(v['state']=='complete' for v in previous['runs'].values())>=489,'Expected the completed 489-run sm-0 recovery')
    override=json.loads((base/'regenerated-ugr-bold-overrides.json').read_text())
    raw=recovery.raw_inputs(corrected_row,recovery.stock.RAW)
    actual_filter=recovery.stock.RAW.parent/'code/fmriprep_config.json'
    core.require(json.loads(actual_filter.read_text())==json.loads(recovery.stock.STOCK.read_text()),'Stock filter changed')
    container=json.loads((RECORDS/'linux2-trust-regeneration-20261027T162428Z-ef53c5/container.json').read_text())
    core.require(recovery.stock.rebuild.sha(recovery.stock.IMAGE)==container['sha256'],'24.1.1 image changed')
    save(record/'container.json',container);save(record/'raw-warpkit-audit.json',raw)
    anchors=json.loads((RECORDS/'linux2-trust-regeneration-20261027T162428Z-ef53c5/input-registry.json').read_text())
    anchor=dict(next(r for r in anchors if r['subject']=='sub-10817'),raw_images=raw['images'])
    case=root/'sub-10817/run-1';work=recovery.WORK/root.name/'sub-10817/run-1'
    runtime=shutil.which('apptainer') or shutil.which('singularity');core.require(runtime and shutil.which('feat_model'),'Container/FSL setup missing')
    if not case.exists():
        anatomy,identities=recovery.stock.stage_inputs(anchor,case,stock=actual_filter,task='ugr')
        work.mkdir(parents=True,exist_ok=False);argv=recovery.build_command(runtime,case,work,'sub-10817',anatomy)
        save(case/'execution-plan.json',dict(argv=argv,anatomical_inputs=identities))
    else:
        saved=json.loads((case/'execution-plan.json').read_text());argv,identities=saved['argv'],saved['anatomical_inputs']
    core.require(identities==[core.identity(i['path'],False) for i in identities],'Cached anatomy changed')
    env=dict(os.environ)
    for prefix in ('APPTAINERENV_','SINGULARITYENV_'):
        env.update({prefix+'TEMPLATEFLOW_HOME':'/opt/templateflow',prefix+'MPLCONFIGDIR':'/opt/mplconfigdir'})
    log=case/'sub-10817.fmriprep.log';state['state']='running_correct_TR_fmriprep';save(record/'state.json',state)
    if not (log.is_file() and 'fMRIPrep finished successfully!' in log.read_text(errors='replace')):
        pre=argv[:argv.index(str(recovery.stock.IMAGE))]+[str(recovery.stock.IMAGE),'python','-c',recovery.stock.ANATOMY_PREFLIGHT,'sub-10817','/rebuild/anatomy-preflight.json'];pre[1]='exec'
        with (case/'anatomy-preflight.log').open('w') as stream:subprocess.run(pre,env=env,stdout=stream,stderr=subprocess.STDOUT,check=True)
        recovery.stock.execute(argv,env,log)
    text=log.read_text(errors='replace');preflight=json.loads((case/'anatomy-preflight.json').read_text())
    core.require(preflight['anatomical_registration_nodes_absent'] and 'ANAT Stage 4: Found pre-computed registrations' in text,'Anatomy reuse not established')
    qc=recovery.bold_qc(corrected_row,case);save(record/'qc'/(KEY+'-bold.json'),qc)
    core.require(identities==[core.identity(i['path'],False) for i in identities],'Anatomy changed during processing')
    used=sbref_selection([dict(row,subject='sub-10817',new_bold=qc['path'],raw_images=raw['images'])])[0]
    core.require(used['status'] in ('same_run','bold_derived_reference_no_sbref'),'Wrong/unknown UGR reference selection')
    core.require(min(qc['coverage'].values())>=.8,'Corrected UGR insula coverage requires review')
    save(record/'preprocessing.json',dict(argv=argv,anatomy_preflight=preflight,anatomical_inputs=identities,reference_selection=used,
        actual_success_evidence=[line for line in text.splitlines() if 'Found pre-computed registrations' in line or 'finished successfully!' in line]))
    override[KEY]=dict(path=qc['path'],identity=qc['identity'],source_record=str(record),timing_correction=CORRECTION,
                       anatomical_reuse_verified=True,qc=str(record/'qc'/(KEY+'-bold.json')))
    override_path=record/'regenerated-ugr-bold-overrides.json';save(override_path,override)
    state['state']='fitting_corrected_sm0_and_matched_sm5_controls';save(record/'state.json',state)
    def sm0_fit():
        item=sm0.plan(row,feat_root,override[KEY]);result=sm0.fit(item,feat_root,record)
        return dict(result,feat=item['target'])
    with ThreadPoolExecutor(max_workers=3) as pool:
        zero=pool.submit(sm0_fit)
        old=pool.submit(controlled_fit,row,qc['path'],CORRECTION['historical_tr_seconds'],root,'same_BOLD_wrong_TR_1p615')
        new=pool.submit(controlled_fit,row,qc['path'],CORRECTION['corrected_tr_seconds'],root,'corrected_TR_1p710')
        sm0_result,old_result,new_result=zero.result(),old.result(),new.result()
    save(record/'matched-TR-statistics-comparison.json',compare_controls(row,old_result,new_result))
    save(record/'corrected-analysis-input.json',dict(row,tr=CORRECTION['corrected_tr_seconds'],historical_tr=CORRECTION['historical_tr_seconds'],
        active_feat=new_result['feat'],source_fsf=str(Path(new_result['feat'])/'design.fsf'),
        filtered_data=str(Path(new_result['feat'])/'filtered_func_data.nii.gz'),native_residual=str(Path(new_result['feat'])/'stats/res4d.nii.gz'),
        new_bold=qc['path'],timing_correction=CORRECTION,analysis_replacement_review_pending=True))
    state['fits']=dict(corrected_sm0=sm0_result,wrong_TR_sm5_control=old_result,corrected_sm5=new_result)
    combined=RECORDS/('linux2-sm-0-feat-'+record.name.removeprefix('linux2-ugr-TR-repair-')+'-TR-repair');combined.mkdir()
    for name in ('qc','errors'):(combined/name).mkdir()
    all_state=dict(state='preflight',runs={},output_root=str(feat_root),source_record=str(source),timing_repair_record=str(record),code_commit=state['code_commit'])
    from types import SimpleNamespace
    sm0.run(SimpleNamespace(source_record=source,jobs=64,execute=True,bold_overrides=override_path),feat_root,combined,all_state)
    try:publish(combined,'Complete corrected-TR UGR sm-0 fit and reuse the 489 completed FEATs')
    except Exception as exc:
        (combined/'publication-error.txt').write_text(str(exc));print('Git publication failed:',exc,flush=True)
    state.update(combined_sm0_record=str(combined),state='TR_repair_complete_QC_review_pending' if all_state['state']=='sm0_FEAT_complete_QC_pending' else 'TR_repair_combined_incomplete')
    core.require(core.identity(qc['path'],False)==qc['identity'],'Source BOLD changed during control fits')
    save(record/'state.json',state)
    mismatches=[x for x in audit if x['status']=='model_TR_mismatch'];unresolved=[x for x in audit if x['status']=='unresolved']
    (record/'report.md').write_text('# UGR timing repair and project timing audit\n\n'
        f'State: {state["state"]}. Corrected {KEY}: 1.615 → 1.710 seconds. '
        'Raw JSON/NIfTI files were not edited. Original task events in seconds and nuisance values were retained; '
        'the HRF/task design was regenerated and required to change. Existing historical fits remain available.\n\n'
        f'Cohort audit: {len(audit)} accepted runs; {len(mismatches)} model/acquisition TR mismatches; '
        f'{len(unresolved)} unresolved timing checks. See cohort-timing-audit.json.\n\n'
        'Trust for sub-10817:\n\n'+''.join(f'- {x["key"]}: acquisition {x.get("tr_seconds")}; model {x.get("modeled_tr_seconds")}; {x["status"]}.\n' for x in trust)+
        '\nSame-BOLD sm-5 fits isolate the effect of the TR change; matched-TR-statistics-comparison.json contains all COPE/ZSTAT contrasts. '
        'The 1.615-second fit is a diagnostic control and must not be promoted as corrected output. '
        'corrected-analysis-input.json identifies the corrected sm-5 candidate.\n\n'
        f'Combined sm-0 completion: {sum(v["state"]=="complete" for v in all_state["runs"].values())}/490. '
        'QC review is pending. Original group statistics, PPI, RSA and existing ICA fits have not been rerun or silently replaced. '
        'The corrected input registry is ready for explicit propagation into those dependent analyses. '
        'Mixed acquisition TRs must remain explicit when preparing later ICA inputs/frequency reports.\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--resume-record',type=Path);args=parser.parse_args()
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:6]
    record=RECORDS/('linux2-ugr-TR-repair-'+stamp);record.mkdir()
    root=Path(json.loads((args.resume_record/'state.json').read_text())['output_root']) if args.resume_record else ROOT/'derivatives/insula_ica'/('ugr-TR-repair-'+stamp)
    core.require(root.resolve().is_relative_to((ROOT/'derivatives/insula_ica').resolve()) and root.name.startswith('ugr-TR-repair-'),'Invalid repair output root')
    root.mkdir(parents=True,exist_ok=True);(record/'qc').mkdir()
    state=dict(state='preflight',output_root=str(root),key=KEY,code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip())
    save(record/'state.json',state);output=ROOT/'code/insula_ica/output';output.mkdir(exist_ok=True)
    with (output/'ugr-anatomy-reuse.lock').open('a') as lock,(output/'sm-0-feat.lock').open('a') as feat_lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);fcntl.flock(feat_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        try:run(record,root,state)
        except Exception as exc:
            state.update(state='failed',error=str(exc));(record/'failure.txt').write_text(traceback.format_exc());save(record/'state.json',state)
            print('ERROR:',exc,flush=True)
            logs=list((root/'sub-10817/run-1').glob('*.log'))+list((root/'sm5-controls').glob('*.log'))
            if logs:
                (record/'error-logs').mkdir()
                for path in logs:
                    with path.open('rb') as stream:
                        stream.seek(max(0,path.stat().st_size-65536));tail=stream.read().decode(errors='replace')
                    (record/'error-logs'/path.name).write_text(tail)
        try:publish(record,'Repair run-specific UGR timing, check Trust/cohort TRs and compare matched task fits')
        except Exception as exc:
            (record/'publication-error.txt').write_text(str(exc));print('Git publication failed:',exc,flush=True)
    print('Final state:',state['state'],'| record:',record,flush=True)
    raise SystemExit(0 if state['state']=='TR_repair_complete_QC_review_pending' else 1)


if __name__=='__main__':main()
