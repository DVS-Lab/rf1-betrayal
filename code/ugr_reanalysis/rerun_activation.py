#!/usr/bin/env python3
"""Reuse valid L1/L2 estimates for separate, original-model UGR group repairs."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import threading
import traceback
import uuid

import nibabel as nib
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'code/insula_ica'))
import insula_ica as core
import status_insula_streams as status
from run_bold_controls import save, publish
from audit_inputs import contrast_diagnostics

RECORDS = ROOT/'code/ugr_reanalysis/records'
POLICY = 'original-UGR-activation-TR1710-valid-run2-for-nonestimable-contrasts-v1'
BAD = {5, 6, 9, 10, 11, 12, 14, 15, 17}
TR_RECORD = ROOT/'code/insula_ica/records/linux2-ugr-TR-repair-20261029T125025Z-5ad843'


def freeze(path, value):
    if path.exists(): core.require(json.loads(path.read_text()) == value, 'Frozen repair inputs changed: '+str(path))
    else: save(path, value)


def original_contrasts(audit):
    rows = core.read_table(audit/'original_UGR_contrast_audit.tsv')
    flagged = [r for r in rows if float(r['cutoff']) == 1e-6 and r['requires_review'] == 'True']
    core.require({r['key'] for r in flagged} == {'sub-11031_ses-01_task-ugr_run-1'} and
                 {int(r['contrast']) for r in flagged} == BAD, 'Original contrast flags differ; review before rerun')
    return BAD


def replace_settings(text, updates):
    """Change only explicitly named FSF settings; preserve all scientific EVs."""
    seen = set(); output = []
    for line in text.splitlines():
        match = re.match(r'^\s*set\s+(\S+)\s+.*$', line)
        if match and match[1] in updates:
            key = match[1]; value = str(updates[key])
            core.require('\n' not in value, 'Multiline FSF setting')
            value = value.replace('\\', '\\\\').replace('"', '\\"')
            line = f'set {key} "{value}"'; seen.add(key)
        output.append(line)
    core.require(seen == set(updates), 'Missing saved FSF setting')
    return '\n'.join(output)+'\n'


def choose_input(original, cope, corrected_feat):
    if '/sub-10817/' in original:
        expected = f'/L1_task-ugr_model-3_type-act_run-1_sm-5.feat/stats/cope{cope}.nii.gz'
        core.require(original.endswith(expected), 'Unexpected sub-10817 group input')
        return corrected_feat/'stats'/f'cope{cope}.nii.gz', 'reuse_completed_TR_1p710_sm5_fit'
    if '/sub-11031/' in original and cope in BAD:
        expected = f'/L2_task-ugr_model-3_type-act_sm-5.gfeat/cope{cope}.feat/stats/cope1.nii.gz'
        core.require(original.endswith(expected), 'Unexpected sub-11031 group input')
        feat = ROOT/'derivatives/fsl/sub-11031/ses-01/L1_task-ugr_model-3_type-act_run-2_sm-5.feat'
        return feat/'stats'/f'cope{cope}.nii.gz', 'use_original_valid_run2_for_this_contrast'
    return Path(original), 'reuse_original_participant_estimate'


def source_file(group, name):
    return next((f/name for f in (group/'cope1.feat', group) if (f/name).is_file()), group/'cope1.feat'/name)


def stage_input(target, folder, reason, reference):
    """Copy small, already-standard-space inputs; FEAT must not touch sources.

    design.lev identifies an already-standard higher-level packet: featregapply
    exits before writing, and find_std reads its stats/mask/bg_image directly.
    No lower-level refit, registration or interpolation occurs here.
    """
    feat=target.parent.parent; name=target.name; number=int(re.fullmatch(r'cope(\d+)\.nii\.gz',name)[1])
    standard=feat/'reg_standard'; higher=(feat/'design.lev').exists()
    if not higher and (standard/'stats'/name).is_file(): location=standard
    else:
        location=feat
        if not higher and reason!='reuse_completed_TR_1p710_sm5_fit':
            transform=feat/'reg/example_func2standard.mat'
            core.require(transform.is_file() and np.allclose(np.loadtxt(transform),np.eye(4),atol=1e-8),
                         'Missing standardized data or verified identity registration: '+str(feat))
    effect=location/'stats'/name; variance=effect.with_name('var'+name)
    mask=location/'mask.nii.gz'
    background=next((location/n for n in ('bg_image.nii.gz','example_func.nii.gz','mean_func.nii.gz') if (location/n).is_file()),None)
    core.require(background is not None,'No standard-space source background: '+str(location))
    df=location/'stats'/f'tdof_t{number}.nii.gz'
    if not df.is_file():df=feat/'stats'/f'tdof_t{number}.nii.gz'
    scalar_df=not df.is_file()
    if scalar_df:df=feat/'stats/dof'
    paths=[effect,variance,mask,background,df]
    paths.extend(feat/n for n in ('design.con','design.lcon') if (feat/n).is_file())
    identities=[core.identity(p) for p in paths]
    for path in paths:
        if path.name.endswith('.nii.gz'):
            im=nib.load(path);core.require(len(im.shape)==3 and core.same_grid(im,reference),'Input packet is not on original group grid: '+str(path))
    signature=dict(source_files=identities,scalar_DOF_expanded_to_map=scalar_df,already_standard_space=True)
    marker=folder/'input-complete.json'
    if marker.is_file():
        old=json.loads(marker.read_text());core.require(old['signature']==signature and all(core.identity(v['path'])==v for v in old['outputs']), 'Staged input/source changed: '+str(folder))
        return folder/'stats'/name,dict(effect=identities[0],variance=identities[1],packet_sources=identities,staged_files=old['outputs'])
    if folder.exists():folder.rename(folder.with_name(folder.name+'.interrupted-'+uuid.uuid4().hex[:6]))
    (folder/'stats').mkdir(parents=True)
    for source,dest in ((effect,folder/'stats'/name),(variance,folder/'stats'/('var'+name)),(mask,folder/'mask.nii.gz'),(background,folder/'bg_image.nii.gz')):
        shutil.copyfile(source,dest)
    dof_out=folder/'stats'/f'tdof_t{number}.nii.gz'
    if scalar_df:
        value=float(df.read_text().strip());core.require(np.isfinite(value) and value>0,'Invalid original residual DOF')
        im=nib.load(effect);dof_image=nib.Nifti1Image(np.full(im.shape,value,dtype='float32'),im.affine,im.header.copy());dof_image.set_data_dtype('float32');nib.save(dof_image,dof_out)
    else:shutil.copyfile(df,dof_out)
    for source in paths[5:]:shutil.copyfile(source,folder/source.name)
    (folder/'design.lev').write_text('2\n')
    outputs=[core.identity(p) for p in folder.rglob('*') if p.is_file()]
    save(marker,dict(signature=signature,outputs=outputs))
    return folder/'stats'/name,dict(effect=identities[0],variance=identities[1],packet_sources=identities,staged_files=outputs)


def plan(audit, output, record):
    audit_state = json.loads((audit/'state.json').read_text())
    core.require(audit_state['UGR_runs_audited'] == 249 and audit_state['UGR_runs_failed'] == 0 and
                 len(audit_state['TR_mismatches']) == 1 and audit_state['TR_mismatches'][0]['key'] == 'sub-10817_ses-01_task-ugr_run-1',
                 'Unexpected original UGR audit state')
    original_contrasts(audit)
    manifest = json.loads((audit/'manifest.json').read_text())
    audited = {str(Path(r['source']).resolve()):r['identity'] for r in manifest if 'identity' in r}
    correction = json.loads((TR_RECORD/'corrected-analysis-input.json').read_text())
    corrected = Path(correction['active_feat'])
    core.require('corrected_TR_1p710.feat' == corrected.name, 'Wrong timing-control fit selected')
    settings = core.fsl_settings(corrected/'design.fsf')
    core.require(np.isclose(float(settings['fmri(tr)']), 1.710) and float(settings['fmri(smooth)']) == 5,
                 'Corrected candidate must have TR=1.710 and smoothing=5')
    core.require(core.sha(corrected/'confoundevs.txt') == correction['confounds_sha256'], 'Corrected nuisance matrix changed')
    original = Path(correction['original_feat'])
    core.require(np.array_equal(core.fsl_matrix(corrected/'design.con'), core.fsl_matrix(original/'design.con')),
                 'Corrected timing fit contrast definitions changed')
    original_con = core.fsl_matrix(ROOT/'derivatives/fsl/sub-11031/ses-01/L1_task-ugr_model-3_type-act_run-2_sm-5.feat/design.con')
    run2_model = ROOT/'derivatives/fsl/sub-11031/ses-01/L1_task-ugr_model-3_type-act_run-2_sm-5.feat/design.mat'
    core.require(not any(r['requires_review'] for r in contrast_diagnostics(core.fsl_matrix(run2_model), original_con)), 'Run2 contrasts need review')
    freeze(record/'reused_L1_provenance.json',dict(correction_record=core.identity(TR_RECORD/'corrected-analysis-input.json'),
        corrected_sm5_model=[core.identity(corrected/name) for name in ('design.fsf','design.con','confoundevs.txt')],
        original_valid_run2_model=[core.identity(run2_model.parent/name) for name in ('design.fsf','design.mat','design.con','confoundevs.txt')],
        original_bad_run1_contrasts=sorted(BAD),policy=POLICY))
    for name, expected in audited.items():
        if '/sub-11031/' in name or '/sub-10817/' in name:
            core.require(core.identity(name) == expected, 'Audited affected model/event source changed: '+name)
    subjects = {'sub-'+s.removeprefix('sub-') for s in (ROOT/'code/sublist_n132.txt').read_text().split()}
    core.require(len(subjects) == 132, 'Unexpected manuscript subject list')
    cases = []; inputs = []
    for case in json.loads((audit/'group_cases.json').read_text()):
        if case['task'] != 'ugr': continue
        group = Path(case['path']); fsf = group/'design.fsf'; text = fsf.read_text(); old = core.fsl_settings(fsf)
        core.require(core.identity(fsf) == audited[str(fsf.resolve())], 'Original group FSF changed since audit')
        cope = int(re.search(r'-cope-(\d+)_', group.name)[1]); variant = re.search(r'-cov-(.+)-flame1$', group.parent.name)[1]
        core.require(cope in (11, 12) and variant in ('ones', 'full', 'covs-vt'), 'Unexpected main activation model')
        name = variant+'-cope-'+str(cope); base = output/'activation'/variant/('cope-'+str(cope)); destination = base.with_suffix('.gfeat')
        replacements = {'fmri(outputdir)':str(base)}; selected = []; cohort = []
        for key, oldpath in sorted(case['saved_input_paths'].items(), key=lambda kv:int(re.search(r'\d+', kv[0])[0])):
            subject = re.search(r'/sub-(\d+)/', oldpath)
            core.require(subject is not None, 'Unknown subject input'); cohort.append('sub-'+subject[1])
            target, reason = choose_input(oldpath, cope, corrected)
            staged, provenance=stage_input(target,output/'inputs'/cohort[-1]/f'cope-{cope}.feat',reason,nib.load(source_file(group,'mask.nii.gz')))
            replacements[key] = str(staged)
            item = dict(case=name, subject=cohort[-1], cope=cope, old_input=oldpath, selected_input=str(target),group_input=str(staged),reason=reason,**provenance)
            selected.append(item); inputs.append(item)
        core.require(len(cohort) == 132 and set(cohort) == subjects and len(set(cohort)) == 132, 'Group cohort/input order invalid')
        core.require(sum(r['reason'] != 'reuse_original_participant_estimate' for r in selected) == 2, 'Expected exactly two repaired group inputs')
        transport = None
        reference = old.get('fmri(regstandard)', '')
        if reference and not Path(reference).is_file() and not Path(reference+'.nii.gz').is_file():
            core.require(Path(reference).name == 'MNI152_T1_2mm_brain', 'Unknown machine-specific reference path')
            local = Path(os.environ['FSLDIR'])/'data/standard/MNI152_T1_2mm_brain.nii.gz'
            core.require(local.is_file(), 'Local standard reference absent')
            replacements['fmri(regstandard)'] = str(local); transport = dict(original=reference, local=core.identity(local))
        fsf_out = record/'fsfs'/(name+'.fsf'); changed = replace_settings(text, replacements)
        freeze(record/'sources'/(name+'.json'), dict(original_fsf=core.identity(fsf), replacements=replacements, reference_path_transport=transport))
        if fsf_out.exists(): core.require(fsf_out.read_text() == changed, 'Generated FSF changed')
        else: fsf_out.write_text(changed)
        shutil.copyfile(fsf, record/'sources'/(name+'-original.fsf'))
        model_ids = [core.identity(source_file(group, name)) for name in ('design.mat', 'design.con', 'design.grp')]
        reference_ids = [core.identity(group/'cope1.feat'/name) for name in ('mask.nii.gz','stats/zstat1.nii.gz','stats/zstat2.nii.gz')]
        core.require(all(audited.get(v['path']) == v for v in model_ids+reference_ids), 'Original group models/maps changed since audit')
        cases.append(dict(name=name, variant=variant, cope=cope, original_group=str(group), output=str(destination), fsf=str(fsf_out),
            signature=dict(policy=POLICY, generated_fsf=core.sha(fsf_out), source_models=model_ids, source_reference_maps=reference_ids, inputs=selected)))
    core.require({(c['variant'], c['cope']) for c in cases} == {(v,k) for v in ('ones','full','covs-vt') for k in (11,12)}, 'All six original main activation models required')
    freeze(record/'plan.json', cases)
    fields=('case','subject','cope','old_input','selected_input','group_input','reason')
    core.write_table(record/'selected_inputs.tsv',[dict({k:r[k] for k in fields},source_effect_sha256=r['effect']['sha256'],source_variance_sha256=r['variance']['sha256']) for r in inputs])
    return cases


def compare_group(case, record):
    source, target = Path(case['original_group']), Path(case['output'])
    for name in ('design.mat','design.con','design.grp'):
        a, b = core.fsl_matrix(source_file(source,name)), core.fsl_matrix(source_file(target,name))
        core.require(a.shape == b.shape and np.allclose(a,b,atol=1e-6,rtol=1e-6), 'Group model/covariates/contrasts changed: '+name)
    reports = []
    for sign in (1,2):
        before = nib.load(source/'cope1.feat/stats'/f'zstat{sign}.nii.gz')
        after = nib.load(target/'cope1.feat/stats'/f'zstat{sign}.nii.gz')
        core.require(core.same_grid(before,after), 'Repaired group-map geometry differs')
        mask = (nib.load(source/'cope1.feat/mask.nii.gz').get_fdata()>0)&(nib.load(target/'cope1.feat/mask.nii.gz').get_fdata()>0)
        old_mask=nib.load(source/'cope1.feat/mask.nii.gz').get_fdata()>0;new_mask=nib.load(target/'cope1.feat/mask.nii.gz').get_fdata()>0
        x,y = before.get_fdata()[mask],after.get_fdata()[mask]
        core.require(mask.any() and np.isfinite(y).all(), 'Invalid repaired group statistics')
        old_cluster = nib.load(source/'cope1.feat'/f'cluster_mask_zstat{sign}.nii.gz').get_fdata()!=0
        new_cluster = nib.load(target/'cope1.feat'/f'cluster_mask_zstat{sign}.nii.gz').get_fdata()!=0
        reports.append(dict(sign=sign, spatial_r=float(np.corrcoef(x,y)[0,1]), rms_difference=float(np.sqrt(np.mean((y-x)**2))),
            common_stat_voxels=int(mask.sum()),old_mask_voxels=int(old_mask.sum()),new_mask_voxels=int(new_mask.sum()),
            group_mask_dice=float(2*mask.sum()/(old_mask.sum()+new_mask.sum())),
            old_thresholded_voxels=int(old_cluster.sum()), new_thresholded_voxels=int(new_cluster.sum()),
            shared_thresholded_voxels=int((old_cluster&new_cluster).sum())))
        for name in (f'stats/zstat{sign}.nii.gz',f'cluster_mask_zstat{sign}.nii.gz',f'thresh_zstat{sign}.nii.gz'):
            path=target/'cope1.feat'/name; shutil.copyfile(path,record/'group_maps'/(case['name']+'-'+path.name))
    save(record/(case['name']+'-comparison.json'),reports)
    return reports


def run_case(case, record, executable):
    folder=Path(case['output']); marker=folder/'repair-complete.json'
    if marker.exists():
        cached=json.loads(marker.read_text())
        core.require(cached['signature']==case['signature'] and all(core.identity(r['path'])==r for r in cached['outputs']), 'Completed repaired group inputs/outputs changed')
        return dict(state='complete', reused_completed_group=True)
    if folder.exists(): folder.rename(folder.with_name(folder.name+'.interrupted-'+uuid.uuid4().hex[:6]))
    folder.parent.mkdir(parents=True,exist_ok=True)
    core.run([str(executable),case['fsf']],record/'commands'/case['name'])
    core.require(all(core.identity(v['path']) == v for item in case['signature']['inputs'] for v in item['packet_sources']+item['staged_files']),
                 'Selected or staged input changed during group fit')
    compare_group(case,record)
    outputs=[core.identity(folder/'cope1.feat'/name) for name in ('stats/zstat1.nii.gz','stats/zstat2.nii.gz','cluster_mask_zstat1.nii.gz','cluster_mask_zstat2.nii.gz','design.mat','design.con','design.grp')]
    save(marker,dict(signature=case['signature'],outputs=outputs))
    return dict(state='complete', reused_completed_group=False)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--jobs',type=int,default=8);p.add_argument('--publish',action='store_true');p.add_argument('--brief',action='store_true');p.add_argument('--new-record',action='store_true');args=p.parse_args()
    RECORDS.mkdir(exist_ok=True);existing=sorted(RECORDS.glob('linux2-ugr-activation-repair-*'))
    if args.brief:
        if not existing: print('No activation repair record.');return
        record=existing[-1];state=json.loads((record/'state.json').read_text());print('Record:',record.name);print('State:',state['state']);print('Groups complete:',sum(r.get('state')=='complete' for r in state['cases'].values()),'/6');print('Errors:',state.get('error','none'))
        processes=status.processes()
        runners=[v for v in (processes or []) if 'rerun_activation.py' in v['args'] and '--brief' not in v['args']]
        print('Runner:', 'unknown' if processes is None else ', '.join(str(v['pid']) for v in runners) or 'not running')
        for key,value in state['cases'].items():
            if value.get('error'):print(key+': '+value['error'])
        if (record/'publication-error.txt').exists():print('Publication error:',(record/'publication-error.txt').read_text())
        pointer=ROOT/'code/ugr_reanalysis/output/activation-launcher.logpath'
        if pointer.exists():print('Launcher errors:',status.log_errors(Path(pointer.read_text().strip())) or 'none found')
        return
    core.require(1<=args.jobs<=16,'Use 1-16 group workers')
    audit=sorted(RECORDS.glob('linux2-original-group-audit-*'))[-1]
    candidates=[r for r in existing if json.loads((r/'state.json').read_text()).get('policy')==POLICY and json.loads((r/'state.json').read_text()).get('source_audit')==str(audit)]
    (ROOT/'code/ugr_reanalysis/output').mkdir(exist_ok=True)
    with (ROOT/'code/ugr_reanalysis/output/activation-repair.lock').open('a') as handle:
        try:fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:print('UGR activation repair already running.');return
        if candidates and not args.new_record:
            record=candidates[-1];state=json.loads((record/'state.json').read_text());output=Path(state['output_root'])
        else:
            stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:6]
            record=RECORDS/('linux2-ugr-activation-repair-'+stamp);record.mkdir();output=ROOT/'derivatives/ugr_reanalysis'/stamp;output.mkdir(parents=True)
            state=dict(policy=POLICY,source_audit=str(audit),output_root=str(output),cases={},state='preparing_original_model_group_repair')
        for directory in ('fsfs','sources','group_maps','commands','errors'):(record/directory).mkdir(exist_ok=True)
        state.pop('error',None);state['last_started_UTC']=datetime.now(timezone.utc).isoformat();save(record/'state.json',state)
        try:
            core.require(output.resolve().is_relative_to((ROOT/'derivatives/ugr_reanalysis').resolve()),'Wrong repair output root')
            cases=plan(audit,output,record);executable=Path(os.environ['FSLDIR'])/'bin/feat';core.require(executable.is_file(),'FSL feat absent');lock=threading.Lock()
            def worker(case):
                try:result=run_case(case,record,executable)
                except Exception as exc:
                    result=dict(state='failed',error=str(exc));(record/'errors'/(case['name']+'.txt')).write_text(traceback.format_exc())
                    for index,path in enumerate((record/'commands'/case['name']/'commands.log',Path(case['output'])/'report.log',Path(case['output'])/'cope1.feat/report.log')):
                        if path.is_file():
                            with path.open('rb') as stream:
                                stream.seek(max(0,path.stat().st_size-24000));tail=stream.read().decode(errors='replace')
                            (record/'errors'/(case['name']+f'-log-tail-{index}.txt')).write_text(tail)
                with lock:
                    state['cases'][case['name']]=result;save(record/'state.json',state);print(case['name']+': '+result['state'],flush=True)
                return result
            pilot=next(c for c in cases if c['variant']=='ones' and c['cope']==11)
            core.require(worker(pilot)['state']=='complete','Main fairness group pilot failed; no expansion')
            state['state']='pilot_complete_expanding_original_group_models';save(record/'state.json',state)
            with ThreadPoolExecutor(max_workers=args.jobs) as pool:
                for future in as_completed([pool.submit(worker,c) for c in cases if c!=pilot]):future.result()
            core.require(len(state['cases'])==6 and all(c['state']=='complete' for c in state['cases'].values()),'Some repaired groups failed')
            state['state']='original_UGR_activation_groups_complete_review_pending';save(record/'state.json',state)
        except Exception as exc:
            state.update(state='UGR_activation_repair_failed',error=str(exc));save(record/'state.json',state);(record/'failure.txt').write_text(traceback.format_exc());print('ERROR:',exc,flush=True)
        (record/'README.md').write_text('# Original-model UGR activation repair\n\n'
            'Original group designs/covariates/thresholds retained. Reuses corrected TR=1.710 sm-5 sub-10817 single-run estimates and original sub-11031 run2 estimates for affected contrasts, '
            'with all other participant estimates unchanged. All N132 subjects retained; no common-slope model substitution. '
            'All group paths point to isolated copies of small 3D effect/variance/DOF/mask/background files, preserving source directories from FEAT registration-cache writes. '
            'Existing standard-space inputs are reused; no resampling occurs. Scalar original L1 DOF is expanded to a constant 3D map when needed. '
            'See plan.json, selected_inputs.tsv, original/generated FSFs and case comparison JSONs. No preprocessing, new L1/L2 fitting, Trust, PPI, seed replacement, or manuscript edits. '
            'Group images remain under the separately versioned derivatives/ugr_reanalysis root; compact review artifacts are in this record. '
            'A source standard-reference path from a different machine may be mapped to the same FSL MNI152 2-mm basename locally; this is recorded and the generated group matrices must match the original.\n')
        if args.publish:
            try:publish(record,'Recompute original-model UGR activation groups with explicit timing and valid-run inputs')
            except Exception as exc:(record/'publication-error.txt').write_text(str(exc));print('Git publication failed:',exc,flush=True)
        print('Final state:',state['state'],'| record:',record,flush=True)
        raise SystemExit(0 if state['state']=='original_UGR_activation_groups_complete_review_pending' else 1)


if __name__=='__main__':main()
