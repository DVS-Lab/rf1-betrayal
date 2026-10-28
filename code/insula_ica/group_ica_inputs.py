"""Frozen accepted input selection and paired residual representations for ICA."""
from collections import Counter
from contextlib import contextmanager
from pathlib import Path
import json
import shutil
import subprocess
import threading

import nibabel as nib
import numpy as np

import insula_ica as core
import feat_residuals as residuals
from rebuild_trust import compare_fit
from run_bold_controls import save


BASE='linux2-trust-regeneration-20261027T162428Z-ef53c5'
STOCK='linux2-trust-stock-filter-20261028T203335Z-c78e1b-resume'
REPAIR='linux2-registration-fix-sub-10834-20261028T151154Z-acdf91'
ACCEPTED_FLAGS={f'sub-{s}_ses-01_task-trust_run-2' for s in ('10720','11326','11433')}
MODES=('joint_bold','native_film')


def select_trust(records):
    """Select corrected run-2 inputs and repaired sub-10834 pair explicitly."""
    baseline=json.loads((records/BASE/'input-registry.json').read_text())
    completed=json.loads((records/STOCK/'state.json').read_text())
    core.require(completed['completed']==108 and completed['failed']==0,'Stock-reference queue incomplete')
    replacement={v['key']:v for v in completed['subjects'].values()}
    repair=json.loads((records/REPAIR/'review-decision.json').read_text())
    core.require(repair['assessment']=='accept_repaired_pair_for_sub_10834_analysis_use','Unreviewed repair')
    repaired={r['key']:r for r in repair['selected_inputs']}
    core.require(len(baseline)==241 and len(replacement)==108 and len(repaired)==2,'Unexpected Trust input counts')
    rows=[]
    for item in baseline:
        key=item['key'];source=records/BASE;row=dict(item)
        if key in replacement:
            row=dict(replacement[key]['row']);source=records/STOCK
            core.require(replacement[key]['reference_selection']['status']=='same_run','Wrong-run reference')
        if key in repaired:
            row['new_bold']=repaired[key]['bold'];source=records/REPAIR
        qc_path=source/'qc'/(key+'.json');qc=json.loads(qc_path.read_text())
        core.require(qc['finite_nonconstant_residuals'] and qc['confounds']['exact_numeric_equality'],'Trust QC failed')
        core.require(not qc['review_flags'] or key in ACCEPTED_FLAGS,'New/unreviewed Trust flag: '+key)
        feat=Path(qc['residual']).parents[1]
        rows.append(dict(key=key,subject=row['subject'].removeprefix('sub-'),task='trust',run=int(row['run']),
            nvol=int(row['nvol']),tr=float(row['tr']),nss=0,original_feat=row['original_feat'],
            active_feat=str(feat),source_fsf=row['original_fsf'],filtered_data=qc['filtered_data'],native_residual=qc['residual'],
            confounds=row['confounds'],confounds_sha256=row['confounds_sha256'],source_qc=str(qc_path),
            original_flags=qc['review_flags'],source_selection_record=str(source),contrast_index=10))
    return rows


def select_inputs(audit):
    records=core.ROOT/'code/insula_ica/records'
    rows=select_trust(records)
    manifest=core.read_table(audit/'manifest.tsv')
    ugr=[r for r in manifest if r['task']=='ugr' and r['status']=='included']
    core.require(len(ugr)==249,'Expected all 249 accepted UGR runs in frozen audit')
    for r in ugr:
        old=Path(r['original_feat'])
        rows.append(dict(key=r['key'],subject=str(r['subject']).removeprefix('sub-'),task='ugr',run=int(r['run']),
            nvol=int(r['nvol']),tr=float(r['tr']),nss=0,original_feat=str(old),active_feat=str(old),
            source_fsf=r.get('source_fsf',str(old/'design.fsf')),filtered_data=r['filtered_data'],native_residual=str(old/'stats/res4d.nii.gz'),
            confounds=r['confound_matrix'],confounds_sha256=core.sha(r['confound_matrix']),
            source_selection_record=str(audit),contrast_index=11,original_flags=[]))
    cohort=set((core.ROOT/'code/sublist_n132.txt').read_text().split())
    core.require(len(cohort)==132 and len(rows)==len({r['key'] for r in rows})==490,'Cohort/run identity mismatch')
    for task in ('trust','ugr'):
        core.require({r['subject'] for r in rows if r['task']==task}==cohort,'Task does not cover N=132')
    core.require(Counter(r['task'] for r in rows)==dict(trust=241,ugr=249),'Task counts changed')
    return sorted(rows,key=lambda r:(r['task'],r['subject'],r['run']))


def inspect_source(row,reference):
    feat=Path(row['active_feat']);historical=Path(row['original_feat'])
    design=core.fsl_matrix(feat/'design.mat');contrasts=core.fsl_matrix(feat/'design.con')
    core.require(len(design)==row['nvol'] and design.shape[1]>0,'Data/model timing mismatch')
    for name in ('design.mat','design.con'):
        a=core.fsl_matrix(historical/name);b=core.fsl_matrix(feat/name)
        core.require(a.shape==b.shape and np.allclose(a,b,rtol=1e-8,atol=1e-8),'Full joint model changed')
    settings=core.fsl_settings(feat/'design.fsf' if (feat/'design.fsf').exists() else Path(row['source_fsf']))
    residuals.validate_model(settings,design,contrasts,row['task'])
    core.require(settings['fmri(prewhiten_yn)']=='1','Unexpected non-FILM scale')
    core.require(int(settings['fmri(ndelete)'])==0,'Nonzero discarded volumes')
    core.require(core.sha(row['confounds'])==row['confounds_sha256'],'Original confounds changed')
    nuisance=core.load_numeric(row['confounds']);fitted=core.load_numeric(feat/'confoundevs.txt')
    core.require(nuisance.shape==fitted.shape and np.array_equal(nuisance,fitted),'Actual fitted nuisance matrix differs')
    core.require(len(nuisance)==len(design),'Nuisance timing differs')
    core.require(float((feat/'stats/dof').read_text())>=37,'Insufficient residual degrees of freedom')
    image=nib.load(row['filtered_data']);mask=nib.load(feat/'mask.nii.gz')
    core.require(core.same_grid(reference,image) and core.same_grid(reference,mask),'Incompatible common MNI grid')
    core.require(image.shape==(*reference.shape,len(design)) and np.isclose(image.header.get_zooms()[3],1.615,atol=1e-5),'Wrong grid/TR/time count')
    pe=[feat/f'stats/pe{k+1}.nii.gz' for k in range(design.shape[1])]
    paths=[Path(row['filtered_data']),feat/'design.mat',feat/'design.con',feat/'mask.nii.gz',feat/'mean_func.nii.gz',
           feat/'stats/sigmasquareds.nii.gz',feat/'stats/dof',feat/'confoundevs.txt',Path(row['confounds']),*pe]
    if Path(row['native_residual']).is_file():paths.append(Path(row['native_residual']))
    core.require(all(p.is_file() for p in paths),'Missing retained fit inputs: '+row['key'])
    identities=[core.identity(p,content=not p.name.endswith('.nii.gz')) for p in paths]
    return design,identities


def verify_frozen_ugr(rows,audit):
    inventory=json.loads((audit/'audit.json').read_text())
    frozen={p['path']:p for p in inventory['inputs']}
    for row in rows:
        if row['task']!='ugr':continue
        old=Path(row['original_feat'])
        for path in (old/'design.mat',old/'design.con',old/'mask.nii.gz',old/'mean_func.nii.gz',
                     Path(row['filtered_data']),Path(row['confounds'])):
            key=str(path.resolve());core.require(key in frozen,'UGR source not in original frozen inventory: '+key)
            item=frozen[key]
            core.require(core.identity(path,content=item['sha256'] is not None)==item,'Frozen UGR source changed: '+key)


def ensure_native(row,root,require_autocorr=False):
    """Reuse native residuals or replay only the retained UGR FILM fit."""
    source=Path(row['native_residual'])
    if source.is_file() and (not require_autocorr or (Path(row['active_feat'])/'stats/threshac1.nii.gz').is_file()):
        return Path(row['active_feat']),False
    core.require(row['task']=='ugr','Accepted Trust residual disappeared')
    original=Path(row['original_feat']);target=root/'ugr_film'/(row['key']+'.feat')
    identity=dict(key=row['key'],design_sha256=core.sha(original/'design.mat'),
                  filtered_identity=core.identity(row['filtered_data'],content=False))
    marker=target/'recovery_complete.json'
    if marker.is_file():
        core.require(json.loads(marker.read_text())==identity and (target/'stats/res4d.nii.gz').is_file(),
                     'Completed UGR recovery inputs changed')
        return target,True
    core.require(not target.exists(),'Partial UGR FILM recovery exists; preserve it for inspection')
    target.mkdir(parents=True)
    argv,log_source=residuals.replay_command(original,target,True)
    for name in ('design.mat','design.con','design.fts','design.fsf','confoundevs.txt','mask.nii.gz','mean_func.nii.gz'):
        if (original/name).is_file():shutil.copy2(original/name,target/name)
    plan=dict(argv=argv,cwd=str(original),historical_command_source=str(log_source),operation='same-input FILM statistics only')
    save(target/'replay.json',plan)
    with (target/'film.log').open('w') as stream:
        subprocess.run(argv,cwd=original,stdout=stream,stderr=subprocess.STDOUT,check=True)
    compare_fit(row,target,strict_historical=True)
    # Both residual representations must use the same joint-fit coefficients.
    design=core.fsl_matrix(original/'design.mat')
    for k in range(design.shape[1]):
        old=nib.load(original/f'stats/pe{k+1}.nii.gz');new=nib.load(target/f'stats/pe{k+1}.nii.gz')
        core.require(core.same_grid(old,new) and np.allclose(old.get_fdata(),new.get_fdata(),rtol=1e-4,atol=1e-4),
                     'UGR recovered coefficients differ from retained fit')
    save(marker,identity)
    return target,True


def joint_bold_residual(data,design,betas):
    """Original-scale residual of an existing joint GLS fit, without refitting."""
    data=np.asarray(data,dtype=np.float64);design=np.asarray(design,dtype=np.float64)
    betas=np.asarray(betas,dtype=np.float64)
    core.require(data.ndim==2 and data.shape[1]==len(design) and betas.shape==(data.shape[0],design.shape[1]),'Residual decomposition dimensions disagree')
    core.require(np.isfinite(data).all() and np.isfinite(design).all() and np.isfinite(betas).all(),'Nonfinite decomposition inputs')
    return (data-data.mean(axis=1,keepdims=True))-betas@design.T


class MemoryGate:
    def __init__(self):
        self.condition=threading.Condition();self.reserved=0

    @contextmanager
    def reserve(self,amount):
        with self.condition:
            while True:
                mem=Path('/proc/meminfo')
                available=next(int(x.split()[1])/1024**2 for x in mem.read_text().splitlines() if x.startswith('MemAvailable:')) if mem.exists() else 32
                if available-self.reserved>=amount+8:
                    self.reserved+=amount;break
                self.condition.wait(timeout=5)
        try:yield
        finally:
            with self.condition:
                self.reserved-=amount;self.condition.notify_all()


def prepare_run(row,root,record,union,reference):
    folder=record/'runs';folder.mkdir(exist_ok=True)
    design,before=inspect_source(row,reference)
    destination=root/'run_arrays';destination.mkdir(exist_ok=True)
    marker=destination/(row['key']+'_complete.json')
    if marker.is_file():
        previous=json.loads(marker.read_text())
        core.require(previous['qc']['source_identities']==before,'Cached residual input changed')
        core.require(all(core.identity(x['path'])==x for x in previous['arrays']),'Cached residual arrays changed')
        native_identity=previous['qc']['native_residual_identity']
        core.require(core.identity(native_identity['path'],content=False)==native_identity,'Native residual source changed')
        save(folder/(row['key']+'.json'),previous['qc'])
        return previous['qc']
    native_feat,replayed=ensure_native(row,root)
    native=nib.load(native_feat/'stats/res4d.nii.gz')
    core.require(core.same_grid(native,reference) and native.shape==(*reference.shape,len(design)),'Native residual geometry mismatch')
    # Temporary whole-volume reads are bounded by the preparation memory gate.
    source=nib.load(row['filtered_data']);data=np.asarray(source.dataobj,dtype=np.float32)[union].astype(float)
    pe=np.column_stack([np.asarray(nib.load(Path(row['active_feat'])/f'stats/pe{k+1}.nii.gz').dataobj)[union]
                        for k in range(design.shape[1])])
    bold=joint_bold_residual(data,design,pe)
    native_data=np.asarray(native.dataobj,dtype=np.float32)
    core.require(np.isfinite(native_data).all(),'Nonfinite native residuals')
    innovation=native_data[union].astype(float);del native_data
    mask=np.asarray(nib.load(Path(row['active_feat'])/'mask.nii.gz').dataobj)>0
    sigma=np.asarray(nib.load(Path(row['active_feat'])/'stats/sigmasquareds.nii.gz').dataobj)
    fitted=(mask & np.isfinite(sigma) & (sigma>0))[union]
    valid=fitted & np.isfinite(bold).all(axis=1) & np.isfinite(innovation).all(axis=1)
    valid &= (bold.std(axis=1)>1e-8)&(innovation.std(axis=1)>1e-8)
    core.require(valid.any(),'No common valid insula residuals')
    array_paths=[destination/(row['key']+'_'+m+'.npy') for m in MODES]+[destination/(row['key']+'_valid.npy')]
    core.require(not any(p.exists() for p in array_paths),'Partial residual arrays exist; preserve them for inspection')
    for mode,values in zip(MODES,(bold,innovation)):
        np.save(destination/(row['key']+'_'+mode+'.npy'),values.astype('float32'))
    np.save(destination/(row['key']+'_valid.npy'),valid)
    indices=np.flatnonzero(valid)[np.linspace(0,int(valid.sum())-1,min(32,int(valid.sum()))).astype(int)]
    def acf1(values):
        return [float(np.corrcoef(v[:-1],v[1:])[0,1]) for v in values[indices]
                if np.std(v[:-1])>1e-12 and np.std(v[1:])>1e-12]
    bold_acf=acf1(bold);native_acf=acf1(innovation)
    q=core.basis(design)
    ordinary=[float(np.linalg.norm(q.T@v)/max(np.linalg.norm(v),1e-12)) for v in bold[indices]]
    qc=dict(key=row['key'],native_feat=str(native_feat),native_recovered_by_film=replayed,
        native_residual_identity=core.identity(native_feat/'stats/res4d.nii.gz',content=False),
        source_identities=before,valid_voxels=int(valid.sum()),
        sampled_voxel_indices=indices.tolist(),median_joint_bold_acf1=float(np.median(bold_acf)) if bold_acf else None,
        median_native_film_acf1=float(np.median(native_acf)) if native_acf else None,
        median_paired_temporal_r=float(np.median([np.corrcoef(a,b)[0,1] for a,b in zip(bold[indices],innovation[indices])])),
        max_ordinary_design_projection_of_joint_GLS_residual=max(ordinary),
        ordinary_projection_note='Diagnostic only: a joint GLS residual is orthogonal under the fitted whitening metric, not necessarily the ordinary Euclidean metric.',
        no_new_regression=True,no_additional_temporal_filter=True,
        representation_note='Arrays contain insula voxels only; whole-brain connectivity remains a separate stage using these same fitted-model definitions.')
    core.require(before==[core.identity(x['path'],content=x['sha256'] is not None) for x in before],'Historical fitted inputs changed')
    save(folder/(row['key']+'.json'),qc)
    save(marker,dict(qc=qc,arrays=[core.identity(p) for p in array_paths]))
    return qc


def crop_geometry(reference,mask):
    indices=np.argwhere(mask);core.require(len(indices)>0,'Empty crop mask')
    lo=indices.min(axis=0);hi=indices.max(axis=0)+1
    slices=tuple(slice(int(a),int(b)) for a,b in zip(lo,hi))
    translation=np.eye(4);translation[:3,3]=lo
    return slices,reference.affine@translation


def balanced_session(arrays):
    """Equal run energy within one participant/task input; keep every frame."""
    blocks=[];weights=[]
    for values in arrays:
        y=np.asarray(values,dtype=float);y=y-y.mean(axis=1,keepdims=True)
        sd=y.std(axis=1,keepdims=True);core.require((sd>1e-8).all(),'Constant ICA input')
        weight=1/np.sqrt(len(arrays)*y.shape[1]);blocks.append((y/sd*weight).astype('float32'));weights.append(float(weight))
    return np.concatenate(blocks,axis=1),weights
