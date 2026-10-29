"""Paired right-insula dim-5 MELODIC-native variance-normalization pilot."""
import json
from pathlib import Path
import traceback

import nibabel as nib
import numpy as np

import run_group_ica as g


AMPLITUDE=1000.0


def prepare_inputs(root,record,mode):
    pilot=root/'native_vn_pilot';folder=pilot/'sessions'/mode/'right';folder.mkdir(parents=True,exist_ok=True)
    (pilot/'masks').mkdir(exist_ok=True)
    for name in ('right.nii.gz','right_crop.nii.gz'):
        im=nib.load(root/'masks'/name);g.freeze_nifti(pilot/'masks'/name,np.asanyarray(im.dataobj),im.affine)
    image=nib.load(root/'reference_crop.nii.gz');g.freeze_nifti(pilot/'reference_crop.nii.gz',np.asanyarray(image.dataobj),image.affine)
    original=json.loads((root/'sessions'/mode/'right/verified-inputs.json').read_text())
    checks=[];listing=[]
    for check in original['inputs']:
        source=check['output'];destination=folder/Path(source['path']).name
        g.core.require(g.core.identity(source['path'])==source,'VN pilot source changed')
        marker=destination.with_name(destination.name+'.complete.json')
        signature=dict(original=source,global_amplitude_multiplier=AMPLITUDE)
        if marker.is_file():
            old=json.loads(marker.read_text());verification=old['verification']
            g.core.require(old['signature']==signature and g.core.identity(destination)==verification['output'],'Scaled VN pilot input changed')
        else:
            g.core.require(not destination.exists(),'Partial VN pilot input exists; preserve it')
            image=nib.load(source['path']);data=np.asarray(image.dataobj,dtype=np.float32)*AMPLITUDE
            g.write_nifti(destination,data,image.affine)
            verification=g.validate_masked_4d(destination,nib.load(pilot/'masks/right_crop.nii.gz'),data.shape[3])
            g.save(marker,dict(signature=signature,verification=verification))
        checks.append(verification);listing.append(str(destination.resolve()))
    (folder/'inputs.txt').write_text('\n'.join(listing)+'\n')
    g.save(folder/'verified-inputs.json',dict(mode=mode,roi='right',mask_sha256=g.core.sha(pilot/'masks/right_crop.nii.gz'),inputs=checks))
    verified=g.verify_ica_inputs(mode,'right',pilot)
    g.save(record/f'native-vn-{mode}-inputs.json',dict(roi='right',mode=mode,inputs=264,global_amplitude_multiplier=AMPLITUDE,
        source_collection_sha256=g.core.sha(root/'sessions'/mode/'right/verified-inputs.json'),verified_inputs_sha256=verified,
        rationale='Global multiplication preserves relative run/voxel weighting and avoids the absolute 0.01 SD floor in classic joint MELODIC variance normalization. No voxelwise mean is added.'))
    return pilot


def run_pilot(root,record,state,masks,slices,reference,baseline):
    from resume_group_ica import retry_fit
    state['native_variance_normalization_pilot']={}
    gate=g.inputs.MemoryGate()
    for mode in g.inputs.MODES:
        try:
            with gate.reserve(2):
                pilot=prepare_inputs(root,record,mode)
                result=retry_fit((mode,'right',5),pilot,masks,slices,reference,record,variance_normalization=True)
            key=f'{mode}/right/dim-5';existing=baseline.get(key)
            if existing:
                a=nib.load(existing['oriented_maps']).get_fdata()[masks['right']]
                b=nib.load(result['oriented_maps']).get_fdata()[masks['right']]
                result['baseline_spatial_comparison']=g.compare_maps(a,b,labels=('manual_SD','native_MELODIC_VN'))
                result['comparison_note']='Same residual representation; baseline manual runwise SD scaling versus additional native MELODIC normalization.'
            g.component_figure(nib.load(result['oriented_maps']).get_fdata(),masks['right'],reference.get_fdata(),reference.affine,
                record/'figures'/f'native_vn_pilot_{mode}_right_dim-5.png',f'Native MELODIC variance-normalization pilot | {mode} | right | 5 components')
        except Exception as exc:
            result=dict(state='failed',error=str(exc))
            (record/'errors'/f'native_vn_pilot_{mode}.txt').write_text(traceback.format_exc())
        state['native_variance_normalization_pilot'][mode]=result;g.save(record/'state.json',state)
        print('Native MELODIC VN pilot '+mode+': '+result['state'],flush=True)
    g.save(record/'native-variance-normalization-pilot.json',state['native_variance_normalization_pilot'])
