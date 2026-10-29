#!/usr/bin/env python3
"""Record exact raw/retained TR evidence for one UGR run; no preprocessing or edits."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import uuid

import nibabel as nib

import insula_ica as core
from run_bold_controls import publish, save


SOURCE=core.ROOT/'code/insula_ica/records/linux2-paired-group-ica-20261028T215022Z-679a3e'
RAW=Path('/ZPOOL/data/projects/rf1-sra-linux2/bids')


def header(path):
    if not path.is_file():return dict(path=str(path),present=False)
    image=nib.load(path);units=image.header.get_xyzt_units()
    tr=float(image.header.get_zooms()[3]) if len(image.shape)==4 else None
    multiplier={'sec':1.,'msec':.001,'usec':.000001}.get(units[1])
    return dict(path=str(path),present=True,shape=list(image.shape),header_TR=tr,time_units=units[1],
                header_TR_seconds=tr*multiplier if tr is not None and multiplier else None,
                identity=core.identity(path,False))


def metadata(path):
    if not path.is_file():return dict(path=str(path),present=False)
    data=json.loads(path.read_text())
    return dict(path=str(path),present=True,sha256=core.sha(path),
                timing={k:data[k] for k in ('RepetitionTime','VolumeTiming','EchoTime','SliceTiming','DelayTime') if k in data})


def collect(row, raw=RAW):
    sub='sub-'+row['subject'].removeprefix('sub-')
    func=raw/sub/'ses-01/func'
    echoes=[]
    for echo in range(1,5):
        name=f'{row["key"]}_echo-{echo}_part-mag_bold'
        echoes.append(dict(echo=echo,image=header(func/(name+'.nii.gz')),json=metadata(func/(name+'.json'))))
    feat=Path(row['active_feat'])
    settings=core.fsl_settings(feat/'design.fsf')
    historical=Path(settings['feat_files(1)'])
    return dict(key=row['key'],expected_TR_seconds=row['tr'],expected_volumes=row['nvol'],
                original_FEAT_TR=settings.get('fmri(tr)'),echoes=echoes,
                retained_filtered_BOLD=header(Path(row['filtered_data'])),
                historical_preproc_BOLD=header(historical),
                historical_preproc_JSON=metadata(historical.with_name(historical.name.removesuffix('.nii.gz')+'.json')),
                original_design=str(feat/'design.fsf'),original_design_sha256=core.sha(feat/'design.fsf'),
                saved_confounds_sha256=core.sha(row['confounds']),
                expected_confounds_sha256=row['confounds_sha256'],
                actions=dict(preprocessing=False,metadata_modified=False,TR_overridden=False))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--key',default='sub-10817_ses-01_task-ugr_run-1')
    parser.add_argument('--publish',action='store_true')
    args=parser.parse_args()
    rows=json.loads((SOURCE/'accepted-inputs.json').read_text())
    selected=[r for r in rows if r['key']==args.key and r['task']=='ugr']
    core.require(len(selected)==1,'Unknown/ambiguous accepted UGR key')
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:6]
    record=core.ROOT/'code/insula_ica/records'/('linux2-ugr-tr-diagnostic-'+stamp)
    record.mkdir()
    evidence=collect(selected[0]);save(record/'evidence.json',evidence)
    (record/'code-commit.txt').write_text(subprocess.check_output(['git','rev-parse','HEAD'],cwd=core.ROOT,text=True))
    print('Run:',args.key,'| expected TR:',evidence['expected_TR_seconds'],'| original FEAT TR:',evidence['original_FEAT_TR'])
    for item in evidence['echoes']:
        im=item['image'];js=item['json']
        print(f'Echo {item["echo"]}: JSON TR={js.get("timing",{}).get("RepetitionTime","missing")}; '
              f'NIfTI TR={im.get("header_TR","missing")} {im.get("time_units","")}; volumes={im.get("shape",[None])[-1]}')
    kept=evidence['retained_filtered_BOLD']
    print('Retained sm-5 BOLD header TR:',kept.get('header_TR','missing'),kept.get('time_units',''))
    print('Historical preprocessing JSON TR:',evidence['historical_preproc_JSON'].get('timing',{}).get('RepetitionTime','missing'))
    print('Saved confounds unchanged:',evidence['saved_confounds_sha256']==evidence['expected_confounds_sha256'])
    if args.publish:
        publish(record,'Record exact TR evidence for the one blocked UGR run without changing inputs')
        print('Diagnostic pushed to insula-ica.')
    print('Saved:',record)


if __name__=='__main__':main()
