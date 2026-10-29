#!/usr/bin/env python3
"""Publish already computed group review volumes skipped by the old ignore rule."""
import argparse
from pathlib import Path
import json

import nibabel as nib
import numpy as np

import insula_ica as core
from run_bold_controls import publish, save


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--record',type=Path);p.add_argument('--publish',action='store_true');args=p.parse_args()
    records=core.ROOT/'code/insula_ica/records'
    if args.record is None:
        choices=sorted(r for r in records.glob('linux2-completed-ICA-exploration-*') if (r/'state.json').is_file())
        core.require(bool(choices),'No completed-ICA exploration record');args.record=choices[-1]
    record=args.record.resolve();core.require(record.is_relative_to(records.resolve()),'Expected a project analysis record')
    state=json.loads((record/'state.json').read_text())
    core.require(state['state']=='exploratory_connectivity_and_task_review_ready','Exploration has not completed')
    expected=[record/'group_maps/wholebrain_common_mask.nii.gz']
    for mode in ('joint_bold','native_film'):
        for roi in core.ROIS:
            expected.extend([record/'primary_maps'/f'{mode}_{roi}_dim-5.nii.gz',record/'group_maps'/f'{mode}_{roi}_pooled_mean.nii.gz'])
    reference=nib.load(expected[0]);identities=[]
    core.require(len(reference.shape)==3,'Coverage artifact must be a 3D mask')
    for path in expected:
        core.require(path.is_file() and path.stat().st_size<50*1024**2,'Missing/oversized group review map: '+str(path))
        image=nib.load(path)
        core.require(core.same_grid(image,reference),'Review map geometry differs')
        core.require((len(image.shape)==3 and path==expected[0]) or (len(image.shape)==4 and image.shape[3]==5),
                     'Only group five-component maps and coverage mask may be published')
        core.require(np.isfinite(np.asanyarray(image.dataobj)).all(),'Nonfinite review map')
        identities.append(core.identity(path))
    save(record/'group_volume_export.json',dict(files=identities,group_derived_maps_only=True,analysis_rerun=False))
    if args.publish:publish(record,'Publish existing group ICA/connectivity volumes and common coverage mask')
    print(f'{len(identities)} existing group review volumes validated'+(' and pushed.' if args.publish else '.'))
    print('No ICA, FEAT, connectivity or task analysis was rerun.')


if __name__=='__main__':main()
