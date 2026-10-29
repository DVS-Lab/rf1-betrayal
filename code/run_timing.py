#!/usr/bin/env python3
"""Validate each run's acquisition TR before FEAT can overwrite its header."""
import argparse
import json
from pathlib import Path
import re
import sys

import nibabel as nib
import numpy as np


def image_tr(path):
    image=nib.load(path)
    if len(image.shape)!=4:raise ValueError('Expected 4D BOLD: '+str(path))
    unit=image.header.get_xyzt_units()[1]
    factor={'sec':1.,'msec':.001,'usec':.000001}.get(unit)
    if factor is None:raise ValueError('Unknown NIfTI time units: '+str(path))
    return float(image.header.get_zooms()[3])*factor,image.shape[3]


def acquisition_timing(folder,key):
    echoes=[]
    for echo in range(1,5):
        stem=f'{key}_echo-{echo}_part-mag_bold'
        image=Path(folder)/(stem+'.nii.gz');sidecar=Path(folder)/(stem+'.json')
        meta=json.loads(sidecar.read_text());tr=float(meta['RepetitionTime'])
        header,frames=image_tr(image)
        if not np.isfinite(tr) or tr<=0 or not np.isclose(tr,header,rtol=0,atol=1e-5):
            raise ValueError(f'{key} echo {echo}: JSON TR={tr}, NIfTI TR={header}; investigate before FEAT')
        if 'SliceTiming' in meta and not all(np.isfinite(x) and 0<=x<tr for x in meta['SliceTiming']):
            raise ValueError('SliceTiming outside the acquisition TR: '+str(sidecar))
        echoes.append(dict(echo=echo,tr_seconds=tr,header_tr_seconds=header,volumes=frames,
                           image=str(image),json=str(sidecar),echo_time_seconds=meta.get('EchoTime'),
                           protocol={k:meta[k] for k in ('ProtocolName','SequenceName','SeriesDescription',
                               'ScanningSequence','SequenceVariant','MultibandAccelerationFactor','SliceThickness') if k in meta}))
    if any(not np.isclose(x['tr_seconds'],echoes[0]['tr_seconds'],rtol=0,atol=1e-5) or
           x['volumes']!=echoes[0]['volumes'] for x in echoes):
        raise ValueError('Echo TRs/time counts disagree: '+key)
    return dict(key=key,tr_seconds=echoes[0]['tr_seconds'],volumes=echoes[0]['volumes'],echoes=echoes)


def fsf_tr(text):
    values=re.findall(r'^set\s+fmri\(tr\)\s+([^\s#]+)',text,re.M)
    if len(values)!=1:raise ValueError('Expected exactly one numeric fmri(tr)')
    return float(values[0].strip('"'))


def stamp(text,tr):
    if not np.isfinite(tr) or tr<=0:raise ValueError('TR must be finite and positive')
    fsf_tr(text)
    return re.sub(r'^set\s+fmri\(tr\)\s+.*$',f'set fmri(tr) {tr:.9g}',text,flags=re.M)


def main():
    parser=argparse.ArgumentParser(description=__doc__);commands=parser.add_subparsers(dest='command',required=True)
    resolve=commands.add_parser('resolve');resolve.add_argument('--raw-func',type=Path,required=True)
    resolve.add_argument('--key',required=True);resolve.add_argument('--bold',type=Path,required=True)
    for name in ('stamp','check-existing'):
        sub=commands.add_parser(name);sub.add_argument('--fsf',type=Path,required=True);sub.add_argument('--tr',type=float,required=True)
    args=parser.parse_args()
    try:
        if args.command=='resolve':
            timing=acquisition_timing(args.raw_func,args.key);tr,volumes=image_tr(args.bold)
            if not np.isclose(tr,timing['tr_seconds'],rtol=0,atol=1e-5) or volumes!=timing['volumes']:
                raise ValueError(f'Pre-FEAT BOLD timing disagrees with raw acquisition: {args.bold}; '
                                 f'BOLD TR={tr}, acquisition TR={timing["tr_seconds"]}')
            print(f'{timing["tr_seconds"]:.9g}')
        elif args.command=='stamp':
            old=args.fsf.read_text();updated=stamp(old,args.tr)
            tmp=args.fsf.with_name(args.fsf.name+'.timing.tmp');tmp.write_text(updated);tmp.replace(args.fsf)
            args.fsf.with_suffix('.timing.json').write_text(json.dumps(dict(tr_seconds=args.tr,
                source='Validated raw echo JSON/NIfTI and pre-FEAT BOLD timing; event files in seconds unchanged.'),indent=2)+'\n')
        elif not np.isclose(fsf_tr(args.fsf.read_text()),args.tr,rtol=0,atol=1e-5):
            raise ValueError('Existing FEAT uses a different TR; preserve it and perform a versioned timing repair: '+str(args.fsf))
    except (OSError,ValueError,KeyError) as exc:
        print('Timing validation failed:',exc,file=sys.stderr);raise SystemExit(1)


if __name__=='__main__':main()
