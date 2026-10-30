#!/usr/bin/env python3
"""Recheck archived physiological EVs without Linux2 or any BOLD reads."""
import argparse
from collections import Counter
import hashlib
import io
import json
from pathlib import Path
import tarfile

import numpy as np
import audit


def review(record):
    rows=json.loads((record/'physiological_EVs.json').read_text())
    manifest=json.loads((record/'source_manifest.json').read_text())
    by_path={r['source']:r for r in manifest if r.get('archive')}
    entries={}
    for row in rows:
        for entry in (row['cached_timecourse_source'],by_path[row['case']+'/design.mat']):
            entries[entry['archive'],entry['member']]=entry
    blobs={}
    for archive in sorted({k[0] for k in entries}):
        with tarfile.open(record/archive,'r|gz') as stream:
            for member in stream:
                key=archive,member.name
                if key not in entries:continue
                data=stream.extractfile(member).read()
                if hashlib.sha256(data).hexdigest()!=entries[key]['sha256']:
                    raise ValueError('Archived source hash mismatch: '+str(key))
                blobs[key]=data
    def read(entry):return blobs[entry['archive'],entry['member']]
    results=[]
    for row in rows:
        # Processing eligibility/column mapping was established from each saved
        # FSF in the initial audit; do not reinterpret unsupported processing.
        if row['model_check']['status'] not in ('matches_up_to_positive_centering_scaling','does_not_match_supported_model'):
            raise ValueError('Unsupported original model check: '+row['case'])
        x=np.loadtxt(io.BytesIO(read(row['cached_timecourse_source']))).ravel()
        entry=by_path[row['case']+'/design.mat']
        matrix=np.loadtxt(io.StringIO(read(entry).decode().split('/Matrix',1)[1]),ndmin=2)
        y=matrix[:,row['model_check']['real_column_zero_based']]
        results.append(dict(case=row['case'],run_key=row['run_key'],seed=row['folder_seed_label'],
            accepted=row['accepted_manuscript_run'],cached_TC_sha256=row['cached_timecourse_source']['sha256'],
            saved_design_sha256=entry['sha256'],raw_text=audit.affine_check(x,y),
            float32_EV_input=audit.affine_check(x.astype(np.float32).astype(float),y),
            model_agreement=audit.model_column_agreement(x,y)))
    destination=record/'review';destination.mkdir(exist_ok=True)
    (destination/'precision-comparison.json').write_text(json.dumps(results,indent=2)+'\n')
    for seed in sorted({r['seed'] for r in results}):
        subset=[r for r in results if r['seed']==seed and r['accepted']]
        print(seed, 'accepted',len(subset),'status',dict(Counter(r['model_agreement']['status'] for r in subset)))
    return results


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('record',type=Path)
    review(parser.parse_args().record)
