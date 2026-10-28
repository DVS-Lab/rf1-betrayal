#!/usr/bin/env python3
"""Publish small progress/reference diagnostics without changing live jobs."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import shutil
import uuid

from run_bold_controls import publish, save
import status_trust_stock_filter as status


PATTERN = re.compile(r'setting before|single.band|sbref|coreg_boldref|finished successfully|'
                     r'ANAT Stage 4|traceback|exception|error',re.I)


def log_evidence(path):
    if not path.is_file():
        return dict(path=str(path),missing=True),''
    lines=path.read_text(errors='replace').splitlines()
    indices=[i for i,line in enumerate(lines) if PATTERN.search(line)]
    selected=sorted(set(indices[:50]+indices[-70:]+list(range(max(0,len(lines)-12),len(lines)))))
    return dict(path=str(path),total_lines=len(lines),matching_lines=len(indices),
                retained_lines=len(selected)), '\n'.join(f'{i+1}: {lines[i]}' for i in selected)+'\n'


def collect(record,destination,snapshot=None):
    state,root,counts,stages,errors=status.inventory(record)
    save(destination/'current-state.json',state)
    if snapshot is not None:
        shutil.copyfile(snapshot,destination/'resource-snapshot.json')
    subjects=[e['subject'] for e in errors if e['subject']!='runner'][:3]
    evidence=[]
    for subject in subjects:
        case=root/subject
        folder=destination/subject;folder.mkdir()
        info,text=log_evidence(case/(subject+'.fmriprep.log'))
        (folder/'reference-log-excerpt.txt').write_text(text)
        plan=record/'cases'/(subject+'.json')
        if plan.is_file():
            shutil.copyfile(plan,folder/'case.json')
        config=case/'fmriprep_config.json'
        if config.is_file():
            shutil.copyfile(config,folder/config.name)
        sidecars=[]
        func=case/'fmriprep-24'/subject/'ses-01/func'
        for source in sorted(func.glob('*boldref.json')):
            if source.stat().st_size<=256*1024:
                shutil.copyfile(source,folder/source.name);sidecars.append(source.name)
        view=case/'bids-view'/subject/'ses-01/func'
        sources=[dict(view=str(p),target=str(p.resolve())) for p in sorted(view.glob('*sbref.nii*'))]
        save(folder/'reference-files.json',sources)
        evidence.append(dict(subject=subject,log=info,sidecars=sidecars,sbrefs=len(sources)))
    save(destination/'diagnostic.json',dict(source_record=str(record),output_root=str(root),
        counts=counts,stages=stages,errors=errors,reference_evidence=evidence,
        analysis_jobs_changed=False,analysis_rerun=False))
    report=(f'# Live Trust rerun diagnostic\n\nSource: {record.name}\n\n'
            f'fMRIPrep successful: {counts.get("fmriprep_finished",0)}/{counts["selected"]}; '
            f'FEAT started: {counts.get("feat_started",0)}; residuals present: {counts.get("residual_present",0)}; '
            f'FEAT/QC complete: {counts["completed_qc"]}; reported failures: {counts["failed"]}.\n\n'
            'Reference log excerpts, unmodified reference sidecars and raw-reference symlink targets are saved '
            'for up to three reported failures. The prior resource snapshot is included when supplied. '
            'Unknown reference selection remains unverified. No jobs were signalled, launched or changed.\n')
    (destination/'README.md').write_text(report)
    return counts


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot',type=Path)
    parser.add_argument('--publish',action='store_true')
    args=parser.parse_args()
    snapshots=sorted((status.ROOT/'code/insula_ica/output').glob('trust-stock-filter-status-*/snapshot.json'))
    snapshot=args.snapshot or (snapshots[-1] if snapshots else None)
    if snapshot is None:
        raise SystemExit('Run status_trust_stock_filter.py first.')
    saved=json.loads(snapshot.read_text());record=Path(saved['record'])
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:6]
    destination=status.ROOT/'code/insula_ica/records'/('linux2-trust-progress-'+stamp)
    destination.mkdir()
    counts=collect(record,destination,snapshot)
    machine=status.resources()
    processes,runners=status.process_table()
    save(destination/'current-resources.json',dict(resources=machine,processes=processes,runner_pids=runners))
    print('Current stages:',json.dumps(counts),flush=True)
    print('Current resources:',json.dumps({k:v for k,v in machine.items() if k!='vmstat'}),flush=True)
    if args.publish:
        publish(destination,'Record live Trust progress and unresolved reference-selection evidence')
    print('Diagnostic record:',destination,flush=True)


if __name__=='__main__':
    main()
