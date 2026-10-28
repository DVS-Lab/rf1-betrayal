#!/usr/bin/env python3
"""Save a compact read-only snapshot of the live stock-filter rerun."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import subprocess


ROOT = Path(__file__).resolve().parents[2]


def descendants(processes, seeds):
    selected = set(seeds)
    while True:
        added = {p['pid'] for p in processes if p['ppid'] in selected} - selected
        if not added:
            return [p for p in processes if p['pid'] in selected]
        selected.update(added)


def inventory(record):
    state = json.loads((record/'state.json').read_text())
    root = Path(state['output_root']) if state.get('output_root') else None
    cases = sorted(root.glob('sub-*')) if root else []
    counts = Counter(selected=state.get('selected_runs',108), staged=len(cases))
    stages = []
    for case in cases:
        log = case/(case.name+'.fmriprep.log')
        text = log.read_text(errors='replace') if log.is_file() else ''
        done = 'fMRIPrep finished successfully!' in text
        feats = list((case/'feat_reconstruction').glob('*.feat'))
        residual = any((p/'stats/res4d.nii.gz').is_file() for p in feats)
        fsfs = list((case/'feat_reconstruction').glob('*.fsf'))
        counts['fmriprep_finished'] += done
        counts['feat_started'] += bool(feats or fsfs)
        counts['residual_present'] += residual
        stages.append(dict(subject=case.name,fmriprep_finished=done,
                           feat_started=bool(feats or fsfs),residual_present=residual))
    finished = state.get('subjects',{})
    counts['completed_qc'] = sum(p.get('state')=='complete_qc_review_pending' for p in finished.values())
    counts['failed'] = sum(p.get('state')=='failed' for p in finished.values())
    counts['qc_reports'] = len(list((record/'qc').glob('*.json')))
    errors = [dict(subject=s,error=p.get('error','')) for s,p in finished.items() if p.get('state')=='failed']
    if state.get('error'):
        errors.append(dict(subject='runner',error=state['error']))
    return state,root,dict(counts),stages,errors


def process_table():
    text = subprocess.check_output(['ps','-u',str(os.getuid()),'-o',
        'pid=,ppid=,pcpu=,pmem=,etime=,stat=,args='],text=True)
    processes=[]
    for line in text.splitlines():
        fields=line.split(None,6)
        if len(fields)!=7:
            continue
        pid,ppid,cpu,mem,elapsed,status,args=fields
        processes.append(dict(pid=int(pid),ppid=int(ppid),cpu=float(cpu),
                              mem_percent=float(mem),elapsed=elapsed,status=status,args=args))
    seeds={p['pid'] for p in processes if p['pid']!=os.getpid() and
           re.search(r'(?:^|\s)(?:\S*/)?rerun_trust_references\.py(?:\s|$)',p['args'])}
    return descendants(processes,seeds),sorted(seeds)


def resources():
    memory={}
    path=Path('/proc/meminfo')
    if path.is_file():
        for line in path.read_text().splitlines():
            key,rest=line.split(':',1)
            if key in ('MemTotal','MemAvailable','SwapTotal','SwapFree'):
                memory[key+'_GiB']=round(int(rest.split()[0])/1024**2,1)
    affinity=len(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else os.cpu_count()
    try:
        sample=subprocess.check_output(['vmstat','1','2'],text=True,timeout=8).strip()
    except (OSError,subprocess.SubprocessError):
        sample='vmstat unavailable'
    return dict(logical_cpus_available=affinity,load_average=list(os.getloadavg()),memory=memory,vmstat=sample)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record',type=Path)
    args=parser.parse_args()
    candidates=sorted(p for p in (ROOT/'code/insula_ica/records').glob('linux2-trust-stock-filter-*')
                      if (p/'state.json').is_file())
    record=args.record or (candidates[-1] if candidates else None)
    if record is None:
        raise SystemExit('No stock-filter rerun state found.')
    state,root,counts,stages,errors=inventory(record)
    processes,runners=process_table()
    machine=resources()
    active=sorted(processes,key=lambda p:p['cpu'],reverse=True)
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    saved=ROOT/'code/insula_ica/output'/('trust-stock-filter-status-'+stamp)
    saved.mkdir(parents=True,exist_ok=False)
    snapshot=dict(record=str(record),output_root=str(root),state=state.get('state'),counts=counts,
                  stages=stages,errors=errors,runner_pids=runners,processes=processes,resources=machine)
    (saved/'snapshot.json').write_text(json.dumps(snapshot,indent=2)+'\n')
    lines=[f'Record: {record.name}',f'State: {state.get("state")}',
           'Runner PID(s): '+(', '.join(map(str,runners)) or 'none'),
           f'Staged subjects: {counts["staged"]}/{counts["selected"]}',
           f'fMRIPrep successful: {counts.get("fmriprep_finished",0)}/{counts["selected"]}',
           f'FEAT started: {counts.get("feat_started",0)}; residuals present: {counts.get("residual_present",0)}',
           f'Finished FEAT + QC: {counts["completed_qc"]}; QC reports: {counts["qc_reports"]}; failed: {counts["failed"]}',
           'Residual presence alone does not establish FEAT/QC completion.',
           f'CPUs available: {machine["logical_cpus_available"]}; load averages: '+', '.join(f'{x:.1f}' for x in machine['load_average']),
           'Memory GiB: '+json.dumps(machine['memory']),
           '\nCPU/I/O sample (last row is the current one-second interval):',machine['vmstat'],
           '\nBusiest pipeline processes (up to 12):']
    for p in active[:12]:
        command=p['args'].replace(str(root),'$OUTPUT')
        lines.append(f'PID {p["pid"]} PPID {p["ppid"]} CPU {p["cpu"]}% MEM {p["mem_percent"]}% elapsed {p["elapsed"]} {p["status"]} | {command[:210]}')
    if not active:
        lines.append('No processes descended from the stock-filter runner found.')
    lines.append('\nRecorded errors (up to 8):')
    lines.extend(f'{p["subject"]}: {p["error"]}' for p in errors[:8])
    if not errors:
        lines.append('None recorded.')
    waiting=[s for s in stages if s['fmriprep_finished'] and not s['feat_started']]
    lines.append('Successful fMRIPrep awaiting FEAT: '+str(len(waiting)))
    lines.append('Full snapshot saved: '+str(saved/'snapshot.json'))
    report='\n'.join(lines)+'\n'
    (saved/'report.txt').write_text(report)
    print(report,end='')


if __name__=='__main__':
    main()
