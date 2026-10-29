#!/usr/bin/env python3
"""Compact, read-only paired ICA progress and error report; no imaging reads."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import subprocess


ROOT=Path(__file__).resolve().parents[2]


def read_json(path):
    return json.loads(path.read_text())


def short(value):
    return ' '.join(str(value).split())[:320]


def process_table():
    try:output=subprocess.check_output(['ps','-u',str(os.getuid()),'-o','pid=,stat=,etime=,args='],text=True,stderr=subprocess.DEVNULL,timeout=10)
    except (OSError,subprocess.SubprocessError):return None
    rows=[]
    for line in output.splitlines():
        parts=line.split(None,3)
        if len(parts)==4 and not parts[1].startswith('Z'):
            rows.append(dict(pid=int(parts[0]),status=parts[1],elapsed=parts[2],args=parts[3]))
    return rows


def latest_log():
    pointer=ROOT/'code/insula_ica/output/paired-group-ica-launch.lock/log'
    if pointer.is_file():
        path=Path(pointer.read_text().strip())
        if path.is_file():return path
    logs=[p for p in Path('/tmp').glob('rf1-paired-group-ica.*') if p.is_file() and p.stat().st_uid==os.getuid()]
    return max(logs,key=lambda p:p.stat().st_mtime_ns) if logs else None


def log_errors(path):
    if path is None:return []
    with path.open('rb') as stream:
        stream.seek(max(0,path.stat().st_size-131072))
        text=stream.read().decode(errors='replace')
    return [short(line) for line in text.splitlines()
        if re.search(r'Traceback|\bERROR:|\b\w*Error:|\bException:|Git publication failed:',line)][-4:]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record',type=Path)
    parser.add_argument('--brief',action='store_true',help='Compact output (also the default)')
    args=parser.parse_args()
    candidates=sorted(p for p in (ROOT/'code/insula_ica/records').glob('linux2-paired-group-ica-*') if (p/'state.json').is_file())
    record=args.record or (candidates[-1] if candidates else None)
    state=read_json(record/'state.json') if record else {}
    root=Path(state['output_root']) if state.get('output_root') else None
    processes=process_table()
    runners=[p for p in processes or [] if Path(p['args'].split()[0]).name.startswith('python') and
        re.search(r'(?:^|\s)(?:\S*/)?(?:run|resume)_group_ica\.py(?:\s|$)',p['args'])]
    active={exe:[p for p in processes or [] if root and str(root) in p['args'] and
        Path(p['args'].split()[0]).name==exe] for exe in ('melodic','film_gls')}
    completed=[];errors=[]
    if root:
        selection=root/'selected_ica_fits.json'
        selected=read_json(selection).get('fits',{}) if selection.is_file() else {}
        for mode in ('joint_bold','native_film'):
            for roi in ('left','right','bilateral'):
                for dim in (3,5,7,0):
                    key=f'{mode}/{roi}/dim-{dim}'
                    chosen=selected.get(key) or state.get('fits',{}).get(key,{})
                    marker=Path(chosen['oriented_maps']).parent/'complete.json' if chosen.get('oriented_maps') else root/'ica'/key/'complete.json'
                    if not marker.is_file():continue
                    try:
                        item=read_json(marker)
                        if not Path(item['oriented_maps']).is_file():raise ValueError('completion marker has no saved map')
                        completed.append(dict(mode=mode,roi=roi,dimension=dim,actual_components=item['actual_components']))
                    except (OSError,ValueError,KeyError) as exc:errors.append(f'{mode}/{roi}/dim-{dim}: {exc}')
    for collection in ('runs','fits'):
        for key,item in state.get(collection,{}).items():
            if item.get('state')=='failed':errors.append(f'{key}: {item.get("error","failed")}')
    if state.get('error'):errors.append('Runner: '+state['error'])
    publication=record/'publication-error.txt' if record else None
    if publication and publication.is_file():errors.append('Git publication: '+publication.read_text())
    log=latest_log();launcher_errors=log_errors(log)
    prepared=sum(r.get('state')=='complete' for r in state.get('runs',{}).values())
    fixed=sum(c['dimension']!=0 for c in completed);automatic=len(completed)-fixed
    status='running: '+', '.join(f'PID {p["pid"]}, elapsed {p["elapsed"]}' for p in runners) if runners else 'not running'
    if processes is None:status='unknown (process inspection unavailable)'
    lines=[f'Record: {record.name if record else "none found"}',f'Runner: {status}',
        f'State: {state.get("state","not recorded")}',f'Paired residual inputs prepared: {prepared}/490',
        f'ICA fits completed: {len(completed)}/24 (fixed: {fixed}/18; automatic: {automatic}/6)',
        f'Active processes: MELODIC {len(active["melodic"])}; FILM {len(active["film_gls"])}' if processes is not None else 'Active process counts unavailable.']
    if automatic:
        lines.append('Automatic component counts: '+ '; '.join(f'{c["mode"]}/{c["roi"]}={c["actual_components"]}'
            for c in completed if c['dimension']==0))
    pilot=state.get('native_variance_normalization_pilot',{})
    if pilot:
        lines.append(f'Native MELODIC VN pilot completed: {sum(v.get("state")=="complete" for v in pilot.values())}/2 (separate from 24 baseline fits)')
        errors.extend(f'VN pilot {mode}: {v.get("error","failed")}' for mode,v in pilot.items() if v.get('state')=='failed')
    lines+=['Recorded errors (up to 8):']+([short(e) for e in errors[:8]] or ['None recorded.'])
    if len(errors)>8:lines.append(f'{len(errors)-8} more errors in the saved snapshot.')
    lines+=['Latest launcher: '+str(log or 'none found'),'Launcher error lines (up to 4):']+(launcher_errors or ['None found in the bounded log tail.'])
    if not runners and state.get('state')=='paired_ICA_complete_review_pending':lines.append('Finished; ICA results are ready for review.')
    elif processes is not None and not runners and state.get('state') not in ('failed','paired_ICA_complete_review_pending'):
        lines.append('No live runner and no completion state; inspect the reported record/log before restarting.')
    saved=ROOT/'code/insula_ica/output'/('paired-ica-status-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
    saved.mkdir(parents=True)
    snapshot=dict(record=str(record) if record else None,state=state,completed=completed,errors=errors,
        runner_processes=runners,active_pipeline_processes=active,launcher_log=str(log) if log else None,launcher_errors=launcher_errors)
    (saved/'snapshot.json').write_text(json.dumps(snapshot,indent=2)+'\n')
    lines.append('Full snapshot: '+str(saved/'snapshot.json'))
    print('\n'.join(lines))


if __name__=='__main__':main()
