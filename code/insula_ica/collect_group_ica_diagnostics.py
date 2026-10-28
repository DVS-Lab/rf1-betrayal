#!/usr/bin/env python3
"""Collect bounded MELODIC logs for review; never start or alter an analysis."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import subprocess
import sys
import uuid


ROOT=Path(__file__).resolve().parents[2]


def tail(path,limit=262144):
    size=path.stat().st_size
    with path.open('rb') as stream:
        stream.seek(max(0,size-limit))
        text=stream.read().decode(errors='replace')
    return text,dict(source=str(path),source_bytes=size,truncated_head=size>limit)


def collect(source,destination):
    state=json.loads((source/'state.json').read_text());root=Path(state['output_root'])
    summary=dict(source_record=str(source),output_root=str(root),state=state['state'],
        preparation=dict(Counter(r['state'] for r in state.get('runs',{}).values())),
        fits=dict(Counter(r['state'] for r in state.get('fits',{}).values())),fit_details=[])
    for mode in ('joint_bold','native_film'):
        for roi in ('left','right','bilateral'):
            for dim in (3,5,7,0):
                key=f'{mode}/{roi}/dim-{dim}';folder=root/'ica'/key
                item=dict(key=key,recorded=state.get('fits',{}).get(key),directory_exists=folder.is_dir(),
                    files=[dict(name=p.name,bytes=p.stat().st_size) for p in sorted(folder.iterdir()) if p.is_file()] if folder.is_dir() else [],logs=[])
                for name,path in [('launcher',folder.parent/f'dim-{dim}.log'),('melodic',folder/'log.txt'),
                                  ('eigenvalues',folder/'eigenvalues_percent')]:
                    if not path.is_file():continue
                    text,identity=tail(path)
                    saved=destination/(key.replace('/','_')+'_'+name+'.txt');saved.write_text(text)
                    item['logs'].append(dict(**identity,saved=saved.name))
                    if name in ('launcher','melodic'):
                        hits=[re.sub(r'\s+',' ',line).strip()[:300] for line in text.splitlines()
                            if re.search(r'converg|error|fail|kill|exception|abort',line,re.I)]
                        item[name+'_diagnostic_lines']=hits[-4:]
                marker=folder/'complete.json'
                if marker.is_file():item['completion']=json.loads(marker.read_text())
                summary['fit_details'].append(item)
    (destination/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    counts=summary['fits'];done=counts.get('complete',0);failed=counts.get('failed',0)
    (destination/'README.md').write_text('# Paired ICA diagnostic collection\n\n'
        f'Source: `{source.name}`. Recorded state: `{state["state"]}`. '
        f'{done}/24 fits completed, {failed} failed, {24-done-failed} have no recorded result. '
        'Bounded logs and file inventories only; no fit, preprocessing, input or completion marker was changed. '
        'See summary.json and the per-fit text logs.\n')
    print(f'Prepared runs: {summary["preparation"].get("complete",0)}/490 | fits complete: {done}/24 | failed: {failed}')
    for item in summary['fit_details']:
        if (item['recorded'] or {}).get('state')!='failed':continue
        hits=item.get('launcher_diagnostic_lines',[])+item.get('melodic_diagnostic_lines',[])
        print(item['key']+': '+(' | '.join(dict.fromkeys(hits[-4:])) if hits else 'No matching diagnostic line; logs saved.'))
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record',type=Path)
    parser.add_argument('--publish',action='store_true')
    args=parser.parse_args()
    records=ROOT/'code/insula_ica/records'
    candidates=sorted(p for p in records.glob('linux2-paired-group-ica-*') if (p/'state.json').is_file())
    source=args.record or (candidates[-1] if candidates else None)
    if source is None:raise SystemExit('No paired ICA record found.')
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:6]
    destination=records/('linux2-paired-ica-diagnostics-'+stamp);destination.mkdir()
    (destination/'code-commit.txt').write_text(subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True))
    collect(source,destination)
    if args.publish:
        logs=ROOT/'code/insula_ica/output';logs.mkdir(exist_ok=True)
        log=logs/('paired-ica-diagnostic-publication-'+stamp+'.txt')
        program='from pathlib import Path; import sys; from run_bold_controls import publish; publish(Path(sys.argv[1]), "Collect paired insula ICA failure logs without rerunning completed work")'
        with log.open('w') as stream:
            result=subprocess.run([sys.executable,'-c',program,str(destination.resolve())],cwd=Path(__file__).parent,
                stdout=stream,stderr=subprocess.STDOUT)
        if result.returncode:raise SystemExit('Publication failed; diagnostics saved. Git details: '+str(log))
        print('Diagnostic record pushed to insula-ica.')
    print('Saved:',destination)


if __name__=='__main__':main()
