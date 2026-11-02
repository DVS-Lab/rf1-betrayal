#!/usr/bin/env python3
"""Compact whole-brain dual-regression progress without participant image reads."""
import argparse
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'code/insula_ica'))
import status_insula_streams as s


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--brief',action='store_true');p.parse_args()
    choices=sorted(p for p in (ROOT/'code/functional_cluster_ica/records').glob('linux2-functional-wholebrain-*') if (p/'state.json').is_file())
    if not choices:print('No whole-brain functional dual-regression record yet.');return
    record=choices[-1];state=s.read(record/'state.json');processes=s.processes()
    runners,_,_=s.pipeline(processes or [],('code/functional_cluster_ica/wholebrain.py',),Path(state['output_root']))
    print('Record:',record.name)
    print('Runner:','unknown' if processes is None else ', '.join(f'PID {r["pid"]}, elapsed {r["elapsed"]}' for r in runners) if runners else 'not running')
    print('State:',state['state'],'| pilots passed:',state.get('pilot_passed',False))
    print('Slots: dual regression',state.get('jobs'),'| randomise',state.get('randomise_jobs'))
    for key,title,total in [('runs','Run-arm dual regressions',980),('groups','Participant-aggregated solutions',8),('randomise','Component randomise jobs',32)]:
        values=state.get(key,{})
        print(title+':',sum(v['state']=='complete' for v in values.values()),'/',total,'; failed:',sum(v['state']=='failed' for v in values.values()),'; reused:',sum(v.get('reused',False) for v in values.values()))
    print('Anatomy + social task plots:','ready' if state.get('anatomy_ready') else 'pending')
    active=[r for r in processes or [] if str(Path(state['output_root'])/'randomise') in r['args'] and 'randomise' in Path(r['args'].split()[0]).name]
    print('Active randomise processes:',len(active) if processes is not None else 'unknown')
    errors=[]
    for section in ('runs','randomise'):errors.extend(k+': '+v.get('error','failed') for k,v in state.get(section,{}).items() if v['state']=='failed')
    if state.get('error'):errors.append(state['error'])
    if (record/'publication-error.txt').is_file():errors.append('Git publication: '+(record/'publication-error.txt').read_text())
    print('Recorded errors:',len(errors),'(up to 6 shown)')
    for error in errors[:6]:print(' ',s.short(error))
    pointer=ROOT/'code/insula_ica/output/functional-wholebrain-launcher.logpath';log=Path(pointer.read_text().strip()) if pointer.is_file() else None
    print('Launcher errors:',' | '.join(s.log_errors(log)) or 'none found' if log else 'log unavailable')
    print('Results:',record)


if __name__=='__main__':main()
