#!/usr/bin/env python3
"""Bounded status of the authorized one-run UGR timing correction."""
import argparse
from pathlib import Path
import status_insula_streams as s


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--brief',action='store_true');p.parse_args()
    choices=sorted(r for r in (s.ROOT/'code/insula_ica/records').glob('linux2-ugr-TR-repair-*') if (r/'state.json').is_file())
    if not choices:print('No UGR timing repair record found.');return
    record=choices[-1];state=s.read(record/'state.json');rows=s.processes()
    runners,_,counts=s.pipeline(rows or [],('repair_ugr_timing.py',),Path(state['output_root']))
    print('Record:',record.name)
    print('Runner:','unknown' if rows is None else ', '.join(f'PID {r["pid"]}, elapsed {r["elapsed"]}' for r in runners) if runners else 'not running')
    print('State:',state['state'],'| corrected TR: 1.710 seconds')
    print('Completed repair/control fits:',sum(v.get('state')=='complete' for v in state.get('fits',{}).values()),'/ 3')
    print('Active FEAT:',counts['feat'],'| FILM:',counts['film_gls'])
    if state.get('combined_sm0_record'):
        combined=s.read(Path(state['combined_sm0_record'])/'state.json')
        print('Combined sm-0 FEAT:',sum(v['state']=='complete' for v in combined['runs'].values()),'/ 490;',combined['state'])
    trust=record/'sub-10817-Trust-timing.json'
    if trust.is_file():
        for item in s.read(trust):print(item['key'],': acquisition TR',item.get('tr_seconds'),'; model',item.get('modeled_tr_seconds'),';',item['status'])
    errors=[state['error']] if state.get('error') else []
    if (record/'publication-error.txt').is_file():errors.append((record/'publication-error.txt').read_text())
    print('Recorded errors:',len(errors))
    for error in errors:print(s.short(error))
    log=s.launcher('ugr-timing-repair','rf1-ugr-timing-repair.*')
    print('Launcher errors:',(' | '.join(s.log_errors(log)) or 'none found') if log else 'log unavailable')
    print('Report:',record/'report.md')


if __name__=='__main__':main()
