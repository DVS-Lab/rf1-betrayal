#!/usr/bin/env python3
"""Compact current completed-ICA exploration status and error summary."""
import argparse
from pathlib import Path
import status_insula_streams as s


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--brief',action='store_true');p.parse_args()
    choices=sorted(p for p in (s.ROOT/'code/insula_ica/records').glob('linux2-completed-ICA-exploration-*') if (p/'state.json').is_file())
    if not choices:print('No completed-ICA exploration record found.');return
    record=choices[-1];state=s.read(record/'state.json');rows=s.processes()
    runners,_,_=s.pipeline(rows or [],('explore_completed_ica.py',),Path(state['output_root']))
    print('Record:',record.name)
    print('Runner:','unknown' if rows is None else ', '.join(f'PID {r["pid"]}, elapsed {r["elapsed"]}' for r in runners) if runners else 'not running')
    print('State:',state['state'])
    print('Connectivity runs completed:',sum(v['state']=='complete' for v in state['runs'].values()),'/ 490')
    print('Figures exported:',len(list((record/'figures').glob('*.png'))))
    if state.get('tests'):print('Task projections:',state['tests'],'; Bonferroni-60 p<.05:',state.get('task_tests_p_bonferroni_60_below_05'))
    errors=[k+': '+v.get('error','failed') for k,v in state['runs'].items() if v['state']=='failed']
    if state.get('error'):errors.append(state['error'])
    errors.extend('Git publication: '+p.read_text() for p in record.glob('*publication-error.txt'))
    print('Recorded errors:',len(errors),'(up to 6 shown)')
    for e in errors[:6]:print(' ',s.short(e))
    log=s.launcher('completed-ICA-exploration','rf1-ICA-exploration.*');tail=s.log_errors(log)
    print('Launcher errors:',' | '.join(tail) if tail else 'none found' if log else 'log unavailable')
    print('Review:',record/'README.md')


if __name__=='__main__':main()
