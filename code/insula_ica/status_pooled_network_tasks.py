#!/usr/bin/env python3
"""Brief pooled-model status; no imaging or large matrix reads."""
import argparse
from pathlib import Path
import status_insula_streams as s


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--brief', action='store_true'); p.parse_args()
    records = sorted((s.ROOT/'code/insula_ica/records').glob('linux2-pooled-network-tasks-*'))
    if not records: print('No pooled UGR task-model record found.'); return
    record = records[-1]; state = s.read(record/'state.json'); rows = s.processes()
    runners, _, _ = s.pipeline(rows or [], ('pooled_network_tasks.py',), Path(state['output_root']))
    print('Record:', record.name)
    print('Runner:', 'unknown' if rows is None else ', '.join(f'PID {p["pid"]}, elapsed {p["elapsed"]}' for p in runners) if runners else 'not running')
    print('State:', state['state'])
    print('Pooled UGR runs:', sum(r['state'] == 'complete' for r in state['runs'].values()), '/ 249; Trust runs reused:', state['reused_Trust_runs'])
    print('Group bar estimates:', state.get('bar_estimates', 0), '/ 150; figures:', len(list((record/'figures').glob('*.png'))), '/ 8')
    errors = [key+': '+item.get('error', 'failed') for key, item in state['runs'].items() if item['state'] == 'failed']
    if state.get('error'): errors.append(state['error'])
    errors.extend('Git publication: '+p.read_text() for p in record.glob('*publication-error.txt'))
    print('Recorded errors:', len(errors), '(up to 5 shown)')
    for error in errors[:5]: print(' ', s.short(error))
    print('Results:', record)


if __name__ == '__main__': main()
