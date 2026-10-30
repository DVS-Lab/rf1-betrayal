#!/usr/bin/env python3
"""Bounded status for direct network timecourse task models."""
import argparse
from pathlib import Path
import status_insula_streams as s


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--brief', action='store_true'); parser.parse_args()
    choices = sorted((s.ROOT/'code/insula_ica/records').glob('linux2-network-timecourses-*'))
    if not choices: print('No network timecourse record found.'); return
    record = choices[-1]; state = s.read(record/'state.json'); rows = s.processes()
    runners, _, _ = s.pipeline(rows or [], ('run_network_timecourse_analysis.py',), Path(state['output_root']))
    print('Record:', record.name)
    print('Runner:', 'unknown' if rows is None else ', '.join(f'PID {r["pid"]}, elapsed {r["elapsed"]}' for r in runners) if runners else 'not running')
    print('State:', state['state'], '| pilot passed:', state['pilot_passed'])
    print('Runs complete:', sum(v['state'] == 'complete' for v in state['runs'].values()), '/ 490')
    print('Standardized OLS fits:', sum(v.get('network_OLS_fits', 0) for v in state['runs'].values()), '/ 2940')
    print('Group bar estimates:', state.get('bar_estimates', 0), '/ 150')
    print('Task-bar figures:', len(list((record/'figures').glob('*five_task_bars.png'))), '/ 8')
    print('Coefficient magnitude review flags:', state.get('coefficient_review_flag_count', 'not audited'),
          '| numerical review required:', state.get('numerical_review_required', 'not audited'))
    errors = [k+': '+v['error'] for k, v in state['runs'].items() if v['state'] == 'failed']
    if state.get('error'): errors.append(state['error'])
    errors.extend('Git publication: '+p.read_text() for p in record.glob('*publication-error.txt'))
    print('Recorded errors:', len(errors), '(up to 6 shown)')
    for error in errors[:6]: print(' ', s.short(error))
    log = s.launcher('network-timecourses', 'rf1-network-timecourses.*'); tail = s.log_errors(log)
    print('Launcher errors:', ' | '.join(tail) if tail else 'none found' if log else 'log unavailable')
    print('Results:', record)


if __name__ == '__main__': main()
