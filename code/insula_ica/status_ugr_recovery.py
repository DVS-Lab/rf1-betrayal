#!/usr/bin/env python3
"""Compact, read-only UGR anatomy reuse and sm-0 FEAT recovery progress."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import status_insula_streams as status


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--brief',action='store_true')
    parser.parse_args()
    records=status.ROOT/'code/insula_ica/records'
    repairs=[p for p in records.glob('linux2-ugr-TR-repair-*') if (p/'state.json').is_file()]
    if repairs:
        from status_ugr_timing_repair import main as repair_status
        repair_status();return
    choices=sorted(p for p in records.glob('linux2-ugr-anatomy-reuse-*') if (p/'state.json').is_file())
    if not choices:
        print('No UGR recovery record found.');return
    record=choices[-1];state=status.read(record/'state.json')
    rows=status.processes();runners,active,counts=status.pipeline(rows or [],('recover_ugr_unsmoothed.py',),Path(state['output_root']))
    errors=[]
    for group in ('subjects','runs','feat_runs'):
        for key,item in state.get(group,{}).items():
            if item.get('error'):errors.append(f'{group}/{key}: {item["error"]}')
    if state.get('error'):errors.append(state['error'])
    for path in (record/'publication-error.txt',):
        if path.is_file():errors.append('Git publication: '+path.read_text())
    audit=record/'raw-warpkit-audit.json'
    if audit.is_file():
        errors.extend(k+': '+v['error'] for k,v in status.read(audit).items() if v.get('error'))
    def done(group):return sum(v.get('state')=='complete' for v in state.get(group,{}).values())
    log=status.launcher('ugr-anatomy-reuse','rf1-ugr-anatomy-reuse.*');tail=status.log_errors(log)
    print('Record:',record.name)
    print('Runner:','unknown (process inspection unavailable)' if rows is None else ', '.join(f'PID {p["pid"]}, elapsed {p["elapsed"]}' for p in runners) if runners else 'not running')
    print('State:',state['state'],'| pilot passed:',state.get('pilot_passed',False))
    print(f'UGR subjects complete: {done("subjects")}/132; fMRIPrep + BOLD QC: {done("runs")}/249; sm-0 FEAT + QC: {done("feat_runs")}/249')
    print('Queue slots: fMRIPrep',state.get('fmriprep_worker_slots','pending'),'; FEAT',state.get('feat_worker_slots','pending'),'; usable CPUs',state.get('usable_cpus','pending'))
    if rows is not None:print(f'Active FEAT {counts["feat"]}; FILM {counts["film_gls"]}; MELODIC {counts["melodic"]}')
    combined=Path(state['combined_sm0_record']) if state.get('combined_sm0_record') else None
    if combined and (combined/'state.json').is_file():
        full=status.read(combined/'state.json');print('Combined Trust + UGR sm-0 FEAT:',sum(v.get('state')=='complete' for v in full['runs'].values()),'/ 490; state:',full['state'])
        if (combined/'publication-error.txt').is_file():errors.append('Combined FEAT Git publication: '+(combined/'publication-error.txt').read_text())
    print(f'Recorded errors: {len(errors)} (up to 6 shown)')
    for error in errors[:6]:print(' ',status.short(error))
    print('Launcher error tail:',' | '.join(tail) if tail else 'none found' if log else 'log unavailable')
    saved=status.ROOT/'code/insula_ica/output'/('ugr-recovery-status-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
    saved.mkdir(parents=True)
    (saved/'snapshot.json').write_text(json.dumps(dict(record=str(record),state=state,errors=errors,runner_processes=runners,
        active_processes=active,launcher_log=str(log) if log else None,launcher_errors=tail),indent=2)+'\n')
    print('Full snapshot:',saved/'snapshot.json')


if __name__=='__main__':main()
