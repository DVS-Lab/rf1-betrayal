#!/usr/bin/env python3
"""Read-only bounded functional-cluster status without image loading."""
import argparse
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'code/insula_ica'))
import status_insula_streams as s


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--brief',action='store_true');p.add_argument('--record',type=Path);args=p.parse_args()
    choices=sorted(p for p in (ROOT/'code/functional_cluster_ica/records').glob('linux2-functional-cluster-*') if (p/'state.json').is_file())
    record=args.record or (choices[-1] if choices else None)
    if record is None:print('No functional-cluster run record yet.');return
    state=json.loads((record/'state.json').read_text());root=Path(state['output_root']);rows=s.processes()
    runners,_,counts=s.pipeline(rows or [],('code/functional_cluster_ica/run.py',),root)
    # Resolve by the actual runner script path, rather than every generic run.py.
    runners=[r for r in rows or [] if Path(r['args'].split()[0]).name.startswith('python') and 'code/functional_cluster_ica/run.py' in r['args']]
    print('Record:',record.name)
    print('Runner:','unknown' if rows is None else ', '.join(f'PID {r["pid"]}, elapsed {r["elapsed"]}' for r in runners) if runners else 'not running')
    print('State:',state['state'])
    for key,title in [('mean_runs','Functional mean task fits'),('residual_runs','Paired residual runs prepared'),('component_runs','Functional IC task fits')]:
        values=state.get(key,{})
        print(title+':',sum(v['state']=='complete' for v in values.values()),'/ 490; failed:',sum(v['state']=='failed' for v in values.values()))
    fits=state.get('fits',{});done=[v for v in fits.values() if v['state']=='complete']
    print('ICA complete:',len(done),'/ 8; fixed:',sum(v['dimensionality']!=0 for v in done),'/ 6; automatic:',sum(v['dimensionality']==0 for v in done),'/ 2')
    for v in done:
        if v['dimensionality']==0:print('Automatic components:',v['mode'],v['actual_components'])
    print('Mean plots:', 'ready' if state.get('mean_plots_ready') else 'pending','| existing anatomical estimates reused:',state.get('anatomical_estimates_reused',465 if state.get('mean_plots_ready') else 0))
    print('Active processes:',', '.join(k+' '+str(v) for k,v in counts.items()))
    errors=[]
    for key in ('mean_runs','residual_runs','fits','component_runs'):
        errors.extend(k+': '+v.get('error','failed') for k,v in state.get(key,{}).items() if v['state']=='failed')
    if state.get('error'):errors.append(state['error'])
    if (record/'publication-error.txt').is_file():errors.append('Git publication: '+(record/'publication-error.txt').read_text())
    print('Recorded errors:',len(errors),'(up to 6 shown)')
    for line in errors[:6]:print(' ',s.short(line))
    pointer=ROOT/'code/insula_ica/output/functional-cluster-launcher.logpath'
    log=Path(pointer.read_text().strip()) if pointer.is_file() else None
    print('Launcher errors:',' | '.join(s.log_errors(log)) or 'none found' if log else 'log unavailable')
    print('Results:',record)


if __name__=='__main__':main()
