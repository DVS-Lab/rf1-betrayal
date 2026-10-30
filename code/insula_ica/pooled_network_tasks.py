#!/usr/bin/env python3
"""Pool high/low UGR pmod columns before standardization; reuse cached timecourses."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import fcntl
import json
from pathlib import Path
import threading
import traceback
import uuid

import numpy as np

import run_network_timecourse_analysis as n
from run_bold_controls import save, publish

POLICY = 'UGR-common-high-low-slopes-pooled-before-zscore-v1'


def pool_design(design, settings):
    titles = ('nonsocial_high_constant', 'nonsocial_high_pmod', 'nonsocial_low_constant',
              'nonsocial_low_pmod', 'social_high_constant', 'social_high_pmod',
              'social_low_constant', 'social_low_pmod')
    for k, title in enumerate(titles, 1):
        n.g.core.require(settings[f'fmri(evtitle{k})'] == title and settings.get(f'fmri(deriv_yn{k})', '0') == '0',
                         'Unexpected UGR EV order/derivatives')
    columns = []; mapping = []
    for old in range(design.shape[1]):
        if old in (3, 7): continue
        source = [old, old+2] if old in (1, 5) else [old]
        columns.append(design[:, source].sum(axis=1))
        mapping.append(dict(new_column=len(columns)-1, original_columns=source,
                            operation='sum before z-scoring' if len(source) == 2 else 'unchanged'))
    pooled = np.column_stack(columns)
    # Columns 1 and 4 carry the pooled nonsocial and social pmods. The four
    # high/low main effects and all RT/missed/nuisance columns remain separate.
    cs = np.zeros((2, pooled.shape[1])); cs[0, 1] = 1; cs[1, 4] = 1
    return pooled, cs, mapping


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--source-record', type=Path)
    p.add_argument('--jobs', type=int, default=16); p.add_argument('--publish', action='store_true'); args = p.parse_args()
    n.g.core.require(1 <= args.jobs <= 32, 'Use 1-32 small cached-model workers')
    source = args.source_record or sorted(n.RECORDS.glob('linux2-network-timecourses-*'))[-1]
    state = json.loads((source/'state.json').read_text()); root = Path(state['output_root'])
    n.g.core.require(state['policy'] == n.POLICY and len(state['runs']) == 490 and
                     all(r['state'] == 'complete' for r in state['runs'].values()), 'Expected complete standardized task source')
    registry = json.loads((source/'accepted-inputs.json').read_text()); ugr = [r for r in registry if r['task'] == 'ugr']
    n.g.core.require(len(ugr) == 249 and len({r['subject'] for r in ugr}) == 132, 'UGR cohort changed')
    lockdir = n.ROOT/'code/insula_ica/output'; lockdir.mkdir(exist_ok=True)
    with (lockdir/'network-pooled-model.lock').open('a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:6]
        record = n.RECORDS/('linux2-pooled-network-tasks-'+stamp); record.mkdir(); (record/'errors').mkdir(); (record/'figures').mkdir()
        output = n.ROOT/'derivatives/insula_ica'/('pooled-network-tasks-'+stamp); output.mkdir(parents=True)
        progress = dict(state='fitting_pooled_UGR_models', policy=POLICY, source_record=str(source),
                        output_root=str(output), runs={}, jobs=args.jobs, reused_Trust_runs=241,
                        no_4D_reads=True, no_connectivity_or_randomise=True)
        save(record/'state.json', progress); lock = threading.Lock()

        def worker(row):
            try:
                folder = root/row['key']; marker = json.loads((folder/'complete.json').read_text())
                expected = {v['path']:v for v in marker['outputs']}
                feat = Path(row['active_feat']); model_path = feat/'design.mat'
                source_models = {v['path']:v for v in marker['signature']['inputs']}
                n.g.core.require(n.g.core.identity(model_path) == source_models[str(model_path.resolve())], 'Saved L1 model changed')
                settings = n.g.core.fsl_settings(feat/'design.fsf')
                design = n.g.core.fsl_matrix(model_path); pooled, cs, mapping = pool_design(design, settings)
                run_output = output/row['key']; run_output.mkdir(); rows = []; inputs = []
                for mode in n.g.inputs.MODES:
                    for roi in n.g.core.ROIS:
                        current = folder/f'{mode}_{roi}'
                        tc_path = current/'standardized_timecourses.npy'; meta_path = current/'standardization.json'
                        for path in (tc_path, meta_path):
                            identity = n.g.core.identity(path)
                            n.g.core.require(identity == expected[str(path.resolve())], 'Cached response/scale changed')
                            inputs.append(identity)
                        meta = json.loads(meta_path.read_text()); zy = np.load(tc_path)
                        # Reconstruct the original centered stage-1 response,
                        # then let the same fitter standardize data/model.
                        tc = zy*np.asarray(meta['response_SD'])+np.asarray(meta['response_mean'])
                        estimates, _ = n.fit_timecourses(tc, pooled, cs, run_output/f'{mode}_{roi}')
                        for j, label in enumerate(n.BAR_LABELS[:2]):
                            for k in range(5):
                                rows.append(dict(mode=mode, roi=roi, subject=row['subject'], task='ugr', run=row['run'],
                                    component=k+1, contrast=label, effect=float(estimates[j, k]), plotted=True,
                                    estimator='common high/low slope; pmods summed before standardization'))
                zx, _, _, active = n.zscore_columns(pooled, allow_constant=True)
                singular = np.linalg.svd(zx[:, active], compute_uv=False)
                condition = float(singular[0]/singular[-1])
                n.g.core.require(condition < 1e6 and all(abs(r['effect']) < 5 for r in rows),
                                 'Pooled model still unstable; no incomplete-cohort summary')
                receipt = dict(key=row['key'], condition=condition, original_model=n.g.core.identity(model_path),
                               mapping=mapping, cached_inputs=inputs, original_columns=design.shape[1], pooled_columns=pooled.shape[1])
                save(record/(row['key']+'.json'), receipt)
                result = dict(state='complete', condition=condition, OLS_fits=6)
            except Exception as exc:
                rows = None; result = dict(state='failed', error=str(exc)); (record/'errors'/(row['key']+'.txt')).write_text(traceback.format_exc())
            with lock:
                progress['runs'][row['key']] = result; save(record/'state.json', progress)
                print(f'Pooled UGR {len(progress["runs"])}/249 | {row["key"]}: {result["state"]}', flush=True)
            return rows

        try:
            with ThreadPoolExecutor(max_workers=args.jobs) as pool:
                results = [f.result() for f in as_completed([pool.submit(worker, row) for row in ugr])]
            n.g.core.require(all(r is not None for r in results), 'Some pooled UGR fits failed; original results preserved')
            effects = [dict(r, component=int(r['component']), run=int(r['run']), effect=float(r['effect']))
                       for r in n.g.core.read_table(source/'run_network_task_effects.tsv') if r['task'] == 'trust']
            effects.extend(r for result in results for r in result)
            participant, group = n.summarize(effects, sorted({r['subject'] for r in registry}))
            for name, values in [('run_network_task_effects.tsv', effects), ('participant_network_task_effects.tsv', participant),
                                 ('group_network_task_effects.tsv', group)]: n.g.core.write_table(record/name, values)
            n.figures(group, record/'figures', note='N=132; all 490 accepted runs retained; pointwise 95% CIs. UGR high/low pmods pooled BEFORE z-scoring; Trust unchanged.')
            progress.update(state='pooled_five_bar_review_ready', bar_estimates=150, participants=132,
                            accepted_runs=490, excluded_runs=0, excluded_participants=0)
            save(record/'state.json', progress)
            (record/'README.md').write_text('All five bars for all five ICs, left/right/bilateral, in both discovery bases. '
                'UGR pmods are summed across high/low within nonsocial/social BEFORE standardization and joint OLS. '
                'The four high/low main effects, RT, missed-trial and all nuisance columns remain. '
                'This assumes a common high/low fairness slope and differs from averaging separately standardized coefficients. '
                'All 249 UGR runs use the same pooled model; 241 Trust runs reuse their unchanged estimates. '
                'All 132 participants and 490 accepted runs are retained. No unidentifiable social-high coefficient is manufactured. '
                'Cached stage-1 timecourses are reused. No 4D BOLD extraction, original FEAT edit, ICA, connectivity or randomise. '
                'Original unpooled results remain available in the source record.\n')
        except Exception as exc:
            progress.update(state='pooled_model_review_failed', error=str(exc)); save(record/'state.json', progress)
            (record/'failure.txt').write_text(traceback.format_exc()); print('ERROR:', exc, flush=True)
        if args.publish:
            try: publish(record, 'Publish pooled-endowment UGR task model and full five-bar figures; retain all subjects and runs')
            except Exception as exc: (record/'publication-error.txt').write_text(str(exc)); print('Git publication failed:', exc, flush=True)
        print('Final state:', progress['state'], '| record:', record, flush=True)
        raise SystemExit(0 if progress['state'] == 'pooled_five_bar_review_ready' else 1)


if __name__ == '__main__': main()
