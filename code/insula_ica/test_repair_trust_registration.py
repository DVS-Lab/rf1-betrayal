"""Guards for the isolated, costly single-participant correction candidate."""
from contextlib import ExitStack
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import repair_trust_registration as r
from diagnose_trust_registration import selected_report, sbref_selection, log_excerpt


class RegistrationRepairTests(unittest.TestCase):
    def test_command_changes_only_reference_selection(self):
        args = ('apptainer', Path('/images/24.sif'), Path('/upstream'), Path('/isolated'), Path('/scratch/fresh'))
        before = r.production_command(*args, r.SUBJECT)
        after = r.candidate_command(*args)
        self.assertEqual(after, before + ['--ignore', 'sbref'])
        self.assertEqual(after[after.index('--participant-label')+1], '10834')
        self.assertIn('/upstream/bids:/bids:ro', after)
        self.assertEqual(after[after.index('--task-id')+1], 'trust')

    def test_rows_restrict_both_runs_without_mutating_old_inputs_or_confounds(self):
        rows = [dict(subject='sub-10834', run=i, new_bold=f'/old/func/run-{i}.nii.gz', confounds=f'/saved/run-{i}.tsv') for i in (1,2)]
        rows += [dict(subject='sub-10402', run=1, new_bold='/old/other.nii.gz')]
        before = json.dumps(rows)
        selected = r.candidate_rows(rows, Path('/fresh'))
        self.assertEqual(json.dumps(rows), before)
        self.assertEqual({x['subject'] for x in selected}, {'sub-10834'})
        self.assertEqual([x['confounds'] for x in selected], ['/saved/run-1.tsv','/saved/run-2.tsv'])
        self.assertTrue(all(x['new_bold'].startswith('/fresh/fmriprep-24/sub-10834/') for x in selected))
        with self.assertRaises(ValueError): r.candidate_rows(rows[:1], Path('/fresh'))

    def test_real_report_names_are_selected_and_generic_warnings_do_not_hide_registration(self):
        for name in ('sub-10834_ses-01_space-MNI152NLin6Asym_T1w.svg',
                     'sub-10834_ses-01_dseg.svg', 'sub-10834_desc-rois_bold.svg',
                     'sub-10834_desc-coreg_bold.svg'):
            self.assertTrue(selected_report(Path(name)), name)
        self.assertFalse(selected_report(Path('sub-10834_desc-carpetplot_bold.svg')))
        with tempfile.TemporaryDirectory() as folder:
            log = Path(folder)/'log.txt'
            log.write_text('Using single-band reference file(s) run-1_sbref.nii.gz.\nBBR fallback rejected\n' +
                           'WARNING getdefaultlocale deprecated\n'*500)
            text, counts = log_excerpt(log)
            self.assertIn('run-1_sbref.nii.gz', text)
            self.assertIn('BBR fallback rejected', text)
            self.assertEqual(counts['matching_lines'], 2)

    def test_successful_preprocessing_reaches_fit_qc_and_retains_original_files(self):
        with tempfile.TemporaryDirectory() as folder, ExitStack() as patches:
            root = Path(folder); source = root/'code/insula_ica/records'/r.SOURCE
            (source/'qc').mkdir(parents=True)
            rows = []
            coverage = {roi:dict(retained_fraction=1.) for roi in ('left','right','bilateral')}
            for run in (1,2):
                bold = root/f'original/run-{run}.nii.gz'; bold.parent.mkdir(exist_ok=True); bold.write_bytes(b'unchanged old BOLD')
                conf = root/f'original/run-{run}.tsv'; conf.write_text('1\n2\n')
                row = dict(subject=r.SUBJECT, run=run, key=f'{r.SUBJECT}_run-{run}', new_bold=str(bold),
                           confounds=str(conf), confounds_sha256=r.core.sha(conf), warpkit=dict(ready=True))
                rows.append(row)
                (source/'qc'/(row['key']+'.json')).write_text(json.dumps(dict(coverage=coverage,historical_frame=dict(spatial_r=.3))))
            original_bytes = json.dumps(rows)
            (source/'input-registry.json').write_text(original_bytes)
            (source/'container.json').write_text(json.dumps(dict(sha256='same-image')))
            record = root/'record'; record.mkdir()
            patches.enter_context(patch.object(r.core, 'ROOT', root))
            patches.enter_context(patch.object(r, 'WORK_BASE', root/'scratch'))
            patches.enter_context(patch.dict(r.os.environ, FSLDIR='/fixture/fsl'))
            patches.enter_context(patch.object(r.shutil, 'which', side_effect=lambda exe: '/usr/bin/'+exe))
            patches.enter_context(patch.object(r.subprocess, 'check_output', return_value='fMRIPrep 24.1.1'))
            patches.enter_context(patch.object(r.subprocess, 'run'))
            patches.enter_context(patch.object(r.rebuild, 'record_container',
                side_effect=lambda out,*args: (out/'container.json').write_text(json.dumps(dict(sha256='same-image')))))
            patches.enter_context(patch.object(r.rebuild, 'verify_outputs'))
            fit = patches.enter_context(patch.object(r.rebuild, 'fit_run'))
            patches.enter_context(patch.object(r, 'collect', return_value=(4,0)))
            patches.enter_context(patch.object(r, 'sbref_selection', return_value=[]))
            patches.enter_context(patch.object(r, 'qc_run', side_effect=lambda row,out:
                dict(coverage=coverage, historical_frame=dict(spatial_r=.99),confounds=dict(exact_numeric_equality=True))))
            patches.enter_context(patch.object(r, 'montage', side_effect=lambda qc,dst,fsldir: dst.write_bytes(b'png fixture')))
            state = {}; r.run(record, state, 'fixture')
            self.assertEqual(state['state'], 'candidate_ready_for_visual_review')
            self.assertFalse(state['automatic_adoption'])
            self.assertEqual(fit.call_count, 2)
            self.assertTrue(all(c.kwargs == dict(strict_historical=False) for c in fit.call_args_list))
            self.assertEqual((source/'input-registry.json').read_text(), original_bytes)
            self.assertTrue(all(Path(row['new_bold']).read_bytes()==b'unchanged old BOLD' for row in rows))
            self.assertEqual(len(list((record/'qc').glob('*.png'))), 2)
            comparisons=json.loads((record/'comparison.json').read_text())
            self.assertTrue(all(c['confounds_exact'] for c in comparisons))

    def test_selected_reference_audit_distinguishes_run_mismatch_and_missing_evidence(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); subject='sub-10834'
            keys=[f'{subject}_ses-01_task-trust_run-{i}' for i in (1,2)]
            rows=[dict(key=k,subject=subject,new_bold=str(root/'fmriprep-24'/subject/'ses-01/func'/(k+'_bold.nii.gz'))) for k in keys]
            (root/(subject+'.fmriprep.log')).write_text('\n'.join(
                f'Report - setting before (/rebuild/fmriprep-24/{subject}/ses-01/func/{k}_part-mag_desc-coreg_boldref.nii.gz) '
                f'and after (/bids/{subject}/ses-01/func/{keys[0]}_echo-1_part-mag_sbref.nii.gz) images' for k in keys))
            rows.append(dict(key='sub-10402_ses-01_task-trust_run-1',subject='sub-10402',
                             new_bold=str(root/'fmriprep-24/sub-10402/ses-01/func/bold.nii.gz')))
            audit=sbref_selection(rows)
            self.assertEqual([r['status'] for r in audit], ['same_run','different_run_reference','unknown_no_retained_log_evidence'])
            self.assertEqual(audit[1]['sbref_run'],'1')


if __name__ == '__main__':
    unittest.main()
