import json
import io
import sys
from contextlib import redirect_stdout
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import status_trust_stock_filter as s


class StatusTests(unittest.TestCase):
    def test_residual_is_separate_from_completed_qc(self):
        with tempfile.TemporaryDirectory() as directory:
            base=Path(directory);record=base/'record';record.mkdir();root=base/'outputs'
            case=root/'sub-10001';case.mkdir(parents=True)
            (case/'sub-10001.fmriprep.log').write_text('fMRIPrep finished successfully!\n')
            stats=case/'feat_reconstruction/run.feat/stats';stats.mkdir(parents=True)
            (stats/'res4d.nii.gz').touch()
            (record/'state.json').write_text(json.dumps(dict(output_root=str(root),selected_runs=108,subjects={})))
            _,_,counts,_,errors=s.inventory(record)
            self.assertEqual(counts['fmriprep_finished'],1)
            self.assertEqual(counts['residual_present'],1)
            self.assertEqual(counts['completed_qc'],0)
            self.assertEqual(errors,[])

    def test_failed_and_pending_are_distinct(self):
        with tempfile.TemporaryDirectory() as directory:
            record=Path(directory)
            (record/'state.json').write_text(json.dumps(dict(subjects={
                'sub-1':dict(state='failed',error='example failure'),
                'sub-2':dict(state='complete_qc_review_pending')})))
            _,_,counts,_,errors=s.inventory(record)
            self.assertEqual(counts['completed_qc'],1)
            self.assertEqual(counts['failed'],1)
            self.assertEqual(errors[0]['error'],'example failure')

    def test_process_tree_includes_grandchildren_and_excludes_other_jobs(self):
        rows=[dict(pid=1,ppid=0),dict(pid=2,ppid=1),dict(pid=3,ppid=2),dict(pid=4,ppid=0)]
        self.assertEqual({p['pid'] for p in s.descendants(rows,{1})},{1,2,3})

    def test_inherited_feat_is_visible_but_not_reported_as_a_scheduler(self):
        rows=('1 0 0.0 0.1 00:01 S python code/insula_ica/resume_trust_stock_filter.py --record example\n'
              '2 999 90.0 0.1 00:02 R bash /fsl/feat /case/run.fsf\n')
        with patch.object(s.subprocess,'check_output',return_value=rows),patch.object(s.os,'getpid',return_value=100):
            processes,runners=s.process_table([2,999])
        self.assertEqual(runners,[1])
        self.assertEqual({p['pid'] for p in processes},{1,2})

    def test_status_distinguishes_running_64_from_success_and_stale_state(self):
        process=dict(pid=12,args='/env/bin/python -u code/insula_ica/resume_trust_stock_filter.py --jobs 64')
        self.assertIn('jobs=64',s.resume_status({},[process],[12]))
        self.assertIn('finished',s.resume_status(dict(state='functional_reruns_complete_qc_pending'),[],[]))
        self.assertIn('not running',s.resume_status(dict(state='running'),[],[]))

    def test_launcher_errors_select_latest_attempt_and_bound_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);old=root/'rf1-trust-stock-resume.old';new=root/'rf1-trust-stock-resume.new'
            old.write_text('ValueError: earlier attempt\n')
            new.write_text('ordinary progress\n'*1000+'ValueError: current attempt\n'*20)
            import os
            os.utime(old,(1,1));os.utime(new,(2,2))
            log,errors=s.launcher_errors(root)
            self.assertEqual(log,str(new));self.assertEqual(len(errors),8)
            self.assertTrue(all('current attempt' in line for line in errors))

    def test_brief_report_shows_64_worker_completion_without_process_or_log_dump(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);record=root/'code/insula_ica/records/linux2-trust-stock-filter-example-resume'
            record.mkdir(parents=True);output=root/'outputs';case=output/'sub-1';case.mkdir(parents=True)
            (case/'sub-1.fmriprep.log').write_text('fMRIPrep finished successfully!\n')
            stats=case/'feat_reconstruction/run.feat/stats';stats.mkdir(parents=True)
            (stats/'res4d.nii.gz').touch();(record/'qc').mkdir();(record/'qc/run.json').write_text('{}')
            state=dict(state='functional_reruns_complete_qc_pending',jobs=64,selected_runs=1,output_root=str(output),
                       subjects={'sub-1':dict(state='complete_qc_review_pending',fmriprep_reused=True,numerical_qc_reused=True)})
            (record/'state.json').write_text(json.dumps(state))
            machine=dict(logical_cpus_available=128,load_average=[0,0,0],memory={},vmstat='sample omitted from brief report')
            captured=io.StringIO()
            with patch.object(s,'ROOT',root),patch.object(sys,'argv',['status','--brief']), \
                 patch.object(s,'process_table',return_value=([],[])),patch.object(s,'resources',return_value=machine), \
                 patch.object(s,'launcher_errors',return_value=(None,[])),redirect_stdout(captured):
                s.main()
            report=captured.getvalue()
            self.assertIn('Worker slots recorded: 64',report)
            self.assertIn('finished; inputs ready for QC review',report)
            self.assertIn('Finished FEAT + QC: 1',report)
            self.assertNotIn('Busiest pipeline processes',report)
            self.assertNotIn(machine['vmstat'],report)


if __name__=='__main__':
    unittest.main()
