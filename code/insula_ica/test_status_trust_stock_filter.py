import json
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


if __name__=='__main__':
    unittest.main()
