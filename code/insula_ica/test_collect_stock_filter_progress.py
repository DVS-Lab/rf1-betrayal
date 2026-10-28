import json
from pathlib import Path
import tempfile
import unittest

import collect_stock_filter_progress as c


class DiagnosticTests(unittest.TestCase):
    def test_relevant_evidence_survives_noisy_logs_and_output_is_bounded(self):
        with tempfile.TemporaryDirectory() as directory:
            log=Path(directory)/'log.txt'
            log.write_text('Report - setting before (/output/coreg_boldref.nii.gz) and after (/raw/run-2_sbref.nii.gz) images\n'
                           +'ordinary node message\n'*4000+'fMRIPrep finished successfully!\n')
            info,text=c.log_evidence(log)
            self.assertIn('run-2_sbref',text)
            self.assertIn('finished successfully',text)
            self.assertLess(info['retained_lines'],140)

    def test_collection_preserves_source_state_and_logs_and_records_unknown_status(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);record=root/'record';record.mkdir();output=root/'output'
            subject='sub-10001';case=output/subject;case.mkdir(parents=True)
            log=case/(subject+'.fmriprep.log');log.write_text('fMRIPrep finished successfully!\n')
            state=dict(output_root=str(output),subjects={subject:dict(state='failed',error='unknown reference')})
            (record/'state.json').write_text(json.dumps(state))
            before=(record/'state.json').read_bytes(),log.read_bytes()
            target=root/'diagnostic';target.mkdir();counts=c.collect(record,target)
            self.assertEqual(counts['failed'],1)
            self.assertEqual(before,((record/'state.json').read_bytes(),log.read_bytes()))
            diagnostic=json.loads((target/'diagnostic.json').read_text())
            self.assertFalse(diagnostic['analysis_jobs_changed'])
            self.assertEqual(diagnostic['errors'][0]['error'],'unknown reference')


if __name__=='__main__':
    unittest.main()
