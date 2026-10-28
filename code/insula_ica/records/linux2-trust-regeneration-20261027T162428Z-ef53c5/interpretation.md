# Interpretation of completed Trust regeneration

Reviewed from GitHub commit `6a499d86` on 2026-10-28.

All 132 participants / 241 accepted Trust runs completed preprocessing and joint-model FEAT preparation. All 241 QC records confirm exact numeric equality of the preserved confounds and finite, nonconstant residual data. No participant or run was excluded.

Twenty runs were flagged by the historical-frame spatial correlation threshold (<0.95). Of those, both sub-10834 runs also have insula coverage below 80%: run 1 left/right/bilateral 61.9%/56.9%/59.3%, run 2 34.5%/37.0%/35.8%. These two runs need priority visual review. The other 18 flagged runs have 100% anatomical insula coverage and are all run 2. Differences from a historical image alone do not establish a registration failure.

No visual QC decisions are assigned from these JSON records. Inspect comparison images and the fMRIPrep reports before deciding whether sub-10834 needs a focused correction. Do not remove participants or reopen exact historical reproduction by default. ICA has not started; the residual-scale decision remains separate.

| Run | Historical frame r | Left coverage | Right coverage | Bilateral coverage |
| --- | ---: | ---: | ---: | ---: |
| sub-10834_ses-01_task-trust_run-2 | 0.3054 | 34.5% | 37.0% | 35.8% |
| sub-10834_ses-01_task-trust_run-1 | 0.3625 | 61.9% | 56.9% | 59.3% |
| sub-10640_ses-01_task-trust_run-2 | 0.9336 | 100.0% | 100.0% | 100.0% |
| sub-10642_ses-01_task-trust_run-2 | 0.8939 | 100.0% | 100.0% | 100.0% |
| sub-10674_ses-01_task-trust_run-2 | 0.9116 | 100.0% | 100.0% | 100.0% |
| sub-10700_ses-01_task-trust_run-2 | 0.8781 | 100.0% | 100.0% | 100.0% |
| sub-10767_ses-01_task-trust_run-2 | 0.9382 | 100.0% | 100.0% | 100.0% |
| sub-10800_ses-01_task-trust_run-2 | 0.9037 | 100.0% | 100.0% | 100.0% |
| sub-10801_ses-01_task-trust_run-2 | 0.9103 | 100.0% | 100.0% | 100.0% |
| sub-10860_ses-01_task-trust_run-2 | 0.9458 | 100.0% | 100.0% | 100.0% |
| sub-10875_ses-01_task-trust_run-2 | 0.9455 | 100.0% | 100.0% | 100.0% |
| sub-10918_ses-01_task-trust_run-2 | 0.9096 | 100.0% | 100.0% | 100.0% |
| sub-10924_ses-01_task-trust_run-2 | 0.9341 | 100.0% | 100.0% | 100.0% |
| sub-10966_ses-01_task-trust_run-2 | 0.9365 | 100.0% | 100.0% | 100.0% |
| sub-10969_ses-01_task-trust_run-2 | 0.8496 | 100.0% | 100.0% | 100.0% |
| sub-11060_ses-01_task-trust_run-2 | 0.9097 | 100.0% | 100.0% | 100.0% |
| sub-11068_ses-01_task-trust_run-2 | 0.9459 | 100.0% | 100.0% | 100.0% |
| sub-11325_ses-01_task-trust_run-2 | 0.9394 | 100.0% | 100.0% | 100.0% |
| sub-11337_ses-01_task-trust_run-2 | 0.9166 | 100.0% | 100.0% | 100.0% |
| sub-11530_ses-01_task-trust_run-2 | 0.7136 | 100.0% | 100.0% | 100.0% |
