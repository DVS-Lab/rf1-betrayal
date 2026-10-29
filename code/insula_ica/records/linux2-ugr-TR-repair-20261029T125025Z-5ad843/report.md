# UGR timing repair and project timing audit

State: TR_repair_complete_QC_review_pending. Corrected sub-10817_ses-01_task-ugr_run-1: 1.615 → 1.710 seconds. Raw JSON/NIfTI files were not edited. Original task events in seconds and nuisance values were retained; the HRF/task design was regenerated and required to change. Existing historical fits remain available.

Cohort audit: 490 accepted runs; 1 model/acquisition TR mismatches; 0 unresolved timing checks. See cohort-timing-audit.json.

Trust for sub-10817:

- sub-10817_ses-01_task-trust_run-1: acquisition 1.615; model 1.615; matches.
- sub-10817_ses-01_task-trust_run-2: acquisition 1.615; model 1.615; matches.

Same-BOLD sm-5 fits isolate the effect of the TR change; matched-TR-statistics-comparison.json contains all COPE/ZSTAT contrasts. The 1.615-second fit is a diagnostic control and must not be promoted as corrected output. corrected-analysis-input.json identifies the corrected sm-5 candidate.

Combined sm-0 completion: 490/490. QC review is pending. Original group statistics, PPI, RSA and existing ICA fits have not been rerun or silently replaced. The corrected input registry is ready for explicit propagation into those dependent analyses. Mixed acquisition TRs must remain explicit when preparing later ICA inputs/frequency reports.
