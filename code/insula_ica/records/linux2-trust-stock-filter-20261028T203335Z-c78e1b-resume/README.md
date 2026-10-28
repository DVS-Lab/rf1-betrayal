# Trust queue resumed with existing outputs

Workers: 64; new fMRIPrep limit: 2. Completed: 108/108; failed: 0.

Existing FEAT children were allowed to finish. Successful fMRIPrep runs were not repeated. Completed FEAT/QC outputs were retained with model, confound and header checks. Reference selection was verified from actual log evidence. Unknown or mismatched selection remains a failure. New stages reserve 28 GiB (fMRIPrep) or 4 GiB (FEAT/QC), plus 12 GiB available-memory headroom. No ICA or exclusions were performed.

## Review of the three flagged runs — 2026-10-28

All three flags concern similarity to the historical reference image (spatial r below .95). All three have 100% left, right and bilateral insula support under the FEAT-mask/finite, nonconstant-residual definition. Correct run-2 SBRef selection is recorded for each.

| Trust run | Left / right / bilateral support | Historical-image r | Montage review |
| --- | --- | --- | --- |
| [sub-10720 run 2](qc/sub-10720_ses-01_task-trust_run-2.png) | 100% / 100% / 100% | .9363 | No visible gross displacement or insula truncation |
| [sub-11326 run 2](qc/sub-11326_ses-01_task-trust_run-2.png) | 100% / 100% / 100% | .9492 | No visible gross displacement or insula truncation |
| [sub-11433 run 2](qc/sub-11433_ses-01_task-trust_run-2.png) | 100% / 100% / 100% | .9239 | No visible gross displacement or insula truncation |

Recommendation: retain these three runs; this review does not indicate another preprocessing rerun or an exclusion. These are historical-image similarity flags, rather than demonstrated registration failures. The available montages compare one BOLD frame with the historical reference at six axial MNI slices. Native BOLD-to-T1 overlays and signal-quality measures were not included in this review, so subtle registration differences or local signal loss are not ruled out. Full mask support alone does not establish signal quality.

See [flagged-run-review.json](flagged-run-review.json) for asset hashes, metrics, review scope and recommendations. The source flags, pipeline state and analysis inputs remain unchanged; this is a recorded review of the three flags.
