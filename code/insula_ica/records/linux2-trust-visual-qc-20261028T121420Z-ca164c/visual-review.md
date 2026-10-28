# Review of the 20 flagged Trust montages

Reviewed all 20 exported PNGs individually on 2026-10-28, from commit `7418f61b`.

**sub-10834 runs 1 and 2: hold these regenerated inputs for registration correction.** Brain position is grossly displaced relative to the atlas insula outlines. Parts of the insula lie outside brain tissue. Bilateral coverage is 59.3% and 35.8%. The retained historical images have substantially more plausible alignment. These observations go beyond ordinary differences in intensity or exact historical reproduction.

**Other 18 runs: no comparable gross misalignment seen in the displayed slices.** Each has 100% insula coverage. The single-frame, six-axial-slice review does not certify all registration stages, motion or full temporal coverage. In particular, correlation alone is not a pass/fail decision; sub-11530 run 2 has substantial historical intensity/anatomy differences without the gross atlas displacement seen in sub-10834.

**Next:** inspect sub-10834 anatomical normalization and each BOLD-to-anatomical registration in the fMRIPrep reports to identify the failed stage. Both runs being affected makes a shared stage worth checking, but the montage does not establish which transform failed. Preserve the current outputs and fit any correction separately. No participant/run exclusions, reruns, or ICA were performed by this review.

| Run | Montage finding |
| --- | --- |
| [sub-10834_ses-01_task-trust_run-2](sub-10834_ses-01_task-trust_run-2.png) | Gross misalignment; hold current inputs |
| [sub-10834_ses-01_task-trust_run-1](sub-10834_ses-01_task-trust_run-1.png) | Gross misalignment; hold current inputs |
| [sub-10640_ses-01_task-trust_run-2](sub-10640_ses-01_task-trust_run-2.png) | No comparable gross misalignment seen in displayed slices |
| [sub-10642_ses-01_task-trust_run-2](sub-10642_ses-01_task-trust_run-2.png) | No comparable gross misalignment seen in displayed slices |
| [sub-10674_ses-01_task-trust_run-2](sub-10674_ses-01_task-trust_run-2.png) | No comparable gross misalignment seen in displayed slices |
| [sub-10700_ses-01_task-trust_run-2](sub-10700_ses-01_task-trust_run-2.png) | No comparable gross misalignment seen in displayed slices |
| [sub-10767_ses-01_task-trust_run-2](sub-10767_ses-01_task-trust_run-2.png) | No comparable gross misalignment seen in displayed slices |
| [sub-10800_ses-01_task-trust_run-2](sub-10800_ses-01_task-trust_run-2.png) | No comparable gross misalignment seen in displayed slices |
| [sub-10801_ses-01_task-trust_run-2](sub-10801_ses-01_task-trust_run-2.png) | No comparable gross misalignment seen in displayed slices |
| [sub-10860_ses-01_task-trust_run-2](sub-10860_ses-01_task-trust_run-2.png) | No comparable gross misalignment seen in displayed slices |
| [sub-10875_ses-01_task-trust_run-2](sub-10875_ses-01_task-trust_run-2.png) | No comparable gross misalignment seen in displayed slices |
| [sub-10918_ses-01_task-trust_run-2](sub-10918_ses-01_task-trust_run-2.png) | No comparable gross misalignment seen in displayed slices |
| [sub-10924_ses-01_task-trust_run-2](sub-10924_ses-01_task-trust_run-2.png) | No comparable gross misalignment seen in displayed slices |
| [sub-10966_ses-01_task-trust_run-2](sub-10966_ses-01_task-trust_run-2.png) | No comparable gross misalignment seen in displayed slices |
| [sub-10969_ses-01_task-trust_run-2](sub-10969_ses-01_task-trust_run-2.png) | No comparable gross misalignment seen in displayed slices |
| [sub-11060_ses-01_task-trust_run-2](sub-11060_ses-01_task-trust_run-2.png) | No comparable gross misalignment seen in displayed slices |
| [sub-11068_ses-01_task-trust_run-2](sub-11068_ses-01_task-trust_run-2.png) | No comparable gross misalignment seen in displayed slices |
| [sub-11325_ses-01_task-trust_run-2](sub-11325_ses-01_task-trust_run-2.png) | No comparable gross misalignment seen in displayed slices |
| [sub-11337_ses-01_task-trust_run-2](sub-11337_ses-01_task-trust_run-2.png) | No comparable gross misalignment seen in displayed slices |
| [sub-11530_ses-01_task-trust_run-2](sub-11530_ses-01_task-trust_run-2.png) | No comparable gross misalignment seen in displayed slices |
