# Independent Trust fMRIPrep controls

Two 24.1.1 runs use separate fresh output/work directories, the same seed, fixed skull-stripping seed and one OpenMP thread. The 25.2.5 case uses the same requested settings. No FEAT, standalone TEDANA, Warpkit fit, or production expansion is run.

These controls measure repeatability under explicit controls. They do not measure the full variation of the earlier unseeded, eight-OpenMP-thread configuration or establish the historical version.

Two fresh unseeded 24.1.1 runs with eight OpenMP threads repeat the sub-10402 pilot settings. They are compared with each other and with the retained pilot BOLD. This is one empirical default-setting pair, not an estimate of the full distribution of possible reruns.

Status: complete_with_errors; 0 input changes detected.

## Historical midpoint frames

| Run | 24 A relative RMS | 24 B relative RMS | 25 A relative RMS |
|---|---:|---:|---:|
| sub-10402_ses-01_task-trust_run-1 | unavailable / different grid | unavailable / different grid | unavailable / different grid |
| sub-10402_ses-01_task-trust_run-2 | unavailable / different grid | unavailable / different grid | unavailable / different grid |
| sub-10478_ses-01_task-trust_run-1 | unavailable / different grid | unavailable / different grid | unavailable / different grid |
| sub-10478_ses-01_task-trust_run-2 | unavailable / different grid | unavailable / different grid | unavailable / different grid |

## All time points within the historical FEAT mask

| Run / comparison | Relative RMS | Voxel-demeaned RMS | Median voxel temporal r | Exact equality |
|---|---:|---:|---:|---|
| sub-10402_ses-01_task-trust_run-1 / v24_repeat | unavailable | unavailable | unavailable | unavailable |
| sub-10402_ses-01_task-trust_run-1 / v25_vs_v24 | unavailable | unavailable | unavailable | unavailable |
| sub-10402_ses-01_task-trust_run-1 / v24_default_repeat | unavailable | unavailable | unavailable | unavailable |
| sub-10402_ses-01_task-trust_run-1 / v24_default_vs_previous_pilot | unavailable | unavailable | unavailable | unavailable |
| sub-10402_ses-01_task-trust_run-2 / v24_repeat | unavailable | unavailable | unavailable | unavailable |
| sub-10402_ses-01_task-trust_run-2 / v25_vs_v24 | unavailable | unavailable | unavailable | unavailable |
| sub-10402_ses-01_task-trust_run-2 / v24_default_repeat | unavailable | unavailable | unavailable | unavailable |
| sub-10402_ses-01_task-trust_run-2 / v24_default_vs_previous_pilot | unavailable | unavailable | unavailable | unavailable |
| sub-10478_ses-01_task-trust_run-1 / v24_repeat | unavailable | unavailable | unavailable | unavailable |
| sub-10478_ses-01_task-trust_run-1 / v25_vs_v24 | unavailable | unavailable | unavailable | unavailable |
| sub-10478_ses-01_task-trust_run-2 / v24_repeat | unavailable | unavailable | unavailable | unavailable |
| sub-10478_ses-01_task-trust_run-2 / v25_vs_v24 | unavailable | unavailable | unavailable | unavailable |

A smaller historical-frame difference means closer agreement for that frame, not proof of the correct historical pipeline or better scientific validity. Grid mismatches are reported without resampling. Numerical equality is reported as evidence, not used to unlock production.

The original 24.1.1 pilot used uncontrolled seeds, no fixed skull-strip seed and eight OpenMP threads; strict historical map equality is not an established acceptance standard for that configuration.

Case exit status and errors are in status.json. Full metrics are in results.json; runtime labels/packages are in containers.json; configs/methods and known versionless outputs are in provenance.json and sources/.

- sub-10402-v24-default-a: name 'case' is not defined
- sub-10402-v24-default-b: name 'case' is not defined
- sub-10402-v25-a: name 'case' is not defined
- sub-10478-v25-a: name 'case' is not defined
- sub-10478-v24-b: name 'case' is not defined
- sub-10402-v24-a: name 'case' is not defined
- sub-10478-v24-a: name 'case' is not defined
- sub-10402-v24-b: name 'case' is not defined
