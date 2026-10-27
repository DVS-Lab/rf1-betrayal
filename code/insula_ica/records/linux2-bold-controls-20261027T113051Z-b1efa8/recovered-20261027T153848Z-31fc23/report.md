# Independent Trust fMRIPrep controls

Two 24.1.1 runs use separate fresh output/work directories, the same seed, fixed skull-stripping seed and one OpenMP thread. The 25.2.5 case uses the same requested settings. No FEAT, standalone TEDANA, Warpkit fit, or production expansion is run.

These controls measure repeatability under explicit controls. They do not measure the full variation of the earlier unseeded, eight-OpenMP-thread configuration or establish the historical version.

Two fresh unseeded 24.1.1 runs with eight OpenMP threads repeat the sub-10402 pilot settings. They are compared with each other and with the retained pilot BOLD. This is one empirical default-setting pair, not an estimate of the full distribution of possible reruns.

Status: complete; 0 input changes detected.

## Historical midpoint frames

| Run | 24 A relative RMS | 24 B relative RMS | 25 A relative RMS |
|---|---:|---:|---:|
| sub-10402_ses-01_task-trust_run-1 | 21.97086% | 21.97086% | 8.41955% |
| sub-10402_ses-01_task-trust_run-2 | 23.59764% | 23.59764% | 10.53213% |
| sub-10478_ses-01_task-trust_run-1 | 7.46964% | 7.46964% | 10.32761% |
| sub-10478_ses-01_task-trust_run-2 | 8.17815% | 8.17815% | 10.39998% |

## All time points within the historical FEAT mask

| Run / comparison | Relative RMS | Voxel-demeaned RMS | Median voxel temporal r | Exact equality |
|---|---:|---:|---:|---|
| sub-10402_ses-01_task-trust_run-1 / v24_repeat | 0.0 | 0.0 | 1.0 | True |
| sub-10402_ses-01_task-trust_run-1 / v25_vs_v24 | 0.22037587563690572 | 0.8531452815433385 | 0.7302743497724706 | False |
| sub-10402_ses-01_task-trust_run-1 / v24_default_repeat | 0.34744052650113433 | 1.2283164623671963 | 0.08898414044397532 | False |
| sub-10402_ses-01_task-trust_run-1 / v24_default_vs_previous_pilot | 0.01999783570004132 | 0.07289876397404121 | 0.9992873329313148 | False |
| sub-10402_ses-01_task-trust_run-2 / v24_repeat | 0.0 | 0.0 | 1.0 | True |
| sub-10402_ses-01_task-trust_run-2 / v25_vs_v24 | 0.22016352659049904 | 0.8635854143777139 | 0.7301031728063365 | False |
| sub-10402_ses-01_task-trust_run-2 / v24_default_repeat | 0.3431211781588327 | 1.2393542451530088 | 0.10591786749568456 | False |
| sub-10402_ses-01_task-trust_run-2 / v24_default_vs_previous_pilot | 0.028696133906918234 | 0.13946721479558466 | 0.9994290459044316 | False |
| sub-10478_ses-01_task-trust_run-1 / v24_repeat | 0.0 | 0.0 | 1.0 | True |
| sub-10478_ses-01_task-trust_run-1 / v25_vs_v24 | 0.08697585763002941 | 0.262926856190319 | 0.9944514361209216 | False |
| sub-10478_ses-01_task-trust_run-2 / v24_repeat | 0.0 | 0.0 | 1.0 | True |
| sub-10478_ses-01_task-trust_run-2 / v25_vs_v24 | 0.08832225565995618 | 0.2594894995039971 | 0.99466496607528 | False |

A smaller historical-frame difference means closer agreement for that frame, not proof of the correct historical pipeline or better scientific validity. Grid mismatches are reported without resampling. Numerical equality is reported as evidence, not used to unlock production.

The original 24.1.1 pilot used uncontrolled seeds, no fixed skull-strip seed and eight OpenMP threads; strict historical map equality is not an established acceptance standard for that configuration.

Case exit status and errors are in status.json. Full metrics are in results.json; runtime labels/packages are in containers.json; configs/methods and known versionless outputs are in provenance.json and sources/.

