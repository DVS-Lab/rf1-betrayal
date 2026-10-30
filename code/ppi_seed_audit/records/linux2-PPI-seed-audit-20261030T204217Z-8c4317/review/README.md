# PPI seed audit review — 2026-10-30

## Which analysis used which seed?

The surviving saved UGR analyses include two complete physiological branches: `AIns` and `AIns-vt`. Each branch has 249 accepted manuscript runs from 132 participants, plus four additional saved runs marked excluded. Both branches feed actual saved N132 group models; neither is merely an unused directory. The plain AIns, AIns-socialagnostic and AInsxAQ group families use plain AIns lower-level models. AIns-vt and AInsxAQ-vt groups use the vt models.

The current Linux2 group wrapper `code/run_L3stats-ugr.sh` explicitly selects `ppi_seed-AIns`. The current UGR PPI extraction scripts select `seed-AIns` lower-level results for COPE11 (fairness) and COPE12 (social > nonsocial). Thus the current main extraction route is plain AIns. Surviving parallel group outputs alone cannot establish which branch supplied every final manuscript figure.

The actual saved physiological input is custom EV12 (zero-based real design column 11), with per-run source filenames ending `_AIns.txt` or `_AIns-vt.txt`. No original Trust PPI physiological cases were found in the scanned saved analysis chains. Trust activation covariates in UGR group models are not evidence of a Trust PPI fit.

## Did the physiological files match the saved designs?

Yes, within FEAT input/output precision. All 506 FEAT cached physiological files are byte-identical to their corresponding external files. This includes all 498 accepted-run cases: 249/249 AIns and 249/249 AIns-vt.

| Accepted branch | Cache equals external bytes | Saved physiological column agreement |
|---|---:|---:|
| AIns | 249/249 | 249/249 after float32 custom-EV input rounding and positive centering/scaling |
| AIns-vt | 249/249 | 249/249 with direct text centering/scaling |

The original audit recorded 253 AIns model mismatches because it compared the high-baseline physiological text directly to the saved model at a relative-residual cutoff of 1e-6. Native FSL source reads shape-2 custom EV values into `float tmpf` using `fscanf` with `%f` before centering. Repeating the comparison after that specific float32 input conversion resolves every case without relaxing the cutoff. All 253 cases of each branch, including the excluded cases, agree with the float32 input representation. The minimum original AIns text/model correlation among accepted runs is 0.99999999995242. These flags do not support an overwritten or different physiological timecourse.

`precision-comparison.json` preserves both comparisons and source SHA256 identities for every case. `review_saved_models.py` verifies archived source hashes and reproduces this review without BOLD reads or Linux2. The original audit outputs remain unchanged; their initial mismatch status should be read with this correction. Future audits include the precision-aware check.

The local FSL 6.0.7.19 reduced-model executable probe aborted with a matrix index error; it did not provide a successful native output validation. Its result is retained in `native-float32-precision-probe.json`. The precision interpretation rests on the explicit native source input reader and successful comparisons to all actual historical saved columns, not on that unsuccessful probe. Six targeted audit tests pass, including float32 input handling and genuine stale/reversed-series detection.

## What remains unproven about the spatial seed?

The seed label and the timecourse actually entered into each model are established. The exact historical binary mask footprint is not. The FSF-referenced original BOLD is absent for every physiological case, and the original live `seed-AIns.nii.gz` and `seed-AIns-vt.nii.gz` masks are also absent. The scripts contain variable-based extraction commands, without a surviving literal execution record that certifies mask bytes. Candidate-mask mean extraction therefore could not be reproduced.

The recovered Git `seed-AIns.nii.gz` has 142 voxels, while the surviving `seed-AIns-clusterthresh.nii.gz` has 127. The recovered `seed-AIns-vt.nii.gz` has one voxel. These are distinct candidates; the PPI audit does not establish that any recovered binary version is the exact historical input. In particular, do not assume that the current 127-voxel mask (or its two anatomical-insula voxels) was the mask used for plain AIns PPI.

No participant images, seed files, physiological inputs, fitted designs or production analyses were modified or rerun by this review.
