# Trust insula coverage across the current sample

Calculated on 2026-10-28 from the 241 QC JSON files in this record (source commit `6a499d86`). The JSON companion records every source SHA256, run, participant summary and threshold calculation.

## Definition and scope

Coverage = atlas-insula voxels inside the reconstructed FEAT mask with finite, nonconstant residuals / all atlas-insula voxels. The atlas is HarvardOxford cortical maxprob-thr25, resampled to the output grid (391 left, 427 right, 818 bilateral voxels). This measures voxel support. It does not establish correct tissue correspondence or adequate signal strength/tSNR. These numbers apply to Trust, not yet to UGR or a final common two-task mask.

## Distribution

| ROI | Mean across runs | Median | Minimum | Runs at 100% |
| --- | ---: | ---: | ---: | ---: |
| left | 99.57% | 100.0% | 34.53% | 239/241 |
| right | 99.56% | 100.0% | 37.00% | 239/241 |
| bilateral | 99.56% | 100.0% | 35.82% | 239/241 |

239/241 runs (99.17%) from 131/132 participants have full coverage in all three ROIs. Both exceptions are sub-10834. The participant-weighted bilateral mean (averaging runs within participant first) is 99.60%; the run-weighted bilateral mean is 99.56%.

| sub-10834 run | Left | Right | Bilateral |
| --- | ---: | ---: | ---: |
| run-1 | 61.9% | 56.9% | 59.3% |
| run-2 | 34.5% | 37.0% | 35.8% |

## Proposed rule sensitivity, not applied exclusions

Require the indicated coverage separately in left and right insula. A bilateral average alone can conceal unilateral loss.

| Minimum in each hemisphere | Runs passing | Participants with at least one passing Trust run | Participants with all accepted Trust runs passing |
| --- | ---: | ---: | ---: |
| 50% | 240/241 | 132/132 | 131/132 |
| 60% | 239/241 | 131/132 | 131/132 |
| 70% | 239/241 | 131/132 | 131/132 |
| 80% | 239/241 | 131/132 | 131/132 |
| 90% | 239/241 | 131/132 | 131/132 |
| 95% | 239/241 | 131/132 | 131/132 |
| 99% | 239/241 | 131/132 | 131/132 |
| 100% | 239/241 | 131/132 | 131/132 |

The 80%, 90% and 95% rules select the same runs. This sample cannot establish which of these thresholds is preferable. The existing config `minimum_roi_coverage=0.8` is a final common-mask gate; per-run regeneration flags do not yet exclude anyone.

## Recommendation for an analysis rule

1. Treat gross registration failure as a separate processing/QC failure: attempt correction before deciding on exclusion. Coverage alone cannot rescue a misplaced ROI.
2. Retain 80% in each hemisphere as a provisional screening threshold, and document any final run eligibility threshold before examining ICA/betrayal effects. Report 90%/95% sensitivity. This is a proposed project rule, not an externally established universal cutoff.
3. Assess preprocessed BOLD signal dropout and tSNR separately from mask support; do not calculate signal-quality tSNR from task-regressed/prewhitened residuals.
4. Apply the same finalized rule to both tasks. If a run remains unusable, exclude it at run level; retain participants only if at least one usable run remains for each required task. Require all three ROIs for a shared left/right/bilateral comparison sample.
5. Quantify the final common spatial mask separately across retained runs/tasks. If sub-10834 were set aside, the remaining 239 Trust runs have 100% common insula support under this metric. UGR and the paired common mask remain to be measured.

No exclusion rules, cohort membership or source images were changed. If both sub-10834 Trust runs remain unusable, an analysis requiring Trust and UGR would have at most N=131 from this cohort.

Published QC guidance distinguishes normalization/coregistration, mask errors and signal loss: [Provins et al., Quality control in functional MRI studies with MRIQC and fMRIPrep](https://www.frontiersin.org/journals/neuroimaging/articles/10.3389/fnimg.2022.1073734/full).
