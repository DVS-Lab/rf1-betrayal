# Interpretation of the linux2 BOLD provenance audit

Reviewed 2026-10-27 EDT. Source results were committed as `5a6d9d43`.
All 92 entries in the supplied `SHA256SUMS` verified, including the terminal
output. The collector reported no read/collection errors. This interpretation
is a later addition, outside that original checksum manifest.

## What the evidence establishes

The BOLD difference is present in all three pilot runs before FEAT prestats.
The historical FEAT midpoint reference was compared to frame 140 of the rebuilt
MNI combined BOLD on the same grid, within the historical FEAT mask:

| Run | Spatial correlation | Relative RMS difference | After best global scaling |
|---|---:|---:|---:|
| sub-10369 Trust 1 | 0.996606 | 4.3241% | 4.3238% |
| sub-10402 Trust 1 | 0.995132 | 5.2011% | 5.1640% |
| sub-10402 Trust 2 | 0.987762 | 8.1693% | 8.1224% |

These discrepancies persist after a fitted global multiplier. They are not
explained by uniform intensity scaling alone. These are single-frame spatial
comparisons, not full-series temporal comparisons. The earlier 5.1591% result
for sub-10402 run 1 used the intersection of historical and rebuilt FEAT masks;
this audit uses the historical mask, accounting for the small metric change.

The included `previous-preprocessing-evidence.json` provides the previously
unshared sub-10369 results. For both completed reconstructed FEAT runs
(sub-10369 run 1 and sub-10402 run 1), the saved nuisance copy equals both the
historical and rebuilt `confoundevs.txt` exactly. Complete designs and contrast
matrices have maximum absolute difference zero, and all FSF settings match
apart from the three relocated paths. Sub-10402 run 2's historical matrix is
preserved exactly but its rebuilt FEAT file/design is still unavailable.
Standalone TEDANA files are not needed to reuse these nuisance values.

All 36 frozen input checks passed across the three runs. JSON comparisons use
SHA256; image comparisons use size/mtime. These establish no detected drift
since the recovery snapshot, not identity with the inputs used years earlier.

The rebuilt environment is fMRIPrep 24.1.1, Nipype 1.8.6, and TemplateFlow
24.2.0. Saved rebuilt settings report multi-echo `curvefit`, slice-time reference
0.5, and no FreeSurfer reconstruction. BOLD summary reports confirm four echoes,
applied slice timing correction, directly measured fieldmap-based distortion
correction, and six-DOF FSL BBR registration.

## What remains unresolved

The historical `fmriprep-24/sub-10369` and `sub-10402` directories are absent
at the recorded upstream location. No historical dataset version, subject
report, or execution config was recovered there. Consequently there are zero
historical/rebuilt execution-config comparison pairs. The archived cluster
script specifies 24.1.1, but it is not an execution record for these old runs.
No exact historical version, seed, combination method, fieldmap input, or
registration setting difference has been demonstrated by this audit.

The individual rebuilt echoes and T1w-space outputs are on different grids
from the historical MNI reference. They were correctly skipped rather than
resampled, so this check does not rule out an echo/combination difference.

## A concrete configuration difference to test

Source comparison with
[`a19c46ef:code/fmriprep24-hpc.sh`](https://github.com/DVS-Lab/rf1-sra-linux2/blob/a19c46ef/code/fmriprep24-hpc.sh)
and its `ses-01.json` shows that the archived command selected all session-01
tasks/runs. The current pilot restricts functional processing to accepted
Trust runs using `--task-id trust` and task/run filters. That is an actual
recovery-versus-archived-script difference, not a demonstrated cause of the BOLD
discrepancy. Current random seeds are recorded; the original seeds are unknown.

Before another cohort rebuild, a controlled single-subject comparison should
separate repeat-run variability from the effect of restoring the archived
functional selection, while holding the container, inputs, resources, and
recorded seeds fixed within the control. This would test a specific hypothesis
without claiming that a fresh rerun can recreate unknown historical seeds.
No such control or production expansion was launched by this audit.
