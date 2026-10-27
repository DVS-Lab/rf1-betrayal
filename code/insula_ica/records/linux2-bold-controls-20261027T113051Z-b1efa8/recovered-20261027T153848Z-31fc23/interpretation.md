# Interpretation of the completed repeat/version controls

Reviewed 2026-10-27 EDT. Numerical results were committed as `63bee268`.
Recovery state is `complete`, all eight preprocessing success markers were
verified, input changes were empty at collection start/end, and all five
original recovery checksum entries verified. This interpretation was added
later and is outside that original checksum manifest.

## Repeatability

The two controlled 24.1.1 runs were exactly equal at every compared voxel and
time point for all four subject/runs: relative RMS 0, maximum absolute difference
0, and voxelwise temporal correlations 1. Comparisons cover the historical
FEAT brain masks; this is not a claim about compressed file bytes or voxels
outside those masks. Both controlled subjects used the same requested seed,
fixed skull-strip seed, and single-thread controls. Saved execution configs
confirmed the controls and identical generated seeds.

The ordinary-setting 24.1.1 pair for sub-10402 did not match:

| Run | Full-series relative RMS difference | Median voxel temporal correlation |
|---|---:|---:|
| Trust 1 | 34.7441% | 0.08898 |
| Trust 2 | 34.3121% | 0.10592 |

Default A was much closer to the earlier retained 24.1.1 pilot: full-series
relative RMS differences 1.9998% and 2.8696%, median voxel temporal correlations
0.99929 and 0.99943. Default B was much farther from the historical frame and
from Default A. Successful preprocessing exit does not establish acceptable
registration or masking. These large differences need QC rather than being
treated as harmless floating-point variation. The affected stage has not been
localized; the comparisons alone do not prove that skull stripping, BOLD
registration, anatomical normalization, or a particular seed is responsible.

## Was 25.2.5 closer to the original?

There is no consistent winner. Relative RMS differences below use the one
retained historical midpoint frame, on the same grid and historical FEAT mask.
Lower is closer; neither the original full 4D BOLD nor original run settings
have been recovered.

| Subject / run | 24.1.1 controlled A (B identical) | 25.2.5 controlled | 24.1.1 ordinary A | 24.1.1 ordinary B |
|---|---:|---:|---:|---:|
| sub-10402 / Trust 1 | 21.9709% | 8.4195% | 5.2750% | 34.4033% |
| sub-10402 / Trust 2 | 23.5976% | 10.5321% | 7.6926% | 35.4935% |
| sub-10478 / Trust 1 | 7.4696% | 10.3276% | Not run | Not run |
| sub-10478 / Trust 2 | 8.1781% | 10.4000% | Not run | Not run |

25.2.5 improves agreement relative to controlled 24.1.1 for sub-10402, but
worsens it for sub-10478. Ordinary 24.1.1 A is closer than either controlled
version for sub-10402. Fixed settings produced repeatable 24.1.1 output but
did not ensure historical agreement. These results do not identify a newer
version as the historical pipeline or justify a cohort-wide version switch.
25.2.5 itself was not repeated in this experiment.

## Next implication

The immediate QC target is sub-10402's inconsistent ordinary-setting repeats
and its controlled output: compare brain masks and anatomical/BOLD registration
and normalization at successive spaces. That can determine where large
differences first appear while retaining the already completed outputs. No
additional preprocessing, model fit, production expansion, or gate change was
performed during this interpretation.
