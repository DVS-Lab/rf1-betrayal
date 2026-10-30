# Design audit: one UGR pmod is redundant with its condition regressor

All 490 accepted run designs were audited. The standardized nonconstant-column
condition number is 19,935,762 for sub-11031 UGR run 1; the next largest over the
other runs is 257.15. The social-high constant and pmod have conditional energy
fractions approximately 3.19e-14 after the other design columns are considered.

The three social-high events all have pmod amplitude -1.458333. Their onsets and
condition belong to the same social-high trials; a constant-amplitude modulator
is a multiple of the condition main effect. Consequently there is no separately
identifiable social-high fairness slope for that run. Increasing a numerical
cutoff can make a finite coefficient appear, but cannot make it identifiable.
The participant's other UGR run and both Trust runs do not have this defect.

The proposed common-high/low model pools the raw convolved high/low pmod columns
within partner type before z-scoring, while retaining four separate condition
main effects and all other task/nuisance columns. Apply the same model to all
249 UGR runs. It assumes a common high/low slope; it is not the previous average
of separately standardized beta coefficients. In this run, social-low variation
identifies the shared social slope. No variation is invented for social-high.
No subject or run need be excluded, and no original FEAT output is edited.
The 241 Trust runs and their partner-specific reciprocity contrasts remain
unchanged. The pooled model can be fit using the cached stage-1 responses only.

## Manuscript cohort and the two input/model issues

The task registry matches code/sublist_n132.txt exactly (132 participants).
It also matches the corrected master QC table's accepted run identities exactly:
241 Trust + 249 UGR = 490, with zero added or missing run keys. No previously
excluded run was added. The input BOLD has been rebuilt where necessary; the
current unsmoothed ICA stream therefore is not a claim of exact historical
image reproduction. The original nuisance identities are preserved.

The two first-level issues identified so far are sub-10817 UGR run 1 acquisition
TR 1.710 seconds versus historical model 1.615, and this sub-11031 UGR run-1
unidentifiable social-high pmod. Their impact on the manuscript's group results
has not been established. The new pooled task model does not automatically
repair or rerun the original paper's voxelwise/group analyses.
