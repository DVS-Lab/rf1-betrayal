# Original seed/group-map audit: findings and activation repair

Reviewed the latest of two equivalent published audits,
`linux2-original-group-audit-20261030T113752Z-f372a9` (8cd43bee).
All 249 original sm-5 UGR run designs were read, with no missing models.
There is one acquisition/model TR mismatch and one near-dependent task design.
The audit collected nine original group cases: UGR fairness and social-minus-
nonsocial contrasts, plus Trust rec-minus-def, each in ones/full/covs-vt models.

## What the seed files establish

The live Linux2 clusterthresh mask has exactly the same image support and
affine as the earlier plotted tracked file, despite different compressed-file
hashes. It is exactly the 127-voxel connected component of the intersection of
the original **ones-model negative UGR COPE11 and negative Trust COPE10** group
cluster masks. Their full intersection is exactly the saved 365-voxel
`ugr-trust_intersection.nii.gz` support. Thus the available clusterthresh mask's
origin is verified against actual group outputs, not inferred only from names.

Those same-sign task-map intersections have these sizes on the original grid:

| Group model | Positive intersection | Negative intersection | Negative-intersection HO25 insula voxels |
|---|---:|---:|---:|
| ones | 560 | 365 | 2 |
| full | 401 | 328 | 2 |
| covs-vt | 0 | 7 | 0 |

The 127-voxel seed has only two HO25 insular-cortex voxels. Its plotted frontal/
opercular location is not a Git-versus-Linux2 data/header discrepancy. These
results do not establish an alternative large anatomical-insula task conjunction.

The live masks directory also contains a 497-voxel anatomical left-insula mask,
`left_insula_maxprob-thr25_anatomical_refgrid.nii.gz`, exactly equal to the saved
`left_insula_maxprob-thr25_dualtask_coverage_mask.nii.gz`. All 497 voxels carry
the HO25 insular label in this comparison. That is an anatomical/coverage ROI,
not the functional task-effect conjunction. The RSA functional intersection
also overlaps anatomical insula by only two voxels.

The file used to extract the paper's actual PPI physiological timecourses is
still unverified. The live directory does not contain `seed-AIns.nii.gz` or
`seed-AIns-vt.nii.gz`, although analysis/extraction scripts reference them.
Verifying the available cluster's group-map origin does not prove those
historical timecourses used the same mask. Do not replace a PPI seed silently.

Numerical receipts and audited source identities are in `verification.json`.

## Exactly what affects UGR

- sub-10817 UGR run 1: historical TR=1.615 s, current acquisition TR=1.710 s.
  It is a single-accepted-run participant for UGR. The existing corrected sm-5
  candidate can feed group analysis directly, as its old single-run estimate did.
  It uses regenerated BOLD, so comparison with historical images includes that
  preprocessing difference. The previously completed same-regenerated-BOLD
  wrong-TR control separately measures the model-timing effect.
- sub-11031 UGR run 1: social-high condition and constant-amplitude pmod are
  near-dependent at the 1e-6 diagnostic cutoff; three equal-amplitude events
  provide no independent social-high fairness slope. Affected original contrasts
  are **5, 6, 9, 10, 11, 12, 14, 15 and 17**, not only the individual pmod.
  This includes both main manuscript UGR contrasts, COPE11 and COPE12. Run 2
  passes contrast diagnostics.
- No other original UGR run was flagged. Both Trust runs for these participants
  have no corresponding issue in the earlier cohort audit.

## Separate original-model activation rerun

The prepared launcher is `bash code/ugr_reanalysis/start_activation.sh`.
It retains the original high/low model, not the new pooled common-slope model.
For the six audited main UGR group fits, each 132-subject input list changes only:

1. sub-10817's direct L1 input: reuse the completed corrected-TR sm-5 candidate.
2. sub-11031's original two-run L2 input for COPE11/12: use the original valid
   run-2 L1 COPE and its matching VARCOPE directly. This follows the same
   single-run input convention already present in the original group analyses.

The remaining 130 participant inputs are unchanged, and N remains 132. No
first-level or participant fixed-effects refit is necessary for these six group
fits because the needed corrected/valid single-run estimates already exist.
Unaffected original-model contrasts retain both sub-11031 runs if analyzed
later. This is contrast-specific run availability defined by the design audit,
not by image findings or task-effect size.

The runner preserves the saved original group EV values, covariates, weights
and thresholds and verifies regenerated design matrices against the originals.
Only output location, isolated input-file paths and, if necessary, a stale
machine-specific standard-reference path are changed. Only two participant
estimates change; the other 130 source estimates are reused. The reference path maps
the same FSL MNI152 2-mm basename to the local installation; it is recorded.

All 132 input paths are redirected into isolated participant/contrast packets
under `<rerun-root>/inputs/sub-*/cope-*.feat/`. These contain copies of small 3D
COPE/VARCOPE, degrees-of-freedom, mask and background images; source standard-space
maps are used when available. Raw first-level maps require an identity transform
and matching grid, except the explicitly corrected candidate already known to
be in MNI space. Existing voxelwise DOF is retained; if only scalar L1 residual
DOF exists, it is expanded to a constant 3D map without estimating new DOF.
Contrast-height metadata is retained where present.

FSL's own `featregapply` exits immediately for the packets' `design.lev` marker,
and `feat5:find_std` reads their local files. This uses FSL's already-standard
input convention, avoiding registration, interpolation or writes into original
source FEAT directories. The behavior was checked in local FSL source. All
packet sources/copies are identity-checked; originals remain read only.

Full group outputs go under
`derivatives/ugr_reanalysis/<timestamp>/activation/<ones|full|covs-vt>/cope-<11|12>.gfeat`.
Code and compact manifests, model snapshots, commands, group maps, comparisons
and errors go under `code/ugr_reanalysis/records/linux2-ugr-activation-repair-*`.
Original paper outputs and source estimates are read only. Completed reruns
resume using validated input/output identities; interrupted group directories
are preserved. A pilot of ones/COPE11 must pass before expanding to five others.

This phase covers the two main activation contrasts across the three saved
model/threshold variants. PPI, other activation contrasts, connectivity,
extractions and manuscript revisions are not marked repaired by this phase.
PPI needs its own design/timing and actual-seed provenance review, followed by
dependent reruns where affected. Changing a seed would require new physiological
timecourses and downstream fits across relevant participants.

Validation: all eight audit/repair tests pass. The six actual saved group FSFs
were exercised in a filesystem fixture with all 132 subject input paths; exactly
two source estimates per fit changed and every covariate/contrast setting was preserved.
All 132 file pointers are isolated from source directories. Packet staging and
scalar DOF expansion were checked against unchanged original source identities.
A native FSL 6.0.7.19 `featregapply` probe on a synthetic staged packet exited
successfully and left all source/packet files unchanged; its receipt is saved in
`native-registration-noop-probe.json`.
Completion/resume and partial-output preservation were tested with simulated
group output. Actual FSL group execution remains gated by the Linux2 pilot.
