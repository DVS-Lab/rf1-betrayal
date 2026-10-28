# sub-10834 registration findings and correction candidate

Reviewed on 2026-10-28. Both coregistration panels show gross mismatch between BOLD and the participant's own T1w anatomy. The T1w image and its WM contours appear broadly consistent in these native-space panels. Therefore the BOLD-to-anatomical stage already fails, before applying the T1w-to-MNI transform. This does not certify the unreviewed baseline normalization panel.

The BOLD reference/SDC report also shows poor reference masking, with much of the superior/anterior brain outside the displayed mask. The reference-mask pathway is a plausible contributor, not an established sole cause. The original SDC report compares the distorted input with the corrected reference; contrast differences alone do not establish a distortion-field error. Do not infer a need to regenerate Warpkit or TEDANA from these images.

## SBRef filter error

The filter we generated overrode SBRef `run` with the full accepted-run list. fMRIPrep 24.1.1 applies these filters as overrides of the current BOLD run's entities, so that list can permit cross-run selection. The retained log explicitly uses run 1's SBRef in run 2's SDC/report path. This is a mistake in our filter, now corrected by inheriting the current run for SBRef matching. It is potentially relevant beyond this participant; a separate extraction of actual selected references from all retained cohort logs accompanies the correction attempt. A different-run reference is a provenance concern, not by itself proof of an unusable registration. Missing log evidence remains unknown.

Primary code reference: [fMRIPrep 24.1.1 get_sbrefs and _select_ref](https://github.com/nipreps/fmriprep/blob/24.1.1/fmriprep/workflows/bold/fit.py). Original completed outputs and frozen records have not been edited.

## Candidate test

Run only sub-10834's two accepted Trust runs with the same 24.1.1 image and production settings plus `--ignore sbref`, generating a BOLD-derived reference. Use fresh output and work directories. Refit the original full FEAT task+nuisance models with the exact preserved confounds. Save before/candidate coverage, historical reference correlations, both montages, and current/candidate registration reports. The collector now recognizes actual normalization and brain-mask report filenames that were missed by its original keyword filter.

This is a candidate correction, not a promise that ignoring SBRef fixes the problem. No candidate is automatically adopted, and no participant is excluded. Visual confirmation of anatomy/insula alignment is required regardless of the numeric coverage result.

PNG snapshots in this record show the two original SVG flicker layers separately. Embedded WebP images were decoded to lossless PNG for SVG rendering; the source SVGs remain unchanged. These snapshots are display conversions, not resampled or registered MRI data.
