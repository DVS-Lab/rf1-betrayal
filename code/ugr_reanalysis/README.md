# Original UGR analysis and functional-seed reassessment

The first six repaired activation groups are complete. See the
[original-versus-corrected review](records/linux2-ugr-activation-repair-20261030T130959Z-0fac9f/review/README.md):
whole-brain patterns remain similar, one small full-model social-minus-nonsocial
cluster loses threshold survival, and anatomical-insula task overlap does not
increase. This does not mark PPI or other activation contrasts repaired.

The original-input audit is complete: see
[findings and the separate activation rerun](records/2026-10-30-original-input-review/README.md).
Run `bash code/ugr_reanalysis/start_activation.sh` to reuse the completed
corrected-TR sm-5 fit and original valid run-2 estimates for the two main UGR
contrasts. It reruns their six saved group models in a fresh versioned
`derivatives/ugr_reanalysis/` directory, preserving N132 and original group
designs/thresholds, without new preprocessing, L1 or L2 fitting. It uses nohup,
automatically resumes validated completed groups and publishes its own record.
Small 3D input packets are staged inside the new output tree so FSL cannot update
registration caches in the original source directories. No 4D data is copied.
Check with the existing Linux2 Python and
`code/ugr_reanalysis/rerun_activation.py --brief`.
PPI and other contrasts remain separately tracked work.

First recover and inspect Jacob's saved group maps, actual contrast definitions,
group inputs and live seed files. The currently tracked clusterthresh mask is
not established as the manuscript's final seed. Historical seed candidates were
recovered separately under
`code/insula_ica/records/2026-10-30-seed-provenance/`; none was promoted into
`masks/` or a production analysis.

## Original-input audit on Linux2

```bash
bash code/ugr_reanalysis/start_audit.sh
```

Uses the existing tedana-26.0.3 Python environment and nohup. It records and
publishes only its own new `code/ugr_reanalysis/records/linux2-original-group-audit-*`
directory. It reads the current project `masks/` and direct `masks/archive/`,
and the original `derivatives/fsl/L3-act` N132 activation outputs for UGR COPE11/12
and Trust COPE10. It captures available saved group designs, input paths,
contrast names, threshold settings and both zstat1/2 unthresholded and
cluster-thresholded maps. Original paths and hashes distinguish copies with the
same basename. No cross-project filesystem scan is performed.

It also audits all 249 accepted **original sm-5** UGR designs, explicitly using
their original manuscript locations rather than the newer ICA sm-0 registry's
image paths. For the two affected participants it copies small original models,
confounds and event files. No participant 4D image is copied, and no FEAT,
preprocessing, L2, L3 or PPI job is launched. Missing originals are reported.

Numerical contrast checks retain original raw-unit contrast definitions when
scaling design columns for SVD. Three cutoffs are reported to distinguish
rounding-induced near dependence from a well-conditioned model. Nonestimable
contrasts and acquisition/model TR disagreements are diagnostic flags requiring
review, not automatic exclusions or model replacements.

Local validation: three targeted numerical tests pass (raw contrast-unit
conversion under a dependent design, constant-column contrast handling and
full-rank contrast estimability). Shell syntax and Python help are checked.
An end-to-end synthetic original-layout probe also audited 249 models, collected
both task group-map cases, preserved saved contrast names, and detected the
deliberately introduced TR mismatch and aliased pmod run exactly.

## Scope of the subsequent repair

The two established issues require first-level review before group rerunning:

1. **sub-10817 UGR run 1:** current raw TR is 1.710 s; the historical model used
   1.615 s. A versioned corrected sm-5 fit and same-regenerated-BOLD wrong-TR
   control already exist in the TR-repair record. Their matched comparison
   isolates model timing from preprocessing differences. That corrected fit has
   not been promoted into the original participant/group analysis.
2. **sub-11031 UGR run 1:** social-high has three identical pmod amplitudes and
   hence no independent social-high slope. The main UGR fairness contrast gives
   this pmod nonzero weight. Audit the original sm-5 design and exact contrast
   definitions before choosing contrast-specific run handling or a common-slope
   sensitivity model. The pooled model that fixed the ICA task plots has not
   repaired the manuscript's voxelwise COPE/VARCOPE or PPI results.

An L3-only rerun would reuse affected first-level estimates. The repair must
propagate corrected estimates and appropriate variance estimates through the
affected participant L2 inputs and relevant UGR group activation contrasts.
The original-model repair and a pooled high/low common-slope sensitivity are
different analyses: pooling changes the model assumptions and must not silently
replace the original paper's high/low effects. Rebuild the two-task intersection
only after verifying matching effects, signs, thresholds and spatial references.
Retain the N132 cohort and document contrast-specific usable runs.

If the functional seed changes, its physiological timecourses and downstream
PPI models also need rerunning; changing group masks alone cannot update those
models. Task/design timing corrections also require regeneration of affected
PPI regressors rather than reuse of stale timecourse/deconvolution caches.

Current L3 wrappers unconditionally delete their target output before rerunning.
Use separately versioned output roots and explicit input substitutions for the
repair; do not execute those wrappers as an in-place manuscript overwrite.
Compare original and repaired unthresholded effects, uncertainty and thresholded
results before revising manuscript statements. Group-level impact is not yet
established.
