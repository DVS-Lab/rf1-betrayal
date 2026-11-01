# Matched 5 mm functional-cluster arm

Same source BOLD, corrected run-specific TR, full temporal model/contrasts, original nuisance matrix, cohort, fixed ICA mask and discovery settings as the sm0 arm. Only FEAT spatial smoothing is changed to 5 mm. Reuse exact verified sm5 fits; rebuild only missing/mismatched first-level fits in this separate root. No fMRIPrep/TEDANA or unsmoothed rerun. Smooth whole-brain BOLD in FEAT before masking; do not smooth an already masked ROI.

Both retained joint BOLD-scale and native FILM residuals; physical 4D masking; 3/5/7/automatic concat ICA with native variance normalization, identical balancing and serialized automatic PCA. Require sm5 support at every frozen sm0 ICA voxel; no erosion, voxel/subject selection or mask substitution. All masks/inputs/settings and candidate rejection reasons are recorded.

Raw FEAT-input task BOLD, mean estimates and anatomical-insula results are reused. New sm5 component maps undergo the same spatial regression on that identical task-preserved input and pooled-UGR/full task+nuisance standardized OLS. Five bars per component, pointwise participant CIs. Spatial-only cross-smoothing matching and paired task differences are descriptive; components are not assumed to correspond by index. This remains exploratory task-selected ROI characterization.
