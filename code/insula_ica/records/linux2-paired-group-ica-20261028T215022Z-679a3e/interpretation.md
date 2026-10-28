# Initial paired ICA attempt: incomplete

Reviewed the published record at commit `311ae4d7` on 2026-10-28.

- All 490 accepted runs completed residual preparation.
- All six ROI/representation collections contain 264 physically masked 4D inputs,
  verified exactly zero outside the corresponding ROI (1,584 files total).
- Common valid support retains all 391 left, 427 right and 818 bilateral atlas
  voxels across all runs and both residual representations.
- Thirteen of 18 fixed-dimensionality fits passed the output/grid/mask/rank checks.
- Five fits returned without `melodic_IC.nii.gz`: joint-BOLD left dim 7, joint-BOLD
  right dims 5/7, and native-FILM right dims 5/7.
- All six bilateral fixed-dimensionality fits passed.
- The fixed-fit gate stopped the driver. None of the six automatic-dimensionality
  fits was started. Component figures and between-scale/dimensionality comparisons
  were not generated because their phase follows all fits.

The recorded state is correctly **failed**, not complete. The error is a missing
MELODIC component-map output; the runtime's parenthetical "possible nonconvergence"
is not an established diagnosis. Actual MELODIC stdout/log.txt files remain on
Linux2 and must be collected before changing optimization settings or explaining
the cause. The GitHub record currently contains Python error reports, not those
underlying logs. No evidence here establishes a memory failure, registration
failure, insufficient dimensionality, or poor component quality.

Use `collect_group_ica_diagnostics.py --publish` on Linux2 to collect bounded
MELODIC logs, eigenvalue reports, completion metadata and file inventories into
a separate diagnostic record. This command changes no analysis inputs/results and
starts no fits. Preserve the 13 completed fits, all residual arrays, recovered UGR
statistics and all masked input files for subsequent diagnosis and resume.
