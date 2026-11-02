# Component anatomy and paired social contrasts

This small-file review reuses the completed `linux2-functional-cluster-20261101T200404Z-74f591`
and `linux2-functional-sm5-20261101T211149Z-7aa636` records. It reads published
component maps and saved participant estimates, with no new participant 4D
reads, FEAT fits, ICAs, or timecourse-model fits. Whole-brain network results
are not present yet: those require the new Linux2 dual-regression launcher.

`component_anatomy.tsv` and `component_atlas_loading_profiles.json` cover all
32 components: sm0/sm5 × BOLD-scale/native FILM × 3/5 components. Atlas fractions
refer to squared loading energy, not exclusive parcels or voxel counts.
Signed loadings are displayed within the fixed 117-voxel discovery mask, on
the native analysis grid and with each map divided by its own absolute peak
for display only. Indices remain local to each solution.

In the unsmoothed five-component BOLD-scale solution, IC2 is weighted toward
inferior frontal pars opercularis and IC3 toward frontal operculum. IC4 and IC5
are ventral frontal orbital patterns (94% and 90% of squared loading energy).
IC1 spans a more dorsal orbital/IFG border. These are fractions of Jacob's
functional cluster, not five anatomical insula parcels. The coarse spatial
differences are visible before whole-brain connectivity is known.

New task figures show four original condition estimates plus two paired
social-specificity contrasts. Friend effects and a friend/stranger average
are retained as secondary table summaries. Negative social-minus-nonsocial
differences are compatible with greater social norm-violation sensitivity
under the original fairness / reciprocate-minus-defect directions. Every
primary paired-contrast pointwise 95% CI in the eight component solutions
includes zero; this is descriptive characterization of a task-selected ROI,
not evidence that social and nonsocial effects are equivalent. No combined
UGR-versus-Trust coefficient or connectivity interaction was fitted.

The anatomical mean reference remains the original left-insula mean; the
functional mean uses the original 127 voxels. Both mean estimates and all task
responses are reused from task-preserved raw input, so the reference bars are
the same in the sm0 and sm5 discovery figures. The original five-bar plots
remain unchanged in the source records.

`synthetic_randomise_validation.json` is a separate native-FSL interface check:
16 simulated observations, 100 permutations, explicit +/- one-sample contrasts,
3D TFCE, finite probability maps and checkpoint reuse. It is not cohort data.

Reproduce this preview with `preview_network_anatomy.py`, passing the two
completed source records, an output directory and `--fsldir`. Production
`start_wholebrain.sh` publishes the same summaries before whole-brain processing.
