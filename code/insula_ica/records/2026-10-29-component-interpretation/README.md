# Interpretation for the paper: heterogeneity, not a five-region claim

## What the current analysis establishes descriptively

The pooled task analysis contains all 132 manuscript participants and all 490
accepted runs. The current figures characterize the five-component unsmoothed
ICA solutions discovered in full task+nuisance residuals. Task responses were
then independently extracted from task-preserved FEAT inputs, with the same
complete task+nuisance model standardized and fit by OLS. UGR now has common
high/low social/nonsocial slopes; Trust retains partner rec-minus-def contrasts.

For the bilateral solution, IC3 has positive Trust rec-minus-def means, strongest
for friend; IC4 has negative means for all three partners. All bilateral UGR
pointwise intervals span zero. The left joint-BOLD IC3 social-UGR interval is
just above zero, but the corresponding native-FILM discovery-basis interval spans
zero. There is not a robust demonstrated shared betrayal/unfairness component
at this stage. These are exploratory same-cohort descriptions; no formal
component x outcome or between-task interaction has been established. Component
signs are fixed by the largest absolute spatial loading, not the task results.

The two task bars are not measurements of the same construct: UGR is an offer
pmod slope; Trust is a difference between two outcome coefficients. Z-scoring
does not turn those into a directly comparable task interaction. A Trust
rec-minus-def difference does not establish absolute activation/deactivation.

## What the ICs look like

The accompanying anatomy figure uses the actual bilateral five-component maps
and a standard MNI T1 underlay. Signed component loadings overlap; the display
threshold is not a statistical or parcel boundary. Sagittal panels are at each
component's left loading peak, coronal panels at its global peak y coordinate.

| Bilateral IC | Descriptive topography | Existing background-connectivity targets |
| --- | --- | --- |
| 1 | Most anterior | Frontal operculum, inferior/middle frontal cortex |
| 2 | Mid/anterior | Frontal/central/parietal operculum, supplementary motor cortex |
| 3 | Posterior and relatively dorsal | Central/parietal operculum, Heschl's gyrus, postcentral cortex, planum temporale |
| 4 | Ventral/inferior mid-to-posterior | Planum polare, cingulate, parietal operculum and posterior parahippocampal regions |
| 5 | Central/inferior, more distributed | Opercular/temporal targets; less clean correspondence across model orders |

These names are anatomical descriptions, not canonical network assignments.
IC4's broad positive fingerprint needs particular scrutiny before calling it
DMN, valuation or a selective limbic network. Existing connectivity values are
descriptive partial regression means, not corrected significant clusters or
formal comparisons between component connectivity maps.

The posterior/sensorimotor association is broadly consistent with prior insular
connectivity parcellation, but does not itself establish an identity for this IC:
[Deen et al. 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3116731/).
The restricted ICA approach seeks overlapping network signals and should be
supported by whole-brain fingerprints and model-order sensitivity, as in
[Leech et al. 2012](https://www.jneurosci.org/content/32/1/215).

## Important interpretation of the spatial intercept

The extraction fits all centered maps plus a spatial intercept. Its map-weight
operator annihilates a spatially uniform ROI signal (verified numerically for
all six primary bases in spatial-intercept-check.json). Such a signal is carried
by the intercept, not these five pattern coefficients. The component timecourses
are not five ordinary regional-average BOLD signals. Therefore weak UGR
component bars do not establish weak insula-wide UGR activation or contradict
the original univariate finding. The next task review should also fit the
insula-mean/intercept timecourse as a reference using the same temporal model.
This reference was not retained in the current cached five-component responses.

## Other dimensionalities: spatial evidence exists, task evidence does not yet

For bilateral joint-BOLD maps:

| Five-component pattern | Three-component match | Seven-component match |
| --- | --- | --- |
| IC1 | IC3, r=0.771 | IC1, r=0.933 |
| IC2 | IC2, r=0.561 | IC2, r=0.946 |
| IC3 | IC1, r=0.784 | IC3, r=0.872; additional relation to IC4, r=0.497 |
| IC4 | No separate one-to-one match; largest r=0.476 | IC5, r=0.950 |
| IC5 | No separate one-to-one match | IC7, r=0.720 and IC6, r=0.596 |

The principal anatomical patterns have promising correspondence at seven
components; the posterior and central/inferior patterns can subdivide. Three
components merge signals that are distinguishable at five. This does not make
five the uniquely correct number. Native-FILM gives a similar bilateral picture;
left-only correspondences are less uniform. Match spatially and align signs
before interpreting task coefficients; numbers are local to each solution.

All 18 corrected unsmoothed fixed fits (3/5/7 x three ROIs x two bases) are
represented in the completed-map snapshot. The current pooled task models were
only run for the five-component bases. No final corrected unsmoothed automatic-
dimensionality record is published in this checkout, so its completion cannot
be claimed. The older smoothed stream reports 19/24 completed fits. The older
matched confounds-only/full-OLS task-removal sensitivity reports 5/6, with the
confounds-only left solution failing to converge. Those older streams do not
substitute for a matched comparison on the corrected current sm-0 inputs.

## Relation to the manuscript's anterior-insula seed

The supplied brief identifies the reviewer concern as treating left anterior
insula as homogeneous. A whole-insula decomposition addresses that concern
only if we establish its relationship to the actual paper seed/ROI. The
strongest positive Trust pattern here is relatively posterior, so it cannot
simply be relabeled anterior-insula validation.

The project wrappers refer to masks/seed-AIns.nii.gz, which is not in this
checkout. A tracked seed-AIns-clusterthresh.nii.gz exists, but its use as the
actual paper PPI seed has not been verified. Do not substitute it silently.
The exact historical seed should be read from Linux2 and overlaid on the ICA
support and maps. This is an unresolved provenance check, not a demonstrated
third error in the paper or proof of the actual seed's overlap.

## How it could help the paper

The defensible contribution would be an exploratory reviewer-response analysis:
within an independent anatomical insula mask, multiple overlapping spatial
signals have distinct background connectivity fingerprints and differing task
responses. That can replace a homogeneous-insula interpretation with a more
precise discussion. The new analysis is not independent validation of the
original sample-defined aIns-vmPFC gPPI, and currently does not establish a shared
betrayal/unfairness mechanism or psychological specificity of a component.

Before a manuscript claim, prioritize:

1. Confirm and overlay the exact original seed/ROI. Include its relation to the
   independent anatomical mask and actual component support.
2. Repeat the same standardized pooled task model for every three- and seven-
   component solution, alongside the ROI-mean/intercept reference. Use spatial
   matches fixed before inspecting task effects; retain all components.
3. Complete a matched task-preserved discovery comparison on the corrected
   sm-0 inputs: nuisance-only versus full-model task removal, holding the fitting
   convention and all other inputs fixed. Extracting raw task timecourses from
   residual-derived maps is not this discovery sensitivity.
4. Establish independent network-template correspondence and corrected group
   connectivity evidence outside the defining ROI before assigning network
   names. Mean maps alone do not demonstrate between-network differences.
5. Keep the two first-level model issues and their impact assessment on the
   original manuscript analyses separate from this exploratory ICA addition.

The bars alone should not carry a new mechanistic or network-specific conclusion.
A supplement or bounded reviewer-response result is a realistic role if the
spatial, task and connectivity sensitivities tell the same story.
