# Task dimensionality review

Reviewed the results published in commit 8c7062f3. All 490 runs completed, every
group estimate includes 132 participants, and the 150 reused five-component
group estimates differ from their source by exactly zero. There are 450
component bar estimates and 15 unique ROI-mean estimates; no recorded run errors
or coefficient-magnitude review flags. No exclusions were added.

## Bilateral five-component IC3

The positive Trust response, especially friend reciprocate minus defect, is
preserved across all three dimensions. Component numbers must follow spatial
correspondence rather than matching labels.

| Discovery basis | Dimension / IC | Spatial r with five-component IC3 | Friend rec-minus-def, standardized contrast [pointwise 95% CI] |
|---|---|---:|---:|
| joint_bold | 3 / IC1 | 0.784 | 0.0201 [0.0086, 0.0317] |
| joint_bold | 5 / IC3 | 1.000 | 0.0283 [0.0166, 0.0401] |
| joint_bold | 7 / IC3 | 0.872 | 0.0274 [0.0155, 0.0393] |
| native_film | 3 / IC3 | 0.870 | 0.0253 [0.0145, 0.0361] |
| native_film | 5 / IC3 | 1.000 | 0.0288 [0.0170, 0.0406] |
| native_film | 7 / IC3 | 0.889 | 0.0274 [0.0155, 0.0392] |

The positive computer contrast is clearer at seven components; the stranger
contrast is less consistently separated from zero. Neither bilateral UGR slope
has a pointwise CI excluding zero in these matched components. This supports a
stable friend contrast more directly than general partner-invariant or
cross-task selectivity. These solutions use the same sample, not independent
replications.

## Bilateral five-component IC4

At seven components this is IC5, with spatial r=0.950 (joint_bold) and 0.933
(native_film). All three Trust group contrasts remain negative. For joint_bold,
computer/friend/stranger estimates change from -0.0106/-0.0150/-0.0099 at five
components to -0.0098/-0.0109/-0.0085 at seven. The seven-component stranger CI
includes zero, whereas computer and friend pointwise CIs remain below zero.
The native_film result has the same qualitative pattern. Thus direction and
spatial identity are retained at seven; the strength of evidence is not
identical for every partner.

At three components IC4 has no separate one-to-one match. Its largest positive
spatial correlation is with joint_bold IC2 (r=0.476), which is assigned to
five-component IC2 (r=0.561), or native_film IC1 (r=0.403), assigned to
five-component IC2 (r=0.794). These coarse components also have negative Trust
group estimates, but mix multiple five-component spatial patterns. They do not
establish preservation of a distinct IC4 subregion at dimension three.

## Mean reference and interpretation

All 15 left/right/bilateral insula-mean pointwise CIs include zero. For bilateral
insula, computer/friend/stranger contrasts are -0.0023/-0.0033/+0.0004, and
nonsocial/social UGR slopes are +0.0063/+0.0062. The means are fitted once per
ROI and reused as the same reference across discovery bases and dimensions.

The weak mean alongside opposing component responses is consistent with
functional heterogeneity within the insula. It does not prove exact cancellation
or a formal between-component difference. Spatial component regression includes
an intercept, so components describe signed spatial variation around the ROI
mean; they are not disjoint parcel-average activation measures. Positive and
negative task contrasts refer to the saved oriented component timecourses.

These results support retaining the IC3 observation and the IC4 observation at
five/seven components as descriptive findings. They do not establish equivalence
across dimensions, corrected significance, a partner interaction, or a
whole-brain connectivity/gPPI result. CIs are pointwise and unadjusted; no new
formal test was run for this review. Component numbering and task patterns in
the separate unilateral decompositions must be assessed within those bases.

## Figures and source tables

- [joint_bold, 3 components, all ROIs and means](figures/joint_bold_dim-3_task_bars_with_mean.png)
- [joint_bold, 7 components, all ROIs and means](figures/joint_bold_dim-7_task_bars_with_mean.png)
- [native_film, 3 components, all ROIs and means](figures/native_film_dim-3_task_bars_with_mean.png)
- [native_film, 7 components, all ROIs and means](figures/native_film_dim-7_task_bars_with_mean.png)
- [Insula-mean responses](figures/insula_mean_task_bars.png)
- [All group estimates](group_task_effects.tsv)
- [Spatial correlations and assignments](spatial_correspondences.json)
- [Spatially matched task comparisons](aligned_task_comparisons.tsv)
- [Five-component reuse verification](reused-five-component-verification.json)
