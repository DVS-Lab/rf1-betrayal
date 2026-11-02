# Completed functional-cluster and smoothing review

This review reads the published `linux2-functional-cluster-20261101T200404Z-74f591` and `linux2-functional-sm5-20261101T211149Z-7aa636` records. It runs no production fit and changes no original results. All 16 planned ICAs and all 490 run-level temporal analyses in each arm completed with no recorded failure/publication errors. N=132 and all accepted runs remain. The smoothed arm reused 242 matching FEATs and ran 248 new matched FEATs. All 465 original anatomical-insula summaries were reused unchanged.

## Masks and models

The functional mean uses the original 127-voxel current cluster candidate. ICA uses 117 common fitted/finite/nonconstant voxels (92.1% of that candidate) across the unsmoothed data. The smoothed arm covers all 127 but retains exactly the same 117-voxel ICA mask for the comparison, with no new erosion or participant/run removal. Two atlas-insula voxels remain in both original and common masks.

Both ICA bases are from the existing full-task-plus-nuisance residuals: BOLD-scale residuals using existing FILM PEs and native voxelwise-prewhitened FILM residuals. Smoothing changes only the FEAT smoothing setting to 5 mm, using the same raw input, task model, nuisance matrix, corrected TRs and FSL version. New discovery maps are projected by simultaneous spatial regression, plus an intercept, onto the same unresidualized masked FEAT input to estimate task timecourses. All temporal columns/responses are z-scored and jointly fitted, with pooled high/low UGR fairness pmods and Trust C/F/S reciprocate-minus-defect contrasts.

## Main observations

The original functional-cluster mean has negative coefficients in all five plotted conditions: UGR nonsocial −0.02155, UGR social −0.02978; Trust computer −0.01739, friend −0.05882, stranger −0.02372. These correspond descriptively to greater response for unfairness and defection, with the largest observed Trust contrast for friend. The anatomical left-insula mean is close to zero in this same model: UGR +0.00507/+0.00589; Trust −0.00101/−0.00637/−0.00183 (see original tables for full precision). This is evidence that the ROIs are not interchangeable in this task model. It is not independent confirmation of task sensitivity in a cluster selected using related task effects, nor a tested difference between the ROIs/partners.

In the five-component BOLD-scale unsmoothed solution, IC3 has an opercular/IFG loading pattern (peak MNI −45.1, 22.5, 5.2), small negative Trust bars and nearly zero UGR bars. It maps closely to smoothed IC2 (r=0.935). It also maps to the three-component solution's IC3 (r=0.955), but its seven-component match is weaker (r=0.800). Thus it is a candidate Trust-weighted pattern, not established task selectivity. Its friend rec-minus-def estimate moves from −0.01088 to −0.00521 with smoothing.

Unsmoothed functional IC4/IC5 are more ventral, with spatial loading energy predominantly in the atlas Frontal Orbital Cortex (93.8% and 89.7%). They show negative nonsocial-UGR and friend-Trust point estimates. Their smoothed matches have r=0.609/0.721. Nonsocial-UGR estimates change from −0.01391/−0.01422 to −0.01259/−0.00967, while friend-Trust estimates attenuate from −0.02032/−0.02095 to −0.01020/−0.00669. This discourages treating the fine five-component task profile as invariant to smoothing. These IC indices refer to this functional ROI, not the earlier anatomical-insula IC3/IC4.

| Requested dimensions | BOLD-scale sm0/sm5 matched absolute r | Native FILM sm0/sm5 matched absolute r |
|---|---|---|
| 3 | 0.865–0.933 | 0.871–0.918 |
| 5 | 0.597–0.935 | 0.716–0.860 |
| 7 | 0.378–0.888 | 0.714–0.916 |

Coarse three-component patterns agree well across smoothing and the two residual representations (within-arm BOLD/FILM matched r=0.970–0.980 for sm0 and 0.993–0.997 for sm5). Finer decompositions are more sensitive, especially the BOLD-scale seven-component solution. This is agreement across these settings, not bootstrap/test–retest validation of distinct parcels.

Automatic model order is 65 components for both unsmoothed bases, versus 39 (BOLD-scale) and 40 (native FILM) after smoothing. These are algorithm-selected model orders within 117 voxels, not evidence for 39–65 meaningful anatomical subregions. No solution was selected for strongest task association, and every component remains in the original outputs.

## Interpretation limits

The original-ROI mean and IC bars are different response estimators. The mean is an average over 127 voxels. The spatial component model uses 117 voxels, includes an intercept and centers its maps, so spatially uniform ROI activity does not enter the component timecourses. Therefore, small component task bars cannot be read as absence of task activity, or their magnitudes as a partition of the mean effect. The exact 117-voxel mean has not been fitted in this record; that is a small useful follow-up if directly comparing common-mode and component task responses.

The discovery ICA is task-regressed background fractionation. A task-preserved/confounds-only discovery ICA of this functional ROI has not been included in these two arms. Thus this analysis does not exclude a task-driven spatial organization of the cluster. No new whole-brain background connectivity, randomise or PPI analysis was run in this functional stream. There is no new component-to-network connectivity result to interpret yet.

The defensible present characterization is a task-sensitive selected functional cluster, with a reproducible coarse residual structure and finer task profiles that depend on decomposition/smoothing. It does not establish five distinct insular mechanisms, prove a common neural code for unfairness and betrayal, or independently validate the original PPI. The exact historical seed-file identity remains unresolved.

## Figures and reproducibility

- `five_component_task_summary.png`: reference means and all five BOLD-scale components; smoothed maps reordered/sign-aligned strictly by spatial correlation. The top panel uses a different y-range from the component panels, with all axes labeled.
- `five_component_map_comparison.png`: matched signed loadings in regional orthogonal views. Same cuts within each pair; no statistical threshold, separately normalized display maxima. These are weighted overlapping patterns, not disjoint parcels.
- `review-summary.json`: source group-table hashes, every fixed solution's atlas-energy anatomy, full residual-scale map correspondence matrices and cross-smoothing summaries.
- `review.py`: regenerates the review from the published records using numpy/scipy/nibabel/matplotlib and standard FSL atlas/template files (`FSLDIR`). No participant 4D image or new Linux2 job is needed.

Error bars and any paired-difference intervals are descriptive pointwise participant CIs without multiple-comparison correction; no significance stars or independent-inference claim are made. No production default, mask, component labels or estimates were altered by this review.
