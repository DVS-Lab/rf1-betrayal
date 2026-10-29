# First review of completed unsmoothed ICA exploration

The exploration completed all 490 run-level connectivity analyses and published component and background-connectivity review figures and tables. Automatic-dimensionality ICA is a separate process; no completion of its
six fits is implied by this exploration record.

## Component anatomy and connectivity

The bilateral five-component joint-BOLD maps show anterior/posterior and
dorsal/ventral differentiation. Peak coordinates are descriptive loading peaks,
not centers of disjoint parcels:

| IC | Peak y, z (mm) | Descriptive connectivity targets outside the seed |
| --- | --- | --- |
| 1 | 22.5, -3.7 | Frontal operculum, inferior/middle frontal cortex |
| 2 | 9.0, -0.7 | Frontal, parietal and central operculum |
| 3 | -9.9, 14.1 | Central/parietal operculum and Heschl's gyrus |
| 4 | -4.5, -12.6 | Planum polare, anterior cingulate and parietal operculum |
| 5 | 0.9, -6.7 | Opercular/temporal targets; negative medial-frontal/temporal means |

These are first anatomical descriptions supported by the component PNGs and
regional coefficient table. They are not independently verified reference-network
assignments or significant connectivity clusters. Background maps pool both tasks
and do not establish task-modulated connectivity. Native-FILM and joint-BOLD
bilateral primary maps match very closely (matched spatial r 0.995–0.999), as do
right maps (0.996–0.997). Left matches are lower (0.891–0.950). These use the same
participants/inputs and are sensitivity comparisons, not independent replication.

## Task analysis withdrawn

The previous task-analysis method was rejected on 2026-10-29. Its tables, figures,
statistics and directional conclusions have been removed. IC3/IC4 task directions
must be re-estimated from network timecourses in task-preserved BOLD. Background
connectivity and component anatomy are retained.

## Coverage and publication follow-up

Common connectivity support contains 67,079 voxels, 54.55% of the union of fitted
masks. This ratio is **not** a fraction of an independent anatomical whole-brain
mask: the union can include noisy/peripheral voxels. The common mask should be
checked against an anatomical reference before interpreting missing targets.
Zeros outside it are missing coverage, not null connectivity.

A nested `code/insula_ica/.gitignore` excluded the intended group NIfTI volumes
while PNGs/tables published successfully. The rule is now narrowly relaxed for
group `primary_maps/` and `group_maps/` artifacts only. Run
`publish_ica_review_maps.py --publish` to validate and publish the existing 13
group-derived volumes without repeating any numerical analysis. Raw/subject
4D imaging remains ignored.
