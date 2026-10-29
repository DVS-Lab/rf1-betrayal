# First review of completed unsmoothed ICA exploration

The exploration completed all 490 run-level connectivity analyses and all 60
planned primary task projections, publishing 26 PNG review figures and derived
tables. Automatic-dimensionality ICA is a separate process; no completion of its
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

## Task projections

Eight of 60 tests have Bonferroni-60 p<0.05; all eight concern Trust's stored
`rec-def` contrast. No UGR offer-pmod projection passed that combined threshold.

- Bilateral IC3 is positive in both representations (corrected p approximately
  0.00014 joint-BOLD and 0.000082 native-FILM).
- Bilateral IC4 is negative in both (0.0061 and 0.0035).
- Right IC3 is positive in both (0.00022 and 0.00018).
- Left IC3 is positive for joint-BOLD (0.0226); left IC4 is negative for native
  FILM (0.0426). Left-side task findings therefore depend more on representation.

Opposing Trust directions are interesting for follow-up, but do not establish a
component-specific psychological mechanism, a significant hemisphere/task
interaction, or validation of the original PPI. The bases were discovered in this
cohort and the tests are exploratory. Different significance decisions are not a
test of a difference. Absolute effect units also differ between task contrasts.

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
