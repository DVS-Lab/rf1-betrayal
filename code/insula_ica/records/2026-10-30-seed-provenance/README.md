# Historical seed candidates and original group-map reassessment

The tracked `seed-AIns-clusterthresh.nii.gz` is a task-intersection cluster,
but its identity as Jacob's final manuscript seed is unestablished. Its name
and inclusion in a saved overlap do not establish the intended paper effect.
The earlier plot's title has been clarified accordingly.

Git history recovered the other AIns/intersection candidates removed during the
June 11, 2026 cleanup (commit 28bc0df6). The last pre-cleanup snapshot is
f70193c5, June 9, 2026. Recovered bytes are kept in this record's `primary_maps/`,
not restored into production `masks/`. `historical-mask-inventory.json` records
original paths, Git blob IDs, hashes, geometry and atlas overlap.
`all-history-mask-versions.json` also checks all five mask-changing commits
reachable in this checkout: there are nine distinct AIns/intersection blobs;
the earlier snapshots provide no additional distinct mask candidate.

| Candidate | Nonzero voxels | Approximate centroid (MNI mm) | HO25 insula overlap |
|---|---:|---|---:|
| Current seed-AIns-clusterthresh | 127 | -43.5, 22.0, -4.7 | 1.6% |
| Archived seed-AIns | 142 | -44.3, 22.3, -4.3 | 1.4% |
| Deleted seed-AIns-vt | 1 | -34.3, 19.8, -12.6 | 0% |
| Archived seed-aIns (lowercase a) | 8 | 7.5, 12.7, -26.0 | 0% |

The 8-voxel mask spans more than one region; its centroid alone is not a parcel
location. None of these candidates provides a convincing anatomical-insula
replacement. All examined historical candidates have matching qform and sform
matrices with MNI codes; that excludes a qform-versus-sform display disagreement,
but does not prove correct original geometry or effect provenance.

Historical scripts also differ: the January 20 intersection script referenced
N128 voxel-thresholded groups, while the current saved script names N132
ones-model cluster masks. The current UGR L3 wrapper selects a full-covariate
model. Recover the actual saved group designs, effects/signs, thresholds and
source masks before deciding which branch produced the paper's seed. Do not
select a candidate merely because it looks more insular.

## UGR repair scope

The established run-TR error in sub-10817 and the unidentifiable social-high
pmod in sub-11031 both need original first-level evaluation before L2/L3
propagation. The TR-repair record already contains a corrected sm-5 candidate
and a matched same-regenerated-BOLD timing control. The COPE11 map comparison
has r=-0.1253 and relative RMS difference=1.4014, demonstrating a large run-level
timing effect in that control. This does not establish the group-level impact
or isolate new preprocessing from the original historical images.

The pooled model that stabilized ICA task plots is a new common-high/low-slope
model, not an automatic repair to the paper's voxelwise estimates or high/low
effects. Audit original sm-5 model contrast estimability, then document the
original-model repair separately from the pooled-model sensitivity. Preserve
variance estimates and propagate changes through participant and group models;
an L3-only rerun would reuse affected estimates. If the seed changes, dependent
PPI regressors and fits must be regenerated as well.

The prepared read-only collector is documented at
[`code/ugr_reanalysis/README.md`](../../../ugr_reanalysis/README.md).
It examines the manuscript's actual project directories and exports small group
maps/designs and live seeds, plus all original UGR design contrast diagnostics.
No rerun has been launched and no production seed has been replaced. Group
effect, sign, threshold and input reviews remain pending Linux2 output.

Recovered file example:

```bash
git show f70193c5:masks/seed-AIns-vt.nii.gz > /tmp/seed-AIns-vt-from-git.nii.gz
```
