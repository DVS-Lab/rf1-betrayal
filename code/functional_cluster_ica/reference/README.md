# Independent resting-network reference

`Smith2009_rsn10_analysis_grid.nii.gz` contains the ten group resting maps
distributed by FMRIB for Smith et al., PNAS 2009, DOI
[10.1073/pnas.0905267106](https://doi.org/10.1073/pnas.0905267106).

Source: [FMRIB BrainMap and resting networks](https://www.fmrib.ox.ac.uk/datasets/brainmap%2Brsns/),
file `PNAS_Smith09_rsn10.nii.gz`. The source and derived SHA256 hashes, dimensions,
affine and labels are pinned in `Smith2009_rsn10_metadata.json`.

The independent templates alone were linearly resampled with
`nibabel.processing.resample_from_to` (order 1) to the existing functional ICA
analysis grid (57 × 70 × 54). No participant data were resampled. Packaging the
derived reference avoids an external download during the Linux2 run.

Correspondence is a descriptive spatial Pearson correlation of each signed
whole-brain group map with every template, using the common brain outside the
original 127-voxel defining ROI plus a 10-mm guard band. It does not establish
exclusive network membership or test a new task effect. Network 9/10 labels
remain broad frontoparietal labels rather than inferred hemisphere names.
