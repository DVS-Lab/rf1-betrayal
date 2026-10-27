# Trust BOLD preprocessing provenance

No preprocessing, model fits, or production jobs were launched.

The original nuisance fit is already preserved; this check does not require standalone TEDANA files.

Historical dataset/subject metadata can have been overwritten. Config comparisons below are candidate evidence, not proof that a config generated the old Trust images.

## Version and execution evidence

Rebuild container record: `{"schema_version": 2, "path": "/ZPOOL/data/tools/fmriprep-24.1.1.simg", "sha256": "6261734ae209de4b103fa051d4c4fe8ea914cb4a4623014fd208454eace0c67f", "version": "24.1.1"}`

### sub-10369

historical: 0 execution configs; 0 methods/summary reports.
Subject HTML version strings: `[]` (subject-level evidence only).
rebuilt: 1 execution configs; 6 methods/summary reports.
Subject HTML version strings: `['24.1.1']` (subject-level evidence only).
- Dataset-level version evidence: `{"GeneratedBy": [{"Name": "fMRIPrep", "Version": "24.1.1", "CodeURL": "https://github.com/nipreps/fmriprep/archive/24.1.1.tar.gz"}]}`
- `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/sub-10369/log/20261026-205055_4addb3a4-5f3b-4526-a716-70013f3ab0a6/fmriprep.toml` — trust_named; environment: `{"cpu_count": 128, "exec_env": "singularity", "free_mem": 109.8, "overcommit_policy": "heuristic", "overcommit_limit": "50%", "nipype_version": "1.8.6", "templateflow_version": "24.2.0", "version": "24.1.1"}`

Report `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/sub-10369/figures/sub-10369_ses-01_desc-about_T1w.html`:

fMRIPrep version: 24.1.1 
 fMRIPrep command: /opt/conda/envs/fmriprep/bin/fmriprep /bids /rebuild/fmriprep-24 participant --participant-label 10369 --task-id trust --bids-filter-file /rebuild/filters/sub-10369.json --stop-on-first-crash --skip-bids-validation --nthreads 14 --mem-mb 28000 --me-output-echos --output-spaces T1w MNI152NLin6Asym --fs-no-reconall --fs-license-file /opts/fs_license.txt -w /scratch 
 Date preprocessed: 2026-10-26 20:51:29 -0400


Report `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/sub-10369/figures/sub-10369_ses-01_desc-summary_T1w.html`:

Subject ID: 10369 
 Structural images: 1 T1-weighted 
 Functional series: 4 
 
 Task: trust (4 runs) 
 
 Standard output spaces: MNI152NLin6Asym, MNI152NLin2009cAsym 
 Non-standard output spaces: T1w 
 FreeSurfer reconstruction: Not run


Report `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/sub-10369/figures/sub-10369_ses-01_task-trust_run-1_part-mag_desc-summary_bold.html`:

Summary 
 
 Original orientation: LAS 
 Repetition time (TR): 1.61s 
 Phase-encoding (PE) direction: Anterior-Posterior 
 Multi-echo EPI sequence: 4 echoes. 
 Slice timing correction: Applied 
 Susceptibility distortion correction: FMB (fieldmap-based) - directly measured B0 map 
 Registration: FSL flirt with boundary-based registration (BBR) metric - 6 dof 
 Non-steady-state volumes: 2


Report `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/logs/CITATION.md`:

Results included in this manuscript come from preprocessing
performed using *fMRIPrep* 24.1.1
(@fmriprep1; @fmriprep2; RRID:SCR_016216),
which is based on *Nipype* 1.8.6
(@nipype1; @nipype2; RRID:SCR_002502).



Preprocessing of B 0 inhomogeneity mappings

: A total of 2 fieldmaps were found available within the input
BIDS structure for this particular subject.
A *B 0 * nonuniformity map (or *fieldmap*) was directly measured with
an MRI scheme designed with that purpose such as SEI (Spiral-Echo Imaging).

Anatomical data preprocessing

: A total of 1 T1-weighted (T1w) images were found within the input
BIDS dataset. The T1w image was corrected for intensity
non-uniformity (INU) with `N4BiasFieldCorrection` [@n4], distributed with ANTs 2.5.3
[@ants, RRID:SCR_004757], and used as T1w-reference throughout the workflow.
The T1w-reference was then skull-stripped with a *Nipype* implementation of
the `antsBrainExtraction.sh` workflow (from ANTs), using OASIS30ANTs
as target template.
Brain tissue segmentation of cerebrospinal fluid (CSF),
white-matter (WM) and gray-matter (GM) was performed on
the brain-extracted T1w using `fast` [FSL (version unknown), RRID:SCR_002823, @fsl_fast].
Volume-based spatial normalization to two standard spaces (MNI152NLin6Asym, MNI152NLin2009cAsym) was performed through
nonlinear registration with `antsRegistration` (ANTs 2.5.3),
using brain-extracted versions of both T1w reference and the T1w template.
The following templates were were selected for spatial normalization
and accessed with *TemplateFlow* [24.2.0, @templateflow]:
*FSL's MNI ICBM 152 non-linear 6th Generation Asymmetric Average Brain Stereotaxic Registration Model* [@mni152nlin6asym, RRID:SCR_002823; TemplateFlow ID: MNI152NLin6Asym], *ICBM 152 Nonlinear Asymmetrical template version 2009c* [@mni152nlin2009casym, RRID:SCR_008796; TemplateFlow ID: MNI152NLin2009cAsym].

Functional data preprocessing

: For each of the 2 BOLD runs found per subject (across all
tasks and sessions), the following preprocessing was performed.
First, a reference volume was generated from the shortest echo of the BOLD run,
using a custom methodology of *fMRIPrep*, for use in head motion correction.
Head-motion parameters with respect to the BOLD reference
(transformation matrices, and six corresponding rotation and translation
parameters) are estimated before any spatiotemporal filtering using
`mcflirt` [FSL , @mcflirt].
The estimated *fieldmap* was then aligned with rigid-registration to the target
EPI (echo-planar imaging) reference run.
The field coefficients were mapped on to the reference EPI using the transform.
The BOLD reference was then co-registered to the T1w reference using
`mri_coreg` (FreeSurfer) followed by `flirt` [FSL , @flirt]
with the boundary-based registration [@bbr] cost-function.
Co-registration was configured with six degrees of freedom.
Several confounding time-series were calculated based on the
*preprocessed BOLD*: framewise displacement (FD), DVARS and
three region-wise global signals.
FD was computed using two formulations following Power (absolute sum of
relative motions, @power_fd_dvars) and Jenkinson (relative root mean square
displacement between affines, @mcflirt).
FD and DVARS are calculated for each functional run, both using their
implementations in *Nipype* [following the definitions by @power_fd_dvars].
The three global signals are extracted within the CSF, the WM, and
the whole-brain masks.
Additionally, a set of physiological regressors were extracted to
allow for component-based noise correction [*CompCor*, @compcor].
Principal components are estimated after high-pass filtering the
*preprocessed BOLD* time-series (using a discrete cosine filter with
128s cut-off) for the two *CompCor* variants: temporal (tCompCor)
and anatomical (aCompCor).
tCompCor components are then calculated from the top 2% variable
voxels within the brain mask.
For aCompCor, three probabilistic masks (CSF, WM and combined CSF+WM)
are generated in anatomical space.
The implementation differs from that of Behzadi et al. in that instead
of eroding the masks by 2 pixels on BOLD space, a mask of pixels that
likely contain a volume fraction of GM is subtracted from the aCompCor masks.
This mask is obtained by thresholding the corresponding partial volume map at 0.05, and it ensures components are not extracted
from voxels containing a minimal fraction of GM.
Finally, these masks are resampled into BOLD space and binarized by
thresholding at 0.99 (as in the original implementation).
Components are also calculated separately within the WM and CSF masks.
For each CompCor decomposition, the *k* components with the largest singular
values are retained, such that the retained components' time series are
sufficient to explain 50 percent of variance across the nuisance mask (CSF,
WM, combined, or temporal). The remaining components are dropped from
consideration.
The head-motion estimates calculated in the correction step were also
placed withi


Report `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/logs/CITATION.html`:

"fMRIPrep citation boilerplate" 
 
 html {
 color: #1a1a1a;
 background-color: #fdfdfd;
 }
 body {
 margin: 0 auto;
 max-width: 36em;
 padding-left: 50px;
 padding-right: 50px;
 padding-top: 50px;
 padding-bottom: 50px;
 hyphens: auto;
 overflow-wrap: break-word;
 text-rendering: optimizeLegibility;
 font-kerning: normal;
 }
 @media (max-width: 600px) {
 body {
 font-size: 0.9em;
 padding: 12px;
 }
 h1 {
 font-size: 1.8em;
 }
 }
 @media print {
 html {
 background-color: white;
 }
 body {
 background-color: transparent;
 color: black;
 font-size: 12pt;
 }
 p, h2, h3 {
 orphans: 3;
 widows: 3;
 }
 h2, h3, h4 {
 page-break-after: avoid;
 }
 }
 p {
 margin: 1em 0;
 }
 a {
 color: #1a1a1a;
 }
 a:visited {
 color: #1a1a1a;
 }
 img {
 max-width: 100%;
 }
 svg {
 height: auto;
 max-width: 100%;
 }
 h1, h2, h3, h4, h5, h6 {
 margin-top: 1.4em;
 }
 h5, h6 {
 font-size: 1em;
 font-style: italic;
 }
 h6 {
 font-weight: normal;
 }
 ol, ul {
 padding-left: 1.7em;
 margin-top: 1em;
 }
 li > ol, li > ul {
 margin-top: 0;
 }
 blockquote {
 margin: 1em 0 1em 1.7em;
 padding-left: 1em;
 border-left: 2px solid #e6e6e6;
 color: #606060;
 }
 code {
 font-family: Menlo, Monaco, Consolas, 'Lucida Console', monospace;
 font-size: 85%;
 margin: 0;
 hyphens: manual;
 }
 pre {
 margin: 1em 0;
 overflow: auto;
 }
 pre code {
 padding: 0;
 overflow: visible;
 overflow-wrap: normal;
 }
 .sourceCode {
 background-color: transparent;
 overflow: visible;
 }
 hr {
 background-color: #1a1a1a;
 border: none;
 height: 1px;
 margin: 1em 0;
 }
 table {
 margin: 1em 0;
 border-collapse: collapse;
 width: 100%;
 overflow-x: auto;
 display: block;
 font-variant-numeric: lining-nums tabular-nums;
 }
 table caption {
 margin-bottom: 0.75em;
 }
 tbody {
 margin-top: 0.5em;
 border-top: 1px solid #1a1a1a;
 border-bottom: 1px solid #1a1a1a;
 }
 th {
 border-top: 1px solid #1a1a1a;
 padding: 0.25em 0.5em 0.25em 0.5em;
 }
 td {
 padding: 0.125em 0.5em 0.25em 0.5em;
 }
 header {
 margin-bottom: 4em;
 text-align: center;
 }
 #TOC li {
 list-style: none;
 }
 #TOC ul {
 padding-left: 1.3em;
 }
 #TOC > ul {
 padding-left: 0;
 }
 #TOC a:not(:hover) {
 text-decoration: none;
 }
 code{white-space: pre-wrap;}
 span.smallcaps{font-variant: small-caps;}
 div.columns{display: flex; gap: min(4vw, 1.5em);}
 div.column{flex: auto; overflow-x: auto;}
 div.hanging-indent{margin-left: 1.5em; text-indent: -1.5em;}
 /* The extra [class] is a hack that increases specificity enough to
 override a similar rule in reveal.js */
 ul.task-list[class]{list-style: none;}
 ul.task-list li input[type="checkbox"] {
 font-size: inherit;
 width: 0.8em;
 margin: 0 0.8em 0.2em -1.6em;
 vertical-align: middle;
 }
 .display.math{display: block; text-align: center; margin: 0.5rem auto;}
 /* CSS for citations */
 div.csl-bib-body { }
 div.csl-entry {
 clear: both;
 margin-bottom: 0em;
 }
 .hanging-indent div.csl-entry {
 margin-left:2em;
 text-indent:-2em;
 }
 div.csl-left-margin {
 min-width:2em;
 float:left;
 }
 div.csl-right-inline {
 margin-left:2em;
 padding-left:1em;
 }
 div.csl-indent {
 margin-left: 2em;
 } 
 
 
 Results included in this manuscript come from preprocessing performed
using fMRIPrep 24.1.1 ( Esteban et al. (2019) ; Esteban et al. (2018) ;
RRID:SCR_016216), which is based on Nipype 1.8.6 ( K. Gorgolewski et al.
(2011) ; K. J.
Gorgolewski et al. (2018) ; RRID:SCR_002502). 
 
 Preprocessing of B 0 inhomogeneity mappings 
 
 A total of 2 fieldmaps were found available within the input BIDS
structure for this particular subject. A B 0 
nonuniformity map (or fieldmap ) was directly measured with an
MRI scheme designed with that purpose such as SEI (Spiral-Echo
Imaging). 
 
 Anatomical data preprocessing 
 
 A total of 1 T1-weighted (T1w) images were found within the input
BIDS dataset. The T1w image was corrected for intensity non-uniformity
(INU) with N4BiasFieldCorrection (Tustison et al. 2010) , distributed with ANTs
2.5.3 (Avants et al. 2008,
RRID:SCR_004757) , and used as T1w-reference throughout the
workflow. The T1w-reference was then skull-stripped with a
 Nipype implementation of the
 antsBrainExtraction.sh workflow (from ANTs), using
OASIS30ANTs as target template. Brain tissue segmentation of
cerebrospinal fluid (CSF), white-matter (WM) and gray-matter (GM) was
performed on the brain-extracted T1w using fast (FSL (version unknown),
RRID:SCR_002823, Zhang, Brady, and Smith 2001) . Volume-based
spatial normalization to two standard spaces (MNI152NLin6Asym,
MNI152NLin2009cAsym) was performed through nonlinear registration with
 antsRegistration (ANTs 2.5.3), using brain-extracted
versions of both T1w reference and the T1w template. The following
templates were were selected for spatial normalization and accessed with
 TemplateFlow (24.2.0, Ciric et al. 2022) : FSL’s
MNI ICBM 152 non-linear 6th Generation Asymmetric Average Brain
Stereotaxic Registration Model [ Evans et al. (2012) ,
RRID:SCR_002823; TemplateFlow ID: MNI152NLin6Asym], ICBM 152
Nonlinear Asymmetrical templ


Report `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/logs/CITATION.bib`:

@article{fmriprep1,
 author = {Esteban, Oscar and Markiewicz, Christopher and Blair, Ross W and Moodie, Craig and Isik, Ayse Ilkay and Erramuzpe Aliaga, Asier and Kent, James and Goncalves, Mathias and DuPre, Elizabeth and Snyder, Madeleine and Oya, Hiroyuki and Ghosh, Satrajit and Wright, Jessey and Durnez, Joke and Poldrack, Russell and Gorgolewski, Krzysztof Jacek},
 title = {{fMRIPrep}: a robust preprocessing pipeline for functional {MRI}},
 volume = {16},
 pages = {111--116},
 year = {2019},
 doi = {10.1038/s41592-018-0235-4},
 journal = {Nature Methods}
}

@article{fmriprep2,
 author = {Esteban, Oscar and Blair, Ross and Markiewicz, Christopher J. and Berleant, Shoshana L. and Moodie, Craig and Ma, Feilong and Isik, Ayse Ilkay and Erramuzpe, Asier and Kent, James D. andGoncalves, Mathias and DuPre, Elizabeth and Sitek, Kevin R. and Gomez, Daniel E. P. and Lurie, Daniel J. and Ye, Zhifang and Poldrack, Russell A. and Gorgolewski, Krzysztof J.},
 title = {fMRIPrep 24.1.1},
 year = 2018,
 doi = {10.5281/zenodo.852659},
 publisher = {Zenodo},
 journal = {Software}
}

@article{templateflow,
 author = {Ciric, R. and Thompson, William H. and Lorenz, R. and Goncalves, M. and MacNicol, E. and Markiewicz, C. J. and Halchenko, Y. O. and Ghosh, S. S. and Gorgolewski, K. J. and Poldrack, R. A. and Esteban, O.},
 title = {{TemplateFlow}: {FAIR}-sharing of multi-scale, multi-species brain models},
 volume = {19},
 pages = {1568--1571},
 year = {2022},
 doi = {10.1038/s41592-022-01681-2},
 journal = {Nature Methods}
}

@article{nipype1,
 author = {Gorgolewski, K. and Burns, C. D. and Madison, C. and Clark, D. and Halchenko, Y. O. and Waskom, M. L. and Ghosh, S.},
 doi = {10.3389/fninf.2011.00013},
 journal = {Frontiers in Neuroinformatics},
 pages = 13,
 shorttitle = {Nipype},
 title = {Nipype: a flexible, lightweight and extensible neuroimaging data processing framework in Python},
 volume = 5,
 year = 2011
}

@article{nipype2,
 author = {Gorgolewski, Krzysztof J. and Esteban, Oscar and Markiewicz, Christopher J. and Ziegler, Erik and Ellis, David Gage and Notter, Michael Philipp and Jarecka, Dorota and Johnson, Hans and Burns, Christopher and Manhães-Savio, Alexandre and Hamalainen, Carlo and Yvernault, Benjamin and Salo, Taylor and Jordan, Kesshi and Goncalves, Mathias and Waskom, Michael and Clark, Daniel and Wong, Jason and Loney, Fred and Modat, Marc and Dewey, Blake E and Madison, Cindee and Visconti di Oleggio Castello, Matteo and Clark, Michael G. and Dayan, Michael and Clark, Dav and Keshavan, Anisha and Pinsard, Basile and Gramfort, Alexandre and Berleant, Shoshana and Nielson, Dylan M. and Bougacha, Salma and Varoquaux, Gael and Cipollini, Ben and Markello, Ross and Rokem, Ariel and Moloney, Brendan and Halchenko, Yaroslav O. and Wassermann , Demian and Hanke, Michael and Horea, Christian and Kaczmarzyk, Jakub and Gilles de Hollander and DuPre, Elizabeth and Gillman, Ashley and Mordom, David and Buchanan, Colin and Tungaraza, Rosalia and Pauli, Wolfgang M. and Iqbal, Shariq and Sikka, Sharad and Mancini, Matteo and Schwartz, Yannick and Malone, Ian B. and Dubois, Mathieu and Frohlich, Caroline and Welch, David and Forbes, Jessica and Kent, James and Watanabe, Aimi and Cumba, Chad and Huntenburg, Julia M. and Kastman, Erik and Nichols, B. Nolan and Eshaghi, Arman and Ginsburg, Daniel and Schaefer, Alexander and Acland, Benjamin and Giavasis, Steven and Kleesiek, Jens and Erickson, Drew and Küttner, René and Haselgrove, Christian and Correa, Carlos and Ghayoor, Ali and Liem, Franz and Millman, Jarrod and Haehn, Daniel and Lai, Jeff and Zhou, Dale and Blair, Ross and Glatard, Tristan and Renfro, Mandy and Liu, Siqi and Kahn, Ari E. and Pérez-García, Fernando and Triplett, William and Lampe, Leonie and Stadler, Jörg and Kong, Xiang-Zhen and Hallquist, Michael and Chetverikov, Andrey and Salvatore, John and Park, Anne and Poldrack, Russell and Craddock, R. Cameron and Inati, Souheil and Hinds, Oliver and Cooper, Gavin and Perkins, L. Nathan and Marina, Ana and Mattfeld, Aaron and Noel, Maxime and Lukas Snoek and Matsubara, K and Cheung, Brian and Rothmei, Simon and Urchs, Sebastian and Durnez, Joke and Mertz, Fred and Geisler, Daniel and Floren, Andrew and Gerhard, Stephan and Sharp, Paul and Molina-Romero, Miguel and Weinstein, Alejandro and Broderick, William and Saase, Victor and Andberg, Sami Kristian and Harms, Robbert and Schlamp, Kai and Arias, Jaime and Papadopoulos Orfanos, Dimitri and Tarbert, Claire and Tambini, Arielle and De La Vega, Alejandro and Nickson, Thomas and Brett, Matthew and Falkiewicz, Marcel and Podranski, Kornelius and Linkersdörfer, Janosch and Flandin, Guillaume and Ort, Eduard and Shachnev, Dmitry and McNamee, Daniel and Davison, Andrew and Varada, Jan and Schwabacher, Isaac and Pellman, John and Perez-Guevara, Martin and Khanuja, Ranjeet and Pannetier, Nicolas and McDermottroe, Conor and Ghosh, Satrajit},
 title = {Nipype},
 year = 2018,
 doi = {10.5281/zenodo.596855},
 publisher


No historical/rebuilt candidate config pair is available. Exact historical execution differences remain unknown.
### sub-10402

historical: 0 execution configs; 0 methods/summary reports.
Subject HTML version strings: `[]` (subject-level evidence only).
rebuilt: 1 execution configs; 7 methods/summary reports.
Subject HTML version strings: `['24.1.1']` (subject-level evidence only).
- Dataset-level version evidence: `{"GeneratedBy": [{"Name": "fMRIPrep", "Version": "24.1.1", "CodeURL": "https://github.com/nipreps/fmriprep/archive/24.1.1.tar.gz"}]}`
- `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/sub-10402/log/20261026-205055_82673113-2774-4213-bf22-c08536e8c6e2/fmriprep.toml` — trust_named; environment: `{"cpu_count": 128, "exec_env": "singularity", "free_mem": 109.7, "overcommit_policy": "heuristic", "overcommit_limit": "50%", "nipype_version": "1.8.6", "templateflow_version": "24.2.0", "version": "24.1.1"}`

Report `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/sub-10402/figures/sub-10402_ses-01_desc-about_T1w.html`:

fMRIPrep version: 24.1.1 
 fMRIPrep command: /opt/conda/envs/fmriprep/bin/fmriprep /bids /rebuild/fmriprep-24 participant --participant-label 10402 --task-id trust --bids-filter-file /rebuild/filters/sub-10402.json --stop-on-first-crash --skip-bids-validation --nthreads 14 --mem-mb 28000 --me-output-echos --output-spaces T1w MNI152NLin6Asym --fs-no-reconall --fs-license-file /opts/fs_license.txt -w /scratch 
 Date preprocessed: 2026-10-26 20:51:32 -0400


Report `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/sub-10402/figures/sub-10402_ses-01_desc-summary_T1w.html`:

Subject ID: 10402 
 Structural images: 1 T1-weighted 
 Functional series: 2 
 
 Task: trust (2 runs) 
 
 Standard output spaces: MNI152NLin6Asym, MNI152NLin2009cAsym 
 Non-standard output spaces: T1w 
 FreeSurfer reconstruction: Not run


Report `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/sub-10402/figures/sub-10402_ses-01_task-trust_run-1_part-mag_desc-summary_bold.html`:

Summary 
 
 Original orientation: LAS 
 Repetition time (TR): 1.61s 
 Phase-encoding (PE) direction: Anterior-Posterior 
 Multi-echo EPI sequence: 4 echoes. 
 Slice timing correction: Applied 
 Susceptibility distortion correction: FMB (fieldmap-based) - directly measured B0 map 
 Registration: FSL flirt with boundary-based registration (BBR) metric - 6 dof 
 Non-steady-state volumes: 1


Report `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/sub-10402/figures/sub-10402_ses-01_task-trust_run-2_part-mag_desc-summary_bold.html`:

Summary 
 
 Original orientation: LAS 
 Repetition time (TR): 1.61s 
 Phase-encoding (PE) direction: Anterior-Posterior 
 Multi-echo EPI sequence: 4 echoes. 
 Slice timing correction: Applied 
 Susceptibility distortion correction: FMB (fieldmap-based) - directly measured B0 map 
 Registration: FSL flirt with boundary-based registration (BBR) metric - 6 dof 
 Non-steady-state volumes: 2


Report `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/logs/CITATION.md`:

Results included in this manuscript come from preprocessing
performed using *fMRIPrep* 24.1.1
(@fmriprep1; @fmriprep2; RRID:SCR_016216),
which is based on *Nipype* 1.8.6
(@nipype1; @nipype2; RRID:SCR_002502).



Preprocessing of B 0 inhomogeneity mappings

: A total of 2 fieldmaps were found available within the input
BIDS structure for this particular subject.
A *B 0 * nonuniformity map (or *fieldmap*) was directly measured with
an MRI scheme designed with that purpose such as SEI (Spiral-Echo Imaging).

Anatomical data preprocessing

: A total of 1 T1-weighted (T1w) images were found within the input
BIDS dataset. The T1w image was corrected for intensity
non-uniformity (INU) with `N4BiasFieldCorrection` [@n4], distributed with ANTs 2.5.3
[@ants, RRID:SCR_004757], and used as T1w-reference throughout the workflow.
The T1w-reference was then skull-stripped with a *Nipype* implementation of
the `antsBrainExtraction.sh` workflow (from ANTs), using OASIS30ANTs
as target template.
Brain tissue segmentation of cerebrospinal fluid (CSF),
white-matter (WM) and gray-matter (GM) was performed on
the brain-extracted T1w using `fast` [FSL (version unknown), RRID:SCR_002823, @fsl_fast].
Volume-based spatial normalization to two standard spaces (MNI152NLin6Asym, MNI152NLin2009cAsym) was performed through
nonlinear registration with `antsRegistration` (ANTs 2.5.3),
using brain-extracted versions of both T1w reference and the T1w template.
The following templates were were selected for spatial normalization
and accessed with *TemplateFlow* [24.2.0, @templateflow]:
*FSL's MNI ICBM 152 non-linear 6th Generation Asymmetric Average Brain Stereotaxic Registration Model* [@mni152nlin6asym, RRID:SCR_002823; TemplateFlow ID: MNI152NLin6Asym], *ICBM 152 Nonlinear Asymmetrical template version 2009c* [@mni152nlin2009casym, RRID:SCR_008796; TemplateFlow ID: MNI152NLin2009cAsym].

Functional data preprocessing

: For each of the 2 BOLD runs found per subject (across all
tasks and sessions), the following preprocessing was performed.
First, a reference volume was generated from the shortest echo of the BOLD run,
using a custom methodology of *fMRIPrep*, for use in head motion correction.
Head-motion parameters with respect to the BOLD reference
(transformation matrices, and six corresponding rotation and translation
parameters) are estimated before any spatiotemporal filtering using
`mcflirt` [FSL , @mcflirt].
The estimated *fieldmap* was then aligned with rigid-registration to the target
EPI (echo-planar imaging) reference run.
The field coefficients were mapped on to the reference EPI using the transform.
The BOLD reference was then co-registered to the T1w reference using
`mri_coreg` (FreeSurfer) followed by `flirt` [FSL , @flirt]
with the boundary-based registration [@bbr] cost-function.
Co-registration was configured with six degrees of freedom.
Several confounding time-series were calculated based on the
*preprocessed BOLD*: framewise displacement (FD), DVARS and
three region-wise global signals.
FD was computed using two formulations following Power (absolute sum of
relative motions, @power_fd_dvars) and Jenkinson (relative root mean square
displacement between affines, @mcflirt).
FD and DVARS are calculated for each functional run, both using their
implementations in *Nipype* [following the definitions by @power_fd_dvars].
The three global signals are extracted within the CSF, the WM, and
the whole-brain masks.
Additionally, a set of physiological regressors were extracted to
allow for component-based noise correction [*CompCor*, @compcor].
Principal components are estimated after high-pass filtering the
*preprocessed BOLD* time-series (using a discrete cosine filter with
128s cut-off) for the two *CompCor* variants: temporal (tCompCor)
and anatomical (aCompCor).
tCompCor components are then calculated from the top 2% variable
voxels within the brain mask.
For aCompCor, three probabilistic masks (CSF, WM and combined CSF+WM)
are generated in anatomical space.
The implementation differs from that of Behzadi et al. in that instead
of eroding the masks by 2 pixels on BOLD space, a mask of pixels that
likely contain a volume fraction of GM is subtracted from the aCompCor masks.
This mask is obtained by thresholding the corresponding partial volume map at 0.05, and it ensures components are not extracted
from voxels containing a minimal fraction of GM.
Finally, these masks are resampled into BOLD space and binarized by
thresholding at 0.99 (as in the original implementation).
Components are also calculated separately within the WM and CSF masks.
For each CompCor decomposition, the *k* components with the largest singular
values are retained, such that the retained components' time series are
sufficient to explain 50 percent of variance across the nuisance mask (CSF,
WM, combined, or temporal). The remaining components are dropped from
consideration.
The head-motion estimates calculated in the correction step were also
placed withi


Report `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/logs/CITATION.html`:

"fMRIPrep citation boilerplate" 
 
 html {
 color: #1a1a1a;
 background-color: #fdfdfd;
 }
 body {
 margin: 0 auto;
 max-width: 36em;
 padding-left: 50px;
 padding-right: 50px;
 padding-top: 50px;
 padding-bottom: 50px;
 hyphens: auto;
 overflow-wrap: break-word;
 text-rendering: optimizeLegibility;
 font-kerning: normal;
 }
 @media (max-width: 600px) {
 body {
 font-size: 0.9em;
 padding: 12px;
 }
 h1 {
 font-size: 1.8em;
 }
 }
 @media print {
 html {
 background-color: white;
 }
 body {
 background-color: transparent;
 color: black;
 font-size: 12pt;
 }
 p, h2, h3 {
 orphans: 3;
 widows: 3;
 }
 h2, h3, h4 {
 page-break-after: avoid;
 }
 }
 p {
 margin: 1em 0;
 }
 a {
 color: #1a1a1a;
 }
 a:visited {
 color: #1a1a1a;
 }
 img {
 max-width: 100%;
 }
 svg {
 height: auto;
 max-width: 100%;
 }
 h1, h2, h3, h4, h5, h6 {
 margin-top: 1.4em;
 }
 h5, h6 {
 font-size: 1em;
 font-style: italic;
 }
 h6 {
 font-weight: normal;
 }
 ol, ul {
 padding-left: 1.7em;
 margin-top: 1em;
 }
 li > ol, li > ul {
 margin-top: 0;
 }
 blockquote {
 margin: 1em 0 1em 1.7em;
 padding-left: 1em;
 border-left: 2px solid #e6e6e6;
 color: #606060;
 }
 code {
 font-family: Menlo, Monaco, Consolas, 'Lucida Console', monospace;
 font-size: 85%;
 margin: 0;
 hyphens: manual;
 }
 pre {
 margin: 1em 0;
 overflow: auto;
 }
 pre code {
 padding: 0;
 overflow: visible;
 overflow-wrap: normal;
 }
 .sourceCode {
 background-color: transparent;
 overflow: visible;
 }
 hr {
 background-color: #1a1a1a;
 border: none;
 height: 1px;
 margin: 1em 0;
 }
 table {
 margin: 1em 0;
 border-collapse: collapse;
 width: 100%;
 overflow-x: auto;
 display: block;
 font-variant-numeric: lining-nums tabular-nums;
 }
 table caption {
 margin-bottom: 0.75em;
 }
 tbody {
 margin-top: 0.5em;
 border-top: 1px solid #1a1a1a;
 border-bottom: 1px solid #1a1a1a;
 }
 th {
 border-top: 1px solid #1a1a1a;
 padding: 0.25em 0.5em 0.25em 0.5em;
 }
 td {
 padding: 0.125em 0.5em 0.25em 0.5em;
 }
 header {
 margin-bottom: 4em;
 text-align: center;
 }
 #TOC li {
 list-style: none;
 }
 #TOC ul {
 padding-left: 1.3em;
 }
 #TOC > ul {
 padding-left: 0;
 }
 #TOC a:not(:hover) {
 text-decoration: none;
 }
 code{white-space: pre-wrap;}
 span.smallcaps{font-variant: small-caps;}
 div.columns{display: flex; gap: min(4vw, 1.5em);}
 div.column{flex: auto; overflow-x: auto;}
 div.hanging-indent{margin-left: 1.5em; text-indent: -1.5em;}
 /* The extra [class] is a hack that increases specificity enough to
 override a similar rule in reveal.js */
 ul.task-list[class]{list-style: none;}
 ul.task-list li input[type="checkbox"] {
 font-size: inherit;
 width: 0.8em;
 margin: 0 0.8em 0.2em -1.6em;
 vertical-align: middle;
 }
 .display.math{display: block; text-align: center; margin: 0.5rem auto;}
 /* CSS for citations */
 div.csl-bib-body { }
 div.csl-entry {
 clear: both;
 margin-bottom: 0em;
 }
 .hanging-indent div.csl-entry {
 margin-left:2em;
 text-indent:-2em;
 }
 div.csl-left-margin {
 min-width:2em;
 float:left;
 }
 div.csl-right-inline {
 margin-left:2em;
 padding-left:1em;
 }
 div.csl-indent {
 margin-left: 2em;
 } 
 
 
 Results included in this manuscript come from preprocessing performed
using fMRIPrep 24.1.1 ( Esteban et al. (2019) ; Esteban et al. (2018) ;
RRID:SCR_016216), which is based on Nipype 1.8.6 ( K. Gorgolewski et al.
(2011) ; K. J.
Gorgolewski et al. (2018) ; RRID:SCR_002502). 
 
 Preprocessing of B 0 inhomogeneity mappings 
 
 A total of 2 fieldmaps were found available within the input BIDS
structure for this particular subject. A B 0 
nonuniformity map (or fieldmap ) was directly measured with an
MRI scheme designed with that purpose such as SEI (Spiral-Echo
Imaging). 
 
 Anatomical data preprocessing 
 
 A total of 1 T1-weighted (T1w) images were found within the input
BIDS dataset. The T1w image was corrected for intensity non-uniformity
(INU) with N4BiasFieldCorrection (Tustison et al. 2010) , distributed with ANTs
2.5.3 (Avants et al. 2008,
RRID:SCR_004757) , and used as T1w-reference throughout the
workflow. The T1w-reference was then skull-stripped with a
 Nipype implementation of the
 antsBrainExtraction.sh workflow (from ANTs), using
OASIS30ANTs as target template. Brain tissue segmentation of
cerebrospinal fluid (CSF), white-matter (WM) and gray-matter (GM) was
performed on the brain-extracted T1w using fast (FSL (version unknown),
RRID:SCR_002823, Zhang, Brady, and Smith 2001) . Volume-based
spatial normalization to two standard spaces (MNI152NLin6Asym,
MNI152NLin2009cAsym) was performed through nonlinear registration with
 antsRegistration (ANTs 2.5.3), using brain-extracted
versions of both T1w reference and the T1w template. The following
templates were were selected for spatial normalization and accessed with
 TemplateFlow (24.2.0, Ciric et al. 2022) : FSL’s
MNI ICBM 152 non-linear 6th Generation Asymmetric Average Brain
Stereotaxic Registration Model [ Evans et al. (2012) ,
RRID:SCR_002823; TemplateFlow ID: MNI152NLin6Asym], ICBM 152
Nonlinear Asymmetrical templ


Report `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/logs/CITATION.bib`:

@article{fmriprep1,
 author = {Esteban, Oscar and Markiewicz, Christopher and Blair, Ross W and Moodie, Craig and Isik, Ayse Ilkay and Erramuzpe Aliaga, Asier and Kent, James and Goncalves, Mathias and DuPre, Elizabeth and Snyder, Madeleine and Oya, Hiroyuki and Ghosh, Satrajit and Wright, Jessey and Durnez, Joke and Poldrack, Russell and Gorgolewski, Krzysztof Jacek},
 title = {{fMRIPrep}: a robust preprocessing pipeline for functional {MRI}},
 volume = {16},
 pages = {111--116},
 year = {2019},
 doi = {10.1038/s41592-018-0235-4},
 journal = {Nature Methods}
}

@article{fmriprep2,
 author = {Esteban, Oscar and Blair, Ross and Markiewicz, Christopher J. and Berleant, Shoshana L. and Moodie, Craig and Ma, Feilong and Isik, Ayse Ilkay and Erramuzpe, Asier and Kent, James D. andGoncalves, Mathias and DuPre, Elizabeth and Sitek, Kevin R. and Gomez, Daniel E. P. and Lurie, Daniel J. and Ye, Zhifang and Poldrack, Russell A. and Gorgolewski, Krzysztof J.},
 title = {fMRIPrep 24.1.1},
 year = 2018,
 doi = {10.5281/zenodo.852659},
 publisher = {Zenodo},
 journal = {Software}
}

@article{templateflow,
 author = {Ciric, R. and Thompson, William H. and Lorenz, R. and Goncalves, M. and MacNicol, E. and Markiewicz, C. J. and Halchenko, Y. O. and Ghosh, S. S. and Gorgolewski, K. J. and Poldrack, R. A. and Esteban, O.},
 title = {{TemplateFlow}: {FAIR}-sharing of multi-scale, multi-species brain models},
 volume = {19},
 pages = {1568--1571},
 year = {2022},
 doi = {10.1038/s41592-022-01681-2},
 journal = {Nature Methods}
}

@article{nipype1,
 author = {Gorgolewski, K. and Burns, C. D. and Madison, C. and Clark, D. and Halchenko, Y. O. and Waskom, M. L. and Ghosh, S.},
 doi = {10.3389/fninf.2011.00013},
 journal = {Frontiers in Neuroinformatics},
 pages = 13,
 shorttitle = {Nipype},
 title = {Nipype: a flexible, lightweight and extensible neuroimaging data processing framework in Python},
 volume = 5,
 year = 2011
}

@article{nipype2,
 author = {Gorgolewski, Krzysztof J. and Esteban, Oscar and Markiewicz, Christopher J. and Ziegler, Erik and Ellis, David Gage and Notter, Michael Philipp and Jarecka, Dorota and Johnson, Hans and Burns, Christopher and Manhães-Savio, Alexandre and Hamalainen, Carlo and Yvernault, Benjamin and Salo, Taylor and Jordan, Kesshi and Goncalves, Mathias and Waskom, Michael and Clark, Daniel and Wong, Jason and Loney, Fred and Modat, Marc and Dewey, Blake E and Madison, Cindee and Visconti di Oleggio Castello, Matteo and Clark, Michael G. and Dayan, Michael and Clark, Dav and Keshavan, Anisha and Pinsard, Basile and Gramfort, Alexandre and Berleant, Shoshana and Nielson, Dylan M. and Bougacha, Salma and Varoquaux, Gael and Cipollini, Ben and Markello, Ross and Rokem, Ariel and Moloney, Brendan and Halchenko, Yaroslav O. and Wassermann , Demian and Hanke, Michael and Horea, Christian and Kaczmarzyk, Jakub and Gilles de Hollander and DuPre, Elizabeth and Gillman, Ashley and Mordom, David and Buchanan, Colin and Tungaraza, Rosalia and Pauli, Wolfgang M. and Iqbal, Shariq and Sikka, Sharad and Mancini, Matteo and Schwartz, Yannick and Malone, Ian B. and Dubois, Mathieu and Frohlich, Caroline and Welch, David and Forbes, Jessica and Kent, James and Watanabe, Aimi and Cumba, Chad and Huntenburg, Julia M. and Kastman, Erik and Nichols, B. Nolan and Eshaghi, Arman and Ginsburg, Daniel and Schaefer, Alexander and Acland, Benjamin and Giavasis, Steven and Kleesiek, Jens and Erickson, Drew and Küttner, René and Haselgrove, Christian and Correa, Carlos and Ghayoor, Ali and Liem, Franz and Millman, Jarrod and Haehn, Daniel and Lai, Jeff and Zhou, Dale and Blair, Ross and Glatard, Tristan and Renfro, Mandy and Liu, Siqi and Kahn, Ari E. and Pérez-García, Fernando and Triplett, William and Lampe, Leonie and Stadler, Jörg and Kong, Xiang-Zhen and Hallquist, Michael and Chetverikov, Andrey and Salvatore, John and Park, Anne and Poldrack, Russell and Craddock, R. Cameron and Inati, Souheil and Hinds, Oliver and Cooper, Gavin and Perkins, L. Nathan and Marina, Ana and Mattfeld, Aaron and Noel, Maxime and Lukas Snoek and Matsubara, K and Cheung, Brian and Rothmei, Simon and Urchs, Sebastian and Durnez, Joke and Mertz, Fred and Geisler, Daniel and Floren, Andrew and Gerhard, Stephan and Sharp, Paul and Molina-Romero, Miguel and Weinstein, Alejandro and Broderick, William and Saase, Victor and Andberg, Sami Kristian and Harms, Robbert and Schlamp, Kai and Arias, Jaime and Papadopoulos Orfanos, Dimitri and Tarbert, Claire and Tambini, Arielle and De La Vega, Alejandro and Nickson, Thomas and Brett, Matthew and Falkiewicz, Marcel and Podranski, Kornelius and Linkersdörfer, Janosch and Flandin, Guillaume and Ort, Eduard and Shachnev, Dmitry and McNamee, Daniel and Davison, Andrew and Varada, Jan and Schwabacher, Isaac and Pellman, John and Perez-Guevara, Martin and Khanuja, Ranjeet and Pannetier, Nicolas and McDermottroe, Conor and Ghosh, Satrajit},
 title = {Nipype},
 year = 2018,
 doi = {10.5281/zenodo.596855},
 publisher


No historical/rebuilt candidate config pair is available. Exact historical execution differences remain unknown.

## Frozen rebuild inputs

These comparisons detect changes since the rebuild snapshot, not changes since the original historical analysis. Image size/mtime matches are not content-hash proof.

### sub-10369_ses-01_task-trust_run-1

12 snapshot checks; 0 missing/changed.

Historical reference vs rebuilt outputs (same midpoint index; historical FEAT mask):

- `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/sub-10369/ses-01/func/sub-10369_ses-01_task-trust_run-1_echo-1_part-mag_desc-preproc_bold.nii.gz`: not compared (different grid or incompatible dimensions).
- `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/sub-10369/ses-01/func/sub-10369_ses-01_task-trust_run-1_echo-2_part-mag_desc-preproc_bold.nii.gz`: not compared (different grid or incompatible dimensions).
- `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/sub-10369/ses-01/func/sub-10369_ses-01_task-trust_run-1_echo-3_part-mag_desc-preproc_bold.nii.gz`: not compared (different grid or incompatible dimensions).
- `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/sub-10369/ses-01/func/sub-10369_ses-01_task-trust_run-1_echo-4_part-mag_desc-preproc_bold.nii.gz`: not compared (different grid or incompatible dimensions).
- `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/sub-10369/ses-01/func/sub-10369_ses-01_task-trust_run-1_part-mag_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz`: r=0.9966056362779995; relative RMS=0.043241450457134636; RMS after global scale=0.043237527636462304
- `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/sub-10369/ses-01/func/sub-10369_ses-01_task-trust_run-1_part-mag_space-T1w_desc-preproc_bold.nii.gz`: not compared (different grid or incompatible dimensions).
### sub-10402_ses-01_task-trust_run-1

12 snapshot checks; 0 missing/changed.

Historical reference vs rebuilt outputs (same midpoint index; historical FEAT mask):

- `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/sub-10402/ses-01/func/sub-10402_ses-01_task-trust_run-1_echo-1_part-mag_desc-preproc_bold.nii.gz`: not compared (different grid or incompatible dimensions).
- `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/sub-10402/ses-01/func/sub-10402_ses-01_task-trust_run-1_echo-2_part-mag_desc-preproc_bold.nii.gz`: not compared (different grid or incompatible dimensions).
- `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/sub-10402/ses-01/func/sub-10402_ses-01_task-trust_run-1_echo-3_part-mag_desc-preproc_bold.nii.gz`: not compared (different grid or incompatible dimensions).
- `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/sub-10402/ses-01/func/sub-10402_ses-01_task-trust_run-1_echo-4_part-mag_desc-preproc_bold.nii.gz`: not compared (different grid or incompatible dimensions).
- `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/sub-10402/ses-01/func/sub-10402_ses-01_task-trust_run-1_part-mag_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz`: r=0.9951317391289521; relative RMS=0.05201085919832398; RMS after global scale=0.05163963250697017
- `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/sub-10402/ses-01/func/sub-10402_ses-01_task-trust_run-1_part-mag_space-T1w_desc-preproc_bold.nii.gz`: not compared (different grid or incompatible dimensions).
### sub-10402_ses-01_task-trust_run-2

12 snapshot checks; 0 missing/changed.

Historical reference vs rebuilt outputs (same midpoint index; historical FEAT mask):

- `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/sub-10402/ses-01/func/sub-10402_ses-01_task-trust_run-2_echo-1_part-mag_desc-preproc_bold.nii.gz`: not compared (different grid or incompatible dimensions).
- `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/sub-10402/ses-01/func/sub-10402_ses-01_task-trust_run-2_echo-2_part-mag_desc-preproc_bold.nii.gz`: not compared (different grid or incompatible dimensions).
- `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/sub-10402/ses-01/func/sub-10402_ses-01_task-trust_run-2_echo-3_part-mag_desc-preproc_bold.nii.gz`: not compared (different grid or incompatible dimensions).
- `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/sub-10402/ses-01/func/sub-10402_ses-01_task-trust_run-2_echo-4_part-mag_desc-preproc_bold.nii.gz`: not compared (different grid or incompatible dimensions).
- `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/sub-10402/ses-01/func/sub-10402_ses-01_task-trust_run-2_part-mag_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz`: r=0.9877619261469892; relative RMS=0.0816934430238724; RMS after global scale=0.08122407075329957
- `/ZPOOL/data/projects/rf1-betrayal/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41/fmriprep-24/sub-10402/ses-01/func/sub-10402_ses-01_task-trust_run-2_part-mag_space-T1w_desc-preproc_bold.nii.gz`: not compared (different grid or incompatible dimensions).

## Errors and limits

0 collection/read errors (see evidence.json).

No config difference is automatically attributed as the cause. Check task/session coverage and timestamps before choosing a controlled rerun. Matching a midpoint frame alone cannot prove full-series equivalence.

Full settings, report text, source identities and snapshot details are in evidence.json and sources/. No images, full nuisance matrices, or participant tables are packaged.
