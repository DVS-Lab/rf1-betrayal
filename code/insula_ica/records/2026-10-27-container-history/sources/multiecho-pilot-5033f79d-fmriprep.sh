# example code for FMRIPREP
# runs FMRIPREP on input subject
# usage: bash run_fmriprep.sh sub
# example: bash run_fmriprep.sh 102

sub=$1

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"
bidsdir=$maindir/bids

# make derivatives folder if it doesn't exist.
# let's keep this out of bids for now
if [ ! -d $maindir/derivatives/fmriprep ]; then
	mkdir -p $maindir/derivatives/fmriprep
fi

scratchdir=/ZPOOL/data/scratch/`whoami`
if [ ! -d $scratchdir ]; then
	mkdir -p $scratchdir
fi



TEMPLATEFLOW_DIR=/ZPOOL/data/tools/templateflow
export APPTAINERENV_TEMPLATEFLOW_HOME=/opt/templateflow
singularity run --cleanenv \
-B ${TEMPLATEFLOW_DIR}:/opt/templateflow \
-B $maindir:/base \
-B /ZPOOL/data/tools/licenses:/opts \
-B $scratchdir:/scratch \
-B /mnt \
/ZPOOL/data/tools/fmriprep-24.1.1.simg \
/base/bids /base/derivatives/fmriprep \
participant --participant_label $sub \
--stop-on-first-crash \
--me-output-echos \
--output-spaces MNI152NLin6Asym \
--bids-filter-file /base/code/fmriprep_config.json \
--fs-no-reconall --fs-license-file /opts/fs_license.txt -w /scratch


# NOTE: the --use-syn-sdc option will only be applied in cases where there is not a valid fmap.
# https://neurostars.org/t/using-use-syn-sdc-with-fieldmap-data/2592/2

