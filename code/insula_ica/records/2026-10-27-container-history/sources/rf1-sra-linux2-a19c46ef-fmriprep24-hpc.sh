#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -N fmriprep-rf1
#PBS -q normal
#PBS -l nodes=1:ppn=28

# load modules and go to workdir   (keep same stack you used in the working script)
module load fsl/6.0.2
source $FSLDIR/etc/fslconf/fsl.sh
module load singularity
cd $PBS_O_WORKDIR

# paths for rf1-sra-linux2
projectname=rf1-sra-linux2
maindir=/gpfs/scratch/tug87422/smithlab-shared/$projectname
scriptdir=$maindir/code
bidsdir=$maindir/bids
logdir=$maindir/logs
mkdir -p $logdir

rm -f $logdir/cmd_fmriprep_${PBS_JOBID}.txt
touch  $logdir/cmd_fmriprep_${PBS_JOBID}.txt

# derivatives + scratch
if [ ! -d $maindir/derivatives ]; then
    mkdir -p $maindir/derivatives
fi
scratchdir=~/scratch/$projectname/fmriprep-24
mkdir -p $scratchdir

# container bind
TEMPLATEFLOW_DIR=/gpfs/scratch/tug87422/smithlab-shared/tools/templateflow
MPLCONFIGDIR_DIR=/gpfs/scratch/tug87422/smithlab-shared/tools/mplconfigdir
export SINGULARITYENV_TEMPLATEFLOW_HOME=/opt/templateflow
export SINGULARITYENV_MPLCONFIGDIR=/opt/mplconfigdir

# subjects comes from qsub -v subjects="XXXXX"
for sub in ${subjects[@]}; do
  # loop sessions; only queue if the session dir exists
  for ses in ses-01 ses-02; do
    if [ -d "$bidsdir/sub-${sub}/${ses}" ]; then
      ses_num=${ses#ses-}
      cfg="/base/code/ses-${ses_num}.json"   # <- per-session filter keeps ANAT/FUNC/FMAP to that session

      echo singularity run --cleanenv \
        -B ${TEMPLATEFLOW_DIR}:/opt/templateflow \
        -B ${MPLCONFIGDIR_DIR}:/opt/mplconfigdir \
        -B $maindir:/base \
        -B ~/work/tools/hpctools/licenses:/opts \
        -B $scratchdir:/scratch \
        ~/work/tools/hpctools/fmriprep-24.1.1.simg \
        /base/bids /base/derivatives/fmriprep-24 \
        participant --participant_label $sub \
        --stop-on-first-crash \
        --skip-bids-validation \
        --nthreads 14 \
        --me-output-echos \
        --output-spaces T1w MNI152NLin6Asym \
        --bids-filter-file $cfg \
        --fs-no-reconall --fs-license-file /opts/fs_license.txt \
        -w /scratch/sub-${sub}_${ses} >> $logdir/cmd_fmriprep_${PBS_JOBID}.txt
    fi
  done
done

torque-launch -p $logdir/chk_fmriprep_${PBS_JOBID}.txt $logdir/cmd_fmriprep_${PBS_JOBID}.txt

