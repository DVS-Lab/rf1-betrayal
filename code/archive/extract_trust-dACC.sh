# directory (borrowed style from L1stats.sh)
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
basedir="$(dirname "$scriptdir")"

# Define input variables
sublist="/gpfs/scratch/tug87422/smithlab-shared/rf1-derrick/code/sublist_DD128.txt"
task="trust"
model="1"
sm="5"
mask="/gpfs/scratch/tug87422/smithlab-shared/rf1-norms/masks/thr_seed-dACC.nii.gz"

# Output directory for the meants text files
outputdir="/gpfs/scratch/tug87422/smithlab-shared/rf1-derrick/derivatives/dACC_trust_meants"
log="/gpfs/scratch/tug87422/smithlab-shared/rf1-derrick/code/dacc_trust_meants_log.txt"
mkdir -p $outputdir

# Loop through each subject
for sub in `cat /gpfs/scratch/tug87422/smithlab-shared/rf1-derrick/code/sublist_DD128.txt`; do
#for sub in 10402; do 
  # Define the path to the filtered functional data
 # /gpfs/scratch/tug87422/smithlab-shared/rf1-sra-trust/derivatives/fsl/sub-10402/L2_task-trust_model-1_type-act_sm-5.gfeat/cope10.feat/stats

  DATA="/gpfs/scratch/tug87422/smithlab-shared/r01-soi/derivatives/fsl/final/sub-${sub}/L2_task-${task}_model-${model}_type-act_sm-${sm}.gfeat/cope10.feat/stats/zstat1.nii.gz"
  DATAC4="/gpfs/scratch/tug87422/smithlab-shared/r01-soi/derivatives/fsl/final/sub-${sub}/L2_task-${task}_model-${model}_type-act_sm-${sm}.gfeat/cope4.feat/stats/zstat1.nii.gz"
  DATAC5="/gpfs/scratch/tug87422/smithlab-shared/r01-soi/derivatives/fsl/final/sub-${sub}/L2_task-${task}_model-${model}_type-act_sm-${sm}.gfeat/cope5.feat/stats/zstat1.nii.gz"
  DATAC8="/gpfs/scratch/tug87422/smithlab-shared/r01-soi/derivatives/fsl/final/sub-${sub}/L2_task-${task}_model-${model}_type-act_sm-${sm}.gfeat/cope8.feat/stats/zstat1.nii.gz"
  DATAC9="/gpfs/scratch/tug87422/smithlab-shared/r01-soi/derivatives/fsl/final/sub-${sub}/L2_task-${task}_model-${model}_type-act_sm-${sm}.gfeat/cope9.feat/stats/zstat1.nii.gz"
  DATArun1="/gpfs/scratch/tug87422/smithlab-shared/r01-soi/derivatives/fsl/final/sub-${sub}/L1_task-${task}_model-${model}_type-act_run-1_sm-${sm}.feat/stats/zstat10.nii.gz"
  DATAC4run1="/gpfs/scratch/tug87422/smithlab-shared/r01-soi/derivatives/fsl/final/sub-${sub}/L1_task-${task}_model-${model}_type-act_run-1_sm-${sm}.feat/stats/zstat4.nii.gz"
  DATAC5run1="/gpfs/scratch/tug87422/smithlab-shared/r01-soi/derivatives/fsl/final/sub-${sub}/L1_task-${task}_model-${model}_type-act_run-1_sm-${sm}.feat/stats/zstat5.nii.gz"
  DATAC8run1="/gpfs/scratch/tug87422/smithlab-shared/r01-soi/derivatives/fsl/final/sub-${sub}/L1_task-${task}_model-${model}_type-act_run-1_sm-${sm}.feat/stats/zstat8.nii.gz"
  DATAC9run1="/gpfs/scratch/tug87422/smithlab-shared/r01-soi/derivatives/fsl/final/sub-${sub}/L1_task-${task}_model-${model}_type-act_run-1_sm-${sm}.feat/stats/zstat9.nii.gz"
  DATArun2="/gpfs/scratch/tug87422/smithlab-shared/r01-soi/derivatives/fsl/final/sub-${sub}/L1_task-${task}_model-${model}_type-act_run-2_sm-${sm}.feat/stats/zstat10.nii.gz"
  DATAC4run2="/gpfs/scratch/tug87422/smithlab-shared/r01-soi/derivatives/fsl/final/sub-${sub}/L1_task-${task}_model-${model}_type-act_run-2_sm-${sm}.feat/stats/zstat4.nii.gz"
  DATAC5run2="/gpfs/scratch/tug87422/smithlab-shared/r01-soi/derivatives/fsl/final/sub-${sub}/L1_task-${task}_model-${model}_type-act_run-2_sm-${sm}.feat/stats/zstat5.nii.gz"
  DATAC8run2="/gpfs/scratch/tug87422/smithlab-shared/r01-soi/derivatives/fsl/final/sub-${sub}/L1_task-${task}_model-${model}_type-act_run-2_sm-${sm}.feat/stats/zstat8.nii.gz"
  DATAC9run2="/gpfs/scratch/tug87422/smithlab-shared/r01-soi/derivatives/fsl/final/sub-${sub}/L1_task-${task}_model-${model}_type-act_run-2_sm-${sm}.feat/stats/zstat9.nii.gz"
  DATAC54="/gpfs/scratch/tug87422/smithlab-shared/r01-soi/derivatives/fsl/final/sub-${sub}/L2_task-${task}_model-${model}_type-act_sm-${sm}.gfeat/cope10.feat/stats/zstat54.nii.gz"
  DATAC98="/gpfs/scratch/tug87422/smithlab-shared/r01-soi/derivatives/fsl/final/sub-${sub}/L2_task-${task}_model-${model}_type-act_sm-${sm}.gfeat/cope10.feat/stats/zstat98.nii.gz"
  DATAC54run1="/gpfs/scratch/tug87422/smithlab-shared/r01-soi/derivatives/fsl/final/sub-${sub}/L1_task-${task}_model-${model}_type-act_run-1_sm-${sm}.feat/stats/zstat54.nii.gz"
  DATAC98run1="/gpfs/scratch/tug87422/smithlab-shared/r01-soi/derivatives/fsl/final/sub-${sub}/L1_task-${task}_model-${model}_type-act_run-1_sm-${sm}.feat/stats/zstat98.nii.gz"
  DATAC54run2="/gpfs/scratch/tug87422/smithlab-shared/r01-soi/derivatives/fsl/final/sub-${sub}/L1_task-${task}_model-${model}_type-act_run-2_sm-${sm}.feat/stats/zstat54.nii.gz"
  DATAC98run2="/gpfs/scratch/tug87422/smithlab-shared/r01-soi/derivatives/fsl/final/sub-${sub}/L1_task-${task}_model-${model}_type-act_run-2_sm-${sm}.feat/stats/zstat98.nii.gz"

  DATASOC="/gpfs/scratch/tug87422/smithlab-shared/r01-soi/derivatives/fsl/final/sub-${sub}/L2_task-${task}_model-${model}_type-act_sm-${sm}.gfeat/cope10.feat/stats/zstat20.nii.gz"
  DATASOCrun1="/gpfs/scratch/tug87422/smithlab-shared/r01-soi/derivatives/fsl/final/sub-${sub}/L1_task-${task}_model-${model}_type-act_run-1_sm-${sm}.feat/stats/zstat20.nii.gz"
  DATASOCrun2="/gpfs/scratch/tug87422/smithlab-shared/r01-soi/derivatives/fsl/final/sub-${sub}/L1_task-${task}_model-${model}_type-act_run-2_sm-${sm}.feat/stats/zstat20.nii.gz"


  # Definec the output file
  
  # Check if the data file exists
  if [ -f "$DATA" ]; then   
     echo "Extracting L2 dACC signal for subject: ${sub}" >> $log
    OUTPUT="${outputdir}/sub-${sub}_L2_dACC_meants.txt" 
    fslmeants -i $DATA -o $OUTPUT -m $mask
    
    #create DATASOC
         
         fslmaths $DATAC5 -sub $DATAC4 $DATAC54
         fslmaths $DATAC9 -sub $DATAC8 $DATAC98
         fslmaths $DATAC98 -sub $DATAC54 $DATASOC
   	
    chmod 777 $OUTPUT
  else
    echo "sub${sub} ${DATA} doesn't exist" >> $log
  fi
  	
   if [ -f "$DATArun1" ]; then
   	 echo "Extracting run1 dACC signal for subject: ${sub}" >> $log
   	 OUTPUT="${outputdir}/sub-${sub}_run-1_dACC_meants.txt"
   	 fslmeants -i $DATArun1 -o $OUTPUT -m $mask

    #create DATASOCrun1
         fslmaths $DATAC5run1 -sub $DATAC4run1 $DATAC54run1
         fslmaths $DATAC9run1 -sub $DATAC8run1 $DATAC98run1
         fslmaths $DATAC98run1 -sub $DATAC54run1 $DATASOCrun1

     chmod 777 $OUTPUT
   else 
   	 echo "WARNING: Run-1 ${DATArun1} not found for subject ${sub}. Skipping." >> $log
   fi
 	if [ -f "$DATArun2" ]; then
   	 echo "Extracting run2 dACC signal for subject: ${sub}" >> $log
   	 OUTPUT="${outputdir}/sub-${sub}_run-2_dACC_meants.txt"
   	 fslmeants -i $DATArun2 -o $OUTPUT -m $mask

   #create DATASOCrun2
      
         fslmaths $DATAC5run2 -sub $DATAC4run2 $DATAC54run2
         fslmaths $DATAC9run2 -sub $DATAC8run2 $DATAC98run2
         fslmaths $DATAC98run2 -sub $DATAC54run2 $DATASOCrun2

     chmod 777 $OUTPUT
   else 
   	 echo "WARNING: Run-2 ${DATArun2} not found for subject ${sub}. Skipping." >> $log
   fi

  # Check if the social data file exists
  if [ -f "$DATASOC" ]; then
     echo "Extracting L2 dACC signal for subject: ${sub}" >> $log
    OUTPUT="${outputdir}/sub-${sub}_L2_dACC_SOC_meants.txt"
    fslmeants -i $DATASOC -o $OUTPUT -m $mask
    chmod 777 $OUTPUT
  else
    echo "sub${sub} ${DATASOC} doesn't exist" >> $log
  fi

   if [ -f "$DATASOCrun1" ]; then
         echo "Extracting run1 dACC signal for subject: ${sub}" >> $log
         OUTPUT="${outputdir}/sub-${sub}_run-1_dACC_SOC_meants.txt"
         fslmeants -i $DATASOCrun1 -o $OUTPUT -m $mask
     chmod 777 $OUTPUT
   else
         echo "WARNING: Run-1 ${DATASOCrun1} not found for subject ${sub}. Skipping." >> $log
   fi
        if [ -f "$DATASOCrun2" ]; then
         echo "Extracting run2 dACC signal for subject: ${sub}" >> $log
         OUTPUT="${outputdir}/sub-${sub}_run-2_dACC_SOC_meants.txt"
         fslmeants -i $DATASOCrun2 -o $OUTPUT -m $mask
     chmod 777 $OUTPUT
   else
         echo "WARNING: Run-2 ${DATASOCrun2} not found for subject ${sub}. Skipping." >> $log
   fi

done

echo "Extraction complete. All output files saved to: $outputdir"

