# directory (borrowed style from L1stats.sh)
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
basedir="$(dirname "$scriptdir")"

# Define input variables
sublist=${basedir}/code/sublist_n132.txt
task="ugr"
model="3"
sm="5"
mask=${basedir}/masks/vmPFC.nii.gz

# Output directory for the meants text files
outputdir=${basedir}/derivatives/extractions/vmPFC_ugr_ppi_meants
log=${basedir}/code/vmPFC_ugr_meants_log.txt
mkdir -p $outputdir

# Loop through each subject
for sub in `cat ${basedir}/code/sublist_n132.txt`; do

  # Define the path to the filtered functional data
  DATA=${basedir}/derivatives/fsl/sub-${sub}/ses-01/L2_task-${task}_model-${model}_type-ppi_seed-AIns_sm-${sm}.gfeat/cope12.feat/stats/zstat1.nii.gz

  DATArun1=${basedir}/derivatives/fsl/sub-${sub}/ses-01/L1_task-${task}_model-${model}_type-ppi_seed-AIns_run-1_sm-${sm}.feat/stats/zstat12.nii.gz
  
  DATArun2=${basedir}/derivatives/fsl/sub-${sub}/ses-01/L1_task-${task}_model-${model}_type-ppi_seed-AIns_run-2_sm-${sm}.feat/stats/zstat12.nii.gz
 

  # Check if the data file exists
  if [ -f "$DATA" ]; then   
     echo "Extracting L2 vmPFC signal for subject: ${sub}" >> $log
    OUTPUT="${outputdir}/sub-${sub}_L2_vmPFC_meants.txt" 
    fslmeants -i $DATA -o $OUTPUT -m $mask
   	
    chmod 777 $OUTPUT
  else
    echo "sub${sub} ${DATA} doesn't exist" >> $log
  fi
  	
   if [ -f "$DATArun1" ]; then
   	 echo "Extracting run1 vmPFC signal for subject: ${sub}" >> $log
   	 OUTPUT="${outputdir}/sub-${sub}_run-1_vmPFC_meants.txt"
   	 fslmeants -i $DATArun1 -o $OUTPUT -m $mask


     chmod 777 $OUTPUT
   else 
   	 echo "WARNING: Run-1 ${DATArun1} not found for subject ${sub}. Skipping." >> $log
   fi
 	if [ -f "$DATArun2" ]; then
   	 echo "Extracting run2 vmPFC signal for subject: ${sub}" >> $log
   	 OUTPUT="${outputdir}/sub-${sub}_run-2_vmPFC_meants.txt"
   	 fslmeants -i $DATArun2 -o $OUTPUT -m $mask
 

     chmod 777 $OUTPUT
   else 
   	 echo "WARNING: Run-2 ${DATArun2} not found for subject ${sub}. Skipping." >> $log
   fi


done

echo "Extraction complete. All output files saved to: $outputdir"

